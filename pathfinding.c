/**
 * @file pathfinding.c
 * @brief Pathfinding and Waypoint Generation for UAV Navigation
 *
 * Implements A* and RRT pathfinding algorithms for UAV navigation.
 * Optimized for ARM Cortex-A76 (Raspberry Pi 5B).
 */

#include "pathfinding.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/* ============================================================================
 * Binary Heap Priority Queue for A*
 * ============================================================================ */

typedef struct {
    astar_node_t** heap;
    int size;
    int capacity;
} priority_queue_t;

static priority_queue_t* pq_create(int capacity) {
    priority_queue_t* pq = (priority_queue_t*)malloc(sizeof(priority_queue_t));
    if (!pq) return NULL;

    pq->heap = (astar_node_t**)malloc(sizeof(astar_node_t*) * capacity);
    if (!pq->heap) {
        free(pq);
        return NULL;
    }

    pq->size = 0;
    pq->capacity = capacity;
    return pq;
}

static void pq_free(priority_queue_t* pq) {
    if (!pq) return;
    free(pq->heap);
    free(pq);
}

static void pq_swap(priority_queue_t* pq, int i, int j) {
    astar_node_t* temp = pq->heap[i];
    pq->heap[i] = pq->heap[j];
    pq->heap[j] = temp;
}

static void pq_heapify_up(priority_queue_t* pq, int idx) {
    while (idx > 0) {
        int parent = (idx - 1) / 2;
        if (pq->heap[idx]->f_score >= pq->heap[parent]->f_score) break;
        pq_swap(pq, idx, parent);
        idx = parent;
    }
}

static void pq_heapify_down(priority_queue_t* pq, int idx) {
    while (1) {
        int left = 2 * idx + 1;
        int right = 2 * idx + 2;
        int smallest = idx;

        if (left < pq->size && pq->heap[left]->f_score < pq->heap[smallest]->f_score)
            smallest = left;
        if (right < pq->size && pq->heap[right]->f_score < pq->heap[smallest]->f_score)
            smallest = right;

        if (smallest == idx) break;

        pq_swap(pq, idx, smallest);
        idx = smallest;
    }
}

static int pq_push(priority_queue_t* pq, astar_node_t* node) {
    if (pq->size >= pq->capacity) return -1;

    pq->heap[pq->size] = node;
    pq_heapify_up(pq, pq->size);
    pq->size++;
    return 0;
}

static astar_node_t* pq_pop(priority_queue_t* pq) {
    if (pq->size == 0) return NULL;

    astar_node_t* result = pq->heap[0];
    pq->size--;

    if (pq->size > 0) {
        pq->heap[0] = pq->heap[pq->size];
        pq_heapify_down(pq, 0);
    }

    return result;
}

/* ============================================================================
 * Geometric Utilities
 * ============================================================================ */

float distance_3d(const position3d_t* p1, const position3d_t* p2) {
    float dx = p2->x - p1->x;
    float dy = p2->y - p1->y;
    float dz = p2->z - p1->z;
    return sqrtf(dx*dx + dy*dy + dz*dz);
}

float calculate_heading(const position3d_t* from, const position3d_t* to) {
    float dx = to->x - from->x;
    float dy = to->y - from->y;
    float heading = atan2f(dy, dx) * 180.0f / M_PI;

    // Convert to compass bearing (0=North, 90=East)
    heading = 90.0f - heading;
    if (heading < 0.0f) heading += 360.0f;
    if (heading >= 360.0f) heading -= 360.0f;

    return heading;
}

/* ============================================================================
 * Obstacle Utilities
 * ============================================================================ */

bool obstacle_collides_point(const obstacle_t* obstacle, const position3d_t* point) {
    if (!obstacle->active) return false;

    float dist = distance_3d(&obstacle->position, point);
    return dist < obstacle->radius;
}

bool obstacle_collides_line(const obstacle_t* obstacle,
                            const position3d_t* start,
                            const position3d_t* end) {
    if (!obstacle->active) return false;

    // Line segment from start to end
    float dx = end->x - start->x;
    float dy = end->y - start->y;
    float dz = end->z - start->z;
    float line_len = sqrtf(dx*dx + dy*dy + dz*dz);

    if (line_len < 1e-6f) {
        return obstacle_collides_point(obstacle, start);
    }

    // Unit direction
    float ux = dx / line_len;
    float uy = dy / line_len;
    float uz = dz / line_len;

    // Project obstacle center onto line
    float sx = obstacle->position.x - start->x;
    float sy = obstacle->position.y - start->y;
    float sz = obstacle->position.z - start->z;

    float projection = sx*ux + sy*uy + sz*uz;
    projection = fmaxf(0.0f, fminf(projection, line_len));

    // Closest point on line
    position3d_t closest = {
        .x = start->x + projection * ux,
        .y = start->y + projection * uy,
        .z = start->z + projection * uz
    };

    float dist = distance_3d(&closest, &obstacle->position);
    return dist < obstacle->radius;
}

/* ============================================================================
 * A* Planner Implementation
 * ============================================================================ */

int astar_planner_init(astar_planner_t* planner,
                       float grid_resolution,
                       float x_min, float x_max,
                       float y_min, float y_max) {
    if (!planner) return -1;

    planner->grid_resolution = grid_resolution;
    planner->x_min = x_min;
    planner->x_max = x_max;
    planner->y_min = y_min;
    planner->y_max = y_max;

    planner->nx = (int)((x_max - x_min) / grid_resolution) + 1;
    planner->ny = (int)((y_max - y_min) / grid_resolution) + 1;

    if (planner->nx > PF_MAX_GRID_DIM || planner->ny > PF_MAX_GRID_DIM) {
        return -1;
    }

    // Allocate obstacle map
    planner->obstacle_map = (bool*)calloc(planner->nx * planner->ny, sizeof(bool));
    if (!planner->obstacle_map) return -1;

    planner->num_obstacles = 0;
    planner->signal_field = NULL;

    return 0;
}

void astar_planner_free(astar_planner_t* planner) {
    if (!planner) return;

    if (planner->obstacle_map) {
        free(planner->obstacle_map);
        planner->obstacle_map = NULL;
    }
}

int astar_add_obstacle(astar_planner_t* planner, const obstacle_t* obstacle) {
    if (!planner || !obstacle) return -1;
    if (planner->num_obstacles >= PF_MAX_OBSTACLES) return -1;

    planner->obstacles[planner->num_obstacles] = *obstacle;
    planner->num_obstacles++;

    // Update obstacle map
    for (int iy = 0; iy < planner->ny; iy++) {
        for (int ix = 0; ix < planner->nx; ix++) {
            float x = planner->x_min + ix * planner->grid_resolution;
            float y = planner->y_min + iy * planner->grid_resolution;
            position3d_t point = {x, y, 0.0f};

            if (obstacle_collides_point(obstacle, &point)) {
                planner->obstacle_map[iy * planner->nx + ix] = true;
            }
        }
    }

    return 0;
}

void astar_set_signal_field(astar_planner_t* planner, void* signal_field) {
    if (!planner) return;
    planner->signal_field = signal_field;
}

static void grid_to_world(const astar_planner_t* planner, int grid_x, int grid_y,
                         float* world_x, float* world_y) {
    *world_x = planner->x_min + grid_x * planner->grid_resolution;
    *world_y = planner->y_min + grid_y * planner->grid_resolution;
}

static void world_to_grid(const astar_planner_t* planner, float world_x, float world_y,
                         int* grid_x, int* grid_y) {
    *grid_x = (int)((world_x - planner->x_min) / planner->grid_resolution);
    *grid_y = (int)((world_y - planner->y_min) / planner->grid_resolution);
}

static float heuristic(const astar_planner_t* planner, int x1, int y1, int x2, int y2) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    return sqrtf((float)(dx*dx + dy*dy)) * planner->grid_resolution;
}

static bool is_valid_grid(const astar_planner_t* planner, int x, int y) {
    return x >= 0 && x < planner->nx && y >= 0 && y < planner->ny;
}

static bool is_line_clear(const astar_planner_t* planner,
                         const position3d_t* start, const position3d_t* end) {
    for (int i = 0; i < planner->num_obstacles; i++) {
        if (obstacle_collides_line(&planner->obstacles[i], start, end)) {
            return false;
        }
    }
    return true;
}

int astar_plan(astar_planner_t* planner,
               const position3d_t* start,
               const position3d_t* goal,
               float weight_distance,
               float weight_signal,
               position3d_t* path,
               int* num_waypoints) {
    if (!planner || !start || !goal || !path || !num_waypoints) return -1;

    // Convert to grid coordinates
    int start_gx, start_gy, goal_gx, goal_gy;
    world_to_grid(planner, start->x, start->y, &start_gx, &start_gy);
    world_to_grid(planner, goal->x, goal->y, &goal_gx, &goal_gy);

    // Validate
    if (!is_valid_grid(planner, start_gx, start_gy) ||
        !is_valid_grid(planner, goal_gx, goal_gy)) {
        return -1;
    }

    int start_idx = start_gy * planner->nx + start_gx;
    int goal_idx = goal_gy * planner->nx + goal_gx;

    if (planner->obstacle_map[start_idx] || planner->obstacle_map[goal_idx]) {
        return -1;
    }

    // Create priority queue
    priority_queue_t* open_set = pq_create(10000);
    if (!open_set) return -1;

    // G-scores (cost from start)
    float* g_scores = (float*)malloc(sizeof(float) * planner->nx * planner->ny);
    bool* closed_set = (bool*)calloc(planner->nx * planner->ny, sizeof(bool));
    astar_node_t** node_map = (astar_node_t**)calloc(planner->nx * planner->ny, sizeof(astar_node_t*));

    if (!g_scores || !closed_set || !node_map) {
        pq_free(open_set);
        free(g_scores);
        free(closed_set);
        free(node_map);
        return -1;
    }

    // Initialize g_scores
    for (int i = 0; i < planner->nx * planner->ny; i++) {
        g_scores[i] = FLT_MAX;
    }

    // Start node
    astar_node_t* start_node = (astar_node_t*)malloc(sizeof(astar_node_t));
    start_node->grid_x = start_gx;
    start_node->grid_y = start_gy;
    start_node->g_score = 0.0f;
    start_node->f_score = heuristic(planner, start_gx, start_gy, goal_gx, goal_gy);
    start_node->parent = NULL;

    g_scores[start_idx] = 0.0f;
    node_map[start_idx] = start_node;
    pq_push(open_set, start_node);

    astar_node_t* goal_node = NULL;
    int nodes_expanded = 0;

    // A* search
    while (open_set->size > 0) {
        astar_node_t* current = pq_pop(open_set);
        nodes_expanded++;

        int current_idx = current->grid_y * planner->nx + current->grid_x;

        if (closed_set[current_idx]) {
            continue;
        }

        // Goal reached
        if (current->grid_x == goal_gx && current->grid_y == goal_gy) {
            goal_node = current;
            break;
        }

        closed_set[current_idx] = true;

        // Explore 8-connected neighbors
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;

                int nx = current->grid_x + dx;
                int ny = current->grid_y + dy;

                if (!is_valid_grid(planner, nx, ny)) continue;

                int neighbor_idx = ny * planner->nx + nx;

                if (closed_set[neighbor_idx]) continue;
                if (planner->obstacle_map[neighbor_idx]) continue;

                // Movement cost
                float move_dist = sqrtf((float)(dx*dx + dy*dy)) * planner->grid_resolution;
                float tentative_g = current->g_score + weight_distance * move_dist;

                if (tentative_g >= g_scores[neighbor_idx]) continue;

                // Create or update neighbor node
                astar_node_t* neighbor = node_map[neighbor_idx];
                if (!neighbor) {
                    neighbor = (astar_node_t*)malloc(sizeof(astar_node_t));
                    neighbor->grid_x = nx;
                    neighbor->grid_y = ny;
                    node_map[neighbor_idx] = neighbor;
                }

                neighbor->g_score = tentative_g;
                neighbor->f_score = tentative_g + heuristic(planner, nx, ny, goal_gx, goal_gy);
                neighbor->parent = current;

                g_scores[neighbor_idx] = tentative_g;
                pq_push(open_set, neighbor);
            }
        }
    }

    int result = -1;

    if (goal_node) {
        // Reconstruct path
        int path_len = 0;
        astar_node_t* node = goal_node;

        while (node) {
            path_len++;
            node = node->parent;
        }

        if (path_len <= PF_MAX_WAYPOINTS) {
            node = goal_node;
            int idx = path_len - 1;

            while (node) {
                float wx, wy;
                grid_to_world(planner, node->grid_x, node->grid_y, &wx, &wy);
                path[idx].x = wx;
                path[idx].y = wy;
                path[idx].z = start->z;
                idx--;
                node = node->parent;
            }

            *num_waypoints = path_len;
            result = 0;
        }
    }

    // Cleanup
    for (int i = 0; i < planner->nx * planner->ny; i++) {
        if (node_map[i]) free(node_map[i]);
    }

    pq_free(open_set);
    free(g_scores);
    free(closed_set);
    free(node_map);

    return result;
}

/* ============================================================================
 * RRT Planner Implementation
 * ============================================================================ */

int rrt_planner_init(rrt_planner_t* planner,
                     float x_min, float x_max,
                     float y_min, float y_max,
                     float step_size,
                     int max_iterations) {
    if (!planner) return -1;

    planner->x_min = x_min;
    planner->x_max = x_max;
    planner->y_min = y_min;
    planner->y_max = y_max;
    planner->step_size = step_size;
    planner->max_iterations = max_iterations;
    planner->num_obstacles = 0;
    planner->num_nodes = 0;

    return 0;
}

void rrt_planner_free(rrt_planner_t* planner) {
    if (!planner) return;
    planner->num_nodes = 0;
}

int rrt_add_obstacle(rrt_planner_t* planner, const obstacle_t* obstacle) {
    if (!planner || !obstacle) return -1;
    if (planner->num_obstacles >= PF_MAX_OBSTACLES) return -1;

    planner->obstacles[planner->num_obstacles] = *obstacle;
    planner->num_obstacles++;
    return 0;
}

static void random_sample(const rrt_planner_t* planner, position3d_t* sample) {
    sample->x = planner->x_min + (float)rand() / RAND_MAX * (planner->x_max - planner->x_min);
    sample->y = planner->y_min + (float)rand() / RAND_MAX * (planner->y_max - planner->y_min);
    sample->z = 0.0f;
}

static int nearest_node(const rrt_planner_t* planner, const position3d_t* point) {
    int nearest_idx = 0;
    float min_dist = FLT_MAX;

    for (int i = 0; i < planner->num_nodes; i++) {
        float dist = distance_3d(&planner->tree[i].position, point);
        if (dist < min_dist) {
            min_dist = dist;
            nearest_idx = i;
        }
    }

    return nearest_idx;
}

static void steer(const position3d_t* from, const position3d_t* to,
                 float step_size, position3d_t* result) {
    float dx = to->x - from->x;
    float dy = to->y - from->y;
    float dist = sqrtf(dx*dx + dy*dy);

    if (dist < step_size) {
        *result = *to;
    } else {
        float scale = step_size / dist;
        result->x = from->x + dx * scale;
        result->y = from->y + dy * scale;
        result->z = 0.0f;
    }
}

static bool is_collision_free_rrt(const rrt_planner_t* planner,
                                   const position3d_t* start,
                                   const position3d_t* end) {
    for (int i = 0; i < planner->num_obstacles; i++) {
        if (obstacle_collides_line(&planner->obstacles[i], start, end)) {
            return false;
        }
    }
    return true;
}

int rrt_plan(rrt_planner_t* planner,
             const position3d_t* start,
             const position3d_t* goal,
             float goal_tolerance,
             position3d_t* path,
             int* num_waypoints) {
    if (!planner || !start || !goal || !path || !num_waypoints) return -1;

    // Initialize tree with start node
    planner->num_nodes = 0;
    planner->tree[0].position = *start;
    planner->tree[0].parent_index = -1;
    planner->num_nodes = 1;

    int goal_idx = -1;

    // RRT main loop
    for (int iter = 0; iter < planner->max_iterations; iter++) {
        position3d_t sample;

        // Sample (10% bias toward goal)
        if ((float)rand() / RAND_MAX < 0.1f) {
            sample = *goal;
        } else {
            random_sample(planner, &sample);
        }

        // Find nearest node
        int nearest_idx = nearest_node(planner, &sample);

        // Steer toward sample
        position3d_t new_pos;
        steer(&planner->tree[nearest_idx].position, &sample, planner->step_size, &new_pos);

        // Check collision
        if (!is_collision_free_rrt(planner, &planner->tree[nearest_idx].position, &new_pos)) {
            continue;
        }

        // Add to tree
        if (planner->num_nodes >= PF_MAX_RRT_NODES) break;

        planner->tree[planner->num_nodes].position = new_pos;
        planner->tree[planner->num_nodes].parent_index = nearest_idx;
        planner->num_nodes++;

        // Check if goal reached
        float dist_to_goal = distance_3d(&new_pos, goal);
        if (dist_to_goal < goal_tolerance) {
            goal_idx = planner->num_nodes - 1;
            break;
        }
    }

    if (goal_idx < 0) {
        return -1;
    }

    // Reconstruct path
    int path_len = 0;
    int current_idx = goal_idx;

    while (current_idx >= 0 && path_len < PF_MAX_WAYPOINTS) {
        path_len++;
        current_idx = planner->tree[current_idx].parent_index;
    }

    current_idx = goal_idx;
    int idx = path_len - 1;

    while (current_idx >= 0 && idx >= 0) {
        path[idx] = planner->tree[current_idx].position;
        idx--;
        current_idx = planner->tree[current_idx].parent_index;
    }

    *num_waypoints = path_len;
    return 0;
}

/* ============================================================================
 * High-Level Path Planner
 * ============================================================================ */

int path_planner_init(path_planner_t* planner,
                      float x_min, float x_max,
                      float y_min, float y_max,
                      pathfinding_algorithm_t algorithm) {
    if (!planner) return -1;

    planner->x_min = x_min;
    planner->x_max = x_max;
    planner->y_min = y_min;
    planner->y_max = y_max;
    planner->algorithm = algorithm;
    planner->altitude = PF_DEFAULT_ALTITUDE;
    planner->num_obstacles = 0;
    planner->signal_field = NULL;

    // Initialize sub-planners
    if (astar_planner_init(&planner->astar, 10.0f, x_min, x_max, y_min, y_max) != 0) {
        return -1;
    }

    if (rrt_planner_init(&planner->rrt, x_min, x_max, y_min, y_max, 50.0f, 5000) != 0) {
        astar_planner_free(&planner->astar);
        return -1;
    }

    return 0;
}

void path_planner_free(path_planner_t* planner) {
    if (!planner) return;

    astar_planner_free(&planner->astar);
    rrt_planner_free(&planner->rrt);
}

int path_planner_add_obstacle(path_planner_t* planner,
                               const position3d_t* position,
                               float radius) {
    if (!planner || !position) return -1;
    if (planner->num_obstacles >= PF_MAX_OBSTACLES) return -1;

    obstacle_t obstacle = {
        .position = *position,
        .radius = radius,
        .active = true
    };

    planner->obstacles[planner->num_obstacles] = obstacle;
    planner->num_obstacles++;

    astar_add_obstacle(&planner->astar, &obstacle);
    rrt_add_obstacle(&planner->rrt, &obstacle);

    return 0;
}

void path_planner_set_signal_field(path_planner_t* planner, void* signal_field) {
    if (!planner) return;
    planner->signal_field = signal_field;
    astar_set_signal_field(&planner->astar, signal_field);
}

int path_planner_plan(path_planner_t* planner,
                      const position3d_t* start,
                      const position3d_t* goal,
                      bool optimize_for_signal,
                      waypoint_t* waypoints,
                      int* num_waypoints) {
    if (!planner || !start || !goal || !waypoints || !num_waypoints) return -1;

    position3d_t raw_path[PF_MAX_WAYPOINTS];
    int raw_count = 0;
    int result = -1;

    // Plan based on algorithm
    if (planner->algorithm == PF_ASTAR) {
        float weight_signal = optimize_for_signal ? 0.1f : 0.0f;
        result = astar_plan(&planner->astar, start, goal, 1.0f, weight_signal,
                           raw_path, &raw_count);
    } else if (planner->algorithm == PF_RRT || planner->algorithm == PF_RRT_STAR) {
        result = rrt_plan(&planner->rrt, start, goal, 50.0f, raw_path, &raw_count);
    }

    if (result != 0 || raw_count == 0) {
        return -1;
    }

    // Convert to waypoints
    for (int i = 0; i < raw_count && i < PF_MAX_WAYPOINTS; i++) {
        waypoints[i].position = raw_path[i];
        waypoints[i].position.z = planner->altitude;
        waypoints[i].speed = PF_DEFAULT_SPEED;
        waypoints[i].signal_strength = -150.0f;
        waypoints[i].action[0] = '\0';

        // Calculate heading
        if (i < raw_count - 1) {
            waypoints[i].heading = calculate_heading(&raw_path[i], &raw_path[i+1]);
        } else {
            waypoints[i].heading = 0.0f;
        }
    }

    *num_waypoints = raw_count;
    return 0;
}

int path_planner_replan(path_planner_t* planner,
                        const position3d_t* current_position,
                        const position3d_t* goal,
                        bool optimize_for_signal,
                        waypoint_t* waypoints,
                        int* num_waypoints) {
    return path_planner_plan(planner, current_position, goal, optimize_for_signal,
                            waypoints, num_waypoints);
}

/* ============================================================================
 * Path Utilities
 * ============================================================================ */

float path_calculate_length(const position3d_t* path, int num_points) {
    if (!path || num_points < 2) return 0.0f;

    float length = 0.0f;
    for (int i = 0; i < num_points - 1; i++) {
        length += distance_3d(&path[i], &path[i+1]);
    }

    return length;
}

int path_smooth(const position3d_t* path,
                int num_points,
                const obstacle_t* obstacles,
                int num_obstacles,
                position3d_t* smoothed_path,
                int* num_smoothed) {
    if (!path || !smoothed_path || !num_smoothed || num_points < 2) return -1;

    // Simple path smoothing: remove unnecessary waypoints
    smoothed_path[0] = path[0];
    int smoothed_count = 1;

    int i = 0;
    while (i < num_points - 1) {
        // Try to skip as far ahead as possible
        bool found = false;
        for (int j = num_points - 1; j > i; j--) {
            // Check if line from i to j is clear
            bool clear = true;
            for (int k = 0; k < num_obstacles; k++) {
                if (obstacle_collides_line(&obstacles[k], &path[i], &path[j])) {
                    clear = false;
                    break;
                }
            }

            if (clear) {
                smoothed_path[smoothed_count++] = path[j];
                i = j;
                found = true;
                break;
            }
        }

        if (!found) {
            i++;
        }
    }

    *num_smoothed = smoothed_count;
    return 0;
}

int path_to_waypoints(const position3d_t* positions,
                      int num_positions,
                      float altitude,
                      float speed,
                      void* signal_field,
                      waypoint_t* waypoints) {
    if (!positions || !waypoints || num_positions <= 0) return -1;

    for (int i = 0; i < num_positions; i++) {
        waypoints[i].position = positions[i];
        waypoints[i].position.z = altitude;
        waypoints[i].speed = speed;
        waypoints[i].signal_strength = -150.0f;
        waypoints[i].action[0] = '\0';

        if (i < num_positions - 1) {
            waypoints[i].heading = calculate_heading(&positions[i], &positions[i+1]);
        } else {
            waypoints[i].heading = 0.0f;
        }
    }

    return 0;
}
