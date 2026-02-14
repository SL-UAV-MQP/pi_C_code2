/**
 * @file pathfinding.h
 * @brief Pathfinding and Waypoint Generation for UAV Navigation
 *
 * Implements A* and RRT pathfinding algorithms for UAV navigation.
 * Integrates with signal field generation and SLAM for optimal path planning.
 *
 * Features:
 * - A* grid-based pathfinding with obstacle avoidance
 * - RRT (Rapidly-exploring Random Tree) sampling-based planning
 * - Signal-aware path optimization
 * - Obstacle representation and collision checking
 * - Path smoothing and waypoint generation
 *
 * Performance Targets:
 * - A* planning: <500ms for 300x300 grid
 * - RRT planning: <1s for complex environments
 * - Path smoothing: <50ms
 */

#ifndef PATHFINDING_H
#define PATHFINDING_H

#include "common_types.h"
#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of obstacles */
#define PF_MAX_OBSTACLES 50

/** Maximum grid dimension (one axis) */
#define PF_MAX_GRID_DIM 500

/** Maximum number of waypoints in path */
#define PF_MAX_WAYPOINTS 500

/** Maximum RRT tree nodes */
#define PF_MAX_RRT_NODES 10000

/** Default UAV altitude (meters) */
#define PF_DEFAULT_ALTITUDE 100.0

/** Default UAV speed (m/s) */
#define PF_DEFAULT_SPEED 5.0

/* ============================================================================
 * Pathfinding Algorithms
 * ============================================================================ */

/**
 * @brief Pathfinding algorithm selection
 */
typedef enum {
    PF_ASTAR = 0,            /**< A* grid-based search */
    PF_RRT,                  /**< RRT sampling-based */
    PF_RRT_STAR,             /**< RRT* (optimized) */
    PF_GRADIENT              /**< Gradient ascent/descent */
} pathfinding_algorithm_t;

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/* position3d_t provided by common_types.h (double precision) */

/**
 * @brief Navigation waypoint
 */
typedef struct {
    position3d_t position;       /**< 3D position [x, y, z] */
    float heading;               /**< Desired heading (degrees) */
    float speed;                 /**< Desired speed (m/s) */
    char action[32];             /**< Special action at waypoint */
    float signal_strength;       /**< Expected signal strength (dBm) */
} waypoint_t;

/**
 * @brief Obstacle representation
 */
typedef struct {
    position3d_t position;       /**< Obstacle center position */
    float radius;                /**< Obstacle radius (meters) */
    bool active;                 /**< Whether obstacle is active */
} obstacle_t;

/**
 * @brief A* node for priority queue
 */
typedef struct astar_node {
    float f_score;               /**< f = g + h (total cost) */
    float g_score;               /**< Cost from start */
    int grid_x;                  /**< Grid X index */
    int grid_y;                  /**< Grid Y index */
    struct astar_node* parent;   /**< Parent node pointer */
} astar_node_t;

/**
 * @brief RRT tree node
 */
typedef struct {
    position3d_t position;       /**< Node position (2D: z=0) */
    int parent_index;            /**< Parent node index (-1 for root) */
} rrt_node_t;

/* ============================================================================
 * A* Planner
 * ============================================================================ */

/**
 * @brief A* pathfinding planner state
 */
typedef struct {
    /* Grid configuration */
    float grid_resolution;       /**< Grid cell size (meters) */
    float x_min, x_max;         /**< X bounds (meters) */
    float y_min, y_max;         /**< Y bounds (meters) */
    int nx, ny;                 /**< Grid dimensions */

    /* Obstacle map */
    bool* obstacle_map;         /**< Obstacle grid [ny x nx] */

    /* Obstacles */
    obstacle_t obstacles[PF_MAX_OBSTACLES];
    int num_obstacles;

    /* Signal field (optional) */
    void* signal_field;         /**< Pointer to signal field generator */
} astar_planner_t;

/**
 * @brief Initialize A* planner
 *
 * @param planner Pointer to planner structure (allocated by caller)
 * @param grid_resolution Grid cell size (meters)
 * @param x_min Minimum X coordinate (meters)
 * @param x_max Maximum X coordinate (meters)
 * @param y_min Minimum Y coordinate (meters)
 * @param y_max Maximum Y coordinate (meters)
 * @return 0 on success, negative error code on failure
 */
int astar_planner_init(astar_planner_t* planner,
                       float grid_resolution,
                       float x_min, float x_max,
                       float y_min, float y_max);

/**
 * @brief Free A* planner resources
 *
 * @param planner Planner to free
 */
void astar_planner_free(astar_planner_t* planner);

/**
 * @brief Add obstacle to environment
 *
 * @param planner Planner instance
 * @param obstacle Obstacle to add
 * @return 0 on success, negative error code on failure
 */
int astar_add_obstacle(astar_planner_t* planner, const obstacle_t* obstacle);

/**
 * @brief Set signal field for signal-aware planning
 *
 * @param planner Planner instance
 * @param signal_field Signal field generator pointer
 */
void astar_set_signal_field(astar_planner_t* planner, void* signal_field);

/**
 * @brief Find path from start to goal using A*
 *
 * Uses A* search with configurable cost weights for distance and signal quality.
 *
 * Performance: <500ms for 300x300 grid
 *
 * @param planner Planner instance
 * @param start Start position [x, y, z]
 * @param goal Goal position [x, y, z]
 * @param weight_distance Weight for distance cost
 * @param weight_signal Weight for signal quality (negative = prefer high signal)
 * @param path Output waypoints (allocated by caller) [PF_MAX_WAYPOINTS]
 * @param num_waypoints Output number of waypoints
 * @return 0 on success, negative error code on failure
 */
int astar_plan(astar_planner_t* planner,
               const position3d_t* start,
               const position3d_t* goal,
               float weight_distance,
               float weight_signal,
               position3d_t* path,
               int* num_waypoints);

/* ============================================================================
 * RRT Planner
 * ============================================================================ */

/**
 * @brief RRT planner state
 */
typedef struct {
    /* Planning bounds */
    float x_min, x_max;
    float y_min, y_max;

    /* RRT parameters */
    float step_size;             /**< Maximum step size (meters) */
    int max_iterations;          /**< Maximum iterations */

    /* Obstacles */
    obstacle_t obstacles[PF_MAX_OBSTACLES];
    int num_obstacles;

    /* Tree nodes */
    rrt_node_t tree[PF_MAX_RRT_NODES];
    int num_nodes;
} rrt_planner_t;

/**
 * @brief Initialize RRT planner
 *
 * @param planner Pointer to planner structure (allocated by caller)
 * @param x_min Minimum X coordinate (meters)
 * @param x_max Maximum X coordinate (meters)
 * @param y_min Minimum Y coordinate (meters)
 * @param y_max Maximum Y coordinate (meters)
 * @param step_size Maximum step size for expansion (meters)
 * @param max_iterations Maximum iterations
 * @return 0 on success, negative error code on failure
 */
int rrt_planner_init(rrt_planner_t* planner,
                     float x_min, float x_max,
                     float y_min, float y_max,
                     float step_size,
                     int max_iterations);

/**
 * @brief Free RRT planner resources
 *
 * @param planner Planner to free
 */
void rrt_planner_free(rrt_planner_t* planner);

/**
 * @brief Add obstacle to environment
 *
 * @param planner Planner instance
 * @param obstacle Obstacle to add
 * @return 0 on success, negative error code on failure
 */
int rrt_add_obstacle(rrt_planner_t* planner, const obstacle_t* obstacle);

/**
 * @brief Find path using RRT
 *
 * Performance: <1s for typical environments
 *
 * @param planner Planner instance
 * @param start Start position [x, y, z]
 * @param goal Goal position [x, y, z]
 * @param goal_tolerance Distance to goal for success (meters)
 * @param path Output waypoints (allocated by caller) [PF_MAX_WAYPOINTS]
 * @param num_waypoints Output number of waypoints
 * @return 0 on success, negative error code on failure
 */
int rrt_plan(rrt_planner_t* planner,
             const position3d_t* start,
             const position3d_t* goal,
             float goal_tolerance,
             position3d_t* path,
             int* num_waypoints);

/* ============================================================================
 * High-Level Path Planner
 * ============================================================================ */

/**
 * @brief High-level path planner integrating multiple algorithms
 */
typedef struct {
    /* Configuration */
    float x_min, x_max;
    float y_min, y_max;
    pathfinding_algorithm_t algorithm;
    float altitude;              /**< Default UAV altitude */

    /* Sub-planners */
    astar_planner_t astar;
    rrt_planner_t rrt;

    /* Obstacles */
    obstacle_t obstacles[PF_MAX_OBSTACLES];
    int num_obstacles;

    /* Signal field */
    void* signal_field;
} path_planner_t;

/**
 * @brief Initialize path planner
 *
 * @param planner Pointer to planner structure (allocated by caller)
 * @param x_min Minimum X bound (meters)
 * @param x_max Maximum X bound (meters)
 * @param y_min Minimum Y bound (meters)
 * @param y_max Maximum Y bound (meters)
 * @param algorithm Pathfinding algorithm to use
 * @return 0 on success, negative error code on failure
 */
int path_planner_init(path_planner_t* planner,
                      float x_min, float x_max,
                      float y_min, float y_max,
                      pathfinding_algorithm_t algorithm);

/**
 * @brief Free path planner resources
 *
 * @param planner Planner to free
 */
void path_planner_free(path_planner_t* planner);

/**
 * @brief Add circular obstacle
 *
 * @param planner Planner instance
 * @param position Obstacle center position
 * @param radius Obstacle radius (meters)
 * @return 0 on success, negative error code on failure
 */
int path_planner_add_obstacle(path_planner_t* planner,
                               const position3d_t* position,
                               float radius);

/**
 * @brief Set signal field for signal-aware planning
 *
 * @param planner Planner instance
 * @param signal_field Signal field generator pointer
 */
void path_planner_set_signal_field(path_planner_t* planner, void* signal_field);

/**
 * @brief Plan path from start to goal
 *
 * @param planner Planner instance
 * @param start Start position [x, y, z]
 * @param goal Goal position [x, y, z]
 * @param optimize_for_signal Whether to optimize for signal strength
 * @param waypoints Output waypoints (allocated by caller) [PF_MAX_WAYPOINTS]
 * @param num_waypoints Output number of waypoints
 * @return 0 on success, negative error code on failure
 */
int path_planner_plan(path_planner_t* planner,
                      const position3d_t* start,
                      const position3d_t* goal,
                      bool optimize_for_signal,
                      waypoint_t* waypoints,
                      int* num_waypoints);

/**
 * @brief Replan path from current position
 *
 * @param planner Planner instance
 * @param current_position Current UAV position
 * @param goal Goal position
 * @param optimize_for_signal Optimize for signal
 * @param waypoints Output waypoints [PF_MAX_WAYPOINTS]
 * @param num_waypoints Output number of waypoints
 * @return 0 on success, negative error code on failure
 */
int path_planner_replan(path_planner_t* planner,
                        const position3d_t* current_position,
                        const position3d_t* goal,
                        bool optimize_for_signal,
                        waypoint_t* waypoints,
                        int* num_waypoints);

/* ============================================================================
 * Obstacle Utilities
 * ============================================================================ */

/**
 * @brief Check if point collides with obstacle
 *
 * @param obstacle Obstacle instance
 * @param point Point to check
 * @return true if collision, false otherwise
 */
bool obstacle_collides_point(const obstacle_t* obstacle, const position3d_t* point);

/**
 * @brief Check if line segment collides with obstacle
 *
 * @param obstacle Obstacle instance
 * @param start Line segment start
 * @param end Line segment end
 * @return true if collision, false otherwise
 */
bool obstacle_collides_line(const obstacle_t* obstacle,
                            const position3d_t* start,
                            const position3d_t* end);

/* ============================================================================
 * Path Utilities
 * ============================================================================ */

/**
 * @brief Calculate path length
 *
 * @param path Array of 3D positions
 * @param num_points Number of points in path
 * @return Total path length (meters)
 */
float path_calculate_length(const position3d_t* path, int num_points);

/**
 * @brief Smooth path by removing unnecessary waypoints
 *
 * @param path Input path
 * @param num_points Input number of points
 * @param obstacles Obstacles for collision checking
 * @param num_obstacles Number of obstacles
 * @param smoothed_path Output smoothed path [PF_MAX_WAYPOINTS]
 * @param num_smoothed Output number of smoothed points
 * @return 0 on success, negative error code on failure
 */
int path_smooth(const position3d_t* path,
                int num_points,
                const obstacle_t* obstacles,
                int num_obstacles,
                position3d_t* smoothed_path,
                int* num_smoothed);

/**
 * @brief Convert position array to waypoint array
 *
 * Adds heading, speed, and signal strength to waypoints.
 *
 * @param positions Input positions
 * @param num_positions Number of positions
 * @param altitude Default altitude (meters)
 * @param speed Default speed (m/s)
 * @param signal_field Signal field (NULL if not available)
 * @param waypoints Output waypoints [num_positions]
 * @return 0 on success, negative error code on failure
 */
int path_to_waypoints(const position3d_t* positions,
                      int num_positions,
                      float altitude,
                      float speed,
                      void* signal_field,
                      waypoint_t* waypoints);

/* ============================================================================
 * Geometric Utilities
 * ============================================================================ */

/**
 * @brief Calculate Euclidean distance between two points
 *
 * @param p1 First point
 * @param p2 Second point
 * @return Distance (meters)
 */
float distance_3d(const position3d_t* p1, const position3d_t* p2);

/**
 * @brief Calculate heading from one point to another
 *
 * @param from Start point
 * @param to End point
 * @return Heading in degrees (0-360, 0=North, 90=East)
 */
float calculate_heading(const position3d_t* from, const position3d_t* to);

#endif /* PATHFINDING_H */
