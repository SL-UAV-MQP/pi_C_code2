/**
 * @file map_fitting.c
 * @brief Map Fitting and Alignment for SLAM Implementation
 *
 * Implements map alignment, loop closure detection, bundle adjustment,
 * and map fusion for SLAM applications.
 *
 * Features:
 * - Procrustes, Horn, Umeyama alignment algorithms
 * - Similarity transformation (rotation, translation, scale)
 * - Loop closure detection for trajectory consistency
 * - Bundle adjustment for map optimization
 * - Map fusion with ground truth landmarks
 */

#include "map_fitting.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

/* ============================================================================
 * Mathematical Constants
 * ============================================================================ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPSILON 1e-10
#define MIN_TIME_SEPARATION 30.0  /* Minimum time between poses for loop closure (seconds) */

/* ============================================================================
 * Internal Helper Functions - Matrix Operations
 * ============================================================================ */

/**
 * @brief Compute centroid of point set
 */
static void compute_centroid(const position3d_t* points, int num_points, position3d_t* centroid) {
    centroid->x = 0.0;
    centroid->y = 0.0;
    centroid->z = 0.0;

    for (int i = 0; i < num_points; i++) {
        centroid->x += points[i].x;
        centroid->y += points[i].y;
        centroid->z += points[i].z;
    }

    centroid->x /= num_points;
    centroid->y /= num_points;
    centroid->z /= num_points;
}

/**
 * @brief Center point set by subtracting centroid
 */
static void center_points(const position3d_t* points, int num_points,
                         const position3d_t* centroid, position3d_t* centered) {
    for (int i = 0; i < num_points; i++) {
        centered[i].x = points[i].x - centroid->x;
        centered[i].y = points[i].y - centroid->y;
        centered[i].z = points[i].z - centroid->z;
    }
}

/**
 * @brief Compute 3x3 cross-covariance matrix H = sum(source_i * target_i^T)
 */
static void compute_cross_covariance(const position3d_t* source, const position3d_t* target,
                                    int num_points, double H[9]) {
    memset(H, 0, 9 * sizeof(double));

    for (int i = 0; i < num_points; i++) {
        H[0] += source[i].x * target[i].x;  // H[0][0]
        H[1] += source[i].x * target[i].y;  // H[0][1]
        H[2] += source[i].x * target[i].z;  // H[0][2]
        H[3] += source[i].y * target[i].x;  // H[1][0]
        H[4] += source[i].y * target[i].y;  // H[1][1]
        H[5] += source[i].y * target[i].z;  // H[1][2]
        H[6] += source[i].z * target[i].x;  // H[2][0]
        H[7] += source[i].z * target[i].y;  // H[2][1]
        H[8] += source[i].z * target[i].z;  // H[2][2]
    }
}

/**
 * @brief Compute sum of squared norms
 */
static double compute_sum_squared_norms(const position3d_t* points, int num_points) {
    double sum = 0.0;
    for (int i = 0; i < num_points; i++) {
        sum += points[i].x * points[i].x +
               points[i].y * points[i].y +
               points[i].z * points[i].z;
    }
    return sum;
}

/**
 * @brief Compute optimal rotation matrix from cross-covariance H using SVD
 *
 * H = U * S * V^T  →  R = V * diag(1,1,det(V*U^T)) * U^T
 * The determinant correction ensures a proper rotation (det=+1, not reflection).
 * Uses LAPACK dgesvd for numerical robustness.
 */
static int compute_rotation_svd(const double H[9], rotation_matrix_t* rotation) {
    // H is stored row-major: H[i*3+j] = H(i,j)
    // LAPACK expects column-major, so we pass H^T which is H with swapped indices
    // Since SVD(H^T) = V * S * U^T, we adjust accordingly

    // Copy H to column-major format for LAPACK (transpose)
    double A[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[j * 3 + i] = H[i * 3 + j];  // column-major
        }
    }

    double U[9], VT[9], S[3];
    double superb[2];

    // SVD: A = U * S * VT  (LAPACK column-major)
    lapack_int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A',
                                      3, 3, A, 3, S, U, 3, VT, 3, superb);

    if (info != 0) {
        // SVD failed → identity rotation
        memset(rotation->data, 0, 9 * sizeof(double));
        rotation->data[0] = rotation->data[4] = rotation->data[8] = 1.0;
        return -1;
    }

    // Compute V = VT^T (column-major)
    double V[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            V[i * 3 + j] = VT[j * 3 + i];
        }
    }

    // Compute det(V * U^T) to ensure proper rotation
    // W = V * U^T (column-major)
    double W[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            W[j * 3 + i] = 0.0;
            for (int k = 0; k < 3; k++) {
                W[j * 3 + i] += V[k * 3 + i] * U[j * 3 + k];
            }
        }
    }

    // Determinant of 3x3 matrix (column-major: W[col*3+row])
    double det = W[0] * (W[4] * W[8] - W[7] * W[5])
               - W[3] * (W[1] * W[8] - W[7] * W[2])
               + W[6] * (W[1] * W[5] - W[4] * W[2]);

    // Correction diagonal: d = diag(1, 1, sign(det))
    double d = (det < 0.0) ? -1.0 : 1.0;

    // R = V * diag(1,1,d) * U^T  →  R(i,j) = sum_k V(i,k) * d_k * U(j,k)
    // In column-major: R[j*3+i] = sum_k V[k*3+i] * d_k * U[k*3+j]
    double R_cm[9];  // column-major
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_cm[j * 3 + i] = 0.0;
            for (int k = 0; k < 3; k++) {
                double dk = (k == 2) ? d : 1.0;
                R_cm[j * 3 + i] += V[k * 3 + i] * dk * U[k * 3 + j];
            }
        }
    }

    // Convert to row-major for rotation->data
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rotation->data[i * 3 + j] = R_cm[j * 3 + i];
        }
    }

    return 0;
}

/**
 * @brief Compute 4x4 symmetric matrix for Horn's method
 */
static void compute_horn_matrix(const double H[9], double N[16]) {
    double Sxx = H[0], Sxy = H[1], Sxz = H[2];
    double Syx = H[3], Syy = H[4], Syz = H[5];
    double Szx = H[6], Szy = H[7], Szz = H[8];

    // N = [
    //   [Sxx+Syy+Szz,  Syz-Szy,      Szx-Sxz,      Sxy-Syx     ],
    //   [Syz-Szy,      Sxx-Syy-Szz,  Sxy+Syx,      Szx+Sxz     ],
    //   [Szx-Sxz,      Sxy+Syx,      -Sxx+Syy-Szz, Syz+Szy     ],
    //   [Sxy-Syx,      Szx+Sxz,      Syz+Szy,      -Sxx-Syy+Szz]
    // ]

    N[0]  = Sxx + Syy + Szz;
    N[1]  = Syz - Szy;
    N[2]  = Szx - Sxz;
    N[3]  = Sxy - Syx;

    N[4]  = Syz - Szy;
    N[5]  = Sxx - Syy - Szz;
    N[6]  = Sxy + Syx;
    N[7]  = Szx + Sxz;

    N[8]  = Szx - Sxz;
    N[9]  = Sxy + Syx;
    N[10] = -Sxx + Syy - Szz;
    N[11] = Syz + Szy;

    N[12] = Sxy - Syx;
    N[13] = Szx + Sxz;
    N[14] = Syz + Szy;
    N[15] = -Sxx - Syy + Szz;
}

/**
 * @brief Power iteration to find largest eigenvector of symmetric 4x4 matrix
 */
static int power_iteration_4x4(const double A[16], double eigenvector[4], int max_iter) {
    // Initialize random vector
    eigenvector[0] = 1.0;
    eigenvector[1] = 0.0;
    eigenvector[2] = 0.0;
    eigenvector[3] = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double temp[4] = {0, 0, 0, 0};

        // temp = A * eigenvector
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                temp[i] += A[i*4 + j] * eigenvector[j];
            }
        }

        // Normalize
        double norm = sqrt(temp[0]*temp[0] + temp[1]*temp[1] +
                          temp[2]*temp[2] + temp[3]*temp[3]);

        if (norm < EPSILON) {
            return -1;
        }

        for (int i = 0; i < 4; i++) {
            eigenvector[i] = temp[i] / norm;
        }
    }

    return 0;
}

/* ============================================================================
 * Utility Functions Implementation
 * ============================================================================ */

double map_calculate_rms_error(const position3d_t* points1,
                                const position3d_t* points2,
                                int num_points) {
    if (num_points <= 0) return 0.0;

    double sum_sq_error = 0.0;
    for (int i = 0; i < num_points; i++) {
        double dx = points1[i].x - points2[i].x;
        double dy = points1[i].y - points2[i].y;
        double dz = points1[i].z - points2[i].z;
        sum_sq_error += dx*dx + dy*dy + dz*dz;
    }

    return sqrt(sum_sq_error / num_points);
}

double map_normalize_angle(double angle) {
    // Normalize to [-180, 180]
    while (angle > 180.0) angle -= 360.0;
    while (angle < -180.0) angle += 360.0;
    return angle;
}

void map_apply_rotation(const rotation_matrix_t* rotation,
                        const position3d_t* position,
                        position3d_t* rotated) {
    const double* R = rotation->data;
    rotated->x = R[0]*position->x + R[1]*position->y + R[2]*position->z;
    rotated->y = R[3]*position->x + R[4]*position->y + R[5]*position->z;
    rotated->z = R[6]*position->x + R[7]*position->y + R[8]*position->z;
}

void map_quaternion_to_rotation(double qw, double qx, double qy, double qz,
                                 rotation_matrix_t* rotation) {
    // Normalize quaternion
    double norm = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
    if (norm < EPSILON) {
        // Identity
        memset(rotation->data, 0, 9 * sizeof(double));
        rotation->data[0] = rotation->data[4] = rotation->data[8] = 1.0;
        return;
    }

    qw /= norm; qx /= norm; qy /= norm; qz /= norm;

    // Convert to rotation matrix
    double qw2 = qw*qw, qx2 = qx*qx, qy2 = qy*qy, qz2 = qz*qz;

    rotation->data[0] = qw2 + qx2 - qy2 - qz2;
    rotation->data[1] = 2*(qx*qy - qw*qz);
    rotation->data[2] = 2*(qx*qz + qw*qy);

    rotation->data[3] = 2*(qx*qy + qw*qz);
    rotation->data[4] = qw2 - qx2 + qy2 - qz2;
    rotation->data[5] = 2*(qy*qz - qw*qx);

    rotation->data[6] = 2*(qx*qz - qw*qy);
    rotation->data[7] = 2*(qy*qz + qw*qx);
    rotation->data[8] = qw2 - qx2 - qy2 + qz2;
}

int map_find_common_landmarks(const slam_map_t* map1,
                               const slam_map_t* map2,
                               char common_ids[][32],
                               int* num_common) {
    if (!map1 || !map2 || !common_ids || !num_common) {
        return -1;
    }

    *num_common = 0;

    for (int i = 0; i < map1->num_landmarks && *num_common < MAP_MAX_LANDMARKS; i++) {
        for (int j = 0; j < map2->num_landmarks; j++) {
            if (strcmp(map1->landmarks[i].id, map2->landmarks[j].id) == 0) {
                strncpy(common_ids[*num_common], map1->landmarks[i].id, 31);
                common_ids[*num_common][31] = '\0';
                (*num_common)++;
                break;
            }
        }
    }

    return 0;
}

/* ============================================================================
 * Map Alignment Implementation
 * ============================================================================ */

int map_alignment_init(map_alignment_t* alignment) {
    if (!alignment) return -1;

    // Initialize identity transform
    memset(&alignment->transform, 0, sizeof(map_transform_t));
    memset(alignment->transform.rotation.data, 0, 9 * sizeof(double));
    alignment->transform.rotation.data[0] = 1.0;
    alignment->transform.rotation.data[4] = 1.0;
    alignment->transform.rotation.data[8] = 1.0;
    alignment->transform.scale = 1.0;
    alignment->transform.alignment_error = 0.0;

    return 0;
}

void map_alignment_free(map_alignment_t* alignment) {
    // No dynamic allocation, nothing to free
    (void)alignment;
}

int map_alignment_procrustes(const position3d_t* source_pts,
                              const position3d_t* target_pts,
                              int num_points,
                              alignment_result_t* result) {
    if (!source_pts || !target_pts || !result || num_points < 3) {
        return -1;
    }

    // Allocate temporary arrays
    position3d_t* source_centered = malloc(num_points * sizeof(position3d_t));
    position3d_t* target_centered = malloc(num_points * sizeof(position3d_t));

    if (!source_centered || !target_centered) {
        free(source_centered);
        free(target_centered);
        return -1;
    }

    // Compute centroids
    position3d_t source_centroid, target_centroid;
    compute_centroid(source_pts, num_points, &source_centroid);
    compute_centroid(target_pts, num_points, &target_centroid);

    // Center point sets
    center_points(source_pts, num_points, &source_centroid, source_centered);
    center_points(target_pts, num_points, &target_centroid, target_centered);

    // Compute cross-covariance matrix H
    double H[9];
    compute_cross_covariance(source_centered, target_centered, num_points, H);

    // Compute rotation via SVD
    compute_rotation_svd(H, &result->transform.rotation);

    // Compute scale
    double source_scale = sqrt(compute_sum_squared_norms(source_centered, num_points));
    double target_scale = sqrt(compute_sum_squared_norms(target_centered, num_points));

    if (source_scale < EPSILON) {
        result->transform.scale = 1.0;
    } else {
        result->transform.scale = target_scale / source_scale;
    }

    // Compute translation: t = target_centroid - s * R * source_centroid
    position3d_t rotated_source_centroid;
    map_apply_rotation(&result->transform.rotation, &source_centroid, &rotated_source_centroid);

    result->transform.translation.x = target_centroid.x -
                                      result->transform.scale * rotated_source_centroid.x;
    result->transform.translation.y = target_centroid.y -
                                      result->transform.scale * rotated_source_centroid.y;
    result->transform.translation.z = target_centroid.z -
                                      result->transform.scale * rotated_source_centroid.z;

    // Compute alignment error
    position3d_t* transformed = malloc(num_points * sizeof(position3d_t));
    if (transformed) {
        for (int i = 0; i < num_points; i++) {
            position3d_t rotated;
            map_apply_rotation(&result->transform.rotation, &source_pts[i], &rotated);
            transformed[i].x = result->transform.scale * rotated.x + result->transform.translation.x;
            transformed[i].y = result->transform.scale * rotated.y + result->transform.translation.y;
            transformed[i].z = result->transform.scale * rotated.z + result->transform.translation.z;
        }

        result->transform.alignment_error = map_calculate_rms_error(transformed, target_pts, num_points);
        free(transformed);
    }

    result->success = true;
    result->method = ALIGN_PROCRUSTES;
    result->num_common_landmarks = num_points;
    result->error_message[0] = '\0';

    free(source_centered);
    free(target_centered);

    return 0;
}

int map_alignment_horn(const position3d_t* source_pts,
                       const position3d_t* target_pts,
                       int num_points,
                       alignment_result_t* result) {
    if (!source_pts || !target_pts || !result || num_points < 3) {
        return -1;
    }

    // Allocate temporary arrays
    position3d_t* source_centered = malloc(num_points * sizeof(position3d_t));
    position3d_t* target_centered = malloc(num_points * sizeof(position3d_t));

    if (!source_centered || !target_centered) {
        free(source_centered);
        free(target_centered);
        return -1;
    }

    // Compute centroids
    position3d_t source_centroid, target_centroid;
    compute_centroid(source_pts, num_points, &source_centroid);
    compute_centroid(target_pts, num_points, &target_centroid);

    // Center point sets
    center_points(source_pts, num_points, &source_centroid, source_centered);
    center_points(target_pts, num_points, &target_centroid, target_centered);

    // Compute cross-covariance matrix H
    double H[9];
    compute_cross_covariance(source_centered, target_centered, num_points, H);

    // Compute Horn's symmetric matrix N
    double N[16];
    compute_horn_matrix(H, N);

    // Find largest eigenvector (quaternion)
    double q[4];
    if (power_iteration_4x4(N, q, 100) != 0) {
        free(source_centered);
        free(target_centered);
        return -1;
    }

    // Convert quaternion to rotation matrix
    map_quaternion_to_rotation(q[0], q[1], q[2], q[3], &result->transform.rotation);

    // Compute scale (similar to Procrustes)
    double source_scale = sqrt(compute_sum_squared_norms(source_centered, num_points));
    double target_scale = sqrt(compute_sum_squared_norms(target_centered, num_points));

    if (source_scale < EPSILON) {
        result->transform.scale = 1.0;
    } else {
        result->transform.scale = target_scale / source_scale;
    }

    // Compute translation
    position3d_t rotated_source_centroid;
    map_apply_rotation(&result->transform.rotation, &source_centroid, &rotated_source_centroid);

    result->transform.translation.x = target_centroid.x -
                                      result->transform.scale * rotated_source_centroid.x;
    result->transform.translation.y = target_centroid.y -
                                      result->transform.scale * rotated_source_centroid.y;
    result->transform.translation.z = target_centroid.z -
                                      result->transform.scale * rotated_source_centroid.z;

    // Compute alignment error
    position3d_t* transformed = malloc(num_points * sizeof(position3d_t));
    if (transformed) {
        for (int i = 0; i < num_points; i++) {
            position3d_t rotated;
            map_apply_rotation(&result->transform.rotation, &source_pts[i], &rotated);
            transformed[i].x = result->transform.scale * rotated.x + result->transform.translation.x;
            transformed[i].y = result->transform.scale * rotated.y + result->transform.translation.y;
            transformed[i].z = result->transform.scale * rotated.z + result->transform.translation.z;
        }

        result->transform.alignment_error = map_calculate_rms_error(transformed, target_pts, num_points);
        free(transformed);
    }

    result->success = true;
    result->method = ALIGN_HORN;
    result->num_common_landmarks = num_points;
    result->error_message[0] = '\0';

    free(source_centered);
    free(target_centered);

    return 0;
}

int map_alignment_umeyama(const position3d_t* source_pts,
                          const position3d_t* target_pts,
                          int num_points,
                          bool estimate_scale,
                          alignment_result_t* result) {
    // Umeyama is similar to Procrustes but handles scale differently
    // For now, delegate to Procrustes
    int ret = map_alignment_procrustes(source_pts, target_pts, num_points, result);

    if (ret == 0) {
        result->method = ALIGN_UMEYAMA;

        if (!estimate_scale) {
            result->transform.scale = 1.0;
        }
    }

    return ret;
}

int map_alignment_align(map_alignment_t* alignment,
                        const slam_map_t* slam_map,
                        const slam_map_t* ground_truth_map,
                        alignment_method_t method,
                        alignment_result_t* result) {
    if (!alignment || !slam_map || !ground_truth_map || !result) {
        return -1;
    }

    // Find common landmarks
    char common_ids[MAP_MAX_LANDMARKS][32];
    int num_common;

    if (map_find_common_landmarks(slam_map, ground_truth_map, common_ids, &num_common) != 0) {
        result->success = false;
        snprintf(result->error_message, sizeof(result->error_message),
                "Failed to find common landmarks");
        return -1;
    }

    if (num_common < 3) {
        result->success = false;
        snprintf(result->error_message, sizeof(result->error_message),
                "Insufficient common landmarks: %d (need >= 3)", num_common);
        return -1;
    }

    // Extract matched points
    position3d_t* source_pts = malloc(num_common * sizeof(position3d_t));
    position3d_t* target_pts = malloc(num_common * sizeof(position3d_t));

    if (!source_pts || !target_pts) {
        free(source_pts);
        free(target_pts);
        result->success = false;
        snprintf(result->error_message, sizeof(result->error_message), "Memory allocation failed");
        return -1;
    }

    for (int i = 0; i < num_common; i++) {
        // Find in slam_map
        for (int j = 0; j < slam_map->num_landmarks; j++) {
            if (strcmp(slam_map->landmarks[j].id, common_ids[i]) == 0) {
                source_pts[i] = slam_map->landmarks[j].position;
                break;
            }
        }

        // Find in ground_truth_map
        for (int j = 0; j < ground_truth_map->num_landmarks; j++) {
            if (strcmp(ground_truth_map->landmarks[j].id, common_ids[i]) == 0) {
                target_pts[i] = ground_truth_map->landmarks[j].position;
                break;
            }
        }
    }

    // Call alignment method
    int ret;
    switch (method) {
        case ALIGN_PROCRUSTES:
            ret = map_alignment_procrustes(source_pts, target_pts, num_common, result);
            break;
        case ALIGN_HORN:
            ret = map_alignment_horn(source_pts, target_pts, num_common, result);
            break;
        case ALIGN_UMEYAMA:
            ret = map_alignment_umeyama(source_pts, target_pts, num_common, true, result);
            break;
        default:
            ret = -1;
            result->success = false;
            snprintf(result->error_message, sizeof(result->error_message), "Unknown alignment method");
    }

    if (ret == 0) {
        // Store transform in alignment context
        alignment->transform = result->transform;
    }

    free(source_pts);
    free(target_pts);

    return ret;
}

int map_alignment_transform_position(const map_alignment_t* alignment,
                                      const position3d_t* slam_position,
                                      position3d_t* transformed) {
    if (!alignment || !slam_position || !transformed) {
        return -1;
    }

    // Apply transformation: transformed = scale * R * slam_position + t
    position3d_t rotated;
    map_apply_rotation(&alignment->transform.rotation, slam_position, &rotated);

    transformed->x = alignment->transform.scale * rotated.x + alignment->transform.translation.x;
    transformed->y = alignment->transform.scale * rotated.y + alignment->transform.translation.y;
    transformed->z = alignment->transform.scale * rotated.z + alignment->transform.translation.z;

    return 0;
}

int map_alignment_transform_map(const map_alignment_t* alignment,
                                 const slam_map_t* slam_map,
                                 slam_map_t* transformed_map) {
    if (!alignment || !slam_map || !transformed_map) {
        return -1;
    }

    transformed_map->num_landmarks = slam_map->num_landmarks;

    for (int i = 0; i < slam_map->num_landmarks; i++) {
        // Copy ID
        strncpy(transformed_map->landmarks[i].id, slam_map->landmarks[i].id, 31);
        transformed_map->landmarks[i].id[31] = '\0';

        // Transform position
        map_alignment_transform_position(alignment,
                                        &slam_map->landmarks[i].position,
                                        &transformed_map->landmarks[i].position);

        // Transform covariance (simplified: copy for now)
        // Proper transformation: Σ' = s^2 * R * Σ * R^T
        memcpy(transformed_map->landmarks[i].uncertainty,
               slam_map->landmarks[i].uncertainty,
               9 * sizeof(double));
    }

    return 0;
}

/* ============================================================================
 * Loop Closure Detection Implementation
 * ============================================================================ */

int loop_closure_detector_init(loop_closure_detector_t* detector,
                                double distance_threshold,
                                double angle_threshold) {
    if (!detector) return -1;

    detector->distance_threshold = distance_threshold;
    detector->angle_threshold = angle_threshold;
    detector->num_poses = 0;

    return 0;
}

void loop_closure_detector_free(loop_closure_detector_t* detector) {
    // No dynamic allocation, nothing to free
    (void)detector;
}

int loop_closure_add_pose(loop_closure_detector_t* detector,
                          double timestamp,
                          const position3d_t* position,
                          double heading) {
    if (!detector || !position) return -1;

    if (detector->num_poses >= MAP_MAX_TRAJECTORY_POSES) {
        // Shift array (remove oldest)
        memmove(&detector->trajectory[0], &detector->trajectory[1],
                (MAP_MAX_TRAJECTORY_POSES - 1) * sizeof(trajectory_pose_t));
        detector->num_poses = MAP_MAX_TRAJECTORY_POSES - 1;
    }

    trajectory_pose_t* pose = &detector->trajectory[detector->num_poses];
    pose->timestamp = timestamp;
    pose->position = *position;
    pose->heading = heading;

    detector->num_poses++;

    return 0;
}

int loop_closure_detect(loop_closure_detector_t* detector,
                        const position3d_t* current_position,
                        double current_heading,
                        loop_closure_t* closures,
                        int* num_closures) {
    if (!detector || !current_position || !closures || !num_closures) {
        return -1;
    }

    *num_closures = 0;

    if (detector->num_poses == 0) {
        return 0;
    }

    double current_time = detector->trajectory[detector->num_poses - 1].timestamp;

    for (int i = 0; i < detector->num_poses - 1; i++) {
        trajectory_pose_t* hist_pose = &detector->trajectory[i];

        // Check time separation
        if (current_time - hist_pose->timestamp < MIN_TIME_SEPARATION) {
            continue;
        }

        // Compute distance
        double dx = current_position->x - hist_pose->position.x;
        double dy = current_position->y - hist_pose->position.y;
        double dz = current_position->z - hist_pose->position.z;
        double distance = sqrt(dx*dx + dy*dy + dz*dz);

        // Check distance threshold
        if (distance > detector->distance_threshold) {
            continue;
        }

        // Check heading difference
        double heading_diff = fabs(map_normalize_angle(current_heading - hist_pose->heading));
        if (heading_diff > detector->angle_threshold) {
            continue;
        }

        // Loop closure detected
        if (*num_closures < MAP_MAX_LANDMARKS) {
            closures[*num_closures].pose_index = i;
            closures[*num_closures].distance = distance;
            (*num_closures)++;
        }
    }

    return 0;
}

/* ============================================================================
 * Map Optimization (Bundle Adjustment) Implementation
 * ============================================================================ */

int map_optimizer_init(map_optimizer_t* optimizer) {
    if (!optimizer) return -1;

    optimizer->max_iterations = 100;
    optimizer->convergence_threshold = 1e-6;

    return 0;
}

void map_optimizer_free(map_optimizer_t* optimizer) {
    // No dynamic allocation, nothing to free
    (void)optimizer;
}

/**
 * @brief Compute reprojection error for a constraint
 */
static double compute_reprojection_error(const position3d_t* landmark_pos,
                                        const measurement_constraint_t* constraint) {
    // Compute expected azimuth/elevation from observer to landmark
    double dx = landmark_pos->x - constraint->observer_pos.x;
    double dy = landmark_pos->y - constraint->observer_pos.y;
    double dz = landmark_pos->z - constraint->observer_pos.z;

    double expected_azimuth = atan2(dy, dx) * 180.0 / M_PI;
    double horizontal_dist = sqrt(dx*dx + dy*dy);
    double expected_elevation = atan2(dz, horizontal_dist) * 180.0 / M_PI;

    // Compute error
    double azimuth_error = map_normalize_angle(constraint->azimuth - expected_azimuth);
    double elevation_error = constraint->elevation - expected_elevation;

    return sqrt(azimuth_error * azimuth_error + elevation_error * elevation_error);
}

int map_optimizer_optimize(map_optimizer_t* optimizer,
                           const slam_map_t* slam_map,
                           const measurement_constraint_t* constraints,
                           int num_constraints,
                           slam_map_t* optimized_map) {
    if (!optimizer || !slam_map || !constraints || !optimized_map || num_constraints <= 0) {
        return -1;
    }

    // Initialize optimized map with slam_map
    memcpy(optimized_map, slam_map, sizeof(slam_map_t));

    // Gauss-Newton optimization
    // For simplicity, implement basic gradient descent
    double learning_rate = 0.01;
    double prev_total_error = 1e10;

    for (int iter = 0; iter < optimizer->max_iterations; iter++) {
        double total_error = 0.0;
        position3d_t gradients[MAP_MAX_LANDMARKS] = {0};

        // Compute total error and gradients
        for (int c = 0; c < num_constraints; c++) {
            const measurement_constraint_t* constraint = &constraints[c];

            // Find landmark
            int landmark_idx = -1;
            for (int i = 0; i < optimized_map->num_landmarks; i++) {
                if (strcmp(optimized_map->landmarks[i].id, constraint->landmark_id) == 0) {
                    landmark_idx = i;
                    break;
                }
            }

            if (landmark_idx < 0) continue;

            landmark_t* landmark = &optimized_map->landmarks[landmark_idx];

            // Compute error
            double error = compute_reprojection_error(&landmark->position, constraint);
            total_error += error * error;

            // Compute numerical gradient (finite difference)
            const double h = 1e-5;

            position3d_t perturbed_pos;

            // Gradient w.r.t. x
            perturbed_pos = landmark->position;
            perturbed_pos.x += h;
            double error_plus = compute_reprojection_error(&perturbed_pos, constraint);
            gradients[landmark_idx].x += (error_plus - error) / h;

            // Gradient w.r.t. y
            perturbed_pos = landmark->position;
            perturbed_pos.y += h;
            error_plus = compute_reprojection_error(&perturbed_pos, constraint);
            gradients[landmark_idx].y += (error_plus - error) / h;

            // Gradient w.r.t. z
            perturbed_pos = landmark->position;
            perturbed_pos.z += h;
            error_plus = compute_reprojection_error(&perturbed_pos, constraint);
            gradients[landmark_idx].z += (error_plus - error) / h;
        }

        // Check convergence
        if (fabs(prev_total_error - total_error) < optimizer->convergence_threshold) {
            break;
        }
        prev_total_error = total_error;

        // Update landmark positions
        for (int i = 0; i < optimized_map->num_landmarks; i++) {
            optimized_map->landmarks[i].position.x -= learning_rate * gradients[i].x;
            optimized_map->landmarks[i].position.y -= learning_rate * gradients[i].y;
            optimized_map->landmarks[i].position.z -= learning_rate * gradients[i].z;
        }
    }

    return 0;
}

/* ============================================================================
 * Map Fusion Implementation
 * ============================================================================ */

int map_fusion_init(map_fusion_t* fusion) {
    if (!fusion) return -1;

    if (map_alignment_init(&fusion->aligner) != 0) return -1;
    if (loop_closure_detector_init(&fusion->loop_detector,
                                    MAP_LOOP_DISTANCE_THRESHOLD,
                                    MAP_LOOP_ANGLE_THRESHOLD) != 0) return -1;
    if (map_optimizer_init(&fusion->optimizer) != 0) return -1;

    return 0;
}

void map_fusion_free(map_fusion_t* fusion) {
    if (!fusion) return;

    map_alignment_free(&fusion->aligner);
    loop_closure_detector_free(&fusion->loop_detector);
    map_optimizer_free(&fusion->optimizer);
}

int map_fusion_fuse(map_fusion_t* fusion,
                    const slam_map_t* slam_map,
                    const slam_map_t* ground_truth_map,
                    const position3d_t* uav_position,
                    map_fusion_result_t* result) {
    if (!fusion || !slam_map || !ground_truth_map || !uav_position || !result) {
        return -1;
    }

    // Step 1: Align SLAM map with ground truth
    alignment_result_t align_result;
    int ret = map_alignment_align(&fusion->aligner, slam_map, ground_truth_map,
                                   ALIGN_PROCRUSTES, &align_result);

    if (ret != 0 || !align_result.success) {
        result->success = false;
        return -1;
    }

    result->transform = align_result.transform;
    result->alignment_error = align_result.transform.alignment_error;

    // Step 2: Transform SLAM map to ground truth frame
    slam_map_t transformed_slam_map;
    map_alignment_transform_map(&fusion->aligner, slam_map, &transformed_slam_map);

    // Step 3: Transform UAV position
    map_alignment_transform_position(&fusion->aligner, uav_position,
                                     &result->transformed_uav_pos);

    // Step 4: Fuse landmarks using weighted average based on uncertainty
    result->fused_map.num_landmarks = 0;

    // Add all ground truth landmarks first
    for (int i = 0; i < ground_truth_map->num_landmarks; i++) {
        result->fused_map.landmarks[result->fused_map.num_landmarks] =
            ground_truth_map->landmarks[i];
        result->fused_map.num_landmarks++;
    }

    // Merge or add transformed SLAM landmarks
    for (int i = 0; i < transformed_slam_map.num_landmarks; i++) {
        landmark_t* slam_lm = &transformed_slam_map.landmarks[i];

        // Check if landmark exists in ground truth
        int gt_idx = -1;
        for (int j = 0; j < ground_truth_map->num_landmarks; j++) {
            if (strcmp(slam_lm->id, ground_truth_map->landmarks[j].id) == 0) {
                gt_idx = j;
                break;
            }
        }

        if (gt_idx >= 0) {
            // Fuse with ground truth using weighted average
            landmark_t* gt_lm = &ground_truth_map->landmarks[gt_idx];
            landmark_t* fused_lm = &result->fused_map.landmarks[gt_idx];

            // Simplified fusion: average positions
            // Proper fusion: use Kalman update with covariances
            fused_lm->position.x = 0.5 * (gt_lm->position.x + slam_lm->position.x);
            fused_lm->position.y = 0.5 * (gt_lm->position.y + slam_lm->position.y);
            fused_lm->position.z = 0.5 * (gt_lm->position.z + slam_lm->position.z);
        } else {
            // Add new SLAM landmark
            if (result->fused_map.num_landmarks < MAP_MAX_LANDMARKS) {
                result->fused_map.landmarks[result->fused_map.num_landmarks] = *slam_lm;
                result->fused_map.num_landmarks++;
            }
        }
    }

    result->success = true;

    return 0;
}
