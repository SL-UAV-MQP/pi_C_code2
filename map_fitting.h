/**
 * @file map_fitting.h
 * @brief Map Fitting and Alignment for SLAM
 *
 * Aligns SLAM-estimated maps with known ground truth tower positions.
 * Performs coordinate transformation, loop closure detection, and map optimization.
 *
 * Features:
 * - Map alignment (Procrustes, Horn, Umeyama methods)
 * - Similarity transformation (rotation, translation, scale)
 * - Loop closure detection for trajectory consistency
 * - Bundle adjustment for map optimization
 * - Map fusion with ground truth landmarks
 *
 * Performance Targets:
 * - Map alignment: <100ms for 10 landmarks
 * - Loop closure detection: <10ms per pose
 * - Bundle adjustment: <500ms for typical map
 */

#ifndef MAP_FITTING_H
#define MAP_FITTING_H

#include "common_types.h"
#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of landmarks in map */
#define MAP_MAX_LANDMARKS 50

/** Maximum trajectory poses for loop closure */
#define MAP_MAX_TRAJECTORY_POSES 1000

/** Maximum constraints for bundle adjustment */
#define MAP_MAX_CONSTRAINTS 500

/** Default loop closure distance threshold (meters) */
#define MAP_LOOP_DISTANCE_THRESHOLD 50.0

/** Default loop closure angle threshold (degrees) */
#define MAP_LOOP_ANGLE_THRESHOLD 30.0

/* ============================================================================
 * Alignment Methods
 * ============================================================================ */

/**
 * @brief Map alignment algorithm selection
 */
typedef enum {
    ALIGN_PROCRUSTES = 0,    /**< Procrustes alignment (similarity transform) */
    ALIGN_HORN,              /**< Horn's absolute orientation (quaternions) */
    ALIGN_UMEYAMA            /**< Umeyama's method */
} alignment_method_t;

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/* position3d_t provided by common_types.h (double precision) */

/**
 * @brief 3x3 rotation matrix
 */
typedef struct {
    double data[9];  /**< Row-major 3x3 matrix */
} rotation_matrix_t;

/**
 * @brief Map alignment transformation
 */
typedef struct {
    rotation_matrix_t rotation;  /**< 3x3 rotation matrix */
    position3d_t translation;    /**< Translation vector */
    double scale;                /**< Scale factor */
    double alignment_error;      /**< RMS alignment error (meters) */
} map_transform_t;

/**
 * @brief Landmark position with ID
 */
typedef struct {
    char id[32];             /**< Landmark identifier */
    position3d_t position;   /**< 3D position */
    double uncertainty[9];   /**< 3x3 covariance matrix (row-major) */
} landmark_t;

/**
 * @brief SLAM map representation
 */
typedef struct {
    landmark_t landmarks[MAP_MAX_LANDMARKS];
    int num_landmarks;
} slam_map_t;

/**
 * @brief Trajectory pose for loop closure
 */
typedef struct {
    double timestamp;        /**< Time (seconds) */
    position3d_t position;   /**< UAV position */
    double heading;          /**< Heading (degrees) */
} trajectory_pose_t;

/**
 * @brief Loop closure detection result
 */
typedef struct {
    int pose_index;          /**< Index of historical pose */
    double distance;         /**< Distance to historical pose (meters) */
} loop_closure_t;

/**
 * @brief Measurement constraint for bundle adjustment
 */
typedef struct {
    char landmark_id[32];    /**< Landmark identifier */
    position3d_t observer_pos; /**< Observer position */
    double azimuth;          /**< Measured azimuth (degrees) */
    double elevation;        /**< Measured elevation (degrees) */
} measurement_constraint_t;

/**
 * @brief Map alignment result
 */
typedef struct {
    bool success;                    /**< Whether alignment succeeded */
    map_transform_t transform;       /**< Transformation parameters */
    alignment_method_t method;       /**< Method used */
    int num_common_landmarks;        /**< Number of matched landmarks */
    char error_message[128];         /**< Error message if failed */
} alignment_result_t;

/* ============================================================================
 * Map Alignment
 * ============================================================================ */

/**
 * @brief Map alignment context
 */
typedef struct {
    map_transform_t transform;       /**< Current transformation */
} map_alignment_t;

/**
 * @brief Initialize map alignment
 *
 * @param alignment Pointer to alignment structure (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_alignment_init(map_alignment_t* alignment);

/**
 * @brief Free map alignment resources
 *
 * @param alignment Alignment to free
 */
void map_alignment_free(map_alignment_t* alignment);

/**
 * @brief Align SLAM map with ground truth positions
 *
 * Finds optimal similarity transformation (rotation, translation, scale)
 * to align SLAM-estimated positions with known ground truth.
 *
 * Performance: <100ms for 10 landmarks
 *
 * @param alignment Alignment instance
 * @param slam_map SLAM-estimated landmark positions
 * @param ground_truth_map Ground truth landmark positions
 * @param method Alignment algorithm to use
 * @param result Output alignment result (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_alignment_align(map_alignment_t* alignment,
                        const slam_map_t* slam_map,
                        const slam_map_t* ground_truth_map,
                        alignment_method_t method,
                        alignment_result_t* result);

/**
 * @brief Procrustes alignment (similarity transform)
 *
 * Finds R, t, s such that: target = s * R * source + t
 * Uses SVD for optimal rotation.
 *
 * @param source_pts Source points (N x 3)
 * @param target_pts Target points (N x 3)
 * @param num_points Number of points
 * @param result Output alignment result
 * @return 0 on success, negative error code on failure
 */
int map_alignment_procrustes(const position3d_t* source_pts,
                              const position3d_t* target_pts,
                              int num_points,
                              alignment_result_t* result);

/**
 * @brief Horn's absolute orientation (quaternion-based)
 *
 * Uses quaternion formulation for rotation estimation.
 *
 * @param source_pts Source points (N x 3)
 * @param target_pts Target points (N x 3)
 * @param num_points Number of points
 * @param result Output alignment result
 * @return 0 on success, negative error code on failure
 */
int map_alignment_horn(const position3d_t* source_pts,
                       const position3d_t* target_pts,
                       int num_points,
                       alignment_result_t* result);

/**
 * @brief Umeyama's method for similarity transformation
 *
 * @param source_pts Source points (N x 3)
 * @param target_pts Target points (N x 3)
 * @param num_points Number of points
 * @param estimate_scale Whether to estimate scale
 * @param result Output alignment result
 * @return 0 on success, negative error code on failure
 */
int map_alignment_umeyama(const position3d_t* source_pts,
                          const position3d_t* target_pts,
                          int num_points,
                          bool estimate_scale,
                          alignment_result_t* result);

/**
 * @brief Transform position from SLAM to ground truth coordinates
 *
 * Applies: transformed = scale * R * slam_position + t
 *
 * @param alignment Alignment instance
 * @param slam_position Position in SLAM frame
 * @param transformed Output position in ground truth frame
 * @return 0 on success, negative error code on failure
 */
int map_alignment_transform_position(const map_alignment_t* alignment,
                                      const position3d_t* slam_position,
                                      position3d_t* transformed);

/**
 * @brief Transform entire SLAM map to ground truth coordinates
 *
 * @param alignment Alignment instance
 * @param slam_map SLAM landmark positions
 * @param transformed_map Output transformed map (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_alignment_transform_map(const map_alignment_t* alignment,
                                 const slam_map_t* slam_map,
                                 slam_map_t* transformed_map);

/* ============================================================================
 * Loop Closure Detection
 * ============================================================================ */

/**
 * @brief Loop closure detector state
 */
typedef struct {
    double distance_threshold;   /**< Distance threshold (meters) */
    double angle_threshold;      /**< Angle threshold (degrees) */
    trajectory_pose_t trajectory[MAP_MAX_TRAJECTORY_POSES];
    int num_poses;               /**< Number of poses in trajectory */
} loop_closure_detector_t;

/**
 * @brief Initialize loop closure detector
 *
 * @param detector Pointer to detector structure (allocated by caller)
 * @param distance_threshold Distance threshold (meters)
 * @param angle_threshold Angle threshold (degrees)
 * @return 0 on success, negative error code on failure
 */
int loop_closure_detector_init(loop_closure_detector_t* detector,
                                double distance_threshold,
                                double angle_threshold);

/**
 * @brief Free loop closure detector resources
 *
 * @param detector Detector to free
 */
void loop_closure_detector_free(loop_closure_detector_t* detector);

/**
 * @brief Add pose to trajectory history
 *
 * @param detector Detector instance
 * @param timestamp Time (seconds)
 * @param position UAV position
 * @param heading UAV heading (degrees)
 * @return 0 on success, negative error code on failure
 */
int loop_closure_add_pose(loop_closure_detector_t* detector,
                          double timestamp,
                          const position3d_t* position,
                          double heading);

/**
 * @brief Detect if current position closes a loop
 *
 * Checks current position against historical poses.
 * Requires minimum time separation (30s) between poses.
 *
 * Performance: <10ms per check
 *
 * @param detector Detector instance
 * @param current_position Current UAV position
 * @param current_heading Current heading (degrees)
 * @param closures Output array of loop closures [MAP_MAX_LANDMARKS]
 * @param num_closures Output number of closures found
 * @return 0 on success, negative error code on failure
 */
int loop_closure_detect(loop_closure_detector_t* detector,
                        const position3d_t* current_position,
                        double current_heading,
                        loop_closure_t* closures,
                        int* num_closures);

/* ============================================================================
 * Map Optimization (Bundle Adjustment)
 * ============================================================================ */

/**
 * @brief Map optimizer context
 */
typedef struct {
    int max_iterations;          /**< Maximum optimization iterations */
    double convergence_threshold; /**< Convergence threshold */
} map_optimizer_t;

/**
 * @brief Initialize map optimizer
 *
 * @param optimizer Pointer to optimizer structure (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_optimizer_init(map_optimizer_t* optimizer);

/**
 * @brief Free map optimizer resources
 *
 * @param optimizer Optimizer to free
 */
void map_optimizer_free(map_optimizer_t* optimizer);

/**
 * @brief Optimize map using least squares bundle adjustment
 *
 * Minimizes reprojection error using Levenberg-Marquardt.
 * Constraints: {landmark_id, observer_pos, azimuth, elevation}
 *
 * Performance: <500ms for typical map
 *
 * @param optimizer Optimizer instance
 * @param slam_map Current SLAM landmark positions
 * @param constraints Measurement constraints
 * @param num_constraints Number of constraints
 * @param optimized_map Output optimized map (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_optimizer_optimize(map_optimizer_t* optimizer,
                           const slam_map_t* slam_map,
                           const measurement_constraint_t* constraints,
                           int num_constraints,
                           slam_map_t* optimized_map);

/* ============================================================================
 * Map Fusion
 * ============================================================================ */

/**
 * @brief Map fusion context
 */
typedef struct {
    map_alignment_t aligner;
    loop_closure_detector_t loop_detector;
    map_optimizer_t optimizer;
} map_fusion_t;

/**
 * @brief Map fusion result
 */
typedef struct {
    bool success;                /**< Whether fusion succeeded */
    slam_map_t fused_map;        /**< Fused landmark positions */
    position3d_t transformed_uav_pos; /**< Transformed UAV position */
    double alignment_error;      /**< Alignment error (meters) */
    map_transform_t transform;   /**< Applied transformation */
} map_fusion_result_t;

/**
 * @brief Initialize map fusion
 *
 * @param fusion Pointer to fusion structure (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_fusion_init(map_fusion_t* fusion);

/**
 * @brief Free map fusion resources
 *
 * @param fusion Fusion to free
 */
void map_fusion_free(map_fusion_t* fusion);

/**
 * @brief Fuse SLAM map with known ground truth landmarks
 *
 * Algorithm:
 * 1. Align SLAM map with ground truth
 * 2. Transform SLAM positions to ground truth frame
 * 3. Fuse using weighted average based on uncertainty
 *
 * @param fusion Fusion instance
 * @param slam_map SLAM landmarks with uncertainties
 * @param ground_truth_map Known landmarks
 * @param uav_position Current UAV position in SLAM frame
 * @param result Output fusion result (allocated by caller)
 * @return 0 on success, negative error code on failure
 */
int map_fusion_fuse(map_fusion_t* fusion,
                    const slam_map_t* slam_map,
                    const slam_map_t* ground_truth_map,
                    const position3d_t* uav_position,
                    map_fusion_result_t* result);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Calculate RMS error between point sets
 *
 * @param points1 First point set
 * @param points2 Second point set
 * @param num_points Number of points
 * @return RMS error (meters)
 */
double map_calculate_rms_error(const position3d_t* points1,
                                const position3d_t* points2,
                                int num_points);

/**
 * @brief Normalize angle to [-180, 180] range
 *
 * @param angle Angle in degrees
 * @return Normalized angle (degrees)
 */
double map_normalize_angle(double angle);

/**
 * @brief Apply rotation matrix to position
 *
 * @param rotation Rotation matrix
 * @param position Input position
 * @param rotated Output rotated position
 */
void map_apply_rotation(const rotation_matrix_t* rotation,
                        const position3d_t* position,
                        position3d_t* rotated);

/**
 * @brief Compute rotation matrix from quaternion
 *
 * @param qw Quaternion w component
 * @param qx Quaternion x component
 * @param qy Quaternion y component
 * @param qz Quaternion z component
 * @param rotation Output rotation matrix
 */
void map_quaternion_to_rotation(double qw, double qx, double qy, double qz,
                                 rotation_matrix_t* rotation);

/**
 * @brief Find common landmarks between two maps
 *
 * @param map1 First map
 * @param map2 Second map
 * @param common_ids Output array of common IDs [MAP_MAX_LANDMARKS]
 * @param num_common Output number of common landmarks
 * @return 0 on success, negative error code on failure
 */
int map_find_common_landmarks(const slam_map_t* map1,
                               const slam_map_t* map2,
                               char common_ids[][32],
                               int* num_common);

#endif /* MAP_FITTING_H */
