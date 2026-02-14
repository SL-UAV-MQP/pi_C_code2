/**
 * @file cfar_2d.h
 * @brief 2D Order-Statistic CFAR Detector API
 *
 * Implements adaptive thresholding for MUSIC AOA validation in
 * multipath environments using OS-CFAR (Order-Statistic CFAR).
 */

#ifndef CFAR_2D_H
#define CFAR_2D_H

#include "common_types.h"
#include <stdbool.h>

/* ============================================================================
 * CFAR Configuration Structure
 * ============================================================================ */

/**
 * @brief CFAR detector configuration
 */
typedef struct {
    int window_size_az;          /**< Window size in azimuth direction */
    int window_size_el;          /**< Window size in elevation direction */
    int guard_size_az;           /**< Guard cells in azimuth */
    int guard_size_el;           /**< Guard cells in elevation */
    double k_factor;             /**< Order statistic selection factor [0-1] */
    double pfa;                  /**< Probability of false alarm */
    double threshold_offset_db;  /**< Additional threshold margin (dB) */
    double min_snr_db;           /**< Minimum SNR for detection (dB) */
    bool use_circular_padding;   /**< Use circular padding for azimuth */
} cfar_config_t;

/**
 * @brief CFAR training mask (precomputed)
 */
typedef struct {
    bool* mask;                  /**< Binary mask [window_height × window_width] */
    int window_height;           /**< Window height (azimuth) */
    int window_width;            /**< Window width (elevation) */
    int n_training_cells;        /**< Number of training cells (sum of mask) */
} cfar_training_mask_t;

/* ============================================================================
 * CFAR Configuration Presets
 * ============================================================================ */

/**
 * @brief Get CMRCM-optimized CFAR configuration
 *
 * Returns validated configuration for CMRCM field environment.
 *
 * @param[out] config  Output configuration structure
 * @return MUSIC_SUCCESS on success
 */
music_status_t cfar_get_cmrcm_config(cfar_config_t* config);

/**
 * @brief Get default CFAR configuration
 *
 * @param[out] config  Output configuration structure
 * @return MUSIC_SUCCESS on success
 */
music_status_t cfar_get_default_config(cfar_config_t* config);

/* ============================================================================
 * CFAR Initialization and Cleanup
 * ============================================================================ */

/**
 * @brief Initialize CFAR training mask
 *
 * Creates binary mask for training cells (excludes CUT and guard cells).
 *
 * Mask layout (X = training cell, G = guard, C = CUT):
 *   X X X X X X X
 *   X X G G G X X
 *   X X G C G X X
 *   X X G G G X X
 *   X X X X X X X
 *
 * @param[in] config  CFAR configuration
 * @param[out] mask  Output training mask (allocated by function)
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 */
music_status_t cfar_create_training_mask(
    const cfar_config_t* config,
    cfar_training_mask_t* mask
);

/**
 * @brief Free CFAR training mask
 *
 * @param[in,out] mask  Mask to free
 */
void cfar_free_training_mask(cfar_training_mask_t* mask);

/**
 * @brief Validate CFAR configuration
 *
 * Checks configuration parameters for validity.
 *
 * Validation rules:
 * - Window sizes must be odd
 * - Guard sizes must be < window_size/2
 * - k_factor must be in (0, 1)
 * - pfa must be positive
 *
 * @param[in] config  Configuration to validate
 * @return MUSIC_SUCCESS if valid, error code otherwise
 */
music_status_t cfar_validate_config(const cfar_config_t* config);

/* ============================================================================
 * Main CFAR Detection
 * ============================================================================ */

/**
 * @brief Apply 2D OS-CFAR detection to spatial spectrum
 *
 * Performs adaptive thresholding using Order-Statistic CFAR.
 * Handles multipath environment with inhomogeneous noise.
 *
 * Algorithm:
 * 1. For each cell (i,j):
 *    a. Extract window around cell
 *    b. Select training cells (exclude CUT and guard)
 *    c. Sort training cells, pick k-th order statistic
 *    d. Compute threshold = k-th value * scale_factor
 *    e. Detect if cell > threshold AND SNR > min_snr
 *
 * OS-CFAR advantages:
 * - Robust to outliers in training cells
 * - Better than CA-CFAR in multipath
 * - Preserves power density information
 *
 * @param[in] spectrum  2D spectrum in dB [height × width]
 * @param[in] height  Number of rows (azimuth bins)
 * @param[in] width  Number of columns (elevation bins)
 * @param[in] config  CFAR configuration
 * @param[in] mask  Precomputed training mask
 * @param[out] result  Detection result (allocated by caller)
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Input spectrum should be in dB
 * @note Output threshold and SNR are in dB
 * @note Use circular padding for azimuth (wraps at ±180°)
 */
music_status_t cfar_detect_2d(
    const double* spectrum,
    int height,
    int width,
    const cfar_config_t* config,
    const cfar_training_mask_t* mask,
    cfar_result_t* result
);

/* ============================================================================
 * CFAR Threshold Computation
 * ============================================================================ */

/**
 * @brief Compute adaptive CFAR threshold map
 *
 * Sliding window computation of OS-CFAR threshold for entire spectrum.
 *
 * @param[in] spectrum  2D spectrum in linear power [height × width]
 * @param[in] height  Number of rows
 * @param[in] width  Number of columns
 * @param[in] config  CFAR configuration
 * @param[in] mask  Training mask
 * @param[out] threshold  Output threshold map in linear power [height × width]
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Input and output are in linear power (not dB)
 * @note Uses k-th order statistic from sorted training cells
 */
music_status_t cfar_compute_threshold(
    const double* spectrum,
    int height,
    int width,
    const cfar_config_t* config,
    const cfar_training_mask_t* mask,
    double* threshold
);

/* ============================================================================
 * CFAR Scale Factor
 * ============================================================================ */

/**
 * @brief Compute CFAR scale factor from Pfa
 *
 * Calculates multiplicative scale factor based on probability of false alarm.
 *
 * Approximation:
 *   α = N_train * (Pfa^(-1/N_train) - 1)
 *
 * Total scale = α * 10^(threshold_offset_db / 10)
 *
 * @param[in] n_training_cells  Number of training cells
 * @param[in] pfa  Probability of false alarm
 * @param[in] threshold_offset_db  Additional offset in dB
 * @param[out] scale_factor  Output scale factor (linear)
 *
 * @return MUSIC_SUCCESS on success
 */
music_status_t cfar_compute_scale_factor(
    int n_training_cells,
    double pfa,
    double threshold_offset_db,
    double* scale_factor
);

/* ============================================================================
 * Noise Floor Estimation
 * ============================================================================ */

/**
 * @brief Estimate noise floor from spectrum using robust statistics
 *
 * Uses lower percentile to avoid signal contamination.
 *
 * @param[in] spectrum_db  Spectrum in dB [height × width]
 * @param[in] size  Total number of elements (height * width)
 * @param[out] noise_floor_db  Estimated noise floor in dB
 *
 * @return MUSIC_SUCCESS on success
 *
 * @note Uses 10th percentile to avoid signal contamination
 */
music_status_t cfar_estimate_noise_floor(
    const double* spectrum_db,
    int size,
    double* noise_floor_db
);

/* ============================================================================
 * Peak Extraction from CFAR Detections
 * ============================================================================ */

/**
 * @brief Extract peak locations from CFAR detection mask
 *
 * Converts binary detection mask to list of peak coordinates and properties.
 *
 * @param[in] spectrum  Original spectrum in dB [height × width]
 * @param[in] detection_mask  Binary detection mask [height × width]
 * @param[in] height  Number of rows
 * @param[in] width  Number of columns
 * @param[in] azimuth_grid  Azimuth angles (degrees) [height]
 * @param[in] elevation_grid  Elevation angles (optional) [width]
 * @param[in] noise_floor_db  Noise floor in dB
 * @param[out] peaks  Output detected sources (allocated by caller)
 *
 * @return MUSIC_SUCCESS on success
 *
 * @note Peaks sorted by power (descending)
 * @note SNR = peak_power - noise_floor
 */
music_status_t cfar_extract_peaks(
    const double* spectrum,
    const bool* detection_mask,
    int height,
    int width,
    const double* azimuth_grid,
    const double* elevation_grid,
    double noise_floor_db,
    detected_sources_t* peaks
);

/* ============================================================================
 * Non-Maximum Suppression
 * ============================================================================ */

/**
 * @brief Apply 2D Non-Maximum Suppression
 *
 * Eliminates duplicate detections by keeping only local maxima.
 *
 * @param[in] detections  Binary detection mask [height × width]
 * @param[in] spectrum  Original spectrum [height × width]
 * @param[in] height  Number of rows
 * @param[in] width  Number of columns
 * @param[in] window_size  NMS window size (default 5)
 * @param[out] suppressed  Output mask after NMS [height × width]
 *
 * @return MUSIC_SUCCESS on success
 *
 * @note Keeps only peaks that are local maxima in window_size neighborhood
 */
music_status_t cfar_apply_nms_2d(
    const bool* detections,
    const double* spectrum,
    int height,
    int width,
    int window_size,
    bool* suppressed
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Sort double array and return indices
 *
 * Helper for order statistic selection.
 *
 * @param[in] data  Input array
 * @param[in] length  Array length
 * @param[out] sorted_data  Sorted array (allocated by caller)
 * @param[out] indices  Original indices (allocated by caller, optional)
 *
 * @return MUSIC_SUCCESS on success
 */
music_status_t sort_with_indices(
    const double* data,
    int length,
    double* sorted_data,
    int* indices
);

/**
 * @brief Compute percentile of array
 *
 * @param[in] data  Input array (need not be sorted)
 * @param[in] length  Array length
 * @param[in] percentile  Percentile value [0-100]
 * @param[out] result  Output percentile value
 *
 * @return MUSIC_SUCCESS on success
 */
music_status_t compute_percentile(
    const double* data,
    int length,
    double percentile,
    double* result
);

/**
 * @brief Print CFAR configuration
 *
 * @param[in] config  Configuration to print
 */
void cfar_print_config(const cfar_config_t* config);

/**
 * @brief Print CFAR detection results
 *
 * @param[in] result  Detection result to print
 */
void cfar_print_results(const cfar_result_t* result);

#endif /* CFAR_2D_H */
