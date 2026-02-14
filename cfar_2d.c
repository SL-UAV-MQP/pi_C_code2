/**
 * @file cfar_2d.c
 * @brief 2D Order-Statistic CFAR Detector Implementation
 */

#include "cfar_2d.h"
#include "music_config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Configuration Presets
 * ============================================================================ */

music_status_t cfar_get_cmrcm_config(cfar_config_t* config) {
    // cmrcm_config static method
    if (!config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    config->window_size_az = CFAR_WINDOW_SIZE_AZ;
    config->window_size_el = CFAR_WINDOW_SIZE_EL;
    config->guard_size_az = CFAR_GUARD_SIZE_AZ;
    config->guard_size_el = CFAR_GUARD_SIZE_EL;
    config->k_factor = CFAR_K_FACTOR;
    config->pfa = CFAR_PFA;
    config->threshold_offset_db = CFAR_THRESHOLD_OFFSET_DB;
    config->min_snr_db = CFAR_MIN_SNR_DB;
    config->use_circular_padding = true;

    return MUSIC_SUCCESS;
}

music_status_t cfar_get_default_config(cfar_config_t* config) {
    if (!config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    config->window_size_az = 9;
    config->window_size_el = 9;
    config->guard_size_az = 3;
    config->guard_size_el = 3;
    config->k_factor = 0.75;
    config->pfa = 1e-4;
    config->threshold_offset_db = 3.0;
    config->min_snr_db = 10.0;
    config->use_circular_padding = false;

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Configuration Validation
 * ============================================================================ */

music_status_t cfar_validate_config(const cfar_config_t* config) {
    // TODO: Implement configuration validation

    if (!config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Check odd window sizes
    if (config->window_size_az % 2 == 0 || config->window_size_el % 2 == 0) {
        fprintf(stderr, "ERROR: Window sizes must be odd\n");
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Check guard sizes
    if (config->guard_size_az >= config->window_size_az / 2 ||
        config->guard_size_el >= config->window_size_el / 2) {
        fprintf(stderr, "ERROR: Guard sizes must be < window_size/2\n");
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Check k_factor
    if (config->k_factor <= 0.0 || config->k_factor >= 1.0) {
        fprintf(stderr, "ERROR: k_factor must be in (0, 1)\n");
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Check pfa
    if (config->pfa <= 0.0 || config->pfa >= 1.0) {
        fprintf(stderr, "ERROR: pfa must be in (0, 1)\n");
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Training Mask Creation
 * ============================================================================ */

music_status_t cfar_create_training_mask(
    const cfar_config_t* config,
    cfar_training_mask_t* mask)
{
    // TODO: Implement training mask creation

    if (!config || !mask) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int h = config->window_size_az;
    int w = config->window_size_el;
    int gh = config->guard_size_az;
    int gw = config->guard_size_el;

    mask->window_height = h;
    mask->window_width = w;
    mask->mask = (bool*)malloc(h * w * sizeof(bool));
    if (!mask->mask) {
        return MUSIC_ERROR_MEMORY;
    }

    // Initialize all to true (training cells)
    for (int i = 0; i < h * w; i++) {
        mask->mask[i] = true;
    }

    // Zero out CUT and guard cells
    int ch = h / 2;  // Center row
    int cw = w / 2;  // Center column

    // FIX: Use row-major indexing consistently (i * w + j)
    // Previously used column-major (i + j * h) here but row-major elsewhere
    for (int i = ch - gh; i <= ch + gh; i++) {
        for (int j = cw - gw; j <= cw + gw; j++) {
            if (i >= 0 && i < h && j >= 0 && j < w) {
                mask->mask[i * w + j] = false;
            }
        }
    }

    // Count training cells
    int count = 0;
    for (int i = 0; i < h * w; i++) {
        if (mask->mask[i]) count++;
    }
    mask->n_training_cells = count;

    if (count < 10) {
        fprintf(stderr, "WARNING: Only %d training cells (consider larger window)\n", count);
    }

    return MUSIC_SUCCESS;
}

void cfar_free_training_mask(cfar_training_mask_t* mask) {
    if (mask && mask->mask) {
        free(mask->mask);
        mask->mask = NULL;
    }
}

/* ============================================================================
 * CFAR Scale Factor Computation
 * ============================================================================ */

music_status_t cfar_compute_scale_factor(
    int n_training_cells,
    double pfa,
    double threshold_offset_db,
    double* scale_factor)
{
    // TODO: Implement scale factor computation

    if (!scale_factor || n_training_cells <= 0) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // NOTE: This is the CA-CFAR (Cell Averaging) scale factor formula.
    // For OS-CFAR (used in detect functions with k-th order statistic),
    // the exact formula differs but this approximation is reasonable
    // when applied as dB offset: threshold = X_k + 10*log10(alpha)
    double alpha = n_training_cells * (pow(pfa, -1.0 / n_training_cells) - 1.0);

    // Apply offset
    double offset_linear = pow(10.0, threshold_offset_db / 10.0);

    *scale_factor = alpha * offset_linear;

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Noise Floor Estimation
 * ============================================================================ */

music_status_t cfar_estimate_noise_floor(
    const double* spectrum_db,
    int size,
    double* noise_floor_db)
{
    // TODO: Implement noise floor estimation using percentile

    if (!spectrum_db || !noise_floor_db || size <= 0) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Copy and sort spectrum
    double* sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) {
        return MUSIC_ERROR_MEMORY;
    }
    memcpy(sorted, spectrum_db, size * sizeof(double));

    // Simple bubble sort (for small arrays)
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (sorted[j] > sorted[j + 1]) {
                double temp = sorted[j];
                sorted[j] = sorted[j + 1];
                sorted[j + 1] = temp;
            }
        }
    }

    // 10th percentile
    int idx = (int)(0.10 * size);
    if (idx >= size) idx = size - 1;
    *noise_floor_db = sorted[idx];

    free(sorted);
    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Main CFAR Detection
 * ============================================================================ */

music_status_t cfar_detect_2d(
    const double* spectrum,
    int height,
    int width,
    const cfar_config_t* config,
    const cfar_training_mask_t* mask,
    cfar_result_t* result)
{

    if (!spectrum || !config || !mask || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (result->height != height || result->width != width) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Validate configuration
    music_status_t status = cfar_validate_config(config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Estimate noise floor
    status = cfar_estimate_noise_floor(spectrum, height * width, &result->noise_floor_db);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Compute scale factor
    double scale_factor;
    status = cfar_compute_scale_factor(
        mask->n_training_cells,
        config->pfa,
        config->threshold_offset_db,
        &scale_factor
    );
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Initialize result arrays
    memset(result->detections, 0, height * width * sizeof(bool));
    result->n_detections = 0;
    result->peak_snr_db = -1000.0;

    // Window half-sizes
    int half_h = config->window_size_az / 2;
    int half_w = config->window_size_el / 2;

    // Allocate temporary buffer for training cells
    double* training_cells = (double*)malloc(mask->n_training_cells * sizeof(double));
    if (!training_cells) {
        return MUSIC_ERROR_MEMORY;
    }

    // main detection loop
    // Sliding window over entire spectrum
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            int cell_idx = row * width + col;
            double cell_value = spectrum[cell_idx];

            // Extract training cells from window
            int n_training = 0;
            for (int i = -half_h; i <= half_h; i++) {
                for (int j = -half_w; j <= half_w; j++) {
                    // Compute window position
                    int win_row = row + i;
                    int win_col = col + j;

                    // Handle boundary conditions
                    if (config->use_circular_padding) {
                        // Circular padding in azimuth (wraps at ±180°)
                        win_row = (win_row + height) % height;
                    } else {
                        if (win_row < 0 || win_row >= height) continue;
                    }
                    if (win_col < 0 || win_col >= width) continue;

                    // Check if this is a training cell (using mask)
                    int mask_row = i + half_h;
                    int mask_col = j + half_w;
                    int mask_idx = mask_row * mask->window_width + mask_col;

                    if (mask->mask[mask_idx]) {
                        // Add to training cells
                        int spec_idx = win_row * width + win_col;
                        training_cells[n_training++] = spectrum[spec_idx];
                    }
                }
            }

            // Ensure we have enough training cells
            if (n_training < 5) {
                result->threshold_db[cell_idx] = cell_value + config->threshold_offset_db;
                result->snr_db[cell_idx] = 0.0;
                continue;
            }

            // order statistic selection
            // Sort training cells
            for (int i = 0; i < n_training - 1; i++) {
                for (int j = 0; j < n_training - i - 1; j++) {
                    if (training_cells[j] > training_cells[j + 1]) {
                        double temp = training_cells[j];
                        training_cells[j] = training_cells[j + 1];
                        training_cells[j + 1] = temp;
                    }
                }
            }

            // Select k-th order statistic
            int k_idx = (int)(config->k_factor * (n_training - 1));
            if (k_idx >= n_training) k_idx = n_training - 1;
            double noise_estimate_db = training_cells[k_idx];

            // Compute adaptive threshold
            // threshold = noise_estimate + 10*log10(scale_factor)
            double threshold_db = noise_estimate_db + 10.0 * log10(scale_factor);
            result->threshold_db[cell_idx] = threshold_db;

            // Compute SNR
            double snr_db = cell_value - noise_estimate_db;
            result->snr_db[cell_idx] = snr_db;

            // Apply detection criteria
            if (cell_value > threshold_db && snr_db > config->min_snr_db) {
                result->detections[cell_idx] = true;
                result->n_detections++;

                if (snr_db > result->peak_snr_db) {
                    result->peak_snr_db = snr_db;
                }
            }
        }
    }

    free(training_cells);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * CFAR Threshold Computation
 * ============================================================================ */

music_status_t cfar_compute_threshold(
    const double* spectrum,
    int height,
    int width,
    const cfar_config_t* config,
    const cfar_training_mask_t* mask,
    double* threshold)
{
    // compute_threshold

    if (!spectrum || !config || !mask || !threshold) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute scale factor
    double scale_factor;
    music_status_t status = cfar_compute_scale_factor(
        mask->n_training_cells,
        config->pfa,
        config->threshold_offset_db,
        &scale_factor
    );
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Window half-sizes
    int half_h = config->window_size_az / 2;
    int half_w = config->window_size_el / 2;

    // Pad spectrum for edge handling
    int padded_height = height + 2 * half_h;
    int padded_width = width + 2 * half_w;
    double* spectrum_padded = (double*)malloc(padded_height * padded_width * sizeof(double));
    if (!spectrum_padded) {
        return MUSIC_ERROR_MEMORY;
    }

    // Apply padding (circular for azimuth, constant for elevation)
    if (config->use_circular_padding) {
        // Circular padding in azimuth (vertical)
        for (int i = 0; i < padded_height; i++) {
            for (int j = 0; j < padded_width; j++) {
                int src_row = (i - half_h + height) % height;
                int src_col = j - half_w;

                // Symmetric padding in elevation (horizontal)
                if (src_col < 0) {
                    src_col = -src_col - 1;
                } else if (src_col >= width) {
                    src_col = 2 * width - src_col - 1;
                }

                if (src_row >= 0 && src_row < height && src_col >= 0 && src_col < width) {
                    spectrum_padded[i * padded_width + j] = spectrum[src_row * width + src_col];
                } else {
                    // Fallback to median value
                    spectrum_padded[i * padded_width + j] = 0.0;
                }
            }
        }
    } else {
        // Constant padding with median value
        // Compute median of original spectrum
        double* temp = (double*)malloc(height * width * sizeof(double));
        if (!temp) {
            free(spectrum_padded);
            return MUSIC_ERROR_MEMORY;
        }
        memcpy(temp, spectrum, height * width * sizeof(double));

        // Simple bubble sort for median
        for (int i = 0; i < height * width - 1; i++) {
            for (int j = 0; j < height * width - i - 1; j++) {
                if (temp[j] > temp[j + 1]) {
                    double t = temp[j];
                    temp[j] = temp[j + 1];
                    temp[j + 1] = t;
                }
            }
        }
        double median_val = temp[(height * width) / 2];
        free(temp);

        // Fill padded array
        for (int i = 0; i < padded_height; i++) {
            for (int j = 0; j < padded_width; j++) {
                int src_row = i - half_h;
                int src_col = j - half_w;

                if (src_row >= 0 && src_row < height && src_col >= 0 && src_col < width) {
                    spectrum_padded[i * padded_width + j] = spectrum[src_row * width + src_col];
                } else {
                    spectrum_padded[i * padded_width + j] = median_val;
                }
            }
        }
    }

    // Allocate temporary buffer for training cells
    double* training_cells = (double*)malloc(mask->n_training_cells * sizeof(double));
    if (!training_cells) {
        free(spectrum_padded);
        return MUSIC_ERROR_MEMORY;
    }

    // sliding window loop
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Extract window from padded spectrum
            int n_training = 0;
            for (int wi = 0; wi < config->window_size_az; wi++) {
                for (int wj = 0; wj < config->window_size_el; wj++) {
                    // Check if this is a training cell
                    int mask_idx = wi * mask->window_width + wj;
                    if (mask->mask[mask_idx]) {
                        // Get value from padded spectrum
                        int padded_row = i + wi;
                        int padded_col = j + wj;
                        int padded_idx = padded_row * padded_width + padded_col;
                        training_cells[n_training++] = spectrum_padded[padded_idx];
                    }
                }
            }

            // Sort training cells
            for (int k = 0; k < n_training - 1; k++) {
                for (int m = 0; m < n_training - k - 1; m++) {
                    if (training_cells[m] > training_cells[m + 1]) {
                        double temp = training_cells[m];
                        training_cells[m] = training_cells[m + 1];
                        training_cells[m + 1] = temp;
                    }
                }
            }

            // Select k-th order statistic
            int k_index = (int)(config->k_factor * (n_training - 1));
            if (k_index >= n_training) k_index = n_training - 1;
            double threshold_val = training_cells[k_index];

            // Apply scaling factor
            threshold[i * width + j] = threshold_val * scale_factor;
        }
    }

    free(training_cells);
    free(spectrum_padded);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Peak Extraction from CFAR Detections
 * ============================================================================ */

music_status_t cfar_extract_peaks(
    const double* spectrum,
    const bool* detection_mask,
    int height,
    int width,
    const double* azimuth_grid,
    const double* elevation_grid,
    double noise_floor_db,
    detected_sources_t* peaks)
{
    // extract_peaks

    if (!spectrum || !detection_mask || !azimuth_grid || !peaks) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Find peak indices
    int* peak_rows = (int*)malloc(height * width * sizeof(int));
    int* peak_cols = (int*)malloc(height * width * sizeof(int));
    double* peak_powers = (double*)malloc(height * width * sizeof(double));

    if (!peak_rows || !peak_cols || !peak_powers) {
        free(peak_rows);
        free(peak_cols);
        free(peak_powers);
        return MUSIC_ERROR_MEMORY;
    }

    int n_peaks = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int idx = i * width + j;
            if (detection_mask[idx]) {
                peak_rows[n_peaks] = i;
                peak_cols[n_peaks] = j;
                peak_powers[n_peaks] = spectrum[idx];
                n_peaks++;
            }
        }
    }

    if (n_peaks == 0) {
        peaks->num_sources = 0;
        free(peak_rows);
        free(peak_cols);
        free(peak_powers);
        return MUSIC_SUCCESS;
    }

    // Sort by power (descending)
    for (int i = 0; i < n_peaks - 1; i++) {
        for (int j = 0; j < n_peaks - i - 1; j++) {
            if (peak_powers[j] < peak_powers[j + 1]) {
                // Swap powers
                double temp_pow = peak_powers[j];
                peak_powers[j] = peak_powers[j + 1];
                peak_powers[j + 1] = temp_pow;

                // Swap rows
                int temp_row = peak_rows[j];
                peak_rows[j] = peak_rows[j + 1];
                peak_rows[j + 1] = temp_row;

                // Swap cols
                int temp_col = peak_cols[j];
                peak_cols[j] = peak_cols[j + 1];
                peak_cols[j + 1] = temp_col;
            }
        }
    }

    // Extract top peaks (up to MAX_SOURCES)
    int num_output = (n_peaks < MAX_SOURCES) ? n_peaks : MAX_SOURCES;
    peaks->num_sources = num_output;

    for (int i = 0; i < num_output; i++) {
        int row = peak_rows[i];
        int col = peak_cols[i];

        peaks->sources[i].azimuth_deg = azimuth_grid[row];

        if (elevation_grid != NULL) {
            peaks->sources[i].elevation_deg = elevation_grid[col];
        } else {
            peaks->sources[i].elevation_deg = 0.0;
        }

        peaks->sources[i].magnitude_db = peak_powers[i];

        double snr = peak_powers[i] - noise_floor_db;
        peaks->sources[i].confidence = (snr > 0.0) ? (snr / 20.0) : 0.0;  // Normalize
        if (peaks->sources[i].confidence > 1.0) peaks->sources[i].confidence = 1.0;

        strcpy(peaks->sources[i].type, "CFAR_peak");
    }

    free(peak_rows);
    free(peak_cols);
    free(peak_powers);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Non-Maximum Suppression
 * ============================================================================ */

music_status_t cfar_apply_nms_2d(
    const bool* detections,
    const double* spectrum,
    int height,
    int width,
    int window_size,
    bool* suppressed)
{
    // apply_nms_2d static method

    if (!detections || !spectrum || !suppressed) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (window_size % 2 == 0) {
        fprintf(stderr, "WARNING: NMS window_size should be odd, adjusting %d -> %d\n",
                window_size, window_size + 1);
        window_size++;
    }

    int half_win = window_size / 2;

    // Find local maxima
    bool* local_maxima = (bool*)malloc(height * width * sizeof(bool));
    if (!local_maxima) {
        return MUSIC_ERROR_MEMORY;
    }

    // For each cell, check if it's a local maximum
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int idx = i * width + j;
            double cell_val = spectrum[idx];
            bool is_max = true;

            // Check neighborhood
            for (int di = -half_win; di <= half_win && is_max; di++) {
                for (int dj = -half_win; dj <= half_win && is_max; dj++) {
                    int ni = i + di;
                    int nj = j + dj;

                    // Circular boundary for azimuth
                    ni = (ni + height) % height;

                    // Skip out of bounds in elevation
                    if (nj < 0 || nj >= width) continue;

                    int n_idx = ni * width + nj;

                    // If any neighbor is strictly greater, not a local max
                    if (spectrum[n_idx] > cell_val) {
                        is_max = false;
                    }
                }
            }

            local_maxima[idx] = is_max;
        }
    }

    // Keep only detections that are also local maxima
    for (int i = 0; i < height * width; i++) {
        suppressed[i] = detections[i] && local_maxima[i];
    }

    free(local_maxima);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

void cfar_print_config(const cfar_config_t* config) {
    if (!config) {
        printf("ERROR: NULL config\n");
        return;
    }

    printf("=== CFAR Configuration ===\n");
    printf("Window size: [%d, %d]\n", config->window_size_az, config->window_size_el);
    printf("Guard size: [%d, %d]\n", config->guard_size_az, config->guard_size_el);
    printf("k-factor: %.2f\n", config->k_factor);
    printf("Pfa: %.2e\n", config->pfa);
    printf("Threshold offset: %.1f dB\n", config->threshold_offset_db);
    printf("Min SNR: %.1f dB\n", config->min_snr_db);
    printf("Circular padding: %s\n", config->use_circular_padding ? "ON" : "OFF");
    printf("\n");
}

void cfar_print_results(const cfar_result_t* result) {
    if (!result) {
        printf("ERROR: NULL result\n");
        return;
    }

    printf("=== CFAR Detection Results ===\n");
    printf("Dimensions: %d × %d\n", result->height, result->width);
    printf("Detections: %d\n", result->n_detections);
    printf("Peak SNR: %.1f dB\n", result->peak_snr_db);
    printf("Noise floor: %.1f dB\n", result->noise_floor_db);
    printf("\n");
}
