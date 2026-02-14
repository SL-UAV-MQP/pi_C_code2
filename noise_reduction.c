/**
 * @file noise_reduction.c
 * @brief Advanced Noise Reduction Algorithms
 *
 * Implements spectral and spatial noise reduction for urban RF environment.
 *
 * Performance targets: Single-channel <50ms, SNR improvement 5-15 dB
 */

#include "noise_reduction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SPEED_OF_LIGHT 299792458.0

/** Default UCA radius for spatial canceller (meters) */
#define NR_UCA_RADIUS 0.176

/** Diagonal loading factor for LCMV when no covariance data available */
#define NR_LCMV_DIAGONAL_LOADING 1e-2

/* ============================================================================
 * Configuration Helpers
 * ============================================================================ */

void noise_reduction_config_default(noise_reduction_config_t* config) {
    config->alpha = NR_ALPHA;
    config->beta = NR_SPECTRAL_FLOOR;
    config->noise_estimation_frames = 10;
    config->wiener_gain_min = NR_WIENER_GAIN_MIN;
    config->num_elements = 6;
    config->null_depth = 30.0;
    config->adaptation_rate = 0.1;
    config->noise_percentile = 10.0;
    config->fft_size = NR_FFT_SIZE;
    config->overlap = NR_OVERLAP;
    config->sample_rate = 61.44e6;
}

/* ============================================================================
 * Noise Floor Estimator
 * ============================================================================ */

music_status_t noise_floor_estimator_init(
    noise_floor_estimator_t* estimator,
    const noise_reduction_config_t* config
) {
    if (!estimator || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    estimator->config = *config;
    estimator->num_bins = config->fft_size / 2;
    estimator->initialized = false;

    estimator->noise_floor = (double*)calloc(estimator->num_bins, sizeof(double));
    if (!estimator->noise_floor) {
        return MUSIC_ERROR_MEMORY;
    }

    return MUSIC_SUCCESS;
}

void noise_floor_estimator_free(noise_floor_estimator_t* estimator) {
    if (estimator && estimator->noise_floor) {
        free(estimator->noise_floor);
        estimator->noise_floor = NULL;
    }
}

music_status_t noise_floor_estimate(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* noise_floor
) {
    // estimate_noise_floor
    if (!estimator || !spectrum || !noise_floor) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Use percentile-based estimation
    double* sorted = (double*)malloc(num_bins * sizeof(double));
    if (!sorted) {
        return MUSIC_ERROR_MEMORY;
    }

    memcpy(sorted, spectrum, num_bins * sizeof(double));

    // Sort for percentile calculation
    for (int i = 1; i < num_bins; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }

    int percentile_idx = (int)(num_bins * estimator->config.noise_percentile / 100.0);
    double percentile_value = sorted[percentile_idx];

    // Fill noise floor with minimum of spectrum and percentile
    for (int i = 0; i < num_bins; i++) {
        noise_floor[i] = (spectrum[i] < percentile_value) ? spectrum[i] : percentile_value;
    }

    free(sorted);
    return MUSIC_SUCCESS;
}

music_status_t noise_floor_update(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins
) {
    if (!estimator || !spectrum) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    double* current_estimate = (double*)malloc(num_bins * sizeof(double));
    if (!current_estimate) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = noise_floor_estimate(estimator, spectrum, num_bins,
                                                  current_estimate);
    if (status != MUSIC_SUCCESS) {
        free(current_estimate);
        return status;
    }

    if (!estimator->initialized) {
        memcpy(estimator->noise_floor, current_estimate, num_bins * sizeof(double));
        estimator->initialized = true;
    } else {
        // Exponential moving average
        double alpha = estimator->config.adaptation_rate;
        for (int i = 0; i < num_bins; i++) {
            estimator->noise_floor[i] = (1.0 - alpha) * estimator->noise_floor[i] +
                                        alpha * current_estimate[i];
        }
    }

    free(current_estimate);
    return MUSIC_SUCCESS;
}

music_status_t noise_floor_get_snr(
    const noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* snr_db
) {
    // get_snr_estimate
    if (!estimator || !spectrum || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (!estimator->initialized) {
        *snr_db = 0.0;
        return MUSIC_SUCCESS;
    }

    double signal_power = 0.0;
    double noise_power = 0.0;

    for (int i = 0; i < num_bins; i++) {
        signal_power += spectrum[i];
        noise_power += estimator->noise_floor[i];
    }

    signal_power /= num_bins;
    noise_power /= num_bins;

    *snr_db = 10.0 * log10((signal_power / (noise_power + 1e-12)) + 1e-12);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Spectral Subtraction
 * ============================================================================ */

music_status_t spectral_subtractor_init(
    spectral_subtractor_t* subtractor,
    const noise_reduction_config_t* config
) {
    if (!subtractor || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    subtractor->config = *config;

    music_status_t status = noise_floor_estimator_init(&subtractor->noise_estimator, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    // Create Hann window
    subtractor->window = (double*)malloc(config->fft_size * sizeof(double));
    if (!subtractor->window) {
        noise_floor_estimator_free(&subtractor->noise_estimator);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < config->fft_size; i++) {
        subtractor->window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (config->fft_size - 1)));
    }

    return MUSIC_SUCCESS;
}

void spectral_subtractor_free(spectral_subtractor_t* subtractor) {
    if (subtractor) {
        noise_floor_estimator_free(&subtractor->noise_estimator);
        if (subtractor->window) {
            free(subtractor->window);
        }
    }
}

music_status_t spectral_subtractor_process(
    spectral_subtractor_t* subtractor,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    if (!subtractor || !iq_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = subtractor->config.fft_size;
    int hop = (int)(N * (1.0 - subtractor->config.overlap));
    int num_frames = (num_samples - N) / hop + 1;
    double alpha = subtractor->config.alpha;
    double beta = subtractor->config.beta;

    if (num_frames < 1 || N <= 0) {
        memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));
        return MUSIC_SUCCESS;
    }

    // Allocate FFTW buffers
    fftw_complex* fft_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* ifft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    double* overlap_buf = (double*)calloc(num_samples * 2, sizeof(double));
    double* window_sum = (double*)calloc(num_samples, sizeof(double));

    if (!fft_in || !fft_out || !ifft_out || !overlap_buf || !window_sum) {
        if (fft_in) fftw_free(fft_in);
        if (fft_out) fftw_free(fft_out);
        if (ifft_out) fftw_free(ifft_out);
        if (overlap_buf) free(overlap_buf);
        if (window_sum) free(window_sum);
        return MUSIC_ERROR_MEMORY;
    }

    fftw_plan plan_fwd = fftw_plan_dft_1d(N, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(N, fft_out, ifft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Noise power estimate (from first frame or estimator)
    double* noise_power = (double*)calloc(N, sizeof(double));
    if (!noise_power) {
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_inv);
        fftw_free(fft_in); fftw_free(fft_out); fftw_free(ifft_out);
        free(overlap_buf); free(window_sum);
        return MUSIC_ERROR_MEMORY;
    }

    // Estimate noise from first frame
    for (int i = 0; i < N; i++) {
        fft_in[i][0] = creal(iq_samples[i]) * subtractor->window[i];
        fft_in[i][1] = cimag(iq_samples[i]) * subtractor->window[i];
    }
    fftw_execute(plan_fwd);
    for (int i = 0; i < N; i++) {
        noise_power[i] = fft_out[i][0] * fft_out[i][0] + fft_out[i][1] * fft_out[i][1];
    }

    // STFT -> spectral subtraction -> ISTFT (overlap-add)
    for (int frame = 0; frame < num_frames; frame++) {
        int offset = frame * hop;

        // Window + FFT
        for (int i = 0; i < N; i++) {
            fft_in[i][0] = creal(iq_samples[offset + i]) * subtractor->window[i];
            fft_in[i][1] = cimag(iq_samples[offset + i]) * subtractor->window[i];
        }
        fftw_execute(plan_fwd);

        // Spectral subtraction: |S_clean|^2 = max(|S|^2 - alpha*N, beta*|S|^2)
        for (int i = 0; i < N; i++) {
            double mag_sq = fft_out[i][0] * fft_out[i][0] + fft_out[i][1] * fft_out[i][1];
            double clean_sq = mag_sq - alpha * noise_power[i];
            double floor_sq = beta * mag_sq;
            if (clean_sq < floor_sq) clean_sq = floor_sq;

            // Apply gain while preserving phase
            double gain = (mag_sq > 1e-20) ? sqrt(clean_sq / mag_sq) : 0.0;
            fft_out[i][0] *= gain;
            fft_out[i][1] *= gain;
        }

        // IFFT
        fftw_execute(plan_inv);

        // Overlap-add (IFFT normalization: 1/N)
        for (int i = 0; i < N; i++) {
            if (offset + i < num_samples) {
                overlap_buf[(offset + i) * 2 + 0] += ifft_out[i][0] / N * subtractor->window[i];
                overlap_buf[(offset + i) * 2 + 1] += ifft_out[i][1] / N * subtractor->window[i];
                window_sum[offset + i] += subtractor->window[i] * subtractor->window[i];
            }
        }
    }

    // Normalize by window overlap sum and write output
    for (int i = 0; i < num_samples; i++) {
        double w = (window_sum[i] > 1e-10) ? window_sum[i] : 1.0;
        output[i] = (overlap_buf[i * 2 + 0] / w) + I * (overlap_buf[i * 2 + 1] / w);
    }

    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(fft_in);
    fftw_free(fft_out);
    fftw_free(ifft_out);
    free(overlap_buf);
    free(window_sum);
    free(noise_power);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Wiener Filter
 * ============================================================================ */

music_status_t wiener_filter_init(
    wiener_filter_t* filter,
    const noise_reduction_config_t* config
) {
    if (!filter || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    filter->config = *config;

    music_status_t status = noise_floor_estimator_init(&filter->noise_estimator, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    filter->window = (double*)malloc(config->fft_size * sizeof(double));
    if (!filter->window) {
        noise_floor_estimator_free(&filter->noise_estimator);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < config->fft_size; i++) {
        filter->window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (config->fft_size - 1)));
    }

    return MUSIC_SUCCESS;
}

void wiener_filter_free(wiener_filter_t* filter) {
    if (filter) {
        noise_floor_estimator_free(&filter->noise_estimator);
        if (filter->window) {
            free(filter->window);
        }
    }
}

music_status_t wiener_compute_gain(
    const double* signal_power,
    const double* noise_power,
    int num_bins,
    double gain_min,
    double* wiener_gain
) {
    // compute_wiener_gain
    if (!signal_power || !noise_power || !wiener_gain) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    for (int i = 0; i < num_bins; i++) {
        double snr = signal_power[i] / (noise_power[i] + 1e-12);
        double gain = snr / (snr + 1.0);
        wiener_gain[i] = (gain > gain_min) ? gain : gain_min;
    }

    return MUSIC_SUCCESS;
}

music_status_t wiener_filter_process(
    wiener_filter_t* filter,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    if (!filter || !iq_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = filter->config.fft_size;
    int hop = (int)(N * (1.0 - filter->config.overlap));
    int num_frames = (num_samples - N) / hop + 1;
    double gain_min = filter->config.wiener_gain_min;

    if (num_frames < 1 || N <= 0) {
        memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));
        return MUSIC_SUCCESS;
    }

    fftw_complex* fft_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* ifft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    double* overlap_buf = (double*)calloc(num_samples * 2, sizeof(double));
    double* window_sum = (double*)calloc(num_samples, sizeof(double));
    double* noise_power = (double*)calloc(N, sizeof(double));

    if (!fft_in || !fft_out || !ifft_out || !overlap_buf || !window_sum || !noise_power) {
        if (fft_in) fftw_free(fft_in);
        if (fft_out) fftw_free(fft_out);
        if (ifft_out) fftw_free(ifft_out);
        if (overlap_buf) free(overlap_buf);
        if (window_sum) free(window_sum);
        if (noise_power) free(noise_power);
        return MUSIC_ERROR_MEMORY;
    }

    fftw_plan plan_fwd = fftw_plan_dft_1d(N, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(N, fft_out, ifft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Estimate noise from first frame
    for (int i = 0; i < N; i++) {
        fft_in[i][0] = creal(iq_samples[i]) * filter->window[i];
        fft_in[i][1] = cimag(iq_samples[i]) * filter->window[i];
    }
    fftw_execute(plan_fwd);
    for (int i = 0; i < N; i++) {
        noise_power[i] = fft_out[i][0] * fft_out[i][0] + fft_out[i][1] * fft_out[i][1];
    }

    // STFT -> Wiener filter -> ISTFT (overlap-add)
    for (int frame = 0; frame < num_frames; frame++) {
        int offset = frame * hop;

        for (int i = 0; i < N; i++) {
            fft_in[i][0] = creal(iq_samples[offset + i]) * filter->window[i];
            fft_in[i][1] = cimag(iq_samples[offset + i]) * filter->window[i];
        }
        fftw_execute(plan_fwd);

        // Wiener gain: G = SNR / (SNR + 1) = (|S|^2 - N) / |S|^2, clamped to gain_min
        for (int i = 0; i < N; i++) {
            double mag_sq = fft_out[i][0] * fft_out[i][0] + fft_out[i][1] * fft_out[i][1];
            double snr = mag_sq / (noise_power[i] + 1e-12);
            double gain = snr / (snr + 1.0);
            if (gain < gain_min) gain = gain_min;
            fft_out[i][0] *= gain;
            fft_out[i][1] *= gain;
        }

        fftw_execute(plan_inv);

        for (int i = 0; i < N; i++) {
            if (offset + i < num_samples) {
                overlap_buf[(offset + i) * 2 + 0] += ifft_out[i][0] / N * filter->window[i];
                overlap_buf[(offset + i) * 2 + 1] += ifft_out[i][1] / N * filter->window[i];
                window_sum[offset + i] += filter->window[i] * filter->window[i];
            }
        }
    }

    for (int i = 0; i < num_samples; i++) {
        double w = (window_sum[i] > 1e-10) ? window_sum[i] : 1.0;
        output[i] = (overlap_buf[i * 2 + 0] / w) + I * (overlap_buf[i * 2 + 1] / w);
    }

    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(fft_in);
    fftw_free(fft_out);
    fftw_free(ifft_out);
    free(overlap_buf);
    free(window_sum);
    free(noise_power);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Spatial Noise Canceller
 * ============================================================================ */

music_status_t spatial_canceller_init(
    spatial_noise_canceller_t* canceller,
    const noise_reduction_config_t* config,
    const double* antenna_positions
) {
    if (!canceller || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    canceller->config = *config;
    canceller->num_elements = config->num_elements;

    canceller->antenna_positions = (double*)malloc(config->num_elements * 3 * sizeof(double));
    if (!canceller->antenna_positions) {
        return MUSIC_ERROR_MEMORY;
    }

    if (antenna_positions) {
        memcpy(canceller->antenna_positions, antenna_positions,
               config->num_elements * 3 * sizeof(double));
    } else {
        // Default UCA positions
        double radius = NR_UCA_RADIUS;
        for (int i = 0; i < config->num_elements; i++) {
            double phi = 2.0 * M_PI * i / config->num_elements;
            canceller->antenna_positions[i * 3 + 0] = radius * cos(phi);  // x
            canceller->antenna_positions[i * 3 + 1] = radius * sin(phi);  // y
            canceller->antenna_positions[i * 3 + 2] = 0.0;                // z
        }
    }

    return MUSIC_SUCCESS;
}

void spatial_canceller_free(spatial_noise_canceller_t* canceller) {
    if (canceller && canceller->antenna_positions) {
        free(canceller->antenna_positions);
        canceller->antenna_positions = NULL;
    }
}

/**
 * @brief Internal: compute UCA steering vector for given azimuth
 *
 * a_m(theta) = exp(j * k * r * cos(theta - phi_m))
 * where phi_m = 2*pi*m/M is the angular position of element m
 *
 * The vector is normalized to unit norm: ||a|| = 1
 */
static void nr_compute_uca_steering_vector(
    int num_elements,
    double radius,
    double frequency,
    double azimuth_deg,
    cdouble_t* steering_vector
) {
    double wavelength = SPEED_OF_LIGHT / frequency;
    double k = 2.0 * M_PI / wavelength;
    double azimuth_rad = azimuth_deg * M_PI / 180.0;

    for (int m = 0; m < num_elements; m++) {
        double phi_m = 2.0 * M_PI * m / num_elements;
        double x_m = radius * cos(phi_m);
        double y_m = radius * sin(phi_m);
        double phase = k * (x_m * cos(azimuth_rad) + y_m * sin(azimuth_rad));
        steering_vector[m] = cexp(I * phase);
    }

    /* Normalize to unit norm */
    double norm = 0.0;
    for (int m = 0; m < num_elements; m++) {
        norm += creal(steering_vector[m] * conj(steering_vector[m]));
    }
    norm = sqrt(norm);
    if (norm > 1e-15) {
        for (int m = 0; m < num_elements; m++) {
            steering_vector[m] /= norm;
        }
    }
}

/**
 * @brief Internal: invert a small complex matrix using LAPACK (zgetrf + zgetri)
 *
 * @param mat  Input matrix (column-major), overwritten on output
 * @param n    Matrix dimension
 * @param inv  Output inverse (column-major), caller-allocated [n*n]
 * @return MUSIC_SUCCESS or error code
 */
static music_status_t nr_invert_complex_matrix(
    const cdouble_t* mat,
    int n,
    cdouble_t* inv
) {
    memcpy(inv, mat, n * n * sizeof(cdouble_t));

    lapack_int* ipiv = (lapack_int*)malloc(n * sizeof(lapack_int));
    if (!ipiv) {
        return MUSIC_ERROR_MEMORY;
    }

    lapack_int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,
                                     (lapack_complex_double*)inv, n, ipiv);
    if (info != 0) {
        free(ipiv);
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    info = LAPACKE_zgetri(LAPACK_COL_MAJOR, n,
                          (lapack_complex_double*)inv, n, ipiv);
    free(ipiv);

    if (info != 0) {
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    return MUSIC_SUCCESS;
}

music_status_t spatial_compute_null_weights(
    const spatial_noise_canceller_t* canceller,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* weights
) {
    /*
     * LCMV beamforming: w = R^{-1} C (C^H R^{-1} C)^{-1} f
     *
     * When no covariance data is available, use diagonal-loaded identity
     * as R, which makes R^{-1} = (1/loading) * I.
     *
     * Constraint matrix C:
     *   Column 0: a(desired_azimuth)   -- unit gain
     *   Column 1..K: a(interference_i) -- null
     *
     * Response vector f: [1, 0, 0, ..., 0]^T
     *
     * If no interferers and desired_azimuth >= 0, falls back to MVDR
     * (single constraint LCMV).
     */

    if (!canceller || !weights) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (frequency <= 0.0) {
        /* Cannot compute steering vectors without valid frequency */
        for (int i = 0; i < canceller->num_elements; i++) {
            weights[i] = 1.0 / canceller->num_elements;
        }
        return MUSIC_SUCCESS;
    }

    int M = canceller->num_elements;
    double radius = NR_UCA_RADIUS;

    /* Determine number of constraints */
    int num_constraints;
    if (desired_azimuth >= 0.0) {
        num_constraints = 1 + num_interferers;
    } else {
        /* No desired direction -- only null interferers if any */
        if (num_interferers <= 0 || !interference_azimuths) {
            /* No constraints at all: uniform weights */
            for (int i = 0; i < M; i++) {
                weights[i] = 1.0 / M;
            }
            return MUSIC_SUCCESS;
        }
        /* Null interferers only -- use first interferer column as "desired"
         * with response 0, but this is unusual. Better: set desired to 0 deg
         * with unit gain and null all interferers. */
        num_constraints = num_interferers;
    }

    if (num_constraints > M) {
        /* Over-determined: too many constraints for array size */
        num_constraints = M - 1;
    }
    if (num_constraints < 1) {
        for (int i = 0; i < M; i++) {
            weights[i] = 1.0 / M;
        }
        return MUSIC_SUCCESS;
    }

    /* -----------------------------------------------------------------
     * Allocate working memory
     * C: M x num_constraints (column-major)
     * f: num_constraints x 1
     * R_inv: M x M
     * R_inv_C: M x num_constraints
     * G: num_constraints x num_constraints  (= C^H R_inv C)
     * G_inv: num_constraints x num_constraints
     * temp_vec: num_constraints x 1  (= G_inv * f)
     * ----------------------------------------------------------------- */
    cdouble_t* C = (cdouble_t*)calloc(M * num_constraints, sizeof(cdouble_t));
    cdouble_t* f = (cdouble_t*)calloc(num_constraints, sizeof(cdouble_t));
    cdouble_t* R_inv = (cdouble_t*)calloc(M * M, sizeof(cdouble_t));
    cdouble_t* R_inv_C = (cdouble_t*)calloc(M * num_constraints, sizeof(cdouble_t));
    cdouble_t* G = (cdouble_t*)calloc(num_constraints * num_constraints, sizeof(cdouble_t));
    cdouble_t* G_inv = (cdouble_t*)calloc(num_constraints * num_constraints, sizeof(cdouble_t));
    cdouble_t* temp_vec = (cdouble_t*)calloc(num_constraints, sizeof(cdouble_t));

    if (!C || !f || !R_inv || !R_inv_C || !G || !G_inv || !temp_vec) {
        free(C); free(f); free(R_inv); free(R_inv_C);
        free(G); free(G_inv); free(temp_vec);
        return MUSIC_ERROR_MEMORY;
    }

    /* -----------------------------------------------------------------
     * Build R_inv = (1/loading) * I  (diagonal-loaded identity inverse)
     * ----------------------------------------------------------------- */
    double loading = NR_LCMV_DIAGONAL_LOADING;
    double inv_loading = 1.0 / loading;
    for (int i = 0; i < M; i++) {
        R_inv[i + i * M] = inv_loading + 0.0 * I;
    }

    /* -----------------------------------------------------------------
     * Build constraint matrix C and response vector f
     * ----------------------------------------------------------------- */
    int col = 0;

    if (desired_azimuth >= 0.0) {
        /* Column 0: desired direction with unit gain */
        nr_compute_uca_steering_vector(M, radius, frequency, desired_azimuth,
                                       &C[col * M]);
        f[col] = 1.0 + 0.0 * I;
        col++;
    }

    /* Interference columns: null response */
    if (interference_azimuths) {
        int max_int = num_constraints - col;
        for (int i = 0; i < max_int && i < num_interferers; i++) {
            nr_compute_uca_steering_vector(M, radius, frequency,
                                           interference_azimuths[i],
                                           &C[col * M]);
            f[col] = 0.0 + 0.0 * I;
            col++;
        }
    }

    /* Actual number of constraints used */
    int Nc = col;

    /* -----------------------------------------------------------------
     * Compute R_inv_C = R_inv * C   (M x Nc)
     * Using cblas_zgemm: R_inv_C = alpha * R_inv * C + beta * R_inv_C
     * ----------------------------------------------------------------- */
    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, Nc, M,
                    &alpha_blas, R_inv, M,
                    C, M,
                    &beta_blas, R_inv_C, M);
    }

    /* -----------------------------------------------------------------
     * Compute G = C^H * R_inv_C   (Nc x Nc)
     * ----------------------------------------------------------------- */
    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
                    Nc, Nc, M,
                    &alpha_blas, C, M,
                    R_inv_C, M,
                    &beta_blas, G, Nc);
    }

    /* -----------------------------------------------------------------
     * Invert G -> G_inv  (small Nc x Nc matrix)
     * ----------------------------------------------------------------- */
    music_status_t status = nr_invert_complex_matrix(G, Nc, G_inv);
    if (status != MUSIC_SUCCESS) {
        /* Fallback: uniform weights on failure */
        for (int i = 0; i < M; i++) {
            weights[i] = 1.0 / M;
        }
        free(C); free(f); free(R_inv); free(R_inv_C);
        free(G); free(G_inv); free(temp_vec);
        return MUSIC_SUCCESS;
    }

    /* -----------------------------------------------------------------
     * Compute temp_vec = G_inv * f   (Nc x 1)
     * ----------------------------------------------------------------- */
    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    Nc, Nc,
                    &alpha_blas, G_inv, Nc,
                    f, 1,
                    &beta_blas, temp_vec, 1);
    }

    /* -----------------------------------------------------------------
     * Compute weights = R_inv_C * temp_vec   (M x 1)
     * w = R_inv * C * (C^H * R_inv * C)^{-1} * f
     * ----------------------------------------------------------------- */
    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    M, Nc,
                    &alpha_blas, R_inv_C, M,
                    temp_vec, 1,
                    &beta_blas, weights, 1);
    }

    /* Cleanup */
    free(C);
    free(f);
    free(R_inv);
    free(R_inv_C);
    free(G);
    free(G_inv);
    free(temp_vec);

    return MUSIC_SUCCESS;
}

music_status_t spatial_estimate_interference(
    const spatial_noise_canceller_t* canceller,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    double frequency,
    int num_sources,
    double* interference_azimuths
) {
    /*
     * Spatial spectrum scanning to estimate interference directions.
     *
     * Algorithm:
     *   1. Estimate spatial covariance matrix R from multichannel data
     *   2. Scan azimuth 0..359 degrees (1-degree resolution)
     *   3. Compute conventional beamformer output power: P(theta) = a^H R a
     *   4. Find the num_sources strongest peaks (with local maximum criterion)
     */

    if (!canceller || !multichannel_samples || !interference_azimuths) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (num_elements < 1 || num_samples < 1 || num_sources < 1) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    if (frequency <= 0.0) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    int M = num_elements;
    double radius = NR_UCA_RADIUS;

    /* -----------------------------------------------------------------
     * Step 1: Estimate covariance matrix R = (1/N) * X * X^H
     * Data layout: multichannel_samples[m + n*M] = element m, sample n
     * For cblas_zgemm we need column-major M x N layout, which is what
     * we already have.
     * ----------------------------------------------------------------- */
    cdouble_t* R = (cdouble_t*)calloc(M * M, sizeof(cdouble_t));
    if (!R) {
        return MUSIC_ERROR_MEMORY;
    }

    {
        cdouble_t alpha_blas = (1.0 / num_samples) + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                    M, M, num_samples,
                    &alpha_blas, multichannel_samples, M,
                    multichannel_samples, M,
                    &beta_blas, R, M);
    }

    /* Add diagonal loading for numerical stability */
    double trace = 0.0;
    for (int i = 0; i < M; i++) {
        trace += creal(R[i + i * M]);
    }
    double loading = NR_LCMV_DIAGONAL_LOADING * trace / M;
    for (int i = 0; i < M; i++) {
        R[i + i * M] += loading;
    }

    /* -----------------------------------------------------------------
     * Step 2: Scan spatial spectrum P(theta) = a(theta)^H * R * a(theta)
     * ----------------------------------------------------------------- */
    int num_scan = 360;
    double* power_spectrum = (double*)malloc(num_scan * sizeof(double));
    cdouble_t* a = (cdouble_t*)malloc(M * sizeof(cdouble_t));
    cdouble_t* Ra = (cdouble_t*)malloc(M * sizeof(cdouble_t));

    if (!power_spectrum || !a || !Ra) {
        free(R);
        if (power_spectrum) free(power_spectrum);
        if (a) free(a);
        if (Ra) free(Ra);
        return MUSIC_ERROR_MEMORY;
    }

    for (int deg = 0; deg < num_scan; deg++) {
        double azimuth = (double)deg;
        nr_compute_uca_steering_vector(M, radius, frequency, azimuth, a);

        /* Ra = R * a */
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    M, M,
                    &alpha_blas, R, M,
                    a, 1,
                    &beta_blas, Ra, 1);

        /* P = a^H * Ra (real part of inner product) */
        cdouble_t p = 0.0;
        for (int m = 0; m < M; m++) {
            p += conj(a[m]) * Ra[m];
        }
        power_spectrum[deg] = creal(p);
    }

    /* -----------------------------------------------------------------
     * Step 3: Find peaks (local maxima) and pick the strongest N
     *
     * A local maximum at index i means:
     *   power_spectrum[i] > power_spectrum[i-1] AND
     *   power_spectrum[i] > power_spectrum[i+1]
     * (with wraparound for circular scan)
     * ----------------------------------------------------------------- */
    /* Find all local maxima */
    int max_peaks = num_scan;
    int* peak_indices = (int*)malloc(max_peaks * sizeof(int));
    double* peak_values = (double*)malloc(max_peaks * sizeof(double));
    int num_peaks = 0;

    if (!peak_indices || !peak_values) {
        free(R); free(power_spectrum); free(a); free(Ra);
        if (peak_indices) free(peak_indices);
        if (peak_values) free(peak_values);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < num_scan; i++) {
        int prev = (i - 1 + num_scan) % num_scan;
        int next = (i + 1) % num_scan;
        if (power_spectrum[i] > power_spectrum[prev] &&
            power_spectrum[i] > power_spectrum[next]) {
            peak_indices[num_peaks] = i;
            peak_values[num_peaks] = power_spectrum[i];
            num_peaks++;
        }
    }

    /* Sort peaks by descending power (simple selection sort, small N) */
    for (int i = 0; i < num_peaks - 1; i++) {
        int max_idx = i;
        for (int j = i + 1; j < num_peaks; j++) {
            if (peak_values[j] > peak_values[max_idx]) {
                max_idx = j;
            }
        }
        if (max_idx != i) {
            double tmp_v = peak_values[i];
            peak_values[i] = peak_values[max_idx];
            peak_values[max_idx] = tmp_v;
            int tmp_i = peak_indices[i];
            peak_indices[i] = peak_indices[max_idx];
            peak_indices[max_idx] = tmp_i;
        }
    }

    /* Fill output: top num_sources peaks */
    int found = (num_peaks < num_sources) ? num_peaks : num_sources;
    for (int i = 0; i < found; i++) {
        interference_azimuths[i] = (double)peak_indices[i];
    }
    /* If fewer peaks than requested, fill remaining with -1 */
    for (int i = found; i < num_sources; i++) {
        interference_azimuths[i] = -1.0;
    }

    /* Cleanup */
    free(R);
    free(power_spectrum);
    free(a);
    free(Ra);
    free(peak_indices);
    free(peak_values);

    return MUSIC_SUCCESS;
}

music_status_t spatial_canceller_process(
    spatial_noise_canceller_t* canceller,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* output
) {
    // process_multichannel
    if (!canceller || !multichannel_samples || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute weights
    cdouble_t* weights = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));
    if (!weights) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = spatial_compute_null_weights(canceller, interference_azimuths,
                                                          num_interferers, frequency,
                                                          desired_azimuth, weights);
    if (status != MUSIC_SUCCESS) {
        free(weights);
        return status;
    }

    // Apply beamforming: y = w^H * x
    for (int n = 0; n < num_samples; n++) {
        output[n] = 0.0;
        for (int m = 0; m < num_elements; m++) {
            output[n] += conj(weights[m]) * multichannel_samples[m + n * num_elements];
        }
    }

    free(weights);
    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Unified Noise Reducer
 * ============================================================================ */

music_status_t noise_reducer_init(
    noise_reducer_t* reducer,
    const noise_reduction_config_t* config
) {
    if (!reducer) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    noise_reduction_config_t default_config;
    if (!config) {
        noise_reduction_config_default(&default_config);
        config = &default_config;
    }

    reducer->config = *config;

    music_status_t status = spectral_subtractor_init(&reducer->spectral_subtractor, config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    status = wiener_filter_init(&reducer->wiener_filter, config);
    if (status != MUSIC_SUCCESS) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        return status;
    }

    status = spatial_canceller_init(&reducer->spatial_canceller, config, NULL);
    if (status != MUSIC_SUCCESS) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        wiener_filter_free(&reducer->wiener_filter);
        return status;
    }

    return MUSIC_SUCCESS;
}

void noise_reducer_free(noise_reducer_t* reducer) {
    if (reducer) {
        spectral_subtractor_free(&reducer->spectral_subtractor);
        wiener_filter_free(&reducer->wiener_filter);
        spatial_canceller_free(&reducer->spatial_canceller);
    }
}

music_status_t noise_reducer_process_single(
    noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    noise_reduction_method_t method,
    cdouble_t* output,
    noise_reduction_metrics_t* metrics
) {
    // process_single_channel
    if (!reducer || !iq_samples || !output || !metrics) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Measure input SNR
    double input_snr;
    noise_reducer_estimate_snr(reducer, iq_samples, num_samples, &input_snr);

    music_status_t status = MUSIC_SUCCESS;

    if (method == NR_SPECTRAL_SUBTRACTION) {
        status = spectral_subtractor_process(&reducer->spectral_subtractor,
                                             iq_samples, num_samples, output);
        strcpy(metrics->method, "SpectralSubtraction");
    } else if (method == NR_WIENER_FILTER) {
        status = wiener_filter_process(&reducer->wiener_filter,
                                       iq_samples, num_samples, output);
        strcpy(metrics->method, "WienerFilter");
    } else {
        memcpy(output, iq_samples, num_samples * sizeof(cdouble_t));
        strcpy(metrics->method, "None");
    }

    // Measure output SNR
    double output_snr;
    noise_reducer_estimate_snr(reducer, output, num_samples, &output_snr);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_ms = (end.tv_sec - start.tv_sec) * 1000.0 +
                       (end.tv_nsec - start.tv_nsec) / 1000000.0;

    // Fill metrics
    metrics->input_snr_db = input_snr;
    metrics->output_snr_db = output_snr;
    metrics->snr_improvement_db = output_snr - input_snr;
    metrics->processing_time_ms = elapsed_ms;
    metrics->samples_processed = num_samples;

    return status;
}

music_status_t noise_reducer_estimate_snr(
    const noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    double* snr_db
) {
    // estimate_snr
    if (!reducer || !iq_samples || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Simple FFT-based SNR estimation
    int fft_size = (num_samples < reducer->config.fft_size) ?
                   num_samples : reducer->config.fft_size;

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);

    if (!in || !out) {
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return MUSIC_ERROR_MEMORY;
    }

    // Copy samples
    for (int i = 0; i < fft_size; i++) {
        in[i][0] = creal(iq_samples[i]);
        in[i][1] = cimag(iq_samples[i]);
    }

    // FFT
    fftw_plan plan = fftw_plan_dft_1d(fft_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Compute power spectrum
    double max_power = 0.0;
    double* power = (double*)malloc(fft_size * sizeof(double));
    double* sorted = (double*)malloc(fft_size * sizeof(double));
    if (!power || !sorted) {
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
        if (power) free(power);
        if (sorted) free(sorted);
        return MUSIC_ERROR_MEMORY;
    }
    for (int i = 0; i < fft_size; i++) {
        power[i] = out[i][0] * out[i][0] + out[i][1] * out[i][1];
        if (power[i] > max_power) {
            max_power = power[i];
        }
    }

    // Estimate noise as median
    memcpy(sorted, power, fft_size * sizeof(double));

    for (int i = 1; i < fft_size; i++) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) {
            sorted[j + 1] = sorted[j];
            j--;
        }
        sorted[j + 1] = key;
    }

    double noise_power = sorted[fft_size / 2];

    *snr_db = 10.0 * log10((max_power / (noise_power + 1e-12)) + 1e-12);

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    free(power);
    free(sorted);

    return MUSIC_SUCCESS;
}

music_status_t create_noise_reducer(
    noise_reducer_t* reducer,
    double sample_rate,
    int fft_size,
    double alpha,
    double adaptation_rate
) {
    // create_noise_reducer
    noise_reduction_config_t config;
    noise_reduction_config_default(&config);

    config.sample_rate = sample_rate;
    config.fft_size = fft_size;
    config.alpha = alpha;
    config.adaptation_rate = adaptation_rate;

    return noise_reducer_init(reducer, &config);
}

music_status_t noise_reducer_process_multi(
    noise_reducer_t* reducer,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    double frequency,
    const double* interference_azimuths,
    int num_interferers,
    double desired_azimuth,
    noise_reduction_method_t temporal_method,
    cdouble_t* output,
    noise_reduction_metrics_t* metrics
) {
    /* Stub: multi-channel noise reduction not yet implemented */
    (void)reducer; (void)multichannel_samples; (void)num_elements;
    (void)num_samples; (void)frequency; (void)interference_azimuths;
    (void)num_interferers; (void)desired_azimuth; (void)temporal_method;
    (void)output; (void)metrics;
    return MUSIC_ERROR_INVALID_CONFIG;
}
