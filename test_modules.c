/**
 * @file test_modules.c
 * @brief Comprehensive unit tests for MQP system modules
 *
 * Tests: noise_reduction, beamforming, covariance, eigendecomp,
 *        cfar_2d, state_machine, lte_detector, p25_detector
 *
 * Build (from c_implementation/):
 *   gcc -Iinclude -o test/test_modules test/test_modules.c \
 *       src/noise_reduction.c src/beamforming.c src/music_uca_6.c \
 *       src/cfar_2d.c src/state_machine.c src/lte_detector.c \
 *       src/p25_detector.c src/common_types.c \
 *       -lfftw3 -lopenblas -lm -lpthread
 */

#include "noise_reduction.h"
#include "beamforming.h"
#include "music_uca_6.h"
#include "music_config.h"
#include "cfar_2d.h"
#include "state_machine.h"
#include "lte_detector.h"
#include "p25_detector.h"
#include "common_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * Test Framework
 * ============================================================================ */

static int test_passed = 0;
static int test_failed = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        if (condition) { \
            test_passed++; \
            printf("  PASS: %s\n", message); \
        } else { \
            test_failed++; \
            printf("  FAIL: %s\n", message); \
        } \
    } while(0)

/* ============================================================================
 * Noise Reduction Tests
 * ============================================================================ */

/**
 * Test 1: Noise floor estimator init, feed known spectrum, verify estimate.
 * Then verify EMA update changes the noise floor.
 */
static void test_noise_floor_estimator(void)
{
    printf("\n[TEST] test_noise_floor_estimator\n");

    noise_reduction_config_t config;
    noise_reduction_config_default(&config);
    config.fft_size = 256;

    noise_floor_estimator_t estimator;
    music_status_t status = noise_floor_estimator_init(&estimator, &config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "noise_floor_estimator_init returns MUSIC_SUCCESS");
    TEST_ASSERT(estimator.initialized == true, "estimator is marked as initialized");

    int num_bins = config.fft_size / 2;
    double* spectrum = (double*)calloc(num_bins, sizeof(double));
    double* noise_floor = (double*)calloc(num_bins, sizeof(double));

    /* Fill spectrum with known flat noise at 0.5 */
    for (int i = 0; i < num_bins; i++) {
        spectrum[i] = 0.5;
    }

    status = noise_floor_estimate(&estimator, spectrum, num_bins, noise_floor);
    TEST_ASSERT(status == MUSIC_SUCCESS, "noise_floor_estimate returns MUSIC_SUCCESS");

    /* Noise floor should be close to the constant input (all bins = 0.5) */
    int close_count = 0;
    for (int i = 0; i < num_bins; i++) {
        if (fabs(noise_floor[i] - 0.5) < 0.3) {
            close_count++;
        }
    }
    TEST_ASSERT(close_count > num_bins / 2,
        "noise floor estimate is close to input for majority of bins");

    /* Save first bin value before update */
    double before_update = estimator.noise_floor ? estimator.noise_floor[0] : -1.0;

    /* Update with a different spectrum (higher noise) */
    for (int i = 0; i < num_bins; i++) {
        spectrum[i] = 5.0;
    }
    status = noise_floor_update(&estimator, spectrum, num_bins);
    TEST_ASSERT(status == MUSIC_SUCCESS, "noise_floor_update returns MUSIC_SUCCESS");

    /* After EMA update with higher values, noise floor should have shifted up */
    double after_update = estimator.noise_floor ? estimator.noise_floor[0] : -1.0;
    TEST_ASSERT(after_update > before_update || before_update < 0.0,
        "EMA update shifts noise floor toward new spectrum");

    noise_floor_estimator_free(&estimator);
    free(spectrum);
    free(noise_floor);
}

/**
 * Test 2: Wiener gain computation.
 * signal_power=[10,20,30], noise_power=[5,5,5]
 * Expected gain: SNR/(SNR+1) = [0.5, 0.75, 0.833], clamped by gain_min.
 */
static void test_wiener_gain(void)
{
    printf("\n[TEST] test_wiener_gain\n");

    double signal_power[3] = {10.0, 20.0, 30.0};
    double noise_power[3]  = {5.0,  5.0,  5.0};
    double gain[3] = {0.0};
    double gain_min = 0.1;

    music_status_t status = wiener_compute_gain(
        signal_power, noise_power, 3, gain_min, gain);
    TEST_ASSERT(status == MUSIC_SUCCESS, "wiener_compute_gain returns MUSIC_SUCCESS");

    /* SNR = (S - N) / N for Wiener, but gain = SNR/(SNR+1) = (S-N)/S
     * Or equivalently: gain = max(1 - N/S, gain_min)
     * Bin 0: SNR = (10-5)/5 = 1.0, gain = 1/(1+1) = 0.5
     * Bin 1: SNR = (20-5)/5 = 3.0, gain = 3/(3+1) = 0.75
     * Bin 2: SNR = (30-5)/5 = 5.0, gain = 5/(5+1) = 0.833
     */
    TEST_ASSERT(fabs(gain[0] - 0.5) < 0.15,
        "gain[0] ~ 0.5 for SNR=1");
    TEST_ASSERT(fabs(gain[1] - 0.75) < 0.15,
        "gain[1] ~ 0.75 for SNR=3");
    TEST_ASSERT(fabs(gain[2] - 0.833) < 0.15,
        "gain[2] ~ 0.833 for SNR=5");

    /* Verify gain_min clamp: if noise > signal, gain should be clamped */
    double low_signal[1] = {1.0};
    double high_noise[1] = {100.0};
    double clamped_gain[1] = {0.0};
    status = wiener_compute_gain(low_signal, high_noise, 1, gain_min, clamped_gain);
    TEST_ASSERT(status == MUSIC_SUCCESS, "wiener_compute_gain with low SNR succeeds");
    TEST_ASSERT(clamped_gain[0] >= gain_min - 1e-9,
        "gain is clamped to gain_min when noise >> signal");
}

/**
 * Test 3: noise_reducer init and free (no crash).
 */
static void test_noise_reducer_init_free(void)
{
    printf("\n[TEST] test_noise_reducer_init_free\n");

    noise_reducer_t reducer;
    memset(&reducer, 0, sizeof(reducer));

    noise_reduction_config_t config;
    noise_reduction_config_default(&config);

    music_status_t status = noise_reducer_init(&reducer, &config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "noise_reducer_init returns MUSIC_SUCCESS");

    /* Free should not crash */
    noise_reducer_free(&reducer);
    TEST_ASSERT(1, "noise_reducer_free completes without crash");
}

/* ============================================================================
 * Beamforming Tests
 * ============================================================================ */

/**
 * Test 4: Init conventional beamformer, compute weights for azimuth=45 deg.
 */
static void test_conventional_beamformer(void)
{
    printf("\n[TEST] test_conventional_beamformer\n");

    beamforming_config_t config;
    beamforming_config_default(&config);

    conventional_beamformer_t bf;
    memset(&bf, 0, sizeof(bf));
    music_status_t status = conventional_beamformer_init(&bf, &config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "conventional_beamformer_init returns MUSIC_SUCCESS");

    cdouble_t weights[BF_NUM_ELEMENTS];
    memset(weights, 0, sizeof(weights));

    status = conventional_compute_weights(&bf, 45.0, 0.0, true, weights);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "conventional_compute_weights(az=45) returns MUSIC_SUCCESS");

    /* Verify 6 weights are non-zero */
    int nonzero_count = 0;
    for (int i = 0; i < BF_NUM_ELEMENTS; i++) {
        if (cabs(weights[i]) > 1e-15) {
            nonzero_count++;
        }
    }
    TEST_ASSERT(nonzero_count == BF_NUM_ELEMENTS,
        "all 6 beamforming weights are non-zero");

    conventional_beamformer_free(&bf);
}

/**
 * Test 5: Compute beam pattern 0-360 deg, verify peak at look direction.
 */
static void test_beam_pattern(void)
{
    printf("\n[TEST] test_beam_pattern\n");

    beamforming_config_t config;
    beamforming_config_default(&config);

    conventional_beamformer_t bf;
    memset(&bf, 0, sizeof(bf));
    music_status_t status = conventional_beamformer_init(&bf, &config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "beamformer init for pattern test");

    double look_az = 90.0;
    beam_pattern_t* pattern = beam_pattern_alloc(361);
    TEST_ASSERT(pattern != NULL, "beam_pattern_alloc returns non-NULL");

    if (pattern != NULL) {
        status = conventional_compute_pattern(
            &bf, look_az, -180.0, 180.0, 1.0, pattern);
        TEST_ASSERT(status == MUSIC_SUCCESS,
            "conventional_compute_pattern returns MUSIC_SUCCESS");

        if (status == MUSIC_SUCCESS && pattern->num_points > 0) {
            /* Find the peak in the pattern */
            int peak_idx = 0;
            double peak_val = -1e30;
            for (int i = 0; i < pattern->num_points; i++) {
                if (pattern->pattern_db[i] > peak_val) {
                    peak_val = pattern->pattern_db[i];
                    peak_idx = i;
                }
            }
            double peak_az = pattern->azimuths[peak_idx];

            /* Peak should be within 5 degrees of look direction */
            double diff = fabs(peak_az - look_az);
            if (diff > 180.0) diff = 360.0 - diff;
            TEST_ASSERT(diff < 5.0,
                "beam pattern peak is within 5 deg of look direction (90 deg)");

            /* Verify pattern has variation (indicating nulls/sidelobes) */
            double min_val = 1e30;
            for (int i = 0; i < pattern->num_points; i++) {
                if (pattern->pattern_db[i] < min_val) {
                    min_val = pattern->pattern_db[i];
                }
            }
            TEST_ASSERT(peak_val - min_val > 3.0,
                "beam pattern has > 3 dB dynamic range (nulls present)");
        }

        beam_pattern_free(pattern);
    }

    conventional_beamformer_free(&bf);
}

/**
 * Test 6: Verify default beamforming config values.
 */
static void test_beamforming_config_default(void)
{
    printf("\n[TEST] test_beamforming_config_default\n");

    beamforming_config_t config;
    memset(&config, 0, sizeof(config));
    beamforming_config_default(&config);

    TEST_ASSERT(config.num_elements == 6,
        "default config num_elements == 6");
    TEST_ASSERT(config.radius > 0.0,
        "default config radius > 0");
    TEST_ASSERT(config.center_frequency > 0.0,
        "default config center_frequency > 0");
    TEST_ASSERT(config.diagonal_loading > 0.0,
        "default config diagonal_loading > 0");
}

/* ============================================================================
 * Covariance Tests
 * ============================================================================ */

/**
 * Test 7: Compute covariance from synthetic 6-channel signal, verify Hermitian.
 */
static void test_covariance_computation(void)
{
    printf("\n[TEST] test_covariance_computation\n");

    int M = 6;
    int N = 256;

    /* Allocate signal matrix: 6 antennas x 256 snapshots */
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    TEST_ASSERT(signals != NULL, "signal_matrix_alloc returns non-NULL");
    if (signals == NULL) return;

    /* Generate synthetic signal: single source at 45 degrees
     * Phase shift for UCA element m: phi_m = k*R*cos(theta - 2*pi*m/M)
     * where theta = 45 deg = pi/4
     */
    double radius = ARRAY_RADIUS_M;
    double k = 2.0 * M_PI / WAVELENGTH_M;
    double theta = 45.0 * M_PI / 180.0;

    for (int n = 0; n < N; n++) {
        /* Base signal: unit complex sinusoid */
        double phase_signal = 2.0 * M_PI * 0.1 * n;
        cdouble_t base_signal = cexp(I * phase_signal);

        for (int m = 0; m < M; m++) {
            double phi_m = 2.0 * M_PI * m / (double)M;
            double phase_shift = k * radius * cos(theta - phi_m);
            cdouble_t steering = cexp(I * phase_shift);
            /* Column-major: element (m, n) at index m + n*M */
            signals->data[m + n * M] = base_signal * steering;
        }
    }

    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    TEST_ASSERT(cov != NULL, "covariance_matrix_alloc returns non-NULL");
    if (cov == NULL) {
        signal_matrix_free(signals);
        return;
    }

    music_status_t status = compute_covariance_matrix(signals, false, cov);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "compute_covariance_matrix returns MUSIC_SUCCESS");

    /* Verify Hermitian property: R[i,j] = conj(R[j,i])
     * Column-major: R[i,j] = data[i + j*M]
     */
    if (status == MUSIC_SUCCESS) {
        int hermitian_ok = 1;
        for (int i = 0; i < M; i++) {
            for (int j = i + 1; j < M; j++) {
                cdouble_t r_ij = cov->data[i + j * M];
                cdouble_t r_ji = cov->data[j + i * M];
                cdouble_t diff = r_ij - conj(r_ji);
                if (cabs(diff) > 1e-6) {
                    hermitian_ok = 0;
                    break;
                }
            }
            if (!hermitian_ok) break;
        }
        TEST_ASSERT(hermitian_ok,
            "covariance matrix is Hermitian: R[i,j] == conj(R[j,i])");

        /* Verify diagonal elements are real and positive */
        int diag_ok = 1;
        for (int i = 0; i < M; i++) {
            cdouble_t diag = cov->data[i + i * M];
            if (fabs(cimag(diag)) > 1e-6 || creal(diag) < 0.0) {
                diag_ok = 0;
                break;
            }
        }
        TEST_ASSERT(diag_ok,
            "diagonal elements are real and non-negative");
    }

    covariance_matrix_free(cov);
    signal_matrix_free(signals);
}

/**
 * Test 8: After covariance computation, verify diagonal loading was applied.
 * Diagonal elements should be >= some minimum (EPSILON * factor).
 */
static void test_diagonal_loading(void)
{
    printf("\n[TEST] test_diagonal_loading\n");

    int M = 6;
    int N = 64;

    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    if (signals == NULL) {
        TEST_ASSERT(0, "signal_matrix_alloc for diagonal loading test");
        return;
    }

    /* Fill with very weak noise to test that diagonal loading prevents singularity */
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            double re = ((double)rand() / RAND_MAX - 0.5) * 1e-8;
            double im = ((double)rand() / RAND_MAX - 0.5) * 1e-8;
            signals->data[m + n * M] = re + I * im;
        }
    }

    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    if (cov == NULL) {
        signal_matrix_free(signals);
        TEST_ASSERT(0, "covariance_matrix_alloc for diagonal loading test");
        return;
    }

    music_status_t status = compute_covariance_matrix(signals, false, cov);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "compute_covariance_matrix (weak signal) returns MUSIC_SUCCESS");

    if (status == MUSIC_SUCCESS) {
        /* Diagonal loading of 1e-10 was added per code review.
         * Diagonal elements should be >= 1e-10 (the loading value).
         */
        int loading_ok = 1;
        double min_diag = 1e30;
        for (int i = 0; i < M; i++) {
            double diag_real = creal(cov->data[i + i * M]);
            if (diag_real < min_diag) min_diag = diag_real;
        }
        /* With 1e-10 diagonal loading, diagonal should be at least ~1e-10 */
        loading_ok = (min_diag >= 1e-12);
        TEST_ASSERT(loading_ok,
            "diagonal elements >= 1e-12 (diagonal loading applied)");
    }

    covariance_matrix_free(cov);
    signal_matrix_free(signals);
}

/* ============================================================================
 * Eigendecomposition Tests
 * ============================================================================ */

/**
 * Test 9: Eigenvalue ordering and non-negativity.
 */
static void test_eigenvalue_ordering(void)
{
    printf("\n[TEST] test_eigenvalue_ordering\n");

    int M = 6;
    int N = 256;

    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    if (signals == NULL) {
        TEST_ASSERT(0, "signal_matrix_alloc for eigendecomp test");
        return;
    }

    /* Generate synthetic signal with noise */
    double radius = ARRAY_RADIUS_M;
    double k = 2.0 * M_PI / WAVELENGTH_M;
    double theta = 60.0 * M_PI / 180.0;

    for (int n = 0; n < N; n++) {
        cdouble_t base = cexp(I * 2.0 * M_PI * 0.05 * n);
        for (int m = 0; m < M; m++) {
            double phi_m = 2.0 * M_PI * m / (double)M;
            double phase = k * radius * cos(theta - phi_m);
            cdouble_t steering = cexp(I * phase);
            /* Signal + noise */
            double noise_re = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            double noise_im = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            signals->data[m + n * M] = base * steering + (noise_re + I * noise_im);
        }
    }

    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    eigendecomp_result_t* eig = eigendecomp_result_alloc(M);

    if (cov == NULL || eig == NULL) {
        TEST_ASSERT(0, "allocation for eigendecomp test");
        if (cov) covariance_matrix_free(cov);
        if (eig) eigendecomp_result_free(eig);
        signal_matrix_free(signals);
        return;
    }

    music_status_t status = compute_covariance_matrix(signals, false, cov);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "compute_covariance_matrix for eigendecomp test");

    if (status == MUSIC_SUCCESS) {
        status = compute_eigendecomposition(cov, eig);
        TEST_ASSERT(status == MUSIC_SUCCESS,
            "compute_eigendecomposition returns MUSIC_SUCCESS");

        if (status == MUSIC_SUCCESS) {
            /* Verify eigenvalues are in descending order */
            int descending_ok = 1;
            for (int i = 0; i < M - 1; i++) {
                if (eig->eigenvalues[i] < eig->eigenvalues[i + 1] - 1e-12) {
                    descending_ok = 0;
                    break;
                }
            }
            TEST_ASSERT(descending_ok,
                "eigenvalues sorted in descending order");

            /* Verify no negative eigenvalues (clamped to 0 per code review fix) */
            int nonneg_ok = 1;
            for (int i = 0; i < M; i++) {
                if (eig->eigenvalues[i] < -1e-12) {
                    nonneg_ok = 0;
                    break;
                }
            }
            TEST_ASSERT(nonneg_ok,
                "all eigenvalues are non-negative (clamped)");

            /* The largest eigenvalue should be significantly bigger (signal subspace) */
            if (M >= 2) {
                double ratio = eig->eigenvalues[0] /
                    (eig->eigenvalues[M - 1] + 1e-30);
                TEST_ASSERT(ratio > 2.0,
                    "largest eigenvalue >> smallest (signal vs noise subspace)");
            }
        }
    }

    eigendecomp_result_free(eig);
    covariance_matrix_free(cov);
    signal_matrix_free(signals);
}

/* ============================================================================
 * CFAR 2D Tests
 * ============================================================================ */

/**
 * Test 10: CFAR mask creation, verify dimensions and guard cells excluded.
 */
static void test_cfar_mask_creation(void)
{
    printf("\n[TEST] test_cfar_mask_creation\n");

    cfar_config_t config;
    music_status_t status = cfar_get_default_config(&config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "cfar_get_default_config returns MUSIC_SUCCESS");

    cfar_training_mask_t mask;
    memset(&mask, 0, sizeof(mask));

    status = cfar_create_training_mask(&config, &mask);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "cfar_create_training_mask returns MUSIC_SUCCESS");

    if (status == MUSIC_SUCCESS) {
        TEST_ASSERT(mask.window_height > 0, "mask window_height > 0");
        TEST_ASSERT(mask.window_width > 0, "mask window_width > 0");
        TEST_ASSERT(mask.n_training_cells > 0, "mask has training cells > 0");
        TEST_ASSERT(mask.mask != NULL, "mask data pointer is not NULL");

        /* The center cell (CUT) should be excluded (false) */
        int center_row = mask.window_height / 2;
        int center_col = mask.window_width / 2;
        bool center_val = mask.mask[center_row * mask.window_width + center_col];
        TEST_ASSERT(center_val == false,
            "center cell (CUT) is excluded from mask (false)");

        /* Training cells should be non-zero count and less than total cells */
        int total_cells = mask.window_height * mask.window_width;
        TEST_ASSERT(mask.n_training_cells < total_cells,
            "training cells < total window cells (guard + CUT excluded)");
    }

    cfar_free_training_mask(&mask);
}

/**
 * Test 11: CFAR detection on spectrum with a known target.
 */
static void test_cfar_detection(void)
{
    printf("\n[TEST] test_cfar_detection\n");

    cfar_config_t config;
    cfar_get_default_config(&config);

    cfar_training_mask_t mask;
    memset(&mask, 0, sizeof(mask));
    music_status_t status = cfar_create_training_mask(&config, &mask);
    if (status != MUSIC_SUCCESS) {
        TEST_ASSERT(0, "cfar_create_training_mask for detection test");
        return;
    }

    /* Create a 36x1 spectrum (azimuth only, 10-deg resolution) with a strong target */
    int height = 36;
    int width = 1;
    double* spectrum = (double*)calloc(height * width, sizeof(double));
    if (spectrum == NULL) {
        TEST_ASSERT(0, "spectrum allocation for CFAR detection test");
        cfar_free_training_mask(&mask);
        return;
    }

    /* Background noise floor at -20 dB, target at 0 dB at index 18 (180 deg) */
    for (int i = 0; i < height * width; i++) {
        spectrum[i] = -20.0;
    }
    spectrum[18] = 0.0;  /* Strong target, 20 dB above noise */

    cfar_result_t* result = cfar_result_alloc(height, width);
    if (result == NULL) {
        TEST_ASSERT(0, "cfar_result_alloc for detection test");
        free(spectrum);
        cfar_free_training_mask(&mask);
        return;
    }

    status = cfar_detect_2d(spectrum, height, width, &config, &mask, result);
    TEST_ASSERT(status == MUSIC_SUCCESS,
        "cfar_detect_2d returns MUSIC_SUCCESS");

    if (status == MUSIC_SUCCESS) {
        TEST_ASSERT(result->n_detections >= 1,
            "CFAR detected at least 1 target in spectrum");

        /* The target cell should be detected */
        if (result->detections != NULL) {
            bool target_detected = result->detections[18 * width + 0];
            TEST_ASSERT(target_detected,
                "target at index 18 is detected by CFAR");
        }
    }

    cfar_result_free(result);
    free(spectrum);
    cfar_free_training_mask(&mask);
}

/* ============================================================================
 * State Machine Tests
 * ============================================================================ */

/**
 * Test 12: State machine init, verify initial state.
 */
static void test_sm_init(void)
{
    printf("\n[TEST] test_sm_init\n");

    state_machine_t sm;
    memset(&sm, 0, sizeof(sm));

    int ret = sm_init(&sm);
    TEST_ASSERT(ret == 0, "sm_init returns 0 (success)");
    TEST_ASSERT(sm.current_state == STATE_HW_INIT,
        "initial state is STATE_HW_INIT");
    TEST_ASSERT(sm.retry_count == 0,
        "initial retry_count is 0");
    TEST_ASSERT(sm.last_error == ERR_NONE,
        "initial last_error is ERR_NONE");
    TEST_ASSERT(sm.error_log_count == 0,
        "initial error_log_count is 0");
}

/**
 * Test 13: State name strings for all states.
 */
static void test_sm_state_names(void)
{
    printf("\n[TEST] test_sm_state_names\n");

    const char* name;

    name = sm_state_name(STATE_HW_INIT);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "STATE_HW_INIT has a non-empty name");

    name = sm_state_name(STATE_SDR_TEST);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "STATE_SDR_TEST has a non-empty name");

    name = sm_state_name(STATE_ALGO_TEST);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "STATE_ALGO_TEST has a non-empty name");

    name = sm_state_name(STATE_LIVE);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "STATE_LIVE has a non-empty name");

    name = sm_state_name(STATE_ERROR_DIAG);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "STATE_ERROR_DIAG has a non-empty name");
}

/**
 * Test 14: Error name strings for selected error codes.
 */
static void test_sm_error_names(void)
{
    printf("\n[TEST] test_sm_error_names\n");

    const char* name;

    name = sm_error_name(ERR_NONE);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "ERR_NONE has a non-empty name");

    name = sm_error_name(ERR_SDR_USB_NOT_FOUND);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "ERR_SDR_USB_NOT_FOUND has a non-empty name");

    name = sm_error_name(ERR_ALGO_MUSIC_FAIL);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "ERR_ALGO_MUSIC_FAIL has a non-empty name");

    name = sm_error_name(ERR_UNKNOWN);
    TEST_ASSERT(name != NULL && strlen(name) > 0,
        "ERR_UNKNOWN has a non-empty name");
}

/**
 * Test 15: Force state to STATE_LIVE, verify sm_is_operational.
 */
static void test_sm_force_state(void)
{
    printf("\n[TEST] test_sm_force_state\n");

    state_machine_t sm;
    memset(&sm, 0, sizeof(sm));
    sm_init(&sm);

    TEST_ASSERT(sm_is_operational(&sm) == false,
        "not operational in STATE_HW_INIT");

    sm_force_state(&sm, STATE_LIVE);
    TEST_ASSERT(sm.current_state == STATE_LIVE,
        "state is STATE_LIVE after force");
    TEST_ASSERT(sm_is_operational(&sm) == true,
        "sm_is_operational returns true in STATE_LIVE");

    sm_force_state(&sm, STATE_ERROR_DIAG);
    TEST_ASSERT(sm.current_state == STATE_ERROR_DIAG,
        "state is STATE_ERROR_DIAG after force");
    TEST_ASSERT(sm_is_operational(&sm) == false,
        "sm_is_operational returns false in STATE_ERROR_DIAG");
}

/**
 * Test 16: Generate report, verify non-empty string.
 */
static void test_sm_generate_report(void)
{
    printf("\n[TEST] test_sm_generate_report\n");

    state_machine_t sm;
    memset(&sm, 0, sizeof(sm));
    sm_init(&sm);

    char buffer[1024];
    memset(buffer, 0, sizeof(buffer));

    sm_generate_report(&sm, buffer, sizeof(buffer));
    TEST_ASSERT(strlen(buffer) > 0,
        "sm_generate_report produces non-empty string");
    TEST_ASSERT(strlen(buffer) < sizeof(buffer) - 1,
        "report does not overflow buffer");
}

/* ============================================================================
 * Protocol Detection Tests
 * ============================================================================ */

/**
 * Test 17: LTE detector init and free.
 */
static void test_lte_detector_init(void)
{
    printf("\n[TEST] test_lte_detector_init\n");

    lte_detector_t detector;
    memset(&detector, 0, sizeof(detector));

    int ret = lte_detector_init(&detector, 30.72e6);
    TEST_ASSERT(ret == 0, "lte_detector_init returns 0 (success)");
    TEST_ASSERT(detector.sample_rate == 30.72e6,
        "LTE detector sample_rate is 30.72 MHz");
    TEST_ASSERT(detector.fft_size > 0,
        "LTE detector fft_size > 0");

    /* Verify PSS sequences were generated (first element of sequence 0 should be non-zero) */
    cdouble_t pss0_first = detector.pss_sequences[0][0];
    TEST_ASSERT(cabs(pss0_first) > 1e-15,
        "PSS sequence 0 has non-zero first element (Zadoff-Chu generated)");

    lte_detector_free(&detector);
    TEST_ASSERT(1, "lte_detector_free completes without crash");
}

/**
 * Test 18: P25 detector init and free.
 */
static void test_p25_detector_init(void)
{
    printf("\n[TEST] test_p25_detector_init\n");

    p25_detector_t detector;
    memset(&detector, 0, sizeof(detector));

    int ret = p25_detector_init(&detector, 1.0e6);
    TEST_ASSERT(ret == 0, "p25_detector_init returns 0 (success)");
    TEST_ASSERT(detector.sample_rate == 1.0e6,
        "P25 detector sample_rate is 1 MHz");
    TEST_ASSERT(detector.symbol_rate == P25_SYMBOL_RATE,
        "P25 detector symbol_rate is 4800");
    TEST_ASSERT(detector.samples_per_symbol > 0,
        "P25 detector samples_per_symbol > 0");

    /* Verify frame sync pattern is loaded (48 bits, should have both 0s and 1s) */
    int sync_ones = 0;
    int sync_zeros = 0;
    for (int i = 0; i < P25_FRAME_SYNC_LENGTH; i++) {
        if (detector.frame_sync[i] == 1) sync_ones++;
        else if (detector.frame_sync[i] == 0) sync_zeros++;
    }
    TEST_ASSERT(sync_ones > 0 && sync_zeros > 0,
        "P25 frame sync pattern contains both 0s and 1s");
    TEST_ASSERT(sync_ones + sync_zeros == P25_FRAME_SYNC_LENGTH,
        "P25 frame sync pattern is exactly 48 bits");

    p25_detector_free(&detector);
    TEST_ASSERT(1, "p25_detector_free completes without crash");
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    printf("==========================================================\n");
    printf("  MQP System Module Unit Tests\n");
    printf("  Modules: noise_reduction, beamforming, covariance,\n");
    printf("           eigendecomp, cfar_2d, state_machine,\n");
    printf("           lte_detector, p25_detector\n");
    printf("==========================================================\n");

    /* Noise Reduction */
    test_noise_floor_estimator();
    test_wiener_gain();
    test_noise_reducer_init_free();

    /* Beamforming */
    test_conventional_beamformer();
    test_beam_pattern();
    test_beamforming_config_default();

    /* Covariance */
    test_covariance_computation();
    test_diagonal_loading();

    /* Eigendecomposition */
    test_eigenvalue_ordering();

    /* CFAR 2D */
    test_cfar_mask_creation();
    test_cfar_detection();

    /* State Machine */
    test_sm_init();
    test_sm_state_names();
    test_sm_error_names();
    test_sm_force_state();
    test_sm_generate_report();

    /* Protocol Detection */
    test_lte_detector_init();
    test_p25_detector_init();

    /* Summary */
    printf("\n==========================================================\n");
    printf("  TEST SUMMARY\n");
    printf("==========================================================\n");
    printf("  Total:  %d\n", test_passed + test_failed);
    printf("  Passed: %d\n", test_passed);
    printf("  Failed: %d\n", test_failed);
    printf("  Result: %s\n", test_failed == 0 ? "ALL PASSED" : "SOME FAILED");
    printf("==========================================================\n");

    return test_failed > 0 ? 1 : 0;
}
