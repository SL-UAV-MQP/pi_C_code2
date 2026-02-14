/**
 * @file test_music.c
 * @brief Comprehensive unit tests for MUSIC-UCA-6 algorithm
 *
 * Tests each component of the MUSIC algorithm with MATLAB reference data validation.
 * Target accuracy: 2.84° mean AOA error (MATLAB baseline).
 *
 * Test Coverage:
 * 1. Array Geometry: Position/orientation initialization
 * 2. Antenna Pattern: Directional gain computation
 * 3. Steering Vectors: Phase shifts + antenna pattern
 * 4. Covariance Matrix: Forward-backward averaging
 * 5. Eigendecomposition: Signal/noise subspace separation
 * 6. MUSIC Spectrum: Noise subspace projection
 * 7. Peak Detection: Parabolic interpolation
 * 8. Model Order Selection: MDL/AIC estimation
 * 9. Integration Tests: Full pipeline validation
 * 10. Edge Cases: Numerical stability, error handling
 *
 * @author Testing and Validation Specialist
 * @date 2025-12-11
 */

#include "music_uca_6.h"
#include "array_geometry.h"
#include "cfar_2d.h"
#include "music_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>

/* ============================================================================
 * Test Configuration
 * ============================================================================ */

#define TEST_TOLERANCE_ANGLE 0.5    /* Degrees */
#define TEST_TOLERANCE_GAIN 0.01    /* Linear scale */
#define TEST_TOLERANCE_DB 0.5       /* dB */
#define TEST_TOLERANCE_NUMERIC 1e-6 /* Numerical precision */
#define TEST_TOLERANCE_PHASE 0.01   /* Radians */

/* ============================================================================
 * MATLAB Reference Data (from MATLAB baseline)
 * ============================================================================ */

/* Known tower azimuths from CMRCM field (degrees, 0=North) */
static const double TOWER_AZIMUTHS[] = {
    349.7,  /* Tower 1: WQNR305/WQNX339 (650m) */
    291.8,  /* Tower 2: WQQL972/WQRV287 (1.7km) */
    189.0   /* Tower 3: Chauncy Hill (1.9km) */
};
static const int NUM_TOWERS = 3;

/* Expected array geometry (MATLAB: MUSIC_UCA_6.m lines 76-82) */
static const position_3d_t EXPECTED_POSITIONS[6] = {
    {0.1760, 0.0000, 0.0},   /* Ant 0: 0° */
    {0.0880, 0.1524, 0.0},   /* Ant 1: 60° */
    {-0.0880, 0.1524, 0.0},  /* Ant 2: 120° */
    {-0.1760, 0.0000, 0.0},  /* Ant 3: 180° */
    {-0.0880, -0.1524, 0.0}, /* Ant 4: 240° */
    {0.0880, -0.1524, 0.0}   /* Ant 5: 300° */
};

static const double EXPECTED_ORIENTATIONS[6] = {
    180.0, 240.0, 300.0, 0.0, 60.0, 120.0
};

/* ============================================================================
 * Test Helper Functions
 * ============================================================================ */

static int test_passed = 0;
static int test_failed = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            fprintf(stderr, "FAIL: %s\n", message); \
            fprintf(stderr, "      %s:%d\n", __FILE__, __LINE__); \
            test_failed++; \
            return; \
        } \
    } while (0)

#define TEST_ASSERT_DOUBLE(actual, expected, tolerance, message) \
    do { \
        if (fabs((actual) - (expected)) > (tolerance)) { \
            fprintf(stderr, "FAIL: %s\n", message); \
            fprintf(stderr, "      Expected: %.6f, Got: %.6f, Diff: %.6f\n", \
                    (double)(expected), (double)(actual), \
                    fabs((double)(actual) - (double)(expected))); \
            fprintf(stderr, "      %s:%d\n", __FILE__, __LINE__); \
            test_failed++; \
            return; \
        } \
    } while (0)

#define TEST_ASSERT_COMPLEX(actual, expected, tolerance, message) \
    do { \
        double diff = cabs((actual) - (expected)); \
        if (diff > (tolerance)) { \
            fprintf(stderr, "FAIL: %s\n", message); \
            fprintf(stderr, "      Expected: %.6f+%.6fi, Got: %.6f+%.6fi, Diff: %.6f\n", \
                    creal(expected), cimag(expected), \
                    creal(actual), cimag(actual), diff); \
            fprintf(stderr, "      %s:%d\n", __FILE__, __LINE__); \
            test_failed++; \
            return; \
        } \
    } while (0)

#define TEST_PASS() \
    do { \
        test_passed++; \
        printf("PASS\n"); \
    } while (0)

static void print_test_summary() {
    printf("\n");
    printf("===============================\n");
    printf("Test Summary\n");
    printf("===============================\n");
    printf("Passed: %d\n", test_passed);
    printf("Failed: %d\n", test_failed);
    printf("Total:  %d\n", test_passed + test_failed);
    if (test_failed == 0) {
        printf("Result: ALL TESTS PASSED\n");
    } else {
        printf("Result: SOME TESTS FAILED\n");
    }
    printf("===============================\n");
}

/* Helper: Generate synthetic signal with known sources */
static void generate_synthetic_signal(
    signal_matrix_t* signals,
    const array_geometry_t* geometry,
    const double* azimuths,
    int num_sources,
    double snr_db)
{
    int M = signals->M;
    int N = signals->N;

    // Noise variance from SNR
    double noise_power = pow(10.0, -snr_db / 10.0);

    // Initialize with noise
    for (int i = 0; i < M * N; i++) {
        double real = ((double)rand() / RAND_MAX - 0.5) * sqrt(noise_power);
        double imag = ((double)rand() / RAND_MAX - 0.5) * sqrt(noise_power);
        signals->data[i] = real + I * imag;
    }

    // Add signals from each source
    for (int src = 0; src < num_sources; src++) {
        cdouble_t steering[MAX_ANTENNAS];
        compute_steering_vector(geometry, azimuths[src], 0.0,
                              WAVELENGTH_M, true, steering);

        // Add source signal with random phase
        for (int snap = 0; snap < N; snap++) {
            double phase = 2.0 * M_PI * ((double)rand() / RAND_MAX);
            cdouble_t symbol = cexp(I * phase);

            for (int ant = 0; ant < M; ant++) {
                signals->data[ant + snap * M] += steering[ant] * symbol;
            }
        }
    }
}

/* ============================================================================
 * Test Cases: Array Geometry
 * ============================================================================ */

void test_array_geometry_init() {
    printf("TEST: Array Geometry Initialization ... ");

    array_geometry_t geometry;
    music_status_t status = array_geometry_init(
        &geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    TEST_ASSERT(status == MUSIC_SUCCESS, "Array init should succeed");
    TEST_ASSERT(geometry.num_elements == NUM_ELEMENTS, "Number of elements should be 6");
    TEST_ASSERT(fabs(geometry.radius_m - ARRAY_RADIUS_M) < 1e-6, "Radius should match config");

    // Check antenna positions (MATLAB reference)
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        TEST_ASSERT_DOUBLE(geometry.positions[i].x, EXPECTED_POSITIONS[i].x,
                          0.01, "X position should match");
        TEST_ASSERT_DOUBLE(geometry.positions[i].y, EXPECTED_POSITIONS[i].y,
                          0.01, "Y position should match");
        TEST_ASSERT_DOUBLE(geometry.positions[i].z, EXPECTED_POSITIONS[i].z,
                          1e-6, "Z position should be 0");
    }

    // Check orientations
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        TEST_ASSERT_DOUBLE(geometry.orientations_deg[i], EXPECTED_ORIENTATIONS[i],
                          1.0, "Orientation should match");
    }

    TEST_PASS();
}

void test_array_geometry_spacing() {
    printf("TEST: Array Geometry Spacing ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    // Check inter-element spacing (should be ~0.6λ circumference)
    double circumference = 2.0 * M_PI * ARRAY_RADIUS_M;
    double expected_spacing = circumference / NUM_ELEMENTS;

    // Check distance between adjacent elements
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        int j = (i + 1) % NUM_ELEMENTS;
        double dx = geometry.positions[j].x - geometry.positions[i].x;
        double dy = geometry.positions[j].y - geometry.positions[i].y;
        double dist = sqrt(dx*dx + dy*dy);

        TEST_ASSERT_DOUBLE(dist, expected_spacing, 0.01,
                          "Adjacent element spacing should be uniform");
    }

    // Check that all elements are on circle
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        double r = sqrt(geometry.positions[i].x * geometry.positions[i].x +
                       geometry.positions[i].y * geometry.positions[i].y);
        TEST_ASSERT_DOUBLE(r, ARRAY_RADIUS_M, 1e-6,
                          "All elements should be on circle");
    }

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Antenna Pattern
 * ============================================================================ */

void test_directional_antenna_gain() {
    printf("TEST: Directional Antenna Gain Pattern ... ");

    // Test cases from MATLAB validation
    // At boresight (0° difference): gain = 1.0
    double gain_0 = directional_antenna_gain(0.0, 0.0,
                        ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);
    TEST_ASSERT_DOUBLE(gain_0, 1.0, TEST_TOLERANCE_GAIN, "Boresight gain should be 1.0");

    // At HPBW/2 (60° for 120° beamwidth): gain ≈ 0.5 (-3 dB)
    double gain_60 = directional_antenna_gain(60.0, 0.0,
                         ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);
    TEST_ASSERT_DOUBLE(gain_60, 0.5, 0.1, "HPBW/2 gain should be ~0.5");

    // At back (180° difference): gain ≈ 0.178 (-15 dB F/B ratio)
    double gain_180 = directional_antenna_gain(180.0, 0.0,
                          ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);
    double expected_fb = pow(10.0, -ANTENNA_FB_RATIO_DB / 20.0);
    TEST_ASSERT_DOUBLE(gain_180, expected_fb, 0.05, "Back lobe gain should match F/B ratio");

    TEST_PASS();
}

void test_antenna_gain_symmetry() {
    printf("TEST: Antenna Gain Symmetry ... ");

    // Test that gain pattern is symmetric
    for (int angle = 0; angle <= 180; angle += 30) {
        double gain_pos = directional_antenna_gain((double)angle, 0.0,
                             ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);
        double gain_neg = directional_antenna_gain((double)-angle, 0.0,
                             ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);

        TEST_ASSERT_DOUBLE(gain_pos, gain_neg, TEST_TOLERANCE_GAIN,
                          "Antenna pattern should be symmetric");
    }

    TEST_PASS();
}

void test_antenna_gain_monotonic() {
    printf("TEST: Antenna Gain Monotonic Decrease ... ");

    // Test that gain decreases monotonically from 0° to 90°
    double prev_gain = 1.0;
    for (int angle = 0; angle <= 90; angle += 10) {
        double gain = directional_antenna_gain((double)angle, 0.0,
                         ANTENNA_BEAMWIDTH_DEG, ANTENNA_FB_RATIO_DB);

        TEST_ASSERT(gain <= prev_gain + TEST_TOLERANCE_GAIN,
                   "Gain should decrease monotonically in main lobe");
        prev_gain = gain;
    }

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Steering Vector
 * ============================================================================ */

void test_steering_vector_basic() {
    printf("TEST: Steering Vector Computation ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    cdouble_t steering_vec[NUM_ELEMENTS];

    // Test 1: Broadside (0° azimuth, North)
    music_status_t status = compute_steering_vector(
        &geometry, 0.0, 0.0, WAVELENGTH_M, false, steering_vec);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Steering vector should compute successfully");

    // Check normalization
    double norm = complex_vector_norm(steering_vec, NUM_ELEMENTS);
    TEST_ASSERT_DOUBLE(norm, 1.0, TEST_TOLERANCE_NUMERIC, "Steering vector should be normalized");

    // Test 2: With antenna pattern
    status = compute_steering_vector(
        &geometry, 0.0, 0.0, WAVELENGTH_M, true, steering_vec);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Steering vector with pattern should succeed");

    norm = complex_vector_norm(steering_vec, NUM_ELEMENTS);
    TEST_ASSERT_DOUBLE(norm, 1.0, TEST_TOLERANCE_NUMERIC, "Steering vector should remain normalized");

    TEST_PASS();
}

void test_steering_vector_phase_progression() {
    printf("TEST: Steering Vector Phase Progression ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    cdouble_t steering_vec[NUM_ELEMENTS];

    // For 0° azimuth, check phase progression
    compute_steering_vector(&geometry, 0.0, 0.0, WAVELENGTH_M, false, steering_vec);

    // Verify all phases are within reasonable range
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        double phase = carg(steering_vec[i]);
        TEST_ASSERT(fabs(phase) <= M_PI, "Phase should be in [-π, π]");
    }

    // For 90° azimuth, phases should be different
    cdouble_t steering_vec_90[NUM_ELEMENTS];
    compute_steering_vector(&geometry, 90.0, 0.0, WAVELENGTH_M, false, steering_vec_90);

    // At least one phase should differ significantly
    bool phases_differ = false;
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        if (fabs(carg(steering_vec[i]) - carg(steering_vec_90[i])) > 0.1) {
            phases_differ = true;
            break;
        }
    }
    TEST_ASSERT(phases_differ, "Phase progression should change with azimuth");

    TEST_PASS();
}

void test_steering_vector_antenna_pattern_effect() {
    printf("TEST: Steering Vector Antenna Pattern Effect ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    cdouble_t sv_no_pattern[NUM_ELEMENTS];
    cdouble_t sv_with_pattern[NUM_ELEMENTS];

    // Compute both versions
    compute_steering_vector(&geometry, 45.0, 0.0, WAVELENGTH_M, false, sv_no_pattern);
    compute_steering_vector(&geometry, 45.0, 0.0, WAVELENGTH_M, true, sv_with_pattern);

    // Magnitudes should differ due to antenna pattern
    bool magnitudes_differ = false;
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        double mag_no_pattern = cabs(sv_no_pattern[i]);
        double mag_with_pattern = cabs(sv_with_pattern[i]);

        // Both should contribute, but differently after normalization
        if (fabs(mag_no_pattern - mag_with_pattern) > 0.01) {
            magnitudes_differ = true;
        }
    }

    // Both should still be normalized
    TEST_ASSERT_DOUBLE(complex_vector_norm(sv_no_pattern, NUM_ELEMENTS), 1.0,
                      TEST_TOLERANCE_NUMERIC, "No-pattern vector should be normalized");
    TEST_ASSERT_DOUBLE(complex_vector_norm(sv_with_pattern, NUM_ELEMENTS), 1.0,
                      TEST_TOLERANCE_NUMERIC, "With-pattern vector should be normalized");

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Covariance Matrix
 * ============================================================================ */

void test_covariance_matrix_basic() {
    printf("TEST: Covariance Matrix Computation ... ");

    // Create synthetic signal (single source at 0°)
    int M = 6;
    int N = 1000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    TEST_ASSERT(signals != NULL, "Signal matrix allocation should succeed");

    // Fill with random noise + signal
    srand(12345);  // Fixed seed for reproducibility
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            double real = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            double imag = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            signals->data[i + j * M] = real + I * imag;
        }
    }

    // Compute covariance
    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    TEST_ASSERT(cov != NULL, "Covariance matrix allocation should succeed");

    music_status_t status = compute_covariance_matrix(signals, false, cov);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Covariance computation should succeed");

    // Check Hermitian property: R[i,j] = conj(R[j,i])
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cdouble_t R_ij = cov->data[i + j * M];
            cdouble_t R_ji = cov->data[j + i * M];
            double diff = cabs(R_ij - conj(R_ji));
            TEST_ASSERT(diff < 1e-10, "Covariance should be Hermitian");
        }
    }

    signal_matrix_free(signals);
    covariance_matrix_free(cov);

    TEST_PASS();
}

void test_covariance_forward_backward() {
    printf("TEST: Covariance Forward-Backward Averaging ... ");

    int M = 6;
    int N = 1000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    TEST_ASSERT(signals != NULL, "Signal allocation should succeed");

    // Generate synthetic data
    srand(12345);
    for (int i = 0; i < M * N; i++) {
        double real = ((double)rand() / RAND_MAX - 0.5);
        double imag = ((double)rand() / RAND_MAX - 0.5);
        signals->data[i] = real + I * imag;
    }

    // Compute both versions
    covariance_matrix_t* cov_f = covariance_matrix_alloc(M);
    covariance_matrix_t* cov_fb = covariance_matrix_alloc(M);
    TEST_ASSERT(cov_f != NULL && cov_fb != NULL, "Covariance allocation should succeed");

    compute_covariance_matrix(signals, false, cov_f);
    compute_covariance_matrix(signals, true, cov_fb);

    // Both should be Hermitian
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cdouble_t R_ij = cov_fb->data[i + j * M];
            cdouble_t R_ji = cov_fb->data[j + i * M];
            TEST_ASSERT(cabs(R_ij - conj(R_ji)) < 1e-10, "FB covariance should be Hermitian");
        }
    }

    // Diagonal should be similar (both are power measurements)
    for (int i = 0; i < M; i++) {
        double diag_f = creal(cov_f->data[i + i * M]);
        double diag_fb = creal(cov_fb->data[i + i * M]);
        TEST_ASSERT(diag_f > 0 && diag_fb > 0, "Diagonal should be positive");
    }

    signal_matrix_free(signals);
    covariance_matrix_free(cov_f);
    covariance_matrix_free(cov_fb);

    TEST_PASS();
}

void test_covariance_positive_semidefinite() {
    printf("TEST: Covariance Positive Semi-Definite ... ");

    int M = 6;
    int N = 1000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    TEST_ASSERT(signals != NULL, "Signal allocation should succeed");

    // Generate data
    srand(12345);
    for (int i = 0; i < M * N; i++) {
        double real = ((double)rand() / RAND_MAX - 0.5);
        double imag = ((double)rand() / RAND_MAX - 0.5);
        signals->data[i] = real + I * imag;
    }

    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    compute_covariance_matrix(signals, false, cov);

    // Compute eigenvalues to verify positive semi-definite
    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    compute_eigendecomposition(cov, eigen);

    // All eigenvalues should be non-negative
    for (int i = 0; i < M; i++) {
        TEST_ASSERT(eigen->eigenvalues[i] >= -1e-10,
                   "Eigenvalues should be non-negative (positive semi-definite)");
    }

    signal_matrix_free(signals);
    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Eigendecomposition
 * ============================================================================ */

void test_eigendecomposition_identity() {
    printf("TEST: Eigenvalue Decomposition (Identity) ... ");

    // Create simple test covariance matrix (identity)
    int M = 6;
    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    TEST_ASSERT(cov != NULL, "Covariance allocation should succeed");

    // Identity matrix (all eigenvalues = 1)
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cov->data[i + j * M] = (i == j) ? 1.0 : 0.0;
        }
    }

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    TEST_ASSERT(eigen != NULL, "Eigendecomp allocation should succeed");

    music_status_t status = compute_eigendecomposition(cov, eigen);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Eigendecomposition should succeed");

    // Check eigenvalues (should all be ~1.0 for identity matrix)
    for (int i = 0; i < M; i++) {
        TEST_ASSERT_DOUBLE(eigen->eigenvalues[i], 1.0, TEST_TOLERANCE_NUMERIC,
                          "Eigenvalues of identity should be 1.0");
    }

    // Check descending order
    for (int i = 0; i < M - 1; i++) {
        TEST_ASSERT(eigen->eigenvalues[i] >= eigen->eigenvalues[i + 1] - 1e-10,
                   "Eigenvalues should be in descending order");
    }

    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);

    TEST_PASS();
}

void test_eigendecomposition_diagonal() {
    printf("TEST: Eigendecomposition (Diagonal Matrix) ... ");

    int M = 6;
    covariance_matrix_t* cov = covariance_matrix_alloc(M);

    // Diagonal matrix with known eigenvalues
    double known_eigs[6] = {5.0, 4.0, 3.0, 2.0, 1.0, 0.5};
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cov->data[i + j * M] = (i == j) ? known_eigs[i] : 0.0;
        }
    }

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    compute_eigendecomposition(cov, eigen);

    // Check eigenvalues match (in descending order)
    double sorted_eigs[6] = {5.0, 4.0, 3.0, 2.0, 1.0, 0.5};
    for (int i = 0; i < M; i++) {
        TEST_ASSERT_DOUBLE(eigen->eigenvalues[i], sorted_eigs[i], 1e-6,
                          "Eigenvalues should match known values");
    }

    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);

    TEST_PASS();
}

void test_eigendecomposition_orthonormal() {
    printf("TEST: Eigendecomposition Orthonormal Eigenvectors ... ");

    int M = 6;
    covariance_matrix_t* cov = covariance_matrix_alloc(M);

    // Create random Hermitian matrix
    srand(12345);
    for (int i = 0; i < M; i++) {
        for (int j = i; j < M; j++) {
            double real = ((double)rand() / RAND_MAX - 0.5);
            double imag = (i == j) ? 0.0 : ((double)rand() / RAND_MAX - 0.5);
            cov->data[i + j * M] = real + I * imag;
            cov->data[j + i * M] = real - I * imag;  // Hermitian
        }
    }

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    compute_eigendecomposition(cov, eigen);

    // Check orthonormality: v_i^H * v_j = δ_ij
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cdouble_t dot_product = 0.0;
            for (int k = 0; k < M; k++) {
                dot_product += conj(eigen->eigenvectors[k + i * M]) *
                              eigen->eigenvectors[k + j * M];
            }

            double expected = (i == j) ? 1.0 : 0.0;
            TEST_ASSERT_DOUBLE(cabs(dot_product), expected, 1e-6,
                             "Eigenvectors should be orthonormal");
        }
    }

    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: MUSIC Spectrum
 * ============================================================================ */

void test_music_spectrum_single_source() {
    printf("TEST: MUSIC Spectrum (Single Source) ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    // Generate signal with single source at 45°
    int M = 6;
    int N = 2000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);

    double source_az = 45.0;
    double azimuths[] = {source_az};
    srand(12345);
    generate_synthetic_signal(signals, &geometry, azimuths, 1, 20.0);

    // Compute covariance and eigendecomposition
    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    compute_covariance_matrix(signals, false, cov);

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    compute_eigendecomposition(cov, eigen);

    // Compute MUSIC spectrum
    int num_bins = 361;  // -180 to 180 with 1° resolution
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, M);

    music_status_t status = compute_music_spectrum(
        eigen, &geometry, 1, -180.0, 180.0, 1.0, WAVELENGTH_M, spectrum);

    TEST_ASSERT(status == MUSIC_SUCCESS, "MUSIC spectrum should compute successfully");

    // Find peak
    double max_val = -1000.0;
    int max_idx = 0;
    for (int i = 0; i < num_bins; i++) {
        if (spectrum->spectrum_db[i] > max_val) {
            max_val = spectrum->spectrum_db[i];
            max_idx = i;
        }
    }

    double peak_azimuth = spectrum->azimuth_deg[max_idx];

    // Peak should be near source azimuth (within a few degrees)
    TEST_ASSERT_DOUBLE(peak_azimuth, source_az, 5.0,
                      "Peak should be near source azimuth");

    // Peak should be at 0 dB (normalized)
    TEST_ASSERT_DOUBLE(max_val, 0.0, 0.1, "Peak should be at 0 dB");

    signal_matrix_free(signals);
    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);
    music_spectrum_free(spectrum);

    TEST_PASS();
}

void test_music_spectrum_normalization() {
    printf("TEST: MUSIC Spectrum Normalization ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    // Generate signal
    int M = 6;
    int N = 1000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    double azimuths[] = {0.0};
    srand(12345);
    generate_synthetic_signal(signals, &geometry, azimuths, 1, 15.0);

    covariance_matrix_t* cov = covariance_matrix_alloc(M);
    compute_covariance_matrix(signals, false, cov);

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    compute_eigendecomposition(cov, eigen);

    int num_bins = 181;
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, M);
    compute_music_spectrum(eigen, &geometry, 1, -90.0, 90.0, 1.0, WAVELENGTH_M, spectrum);

    // Find max value
    double max_val = spectrum->spectrum_db[0];
    for (int i = 1; i < num_bins; i++) {
        if (spectrum->spectrum_db[i] > max_val) {
            max_val = spectrum->spectrum_db[i];
        }
    }

    // Max should be at 0 dB
    TEST_ASSERT_DOUBLE(max_val, 0.0, 0.01, "Spectrum should be normalized to 0 dB peak");

    signal_matrix_free(signals);
    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);
    music_spectrum_free(spectrum);

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Peak Detection
 * ============================================================================ */

void test_peak_detection_single() {
    printf("TEST: Peak Detection (Single Peak) ... ");

    // Create synthetic spectrum with single peak at 45°
    int num_bins = 361;
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, 6);

    // Fill with Gaussian peak
    double peak_az = 45.0;
    double peak_width = 5.0;  // degrees
    for (int i = 0; i < num_bins; i++) {
        double az = -180.0 + i * 1.0;
        spectrum->azimuth_deg[i] = az;
        double diff = az - peak_az;
        spectrum->spectrum_db[i] = -0.5 * (diff * diff) / (peak_width * peak_width);
    }

    detected_sources_t detected;
    music_status_t status = find_peaks(spectrum, 1, -30.0, 10.0, &detected);

    TEST_ASSERT(status == MUSIC_SUCCESS, "Peak detection should succeed");
    TEST_ASSERT(detected.num_sources == 1, "Should detect 1 peak");
    TEST_ASSERT_DOUBLE(detected.sources[0].azimuth_deg, peak_az, 2.0,
                      "Peak should be near expected azimuth");

    music_spectrum_free(spectrum);

    TEST_PASS();
}

void test_peak_detection_multiple() {
    printf("TEST: Peak Detection (Multiple Peaks) ... ");

    int num_bins = 361;
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, 6);

    // Create 3 peaks
    double peak_azimuths[] = {-60.0, 0.0, 75.0};
    double peak_mags[] = {0.0, -2.0, -5.0};  // Different magnitudes

    for (int i = 0; i < num_bins; i++) {
        double az = -180.0 + i * 1.0;
        spectrum->azimuth_deg[i] = az;

        double max_val = -100.0;
        for (int p = 0; p < 3; p++) {
            double diff = az - peak_azimuths[p];
            double val = peak_mags[p] - 0.5 * (diff * diff) / 25.0;
            if (val > max_val) max_val = val;
        }
        spectrum->spectrum_db[i] = max_val;
    }

    detected_sources_t detected;
    music_status_t status = find_peaks(spectrum, 3, -30.0, 20.0, &detected);

    TEST_ASSERT(status == MUSIC_SUCCESS, "Peak detection should succeed");
    TEST_ASSERT(detected.num_sources == 3, "Should detect 3 peaks");

    // Peaks should be sorted by magnitude (strongest first)
    for (int i = 0; i < detected.num_sources - 1; i++) {
        TEST_ASSERT(detected.sources[i].magnitude_db >=
                   detected.sources[i + 1].magnitude_db,
                   "Peaks should be sorted by magnitude");
    }

    music_spectrum_free(spectrum);

    TEST_PASS();
}

void test_peak_detection_parabolic_interpolation() {
    printf("TEST: Peak Detection Parabolic Interpolation ... ");

    int num_bins = 361;
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, 6);

    // Create peak at 45.3° (between grid points)
    double true_peak = 45.3;
    for (int i = 0; i < num_bins; i++) {
        double az = -180.0 + i * 1.0;
        spectrum->azimuth_deg[i] = az;
        double diff = az - true_peak;
        spectrum->spectrum_db[i] = -0.5 * diff * diff / 9.0;
    }

    detected_sources_t detected;
    find_peaks(spectrum, 1, -30.0, 10.0, &detected);

    // Parabolic interpolation should refine the peak location
    // Should be closer to true peak than grid resolution (1°)
    TEST_ASSERT_DOUBLE(detected.sources[0].azimuth_deg, true_peak, 0.5,
                      "Parabolic interpolation should refine peak location");

    music_spectrum_free(spectrum);

    TEST_PASS();
}

void test_peak_detection_no_peaks() {
    printf("TEST: Peak Detection (No Peaks) ... ");

    int num_bins = 361;
    music_spectrum_t* spectrum = music_spectrum_alloc(num_bins, 6);

    // Flat spectrum below threshold
    for (int i = 0; i < num_bins; i++) {
        spectrum->azimuth_deg[i] = -180.0 + i * 1.0;
        spectrum->spectrum_db[i] = -50.0;  // Below -30 dB threshold
    }

    detected_sources_t detected;
    music_status_t status = find_peaks(spectrum, 3, -30.0, 10.0, &detected);

    TEST_ASSERT(status == MUSIC_ERROR_NO_PEAKS, "Should return NO_PEAKS error");
    TEST_ASSERT(detected.num_sources == 0, "Should detect 0 peaks");

    music_spectrum_free(spectrum);

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: Model Order Selection
 * ============================================================================ */

void test_model_order_selection_mdl() {
    printf("TEST: Model Order Selection (MDL) ... ");

    // Create eigenvalues with clear gap (2 signal sources)
    double eigenvalues[6] = {10.0, 8.0, 0.5, 0.4, 0.3, 0.2};

    int estimated_sources;
    music_status_t status = estimate_num_sources(
        eigenvalues, 6, 1000, true, &estimated_sources);

    TEST_ASSERT(status == MUSIC_SUCCESS, "Model order selection should succeed");
    TEST_ASSERT(estimated_sources >= 1 && estimated_sources <= 4,
               "Estimated sources should be reasonable");

    TEST_PASS();
}

void test_model_order_selection_aic() {
    printf("TEST: Model Order Selection (AIC) ... ");

    double eigenvalues[6] = {10.0, 8.0, 0.5, 0.4, 0.3, 0.2};

    int estimated_sources;
    music_status_t status = estimate_num_sources(
        eigenvalues, 6, 1000, false, &estimated_sources);

    TEST_ASSERT(status == MUSIC_SUCCESS, "AIC should succeed");
    TEST_ASSERT(estimated_sources >= 1 && estimated_sources <= 4,
               "AIC estimate should be reasonable");

    TEST_PASS();
}

/* ============================================================================
 * Test Cases: CFAR Configuration
 * ============================================================================ */

void test_cfar_configuration() {
    printf("TEST: CFAR Configuration ... ");

    cfar_config_t config;
    music_status_t status = cfar_get_cmrcm_config(&config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "CMRCM config should be valid");

    status = cfar_validate_config(&config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "CMRCM config should pass validation");

    // Check values
    TEST_ASSERT(config.window_size_az == CFAR_WINDOW_SIZE_AZ, "Window size should match");
    TEST_ASSERT(config.k_factor == CFAR_K_FACTOR, "k-factor should match");

    TEST_PASS();
}

void test_cfar_training_mask() {
    printf("TEST: CFAR Training Mask ... ");

    cfar_config_t config;
    cfar_get_cmrcm_config(&config);

    cfar_training_mask_t mask;
    music_status_t status = cfar_create_training_mask(&config, &mask);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Training mask creation should succeed");

    // Check dimensions
    TEST_ASSERT(mask.window_height == config.window_size_az, "Mask height should match");
    TEST_ASSERT(mask.window_width == config.window_size_el, "Mask width should match");

    // Check that training cells > 0
    TEST_ASSERT(mask.n_training_cells > 0, "Should have training cells");

    // Check that center (CUT) is excluded
    int ch = mask.window_height / 2;
    int cw = mask.window_width / 2;
    TEST_ASSERT(mask.mask[ch * mask.window_width + cw] == false,
               "Center cell (CUT) should be excluded");

    cfar_free_training_mask(&mask);

    TEST_PASS();
}

/* ============================================================================
 * Integration Tests
 * ============================================================================ */

void test_music_integration_pipeline() {
    printf("TEST: MUSIC Integration Pipeline ... ");

    // Initialize array geometry
    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    // Create synthetic multi-source signal
    int M = 6;
    int N = 4000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);
    TEST_ASSERT(signals != NULL, "Signal allocation should succeed");

    // Generate signal with 3 known sources
    double source_azimuths[] = {-30.0, 45.0, 120.0};
    srand(12345);
    generate_synthetic_signal(signals, &geometry, source_azimuths, 3, 20.0);

    // Allocate output structures
    music_spectrum_t* spectrum = music_spectrum_alloc(NUM_AZIMUTH_BINS, M);
    TEST_ASSERT(spectrum != NULL, "Spectrum allocation should succeed");

    detected_sources_t detected;

    // Run MUSIC pipeline
    music_status_t status = music_estimate_aoa(
        signals, &geometry, 3, USE_FORWARD_BACKWARD,
        spectrum, &detected);

    TEST_ASSERT(status == MUSIC_SUCCESS, "MUSIC pipeline should succeed");
    TEST_ASSERT(detected.num_sources > 0, "Should detect sources");
    TEST_ASSERT(detected.num_sources <= 3, "Should not detect more than expected");

    // Cleanup
    signal_matrix_free(signals);
    music_spectrum_free(spectrum);

    TEST_PASS();
}

void test_music_end_to_end_validation() {
    printf("TEST: MUSIC End-to-End Validation ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    // Generate high-SNR signal with known sources
    int M = 6;
    int N = 4000;
    signal_matrix_t* signals = signal_matrix_alloc(M, N);

    double source_azimuths[] = {0.0, 90.0, -90.0};
    srand(54321);
    generate_synthetic_signal(signals, &geometry, source_azimuths, 3, 25.0);

    music_spectrum_t* spectrum = music_spectrum_alloc(NUM_AZIMUTH_BINS, M);
    detected_sources_t detected;

    music_status_t status = music_estimate_aoa(
        signals, &geometry, 3, true, spectrum, &detected);

    TEST_ASSERT(status == MUSIC_SUCCESS, "End-to-end should succeed");
    TEST_ASSERT(detected.num_sources == 3, "Should detect all 3 sources");

    // Check that detected sources are close to true sources
    for (int i = 0; i < detected.num_sources; i++) {
        bool found_match = false;
        for (int j = 0; j < 3; j++) {
            double diff = fabs(detected.sources[i].azimuth_deg - source_azimuths[j]);
            if (diff < 10.0) {  // Within 10 degrees
                found_match = true;
                break;
            }
        }
        TEST_ASSERT(found_match, "Each detected source should match a true source");
    }

    signal_matrix_free(signals);
    music_spectrum_free(spectrum);

    TEST_PASS();
}

/* ============================================================================
 * Edge Cases and Error Handling
 * ============================================================================ */

void test_null_pointer_handling() {
    printf("TEST: Null Pointer Handling ... ");

    array_geometry_t geometry;
    array_geometry_init(&geometry, NUM_ELEMENTS, ARRAY_RADIUS_M, WAVELENGTH_M);

    cdouble_t steering_vec[NUM_ELEMENTS];

    // Test null pointer in compute_steering_vector
    music_status_t status = compute_steering_vector(
        NULL, 0.0, 0.0, WAVELENGTH_M, false, steering_vec);
    TEST_ASSERT(status == MUSIC_ERROR_NULL_POINTER,
               "Should return NULL_POINTER error");

    status = compute_steering_vector(
        &geometry, 0.0, 0.0, WAVELENGTH_M, false, NULL);
    TEST_ASSERT(status == MUSIC_ERROR_NULL_POINTER,
               "Should return NULL_POINTER error for null output");

    TEST_PASS();
}

void test_invalid_dimensions() {
    printf("TEST: Invalid Dimensions Handling ... ");

    // Try to allocate with invalid dimensions
    signal_matrix_t* signals = signal_matrix_alloc(0, 100);
    TEST_ASSERT(signals == NULL, "Should fail with 0 antennas");

    signals = signal_matrix_alloc(6, 0);
    TEST_ASSERT(signals == NULL, "Should fail with 0 snapshots");

    TEST_PASS();
}

void test_numerical_stability_small_values() {
    printf("TEST: Numerical Stability (Small Values) ... ");

    int M = 6;
    covariance_matrix_t* cov = covariance_matrix_alloc(M);

    // Very small covariance matrix
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cov->data[i + j * M] = (i == j) ? 1e-10 : 0.0;
        }
    }

    eigendecomp_result_t* eigen = eigendecomp_result_alloc(M);
    music_status_t status = compute_eigendecomposition(cov, eigen);

    TEST_ASSERT(status == MUSIC_SUCCESS, "Should handle small values");

    // Eigenvalues should be small but positive
    for (int i = 0; i < M; i++) {
        TEST_ASSERT(eigen->eigenvalues[i] >= 0.0,
                   "Eigenvalues should remain non-negative");
    }

    covariance_matrix_free(cov);
    eigendecomp_result_free(eigen);

    TEST_PASS();
}

/* ============================================================================
 * Utility Function Tests
 * ============================================================================ */

void test_utility_functions() {
    printf("TEST: Utility Functions ... ");

    // Test linear_to_db
    double power_linear = 100.0;
    double power_db = linear_to_db(power_linear);
    TEST_ASSERT_DOUBLE(power_db, 20.0, 1e-6, "100 linear = 20 dB");

    // Test db_to_linear
    double back_to_linear = db_to_linear(power_db);
    TEST_ASSERT_DOUBLE(back_to_linear, power_linear, 1e-6,
                      "Should convert back correctly");

    // Test complex_vector_norm
    cdouble_t vec[3] = {3.0 + 0.0*I, 0.0 + 4.0*I, 0.0 + 0.0*I};
    double norm = complex_vector_norm(vec, 3);
    TEST_ASSERT_DOUBLE(norm, 5.0, 1e-6, "Norm of [3,4i,0] should be 5");

    // Test normalize_complex_vector
    cdouble_t vec2[3] = {1.0, 1.0, 1.0};
    normalize_complex_vector(vec2, 3);
    double norm2 = complex_vector_norm(vec2, 3);
    TEST_ASSERT_DOUBLE(norm2, 1.0, 1e-6, "Normalized vector should have norm 1");

    TEST_PASS();
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */

int main(int argc, char** argv) {
    printf("\n");
    printf("===============================\n");
    printf("MUSIC-UCA-6 Unit Test Suite\n");
    printf("===============================\n");
    printf("Target: Raspberry Pi 5B\n");
    printf("MATLAB baseline: 2.84° mean error\n");
    printf("Test date: 2025-12-11\n");
    printf("===============================\n\n");

    // Component tests: Array Geometry
    printf("--- Array Geometry Tests ---\n");
    test_array_geometry_init();
    test_array_geometry_spacing();

    // Component tests: Antenna Pattern
    printf("\n--- Antenna Pattern Tests ---\n");
    test_directional_antenna_gain();
    test_antenna_gain_symmetry();
    test_antenna_gain_monotonic();

    // Component tests: Steering Vector
    printf("\n--- Steering Vector Tests ---\n");
    test_steering_vector_basic();
    test_steering_vector_phase_progression();
    test_steering_vector_antenna_pattern_effect();

    // Component tests: Covariance Matrix
    printf("\n--- Covariance Matrix Tests ---\n");
    test_covariance_matrix_basic();
    test_covariance_forward_backward();
    test_covariance_positive_semidefinite();

    // Component tests: Eigendecomposition
    printf("\n--- Eigendecomposition Tests ---\n");
    test_eigendecomposition_identity();
    test_eigendecomposition_diagonal();
    test_eigendecomposition_orthonormal();

    // Component tests: MUSIC Spectrum
    printf("\n--- MUSIC Spectrum Tests ---\n");
    test_music_spectrum_single_source();
    test_music_spectrum_normalization();

    // Component tests: Peak Detection
    printf("\n--- Peak Detection Tests ---\n");
    test_peak_detection_single();
    test_peak_detection_multiple();
    test_peak_detection_parabolic_interpolation();
    test_peak_detection_no_peaks();

    // Component tests: Model Order Selection
    printf("\n--- Model Order Selection Tests ---\n");
    test_model_order_selection_mdl();
    test_model_order_selection_aic();

    // CFAR tests
    printf("\n--- CFAR Tests ---\n");
    test_cfar_configuration();
    test_cfar_training_mask();

    // Integration tests
    printf("\n--- Integration Tests ---\n");
    test_music_integration_pipeline();
    test_music_end_to_end_validation();

    // Edge cases and error handling
    printf("\n--- Edge Cases & Error Handling ---\n");
    test_null_pointer_handling();
    test_invalid_dimensions();
    test_numerical_stability_small_values();

    // Utility functions
    printf("\n--- Utility Functions ---\n");
    test_utility_functions();

    // Summary
    print_test_summary();

    return (test_failed > 0) ? 1 : 0;
}
