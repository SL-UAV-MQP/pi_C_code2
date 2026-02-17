/**
 * @file beamforming.c
 * @brief Beamforming Algorithms for UCA 6-Element Array
 *
 * Implements conventional (Delay-and-Sum) and adaptive (MVDR/Capon)
 * beamforming for spatial filtering and interference rejection.
 *
 * Performance targets: Weight computation <1ms, Beamforming <5ms
 */

#include "beamforming.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>
#include <lapacke.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Speed of light */
#define SPEED_OF_LIGHT 299792458.0

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Compute UCA steering vector
 */
static void compute_uca_steering_vector(
    int num_elements,
    double radius,
    double frequency,
    double azimuth_deg,
    cdouble_t* steering_vector
) {
    double wavelength = SPEED_OF_LIGHT / frequency;
    double k = 2.0 * M_PI / wavelength;
    double azimuth_rad = azimuth_deg * M_PI / 180.0;

    // UCA element positions
    for (int m = 0; m < num_elements; m++) {
        double phi_m = 2.0 * M_PI * m / num_elements;
        double x_m = radius * cos(phi_m);
        double y_m = radius * sin(phi_m);

        // Phase shift for plane wave from azimuth direction
        double phase = k * (x_m * cos(azimuth_rad) + y_m * sin(azimuth_rad));
        steering_vector[m] = cexp(I * phase);
    }

    // Normalize
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
 * @brief Compute Hermitian matrix product: R = (1/N) * X * X^H
 */
static void compute_covariance_matrix(
    const cdouble_t* signals,
    int M,
    int N,
    cdouble_t* R
) {
    // R = (1/N) * signals * signals^H
    // signals is M x N (column-major)
    // R is M x M (column-major)

    // cblas_zgemm expects const void* pointing to double[2] (complex double)
    cdouble_t alpha = (1.0 / N) + 0.0 * I;
    cdouble_t beta = 0.0 + 0.0 * I;

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                M, M, N,
                &alpha, signals, M,
                signals, M,
                &beta, R, M);
}

/**
 * @brief Invert complex matrix using LAPACK
 */
static music_status_t invert_matrix(
    cdouble_t* matrix,
    int n,
    cdouble_t* inverse
) {
    // Copy matrix to inverse buffer
    memcpy(inverse, matrix, n * n * sizeof(cdouble_t));

    // LU factorization
    lapack_int* ipiv = (lapack_int*)malloc(n * sizeof(lapack_int));
    if (!ipiv) {
        return MUSIC_ERROR_MEMORY;
    }

    lapack_int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,
                                     (lapack_complex_double*)inverse, n, ipiv);

    if (info != 0) {
        free(ipiv);
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    // Matrix inversion
    info = LAPACKE_zgetri(LAPACK_COL_MAJOR, n,
                          (lapack_complex_double*)inverse, n, ipiv);

    free(ipiv);

    if (info != 0) {
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Configuration Helpers
 * ============================================================================ */

void beamforming_config_default(beamforming_config_t* config) {
    config->num_elements = BF_NUM_ELEMENTS;
    config->radius = BF_UCA_RADIUS;
    config->center_frequency = BF_CENTER_FREQUENCY;
    config->sample_rate = 61.44e6;
    config->antenna_beamwidth = 60.0;
    config->antenna_fb_ratio = 15.0;
    config->diagonal_loading = BF_DIAGONAL_LOADING;
    config->use_forward_backward = true;
}

/* ============================================================================
 * Conventional Beamforming
 * ============================================================================ */

music_status_t conventional_beamformer_init(
    conventional_beamformer_t* bf,
    const beamforming_config_t* config
) {
    if (!bf || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    bf->config = *config;
    bf->steering_calculator = NULL;  // Would initialize UCASteeringVector here

    return MUSIC_SUCCESS;
}

void conventional_beamformer_free(conventional_beamformer_t* bf) {
    if (bf && bf->steering_calculator) {
        // Free steering calculator if allocated
        bf->steering_calculator = NULL;
    }
}

music_status_t conventional_compute_weights(
    const conventional_beamformer_t* bf,
    double azimuth,
    double elevation,
    bool normalize,
    cdouble_t* weights
) {
    // compute_weights
    if (!bf || !weights) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Get steering vector for look direction
    cdouble_t* steering = (cdouble_t*)malloc(bf->config.num_elements * sizeof(cdouble_t));
    if (!steering) {
        return MUSIC_ERROR_MEMORY;
    }

    compute_uca_steering_vector(bf->config.num_elements, bf->config.radius,
                                bf->config.center_frequency, azimuth, steering);

    // Conventional weights: w = a* (conjugate)
    for (int i = 0; i < bf->config.num_elements; i++) {
        weights[i] = conj(steering[i]);
    }

    // Normalize if requested
    if (normalize) {
        double norm = 0.0;
        for (int i = 0; i < bf->config.num_elements; i++) {
            norm += creal(weights[i] * conj(weights[i]));
        }
        norm = sqrt(norm);

        for (int i = 0; i < bf->config.num_elements; i++) {
            weights[i] /= norm;
        }
    }

    free(steering);
    return MUSIC_SUCCESS;
}

music_status_t conventional_apply_beamformer(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    const cdouble_t* weights,
    cdouble_t* output
) {
    // apply_beamformer
    // Output: y = w^H * x
    if (!signals || !weights || !output) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Manual beamforming: y[n] = w^H * x[n] = sum_m conj(weights[m]) * signals[m + n*M]
    // Note: cblas_zgemv with CblasConjTrans would double-conjugate the weights,
    // so we use a direct loop instead (M=6 elements, negligible overhead).
    for (int n = 0; n < num_snapshots; n++) {
        cdouble_t sum = 0.0;
        for (int m = 0; m < num_elements; m++) {
            sum += conj(weights[m]) * signals[m + n * num_elements];
        }
        output[n] = sum;
    }

    return MUSIC_SUCCESS;
}

music_status_t conventional_steer_beam(
    const conventional_beamformer_t* bf,
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    double azimuth,
    double elevation,
    cdouble_t* output
) {
    // steer_beam
    cdouble_t* weights = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));
    if (!weights) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = conventional_compute_weights(bf, azimuth, elevation,
                                                          true, weights);
    if (status != MUSIC_SUCCESS) {
        free(weights);
        return status;
    }

    status = conventional_apply_beamformer(signals, num_elements, num_snapshots,
                                           weights, output);

    free(weights);
    return status;
}

music_status_t conventional_compute_pattern(
    const conventional_beamformer_t* bf,
    double azimuth_look,
    double azimuth_min,
    double azimuth_max,
    double resolution,
    beam_pattern_t* pattern
) {
    // compute_beam_pattern
    if (!bf || !pattern) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute weights for look direction
    cdouble_t* weights = (cdouble_t*)malloc(bf->config.num_elements * sizeof(cdouble_t));
    if (!weights) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = conventional_compute_weights(bf, azimuth_look, 0.0,
                                                          true, weights);
    if (status != MUSIC_SUCCESS) {
        free(weights);
        return status;
    }

    // Create azimuth grid
    int num_points = (int)((azimuth_max - azimuth_min) / resolution) + 1;
    if (num_points > pattern->num_points) {
        free(weights);
        return MUSIC_ERROR_INVALID_SIZE;
    }

    cdouble_t* steering = (cdouble_t*)malloc(bf->config.num_elements * sizeof(cdouble_t));
    if (!steering) {
        free(weights);
        return MUSIC_ERROR_MEMORY;
    }

    double max_response = 0.0;

    // Compute response for each direction
    for (int i = 0; i < num_points; i++) {
        double az = azimuth_min + i * resolution;
        pattern->azimuths[i] = az;

        // Get steering vector
        compute_uca_steering_vector(bf->config.num_elements, bf->config.radius,
                                    bf->config.center_frequency, az, steering);

        // Beam response: |w^H * a|^2
        cdouble_t response = 0.0;
        for (int j = 0; j < bf->config.num_elements; j++) {
            response += conj(weights[j]) * steering[j];
        }

        double response_mag = cabs(response);
        double response_power = response_mag * response_mag;

        if (response_power > max_response) {
            max_response = response_power;
        }

        pattern->pattern_db[i] = response_power;
    }

    // Normalize to dB
    for (int i = 0; i < num_points; i++) {
        pattern->pattern_db[i] = 10.0 * log10(pattern->pattern_db[i] / (max_response + 1e-12));
    }

    pattern->num_points = num_points;
    pattern->look_direction = azimuth_look;

    free(weights);
    free(steering);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Adaptive Beamforming (MVDR/Capon)
 * ============================================================================ */

music_status_t adaptive_beamformer_init(
    adaptive_beamformer_t* bf,
    const beamforming_config_t* config
) {
    if (!bf || !config) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    bf->config = *config;
    bf->steering_calculator = NULL;

    return MUSIC_SUCCESS;
}

void adaptive_beamformer_free(adaptive_beamformer_t* bf) {
    if (bf && bf->steering_calculator) {
        bf->steering_calculator = NULL;
    }
}

music_status_t adaptive_estimate_covariance(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    bool use_forward_backward,
    cdouble_t* R
) {
    // estimate_covariance
    if (!signals || !R) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Forward covariance: R_forward = (1/N) * X * X^H
    compute_covariance_matrix(signals, num_elements, num_snapshots, R);

    if (use_forward_backward) {
        // Backward covariance using exchange matrix J
        cdouble_t* R_backward = (cdouble_t*)malloc(num_elements * num_elements *
                                                    sizeof(cdouble_t));
        if (!R_backward) {
            return MUSIC_ERROR_MEMORY;
        }

        // Exchange matrix operation: J @ R_forward.conj() @ J
        // R is column-major (built by cblas_zgemm): element (i,j) at index [i + j*M]
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < num_elements; j++) {
                int ii = num_elements - 1 - i;
                int jj = num_elements - 1 - j;
                R_backward[i + j * num_elements] = conj(R[ii + jj * num_elements]);
            }
        }

        // Average: R = 0.5 * (R_forward + R_backward)
        for (int i = 0; i < num_elements * num_elements; i++) {
            R[i] = 0.5 * (R[i] + R_backward[i]);
        }

        free(R_backward);
    }

    return MUSIC_SUCCESS;
}

music_status_t adaptive_compute_mvdr_weights(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double azimuth,
    double elevation,
    double diagonal_loading,
    cdouble_t* weights
) {
    // compute_mvdr_weights
    // MVDR: w = (R + eI)^(-1) * a / (a^H * (R + eI)^(-1) * a)
    if (!bf || !R || !weights) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (diagonal_loading <= 0.0) {
        diagonal_loading = bf->config.diagonal_loading;
    }

    // Apply diagonal loading: R_loaded = R + e * tr(R)/M * I
    cdouble_t* R_loaded = (cdouble_t*)malloc(num_elements * num_elements *
                                              sizeof(cdouble_t));
    if (!R_loaded) {
        return MUSIC_ERROR_MEMORY;
    }

    // Compute trace
    double trace = 0.0;
    for (int i = 0; i < num_elements; i++) {
        trace += creal(R[i * num_elements + i]);
    }

    double loading = diagonal_loading * trace / num_elements;

    // Add loading to diagonal
    for (int i = 0; i < num_elements * num_elements; i++) {
        R_loaded[i] = R[i];
    }
    for (int i = 0; i < num_elements; i++) {
        R_loaded[i * num_elements + i] += loading;
    }

    // Invert R_loaded
    cdouble_t* R_inv = (cdouble_t*)malloc(num_elements * num_elements *
                                           sizeof(cdouble_t));
    if (!R_inv) {
        free(R_loaded);
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = invert_matrix(R_loaded, num_elements, R_inv);
    if (status != MUSIC_SUCCESS) {
        free(R_loaded);
        free(R_inv);
        return status;
    }

    // Get steering vector
    cdouble_t* a = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));
    if (!a) {
        free(R_loaded);
        free(R_inv);
        return MUSIC_ERROR_MEMORY;
    }

    compute_uca_steering_vector(num_elements, bf->config.radius,
                                bf->config.center_frequency, azimuth, a);

    // Compute R_inv * a
    cdouble_t* R_inv_a = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));
    if (!R_inv_a) {
        free(R_loaded);
        free(R_inv);
        free(a);
        return MUSIC_ERROR_MEMORY;
    }

    cdouble_t alpha = 1.0;
    cdouble_t beta = 0.0;
    cblas_zgemv(CblasColMajor, CblasNoTrans,
                num_elements, num_elements,
                &alpha, R_inv, num_elements,
                a, 1,
                &beta, R_inv_a, 1);

    // Compute a^H * R_inv * a
    cdouble_t denominator = 0.0;
    for (int i = 0; i < num_elements; i++) {
        denominator += conj(a[i]) * R_inv_a[i];
    }

    // MVDR weights: w = R_inv_a / denominator
    for (int i = 0; i < num_elements; i++) {
        weights[i] = R_inv_a[i] / denominator;
    }

    // Cleanup
    free(R_loaded);
    free(R_inv);
    free(a);
    free(R_inv_a);

    return MUSIC_SUCCESS;
}

music_status_t adaptive_compute_mpdr_weights(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double desired_azimuth,
    const double* interference_azimuths,
    int num_interferers,
    cdouble_t* weights
) {
    /*
     * LCMV beamforming: w = R^{-1} C (C^H R^{-1} C)^{-1} f
     *
     * Constraint matrix C (num_elements x num_constraints, column-major):
     *   Column 0: a(desired_azimuth)   -- unit gain
     *   Column 1..K: a(interference_i) -- null response
     *
     * Response vector f: [1, 0, 0, ..., 0]^T
     *
     * R is the spatial covariance matrix with diagonal loading applied.
     */
    if (!bf || !R || !weights) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = num_elements;
    int num_constraints = 1 + num_interferers;

    /* Clamp constraints to M-1 to avoid over-determined system */
    if (num_constraints > M) {
        num_constraints = M;
    }
    int actual_interferers = num_constraints - 1;

    /* Build constraint matrix C and response vector f */
    cdouble_t* C = (cdouble_t*)calloc(M * num_constraints, sizeof(cdouble_t));
    cdouble_t* f = (cdouble_t*)calloc(num_constraints, sizeof(cdouble_t));

    if (!C || !f) {
        if (C) free(C);
        if (f) free(f);
        return MUSIC_ERROR_MEMORY;
    }

    /* First constraint: desired direction (unit gain) */
    compute_uca_steering_vector(M, bf->config.radius,
                                bf->config.center_frequency, desired_azimuth,
                                &C[0]);  /* Column 0 */
    f[0] = 1.0 + 0.0 * I;

    /* Interference constraints (nulls) */
    if (interference_azimuths) {
        for (int i = 0; i < actual_interferers; i++) {
            compute_uca_steering_vector(M, bf->config.radius,
                                        bf->config.center_frequency,
                                        interference_azimuths[i],
                                        &C[(i + 1) * M]);
            f[i + 1] = 0.0 + 0.0 * I;
        }
    }

    int Nc = 1 + actual_interferers;

    /* -----------------------------------------------------------------
     * Apply diagonal loading to R: R_loaded = R + eps * tr(R)/M * I
     * ----------------------------------------------------------------- */
    cdouble_t* R_loaded = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
    if (!R_loaded) {
        free(C); free(f);
        return MUSIC_ERROR_MEMORY;
    }

    double trace = 0.0;
    for (int i = 0; i < M; i++) {
        trace += creal(R[i * M + i]);
    }
    double loading = bf->config.diagonal_loading * trace / M;

    memcpy(R_loaded, R, M * M * sizeof(cdouble_t));
    for (int i = 0; i < M; i++) {
        R_loaded[i * M + i] += loading;
    }

    /* -----------------------------------------------------------------
     * Invert R_loaded -> R_inv  (M x M)
     * ----------------------------------------------------------------- */
    cdouble_t* R_inv = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
    if (!R_inv) {
        free(C); free(f); free(R_loaded);
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = invert_matrix(R_loaded, M, R_inv);
    if (status != MUSIC_SUCCESS) {
        /* Fallback to MVDR if inversion fails */
        free(C); free(f); free(R_loaded); free(R_inv);
        return adaptive_compute_mvdr_weights(bf, R, M, desired_azimuth, 0.0,
                                             bf->config.diagonal_loading, weights);
    }

    /* -----------------------------------------------------------------
     * Compute R_inv_C = R_inv * C   (M x Nc)
     * ----------------------------------------------------------------- */
    cdouble_t* R_inv_C = (cdouble_t*)calloc(M * Nc, sizeof(cdouble_t));
    if (!R_inv_C) {
        free(C); free(f); free(R_loaded); free(R_inv);
        return MUSIC_ERROR_MEMORY;
    }

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
     * Compute G = C^H * R_inv * C = C^H * R_inv_C   (Nc x Nc)
     * ----------------------------------------------------------------- */
    cdouble_t* G = (cdouble_t*)calloc(Nc * Nc, sizeof(cdouble_t));
    if (!G) {
        free(C); free(f); free(R_loaded); free(R_inv); free(R_inv_C);
        return MUSIC_ERROR_MEMORY;
    }

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
     * Invert G -> G_inv   (Nc x Nc, small matrix)
     * ----------------------------------------------------------------- */
    cdouble_t* G_inv = (cdouble_t*)calloc(Nc * Nc, sizeof(cdouble_t));
    if (!G_inv) {
        free(C); free(f); free(R_loaded); free(R_inv); free(R_inv_C); free(G);
        return MUSIC_ERROR_MEMORY;
    }

    status = invert_matrix(G, Nc, G_inv);
    if (status != MUSIC_SUCCESS) {
        /* Fallback to MVDR if small matrix inversion fails */
        free(C); free(f); free(R_loaded); free(R_inv);
        free(R_inv_C); free(G); free(G_inv);
        return adaptive_compute_mvdr_weights(bf, R, M, desired_azimuth, 0.0,
                                             bf->config.diagonal_loading, weights);
    }

    /* -----------------------------------------------------------------
     * Compute temp = G_inv * f   (Nc x 1)
     * ----------------------------------------------------------------- */
    cdouble_t* temp = (cdouble_t*)calloc(Nc, sizeof(cdouble_t));
    if (!temp) {
        free(C); free(f); free(R_loaded); free(R_inv);
        free(R_inv_C); free(G); free(G_inv);
        return MUSIC_ERROR_MEMORY;
    }

    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    Nc, Nc,
                    &alpha_blas, G_inv, Nc,
                    f, 1,
                    &beta_blas, temp, 1);
    }

    /* -----------------------------------------------------------------
     * Compute weights = R_inv_C * temp   (M x 1)
     * w = R^{-1} * C * (C^H * R^{-1} * C)^{-1} * f
     * ----------------------------------------------------------------- */
    {
        cdouble_t alpha_blas = 1.0 + 0.0 * I;
        cdouble_t beta_blas = 0.0 + 0.0 * I;
        cblas_zgemv(CblasColMajor, CblasNoTrans,
                    M, Nc,
                    &alpha_blas, R_inv_C, M,
                    temp, 1,
                    &beta_blas, weights, 1);
    }

    /* Cleanup */
    free(C);
    free(f);
    free(R_loaded);
    free(R_inv);
    free(R_inv_C);
    free(G);
    free(G_inv);
    free(temp);

    return MUSIC_SUCCESS;
}

music_status_t adaptive_apply_beamformer(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    const cdouble_t* weights,
    cdouble_t* output
) {
    // apply_beamformer - same as conventional
    return conventional_apply_beamformer(signals, num_elements, num_snapshots,
                                         weights, output);
}

music_status_t adaptive_steer_beam(
    const adaptive_beamformer_t* bf,
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    double azimuth,
    double elevation,
    beamforming_result_t* result
) {
    // steer_adaptive_beam
    if (!bf || !signals || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Estimate covariance
    cdouble_t* R = (cdouble_t*)malloc(num_elements * num_elements * sizeof(cdouble_t));
    if (!R) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = adaptive_estimate_covariance(signals, num_elements,
                                                          num_snapshots,
                                                          bf->config.use_forward_backward, R);
    if (status != MUSIC_SUCCESS) {
        free(R);
        return status;
    }

    // Compute MVDR weights
    status = adaptive_compute_mvdr_weights(bf, R, num_elements, azimuth, elevation,
                                           0.0, result->weights);
    if (status != MUSIC_SUCCESS) {
        free(R);
        return status;
    }

    // Apply beamformer
    status = adaptive_apply_beamformer(signals, num_elements, num_snapshots,
                                       result->weights, result->output);

    free(R);

    strcpy(result->method, "MVDR");
    result->num_samples = num_snapshots;

    return status;
}

music_status_t adaptive_compute_pattern(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double azimuth_look,
    double azimuth_min,
    double azimuth_max,
    double resolution,
    beam_pattern_t* pattern
) {
    /* Stub: adaptive beam pattern computation not yet implemented */
    (void)bf; (void)R; (void)num_elements; (void)azimuth_look;
    (void)azimuth_min; (void)azimuth_max; (void)resolution; (void)pattern;
    return MUSIC_ERROR_INVALID_CONFIG;
}

/* ============================================================================
 * Unified Beamforming Processor
 * ============================================================================ */

music_status_t beamforming_processor_init(
    beamforming_processor_t* processor,
    const beamforming_config_t* config
) {
    if (!processor) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (config) {
        processor->config = *config;
    } else {
        beamforming_config_default(&processor->config);
    }

    music_status_t status = conventional_beamformer_init(&processor->conventional,
                                                          &processor->config);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    status = adaptive_beamformer_init(&processor->adaptive, &processor->config);
    if (status != MUSIC_SUCCESS) {
        conventional_beamformer_free(&processor->conventional);
        return status;
    }

    return MUSIC_SUCCESS;
}

void beamforming_processor_free(beamforming_processor_t* processor) {
    if (processor) {
        conventional_beamformer_free(&processor->conventional);
        adaptive_beamformer_free(&processor->adaptive);
    }
}

music_status_t beamforming_process_signals(
    beamforming_processor_t* processor,
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    beamforming_method_t method,
    double desired_azimuth,
    const double* interference_azimuths,
    int num_interferers,
    beamforming_result_t* result
) {
    // process_signals
    if (!processor || !signals || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    music_status_t status = MUSIC_SUCCESS;

    if (method == BF_METHOD_CONVENTIONAL) {
        strcpy(result->method, "Conventional");
        status = conventional_steer_beam(&processor->conventional, signals,
                                         num_elements, num_snapshots,
                                         desired_azimuth, 0.0, result->output);
    } else if (method == BF_METHOD_ADAPTIVE) {
        status = adaptive_steer_beam(&processor->adaptive, signals,
                                     num_elements, num_snapshots,
                                     desired_azimuth, 0.0, result);
    } else if (method == BF_METHOD_MPDR && interference_azimuths && num_interferers > 0) {
        strcpy(result->method, "MPDR");

        /* Estimate covariance for LCMV */
        cdouble_t* R = (cdouble_t*)malloc(num_elements * num_elements * sizeof(cdouble_t));
        if (!R) {
            return MUSIC_ERROR_MEMORY;
        }

        status = adaptive_estimate_covariance(signals, num_elements, num_snapshots,
                                               processor->config.use_forward_backward, R);
        if (status != MUSIC_SUCCESS) {
            free(R);
            return status;
        }

        /* Compute LCMV weights */
        status = adaptive_compute_mpdr_weights(&processor->adaptive, R, num_elements,
                                                desired_azimuth, interference_azimuths,
                                                num_interferers, result->weights);
        if (status != MUSIC_SUCCESS) {
            free(R);
            return status;
        }

        /* Apply beamformer */
        status = adaptive_apply_beamformer(signals, num_elements, num_snapshots,
                                           result->weights, result->output);
        free(R);
    } else if (method == BF_METHOD_MPDR) {
        /* BUG-18 fix: MPDR without interference data - fall back to MVDR */
        status = adaptive_steer_beam(&processor->adaptive, signals,
                                     num_elements, num_snapshots,
                                     desired_azimuth, 0.0, result);
    } else {
        /* Unknown method */
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    result->num_samples = num_snapshots;

    // Estimate SINR
    beamforming_estimate_sinr(signals, result->output, result->weights,
                              num_elements, num_snapshots, &result->sinr_db);

    return status;
}

music_status_t beamforming_estimate_sinr(
    const cdouble_t* input_signals,
    const cdouble_t* output_signal,
    const cdouble_t* weights,
    int num_elements,
    int num_snapshots,
    double* sinr_db
) {
    // estimate_sinr
    if (!input_signals || !output_signal || !weights || !sinr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Input power (average across elements)
    double input_power = 0.0;
    for (int m = 0; m < num_elements; m++) {
        for (int n = 0; n < num_snapshots; n++) {
            cdouble_t val = input_signals[m + n * num_elements];
            input_power += cabs(val) * cabs(val);
        }
    }
    input_power /= (num_elements * num_snapshots);

    // Output power
    double output_power = 0.0;
    for (int n = 0; n < num_snapshots; n++) {
        output_power += cabs(output_signal[n]) * cabs(output_signal[n]);
    }
    output_power /= num_snapshots;

    // Array gain (power ratio)
    double array_gain = output_power / ((input_power / num_elements) + 1e-12);

    *sinr_db = 10.0 * log10(array_gain + 1e-12);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Memory Management
 * ============================================================================ */

beam_pattern_t* beam_pattern_alloc(int num_points) {
    beam_pattern_t* pattern = (beam_pattern_t*)malloc(sizeof(beam_pattern_t));
    if (!pattern) {
        return NULL;
    }

    pattern->num_points = num_points;
    pattern->azimuths = (double*)malloc(num_points * sizeof(double));
    pattern->pattern_db = (double*)malloc(num_points * sizeof(double));

    if (!pattern->azimuths || !pattern->pattern_db) {
        beam_pattern_free(pattern);
        return NULL;
    }

    return pattern;
}

void beam_pattern_free(beam_pattern_t* pattern) {
    if (pattern) {
        if (pattern->azimuths) free(pattern->azimuths);
        if (pattern->pattern_db) free(pattern->pattern_db);
        free(pattern);
    }
}

beamforming_result_t* beamforming_result_alloc(int num_elements, int num_samples) {
    beamforming_result_t* result = (beamforming_result_t*)malloc(sizeof(beamforming_result_t));
    if (!result) {
        return NULL;
    }

    result->num_samples = num_samples;
    result->output = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
    result->weights = (cdouble_t*)malloc(num_elements * sizeof(cdouble_t));

    if (!result->output || !result->weights) {
        beamforming_result_free(result);
        return NULL;
    }

    /* BUG-20 fix: initialize method and sinr to avoid uninitialized reads */
    result->method[0] = '\0';
    result->sinr_db = 0.0;

    return result;
}

void beamforming_result_free(beamforming_result_t* result) {
    if (result) {
        if (result->output) free(result->output);
        if (result->weights) free(result->weights);
        free(result);
    }
}

music_status_t create_beamformer(
    beamforming_processor_t* processor,
    double frequency,
    double sample_rate,
    double diagonal_loading
) {
    // create_beamformer
    beamforming_config_t config;
    beamforming_config_default(&config);

    config.center_frequency = frequency;
    config.sample_rate = sample_rate;
    config.diagonal_loading = diagonal_loading;

    return beamforming_processor_init(processor, &config);
}
