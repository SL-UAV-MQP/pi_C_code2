/**
 * @file spectrum.c
 * @brief MUSIC pseudo-spectrum computation
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <math.h>
#include <complex.h>
#include <string.h>

/* BLAS function declarations */
/* ZGEMV: Matrix-vector multiply: y := alpha*A*x + beta*y */
extern void zgemv_(
    const char* trans,
    const int* m,
    const int* n,
    const cdouble_t* alpha,
    const cdouble_t* A,
    const int* lda,
    const cdouble_t* x,
    const int* incx,
    const cdouble_t* beta,
    cdouble_t* y,
    const int* incy
);

/* ============================================================================
 * MUSIC Spectrum Computation
 * ============================================================================ */

music_status_t compute_music_spectrum(
    const eigendecomp_result_t* eigendecomp,
    const array_geometry_t* geometry,
    int num_sources,
    double azimuth_min_deg,
    double azimuth_max_deg,
    double azimuth_resolution_deg,
    double wavelength_m,
    music_spectrum_t* spectrum)
{
    // Computes MUSIC pseudo-spectrum using noise subspace projection

    if (!eigendecomp || !geometry || !spectrum) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = eigendecomp->M;
    if (num_sources >= M || num_sources < 1) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Calculate number of azimuth bins
    int num_bins = (int)((azimuth_max_deg - azimuth_min_deg) / azimuth_resolution_deg) + 1;
    if (num_bins != spectrum->num_azimuth_bins) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Extract noise subspace (eigenvectors corresponding to M-K smallest eigenvalues)
    int noise_dim = M - num_sources;
    cdouble_t* noise_subspace = eigendecomp->eigenvectors + num_sources * M;  // Points to column num_sources

    // Allocate temporary steering vector
    cdouble_t* steering_vec = (cdouble_t*)malloc(M * sizeof(cdouble_t));
    cdouble_t* projection = (cdouble_t*)malloc(noise_dim * sizeof(cdouble_t));
    if (!steering_vec || !projection) {
        if (steering_vec) free(steering_vec);
        if (projection) free(projection);
        return MUSIC_ERROR_MEMORY;
    }

    // Compute spectrum for each azimuth angle
    double* spectrum_linear = (double*)malloc(num_bins * sizeof(double));
    if (!spectrum_linear) {
        free(steering_vec);
        free(projection);
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < num_bins; i++) {
        double azimuth = azimuth_min_deg + i * azimuth_resolution_deg;
        spectrum->azimuth_deg[i] = azimuth;

        // Compute steering vector
        music_status_t status = compute_steering_vector(
            geometry, azimuth, 0.0, wavelength_m, true, steering_vec);
        if (status != MUSIC_SUCCESS) {
            free(steering_vec);
            free(projection);
            free(spectrum_linear);
            return status;
        }

        // Project onto noise subspace: projection = E_n^H * a
        // Use BLAS ZGEMV: y = alpha * A^H * x + beta * y
        const char trans = 'C';  // Conjugate transpose
        cdouble_t alpha = 1.0 + 0.0*I;
        cdouble_t beta = 0.0 + 0.0*I;
        int inc = 1;

        zgemv_(&trans, &M, &noise_dim, &alpha, noise_subspace, &M,
               steering_vec, &inc, &beta, projection, &inc);

        // Compute ||projection||^2
        double norm_sq = 0.0;
        for (int j = 0; j < noise_dim; j++) {
            double mag = cabs(projection[j]);
            norm_sq += mag * mag;
        }

        // MUSIC spectrum: P(θ) = 1 / ||projection||^2
        spectrum_linear[i] = 1.0 / (norm_sq + EPSILON);
    }

    // Find maximum for normalization
    double max_val = spectrum_linear[0];
    for (int i = 1; i < num_bins; i++) {
        if (spectrum_linear[i] > max_val) {
            max_val = spectrum_linear[i];
        }
    }

    // Convert to dB (normalized)
    for (int i = 0; i < num_bins; i++) {
        spectrum->spectrum_db[i] = 10.0 * log10(spectrum_linear[i] / max_val + EPSILON);
    }

    // Copy eigenvalues to spectrum result
    spectrum->num_sources = num_sources;
    memcpy(spectrum->eigenvalues, eigendecomp->eigenvalues,
           spectrum->num_eigenvalues * sizeof(double));

    // Cleanup
    free(steering_vec);
    free(projection);
    free(spectrum_linear);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

double linear_to_db(double linear_power) {
    return 10.0 * log10(linear_power + EPSILON);
}

double db_to_linear(double db_power) {
    return pow(10.0, db_power / 10.0);
}

void print_spectrum_stats(const music_spectrum_t* spectrum) {
    if (!spectrum) {
        printf("ERROR: NULL spectrum\n");
        return;
    }

    // Find max and min
    double max_val = spectrum->spectrum_db[0];
    double min_val = spectrum->spectrum_db[0];
    int max_idx = 0;

    for (int i = 1; i < spectrum->num_azimuth_bins; i++) {
        if (spectrum->spectrum_db[i] > max_val) {
            max_val = spectrum->spectrum_db[i];
            max_idx = i;
        }
        if (spectrum->spectrum_db[i] < min_val) {
            min_val = spectrum->spectrum_db[i];
        }
    }

    printf("=== MUSIC Spectrum Statistics ===\n");
    printf("Azimuth bins: %d\n", spectrum->num_azimuth_bins);
    printf("Azimuth range: [%.1f°, %.1f°]\n",
           spectrum->azimuth_deg[0],
           spectrum->azimuth_deg[spectrum->num_azimuth_bins - 1]);
    printf("Resolution: %.2f°\n",
           spectrum->azimuth_deg[1] - spectrum->azimuth_deg[0]);
    printf("Spectrum range: [%.1f, %.1f] dB\n", min_val, max_val);
    printf("Peak at %.1f° (%.1f dB)\n",
           spectrum->azimuth_deg[max_idx], max_val);
    printf("Number of sources: %d\n", spectrum->num_sources);
    printf("\n");
}
