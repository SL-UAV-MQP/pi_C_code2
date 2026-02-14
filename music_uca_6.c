/**
 * @file music_uca_6.c
 * @brief Main MUSIC-UCA-6 algorithm pipeline
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ============================================================================
 * Main MUSIC AOA Estimation Pipeline
 * ============================================================================ */

music_status_t music_estimate_aoa(
    const signal_matrix_t* signals,
    const array_geometry_t* geometry,
    int num_sources,
    bool use_fb_averaging,
    music_spectrum_t* spectrum,
    detected_sources_t* detected)
{
    // estimate_aoa
    // End-to-end MUSIC pipeline with profiling

    if (!signals || !geometry || !spectrum || !detected) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Validate input dimensions
    if (signals->M != geometry->num_elements) {
        fprintf(stderr, "ERROR: Expected %d antennas, got %d\n",
                geometry->num_elements, signals->M);
        return MUSIC_ERROR_INVALID_SIZE;
    }

    if (num_sources < 1 || num_sources >= signals->M) {
        fprintf(stderr, "ERROR: Invalid num_sources=%d (must be 1 to %d)\n",
                num_sources, signals->M - 1);
        return MUSIC_ERROR_INVALID_CONFIG;
    }

#if MUSIC_PROFILE
    clock_t t_start = clock();
    clock_t t_stage;
#endif

    if (MUSIC_VERBOSE) {
        printf("Starting MUSIC AOA estimation:\n");
        printf("  Antennas: %d\n", signals->M);
        printf("  Snapshots: %d\n", signals->N);
        printf("  Sources: %d\n", num_sources);
        printf("  Forward-backward averaging: %s\n", use_fb_averaging ? "ON" : "OFF");
        printf("\n");
    }

    // Step 1: Compute covariance matrix
    covariance_matrix_t* covariance = covariance_matrix_alloc(signals->M);
    if (!covariance) {
        return MUSIC_ERROR_MEMORY;
    }

#if MUSIC_PROFILE
    t_stage = clock();
#endif

    music_status_t status = compute_covariance_matrix(signals, use_fb_averaging, covariance);
    if (status != MUSIC_SUCCESS) {
        covariance_matrix_free(covariance);
        return status;
    }

#if MUSIC_PROFILE
    double t_cov = (double)(clock() - t_stage) / CLOCKS_PER_SEC * 1000.0;
    printf("[PROFILE] Covariance: %.2f ms\n", t_cov);
#endif

    // Step 2: Eigendecomposition
    eigendecomp_result_t* eigendecomp = eigendecomp_result_alloc(signals->M);
    if (!eigendecomp) {
        covariance_matrix_free(covariance);
        return MUSIC_ERROR_MEMORY;
    }

#if MUSIC_PROFILE
    t_stage = clock();
#endif

    status = compute_eigendecomposition(covariance, eigendecomp);
    covariance_matrix_free(covariance);  // Free after use
    if (status != MUSIC_SUCCESS) {
        eigendecomp_result_free(eigendecomp);
        return status;
    }

#if MUSIC_PROFILE
    double t_eigen = (double)(clock() - t_stage) / CLOCKS_PER_SEC * 1000.0;
    printf("[PROFILE] Eigendecomposition: %.2f ms\n", t_eigen);
#endif

    // Step 3: Compute MUSIC spectrum
#if MUSIC_PROFILE
    t_stage = clock();
#endif

    status = compute_music_spectrum(
        eigendecomp,
        geometry,
        num_sources,
        AZIMUTH_MIN_DEG,
        AZIMUTH_MAX_DEG,
        AZIMUTH_RESOLUTION_DEG,
        geometry->wavelength_m,
        spectrum
    );

    if (status != MUSIC_SUCCESS) {
        eigendecomp_result_free(eigendecomp);
        return status;
    }

#if MUSIC_PROFILE
    double t_spectrum = (double)(clock() - t_stage) / CLOCKS_PER_SEC * 1000.0;
    printf("[PROFILE] MUSIC spectrum: %.2f ms (%d bins)\n",
           t_spectrum, spectrum->num_azimuth_bins);
#endif

    // Step 4: Find peaks
#if MUSIC_PROFILE
    t_stage = clock();
#endif

    status = find_peaks(
        spectrum,
        num_sources,
        MIN_PEAK_HEIGHT_DB,
        MIN_PEAK_SEPARATION_DEG,
        detected
    );

    eigendecomp_result_free(eigendecomp);  // Free after use

    if (status != MUSIC_SUCCESS) {
        if (status == MUSIC_ERROR_NO_PEAKS && MUSIC_VERBOSE) {
            fprintf(stderr, "WARNING: No peaks found in MUSIC spectrum\n");
        }
        return status;
    }

#if MUSIC_PROFILE
    double t_peaks = (double)(clock() - t_stage) / CLOCKS_PER_SEC * 1000.0;
    double t_total = (double)(clock() - t_start) / CLOCKS_PER_SEC * 1000.0;
    printf("[PROFILE] Peak detection: %.2f ms\n", t_peaks);
    printf("[PROFILE] TOTAL: %.2f ms\n\n", t_total);

    if (t_total > MAX_PROCESSING_TIME_MS) {
        fprintf(stderr, "WARNING: Processing time %.2f ms exceeds target %d ms\n",
                t_total, MAX_PROCESSING_TIME_MS);
    }
#endif

    if (MUSIC_VERBOSE) {
        printf("MUSIC AOA estimation complete:\n");
        printf("  Detected %d sources\n\n", detected->num_sources);
        print_detected_sources(detected);
    }

    return MUSIC_SUCCESS;
}
