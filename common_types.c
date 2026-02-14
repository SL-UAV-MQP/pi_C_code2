/**
 * @file common_types.c
 * @brief Memory management for common data structures
 *
 * Implements allocation and deallocation functions for MUSIC algorithm
 * data structures.
 */

#include "common_types.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Signal Matrix Memory Management
 * ============================================================================ */

signal_matrix_t* signal_matrix_alloc(int M, int N) {
    if (M <= 0 || N <= 0 || M > MAX_ANTENNAS || N > MAX_SNAPSHOTS) {
        return NULL;
    }

    signal_matrix_t* mat = (signal_matrix_t*)malloc(sizeof(signal_matrix_t));
    if (!mat) {
        return NULL;
    }

    mat->M = M;
    mat->N = N;
    mat->data = (cdouble_t*)calloc(M * N, sizeof(cdouble_t));

    if (!mat->data) {
        free(mat);
        return NULL;
    }

    return mat;
}

void signal_matrix_free(signal_matrix_t* mat) {
    if (mat) {
        if (mat->data) {
            free(mat->data);
        }
        free(mat);
    }
}

/* ============================================================================
 * Covariance Matrix Memory Management
 * ============================================================================ */

covariance_matrix_t* covariance_matrix_alloc(int M) {
    if (M <= 0 || M > MAX_ANTENNAS) {
        return NULL;
    }

    covariance_matrix_t* mat = (covariance_matrix_t*)malloc(sizeof(covariance_matrix_t));
    if (!mat) {
        return NULL;
    }

    mat->M = M;
    mat->data = (cdouble_t*)calloc(M * M, sizeof(cdouble_t));

    if (!mat->data) {
        free(mat);
        return NULL;
    }

    return mat;
}

void covariance_matrix_free(covariance_matrix_t* mat) {
    if (mat) {
        if (mat->data) {
            free(mat->data);
        }
        free(mat);
    }
}

/* ============================================================================
 * Eigendecomposition Result Memory Management
 * ============================================================================ */

eigendecomp_result_t* eigendecomp_result_alloc(int M) {
    if (M <= 0 || M > MAX_ANTENNAS) {
        return NULL;
    }

    eigendecomp_result_t* result = (eigendecomp_result_t*)malloc(sizeof(eigendecomp_result_t));
    if (!result) {
        return NULL;
    }

    result->M = M;
    result->eigenvalues = (double*)calloc(M, sizeof(double));
    result->eigenvectors = (cdouble_t*)calloc(M * M, sizeof(cdouble_t));

    if (!result->eigenvalues || !result->eigenvectors) {
        if (result->eigenvalues) free(result->eigenvalues);
        if (result->eigenvectors) free(result->eigenvectors);
        free(result);
        return NULL;
    }

    return result;
}

void eigendecomp_result_free(eigendecomp_result_t* result) {
    if (result) {
        if (result->eigenvalues) {
            free(result->eigenvalues);
        }
        if (result->eigenvectors) {
            free(result->eigenvectors);
        }
        free(result);
    }
}

/* ============================================================================
 * MUSIC Spectrum Memory Management
 * ============================================================================ */

music_spectrum_t* music_spectrum_alloc(int num_azimuth_bins, int num_eigenvalues) {
    if (num_azimuth_bins <= 0 || num_eigenvalues <= 0 ||
        num_azimuth_bins > MAX_AZIMUTH_BINS || num_eigenvalues > MAX_ANTENNAS) {
        return NULL;
    }

    music_spectrum_t* spectrum = (music_spectrum_t*)malloc(sizeof(music_spectrum_t));
    if (!spectrum) {
        return NULL;
    }

    spectrum->num_azimuth_bins = num_azimuth_bins;
    spectrum->num_eigenvalues = num_eigenvalues;
    spectrum->num_sources = 0;

    spectrum->azimuth_deg = (double*)calloc(num_azimuth_bins, sizeof(double));
    spectrum->spectrum_db = (double*)calloc(num_azimuth_bins, sizeof(double));
    spectrum->eigenvalues = (double*)calloc(num_eigenvalues, sizeof(double));

    if (!spectrum->azimuth_deg || !spectrum->spectrum_db || !spectrum->eigenvalues) {
        if (spectrum->azimuth_deg) free(spectrum->azimuth_deg);
        if (spectrum->spectrum_db) free(spectrum->spectrum_db);
        if (spectrum->eigenvalues) free(spectrum->eigenvalues);
        free(spectrum);
        return NULL;
    }

    return spectrum;
}

void music_spectrum_free(music_spectrum_t* spectrum) {
    if (spectrum) {
        if (spectrum->azimuth_deg) free(spectrum->azimuth_deg);
        if (spectrum->spectrum_db) free(spectrum->spectrum_db);
        if (spectrum->eigenvalues) free(spectrum->eigenvalues);
        free(spectrum);
    }
}

/* ============================================================================
 * CFAR Result Memory Management
 * ============================================================================ */

cfar_result_t* cfar_result_alloc(int height, int width) {
    if (height <= 0 || width <= 0) {
        return NULL;
    }

    cfar_result_t* result = (cfar_result_t*)malloc(sizeof(cfar_result_t));
    if (!result) {
        return NULL;
    }

    result->height = height;
    result->width = width;
    result->n_detections = 0;
    result->peak_snr_db = 0.0;
    result->noise_floor_db = 0.0;

    int size = height * width;
    result->detections = (bool*)calloc(size, sizeof(bool));
    result->threshold_db = (double*)calloc(size, sizeof(double));
    result->snr_db = (double*)calloc(size, sizeof(double));

    if (!result->detections || !result->threshold_db || !result->snr_db) {
        if (result->detections) free(result->detections);
        if (result->threshold_db) free(result->threshold_db);
        if (result->snr_db) free(result->snr_db);
        free(result);
        return NULL;
    }

    return result;
}

void cfar_result_free(cfar_result_t* result) {
    if (result) {
        if (result->detections) free(result->detections);
        if (result->threshold_db) free(result->threshold_db);
        if (result->snr_db) free(result->snr_db);
        free(result);
    }
}
