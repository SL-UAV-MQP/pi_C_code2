/**
 * @file steering_vector.c
 * @brief UCA steering vector computation
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <math.h>
#include <complex.h>
#include <string.h>

/* ============================================================================
 * Steering Vector Computation
 * ============================================================================ */

music_status_t compute_steering_vector(
    const array_geometry_t* geometry,
    double azimuth_deg,
    double elevation_deg,
    double wavelength_m,
    bool include_antenna_pattern,
    cdouble_t* steering_vec)
{
    // Computes phase-shifted steering vector with optional antenna pattern

    if (!geometry || !steering_vec) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Convert to radians
    double az_rad = azimuth_deg * M_PI / 180.0;
    double el_rad = elevation_deg * M_PI / 180.0;

    // Wave number
    double k = 2.0 * M_PI / wavelength_m;

    // DOA unit vector
    double doa[3];
    doa[0] = cos(el_rad) * cos(az_rad);  // X (East)
    doa[1] = cos(el_rad) * sin(az_rad);  // Y (North)
    doa[2] = sin(el_rad);                 // Z (Up)

    // Compute steering vector with phase shifts
    for (int m = 0; m < geometry->num_elements; m++) {
        // Dot product: r_m · doa
        double dot_product =
            geometry->positions[m].x * doa[0] +
            geometry->positions[m].y * doa[1] +
            geometry->positions[m].z * doa[2];

        // Phase shift
        double phase = k * dot_product;

        // Base steering vector element
        steering_vec[m] = cexp(I * phase);

        // Apply directional antenna gain
        if (include_antenna_pattern) {
            double gain = directional_antenna_gain(
                azimuth_deg,
                geometry->orientations_deg[m],
                ANTENNA_BEAMWIDTH_DEG,
                ANTENNA_FB_RATIO_DB
            );
            steering_vec[m] *= gain;
        }
    }

    // Normalize
    music_status_t status = normalize_complex_vector(steering_vec, geometry->num_elements);
    if (status != MUSIC_SUCCESS) {
        return status;
    }

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Complex Vector Utilities
 * ============================================================================ */

double complex_vector_norm(const cdouble_t* vec, int length) {
    double sum_sq = 0.0;
    for (int i = 0; i < length; i++) {
        double magnitude = cabs(vec[i]);
        sum_sq += magnitude * magnitude;
    }
    return sqrt(sum_sq);
}

music_status_t normalize_complex_vector(cdouble_t* vec, int length) {
    if (!vec || length <= 0) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    double norm = complex_vector_norm(vec, length);

    if (norm < EPSILON) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    for (int i = 0; i < length; i++) {
        vec[i] /= norm;
    }

    return MUSIC_SUCCESS;
}
