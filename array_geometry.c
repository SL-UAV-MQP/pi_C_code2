/**
 * @file array_geometry.c
 * @brief UCA antenna array geometry implementation
 */

#include "array_geometry.h"
#include "music_config.h"
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * Array Geometry Initialization
 * ============================================================================ */

music_status_t array_geometry_init(
    array_geometry_t* geometry,
    int num_elements,
    double radius_m,
    double wavelength_m)
{
    // Validates inputs and computes UCA antenna positions and orientations

    if (!geometry) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (num_elements != NUM_ELEMENTS) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    if (radius_m <= 0 || wavelength_m <= 0) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    // Initialize structure fields
    geometry->num_elements = num_elements;
    geometry->radius_m = radius_m;
    geometry->wavelength_m = wavelength_m;

    // Compute antenna positions and orientations
    for (int i = 0; i < num_elements; i++) {
        double angle_rad = 2.0 * M_PI * i / num_elements;

        // Position (planar UCA)
        geometry->positions[i].x = radius_m * cos(angle_rad);
        geometry->positions[i].y = radius_m * sin(angle_rad);
        geometry->positions[i].z = 0.0;

        // Orientation (pointing toward center)
        double orientation = fmod(angle_rad * 180.0 / M_PI + 180.0, 360.0);
        geometry->orientations_deg[i] = orientation;
    }

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Accessor Functions
 * ============================================================================ */

music_status_t array_geometry_get_position(
    const array_geometry_t* geometry,
    int antenna_idx,
    position_3d_t* position)
{
    if (!geometry || !position) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (antenna_idx < 0 || antenna_idx >= geometry->num_elements) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    *position = geometry->positions[antenna_idx];
    return MUSIC_SUCCESS;
}

music_status_t array_geometry_get_orientation(
    const array_geometry_t* geometry,
    int antenna_idx,
    double* orientation_deg)
{
    if (!geometry || !orientation_deg) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (antenna_idx < 0 || antenna_idx >= geometry->num_elements) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    *orientation_deg = geometry->orientations_deg[antenna_idx];
    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Directional Antenna Gain Pattern
 * ============================================================================ */

double directional_antenna_gain(
    double azimuth_deg,
    double antenna_pointing_deg,
    double beamwidth_deg,
    double fb_ratio_db)
{
    // Gaussian beam pattern with front-to-back ratio

    double angle_diff = fabs(azimuth_deg - antenna_pointing_deg);

    // Wrap to [-180, 180]
    if (angle_diff > 180.0) {
        angle_diff = 360.0 - angle_diff;
    }

    // Gaussian sigma from HPBW
    double sigma = beamwidth_deg / (2.0 * sqrt(2.0 * log(2.0)));
    double gain;

    if (fabs(angle_diff) > 90.0) {
        // Back lobe
        double fb_linear = pow(10.0, -fb_ratio_db / 20.0);
        double angle_from_back = angle_diff - 180.0;
        gain = fb_linear * exp(-(angle_from_back * angle_from_back) / (2.0 * sigma * sigma));
    } else {
        // Main lobe
        gain = exp(-(angle_diff * angle_diff) / (2.0 * sigma * sigma));
    }

    return gain;
}

/* ============================================================================
 * Angle Utilities
 * ============================================================================ */

double angle_wrap_180(double angle_deg) {
    angle_deg = fmod(angle_deg + 180.0, 360.0);
    if (angle_deg < 0) {
        angle_deg += 360.0;
    }
    return angle_deg - 180.0;
}

double angle_difference(double angle1_deg, double angle2_deg) {
    double diff = fabs(angle1_deg - angle2_deg);
    if (diff > 180.0) {
        diff = 360.0 - diff;
    }
    return diff;
}

/* ============================================================================
 * Debug Helpers
 * ============================================================================ */

void array_geometry_print(const array_geometry_t* geometry) {
    if (!geometry) {
        printf("ERROR: NULL geometry\n");
        return;
    }

    printf("=== UCA Array Geometry ===\n");
    printf("Elements: %d\n", geometry->num_elements);
    printf("Radius: %.4f m (%.2f λ)\n",
           geometry->radius_m,
           geometry->radius_m / geometry->wavelength_m);
    printf("Wavelength: %.4f m\n\n", geometry->wavelength_m);

    printf("Antenna Positions and Orientations:\n");
    for (int i = 0; i < geometry->num_elements; i++) {
        printf("  Ant %d: pos=(%.4f, %.4f, %.4f) m, pointing=%.1f°\n",
               i,
               geometry->positions[i].x,
               geometry->positions[i].y,
               geometry->positions[i].z,
               geometry->orientations_deg[i]);
    }
    printf("\n");
}
