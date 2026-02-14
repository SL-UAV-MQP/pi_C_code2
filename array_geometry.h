/**
 * @file array_geometry.h
 * @brief UCA antenna array geometry calculations
 *
 * Computes 6-element Uniform Circular Array positions and orientations.
 * Maps to MATLAB compute_antenna_geometry() method.
 */

#ifndef ARRAY_GEOMETRY_H
#define ARRAY_GEOMETRY_H

#include "common_types.h"

/**
 * @brief Initialize UCA antenna array geometry
 *
 * Computes antenna positions in planar circular array (X-Y plane)
 * and orientation angles (pointing toward center).
 *
 * MATLAB: MUSIC_UCA_6.m lines 69-83 (compute_antenna_geometry)
 *
 * Array layout (view from above, looking down at Z=0 plane):
 *        Y (North)
 *        ^
 *        |
 *   Ant2 | Ant1
 *        |
 * -------+-------> X (East)
 *        |
 *   Ant3 | Ant6
 *        |
 *   Ant4   Ant5
 *
 * Antenna numbering: 1-6 (clockwise from East, starting at 0°)
 * Orientations: Each antenna points toward center (+180° from position angle)
 *
 * @param[out] geometry  Pointer to array_geometry_t structure to initialize
 * @param[in]  num_elements  Number of antenna elements (must be 6)
 * @param[in]  radius_m      Array radius in meters
 * @param[in]  wavelength_m  Wavelength at center frequency in meters
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Antenna positions are in meters, Z=0 (planar array)
 * @note Orientations are in degrees, 0° = North, 90° = East
 */
music_status_t array_geometry_init(
    array_geometry_t* geometry,
    int num_elements,
    double radius_m,
    double wavelength_m
);

/**
 * @brief Get antenna position by index
 *
 * @param[in] geometry  Pointer to initialized geometry
 * @param[in] antenna_idx  Antenna index (0-5)
 * @param[out] position  Output position vector
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 */
music_status_t array_geometry_get_position(
    const array_geometry_t* geometry,
    int antenna_idx,
    position_3d_t* position
);

/**
 * @brief Get antenna orientation by index
 *
 * @param[in] geometry  Pointer to initialized geometry
 * @param[in] antenna_idx  Antenna index (0-5)
 * @param[out] orientation_deg  Output orientation in degrees
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 */
music_status_t array_geometry_get_orientation(
    const array_geometry_t* geometry,
    int antenna_idx,
    double* orientation_deg
);

/**
 * @brief Print array geometry information
 *
 * Debug helper to display antenna positions and orientations.
 *
 * @param[in] geometry  Pointer to initialized geometry
 */
void array_geometry_print(const array_geometry_t* geometry);

/**
 * @brief Compute directional antenna gain pattern
 *
 * Calculates antenna gain using Gaussian beam pattern with front-to-back ratio.
 * Used to apply directional antenna characteristics to steering vectors.
 *
 * MATLAB: MUSIC_UCA_6.m lines 85-114 (directional_antenna_gain)
 *
 * Pattern model:
 * - Main lobe: Gaussian with HPBW = 120°
 * - Back lobe: Attenuated by 15 dB F/B ratio
 *
 * @param[in] azimuth_deg  Signal azimuth angle (degrees, -180 to 180)
 * @param[in] antenna_pointing_deg  Antenna pointing direction (degrees)
 * @param[in] beamwidth_deg  Half-power beamwidth (degrees)
 * @param[in] fb_ratio_db  Front-to-back ratio (dB)
 *
 * @return Antenna gain in linear scale (0-1)
 *
 * @note Returns values in [0, 1] range (linear scale, not dB)
 * @note 0° angle difference = maximum gain (1.0)
 * @note >90° angle difference = back lobe (attenuated by F/B ratio)
 */
double directional_antenna_gain(
    double azimuth_deg,
    double antenna_pointing_deg,
    double beamwidth_deg,
    double fb_ratio_db
);

/**
 * @brief Wrap angle to [-180, 180] range
 *
 * Helper function for angle difference calculations.
 *
 * @param[in] angle_deg  Input angle in degrees
 * @return Wrapped angle in [-180, 180] range
 */
double angle_wrap_180(double angle_deg);

/**
 * @brief Calculate angular difference with wrap-around
 *
 * Computes minimum angular distance between two angles,
 * accounting for 0°/360° wrap-around.
 *
 * @param[in] angle1_deg  First angle in degrees
 * @param[in] angle2_deg  Second angle in degrees
 * @return Angular difference in [0, 180] range
 */
double angle_difference(double angle1_deg, double angle2_deg);

#endif /* ARRAY_GEOMETRY_H */
