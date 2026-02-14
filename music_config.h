/**
 * @file music_config.h
 * @brief Configuration constants for MUSIC-UCA-6 algorithm
 *
 * Contains all validated parameters from MATLAB simulation (MQP_CMRCM_Simulation.m)
 * These values achieved 2.84° mean AOA error in field testing.
 */

#ifndef MUSIC_CONFIG_H
#define MUSIC_CONFIG_H

#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * Physical Array Configuration
 * ============================================================================ */

/** Number of antenna elements in UCA */
#define NUM_ELEMENTS 6

/** Array radius in meters (0.6λ at 1050 MHz) */
#define ARRAY_RADIUS_M 0.176

/** Center frequency in Hz */
#define CENTER_FREQ_HZ 1050000000.0  /* 1050 MHz */

/** Speed of light in m/s */
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792458.0
#endif

/** Wavelength in meters (calculated: c / f) */
#define WAVELENGTH_M (SPEED_OF_LIGHT / CENTER_FREQ_HZ)  /* ~0.2855 m */

/** Wave number k = 2π/λ */
#define WAVE_NUMBER (2.0 * M_PI / WAVELENGTH_M)  /* ~22.0 rad/m */

/* ============================================================================
 * Directional Antenna Pattern Parameters
 * Optimized for 6-element UCA with parabolic reflectors
 * ============================================================================ */

/** Antenna Half-Power Beamwidth (HPBW) in degrees */
#define ANTENNA_BEAMWIDTH_DEG 120.0

/** Front-to-Back ratio in dB */
#define ANTENNA_FB_RATIO_DB 15.0

/** Gaussian beam pattern sigma (derived from HPBW) */
#define ANTENNA_SIGMA_DEG (ANTENNA_BEAMWIDTH_DEG / (2.0 * sqrt(2.0 * log(2.0))))  /* ~50.9° */

/** Front-to-back ratio in linear scale */
#define ANTENNA_FB_RATIO_LINEAR (pow(10.0, -ANTENNA_FB_RATIO_DB / 20.0))  /* ~0.178 */

/* ============================================================================
 * MUSIC Algorithm Parameters
 * Validated values from field testing
 * ============================================================================ */

/** Number of signal sources (increased to 4 for ghost/multipath handling) */
#define NUM_SOURCES 4

/** Auto-estimate number of sources using MDL/AIC */
#define AUTO_ESTIMATE_SOURCES false

/** Use forward-backward averaging for coherent sources */
#define USE_FORWARD_BACKWARD true

/* ============================================================================
 * Spectrum Computation Parameters
 * ============================================================================ */

/** Azimuth scan range minimum (degrees) */
#define AZIMUTH_MIN_DEG -180.0

/** Azimuth scan range maximum (degrees) */
#define AZIMUTH_MAX_DEG 180.0

/** Azimuth angular resolution (degrees) */
#define AZIMUTH_RESOLUTION_DEG 0.2  /* 0.2° for high accuracy, 1800 bins */

/** Number of azimuth bins (calculated) */
#define NUM_AZIMUTH_BINS ((int)((AZIMUTH_MAX_DEG - AZIMUTH_MIN_DEG) / AZIMUTH_RESOLUTION_DEG + 1))

/* ============================================================================
 * Peak Detection Parameters
 * Optimized to prevent ghost peaks (e.g., 180° vs 189°)
 * ============================================================================ */

/** Minimum peak height in dB (lowered for weak signals) */
#define MIN_PEAK_HEIGHT_DB -30.0

/** Minimum peak separation in degrees (increased to prevent ghost peaks) */
#define MIN_PEAK_SEPARATION_DEG 30.0

/** Number of peaks to detect (default: NUM_SOURCES) */
#define NUM_PEAKS_TO_DETECT NUM_SOURCES

/* ============================================================================
 * Signal Processing Parameters
 * ============================================================================ */

/** Number of snapshots for covariance estimation */
#define NUM_SNAPSHOTS 4000

/** SDR sample rate (Hz) - 30.72 MSPS for LTE standard, SDR actual: 20 MSPS */
#define SAMPLE_RATE_HZ 20000000.0  /* 20 MSPS (ADALM-PLUTO max USB 2.0) */

/* ============================================================================
 * Signal Strength Boost Factors
 * Optimized for 4-source subspace (balanced tower strengths)
 * ============================================================================ */

/** Tower 1 boost factor (closest tower, 650m) */
#define TOWER1_BOOST_FACTOR 1.25

/** Tower 2 boost factor (Juniper Hill, 1.7km) */
#define TOWER2_BOOST_FACTOR 2.70

/** Tower 3 boost factor (Chauncy Hill, 1.9km) */
#define TOWER3_BOOST_FACTOR 3.20

/* ============================================================================
 * CFAR Detection Parameters
 * Optimized for CMRCM multipath environment
 * ============================================================================ */

/** CFAR window size [azimuth, elevation] */
#define CFAR_WINDOW_SIZE_AZ 11
#define CFAR_WINDOW_SIZE_EL 7

/** CFAR guard cell size [azimuth, elevation] */
#define CFAR_GUARD_SIZE_AZ 3
#define CFAR_GUARD_SIZE_EL 2

/** CFAR order statistic k-factor (0-1) */
#define CFAR_K_FACTOR 0.7

/** CFAR probability of false alarm */
#define CFAR_PFA 1e-4

/** CFAR threshold offset in dB */
#define CFAR_THRESHOLD_OFFSET_DB 2.5

/** CFAR minimum SNR for detection (dB) */
#define CFAR_MIN_SNR_DB 10.0

/* ============================================================================
 * Real-Time Performance Constraints
 * Target: Raspberry Pi 5B (ARM Cortex-A76, 4 cores @ 2.4 GHz)
 * ============================================================================ */

/** Maximum processing time per frame (milliseconds) */
#define MAX_PROCESSING_TIME_MS 100

/** Target update rate (Hz) */
#define TARGET_UPDATE_RATE_HZ 10

/** Enable ARM NEON SIMD optimizations */
#define USE_ARM_NEON 1

/** Number of threads for parallel processing */
#define NUM_THREADS 4

/* ============================================================================
 * Numerical Stability Constants
 * ============================================================================ */

/** Epsilon for division by zero prevention */
#ifndef EPSILON
#define EPSILON 1e-12
#endif

/** Minimum eigenvalue threshold (for noise subspace selection) */
#define MIN_EIGENVALUE 1e-10

/** Maximum condition number for covariance matrix */
#define MAX_CONDITION_NUMBER 1e6

/* ============================================================================
 * Debug and Logging Configuration
 * ============================================================================ */

/** Enable verbose logging */
#define MUSIC_VERBOSE 0

/** Enable performance profiling (set 1 with DEBUG=1 build) */
#define MUSIC_PROFILE 0

/** Enable MATLAB comparison mode (output intermediate results) */
#define MATLAB_COMPARE_MODE 0

#endif /* MUSIC_CONFIG_H */
