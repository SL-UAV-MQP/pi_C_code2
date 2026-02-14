/**
 * @file common_types.h
 * @brief Common data types and structures for MUSIC-UCA-6 algorithm
 *
 * Defines shared data structures used across all MUSIC algorithm components.
 * Maps to MATLAB class properties and function return values.
 */

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <complex.h>
#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of antenna elements in UCA */
#define MAX_ANTENNAS 6

/** Maximum number of sources to detect */
#define MAX_SOURCES 4

/** Maximum number of azimuth angles in spectrum */
#define MAX_AZIMUTH_BINS 3600

/** Maximum number of snapshots for covariance estimation */
#define MAX_SNAPSHOTS 8192

/* ============================================================================
 * Return Codes
 * ============================================================================ */

/** Operation return codes */
typedef enum {
    MUSIC_SUCCESS = 0,           /**< Operation successful */
    MUSIC_ERROR_NULL_POINTER,    /**< Null pointer argument */
    MUSIC_ERROR_INVALID_SIZE,    /**< Invalid dimension */
    MUSIC_ERROR_INVALID_CONFIG,  /**< Invalid configuration */
    MUSIC_ERROR_LAPACK_FAILED,   /**< LAPACK operation failed */
    MUSIC_ERROR_NO_PEAKS,        /**< No peaks detected in spectrum */
    MUSIC_ERROR_MEMORY,          /**< Memory allocation failed */
} music_status_t;

/* ============================================================================
 * Complex Number Helpers
 * ============================================================================ */

/** Complex float (single precision) */
typedef float complex cfloat_t;

/** Complex double (double precision) */
typedef double complex cdouble_t;

/* ============================================================================
 * Array Geometry Structures
 * ============================================================================ */

/**
 * @brief 3D position vector [x, y, z] in meters (double precision)
 */
typedef struct {
    double x;  /**< X coordinate (East) */
    double y;  /**< Y coordinate (North) */
    double z;  /**< Z coordinate (Up) */
} position_3d_t;

/** Alias for cross-module compatibility (pathfinding, map_fitting) */
typedef position_3d_t position3d_t;

/**
 * @brief UCA antenna array geometry

 */
typedef struct {
    int num_elements;                          /**< Number of antenna elements (6) */
    double radius_m;                           /**< Array radius in meters */
    double wavelength_m;                       /**< Wavelength at center frequency */
    position_3d_t positions[MAX_ANTENNAS];     /**< Antenna positions in meters */
    double orientations_deg[MAX_ANTENNAS];     /**< Antenna pointing directions (degrees) */
} array_geometry_t;

/* ============================================================================
 * Signal Processing Structures
 * ============================================================================ */

/**
 * @brief Complex signal matrix (M antennas × N snapshots)
 */
typedef struct {
    int M;                           /**< Number of antennas */
    int N;                           /**< Number of snapshots */
    cdouble_t* data;                 /**< Flattened column-major data [M*N] */
} signal_matrix_t;

/**
 * @brief Covariance matrix (M × M complex Hermitian)
 */
typedef struct {
    int M;                           /**< Matrix dimension */
    cdouble_t* data;                 /**< Flattened column-major data [M*M] */
} covariance_matrix_t;

/**
 * @brief Eigendecomposition result
 */
typedef struct {
    int M;                           /**< Matrix dimension */
    double* eigenvalues;             /**< Real eigenvalues (descending order) [M] */
    cdouble_t* eigenvectors;         /**< Complex eigenvectors column-major [M*M] */
} eigendecomp_result_t;

/* ============================================================================
 * MUSIC Spectrum Structures
 * ============================================================================ */

/**
 * @brief MUSIC spectrum result
 */
typedef struct {
    int num_azimuth_bins;            /**< Number of azimuth angle bins */
    double* azimuth_deg;             /**< Azimuth angles in degrees */
    double* spectrum_db;             /**< MUSIC pseudo-spectrum in dB */
    double* eigenvalues;             /**< Copy of eigenvalues for reference */
    int num_eigenvalues;             /**< Number of eigenvalues */
    int num_sources;                 /**< Number of sources used */
} music_spectrum_t;

/* ============================================================================
 * Peak Detection Structures
 * ============================================================================ */

/**
 * @brief Detected source information
 */
typedef struct {
    double azimuth_deg;              /**< Azimuth angle (degrees) */
    double elevation_deg;            /**< Elevation angle (degrees, always 0 for UCA) */
    double magnitude_db;             /**< Peak magnitude in dB */
    double confidence;               /**< Confidence score [0-1] */
    char type[32];                   /**< Peak type ("MUSIC_peak") */
} detected_source_t;

/**
 * @brief Array of detected sources
 */
typedef struct {
    int num_sources;                           /**< Number of detected sources */
    detected_source_t sources[MAX_SOURCES];    /**< Detected source array */
} detected_sources_t;

/* ============================================================================
 * CFAR Detection Structures
 * ============================================================================ */

/**
 * @brief 2D CFAR detection result
 */
typedef struct {
    int height;                      /**< Number of rows (azimuth bins) */
    int width;                       /**< Number of columns (elevation bins) */
    bool* detections;                /**< Binary detection mask [height*width] */
    double* threshold_db;            /**< Adaptive threshold in dB [height*width] */
    double* snr_db;                  /**< SNR at each position in dB [height*width] */
    double noise_floor_db;           /**< Estimated noise floor in dB */
    int n_detections;                /**< Total number of detections */
    double peak_snr_db;              /**< Maximum SNR in dB */
} cfar_result_t;

/* ============================================================================
 * Memory Management Helpers
 * ============================================================================ */

/**
 * @brief Allocate signal matrix
 * @param M Number of antennas
 * @param N Number of snapshots
 * @return Pointer to allocated matrix or NULL on failure
 */
signal_matrix_t* signal_matrix_alloc(int M, int N);

/**
 * @brief Free signal matrix
 * @param mat Matrix to free
 */
void signal_matrix_free(signal_matrix_t* mat);

/**
 * @brief Allocate covariance matrix
 * @param M Matrix dimension
 * @return Pointer to allocated matrix or NULL on failure
 */
covariance_matrix_t* covariance_matrix_alloc(int M);

/**
 * @brief Free covariance matrix
 * @param mat Matrix to free
 */
void covariance_matrix_free(covariance_matrix_t* mat);

/**
 * @brief Allocate eigendecomposition result
 * @param M Matrix dimension
 * @return Pointer to allocated result or NULL on failure
 */
eigendecomp_result_t* eigendecomp_result_alloc(int M);

/**
 * @brief Free eigendecomposition result
 * @param result Result to free
 */
void eigendecomp_result_free(eigendecomp_result_t* result);

/**
 * @brief Allocate MUSIC spectrum result
 * @param num_azimuth_bins Number of azimuth bins
 * @param num_eigenvalues Number of eigenvalues
 * @return Pointer to allocated spectrum or NULL on failure
 */
music_spectrum_t* music_spectrum_alloc(int num_azimuth_bins, int num_eigenvalues);

/**
 * @brief Free MUSIC spectrum result
 * @param spectrum Spectrum to free
 */
void music_spectrum_free(music_spectrum_t* spectrum);

/**
 * @brief Allocate CFAR result
 * @param height Number of rows
 * @param width Number of columns
 * @return Pointer to allocated result or NULL on failure
 */
cfar_result_t* cfar_result_alloc(int height, int width);

/**
 * @brief Free CFAR result
 * @param result Result to free
 */
void cfar_result_free(cfar_result_t* result);

#endif /* COMMON_TYPES_H */
