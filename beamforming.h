/**
 * @file beamforming.h
 * @brief Beamforming Algorithms for UCA 6-Element Array
 *
 * Implements conventional (Delay-and-Sum) and adaptive (MVDR/Capon)
 * beamforming for spatial filtering and interference rejection.
 *
 * Features:
 * - Conventional beamforming with steering
 * - Adaptive MVDR (Minimum Variance Distortionless Response) beamforming
 * - MPDR (Multiple Point Distortionless Response) with nulls
 * - Covariance estimation with forward-backward averaging
 * - Beam pattern computation
 * - Integration with MUSIC AOA estimator
 *
 * Performance Targets:
 * - Weight computation: <1ms
 * - Beamforming application: <5ms for 10k snapshots
 * - Covariance estimation: <10ms
 */

#ifndef BEAMFORMING_H
#define BEAMFORMING_H

#include "common_types.h"
#include <complex.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Default number of antenna elements (UCA) */
#define BF_NUM_ELEMENTS 6

/** Default UCA radius in meters (0.6λ at 1050 MHz) */
#define BF_UCA_RADIUS 0.176

/** Default center frequency (1050 MHz for 5G/LTE) */
#define BF_CENTER_FREQUENCY 1050e6

/** Default diagonal loading factor for MVDR */
#define BF_DIAGONAL_LOADING 0.01

/** Maximum number of interference nulls */
#define BF_MAX_NULLS 5

/** Maximum azimuth scan points for beam pattern */
#define BF_MAX_PATTERN_POINTS 720

/* ============================================================================
 * Beamforming Configuration
 * ============================================================================ */

/**
 * @brief Beamforming configuration parameters
 */
typedef struct {
    int num_elements;             /**< Number of antenna elements */
    double radius;                /**< UCA radius (meters) */
    double center_frequency;      /**< Center frequency (Hz) */
    double sample_rate;           /**< Sample rate (Hz) */

    /* Antenna characteristics */
    double antenna_beamwidth;     /**< Antenna HPBW (degrees) */
    double antenna_fb_ratio;      /**< Front-to-back ratio (dB) */

    /* Beamforming parameters */
    double diagonal_loading;      /**< MVDR regularization factor */
    bool use_forward_backward;    /**< Use forward-backward averaging */
} beamforming_config_t;

/* ============================================================================
 * Beamforming Result Structures
 * ============================================================================ */

/**
 * @brief Beamforming weights (complex)
 */
typedef struct {
    int num_elements;             /**< Number of weights */
    cdouble_t* weights;           /**< Complex beamforming weights [num_elements] */
    double sinr_db;               /**< Estimated SINR (dB) */
} beamforming_weights_t;

/**
 * @brief Beamforming output result
 */
typedef struct {
    int num_samples;              /**< Number of output samples */
    cdouble_t* output;            /**< Beamformed output signal [num_samples] */
    cdouble_t* weights;           /**< Applied weights [num_elements] */
    double sinr_db;               /**< Estimated SINR (dB) */
    char method[32];              /**< Beamforming method used */
} beamforming_result_t;

/**
 * @brief Beam pattern result
 */
typedef struct {
    int num_points;               /**< Number of azimuth points */
    double* azimuths;             /**< Azimuth angles (degrees) [num_points] */
    double* pattern_db;           /**< Beam pattern (dB, normalized) [num_points] */
    double look_direction;        /**< Look direction (degrees) */
    double beamwidth_3db;         /**< 3dB beamwidth (degrees) */
} beam_pattern_t;

/* covariance_matrix_t provided by common_types.h (field: M, not size) */

/* ============================================================================
 * Conventional Beamforming
 * ============================================================================ */

/**
 * @brief Conventional beamformer context
 */
typedef struct {
    beamforming_config_t config;  /**< Configuration */
    void* steering_calculator;    /**< UCA steering vector calculator */
} conventional_beamformer_t;

/**
 * @brief Initialize conventional beamformer
 *
 * @param bf Pointer to beamformer structure (allocated by caller)
 * @param config Beamforming configuration
 * @return MUSIC_SUCCESS or error code
 */
music_status_t conventional_beamformer_init(
    conventional_beamformer_t* bf,
    const beamforming_config_t* config
);

/**
 * @brief Free conventional beamformer resources
 *
 * @param bf Beamformer to free
 */
void conventional_beamformer_free(conventional_beamformer_t* bf);

/**
 * @brief Compute conventional beamforming weights
 *
 * Weights are conjugate of steering vector: w = a*
 *
 * @param bf Beamformer instance
 * @param azimuth Look direction azimuth (degrees)
 * @param elevation Look direction elevation (degrees)
 * @param normalize Normalize weights to unit norm
 * @param weights Output weights (allocated by caller) [num_elements]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t conventional_compute_weights(
    const conventional_beamformer_t* bf,
    double azimuth,
    double elevation,
    bool normalize,
    cdouble_t* weights
);

/**
 * @brief Apply beamforming weights to signals
 *
 * Output: y = w^H * x
 *
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of antenna elements
 * @param num_snapshots Number of time snapshots
 * @param weights Beamforming weights [num_elements]
 * @param output Output beamformed signal [num_snapshots]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t conventional_apply_beamformer(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    const cdouble_t* weights,
    cdouble_t* output
);

/**
 * @brief Steer beam towards specified direction
 *
 * @param bf Beamformer instance
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of elements
 * @param num_snapshots Number of snapshots
 * @param azimuth Desired azimuth (degrees)
 * @param elevation Desired elevation (degrees)
 * @param output Output beamformed signal [num_snapshots]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t conventional_steer_beam(
    const conventional_beamformer_t* bf,
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    double azimuth,
    double elevation,
    cdouble_t* output
);

/**
 * @brief Compute conventional beam pattern
 *
 * @param bf Beamformer instance
 * @param azimuth_look Look direction (degrees)
 * @param azimuth_min Minimum azimuth (degrees)
 * @param azimuth_max Maximum azimuth (degrees)
 * @param resolution Angular resolution (degrees)
 * @param pattern Output beam pattern (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t conventional_compute_pattern(
    const conventional_beamformer_t* bf,
    double azimuth_look,
    double azimuth_min,
    double azimuth_max,
    double resolution,
    beam_pattern_t* pattern
);

/* ============================================================================
 * Adaptive Beamforming (MVDR/Capon)
 * ============================================================================ */

/**
 * @brief Adaptive beamformer context
 */
typedef struct {
    beamforming_config_t config;  /**< Configuration */
    void* steering_calculator;    /**< UCA steering vector calculator */
} adaptive_beamformer_t;

/**
 * @brief Initialize adaptive beamformer
 *
 * @param bf Pointer to beamformer structure (allocated by caller)
 * @param config Beamforming configuration
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_beamformer_init(
    adaptive_beamformer_t* bf,
    const beamforming_config_t* config
);

/**
 * @brief Free adaptive beamformer resources
 *
 * @param bf Beamformer to free
 */
void adaptive_beamformer_free(adaptive_beamformer_t* bf);

/**
 * @brief Estimate spatial covariance matrix
 *
 * R = (1/N) * X * X^H
 *
 * Optionally uses forward-backward averaging for enhanced performance.
 *
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of antenna elements
 * @param num_snapshots Number of snapshots
 * @param use_forward_backward Use forward-backward averaging
 * @param R Output covariance matrix (allocated by caller) [M x M]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_estimate_covariance(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    bool use_forward_backward,
    cdouble_t* R
);

/**
 * @brief Compute MVDR beamforming weights
 *
 * MVDR solution: w = (R + εI)^(-1) * a / (a^H * (R + εI)^(-1) * a)
 *
 * where:
 * - R: covariance matrix
 * - a: steering vector (look direction)
 * - ε: diagonal loading factor
 *
 * Performance target: <5ms
 *
 * @param bf Beamformer instance
 * @param R Covariance matrix [M x M]
 * @param num_elements Number of elements
 * @param azimuth Look direction azimuth (degrees)
 * @param elevation Look direction elevation (degrees)
 * @param diagonal_loading Regularization factor (0 for config default)
 * @param weights Output MVDR weights [num_elements]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_compute_mvdr_weights(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double azimuth,
    double elevation,
    double diagonal_loading,
    cdouble_t* weights
);

/**
 * @brief Compute MPDR weights with explicit nulls
 *
 * Multiple Point Distortionless Response:
 * - Unit gain in desired direction
 * - Nulls in specified interference directions
 *
 * @param bf Beamformer instance
 * @param R Covariance matrix [M x M]
 * @param num_elements Number of elements
 * @param desired_azimuth Desired signal azimuth (degrees)
 * @param interference_azimuths Interference azimuths (degrees)
 * @param num_interferers Number of interferers
 * @param weights Output MPDR weights [num_elements]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_compute_mpdr_weights(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double desired_azimuth,
    const double* interference_azimuths,
    int num_interferers,
    cdouble_t* weights
);

/**
 * @brief Apply adaptive beamforming weights
 *
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of elements
 * @param num_snapshots Number of snapshots
 * @param weights Beamforming weights [num_elements]
 * @param output Output signal [num_snapshots]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_apply_beamformer(
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    const cdouble_t* weights,
    cdouble_t* output
);

/**
 * @brief Adaptively steer beam with MVDR
 *
 * @param bf Beamformer instance
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of elements
 * @param num_snapshots Number of snapshots
 * @param azimuth Desired azimuth (degrees)
 * @param elevation Desired elevation (degrees)
 * @param result Output beamforming result (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_steer_beam(
    const adaptive_beamformer_t* bf,
    const cdouble_t* signals,
    int num_elements,
    int num_snapshots,
    double azimuth,
    double elevation,
    beamforming_result_t* result
);

/**
 * @brief Compute adaptive beam pattern
 *
 * @param bf Beamformer instance
 * @param R Covariance matrix [M x M]
 * @param num_elements Number of elements
 * @param azimuth_look Look direction (degrees)
 * @param azimuth_min Minimum azimuth (degrees)
 * @param azimuth_max Maximum azimuth (degrees)
 * @param resolution Angular resolution (degrees)
 * @param pattern Output beam pattern (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t adaptive_compute_pattern(
    const adaptive_beamformer_t* bf,
    const cdouble_t* R,
    int num_elements,
    double azimuth_look,
    double azimuth_min,
    double azimuth_max,
    double resolution,
    beam_pattern_t* pattern
);

/* ============================================================================
 * Unified Beamforming Processor
 * ============================================================================ */

/**
 * @brief Beamforming method selection
 */
typedef enum {
    BF_METHOD_CONVENTIONAL = 0,   /**< Conventional delay-and-sum */
    BF_METHOD_ADAPTIVE,           /**< Adaptive MVDR */
    BF_METHOD_MPDR                /**< MPDR with nulls */
} beamforming_method_t;

/**
 * @brief Unified beamforming processor
 */
typedef struct {
    beamforming_config_t config;
    conventional_beamformer_t conventional;
    adaptive_beamformer_t adaptive;
} beamforming_processor_t;

/**
 * @brief Initialize beamforming processor
 *
 * @param processor Pointer to processor structure (allocated by caller)
 * @param config Configuration (NULL for defaults)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t beamforming_processor_init(
    beamforming_processor_t* processor,
    const beamforming_config_t* config
);

/**
 * @brief Free beamforming processor resources
 *
 * @param processor Processor to free
 */
void beamforming_processor_free(beamforming_processor_t* processor);

/**
 * @brief Process signals with beamforming
 *
 * @param processor Beamforming processor
 * @param signals Input signals (num_elements x num_snapshots)
 * @param num_elements Number of elements
 * @param num_snapshots Number of snapshots
 * @param method Beamforming method
 * @param desired_azimuth Desired signal azimuth (degrees)
 * @param interference_azimuths Interference azimuths (degrees, for MPDR)
 * @param num_interferers Number of interferers
 * @param result Output beamforming result (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
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
);

/**
 * @brief Estimate SINR from beamformer output
 *
 * @param input_signals Input array signals (num_elements x num_snapshots)
 * @param output_signal Beamformed output [num_snapshots]
 * @param weights Beamforming weights [num_elements]
 * @param num_elements Number of elements
 * @param num_snapshots Number of snapshots
 * @param sinr_db Output SINR in dB
 * @return MUSIC_SUCCESS or error code
 */
music_status_t beamforming_estimate_sinr(
    const cdouble_t* input_signals,
    const cdouble_t* output_signal,
    const cdouble_t* weights,
    int num_elements,
    int num_snapshots,
    double* sinr_db
);

/* ============================================================================
 * Memory Management
 * ============================================================================ */

/**
 * @brief Allocate beam pattern structure
 *
 * @param num_points Number of azimuth points
 * @return Pointer to allocated pattern or NULL on failure
 */
beam_pattern_t* beam_pattern_alloc(int num_points);

/**
 * @brief Free beam pattern
 *
 * @param pattern Pattern to free
 */
void beam_pattern_free(beam_pattern_t* pattern);

/**
 * @brief Allocate beamforming result structure
 *
 * @param num_elements Number of antenna elements
 * @param num_samples Number of output samples
 * @return Pointer to allocated result or NULL on failure
 */
beamforming_result_t* beamforming_result_alloc(int num_elements, int num_samples);

/**
 * @brief Free beamforming result
 *
 * @param result Result to free
 */
void beamforming_result_free(beamforming_result_t* result);

/* covariance_matrix_alloc / covariance_matrix_free provided by common_types.h */

/* ============================================================================
 * Factory Functions
 * ============================================================================ */

/**
 * @brief Create default beamforming configuration
 *
 * @param config Output configuration structure
 */
void beamforming_config_default(beamforming_config_t* config);

/**
 * @brief Create beamforming processor with custom parameters
 *
 * @param processor Processor structure (allocated by caller)
 * @param frequency Center frequency (Hz)
 * @param sample_rate Sample rate (Hz)
 * @param diagonal_loading MVDR regularization factor
 * @return MUSIC_SUCCESS or error code
 */
music_status_t create_beamformer(
    beamforming_processor_t* processor,
    double frequency,
    double sample_rate,
    double diagonal_loading
);

#endif /* BEAMFORMING_H */
