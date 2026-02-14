/**
 * @file noise_reduction.h
 * @brief Advanced Noise Reduction Algorithms
 *
 * Implements spectral and spatial noise reduction for urban RF environment.
 * Designed for weak signals (-60 to -90 dBm) in multipath/interference conditions.
 *
 * Features:
 * - Spectral subtraction with over-subtraction and spectral floor
 * - Wiener filtering with adaptive noise estimation
 * - Spatial noise cancellation (null steering)
 * - Adaptive noise floor estimation with percentile-based method
 * - SNR estimation and improvement measurement
 *
 * Performance Targets:
 * - Single-channel processing: <50ms (real-time capable)
 * - SNR improvement: 5-15 dB typical
 * - Multi-channel spatial filtering: <100ms
 */

#ifndef NOISE_REDUCTION_H
#define NOISE_REDUCTION_H

#include "common_types.h"
#include <complex.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Default FFT size for STFT processing */
#define NR_FFT_SIZE 2048

/** Default overlap for STFT (0.75 = 75%) */
#define NR_OVERLAP 0.75

/** Maximum number of noise estimation frames */
#define NR_MAX_HISTORY_FRAMES 20

/** Default spectral floor parameter */
#define NR_SPECTRAL_FLOOR 0.01

/** Default over-subtraction factor */
#define NR_ALPHA 2.0

/** Minimum Wiener gain (to prevent over-suppression) */
#define NR_WIENER_GAIN_MIN 0.1

/* ============================================================================
 * Noise Reduction Methods
 * ============================================================================ */

/**
 * @brief Noise reduction algorithm selection
 */
typedef enum {
    NR_SPECTRAL_SUBTRACTION = 0,  /**< Spectral subtraction method */
    NR_WIENER_FILTER,             /**< Wiener filtering */
    NR_SPATIAL_CANCELLATION,      /**< Spatial null steering */
    NR_ADAPTIVE_ESTIMATION,       /**< Adaptive noise floor estimation */
    NR_COMBINED                   /**< Combined spatial + temporal */
} noise_reduction_method_t;

/* ============================================================================
 * Configuration Structures
 * ============================================================================ */

/**
 * @brief Noise reduction configuration
 */
typedef struct {
    /* Spectral subtraction parameters */
    double alpha;                      /**< Over-subtraction factor (1.0-3.0) */
    double beta;                       /**< Spectral floor parameter */

    /* Wiener filter parameters */
    int noise_estimation_frames;       /**< Number of frames for noise estimation */
    double wiener_gain_min;            /**< Minimum Wiener gain */

    /* Spatial filtering parameters */
    int num_elements;                  /**< Number of array elements */
    double null_depth;                 /**< Null depth (dB) */

    /* Adaptive estimation */
    double adaptation_rate;            /**< Noise floor adaptation rate (0-1) */
    double noise_percentile;           /**< Percentile for noise estimation */

    /* General parameters */
    int fft_size;                      /**< FFT size */
    double overlap;                    /**< STFT overlap ratio (0.0-1.0) */
    double sample_rate;                /**< Sample rate (Hz) */
} noise_reduction_config_t;

/* ============================================================================
 * Performance Metrics
 * ============================================================================ */

/**
 * @brief Noise reduction performance metrics
 */
typedef struct {
    char method[32];                   /**< Method name */
    double input_snr_db;               /**< Input SNR (dB) */
    double output_snr_db;              /**< Output SNR (dB) */
    double snr_improvement_db;         /**< SNR improvement (dB) */
    double input_power_dbm;            /**< Input power (dBm) */
    double output_power_dbm;           /**< Output power (dBm) */
    double processing_time_ms;         /**< Processing time (ms) */
    int samples_processed;             /**< Number of samples processed */
} noise_reduction_metrics_t;

/* ============================================================================
 * Noise Floor Estimator
 * ============================================================================ */

/**
 * @brief Noise floor estimator state
 */
typedef struct {
    noise_reduction_config_t config;   /**< Configuration */
    double* noise_floor;               /**< Current noise floor estimate [fft_size/2] */
    int num_bins;                      /**< Number of frequency bins */
    bool initialized;                  /**< Whether estimator is initialized */
} noise_floor_estimator_t;

/**
 * @brief Initialize noise floor estimator
 *
 * @param estimator Pointer to estimator structure (allocated by caller)
 * @param config Noise reduction configuration
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_floor_estimator_init(
    noise_floor_estimator_t* estimator,
    const noise_reduction_config_t* config
);

/**
 * @brief Free noise floor estimator resources
 *
 * @param estimator Estimator to free
 */
void noise_floor_estimator_free(noise_floor_estimator_t* estimator);

/**
 * @brief Estimate noise floor from power spectrum
 *
 * Uses percentile-based estimation with median filtering.
 *
 * @param estimator Estimator instance
 * @param spectrum Power spectrum (linear scale) [num_bins]
 * @param num_bins Number of frequency bins
 * @param noise_floor Output noise floor estimate [num_bins]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_floor_estimate(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* noise_floor
);

/**
 * @brief Update noise floor adaptively
 *
 * Uses exponential moving average for adaptation.
 *
 * @param estimator Estimator instance
 * @param spectrum Current power spectrum [num_bins]
 * @param num_bins Number of bins
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_floor_update(
    noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins
);

/**
 * @brief Estimate overall SNR
 *
 * @param estimator Estimator instance
 * @param spectrum Power spectrum [num_bins]
 * @param num_bins Number of bins
 * @param snr_db Output SNR in dB
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_floor_get_snr(
    const noise_floor_estimator_t* estimator,
    const double* spectrum,
    int num_bins,
    double* snr_db
);

/* ============================================================================
 * Spectral Subtraction
 * ============================================================================ */

/**
 * @brief Spectral subtractor state
 */
typedef struct {
    noise_reduction_config_t config;
    noise_floor_estimator_t noise_estimator;
    double* window;                    /**< Window function [fft_size] */
} spectral_subtractor_t;

/**
 * @brief Initialize spectral subtractor
 *
 * @param subtractor Pointer to subtractor structure (allocated by caller)
 * @param config Noise reduction configuration
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spectral_subtractor_init(
    spectral_subtractor_t* subtractor,
    const noise_reduction_config_t* config
);

/**
 * @brief Free spectral subtractor resources
 *
 * @param subtractor Subtractor to free
 */
void spectral_subtractor_free(spectral_subtractor_t* subtractor);

/**
 * @brief Apply spectral subtraction to IQ samples
 *
 * Algorithm:
 * 1. Compute STFT
 * 2. Estimate noise floor for each frame
 * 3. Apply spectral subtraction: S_clean = sqrt(max(|S|^2 - α*N, β*|S|^2))
 * 4. Inverse STFT
 *
 * Performance target: <30ms for 10ms of data
 *
 * @param subtractor Subtractor instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param output Cleaned IQ samples [num_samples] (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spectral_subtractor_process(
    spectral_subtractor_t* subtractor,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
);

/* ============================================================================
 * Wiener Filter
 * ============================================================================ */

/**
 * @brief Wiener filter state
 */
typedef struct {
    noise_reduction_config_t config;
    noise_floor_estimator_t noise_estimator;
    double* window;                    /**< Window function [fft_size] */
} wiener_filter_t;

/**
 * @brief Initialize Wiener filter
 *
 * @param filter Pointer to filter structure (allocated by caller)
 * @param config Noise reduction configuration
 * @return MUSIC_SUCCESS or error code
 */
music_status_t wiener_filter_init(
    wiener_filter_t* filter,
    const noise_reduction_config_t* config
);

/**
 * @brief Free Wiener filter resources
 *
 * @param filter Filter to free
 */
void wiener_filter_free(wiener_filter_t* filter);

/**
 * @brief Compute Wiener filter gain
 *
 * Wiener gain: G = SNR / (SNR + 1) = (S - N) / S
 *
 * @param signal_power Signal power spectrum [num_bins]
 * @param noise_power Noise power spectrum [num_bins]
 * @param num_bins Number of frequency bins
 * @param gain_min Minimum gain threshold
 * @param wiener_gain Output Wiener gain [num_bins]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t wiener_compute_gain(
    const double* signal_power,
    const double* noise_power,
    int num_bins,
    double gain_min,
    double* wiener_gain
);

/**
 * @brief Apply Wiener filtering
 *
 * Performance target: <30ms for 10ms of data
 *
 * @param filter Filter instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param output Filtered IQ samples [num_samples] (allocated by caller)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t wiener_filter_process(
    wiener_filter_t* filter,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
);

/* ============================================================================
 * Spatial Noise Canceller
 * ============================================================================ */

/**
 * @brief Spatial noise canceller (array-based null steering)
 */
typedef struct {
    noise_reduction_config_t config;
    double* antenna_positions;         /**< Antenna positions [num_elements x 3] */
    int num_elements;                  /**< Number of array elements */
} spatial_noise_canceller_t;

/**
 * @brief Initialize spatial noise canceller
 *
 * @param canceller Pointer to canceller structure (allocated by caller)
 * @param config Noise reduction configuration
 * @param antenna_positions Antenna positions [num_elements x 3] (NULL for default UCA)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spatial_canceller_init(
    spatial_noise_canceller_t* canceller,
    const noise_reduction_config_t* config,
    const double* antenna_positions
);

/**
 * @brief Free spatial canceller resources
 *
 * @param canceller Canceller to free
 */
void spatial_canceller_free(spatial_noise_canceller_t* canceller);

/**
 * @brief Compute null steering weights using LCMV
 *
 * Linearly Constrained Minimum Variance (LCMV) beamforming:
 * - Unit gain in desired direction (optional)
 * - Nulls in interference directions
 *
 * @param canceller Canceller instance
 * @param interference_azimuths Interference directions (degrees)
 * @param num_interferers Number of interference sources
 * @param frequency Operating frequency (Hz)
 * @param desired_azimuth Desired signal azimuth (degrees, <0 for omni)
 * @param weights Output complex weights [num_elements]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spatial_compute_null_weights(
    const spatial_noise_canceller_t* canceller,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* weights
);

/**
 * @brief Apply spatial noise cancellation to multichannel data
 *
 * @param canceller Canceller instance
 * @param multichannel_samples Input (num_elements x num_samples)
 * @param num_elements Number of channels
 * @param num_samples Number of samples
 * @param interference_azimuths Interference directions (degrees)
 * @param num_interferers Number of interferers
 * @param frequency Operating frequency (Hz)
 * @param desired_azimuth Desired signal direction (degrees, <0 for omni)
 * @param output Output filtered signal [num_samples]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spatial_canceller_process(
    spatial_noise_canceller_t* canceller,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    const double* interference_azimuths,
    int num_interferers,
    double frequency,
    double desired_azimuth,
    cdouble_t* output
);

/**
 * @brief Estimate interference directions from multichannel data
 *
 * Uses spatial spectrum scanning to find interference sources.
 *
 * @param canceller Canceller instance
 * @param multichannel_samples Input (num_elements x num_samples)
 * @param num_elements Number of elements
 * @param num_samples Number of samples
 * @param frequency Operating frequency (Hz)
 * @param num_sources Number of sources to find
 * @param interference_azimuths Output interference azimuths [num_sources]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t spatial_estimate_interference(
    const spatial_noise_canceller_t* canceller,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    double frequency,
    int num_sources,
    double* interference_azimuths
);

/* ============================================================================
 * Unified Noise Reducer
 * ============================================================================ */

/**
 * @brief Unified noise reducer combining all methods
 */
typedef struct {
    noise_reduction_config_t config;
    spectral_subtractor_t spectral_subtractor;
    wiener_filter_t wiener_filter;
    spatial_noise_canceller_t spatial_canceller;
} noise_reducer_t;

/**
 * @brief Initialize noise reducer
 *
 * @param reducer Pointer to reducer structure (allocated by caller)
 * @param config Configuration (NULL for defaults)
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_reducer_init(
    noise_reducer_t* reducer,
    const noise_reduction_config_t* config
);

/**
 * @brief Free noise reducer resources
 *
 * @param reducer Reducer to free
 */
void noise_reducer_free(noise_reducer_t* reducer);

/**
 * @brief Process single-channel IQ data
 *
 * @param reducer Noise reducer instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param method Noise reduction method
 * @param output Cleaned samples [num_samples] (allocated by caller)
 * @param metrics Output performance metrics
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_reducer_process_single(
    noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    noise_reduction_method_t method,
    cdouble_t* output,
    noise_reduction_metrics_t* metrics
);

/**
 * @brief Process multichannel IQ data with spatial + temporal filtering
 *
 * @param reducer Noise reducer instance
 * @param multichannel_samples Input (num_elements x num_samples)
 * @param num_elements Number of channels
 * @param num_samples Number of samples
 * @param frequency Operating frequency (Hz)
 * @param interference_azimuths Interference directions (NULL for auto-detect)
 * @param num_interferers Number of interferers (0 for auto-detect)
 * @param desired_azimuth Desired signal direction (degrees, <0 for omni)
 * @param temporal_method Temporal filtering method
 * @param output Cleaned signal [num_samples] (allocated by caller)
 * @param metrics Output performance metrics
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_reducer_process_multi(
    noise_reducer_t* reducer,
    const cdouble_t* multichannel_samples,
    int num_elements,
    int num_samples,
    double frequency,
    const double* interference_azimuths,
    int num_interferers,
    double desired_azimuth,
    noise_reduction_method_t temporal_method,
    cdouble_t* output,
    noise_reduction_metrics_t* metrics
);

/**
 * @brief Estimate SNR of IQ samples
 *
 * @param reducer Reducer instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param snr_db Output SNR in dB
 * @return MUSIC_SUCCESS or error code
 */
music_status_t noise_reducer_estimate_snr(
    const noise_reducer_t* reducer,
    const cdouble_t* iq_samples,
    int num_samples,
    double* snr_db
);

/* ============================================================================
 * Factory Functions
 * ============================================================================ */

/**
 * @brief Create default noise reduction configuration
 *
 * @param config Output configuration structure
 */
void noise_reduction_config_default(noise_reduction_config_t* config);

/**
 * @brief Create noise reducer with custom parameters
 *
 * @param reducer Reducer structure (allocated by caller)
 * @param sample_rate Sample rate (Hz)
 * @param fft_size FFT size
 * @param alpha Over-subtraction factor
 * @param adaptation_rate Noise floor adaptation rate
 * @return MUSIC_SUCCESS or error code
 */
music_status_t create_noise_reducer(
    noise_reducer_t* reducer,
    double sample_rate,
    int fft_size,
    double alpha,
    double adaptation_rate
);

#endif /* NOISE_REDUCTION_H */
