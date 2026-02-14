/**
 * @file signal_processor.h
 * @brief Signal preprocessing for SDR IQ data
 *
 * Implements FFT, PSD, spectrogram, band normalization, filtering,
 * peak detection, downsampling, and SNR estimation for real-time
 * signal processing on Raspberry Pi.
 */

#ifndef SIGNAL_PROCESSOR_H
#define SIGNAL_PROCESSOR_H

#include "common_types.h"
#include <complex.h>
#include <stdbool.h>
#include <fftw3.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum FFT size supported */
#define MAX_FFT_SIZE 8192

/** Maximum number of peaks to detect */
#define MAX_PEAKS 100

/** Maximum filter order */
#define MAX_FILTER_ORDER 20

/** Epsilon for numerical stability - guarded to avoid conflict with music_config.h */
#ifndef EPSILON
#define EPSILON 1e-12
#endif

/* ============================================================================
 * Window Function Types
 * ============================================================================ */

/** Window function types for FFT/PSD */
typedef enum {
    WINDOW_RECTANGULAR = 0,  /**< Rectangular (no window) */
    WINDOW_HAMMING,          /**< Hamming window */
    WINDOW_HANN,             /**< Hann window */
    WINDOW_BLACKMAN,         /**< Blackman window */
    WINDOW_KAISER            /**< Kaiser window (beta=8.6) */
} window_type_t;

/* ============================================================================
 * Normalization Methods
 * ============================================================================ */

/** Band normalization methods */
typedef enum {
    NORMALIZE_PSD = 0,   /**< Normalize based on PSD */
    NORMALIZE_MEAN,      /**< Zero mean, unit variance */
    NORMALIZE_MAX        /**< Normalize to max amplitude */
} normalize_method_t;

/* ============================================================================
 * Signal Processor Configuration
 * ============================================================================ */

/**
 * @brief Signal processor configuration and state
 */
typedef struct {
    double sample_rate;          /**< Sample rate in Hz */
    int fft_size;                /**< FFT size (number of points) */
    window_type_t window_type;   /**< Window function type */
    double overlap;              /**< Overlap ratio (0.0 to 1.0) */
    int hop_size;                /**< Hop size for STFT */
    double* window;              /**< Precomputed window function [fft_size] */
    /* Cached FFTW plan and buffers (created once, reused) */
    fftw_complex* fftw_in;       /**< FFTW input buffer [fft_size] */
    fftw_complex* fftw_out;      /**< FFTW output buffer [fft_size] */
    fftw_plan fftw_forward;      /**< Cached forward FFT plan */
} signal_processor_t;

/* ============================================================================
 * FFT/PSD Result Structures
 * ============================================================================ */

/**
 * @brief FFT result (frequency spectrum)
 */
typedef struct {
    int num_bins;           /**< Number of frequency bins */
    double* frequencies;    /**< Frequency array in Hz [num_bins] */
    double* magnitude_db;   /**< Magnitude spectrum in dB [num_bins] */
} fft_result_t;

/**
 * @brief PSD result (Power Spectral Density)
 */
typedef struct {
    int num_bins;           /**< Number of frequency bins */
    double* frequencies;    /**< Frequency array in Hz [num_bins] */
    double* psd_db;         /**< PSD in dB/Hz [num_bins] */
} psd_result_t;

/**
 * @brief Spectrogram result
 */
typedef struct {
    int num_time_bins;      /**< Number of time bins */
    int num_freq_bins;      /**< Number of frequency bins */
    double* times;          /**< Time array in seconds [num_time_bins] */
    double* frequencies;    /**< Frequency array in Hz [num_freq_bins] */
    double* spectrogram_db; /**< Spectrogram in dB [num_freq_bins * num_time_bins] (column-major) */
} spectrogram_result_t;

/**
 * @brief Peak detection result
 */
typedef struct {
    int num_peaks;              /**< Number of detected peaks */
    int indices[MAX_PEAKS];     /**< Peak indices in spectrum */
    double frequencies[MAX_PEAKS];  /**< Peak frequencies in Hz */
    double magnitudes[MAX_PEAKS];   /**< Peak magnitudes in dB */
    double prominences[MAX_PEAKS];  /**< Peak prominences */
    double widths[MAX_PEAKS];       /**< Peak widths in bins */
} peak_result_t;

/**
 * @brief IQ data buffer (complex samples)
 */
typedef struct {
    int length;              /**< Number of samples */
    cdouble_t* samples;      /**< Complex IQ samples [length] */
} iq_buffer_t;

/**
 * @brief Butterworth filter coefficients
 */
typedef struct {
    int order;               /**< Filter order */
    double b[MAX_FILTER_ORDER + 1];  /**< Numerator coefficients */
    double a[MAX_FILTER_ORDER + 1];  /**< Denominator coefficients */
} filter_coeffs_t;

/* ============================================================================
 * Initialization and Cleanup
 * ============================================================================ */

/**
 * @brief Initialize signal processor
 *
 * @param sample_rate Sample rate in Hz
 * @param fft_size FFT size (must be power of 2, <= MAX_FFT_SIZE)
 * @param window_type Window function type
 * @param overlap Overlap ratio for STFT (0.0 to 1.0)
 * @return Pointer to allocated processor or NULL on failure
 */
signal_processor_t* signal_processor_init(
    double sample_rate,
    int fft_size,
    window_type_t window_type,
    double overlap
);

/**
 * @brief Free signal processor resources
 *
 * @param processor Processor to free
 */
void signal_processor_free(signal_processor_t* processor);

/* ============================================================================
 * FFT and Spectral Analysis
 * ============================================================================ */

/**
 * @brief Compute FFT of IQ samples
 *
 * Target performance: <5ms for 2048-point FFT
 *
 * @param processor Signal processor instance
 * @param iq_samples Complex IQ samples (length >= fft_size)
 * @param num_samples Number of samples in buffer
 * @param result Pointer to allocated FFT result
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_compute_fft(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    fft_result_t* result
);

/**
 * @brief Compute Power Spectral Density using Welch's method
 *
 * Target performance: <10ms
 *
 * @param processor Signal processor instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param result Pointer to allocated PSD result
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_compute_psd(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    psd_result_t* result
);

/**
 * @brief Compute spectrogram using STFT
 *
 * @param processor Signal processor instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param result Pointer to allocated spectrogram result
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_compute_spectrogram(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    spectrogram_result_t* result
);

/* ============================================================================
 * Normalization and Filtering
 * ============================================================================ */

/**
 * @brief Normalize signal power across frequency band
 *
 * @param processor Signal processor instance
 * @param iq_samples Complex IQ samples (modified in-place)
 * @param num_samples Number of samples
 * @param method Normalization method
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_normalize_band(
    const signal_processor_t* processor,
    cdouble_t* iq_samples,
    int num_samples,
    normalize_method_t method
);

/**
 * @brief Design Butterworth bandpass filter
 *
 * @param sample_rate Sample rate in Hz
 * @param low_freq Low cutoff frequency (Hz)
 * @param high_freq High cutoff frequency (Hz)
 * @param order Filter order
 * @param coeffs Output filter coefficients
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_design_bandpass(
    double sample_rate,
    double low_freq,
    double high_freq,
    int order,
    filter_coeffs_t* coeffs
);

/**
 * @brief Apply bandpass filter to IQ samples
 *
 * Target performance: <5ms
 *
 * @param coeffs Filter coefficients
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param output Output filtered samples [num_samples]
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_apply_filter(
    const filter_coeffs_t* coeffs,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
);

/* ============================================================================
 * Peak Detection and Analysis
 * ============================================================================ */

/**
 * @brief Detect peaks in spectrum
 *
 * @param frequencies Frequency array in Hz
 * @param spectrum Magnitude spectrum in dB
 * @param num_bins Number of frequency bins
 * @param height Minimum peak height (or NAN for auto)
 * @param threshold Threshold above noise floor (or NAN for auto)
 * @param distance Minimum distance between peaks (bins)
 * @param result Pointer to peak detection result
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_detect_peaks(
    const double* frequencies,
    const double* spectrum,
    int num_bins,
    double height,
    double threshold,
    int distance,
    peak_result_t* result
);

/**
 * @brief Downsample IQ data with anti-aliasing filter
 *
 * @param iq_samples Input complex IQ samples
 * @param num_samples Number of input samples
 * @param factor Downsampling factor
 * @param output Output downsampled samples [num_samples/factor]
 * @param output_length Pointer to store output length
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_downsample(
    const cdouble_t* iq_samples,
    int num_samples,
    int factor,
    cdouble_t* output,
    int* output_length
);

/**
 * @brief Estimate Signal-to-Noise Ratio
 *
 * @param processor Signal processor instance
 * @param iq_samples Complex IQ samples
 * @param num_samples Number of samples
 * @param signal_band_low Low frequency of signal band (Hz)
 * @param signal_band_high High frequency of signal band (Hz)
 * @param snr_db Output SNR in dB
 * @return MUSIC_SUCCESS or error code
 */
music_status_t signal_processor_estimate_snr(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    double signal_band_low,
    double signal_band_high,
    double* snr_db
);

/* ============================================================================
 * Memory Management Helpers
 * ============================================================================ */

/**
 * @brief Allocate FFT result structure
 *
 * @param num_bins Number of frequency bins
 * @return Pointer to allocated result or NULL on failure
 */
fft_result_t* fft_result_alloc(int num_bins);

/**
 * @brief Free FFT result
 *
 * @param result Result to free
 */
void fft_result_free(fft_result_t* result);

/**
 * @brief Allocate PSD result structure
 *
 * @param num_bins Number of frequency bins
 * @return Pointer to allocated result or NULL on failure
 */
psd_result_t* psd_result_alloc(int num_bins);

/**
 * @brief Free PSD result
 *
 * @param result Result to free
 */
void psd_result_free(psd_result_t* result);

/**
 * @brief Allocate spectrogram result structure
 *
 * @param num_time_bins Number of time bins
 * @param num_freq_bins Number of frequency bins
 * @return Pointer to allocated result or NULL on failure
 */
spectrogram_result_t* spectrogram_result_alloc(int num_time_bins, int num_freq_bins);

/**
 * @brief Free spectrogram result
 *
 * @param result Result to free
 */
void spectrogram_result_free(spectrogram_result_t* result);

/**
 * @brief Allocate IQ buffer
 *
 * @param length Number of samples
 * @return Pointer to allocated buffer or NULL on failure
 */
iq_buffer_t* iq_buffer_alloc(int length);

/**
 * @brief Free IQ buffer
 *
 * @param buffer Buffer to free
 */
void iq_buffer_free(iq_buffer_t* buffer);

#endif /* SIGNAL_PROCESSOR_H */
