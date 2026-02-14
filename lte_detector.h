/**
 * @file lte_detector.h
 * @brief LTE/5G signal detection and Cell ID extraction
 *
 * Implements PSS/SSS detection and Physical Cell ID extraction for LTE signals.
 * Based on 3GPP TS 36.211 specifications.
 *
 */

#ifndef LTE_DETECTOR_H
#define LTE_DETECTOR_H

#include "common_types.h"
#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of PSS sequences (N_ID_2 = 0, 1, 2) */
#define LTE_NUM_PSS_SEQUENCES 3

/** Maximum number of SSS sequences (N_ID_1 = 0 to 167) */
#define LTE_NUM_SSS_SEQUENCES 168

/** PSS/SSS sequence length (subcarriers) */
#define LTE_SYNC_SEQ_LENGTH 62

/** Maximum Physical Cell ID (0-503) */
#define LTE_MAX_CELL_ID 503

/** Maximum number of detections per frame */
#define LTE_MAX_DETECTIONS 10

/** Maximum IQ samples per detection call */
#define LTE_MAX_IQ_SAMPLES (30720 * 100)  // 100ms at 30.72 MSPS

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/* cdouble_t provided by common_types.h */

/**
 * @brief PSS detection result
 */
typedef struct {
    int n_id_2;                    /**< PSS identifier (0, 1, or 2) */
    int timing_offset;             /**< Sample offset of PSS peak */
    double correlation;            /**< Normalized correlation value */
    char type[8];                  /**< Detection type "PSS" */
} pss_detection_t;

/**
 * @brief SSS detection result
 */
typedef struct {
    int n_id_1;                    /**< SSS identifier (0-167) */
    int n_id_2;                    /**< PSS identifier (0-2) */
    double correlation;            /**< Normalized correlation value */
    char type[8];                  /**< Detection type "SSS" */
} sss_detection_t;

/**
 * @brief LTE Cell detection result
 */
typedef struct {
    int cell_id;                   /**< Physical Cell ID (0-503) */
    int n_id_1;                    /**< SSS identifier (0-167) */
    int n_id_2;                    /**< PSS identifier (0-2) */
    int timing_offset;             /**< Sample offset */
    double pss_correlation;        /**< PSS correlation strength */
    double sss_correlation;        /**< SSS correlation strength */
    char signal_type[8];           /**< "LTE" or "5GNR" */
} lte_cell_detection_t;

/**
 * @brief LTE detector state
 */
typedef struct {
    double sample_rate;                                  /**< Sample rate in Hz */
    int fft_size;                                        /**< FFT size */
    double subcarrier_spacing;                           /**< 15 kHz for LTE */

    // Pre-generated PSS sequences (Zadoff-Chu)
    cdouble_t pss_sequences[LTE_NUM_PSS_SEQUENCES][LTE_SYNC_SEQ_LENGTH];

    // Pre-generated SSS sequences (M-sequences)
    cdouble_t sss_sequences[LTE_NUM_SSS_SEQUENCES][LTE_SYNC_SEQ_LENGTH];

    // Working buffers for correlation
    cdouble_t* correlation_buffer;
    double* correlation_magnitude;
    int correlation_buffer_size;
} lte_detector_t;

/* ============================================================================
 * Function Declarations
 * ============================================================================ */

/**
 * @brief Initialize LTE detector
 *
 * @param detector Pointer to detector structure
 * @param sample_rate Sample rate in Hz (default: 30.72e6)
 * @return 0 on success, -1 on failure
 */
int lte_detector_init(lte_detector_t* detector, double sample_rate);

/**
 * @brief Free LTE detector resources
 *
 * @param detector Pointer to detector structure
 */
void lte_detector_free(lte_detector_t* detector);

/**
 * @brief Generate PSS sequences (Zadoff-Chu)
 *
 * Generates 3 PSS sequences for N_ID_2 = 0, 1, 2
 * Based on Zadoff-Chu sequences with roots u = 25, 29, 34
 *
 * @param detector Pointer to detector structure
 * @return 0 on success, -1 on failure
 */
int lte_generate_pss_sequences(lte_detector_t* detector);

/**
 * @brief Generate SSS sequences (M-sequences)
 *
 * Generates 168 SSS sequences for N_ID_1 = 0 to 167
 * Based on two M-sequences with different taps
 *
 * @param detector Pointer to detector structure
 * @return 0 on success, -1 on failure
 */
int lte_generate_sss_sequences(lte_detector_t* detector);

/**
 * @brief Detect PSS in IQ samples
 *
 * Performs cross-correlation with all 3 PSS sequences and finds peaks
 * above threshold. Target: <20ms processing time.
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param threshold Detection threshold (0.0 to 1.0, default 0.6)
 * @param detections Output array for detections
 * @param max_detections Maximum detections to return
 * @return Number of detections found (0 to max_detections)
 */
int lte_detect_pss(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    double threshold,
    pss_detection_t* detections,
    int max_detections
);

/**
 * @brief Detect SSS given PSS timing
 *
 * SSS is located in the symbol before PSS. Correlates with all 168
 * SSS sequences to find best match.
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param pss_timing PSS timing offset from PSS detection
 * @param n_id_2 N_ID_2 from PSS (0, 1, or 2)
 * @param threshold Detection threshold (default 0.5)
 * @param sss_det Output SSS detection result
 * @return 0 if SSS detected, -1 if not found
 */
int lte_detect_sss(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int pss_timing,
    int n_id_2,
    double threshold,
    sss_detection_t* sss_det
);

/**
 * @brief Calculate Physical Cell ID from N_ID_1 and N_ID_2
 *
 * Cell ID = 3 * N_ID_1 + N_ID_2
 *
 * @param n_id_1 SSS identifier (0-167)
 * @param n_id_2 PSS identifier (0-2)
 * @return Physical Cell ID (0-503)
 */
int lte_extract_cell_id(int n_id_1, int n_id_2);

/**
 * @brief Full LTE signal detection pipeline
 *
 * Complete detection: PSS -> SSS -> Cell ID extraction
 * Target: <30ms total processing time.
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param cells Output array for detected cells
 * @param max_cells Maximum cells to detect
 * @return Number of cells detected (0 to max_cells)
 */
int lte_detect_signal(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    lte_cell_detection_t* cells,
    int max_cells
);

/**
 * @brief Estimate carrier frequency offset using PSS
 *
 * Uses phase difference between consecutive PSS symbols (5ms apart)
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param pss_timing PSS timing offset
 * @param n_id_2 PSS identifier
 * @return Frequency offset in Hz
 */
double lte_estimate_frequency_offset(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int pss_timing,
    int n_id_2
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Compute cross-correlation between two complex sequences
 *
 * @param signal Input signal
 * @param signal_len Signal length
 * @param reference Reference sequence
 * @param ref_len Reference length
 * @param output Output correlation (size: signal_len - ref_len + 1)
 */
void lte_cross_correlate(
    const cdouble_t* signal,
    int signal_len,
    const cdouble_t* reference,
    int ref_len,
    cdouble_t* output
);

/**
 * @brief Find peaks in a magnitude array
 *
 * @param data Input magnitude data
 * @param len Data length
 * @param threshold Minimum peak height
 * @param min_distance Minimum samples between peaks
 * @param peaks Output array of peak indices
 * @param max_peaks Maximum peaks to find
 * @return Number of peaks found
 */
int lte_find_peaks(
    const double* data,
    int len,
    double threshold,
    int min_distance,
    int* peaks,
    int max_peaks
);

#endif /* LTE_DETECTOR_H */
