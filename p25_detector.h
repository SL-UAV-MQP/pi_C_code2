/**
 * @file p25_detector.h
 * @brief P25 (Project 25) Land Mobile Radio signal detection
 *
 * Implements detection and decoding of P25 Phase 1 and Phase 2 signals.
 */

#ifndef P25_DETECTOR_H
#define P25_DETECTOR_H

#include "common_types.h"
#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** P25 symbol rate (symbols per second) */
#define P25_SYMBOL_RATE 4800

/** P25 Frame Sync pattern length (bits) */
#define P25_FRAME_SYNC_LENGTH 48

/** Maximum number of P25 detections per call */
#define P25_MAX_DETECTIONS 10

/** Maximum IQ samples per detection call */
#define P25_MAX_IQ_SAMPLES (1000000 * 1)  // 1 second at 1 MSPS

/** Maximum bit buffer size */
#define P25_MAX_BIT_BUFFER (100000)

/** Network Identifier (NID) length in bits */
#define P25_NID_LENGTH 64

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/* cdouble_t provided by common_types.h */

/**
 * @brief P25 Network Identifier (NID) structure
 */
typedef struct {
    uint16_t nac;              /**< Network Access Code (12 bits) */
    uint8_t duid;              /**< Data Unit ID (4 bits) */
    int position;              /**< Bit position in stream */
} p25_nid_t;

/**
 * @brief P25 detection result
 */
typedef struct {
    char signal_type[8];       /**< "P25" */
    uint16_t nac;              /**< Network Access Code */
    uint8_t duid;              /**< Data Unit ID */
    int timing;                /**< Frame sync bit position */
    char modulation[8];        /**< "C4FM" */
} p25_detection_t;

/**
 * @brief P25 detector state
 */
typedef struct {
    double sample_rate;                        /**< Sample rate in Hz */
    int symbol_rate;                           /**< Symbol rate (4800 sps) */
    int samples_per_symbol;                    /**< Samples per symbol */

    // Frame sync pattern (48 bits)
    int frame_sync[P25_FRAME_SYNC_LENGTH];

    // Working buffers
    int* symbol_buffer;
    int symbol_buffer_size;
    int* bit_buffer;
    int bit_buffer_size;
} p25_detector_t;

/* ============================================================================
 * Function Declarations
 * ============================================================================ */

/**
 * @brief Initialize P25 detector
 *
 * @param detector Pointer to detector structure
 * @param sample_rate Sample rate in Hz (default: 1e6)
 * @return 0 on success, -1 on failure
 */
int p25_detector_init(p25_detector_t* detector, double sample_rate);

/**
 * @brief Free P25 detector resources
 *
 * @param detector Pointer to detector structure
 */
void p25_detector_free(p25_detector_t* detector);

/**
 * @brief Demodulate C4FM (4-level FSK) signal
 *
 * Converts IQ samples to symbol stream (dibits: 0, 1, 2, 3)
 * Target: <10ms processing time.
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param symbols Output symbol buffer
 * @param max_symbols Maximum symbols to output
 * @return Number of symbols demodulated
 */
int p25_demodulate_c4fm(
    p25_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int* symbols,
    int max_symbols
);

/**
 * @brief Convert dibits to bit stream
 *
 * @param symbols Input symbol array (dibits 0-3)
 * @param num_symbols Number of symbols
 * @param bits Output bit array
 * @return Number of bits generated (2 * num_symbols)
 */
int p25_symbols_to_bits(
    const int* symbols,
    int num_symbols,
    int* bits
);

/**
 * @brief Detect P25 frame sync pattern
 *
 * Uses Hamming distance to find sync pattern with error tolerance.
 * Target: <15ms processing time.
 *
 * @param detector Pointer to detector structure
 * @param bits Bit stream
 * @param num_bits Number of bits
 * @param threshold Maximum bit errors allowed (default: 4)
 * @param sync_positions Output array of sync positions
 * @param max_positions Maximum positions to find
 * @return Number of sync positions found
 */
int p25_detect_frame_sync(
    p25_detector_t* detector,
    const int* bits,
    int num_bits,
    int threshold,
    int* sync_positions,
    int max_positions
);

/**
 * @brief Decode Network Identifier (NID) from P25 frame
 *
 * Extracts NAC and DUID from bits following frame sync.
 *
 * @param detector Pointer to detector structure
 * @param bits Bit stream
 * @param num_bits Total number of bits
 * @param sync_pos Frame sync position
 * @param nid Output NID structure
 * @return 0 on success, -1 if insufficient bits
 */
int p25_decode_nid(
    p25_detector_t* detector,
    const int* bits,
    int num_bits,
    int sync_pos,
    p25_nid_t* nid
);

/**
 * @brief Full P25 signal detection pipeline
 *
 * Complete detection: Demodulate -> Frame sync -> NID decode
 * Target: <15ms total processing time.
 *
 * @param detector Pointer to detector structure
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @param detections Output array for detections
 * @param max_detections Maximum detections to return
 * @return Number of P25 frames detected
 */
int p25_detect_signal(
    p25_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    p25_detection_t* detections,
    int max_detections
);

/**
 * @brief Estimate P25 signal strength
 *
 * @param iq_samples Complex IQ sample buffer
 * @param num_samples Number of samples
 * @return Signal strength in dB
 */
double p25_estimate_signal_strength(
    const cdouble_t* iq_samples,
    int num_samples
);

#endif /* P25_DETECTOR_H */
