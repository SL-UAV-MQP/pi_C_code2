/**
 * @file p25_detector.c
 * @brief P25 (Project 25) Land Mobile Radio signal detection
 *
 * Implements detection and decoding of P25 Phase 1 signals.
 *
 * Performance target: Demodulation <10ms, Detection <15ms
 */

#include "p25_detector.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Unwrap phase angle
 */
static void unwrap_phase(const double* phase_in, int length, double* phase_out) {
    phase_out[0] = phase_in[0];

    for (int i = 1; i < length; i++) {
        double diff = phase_in[i] - phase_in[i-1];

        // Unwrap when difference > pi
        while (diff > M_PI) {
            diff -= 2.0 * M_PI;
        }
        while (diff < -M_PI) {
            diff += 2.0 * M_PI;
        }

        phase_out[i] = phase_out[i-1] + diff;
    }
}

/**
 * @brief Compute instantaneous frequency from IQ samples
 */
static void compute_instantaneous_frequency(
    const cdouble_t* iq_samples,
    int num_samples,
    double sample_rate,
    double* inst_freq
) {
    // Calculate phase
    double* phase = (double*)malloc(num_samples * sizeof(double));
    double* phase_unwrapped = (double*)malloc(num_samples * sizeof(double));

    for (int i = 0; i < num_samples; i++) {
        phase[i] = carg(iq_samples[i]);
    }

    // Unwrap phase
    unwrap_phase(phase, num_samples, phase_unwrapped);

    // Compute instantaneous frequency: f = diff(phase) / (2*pi) * fs
    for (int i = 0; i < num_samples - 1; i++) {
        inst_freq[i] = (phase_unwrapped[i+1] - phase_unwrapped[i]) /
                      (2.0 * M_PI) * sample_rate;
    }
    inst_freq[num_samples - 1] = inst_freq[num_samples - 2];

    free(phase);
    free(phase_unwrapped);
}

/**
 * @brief Calculate Hamming distance between bit sequences
 */
static int hamming_distance(const int* bits1, const int* bits2, int length) {
    int distance = 0;
    for (int i = 0; i < length; i++) {
        if (bits1[i] != bits2[i]) {
            distance++;
        }
    }
    return distance;
}

/* ============================================================================
 * Initialization and Cleanup
 * ============================================================================ */

int p25_detector_init(p25_detector_t* detector, double sample_rate) {
    if (!detector) {
        return -1;
    }

    detector->sample_rate = sample_rate;
    detector->symbol_rate = P25_SYMBOL_RATE;
    detector->samples_per_symbol = (int)(sample_rate / P25_SYMBOL_RATE);

    // FIX: P25 Frame Sync pattern (48 bits) per TIA-102.BAAB standard
    // Hex: 0x5575F5FF77FF
    // Binary: 0101 0101 0111 0101 1111 0101 1111 1111 0111 0111 1111 1111
    int sync_pattern[P25_FRAME_SYNC_LENGTH] = {
        0,1,0,1, 0,1,0,1, 0,1,1,1, 0,1,0,1,  /* 0x5575 */
        1,1,1,1, 0,1,0,1, 1,1,1,1, 1,1,1,1,  /* 0xF5FF */
        0,1,1,1, 0,1,1,1, 1,1,1,1, 1,1,1,1   /* 0x77FF */
    };

    memcpy(detector->frame_sync, sync_pattern, P25_FRAME_SYNC_LENGTH * sizeof(int));

    // Allocate working buffers
    detector->symbol_buffer_size = P25_MAX_IQ_SAMPLES / detector->samples_per_symbol;
    detector->symbol_buffer = (int*)malloc(detector->symbol_buffer_size * sizeof(int));

    detector->bit_buffer_size = P25_MAX_BIT_BUFFER;
    detector->bit_buffer = (int*)malloc(detector->bit_buffer_size * sizeof(int));

    if (!detector->symbol_buffer || !detector->bit_buffer) {
        p25_detector_free(detector);
        return -1;
    }

    return 0;
}

void p25_detector_free(p25_detector_t* detector) {
    if (detector) {
        if (detector->symbol_buffer) {
            free(detector->symbol_buffer);
            detector->symbol_buffer = NULL;
        }
        if (detector->bit_buffer) {
            free(detector->bit_buffer);
            detector->bit_buffer = NULL;
        }
    }
}

/* ============================================================================
 * Demodulation and Symbol Processing
 * ============================================================================ */

int p25_demodulate_c4fm(
    p25_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int* symbols,
    int max_symbols
) {
    // demodulate_c4fm
    if (!detector || !iq_samples || !symbols) {
        return 0;
    }

    // Calculate instantaneous frequency
    double* inst_freq = (double*)malloc(num_samples * sizeof(double));
    if (!inst_freq) {
        return 0;
    }

    compute_instantaneous_frequency(iq_samples, num_samples,
                                    detector->sample_rate, inst_freq);

    // Resample to symbol rate
    int num_symbols = num_samples / detector->samples_per_symbol;
    if (num_symbols > max_symbols) {
        num_symbols = max_symbols;
    }

    // Average frequency over symbol period and map to symbols
    for (int i = 0; i < num_symbols; i++) {
        int start_idx = i * detector->samples_per_symbol;
        int end_idx = start_idx + detector->samples_per_symbol;

        if (end_idx > num_samples) {
            break;
        }

        // Average frequency over symbol period
        double avg_freq = 0.0;
        for (int j = start_idx; j < end_idx; j++) {
            avg_freq += inst_freq[j];
        }
        avg_freq /= detector->samples_per_symbol;

        // FIX: P25 C4FM symbol mapping with Gray coding per TIA-102.BAAA
        // Deviation levels: +1800 Hz, +600 Hz, -600 Hz, -1800 Hz
        // Gray-coded dibits:  01 (1),   00 (0),   10 (2),   11 (3)
        if (avg_freq > 1200.0) {
            symbols[i] = 1;       /* +1800 Hz → dibit 01 */
        } else if (avg_freq > 0.0) {
            symbols[i] = 0;       /* +600 Hz  → dibit 00 */
        } else if (avg_freq > -1200.0) {
            symbols[i] = 2;       /* -600 Hz  → dibit 10 */
        } else {
            symbols[i] = 3;       /* -1800 Hz → dibit 11 */
        }
    }

    free(inst_freq);
    return num_symbols;
}

int p25_symbols_to_bits(
    const int* symbols,
    int num_symbols,
    int* bits
) {
    // symbols_to_bits
    if (!symbols || !bits) {
        return 0;
    }

    // Convert dibits to bit stream
    for (int i = 0; i < num_symbols; i++) {
        bits[2*i] = (symbols[i] >> 1) & 1;      // MSB
        bits[2*i + 1] = symbols[i] & 1;         // LSB
    }

    return num_symbols * 2;
}

/* ============================================================================
 * Frame Sync Detection
 * ============================================================================ */

int p25_detect_frame_sync(
    p25_detector_t* detector,
    const int* bits,
    int num_bits,
    int threshold,
    int* sync_positions,
    int max_positions
) {
    // detect_frame_sync
    if (!detector || !bits || !sync_positions) {
        return 0;
    }

    int num_positions = 0;
    int sync_length = P25_FRAME_SYNC_LENGTH;

    // Scan through bit stream looking for frame sync pattern
    for (int i = 0; i < num_bits - sync_length && num_positions < max_positions; i++) {
        // Calculate Hamming distance
        int hamming_dist = hamming_distance(&bits[i], detector->frame_sync, sync_length);

        // Check if within threshold
        if (hamming_dist <= threshold) {
            sync_positions[num_positions] = i;
            num_positions++;
        }
    }

    return num_positions;
}

/* ============================================================================
 * NID Decoding
 * ============================================================================ */

int p25_decode_nid(
    p25_detector_t* detector,
    const int* bits,
    int num_bits,
    int sync_pos,
    p25_nid_t* nid
) {
    // decode_nid
    if (!detector || !bits || !nid) {
        return -1;
    }

    // NID comes after frame sync
    int nid_start = sync_pos + P25_FRAME_SYNC_LENGTH;
    int nid_length = P25_NID_LENGTH;

    if (nid_start + nid_length > num_bits) {
        return -1;
    }

    const int* nid_bits = &bits[nid_start];

    // Extract NAC (Network Access Code) - first 12 bits
    uint16_t nac = 0;
    for (int i = 0; i < 12; i++) {
        nac = (nac << 1) | nid_bits[i];
    }

    // Extract DUID (Data Unit ID) - 4 bits at position 12-15
    uint8_t duid = 0;
    for (int i = 12; i < 16; i++) {
        duid = (duid << 1) | nid_bits[i];
    }

    nid->nac = nac;
    nid->duid = duid;
    nid->position = sync_pos;

    return 0;
}

/* ============================================================================
 * Full Detection Pipeline
 * ============================================================================ */

int p25_detect_signal(
    p25_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    p25_detection_t* detections,
    int max_detections
) {
    // detect_p25_signal
    if (!detector || !iq_samples || !detections) {
        return 0;
    }

    int num_detections = 0;

    // Step 1: Demodulate C4FM
    int num_symbols = p25_demodulate_c4fm(detector, iq_samples, num_samples,
                                          detector->symbol_buffer,
                                          detector->symbol_buffer_size);

    if (num_symbols == 0) {
        return 0;
    }

    // Step 2: Convert to bits
    int num_bits = p25_symbols_to_bits(detector->symbol_buffer, num_symbols,
                                       detector->bit_buffer);

    if (num_bits == 0) {
        return 0;
    }

    // Step 3: Find frame sync
    int* sync_positions = (int*)malloc(P25_MAX_DETECTIONS * sizeof(int));
    if (!sync_positions) {
        return 0;
    }

    int num_syncs = p25_detect_frame_sync(detector, detector->bit_buffer, num_bits,
                                          4, sync_positions, P25_MAX_DETECTIONS);

    // Step 4: Decode NIDs
    for (int i = 0; i < num_syncs && num_detections < max_detections; i++) {
        p25_nid_t nid;

        int nid_result = p25_decode_nid(detector, detector->bit_buffer, num_bits,
                                        sync_positions[i], &nid);

        if (nid_result == 0) {
            strcpy(detections[num_detections].signal_type, "P25");
            detections[num_detections].nac = nid.nac;
            detections[num_detections].duid = nid.duid;
            detections[num_detections].timing = nid.position;
            strcpy(detections[num_detections].modulation, "C4FM");

            num_detections++;
        }
    }

    free(sync_positions);

    return num_detections;
}

/* ============================================================================
 * Signal Strength Estimation
 * ============================================================================ */

double p25_estimate_signal_strength(
    const cdouble_t* iq_samples,
    int num_samples
) {
    // estimate_signal_strength
    if (!iq_samples || num_samples == 0) {
        return -100.0;  // Very weak signal
    }

    // Calculate average power
    double power = 0.0;
    for (int i = 0; i < num_samples; i++) {
        double mag = cabs(iq_samples[i]);
        power += mag * mag;
    }
    power /= num_samples;

    // Convert to dB
    double power_db = 10.0 * log10(power + 1e-12);

    return power_db;
}
