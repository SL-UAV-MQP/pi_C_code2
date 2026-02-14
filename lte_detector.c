/**
 * @file lte_detector.c
 * @brief LTE/5G signal detection and Cell ID extraction
 *
 * Implements PSS/SSS detection and Physical Cell ID extraction for LTE signals.
 * Based on 3GPP TS 36.211 specifications.
 *
 * 3GPP TS 36.211 Section 6.11: Synchronization Signals
 * =====================================================
 *
 * LTE uses two synchronization signals for cell search:
 * 1. PSS (Primary Synchronization Signal) - Section 6.11.1
 *    - Based on Zadoff-Chu sequences
 *    - 3 possible sequences for N_ID_2 = 0, 1, 2
 *    - Root indices: 25, 29, 34
 *
 * 2. SSS (Secondary Synchronization Signal) - Section 6.11.2
 *    - Based on maximum-length sequences (m-sequences)
 *    - 168 possible sequences for N_ID_1 = 0, 1, ..., 167
 *    - Combined with PSS to identify 504 unique cells
 *
 * Physical Cell ID Calculation:
 *   N_ID_cell = 3 * N_ID_1 + N_ID_2
 *   Range: 0 to 503 (504 unique cell identities)
 *
 * SSS M-Sequence Generation (Section 6.11.2.1):
 *   - Characteristic polynomial: x^5 + x^2 + 1
 *   - Sequence length: 31
 *   - Two indices m0, m1 derived from N_ID_1
 *   - Scrambling with c0, c1 (based on N_ID_2)
 *   - Additional scrambling with z0, z1 (based on m0, m1)
 *
 * Performance target: PSS detection <20ms, Full detection <30ms
 * @updated 2025-12-11 - Fixed SSS generation per 3GPP TS 36.211 Section 6.11.2
 */

#include "lte_detector.h"
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
 * @brief Generate Zadoff-Chu sequence for PSS
 * generate_pss_sequences
 */
static void generate_zadoff_chu(int root_index, int length, cdouble_t* sequence) {
    // FIX: Per 3GPP TS 36.211 Section 6.11.1.1, PSS is generated from
    // Zadoff-Chu sequence of length N_ZC = 63 mapped to 62 subcarriers.
    // The DC subcarrier (n = 31 in 0-indexed Zadoff-Chu, mapped to k=0)
    // is excluded. The mapping is:
    //   d(n) = ZC(n)     for n = 0, 1, ..., 30
    //   d(n) = ZC(n+1)   for n = 31, 32, ..., 61
    for (int n = 0; n < length; n++) {
        int zc_n = (n < 31) ? n : (n + 1);  /* Skip n=31 (DC subcarrier) */
        double phase = -M_PI * root_index * zc_n * (zc_n + 1) / 63.0;
        sequence[n] = cexp(I * phase);
    }
}

/**
 * @brief Generate M-sequence for SSS according to 3GPP TS 36.211 Section 6.11.2
 *
 * Generates a maximum-length sequence (m-sequence) of length 31 using the
 * characteristic polynomial x^5 + x^2 + 1.
 *
 * The m-sequence is defined by:
 *   x(i+5) = (x(i+2) + x(i)) mod 2
 *
 * with initial state [0 0 0 0 1] (all zeros except the last position).
 *
 * @param sequence Output array of length 31 containing the m-sequence (0 or 1)
 */
static void generate_m_sequence_base(int* sequence) {
    // Initial state: x[0..4] = [0, 0, 0, 0, 1]
    int x[5] = {0, 0, 0, 0, 1};

    // Generate 31 bits using the characteristic polynomial x^5 + x^2 + 1
    for (int i = 0; i < 31; i++) {
        sequence[i] = x[0];

        // Feedback: x(i+5) = (x(i+2) + x(i)) mod 2
        int feedback = (x[2] + x[0]) % 2;

        // Shift register: move all bits left
        for (int j = 0; j < 4; j++) {
            x[j] = x[j + 1];
        }
        x[4] = feedback;
    }
}

/**
 * @brief Generate cyclic shift of m-sequence
 *
 * @param m_seq Base m-sequence of length 31
 * @param shift Cyclic shift amount
 * @param output Output shifted sequence of length 31
 */
static void cyclic_shift_m_sequence(const int* m_seq, int shift, int* output) {
    for (int i = 0; i < 31; i++) {
        output[i] = m_seq[(i + shift) % 31];
    }
}

/**
 * @brief Generate scrambling sequence c for SSS (3GPP TS 36.211 Section 6.11.2.1)
 *
 * The scrambling sequences c0 and c1 are based on cyclic shifts of the m-sequence:
 *   c0(n) = x((n + N_ID_2) mod 31)
 *   c1(n) = x((n + N_ID_2 + 3) mod 31)
 *
 * @param n_id_2 Physical layer cell identity group (0, 1, or 2)
 * @param c0 Output scrambling sequence c0 (length 31)
 * @param c1 Output scrambling sequence c1 (length 31)
 */
static void generate_sss_scrambling(int n_id_2, int* c0, int* c1) {
    int x[31];
    generate_m_sequence_base(x);

    // c0(n) = x((n + N_ID_2) mod 31)
    cyclic_shift_m_sequence(x, n_id_2, c0);

    // c1(n) = x((n + N_ID_2 + 3) mod 31)
    cyclic_shift_m_sequence(x, (n_id_2 + 3) % 31, c1);
}

/**
 * @brief Generate z scrambling sequence for SSS (3GPP TS 36.211 Section 6.11.2.1)
 *
 * The z sequences z0 and z1 are defined as:
 *   z0(n) = x((n + (m0 mod 8)) mod 31)
 *   z1(n) = x((n + (m1 mod 8)) mod 31)
 *
 * @param m0 Parameter derived from N_ID_1
 * @param m1 Parameter derived from N_ID_1
 * @param z0 Output scrambling sequence z0 (length 31)
 * @param z1 Output scrambling sequence z1 (length 31)
 */
static void generate_sss_z_sequences(int m0, int m1, int* z0, int* z1) {
    int x[31];
    generate_m_sequence_base(x);

    // z0(n) = x((n + (m0 mod 8)) mod 31)
    cyclic_shift_m_sequence(x, m0 % 8, z0);

    // z1(n) = x((n + (m1 mod 8)) mod 31)
    cyclic_shift_m_sequence(x, m1 % 8, z1);
}

/**
 * @brief Generate SSS sequence with N_ID_2 scrambling applied
 *
 * This function generates a complete SSS sequence for a given N_ID_1 and N_ID_2.
 * Used during SSS detection to correlate against received signal.
 *
 * @param n_id_1 Physical layer cell-identity group (0 to 167)
 * @param n_id_2 Physical layer identity within group (0, 1, 2)
 * @param sss_seq Output SSS sequence (length 62)
 */
static void generate_sss_with_scrambling(int n_id_1, int n_id_2, double* sss_seq) {
    // Generate base m-sequence
    int x[31];
    generate_m_sequence_base(x);

    // Calculate m0 and m1 according to 3GPP TS 36.211 Table 6.11.2.1-1
    int q_prime = n_id_1 / 30;
    int q = (n_id_1 + q_prime * (q_prime + 1) / 2) / 30;
    int m_prime = n_id_1 + q * (q + 1) / 2;
    int m0 = m_prime % 31;
    int m1 = (m0 + (m_prime / 31) + 1) % 31;

    // Generate s0 and s1 sequences
    int s0_m0[31], s1_m1[31];
    cyclic_shift_m_sequence(x, m0, s0_m0);
    cyclic_shift_m_sequence(x, m1, s1_m1);

    // Generate z sequences
    int z0_m0[31], z1_m1[31];
    generate_sss_z_sequences(m0, m1, z0_m0, z1_m1);

    // Generate c sequences with actual N_ID_2
    int c0[31], c1[31];
    generate_sss_scrambling(n_id_2, c0, c1);

    // Generate SSS sequence for subframe 0
    for (int n = 0; n < 31; n++) {
        double s0_val = 1.0 - 2.0 * s0_m0[n];
        double s1_val = 1.0 - 2.0 * s1_m1[n];
        double c0_val = 1.0 - 2.0 * c0[n];
        double c1_val = 1.0 - 2.0 * c1[n];
        double z0_val = 1.0 - 2.0 * z0_m0[n];

        sss_seq[2*n]   = s0_val * c0_val;
        sss_seq[2*n+1] = s1_val * c1_val * z0_val;
    }
}

/**
 * @brief Cross-correlate complex sequences
 */
static void cross_correlate_complex(
    const cdouble_t* signal,
    int signal_len,
    const cdouble_t* reference,
    int ref_len,
    cdouble_t* output
) {
    int output_len = signal_len - ref_len + 1;

    for (int lag = 0; lag < output_len; lag++) {
        output[lag] = 0.0;
        for (int i = 0; i < ref_len; i++) {
            output[lag] += signal[lag + i] * conj(reference[i]);
        }
    }
}

/* ============================================================================
 * Initialization and Cleanup
 * ============================================================================ */

int lte_detector_init(lte_detector_t* detector, double sample_rate) {
    if (!detector) {
        return -1;
    }

    detector->sample_rate = sample_rate;
    detector->fft_size = 2048;
    detector->subcarrier_spacing = 15e3;  // 15 kHz for LTE

    // Generate PSS sequences 
    if (lte_generate_pss_sequences(detector) != 0) {
        return -1;
    }

    // Generate SSS sequences
    if (lte_generate_sss_sequences(detector) != 0) {
        return -1;
    }

    // Allocate correlation buffers
    detector->correlation_buffer_size = LTE_MAX_IQ_SAMPLES;
    detector->correlation_buffer = (cdouble_t*)malloc(
        detector->correlation_buffer_size * sizeof(cdouble_t));
    detector->correlation_magnitude = (double*)malloc(
        detector->correlation_buffer_size * sizeof(double));

    if (!detector->correlation_buffer || !detector->correlation_magnitude) {
        lte_detector_free(detector);
        return -1;
    }

    return 0;
}

void lte_detector_free(lte_detector_t* detector) {
    if (detector) {
        if (detector->correlation_buffer) {
            free(detector->correlation_buffer);
            detector->correlation_buffer = NULL;
        }
        if (detector->correlation_magnitude) {
            free(detector->correlation_magnitude);
            detector->correlation_magnitude = NULL;
        }
    }
}

/* ============================================================================
 * Sequence Generation
 * ============================================================================ */

/**
 * 3GPP TS 36.211 Section 6.11.2: SSS Sequence Generation
 * ========================================================
 *
 * The SSS is constructed from interleaved sequences in the frequency domain:
 *
 * Subframe 0: d(2n)   = s0_m0(n) * c0(n)
 *            d(2n+1) = s1_m1(n) * c1(n) * z0_m0(n)
 *
 * Subframe 5: d(2n)   = s1_m1(n) * c0(n)
 *            d(2n+1) = s0_m0(n) * c1(n) * z1_m1(n)
 *
 * where n = 0, 1, ..., 30 (total 62 subcarriers)
 *
 * Sequence Definitions:
 * ---------------------
 * 1. Base m-sequence x(i):
 *    - Characteristic polynomial: x^5 + x^2 + 1
 *    - Recurrence: x(i+5) = (x(i+2) + x(i)) mod 2
 *    - Initial state: [0 0 0 0 1]
 *    - Length: 31
 *
 * 2. Indices m0, m1 (Table 6.11.2.1-1):
 *    - q' = floor(N_ID_1 / 30)
 *    - q  = floor((N_ID_1 + q'*(q'+1)/2) / 30)
 *    - m' = N_ID_1 + q*(q+1)/2
 *    - m0 = m' mod 31
 *    - m1 = (m0 + floor(m'/31) + 1) mod 31
 *
 * 3. Sequences s0 and s1:
 *    - s0_m0(n) = 1 - 2*x((n + m0) mod 31)    (cyclic shift by m0)
 *    - s1_m1(n) = 1 - 2*x((n + m1) mod 31)    (cyclic shift by m1)
 *
 * 4. Scrambling sequences c0, c1:
 *    - c0(n) = 1 - 2*x((n + N_ID_2) mod 31)
 *    - c1(n) = 1 - 2*x((n + N_ID_2 + 3) mod 31)
 *
 * 5. Scrambling sequences z0, z1:
 *    - z0_m0(n) = 1 - 2*x((n + (m0 mod 8)) mod 31)
 *    - z1_m1(n) = 1 - 2*x((n + (m1 mod 8)) mod 31)
 *
 * BPSK Mapping:
 * -------------
 * Binary 0 -> +1
 * Binary 1 -> -1
 * Formula: s = 1 - 2*bit
 *
 * References:
 * -----------
 * - 3GPP TS 36.211 V13.2.0 Section 6.11.2
 * - MATLAB lteSSS function
 * - srsRAN SSS generation implementation
 */

int lte_generate_pss_sequences(lte_detector_t* detector) {
    // generate_pss_sequences
    if (!detector) {
        return -1;
    }

    // PSS root indices for N_ID_2 = 0, 1, 2
    int roots[LTE_NUM_PSS_SEQUENCES] = {25, 29, 34};

    for (int n_id_2 = 0; n_id_2 < LTE_NUM_PSS_SEQUENCES; n_id_2++) {
        generate_zadoff_chu(roots[n_id_2], LTE_SYNC_SEQ_LENGTH,
                           detector->pss_sequences[n_id_2]);
    }

    return 0;
}

int lte_generate_sss_sequences(lte_detector_t* detector) {
    /**
     * LTE SSS Generation according to 3GPP TS 36.211 Section 6.11.2
     *
     * The SSS is constructed from two sequences that are generated based on:
     * - Physical layer cell identity (N_ID_cell = 3*N_ID_1 + N_ID_2)
     * - N_ID_1: Physical layer cell-identity group (0 to 167)
     * - N_ID_2: Physical layer identity within group (0, 1, 2)
     *
     * Algorithm:
     * 1. Generate base m-sequence x of length 31
     * 2. For each N_ID_1, compute indices m0 and m1
     * 3. Generate sequences s0 and s1 from cyclic shifts of x
     * 4. Apply scrambling with c0, c1, z0, z1 sequences
     * 5. Construct SSS in subframe 0 and 5 (different scrambling)
     *
     * SSS Formula (3GPP TS 36.211 Section 6.11.2.1):
     *   Subframe 0: d(2n)   = s0_m0(n) * c0(n)
     *              d(2n+1) = s1_m1(n) * c1(n) * z0_m0(n)
     *
     *   Subframe 5: d(2n)   = s1_m1(n) * c0(n)
     *              d(2n+1) = s0_m0(n) * c1(n) * z1_m1(n)
     *
     * where n = 0, 1, ..., 30 (31 positions)
     */

    if (!detector) {
        return -1;
    }

    // Generate base m-sequence
    int x[31];
    generate_m_sequence_base(x);

    // Generate SSS for each N_ID_1 (0 to 167)
    for (int n_id_1 = 0; n_id_1 < LTE_NUM_SSS_SEQUENCES; n_id_1++) {
        // Calculate m0 and m1 indices according to 3GPP TS 36.211 Table 6.11.2.1-1
        // q_prime = floor(N_ID_1 / 30)
        int q_prime = n_id_1 / 30;

        // q = floor((N_ID_1 + q_prime * (q_prime + 1) / 2) / 30)
        int q = (n_id_1 + q_prime * (q_prime + 1) / 2) / 30;

        // m_prime = N_ID_1 + q * (q + 1) / 2
        int m_prime = n_id_1 + q * (q + 1) / 2;

        // m0 = m' mod 31
        int m0 = m_prime % 31;

        // m1 = (m0 + floor(m'/31) + 1) mod 31
        int m1 = (m0 + (m_prime / 31) + 1) % 31;

        // Generate s0 and s1 sequences (cyclic shifts of base m-sequence)
        int s0_m0[31], s1_m1[31];
        cyclic_shift_m_sequence(x, m0, s0_m0);
        cyclic_shift_m_sequence(x, m1, s1_m1);

        // Generate z0 and z1 scrambling sequences
        int z0_m0[31], z1_m1[31];
        generate_sss_z_sequences(m0, m1, z0_m0, z1_m1);

        // Generate scrambling sequences c0 and c1 for each N_ID_2
        // For SSS sequence storage, we use N_ID_2 = 0 as reference
        // (actual N_ID_2 scrambling will be applied during detection)
        int c0[31], c1[31];
        generate_sss_scrambling(0, c0, c1);  // N_ID_2 = 0 reference

        // Generate SSS sequence for subframe 0
        // d(2n)   = s0_m0(n) * c0(n)
        // d(2n+1) = s1_m1(n) * c1(n) * z0_m0(n)
        for (int n = 0; n < 31 && (2*n+1) < LTE_SYNC_SEQ_LENGTH; n++) {
            // Convert binary (0,1) to BPSK (-1,+1): s = 1 - 2*bit
            double s0_val = 1.0 - 2.0 * s0_m0[n];
            double s1_val = 1.0 - 2.0 * s1_m1[n];
            double c0_val = 1.0 - 2.0 * c0[n];
            double c1_val = 1.0 - 2.0 * c1[n];
            double z0_val = 1.0 - 2.0 * z0_m0[n];

            // Subframe 0 SSS
            detector->sss_sequences[n_id_1][2*n]   = s0_val * c0_val;
            detector->sss_sequences[n_id_1][2*n+1] = s1_val * c1_val * z0_val;
        }

        // Note: In practice, subframe 5 uses different interleaving:
        // d(2n)   = s1_m1(n) * c0(n)
        // d(2n+1) = s0_m0(n) * c1(n) * z1_m1(n)
        // For detection, we use subframe 0 pattern as reference
    }

    return 0;
}

/* ============================================================================
 * PSS/SSS Detection
 * ============================================================================ */

int lte_detect_pss(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    double threshold,
    pss_detection_t* detections,
    int max_detections
) {
    // detect_pss
    if (!detector || !iq_samples || !detections) {
        return 0;
    }

    int num_detections = 0;

    // Correlate with each PSS sequence
    for (int n_id_2 = 0; n_id_2 < LTE_NUM_PSS_SEQUENCES && num_detections < max_detections; n_id_2++) {
        int corr_len = num_samples - LTE_SYNC_SEQ_LENGTH + 1;
        if (corr_len <= 0) continue;

        // Cross-correlation
        cross_correlate_complex(iq_samples, num_samples,
                               detector->pss_sequences[n_id_2],
                               LTE_SYNC_SEQ_LENGTH,
                               detector->correlation_buffer);

        // Normalize 
        double pss_norm = 0.0;
        for (int i = 0; i < LTE_SYNC_SEQ_LENGTH; i++) {
            pss_norm += cabs(detector->pss_sequences[n_id_2][i]) *
                       cabs(detector->pss_sequences[n_id_2][i]);
        }
        pss_norm = sqrt(pss_norm);

        for (int i = 0; i < corr_len; i++) {
            detector->correlation_magnitude[i] = cabs(detector->correlation_buffer[i]) /
                                                (pss_norm * sqrt(LTE_SYNC_SEQ_LENGTH));
        }

        // Find peaks above threshold 
        int* peaks;
        int num_peaks = lte_find_peaks(detector->correlation_magnitude, corr_len,
                                       threshold, 100, NULL, 0);

        if (num_peaks > 0) {
            peaks = (int*)malloc(num_peaks * sizeof(int));
            lte_find_peaks(detector->correlation_magnitude, corr_len,
                          threshold, 100, peaks, num_peaks);

            for (int i = 0; i < num_peaks && num_detections < max_detections; i++) {
                detections[num_detections].n_id_2 = n_id_2;
                detections[num_detections].timing_offset = peaks[i];
                detections[num_detections].correlation =
                    detector->correlation_magnitude[peaks[i]];
                strcpy(detections[num_detections].type, "PSS");
                num_detections++;
            }

            free(peaks);
        }
    }

    return num_detections;
}

int lte_detect_sss(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int pss_timing,
    int n_id_2,
    double threshold,
    sss_detection_t* sss_det
) {
    /**
     * SSS Detection with proper N_ID_2 scrambling
     *
     * According to 3GPP TS 36.211 Section 6.11.2:
     * - SSS is located in the same subframe as PSS
     * - In FDD: SSS is in subframe 0 and 5 (slots 0 and 10)
     * - SSS is mapped to the same OFDM symbol as PSS in adjacent slots
     *
     * Detection algorithm:
     * 1. Extract SSS region from IQ samples (before PSS)
     * 2. For each N_ID_1 (0 to 167):
     *    a. Generate SSS with proper N_ID_2 scrambling
     *    b. Correlate with received signal
     *    c. Normalize correlation
     * 3. Select N_ID_1 with highest correlation
     * 4. Check if correlation exceeds threshold
     */
    if (!detector || !iq_samples || !sss_det) {
        return -1;
    }

    // SSS is located before PSS
    int sss_start = (pss_timing - 100 > 0) ? pss_timing - 100 : 0;
    int sss_end = pss_timing;
    int sss_len = sss_end - sss_start;

    if (sss_len < LTE_SYNC_SEQ_LENGTH) {
        return -1;
    }

    // Correlate with all SSS sequences 
    double best_correlation = 0.0;
    int best_n_id_1 = -1;

    // Generate SSS with proper N_ID_2 scrambling for each N_ID_1
    double sss_reference[62];

    for (int n_id_1 = 0; n_id_1 < LTE_NUM_SSS_SEQUENCES; n_id_1++) {
        // Generate SSS with correct N_ID_2 scrambling
        generate_sss_with_scrambling(n_id_1, n_id_2, sss_reference);

        cdouble_t corr_sum = 0.0;

        for (int i = 0; i < LTE_SYNC_SEQ_LENGTH && (sss_start + i) < num_samples; i++) {
            // Correlate with real-valued SSS sequence
            corr_sum += iq_samples[sss_start + i] * sss_reference[i];
        }

        double corr_mag = cabs(corr_sum);

        // Normalize by SSS energy
        double sss_norm = 0.0;
        for (int i = 0; i < LTE_SYNC_SEQ_LENGTH; i++) {
            sss_norm += sss_reference[i] * sss_reference[i];
        }
        sss_norm = sqrt(sss_norm);

        double corr_normalized = corr_mag / (sss_norm * sqrt(LTE_SYNC_SEQ_LENGTH));

        if (corr_normalized > best_correlation) {
            best_correlation = corr_normalized;
            best_n_id_1 = n_id_1;
        }
    }

    // Check threshold
    if (best_correlation > threshold) {
        sss_det->n_id_1 = best_n_id_1;
        sss_det->n_id_2 = n_id_2;
        sss_det->correlation = best_correlation;
        strcpy(sss_det->type, "SSS");
        return 0;
    }

    return -1;
}

/* ============================================================================
 * Cell ID Extraction and Full Detection
 * ============================================================================ */

int lte_extract_cell_id(int n_id_1, int n_id_2) {
    // extract_cell_id
    // Cell ID = 3 * N_ID_1 + N_ID_2
    return 3 * n_id_1 + n_id_2;
}

int lte_detect_signal(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    lte_cell_detection_t* cells,
    int max_cells
) {
    // detect_lte_signal
    if (!detector || !iq_samples || !cells) {
        return 0;
    }

    int num_cells = 0;

    // Step 1: Detect PSS
    pss_detection_t* pss_detections = (pss_detection_t*)malloc(
        LTE_MAX_DETECTIONS * sizeof(pss_detection_t));

    if (!pss_detections) {
        return 0;
    }

    int num_pss = lte_detect_pss(detector, iq_samples, num_samples, 0.6,
                                 pss_detections, LTE_MAX_DETECTIONS);

    // Step 2: For each PSS, detect SSS
    for (int i = 0; i < num_pss && num_cells < max_cells; i++) {
        sss_detection_t sss_det;

        int sss_result = lte_detect_sss(detector, iq_samples, num_samples,
                                        pss_detections[i].timing_offset,
                                        pss_detections[i].n_id_2,
                                        0.5, &sss_det);

        if (sss_result == 0) {
            // Calculate Cell ID
            int cell_id = lte_extract_cell_id(sss_det.n_id_1, sss_det.n_id_2);

            cells[num_cells].cell_id = cell_id;
            cells[num_cells].n_id_1 = sss_det.n_id_1;
            cells[num_cells].n_id_2 = sss_det.n_id_2;
            cells[num_cells].timing_offset = pss_detections[i].timing_offset;
            cells[num_cells].pss_correlation = pss_detections[i].correlation;
            cells[num_cells].sss_correlation = sss_det.correlation;
            strcpy(cells[num_cells].signal_type, "LTE");

            num_cells++;
        }
    }

    free(pss_detections);

    return num_cells;
}

double lte_estimate_frequency_offset(
    lte_detector_t* detector,
    const cdouble_t* iq_samples,
    int num_samples,
    int pss_timing,
    int n_id_2
) {
    // estimate_frequency_offset
    if (!detector || !iq_samples) {
        return 0.0;
    }

    // Extract two consecutive PSS symbols (5ms apart in LTE)
    int symbol_spacing = (int)(0.005 * detector->sample_rate);
    int pss1_start = pss_timing;
    int pss2_start = pss_timing + symbol_spacing;

    if (pss2_start + LTE_SYNC_SEQ_LENGTH >= num_samples) {
        return 0.0;
    }

    // Calculate phase difference
    cdouble_t phase_sum = 0.0;
    for (int i = 0; i < LTE_SYNC_SEQ_LENGTH; i++) {
        cdouble_t pss1 = iq_samples[pss1_start + i];
        cdouble_t pss2 = iq_samples[pss2_start + i];
        phase_sum += pss2 * conj(pss1);
    }

    double phase_diff = carg(phase_sum);

    // Convert to frequency offset
    double freq_offset = phase_diff / (2.0 * M_PI * 0.005);

    return freq_offset;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

void lte_cross_correlate(
    const cdouble_t* signal,
    int signal_len,
    const cdouble_t* reference,
    int ref_len,
    cdouble_t* output
) {
    cross_correlate_complex(signal, signal_len, reference, ref_len, output);
}

int lte_find_peaks(
    const double* data,
    int len,
    double threshold,
    int min_distance,
    int* peaks,
    int max_peaks
) {
    if (!data) {
        return 0;
    }

    int num_peaks = 0;
    int* temp_peaks = (int*)malloc(len * sizeof(int));
    if (!temp_peaks) {
        return 0;
    }

    // Find local maxima above threshold
    for (int i = 1; i < len - 1; i++) {
        if (data[i] > data[i-1] && data[i] > data[i+1] && data[i] >= threshold) {
            temp_peaks[num_peaks++] = i;
        }
    }

    // Apply distance constraint
    if (peaks && max_peaks > 0) {
        int out_idx = 0;
        for (int i = 0; i < num_peaks && out_idx < max_peaks; i++) {
            bool too_close = false;
            for (int j = 0; j < out_idx; j++) {
                if (abs(temp_peaks[i] - peaks[j]) < min_distance) {
                    too_close = true;
                    break;
                }
            }
            if (!too_close) {
                peaks[out_idx++] = temp_peaks[i];
            }
        }
        num_peaks = out_idx;
    }

    free(temp_peaks);
    return num_peaks;
}
