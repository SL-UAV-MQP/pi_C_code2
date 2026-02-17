/**
 * @file signal_processor.c
 * @brief Signal preprocessing for SDR IQ data
 *
 * Implements FFT, PSD, spectrogram, band normalization, filtering,
 * peak detection, downsampling, and SNR estimation for real-time
 * signal processing on Raspberry Pi.
 *
 * Performance targets: FFT <5ms, PSD <10ms, filtering <5ms
 */

#include "signal_processor.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Create window function
 */
static void create_window(window_type_t type, int size, double* window) {
    /* BUG-04 fix: guard against size <= 1 to avoid division by zero */
    if (size <= 1) {
        if (size == 1) window[0] = 1.0;
        return;
    }

    switch (type) {
        case WINDOW_HAMMING:
            for (int n = 0; n < size; n++) {
                window[n] = 0.54 - 0.46 * cos(2.0 * M_PI * n / (size - 1));
            }
            break;

        case WINDOW_HANN:
            for (int n = 0; n < size; n++) {
                window[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / (size - 1)));
            }
            break;

        case WINDOW_BLACKMAN:
            for (int n = 0; n < size; n++) {
                double a0 = 0.42;
                double a1 = 0.5;
                double a2 = 0.08;
                window[n] = a0 - a1 * cos(2.0 * M_PI * n / (size - 1)) +
                           a2 * cos(4.0 * M_PI * n / (size - 1));
            }
            break;

        case WINDOW_KAISER:
            // Kaiser window with beta=8.6
            // I0(x) = sum_{k=0}^{N} ((x/2)^k / k!)^2
            {
                double beta = 8.6;
                double alpha = (size - 1) / 2.0;
                // Precompute I0(beta) denominator
                double i0_beta = 1.0;
                {
                    double term = 1.0;
                    double half_beta = beta / 2.0;
                    for (int k = 1; k <= 25; k++) {
                        term *= (half_beta / k);
                        i0_beta += term * term;
                    }
                }
                for (int n = 0; n < size; n++) {
                    double x = (n - alpha) / alpha;
                    double arg = beta * sqrt(1.0 - x * x);
                    // Compute I0(arg) using 25-term series
                    double i0_arg = 1.0;
                    double term = 1.0;
                    double half_arg = arg / 2.0;
                    for (int k = 1; k <= 25; k++) {
                        term *= (half_arg / k);
                        i0_arg += term * term;
                    }
                    window[n] = i0_arg / i0_beta;
                }
            }
            break;

        case WINDOW_RECTANGULAR:
        default:
            for (int n = 0; n < size; n++) {
                window[n] = 1.0;
            }
            break;
    }
}

/**
 * @brief Compute magnitude in dB from complex value
 */
static inline double magnitude_db(cdouble_t z) {
    double mag = cabs(z);
    return 20.0 * log10(mag + EPSILON);
}

/**
 * @brief Compute power in dB from complex value
 */
static inline double power_db(cdouble_t z) {
    double power = creal(z) * creal(z) + cimag(z) * cimag(z);
    return 10.0 * log10(power + EPSILON);
}

/* ============================================================================
 * Initialization and Cleanup
 * ============================================================================ */

signal_processor_t* signal_processor_init(
    double sample_rate,
    int fft_size,
    window_type_t window_type,
    double overlap
) {
    // Validate inputs
    if (sample_rate <= 0 || fft_size <= 0 || fft_size > MAX_FFT_SIZE) {
        return NULL;
    }

    if (overlap < 0.0 || overlap >= 1.0) {
        return NULL;
    }

    // Check if fft_size is power of 2
    if ((fft_size & (fft_size - 1)) != 0) {
        return NULL;
    }

    // Allocate processor structure
    signal_processor_t* processor = (signal_processor_t*)malloc(sizeof(signal_processor_t));
    if (!processor) {
        return NULL;
    }

    // Initialize parameters
    processor->sample_rate = sample_rate;
    processor->fft_size = fft_size;
    processor->window_type = window_type;
    processor->overlap = overlap;
    processor->hop_size = (int)(fft_size * (1.0 - overlap));

    // Allocate and create window
    processor->window = (double*)malloc(fft_size * sizeof(double));
    if (!processor->window) {
        free(processor);
        return NULL;
    }

    create_window(window_type, fft_size, processor->window);

    // Create cached FFTW plan and buffers (created once, reused for all FFT calls)
    processor->fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);
    processor->fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_size);
    if (!processor->fftw_in || !processor->fftw_out) {
        if (processor->fftw_in) fftw_free(processor->fftw_in);
        if (processor->fftw_out) fftw_free(processor->fftw_out);
        free(processor->window);
        free(processor);
        return NULL;
    }

    processor->fftw_forward = fftw_plan_dft_1d(fft_size, processor->fftw_in,
                                                processor->fftw_out,
                                                FFTW_FORWARD, FFTW_MEASURE);

    if (!processor->fftw_forward) {
        fftw_free(processor->fftw_in);
        fftw_free(processor->fftw_out);
        free(processor->window);
        free(processor);
        return NULL;
    }

    return processor;
}

void signal_processor_free(signal_processor_t* processor) {
    if (processor) {
        if (processor->fftw_forward) {
            fftw_destroy_plan(processor->fftw_forward);
        }
        if (processor->fftw_in) fftw_free(processor->fftw_in);
        if (processor->fftw_out) fftw_free(processor->fftw_out);
        if (processor->window) {
            free(processor->window);
        }
        free(processor);
    }
}

/* ============================================================================
 * FFT and Spectral Analysis
 * ============================================================================ */

music_status_t signal_processor_compute_fft(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    fft_result_t* result
) {
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (num_samples < processor->fft_size) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    int N = processor->fft_size;
    fftw_complex* in = processor->fftw_in;
    fftw_complex* out = processor->fftw_out;

    // Apply window and copy to cached input buffer
    for (int i = 0; i < N; i++) {
        double window_val = processor->window[i];
        in[i][0] = creal(iq_samples[i]) * window_val;
        in[i][1] = cimag(iq_samples[i]) * window_val;
    }

    // Execute cached FFT plan
    fftw_execute(processor->fftw_forward);

    // Compute frequency axis and magnitude
    double df = processor->sample_rate / N;
    for (int i = 0; i < N; i++) {
        int k = (i + N/2) % N;
        result->frequencies[i] = (k - N/2) * df;
        cdouble_t val = out[k][0] + I * out[k][1];
        result->magnitude_db[i] = magnitude_db(val);
    }

    result->num_bins = N;

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_compute_psd(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    psd_result_t* result
) {
    // compute_psd using Welch's method
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = processor->fft_size;
    int hop = processor->hop_size;
    int num_segments = (num_samples - N) / hop + 1;

    if (num_segments < 1) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Allocate working buffer for accumulation
    double* psd_accumulator = (double*)calloc(N, sizeof(double));
    if (!psd_accumulator) {
        return MUSIC_ERROR_MEMORY;
    }

    fftw_complex* in = processor->fftw_in;
    fftw_complex* out = processor->fftw_out;

    // Window normalization factor
    double window_norm = 0.0;
    for (int i = 0; i < N; i++) {
        window_norm += processor->window[i] * processor->window[i];
    }

    // Process each segment (Welch's method)
    for (int seg = 0; seg < num_segments; seg++) {
        int offset = seg * hop;

        // Apply window
        for (int i = 0; i < N; i++) {
            double window_val = processor->window[i];
            in[i][0] = creal(iq_samples[offset + i]) * window_val;
            in[i][1] = cimag(iq_samples[offset + i]) * window_val;
        }

        // Execute cached plan
        fftw_execute(processor->fftw_forward);

        // Accumulate power spectrum
        for (int i = 0; i < N; i++) {
            double power = out[i][0] * out[i][0] + out[i][1] * out[i][1];
            psd_accumulator[i] += power;
        }
    }

    // Average and normalize
    double scale = 1.0 / (num_segments * window_norm * processor->sample_rate);
    double df = processor->sample_rate / N;

    for (int i = 0; i < N; i++) {
        int k = (i + N/2) % N;
        result->frequencies[i] = (k - N/2) * df;
        double psd_linear = psd_accumulator[k] * scale;
        result->psd_db[i] = 10.0 * log10(psd_linear + EPSILON);
    }

    result->num_bins = N;

    free(psd_accumulator);

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_compute_spectrogram(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    spectrogram_result_t* result
) {
    // compute_spectrogram using STFT
    if (!processor || !iq_samples || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int N = processor->fft_size;
    int hop = processor->hop_size;
    int num_time_bins = (num_samples - N) / hop + 1;

    if (num_time_bins < 1) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    fftw_complex* in = processor->fftw_in;
    fftw_complex* out = processor->fftw_out;

    // Compute time and frequency axes
    double dt = (double)hop / processor->sample_rate;
    double df = processor->sample_rate / N;

    for (int t = 0; t < num_time_bins; t++) {
        result->times[t] = t * dt;
    }

    for (int f = 0; f < N; f++) {
        int k = (f + N/2) % N;
        result->frequencies[f] = (k - N/2) * df;
    }

    // Compute STFT using cached plan
    for (int t = 0; t < num_time_bins; t++) {
        int offset = t * hop;

        // Apply window
        for (int i = 0; i < N; i++) {
            double window_val = processor->window[i];
            in[i][0] = creal(iq_samples[offset + i]) * window_val;
            in[i][1] = cimag(iq_samples[offset + i]) * window_val;
        }

        // Execute cached plan
        fftw_execute(processor->fftw_forward);

        // Store magnitude in dB (column-major order)
        for (int f = 0; f < N; f++) {
            int k = (f + N/2) % N;
            cdouble_t val = out[k][0] + I * out[k][1];
            result->spectrogram_db[f * num_time_bins + t] = magnitude_db(val);
        }
    }

    result->num_time_bins = num_time_bins;
    result->num_freq_bins = N;

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Normalization and Filtering
 * ============================================================================ */

music_status_t signal_processor_normalize_band(
    const signal_processor_t* processor,
    cdouble_t* iq_samples,
    int num_samples,
    normalize_method_t method
) {
    // normalize_band
    if (!processor || !iq_samples) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (method == NORMALIZE_PSD) {
        // Normalize based on PSD
        psd_result_t* psd = psd_result_alloc(processor->fft_size);
        if (!psd) {
            return MUSIC_ERROR_MEMORY;
        }

        music_status_t status = signal_processor_compute_psd(processor, iq_samples,
                                                             num_samples, psd);
        if (status != MUSIC_SUCCESS) {
            psd_result_free(psd);
            return status;
        }

        // Compute mean power
        double mean_power_db = 0.0;
        for (int i = 0; i < psd->num_bins; i++) {
            mean_power_db += psd->psd_db[i];
        }
        mean_power_db /= psd->num_bins;
        double mean_power_linear = pow(10.0, mean_power_db / 10.0);

        // Current signal power
        double current_power = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            current_power += mag * mag;
        }
        current_power /= num_samples;

        // Scale factor
        double scale_factor = sqrt(mean_power_linear / (current_power + EPSILON));

        // Apply scaling
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] *= scale_factor;
        }

        psd_result_free(psd);

    } else if (method == NORMALIZE_MEAN) {
        // Zero mean, unit variance
        cdouble_t mean = 0.0;
        for (int i = 0; i < num_samples; i++) {
            mean += iq_samples[i];
        }
        mean /= num_samples;

        // Remove mean
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] -= mean;
        }

        // Compute standard deviation
        double variance = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            variance += mag * mag;
        }
        double std = sqrt(variance / num_samples);

        // Normalize to unit variance
        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] /= (std + EPSILON);
        }

    } else if (method == NORMALIZE_MAX) {
        // Normalize to max amplitude
        double max_amp = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double mag = cabs(iq_samples[i]);
            if (mag > max_amp) {
                max_amp = mag;
            }
        }

        for (int i = 0; i < num_samples; i++) {
            iq_samples[i] /= (max_amp + EPSILON);
        }
    }

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_design_bandpass(
    double sample_rate,
    double low_freq,
    double high_freq,
    int order,
    filter_coeffs_t* coeffs
) {
    if (!coeffs || order <= 0) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    /* Limit to order 4 lowpass prototype => 8th order bandpass => 4 biquads.
     * b[]/a[] have 21 elements each, so max 7 biquad sections (7*3=21).
     * order field = number of biquad sections. */
    if (order > 4) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    double nyquist = sample_rate / 2.0;
    double low_norm = low_freq / nyquist;
    double high_norm = high_freq / nyquist;

    if (low_norm < 0.001) low_norm = 0.001;
    if (high_norm > 0.999) high_norm = 0.999;
    if (low_norm >= high_norm) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    /* ----------------------------------------------------------------
     * Step 1: Pre-warp digital frequencies to analog domain
     * ---------------------------------------------------------------- */
    double wL = 2.0 * sample_rate * tan(M_PI * low_freq / sample_rate);
    double wH = 2.0 * sample_rate * tan(M_PI * high_freq / sample_rate);
    double w0 = sqrt(wL * wH);          /* center frequency */
    double BW = wH - wL;                /* bandwidth */

    /* ----------------------------------------------------------------
     * Step 2: Analog Butterworth lowpass prototype poles (unit circle)
     * For order N, poles are at s_k = exp(j * pi * (2k+N+1)/(2N)),
     * k = 0..N-1. We only take the left-half-plane poles.
     * ---------------------------------------------------------------- */
    int N = order;   /* lowpass prototype order */
    cdouble_t lp_poles[4]; /* max order 4 */

    for (int k = 0; k < N; k++) {
        double angle = M_PI * (2.0*k + N + 1.0) / (2.0 * N);
        lp_poles[k] = cos(angle) + I * sin(angle);
    }

    /* ----------------------------------------------------------------
     * Step 3: Lowpass-to-bandpass transformation
     * Each LP pole p_k maps to two BP poles:
     *   s = (p_k * BW/2) +/- sqrt((p_k*BW/2)^2 - w0^2)
     * This doubles the number of poles: 2N total.
     * ---------------------------------------------------------------- */
    int num_bp_poles = 2 * N;
    cdouble_t bp_poles[8]; /* max 2*4=8 */

    for (int k = 0; k < N; k++) {
        cdouble_t half_bw_pk = lp_poles[k] * (BW / 2.0);
        cdouble_t disc = csqrt(half_bw_pk * half_bw_pk - w0 * w0);
        bp_poles[2*k]     = half_bw_pk + disc;
        bp_poles[2*k + 1] = half_bw_pk - disc;
    }

    /* Bandpass numerator: s^N (N zeros at origin) */

    /* ----------------------------------------------------------------
     * Step 4: Bilinear transform: s = 2*fs*(z-1)/(z+1)
     * Each analog pole p maps to digital pole: z = (1 + p/(2*fs)) / (1 - p/(2*fs))
     * Each zero at s=0 maps to z = -1
     * ---------------------------------------------------------------- */
    cdouble_t dz_poles[8];
    double fs2 = 2.0 * sample_rate;

    for (int k = 0; k < num_bp_poles; k++) {
        cdouble_t sp = bp_poles[k] / fs2;
        dz_poles[k] = (1.0 + sp) / (1.0 - sp);
    }

    /* Bilinear transform maps:
     * s=0 -> z=+1,  s=inf -> z=-1.
     * So the N zeros at s=0 map to N zeros at z=+1,
     * and N implicit zeros at s=inf map to N zeros at z=-1.
     * Each biquad section gets one zero at z=+1 and one at z=-1,
     * giving numerator (1 - z^-2) = [1, 0, -1]. */

    /* Compute digital center frequency for gain normalization */
    double fc = (low_freq + high_freq) / 2.0;
    double wc = 2.0 * M_PI * fc / sample_rate;

    /* ----------------------------------------------------------------
     * Step 5: Factor into cascaded biquad (SOS) sections
     * Pair conjugate poles, pair conjugate zeros, form biquad sections.
     * For bandpass: we have N zeros at z=+1 and N at z=-1, and 2N poles
     * in conjugate pairs. We form N biquad sections.
     * ---------------------------------------------------------------- */

    /* Instead of factoring the already-built polynomial, we build SOS
     * directly from the poles and zeros for better numerical accuracy. */

    /* Pair up poles into conjugate pairs */
    /* Digital poles come from conjugate analog pairs, so
     * bp_poles[2k] and bp_poles[2k+1] come from the same LP pole.
     * Their bilinear images are also conjugate pairs if the LP pole
     * is complex (which it is for N>1), or real pairs for N=1. */

    /* BUG-01 fix: Pair poles into TRUE conjugate pairs.
     * For bandpass, 2N digital poles must be grouped into N conjugate pairs
     * (p, conj(p)) so that each biquad has REAL coefficients.
     * Previous code paired poles [2k, 2k+1] from the same LP pole, but
     * those are NOT conjugates for order > 1.
     *
     * Strategy: mark each pole as used, then for each unused pole find
     * its conjugate among the remaining poles. */
    bool used[8] = {false};  /* max 2N = 2*4 = 8 poles */
    int pair_idx[4][2];
    int npairs = 0;

    for (int i = 0; i < 2*N; i++) {
        if (used[i]) continue;

        /* Find conjugate partner: look for pole closest to conj(dz_poles[i]) */
        cdouble_t target = conj(dz_poles[i]);
        int best_j = -1;
        double best_dist = 1e30;

        for (int j = i + 1; j < 2*N; j++) {
            if (used[j]) continue;
            double dist = cabs(dz_poles[j] - target);
            if (dist < best_dist) {
                best_dist = dist;
                best_j = j;
            }
        }

        if (best_j >= 0 && best_dist < 1e-8) {
            /* Found conjugate pair */
            pair_idx[npairs][0] = i;
            pair_idx[npairs][1] = best_j;
            used[i] = true;
            used[best_j] = true;
            npairs++;
        } else {
            /* Real pole or no match found - pair with itself */
            pair_idx[npairs][0] = i;
            pair_idx[npairs][1] = i;
            used[i] = true;
            npairs++;
        }
    }

    /* Each biquad section:
     * Numerator: one factor of (1 - z^-1)(1 + z^-1) = (1 - z^-2) = [1, 0, -1]
     * Denominator: (1 - p*z^-1)(1 - conj(p)*z^-1) where p, conj(p) are paired
     *            = [1, -2*Re(p), |p|^2]
     */
    for (int s = 0; s < npairs; s++) {
        int i0 = pair_idx[s][0];
        int i1 = pair_idx[s][1];

        /* Denominator from conjugate pair → guaranteed real coefficients */
        double a0 = 1.0;
        double a1, a2;
        if (i0 == i1) {
            /* Real pole (self-conjugate) */
            double p = creal(dz_poles[i0]);
            a1 = -2.0 * p;
            a2 = p * p;
        } else {
            /* Conjugate pair: a1 = -(p + conj(p)) = -2*Re(p), a2 = |p|^2 */
            a1 = -2.0 * creal(dz_poles[i0]);
            a2 = creal(dz_poles[i0]) * creal(dz_poles[i0])
               + cimag(dz_poles[i0]) * cimag(dz_poles[i0]);
        }

        /* Numerator: [1, 0, -1] for each section */
        double b0_sec = 1.0;
        double b1_sec = 0.0;
        double b2_sec = -1.0;

        coeffs->b[s*3 + 0] = b0_sec;
        coeffs->b[s*3 + 1] = b1_sec;
        coeffs->b[s*3 + 2] = b2_sec;
        coeffs->a[s*3 + 0] = a0;
        coeffs->a[s*3 + 1] = a1;
        coeffs->a[s*3 + 2] = a2;
    }

    /* Apply overall gain: distribute across first section */
    /* Recompute gain using SOS form evaluated at center frequency */
    cdouble_t total_resp = 1.0 + 0.0*I;
    cdouble_t ejwn1 = cos(-wc) + I * sin(-wc);    /* z^-1 */
    cdouble_t ejwn2 = cos(-2*wc) + I * sin(-2*wc); /* z^-2 */

    for (int s = 0; s < npairs; s++) {
        cdouble_t num_s = coeffs->b[s*3+0]
                        + coeffs->b[s*3+1] * ejwn1
                        + coeffs->b[s*3+2] * ejwn2;
        cdouble_t den_s = coeffs->a[s*3+0]
                        + coeffs->a[s*3+1] * ejwn1
                        + coeffs->a[s*3+2] * ejwn2;
        if (cabs(den_s) < 1e-12) {
            total_resp *= num_s / (den_s + 1e-12);
        } else {
            total_resp *= num_s / den_s;
        }
    }

    double sos_gain = 1.0 / (cabs(total_resp) + EPSILON);

    /* Apply gain to first section numerator */
    coeffs->b[0] *= sos_gain;
    coeffs->b[1] *= sos_gain;
    coeffs->b[2] *= sos_gain;

    /* Store number of biquad sections as 'order' */
    coeffs->order = npairs;

    return MUSIC_SUCCESS;
}

/**
 * @brief Apply a single biquad section using Direct Form II Transposed
 *
 * Processes real-valued samples through one second-order section.
 * Uses two state variables (w1, w2) for the transposed structure.
 *
 * @param b Numerator coefficients [b0, b1, b2]
 * @param a Denominator coefficients [1, a1, a2]
 * @param x Input samples
 * @param y Output samples
 * @param n Number of samples
 * @param w1 State variable 1 (in/out)
 * @param w2 State variable 2 (in/out)
 */
static void biquad_filter(const double* b, const double* a,
                           const double* x, double* y, int n,
                           double* w1, double* w2) {
    double s1 = *w1, s2 = *w2;
    for (int i = 0; i < n; i++) {
        double out = b[0] * x[i] + s1;
        s1 = b[1] * x[i] - a[1] * out + s2;
        s2 = b[2] * x[i] - a[2] * out;
        y[i] = out;
    }
    *w1 = s1;
    *w2 = s2;
}

/**
 * @brief Apply cascaded biquad SOS filter to a real-valued signal
 *
 * @param coeffs Filter coefficients in SOS format
 * @param x Input signal
 * @param y Output signal
 * @param n Number of samples
 */
static int apply_sos_filter(const filter_coeffs_t* coeffs,
                             const double* x, double* y, int n) {
    int num_sections = coeffs->order;

    /* Copy input to output as starting point */
    for (int i = 0; i < n; i++) y[i] = x[i];

    /* Apply each biquad section in cascade */
    /* Use a temp buffer to avoid in-place issues */
    double* temp = (double*)malloc(n * sizeof(double));
    if (!temp) return -1;  /* BUG-03 fix: return error instead of silent fail */

    for (int s = 0; s < num_sections; s++) {
        double w1 = 0.0, w2 = 0.0;
        biquad_filter(&coeffs->b[s*3], &coeffs->a[s*3],
                      y, temp, n, &w1, &w2);
        for (int i = 0; i < n; i++) y[i] = temp[i];
    }

    free(temp);
    return 0;
}

music_status_t signal_processor_apply_filter(
    const filter_coeffs_t* coeffs,
    const cdouble_t* iq_samples,
    int num_samples,
    cdouble_t* output
) {
    if (!coeffs || !iq_samples || !output || num_samples <= 0) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    if (coeffs->order <= 0 || coeffs->order > 7) {
        return MUSIC_ERROR_INVALID_CONFIG;
    }

    /* ----------------------------------------------------------------
     * Zero-phase filtering (filtfilt equivalent)
     *
     * 1. Reflect-pad signal at boundaries to reduce transients
     * 2. Forward pass: apply cascaded biquad IIR filter
     * 3. Time-reverse the result
     * 4. Backward pass: apply same filter again
     * 5. Time-reverse to get final zero-phase output
     *
     * Real and imaginary parts are filtered independently.
     * ---------------------------------------------------------------- */

    /* Reflection padding length: 3 * (number_of_biquads) to settle transients */
    int pad_len = 3 * coeffs->order;
    if (pad_len > num_samples - 1) {
        pad_len = num_samples - 1;
    }

    int padded_len = num_samples + 2 * pad_len;

    /* Allocate working buffers for real and imaginary parts */
    double* re_in  = (double*)malloc(padded_len * sizeof(double));
    double* im_in  = (double*)malloc(padded_len * sizeof(double));
    double* re_out = (double*)malloc(padded_len * sizeof(double));
    double* im_out = (double*)malloc(padded_len * sizeof(double));

    if (!re_in || !im_in || !re_out || !im_out) {
        free(re_in); free(im_in); free(re_out); free(im_out);
        return MUSIC_ERROR_MEMORY;
    }

    /* Separate real and imaginary parts with reflection padding */
    /* Front reflection: mirror the first pad_len samples about index 0 */
    double re0 = creal(iq_samples[0]);
    double im0 = cimag(iq_samples[0]);
    for (int i = 0; i < pad_len; i++) {
        re_in[i] = 2.0 * re0 - creal(iq_samples[pad_len - i]);
        im_in[i] = 2.0 * im0 - cimag(iq_samples[pad_len - i]);
    }

    /* Original signal */
    for (int i = 0; i < num_samples; i++) {
        re_in[pad_len + i] = creal(iq_samples[i]);
        im_in[pad_len + i] = cimag(iq_samples[i]);
    }

    /* Back reflection: mirror the last pad_len samples about last index */
    double reN = creal(iq_samples[num_samples - 1]);
    double imN = cimag(iq_samples[num_samples - 1]);
    for (int i = 0; i < pad_len; i++) {
        re_in[pad_len + num_samples + i] = 2.0 * reN
            - creal(iq_samples[num_samples - 2 - i]);
        im_in[pad_len + num_samples + i] = 2.0 * imN
            - cimag(iq_samples[num_samples - 2 - i]);
    }

    /* Forward pass on real part */
    apply_sos_filter(coeffs, re_in, re_out, padded_len);

    /* Reverse real part */
    for (int i = 0; i < padded_len / 2; i++) {
        double tmp = re_out[i];
        re_out[i] = re_out[padded_len - 1 - i];
        re_out[padded_len - 1 - i] = tmp;
    }

    /* Backward pass on real part */
    apply_sos_filter(coeffs, re_out, re_in, padded_len);

    /* Reverse again */
    for (int i = 0; i < padded_len / 2; i++) {
        double tmp = re_in[i];
        re_in[i] = re_in[padded_len - 1 - i];
        re_in[padded_len - 1 - i] = tmp;
    }

    /* Forward pass on imaginary part */
    apply_sos_filter(coeffs, im_in, im_out, padded_len);

    /* Reverse imaginary part */
    for (int i = 0; i < padded_len / 2; i++) {
        double tmp = im_out[i];
        im_out[i] = im_out[padded_len - 1 - i];
        im_out[padded_len - 1 - i] = tmp;
    }

    /* Backward pass on imaginary part */
    apply_sos_filter(coeffs, im_out, im_in, padded_len);

    /* Reverse again */
    for (int i = 0; i < padded_len / 2; i++) {
        double tmp = im_in[i];
        im_in[i] = im_in[padded_len - 1 - i];
        im_in[padded_len - 1 - i] = tmp;
    }

    /* Extract the non-padded portion and recombine to complex */
    for (int i = 0; i < num_samples; i++) {
        output[i] = re_in[pad_len + i] + I * im_in[pad_len + i];
    }

    free(re_in);
    free(im_in);
    free(re_out);
    free(im_out);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Peak Detection and Analysis
 * ============================================================================ */

music_status_t signal_processor_detect_peaks(
    const double* frequencies,
    const double* spectrum,
    int num_bins,
    double height,
    double threshold,
    int distance,
    peak_result_t* result
) {
    // detect_peaks using scipy.signal.find_peaks
    if (!frequencies || !spectrum || !result) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Auto-detect threshold if needed
    if (isnan(threshold)) {
        // Compute median as noise floor
        double* sorted = (double*)malloc(num_bins * sizeof(double));
        if (!sorted) return MUSIC_ERROR_MEMORY;  /* BUG-08 fix */
        memcpy(sorted, spectrum, num_bins * sizeof(double));

        // Simple insertion sort for median
        for (int i = 1; i < num_bins; i++) {
            double key = sorted[i];
            int j = i - 1;
            while (j >= 0 && sorted[j] > key) {
                sorted[j + 1] = sorted[j];
                j--;
            }
            sorted[j + 1] = key;
        }

        double noise_floor = sorted[num_bins / 2];
        threshold = noise_floor + 10.0;  // 10 dB above noise floor

        free(sorted);
    }

    if (isnan(height)) {
        height = threshold;
    }

    // Find peaks
    int num_peaks = 0;
    for (int i = 1; i < num_bins - 1 && num_peaks < MAX_PEAKS; i++) {
        // Check if local maximum
        if (spectrum[i] > spectrum[i-1] && spectrum[i] > spectrum[i+1]) {
            // Check height threshold
            if (spectrum[i] >= height) {
                // Check distance from previous peaks
                bool too_close = false;
                for (int j = 0; j < num_peaks; j++) {
                    if (abs(i - result->indices[j]) < distance) {
                        too_close = true;
                        break;
                    }
                }

                if (!too_close) {
                    result->indices[num_peaks] = i;
                    result->frequencies[num_peaks] = frequencies[i];
                    result->magnitudes[num_peaks] = spectrum[i];

                    // Compute prominence (simplified)
                    result->prominences[num_peaks] = spectrum[i] - threshold;

                    // Compute width (simplified)
                    result->widths[num_peaks] = 1.0;

                    num_peaks++;
                }
            }
        }
    }

    result->num_peaks = num_peaks;

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_downsample(
    const cdouble_t* iq_samples,
    int num_samples,
    int factor,
    cdouble_t* output,
    int* output_length
) {
    // downsample using scipy.signal.decimate
    if (!iq_samples || !output || !output_length || factor <= 0) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Simplified decimation (no anti-aliasing filter for now)
    int out_len = num_samples / factor;

    for (int i = 0; i < out_len; i++) {
        output[i] = iq_samples[i * factor];
    }

    *output_length = out_len;

    return MUSIC_SUCCESS;
}

music_status_t signal_processor_estimate_snr(
    const signal_processor_t* processor,
    const cdouble_t* iq_samples,
    int num_samples,
    double signal_band_low,
    double signal_band_high,
    double* snr_db
) {
    // estimate_snr
    if (!processor || !iq_samples || !snr_db) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    // Compute PSD
    psd_result_t* psd = psd_result_alloc(processor->fft_size);
    if (!psd) {
        return MUSIC_ERROR_MEMORY;
    }

    music_status_t status = signal_processor_compute_psd(processor, iq_samples,
                                                         num_samples, psd);
    if (status != MUSIC_SUCCESS) {
        psd_result_free(psd);
        return status;
    }

    // Find signal band
    double signal_power_sum = 0.0;
    int signal_count = 0;
    double noise_power_sum = 0.0;
    int noise_count = 0;

    for (int i = 0; i < psd->num_bins; i++) {
        double freq = psd->frequencies[i];
        double power_linear = pow(10.0, psd->psd_db[i] / 10.0);

        if (freq >= signal_band_low && freq <= signal_band_high) {
            signal_power_sum += power_linear;
            signal_count++;
        } else {
            noise_power_sum += power_linear;
            noise_count++;
        }
    }

    double signal_power = signal_power_sum / (signal_count + EPSILON);
    double noise_power = noise_power_sum / (noise_count + EPSILON);

    *snr_db = 10.0 * log10((signal_power / (noise_power + EPSILON)) + EPSILON);

    psd_result_free(psd);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Memory Management Helpers
 * ============================================================================ */

fft_result_t* fft_result_alloc(int num_bins) {
    fft_result_t* result = (fft_result_t*)malloc(sizeof(fft_result_t));
    if (!result) {
        return NULL;
    }

    result->num_bins = num_bins;
    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitude_db = (double*)malloc(num_bins * sizeof(double));

    if (!result->frequencies || !result->magnitude_db) {
        fft_result_free(result);
        return NULL;
    }

    return result;
}

void fft_result_free(fft_result_t* result) {
    if (result) {
        if (result->frequencies) free(result->frequencies);
        if (result->magnitude_db) free(result->magnitude_db);
        free(result);
    }
}

psd_result_t* psd_result_alloc(int num_bins) {
    psd_result_t* result = (psd_result_t*)malloc(sizeof(psd_result_t));
    if (!result) {
        return NULL;
    }

    result->num_bins = num_bins;
    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->psd_db = (double*)malloc(num_bins * sizeof(double));

    if (!result->frequencies || !result->psd_db) {
        psd_result_free(result);
        return NULL;
    }

    return result;
}

void psd_result_free(psd_result_t* result) {
    if (result) {
        if (result->frequencies) free(result->frequencies);
        if (result->psd_db) free(result->psd_db);
        free(result);
    }
}

spectrogram_result_t* spectrogram_result_alloc(int num_time_bins, int num_freq_bins) {
    spectrogram_result_t* result = (spectrogram_result_t*)malloc(sizeof(spectrogram_result_t));
    if (!result) {
        return NULL;
    }

    result->num_time_bins = num_time_bins;
    result->num_freq_bins = num_freq_bins;
    result->times = (double*)malloc(num_time_bins * sizeof(double));
    result->frequencies = (double*)malloc(num_freq_bins * sizeof(double));
    result->spectrogram_db = (double*)malloc(num_time_bins * num_freq_bins * sizeof(double));

    if (!result->times || !result->frequencies || !result->spectrogram_db) {
        spectrogram_result_free(result);
        return NULL;
    }

    return result;
}

void spectrogram_result_free(spectrogram_result_t* result) {
    if (result) {
        if (result->times) free(result->times);
        if (result->frequencies) free(result->frequencies);
        if (result->spectrogram_db) free(result->spectrogram_db);
        free(result);
    }
}

iq_buffer_t* iq_buffer_alloc(int length) {
    iq_buffer_t* buffer = (iq_buffer_t*)malloc(sizeof(iq_buffer_t));
    if (!buffer) {
        return NULL;
    }

    buffer->length = length;
    buffer->samples = (cdouble_t*)malloc(length * sizeof(cdouble_t));

    if (!buffer->samples) {
        free(buffer);
        return NULL;
    }

    return buffer;
}

void iq_buffer_free(iq_buffer_t* buffer) {
    if (buffer) {
        if (buffer->samples) free(buffer->samples);
        free(buffer);
    }
}
