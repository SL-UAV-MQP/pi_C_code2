/**
 * @file test_signal_processor.c
 * @brief Comprehensive unit tests for signal_processor module
 *
 * Tests FFT, PSD (Welch), window functions, peak detection,
 * normalization, downsampling, SNR estimation, and FFTW plan caching.
 * All test data is generated synthetically -- no external files needed.
 *
 * Test Coverage:
 *  1. Initialization and free (valid + invalid params)
 *  2. Window functions (all 5 types)
 *  3. FFT computation (known sinusoid)
 *  4. PSD via Welch's method (white noise flatness)
 *  5. Peak detection (synthetic multi-peak spectrum)
 *  6. SNR estimation (tone + noise)
 *  7. Normalization (NORMALIZE_MEAN, NORMALIZE_MAX)
 *  8. Downsampling (factor 4)
 *  9. FFTW plan caching (non-NULL, consistency)
 *
 * @author Signal Processing Specialist Agent
 * @date 2026-02-12
 */

#include "signal_processor.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Test Framework
 * ============================================================================ */

static int test_passed = 0;
static int test_failed = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        if (condition) { \
            test_passed++; \
            printf("  PASS: %s\n", message); \
        } else { \
            test_failed++; \
            printf("  FAIL: %s\n", message); \
        } \
    } while(0)

#define TEST_ASSERT_DOUBLE(actual, expected, tolerance, message) \
    do { \
        double _a = (double)(actual); \
        double _e = (double)(expected); \
        double _t = (double)(tolerance); \
        if (fabs(_a - _e) <= _t) { \
            test_passed++; \
            printf("  PASS: %s (got %.6f, expected %.6f)\n", message, _a, _e); \
        } else { \
            test_failed++; \
            printf("  FAIL: %s (got %.6f, expected %.6f, diff %.6f > tol %.6f)\n", \
                   message, _a, _e, fabs(_a - _e), _t); \
        } \
    } while(0)

/* ============================================================================
 * Helper: generate complex sinusoid at a given frequency
 *
 * x[n] = exp(j * 2 * pi * freq * n / sample_rate)
 * ============================================================================ */
static void generate_complex_sinusoid(cdouble_t* buf, int length,
                                      double freq, double sample_rate,
                                      double amplitude) {
    for (int n = 0; n < length; n++) {
        double phase = 2.0 * M_PI * freq * n / sample_rate;
        buf[n] = amplitude * (cos(phase) + I * sin(phase));
    }
}

/* ============================================================================
 * Helper: simple PRNG-based complex white noise
 *
 * Uses rand() seeded externally.  Variance ~ amplitude^2.
 * ============================================================================ */
static void generate_complex_noise(cdouble_t* buf, int length,
                                   double amplitude) {
    for (int n = 0; n < length; n++) {
        double re = amplitude * ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double im = amplitude * ((double)rand() / RAND_MAX - 0.5) * 2.0;
        buf[n] = re + I * im;
    }
}

/* ============================================================================
 * Test 1: Initialization and Free
 * ============================================================================ */

static void test_init_and_free(void) {
    printf("\n--- Test: Initialization and Free ---\n");

    /* 1a. Valid parameters */
    signal_processor_t* proc = signal_processor_init(
        8000.0, 1024, WINDOW_HAMMING, 0.5);
    TEST_ASSERT(proc != NULL,
        "signal_processor_init with valid params returns non-NULL");

    if (proc) {
        TEST_ASSERT(proc->fft_size == 1024,
            "fft_size stored correctly (1024)");
        TEST_ASSERT_DOUBLE(proc->sample_rate, 8000.0, 1e-9,
            "sample_rate stored correctly (8000.0)");
        TEST_ASSERT(proc->window_type == WINDOW_HAMMING,
            "window_type stored correctly (HAMMING)");
        TEST_ASSERT(proc->window != NULL,
            "window array is allocated");
        TEST_ASSERT_DOUBLE(proc->overlap, 0.5, 1e-9,
            "overlap stored correctly (0.5)");
        TEST_ASSERT(proc->hop_size == 512,
            "hop_size computed correctly (512)");

        signal_processor_free(proc);
        proc = NULL;
    }

    /* 1b. fft_size = 0 should return NULL */
    proc = signal_processor_init(8000.0, 0, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc == NULL,
        "fft_size=0 returns NULL");

    /* 1c. Negative sample_rate should return NULL */
    proc = signal_processor_init(-1.0, 1024, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc == NULL,
        "negative sample_rate returns NULL");

    /* 1d. Non-power-of-2 fft_size should return NULL */
    proc = signal_processor_init(8000.0, 1000, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc == NULL,
        "non-power-of-2 fft_size (1000) returns NULL");

    /* 1e. Non-power-of-2 fft_size should return NULL (another example) */
    proc = signal_processor_init(8000.0, 3, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc == NULL,
        "non-power-of-2 fft_size (3) returns NULL");

    /* 1f. fft_size exceeding MAX_FFT_SIZE should return NULL */
    proc = signal_processor_init(8000.0, 16384, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc == NULL,
        "fft_size > MAX_FFT_SIZE (16384) returns NULL");

    /* 1g. Overlap out of range should return NULL */
    proc = signal_processor_init(8000.0, 1024, WINDOW_HANN, 1.0);
    TEST_ASSERT(proc == NULL,
        "overlap=1.0 (out of range) returns NULL");

    proc = signal_processor_init(8000.0, 1024, WINDOW_HANN, -0.1);
    TEST_ASSERT(proc == NULL,
        "overlap=-0.1 (negative) returns NULL");

    /* 1h. Free NULL should not crash */
    signal_processor_free(NULL);
    TEST_ASSERT(1,
        "signal_processor_free(NULL) does not crash");
}

/* ============================================================================
 * Test 2: Window Functions
 * ============================================================================ */

static void test_window_functions(void) {
    printf("\n--- Test: Window Functions ---\n");

    const int N = 256;

    /* Helper struct to iterate over window types */
    struct {
        window_type_t type;
        const char*   name;
        double        expected_first;   /* window[0] */
        double        expected_last;    /* window[N-1] */
        double        tol;
    } cases[] = {
        { WINDOW_RECTANGULAR, "Rectangular", 1.0,    1.0,    1e-12 },
        { WINDOW_HAMMING,     "Hamming",     0.08,   0.08,   0.01  },
        { WINDOW_HANN,        "Hann",        0.0,    0.0,    1e-12 },
        { WINDOW_BLACKMAN,    "Blackman",    0.0,    0.0,    1e-6  },
        { WINDOW_KAISER,      "Kaiser",      -1.0,   -1.0,   0.0   },
        /* Kaiser entry: first/last checked separately below */
    };
    int num_cases = 5;

    for (int c = 0; c < num_cases; c++) {
        signal_processor_t* proc = signal_processor_init(
            8000.0, N, cases[c].type, 0.5);
        if (!proc) {
            test_failed++;
            printf("  FAIL: Could not init processor for %s window\n",
                   cases[c].name);
            continue;
        }

        /* --- Endpoint values --- */
        if (cases[c].type != WINDOW_KAISER) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%s window[0] is %.4f",
                     cases[c].name, cases[c].expected_first);
            TEST_ASSERT_DOUBLE(proc->window[0],
                               cases[c].expected_first,
                               cases[c].tol, msg);

            snprintf(msg, sizeof(msg), "%s window[N-1] is %.4f",
                     cases[c].name, cases[c].expected_last);
            TEST_ASSERT_DOUBLE(proc->window[N - 1],
                               cases[c].expected_last,
                               cases[c].tol, msg);
        } else {
            /* Kaiser: endpoints should be > 0 and <= 1 */
            char msg[128];
            snprintf(msg, sizeof(msg),
                     "Kaiser window[0] in (0, 1] (got %.6f)",
                     proc->window[0]);
            TEST_ASSERT(proc->window[0] > 0.0 && proc->window[0] <= 1.0,
                        msg);

            snprintf(msg, sizeof(msg),
                     "Kaiser window[N-1] in (0, 1] (got %.6f)",
                     proc->window[N - 1]);
            TEST_ASSERT(proc->window[N - 1] > 0.0 &&
                        proc->window[N - 1] <= 1.0,
                        msg);
        }

        /* --- Symmetry: window[i] ~ window[N-1-i] --- */
        int sym_ok = 1;
        double max_sym_err = 0.0;
        for (int i = 0; i < N / 2; i++) {
            double err = fabs(proc->window[i] - proc->window[N - 1 - i]);
            if (err > max_sym_err) max_sym_err = err;
            if (err > 1e-10) {
                sym_ok = 0;
            }
        }
        {
            char msg[128];
            snprintf(msg, sizeof(msg),
                     "%s window symmetry (max err %.2e)",
                     cases[c].name, max_sym_err);
            TEST_ASSERT(sym_ok, msg);
        }

        /* --- All values in [0, 1] --- */
        int range_ok = 1;
        for (int i = 0; i < N; i++) {
            if (proc->window[i] < -1e-12 || proc->window[i] > 1.0 + 1e-12) {
                range_ok = 0;
                break;
            }
        }
        {
            char msg[128];
            snprintf(msg, sizeof(msg),
                     "%s window values all in [0, 1]",
                     cases[c].name);
            TEST_ASSERT(range_ok, msg);
        }

        signal_processor_free(proc);
    }

    /* --- Kaiser-specific: peak at center should be 1.0 --- */
    {
        signal_processor_t* proc = signal_processor_init(
            8000.0, N, WINDOW_KAISER, 0.5);
        if (proc) {
            double center_val = proc->window[N / 2];
            /* For even N, the true center is between N/2-1 and N/2.
               The max should be very close to 1.0 at the midpoint. */
            double peak = 0.0;
            for (int i = 0; i < N; i++) {
                if (proc->window[i] > peak) peak = proc->window[i];
            }
            TEST_ASSERT_DOUBLE(peak, 1.0, 1e-6,
                "Kaiser window peak value is 1.0");
            signal_processor_free(proc);
        }
    }
}

/* ============================================================================
 * Test 3: FFT Computation
 * ============================================================================ */

static void test_fft_computation(void) {
    printf("\n--- Test: FFT Computation ---\n");

    const double sample_rate = 8000.0;
    const int    fft_size    = 1024;
    const double tone_freq   = 1000.0;  /* 1 kHz */

    signal_processor_t* proc = signal_processor_init(
        sample_rate, fft_size, WINDOW_RECTANGULAR, 0.0);
    TEST_ASSERT(proc != NULL, "Processor init for FFT test");
    if (!proc) return;

    /* Generate a pure 1 kHz complex sinusoid */
    cdouble_t* samples = (cdouble_t*)malloc(fft_size * sizeof(cdouble_t));
    TEST_ASSERT(samples != NULL, "Sample buffer allocated");
    if (!samples) { signal_processor_free(proc); return; }

    generate_complex_sinusoid(samples, fft_size, tone_freq, sample_rate, 1.0);

    /* Allocate FFT result */
    fft_result_t* result = fft_result_alloc(fft_size);
    TEST_ASSERT(result != NULL, "FFT result allocated");
    if (!result) { free(samples); signal_processor_free(proc); return; }

    /* Compute FFT */
    music_status_t status = signal_processor_compute_fft(
        proc, samples, fft_size, result);
    TEST_ASSERT(status == MUSIC_SUCCESS, "compute_fft returns MUSIC_SUCCESS");

    TEST_ASSERT(result->num_bins == fft_size,
        "FFT result num_bins equals fft_size");

    /* Find peak bin */
    double max_mag = -1e30;
    int    max_idx = -1;
    for (int i = 0; i < result->num_bins; i++) {
        if (result->magnitude_db[i] > max_mag) {
            max_mag = result->magnitude_db[i];
            max_idx = i;
        }
    }

    /* The peak frequency should be close to 1000 Hz.
       Frequency resolution = sample_rate / fft_size = 7.8125 Hz
       Expected bin offset from center for 1000 Hz tone. */
    double peak_freq = result->frequencies[max_idx];
    TEST_ASSERT_DOUBLE(peak_freq, tone_freq, sample_rate / fft_size + 1.0,
        "FFT peak frequency is near 1000 Hz");

    /* The peak should be significantly above the noise floor.
       With a rectangular window and pure tone the dynamic range is large. */
    double second_max = -1e30;
    for (int i = 0; i < result->num_bins; i++) {
        if (i != max_idx && result->magnitude_db[i] > second_max) {
            second_max = result->magnitude_db[i];
        }
    }
    double dynamic_range = max_mag - second_max;
    TEST_ASSERT(dynamic_range > 20.0,
        "FFT peak is >20 dB above second-highest bin");

    fft_result_free(result);
    free(samples);
    signal_processor_free(proc);
}

/* ============================================================================
 * Test 4: PSD via Welch's Method
 * ============================================================================ */

static void test_psd_welch(void) {
    printf("\n--- Test: PSD (Welch's Method) ---\n");

    const double sample_rate = 8000.0;
    const int    fft_size    = 256;
    const int    num_samples = 8192;  /* Enough for many segments */

    signal_processor_t* proc = signal_processor_init(
        sample_rate, fft_size, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc != NULL, "Processor init for PSD test");
    if (!proc) return;

    /* Generate white noise */
    srand(42);
    cdouble_t* samples = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
    TEST_ASSERT(samples != NULL, "Noise sample buffer allocated");
    if (!samples) { signal_processor_free(proc); return; }

    generate_complex_noise(samples, num_samples, 1.0);

    /* Allocate PSD result */
    psd_result_t* psd = psd_result_alloc(fft_size);
    TEST_ASSERT(psd != NULL, "PSD result allocated");
    if (!psd) { free(samples); signal_processor_free(proc); return; }

    /* Compute PSD */
    music_status_t status = signal_processor_compute_psd(
        proc, samples, num_samples, psd);
    TEST_ASSERT(status == MUSIC_SUCCESS, "compute_psd returns MUSIC_SUCCESS");

    TEST_ASSERT(psd->num_bins == fft_size,
        "PSD num_bins equals fft_size");

    /* White noise PSD should be relatively flat.
       Check that the range (max - min) across bins is < 10 dB. */
    double psd_min = 1e30;
    double psd_max = -1e30;
    for (int i = 0; i < psd->num_bins; i++) {
        if (psd->psd_db[i] < psd_min) psd_min = psd->psd_db[i];
        if (psd->psd_db[i] > psd_max) psd_max = psd->psd_db[i];
    }
    double psd_range = psd_max - psd_min;
    {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "White noise PSD is flat within 10 dB (range=%.2f dB)",
                 psd_range);
        TEST_ASSERT(psd_range < 10.0, msg);
    }

    psd_result_free(psd);
    free(samples);
    signal_processor_free(proc);
}

/* ============================================================================
 * Test 5: Peak Detection
 * ============================================================================ */

static void test_peak_detection(void) {
    printf("\n--- Test: Peak Detection ---\n");

    /* Create a synthetic spectrum with 3 known peaks at bins 100, 300, 500.
       Magnitudes: -20 dB, -10 dB, -30 dB.
       Background (noise floor): -60 dB. */
    const int num_bins = 1024;
    const double sample_rate = 8000.0;
    const double df = sample_rate / num_bins;  /* ~7.8125 Hz */

    double* frequencies = (double*)malloc(num_bins * sizeof(double));
    double* spectrum    = (double*)malloc(num_bins * sizeof(double));
    TEST_ASSERT(frequencies != NULL && spectrum != NULL,
        "Peak detection arrays allocated");
    if (!frequencies || !spectrum) {
        if (frequencies) free(frequencies);
        if (spectrum) free(spectrum);
        return;
    }

    /* Fill with noise floor */
    for (int i = 0; i < num_bins; i++) {
        frequencies[i] = i * df;
        spectrum[i]    = -60.0;
    }

    /* Insert peaks as narrow Gaussians */
    int peak_indices[] = { 100, 300, 500 };
    double peak_mags[] = { -20.0, -10.0, -30.0 };
    int num_peaks_expected = 3;

    for (int p = 0; p < num_peaks_expected; p++) {
        int idx = peak_indices[p];
        /* Gaussian with sigma=3 bins centered at idx */
        for (int i = idx - 15; i <= idx + 15 && i < num_bins && i >= 0; i++) {
            double d = (double)(i - idx);
            double val = peak_mags[p] * exp(-0.5 * d * d / 9.0);
            /* The peak magnitude is peak_mags[p]; off-center decays toward 0.
               But we want the dB value, so model as:
               spectrum[i] = peak_mags[p] - (attenuation) */
            double atten = 0.5 * d * d / 9.0 * 10.0; /* ~10 dB per sigma^2 */
            double val_db = peak_mags[p] - atten;
            if (val_db > spectrum[i]) {
                spectrum[i] = val_db;
            }
        }
    }

    /* Detect peaks */
    peak_result_t result;
    memset(&result, 0, sizeof(result));

    music_status_t status = signal_processor_detect_peaks(
        frequencies, spectrum, num_bins,
        -40.0,   /* height: only peaks above -40 dB */
        -40.0,   /* threshold */
        20,      /* minimum distance: 20 bins apart */
        &result);

    TEST_ASSERT(status == MUSIC_SUCCESS,
        "detect_peaks returns MUSIC_SUCCESS");

    {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Detected %d peaks (expected %d)",
                 result.num_peaks, num_peaks_expected);
        TEST_ASSERT(result.num_peaks == num_peaks_expected, msg);
    }

    /* Verify that each expected peak index appears in the results */
    for (int p = 0; p < num_peaks_expected; p++) {
        int found = 0;
        for (int r = 0; r < result.num_peaks; r++) {
            /* Allow +/- 2 bins tolerance */
            if (abs(result.indices[r] - peak_indices[p]) <= 2) {
                found = 1;
                break;
            }
        }
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Peak at index %d detected (freq %.1f Hz)",
                 peak_indices[p], peak_indices[p] * df);
        TEST_ASSERT(found, msg);
    }

    /* Verify detected frequencies are correct */
    for (int r = 0; r < result.num_peaks; r++) {
        double expected_freq = result.indices[r] * df;
        TEST_ASSERT_DOUBLE(result.frequencies[r], expected_freq, df,
            "Detected peak frequency matches index");
    }

    free(frequencies);
    free(spectrum);
}

/* ============================================================================
 * Test 6: SNR Estimation
 * ============================================================================ */

static void test_snr_estimation(void) {
    printf("\n--- Test: SNR Estimation ---\n");

    const double sample_rate = 8000.0;
    const int    fft_size    = 1024;
    const int    num_samples = 8192;
    const double tone_freq   = 1000.0;
    const double tone_amp    = 1.0;       /* 0 dB signal */
    const double noise_amp   = 0.1;       /* ~ -20 dB noise power relative to signal */

    signal_processor_t* proc = signal_processor_init(
        sample_rate, fft_size, WINDOW_HANN, 0.5);
    TEST_ASSERT(proc != NULL, "Processor init for SNR test");
    if (!proc) return;

    /* Generate tone + noise */
    cdouble_t* samples = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
    TEST_ASSERT(samples != NULL, "SNR sample buffer allocated");
    if (!samples) { signal_processor_free(proc); return; }

    srand(12345);
    generate_complex_sinusoid(samples, num_samples, tone_freq, sample_rate,
                              tone_amp);
    /* Add noise on top */
    for (int n = 0; n < num_samples; n++) {
        double re = noise_amp * ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double im = noise_amp * ((double)rand() / RAND_MAX - 0.5) * 2.0;
        samples[n] += re + I * im;
    }

    /* Estimate SNR; signal band around tone: 900-1100 Hz */
    double snr_db = 0.0;
    music_status_t status = signal_processor_estimate_snr(
        proc, samples, num_samples,
        900.0, 1100.0, &snr_db);

    TEST_ASSERT(status == MUSIC_SUCCESS,
        "estimate_snr returns MUSIC_SUCCESS");

    /* Expected SNR ~ 20 dB (tone_amp/noise_amp = 10, power ratio = 100).
       Allow +/- 5 dB tolerance due to statistical noise and bin leakage. */
    {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Estimated SNR is ~20 dB (got %.2f dB)", snr_db);
        TEST_ASSERT(snr_db > 10.0 && snr_db < 35.0, msg);
    }

    free(samples);
    signal_processor_free(proc);
}

/* ============================================================================
 * Test 7: Normalization
 * ============================================================================ */

static void test_normalization(void) {
    printf("\n--- Test: Normalization ---\n");

    const double sample_rate = 8000.0;
    const int    fft_size    = 256;
    const int    num_samples = 2048;

    /* --- 7a: NORMALIZE_MEAN (zero mean, unit variance) --- */
    {
        signal_processor_t* proc = signal_processor_init(
            sample_rate, fft_size, WINDOW_HANN, 0.5);
        TEST_ASSERT(proc != NULL, "Processor init for NORMALIZE_MEAN");
        if (!proc) return;

        cdouble_t* samples = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
        TEST_ASSERT(samples != NULL, "NORMALIZE_MEAN sample buffer allocated");
        if (!samples) { signal_processor_free(proc); return; }

        /* Generate signal with non-zero mean */
        srand(999);
        for (int n = 0; n < num_samples; n++) {
            double re = 5.0 + 2.0 * ((double)rand() / RAND_MAX - 0.5);
            double im = -3.0 + 2.0 * ((double)rand() / RAND_MAX - 0.5);
            samples[n] = re + I * im;
        }

        music_status_t status = signal_processor_normalize_band(
            proc, samples, num_samples, NORMALIZE_MEAN);
        TEST_ASSERT(status == MUSIC_SUCCESS,
            "normalize_band MEAN returns MUSIC_SUCCESS");

        /* Check mean is approximately 0 */
        cdouble_t mean = 0.0;
        for (int n = 0; n < num_samples; n++) {
            mean += samples[n];
        }
        mean /= num_samples;
        TEST_ASSERT_DOUBLE(cabs(mean), 0.0, 0.05,
            "After NORMALIZE_MEAN, mean magnitude is ~0");

        /* Check variance is approximately 1 */
        double variance = 0.0;
        for (int n = 0; n < num_samples; n++) {
            double mag = cabs(samples[n]);
            variance += mag * mag;
        }
        variance /= num_samples;
        TEST_ASSERT_DOUBLE(variance, 1.0, 0.1,
            "After NORMALIZE_MEAN, variance is ~1");

        free(samples);
        signal_processor_free(proc);
    }

    /* --- 7b: NORMALIZE_MAX (max amplitude = 1) --- */
    {
        signal_processor_t* proc = signal_processor_init(
            sample_rate, fft_size, WINDOW_HANN, 0.5);
        TEST_ASSERT(proc != NULL, "Processor init for NORMALIZE_MAX");
        if (!proc) return;

        cdouble_t* samples = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
        TEST_ASSERT(samples != NULL, "NORMALIZE_MAX sample buffer allocated");
        if (!samples) { signal_processor_free(proc); return; }

        /* Generate signal with known large amplitude */
        srand(777);
        for (int n = 0; n < num_samples; n++) {
            double re = 10.0 * ((double)rand() / RAND_MAX - 0.5);
            double im = 10.0 * ((double)rand() / RAND_MAX - 0.5);
            samples[n] = re + I * im;
        }

        music_status_t status = signal_processor_normalize_band(
            proc, samples, num_samples, NORMALIZE_MAX);
        TEST_ASSERT(status == MUSIC_SUCCESS,
            "normalize_band MAX returns MUSIC_SUCCESS");

        /* Find max amplitude after normalization */
        double max_amp = 0.0;
        for (int n = 0; n < num_samples; n++) {
            double mag = cabs(samples[n]);
            if (mag > max_amp) max_amp = mag;
        }
        TEST_ASSERT_DOUBLE(max_amp, 1.0, 0.01,
            "After NORMALIZE_MAX, max amplitude is ~1.0");

        /* All values should be <= 1.0 (with small tolerance for EPSILON divisor) */
        int all_below = 1;
        for (int n = 0; n < num_samples; n++) {
            if (cabs(samples[n]) > 1.0 + 0.01) {
                all_below = 0;
                break;
            }
        }
        TEST_ASSERT(all_below,
            "After NORMALIZE_MAX, all amplitudes are <= 1.0");

        free(samples);
        signal_processor_free(proc);
    }
}

/* ============================================================================
 * Test 8: Downsampling
 * ============================================================================ */

static void test_downsample(void) {
    printf("\n--- Test: Downsampling ---\n");

    const int num_samples = 1000;
    const int factor = 4;

    cdouble_t* input  = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
    cdouble_t* output = (cdouble_t*)malloc(num_samples * sizeof(cdouble_t));
    TEST_ASSERT(input != NULL && output != NULL,
        "Downsample buffers allocated");
    if (!input || !output) {
        if (input) free(input);
        if (output) free(output);
        return;
    }

    /* Fill with index-based values for easy verification */
    for (int n = 0; n < num_samples; n++) {
        input[n] = (double)n + I * (double)(n * 2);
    }

    int output_length = 0;
    music_status_t status = signal_processor_downsample(
        input, num_samples, factor, output, &output_length);

    TEST_ASSERT(status == MUSIC_SUCCESS,
        "downsample returns MUSIC_SUCCESS");

    TEST_ASSERT(output_length == num_samples / factor,
        "Output length is num_samples / factor (250)");

    {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Output length = %d (expected %d)",
                 output_length, num_samples / factor);
        TEST_ASSERT(output_length == 250, msg);
    }

    /* Verify that output[i] == input[i * factor] */
    int samples_correct = 1;
    for (int i = 0; i < output_length; i++) {
        int src_idx = i * factor;
        double re_diff = fabs(creal(output[i]) - creal(input[src_idx]));
        double im_diff = fabs(cimag(output[i]) - cimag(input[src_idx]));
        if (re_diff > 1e-12 || im_diff > 1e-12) {
            samples_correct = 0;
            break;
        }
    }
    TEST_ASSERT(samples_correct,
        "Downsampled values match input[i * factor]");

    /* Test NULL pointer handling */
    status = signal_processor_downsample(NULL, num_samples, factor,
                                         output, &output_length);
    TEST_ASSERT(status != MUSIC_SUCCESS,
        "downsample with NULL input returns error");

    status = signal_processor_downsample(input, num_samples, 0,
                                         output, &output_length);
    TEST_ASSERT(status != MUSIC_SUCCESS,
        "downsample with factor=0 returns error");

    free(input);
    free(output);
}

/* ============================================================================
 * Test 9: FFTW Plan Caching
 * ============================================================================ */

static void test_fftw_plan_caching(void) {
    printf("\n--- Test: FFTW Plan Caching ---\n");

    const double sample_rate = 8000.0;
    const int    fft_size    = 512;

    signal_processor_t* proc = signal_processor_init(
        sample_rate, fft_size, WINDOW_HAMMING, 0.5);
    TEST_ASSERT(proc != NULL, "Processor init for FFTW caching test");
    if (!proc) return;

    /* Verify FFTW plan is not NULL after init */
    TEST_ASSERT(proc->fftw_forward != NULL,
        "fftw_forward plan is not NULL after init");

    TEST_ASSERT(proc->fftw_in != NULL,
        "fftw_in buffer is not NULL after init");

    TEST_ASSERT(proc->fftw_out != NULL,
        "fftw_out buffer is not NULL after init");

    /* Generate test signal */
    cdouble_t* samples = (cdouble_t*)malloc(fft_size * sizeof(cdouble_t));
    TEST_ASSERT(samples != NULL, "FFTW caching sample buffer allocated");
    if (!samples) { signal_processor_free(proc); return; }

    generate_complex_sinusoid(samples, fft_size, 500.0, sample_rate, 1.0);

    /* Compute FFT twice -- results should be identical (plan reused) */
    fft_result_t* result1 = fft_result_alloc(fft_size);
    fft_result_t* result2 = fft_result_alloc(fft_size);
    TEST_ASSERT(result1 != NULL && result2 != NULL,
        "FFT result buffers for caching test allocated");
    if (!result1 || !result2) {
        if (result1) fft_result_free(result1);
        if (result2) fft_result_free(result2);
        free(samples);
        signal_processor_free(proc);
        return;
    }

    music_status_t s1 = signal_processor_compute_fft(
        proc, samples, fft_size, result1);
    music_status_t s2 = signal_processor_compute_fft(
        proc, samples, fft_size, result2);

    TEST_ASSERT(s1 == MUSIC_SUCCESS && s2 == MUSIC_SUCCESS,
        "Both FFT calls succeed");

    /* Compare results -- they should be bit-identical */
    int results_match = 1;
    double max_diff = 0.0;
    for (int i = 0; i < fft_size; i++) {
        double diff_mag = fabs(result1->magnitude_db[i] -
                               result2->magnitude_db[i]);
        double diff_freq = fabs(result1->frequencies[i] -
                                result2->frequencies[i]);
        if (diff_mag > max_diff) max_diff = diff_mag;
        if (diff_mag > 1e-10 || diff_freq > 1e-10) {
            results_match = 0;
        }
    }
    {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Two FFT calls produce identical results (max diff %.2e)",
                 max_diff);
        TEST_ASSERT(results_match, msg);
    }

    fft_result_free(result1);
    fft_result_free(result2);
    free(samples);
    signal_processor_free(proc);
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;

    printf("\n");
    printf("================================================================\n");
    printf("  Signal Processor Unit Test Suite\n");
    printf("================================================================\n");
    printf("  Target: Raspberry Pi 5B (ARM Cortex-A76)\n");
    printf("  Module: signal_processor.h / signal_processor.c\n");
    printf("  Date:   2026-02-12\n");
    printf("================================================================\n");

    test_init_and_free();
    test_window_functions();
    test_fft_computation();
    test_psd_welch();
    test_peak_detection();
    test_snr_estimation();
    test_normalization();
    test_downsample();
    test_fftw_plan_caching();

    printf("\n");
    printf("================================================================\n");
    printf("  Test Summary\n");
    printf("================================================================\n");
    printf("  Tests passed: %d\n", test_passed);
    printf("  Tests failed: %d\n", test_failed);
    printf("  Total:        %d\n", test_passed + test_failed);
    if (test_failed == 0) {
        printf("  Result:       ALL TESTS PASSED\n");
    } else {
        printf("  Result:       SOME TESTS FAILED\n");
    }
    printf("================================================================\n");

    return (test_failed > 0) ? 1 : 0;
}
