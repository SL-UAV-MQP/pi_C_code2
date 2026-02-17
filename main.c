/**
 * @file main.c
 * @brief MQP UAV Localization - Real-Time Pipeline Application
 *
 * Complete real-time pipeline:
 *   SDR Capture → Bandpass Filter → Noise Reduction → MUSIC AOA
 *   → Protocol Detection (LTE/P25) → Cell ID → Triangulation → Position
 *
 * Usage:
 *   ./mqp_localize [options]
 *     -s <site>     Site: "wpi" or "cmrcm" (default: wpi)
 *     -f <freq_hz>  Center frequency in Hz (default: from site config)
 *     -b <band>     Band index (0-5 for WPI, default: 0 = LTE Band 13)
 *     -g <gain_db>  SDR gain in dB (0-76, default: 60)
 *     -n <count>    Number of frames to process (0 = infinite, default: 0)
 *     -l <path>     Log file base path (default: /tmp/mqp_log)
 *     -d            Debug mode (verbose output)
 *     -h            Show help
 *
 * State Machine Flow:
 *   HW_INIT → SDR_TEST → ALGO_TEST → LIVE → (loop until Ctrl+C)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <stdbool.h>
#include <getopt.h>

/* MQP modules */
#include "common_types.h"
#include "music_uca_6.h"
#include "array_geometry.h"
#include "signal_processor.h"
#include "noise_reduction.h"
#include "beamforming.h"
#include "lte_detector.h"
#include "p25_detector.h"
#include "cell_tower_db.h"
#include "aoa_triangulation.h"
#include "state_machine.h"
#include "diagnostic_logger.h"
#include "wpi_test_config.h"
#include "cfar_2d.h"

#ifdef USE_SDR
#include "sdr_pluto.h"
#endif

/* ============================================================================
 * Configuration
 * ============================================================================ */

/** Number of snapshots per MUSIC frame */
#define PIPELINE_NUM_SNAPSHOTS    4000

/** MUSIC azimuth scan resolution */
#define PIPELINE_AZ_RESOLUTION    0.5

/** Real-time target (ms) */
#define PIPELINE_REALTIME_MS      200.0

/** SDR sample rate */
#define PIPELINE_SAMPLE_RATE      20000000

/** FFT size for signal processor */
#define PIPELINE_FFT_SIZE         2048

/** Minimum AOA bearings for triangulation */
#define MIN_BEARINGS_FOR_FIX      2

/** Maximum bearing history */
#define MAX_BEARING_HISTORY       20

/** Bearing timeout (seconds) - discard bearings older than this */
#define BEARING_TIMEOUT_S         10.0

/* ============================================================================
 * Pipeline Context
 * ============================================================================ */

/** Command-line options */
typedef struct {
    test_site_t site;
    int band_index;
    uint64_t freq_hz;
    int gain_db;
    int max_frames;           /* 0 = infinite */
    char log_path[256];
    bool debug;
} cli_options_t;

/** Pre-allocated pipeline buffers (MA-2/MA-6 fix: avoid per-frame malloc) */
typedef struct {
    signal_matrix_t* signals;
    cdouble_t* filtered;      /* PIPELINE_NUM_SNAPSHOTS */
    cdouble_t* ch_data;       /* PIPELINE_NUM_SNAPSHOTS, reused per channel */
    cdouble_t* nr_output;     /* PIPELINE_NUM_SNAPSHOTS */
    cdouble_t* ch0_data;      /* PIPELINE_NUM_SNAPSHOTS, for protocol detection */
    filter_coeffs_t bp_filter;
    bool filter_valid;
    bool allocated;
} pipeline_buffers_t;

/** Pipeline runtime context */
typedef struct {
    /* Configuration */
    cli_options_t opts;

    /* State machine */
    state_machine_t sm;

    /* Site configuration */
    wpi_test_config_t wpi_config;
    cell_tower_database_t cmrcm_tower_db;   /* owned DB for CMRCM site */
    cell_tower_database_t* tower_db;        /* MA-1 fix: pointer to active DB */

    /* Pre-allocated buffers */
    pipeline_buffers_t buffers;

    /* Signal processing */
    signal_processor_t* sig_proc;
    noise_reducer_t noise_reducer;
    beamforming_processor_t beamformer;

    /* Protocol detectors */
    lte_detector_t lte_det;
    p25_detector_t p25_det;

    /* Array geometry */
    array_geometry_t geometry;

    /* Diagnostic logger */
    diag_logger_t* logger;
    diag_logger_t cmrcm_logger;  /* owned logger for CMRCM site */
    bool logger_ready;

    /* SDR */
#ifdef USE_SDR
    sdr_pluto_context_t* sdr_ctx;
#endif

    /* Bearing accumulator for triangulation */
    aoa_bearing_t bearings[MAX_BEARING_HISTORY];
    int num_bearings;

    /* ENU origin (site reference) */
    enu_origin_t enu_origin;

    /* Statistics */
    uint64_t frames_processed;
    uint64_t frames_with_detections;
    double total_pipeline_ms;

    /* Running flag */
    volatile bool running;
} pipeline_context_t;

/* Global for signal handler */
static pipeline_context_t* g_ctx = NULL;

/* ============================================================================
 * Signal Handler (Ctrl+C graceful shutdown)
 * ============================================================================ */

static void signal_handler(int sig) {
    (void)sig;
    if (g_ctx) {
        g_ctx->running = false;
    }
    /* MA-4 fix: printf is not async-signal-safe, message printed after main loop */
}

/* ============================================================================
 * Timing Helpers
 * ============================================================================ */

static double get_time_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

/* ============================================================================
 * CLI Parsing
 * ============================================================================ */

static void print_usage(const char* prog) {
    printf("MQP UAV Localization System\n");
    printf("===========================\n\n");
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  -s <site>     Site: \"wpi\" or \"cmrcm\" (default: wpi)\n");
    printf("  -b <band>     Band index (0-5 for WPI, default: 0)\n");
    printf("  -f <freq_hz>  Override center frequency (Hz)\n");
    printf("  -g <gain_db>  SDR gain in dB (0-76, default: 60)\n");
    printf("  -n <count>    Number of frames (0 = infinite, default: 0)\n");
    printf("  -l <path>     Log file base path (default: /tmp/mqp_log)\n");
    printf("  -d            Debug mode\n");
    printf("  -h            Show this help\n\n");
    printf("WPI Band Index:\n");
    printf("  0: LTE Band 13 (751 MHz, Verizon)\n");
    printf("  1: LTE Band 5  (881 MHz, AT&T/VZW)\n");
    printf("  2: LTE Band 2  (1960 MHz, T-Mobile)\n");
    printf("  3: LTE Band 4  (2132 MHz, VZW/AT&T)\n");
    printf("  4: P25 Worcester PD 1 (460.025 MHz)\n");
    printf("  5: P25 Worcester PD 2 (460.125 MHz)\n\n");
    printf("Example:\n");
    printf("  %s -s wpi -b 0 -g 60 -l /tmp/wpi_test_001\n", prog);
    printf("  %s -s wpi -b 4 -d   # P25 debug mode\n\n", prog);
}

static int parse_options(int argc, char* argv[], cli_options_t* opts) {
    /* Defaults */
    opts->site = SITE_WPI_CAMPUS;
    opts->band_index = 0;
    opts->freq_hz = 0;  /* 0 = from site config */
    opts->gain_db = 60;
    opts->max_frames = 0;
    strncpy(opts->log_path, "/tmp/mqp_log", sizeof(opts->log_path) - 1);
    opts->debug = false;

    int opt;
    while ((opt = getopt(argc, argv, "s:b:f:g:n:l:dh")) != -1) {
        switch (opt) {
            case 's':
                if (strcmp(optarg, "cmrcm") == 0 || strcmp(optarg, "CMRCM") == 0) {
                    opts->site = SITE_CMRCM;
                } else if (strcmp(optarg, "wpi") == 0 || strcmp(optarg, "WPI") == 0) {
                    opts->site = SITE_WPI_CAMPUS;
                } else {
                    fprintf(stderr, "Unknown site: %s (use 'wpi' or 'cmrcm')\n", optarg);
                    return -1;
                }
                break;
            case 'b':
                opts->band_index = atoi(optarg);
                break;
            case 'f':
                opts->freq_hz = (uint64_t)atof(optarg);
                break;
            case 'g':
                opts->gain_db = atoi(optarg);
                if (opts->gain_db < 0 || opts->gain_db > 76) {
                    fprintf(stderr, "Gain must be 0-76 dB\n");
                    return -1;
                }
                break;
            case 'n':
                opts->max_frames = atoi(optarg);
                break;
            case 'l':
                strncpy(opts->log_path, optarg, sizeof(opts->log_path) - 1);
                break;
            case 'd':
                opts->debug = true;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(0);
            default:
                print_usage(argv[0]);
                return -1;
        }
    }

    return 0;
}

/* ============================================================================
 * Pipeline Initialization
 * ============================================================================ */

static int pipeline_init(pipeline_context_t* ctx) {
    printf("\n");
    printf("========================================\n");
    printf("  MQP UAV Localization System v1.0\n");
    printf("  Urban Signals of Opportunity\n");
    printf("========================================\n\n");

    /* Initialize state machine */
    sm_init(&ctx->sm);

    /* Initialize site configuration */
    if (ctx->opts.site == SITE_WPI_CAMPUS) {
        printf("[INIT] Site: WPI Campus (Worcester, MA)\n");
        if (wpi_test_config_init(&ctx->wpi_config, ctx->opts.log_path) != 0) {
            fprintf(stderr, "[INIT] Failed to initialize WPI config\n");
            return -1;
        }
        ctx->tower_db = &ctx->wpi_config.tower_db;  /* MA-1 fix: pointer, not copy */
        ctx->logger = &ctx->wpi_config.logger;
        ctx->logger_ready = ctx->wpi_config.logger_initialized;

        /* Get center frequency from band config */
        if (ctx->opts.freq_hz == 0) {
            ctx->opts.freq_hz = (uint64_t)wpi_get_sdr_tune_freq(
                &ctx->wpi_config, ctx->opts.band_index);
            if (ctx->opts.freq_hz == 0) {
                fprintf(stderr, "[INIT] Invalid band index %d\n", ctx->opts.band_index);
                return -1;
            }
        }

        /* Set ENU origin to WPI campus center */
        enu_set_origin(&ctx->enu_origin,
                       ctx->tower_db->ref_lat,
                       ctx->tower_db->ref_lon,
                       150.0);  /* ~150m ASL */
    } else {
        printf("[INIT] Site: CMRCM Field\n");
        if (cell_tower_db_init_site(&ctx->cmrcm_tower_db, SITE_CMRCM) != 0) {
            fprintf(stderr, "[INIT] Failed to initialize CMRCM tower DB\n");
            return -1;
        }
        ctx->tower_db = &ctx->cmrcm_tower_db;  /* MA-1 fix: pointer */

        /* Initialize logger directly */
        if (diag_logger_init(&ctx->cmrcm_logger, ctx->opts.log_path, DIAG_LEVEL_INFO) == 0) {
            ctx->logger = &ctx->cmrcm_logger;
            ctx->logger_ready = true;
        }

        if (ctx->opts.freq_hz == 0) {
            ctx->opts.freq_hz = 751000000;  /* Default: LTE Band 13 */
        }

        enu_set_origin(&ctx->enu_origin,
                       ctx->tower_db->ref_lat,
                       ctx->tower_db->ref_lon,
                       0.0);
    }

    printf("[INIT] Frequency: %.3f MHz\n", ctx->opts.freq_hz / 1e6);
    printf("[INIT] Gain: %d dB\n", ctx->opts.gain_db);
    printf("[INIT] Log path: %s\n", ctx->opts.log_path);

    /* Initialize array geometry (6-element UCA) */
    double wavelength = 3e8 / (double)ctx->opts.freq_hz;
    double radius = 0.6 * wavelength;  /* 0.6λ standard */
    if (array_geometry_init(&ctx->geometry, MAX_ANTENNAS, radius, wavelength) != MUSIC_SUCCESS) {
        fprintf(stderr, "[INIT] Failed to initialize array geometry\n");
        return -1;
    }
    printf("[INIT] Array: 6-element UCA, r=%.3f m, λ=%.3f m\n", radius, wavelength);

    /* Initialize signal processor */
    ctx->sig_proc = signal_processor_init(
        (double)PIPELINE_SAMPLE_RATE,
        PIPELINE_FFT_SIZE,
        WINDOW_HAMMING,
        0.5
    );
    if (!ctx->sig_proc) {
        fprintf(stderr, "[INIT] Failed to initialize signal processor\n");
        return -1;
    }

    /* Initialize noise reducer */
    if (create_noise_reducer(&ctx->noise_reducer,
                             (double)PIPELINE_SAMPLE_RATE,
                             NR_FFT_SIZE,
                             NR_ALPHA,
                             0.05) != MUSIC_SUCCESS) {
        fprintf(stderr, "[INIT] Failed to initialize noise reducer\n");
        return -1;
    }

    /* Initialize beamformer */
    if (create_beamformer(&ctx->beamformer,
                          (double)ctx->opts.freq_hz,
                          (double)PIPELINE_SAMPLE_RATE,
                          BF_DIAGONAL_LOADING) != MUSIC_SUCCESS) {
        fprintf(stderr, "[INIT] Failed to initialize beamformer\n");
        return -1;
    }

    /* Initialize protocol detectors */
    if (lte_detector_init(&ctx->lte_det, (double)PIPELINE_SAMPLE_RATE) != 0) {
        fprintf(stderr, "[INIT] Failed to initialize LTE detector\n");
        return -1;
    }

    if (p25_detector_init(&ctx->p25_det, (double)PIPELINE_SAMPLE_RATE) != 0) {
        fprintf(stderr, "[INIT] Failed to initialize P25 detector\n");
        return -1;
    }

    /* Initialize SDR */
#ifdef USE_SDR
    printf("[INIT] Initializing 3x ADALM-PLUTO SDR...\n");
    ctx->sdr_ctx = sdr_pluto_init();
    if (!ctx->sdr_ctx) {
        fprintf(stderr, "[INIT] Failed to initialize SDR context\n");
        return -1;
    }

    sdr_status_t sdr_ret = sdr_pluto_connect(ctx->sdr_ctx);
    if (sdr_ret != SDR_SUCCESS) {
        fprintf(stderr, "[INIT] SDR connect failed: %s\n",
                sdr_pluto_get_error_string(sdr_ret));
        return -1;
    }
    printf("[INIT] SDR: %d devices connected\n", ctx->sdr_ctx->num_connected);

    /* Configure SDR */
    sdr_config_t sdr_cfg;
    sdr_pluto_get_default_config(&sdr_cfg);
    sdr_cfg.center_freq_hz = ctx->opts.freq_hz;
    sdr_cfg.sample_rate_sps = PIPELINE_SAMPLE_RATE;
    sdr_cfg.gain_mode = GAIN_MODE_MANUAL;
    sdr_cfg.gain_db = ctx->opts.gain_db;
    sdr_cfg.buffer_size = PIPELINE_NUM_SNAPSHOTS;

    sdr_ret = sdr_pluto_configure(ctx->sdr_ctx, &sdr_cfg);
    if (sdr_ret != SDR_SUCCESS) {
        fprintf(stderr, "[INIT] SDR configure failed: %s\n",
                sdr_pluto_get_error_string(sdr_ret));
        return -1;
    }
    printf("[INIT] SDR configured: %.1f MHz, %d MSPS, %d dB gain\n",
           ctx->opts.freq_hz / 1e6,
           PIPELINE_SAMPLE_RATE / 1000000,
           ctx->opts.gain_db);
#else
    printf("[INIT] SDR: Disabled (compile with USE_SDR=1 for hardware)\n");
#endif

    /* MA-2/MA-6 fix: pre-allocate pipeline buffers */
    ctx->buffers.signals = signal_matrix_alloc(MAX_ANTENNAS, PIPELINE_NUM_SNAPSHOTS);
    ctx->buffers.filtered = (cdouble_t*)malloc(PIPELINE_NUM_SNAPSHOTS * sizeof(cdouble_t));
    ctx->buffers.ch_data = (cdouble_t*)malloc(PIPELINE_NUM_SNAPSHOTS * sizeof(cdouble_t));
    ctx->buffers.nr_output = (cdouble_t*)malloc(PIPELINE_NUM_SNAPSHOTS * sizeof(cdouble_t));
    ctx->buffers.ch0_data = (cdouble_t*)malloc(PIPELINE_NUM_SNAPSHOTS * sizeof(cdouble_t));
    if (!ctx->buffers.signals || !ctx->buffers.filtered || !ctx->buffers.ch_data ||
        !ctx->buffers.nr_output || !ctx->buffers.ch0_data) {
        fprintf(stderr, "[INIT] Failed to pre-allocate pipeline buffers\n");
        return -1;
    }
    ctx->buffers.allocated = true;

    /* MA-6 fix: design bandpass filter once (constant parameters) */
    double bw_half = 5e6;
    double low_cut = 1e6;
    double high_cut = bw_half;
    if (signal_processor_design_bandpass(
            (double)PIPELINE_SAMPLE_RATE,
            low_cut, high_cut, 4, &ctx->buffers.bp_filter) == MUSIC_SUCCESS) {
        ctx->buffers.filter_valid = true;
    } else {
        fprintf(stderr, "[INIT] WARNING: Bandpass filter design failed\n");
        ctx->buffers.filter_valid = false;
    }

    /* Bearing accumulator */
    ctx->num_bearings = 0;

    /* Stats */
    ctx->frames_processed = 0;
    ctx->frames_with_detections = 0;
    ctx->total_pipeline_ms = 0;

    ctx->running = true;
    g_ctx = ctx;

    printf("[INIT] Initialization complete.\n\n");
    return 0;
}

/* ============================================================================
 * Pipeline Cleanup
 * ============================================================================ */

static void pipeline_cleanup(pipeline_context_t* ctx) {
    printf("\n[SHUTDOWN] Cleaning up...\n");

    /* Print statistics */
    printf("\n========================================\n");
    printf("  Session Statistics\n");
    printf("========================================\n");
    printf("  Frames processed:    %lu\n", (unsigned long)ctx->frames_processed);
    printf("  Frames w/detections: %lu\n", (unsigned long)ctx->frames_with_detections);
    if (ctx->frames_processed > 0) {
        printf("  Avg pipeline time:   %.1f ms\n",
               ctx->total_pipeline_ms / ctx->frames_processed);
        printf("  Avg rate:            %.1f Hz\n",
               1000.0 / (ctx->total_pipeline_ms / ctx->frames_processed));
    }
    printf("========================================\n");

    /* Print diagnostic summary */
    if (ctx->logger_ready) {
        diag_print_summary(ctx->logger);
    }

    /* Free pre-allocated buffers */
    if (ctx->buffers.allocated) {
        if (ctx->buffers.signals) signal_matrix_free(ctx->buffers.signals);
        free(ctx->buffers.filtered);
        free(ctx->buffers.ch_data);
        free(ctx->buffers.nr_output);
        free(ctx->buffers.ch0_data);
        ctx->buffers.allocated = false;
    }

    /* Free resources */
    if (ctx->sig_proc) {
        signal_processor_free(ctx->sig_proc);
    }
    noise_reducer_free(&ctx->noise_reducer);
    beamforming_processor_free(&ctx->beamformer);
    lte_detector_free(&ctx->lte_det);
    p25_detector_free(&ctx->p25_det);

#ifdef USE_SDR
    if (ctx->sdr_ctx) {
        sdr_pluto_disconnect(ctx->sdr_ctx);
        sdr_pluto_free(ctx->sdr_ctx);
    }
#endif

    /* Close logger & WPI config */
    if (ctx->opts.site == SITE_WPI_CAMPUS) {
        wpi_test_config_free(&ctx->wpi_config);
    } else if (ctx->logger_ready) {
        diag_logger_close(&ctx->cmrcm_logger);
    }

    printf("[SHUTDOWN] Done.\n\n");
}

/* ============================================================================
 * Core Pipeline: Single Frame Processing
 * ============================================================================ */

/**
 * @brief Process one frame of SDR data through the full pipeline
 *
 * Pipeline stages with timing:
 *   1. SDR Capture        (~1 ms)
 *   2. Signal Processing  (~5 ms)   - Bandpass filter + normalization
 *   3. Noise Reduction    (~20 ms)  - Spectral subtraction
 *   4. MUSIC AOA          (~60 ms)  - Covariance → eigendecomp → spectrum → peaks
 *   5. Protocol Detection (~30 ms)  - LTE PSS/SSS or P25 C4FM
 *   6. Triangulation      (~1 ms)   - If enough bearings accumulated
 *   Total target: <200 ms (5 Hz real-time)
 */
static int pipeline_process_frame(pipeline_context_t* ctx) {
    double t_start = get_time_s();
    double t_stage;

    diag_timing_entry_t timing;
    memset(&timing, 0, sizeof(timing));
    timing.timestamp_s = ctx->logger_ready ? diag_get_elapsed_s(ctx->logger) : 0;

    /* ---- Stage 1: SDR Capture ---- */
    t_stage = get_time_s();

    /* Use pre-allocated signal matrix (MA-2 fix) */
    signal_matrix_t* signals = ctx->buffers.signals;

#ifdef USE_SDR
    /* Synchronize and capture from 3x PLUTOs */
    sdr_pluto_synchronize_devices(ctx->sdr_ctx);

    sdr_capture_result_t* capture = sdr_capture_result_alloc(
        PIPELINE_NUM_SNAPSHOTS, NUM_RX_CHANNELS);
    if (!capture) {
        return -1;
    }

    sdr_status_t cap_ret = sdr_pluto_capture(
        ctx->sdr_ctx, PIPELINE_NUM_SNAPSHOTS, capture);
    if (cap_ret != SDR_SUCCESS) {
        fprintf(stderr, "[PIPE] SDR capture failed: %s\n",
                sdr_pluto_get_error_string(cap_ret));
        sdr_capture_result_free(capture);
        return -1;
    }

    /* Convert to signal matrix */
    sdr_pluto_convert_to_signal_matrix(capture, signals);
    sdr_capture_result_free(capture);
#else
    /* No SDR - generate synthetic test signal for validation */
    double freq_hz = (double)ctx->opts.freq_hz;
    double wavelength = 3e8 / freq_hz;
    double sample_rate = (double)PIPELINE_SAMPLE_RATE;

    /* Simulate signal from first tower's azimuth */
    double test_azimuth_deg = 45.0;  /* default */
    if (ctx->tower_db->num_towers > 0) {
        test_azimuth_deg = ctx->tower_db->towers[0].azimuth_from_ref;
    }
    /* Generate steering vector for test signal */
    cdouble_t steer[MAX_ANTENNAS];
    compute_steering_vector(&ctx->geometry, test_azimuth_deg, 0.0,
                            wavelength, true, steer);

    /* Generate snapshots: signal + noise */
    for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
        double t = n / sample_rate;
        cdouble_t sig = cexp(I * 2.0 * M_PI * 1000.0 * t);  /* 1 kHz tone */
        for (int m = 0; m < MAX_ANTENNAS; m++) {
            /* Signal + white noise */
            double noise_re = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            double noise_im = ((double)rand() / RAND_MAX - 0.5) * 0.1;
            signals->data[n * MAX_ANTENNAS + m] =
                steer[m] * sig + (noise_re + I * noise_im);
        }
    }
#endif

    timing.sdr_capture_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Stage 2: Signal Processing (per channel) ---- */
    t_stage = get_time_s();

    /* MA-6 fix: use pre-designed filter; MA-2 fix: use pre-allocated buffers */
    if (ctx->buffers.filter_valid) {
        for (int ch = 0; ch < MAX_ANTENNAS; ch++) {
            /* Extract channel data into pre-allocated buffer */
            for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
                ctx->buffers.ch_data[n] = signals->data[n * MAX_ANTENNAS + ch];
            }

            if (signal_processor_apply_filter(&ctx->buffers.bp_filter,
                    ctx->buffers.ch_data, PIPELINE_NUM_SNAPSHOTS,
                    ctx->buffers.filtered) == MUSIC_SUCCESS) {
                for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
                    signals->data[n * MAX_ANTENNAS + ch] = ctx->buffers.filtered[n];
                }
            }
        }
    }

    timing.signal_proc_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Stage 3: Noise Reduction ---- */
    t_stage = get_time_s();

    /* Apply spectral subtraction per channel (MA-2 fix: pre-allocated buffers) */
    noise_reduction_metrics_t nr_metrics;

    for (int ch = 0; ch < MAX_ANTENNAS; ch++) {
        for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
            ctx->buffers.ch_data[n] = signals->data[n * MAX_ANTENNAS + ch];
        }

        if (noise_reducer_process_single(&ctx->noise_reducer,
                ctx->buffers.ch_data, PIPELINE_NUM_SNAPSHOTS,
                NR_SPECTRAL_SUBTRACTION,
                ctx->buffers.nr_output, &nr_metrics) == MUSIC_SUCCESS) {
            for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
                signals->data[n * MAX_ANTENNAS + ch] = ctx->buffers.nr_output[n];
            }
        }
    }

    timing.noise_reduce_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Stage 4: MUSIC AOA Estimation ---- */
    t_stage = get_time_s();

    /* Estimate number of sources */
    int num_sources = 1;  /* Conservative default */

    /* Allocate MUSIC structures */
    int num_az_bins = (int)(360.0 / PIPELINE_AZ_RESOLUTION);
    music_spectrum_t* spectrum = music_spectrum_alloc(num_az_bins, MAX_ANTENNAS);
    detected_sources_t detected;
    memset(&detected, 0, sizeof(detected));

    bool aoa_success = false;
    if (spectrum) {
        if (music_estimate_aoa(signals, &ctx->geometry,
                               num_sources, true,
                               spectrum, &detected) == MUSIC_SUCCESS) {
            aoa_success = true;
        }
    }

    timing.music_aoa_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Stage 5: Protocol Detection ---- */
    t_stage = get_time_s();

    /* Use channel 0 for protocol detection (MA-2 fix: pre-allocated buffer) */
    int detected_cell_id = -1;
    char detected_type[16] = "UNKNOWN";

    for (int n = 0; n < PIPELINE_NUM_SNAPSHOTS; n++) {
        ctx->buffers.ch0_data[n] = signals->data[n * MAX_ANTENNAS + 0];
    }

    /* Try LTE detection first */
    lte_cell_detection_t lte_cells[LTE_MAX_DETECTIONS];
    int n_lte = lte_detect_signal(&ctx->lte_det, ctx->buffers.ch0_data,
                                   PIPELINE_NUM_SNAPSHOTS,
                                   lte_cells, LTE_MAX_DETECTIONS);
    if (n_lte > 0) {
        detected_cell_id = lte_cells[0].cell_id;
        strncpy(detected_type, "LTE", sizeof(detected_type) - 1);
    }

    /* Try P25 detection */
    if (n_lte == 0) {
        p25_detection_t p25_dets[P25_MAX_DETECTIONS];
        int n_p25 = p25_detect_signal(&ctx->p25_det, ctx->buffers.ch0_data,
                                       PIPELINE_NUM_SNAPSHOTS,
                                       p25_dets, P25_MAX_DETECTIONS);
        if (n_p25 > 0) {
            strncpy(detected_type, "P25", sizeof(detected_type) - 1);
        }
    }

    timing.protocol_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Stage 6: Triangulation ---- */
    t_stage = get_time_s();

    aoa_result_t tri_result;
    memset(&tri_result, 0, sizeof(tri_result));
    wgs84_coord_t est_position;
    memset(&est_position, 0, sizeof(est_position));
    bool position_valid = false;

    if (aoa_success && detected.num_sources > 0) {
        double now = get_time_s();

        /* Add new bearing to accumulator */
        for (int s = 0; s < detected.num_sources && ctx->num_bearings < MAX_BEARING_HISTORY; s++) {
            /* Look up tower by Cell ID or azimuth to get observer position */
            const cell_tower_t* tower = NULL;

            if (detected_cell_id >= 0) {
                tower = cell_tower_db_get_by_cell_id(ctx->tower_db, detected_cell_id);
            }
            if (!tower) {
                tower = cell_tower_db_get_by_azimuth(ctx->tower_db,
                    detected.sources[s].azimuth_deg, 15.0);
            }

            if (tower) {
                /* Tower position in ENU */
                wgs84_coord_t tower_wgs84 = {tower->latitude, tower->longitude, 0};
                double t_east, t_north, t_up;
                wgs84_to_enu(&ctx->enu_origin, &tower_wgs84,
                             &t_east, &t_north, &t_up);

                /* Bearing from tower to us = reverse of our AOA */
                double reverse_az = fmod(detected.sources[s].azimuth_deg + 180.0, 360.0);

                aoa_bearing_t* b = &ctx->bearings[ctx->num_bearings];
                b->observer_x = t_east;
                b->observer_y = t_north;
                b->azimuth_deg = reverse_az;
                b->confidence = detected.sources[s].confidence;
                b->timestamp = now;
                ctx->num_bearings++;
            }
        }

        /* Remove stale bearings */
        int write_idx = 0;
        for (int i = 0; i < ctx->num_bearings; i++) {
            if (now - ctx->bearings[i].timestamp < BEARING_TIMEOUT_S) {
                if (write_idx != i) {
                    ctx->bearings[write_idx] = ctx->bearings[i];
                }
                write_idx++;
            }
        }
        ctx->num_bearings = write_idx;

        /* Triangulate if enough bearings */
        if (ctx->num_bearings >= MIN_BEARINGS_FOR_FIX) {
            if (aoa_triangulate(ctx->bearings, ctx->num_bearings,
                                &tri_result) == 0 && tri_result.success) {
                /* Convert ENU result to WGS84 */
                enu_to_wgs84(&ctx->enu_origin,
                             tri_result.est_x, tri_result.est_y, 0.0,
                             &est_position);
                position_valid = true;
            }
        }
    }

    timing.triangulation_ms = (get_time_s() - t_stage) * 1000.0;

    /* ---- Timing Summary ---- */
    double t_total_ms = (get_time_s() - t_start) * 1000.0;
    timing.total_pipeline_ms = t_total_ms;
    timing.met_realtime = (t_total_ms <= PIPELINE_REALTIME_MS);

    ctx->frames_processed++;
    ctx->total_pipeline_ms += t_total_ms;
    if (aoa_success && detected.num_sources > 0) {
        ctx->frames_with_detections++;
    }

    /* ---- Console Output ---- */
    printf("\r[%06lu] ", (unsigned long)ctx->frames_processed);

    if (aoa_success && detected.num_sources > 0) {
        printf("AOA: ");
        for (int s = 0; s < detected.num_sources; s++) {
            printf("%.1f° ", detected.sources[s].azimuth_deg);
        }
        printf("(%s", detected_type);
        if (detected_cell_id >= 0) printf(" CID=%d", detected_cell_id);
        printf(") ");
    } else {
        printf("No signal ");
    }

    if (position_valid) {
        printf("POS: %.6f, %.6f (GDOP=%.1f, err=%.1fdeg) ",
               est_position.latitude_deg,
               est_position.longitude_deg,
               tri_result.gdop,
               tri_result.residual_deg);
    } else if (ctx->num_bearings > 0) {
        printf("[%d bearings] ", ctx->num_bearings);
    }

    printf("| %.0fms %s",
           t_total_ms,
           timing.met_realtime ? "OK" : "SLOW");

    if (ctx->opts.debug) {
        printf("\n  SDR:%.1f SIG:%.1f NR:%.1f MUSIC:%.1f PROTO:%.1f TRI:%.1f",
               timing.sdr_capture_ms, timing.signal_proc_ms,
               timing.noise_reduce_ms, timing.music_aoa_ms,
               timing.protocol_ms, timing.triangulation_ms);
    }
    printf("    \n");
    fflush(stdout);

    /* ---- Diagnostic Logging ---- */
    if (ctx->logger_ready) {
        /* Log timing */
        diag_log_pipeline_timing(ctx->logger, &timing);

        /* Log AOA */
        if (aoa_success && detected.num_sources > 0) {
            for (int s = 0; s < detected.num_sources; s++) {
                diag_aoa_entry_t aoa_entry;
                aoa_entry.timestamp_s = timing.timestamp_s;
                aoa_entry.azimuth_deg = detected.sources[s].azimuth_deg;
                aoa_entry.elevation_deg = 0.0;
                aoa_entry.confidence = detected.sources[s].confidence;
                aoa_entry.num_sources = detected.num_sources;
                aoa_entry.music_spectrum_peak = detected.sources[s].magnitude_db;

                /* Check against expected azimuth */
                const cell_tower_t* closest = cell_tower_db_get_by_azimuth(
                    ctx->tower_db, detected.sources[s].azimuth_deg, 20.0);
                if (closest) {
                    aoa_entry.expected_az = closest->azimuth_from_ref;
                    aoa_entry.az_error = detected.sources[s].azimuth_deg -
                                         closest->azimuth_from_ref;
                } else {
                    aoa_entry.expected_az = -1.0;
                    aoa_entry.az_error = 0.0;
                }

                diag_log_aoa(ctx->logger, &aoa_entry);
            }
        }

        /* Log position */
        if (position_valid) {
            diag_position_entry_t pos_entry;
            pos_entry.timestamp_s = timing.timestamp_s;
            pos_entry.est_lat = est_position.latitude_deg;
            pos_entry.est_lon = est_position.longitude_deg;
            pos_entry.est_alt = est_position.altitude_m;
            pos_entry.gdop = tri_result.gdop;
            pos_entry.error_m = tri_result.residual_deg;
            pos_entry.num_bearings = tri_result.num_bearings_used;
            diag_log_position(ctx->logger, &pos_entry);
        }

        /* Log signal */
        if (detected_cell_id >= 0 || strcmp(detected_type, "P25") == 0) {
            diag_signal_entry_t sig_entry;
            memset(&sig_entry, 0, sizeof(sig_entry));
            sig_entry.timestamp_s = timing.timestamp_s;
            sig_entry.frequency_hz = (double)ctx->opts.freq_hz;
            sig_entry.cell_id = detected_cell_id;
            strncpy(sig_entry.signal_type, detected_type, sizeof(sig_entry.signal_type) - 1);

            const cell_tower_t* sig_tower = cell_tower_db_get_by_cell_id(
                ctx->tower_db, detected_cell_id);
            if (sig_tower) {
                strncpy(sig_entry.tower_name, sig_tower->name,
                        sizeof(sig_entry.tower_name) - 1);
            }
            diag_log_signal(ctx->logger, &sig_entry);
        }
    }

    /* Free MUSIC spectrum (signals is pre-allocated, not freed here) */
    if (spectrum) music_spectrum_free(spectrum);

    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

int main(int argc, char* argv[]) {
    pipeline_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));

    /* Parse command-line options */
    if (parse_options(argc, argv, &ctx.opts) != 0) {
        return 1;
    }

    /* Install signal handler */
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    /* Initialize pipeline */
    if (pipeline_init(&ctx) != 0) {
        fprintf(stderr, "[MAIN] Pipeline initialization failed\n");
        pipeline_cleanup(&ctx);
        return 1;
    }

    /* Run state machine through init phases */
    printf("[MAIN] Running state machine startup sequence...\n\n");

    /* Step through HW_INIT → SDR_TEST → ALGO_TEST */
    while (ctx.running && !sm_is_operational(&ctx.sm)) {
        system_state_t state = sm_step(&ctx.sm);

        if (state == STATE_ERROR_DIAG) {
            fprintf(stderr, "\n[MAIN] State machine entered ERROR state.\n");
            sm_print_status(&ctx.sm);

            /* In debug mode, continue to LIVE anyway for testing */
            if (ctx.opts.debug) {
                printf("[MAIN] Debug mode: forcing LIVE state...\n");
                sm_force_state(&ctx.sm, STATE_LIVE);
            } else {
                pipeline_cleanup(&ctx);
                return 1;
            }
        }
    }

    printf("\n[MAIN] System operational. Starting real-time pipeline...\n");
    printf("[MAIN] Press Ctrl+C to stop.\n\n");

    /* Print tower database for reference */
    cell_tower_db_print_summary(ctx.tower_db);
    printf("\n");

    /* ---- Main Loop: LIVE State ---- */
    while (ctx.running) {
        /* Check frame limit */
        if (ctx.opts.max_frames > 0 &&
            (int)ctx.frames_processed >= ctx.opts.max_frames) {
            printf("\n[MAIN] Reached frame limit (%d)\n", ctx.opts.max_frames);
            break;
        }

        /* Process one frame */
        int ret = pipeline_process_frame(&ctx);
        if (ret != 0) {
            /* Log error but continue */
            ctx.sm.live_errors++;
            if (ctx.logger_ready) {
                diag_log_error(ctx.logger, DIAG_LEVEL_WARN,
                               "PIPELINE", "Frame processing failed");
            }
        }

        /* Update state machine stats */
        ctx.sm.live_frames_processed = ctx.frames_processed;
        if (ctx.frames_processed > 0) {
            ctx.sm.live_avg_latency_ms =
                ctx.total_pipeline_ms / ctx.frames_processed;
        }
    }

    /* Shutdown (print message deferred from signal handler for async-signal safety) */
    printf("\n[MAIN] Shutting down...\n");
    pipeline_cleanup(&ctx);

    return 0;
}
