/**
 * @file test_wpi_integration.c
 * @brief Integration tests for WPI campus testing features
 *
 * Tests: site-configurable tower DB, Butterworth filter, LCMV weights,
 *        diagnostic logger, WPI test config, and campus self-test.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "common_types.h"
#include "cell_tower_db.h"
#include "signal_processor.h"
#include "beamforming.h"
#include "noise_reduction.h"
#include "diagnostic_logger.h"
#include "wpi_test_config.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST_ASSERT(cond, msg) do { \
    tests_run++; \
    if (cond) { tests_passed++; printf("  PASS: %s\n", msg); } \
    else { printf("  FAIL: %s\n", msg); } \
} while(0)

/* ============================================================================
 * Test 1: WPI Tower Database
 * ============================================================================ */

static void test_wpi_tower_db(void) {
    printf("\n--- Test: WPI Tower Database ---\n");

    cell_tower_database_t db;

    /* Test CMRCM init (backward compatibility) */
    int ret = cell_tower_db_init(&db);
    TEST_ASSERT(ret == 0, "CMRCM init succeeds");
    TEST_ASSERT(db.active_site == SITE_CMRCM, "Default site is CMRCM");
    TEST_ASSERT(db.num_towers == 3, "CMRCM has 3 towers");

    /* Test WPI campus init */
    ret = cell_tower_db_init_site(&db, SITE_WPI_CAMPUS);
    TEST_ASSERT(ret == 0, "WPI init succeeds");
    TEST_ASSERT(db.active_site == SITE_WPI_CAMPUS, "Active site is WPI");
    TEST_ASSERT(db.num_towers == 7, "WPI has 7 towers");
    TEST_ASSERT(strcmp(db.site_name, "WPI Campus") == 0, "Site name is WPI Campus");

    /* Test reference coordinates */
    TEST_ASSERT(fabs(db.ref_lat - 42.274580) < 0.001, "WPI lat correct");
    TEST_ASSERT(fabs(db.ref_lon - (-71.806520)) < 0.001, "WPI lon correct");

    /* Test tower lookup by name */
    const cell_tower_t* tower = cell_tower_db_get_by_name(&db, "Park_Ave_ATT");
    TEST_ASSERT(tower != NULL, "Find Park_Ave_ATT tower");
    if (tower) {
        TEST_ASSERT(strcmp(tower->carrier, "AT&T") == 0, "Park Ave carrier is AT&T");
        TEST_ASSERT(tower->num_freq_bands == 4, "Park Ave has 4 bands");
    }

    /* Test distance computation */
    double dist = cell_tower_db_calculate_distance(&db, 42.274200, -71.808800);
    TEST_ASSERT(dist < 500.0, "On-campus tower is within 500m");
    TEST_ASSERT(dist > 0.0, "Distance is positive");

    /* Test azimuth computation */
    double az = cell_tower_db_calculate_azimuth(&db, 42.278400, -71.802100);
    TEST_ASSERT(az >= 0.0 && az <= 360.0, "Azimuth in valid range");

    /* Test expected azimuths */
    double azimuths[MAX_TOWERS];
    int n = cell_tower_db_get_expected_azimuths(&db, azimuths, MAX_TOWERS);
    TEST_ASSERT(n == 7, "Got 7 expected azimuths");

    /* Test custom site init */
    ret = cell_tower_db_init_custom(&db, "Custom Test", 42.0, -71.0);
    TEST_ASSERT(ret == 0, "Custom site init succeeds");
    TEST_ASSERT(db.num_towers == 0, "Custom site starts empty");

    /* Test site getter */
    cell_tower_db_init_site(&db, SITE_WPI_CAMPUS);
    TEST_ASSERT(cell_tower_db_get_site(&db) == SITE_WPI_CAMPUS, "get_site returns WPI");

    /* Test azimuth-based lookup */
    /* Park Ave ATT is ~NE of WPI, so azimuth should be roughly 20-80 deg */
    const cell_tower_t* az_tower = cell_tower_db_get_by_azimuth(&db, az, 15.0);
    TEST_ASSERT(az_tower != NULL, "Find tower by azimuth");
}

/* ============================================================================
 * Test 2: Butterworth Bandpass Filter
 * ============================================================================ */

static void test_butterworth_filter(void) {
    printf("\n--- Test: Butterworth Bandpass Filter ---\n");

    filter_coeffs_t coeffs;
    memset(&coeffs, 0, sizeof(coeffs));

    /* Design bandpass filter: 800-900 MHz at 61.44 MHz sample rate
     * Note: These are baseband frequencies since SDR centers at carrier */
    double fs = 61.44e6;
    double f_low = 1e6;   /* 1 MHz from center */
    double f_high = 5e6;  /* 5 MHz from center */

    music_status_t status = signal_processor_design_bandpass(fs, f_low, f_high, 2, &coeffs);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Bandpass design succeeds");
    TEST_ASSERT(coeffs.order > 0, "Filter order > 0 (SOS sections)");
    TEST_ASSERT(coeffs.order <= 7, "Filter order within bounds");

    /* Verify coefficients are not identity */
    bool not_identity = false;
    for (int i = 0; i < coeffs.order * 3; i++) {
        if (fabs(coeffs.a[i]) > 1e-10 && i != 0) {
            not_identity = true;
            break;
        }
    }
    TEST_ASSERT(not_identity, "Filter coefficients are not identity");

    /* Test filter application on synthetic signal */
    int N = 4096;
    cdouble_t* input = (cdouble_t*)malloc(N * sizeof(cdouble_t));
    cdouble_t* output = (cdouble_t*)malloc(N * sizeof(cdouble_t));

    /* Generate tone at 3 MHz (within passband) + tone at 15 MHz (stopband) */
    for (int i = 0; i < N; i++) {
        double t = (double)i / fs;
        double sig_pass = cos(2.0 * M_PI * 3e6 * t);
        double sig_stop = cos(2.0 * M_PI * 15e6 * t);
        input[i] = (sig_pass + sig_stop) + 0.0 * I;
    }

    status = signal_processor_apply_filter(&coeffs, input, N, output);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Filter application succeeds");

    /* Check that passband signal is preserved and stopband is attenuated */
    double power_in_pass = 0.0, power_in_stop = 0.0;
    double power_out = 0.0;
    for (int i = N/4; i < 3*N/4; i++) {  /* Skip transients at edges */
        power_out += cabs(output[i]) * cabs(output[i]);
        power_in_pass += cos(2.0 * M_PI * 3e6 * (double)i / fs) *
                         cos(2.0 * M_PI * 3e6 * (double)i / fs);
        power_in_stop += cos(2.0 * M_PI * 15e6 * (double)i / fs) *
                         cos(2.0 * M_PI * 15e6 * (double)i / fs);
    }

    /* Output should have less power than combined input (stopband removed) */
    double total_in = power_in_pass + power_in_stop;
    TEST_ASSERT(power_out < total_in * 0.9, "Filter attenuates out-of-band signal");
    TEST_ASSERT(power_out > 0.01, "Filter preserves in-band signal");

    /* Test edge cases */
    status = signal_processor_design_bandpass(fs, f_low, f_high, 0, &coeffs);
    TEST_ASSERT(status != MUSIC_SUCCESS, "Rejects order 0");

    status = signal_processor_design_bandpass(fs, 10e6, 5e6, 2, &coeffs);
    TEST_ASSERT(status != MUSIC_SUCCESS, "Rejects low > high");

    free(input);
    free(output);
}

/* ============================================================================
 * Test 3: LCMV Beamforming (via MPDR)
 * ============================================================================ */

static void test_lcmv_beamforming(void) {
    printf("\n--- Test: LCMV Beamforming ---\n");

    /* Initialize adaptive beamformer */
    adaptive_beamformer_t bf;
    beamforming_config_t config;
    beamforming_config_default(&config);

    music_status_t status = adaptive_beamformer_init(&bf, &config);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Adaptive beamformer init");

    int M = config.num_elements;
    int N = 256;

    /* Create synthetic multichannel data: signal from 45 deg + interference from 180 deg */
    cdouble_t* signals = (cdouble_t*)calloc(M * N, sizeof(cdouble_t));

    double wavelength = 299792458.0 / config.center_frequency;
    double k = 2.0 * M_PI / wavelength;

    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            double phi_m = 2.0 * M_PI * m / M;
            double x_m = config.radius * cos(phi_m);
            double y_m = config.radius * sin(phi_m);

            /* Signal at 45 deg */
            double phase_sig = k * (x_m * cos(45.0 * M_PI / 180.0) +
                                    y_m * sin(45.0 * M_PI / 180.0));
            /* Interference at 180 deg */
            double phase_int = k * (x_m * cos(180.0 * M_PI / 180.0) +
                                    y_m * sin(180.0 * M_PI / 180.0));

            signals[m + n * M] = cexp(I * phase_sig) + 2.0 * cexp(I * phase_int) +
                                 0.1 * ((double)rand() / RAND_MAX - 0.5);
        }
    }

    /* Estimate covariance */
    cdouble_t* R = (cdouble_t*)calloc(M * M, sizeof(cdouble_t));
    status = adaptive_estimate_covariance(signals, M, N, true, R);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Covariance estimation");

    /* Compute MPDR/LCMV weights: unit gain at 45 deg, null at 180 deg */
    cdouble_t* weights = (cdouble_t*)calloc(M, sizeof(cdouble_t));
    double interference_az[] = {180.0};

    status = adaptive_compute_mpdr_weights(&bf, R, M, 45.0, interference_az, 1, weights);
    TEST_ASSERT(status == MUSIC_SUCCESS, "MPDR/LCMV weight computation");

    /* Verify weights are non-zero */
    double weight_norm = 0.0;
    for (int m = 0; m < M; m++) {
        weight_norm += cabs(weights[m]) * cabs(weights[m]);
    }
    TEST_ASSERT(weight_norm > 1e-10, "LCMV weights are non-zero");

    /* TW-1 fix: verify LCMV null-steering property */
    /* Compute beam response at desired direction (45 deg) */
    cdouble_t response_desired = 0.0;
    cdouble_t response_interf = 0.0;
    for (int m = 0; m < M; m++) {
        double phi_m = 2.0 * M_PI * m / M;
        double x_m = config.radius * cos(phi_m);
        double y_m = config.radius * sin(phi_m);

        /* Steering vector at desired direction (45 deg) */
        double phase_d = k * (x_m * cos(45.0 * M_PI / 180.0) +
                              y_m * sin(45.0 * M_PI / 180.0));
        response_desired += conj(weights[m]) * cexp(I * phase_d);

        /* Steering vector at interference direction (180 deg) */
        double phase_i = k * (x_m * cos(180.0 * M_PI / 180.0) +
                              y_m * sin(180.0 * M_PI / 180.0));
        response_interf += conj(weights[m]) * cexp(I * phase_i);
    }
    TEST_ASSERT(cabs(response_desired) > 0.3, "LCMV: desired direction gain > 0.3");
    TEST_ASSERT(cabs(response_interf) < cabs(response_desired),
                "LCMV: interference response < desired response");

    adaptive_beamformer_free(&bf);
    free(signals);
    free(R);
    free(weights);
}

/* ============================================================================
 * Test 4: Spatial LCMV Noise Canceller
 * ============================================================================ */

static void test_lcmv_spatial(void) {
    printf("\n--- Test: LCMV Spatial Noise Canceller ---\n");

    noise_reduction_config_t config;
    noise_reduction_config_default(&config);

    spatial_noise_canceller_t canceller;
    music_status_t status = spatial_canceller_init(&canceller, &config, NULL);
    TEST_ASSERT(status == MUSIC_SUCCESS, "Spatial canceller init");

    /* Compute LCMV weights: desired at 90 deg, interference at 270 deg */
    cdouble_t weights[6];
    double interference[] = {270.0};

    status = spatial_compute_null_weights(&canceller, interference, 1,
                                          1050e6, 90.0, weights);
    TEST_ASSERT(status == MUSIC_SUCCESS, "LCMV null weight computation");

    /* Verify weights are non-uniform (not placeholder anymore) */
    bool non_uniform = false;
    for (int i = 1; i < 6; i++) {
        if (cabs(weights[i] - weights[0]) > 1e-10) {
            non_uniform = true;
            break;
        }
    }
    TEST_ASSERT(non_uniform, "LCMV weights are non-uniform (not placeholder)");

    /* TW-1 fix: verify spatial null-steering property */
    double freq_hz = 1050e6;
    double wl = 299792458.0 / freq_hz;
    double kr = 2.0 * M_PI / wl;
    double sp_radius = 0.6 * wl;  /* matches NR_UCA_RADIUS default */

    cdouble_t resp_desired = 0.0;
    cdouble_t resp_interf = 0.0;
    for (int m = 0; m < 6; m++) {
        double phi_m = 2.0 * M_PI * m / 6.0;
        double x_m = sp_radius * cos(phi_m);
        double y_m = sp_radius * sin(phi_m);

        double phase_d = kr * (x_m * cos(90.0 * M_PI / 180.0) +
                               y_m * sin(90.0 * M_PI / 180.0));
        resp_desired += conj(weights[m]) * cexp(I * phase_d);

        double phase_i = kr * (x_m * cos(270.0 * M_PI / 180.0) +
                               y_m * sin(270.0 * M_PI / 180.0));
        resp_interf += conj(weights[m]) * cexp(I * phase_i);
    }
    TEST_ASSERT(cabs(resp_interf) < cabs(resp_desired),
                "Spatial LCMV: null direction response < desired response");

    spatial_canceller_free(&canceller);
}

/* ============================================================================
 * Test 5: Diagnostic Logger
 * ============================================================================ */

static void test_diagnostic_logger(void) {
    printf("\n--- Test: Diagnostic Logger ---\n");

    diag_logger_t logger;

    /* Init with temp path */
#ifdef _WIN32
    const char* base = "C:\\Users\\Public\\mqp_test_log";
#else
    const char* base = "/tmp/mqp_test_log";
#endif

    int ret = diag_logger_init(&logger, base, DIAG_LEVEL_DEBUG);
    TEST_ASSERT(ret == 0, "Logger init succeeds");
    TEST_ASSERT(logger.enabled, "Logger is enabled");

    /* Log signal entry */
    diag_signal_entry_t sig;
    memset(&sig, 0, sizeof(sig));
    sig.timestamp_s = diag_get_elapsed_s(&logger);
    sig.frequency_hz = 751e6;
    sig.snr_db = 15.0;
    sig.rssi_dbm = -85.0;
    sig.cell_id = 123;
    strcpy(sig.signal_type, "LTE");
    strcpy(sig.tower_name, "Park_Ave_ATT");

    ret = diag_log_signal(&logger, &sig);
    TEST_ASSERT(ret == 0, "Signal log write");

    /* Log AOA entry */
    diag_aoa_entry_t aoa;
    memset(&aoa, 0, sizeof(aoa));
    aoa.timestamp_s = diag_get_elapsed_s(&logger);
    aoa.azimuth_deg = 45.0;
    aoa.elevation_deg = 5.0;
    aoa.confidence = 0.85;
    aoa.expected_az = 47.0;
    aoa.az_error = -2.0;
    aoa.num_sources = 3;

    ret = diag_log_aoa(&logger, &aoa);
    TEST_ASSERT(ret == 0, "AOA log write");

    /* Log timing entry */
    diag_timing_entry_t timing;
    memset(&timing, 0, sizeof(timing));
    timing.timestamp_s = diag_get_elapsed_s(&logger);
    timing.sdr_capture_ms = 25.0;
    timing.signal_proc_ms = 15.0;
    timing.noise_reduce_ms = 30.0;
    timing.music_aoa_ms = 55.0;
    timing.protocol_ms = 20.0;
    timing.triangulation_ms = 5.0;
    timing.total_pipeline_ms = 150.0;
    timing.met_realtime = true;

    ret = diag_log_pipeline_timing(&logger, &timing);
    TEST_ASSERT(ret == 0, "Timing log write");

    /* Log error */
    ret = diag_log_error(&logger, DIAG_LEVEL_WARN, "TEST", "Test warning message");
    TEST_ASSERT(ret == 0, "Error log write");

    /* Check counters */
    TEST_ASSERT(logger.signal_count == 1, "Signal count == 1");
    TEST_ASSERT(logger.aoa_count == 1, "AOA count == 1");
    TEST_ASSERT(logger.timing_count == 1, "Timing count == 1");

    /* Get elapsed time */
    double elapsed = diag_get_elapsed_s(&logger);
    TEST_ASSERT(elapsed >= 0.0, "Elapsed time >= 0");

    diag_logger_close(&logger);
    TEST_ASSERT(!logger.enabled, "Logger disabled after close");
}

/* ============================================================================
 * Test 6: WPI Test Configuration
 * ============================================================================ */

static void test_wpi_config(void) {
    printf("\n--- Test: WPI Test Configuration ---\n");

    wpi_test_config_t config;

#ifdef _WIN32
    const char* log_path = "C:\\Users\\Public\\mqp_wpi_test";
#else
    const char* log_path = "/tmp/mqp_wpi_test";
#endif

    int ret = wpi_test_config_init(&config, log_path);
    TEST_ASSERT(ret == 0, "WPI config init succeeds");
    TEST_ASSERT(config.tower_db.active_site == SITE_WPI_CAMPUS, "Tower DB is WPI");
    TEST_ASSERT(config.num_zones == 3, "3 test zones defined");
    TEST_ASSERT(config.num_bands == 6, "6 test bands defined (4 LTE + 2 P25)");

    /* Check zone coordinates are near WPI */
    for (int i = 0; i < config.num_zones; i++) {
        double dist = cell_tower_db_calculate_distance(&config.tower_db,
                        config.zones[i].lat, config.zones[i].lon);
        TEST_ASSERT(dist < 1000.0, "Test zone is within 1km of WPI");
    }

    /* Check frequency bands */
    TEST_ASSERT(config.bands[0].center_freq_hz > 700e6, "Band 0 freq > 700 MHz");
    TEST_ASSERT(strcmp(config.bands[0].signal_type, "LTE") == 0, "Band 0 is LTE");
    TEST_ASSERT(strcmp(config.bands[4].signal_type, "P25") == 0, "Band 4 is P25");

    /* Test SDR tune frequency */
    double tune = wpi_get_sdr_tune_freq(&config, 0);
    TEST_ASSERT(tune > 0.0, "SDR tune freq > 0");
    tune = wpi_get_sdr_tune_freq(&config, -1);
    TEST_ASSERT(tune == 0.0, "Invalid band returns 0");

    /* Run self-test */
    ret = wpi_campus_self_test(&config);
    TEST_ASSERT(ret == 0, "Campus self-test passes");

    wpi_test_config_free(&config);
}

/* ============================================================================
 * Main
 * ============================================================================ */

/* ============================================================================
 * Test 7: Logger Error Handling (TW-3)
 * ============================================================================ */

static void test_logger_error_handling(void) {
    printf("\n--- Test: Logger Error Handling ---\n");

    /* Test NULL logger */
    int ret = diag_log_signal(NULL, NULL);
    TEST_ASSERT(ret != 0, "NULL logger returns error");

    /* Test init with invalid path */
    diag_logger_t bad_logger;
    ret = diag_logger_init(&bad_logger, "/nonexistent/deeply/nested/path/log",
                           DIAG_LEVEL_DEBUG);
    /* May succeed or fail depending on OS - just verify no crash */
    if (ret == 0) {
        diag_logger_close(&bad_logger);
    }
    TEST_ASSERT(1, "Invalid path init does not crash");

    /* Test logging to closed logger */
    diag_logger_t closed_logger;
    memset(&closed_logger, 0, sizeof(closed_logger));
    closed_logger.enabled = false;
    diag_signal_entry_t dummy_sig;
    memset(&dummy_sig, 0, sizeof(dummy_sig));
    ret = diag_log_signal(&closed_logger, &dummy_sig);
    TEST_ASSERT(ret == 0 || ret == -1, "Logging to disabled logger does not crash");
}

int main(void) {
    srand(42);  /* TW-2 fix: deterministic random for reproducible tests */

    printf("==============================================\n");
    printf(" MQP WPI Campus Integration Tests\n");
    printf("==============================================\n");

    test_wpi_tower_db();
    test_butterworth_filter();
    test_lcmv_beamforming();
    test_lcmv_spatial();
    test_diagnostic_logger();
    test_wpi_config();
    test_logger_error_handling();

    printf("\n==============================================\n");
    printf(" Results: %d/%d passed\n", tests_passed, tests_run);
    printf("==============================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
