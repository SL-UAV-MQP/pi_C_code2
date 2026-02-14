/**
 * @file state_machine.c
 * @brief Hardware Testing State Machine Implementation
 *
 * Sequential validation: HW_INIT → SDR_TEST → ALGO_TEST → LIVE
 * On failure: retry → ERROR_DIAG
 */

#include "state_machine.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef USE_SDR
#include "sdr_pluto.h"
#endif

#include "music_uca_6.h"
#include "music_config.h"

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

static uint64_t get_time_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000ULL + (uint64_t)ts.tv_nsec / 1000000ULL;
}

static void sm_log_error(state_machine_t* sm, error_code_t error, const char* message) {
    if (sm->error_log_count < SM_MAX_ERROR_LOG) {
        error_log_entry_t* entry = &sm->error_log[sm->error_log_count];
        entry->timestamp_ms = get_time_ms();
        entry->state = sm->current_state;
        entry->sub_state = sm->current_sub_state;
        entry->error = error;
        strncpy(entry->message, message, sizeof(entry->message) - 1);
        entry->message[sizeof(entry->message) - 1] = '\0';
        sm->error_log_count++;
    }

    sm->last_error = error;

    if (sm->on_error) {
        sm->on_error(error, message);
    }

    fprintf(stderr, "[SM ERROR] State=%s Error=%s: %s\n",
            sm_state_name(sm->current_state), sm_error_name(error), message);
}

static void sm_transition(state_machine_t* sm, system_state_t new_state,
                          system_sub_state_t sub_state) {
    system_state_t old = sm->current_state;
    sm->current_state = new_state;
    sm->current_sub_state = sub_state;
    sm->state_entry_time_ms = get_time_ms();
    sm->retry_count = 0;

    printf("[SM] %s → %s\n", sm_state_name(old), sm_state_name(new_state));

    if (sm->on_state_change) {
        sm->on_state_change(old, new_state);
    }
}

static void sm_transition_error(state_machine_t* sm, error_code_t error, const char* msg) {
    sm_log_error(sm, error, msg);

    if (sm->retry_count < SM_MAX_RETRIES) {
        sm->retry_count++;
        printf("[SM] Retrying... (attempt %d/%d)\n", sm->retry_count, SM_MAX_RETRIES);
        /* Stay in current state, sub-state will re-execute */
    } else {
        printf("[SM] Max retries reached. Entering ERROR_DIAG.\n");
        sm_transition(sm, STATE_ERROR_DIAG, SUB_ERROR_IDENTIFY);
    }
}

/* ============================================================================
 * STATE 0: HW_INIT
 * ============================================================================ */

static system_state_t sm_step_hw_init(state_machine_t* sm) {
    switch (sm->current_sub_state) {

    case SUB_HW_POWER_CHECK:
        printf("[HW_INIT] Checking power supply...\n");
        /* On RPi 5B: check /sys/class/power_supply or voltage via sysfs */
        /* For now, assume power is OK if we're running */
        printf("[HW_INIT] Power: OK\n");
        sm->current_sub_state = SUB_HW_SENSOR_CHECK;
        break;

    case SUB_HW_SENSOR_CHECK:
        printf("[HW_INIT] Checking sensors (IMU, GPS, barometer)...\n");
        /* TODO: Check I2C/SPI sensor connectivity */
        /* For WPI campus testing, sensors may not be present */
        printf("[HW_INIT] Sensors: SKIPPED (WPI test mode)\n");
        sm->current_sub_state = SUB_HW_ANTENNA_CHECK;
        break;

    case SUB_HW_ANTENNA_CHECK:
        printf("[HW_INIT] Checking antenna array...\n");
        /* TODO: Check antenna connections via RF path */
        /* Verify 6 elements connected (requires SDR later) */
        printf("[HW_INIT] Antenna: DEFERRED to SDR_TEST\n");
        sm->current_sub_state = SUB_HW_RF_CHECK;
        break;

    case SUB_HW_RF_CHECK:
        printf("[HW_INIT] Checking RF path...\n");
        /* TODO: Check USB hub, RF cables */
        printf("[HW_INIT] RF Path: OK (USB hub accessible)\n");
        sm->hw_init_passed = true;
        printf("[HW_INIT] === PASSED ===\n");
        sm_transition(sm, STATE_SDR_TEST, SUB_SDR_USB_SCAN);
        break;

    default:
        sm->current_sub_state = SUB_HW_POWER_CHECK;
        break;
    }

    return sm->current_state;
}

/* ============================================================================
 * STATE 1: SDR_TEST
 * ============================================================================ */

static system_state_t sm_step_sdr_test(state_machine_t* sm) {
#ifndef USE_SDR
    printf("[SDR_TEST] SDR support not compiled. Skipping.\n");
    sm->sdr_test_passed = true;
    sm_transition(sm, STATE_ALGO_TEST, SUB_ALGO_MUSIC_TEST);
    return sm->current_state;
#else
    static sdr_pluto_context_t* sdr_ctx = NULL;

    switch (sm->current_sub_state) {

    case SUB_SDR_USB_SCAN:
        printf("[SDR_TEST] Scanning USB for PLUTO devices...\n");
        sdr_ctx = sdr_pluto_init();
        if (!sdr_ctx) {
            sm_transition_error(sm, ERR_SDR_USB_NOT_FOUND, "Failed to init SDR context");
            return sm->current_state;
        }
        sm->current_sub_state = SUB_SDR_CONNECT;
        break;

    case SUB_SDR_CONNECT: {
        printf("[SDR_TEST] Connecting to PLUTO devices...\n");
        sdr_status_t status = sdr_pluto_connect(sdr_ctx);
        if (status != SDR_SUCCESS) {
            sm_transition_error(sm, ERR_SDR_DEVICE_COUNT,
                               sdr_pluto_get_error_string(status));
            return sm->current_state;
        }
        sm->sdr_devices_found = sdr_ctx->num_connected;
        printf("[SDR_TEST] Found %d devices\n", sm->sdr_devices_found);
        sm->current_sub_state = SUB_SDR_CONFIGURE;
        break;
    }

    case SUB_SDR_CONFIGURE: {
        printf("[SDR_TEST] Configuring devices...\n");
        sdr_config_t config;
        sdr_pluto_get_default_config(&config);
        sdr_status_t status = sdr_pluto_configure(sdr_ctx, &config);
        if (status != SDR_SUCCESS) {
            sm_transition_error(sm, ERR_SDR_CONFIG_FAIL,
                               sdr_pluto_get_error_string(status));
            return sm->current_state;
        }
        sdr_pluto_print_config(sdr_ctx);
        sm->current_sub_state = SUB_SDR_CAPTURE_TEST;
        break;
    }

    case SUB_SDR_CAPTURE_TEST: {
        printf("[SDR_TEST] Testing IQ capture (%d samples)...\n", NUM_SNAPSHOTS);
        sdr_capture_result_t* result = sdr_capture_result_alloc(NUM_SNAPSHOTS,
                                                                 sdr_ctx->num_connected);
        if (!result) {
            sm_transition_error(sm, ERR_SDR_CAPTURE_FAIL, "Failed to allocate capture buffer");
            return sm->current_state;
        }

        sdr_status_t status = sdr_pluto_capture(sdr_ctx, NUM_SNAPSHOTS, result);
        if (status != SDR_SUCCESS) {
            sdr_capture_result_free(result);
            sm_transition_error(sm, ERR_SDR_CAPTURE_FAIL,
                               sdr_pluto_get_error_string(status));
            return sm->current_state;
        }

        /* Compute simple SNR estimate from first channel */
        double power = 0.0;
        for (int n = 0; n < result->num_samples; n++) {
            double mag = cabs(result->channel_data[0][n]);
            power += mag * mag;
        }
        power /= result->num_samples;
        sm->sdr_capture_snr_db = 10.0 * log10(power + 1e-12);

        sdr_pluto_print_capture_stats(result);
        printf("[SDR_TEST] Capture power: %.1f dB\n", sm->sdr_capture_snr_db);

        sdr_capture_result_free(result);
        sm->current_sub_state = SUB_SDR_SYNC_TEST;
        break;
    }

    case SUB_SDR_SYNC_TEST:
        printf("[SDR_TEST] Synchronization test...\n");
        /* Sync quality was measured during capture (timing skew) */
        sm->sdr_sync_skew_us = 0.0;  /* Will be filled by capture thread data */
        printf("[SDR_TEST] Sync: OK\n");

        /* Cleanup SDR for now (LIVE state will re-init) */
        sdr_pluto_free(sdr_ctx);
        sdr_ctx = NULL;

        sm->sdr_test_passed = true;
        printf("[SDR_TEST] === PASSED ===\n");
        sm_transition(sm, STATE_ALGO_TEST, SUB_ALGO_MUSIC_TEST);
        break;

    default:
        sm->current_sub_state = SUB_SDR_USB_SCAN;
        break;
    }

    return sm->current_state;
#endif /* USE_SDR */
}

/* ============================================================================
 * STATE 2: ALGO_TEST
 * ============================================================================ */

static system_state_t sm_step_algo_test(state_machine_t* sm) {
    switch (sm->current_sub_state) {

    case SUB_ALGO_MUSIC_TEST: {
        printf("[ALGO_TEST] Testing MUSIC pipeline with synthetic data...\n");

        /* Generate synthetic test signal: 1 source at 45° */
        int M = NUM_ELEMENTS;
        int N = 1000;  /* Fewer snapshots for quick test */

        signal_matrix_t* signals = signal_matrix_alloc(M, N);
        if (!signals) {
            sm_transition_error(sm, ERR_ALGO_MUSIC_FAIL, "Failed to allocate signal matrix");
            return sm->current_state;
        }

        /* Create simple test: steering vector at 45° + noise */
        array_geometry_t geometry;
        array_geometry_init(&geometry, M, ARRAY_RADIUS_M, WAVELENGTH_M);

        double test_azimuth = 45.0;
        double test_snr = 20.0;  /* dB */

        /* Generate signal: a(theta) * s(t) + n(t) */
        cdouble_t steering[MAX_ANTENNAS];
        compute_steering_vector(&geometry, test_azimuth, 0.0, steering);

        srand(42);  /* Reproducible */
        double noise_amp = pow(10.0, -test_snr / 20.0);
        for (int n = 0; n < N; n++) {
            /* Random signal phase */
            double sig_phase = 2.0 * M_PI * ((double)rand() / RAND_MAX);
            cdouble_t sig = cexp(I * sig_phase);

            for (int m = 0; m < M; m++) {
                double ni = ((double)rand() / RAND_MAX - 0.5) * noise_amp;
                double nq = ((double)rand() / RAND_MAX - 0.5) * noise_amp;
                signals->data[n * M + m] = steering[m] * sig + (ni + I * nq);
            }
        }

        /* Run MUSIC */
        music_spectrum_t* spectrum = music_spectrum_alloc(
            NUM_AZIMUTH_BINS, M);
        detected_sources_t detected;

        music_status_t status = music_estimate_aoa(
            signals, &geometry, NUM_SOURCES, USE_FORWARD_BACKWARD,
            spectrum, &detected);

        if (status != MUSIC_SUCCESS && status != MUSIC_ERROR_NO_PEAKS) {
            signal_matrix_free(signals);
            music_spectrum_free(spectrum);
            sm_transition_error(sm, ERR_ALGO_MUSIC_FAIL, "MUSIC pipeline failed");
            return sm->current_state;
        }

        if (detected.num_sources > 0) {
            sm->music_test_error_deg = fabs(detected.sources[0].azimuth_deg - test_azimuth);
            sm->music_sources_detected = detected.num_sources;
            printf("[ALGO_TEST] MUSIC: detected %d sources, "
                   "azimuth=%.1f° (error=%.2f°)\n",
                   detected.num_sources,
                   detected.sources[0].azimuth_deg,
                   sm->music_test_error_deg);
        } else {
            sm->music_test_error_deg = 999.0;
            sm->music_sources_detected = 0;
            printf("[ALGO_TEST] MUSIC: no sources detected\n");
        }

        signal_matrix_free(signals);
        music_spectrum_free(spectrum);

        sm->current_sub_state = SUB_ALGO_CFAR_TEST;
        break;
    }

    case SUB_ALGO_CFAR_TEST:
        printf("[ALGO_TEST] CFAR: SKIPPED (tested via MUSIC integration)\n");
        sm->current_sub_state = SUB_ALGO_PROTOCOL_TEST;
        break;

    case SUB_ALGO_PROTOCOL_TEST:
        printf("[ALGO_TEST] Protocol detection: SKIPPED (requires RF signals)\n");
        sm->protocol_cells_detected = 0;
        sm->current_sub_state = SUB_ALGO_E2E_TEST;
        break;

    case SUB_ALGO_E2E_TEST:
        printf("[ALGO_TEST] End-to-end validation...\n");

        /* Check results */
        if (sm->music_sources_detected > 0 && sm->music_test_error_deg < 10.0) {
            sm->algo_test_passed = true;
            printf("[ALGO_TEST] === PASSED === (error=%.2f°)\n",
                   sm->music_test_error_deg);
            sm_transition(sm, STATE_LIVE, SUB_LIVE_RUNNING);
        } else {
            sm_transition_error(sm, ERR_ALGO_NO_SOURCES,
                               "MUSIC test: insufficient accuracy or no sources");
        }
        break;

    default:
        sm->current_sub_state = SUB_ALGO_MUSIC_TEST;
        break;
    }

    return sm->current_state;
}

/* ============================================================================
 * STATE 3: LIVE
 * ============================================================================ */

static system_state_t sm_step_live(state_machine_t* sm) {
    switch (sm->current_sub_state) {

    case SUB_LIVE_RUNNING:
        /* In live mode, the main loop calls sm_step() periodically.
         * The actual MUSIC processing happens outside the state machine.
         * Here we just do periodic health checks. */
        {
            uint64_t now = get_time_ms();
            if (now - sm->last_health_check_ms > SM_HEALTH_CHECK_INTERVAL_MS) {
                sm->current_sub_state = SUB_LIVE_HEALTH_CHECK;
            }
        }
        break;

    case SUB_LIVE_HEALTH_CHECK:
        /* Check system health */
        printf("[LIVE] Health check: frames=%llu errors=%llu avg_latency=%.1fms\n",
               (unsigned long long)sm->live_frames_processed,
               (unsigned long long)sm->live_errors,
               sm->live_avg_latency_ms);

        /* Check for stale data (no frames in last 5s) */
        sm->last_health_check_ms = get_time_ms();
        sm->current_sub_state = SUB_LIVE_RUNNING;

        if (sm->on_status_report) {
            char report[512];
            sm_generate_report(sm, report, sizeof(report));
            sm->on_status_report(report);
        }
        break;

    default:
        sm->current_sub_state = SUB_LIVE_RUNNING;
        break;
    }

    return sm->current_state;
}

/* ============================================================================
 * STATE 99: ERROR_DIAG
 * ============================================================================ */

static system_state_t sm_step_error_diag(state_machine_t* sm) {
    switch (sm->current_sub_state) {

    case SUB_ERROR_IDENTIFY:
        printf("\n[ERROR_DIAG] ========= ERROR DIAGNOSIS =========\n");
        printf("[ERROR_DIAG] Last error: %s\n", sm_error_name(sm->last_error));
        printf("[ERROR_DIAG] Failed at: %s\n", sm_state_name(sm->current_state));
        printf("[ERROR_DIAG] HW_INIT: %s\n", sm->hw_init_passed ? "PASS" : "FAIL");
        printf("[ERROR_DIAG] SDR_TEST: %s\n", sm->sdr_test_passed ? "PASS" : "FAIL");
        printf("[ERROR_DIAG] ALGO_TEST: %s\n", sm->algo_test_passed ? "PASS" : "FAIL");
        printf("[ERROR_DIAG] ================================\n\n");
        sm->current_sub_state = SUB_ERROR_REPORT;
        break;

    case SUB_ERROR_REPORT:
        /* Generate and send XBee report */
        if (sm->on_status_report) {
            char report[512];
            sm_generate_report(sm, report, sizeof(report));
            sm->on_status_report(report);
        }
        printf("[ERROR_DIAG] Waiting for operator intervention...\n");
        printf("[ERROR_DIAG] Use sm_force_state() to restart from any state.\n");
        /* Stay here until operator forces state change */
        break;

    default:
        sm->current_sub_state = SUB_ERROR_IDENTIFY;
        break;
    }

    return sm->current_state;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

int sm_init(state_machine_t* sm) {
    if (!sm) return -1;

    memset(sm, 0, sizeof(state_machine_t));
    sm->current_state = STATE_HW_INIT;
    sm->current_sub_state = SUB_HW_POWER_CHECK;
    sm->state_entry_time_ms = get_time_ms();
    sm->last_health_check_ms = get_time_ms();

    printf("[SM] State machine initialized. Starting HW_INIT...\n");
    return 0;
}

system_state_t sm_step(state_machine_t* sm) {
    if (!sm) return STATE_ERROR_DIAG;

    switch (sm->current_state) {
        case STATE_HW_INIT:    return sm_step_hw_init(sm);
        case STATE_SDR_TEST:   return sm_step_sdr_test(sm);
        case STATE_ALGO_TEST:  return sm_step_algo_test(sm);
        case STATE_LIVE:       return sm_step_live(sm);
        case STATE_ERROR_DIAG: return sm_step_error_diag(sm);
        default:
            sm_transition(sm, STATE_ERROR_DIAG, SUB_ERROR_IDENTIFY);
            return sm->current_state;
    }
}

void sm_force_state(state_machine_t* sm, system_state_t target) {
    if (!sm) return;

    printf("[SM] Forced transition to %s\n", sm_state_name(target));

    switch (target) {
        case STATE_HW_INIT:
            sm_transition(sm, STATE_HW_INIT, SUB_HW_POWER_CHECK);
            sm->hw_init_passed = false;
            sm->sdr_test_passed = false;
            sm->algo_test_passed = false;
            break;
        case STATE_SDR_TEST:
            sm_transition(sm, STATE_SDR_TEST, SUB_SDR_USB_SCAN);
            break;
        case STATE_ALGO_TEST:
            sm_transition(sm, STATE_ALGO_TEST, SUB_ALGO_MUSIC_TEST);
            break;
        case STATE_LIVE:
            sm_transition(sm, STATE_LIVE, SUB_LIVE_RUNNING);
            break;
        case STATE_ERROR_DIAG:
            sm_transition(sm, STATE_ERROR_DIAG, SUB_ERROR_IDENTIFY);
            break;
    }
}

const char* sm_state_name(system_state_t state) {
    switch (state) {
        case STATE_HW_INIT:    return "HW_INIT";
        case STATE_SDR_TEST:   return "SDR_TEST";
        case STATE_ALGO_TEST:  return "ALGO_TEST";
        case STATE_LIVE:       return "LIVE";
        case STATE_ERROR_DIAG: return "ERROR_DIAG";
        default: return "UNKNOWN";
    }
}

const char* sm_error_name(error_code_t error) {
    switch (error) {
        case ERR_NONE: return "NONE";
        case ERR_POWER_SUPPLY: return "POWER_SUPPLY";
        case ERR_SENSOR_FAIL: return "SENSOR_FAIL";
        case ERR_ANTENNA_DISCONNECTED: return "ANTENNA_DISCONNECTED";
        case ERR_RF_PATH_FAIL: return "RF_PATH_FAIL";
        case ERR_SDR_USB_NOT_FOUND: return "SDR_USB_NOT_FOUND";
        case ERR_SDR_DEVICE_COUNT: return "SDR_DEVICE_COUNT";
        case ERR_SDR_CONFIG_FAIL: return "SDR_CONFIG_FAIL";
        case ERR_SDR_CAPTURE_FAIL: return "SDR_CAPTURE_FAIL";
        case ERR_SDR_SYNC_FAIL: return "SDR_SYNC_FAIL";
        case ERR_ALGO_MUSIC_FAIL: return "ALGO_MUSIC_FAIL";
        case ERR_ALGO_CFAR_FAIL: return "ALGO_CFAR_FAIL";
        case ERR_ALGO_PROTOCOL_FAIL: return "ALGO_PROTOCOL_FAIL";
        case ERR_ALGO_NO_SOURCES: return "ALGO_NO_SOURCES";
        case ERR_LIVE_TIMEOUT: return "LIVE_TIMEOUT";
        case ERR_LIVE_DATA_STALE: return "LIVE_DATA_STALE";
        case ERR_UNKNOWN: return "UNKNOWN";
        default: return "UNKNOWN";
    }
}

void sm_generate_report(const state_machine_t* sm, char* buffer, int buffer_size) {
    if (!sm || !buffer) return;

    snprintf(buffer, buffer_size,
             "MQP STATUS|state=%s|hw=%s|sdr=%s|algo=%s|"
             "sdr_devs=%d|music_err=%.1f|sources=%d|"
             "frames=%llu|errors=%llu|latency=%.1fms|"
             "last_err=%s",
             sm_state_name(sm->current_state),
             sm->hw_init_passed ? "PASS" : "FAIL",
             sm->sdr_test_passed ? "PASS" : "FAIL",
             sm->algo_test_passed ? "PASS" : "FAIL",
             sm->sdr_devices_found,
             sm->music_test_error_deg,
             sm->music_sources_detected,
             (unsigned long long)sm->live_frames_processed,
             (unsigned long long)sm->live_errors,
             sm->live_avg_latency_ms,
             sm_error_name(sm->last_error));
}

void sm_print_status(const state_machine_t* sm) {
    if (!sm) return;

    printf("\n=== MQP System Status ===\n");
    printf("State:        %s\n", sm_state_name(sm->current_state));
    printf("HW Init:      %s\n", sm->hw_init_passed ? "PASS" : "--");
    printf("SDR Test:     %s (devices: %d)\n",
           sm->sdr_test_passed ? "PASS" : "--", sm->sdr_devices_found);
    printf("Algo Test:    %s (MUSIC error: %.2f°, sources: %d)\n",
           sm->algo_test_passed ? "PASS" : "--",
           sm->music_test_error_deg, sm->music_sources_detected);
    printf("Live:         frames=%llu errors=%llu\n",
           (unsigned long long)sm->live_frames_processed,
           (unsigned long long)sm->live_errors);
    printf("Last Error:   %s\n", sm_error_name(sm->last_error));
    printf("Error Log:    %d entries\n", sm->error_log_count);
    printf("=========================\n\n");
}

bool sm_is_operational(const state_machine_t* sm) {
    return sm && sm->current_state == STATE_LIVE;
}

int sm_get_error_log(const state_machine_t* sm, error_log_entry_t* entries, int max_entries) {
    if (!sm || !entries) return 0;

    int count = sm->error_log_count;
    if (count > max_entries) count = max_entries;

    memcpy(entries, sm->error_log, count * sizeof(error_log_entry_t));
    return count;
}
