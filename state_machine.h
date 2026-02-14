/**
 * @file state_machine.h
 * @brief Hardware Testing State Machine for MQP System
 *
 * Implements sequential hardware validation:
 *   STATE_HW_INIT  → Sensor/Antenna/RF checks
 *   STATE_SDR_TEST → USB discovery, configuration, IQ capture, sync
 *   STATE_ALGO_TEST → MUSIC pipeline, protocol detection, CFAR, end-to-end
 *   STATE_LIVE      → Continuous operation, health monitoring, error recovery
 *   STATE_ERROR_DIAG → Diagnostic mode with XBee reporting
 *
 * On failure: retry once, if still failing → STATE_ERROR_DIAG
 * Designed for WPI campus pre-testing before CMRCM field deployment.
 */

#ifndef STATE_MACHINE_H
#define STATE_MACHINE_H

#include "common_types.h"
#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * State Definitions
 * ============================================================================ */

/** System states */
typedef enum {
    STATE_HW_INIT    = 0,    /**< Hardware initialization & sensor checks */
    STATE_SDR_TEST   = 1,    /**< SDR device connection & capture test */
    STATE_ALGO_TEST  = 2,    /**< Algorithm pipeline validation */
    STATE_LIVE       = 3,    /**< Live continuous operation */
    STATE_ERROR_DIAG = 99    /**< Error diagnostics & reporting */
} system_state_t;

/** Sub-state for detailed progress tracking */
typedef enum {
    /* HW_INIT sub-states */
    SUB_HW_POWER_CHECK = 0,
    SUB_HW_SENSOR_CHECK,
    SUB_HW_ANTENNA_CHECK,
    SUB_HW_RF_CHECK,

    /* SDR_TEST sub-states */
    SUB_SDR_USB_SCAN = 10,
    SUB_SDR_CONNECT,
    SUB_SDR_CONFIGURE,
    SUB_SDR_CAPTURE_TEST,
    SUB_SDR_SYNC_TEST,

    /* ALGO_TEST sub-states */
    SUB_ALGO_MUSIC_TEST = 20,
    SUB_ALGO_CFAR_TEST,
    SUB_ALGO_PROTOCOL_TEST,
    SUB_ALGO_E2E_TEST,

    /* LIVE sub-states */
    SUB_LIVE_RUNNING = 30,
    SUB_LIVE_HEALTH_CHECK,

    /* ERROR sub-states */
    SUB_ERROR_IDENTIFY = 90,
    SUB_ERROR_REPORT,
    SUB_ERROR_RETRY
} system_sub_state_t;

/** Error codes for diagnostics */
typedef enum {
    ERR_NONE = 0,
    ERR_POWER_SUPPLY,
    ERR_SENSOR_FAIL,
    ERR_ANTENNA_DISCONNECTED,
    ERR_RF_PATH_FAIL,
    ERR_SDR_USB_NOT_FOUND,
    ERR_SDR_DEVICE_COUNT,
    ERR_SDR_CONFIG_FAIL,
    ERR_SDR_CAPTURE_FAIL,
    ERR_SDR_SYNC_FAIL,
    ERR_ALGO_MUSIC_FAIL,
    ERR_ALGO_CFAR_FAIL,
    ERR_ALGO_PROTOCOL_FAIL,
    ERR_ALGO_NO_SOURCES,
    ERR_LIVE_TIMEOUT,
    ERR_LIVE_DATA_STALE,
    ERR_UNKNOWN
} error_code_t;

/* ============================================================================
 * State Machine Context
 * ============================================================================ */

/** Maximum retry attempts per state */
#define SM_MAX_RETRIES 2

/** Health check interval in milliseconds */
#define SM_HEALTH_CHECK_INTERVAL_MS 5000

/** Maximum error log entries */
#define SM_MAX_ERROR_LOG 32

/** Error log entry */
typedef struct {
    uint64_t timestamp_ms;       /**< Error timestamp */
    system_state_t state;        /**< State when error occurred */
    system_sub_state_t sub_state; /**< Sub-state when error occurred */
    error_code_t error;          /**< Error code */
    char message[128];           /**< Human-readable message */
} error_log_entry_t;

/** State machine context */
typedef struct {
    system_state_t current_state;
    system_sub_state_t current_sub_state;

    /* Retry tracking */
    int retry_count;
    error_code_t last_error;

    /* Timing */
    uint64_t state_entry_time_ms;
    uint64_t last_health_check_ms;

    /* Test results */
    bool hw_init_passed;
    bool sdr_test_passed;
    bool algo_test_passed;

    /* SDR test results */
    int sdr_devices_found;
    double sdr_capture_snr_db;
    double sdr_sync_skew_us;

    /* Algorithm test results */
    double music_test_error_deg;
    int music_sources_detected;
    int protocol_cells_detected;

    /* Live operation stats */
    uint64_t live_frames_processed;
    uint64_t live_errors;
    double live_avg_latency_ms;

    /* Error log */
    error_log_entry_t error_log[SM_MAX_ERROR_LOG];
    int error_log_count;

    /* Callbacks (optional, for XBee reporting) */
    void (*on_state_change)(system_state_t old_state, system_state_t new_state);
    void (*on_error)(error_code_t error, const char* message);
    void (*on_status_report)(const char* report);
} state_machine_t;

/* ============================================================================
 * API Functions
 * ============================================================================ */

/**
 * @brief Initialize state machine
 * @param sm Pointer to state machine context (allocated by caller)
 * @return 0 on success
 */
int sm_init(state_machine_t* sm);

/**
 * @brief Run one step of the state machine
 *
 * Call this in a loop. Each call advances the state machine by one sub-step.
 * Returns the current state after the step.
 *
 * @param sm State machine context
 * @return Current system state after step
 */
system_state_t sm_step(state_machine_t* sm);

/**
 * @brief Force transition to a specific state
 * @param sm State machine context
 * @param target Target state
 */
void sm_force_state(state_machine_t* sm, system_state_t target);

/**
 * @brief Get human-readable state name
 * @param state System state
 * @return State name string
 */
const char* sm_state_name(system_state_t state);

/**
 * @brief Get human-readable error name
 * @param error Error code
 * @return Error name string
 */
const char* sm_error_name(error_code_t error);

/**
 * @brief Generate status report string
 * @param sm State machine context
 * @param buffer Output buffer
 * @param buffer_size Buffer size
 */
void sm_generate_report(const state_machine_t* sm, char* buffer, int buffer_size);

/**
 * @brief Print current state machine status to stdout
 * @param sm State machine context
 */
void sm_print_status(const state_machine_t* sm);

/**
 * @brief Check if system is in operational (LIVE) state
 * @param sm State machine context
 * @return true if in LIVE state
 */
bool sm_is_operational(const state_machine_t* sm);

/**
 * @brief Get error log
 * @param sm State machine context
 * @param entries Output array
 * @param max_entries Maximum entries to return
 * @return Number of entries returned
 */
int sm_get_error_log(const state_machine_t* sm, error_log_entry_t* entries, int max_entries);

#endif /* STATE_MACHINE_H */
