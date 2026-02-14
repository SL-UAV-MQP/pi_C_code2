/**
 * @file diagnostic_logger.h
 * @brief Field diagnostic logger for real-time debugging
 *
 * Provides timestamped CSV logging for signal metrics, pipeline timing,
 * and error capture during field testing. Essential for WPI campus
 * debugging where real-time terminal access is limited during drone flights.
 *
 * Usage:
 *   diag_logger_t logger;
 *   diag_logger_init(&logger, "/tmp/mqp_flight_001.csv", DIAG_LEVEL_INFO);
 *   diag_log_signal(&logger, &signal_entry);
 *   diag_log_aoa(&logger, &aoa_entry);
 *   diag_log_pipeline_timing(&logger, &timing);
 *   diag_logger_close(&logger);
 */

#ifndef DIAGNOSTIC_LOGGER_H
#define DIAGNOSTIC_LOGGER_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>

/* ============================================================================
 * Configuration
 * ============================================================================ */

/** Maximum log file path length */
#define DIAG_MAX_PATH 256

/** Maximum log message length */
#define DIAG_MAX_MSG 512

/** Maximum number of pipeline stages to track */
#define DIAG_MAX_STAGES 16

/** Ring buffer size for in-memory log (last N entries) */
#define DIAG_RING_SIZE 256

/* ============================================================================
 * Log Levels
 * ============================================================================ */

typedef enum {
    DIAG_LEVEL_ERROR = 0,   /**< Errors only */
    DIAG_LEVEL_WARN  = 1,   /**< Warnings + errors */
    DIAG_LEVEL_INFO  = 2,   /**< Normal operation info */
    DIAG_LEVEL_DEBUG = 3,   /**< Verbose debug output */
    DIAG_LEVEL_TRACE = 4    /**< Trace-level (every sample) */
} diag_level_t;

/* ============================================================================
 * Log Entry Types
 * ============================================================================ */

/** Signal detection log entry */
typedef struct {
    double timestamp_s;       /**< Seconds since logger start */
    double frequency_hz;      /**< Detected signal frequency */
    double snr_db;            /**< Signal-to-noise ratio */
    double rssi_dbm;          /**< Received signal strength */
    int cell_id;              /**< LTE Physical Cell ID (-1 if N/A) */
    char signal_type[16];     /**< "LTE", "P25", "FM", "UNKNOWN" */
    char tower_name[64];      /**< Matched tower name (empty if none) */
} diag_signal_entry_t;

/** AOA estimation log entry */
typedef struct {
    double timestamp_s;       /**< Seconds since logger start */
    double azimuth_deg;       /**< Estimated azimuth (degrees) */
    double elevation_deg;     /**< Estimated elevation (degrees) */
    double confidence;        /**< Estimation confidence (0-1) */
    double expected_az;       /**< Expected azimuth from tower DB */
    double az_error;          /**< Azimuth error (estimated - expected) */
    int num_sources;          /**< Number of MUSIC sources detected */
    double music_spectrum_peak; /**< Peak value of MUSIC spectrum */
} diag_aoa_entry_t;

/** Pipeline timing entry */
typedef struct {
    double timestamp_s;       /**< Seconds since logger start */
    double sdr_capture_ms;    /**< SDR capture time */
    double signal_proc_ms;    /**< Signal processing time */
    double noise_reduce_ms;   /**< Noise reduction time */
    double music_aoa_ms;      /**< MUSIC AOA estimation time */
    double protocol_ms;       /**< Protocol detection time */
    double triangulation_ms;  /**< Triangulation time */
    double total_pipeline_ms; /**< Total pipeline time */
    bool met_realtime;        /**< Did pipeline meet <200ms target? */
} diag_timing_entry_t;

/** Position estimate log entry */
typedef struct {
    double timestamp_s;
    double est_lat;           /**< Estimated latitude */
    double est_lon;           /**< Estimated longitude */
    double est_alt;           /**< Estimated altitude (m) */
    double gdop;              /**< Geometric dilution of precision */
    double error_m;           /**< Error from known position (m, if available) */
    int num_bearings;         /**< Number of bearings used */
} diag_position_entry_t;

/* ============================================================================
 * Logger State
 * ============================================================================ */

typedef struct {
    /* File output */
    FILE* signal_file;        /**< Signal detection CSV */
    FILE* aoa_file;           /**< AOA estimation CSV */
    FILE* timing_file;        /**< Pipeline timing CSV */
    FILE* position_file;      /**< Position estimate CSV */
    FILE* error_file;         /**< Error/warning log */

    /* Configuration */
    diag_level_t level;       /**< Current log level */
    char base_path[DIAG_MAX_PATH]; /**< Base path for log files */
    bool enabled;             /**< Logger enabled flag */

    /* Timing */
    struct timespec start_time; /**< Logger start time */

    /* Statistics */
    uint32_t signal_count;
    uint32_t aoa_count;
    uint32_t timing_count;
    uint32_t position_count;
    uint32_t error_count;
    uint32_t warn_count;
} diag_logger_t;

/* ============================================================================
 * Initialization / Cleanup
 * ============================================================================ */

/**
 * @brief Initialize diagnostic logger
 *
 * Creates CSV files at base_path with suffixes:
 *   _signal.csv, _aoa.csv, _timing.csv, _position.csv, _errors.log
 *
 * @param logger Logger instance
 * @param base_path Base file path (e.g., "/tmp/mqp_flight_001")
 * @param level Minimum log level
 * @return 0 on success, -1 on failure
 */
int diag_logger_init(diag_logger_t* logger, const char* base_path, diag_level_t level);

/**
 * @brief Close logger and flush all files
 */
void diag_logger_close(diag_logger_t* logger);

/**
 * @brief Enable/disable logging at runtime
 */
void diag_logger_set_enabled(diag_logger_t* logger, bool enabled);

/**
 * @brief Change log level at runtime
 */
void diag_logger_set_level(diag_logger_t* logger, diag_level_t level);

/* ============================================================================
 * Logging Functions
 * ============================================================================ */

/**
 * @brief Log signal detection event
 */
int diag_log_signal(diag_logger_t* logger, const diag_signal_entry_t* entry);

/**
 * @brief Log AOA estimation result
 */
int diag_log_aoa(diag_logger_t* logger, const diag_aoa_entry_t* entry);

/**
 * @brief Log pipeline timing
 */
int diag_log_pipeline_timing(diag_logger_t* logger, const diag_timing_entry_t* entry);

/**
 * @brief Log position estimate
 */
int diag_log_position(diag_logger_t* logger, const diag_position_entry_t* entry);

/**
 * @brief Log error/warning message
 */
int diag_log_error(diag_logger_t* logger, diag_level_t level, const char* module, const char* msg);

/**
 * @brief Get elapsed time since logger start (seconds)
 */
double diag_get_elapsed_s(const diag_logger_t* logger);

/**
 * @brief Print summary statistics to stdout
 */
void diag_print_summary(const diag_logger_t* logger);

#endif /* DIAGNOSTIC_LOGGER_H */
