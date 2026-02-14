/**
 * @file diagnostic_logger.c
 * @brief Field diagnostic logger for real-time debugging
 *
 * Writes timestamped CSV files for post-flight analysis.
 * All logging functions are non-blocking and safe for real-time use.
 */

#include "diagnostic_logger.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

static double get_elapsed(const diag_logger_t* logger) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - logger->start_time.tv_sec) +
           (now.tv_nsec - logger->start_time.tv_nsec) / 1e9;
}

static FILE* open_csv(const char* base_path, const char* suffix, const char* header) {
    char path[DIAG_MAX_PATH + 32];
    snprintf(path, sizeof(path), "%s%s", base_path, suffix);

    FILE* f = fopen(path, "w");
    if (f && header) {
        fprintf(f, "%s\n", header);
        fflush(f);
    }
    return f;
}

/* ============================================================================
 * Initialization / Cleanup
 * ============================================================================ */

int diag_logger_init(diag_logger_t* logger, const char* base_path, diag_level_t level) {
    if (!logger || !base_path) {
        return -1;
    }

    memset(logger, 0, sizeof(diag_logger_t));
    strncpy(logger->base_path, base_path, DIAG_MAX_PATH - 1);
    logger->level = level;
    logger->enabled = true;
    clock_gettime(CLOCK_MONOTONIC, &logger->start_time);

    /* Open CSV files with headers */
    logger->signal_file = open_csv(base_path, "_signal.csv",
        "timestamp_s,frequency_hz,snr_db,rssi_dbm,cell_id,signal_type,tower_name");

    logger->aoa_file = open_csv(base_path, "_aoa.csv",
        "timestamp_s,azimuth_deg,elevation_deg,confidence,expected_az,az_error,num_sources,spectrum_peak");

    logger->timing_file = open_csv(base_path, "_timing.csv",
        "timestamp_s,sdr_capture_ms,signal_proc_ms,noise_reduce_ms,music_aoa_ms,protocol_ms,triangulation_ms,total_ms,realtime_ok");

    logger->position_file = open_csv(base_path, "_position.csv",
        "timestamp_s,est_lat,est_lon,est_alt,gdop,error_m,num_bearings");

    logger->error_file = open_csv(base_path, "_errors.log",
        "timestamp_s,level,module,message");

    /* Check if at least error file opened */
    if (!logger->error_file) {
        diag_logger_close(logger);
        return -1;
    }

    /* Log startup */
    diag_log_error(logger, DIAG_LEVEL_INFO, "INIT",
                   "Diagnostic logger started");

    return 0;
}

void diag_logger_close(diag_logger_t* logger) {
    if (!logger) return;

    diag_log_error(logger, DIAG_LEVEL_INFO, "SHUTDOWN",
                   "Logger closing - printing summary");

    if (logger->signal_file)   fclose(logger->signal_file);
    if (logger->aoa_file)      fclose(logger->aoa_file);
    if (logger->timing_file)   fclose(logger->timing_file);
    if (logger->position_file) fclose(logger->position_file);
    if (logger->error_file)    fclose(logger->error_file);

    logger->signal_file = NULL;
    logger->aoa_file = NULL;
    logger->timing_file = NULL;
    logger->position_file = NULL;
    logger->error_file = NULL;
    logger->enabled = false;
}

void diag_logger_set_enabled(diag_logger_t* logger, bool enabled) {
    if (logger) logger->enabled = enabled;
}

void diag_logger_set_level(diag_logger_t* logger, diag_level_t level) {
    if (logger) logger->level = level;
}

/* ============================================================================
 * Logging Functions
 * ============================================================================ */

int diag_log_signal(diag_logger_t* logger, const diag_signal_entry_t* entry) {
    if (!logger || !logger->enabled || !logger->signal_file || !entry) {
        return -1;
    }

    if (logger->level < DIAG_LEVEL_INFO) return 0;

    fprintf(logger->signal_file,
            "%.6f,%.1f,%.2f,%.2f,%d,%s,%s\n",
            entry->timestamp_s,
            entry->frequency_hz,
            entry->snr_db,
            entry->rssi_dbm,
            entry->cell_id,
            entry->signal_type,
            entry->tower_name);

    /* Flush periodically (every 10 entries) for crash safety */
    logger->signal_count++;
    if (logger->signal_count % 10 == 0) {
        fflush(logger->signal_file);
    }

    return 0;
}

int diag_log_aoa(diag_logger_t* logger, const diag_aoa_entry_t* entry) {
    if (!logger || !logger->enabled || !logger->aoa_file || !entry) {
        return -1;
    }

    if (logger->level < DIAG_LEVEL_INFO) return 0;

    fprintf(logger->aoa_file,
            "%.6f,%.2f,%.2f,%.4f,%.2f,%.2f,%d,%.4f\n",
            entry->timestamp_s,
            entry->azimuth_deg,
            entry->elevation_deg,
            entry->confidence,
            entry->expected_az,
            entry->az_error,
            entry->num_sources,
            entry->music_spectrum_peak);

    logger->aoa_count++;
    if (logger->aoa_count % 10 == 0) {
        fflush(logger->aoa_file);
    }

    return 0;
}

int diag_log_pipeline_timing(diag_logger_t* logger, const diag_timing_entry_t* entry) {
    if (!logger || !logger->enabled || !logger->timing_file || !entry) {
        return -1;
    }

    if (logger->level < DIAG_LEVEL_INFO) return 0;

    fprintf(logger->timing_file,
            "%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d\n",
            entry->timestamp_s,
            entry->sdr_capture_ms,
            entry->signal_proc_ms,
            entry->noise_reduce_ms,
            entry->music_aoa_ms,
            entry->protocol_ms,
            entry->triangulation_ms,
            entry->total_pipeline_ms,
            entry->met_realtime ? 1 : 0);

    logger->timing_count++;
    if (logger->timing_count % 5 == 0) {
        fflush(logger->timing_file);
    }

    /* Warn if real-time target missed */
    if (!entry->met_realtime) {
        char msg[DIAG_MAX_MSG];
        snprintf(msg, sizeof(msg),
                 "Pipeline missed real-time target: %.1f ms (limit 200 ms)",
                 entry->total_pipeline_ms);
        diag_log_error(logger, DIAG_LEVEL_WARN, "TIMING", msg);
    }

    return 0;
}

int diag_log_position(diag_logger_t* logger, const diag_position_entry_t* entry) {
    if (!logger || !logger->enabled || !logger->position_file || !entry) {
        return -1;
    }

    if (logger->level < DIAG_LEVEL_INFO) return 0;

    fprintf(logger->position_file,
            "%.6f,%.8f,%.8f,%.2f,%.4f,%.2f,%d\n",
            entry->timestamp_s,
            entry->est_lat,
            entry->est_lon,
            entry->est_alt,
            entry->gdop,
            entry->error_m,
            entry->num_bearings);

    logger->position_count++;
    if (logger->position_count % 5 == 0) {
        fflush(logger->position_file);
    }

    return 0;
}

int diag_log_error(diag_logger_t* logger, diag_level_t level, const char* module, const char* msg) {
    if (!logger || !logger->error_file || !module || !msg) {
        return -1;
    }

    if (level > logger->level) return 0;

    const char* level_str;
    switch (level) {
        case DIAG_LEVEL_ERROR: level_str = "ERROR"; break;
        case DIAG_LEVEL_WARN:  level_str = "WARN";  break;
        case DIAG_LEVEL_INFO:  level_str = "INFO";  break;
        case DIAG_LEVEL_DEBUG: level_str = "DEBUG"; break;
        case DIAG_LEVEL_TRACE: level_str = "TRACE"; break;
        default:               level_str = "???";   break;
    }

    double elapsed = get_elapsed(logger);

    fprintf(logger->error_file, "%.6f,%s,%s,%s\n",
            elapsed, level_str, module, msg);

    /* Always flush errors and warnings immediately */
    if (level <= DIAG_LEVEL_WARN) {
        fflush(logger->error_file);
    }

    if (level == DIAG_LEVEL_ERROR) logger->error_count++;
    if (level == DIAG_LEVEL_WARN)  logger->warn_count++;

    return 0;
}

double diag_get_elapsed_s(const diag_logger_t* logger) {
    if (!logger) return 0.0;
    return get_elapsed(logger);
}

void diag_print_summary(const diag_logger_t* logger) {
    if (!logger) return;

    double elapsed = get_elapsed(logger);

    printf("\n=== Diagnostic Logger Summary ===\n");
    printf("Run time:         %.1f seconds\n", elapsed);
    printf("Signal entries:   %u\n", logger->signal_count);
    printf("AOA entries:      %u\n", logger->aoa_count);
    printf("Timing entries:   %u\n", logger->timing_count);
    printf("Position entries: %u\n", logger->position_count);
    printf("Errors:           %u\n", logger->error_count);
    printf("Warnings:         %u\n", logger->warn_count);

    if (logger->timing_count > 0) {
        printf("Avg pipeline rate: %.1f Hz\n",
               logger->timing_count / (elapsed + 1e-9));
    }

    printf("Log files: %s_*.csv\n", logger->base_path);
    printf("================================\n\n");
}
