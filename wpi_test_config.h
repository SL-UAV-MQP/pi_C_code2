/**
 * @file wpi_test_config.h
 * @brief WPI Campus test configuration for pre-deployment validation
 *
 * Defines test parameters, flight zones, expected signal characteristics,
 * and self-test routines for WPI campus testing before CMRCM deployment.
 *
 * WPI Campus Overview:
 *   Location: Worcester, MA (42.2746°N, 71.8065°W)
 *   Elevation: ~150m ASL (hilly terrain)
 *   Campus size: ~32 acres
 *
 * Recommended test zones (open areas for drone flight):
 *   1. WPI Quad (42.2743, -71.8070) - central open area
 *   2. Park Ave field (42.2760, -71.8035) - sports fields, best LOS
 *   3. Parking lot D (42.2730, -71.8075) - flat, fewer obstacles
 *
 * Worcester P25 Public Safety (for P25 detector testing):
 *   Worcester PD: 460.025, 460.125, 460.225, 460.475 MHz
 *   Worcester FD: 460.525, 460.575 MHz
 *   MA State Police: 154.920 MHz (VHF - out of PLUTO range at 850 MHz config)
 *
 * LTE bands available (for LTE detector + MUSIC AOA testing):
 *   Band 13 (Verizon):  746-756 MHz DL
 *   Band 5  (AT&T/VZW): 869-894 MHz DL
 *   Band 2  (TMO/ATT):  1930-1990 MHz DL
 *   Band 4  (VZW/ATT):  2110-2155 MHz DL
 */

#ifndef WPI_TEST_CONFIG_H
#define WPI_TEST_CONFIG_H

#include "cell_tower_db.h"
#include "diagnostic_logger.h"
#include <stdbool.h>

/* ============================================================================
 * WPI Campus Constants
 * ============================================================================ */

/* GPS coordinates for test zones */
#define WPI_QUAD_LAT       42.274300
#define WPI_QUAD_LON      -71.807000
#define WPI_PARK_AVE_LAT   42.276000
#define WPI_PARK_AVE_LON  -71.803500
#define WPI_PARKING_D_LAT  42.273000
#define WPI_PARKING_D_LON -71.807500

/* Maximum test altitude AGL (meters) - conservative for campus */
#define WPI_MAX_ALT_AGL    30.0

/* Expected SNR ranges for WPI campus (dB) */
#define WPI_LTE_SNR_MIN     5.0    /* Minimum expected LTE SNR */
#define WPI_LTE_SNR_MAX    35.0    /* Maximum expected LTE SNR */
#define WPI_P25_SNR_MIN     0.0    /* P25 may be weak on campus */
#define WPI_P25_SNR_MAX    20.0

/* Worcester P25 frequencies (MHz) */
#define WPI_P25_FREQ_1  460.025e6
#define WPI_P25_FREQ_2  460.125e6
#define WPI_P25_FREQ_3  460.225e6
#define WPI_P25_FREQ_4  460.475e6
#define WPI_P25_FREQ_5  460.525e6
#define WPI_P25_FREQ_6  460.575e6

/* LTE downlink center frequencies for PLUTO tuning (Hz) */
#define WPI_LTE_BAND13_CENTER  751e6    /* Verizon 700 MHz */
#define WPI_LTE_BAND5_CENTER   881.5e6  /* AT&T/VZW 850 MHz */
#define WPI_LTE_BAND2_CENTER   1960e6   /* T-Mobile/AT&T 1900 MHz */
#define WPI_LTE_BAND4_CENTER   2132.5e6 /* VZW/AT&T 2100 MHz */

/* ============================================================================
 * Test Configuration Structure
 * ============================================================================ */

/** Test flight zone */
typedef struct {
    char name[64];
    double lat;
    double lon;
    double recommended_alt_m;  /**< Recommended test altitude AGL */
    double clear_radius_m;     /**< Obstacle-free radius */
} test_zone_t;

/** Test frequency band */
typedef struct {
    char name[32];
    double center_freq_hz;
    double bandwidth_hz;
    char signal_type[16];      /**< "LTE" or "P25" */
    double expected_snr_min;
    double expected_snr_max;
} test_band_t;

/** WPI test configuration */
typedef struct {
    /* Site database */
    cell_tower_database_t tower_db;

    /* Test zones */
    test_zone_t zones[4];
    int num_zones;

    /* Test frequency bands */
    test_band_t bands[10];
    int num_bands;

    /* Logger */
    diag_logger_t logger;
    bool logger_initialized;

    /* Self-test results */
    bool hw_test_passed;
    bool sdr_test_passed;
    bool signal_detected;
} wpi_test_config_t;

/* ============================================================================
 * API Functions
 * ============================================================================ */

/**
 * @brief Initialize WPI campus test configuration
 *
 * Sets up tower database, test zones, frequency bands, and logger.
 *
 * @param config Test configuration (allocated by caller)
 * @param log_base_path Base path for log files (e.g., "/tmp/wpi_test_001")
 * @return 0 on success, -1 on failure
 */
int wpi_test_config_init(wpi_test_config_t* config, const char* log_base_path);

/**
 * @brief Clean up test configuration
 */
void wpi_test_config_free(wpi_test_config_t* config);

/**
 * @brief Run pre-flight self-test on WPI campus
 *
 * Checks:
 *   1. SDR hardware connectivity (if USE_SDR)
 *   2. GPS reference point validity
 *   3. Tower database loaded correctly
 *   4. Logger writing to disk
 *   5. Algorithm pipeline smoke test (synthetic signal)
 *
 * @param config Test configuration
 * @return 0 if all tests pass, -1 if any fail
 */
int wpi_campus_self_test(wpi_test_config_t* config);

/**
 * @brief Print test configuration summary
 */
void wpi_test_config_print(const wpi_test_config_t* config);

/**
 * @brief Get recommended PLUTO SDR tuning frequency for a test band
 *
 * @param config Test configuration
 * @param band_index Index into bands array
 * @return Center frequency in Hz, or 0 on error
 */
double wpi_get_sdr_tune_freq(const wpi_test_config_t* config, int band_index);

#endif /* WPI_TEST_CONFIG_H */
