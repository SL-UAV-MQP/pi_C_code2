/**
 * @file wpi_test_config.c
 * @brief WPI Campus test configuration implementation
 */

#include "wpi_test_config.h"
#include "state_machine.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Initialization
 * ============================================================================ */

int wpi_test_config_init(wpi_test_config_t* config, const char* log_base_path) {
    if (!config) return -1;

    memset(config, 0, sizeof(wpi_test_config_t));

    /* Initialize tower database for WPI campus */
    if (cell_tower_db_init_site(&config->tower_db, SITE_WPI_CAMPUS) != 0) {
        return -1;
    }

    /* Define test zones */
    config->num_zones = 3;

    strncpy(config->zones[0].name, "WPI Quad (Central)", sizeof(config->zones[0].name) - 1);
    config->zones[0].name[sizeof(config->zones[0].name) - 1] = '\0';
    config->zones[0].lat = WPI_QUAD_LAT;
    config->zones[0].lon = WPI_QUAD_LON;
    config->zones[0].recommended_alt_m = 15.0;
    config->zones[0].clear_radius_m = 30.0;

    strncpy(config->zones[1].name, "Park Ave Fields", sizeof(config->zones[1].name) - 1);
    config->zones[1].name[sizeof(config->zones[1].name) - 1] = '\0';
    config->zones[1].lat = WPI_PARK_AVE_LAT;
    config->zones[1].lon = WPI_PARK_AVE_LON;
    config->zones[1].recommended_alt_m = 25.0;
    config->zones[1].clear_radius_m = 80.0;

    strncpy(config->zones[2].name, "Parking Lot D", sizeof(config->zones[2].name) - 1);
    config->zones[2].name[sizeof(config->zones[2].name) - 1] = '\0';
    config->zones[2].lat = WPI_PARKING_D_LAT;
    config->zones[2].lon = WPI_PARKING_D_LON;
    config->zones[2].recommended_alt_m = 20.0;
    config->zones[2].clear_radius_m = 40.0;

    /* Define test frequency bands */
    config->num_bands = 6;

    /* LTE bands */
    strcpy(config->bands[0].name, "LTE Band 13 (VZW)");
    config->bands[0].center_freq_hz = WPI_LTE_BAND13_CENTER;
    config->bands[0].bandwidth_hz = 10e6;
    strcpy(config->bands[0].signal_type, "LTE");
    config->bands[0].expected_snr_min = WPI_LTE_SNR_MIN;
    config->bands[0].expected_snr_max = WPI_LTE_SNR_MAX;

    strcpy(config->bands[1].name, "LTE Band 5 (850)");
    config->bands[1].center_freq_hz = WPI_LTE_BAND5_CENTER;
    config->bands[1].bandwidth_hz = 10e6;
    strcpy(config->bands[1].signal_type, "LTE");
    config->bands[1].expected_snr_min = WPI_LTE_SNR_MIN;
    config->bands[1].expected_snr_max = WPI_LTE_SNR_MAX;

    strcpy(config->bands[2].name, "LTE Band 2 (1900)");
    config->bands[2].center_freq_hz = WPI_LTE_BAND2_CENTER;
    config->bands[2].bandwidth_hz = 20e6;
    strcpy(config->bands[2].signal_type, "LTE");
    config->bands[2].expected_snr_min = WPI_LTE_SNR_MIN;
    config->bands[2].expected_snr_max = WPI_LTE_SNR_MAX;

    strcpy(config->bands[3].name, "LTE Band 4 (2100)");
    config->bands[3].center_freq_hz = WPI_LTE_BAND4_CENTER;
    config->bands[3].bandwidth_hz = 20e6;
    strcpy(config->bands[3].signal_type, "LTE");
    config->bands[3].expected_snr_min = WPI_LTE_SNR_MIN;
    config->bands[3].expected_snr_max = WPI_LTE_SNR_MAX;

    /* P25 bands (Worcester PD) */
    strcpy(config->bands[4].name, "P25 Worcester PD 1");
    config->bands[4].center_freq_hz = WPI_P25_FREQ_1;
    config->bands[4].bandwidth_hz = 12500;
    strcpy(config->bands[4].signal_type, "P25");
    config->bands[4].expected_snr_min = WPI_P25_SNR_MIN;
    config->bands[4].expected_snr_max = WPI_P25_SNR_MAX;

    strcpy(config->bands[5].name, "P25 Worcester PD 2");
    config->bands[5].center_freq_hz = WPI_P25_FREQ_2;
    config->bands[5].bandwidth_hz = 12500;
    strcpy(config->bands[5].signal_type, "P25");
    config->bands[5].expected_snr_min = WPI_P25_SNR_MIN;
    config->bands[5].expected_snr_max = WPI_P25_SNR_MAX;

    /* Initialize logger */
    if (log_base_path) {
        if (diag_logger_init(&config->logger, log_base_path, DIAG_LEVEL_INFO) == 0) {
            config->logger_initialized = true;
        }
    }

    return 0;
}

void wpi_test_config_free(wpi_test_config_t* config) {
    if (!config) return;

    if (config->logger_initialized) {
        diag_print_summary(&config->logger);
        diag_logger_close(&config->logger);
    }

    /* WC-4 fix: zero out struct to prevent stale data access */
    memset(config, 0, sizeof(*config));
}

/* ============================================================================
 * Self-Test
 * ============================================================================ */

int wpi_campus_self_test(wpi_test_config_t* config) {
    if (!config) return -1;

    int pass_count = 0;
    int total_tests = 0;

    printf("\n=== WPI Campus Pre-Flight Self-Test ===\n\n");

    /* Test 1: Tower database */
    total_tests++;
    printf("[TEST 1] Tower database... ");
    if (config->tower_db.num_towers > 0 &&
        config->tower_db.active_site == SITE_WPI_CAMPUS) {
        printf("PASS (%d towers loaded)\n", config->tower_db.num_towers);
        pass_count++;
    } else {
        printf("FAIL (no towers or wrong site)\n");
    }

    /* Test 2: Tower azimuth sanity check */
    total_tests++;
    printf("[TEST 2] Tower azimuths... ");
    {
        double azimuths[MAX_TOWERS];
        int n = cell_tower_db_get_expected_azimuths(&config->tower_db,
                                                      azimuths, MAX_TOWERS);
        bool all_valid = true;
        for (int i = 0; i < n; i++) {
            if (azimuths[i] < 0.0 || azimuths[i] > 360.0) {
                all_valid = false;
                break;
            }
        }
        if (all_valid && n > 0) {
            printf("PASS (%d azimuths in [0,360])\n", n);
            pass_count++;
        } else {
            printf("FAIL\n");
        }
    }

    /* Test 3: Test zones within WPI campus bounds */
    total_tests++;
    printf("[TEST 3] Test zones... ");
    {
        bool zones_valid = true;
        for (int i = 0; i < config->num_zones; i++) {
            double dist = cell_tower_db_calculate_distance(&config->tower_db,
                            config->zones[i].lat, config->zones[i].lon);
            if (dist > 2000.0) {  /* More than 2km from reference = suspicious */
                zones_valid = false;
                printf("\n  Zone '%s' is %.0fm from reference (too far)\n",
                       config->zones[i].name, dist);
            }
        }
        if (zones_valid) {
            printf("PASS (%d zones configured)\n", config->num_zones);
            pass_count++;
        } else {
            printf("FAIL\n");
        }
    }

    /* Test 4: Frequency band sanity */
    total_tests++;
    printf("[TEST 4] Frequency bands... ");
    {
        bool bands_valid = true;
        for (int i = 0; i < config->num_bands; i++) {
            if (config->bands[i].center_freq_hz < 100e6 ||
                config->bands[i].center_freq_hz > 6e9) {
                bands_valid = false;
            }
        }
        if (bands_valid && config->num_bands > 0) {
            /* WC-2 fix: count band types dynamically */
            int lte_count = 0, p25_count = 0;
            for (int i = 0; i < config->num_bands; i++) {
                if (strncmp(config->bands[i].signal_type, "P25", 3) == 0) {
                    p25_count++;
                } else {
                    lte_count++;
                }
            }
            printf("PASS (%d bands: %d LTE, %d P25)\n",
                   config->num_bands, lte_count, p25_count);
            pass_count++;
        } else {
            printf("FAIL\n");
        }
    }

    /* Test 5: Logger */
    total_tests++;
    printf("[TEST 5] Diagnostic logger... ");
    if (config->logger_initialized) {
        diag_signal_entry_t test_entry;
        memset(&test_entry, 0, sizeof(test_entry));
        test_entry.timestamp_s = diag_get_elapsed_s(&config->logger);
        test_entry.frequency_hz = 751e6;
        test_entry.snr_db = 10.0;
        test_entry.cell_id = -1;
        strcpy(test_entry.signal_type, "TEST");
        strcpy(test_entry.tower_name, "self_test");

        if (diag_log_signal(&config->logger, &test_entry) == 0) {
            printf("PASS (writing to %s_*.csv)\n", config->logger.base_path);
            pass_count++;
        } else {
            printf("FAIL (write error)\n");
        }
    } else {
        printf("SKIP (no log path configured)\n");
        pass_count++;  /* Not a failure, just not configured */
    }

    /* Summary */
    printf("\n=== Self-Test Result: %d/%d passed ===\n\n",
           pass_count, total_tests);

    config->hw_test_passed = (pass_count == total_tests);

    if (config->logger_initialized) {
        char msg[DIAG_MAX_MSG];
        snprintf(msg, sizeof(msg), "Self-test: %d/%d passed", pass_count, total_tests);
        diag_log_error(&config->logger,
                       config->hw_test_passed ? DIAG_LEVEL_INFO : DIAG_LEVEL_ERROR,
                       "SELF_TEST", msg);
    }

    return config->hw_test_passed ? 0 : -1;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

void wpi_test_config_print(const wpi_test_config_t* config) {
    if (!config) return;

    printf("\n=== WPI Campus Test Configuration ===\n");
    printf("Site: %s (%.6f, %.6f)\n",
           config->tower_db.site_name,
           config->tower_db.ref_lat,
           config->tower_db.ref_lon);

    printf("\nTest Zones:\n");
    for (int i = 0; i < config->num_zones; i++) {
        printf("  %d. %s\n", i + 1, config->zones[i].name);
        printf("     GPS: %.6f, %.6f\n", config->zones[i].lat, config->zones[i].lon);
        printf("     Alt: %.0fm AGL, Clear: %.0fm radius\n",
               config->zones[i].recommended_alt_m,
               config->zones[i].clear_radius_m);
    }

    printf("\nTest Frequency Bands:\n");
    for (int i = 0; i < config->num_bands; i++) {
        printf("  %d. %s\n", i + 1, config->bands[i].name);
        printf("     Center: %.3f MHz, BW: %.3f MHz (%s)\n",
               config->bands[i].center_freq_hz / 1e6,
               config->bands[i].bandwidth_hz / 1e6,
               config->bands[i].signal_type);
        printf("     Expected SNR: %.0f - %.0f dB\n",
               config->bands[i].expected_snr_min,
               config->bands[i].expected_snr_max);
    }

    printf("\nTower Database:\n");
    cell_tower_db_print_summary(&config->tower_db);
}

double wpi_get_sdr_tune_freq(const wpi_test_config_t* config, int band_index) {
    if (!config || band_index < 0 || band_index >= config->num_bands) {
        return 0.0;
    }
    return config->bands[band_index].center_freq_hz;
}
