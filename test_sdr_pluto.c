/**
 * @file test_sdr_pluto.c
 * @brief Test program for ADALM-PLUTO SDR interface
 *
 * Demonstrates complete workflow:
 * 1. Initialize and connect to 3x PLUTO SDRs
 * 2. Configure for MUSIC processing (1050 MHz, 20 MSPS)
 * 3. Capture synchronized IQ samples
 * 4. Convert to signal_matrix_t format
 * 5. Run MUSIC AOA estimation
 *
 * @author SDR Specialist Agent
 * @date 2025-12-11
 */

#include "sdr_pluto.h"
#include "music_uca_6.h"
#include "array_geometry.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * @brief Test 1: Basic SDR initialization and connection
 */
int test_basic_connection(void) {
    printf("\n========================================\n");
    printf("Test 1: Basic SDR Connection\n");
    printf("========================================\n");

    // Initialize SDR context
    sdr_pluto_context_t* sdr_ctx = sdr_pluto_init();
    if (!sdr_ctx) {
        fprintf(stderr, "FAIL: Failed to initialize SDR context\n");
        return -1;
    }
    printf("PASS: SDR context initialized\n");

    // Connect to PLUTO devices
    sdr_status_t status = sdr_pluto_connect(sdr_ctx);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to connect to PLUTO devices: %s\n",
                sdr_pluto_get_error_string(status));
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Connected to %d PLUTO device(s)\n", sdr_ctx->num_connected);

    // Print device info
    for (int i = 0; i < sdr_ctx->num_connected; i++) {
        char info[256];
        status = sdr_pluto_get_device_info(&sdr_ctx->devices[i], info, sizeof(info));
        if (status == SDR_SUCCESS) {
            printf("      %s\n", info);
        }
    }

    // Cleanup
    sdr_pluto_free(sdr_ctx);
    printf("PASS: SDR cleanup successful\n");

    return 0;
}

/**
 * @brief Test 2: SDR configuration
 */
int test_configuration(void) {
    printf("\n========================================\n");
    printf("Test 2: SDR Configuration\n");
    printf("========================================\n");

    sdr_pluto_context_t* sdr_ctx = sdr_pluto_init();
    if (!sdr_ctx) {
        fprintf(stderr, "FAIL: Failed to initialize SDR context\n");
        return -1;
    }

    sdr_status_t status = sdr_pluto_connect(sdr_ctx);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to connect to PLUTO devices\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    // Get default configuration
    sdr_config_t config;
    status = sdr_pluto_get_default_config(&config);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to get default config\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Default configuration retrieved\n");

    // Customize for MUSIC processing
    config.center_freq_hz = 1050000000ULL;  // 1050 MHz (LTE Band 4)
    config.sample_rate_sps = 20000000ULL;   // 20 MSPS
    config.rf_bandwidth_hz = 30000000ULL;   // 30 MHz
    config.gain_mode = GAIN_MODE_MANUAL;
    config.gain_db = 60;                    // 60 dB gain
    config.buffer_size = 4000;              // 4000 snapshots for MUSIC

    // Apply configuration
    status = sdr_pluto_configure(sdr_ctx, &config);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to configure SDR: %s\n",
                sdr_pluto_get_error_string(status));
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: SDR configuration applied\n");

    // Print configuration
    sdr_pluto_print_config(sdr_ctx);

    // Test parameter changes
    printf("\nTesting parameter updates...\n");

    status = sdr_pluto_set_frequency(sdr_ctx, 900000000ULL);  // 900 MHz
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to set frequency\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Frequency updated to 900 MHz\n");

    status = sdr_pluto_set_gain(sdr_ctx, 50);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to set gain\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Gain updated to 50 dB\n");

    sdr_pluto_free(sdr_ctx);
    return 0;
}

/**
 * @brief Test 3: IQ data capture
 */
int test_capture(void) {
    printf("\n========================================\n");
    printf("Test 3: IQ Data Capture\n");
    printf("========================================\n");

    sdr_pluto_context_t* sdr_ctx = sdr_pluto_init();
    if (!sdr_ctx) {
        fprintf(stderr, "FAIL: Failed to initialize SDR context\n");
        return -1;
    }

    sdr_status_t status = sdr_pluto_connect(sdr_ctx);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to connect to PLUTO devices\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    // Configure SDR
    sdr_config_t config;
    sdr_pluto_get_default_config(&config);
    config.buffer_size = 1000;  // Small test capture
    status = sdr_pluto_configure(sdr_ctx, &config);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to configure SDR\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    // Allocate capture result
    int num_samples = 1000;
    int num_channels = sdr_ctx->num_connected * 2;
    sdr_capture_result_t* capture = sdr_capture_result_alloc(num_samples, num_channels);
    if (!capture) {
        fprintf(stderr, "FAIL: Failed to allocate capture result\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Capture buffer allocated (%d samples x %d channels)\n",
           num_samples, num_channels);

    // Capture IQ data
    status = sdr_pluto_capture(sdr_ctx, num_samples, capture);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to capture IQ data: %s\n",
                sdr_pluto_get_error_string(status));
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Captured %d samples\n", num_samples);

    // Print capture statistics
    sdr_pluto_print_capture_stats(capture);

    // Verify data (check for non-zero samples)
    int non_zero_count = 0;
    for (int ch = 0; ch < capture->num_channels; ch++) {
        for (int n = 0; n < 10 && n < num_samples; n++) {
            cdouble_t sample = capture->channel_data[ch][n];
            if (cabs(sample) > 0.001) {
                non_zero_count++;
            }
        }
    }

    if (non_zero_count > 0) {
        printf("PASS: Captured data contains %d non-zero samples\n", non_zero_count);
    } else {
        fprintf(stderr, "WARNING: All samples appear to be zero (no signal detected)\n");
    }

    // Print first few samples from channel 0
    printf("\nFirst 5 samples from channel 0:\n");
    for (int n = 0; n < 5 && n < num_samples; n++) {
        cdouble_t sample = capture->channel_data[0][n];
        printf("  [%d] = %.6f + j%.6f (mag=%.6f)\n",
               n, creal(sample), cimag(sample), cabs(sample));
    }

    sdr_capture_result_free(capture);
    sdr_pluto_free(sdr_ctx);
    return 0;
}

/**
 * @brief Test 4: Complete MUSIC workflow
 */
int test_music_workflow(void) {
    printf("\n========================================\n");
    printf("Test 4: Complete MUSIC Workflow\n");
    printf("========================================\n");

    // Initialize SDR
    sdr_pluto_context_t* sdr_ctx = sdr_pluto_init();
    if (!sdr_ctx) {
        fprintf(stderr, "FAIL: Failed to initialize SDR context\n");
        return -1;
    }

    sdr_status_t status = sdr_pluto_connect(sdr_ctx);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to connect to PLUTO devices\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    // Configure for MUSIC processing
    sdr_config_t config;
    sdr_pluto_get_default_config(&config);
    config.center_freq_hz = 1050000000ULL;
    config.sample_rate_sps = 20000000ULL;
    config.buffer_size = 4000;  // NUM_SNAPSHOTS
    status = sdr_pluto_configure(sdr_ctx, &config);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to configure SDR\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: SDR configured for MUSIC\n");

    // Capture IQ data
    int num_samples = 4000;
    int num_channels = 6;  // Required for MUSIC-UCA-6
    sdr_capture_result_t* capture = sdr_capture_result_alloc(num_samples, num_channels);
    if (!capture) {
        fprintf(stderr, "FAIL: Failed to allocate capture buffer\n");
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    status = sdr_pluto_capture(sdr_ctx, num_samples, capture);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to capture IQ data\n");
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Captured %d snapshots from %d channels\n",
           num_samples, num_channels);

    // Convert to signal matrix
    signal_matrix_t* signal_matrix = signal_matrix_alloc(num_channels, num_samples);
    if (!signal_matrix) {
        fprintf(stderr, "FAIL: Failed to allocate signal matrix\n");
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    status = sdr_pluto_convert_to_signal_matrix(capture, signal_matrix);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to convert to signal matrix\n");
        signal_matrix_free(signal_matrix);
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    printf("PASS: Converted to signal matrix (M=%d, N=%d)\n",
           signal_matrix->M, signal_matrix->N);

    // Create array geometry (UCA with 6 elements)
    array_geometry_t geometry_storage;
    double wavelength = 299792458.0 / 1050000000.0;  /* c / freq */
    music_status_t geom_status = array_geometry_init(
        &geometry_storage, 6, 0.176, wavelength);
    if (geom_status != MUSIC_SUCCESS) {
        fprintf(stderr, "FAIL: Failed to create array geometry\n");
        signal_matrix_free(signal_matrix);
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }
    array_geometry_t* geometry = &geometry_storage;
    printf("PASS: Created UCA geometry\n");

    // Allocate MUSIC spectrum and detection results
    int num_azimuth_bins = 1800;  // 360° / 0.2°
    music_spectrum_t* spectrum = music_spectrum_alloc(num_azimuth_bins, num_channels);
    detected_sources_t* detected = (detected_sources_t*)calloc(1, sizeof(detected_sources_t));

    if (!spectrum || !detected) {
        fprintf(stderr, "FAIL: Failed to allocate MUSIC results\n");
        if (spectrum) music_spectrum_free(spectrum);
        if (detected) free(detected);
        /* geometry is stack-allocated, no free needed */
        signal_matrix_free(signal_matrix);
        sdr_capture_result_free(capture);
        sdr_pluto_free(sdr_ctx);
        return -1;
    }

    // Run MUSIC AOA estimation
    printf("\nRunning MUSIC algorithm...\n");
    music_status_t music_status = music_estimate_aoa(
        signal_matrix,
        geometry,
        4,              // num_sources
        true,           // use_fb_averaging
        spectrum,
        detected
    );

    if (music_status != MUSIC_SUCCESS) {
        fprintf(stderr, "FAIL: MUSIC estimation failed (code %d)\n", music_status);
    } else {
        printf("PASS: MUSIC estimation successful\n");
        printf("\nDetected sources:\n");
        print_detected_sources(detected);
    }

    // Cleanup
    music_spectrum_free(spectrum);
    free(detected);
    /* geometry is stack-allocated, no free needed */
    signal_matrix_free(signal_matrix);
    sdr_capture_result_free(capture);
    sdr_pluto_free(sdr_ctx);

    return (music_status == MUSIC_SUCCESS) ? 0 : -1;
}

/**
 * @brief Main test runner
 */
int main(int argc, char* argv[]) {
    printf("\n");
    printf("========================================\n");
    printf("ADALM-PLUTO SDR Test Suite\n");
    printf("========================================\n");
    printf("\n");
    printf("This test suite requires:\n");
    printf("  - 3x ADALM-PLUTO SDR connected via USB\n");
    printf("  - LibIIO installed and configured\n");
    printf("  - Antenna connected to RX ports\n");
    printf("\n");

    int test_selection = 0;
    if (argc > 1) {
        test_selection = atoi(argv[1]);
    }

    int result = 0;

    if (test_selection == 0 || test_selection == 1) {
        result |= test_basic_connection();
    }

    if (test_selection == 0 || test_selection == 2) {
        result |= test_configuration();
    }

    if (test_selection == 0 || test_selection == 3) {
        result |= test_capture();
    }

    if (test_selection == 0 || test_selection == 4) {
        result |= test_music_workflow();
    }

    printf("\n========================================\n");
    if (result == 0) {
        printf("All tests PASSED\n");
    } else {
        printf("Some tests FAILED\n");
    }
    printf("========================================\n\n");

    return result;
}
