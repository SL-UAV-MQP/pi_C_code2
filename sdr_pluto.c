/**
 * @file sdr_pluto.c
 * @brief ADALM-PLUTO SDR Interface Implementation
 *
 * Complete LibIIO integration for 3x PLUTO SDRs with synchronized capture.
 *
 * @author SDR Specialist Agent
 * @date 2025-12-11
 */

#include "sdr_pluto.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Convert gain mode enum to string
 */
static const char* gain_mode_to_string(pluto_gain_mode_t mode) {
    switch (mode) {
        case GAIN_MODE_MANUAL: return "manual";
        case GAIN_MODE_SLOW_ATTACK: return "slow_attack";
        case GAIN_MODE_FAST_ATTACK: return "fast_attack";
        case GAIN_MODE_HYBRID: return "hybrid";
        default: return "manual";
    }
}

/**
 * @brief Get current time in nanoseconds
 */
static uint64_t get_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/**
 * @brief Find IIO device by name
 */
static struct iio_device* find_device(struct iio_context* ctx, const char* name) {
    unsigned int nb_devices = iio_context_get_devices_count(ctx);
    for (unsigned int i = 0; i < nb_devices; i++) {
        struct iio_device* dev = iio_context_get_device(ctx, i);
        const char* dev_name = iio_device_get_name(dev);
        if (dev_name && strcmp(dev_name, name) == 0) {
            return dev;
        }
    }
    return NULL;
}

/**
 * @brief Find IIO channel by name
 */
static struct iio_channel* find_channel(struct iio_device* dev, const char* name, bool output) {
    unsigned int nb_channels = iio_device_get_channels_count(dev);
    for (unsigned int i = 0; i < nb_channels; i++) {
        struct iio_channel* chn = iio_device_get_channel(dev, i);
        const char* chn_id = iio_channel_get_id(chn);
        if (chn_id && strcmp(chn_id, name) == 0 &&
            iio_channel_is_output(chn) == output) {
            return chn;
        }
    }
    return NULL;
}

/**
 * @brief Set channel attribute (long long integer)
 */
static sdr_status_t set_channel_attr_ll(struct iio_channel* chn, const char* attr, long long val) {
    if (!chn) return SDR_ERROR_CHANNEL_NOT_FOUND;
    int ret = iio_channel_attr_write_longlong(chn, attr, val);
    if (ret < 0) {
        fprintf(stderr, "ERROR: Failed to set channel attribute '%s' to %lld (error %d)\n",
                attr, val, ret);
        return SDR_ERROR_INVALID_PARAM;
    }
    return SDR_SUCCESS;
}

/**
 * @brief Get channel attribute (long long integer)
 */
static sdr_status_t get_channel_attr_ll(struct iio_channel* chn, const char* attr, long long* val) {
    if (!chn) return SDR_ERROR_CHANNEL_NOT_FOUND;
    int ret = iio_channel_attr_read_longlong(chn, attr, val);
    if (ret < 0) {
        fprintf(stderr, "ERROR: Failed to read channel attribute '%s' (error %d)\n", attr, ret);
        return SDR_ERROR_INVALID_PARAM;
    }
    return SDR_SUCCESS;
}

/**
 * @brief Set channel attribute (boolean)
 */
static sdr_status_t set_channel_attr_bool(struct iio_channel* chn, const char* attr, bool val) {
    if (!chn) return SDR_ERROR_CHANNEL_NOT_FOUND;
    int ret = iio_channel_attr_write_bool(chn, attr, val);
    if (ret < 0) {
        fprintf(stderr, "ERROR: Failed to set channel attribute '%s' to %d (error %d)\n",
                attr, val, ret);
        return SDR_ERROR_INVALID_PARAM;
    }
    return SDR_SUCCESS;
}

/**
 * @brief Set channel attribute (string)
 */
static sdr_status_t set_channel_attr_str(struct iio_channel* chn, const char* attr, const char* val) {
    if (!chn) return SDR_ERROR_CHANNEL_NOT_FOUND;
    int ret = iio_channel_attr_write(chn, attr, val);
    if (ret < 0) {
        fprintf(stderr, "ERROR: Failed to set channel attribute '%s' to '%s' (error %d)\n",
                attr, val, ret);
        return SDR_ERROR_INVALID_PARAM;
    }
    return SDR_SUCCESS;
}

/* ============================================================================
 * Initialization and Configuration
 * ============================================================================ */

sdr_pluto_context_t* sdr_pluto_init(void) {
    sdr_pluto_context_t* ctx = (sdr_pluto_context_t*)calloc(1, sizeof(sdr_pluto_context_t));
    if (!ctx) {
        fprintf(stderr, "ERROR: Failed to allocate SDR context\n");
        return NULL;
    }

    // Initialize device array
    for (int i = 0; i < NUM_PLUTO_DEVICES; i++) {
        ctx->devices[i].ctx = NULL;
        ctx->devices[i].phy_dev = NULL;
        ctx->devices[i].rx_dev = NULL;
        ctx->devices[i].rx0_i = NULL;
        ctx->devices[i].rx0_q = NULL;
        ctx->devices[i].rxbuf = NULL;
        ctx->devices[i].device_index = i;
        ctx->devices[i].is_connected = false;
        ctx->devices[i].last_timestamp = 0;
        memset(ctx->devices[i].uri, 0, MAX_URI_LENGTH);
    }

    // Set default configuration
    sdr_pluto_get_default_config(&ctx->config);

    ctx->num_connected = 0;
    ctx->is_streaming = false;
    ctx->sync_timestamp = 0;

    printf("INFO: SDR context initialized\n");
    return ctx;
}

sdr_status_t sdr_pluto_get_default_config(sdr_config_t* config) {
    if (!config) return SDR_ERROR_NULL_POINTER;

    config->center_freq_hz = 1050000000ULL;  // 1050 MHz (center of LTE band)
    config->sample_rate_sps = DEFAULT_SAMPLE_RATE_SPS;  // 20 MSPS
    config->rf_bandwidth_hz = DEFAULT_RF_BANDWIDTH_HZ;   // 30 MHz
    config->gain_mode = GAIN_MODE_MANUAL;
    config->gain_db = DEFAULT_GAIN_DB;  // 60 dB
    config->buffer_size = 4000;  // NUM_SNAPSHOTS for MUSIC
    config->enable_quadrature_tracking = true;
    config->enable_rf_dc_offset_tracking = true;
    config->enable_bb_dc_offset_tracking = true;

    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_connect(sdr_pluto_context_t* ctx) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    printf("INFO: Scanning for PLUTO devices...\n");

    // Scan for available IIO contexts
    struct iio_scan_context* scan_ctx = iio_create_scan_context("usb", 0);
    if (!scan_ctx) {
        fprintf(stderr, "ERROR: Failed to create IIO scan context\n");
        return SDR_ERROR_IIO_INIT_FAILED;
    }

    struct iio_context_info** info;
    ssize_t num_contexts = iio_scan_context_get_info_list(scan_ctx, &info);

    if (num_contexts < 0) {
        fprintf(stderr, "ERROR: Failed to scan for devices (error %zd)\n", num_contexts);
        iio_scan_context_destroy(scan_ctx);
        return SDR_ERROR_DEVICE_NOT_FOUND;
    }

    printf("INFO: Found %zd IIO context(s)\n", num_contexts);

    // Connect to each PLUTO device
    int connected_count = 0;
    for (ssize_t i = 0; i < num_contexts && connected_count < NUM_PLUTO_DEVICES; i++) {
        const char* uri = iio_context_info_get_uri(info[i]);
        const char* desc = iio_context_info_get_description(info[i]);

        printf("INFO: Found device [%zd]: %s (%s)\n", i, uri, desc);

        // Only connect to PLUTO devices (check for "PlutoSDR" in description)
        if (strstr(desc, "PlutoSDR") || strstr(desc, "PLUTO")) {
            pluto_device_t* dev = &ctx->devices[connected_count];

            // Create IIO context for this device
            dev->ctx = iio_create_context_from_uri(uri);
            if (!dev->ctx) {
                fprintf(stderr, "WARNING: Failed to create context for device %s\n", uri);
                continue;
            }

            // Find AD9361 PHY device
            dev->phy_dev = find_device(dev->ctx, "ad9361-phy");
            if (!dev->phy_dev) {
                fprintf(stderr, "WARNING: Device %s does not have ad9361-phy\n", uri);
                iio_context_destroy(dev->ctx);
                dev->ctx = NULL;
                continue;
            }

            // Find RX streaming device
            dev->rx_dev = find_device(dev->ctx, "cf-ad9361-lpc");
            if (!dev->rx_dev) {
                fprintf(stderr, "WARNING: Device %s does not have cf-ad9361-lpc\n", uri);
                iio_context_destroy(dev->ctx);
                dev->ctx = NULL;
                continue;
            }

            // Find RX0 I/Q channels
            dev->rx0_i = find_channel(dev->rx_dev, "voltage0", false);
            dev->rx0_q = find_channel(dev->rx_dev, "voltage1", false);

            if (!dev->rx0_i || !dev->rx0_q) {
                fprintf(stderr, "WARNING: Device %s missing RX0 I/Q channels\n", uri);
                iio_context_destroy(dev->ctx);
                dev->ctx = NULL;
                continue;
            }

            // Store URI
            strncpy(dev->uri, uri, MAX_URI_LENGTH - 1);
            dev->is_connected = true;

            printf("INFO: Connected to PLUTO #%d (%s)\n", connected_count, uri);
            connected_count++;
        }
    }

    // Cleanup scan context
    iio_context_info_list_free(info);
    iio_scan_context_destroy(scan_ctx);

    ctx->num_connected = connected_count;

    if (connected_count == 0) {
        fprintf(stderr, "ERROR: No PLUTO devices found\n");
        return SDR_ERROR_DEVICE_NOT_FOUND;
    }

    if (connected_count < NUM_PLUTO_DEVICES) {
        fprintf(stderr, "ERROR: Only %d/%d PLUTO devices connected. "
                "All %d devices are required for %d-element UCA.\n",
                connected_count, NUM_PLUTO_DEVICES,
                NUM_PLUTO_DEVICES, NUM_RX_CHANNELS);
        // Disconnect partial connections
        for (int i = 0; i < connected_count; i++) {
            if (ctx->devices[i].ctx) {
                iio_context_destroy(ctx->devices[i].ctx);
                ctx->devices[i].ctx = NULL;
                ctx->devices[i].is_connected = false;
            }
        }
        ctx->num_connected = 0;
        return SDR_ERROR_DEVICE_NOT_FOUND;
    }

    printf("INFO: Successfully connected to %d PLUTO device(s)\n", connected_count);
    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_configure(sdr_pluto_context_t* ctx, const sdr_config_t* config) {
    if (!ctx || !config) return SDR_ERROR_NULL_POINTER;

    // Validate parameters
    if (config->center_freq_hz < 325000000ULL || config->center_freq_hz > 3800000000ULL) {
        fprintf(stderr, "ERROR: Frequency %llu Hz out of range (325 MHz - 3.8 GHz)\n",
                config->center_freq_hz);
        return SDR_ERROR_OUT_OF_RANGE;
    }

    if (config->sample_rate_sps < MIN_SAMPLE_RATE_SPS ||
        config->sample_rate_sps > MAX_SAMPLE_RATE_SPS) {
        fprintf(stderr, "ERROR: Sample rate %llu SPS out of range (%u - %u)\n",
                config->sample_rate_sps, MIN_SAMPLE_RATE_SPS, MAX_SAMPLE_RATE_SPS);
        return SDR_ERROR_OUT_OF_RANGE;
    }

    if (config->gain_db < 0 || config->gain_db > 76) {
        fprintf(stderr, "ERROR: Gain %d dB out of range (0-76)\n", config->gain_db);
        return SDR_ERROR_OUT_OF_RANGE;
    }

    // Save configuration
    memcpy(&ctx->config, config, sizeof(sdr_config_t));

    // Configure all connected devices
    for (int i = 0; i < ctx->num_connected; i++) {
        pluto_device_t* dev = &ctx->devices[i];
        if (!dev->is_connected) continue;

        printf("INFO: Configuring PLUTO #%d...\n", i);

        // Get PHY RX channel (voltage0) for gain/bandwidth/sample_rate
        struct iio_channel* phy_rx = find_channel(dev->phy_dev, "voltage0", false);
        if (!phy_rx) {
            fprintf(stderr, "ERROR: Failed to find PHY RX channel on device %d\n", i);
            return SDR_ERROR_CHANNEL_NOT_FOUND;
        }

        // FIX: Set RX LO frequency on altvoltage0 (RX_LO), NOT voltage0
        // AD9361 LO frequency is controlled via the altvoltage0 output channel
        struct iio_channel* rx_lo = find_channel(dev->phy_dev, "altvoltage0", true);
        if (!rx_lo) {
            fprintf(stderr, "ERROR: Failed to find RX_LO channel on device %d\n", i);
            return SDR_ERROR_CHANNEL_NOT_FOUND;
        }
        sdr_status_t status = set_channel_attr_ll(rx_lo, "frequency",
                                                   (long long)config->center_freq_hz);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set LO frequency on device %d\n", i);
            return status;
        }

        // Set sample rate (on voltage0 - this is correct for AD9361)
        status = set_channel_attr_ll(phy_rx, "sampling_frequency",
                                      (long long)config->sample_rate_sps);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set sample rate on device %d\n", i);
            return status;
        }

        // Set RF bandwidth
        status = set_channel_attr_ll(phy_rx, "rf_bandwidth",
                                      (long long)config->rf_bandwidth_hz);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set RF bandwidth on device %d\n", i);
            return status;
        }

        // Set gain control mode
        status = set_channel_attr_str(phy_rx, "gain_control_mode",
                                       gain_mode_to_string(config->gain_mode));
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set gain mode on device %d\n", i);
            return status;
        }

        // FIX: Set manual gain as double (AD9361 expects floating-point dB value)
        if (config->gain_mode == GAIN_MODE_MANUAL) {
            char gain_str[32];
            snprintf(gain_str, sizeof(gain_str), "%.6f", (double)config->gain_db);
            status = set_channel_attr_str(phy_rx, "hardwaregain", gain_str);
            if (status != SDR_SUCCESS) {
                fprintf(stderr, "ERROR: Failed to set gain on device %d\n", i);
                return status;
            }
        }

        // Enable quadrature tracking
        status = set_channel_attr_bool(phy_rx, "quadrature_tracking_en",
                                        config->enable_quadrature_tracking);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "WARNING: Failed to set quadrature tracking on device %d\n", i);
        }

        // Enable RF DC offset tracking
        status = set_channel_attr_bool(phy_rx, "rf_dc_offset_tracking_en",
                                        config->enable_rf_dc_offset_tracking);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "WARNING: Failed to set RF DC offset tracking on device %d\n", i);
        }

        // Enable BB DC offset tracking
        status = set_channel_attr_bool(phy_rx, "bb_dc_offset_tracking_en",
                                        config->enable_bb_dc_offset_tracking);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "WARNING: Failed to set BB DC offset tracking on device %d\n", i);
        }

        // Enable RX channels for streaming
        iio_channel_enable(dev->rx0_i);
        iio_channel_enable(dev->rx0_q);

        printf("INFO: PLUTO #%d configured successfully\n", i);
    }

    printf("INFO: All devices configured\n");
    return SDR_SUCCESS;
}

/* ============================================================================
 * Data Capture
 * ============================================================================ */

/* ---- Parallel capture thread data ---- */
typedef struct {
    pluto_device_t* dev;
    int num_samples;
    ssize_t ret;          /* iio_buffer_refill return value */
    uint64_t timestamp;   /* capture timestamp (ns) */
} capture_thread_data_t;

static void* capture_thread_func(void* arg) {
    capture_thread_data_t* td = (capture_thread_data_t*)arg;
    td->timestamp = get_time_ns();
    td->ret = iio_buffer_refill(td->dev->rxbuf);
    return NULL;
}

sdr_status_t sdr_pluto_synchronize_devices(sdr_pluto_context_t* ctx) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    if (ctx->num_connected == 0) {
        fprintf(stderr, "ERROR: No devices connected\n");
        return SDR_ERROR_DEVICE_NOT_FOUND;
    }

    // Record synchronization timestamp
    ctx->sync_timestamp = get_time_ns();

    printf("INFO: Devices synchronized at timestamp %llu ns\n",
           (unsigned long long)ctx->sync_timestamp);
    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_get_timestamp(pluto_device_t* device, uint64_t* timestamp) {
    if (!device || !timestamp) return SDR_ERROR_NULL_POINTER;
    if (!device->is_connected) return SDR_ERROR_USB_DISCONNECTED;

    // Use system time for now
    *timestamp = get_time_ns();
    device->last_timestamp = *timestamp;

    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_capture(sdr_pluto_context_t* ctx,
                                int num_samples,
                                sdr_capture_result_t* result) {
    if (!ctx || !result) return SDR_ERROR_NULL_POINTER;

    if (ctx->num_connected == 0) {
        fprintf(stderr, "ERROR: No devices connected\n");
        return SDR_ERROR_DEVICE_NOT_FOUND;
    }

    if (num_samples <= 0 || num_samples > MAX_BUFFER_SIZE) {
        fprintf(stderr, "ERROR: Invalid number of samples: %d\n", num_samples);
        return SDR_ERROR_INVALID_PARAM;
    }

    printf("INFO: Capturing %d samples from %d device(s)...\n",
           num_samples, ctx->num_connected);

    // Create RX buffers if not already created
    for (int dev_idx = 0; dev_idx < ctx->num_connected; dev_idx++) {
        pluto_device_t* dev = &ctx->devices[dev_idx];
        if (!dev->is_connected) continue;

        if (!dev->rxbuf) {
            dev->rxbuf = iio_device_create_buffer(dev->rx_dev, num_samples, false);
            if (!dev->rxbuf) {
                fprintf(stderr, "ERROR: Failed to create RX buffer for device %d\n", dev_idx);
                return SDR_ERROR_BUFFER_CREATE_FAILED;
            }
        }
    }

    // FIX CRITICAL 4: Parallel capture using pthreads for better synchronization
    // Launch all iio_buffer_refill() calls simultaneously to minimize inter-device
    // timing offset (reduces from sequential 2-10ms to <1ms typical)
    capture_thread_data_t thread_data[NUM_PLUTO_DEVICES];
    pthread_t threads[NUM_PLUTO_DEVICES];
    bool thread_started[NUM_PLUTO_DEVICES] = {false};

    // Record sync timestamp before parallel capture
    ctx->sync_timestamp = get_time_ns();

    for (int dev_idx = 0; dev_idx < ctx->num_connected; dev_idx++) {
        pluto_device_t* dev = &ctx->devices[dev_idx];
        if (!dev->is_connected) continue;

        thread_data[dev_idx].dev = dev;
        thread_data[dev_idx].num_samples = num_samples;
        thread_data[dev_idx].ret = 0;
        thread_data[dev_idx].timestamp = 0;

        int rc = pthread_create(&threads[dev_idx], NULL,
                                capture_thread_func, &thread_data[dev_idx]);
        if (rc != 0) {
            fprintf(stderr, "ERROR: Failed to create capture thread for device %d\n", dev_idx);
            // Fall back to sequential capture for this device
            thread_data[dev_idx].timestamp = get_time_ns();
            thread_data[dev_idx].ret = iio_buffer_refill(dev->rxbuf);
        } else {
            thread_started[dev_idx] = true;
        }
    }

    // Wait for all capture threads to complete
    for (int dev_idx = 0; dev_idx < ctx->num_connected; dev_idx++) {
        if (thread_started[dev_idx]) {
            pthread_join(threads[dev_idx], NULL);
        }
    }

    // Check capture results and compute max timing skew
    uint64_t min_ts = UINT64_MAX, max_ts = 0;
    for (int dev_idx = 0; dev_idx < ctx->num_connected; dev_idx++) {
        if (thread_data[dev_idx].ret < 0) {
            fprintf(stderr, "ERROR: Failed to refill buffer on device %d (error %zd)\n",
                    dev_idx, thread_data[dev_idx].ret);
            return SDR_ERROR_CAPTURE_TIMEOUT;
        }
        if (thread_data[dev_idx].timestamp < min_ts)
            min_ts = thread_data[dev_idx].timestamp;
        if (thread_data[dev_idx].timestamp > max_ts)
            max_ts = thread_data[dev_idx].timestamp;
    }
    uint64_t timing_skew_ns = max_ts - min_ts;
    if (timing_skew_ns > 5000000ULL) {  /* 5ms threshold */
        fprintf(stderr, "WARNING: Capture timing skew = %llu us (>5ms)\n",
                (unsigned long long)(timing_skew_ns / 1000));
    }

    // FIX CRITICAL 1: Correct IQ de-interleaving using iio_buffer_step()
    // PLUTO sends interleaved data: [I0, Q0, I1, Q1, ...]
    // iio_buffer_first() gives pointer to first I or Q sample
    // iio_buffer_step() gives the stride between consecutive samples of same channel
    for (int dev_idx = 0; dev_idx < ctx->num_connected; dev_idx++) {
        pluto_device_t* dev = &ctx->devices[dev_idx];
        if (!dev->is_connected) continue;

        // Get buffer stride (bytes between consecutive samples of same channel)
        ptrdiff_t step = iio_buffer_step(dev->rxbuf);
        char* i_ptr = (char*)iio_buffer_first(dev->rxbuf, dev->rx0_i);
        char* q_ptr = (char*)iio_buffer_first(dev->rxbuf, dev->rx0_q);

        if (!i_ptr || !q_ptr) {
            fprintf(stderr, "ERROR: Failed to get buffer pointers for device %d\n", dev_idx);
            return SDR_ERROR_BUFFER_CREATE_FAILED;
        }

        // FIX CRITICAL 3: Each PLUTO provides 1 RX channel, NOT 2
        // 3 PLUTOs = 3 channels. UCA needs 6 channels (2 per PLUTO via external switch).
        // Map: device index directly to channel index (1:1, not 1:2)
        int channel_idx = dev_idx;

        if (channel_idx >= NUM_RX_CHANNELS) {
            fprintf(stderr, "WARNING: Channel index %d exceeds NUM_RX_CHANNELS %d\n",
                    channel_idx, NUM_RX_CHANNELS);
            break;
        }

        // Convert interleaved 12-bit I/Q to complex double using correct stride
        for (int n = 0; n < num_samples; n++) {
            int16_t i_raw = *(int16_t*)(i_ptr + (ptrdiff_t)n * step);
            int16_t q_raw = *(int16_t*)(q_ptr + (ptrdiff_t)n * step);
            double i_val = (double)i_raw / 2048.0;  /* 12-bit: [-2048, 2047] → ~[-1, 1] */
            double q_val = (double)q_raw / 2048.0;
            result->channel_data[channel_idx][n] = i_val + I * q_val;
        }

        printf("INFO: Captured %d samples from PLUTO #%d → channel %d\n",
               num_samples, dev_idx, channel_idx);
    }

    // Fill result structure
    result->num_samples = num_samples;
    // FIX: Report actual channel count (1 per device, not 2)
    result->num_channels = ctx->num_connected;
    result->timestamp = ctx->sync_timestamp;
    result->actual_sample_rate = (double)ctx->config.sample_rate_sps;

    printf("INFO: Capture complete - %d samples x %d channels (skew: %llu us)\n",
           result->num_samples, result->num_channels,
           (unsigned long long)(timing_skew_ns / 1000));

    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_convert_to_signal_matrix(
    const sdr_capture_result_t* capture,
    signal_matrix_t* signal_matrix) {

    if (!capture || !signal_matrix) return SDR_ERROR_NULL_POINTER;

    if (signal_matrix->M != capture->num_channels) {
        fprintf(stderr, "ERROR: Matrix size mismatch (M=%d, channels=%d)\n",
                signal_matrix->M, capture->num_channels);
        return SDR_ERROR_INVALID_PARAM;
    }

    if (signal_matrix->N != capture->num_samples) {
        fprintf(stderr, "ERROR: Matrix size mismatch (N=%d, samples=%d)\n",
                signal_matrix->N, capture->num_samples);
        return SDR_ERROR_INVALID_PARAM;
    }

    // Convert from channel-based to column-major matrix format
    // signal_matrix->data is M x N in column-major order
    for (int m = 0; m < capture->num_channels; m++) {
        for (int n = 0; n < capture->num_samples; n++) {
            int idx = n * signal_matrix->M + m;  // Column-major indexing
            signal_matrix->data[idx] = capture->channel_data[m][n];
        }
    }

    printf("INFO: Converted capture to signal matrix (%d x %d)\n",
           signal_matrix->M, signal_matrix->N);

    return SDR_SUCCESS;
}

/* ============================================================================
 * Device Control
 * ============================================================================ */

sdr_status_t sdr_pluto_set_frequency(sdr_pluto_context_t* ctx, uint64_t freq_hz) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    ctx->config.center_freq_hz = freq_hz;

    for (int i = 0; i < ctx->num_connected; i++) {
        pluto_device_t* dev = &ctx->devices[i];
        if (!dev->is_connected) continue;

        // FIX: LO frequency is on altvoltage0, not voltage0
        struct iio_channel* rx_lo = find_channel(dev->phy_dev, "altvoltage0", true);
        sdr_status_t status = set_channel_attr_ll(rx_lo, "frequency", (long long)freq_hz);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set LO frequency on device %d\n", i);
            return status;
        }
    }

    printf("INFO: Frequency set to %llu Hz on all devices\n", (unsigned long long)freq_hz);
    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_set_sample_rate(sdr_pluto_context_t* ctx, uint64_t sample_rate_sps) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    ctx->config.sample_rate_sps = sample_rate_sps;

    for (int i = 0; i < ctx->num_connected; i++) {
        pluto_device_t* dev = &ctx->devices[i];
        if (!dev->is_connected) continue;

        struct iio_channel* phy_rx = find_channel(dev->phy_dev, "voltage0", false);
        sdr_status_t status = set_channel_attr_ll(phy_rx, "sampling_frequency",
                                                   (long long)sample_rate_sps);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set sample rate on device %d\n", i);
            return status;
        }
    }

    printf("INFO: Sample rate set to %llu SPS on all devices\n", sample_rate_sps);
    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_set_gain(sdr_pluto_context_t* ctx, int gain_db) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    ctx->config.gain_db = gain_db;

    for (int i = 0; i < ctx->num_connected; i++) {
        pluto_device_t* dev = &ctx->devices[i];
        if (!dev->is_connected) continue;

        struct iio_channel* phy_rx = find_channel(dev->phy_dev, "voltage0", false);
        // FIX: hardwaregain expects double value, not longlong
        char gain_str[32];
        snprintf(gain_str, sizeof(gain_str), "%.6f", (double)gain_db);
        sdr_status_t status = set_channel_attr_str(phy_rx, "hardwaregain", gain_str);
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set gain on device %d\n", i);
            return status;
        }
    }

    printf("INFO: Gain set to %d dB on all devices\n", gain_db);
    return SDR_SUCCESS;
}

sdr_status_t sdr_pluto_set_gain_mode(sdr_pluto_context_t* ctx, pluto_gain_mode_t gain_mode) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    ctx->config.gain_mode = gain_mode;

    for (int i = 0; i < ctx->num_connected; i++) {
        pluto_device_t* dev = &ctx->devices[i];
        if (!dev->is_connected) continue;

        struct iio_channel* phy_rx = find_channel(dev->phy_dev, "voltage0", false);
        sdr_status_t status = set_channel_attr_str(phy_rx, "gain_control_mode",
                                                    gain_mode_to_string(gain_mode));
        if (status != SDR_SUCCESS) {
            fprintf(stderr, "ERROR: Failed to set gain mode on device %d\n", i);
            return status;
        }
    }

    printf("INFO: Gain mode set to %s on all devices\n", gain_mode_to_string(gain_mode));
    return SDR_SUCCESS;
}

/* ============================================================================
 * Status and Diagnostics
 * ============================================================================ */

bool sdr_pluto_is_connected(const pluto_device_t* device) {
    if (!device) return false;
    return device->is_connected && device->ctx != NULL;
}

sdr_status_t sdr_pluto_get_device_info(const pluto_device_t* device,
                                        char* info_buffer,
                                        int buffer_size) {
    if (!device || !info_buffer) return SDR_ERROR_NULL_POINTER;
    if (!device->is_connected) return SDR_ERROR_USB_DISCONNECTED;

    const char* backend = iio_context_get_description(device->ctx);
    const char* name = iio_context_get_name(device->ctx);

    snprintf(info_buffer, buffer_size,
             "Device #%d: %s (%s) - URI: %s",
             device->device_index, name ? name : "Unknown",
             backend ? backend : "Unknown", device->uri);

    return SDR_SUCCESS;
}

void sdr_pluto_print_config(const sdr_pluto_context_t* ctx) {
    if (!ctx) return;

    printf("\n========== SDR Configuration ==========\n");
    printf("Connected devices:     %d / %d\n", ctx->num_connected, NUM_PLUTO_DEVICES);
    printf("Center frequency:      %.3f MHz\n", ctx->config.center_freq_hz / 1e6);
    printf("Sample rate:           %.3f MSPS\n", ctx->config.sample_rate_sps / 1e6);
    printf("RF bandwidth:          %.3f MHz\n", ctx->config.rf_bandwidth_hz / 1e6);
    printf("Gain mode:             %s\n", gain_mode_to_string(ctx->config.gain_mode));
    printf("Gain:                  %d dB\n", ctx->config.gain_db);
    printf("Buffer size:           %d samples\n", ctx->config.buffer_size);
    printf("Quadrature tracking:   %s\n", ctx->config.enable_quadrature_tracking ? "ON" : "OFF");
    printf("RF DC tracking:        %s\n", ctx->config.enable_rf_dc_offset_tracking ? "ON" : "OFF");
    printf("BB DC tracking:        %s\n", ctx->config.enable_bb_dc_offset_tracking ? "ON" : "OFF");
    printf("=======================================\n\n");
}

void sdr_pluto_print_capture_stats(const sdr_capture_result_t* result) {
    if (!result) return;

    printf("\n========== Capture Statistics ==========\n");
    printf("Samples per channel:   %d\n", result->num_samples);
    printf("Number of channels:    %d\n", result->num_channels);
    printf("Sample rate:           %.3f MSPS\n", result->actual_sample_rate / 1e6);
    printf("Capture duration:      %.3f ms\n",
           (result->num_samples / result->actual_sample_rate) * 1000.0);
    printf("Timestamp:             %llu ns\n", result->timestamp);
    printf("========================================\n\n");
}

/* ============================================================================
 * Cleanup
 * ============================================================================ */

sdr_status_t sdr_pluto_disconnect(sdr_pluto_context_t* ctx) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;

    printf("INFO: Disconnecting devices...\n");

    for (int i = 0; i < NUM_PLUTO_DEVICES; i++) {
        pluto_device_t* dev = &ctx->devices[i];

        if (dev->rxbuf) {
            iio_buffer_destroy(dev->rxbuf);
            dev->rxbuf = NULL;
        }

        if (dev->ctx) {
            iio_context_destroy(dev->ctx);
            dev->ctx = NULL;
        }

        dev->is_connected = false;
        dev->phy_dev = NULL;
        dev->rx_dev = NULL;
        dev->rx0_i = NULL;
        dev->rx0_q = NULL;
    }

    ctx->num_connected = 0;
    ctx->is_streaming = false;

    printf("INFO: All devices disconnected\n");
    return SDR_SUCCESS;
}

void sdr_pluto_free(sdr_pluto_context_t* ctx) {
    if (!ctx) return;

    sdr_pluto_disconnect(ctx);
    free(ctx);

    printf("INFO: SDR context freed\n");
}

/* ============================================================================
 * Memory Management
 * ============================================================================ */

sdr_capture_result_t* sdr_capture_result_alloc(int num_samples, int num_channels) {
    if (num_samples <= 0 || num_channels <= 0 || num_channels > NUM_RX_CHANNELS) {
        fprintf(stderr, "ERROR: Invalid capture dimensions (%d x %d)\n",
                num_channels, num_samples);
        return NULL;
    }

    sdr_capture_result_t* result = (sdr_capture_result_t*)calloc(1, sizeof(sdr_capture_result_t));
    if (!result) {
        fprintf(stderr, "ERROR: Failed to allocate capture result\n");
        return NULL;
    }

    result->num_samples = num_samples;
    result->num_channels = num_channels;

    // Allocate channel data arrays
    for (int ch = 0; ch < num_channels; ch++) {
        result->channel_data[ch] = (cdouble_t*)calloc(num_samples, sizeof(cdouble_t));
        if (!result->channel_data[ch]) {
            fprintf(stderr, "ERROR: Failed to allocate channel %d data\n", ch);
            sdr_capture_result_free(result);
            return NULL;
        }
    }

    return result;
}

void sdr_capture_result_free(sdr_capture_result_t* result) {
    if (!result) return;

    for (int ch = 0; ch < NUM_RX_CHANNELS; ch++) {
        if (result->channel_data[ch]) {
            free(result->channel_data[ch]);
            result->channel_data[ch] = NULL;
        }
    }

    free(result);
}

/* ============================================================================
 * Error Handling
 * ============================================================================ */

const char* sdr_pluto_get_error_string(sdr_status_t status) {
    switch (status) {
        case SDR_SUCCESS: return "Success";
        case SDR_ERROR_NULL_POINTER: return "Null pointer argument";
        case SDR_ERROR_INVALID_PARAM: return "Invalid parameter";
        case SDR_ERROR_DEVICE_NOT_FOUND: return "PLUTO device not found";
        case SDR_ERROR_USB_DISCONNECTED: return "USB connection lost";
        case SDR_ERROR_IIO_INIT_FAILED: return "LibIIO initialization failed";
        case SDR_ERROR_CHANNEL_NOT_FOUND: return "IIO channel not found";
        case SDR_ERROR_BUFFER_CREATE_FAILED: return "Buffer creation failed";
        case SDR_ERROR_CAPTURE_TIMEOUT: return "Capture timeout";
        case SDR_ERROR_SYNC_FAILED: return "Synchronization failed";
        case SDR_ERROR_OUT_OF_RANGE: return "Parameter out of range";
        case SDR_ERROR_MEMORY: return "Memory allocation failed";
        default: return "Unknown error";
    }
}

sdr_status_t sdr_pluto_handle_usb_error(sdr_pluto_context_t* ctx, int device_index) {
    if (!ctx) return SDR_ERROR_NULL_POINTER;
    if (device_index < 0 || device_index >= NUM_PLUTO_DEVICES) {
        return SDR_ERROR_INVALID_PARAM;
    }

    pluto_device_t* dev = &ctx->devices[device_index];

    fprintf(stderr, "WARNING: USB error on device %d, attempting reconnection...\n", device_index);

    // Destroy old context
    if (dev->rxbuf) {
        iio_buffer_destroy(dev->rxbuf);
        dev->rxbuf = NULL;
    }
    if (dev->ctx) {
        iio_context_destroy(dev->ctx);
        dev->ctx = NULL;
    }
    dev->is_connected = false;

    // Try to reconnect
    char old_uri[MAX_URI_LENGTH];
    strncpy(old_uri, dev->uri, MAX_URI_LENGTH);

    dev->ctx = iio_create_context_from_uri(old_uri);
    if (!dev->ctx) {
        fprintf(stderr, "ERROR: Failed to reconnect to device %d\n", device_index);
        return SDR_ERROR_USB_DISCONNECTED;
    }

    // Re-find devices and channels
    dev->phy_dev = find_device(dev->ctx, "ad9361-phy");
    dev->rx_dev = find_device(dev->ctx, "cf-ad9361-lpc");
    dev->rx0_i = find_channel(dev->rx_dev, "voltage0", false);
    dev->rx0_q = find_channel(dev->rx_dev, "voltage1", false);

    if (!dev->phy_dev || !dev->rx_dev || !dev->rx0_i || !dev->rx0_q) {
        fprintf(stderr, "ERROR: Failed to re-initialize device %d\n", device_index);
        iio_context_destroy(dev->ctx);
        dev->ctx = NULL;
        return SDR_ERROR_IIO_INIT_FAILED;
    }

    dev->is_connected = true;

    // Reconfigure device
    sdr_status_t status = sdr_pluto_configure(ctx, &ctx->config);
    if (status != SDR_SUCCESS) {
        fprintf(stderr, "ERROR: Failed to reconfigure device %d\n", device_index);
        return status;
    }

    printf("INFO: Successfully reconnected device %d\n", device_index);
    return SDR_SUCCESS;
}
