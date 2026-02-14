/**
 * @file sdr_pluto.h
 * @brief ADALM-PLUTO SDR Interface using LibIIO
 *
 * Complete LibIIO integration for controlling 3x ADALM-PLUTO SDRs
 * to capture synchronized IQ samples for MUSIC DOA estimation.
 *
 * Hardware:
 * - 3x ADALM-PLUTO SDR (AD9363 RF Transceiver)
 * - USB 2.0/3.0 interface
 * - Frequency range: 325 MHz - 3.8 GHz
 * - Sample rate: up to 61.44 MSPS
 *
 * Target configuration:
 * - Frequency: 700-1400 MHz (LTE/P25 bands)
 * - Sample rate: 20 MSPS (configurable)
 * - 6 channels total (2 channels per PLUTO: RX1_I/Q)
 * - Synchronized capture across all 3 devices
 */

#ifndef SDR_PLUTO_H
#define SDR_PLUTO_H

#include "common_types.h"
#include <iio.h>
#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Number of PLUTO SDR devices */
#define NUM_PLUTO_DEVICES 3

/** Total number of receive channels (2 per PLUTO for I/Q) */
#define NUM_RX_CHANNELS 6

/** Maximum sample rate in samples per second */
#define MAX_SAMPLE_RATE_SPS 61440000

/** Minimum sample rate in samples per second */
#define MIN_SAMPLE_RATE_SPS 2084000

/** Default sample rate (20 MSPS for MUSIC processing) */
#define DEFAULT_SAMPLE_RATE_SPS 20000000

/** Default buffer size in samples (for covariance estimation) */
#define DEFAULT_BUFFER_SIZE 8192

/** Maximum buffer size in samples */
#define MAX_BUFFER_SIZE 16384

/** Default RF bandwidth (Hz) - typically 1.5x sample rate */
#define DEFAULT_RF_BANDWIDTH_HZ 30000000

/** Default RX gain mode: "manual", "slow_attack", "fast_attack", "hybrid" */
#define DEFAULT_GAIN_MODE "manual"

/** Default manual gain in dB (0 to 76 dB) */
#define DEFAULT_GAIN_DB 60

/** USB connection timeout in milliseconds */
#define USB_TIMEOUT_MS 5000

/** Maximum device URI length */
#define MAX_URI_LENGTH 256

/* ============================================================================
 * Gain Control Modes
 * ============================================================================ */

/** Gain control modes for PLUTO */
typedef enum {
    GAIN_MODE_MANUAL = 0,      /**< Manual gain control */
    GAIN_MODE_SLOW_ATTACK,     /**< AGC slow attack */
    GAIN_MODE_FAST_ATTACK,     /**< AGC fast attack */
    GAIN_MODE_HYBRID           /**< AGC hybrid mode */
} pluto_gain_mode_t;

/* ============================================================================
 * Error Codes
 * ============================================================================ */

/** SDR operation return codes */
typedef enum {
    SDR_SUCCESS = 0,                /**< Operation successful */
    SDR_ERROR_NULL_POINTER,         /**< Null pointer argument */
    SDR_ERROR_INVALID_PARAM,        /**< Invalid parameter */
    SDR_ERROR_DEVICE_NOT_FOUND,     /**< PLUTO device not found */
    SDR_ERROR_USB_DISCONNECTED,     /**< USB connection lost */
    SDR_ERROR_IIO_INIT_FAILED,      /**< LibIIO initialization failed */
    SDR_ERROR_CHANNEL_NOT_FOUND,    /**< IIO channel not found */
    SDR_ERROR_BUFFER_CREATE_FAILED, /**< Buffer creation failed */
    SDR_ERROR_CAPTURE_TIMEOUT,      /**< Capture timeout */
    SDR_ERROR_SYNC_FAILED,          /**< Synchronization failed */
    SDR_ERROR_OUT_OF_RANGE,         /**< Parameter out of valid range */
    SDR_ERROR_MEMORY                /**< Memory allocation failed */
} sdr_status_t;

/* ============================================================================
 * Configuration Structures
 * ============================================================================ */

/**
 * @brief SDR device configuration parameters
 */
typedef struct {
    uint64_t center_freq_hz;        /**< Center frequency in Hz */
    uint64_t sample_rate_sps;       /**< Sample rate in samples/sec */
    uint64_t rf_bandwidth_hz;       /**< RF bandwidth in Hz */
    pluto_gain_mode_t gain_mode;    /**< Gain control mode */
    int gain_db;                    /**< Manual gain in dB (0-76) */
    int buffer_size;                /**< Buffer size in samples */
    bool enable_quadrature_tracking; /**< Enable quadrature tracking */
    bool enable_rf_dc_offset_tracking; /**< Enable RF DC offset tracking */
    bool enable_bb_dc_offset_tracking; /**< Enable baseband DC offset tracking */
} sdr_config_t;

/**
 * @brief Single PLUTO SDR device context
 */
typedef struct {
    struct iio_context* ctx;        /**< IIO context */
    struct iio_device* phy_dev;     /**< AD9361 PHY device (config) */
    struct iio_device* rx_dev;      /**< RX streaming device */
    struct iio_channel* rx0_i;      /**< RX0 I channel */
    struct iio_channel* rx0_q;      /**< RX0 Q channel */
    struct iio_buffer* rxbuf;       /**< RX buffer */

    char uri[MAX_URI_LENGTH];       /**< Device URI (e.g., "usb:1.2.5") */
    int device_index;               /**< Device index (0-2) */
    bool is_connected;              /**< Connection status */
    uint64_t last_timestamp;        /**< Last sample timestamp */
} pluto_device_t;

/**
 * @brief Multi-SDR context for 3x PLUTO devices
 */
typedef struct {
    pluto_device_t devices[NUM_PLUTO_DEVICES]; /**< Array of PLUTO devices */
    sdr_config_t config;            /**< Shared configuration */
    int num_connected;              /**< Number of connected devices */
    bool is_streaming;              /**< Streaming status */
    uint64_t sync_timestamp;        /**< Synchronization timestamp */
} sdr_pluto_context_t;

/**
 * @brief Captured IQ data from all 6 channels
 */
typedef struct {
    int num_samples;                /**< Number of samples per channel */
    int num_channels;               /**< Number of channels (6) */
    cdouble_t* channel_data[NUM_RX_CHANNELS]; /**< IQ data per channel */
    uint64_t timestamp;             /**< Capture timestamp (nanoseconds) */
    double actual_sample_rate;      /**< Actual sample rate achieved */
} sdr_capture_result_t;

/* ============================================================================
 * Initialization and Configuration
 * ============================================================================ */

/**
 * @brief Initialize PLUTO SDR context
 *
 * Allocates and initializes the multi-SDR context structure.
 * Does not connect to devices yet - use sdr_pluto_connect().
 *
 * @return Pointer to allocated context or NULL on failure
 */
sdr_pluto_context_t* sdr_pluto_init(void);

/**
 * @brief Connect to all PLUTO devices
 *
 * Scans for available PLUTO devices via USB and establishes connections.
 * Attempts to connect to exactly 3 PLUTO SDRs.
 *
 * USB URIs are automatically detected:
 * - usb:1.2.5 (PLUTO #1)
 * - usb:1.3.5 (PLUTO #2)
 * - usb:1.4.5 (PLUTO #3)
 *
 * @param ctx SDR context
 * @return SDR_SUCCESS or error code
 *
 * @note Requires all 3 PLUTOs to be connected via USB
 * @note Devices are identified by USB bus/port numbers
 */
sdr_status_t sdr_pluto_connect(sdr_pluto_context_t* ctx);

/**
 * @brief Configure SDR parameters
 *
 * Sets center frequency, sample rate, gain, bandwidth, etc.
 * Applied to all connected PLUTO devices.
 *
 * @param ctx SDR context
 * @param config Configuration parameters
 * @return SDR_SUCCESS or error code
 *
 * @note Validates all parameters before applying
 * @note Automatically calculates RF bandwidth if set to 0
 */
sdr_status_t sdr_pluto_configure(
    sdr_pluto_context_t* ctx,
    const sdr_config_t* config
);

/**
 * @brief Get default configuration
 *
 * Returns a default configuration optimized for MUSIC DOA estimation:
 * - Center frequency: 1050 MHz
 * - Sample rate: 20 MSPS
 * - Gain: 60 dB manual
 * - Buffer size: 4000 samples
 *
 * @param config Output configuration structure
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_get_default_config(sdr_config_t* config);

/* ============================================================================
 * Data Capture
 * ============================================================================ */

/**
 * @brief Capture synchronized IQ samples from all 6 channels
 *
 * Captures N snapshots from all 3 PLUTOs simultaneously.
 * Returns complex IQ data organized by channel.
 *
 * Capture sequence:
 * 1. Synchronize device timestamps
 * 2. Trigger simultaneous capture
 * 3. Read IQ buffers from all devices
 * 4. Verify timestamp alignment (<1 sample drift)
 * 5. Convert to complex double format
 *
 * @param ctx SDR context
 * @param num_samples Number of samples to capture per channel
 * @param result Output capture result (allocated by caller)
 * @return SDR_SUCCESS or error code
 *
 * @note num_samples should match NUM_SNAPSHOTS (4000 for MUSIC)
 * @note result->channel_data must be pre-allocated
 * @note Capture time: ~200μs per snapshot at 20 MSPS
 */
sdr_status_t sdr_pluto_capture(
    sdr_pluto_context_t* ctx,
    int num_samples,
    sdr_capture_result_t* result
);

/**
 * @brief Convert captured IQ data to signal_matrix_t format
 *
 * Converts raw SDR capture to MUSIC-compatible signal matrix.
 * Output is M=6 antennas x N snapshots (column-major).
 *
 * @param capture Captured IQ data from SDR
 * @param signal_matrix Output signal matrix (allocated by caller)
 * @return SDR_SUCCESS or error code
 *
 * @note Output matrix must be pre-allocated with M=6, N=num_samples
 * @note Data is copied and converted to complex double
 */
sdr_status_t sdr_pluto_convert_to_signal_matrix(
    const sdr_capture_result_t* capture,
    signal_matrix_t* signal_matrix
);

/* ============================================================================
 * Synchronization
 * ============================================================================ */

/**
 * @brief Synchronize timestamps across all 3 PLUTOs
 *
 * Ensures all devices start capturing at the same time.
 * Uses hardware timestamp registers for sub-sample accuracy.
 *
 * Synchronization method:
 * 1. Read current timestamp from each device
 * 2. Calculate average timestamp
 * 3. Set trigger time to average + 1ms (allows setup time)
 * 4. Verify all devices triggered within tolerance
 *
 * @param ctx SDR context
 * @return SDR_SUCCESS or error code
 *
 * @note Target synchronization accuracy: <1 sample (50ns at 20 MSPS)
 * @note Must be called before each capture sequence
 */
sdr_status_t sdr_pluto_synchronize_devices(sdr_pluto_context_t* ctx);

/**
 * @brief Get timestamp from specific device
 *
 * Reads hardware timestamp counter from PLUTO.
 *
 * @param device PLUTO device
 * @param timestamp Output timestamp in nanoseconds
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_get_timestamp(
    pluto_device_t* device,
    uint64_t* timestamp
);

/* ============================================================================
 * Device Control
 * ============================================================================ */

/**
 * @brief Set center frequency for all devices
 *
 * @param ctx SDR context
 * @param freq_hz Center frequency in Hz (325 MHz - 3.8 GHz)
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_set_frequency(
    sdr_pluto_context_t* ctx,
    uint64_t freq_hz
);

/**
 * @brief Set sample rate for all devices
 *
 * @param ctx SDR context
 * @param sample_rate_sps Sample rate in samples/sec
 * @return SDR_SUCCESS or error code
 *
 * @note Valid range: 2.084 - 61.44 MSPS
 * @note RF bandwidth is automatically adjusted
 */
sdr_status_t sdr_pluto_set_sample_rate(
    sdr_pluto_context_t* ctx,
    uint64_t sample_rate_sps
);

/**
 * @brief Set gain for all devices
 *
 * @param ctx SDR context
 * @param gain_db Gain in dB (0-76 for manual mode)
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_set_gain(
    sdr_pluto_context_t* ctx,
    int gain_db
);

/**
 * @brief Set gain mode for all devices
 *
 * @param ctx SDR context
 * @param gain_mode Gain control mode
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_set_gain_mode(
    sdr_pluto_context_t* ctx,
    pluto_gain_mode_t gain_mode
);

/* ============================================================================
 * Status and Diagnostics
 * ============================================================================ */

/**
 * @brief Check if device is connected
 *
 * @param device PLUTO device
 * @return true if connected, false otherwise
 */
bool sdr_pluto_is_connected(const pluto_device_t* device);

/**
 * @brief Get device information
 *
 * Retrieves firmware version, serial number, etc.
 *
 * @param device PLUTO device
 * @param info_buffer Output buffer for device info string
 * @param buffer_size Size of info_buffer
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_get_device_info(
    const pluto_device_t* device,
    char* info_buffer,
    int buffer_size
);

/**
 * @brief Print SDR configuration
 *
 * Displays current settings for all devices.
 *
 * @param ctx SDR context
 */
void sdr_pluto_print_config(const sdr_pluto_context_t* ctx);

/**
 * @brief Print capture statistics
 *
 * Displays capture timing, sample rate, buffer fill, etc.
 *
 * @param result Capture result
 */
void sdr_pluto_print_capture_stats(const sdr_capture_result_t* result);

/* ============================================================================
 * Cleanup
 * ============================================================================ */

/**
 * @brief Disconnect all devices
 *
 * Stops streaming, destroys buffers, closes IIO contexts.
 *
 * @param ctx SDR context
 * @return SDR_SUCCESS or error code
 */
sdr_status_t sdr_pluto_disconnect(sdr_pluto_context_t* ctx);

/**
 * @brief Free SDR context
 *
 * Releases all resources. Calls sdr_pluto_disconnect() if needed.
 *
 * @param ctx SDR context to free
 */
void sdr_pluto_free(sdr_pluto_context_t* ctx);

/* ============================================================================
 * Memory Management Helpers
 * ============================================================================ */

/**
 * @brief Allocate capture result structure
 *
 * @param num_samples Number of samples per channel
 * @param num_channels Number of channels (6)
 * @return Pointer to allocated result or NULL on failure
 */
sdr_capture_result_t* sdr_capture_result_alloc(int num_samples, int num_channels);

/**
 * @brief Free capture result
 *
 * @param result Result to free
 */
void sdr_capture_result_free(sdr_capture_result_t* result);

/* ============================================================================
 * Error Handling
 * ============================================================================ */

/**
 * @brief Get error message string
 *
 * @param status Error code
 * @return Human-readable error message
 */
const char* sdr_pluto_get_error_string(sdr_status_t status);

/**
 * @brief Handle USB disconnection error
 *
 * Attempts to reconnect to device.
 *
 * @param ctx SDR context
 * @param device_index Index of disconnected device
 * @return SDR_SUCCESS if reconnected, error code otherwise
 */
sdr_status_t sdr_pluto_handle_usb_error(
    sdr_pluto_context_t* ctx,
    int device_index
);

#endif /* SDR_PLUTO_H */
