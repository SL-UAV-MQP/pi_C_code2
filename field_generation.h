/**
 * @file field_generation.h
 * @brief Signal Strength Field Generation
 *
 * Generates predicted signal strength heatmaps for path planning.
 * Implements multiple path-loss propagation models (Free Space, Log-Distance,
 * Hata Urban, Two-Ray Ground) for RF coverage prediction.
 *
 * Features:
 * - 2D signal field computation on regular grid
 * - Multiple propagation models
 * - Signal gradient computation
 * - Coverage hole detection
 * - Multi-source field fusion
 *
 * Performance Target: <100ms for 300x300 grid
 */

#ifndef FIELD_GENERATION_H
#define FIELD_GENERATION_H

#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of signal sources */
#define MAX_SIGNAL_SOURCES 16

/** Maximum grid dimension (one axis) */
#define MAX_GRID_DIM 500

/** Default path loss exponent for urban environment */
#define DEFAULT_PATH_LOSS_EXPONENT 3.5

/** Speed of light (m/s) - guarded to avoid conflict with music_config.h */
#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792458.0
#endif

/* ============================================================================
 * Enums and Type Definitions
 * ============================================================================ */

/**
 * @brief Path loss propagation models
 */
typedef enum {
    PATH_LOSS_FREE_SPACE = 0,    /**< Free space path loss (FSPL) */
    PATH_LOSS_LOG_DISTANCE,      /**< Log-distance path loss model */
    PATH_LOSS_HATA_URBAN,        /**< Hata urban propagation model */
    PATH_LOSS_TWO_RAY            /**< Two-ray ground reflection model */
} path_loss_model_t;

/**
 * @brief City type for Hata model
 */
typedef enum {
    CITY_TYPE_SMALL = 0,
    CITY_TYPE_MEDIUM,
    CITY_TYPE_LARGE
} city_type_t;

/**
 * @brief 3D position vector
 */
typedef struct {
    float x;  /**< X coordinate (meters, East) */
    float y;  /**< Y coordinate (meters, North) */
    float z;  /**< Z coordinate (meters, Up) */
} point3d_t;

/**
 * @brief Signal source (cell tower or transmitter)
 */
typedef struct {
    point3d_t position;           /**< Source position in ENU coordinates */
    float transmit_power_dbm;     /**< Transmit power (dBm) */
    float frequency_mhz;          /**< Carrier frequency (MHz) */
    float antenna_gain_dbi;       /**< Antenna gain (dBi) */
    char source_id[32];           /**< Source identifier */
    bool active;                  /**< Whether source is active */
} signal_source_t;

/**
 * @brief 2D signal strength field
 */
typedef struct {
    int width;                    /**< Grid width (number of x cells) */
    int height;                   /**< Grid height (number of y cells) */
    float resolution;             /**< Grid resolution (meters per cell) */
    float x_min, x_max;          /**< X bounds (meters) */
    float y_min, y_max;          /**< Y bounds (meters) */
    float altitude;               /**< UAV altitude above ground (meters) */
    float* signal_strength;       /**< Signal strength grid (dBm) [height x width] */
    float* gradient_x;            /**< Signal gradient X component (dBm/m) */
    float* gradient_y;            /**< Signal gradient Y component (dBm/m) */
} signal_field_t;

/**
 * @brief Signal field generator context
 */
typedef struct {
    signal_field_t field;                        /**< Signal field data */
    signal_source_t sources[MAX_SIGNAL_SOURCES]; /**< Signal sources */
    int num_sources;                             /**< Number of active sources */
    path_loss_model_t model;                     /**< Path loss model */
    float path_loss_exponent;                    /**< Path loss exponent (log-distance) */
    float shadowing_std_db;                      /**< Log-normal shadowing std dev */

    /* Hata model parameters */
    float base_height_m;                         /**< Base station height (meters) */
    float mobile_height_m;                       /**< Mobile antenna height (meters) */
    city_type_t city_type;                       /**< City type for Hata model */
} field_generator_t;

/**
 * @brief Coverage hole detection result
 */
typedef struct {
    int num_holes;                /**< Number of coverage holes found */
    point3d_t* positions;         /**< Positions of coverage holes (x,y,0) */
    float* signal_strength;       /**< Signal strength at each hole (dBm) */
} coverage_holes_t;

/* ============================================================================
 * Core Functions
 * ============================================================================ */

/**
 * @brief Initialize field generator
 *
 * @param gen Pointer to generator structure (allocated by caller)
 * @param x_min Minimum X coordinate (meters)
 * @param x_max Maximum X coordinate (meters)
 * @param y_min Minimum Y coordinate (meters)
 * @param y_max Maximum Y coordinate (meters)
 * @param resolution Grid resolution (meters per cell)
 * @param altitude UAV altitude above ground (meters)
 * @param model Path loss model to use
 * @return 0 on success, negative error code on failure
 */
int field_generator_init(field_generator_t* gen,
                        float x_min, float x_max,
                        float y_min, float y_max,
                        float resolution,
                        float altitude,
                        path_loss_model_t model);

/**
 * @brief Free field generator resources
 *
 * @param gen Pointer to generator structure
 */
void field_generator_free(field_generator_t* gen);

/**
 * @brief Add signal source to field
 *
 * @param gen Field generator
 * @param source Signal source to add
 * @return 0 on success, negative error code on failure
 */
int field_generator_add_source(field_generator_t* gen,
                               const signal_source_t* source);

/**
 * @brief Remove all signal sources
 *
 * @param gen Field generator
 */
void field_generator_clear_sources(field_generator_t* gen);

/**
 * @brief Compute signal strength field from all sources
 *
 * Computes received signal strength at each grid point based on
 * transmitter positions, powers, and selected propagation model.
 *
 * Performance: <100ms for 300x300 grid on Raspberry Pi 5B
 *
 * @param gen Field generator
 * @return 0 on success, negative error code on failure
 */
int field_generator_compute(field_generator_t* gen);

/**
 * @brief Compute signal gradient field (for gradient-based navigation)
 *
 * Uses central differences to compute signal gradient.
 * Must be called after field_generator_compute().
 *
 * @param gen Field generator
 * @return 0 on success, negative error code on failure
 */
int field_generator_compute_gradient(field_generator_t* gen);

/**
 * @brief Get signal strength at a specific position (bilinear interpolation)
 *
 * @param gen Field generator
 * @param x X coordinate (meters)
 * @param y Y coordinate (meters)
 * @return Signal strength (dBm), or -150.0 if out of bounds
 */
float field_generator_get_signal(const field_generator_t* gen, float x, float y);

/**
 * @brief Get signal gradient at position (bilinear interpolation)
 *
 * @param gen Field generator
 * @param x X coordinate (meters)
 * @param y Y coordinate (meters)
 * @param grad_x Output: gradient X component (dBm/m)
 * @param grad_y Output: gradient Y component (dBm/m)
 * @return 0 on success, negative error code on failure
 */
int field_generator_get_gradient(const field_generator_t* gen,
                                 float x, float y,
                                 float* grad_x, float* grad_y);

/**
 * @brief Find coverage holes (regions with signal below threshold)
 *
 * @param gen Field generator
 * @param threshold_dbm Signal threshold (dBm)
 * @param result Output: coverage holes (caller must free)
 * @return 0 on success, negative error code on failure
 */
int field_generator_find_coverage_holes(const field_generator_t* gen,
                                        float threshold_dbm,
                                        coverage_holes_t* result);

/**
 * @brief Free coverage holes result
 *
 * @param holes Coverage holes structure
 */
void coverage_holes_free(coverage_holes_t* holes);

/* ============================================================================
 * Path Loss Model Functions
 * ============================================================================ */

/**
 * @brief Free space path loss (FSPL)
 *
 * FSPL(dB) = 20*log10(d_km) + 20*log10(f_MHz) + 32.45
 *
 * @param distance_m Distance (meters)
 * @param frequency_mhz Frequency (MHz)
 * @return Path loss (dB)
 */
float path_loss_free_space(float distance_m, float frequency_mhz);

/**
 * @brief Log-distance path loss model
 *
 * PL(d) = PL(d0) + 10*n*log10(d/d0)
 *
 * @param distance_m Distance (meters)
 * @param path_loss_exponent Path loss exponent (2-5, typical urban: 3.5)
 * @param reference_distance Reference distance (meters)
 * @param reference_loss_db Path loss at reference distance (dB)
 * @return Path loss (dB)
 */
float path_loss_log_distance(float distance_m,
                             float path_loss_exponent,
                             float reference_distance,
                             float reference_loss_db);

/**
 * @brief Hata urban propagation model
 *
 * Valid ranges:
 * - Distance: 1-20 km
 * - Frequency: 150-1500 MHz
 * - Base height: 30-200 m
 * - Mobile height: 1-10 m
 *
 * @param distance_m Distance (meters)
 * @param frequency_mhz Frequency (MHz)
 * @param base_height_m Base station height (meters)
 * @param mobile_height_m Mobile antenna height (meters)
 * @param city_type City type (small/medium/large)
 * @return Path loss (dB)
 */
float path_loss_hata_urban(float distance_m,
                           float frequency_mhz,
                           float base_height_m,
                           float mobile_height_m,
                           city_type_t city_type);

/**
 * @brief Two-ray ground reflection model
 *
 * For distances beyond crossover distance, uses two-ray model.
 * Otherwise falls back to free space model.
 *
 * @param distance_m Distance (meters)
 * @param tx_height_m Transmitter height (meters)
 * @param rx_height_m Receiver height (meters)
 * @param frequency_mhz Frequency (MHz)
 * @return Path loss (dB)
 */
float path_loss_two_ray_ground(float distance_m,
                               float tx_height_m,
                               float rx_height_m,
                               float frequency_mhz);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Create CMRCM test field with known towers
 *
 * Creates a field generator configured for the CMRCM test area
 * with known cell tower positions.
 *
 * @param gen Field generator (allocated by caller)
 * @param resolution Grid resolution (meters)
 * @param altitude UAV altitude (meters)
 * @return 0 on success, negative error code on failure
 */
int create_cmrcm_field(field_generator_t* gen,
                      float resolution,
                      float altitude);

/**
 * @brief Get field statistics
 *
 * @param gen Field generator
 * @param min_dbm Output: minimum signal strength
 * @param max_dbm Output: maximum signal strength
 * @param mean_dbm Output: mean signal strength
 * @return 0 on success, negative error code on failure
 */
int field_generator_get_stats(const field_generator_t* gen,
                              float* min_dbm,
                              float* max_dbm,
                              float* mean_dbm);

#endif /* FIELD_GENERATION_H */
