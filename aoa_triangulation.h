/**
 * @file aoa_triangulation.h
 * @brief AOA-Based Triangulation Solver
 *
 * Computes transmitter position from multiple Angle-of-Arrival (AOA)
 * measurements taken at known observer positions.
 *
 * Features:
 * - 2-bearing intersection (analytic)
 * - Multi-bearing WLS (Weighted Least Squares)
 * - GDOP (Geometric Dilution of Precision)
 * - ENU ↔ WGS84 (Lat/Lon) coordinate conversion
 *
 * Usage:
 *   1. Collect AOA measurements from multiple positions
 *   2. Call aoa_triangulate() to estimate transmitter position
 *   3. Use aoa_compute_gdop() to assess geometry quality
 *   4. Convert result to Lat/Lon with enu_to_wgs84()
 */

#ifndef AOA_TRIANGULATION_H
#define AOA_TRIANGULATION_H

#include "common_types.h"
#include <stdbool.h>

/* ============================================================================
 * Configuration
 * ============================================================================ */

/** Maximum number of bearing measurements */
#define AOA_MAX_BEARINGS 20

/** Minimum bearings required for triangulation */
#define AOA_MIN_BEARINGS 2

/** Maximum angular separation for "parallel" bearings (degrees) */
#define AOA_PARALLEL_THRESHOLD_DEG 5.0

/** Earth's semi-major axis (WGS84) in meters */
#define WGS84_A 6378137.0

/** Earth's flattening (WGS84) */
#define WGS84_F (1.0 / 298.257223563)

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/**
 * @brief Single AOA bearing measurement
 */
typedef struct {
    double observer_x;       /**< Observer X position (East, meters) in ENU */
    double observer_y;       /**< Observer Y position (North, meters) in ENU */
    double azimuth_deg;      /**< Measured azimuth (degrees, 0=North, CW positive) */
    double confidence;       /**< Measurement confidence [0-1] (used as weight) */
    double timestamp;        /**< Measurement timestamp (seconds) */
} aoa_bearing_t;

/**
 * @brief Triangulation result
 */
typedef struct {
    bool success;            /**< Whether triangulation converged */
    double est_x;            /**< Estimated X (East, meters) */
    double est_y;            /**< Estimated Y (North, meters) */
    double error_ellipse_a;  /**< Error ellipse semi-major axis (meters) */
    double error_ellipse_b;  /**< Error ellipse semi-minor axis (meters) */
    double error_ellipse_theta; /**< Error ellipse orientation (degrees) */
    double gdop;             /**< Geometric Dilution of Precision */
    double residual_deg;     /**< Mean bearing residual (degrees) */
    int num_bearings_used;   /**< Number of bearings used */
} aoa_result_t;

/**
 * @brief WGS84 geographic coordinate
 */
typedef struct {
    double latitude_deg;     /**< Latitude (degrees, positive = North) */
    double longitude_deg;    /**< Longitude (degrees, positive = East) */
    double altitude_m;       /**< Altitude above WGS84 ellipsoid (meters) */
} wgs84_coord_t;

/**
 * @brief ENU reference origin
 */
typedef struct {
    double ref_lat_deg;      /**< Reference latitude (degrees) */
    double ref_lon_deg;      /**< Reference longitude (degrees) */
    double ref_alt_m;        /**< Reference altitude (meters) */
} enu_origin_t;

/* ============================================================================
 * Triangulation Functions
 * ============================================================================ */

/**
 * @brief Triangulate transmitter position from multiple bearings
 *
 * For 2 bearings: uses analytic line intersection
 * For 3+ bearings: uses Weighted Least Squares (WLS)
 *
 * @param bearings Array of bearing measurements
 * @param num_bearings Number of bearings (>= 2)
 * @param result Output triangulation result
 * @return 0 on success, negative on failure
 */
int aoa_triangulate(const aoa_bearing_t* bearings, int num_bearings,
                    aoa_result_t* result);

/**
 * @brief Two-bearing intersection (analytic solution)
 *
 * Finds intersection of two bearing lines.
 * Returns failure if bearings are nearly parallel.
 *
 * @param b1 First bearing
 * @param b2 Second bearing
 * @param result Output result
 * @return 0 on success, -1 if parallel
 */
int aoa_intersect_two(const aoa_bearing_t* b1, const aoa_bearing_t* b2,
                      aoa_result_t* result);

/**
 * @brief Multi-bearing Weighted Least Squares solver
 *
 * Minimizes weighted sum of squared bearing residuals:
 *   min_p Σ w_i * (θ_measured_i - θ_predicted_i(p))²
 *
 * Uses iterative Gauss-Newton with initial estimate from
 * best pair intersection.
 *
 * @param bearings Array of bearing measurements
 * @param num_bearings Number of bearings (>= 3)
 * @param result Output result
 * @return 0 on success, negative on failure
 */
int aoa_wls_solve(const aoa_bearing_t* bearings, int num_bearings,
                  aoa_result_t* result);

/**
 * @brief Compute GDOP (Geometric Dilution of Precision)
 *
 * Lower GDOP = better geometry for position estimation.
 * GDOP < 2: excellent, 2-5: good, 5-10: moderate, >10: poor
 *
 * @param bearings Array of bearing measurements
 * @param num_bearings Number of bearings
 * @param est_x Estimated X position
 * @param est_y Estimated Y position
 * @return GDOP value (>0), or -1.0 on error
 */
double aoa_compute_gdop(const aoa_bearing_t* bearings, int num_bearings,
                        double est_x, double est_y);

/**
 * @brief Compute bearing from observer to target
 *
 * @param obs_x Observer X (East, meters)
 * @param obs_y Observer Y (North, meters)
 * @param tgt_x Target X (East, meters)
 * @param tgt_y Target Y (North, meters)
 * @return Bearing in degrees (0=North, CW positive, [0, 360))
 */
double aoa_compute_bearing(double obs_x, double obs_y,
                           double tgt_x, double tgt_y);

/* ============================================================================
 * Coordinate Conversion Functions
 * ============================================================================ */

/**
 * @brief Convert ENU (East-North-Up) to WGS84 (Lat/Lon/Alt)
 *
 * @param origin ENU reference origin
 * @param east East coordinate (meters)
 * @param north North coordinate (meters)
 * @param up Up coordinate (meters)
 * @param wgs84 Output WGS84 coordinate
 */
void enu_to_wgs84(const enu_origin_t* origin,
                  double east, double north, double up,
                  wgs84_coord_t* wgs84);

/**
 * @brief Convert WGS84 (Lat/Lon/Alt) to ENU (East-North-Up)
 *
 * @param origin ENU reference origin
 * @param wgs84 Input WGS84 coordinate
 * @param east Output East coordinate (meters)
 * @param north Output North coordinate (meters)
 * @param up Output Up coordinate (meters)
 */
void wgs84_to_enu(const enu_origin_t* origin,
                  const wgs84_coord_t* wgs84,
                  double* east, double* north, double* up);

/**
 * @brief Set ENU origin from WGS84 coordinate
 *
 * @param origin Output origin structure
 * @param lat_deg Latitude (degrees)
 * @param lon_deg Longitude (degrees)
 * @param alt_m Altitude (meters)
 */
void enu_set_origin(enu_origin_t* origin,
                    double lat_deg, double lon_deg, double alt_m);

#endif /* AOA_TRIANGULATION_H */
