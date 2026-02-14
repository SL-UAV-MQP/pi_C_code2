/**
 * @file aoa_triangulation.c
 * @brief AOA-Based Triangulation Solver Implementation
 *
 * Implements bearing→position estimation using line intersection
 * and Weighted Least Squares (WLS) optimization.
 */

#include "aoa_triangulation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD(d) ((d) * M_PI / 180.0)
#define RAD2DEG(r) ((r) * 180.0 / M_PI)

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

/**
 * @brief Normalize angle to [0, 360) degrees
 */
static double normalize_angle_360(double deg) {
    deg = fmod(deg, 360.0);
    if (deg < 0.0) deg += 360.0;
    return deg;
}

/**
 * @brief Normalize angle difference to [-180, 180) degrees
 */
static double angle_diff(double a, double b) {
    double d = fmod(a - b + 180.0, 360.0);
    if (d < 0.0) d += 360.0;
    return d - 180.0;
}

/* ============================================================================
 * Bearing Computation
 * ============================================================================ */

double aoa_compute_bearing(double obs_x, double obs_y,
                           double tgt_x, double tgt_y) {
    double dx = tgt_x - obs_x;  /* East */
    double dy = tgt_y - obs_y;  /* North */

    /* atan2(East, North) gives bearing from North, CW positive */
    double bearing = RAD2DEG(atan2(dx, dy));
    return normalize_angle_360(bearing);
}

/* ============================================================================
 * Two-Bearing Intersection
 * ============================================================================ */

int aoa_intersect_two(const aoa_bearing_t* b1, const aoa_bearing_t* b2,
                      aoa_result_t* result) {
    if (!b1 || !b2 || !result) return -1;

    memset(result, 0, sizeof(aoa_result_t));

    /* Convert azimuths to unit direction vectors (North=+Y, East=+X) */
    double az1 = DEG2RAD(b1->azimuth_deg);
    double az2 = DEG2RAD(b2->azimuth_deg);

    /* Direction vectors: azimuth measured CW from North */
    double dx1 = sin(az1), dy1 = cos(az1);
    double dx2 = sin(az2), dy2 = cos(az2);

    /* Check for parallel bearings */
    double cross = dx1 * dy2 - dy1 * dx2;
    if (fabs(cross) < sin(DEG2RAD(AOA_PARALLEL_THRESHOLD_DEG))) {
        result->success = false;
        return -1;
    }

    /* Line intersection:
     * P1 + t1 * d1 = P2 + t2 * d2
     *
     * Solving for t1:
     * t1 = ((P2 - P1) × d2) / (d1 × d2)
     *
     * where × denotes 2D cross product: a × b = ax*by - ay*bx
     */
    double dpx = b2->observer_x - b1->observer_x;
    double dpy = b2->observer_y - b1->observer_y;

    double t1 = (dpx * dy2 - dpy * dx2) / cross;

    result->est_x = b1->observer_x + t1 * dx1;
    result->est_y = b1->observer_y + t1 * dy1;
    result->success = true;
    result->num_bearings_used = 2;

    /* Check if intersection is in front of both observers (t > 0) */
    double t2 = (dpx * dy1 - dpy * dx1) / cross;
    if (t1 < 0.0 || t2 < 0.0) {
        /* Intersection behind observer(s) - still report but flag */
        result->confidence = 0.0;
    }

    /* Simple error estimate based on intersection angle */
    double angle_sep = fabs(angle_diff(b1->azimuth_deg, b2->azimuth_deg));
    result->gdop = 1.0 / sin(DEG2RAD(fmin(angle_sep, 180.0 - angle_sep)));

    return 0;
}

/* ============================================================================
 * Multi-Bearing WLS Solver
 * ============================================================================ */

int aoa_wls_solve(const aoa_bearing_t* bearings, int num_bearings,
                  aoa_result_t* result) {
    if (!bearings || !result || num_bearings < 2) return -1;

    memset(result, 0, sizeof(aoa_result_t));

    /* ---- Initial estimate from best pair intersection ---- */
    double best_sep = 0.0;
    int best_i = 0, best_j = 1;

    for (int i = 0; i < num_bearings; i++) {
        for (int j = i + 1; j < num_bearings; j++) {
            double sep = fabs(angle_diff(bearings[i].azimuth_deg,
                                          bearings[j].azimuth_deg));
            /* Best pair is closest to 90° separation */
            double quality = fabs(sep - 90.0);
            if (i == 0 && j == 1) {
                best_sep = quality;
            } else if (quality < best_sep) {
                best_sep = quality;
                best_i = i;
                best_j = j;
            }
        }
    }

    aoa_result_t init_result;
    if (aoa_intersect_two(&bearings[best_i], &bearings[best_j], &init_result) != 0) {
        result->success = false;
        return -1;
    }

    double est_x = init_result.est_x;
    double est_y = init_result.est_y;

    /* ---- Gauss-Newton WLS iteration ---- */
    /* Minimize: Σ w_i * (θ_measured_i - atan2(x-xi, y-yi))²
     *
     * Linearize around current estimate:
     *   θ_pred = atan2(est_x - obs_x, est_y - obs_y)
     *   dθ/dx = (est_y - obs_y) / r²
     *   dθ/dy = -(est_x - obs_x) / r²
     *   where r² = (est_x-obs_x)² + (est_y-obs_y)²
     *
     * Normal equations: (J^T W J) Δp = J^T W Δθ
     */
    const int MAX_ITER = 20;
    const double CONVERGE_M = 0.01;  /* 1cm convergence threshold */

    for (int iter = 0; iter < MAX_ITER; iter++) {
        /* Build normal equations (2x2 system) */
        double JtWJ[4] = {0};   /* [a b; c d] flattened */
        double JtWr[2] = {0};   /* right-hand side */
        double total_residual = 0.0;

        for (int i = 0; i < num_bearings; i++) {
            double dx = est_x - bearings[i].observer_x;
            double dy = est_y - bearings[i].observer_y;
            double r2 = dx * dx + dy * dy;

            if (r2 < 1.0) r2 = 1.0;  /* Guard against observer-on-target */

            /* Predicted bearing */
            double pred_deg = aoa_compute_bearing(
                bearings[i].observer_x, bearings[i].observer_y, est_x, est_y);

            /* Bearing residual (measured - predicted) */
            double residual = DEG2RAD(angle_diff(bearings[i].azimuth_deg, pred_deg));

            /* Jacobian: dθ/d(est_x), dθ/d(est_y)
             * θ = atan2(dx, dy), where dx=East, dy=North
             * dθ/d(est_x) = dy / r²  (in radians)
             * dθ/d(est_y) = -dx / r²
             */
            double j_x = dy / r2;
            double j_y = -dx / r2;

            /* Weight */
            double w = bearings[i].confidence;
            if (w <= 0.0) w = 1.0;

            /* Accumulate normal equations */
            JtWJ[0] += w * j_x * j_x;
            JtWJ[1] += w * j_x * j_y;
            JtWJ[2] += w * j_y * j_x;
            JtWJ[3] += w * j_y * j_y;

            JtWr[0] += w * j_x * residual;
            JtWr[1] += w * j_y * residual;

            total_residual += fabs(RAD2DEG(residual));
        }

        /* Solve 2x2 system: JtWJ * delta = JtWr */
        double det = JtWJ[0] * JtWJ[3] - JtWJ[1] * JtWJ[2];
        if (fabs(det) < 1e-20) {
            /* Singular - geometry too poor */
            result->success = false;
            return -1;
        }

        double delta_x = (JtWJ[3] * JtWr[0] - JtWJ[1] * JtWr[1]) / det;
        double delta_y = (-JtWJ[2] * JtWr[0] + JtWJ[0] * JtWr[1]) / det;

        est_x += delta_x;
        est_y += delta_y;

        /* Check convergence */
        double step = sqrt(delta_x * delta_x + delta_y * delta_y);
        if (step < CONVERGE_M) {
            break;
        }
    }

    /* ---- Compute final statistics ---- */
    result->est_x = est_x;
    result->est_y = est_y;
    result->success = true;
    result->num_bearings_used = num_bearings;

    /* Mean bearing residual */
    double sum_residual = 0.0;
    for (int i = 0; i < num_bearings; i++) {
        double pred = aoa_compute_bearing(
            bearings[i].observer_x, bearings[i].observer_y, est_x, est_y);
        sum_residual += fabs(angle_diff(bearings[i].azimuth_deg, pred));
    }
    result->residual_deg = sum_residual / num_bearings;

    /* GDOP */
    result->gdop = aoa_compute_gdop(bearings, num_bearings, est_x, est_y);

    /* Error ellipse from (JtWJ)^-1 covariance */
    /* Recompute JtWJ at final estimate */
    double C[4] = {0};
    for (int i = 0; i < num_bearings; i++) {
        double dx = est_x - bearings[i].observer_x;
        double dy = est_y - bearings[i].observer_y;
        double r2 = dx * dx + dy * dy;
        if (r2 < 1.0) r2 = 1.0;

        double j_x = dy / r2;
        double j_y = -dx / r2;
        double w = bearings[i].confidence;
        if (w <= 0.0) w = 1.0;

        C[0] += w * j_x * j_x;
        C[1] += w * j_x * j_y;
        C[2] += w * j_y * j_x;
        C[3] += w * j_y * j_y;
    }

    double det = C[0] * C[3] - C[1] * C[2];
    if (fabs(det) > 1e-20) {
        /* Covariance = (JtWJ)^-1 */
        double cov_xx = C[3] / det;
        double cov_yy = C[0] / det;
        double cov_xy = -C[1] / det;

        /* Error ellipse from eigenvalues of covariance */
        double trace = cov_xx + cov_yy;
        double discr = sqrt(fmax(0.0,
            (cov_xx - cov_yy) * (cov_xx - cov_yy) + 4.0 * cov_xy * cov_xy));
        double lambda1 = (trace + discr) / 2.0;
        double lambda2 = (trace - discr) / 2.0;

        result->error_ellipse_a = sqrt(fmax(0.0, lambda1));
        result->error_ellipse_b = sqrt(fmax(0.0, lambda2));
        result->error_ellipse_theta = 0.5 * RAD2DEG(atan2(2.0 * cov_xy,
                                                            cov_xx - cov_yy));
    }

    return 0;
}

/* ============================================================================
 * Main Triangulation API
 * ============================================================================ */

int aoa_triangulate(const aoa_bearing_t* bearings, int num_bearings,
                    aoa_result_t* result) {
    if (!bearings || !result) return -1;
    if (num_bearings < AOA_MIN_BEARINGS) {
        fprintf(stderr, "ERROR: Need at least %d bearings (got %d)\n",
                AOA_MIN_BEARINGS, num_bearings);
        return -1;
    }

    memset(result, 0, sizeof(aoa_result_t));

    if (num_bearings == 2) {
        return aoa_intersect_two(&bearings[0], &bearings[1], result);
    } else {
        return aoa_wls_solve(bearings, num_bearings, result);
    }
}

/* ============================================================================
 * GDOP Computation
 * ============================================================================ */

double aoa_compute_gdop(const aoa_bearing_t* bearings, int num_bearings,
                        double est_x, double est_y) {
    if (!bearings || num_bearings < 2) return -1.0;

    /* GDOP = sqrt(trace((H^T H)^-1))
     * where H is the Jacobian matrix of bearing measurements.
     *
     * H[i,0] = cos(θ_i) / r_i  (direction cosine to target)
     * H[i,1] = sin(θ_i) / r_i
     *
     * For AOA-only (no range):
     * H[i,0] = -sin(θ_i)  (unit vector perpendicular to bearing)
     * H[i,1] = cos(θ_i)
     *
     * Actually for 2D AOA: use bearing geometry matrix
     */
    double HtH[4] = {0};

    for (int i = 0; i < num_bearings; i++) {
        double dx = est_x - bearings[i].observer_x;
        double dy = est_y - bearings[i].observer_y;
        double r = sqrt(dx * dx + dy * dy);

        if (r < 1.0) continue;

        /* Bearing direction (unit vector from observer to target) */
        double ux = dx / r;
        double uy = dy / r;

        /* Perpendicular (unit normal to bearing line) */
        double nx = -uy;
        double ny = ux;

        /* H row: [nx/r, ny/r] — sensitivity of bearing to position */
        double hx = nx / r;
        double hy = ny / r;

        HtH[0] += hx * hx;
        HtH[1] += hx * hy;
        HtH[2] += hy * hx;
        HtH[3] += hy * hy;
    }

    /* GDOP = sqrt(trace((HtH)^-1)) */
    double det = HtH[0] * HtH[3] - HtH[1] * HtH[2];
    if (fabs(det) < 1e-20) return 999.0;

    double inv_trace = (HtH[3] + HtH[0]) / det;
    if (inv_trace < 0.0) return 999.0;

    return sqrt(inv_trace);
}

/* ============================================================================
 * Coordinate Conversion (ENU ↔ WGS84)
 * ============================================================================ */

void enu_set_origin(enu_origin_t* origin,
                    double lat_deg, double lon_deg, double alt_m) {
    if (!origin) return;
    origin->ref_lat_deg = lat_deg;
    origin->ref_lon_deg = lon_deg;
    origin->ref_alt_m = alt_m;
}

void enu_to_wgs84(const enu_origin_t* origin,
                  double east, double north, double up,
                  wgs84_coord_t* wgs84) {
    if (!origin || !wgs84) return;

    double lat0 = DEG2RAD(origin->ref_lat_deg);
    double lon0 = DEG2RAD(origin->ref_lon_deg);

    /* WGS84 ellipsoid radii of curvature */
    double e2 = 2.0 * WGS84_F - WGS84_F * WGS84_F;
    double sin_lat = sin(lat0);
    double N = WGS84_A / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double M = WGS84_A * (1.0 - e2) / pow(1.0 - e2 * sin_lat * sin_lat, 1.5);

    /* Approximate conversion (valid for small displacements) */
    double dlat = north / M;
    double dlon = east / (N * cos(lat0));

    wgs84->latitude_deg = origin->ref_lat_deg + RAD2DEG(dlat);
    wgs84->longitude_deg = origin->ref_lon_deg + RAD2DEG(dlon);
    wgs84->altitude_m = origin->ref_alt_m + up;
}

void wgs84_to_enu(const enu_origin_t* origin,
                  const wgs84_coord_t* wgs84,
                  double* east, double* north, double* up) {
    if (!origin || !wgs84 || !east || !north || !up) return;

    double lat0 = DEG2RAD(origin->ref_lat_deg);
    double lon0 = DEG2RAD(origin->ref_lon_deg);
    double lat = DEG2RAD(wgs84->latitude_deg);
    double lon = DEG2RAD(wgs84->longitude_deg);

    /* WGS84 ellipsoid radii of curvature */
    double e2 = 2.0 * WGS84_F - WGS84_F * WGS84_F;
    double sin_lat = sin(lat0);
    double N = WGS84_A / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double M = WGS84_A * (1.0 - e2) / pow(1.0 - e2 * sin_lat * sin_lat, 1.5);

    *north = (lat - lat0) * M;
    *east = (lon - lon0) * N * cos(lat0);
    *up = wgs84->altitude_m - origin->ref_alt_m;
}
