/**
 * @file test_aoa_triangulation.c
 * @brief Comprehensive unit tests for AOA triangulation module
 *
 * Tests cover:
 *   1. Two-bearing intersection (analytic)
 *   2. Multi-bearing WLS solver
 *   3. GDOP computation
 *   4. Bearing computation
 *   5. ENU <-> WGS84 roundtrip
 *   6. Parallel bearing rejection
 *   7. aoa_triangulate() dispatcher
 *   8. Error ellipse properties
 *
 * All test data is synthetic. No external data files required.
 */

#include "aoa_triangulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Test Framework
 * ============================================================================ */

static int test_passed = 0;
static int test_failed = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        if (condition) { \
            test_passed++; \
            printf("  PASS: %s\n", message); \
        } else { \
            test_failed++; \
            printf("  FAIL: %s\n", message); \
        } \
    } while(0)

/* ============================================================================
 * Test 1: Two-Bearing Intersection
 * ============================================================================ */

static void test_two_bearing_intersection(void) {
    printf("\n--- Test 1: Two-Bearing Intersection ---\n");

    /* Observer 1 at origin, Observer 2 at (100, 0).
     * Target at (50, 50).
     *
     * From (0,0) to (50,50):
     *   dx=50 (East), dy=50 (North)
     *   bearing = atan2(50, 50) = 45 degrees (NE)
     *
     * From (100,0) to (50,50):
     *   dx=-50 (West), dy=50 (North)
     *   bearing = atan2(-50, 50) = -45 => 315 degrees (NW)
     */
    aoa_bearing_t b1 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 45.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_bearing_t b2 = {
        .observer_x = 100.0,
        .observer_y = 0.0,
        .azimuth_deg = 315.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_result_t result;
    int ret = aoa_intersect_two(&b1, &b2, &result);

    TEST_ASSERT(ret == 0, "intersect_two returns success");
    TEST_ASSERT(result.success == true, "result.success is true");
    TEST_ASSERT(fabs(result.est_x - 50.0) < 1.0,
                "est_x within 1m of 50.0");
    TEST_ASSERT(fabs(result.est_y - 50.0) < 1.0,
                "est_y within 1m of 50.0");
    TEST_ASSERT(result.num_bearings_used == 2,
                "num_bearings_used == 2");

    printf("  (est_x=%.2f, est_y=%.2f)\n", result.est_x, result.est_y);

    /* Parallel bearings should fail */
    aoa_bearing_t b3 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 90.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_bearing_t b4 = {
        .observer_x = 0.0,
        .observer_y = 100.0,
        .azimuth_deg = 90.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    ret = aoa_intersect_two(&b3, &b4, &result);
    TEST_ASSERT(ret != 0, "parallel bearings (both 90 deg) return failure");
    TEST_ASSERT(result.success == false,
                "parallel bearings: result.success is false");
}

/* ============================================================================
 * Test 2: WLS Multi-Bearing Solver
 * ============================================================================ */

static void test_wls_solve(void) {
    printf("\n--- Test 2: WLS Multi-Bearing Solver ---\n");

    /* 4 observers at corners of a 100m square. Target at (50, 50).
     * Compute true bearings, add small noise (up to +/-2 degrees). */
    double target_x = 50.0, target_y = 50.0;

    double obs_positions[4][2] = {
        {  0.0,   0.0},
        {100.0,   0.0},
        {100.0, 100.0},
        {  0.0, 100.0}
    };

    /* Fixed noise offsets within +/-2 degrees for reproducibility */
    double noise_deg[4] = { 1.2, -0.8, 1.5, -1.1 };

    aoa_bearing_t bearings[4];
    for (int i = 0; i < 4; i++) {
        double true_bearing = aoa_compute_bearing(
            obs_positions[i][0], obs_positions[i][1],
            target_x, target_y);

        bearings[i].observer_x = obs_positions[i][0];
        bearings[i].observer_y = obs_positions[i][1];
        bearings[i].azimuth_deg = true_bearing + noise_deg[i];
        bearings[i].confidence = 1.0;
        bearings[i].timestamp = (double)i;
    }

    aoa_result_t result;
    int ret = aoa_wls_solve(bearings, 4, &result);

    TEST_ASSERT(ret == 0, "wls_solve returns success");
    TEST_ASSERT(result.success == true, "result.success is true");

    double error = sqrt((result.est_x - target_x) * (result.est_x - target_x) +
                        (result.est_y - target_y) * (result.est_y - target_y));

    TEST_ASSERT(error < 5.0,
                "WLS position error < 5m with +/-2 deg noise");
    TEST_ASSERT(result.num_bearings_used == 4,
                "num_bearings_used == 4");
    TEST_ASSERT(result.residual_deg >= 0.0,
                "residual_deg >= 0");

    printf("  (est_x=%.2f, est_y=%.2f, error=%.2fm, residual=%.2f deg)\n",
           result.est_x, result.est_y, error, result.residual_deg);
}

/* ============================================================================
 * Test 3: GDOP Computation
 * ============================================================================ */

static void test_gdop_computation(void) {
    printf("\n--- Test 3: GDOP Computation ---\n");

    double target_x = 50.0, target_y = 50.0;

    /* Good geometry: 4 observers surrounding target at 100m square corners */
    aoa_bearing_t good_geo[4];
    double good_obs[4][2] = {
        {  0.0,   0.0},
        {100.0,   0.0},
        {100.0, 100.0},
        {  0.0, 100.0}
    };

    for (int i = 0; i < 4; i++) {
        good_geo[i].observer_x = good_obs[i][0];
        good_geo[i].observer_y = good_obs[i][1];
        good_geo[i].azimuth_deg = aoa_compute_bearing(
            good_obs[i][0], good_obs[i][1], target_x, target_y);
        good_geo[i].confidence = 1.0;
        good_geo[i].timestamp = 0.0;
    }

    double gdop_good = aoa_compute_gdop(good_geo, 4, target_x, target_y);
    TEST_ASSERT(gdop_good > 0.0, "good geometry GDOP > 0");
    TEST_ASSERT(gdop_good < 5.0, "good geometry GDOP < 5 (well-surrounded)");
    printf("  Good geometry GDOP = %.2f\n", gdop_good);

    /* Poor geometry: all observers on one side (collinear, same y) */
    aoa_bearing_t poor_geo[4];
    double poor_obs[4][2] = {
        {  0.0, -100.0},
        { 10.0, -100.0},
        { 20.0, -100.0},
        { 30.0, -100.0}
    };

    for (int i = 0; i < 4; i++) {
        poor_geo[i].observer_x = poor_obs[i][0];
        poor_geo[i].observer_y = poor_obs[i][1];
        poor_geo[i].azimuth_deg = aoa_compute_bearing(
            poor_obs[i][0], poor_obs[i][1], target_x, target_y);
        poor_geo[i].confidence = 1.0;
        poor_geo[i].timestamp = 0.0;
    }

    double gdop_poor = aoa_compute_gdop(poor_geo, 4, target_x, target_y);
    TEST_ASSERT(gdop_poor > 10.0,
                "poor geometry GDOP > 10 (all observers on same side)");
    printf("  Poor geometry GDOP = %.2f\n", gdop_poor);

    TEST_ASSERT(gdop_poor > gdop_good,
                "poor GDOP > good GDOP");
}

/* ============================================================================
 * Test 4: Bearing Computation
 * ============================================================================ */

static void test_compute_bearing(void) {
    printf("\n--- Test 4: Bearing Computation ---\n");

    /* Convention: 0 = North (+Y), 90 = East (+X), 180 = South, 270 = West */

    /* (0,0) -> (100, 0): due East = 90 degrees */
    double b_east = aoa_compute_bearing(0.0, 0.0, 100.0, 0.0);
    TEST_ASSERT(fabs(b_east - 90.0) < 0.01,
                "(0,0)->(100,0) = 90 deg (East)");
    printf("  East bearing = %.4f deg\n", b_east);

    /* (0,0) -> (0, 100): due North = 0 degrees */
    double b_north = aoa_compute_bearing(0.0, 0.0, 0.0, 100.0);
    TEST_ASSERT(fabs(b_north - 0.0) < 0.01,
                "(0,0)->(0,100) = 0 deg (North)");
    printf("  North bearing = %.4f deg\n", b_north);

    /* (0,0) -> (-100, 0): due West = 270 degrees */
    double b_west = aoa_compute_bearing(0.0, 0.0, -100.0, 0.0);
    TEST_ASSERT(fabs(b_west - 270.0) < 0.01,
                "(0,0)->(-100,0) = 270 deg (West)");
    printf("  West bearing = %.4f deg\n", b_west);

    /* (0,0) -> (0, -100): due South = 180 degrees */
    double b_south = aoa_compute_bearing(0.0, 0.0, 0.0, -100.0);
    TEST_ASSERT(fabs(b_south - 180.0) < 0.01,
                "(0,0)->(0,-100) = 180 deg (South)");
    printf("  South bearing = %.4f deg\n", b_south);

    /* (0,0) -> (100, 100): NE = 45 degrees */
    double b_ne = aoa_compute_bearing(0.0, 0.0, 100.0, 100.0);
    TEST_ASSERT(fabs(b_ne - 45.0) < 0.01,
                "(0,0)->(100,100) = 45 deg (NE)");
    printf("  NE bearing = %.4f deg\n", b_ne);
}

/* ============================================================================
 * Test 5: ENU <-> WGS84 Roundtrip
 * ============================================================================ */

static void test_enu_wgs84_roundtrip(void) {
    printf("\n--- Test 5: ENU <-> WGS84 Roundtrip ---\n");

    /* Set origin at WPI campus: 42.2746 N, 71.8063 W */
    enu_origin_t origin;
    enu_set_origin(&origin, 42.2746, -71.8063, 200.0);

    /* Original ENU coordinates */
    double east_in = 100.0;
    double north_in = 200.0;
    double up_in = 50.0;

    /* ENU -> WGS84 */
    wgs84_coord_t wgs84;
    enu_to_wgs84(&origin, east_in, north_in, up_in, &wgs84);

    printf("  WGS84: lat=%.8f, lon=%.8f, alt=%.2f\n",
           wgs84.latitude_deg, wgs84.longitude_deg, wgs84.altitude_m);

    /* Sanity check on WGS84: latitude should have increased (moved North),
     * longitude should have increased (moved East, since negative longitude
     * becomes less negative) */
    TEST_ASSERT(wgs84.latitude_deg > origin.ref_lat_deg,
                "latitude increased (moved North)");
    TEST_ASSERT(wgs84.longitude_deg > origin.ref_lon_deg,
                "longitude increased (moved East)");
    TEST_ASSERT(fabs(wgs84.altitude_m - 250.0) < 0.01,
                "altitude = ref_alt + up = 250m");

    /* WGS84 -> ENU (roundtrip) */
    double east_out, north_out, up_out;
    wgs84_to_enu(&origin, &wgs84, &east_out, &north_out, &up_out);

    printf("  Roundtrip: east=%.6f, north=%.6f, up=%.6f\n",
           east_out, north_out, up_out);

    double err_east = fabs(east_out - east_in);
    double err_north = fabs(north_out - north_in);
    double err_up = fabs(up_out - up_in);

    TEST_ASSERT(err_east < 0.01,
                "East roundtrip error < 0.01m");
    TEST_ASSERT(err_north < 0.01,
                "North roundtrip error < 0.01m");
    TEST_ASSERT(err_up < 0.01,
                "Up roundtrip error < 0.01m");

    printf("  Errors: dE=%.6fm, dN=%.6fm, dU=%.6fm\n",
           err_east, err_north, err_up);
}

/* ============================================================================
 * Test 6: Parallel Bearing Rejection
 * ============================================================================ */

static void test_parallel_bearing_rejection(void) {
    printf("\n--- Test 6: Parallel Bearing Rejection ---\n");

    /* Two observers with nearly identical bearing direction.
     * Angular separation < AOA_PARALLEL_THRESHOLD_DEG (5 degrees).
     * This should be rejected by aoa_intersect_two(). */

    /* Observer 1 at (0, 0) bearing 45 degrees */
    aoa_bearing_t b1 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 45.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    /* Observer 2 at (10, 0) bearing 47 degrees (only 2 degrees different) */
    aoa_bearing_t b2 = {
        .observer_x = 10.0,
        .observer_y = 0.0,
        .azimuth_deg = 47.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_result_t result;
    int ret = aoa_intersect_two(&b1, &b2, &result);

    TEST_ASSERT(ret != 0,
                "nearly parallel bearings (2 deg diff) return failure");
    TEST_ASSERT(result.success == false,
                "nearly parallel: result.success is false");

    /* Exactly parallel (same azimuth) */
    aoa_bearing_t b3 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 180.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_bearing_t b4 = {
        .observer_x = 50.0,
        .observer_y = 0.0,
        .azimuth_deg = 180.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    ret = aoa_intersect_two(&b3, &b4, &result);
    TEST_ASSERT(ret != 0,
                "exactly parallel bearings (both 180 deg) return failure");

    /* Anti-parallel (0 and 180 degrees) should also fail (cross product ~ 0) */
    aoa_bearing_t b5 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 0.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_bearing_t b6 = {
        .observer_x = 50.0,
        .observer_y = 0.0,
        .azimuth_deg = 180.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    ret = aoa_intersect_two(&b5, &b6, &result);
    TEST_ASSERT(ret != 0,
                "anti-parallel bearings (0 and 180 deg) return failure");

    /* Just barely above threshold should succeed */
    aoa_bearing_t b7 = {
        .observer_x = 0.0,
        .observer_y = 0.0,
        .azimuth_deg = 0.0,
        .confidence = 1.0,
        .timestamp = 0.0
    };

    aoa_bearing_t b8 = {
        .observer_x = 100.0,
        .observer_y = 0.0,
        .azimuth_deg = 350.0,  /* 10 degrees difference from 0 => above 5 deg threshold */
        .confidence = 1.0,
        .timestamp = 0.0
    };

    ret = aoa_intersect_two(&b7, &b8, &result);
    TEST_ASSERT(ret == 0,
                "10 deg separation (above threshold) returns success");
}

/* ============================================================================
 * Test 7: aoa_triangulate() Dispatcher
 * ============================================================================ */

static void test_triangulate_dispatcher(void) {
    printf("\n--- Test 7: aoa_triangulate() Dispatcher ---\n");

    double target_x = 50.0, target_y = 50.0;

    /* 2 bearings: should call intersect_two, num_bearings_used == 2 */
    aoa_bearing_t two_bearings[2] = {
        { .observer_x =   0.0, .observer_y = 0.0,
          .azimuth_deg = 45.0, .confidence = 1.0, .timestamp = 0.0 },
        { .observer_x = 100.0, .observer_y = 0.0,
          .azimuth_deg = 315.0, .confidence = 1.0, .timestamp = 0.0 }
    };

    aoa_result_t result2;
    int ret2 = aoa_triangulate(two_bearings, 2, &result2);

    TEST_ASSERT(ret2 == 0,
                "triangulate with 2 bearings returns success");
    TEST_ASSERT(result2.num_bearings_used == 2,
                "2 bearings: num_bearings_used == 2");
    TEST_ASSERT(result2.success == true,
                "2 bearings: result.success is true");

    printf("  2 bearings: est=(%.2f, %.2f)\n", result2.est_x, result2.est_y);

    /* 4 bearings: should call wls_solve, num_bearings_used == 4 */
    aoa_bearing_t four_bearings[4];
    double obs4[4][2] = {
        {  0.0,   0.0},
        {100.0,   0.0},
        {100.0, 100.0},
        {  0.0, 100.0}
    };

    for (int i = 0; i < 4; i++) {
        four_bearings[i].observer_x = obs4[i][0];
        four_bearings[i].observer_y = obs4[i][1];
        four_bearings[i].azimuth_deg = aoa_compute_bearing(
            obs4[i][0], obs4[i][1], target_x, target_y);
        four_bearings[i].confidence = 1.0;
        four_bearings[i].timestamp = (double)i;
    }

    aoa_result_t result4;
    int ret4 = aoa_triangulate(four_bearings, 4, &result4);

    TEST_ASSERT(ret4 == 0,
                "triangulate with 4 bearings returns success");
    TEST_ASSERT(result4.num_bearings_used == 4,
                "4 bearings: num_bearings_used == 4");
    TEST_ASSERT(result4.success == true,
                "4 bearings: result.success is true");

    double err4 = sqrt((result4.est_x - target_x) * (result4.est_x - target_x) +
                       (result4.est_y - target_y) * (result4.est_y - target_y));
    TEST_ASSERT(err4 < 1.0,
                "4 bearings (perfect data): position error < 1m");

    printf("  4 bearings: est=(%.2f, %.2f), error=%.4fm\n",
           result4.est_x, result4.est_y, err4);

    /* 3 bearings: should also call wls_solve, num_bearings_used == 3 */
    aoa_bearing_t three_bearings[3];
    for (int i = 0; i < 3; i++) {
        three_bearings[i] = four_bearings[i];
    }

    aoa_result_t result3;
    int ret3 = aoa_triangulate(three_bearings, 3, &result3);

    TEST_ASSERT(ret3 == 0,
                "triangulate with 3 bearings returns success");
    TEST_ASSERT(result3.num_bearings_used == 3,
                "3 bearings: num_bearings_used == 3");

    printf("  3 bearings: est=(%.2f, %.2f)\n", result3.est_x, result3.est_y);

    /* 1 bearing: should fail (below AOA_MIN_BEARINGS) */
    aoa_bearing_t one_bearing[1] = {
        { .observer_x = 0.0, .observer_y = 0.0,
          .azimuth_deg = 45.0, .confidence = 1.0, .timestamp = 0.0 }
    };

    aoa_result_t result1;
    int ret1 = aoa_triangulate(one_bearing, 1, &result1);

    TEST_ASSERT(ret1 != 0,
                "triangulate with 1 bearing returns failure");

    /* NULL inputs: should fail gracefully */
    int ret_null = aoa_triangulate(NULL, 2, &result1);
    TEST_ASSERT(ret_null != 0,
                "triangulate with NULL bearings returns failure");
}

/* ============================================================================
 * Test 8: Error Ellipse Properties
 * ============================================================================ */

static void test_error_ellipse(void) {
    printf("\n--- Test 8: Error Ellipse Properties ---\n");

    /* Run WLS with known geometry and verify error ellipse constraints:
     *   - semi-major axis a > 0
     *   - semi-minor axis b > 0
     *   - a >= b (by definition) */

    double target_x = 50.0, target_y = 50.0;
    double obs[4][2] = {
        {  0.0,   0.0},
        {100.0,   0.0},
        {100.0, 100.0},
        {  0.0, 100.0}
    };

    /* Add small noise so the error ellipse is meaningful */
    double noise[4] = { 0.5, -0.3, 0.7, -0.4 };

    aoa_bearing_t bearings[4];
    for (int i = 0; i < 4; i++) {
        bearings[i].observer_x = obs[i][0];
        bearings[i].observer_y = obs[i][1];
        bearings[i].azimuth_deg = aoa_compute_bearing(
            obs[i][0], obs[i][1], target_x, target_y) + noise[i];
        bearings[i].confidence = 1.0;
        bearings[i].timestamp = (double)i;
    }

    aoa_result_t result;
    int ret = aoa_wls_solve(bearings, 4, &result);

    TEST_ASSERT(ret == 0, "WLS solve returns success");
    TEST_ASSERT(result.success == true, "result.success is true");

    TEST_ASSERT(result.error_ellipse_a > 0.0,
                "error_ellipse_a > 0 (semi-major axis positive)");
    TEST_ASSERT(result.error_ellipse_b > 0.0,
                "error_ellipse_b > 0 (semi-minor axis positive)");
    TEST_ASSERT(result.error_ellipse_a >= result.error_ellipse_b,
                "error_ellipse_a >= error_ellipse_b (a is semi-major)");

    printf("  Error ellipse: a=%.4f m, b=%.4f m, theta=%.2f deg\n",
           result.error_ellipse_a, result.error_ellipse_b,
           result.error_ellipse_theta);

    /* GDOP should be populated and positive */
    TEST_ASSERT(result.gdop > 0.0,
                "GDOP > 0 in WLS result");

    printf("  GDOP=%.4f, residual=%.4f deg\n",
           result.gdop, result.residual_deg);

    /* With asymmetric observer placement, the error ellipse should
     * be elongated (a significantly larger than b).
     * Place all observers on one side for an elongated ellipse. */
    double asym_obs[4][2] = {
        {-100.0, 0.0},
        { -80.0, 0.0},
        { -60.0, 0.0},
        { -40.0, 0.0}
    };

    aoa_bearing_t asym_bearings[4];
    for (int i = 0; i < 4; i++) {
        asym_bearings[i].observer_x = asym_obs[i][0];
        asym_bearings[i].observer_y = asym_obs[i][1];
        asym_bearings[i].azimuth_deg = aoa_compute_bearing(
            asym_obs[i][0], asym_obs[i][1], target_x, target_y) + noise[i];
        asym_bearings[i].confidence = 1.0;
        asym_bearings[i].timestamp = (double)i;
    }

    aoa_result_t asym_result;
    ret = aoa_wls_solve(asym_bearings, 4, &asym_result);

    if (ret == 0 && asym_result.success) {
        TEST_ASSERT(asym_result.error_ellipse_a > result.error_ellipse_a,
                    "asymmetric geometry has larger semi-major axis");
        printf("  Asymmetric ellipse: a=%.4f m, b=%.4f m\n",
               asym_result.error_ellipse_a, asym_result.error_ellipse_b);
    } else {
        printf("  (asymmetric WLS did not converge - skipping comparison)\n");
    }
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("==============================================\n");
    printf("  AOA Triangulation Unit Tests\n");
    printf("==============================================\n");

    test_two_bearing_intersection();
    test_wls_solve();
    test_gdop_computation();
    test_compute_bearing();
    test_enu_wgs84_roundtrip();
    test_parallel_bearing_rejection();
    test_triangulate_dispatcher();
    test_error_ellipse();

    printf("\n==============================================\n");
    printf("  Summary: %d passed, %d failed, %d total\n",
           test_passed, test_failed, test_passed + test_failed);
    printf("==============================================\n");

    return (test_failed > 0) ? 1 : 0;
}
