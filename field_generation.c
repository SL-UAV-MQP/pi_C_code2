/**
 * @file field_generation.c
 * @brief Signal Strength Field Generation
 *
 * Generates predicted signal strength heatmaps for path planning.
 * Implements multiple path-loss propagation models.
 */

#include "field_generation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Path Loss Model Functions
 * ============================================================================ */

float path_loss_free_space(float distance_m, float frequency_mhz) {
    if (distance_m < 1.0f) {
        distance_m = 1.0f;
    }

    // FSPL = 20*log10(d_km) + 20*log10(f_MHz) + 32.45
    float distance_km = distance_m / 1000.0f;
    float fspl = 20.0f * log10f(distance_km) + 20.0f * log10f(frequency_mhz) + 32.45f;
    return fspl;
}

float path_loss_log_distance(float distance_m,
                             float path_loss_exponent,
                             float reference_distance,
                             float reference_loss_db) {
    if (distance_m < reference_distance) {
        distance_m = reference_distance;
    }

    // PL(d) = PL(d0) + 10*n*log10(d/d0)
    float pl = reference_loss_db + 10.0f * path_loss_exponent * log10f(distance_m / reference_distance);
    return pl;
}

float path_loss_hata_urban(float distance_m,
                           float frequency_mhz,
                           float base_height_m,
                           float mobile_height_m,
                           city_type_t city_type) {
    float distance_km = distance_m / 1000.0f;

    // Clamp to valid ranges
    if (distance_km < 1.0f) distance_km = 1.0f;
    if (distance_km > 20.0f) distance_km = 20.0f;
    if (frequency_mhz < 150.0f) frequency_mhz = 150.0f;
    if (frequency_mhz > 1500.0f) frequency_mhz = 1500.0f;

    // Mobile antenna height correction factor
    float C_h;
    if (city_type == CITY_TYPE_LARGE && frequency_mhz > 400.0f) {
        C_h = 3.2f * powf(log10f(11.75f * mobile_height_m), 2.0f) - 4.97f;
    } else {
        C_h = 0.8f + (1.1f * log10f(frequency_mhz) - 0.7f) * mobile_height_m
              - 1.56f * log10f(frequency_mhz);
    }

    // Hata model
    float pl = (69.55f + 26.16f * log10f(frequency_mhz)
                - 13.82f * log10f(base_height_m)
                - C_h
                + (44.9f - 6.55f * log10f(base_height_m)) * log10f(distance_km));

    return pl;
}

float path_loss_two_ray_ground(float distance_m,
                               float tx_height_m,
                               float rx_height_m,
                               float frequency_mhz) {
    if (distance_m < 1.0f) {
        distance_m = 1.0f;
    }

    // Wavelength
    float wavelength = SPEED_OF_LIGHT / (frequency_mhz * 1e6f);

    // Crossover distance (where two-ray model becomes valid)
    float crossover = (4.0f * M_PI * tx_height_m * rx_height_m) / wavelength;

    if (distance_m < crossover) {
        // Use free space for short distances
        return path_loss_free_space(distance_m, frequency_mhz);
    }

    // Two-ray model: PL = 40*log10(d) - 20*log10(ht) - 20*log10(hr)
    float pl = 40.0f * log10f(distance_m) - 20.0f * log10f(tx_height_m) - 20.0f * log10f(rx_height_m);
    return pl;
}

/* ============================================================================
 * Core Functions
 * ============================================================================ */

int field_generator_init(field_generator_t* gen,
                        float x_min, float x_max,
                        float y_min, float y_max,
                        float resolution,
                        float altitude,
                        path_loss_model_t model) {
    if (!gen) return -1;
    if (resolution <= 0.0f) return -1;

    gen->field.x_min = x_min;
    gen->field.x_max = x_max;
    gen->field.y_min = y_min;
    gen->field.y_max = y_max;
    gen->field.resolution = resolution;
    gen->field.altitude = altitude;

    // Calculate grid dimensions
    gen->field.width = (int)((x_max - x_min) / resolution) + 1;
    gen->field.height = (int)((y_max - y_min) / resolution) + 1;

    if (gen->field.width > MAX_GRID_DIM || gen->field.height > MAX_GRID_DIM) {
        return -1;
    }

    // Allocate grids
    int grid_size = gen->field.width * gen->field.height;
    gen->field.signal_strength = (float*)malloc(sizeof(float) * grid_size);
    gen->field.gradient_x = (float*)malloc(sizeof(float) * grid_size);
    gen->field.gradient_y = (float*)malloc(sizeof(float) * grid_size);

    if (!gen->field.signal_strength || !gen->field.gradient_x || !gen->field.gradient_y) {
        free(gen->field.signal_strength);
        free(gen->field.gradient_x);
        free(gen->field.gradient_y);
        return -1;
    }

    // Initialize to low signal
    for (int i = 0; i < grid_size; i++) {
        gen->field.signal_strength[i] = -150.0f;
        gen->field.gradient_x[i] = 0.0f;
        gen->field.gradient_y[i] = 0.0f;
    }

    gen->num_sources = 0;
    gen->model = model;
    gen->path_loss_exponent = DEFAULT_PATH_LOSS_EXPONENT;
    gen->shadowing_std_db = 8.0f;
    gen->base_height_m = 30.0f;
    gen->mobile_height_m = 1.5f;
    gen->city_type = CITY_TYPE_MEDIUM;

    return 0;
}

void field_generator_free(field_generator_t* gen) {
    if (!gen) return;

    if (gen->field.signal_strength) {
        free(gen->field.signal_strength);
        gen->field.signal_strength = NULL;
    }

    if (gen->field.gradient_x) {
        free(gen->field.gradient_x);
        gen->field.gradient_x = NULL;
    }

    if (gen->field.gradient_y) {
        free(gen->field.gradient_y);
        gen->field.gradient_y = NULL;
    }
}

int field_generator_add_source(field_generator_t* gen,
                               const signal_source_t* source) {
    if (!gen || !source) return -1;
    if (gen->num_sources >= MAX_SIGNAL_SOURCES) return -1;

    gen->sources[gen->num_sources] = *source;
    gen->num_sources++;

    return 0;
}

void field_generator_clear_sources(field_generator_t* gen) {
    if (!gen) return;
    gen->num_sources = 0;
}

int field_generator_compute(field_generator_t* gen) {
    if (!gen) return -1;
    if (gen->num_sources == 0) return -1;

    int grid_size = gen->field.width * gen->field.height;

    // Initialize field to low values
    for (int i = 0; i < grid_size; i++) {
        gen->field.signal_strength[i] = -150.0f;
    }

    // Compute contribution from each source
    for (int src_idx = 0; src_idx < gen->num_sources; src_idx++) {
        const signal_source_t* source = &gen->sources[src_idx];
        if (!source->active) continue;

        // Source position
        float tx_x = source->position.x;
        float tx_y = source->position.y;
        float tx_z = source->position.z;

        // For each grid point
        for (int iy = 0; iy < gen->field.height; iy++) {
            for (int ix = 0; ix < gen->field.width; ix++) {
                int idx = iy * gen->field.width + ix;

                // Receiver position
                float rx_x = gen->field.x_min + ix * gen->field.resolution;
                float rx_y = gen->field.y_min + iy * gen->field.resolution;
                float rx_z = gen->field.altitude;

                // 3D distance
                float dx = rx_x - tx_x;
                float dy = rx_y - tx_y;
                float dz = rx_z - tx_z;
                float distance = sqrtf(dx*dx + dy*dy + dz*dz);

                // Compute path loss
                float path_loss;
                switch (gen->model) {
                    case PATH_LOSS_FREE_SPACE:
                        path_loss = path_loss_free_space(distance, source->frequency_mhz);
                        break;

                    case PATH_LOSS_LOG_DISTANCE:
                        path_loss = path_loss_log_distance(distance, gen->path_loss_exponent, 1.0f, 40.0f);
                        break;

                    case PATH_LOSS_HATA_URBAN:
                        path_loss = path_loss_hata_urban(distance, source->frequency_mhz,
                                                         tx_z, rx_z, gen->city_type);
                        break;

                    case PATH_LOSS_TWO_RAY:
                        path_loss = path_loss_two_ray_ground(distance, tx_z, rx_z, source->frequency_mhz);
                        break;

                    default:
                        path_loss = path_loss_log_distance(distance, gen->path_loss_exponent, 1.0f, 40.0f);
                        break;
                }

                // Received signal strength: Pr = Pt + Gt - PL
                float received_power = source->transmit_power_dbm + source->antenna_gain_dbi - path_loss;

                // Take maximum across sources (assumes different frequencies)
                if (received_power > gen->field.signal_strength[idx]) {
                    gen->field.signal_strength[idx] = received_power;
                }
            }
        }
    }

    return 0;
}

int field_generator_compute_gradient(field_generator_t* gen) {
    if (!gen) return -1;
    if (!gen->field.signal_strength) return -1;

    int width = gen->field.width;
    int height = gen->field.height;
    float res = gen->field.resolution;

    // Compute gradients using central differences
    for (int iy = 0; iy < height; iy++) {
        for (int ix = 0; ix < width; ix++) {
            int idx = iy * width + ix;

            // Gradient X (forward/backward/central difference)
            if (ix == 0) {
                gen->field.gradient_x[idx] = (gen->field.signal_strength[idx + 1]
                                             - gen->field.signal_strength[idx]) / res;
            } else if (ix == width - 1) {
                gen->field.gradient_x[idx] = (gen->field.signal_strength[idx]
                                             - gen->field.signal_strength[idx - 1]) / res;
            } else {
                gen->field.gradient_x[idx] = (gen->field.signal_strength[idx + 1]
                                             - gen->field.signal_strength[idx - 1]) / (2.0f * res);
            }

            // Gradient Y
            if (iy == 0) {
                gen->field.gradient_y[idx] = (gen->field.signal_strength[idx + width]
                                             - gen->field.signal_strength[idx]) / res;
            } else if (iy == height - 1) {
                gen->field.gradient_y[idx] = (gen->field.signal_strength[idx]
                                             - gen->field.signal_strength[idx - width]) / res;
            } else {
                gen->field.gradient_y[idx] = (gen->field.signal_strength[idx + width]
                                             - gen->field.signal_strength[idx - width]) / (2.0f * res);
            }
        }
    }

    return 0;
}

float field_generator_get_signal(const field_generator_t* gen, float x, float y) {
    if (!gen || !gen->field.signal_strength) return -150.0f;

    // Check bounds
    if (x < gen->field.x_min || x > gen->field.x_max ||
        y < gen->field.y_min || y > gen->field.y_max) {
        return -150.0f;
    }

    // Bilinear interpolation
    float x_norm = (x - gen->field.x_min) / gen->field.resolution;
    float y_norm = (y - gen->field.y_min) / gen->field.resolution;

    int x_i = (int)floorf(x_norm);
    int y_i = (int)floorf(y_norm);

    if (x_i >= gen->field.width - 1 || y_i >= gen->field.height - 1) {
        int xi = (x_i < gen->field.width) ? x_i : gen->field.width - 1;
        int yi = (y_i < gen->field.height) ? y_i : gen->field.height - 1;
        return gen->field.signal_strength[yi * gen->field.width + xi];
    }

    // Interpolation weights
    float wx = x_norm - x_i;
    float wy = y_norm - y_i;

    // Four corner values
    float s00 = gen->field.signal_strength[y_i * gen->field.width + x_i];
    float s10 = gen->field.signal_strength[y_i * gen->field.width + x_i + 1];
    float s01 = gen->field.signal_strength[(y_i + 1) * gen->field.width + x_i];
    float s11 = gen->field.signal_strength[(y_i + 1) * gen->field.width + x_i + 1];

    // Bilinear interpolation
    float signal = (1.0f - wx) * (1.0f - wy) * s00
                 + wx * (1.0f - wy) * s10
                 + (1.0f - wx) * wy * s01
                 + wx * wy * s11;

    return signal;
}

int field_generator_get_gradient(const field_generator_t* gen,
                                 float x, float y,
                                 float* grad_x, float* grad_y) {
    if (!gen || !gen->field.gradient_x || !gen->field.gradient_y) return -1;
    if (!grad_x || !grad_y) return -1;

    // Check bounds
    if (x < gen->field.x_min || x > gen->field.x_max ||
        y < gen->field.y_min || y > gen->field.y_max) {
        *grad_x = 0.0f;
        *grad_y = 0.0f;
        return -1;
    }

    // Bilinear interpolation for gradients
    float x_norm = (x - gen->field.x_min) / gen->field.resolution;
    float y_norm = (y - gen->field.y_min) / gen->field.resolution;

    int x_i = (int)floorf(x_norm);
    int y_i = (int)floorf(y_norm);

    if (x_i >= gen->field.width - 1 || y_i >= gen->field.height - 1) {
        int xi = (x_i < gen->field.width) ? x_i : gen->field.width - 1;
        int yi = (y_i < gen->field.height) ? y_i : gen->field.height - 1;
        int idx = yi * gen->field.width + xi;
        *grad_x = gen->field.gradient_x[idx];
        *grad_y = gen->field.gradient_y[idx];
        return 0;
    }

    float wx = x_norm - x_i;
    float wy = y_norm - y_i;

    // Interpolate gradient_x
    int idx00 = y_i * gen->field.width + x_i;
    int idx10 = y_i * gen->field.width + x_i + 1;
    int idx01 = (y_i + 1) * gen->field.width + x_i;
    int idx11 = (y_i + 1) * gen->field.width + x_i + 1;

    *grad_x = (1.0f - wx) * (1.0f - wy) * gen->field.gradient_x[idx00]
            + wx * (1.0f - wy) * gen->field.gradient_x[idx10]
            + (1.0f - wx) * wy * gen->field.gradient_x[idx01]
            + wx * wy * gen->field.gradient_x[idx11];

    *grad_y = (1.0f - wx) * (1.0f - wy) * gen->field.gradient_y[idx00]
            + wx * (1.0f - wy) * gen->field.gradient_y[idx10]
            + (1.0f - wx) * wy * gen->field.gradient_y[idx01]
            + wx * wy * gen->field.gradient_y[idx11];

    return 0;
}

int field_generator_find_coverage_holes(const field_generator_t* gen,
                                        float threshold_dbm,
                                        coverage_holes_t* result) {
    if (!gen || !result) return -1;

    // Count holes
    int hole_count = 0;
    for (int i = 0; i < gen->field.width * gen->field.height; i++) {
        if (gen->field.signal_strength[i] < threshold_dbm) {
            hole_count++;
        }
    }

    result->num_holes = hole_count;

    if (hole_count == 0) {
        result->positions = NULL;
        result->signal_strength = NULL;
        return 0;
    }

    // Allocate arrays
    result->positions = (point3d_t*)malloc(sizeof(point3d_t) * hole_count);
    result->signal_strength = (float*)malloc(sizeof(float) * hole_count);

    if (!result->positions || !result->signal_strength) {
        free(result->positions);
        free(result->signal_strength);
        return -1;
    }

    // Fill arrays
    int idx = 0;
    for (int iy = 0; iy < gen->field.height; iy++) {
        for (int ix = 0; ix < gen->field.width; ix++) {
            int grid_idx = iy * gen->field.width + ix;
            if (gen->field.signal_strength[grid_idx] < threshold_dbm) {
                result->positions[idx].x = gen->field.x_min + ix * gen->field.resolution;
                result->positions[idx].y = gen->field.y_min + iy * gen->field.resolution;
                result->positions[idx].z = 0.0f;
                result->signal_strength[idx] = gen->field.signal_strength[grid_idx];
                idx++;
            }
        }
    }

    return 0;
}

void coverage_holes_free(coverage_holes_t* holes) {
    if (!holes) return;

    if (holes->positions) {
        free(holes->positions);
        holes->positions = NULL;
    }

    if (holes->signal_strength) {
        free(holes->signal_strength);
        holes->signal_strength = NULL;
    }

    holes->num_holes = 0;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

int create_cmrcm_field(field_generator_t* gen,
                      float resolution,
                      float altitude) {
    if (!gen) return -1;

    // CMRCM field bounds (approx 3km x 3km)
    float bounds_min = -1500.0f;
    float bounds_max = 1500.0f;

    if (field_generator_init(gen, bounds_min, bounds_max, bounds_min, bounds_max,
                            resolution, altitude, PATH_LOSS_LOG_DISTANCE) != 0) {
        return -1;
    }

    // Add known cell towers (from CellMapper database)

    // Tower 1: 650m north
    signal_source_t tower1 = {
        .position = {0.0f, 650.0f, 50.0f},
        .transmit_power_dbm = 43.0f,
        .frequency_mhz = 1900.0f,
        .antenna_gain_dbi = 17.0f,
        .active = true
    };
    strncpy(tower1.source_id, "WQNR305_WQNX339", sizeof(tower1.source_id) - 1);
    field_generator_add_source(gen, &tower1);

    // Tower 2: 1.7km NW
    signal_source_t tower2 = {
        .position = {-1200.0f, 900.0f, 60.0f},
        .transmit_power_dbm = 43.0f,
        .frequency_mhz = 850.0f,
        .antenna_gain_dbi = 17.0f,
        .active = true
    };
    strncpy(tower2.source_id, "WQQL972_WQRV287", sizeof(tower2.source_id) - 1);
    field_generator_add_source(gen, &tower2);

    // Tower 3: 1.9km SE
    signal_source_t tower3 = {
        .position = {800.0f, -1700.0f, 55.0f},
        .transmit_power_dbm = 43.0f,
        .frequency_mhz = 2100.0f,
        .antenna_gain_dbi = 17.0f,
        .active = true
    };
    strncpy(tower3.source_id, "Chauncy_Hill", sizeof(tower3.source_id) - 1);
    field_generator_add_source(gen, &tower3);

    return 0;
}

int field_generator_get_stats(const field_generator_t* gen,
                              float* min_dbm,
                              float* max_dbm,
                              float* mean_dbm) {
    if (!gen || !min_dbm || !max_dbm || !mean_dbm) return -1;
    if (!gen->field.signal_strength) return -1;

    int grid_size = gen->field.width * gen->field.height;
    float min_val = gen->field.signal_strength[0];
    float max_val = gen->field.signal_strength[0];
    double sum = 0.0;

    for (int i = 0; i < grid_size; i++) {
        float val = gen->field.signal_strength[i];
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
        sum += val;
    }

    *min_dbm = min_val;
    *max_dbm = max_val;
    *mean_dbm = (float)(sum / grid_size);

    return 0;
}
