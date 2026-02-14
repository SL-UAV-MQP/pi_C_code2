/**
 * @file peak_detection.c
 * @brief Peak detection in MUSIC spectrum
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Peak Detection with Parabolic Interpolation
 * ============================================================================ */

music_status_t find_peaks(
    const music_spectrum_t* spectrum,
    int num_peaks,
    double min_peak_height_db,
    double min_peak_separation_deg,
    detected_sources_t* detected)
{
    // find_peaks
    // Detects peaks with parabolic interpolation for sub-degree accuracy

    if (!spectrum || !detected) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int num_bins = spectrum->num_azimuth_bins;
    if (num_bins < 3) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Calculate resolution and window size
    double resolution = spectrum->azimuth_deg[1] - spectrum->azimuth_deg[0];
    int window_half = (int)(min_peak_separation_deg / resolution / 2.0);
    if (window_half < 1) window_half = 1;

    // Find local maxima
    bool* is_local_max = (bool*)calloc(num_bins, sizeof(bool));
    if (!is_local_max) {
        return MUSIC_ERROR_MEMORY;
    }

    for (int i = 0; i < num_bins; i++) {
        // Check if current point is maximum in its window
        bool is_max = true;
        double current_val = spectrum->spectrum_db[i];

        for (int j = -window_half; j <= window_half; j++) {
            int idx = i + j;
            if (idx < 0 || idx >= num_bins || idx == i) continue;

            if (spectrum->spectrum_db[idx] > current_val) {
                is_max = false;
                break;
            }
        }

        is_local_max[i] = is_max && (current_val > min_peak_height_db);
    }

    // Count peaks and create list
    int n_peaks = 0;
    for (int i = 0; i < num_bins; i++) {
        if (is_local_max[i]) n_peaks++;
    }

    if (n_peaks == 0) {
        free(is_local_max);
        detected->num_sources = 0;
        return MUSIC_ERROR_NO_PEAKS;
    }

    // Extract peak indices and values
    int* peak_indices = (int*)malloc(n_peaks * sizeof(int));
    double* peak_values = (double*)malloc(n_peaks * sizeof(double));
    if (!peak_indices || !peak_values) {
        free(is_local_max);
        if (peak_indices) free(peak_indices);
        if (peak_values) free(peak_values);
        return MUSIC_ERROR_MEMORY;
    }

    int peak_count = 0;
    for (int i = 0; i < num_bins; i++) {
        if (is_local_max[i]) {
            peak_indices[peak_count] = i;
            peak_values[peak_count] = spectrum->spectrum_db[i];
            peak_count++;
        }
    }

    // Sort peaks by magnitude (descending) using bubble sort (simple for small N)
    for (int i = 0; i < n_peaks - 1; i++) {
        for (int j = 0; j < n_peaks - i - 1; j++) {
            if (peak_values[j] < peak_values[j + 1]) {
                // Swap values
                double temp_val = peak_values[j];
                peak_values[j] = peak_values[j + 1];
                peak_values[j + 1] = temp_val;
                // Swap indices
                int temp_idx = peak_indices[j];
                peak_indices[j] = peak_indices[j + 1];
                peak_indices[j + 1] = temp_idx;
            }
        }
    }

    // Extract top N peaks with parabolic interpolation
    int n_output = (num_peaks < n_peaks) ? num_peaks : n_peaks;
    if (n_output > MAX_SOURCES) n_output = MAX_SOURCES;

    detected->num_sources = n_output;

    for (int i = 0; i < n_output; i++) {
        int peak_idx = peak_indices[i];
        double azimuth = spectrum->azimuth_deg[peak_idx];
        double magnitude = spectrum->spectrum_db[peak_idx];

        // Parabolic interpolation for sub-degree accuracy
        if (peak_idx > 0 && peak_idx < num_bins - 1) {
            double y1 = spectrum->spectrum_db[peak_idx - 1];
            double y2 = spectrum->spectrum_db[peak_idx];
            double y3 = spectrum->spectrum_db[peak_idx + 1];

            // FIX: Proper denominator guard instead of adding EPSILON
            // Adding EPSILON to a near-zero denominator produces extreme offsets.
            // Instead: if denominator is too small (flat region), skip interpolation.
            double denominator = y1 - 2.0 * y2 + y3;
            if (fabs(denominator) > EPSILON) {
                double delta = 0.5 * (y1 - y3) / denominator;
                // Clamp delta to [-0.5, 0.5] bins (physical constraint)
                if (delta > 0.5) delta = 0.5;
                if (delta < -0.5) delta = -0.5;
                azimuth += delta * resolution;
            }
            // If denominator ≈ 0 (flat region), keep grid-point azimuth
        }

        // Fill detected source structure
        detected->sources[i].azimuth_deg = azimuth;
        detected->sources[i].elevation_deg = 0.0;
        detected->sources[i].magnitude_db = magnitude;

        // Confidence: normalized to [0, 1] based on threshold
        detected->sources[i].confidence =
            (magnitude - min_peak_height_db) / (0.0 - min_peak_height_db);
        if (detected->sources[i].confidence > 1.0) {
            detected->sources[i].confidence = 1.0;
        }

        strncpy(detected->sources[i].type, "MUSIC_peak", 31);
        detected->sources[i].type[31] = '\0';
    }

    // Cleanup
    free(is_local_max);
    free(peak_indices);
    free(peak_values);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Debug Helpers
 * ============================================================================ */

void print_detected_sources(const detected_sources_t* detected) {
    if (!detected) {
        printf("ERROR: NULL detected sources\n");
        return;
    }

    printf("=== Detected Sources ===\n");
    printf("Number of sources: %d\n\n", detected->num_sources);

    for (int i = 0; i < detected->num_sources; i++) {
        const detected_source_t* src = &detected->sources[i];
        printf("Source %d:\n", i + 1);
        printf("  Azimuth: %.2f°\n", src->azimuth_deg);
        printf("  Elevation: %.2f°\n", src->elevation_deg);
        printf("  Magnitude: %.1f dB\n", src->magnitude_db);
        printf("  Confidence: %.3f\n", src->confidence);
        printf("  Type: %s\n", src->type);
        printf("\n");
    }
}
