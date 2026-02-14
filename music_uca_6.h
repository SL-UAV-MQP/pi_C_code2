/**
 * @file music_uca_6.h
 * @brief Main MUSIC-UCA-6 algorithm API
 *
 * Complete MUSIC (Multiple Signal Classification) AOA estimation
 * for 6-element Uniform Circular Array with directional antennas.
 */

#ifndef MUSIC_UCA_6_H
#define MUSIC_UCA_6_H

#include "common_types.h"
#include "array_geometry.h"

/* ============================================================================
 * Main MUSIC Pipeline
 * ============================================================================ */

/**
 * @brief Complete MUSIC AOA estimation pipeline
 *
 * End-to-end processing from IQ samples to detected source angles.
 * Combines all steps: covariance → eigendecomposition → spectrum → peaks.
 *
 * Processing pipeline:
 * 1. Compute spatial covariance matrix (with optional FB averaging)
 * 2. Eigendecomposition to separate signal/noise subspaces
 * 3. Compute MUSIC pseudo-spectrum over azimuth grid
 * 4. Find peaks with parabolic interpolation
 *
 * @param[in] signals  Signal matrix (6 × N snapshots) in column-major
 * @param[in] geometry  Array geometry structure
 * @param[in] num_sources  Number of sources to detect
 * @param[in] use_fb_averaging  Enable forward-backward averaging
 * @param[out] spectrum  Output MUSIC spectrum (allocated by caller)
 * @param[out] detected  Output detected sources (allocated by caller)
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Input signals must be 6 × N complex matrix (column-major)
 * @note Output structures must be pre-allocated by caller
 * @note Processing time target: <100ms on Raspberry Pi 5B
 */
music_status_t music_estimate_aoa(
    const signal_matrix_t* signals,
    const array_geometry_t* geometry,
    int num_sources,
    bool use_fb_averaging,
    music_spectrum_t* spectrum,
    detected_sources_t* detected
);

/* ============================================================================
 * Steering Vector Computation
 * ============================================================================ */

/**
 * @brief Compute UCA steering vector with directional antenna pattern
 *
 * Calculates phase shifts for each antenna element based on DOA,
 * then applies directional antenna gain pattern.
 *
 * Mathematics:
 * - Phase shift: φ_m = k · (r_m · d)  where k = 2π/λ
 * - r_m = antenna position vector
 * - d = unit DOA vector [cos(el)cos(az), cos(el)sin(az), sin(el)]
 * - Base: a_m = exp(jφ_m)
 * - With pattern: a_m = a_m * G_m(az)
 * - Normalized: a = a / ||a||
 *
 * @param[in] geometry  Array geometry structure
 * @param[in] azimuth_deg  Azimuth angle in degrees (-180 to 180)
 * @param[in] elevation_deg  Elevation angle in degrees (0 for planar)
 * @param[in] wavelength_m  Wavelength in meters
 * @param[in] include_antenna_pattern  Apply directional gain (true/false)
 * @param[out] steering_vec  Output steering vector [6 × 1] complex
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Output vector is normalized: ||a|| = 1
 * @note Elevation is typically 0° for planar UCA
 */
music_status_t compute_steering_vector(
    const array_geometry_t* geometry,
    double azimuth_deg,
    double elevation_deg,
    double wavelength_m,
    bool include_antenna_pattern,
    cdouble_t* steering_vec  /* [6] */
);

/* ============================================================================
 * Covariance Matrix Computation
 * ============================================================================ */

/**
 * @brief Compute spatial covariance matrix with optional FB averaging
 *
 * Estimates covariance from snapshot data using sample averaging.
 * Optional forward-backward averaging improves estimation for coherent sources.
 *
 * Forward covariance:
 *   R_f = (1/N) * X * X^H
 *
 * Backward covariance (if enabled):
 *   R_b = J * conj(R_f) * J  where J = exchange matrix
 *   R = 0.5 * (R_f + R_b)
 *
 * @param[in] signals  Signal matrix (M × N) in column-major
 * @param[in] use_forward_backward  Enable FB averaging
 * @param[out] covariance  Output covariance matrix (M × M) in column-major
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Uses BLAS ZHERK/ZGEMM for efficient computation
 * @note Output is Hermitian positive semi-definite
 * @note FB averaging helps with multipath/coherent sources
 */
music_status_t compute_covariance_matrix(
    const signal_matrix_t* signals,
    bool use_forward_backward,
    covariance_matrix_t* covariance
);

/* ============================================================================
 * Eigendecomposition
 * ============================================================================ */

/**
 * @brief Eigendecomposition of Hermitian covariance matrix
 *
 * Computes eigenvalues and eigenvectors, sorted in descending order.
 * Separates signal subspace (large eigenvalues) from noise subspace.
 *
 * Decomposition:
 *   R = V * Λ * V^H
 * where:
 *   Λ = diag(λ_1, λ_2, ..., λ_M)  (descending)
 *   V = [v_1, v_2, ..., v_M]  (eigenvectors)
 *
 * Signal subspace: V_s = [v_1, ..., v_K]  (K largest)
 * Noise subspace:  V_n = [v_{K+1}, ..., v_M]
 *
 * @param[in] covariance  Input covariance matrix (M × M)
 * @param[out] eigendecomp  Output eigenvalues and eigenvectors
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Uses LAPACK ZHEEV for Hermitian eigenvalue problem
 * @note Eigenvalues are real, sorted descending
 * @note Eigenvectors are complex, orthonormal columns
 */
music_status_t compute_eigendecomposition(
    const covariance_matrix_t* covariance,
    eigendecomp_result_t* eigendecomp
);

/* ============================================================================
 * MUSIC Spectrum Computation
 * ============================================================================ */

/**
 * @brief Compute MUSIC angular pseudo-spectrum
 *
 * Scans azimuth angles and computes MUSIC metric using noise subspace.
 * Higher values indicate likely source directions.
 *
 * MUSIC spectrum:
 *   P_MUSIC(θ) = 1 / ||E_n^H * a(θ)||^2
 * where:
 *   E_n = noise subspace eigenvectors
 *   a(θ) = steering vector at angle θ
 *
 * In dB: P_dB(θ) = 10 * log10(P_MUSIC(θ) / max(P_MUSIC))
 *
 * @param[in] eigendecomp  Eigendecomposition result
 * @param[in] geometry  Array geometry
 * @param[in] num_sources  Number of signal sources (K)
 * @param[in] azimuth_min_deg  Minimum azimuth to scan (degrees)
 * @param[in] azimuth_max_deg  Maximum azimuth to scan (degrees)
 * @param[in] azimuth_resolution_deg  Angular step size (degrees)
 * @param[in] wavelength_m  Wavelength in meters
 * @param[out] spectrum  Output MUSIC spectrum
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note Noise subspace = eigenvectors corresponding to M-K smallest eigenvalues
 * @note Azimuth range typically [-180, 180] with 0.2° resolution
 * @note Output normalized to [0 dB, -X dB] range
 */
music_status_t compute_music_spectrum(
    const eigendecomp_result_t* eigendecomp,
    const array_geometry_t* geometry,
    int num_sources,
    double azimuth_min_deg,
    double azimuth_max_deg,
    double azimuth_resolution_deg,
    double wavelength_m,
    music_spectrum_t* spectrum
);

/* ============================================================================
 * Peak Detection
 * ============================================================================ */

/**
 * @brief Find peaks in MUSIC spectrum
 *
 * Detects local maxima above threshold with minimum separation.
 * Uses parabolic interpolation for sub-degree accuracy.
 *
 * Algorithm:
 * 1. Find local maxima using sliding window
 * 2. Filter by minimum height (-30 dB threshold)
 * 3. Enforce minimum separation (30° to prevent ghosts)
 * 4. Sort by magnitude, select top N
 * 5. Parabolic interpolation for refined angle
 *
 * Parabolic interpolation:
 *   δ = 0.5 * (y[i-1] - y[i+1]) / (y[i-1] - 2*y[i] + y[i+1])
 *   θ_refined = θ[i] + δ * Δθ
 *
 * @param[in] spectrum  MUSIC spectrum result
 * @param[in] num_peaks  Number of peaks to detect
 * @param[in] min_peak_height_db  Minimum peak height in dB (-30 dB)
 * @param[in] min_peak_separation_deg  Minimum separation (30°)
 * @param[out] detected  Output detected sources
 *
 * @return MUSIC_SUCCESS on success, MUSIC_ERROR_NO_PEAKS if no peaks found
 *
 * @note Peaks sorted by magnitude (strongest first)
 * @note Confidence score: (peak - threshold) / (0 - threshold)
 * @note Type field set to "MUSIC_peak"
 */
music_status_t find_peaks(
    const music_spectrum_t* spectrum,
    int num_peaks,
    double min_peak_height_db,
    double min_peak_separation_deg,
    detected_sources_t* detected
);

/* ============================================================================
 * Source Estimation (Model Order Selection)
 * ============================================================================ */

/**
 * @brief Estimate number of sources using MDL or AIC
 *
 * Information theoretic criteria for model order selection.
 *
 * MDL (Minimum Description Length):
 *   MDL(k) = -N*(M-k)*log(GM/AM) + 0.5*k*(2M-k)*log(N)
 * where:
 *   GM = geometric mean of noise eigenvalues
 *   AM = arithmetic mean of noise eigenvalues
 *   N = number of snapshots
 *   M = number of antennas
 *   k = assumed number of sources
 *
 * Estimate: K = argmin MDL(k)
 *
 * @param[in] eigenvalues  Array of eigenvalues (descending)
 * @param[in] num_eigenvalues  Number of eigenvalues (M)
 * @param[in] num_snapshots  Number of snapshots used (N)
 * @param[in] use_mdl  true for MDL, false for AIC
 * @param[out] estimated_sources  Estimated number of sources
 *
 * @return MUSIC_SUCCESS on success, error code otherwise
 *
 * @note MDL is more conservative (tends to underestimate)
 * @note AIC is more aggressive (tends to overestimate)
 * @note For CMRCM field, manual K=4 performed better than auto-estimate
 */
music_status_t estimate_num_sources(
    const double* eigenvalues,
    int num_eigenvalues,
    int num_snapshots,
    bool use_mdl,
    int* estimated_sources
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Convert linear power to dB
 * @param linear_power Linear power value
 * @return Power in dB (10*log10(power))
 */
double linear_to_db(double linear_power);

/**
 * @brief Convert dB to linear power
 * @param db_power Power in dB
 * @return Linear power (10^(db/10))
 */
double db_to_linear(double db_power);

/**
 * @brief Compute complex vector norm (L2)
 * @param vec Complex vector
 * @param length Vector length
 * @return L2 norm ||vec||_2
 */
double complex_vector_norm(const cdouble_t* vec, int length);

/**
 * @brief Normalize complex vector in-place
 * @param vec Complex vector to normalize
 * @param length Vector length
 * @return MUSIC_SUCCESS on success
 */
music_status_t normalize_complex_vector(cdouble_t* vec, int length);

/**
 * @brief Print detected sources
 * @param detected Detected sources structure
 */
void print_detected_sources(const detected_sources_t* detected);

/**
 * @brief Print MUSIC spectrum statistics
 * @param spectrum MUSIC spectrum structure
 */
void print_spectrum_stats(const music_spectrum_t* spectrum);

#endif /* MUSIC_UCA_6_H */
