/**
 * @file eigendecomp.c
 * @brief Eigenvalue decomposition wrapper for LAPACK
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <stdlib.h>
#include <string.h>
#include <complex.h>

/* LAPACK function declaration (external linkage) */
/* ZHEEV: Computes eigenvalues/eigenvectors of complex Hermitian matrix */
extern void zheev_(
    const char* jobz,      /* 'N' = eigenvalues only, 'V' = eigenvectors too */
    const char* uplo,      /* 'U' = upper triangle, 'L' = lower triangle */
    const int* n,          /* Order of matrix A */
    cdouble_t* A,          /* Input matrix (destroyed), output eigenvectors */
    const int* lda,        /* Leading dimension of A */
    double* W,             /* Output eigenvalues (ascending order) */
    cdouble_t* work,       /* Workspace */
    const int* lwork,      /* Size of workspace (-1 for query) */
    double* rwork,         /* Real workspace */
    int* info              /* Status: 0=success, <0=illegal arg, >0=failed */
);

/* ============================================================================
 * Eigendecomposition
 * ============================================================================ */

music_status_t compute_eigendecomposition(
    const covariance_matrix_t* covariance,
    eigendecomp_result_t* eigendecomp)
{
    // TODO: Implement eigendecomposition using LAPACK ZHEEV
    // eigendecomposition
    //
    // Implementation steps:
    //
    // 1. Copy covariance matrix (ZHEEV destroys input)
    //    [V, D] = eig(R)
    //
    // 2. Call LAPACK ZHEEV to compute eigenvalues and eigenvectors
    //    Parameters:
    //    - jobz = 'V' (compute eigenvectors)
    //    - uplo = 'U' (upper triangular)
    //    - n = M (matrix dimension)
    //    - A = copy of covariance matrix (destroyed, replaced with eigenvectors)
    //    - W = output eigenvalues (in ascending order)
    //
    // 3. Sort eigenvalues in descending order
    //    [eigenvalues, idx] = sort(diag(D), 'descend')
    //    LAPACK returns ascending order, so reverse array
    //
    // 4. Reorder eigenvectors to match sorted eigenvalues
    //    eigenvectors = V(:, idx)
    //
    // Expected behavior:
    // - First K eigenvalues are large (signal subspace)
    // - Last M-K eigenvalues are small and similar (noise subspace)
    // - For CMRCM with K=4 sources: expect eigenvalue drop after index 4
    //
    // Workspace sizing:
    // - Query optimal workspace size by calling with lwork=-1
    // - Allocate work array of optimal size
    // - Call again with actual computation

    if (!covariance || !eigendecomp) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = covariance->M;
    if (eigendecomp->M != M) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Copy covariance matrix (ZHEEV destroys input)
    cdouble_t* A_work = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
    if (!A_work) {
        return MUSIC_ERROR_MEMORY;
    }
    memcpy(A_work, covariance->data, M * M * sizeof(cdouble_t));

    // Allocate workspace
    const char jobz = 'V';  // Compute eigenvectors
    const char uplo = 'U';  // Upper triangular
    int lwork = -1;         // Query workspace size
    cdouble_t work_query;
    double* rwork = (double*)malloc((3 * M - 2) * sizeof(double));
    int info;

    // Query optimal workspace size
    zheev_(&jobz, &uplo, &M, A_work, &M, eigendecomp->eigenvalues,
           &work_query, &lwork, rwork, &info);

    if (info != 0) {
        free(A_work);
        free(rwork);
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    // Allocate optimal workspace
    lwork = (int)creal(work_query);
    cdouble_t* work = (cdouble_t*)malloc(lwork * sizeof(cdouble_t));
    if (!work) {
        free(A_work);
        free(rwork);
        return MUSIC_ERROR_MEMORY;
    }

    // Perform eigendecomposition
    zheev_(&jobz, &uplo, &M, A_work, &M, eigendecomp->eigenvalues,
           work, &lwork, rwork, &info);

    if (info != 0) {
        free(A_work);
        free(work);
        free(rwork);
        return MUSIC_ERROR_LAPACK_FAILED;
    }

    // FIX: Clamp negative eigenvalues to zero
    // Hermitian matrices should have non-negative eigenvalues, but numerical
    // errors can produce small negative values. Clamp to 0.
    for (int i = 0; i < M; i++) {
        if (eigendecomp->eigenvalues[i] < 0.0) {
            eigendecomp->eigenvalues[i] = 0.0;
        }
    }

    // ZHEEV returns eigenvalues in ascending order - reverse to descending
    for (int i = 0; i < M / 2; i++) {
        double temp = eigendecomp->eigenvalues[i];
        eigendecomp->eigenvalues[i] = eigendecomp->eigenvalues[M - 1 - i];
        eigendecomp->eigenvalues[M - 1 - i] = temp;
    }

    // Copy and reorder eigenvectors (columns) to match descending eigenvalues
    for (int j = 0; j < M; j++) {
        int j_src = M - 1 - j;  // Reverse column order
        for (int i = 0; i < M; i++) {
            eigendecomp->eigenvectors[i + j * M] = A_work[i + j_src * M];
        }
    }

    // Cleanup
    free(A_work);
    free(work);
    free(rwork);

    return MUSIC_SUCCESS;
}

/* ============================================================================
 * Source Number Estimation
 * ============================================================================ */

music_status_t estimate_num_sources(
    const double* eigenvalues,
    int num_eigenvalues,
    int num_snapshots,
    bool use_mdl,
    int* estimated_sources)
{
    // TODO: Implement MDL/AIC model order selection
    // estimate_num_sources
    //
    // Implementation steps:
    //
    // MDL criterion (Minimum Description Length):
    // for k = 0:M-1
    //     noise_eigs = eigenvalues(k+1:end)
    //     geo_mean = exp(mean(log(noise_eigs + 1e-12)))
    //     arith_mean = mean(noise_eigs)
    //     mdl(k+1) = -N*(M-k)*log(geo_mean/arith_mean) + 0.5*k*(2*M-k)*log(N)
    // end
    // [~, num_src] = min(mdl)
    // num_src = num_src - 1  % Adjust for 0-indexing
    //
    // AIC criterion (Akaike Information Criterion):
    // Similar to MDL but with different penalty term:
    // aic(k+1) = -2*N*(M-k)*log(geo_mean/arith_mean) + 2*k*(2*M-k)
    //
    // Interpretation:
    // - Geometric mean / arithmetic mean ratio measures eigenvalue uniformity
    // - For noise-only eigenvalues, this ratio approaches 1 (log → 0)
    // - For mixed signal+noise, ratio < 1 (log < 0, contributes to criterion)
    // - Penalty term discourages overestimation
    // - MDL more conservative than AIC
    //
    // Note: For CMRCM field, manual K=4 worked better than auto-estimation

    if (!eigenvalues || !estimated_sources) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = num_eigenvalues;
    int N = num_snapshots;

    double* criterion = (double*)malloc(M * sizeof(double));
    if (!criterion) {
        return MUSIC_ERROR_MEMORY;
    }

    // Compute MDL or AIC for each k
    for (int k = 0; k < M; k++) {
        // Extract noise eigenvalues (indices k to M-1)
        int n_noise = M - k;
        if (n_noise == 0) {
            criterion[k] = INFINITY;
            continue;
        }

        // Compute geometric and arithmetic means
        double log_sum = 0.0;
        double arith_sum = 0.0;
        for (int i = k; i < M; i++) {
            log_sum += log(eigenvalues[i] + EPSILON);
            arith_sum += eigenvalues[i];
        }
        double geo_mean = exp(log_sum / n_noise);
        double arith_mean = arith_sum / n_noise;

        // Compute criterion
        if (use_mdl) {
            // MDL
            criterion[k] = -N * n_noise * log(geo_mean / arith_mean) +
                          0.5 * k * (2 * M - k) * log(N);
        } else {
            // AIC
            criterion[k] = -2 * N * n_noise * log(geo_mean / arith_mean) +
                          2 * k * (2 * M - k);
        }
    }

    // Find minimum
    int min_idx = 0;
    double min_val = criterion[0];
    for (int k = 1; k < M; k++) {
        if (criterion[k] < min_val) {
            min_val = criterion[k];
            min_idx = k;
        }
    }

    *estimated_sources = min_idx;

    free(criterion);
    return MUSIC_SUCCESS;
}
