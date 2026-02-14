/**
 * @file covariance.c
 * @brief Spatial covariance matrix computation
 */

#include "music_uca_6.h"
#include "music_config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* BLAS/LAPACK function declarations (external linkage) */
/* ZHERK: Hermitian rank-k update: C := alpha*A*A^H + beta*C */
extern void zherk_(
    const char* uplo,    /* 'U' or 'L' */
    const char* trans,   /* 'N' or 'C' (conjugate transpose) */
    const int* n,        /* Order of matrix C */
    const int* k,        /* Number of columns of A */
    const double* alpha, /* Scalar alpha */
    const cdouble_t* A,  /* Matrix A */
    const int* lda,      /* Leading dimension of A */
    const double* beta,  /* Scalar beta */
    cdouble_t* C,        /* Output matrix C */
    const int* ldc       /* Leading dimension of C */
);

/* ============================================================================
 * Covariance Matrix Computation
 * ============================================================================ */

music_status_t compute_covariance_matrix(
    const signal_matrix_t* signals,
    bool use_forward_backward,
    covariance_matrix_t* covariance)
{
    // compute_covariance_matrix
    // Computes spatial covariance using BLAS ZHERK with optional FB averaging

    if (!signals || !covariance) {
        return MUSIC_ERROR_NULL_POINTER;
    }

    int M = signals->M;
    int N = signals->N;

    if (covariance->M != M) {
        return MUSIC_ERROR_INVALID_SIZE;
    }

    // Compute forward covariance using BLAS ZHERK
    // R_forward = (1/N) * X * X^H
    const char uplo = 'U';    // Upper triangular
    const char trans = 'N';   // No transpose
    double alpha = 1.0 / N;   // Scaling factor
    double beta = 0.0;        // Initialize output

    zherk_(&uplo, &trans, &M, &N, &alpha, signals->data, &M, &beta, covariance->data, &M);

    // Fill lower triangular part (ZHERK only computes upper)
    for (int j = 0; j < M; j++) {
        for (int i = j + 1; i < M; i++) {
            // R[i,j] = conj(R[j,i])
            covariance->data[i + j * M] = conj(covariance->data[j + i * M]);
        }
    }

    if (use_forward_backward) {
        // Phase-mode transformation for UCA:
        // Transform to beamspace using DFT matrix F, where F[m,n] = (1/sqrt(M)) * exp(-j*2*pi*m*n/M)
        // In beamspace, UCA manifold has Vandermonde structure → FB averaging is valid.
        // R_beamspace = F * R * F^H, then apply FB averaging, then R = F^H * R_bs * F

        // Step 1: Build DFT matrix F (M x M)
        cdouble_t* F = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
        cdouble_t* FH = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
        cdouble_t* temp = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
        cdouble_t* R_bs = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));

        if (!F || !FH || !temp || !R_bs) {
            free(F); free(FH); free(temp); free(R_bs);
            return MUSIC_ERROR_MEMORY;
        }

        double inv_sqrt_M = 1.0 / sqrt((double)M);
        for (int m = 0; m < M; m++) {
            for (int n = 0; n < M; n++) {
                double phase = -2.0 * M_PI * m * n / M;
                F[m + n * M] = inv_sqrt_M * (cos(phase) + I * sin(phase));
                FH[n + m * M] = conj(F[m + n * M]);
            }
        }

        // Step 2: R_bs = F * R * F^H
        // temp = R * F^H (M x M) x (M x M)
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                temp[i + j * M] = 0.0;
                for (int k = 0; k < M; k++) {
                    temp[i + j * M] += covariance->data[i + k * M] * FH[k + j * M];
                }
            }
        }
        // R_bs = F * temp
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                R_bs[i + j * M] = 0.0;
                for (int k = 0; k < M; k++) {
                    R_bs[i + j * M] += F[i + k * M] * temp[k + j * M];
                }
            }
        }

        // Step 3: FB averaging in beamspace (now valid due to Vandermonde structure)
        cdouble_t* R_backward = (cdouble_t*)malloc(M * M * sizeof(cdouble_t));
        if (!R_backward) {
            free(F); free(FH); free(temp); free(R_bs);
            return MUSIC_ERROR_MEMORY;
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                int i_flip = M - 1 - i;
                int j_flip = M - 1 - j;
                R_backward[i + j * M] = conj(R_bs[i_flip + j_flip * M]);
            }
        }
        for (int i = 0; i < M * M; i++) {
            R_bs[i] = 0.5 * (R_bs[i] + R_backward[i]);
        }

        // Step 4: Transform back: R = F^H * R_bs * F
        // temp = R_bs * F
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                temp[i + j * M] = 0.0;
                for (int k = 0; k < M; k++) {
                    temp[i + j * M] += R_bs[i + k * M] * F[k + j * M];
                }
            }
        }
        // covariance = F^H * temp
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                covariance->data[i + j * M] = 0.0;
                for (int k = 0; k < M; k++) {
                    covariance->data[i + j * M] += FH[i + k * M] * temp[k + j * M];
                }
            }
        }

        free(F);
        free(FH);
        free(temp);
        free(R_bs);
        free(R_backward);
    }

    // FIX: Add diagonal loading for numerical stability
    // Prevents ill-conditioning when noise eigenvalues are very small
    double diag_load = EPSILON * 100.0;  // 1e-10: small enough to not bias, large enough to stabilize
    for (int i = 0; i < M; i++) {
        covariance->data[i + i * M] += diag_load;
    }

    return MUSIC_SUCCESS;
}
