/* cblas_example.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cblas.h"

int main ( )
{
    CBLAS_LAYOUT Layout;
    CBLAS_UPLO Uplo;
    CBLAS_TRANSPOSE Trans;

    int m, n;
    double alpha, beta;
    double *A, *C, *Mean, *M_1, *norm;
    int lda, ldc, i, j;

    Layout = CblasColMajor;
    Uplo = CblasLower;
    Trans = CblasConjTrans;

    m = 10; /* Size of Column ( the number of rows ) */
    n = 4; /* Size of Row ( the number of columns ) */
    alpha = 1;
    beta = 0;
    lda = 10; /* Leading dimension of 10 * 4 matrix is 10 */
    ldc = 4;

    A = (double *)malloc(sizeof(double)*m*n);
    Mean = (double *)malloc(sizeof(double)*1*n);
    M_1 = (double *)malloc(sizeof(double)*m*1);
    C = (double *)malloc(sizeof(double)*n*n);
    norm = (double *)malloc(sizeof(double)*1*n);

    FILE* fp = fopen("matrix.bin", "r");
    fread(A, sizeof(double), m*n, fp);

    // Mean: 1 * n matrix 
    for (i = 0; i < n; ++i) {
        double sum = 0.0;
        for (j = 0; j < m; ++j) {
            sum += A[i*m+j];
        }
        Mean[i] = sum / m;
    }

//    // print matrix Mean
//    printf("======= Mean =======\n");
//    for (i = 0; i < n; ++i) {
//        printf("%f ", Mean[i]);
//    }
//    printf("\n");

    // M_1: m * 1 matrix
    for (i = 0; i < m; ++i) {
        M_1[i] = -1.0;
    }

    // A = A - M_1 * Mean, i.e. A = A - u
    cblas_dger(Layout, m, n, alpha, M_1, 1, Mean, 1, A, m);

//    // print matrix (A-u)
//    printf("\n======= A-u =======\n");
//    for (i = 0; i < m; ++i) {
//        for (j = 0; j < n; ++j) {
//            printf("%f ", A[i+j*m]);
//        }
//        printf("\n");
//    }

    // compute column norm
    for (i = 0; i < n; ++i) {
        norm[i] = cblas_dnrm2(m, &A[i*m], 1);
        cblas_dscal(m, 1/norm[i], &A[i*m], 1);
    }

//    // print the norm array
//    printf("\n======= column norm =======\n");
//    for (i = 0; i < n; ++i) {
//        printf("%f ", norm[i]);
//    }
//    printf("\n");

//    // print the normalized matrix (A-u)
//    printf("\n======= (A-u)/norm =======\n");
//    for (i = 0; i < m; ++i) {
//        for (j = 0; j < n; ++j) {
//            printf("%f ", A[i+j*m]);
//        }
//        printf("\n");
//    }


    // correlation matrix C
    cblas_dsyrk(Layout, Uplo, Trans, n, m, alpha, A, lda, beta, C, ldc);
    printf("\n======= Correlation =======\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%f ", C[i+j*n]);
        }
        printf("\n");
    }
    printf("\n");


    free(norm);
    free(C);
    free(M_1);
    free(Mean);
    free(A);
    fclose(fp);
    return 0;
}
