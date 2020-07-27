/*
 *
 *  Created on: April 17, 2019
 *      Author: Grey Ballard (ballard@wfu.edu)
 *              Binfeng Xu (xub16@wfu.edu)
 */

#include "Tucker.hpp"
#include "Tucker_IO_Util.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "assert.h"
using namespace std;
extern "C" {
    //void dsyrk_(char*, char*, int*, int*,
    //double*, double*, int*, double*, double*,
    //int*);
    int dger_(int*, int*, double*,
	double*, int*, double*, int*,
	double*, int*);
    double dnrm2_(int*, double*, int*);
   // void dscal_(int*, double*, double*, int*);
}

int main(int argc, char* argv[])
{

  //
  // Get the name of the input file
  //
  std::string paramfn1 = Tucker::parseString(argc, (const char**)argv,
      "--metadata", "metadata.txt");
  std::string paramfn2 = Tucker::parseString(argc, (const char**)argv,
      "--matrix", "matrix.bin");


  //////////////////////////////////////

   // CBLAS_LAYOUT Layout;
   // CBLAS_UPLO Uplo;
   // CBLAS_TRANSPOSE Trans;

    int m, n, one = 1;
    double alpha, beta, tmp;
    double *A, *C, *Mean, *M_1, *norm;
    int lda, ldc, i, j;
    char charL, charC;
   // Layout = CblasColMajor;
   // Uplo = CblasLower;
   // Trans = CblasConjTrans;

    //int* rowcol;
    
    //Read in metadata
   // rowcol = (int *)malloc(sizeof(int)*2);
   // FILE* fp1 = fopen((const char*) paramfn1.c_str(),'r');
   // fread(rowcol, sizeof(i''nt), 2, fp1);
    int *rowcol = new int[2];
    ifstream metaRead;
    metaRead.open(paramfn1);
    metaRead >> rowcol[0];
    metaRead >> rowcol[1];
    //std::vector<std::string> fileAsString = Tucker::getFileAsStrings(paramfn);
    m = rowcol[0]; /* Size of Column ( the number of rows ) */
    n = rowcol[1]; /* Size of Row ( the number of columns ) */
    alpha = 1.0;
    beta = 0.0;
    lda = rowcol[0]; /* Leading dimension of 10 * 4 matrix is 10 */
    ldc = rowcol[1];

    
    A = (double *)malloc(sizeof(double)*m*n);
    Mean = (double *)malloc(sizeof(double)*1*n);
    M_1 = (double *)malloc(sizeof(double)*m*1);
    C = (double *)malloc(sizeof(double)*n*n);
    norm = (double *)malloc(sizeof(double)*1*n);

    ifstream matrixRead;
    matrixRead.open(paramfn2);
    matrixRead.read((char*) A,sizeof(double)*m*n);
    

    // Mean: 1 * n matrix
    for (i = 0; i < n; ++i) {
        double sum = 0.0;
        for (j = 0; j < m; ++j) {
            sum += A[i*m+j];
        }
        Mean[i] = sum / m;
    }

    // print matrix Mean
    printf("======= Mean =======\n");
    for (i = 0; i < n; ++i) {
        printf("%f ", Mean[i]);
	cout << Mean[i] << " ";
    }
    printf("\n");

    // M_1: m * 1 matrix
    for (i = 0; i < m; ++i) {
        M_1[i] = -1.0;
    }

    // A = A - M_1 * Mean, i.e. A = A - u
    dger_(&m, &n, &alpha, M_1, &one, Mean, &one, A,&m);

    // print matrix (A-u)
    printf("\n======= A-u =======\n");
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            //printf("%f ", A[i+j*m]);
        }
        //printf("\n");
    }

    // compute column norm

    for (i = 0; i < n; ++i) {
        norm[i] = dnrm2_(&m, &A[i*m], &one);
	tmp = 1.0 / norm[i];
        Tucker::dscal_(&m,&tmp, &A[i*m], &one);
    }

    // print the norm array
    printf("\n======= column norm =======\n");
    for (i = 0; i < n; ++i) {
        printf("%f ", norm[i]);
    }
    printf("\n");

    // print the normalized matrix (A-u)
    printf("\n======= (A-u)/norm =======\n");
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
           // printf("%f ", A[i+j*m]);
        }
       // printf("\n");
    }


    // correlation matrix C
    charL = 'L';
    charC = 'C';
    Tucker::dsyrk_(&charL,&charC, &n, &m, &alpha, A, &lda, &beta, C, &ldc);
    
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
    
  /////////////////////
  // Perform STHOSVD //
  /////////////////////
//GB: we'll eventually  use this code to compute the Tucker model
  /*  const Tucker::TuckerTensor* solution;
  solution = Tucker::STHOSVD(X, tol);
  //solution = Tucker::STHOSVD(X, R_dims);

  double xnorm = std::sqrt(X->norm2());
  double gnorm = std::sqrt(solution->G->norm2());
  std::cout << "Norm of input tensor: " << xnorm << std::endl;
  std::cout << "Norm of core tensor: " << gnorm << std::endl;
*/

  //
  // Free memory
  //
  //Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  //if(R_dims) Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(R_dims);

  //Tucker::MemoryManager::printMaxMemUsage();

  return EXIT_SUCCESS;

}
