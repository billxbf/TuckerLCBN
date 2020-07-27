/** \copyright
 * Copyright (2016) Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 * \n\n
 * BSD 2-Clause License
 * \n\n
 * Copyright (c) 2016, Sandia Corporation
 * All rights reserved.
 * \n\n
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * \n\n
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * \n\n
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * .
 * \n
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @file
 * \brief Contains functions relevant to computing a %Tucker decomposition.
 *
 * @author Alicia Klinvex
 */

#include "Tucker.hpp"
#include <algorithm>
#include <limits>
#include <fstream>
#include <cmath>
#include <cassert>
#include <chrono>
#include <iomanip>

/** \namespace Tucker \brief Contains the data structures and functions
 * necessary for a sequential tucker decomposition
 */
namespace Tucker {

/** \example Tucker_gram_test_file.cpp
 * \example Tucker_gram_test_nofile.cpp
 */
Matrix* computeGram(const Tensor* Y, const int n)
{
  if(Y == 0) {
    throw std::runtime_error("Tucker::computeGram(const Tensor* Y, const int n): Y is a null pointer");
  }
  if(Y->getNumElements() == 0) {
    throw std::runtime_error("Tucker::computeGram(const Tensor* Y, const int n): Y->getNumElements() == 0");
  }
  if(n < 0 || n >= Y->N()) {
    std::ostringstream oss;
    oss << "Tucker::computeGram(const Tensor* Y, const int n): n = "
        << n << " is not in the range [0," << Y->N() << ")";
    throw std::runtime_error(oss.str());
  }

  // Allocate memory for the Gram matrix
  int nrows = Y->size(n);
  Matrix* S = MemoryManager::safe_new<Matrix>(nrows,nrows);

  computeGram(Y, n, S->data(), nrows);

  return S;
}


void computeGram(const Tensor* Y, const int n, double* gram,
    const int stride)
{
  if(Y == 0) {
    throw std::runtime_error("Tucker::computeGram(const Tensor* Y, const int n, double* gram, const int stride): Y is a null pointer");
  }
  if(gram == 0) {
    throw std::runtime_error("Tucker::computeGram(const Tensor* Y, const int n, double* gram, const int stride): gram is a null pointer");
  }
  if(Y->getNumElements() == 0) {
    throw std::runtime_error("Tucker::computeGram(const Tensor* Y, const int n, double* gram, const int stride): Y->getNumElements() == 0");
  }
  if(stride < 1) {
    std::ostringstream oss;
    oss << "Tucker::computeGram(const Tensor* Y, const int n, double* gram, "
        << "const int stride): stride = " << stride << " < 1";
    throw std::runtime_error(oss.str());
  }
  if(n < 0 || n >= Y->N()) {
    std::ostringstream oss;
    oss << "Tucker::computeGram(const Tensor* Y, const int n, double* gram, "
        << "const int stride): n = " << n << " is not in the range [0,"
        << Y->N() << ")";
    throw std::runtime_error(oss.str());
  }

  int nrows = Y->size(n);

  // n = 0 is a special case
  // Y_0 is stored column major
  if(n == 0)
  {
    // Compute number of columns of Y_n
    // Technically, we could divide the total number of entries by n,
    // but that seems like a bad decision
    int ncols =1;
    for(int i=0; i<Y->N(); i++) {
      if(i != n) {
        ncols *= Y->size(i);
      }
    }

    // Call symmetric rank-k update
    // call dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    // C := alpha*A*A' + beta*C
    char uplo = 'U';
    char trans = 'N';
    double alpha = 1;
    double beta = 0;
    dsyrk_(&uplo, &trans, &nrows, &ncols, &alpha,
        Y->data(), &nrows, &beta, gram, &stride);
  }
  else
  {
    int ncols = 1;
    int nmats = 1;

    // Count the number of columns
    for(int i=0; i<n; i++) {
      ncols *= Y->size(i);
    }

    // Count the number of matrices
    for(int i=n+1; i<Y->N(); i++) {
      nmats *= Y->size(i);
    }

    // For each matrix...
    for(int i=0; i<nmats; i++) {
      // Call symmetric rank-k update
      // call dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      // C := alpha*A*A' + beta*C
      char uplo = 'U';
      char trans = 'T';
      double alpha = 1;
      double beta;
      if(i==0)
        beta = 0;
      else
        beta = 1;
      dsyrk_(&uplo, &trans, &nrows, &ncols, &alpha,
          Y->data()+i*nrows*ncols, &ncols, &beta,
          gram, &stride);
    }
  }
}

/// \example Tucker_eig_test.cpp
void computeEigenpairs(Matrix* G, double*& eigenvalues,
    const bool flipSign)
{
  if(G == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, const bool flipSign): G is a null pointer");
  }
  if(G->getNumElements() == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, const bool flipSign): G has no entries");
  }

  int nrows = G->nrows();
  eigenvalues = MemoryManager::safe_new_array<double>(nrows);

  // Compute the leading eigenvectors of S
  // call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
  char jobz = 'V';
  char uplo = 'U';
  int lwork = 8*nrows;
  double* work = MemoryManager::safe_new_array<double>(lwork);
  int info;
  dsyev_(&jobz, &uplo, &nrows, G->data(), &nrows,
      eigenvalues, work, &lwork, &info);

  // Check the error code
  if(info != 0)
    std::cerr << "Error: invalid error code returned by dsyev (" << info << ")\n";

  // The user will expect the eigenvalues to be sorted in descending order
  // LAPACK gives us the eigenvalues in ascending order
  for(int esubs=0; esubs<nrows-esubs-1; esubs++) {
    double temp = eigenvalues[esubs];
    eigenvalues[esubs] = eigenvalues[nrows-esubs-1];
    eigenvalues[nrows-esubs-1] = temp;
  }

  // Sort the eigenvectors too
  double* Gptr = G->data();
  const int ONE = 1;
  for(int esubs=0; esubs<nrows-esubs-1; esubs++) {
    dswap_(&nrows, Gptr+esubs*nrows, &ONE,
        Gptr+(nrows-esubs-1)*nrows, &ONE);
  }

  // Flip the sign if necessary
  if(flipSign)
  {
    for(int c=0; c<nrows; c++)
    {
      int maxIndex=0;
      double maxVal = std::abs(Gptr[c*nrows]);
      for(int r=1; r<nrows; r++)
      {
        double testVal = std::abs(Gptr[c*nrows+r]);
        if(testVal > maxVal) {
          maxIndex = r;
          maxVal = testVal;
        }
      }

      if(Gptr[c*nrows+maxIndex] < 0) {
        const double NEGONE = -1;
        dscal_(&nrows, &NEGONE, Gptr+c*nrows, &ONE);
      }
    }
  }

  MemoryManager::safe_delete_array<double>(work,lwork);
}

void computeEigenpairs(Matrix* G, double*& eigenvalues,
    Matrix*& eigenvectors, const int numEvecs, const bool flipSign)
{
  if(G == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, Matrix*& eigenvectors, const int numEvecs, const bool flipSign): G is a null pointer");
  }
  if(G->getNumElements() == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, Matrix*& eigenvectors, const int numEvecs, const bool flipSign): G has no entries");
  }
  if(numEvecs < 1) {
    std::ostringstream oss;
    oss << "Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, "
        << "Matrix*& eigenvectors, const int numEvecs, "
        << "const bool flipSign): numEvecs = " << numEvecs << " < 1";
    throw std::runtime_error(oss.str());
  }
  if(numEvecs > G->nrows()) {
    std::ostringstream oss;
    oss << "Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, "
        << "Matrix*& eigenvectors, const int numEvecs, "
        << "const bool flipSign): numEvecs = " << numEvecs
        << " > G->nrows() = " << G->nrows();
    throw std::runtime_error(oss.str());
  }

  computeEigenpairs(G, eigenvalues, flipSign);

  // Allocate memory for eigenvectors
  int numRows = G->nrows();
  eigenvectors = MemoryManager::safe_new<Matrix>(numRows,numEvecs);

  // Copy appropriate eigenvectors
  int nToCopy = numRows*numEvecs;
  const int ONE = 1;
  dcopy_(&nToCopy, G->data(), &ONE, eigenvectors->data(), &ONE);
}


void computeEigenpairs(Matrix* G, double*& eigenvalues,
    Matrix*& eigenvectors, const double thresh, const bool flipSign)
{
  if(G == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, Matrix*& eigenvectors, const int numEvecs, const bool flipSign): G is a null pointer");
  }
  if(G->getNumElements() == 0) {
    throw std::runtime_error("Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, Matrix*& eigenvectors, const int numEvecs, const bool flipSign): G has no entries");
  }
  if(thresh < 0) {
    std::ostringstream oss;
    oss << "Tucker::computeEigenpairs(Matrix* G, double*& eigenvalues, "
        << "Matrix*& eigenvectors, const int numEvecs, "
        << "const bool flipSign): thresh = " << thresh << " < 0";
    throw std::runtime_error(oss.str());
  }

  computeEigenpairs(G, eigenvalues, flipSign);

  // Compute number of things to copy
  int nrows = G->nrows();
  int numEvecs=nrows;
  double sum = 0;
  for(int i=nrows-1; i>=0; i--) {
    sum += eigenvalues[i];
    if(sum > thresh) {
      break;
    }
    numEvecs--;
  }

  // Allocate memory for eigenvectors
  int numRows = G->nrows();
  eigenvectors = MemoryManager::safe_new<Matrix>(numRows,numEvecs);

  // Copy appropriate eigenvectors
  int nToCopy = numRows*numEvecs;
  const int ONE = 1;
  dcopy_(&nToCopy, G->data(), &ONE, eigenvectors->data(), &ONE);
}


const struct TuckerTensor* STHOSVD(const Tensor* X,
    const double epsilon, bool flipSign)
{
  if(X == 0) {
    throw std::runtime_error("Tucker::STHOSVD(const Tensor* X, const double epsilon, bool flipSign): X is a null pointer");
  }
  if(X->getNumElements() == 0) {
    throw std::runtime_error("Tucker::STHOSVD(const Tensor* X, const double epsilon, bool flipSign): X has no entries");
  }
  if(epsilon < 0) {
    std::ostringstream oss;
    oss << "Tucker::STHOSVD(const Tensor* const X, const double epsilon, "
        << "bool flipSign): epsilon = " << epsilon << " < 0";
    throw std::runtime_error(oss.str());
  }

  int ndims = X->N();

  // Create a struct to store the factorization
  struct TuckerTensor* factorization = MemoryManager::safe_new<struct TuckerTensor>(ndims);
  factorization->total_timer_.start();

  // Compute the threshold
  double tensorNorm = X->norm2();
  double thresh = epsilon*epsilon*tensorNorm/X->N();
  std::cout << "\tAutoST-HOSVD::Tensor Norm: "
      << std::sqrt(tensorNorm) << "...\n";
  std::cout << "\tAutoST-HOSVD::Relative Threshold: "
      << thresh << "...\n";

  const Tensor* Y = X;

  // For each dimension...
  for(int n=0; n<X->N(); n++)
  {
    // Compute the Gram matrix
    // S = Y_n*Y_n'
    std::cout << "\tAutoST-HOSVD::Starting Gram(" << n << ")...\n";
    factorization->gram_timer_[n].start();
    Matrix* S = computeGram(Y,n);
    factorization->gram_timer_[n].stop();
    std::cout << "\tAutoST-HOSVD::Gram(" << n << ") time: "
        << factorization->gram_timer_[n].duration() << "s\n";

    // Compute the leading eigenvectors of S
    // call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    std::cout << "\tAutoST-HOSVD::Starting Evecs(" << n << ")...\n";
    factorization->eigen_timer_[n].start();
    computeEigenpairs(S, factorization->eigenvalues[n],
        factorization->U[n], thresh, flipSign);
    factorization->eigen_timer_[n].stop();
    std::cout << "\tAutoST-HOSVD::EVECS(" << n << ") time: "
        << factorization->eigen_timer_[n].duration() << "s\n";

    // Free the Gram matrix
    MemoryManager::safe_delete<Matrix>(S);

    // Perform the tensor times matrix multiplication
    std::cout << "\tAutoST-HOSVD::Starting TTM(" << n << ")...\n";
    factorization->ttm_timer_[n].start();
    Tensor* temp = ttm(Y,n,factorization->U[n],true);
    factorization->ttm_timer_[n].stop();
    std::cout << "\tAutoST-HOSVD::TTM(" << n << ") time: "
        << factorization->ttm_timer_[n].duration() << "s\n";
    if(n > 0) {
      MemoryManager::safe_delete<const Tensor>(Y);
    }
    Y = temp;

    size_t nnz = Y->getNumElements();
    std::cout << "Local tensor size after STHOSVD iteration "
        << n << ": " << Y->size() << ", or ";
    Tucker::printBytes(nnz*sizeof(double));
  }

  factorization->G = const_cast<Tensor*>(Y);
  factorization->total_timer_.stop();
  return factorization;
}


const struct TuckerTensor* STHOSVD(const Tensor* X,
    const SizeArray* reducedI, bool flipSign)
{
  if(X == 0) {
    throw std::runtime_error("Tucker::STHOSVD(const Tensor* X, const SizeArray* reducedI, bool flipSign): X is a null pointer");
  }
  if(X->getNumElements() == 0) {
    throw std::runtime_error("Tucker::STHOSVD(const Tensor* X, const SizeArray* reducedI, bool flipSign): X has no entries");
  }
  if(X->N() != reducedI->size()) {
    std::ostringstream oss;
    oss << "Tucker::STHOSVD(const Tensor* X, const SizeArray* reducedI, "
        << "bool flipSign): X->N() = " << X->N()
        << " != reducedI->size() = " << reducedI->size();
    throw std::runtime_error(oss.str());
  }
  for(int i=0; i<reducedI->size(); i++) {
    if((*reducedI)[i] <= 0) {
      std::ostringstream oss;
      oss << "Tucker::STHOSVD(const Tensor* X, const SizeArray* reducedI, "
      << "bool flipSign): reducedI[" << i << "] = " << (*reducedI)[i]
      << " <= 0";
      throw std::runtime_error(oss.str());
    }
    if((*reducedI)[i] > X->size(i)) {
      std::ostringstream oss;
      oss << "Tucker::STHOSVD(const Tensor* X, const SizeArray* reducedI, "
      << "bool flipSign): reducedI[" << i << "] = " << (*reducedI)[i]
      << " > X->size(" << i << ") = " << X->size(i);
      throw std::runtime_error(oss.str());
    }
  }

  // Create a struct to store the factorization
  struct TuckerTensor* factorization = MemoryManager::safe_new<struct TuckerTensor>(X->N());
  factorization->total_timer_.start();

  const Tensor* Y = X;

  // For each dimension...
  for(int n=0; n<X->N(); n++)
  {
    // Compute the Gram matrix
    // S = Y_n*Y_n'
    factorization->gram_timer_[n].start();
    Matrix* S = computeGram(Y,n);
    factorization->gram_timer_[n].stop();

    // Compute the leading eigenvectors of S
    // call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    factorization->eigen_timer_[n].start();
    computeEigenpairs(S, factorization->eigenvalues[n],
        factorization->U[n], (*reducedI)[n], flipSign);
    factorization->eigen_timer_[n].stop();

    MemoryManager::safe_delete<Matrix>(S);

    // Perform the tensor times matrix multiplication
    factorization->ttm_timer_[n].start();
    Tensor* temp = ttm(Y,n,factorization->U[n],true);
    factorization->ttm_timer_[n].stop();
    if(n > 0) {
      MemoryManager::safe_delete<const Tensor>(Y);
    }
    Y = temp;
  }

  factorization->G = const_cast<Tensor*>(Y);
  factorization->total_timer_.stop();
  return factorization;
}


Tensor* ttm(const Tensor* X, const int n,
    const Matrix* U, bool Utransp)
{
  if(X == 0) {
    throw std::runtime_error("Tucker::ttm(const Tensor* X, const int n, const Matrix* U, bool Utransp): X is a null pointer");
  }
  if(X->getNumElements() == 0) {
    throw std::runtime_error("Tucker::ttm(const Tensor* X, const int n, const Matrix* U, bool Utransp): X has no entries");
  }
  if(U == 0) {
    throw std::runtime_error("Tucker::ttm(const Tensor* X, const int n, const Matrix* U, bool Utransp): U is a null pointer");
  }
  if(U->getNumElements() == 0) {
    throw std::runtime_error("Tucker::ttm(const Tensor* X, const int n, const Matrix* U, bool Utransp): U has no entries");
  }
  if(n < 0 || n >= X->N()) {
    std::ostringstream oss;
    oss << "Tucker::ttm(const Tensor* X, const int n, const Matrix* U, "
        << "bool Utransp): n = " << n << " is not in the range [0,"
        << X->N() << ")";
    throw std::runtime_error(oss.str());
  }
  if(!Utransp && U->ncols() != X->size(n)) {
    std::ostringstream oss;
    // TODO: amk Oct 17 2016 Finish adding exceptions to this file
  }

  // Compute the number of rows for the resulting "matrix"
  int nrows;
  if(Utransp)
    nrows = U->ncols();
  else
    nrows = U->nrows();

  // Allocate space for the new tensor
  SizeArray I(X->N());
  for(int i=0; i<I.size(); i++) {
    if(i != n) {
      I[i] = X->size(i);
    }
    else {
      I[i] = nrows;
    }
  }
  Tensor* Y = MemoryManager::safe_new<Tensor>(I);

  // Call TTM
  ttm(X, n, U, Y, Utransp);

  // Return the tensor
  return Y;
}


Tensor* ttm(const Tensor* const X, const int n,
    const double* const Uptr, const int dimU,
    const int strideU, bool Utransp)
{
  // Allocate space for the new tensor
  SizeArray I(X->N());
  for(int i=0; i<I.size(); i++) {
    if(i != n) {
      I[i] = X->size(i);
    }
    else {
      I[i] = dimU;
    }
  }
  Tensor* Y = MemoryManager::safe_new<Tensor>(I);

  // Call TTM
  ttm(X, n, Uptr, strideU, Y, Utransp);

  // Return the tensor
  return Y;
}


void ttm(const Tensor* const X, const int n,
    const Matrix* const U, Tensor* Y, bool Utransp)
{
  // Check that the input is valid
  assert(U != 0);
  if(Utransp) {
    assert(U->nrows() == X->size(n));
    assert(U->ncols() == Y->size(n));
  }
  else {
    assert(U->ncols() == X->size(n));
    assert(U->nrows() == Y->size(n));
  }
  ttm(X, n, U->data(), U->nrows(), Y, Utransp);
}


void ttm(const Tensor* const X, const int n,
    const double* const Uptr, const int strideU,
    Tensor* Y, bool Utransp)
{
  // Check that the input is valid
  assert(X != 0);
  assert(Uptr != 0);
  assert(Y != 0);
  assert(n >= 0 && n < X->N());
  for(int i=0; i<X->N(); i++) {
    if(i != n) {
      assert(X->size(i) == Y->size(i));
    }
  }

  // Obtain the number of rows and columns of U
  int Unrows, Uncols;
  if(Utransp) {
    Unrows = X->size(n);
    Uncols = Y->size(n);
  }
  else {
    Uncols = X->size(n);
    Unrows = Y->size(n);
  }

  // n = 0 is a special case
  // Y_0 is stored column major
  if(n == 0)
  {
    // Compute number of columns of Y_n
    // Technically, we could divide the total number of entries by n,
    // but that seems like a bad decision
    size_t ncols = X->size().prod(1,X->N()-1);

    if(ncols > std::numeric_limits<int>::max()) {
      std::ostringstream oss;
      oss << "Error in Tucker::ttm: " << ncols
          << " is larger than std::numeric_limits<int>::max() ("
          << std::numeric_limits<int>::max() << ")";
      throw std::runtime_error(oss.str());
    }

    // Call matrix matrix multiply
    // call dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    // C := alpha*op( A )*op( B ) + beta*C
    // A, B and C are matrices, with op( A ) an m by k matrix,
    // op( B ) a k by n matrix and C an m by n matrix.
    char transa;
    char transb = 'N';
    int m = Y->size(n);
    int blas_n = (int)ncols;
    int k = X->size(n);
    int lda = strideU;
    int ldb = k;
    int ldc = m;
    double alpha = 1;
    double beta = 0;

    if(Utransp) {
      transa = 'T';
    } else {
      transa = 'N';
    }
    dgemm_(&transa, &transb, &m, &blas_n, &k, &alpha, Uptr,
        &lda, X->data(), &ldb, &beta, Y->data(), &ldc);
  }
  else
  {
    // Count the number of columns
    size_t ncols = X->size().prod(0,n-1);

    // Count the number of matrices
    size_t nmats = X->size().prod(n+1,X->N()-1,1);

    if(ncols > std::numeric_limits<int>::max()) {
      std::ostringstream oss;
      oss << "Error in Tucker::ttm: " << ncols
          << " is larger than std::numeric_limits<int>::max() ("
          << std::numeric_limits<int>::max() << ")";
      throw std::runtime_error(oss.str());
    }

    // For each matrix...
    for(size_t i=0; i<nmats; i++) {
      // Call matrix matrix multiply
      // call dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      // C := alpha*op( A )*op( B ) + beta*C
      // A, B and C are matrices, with op( A ) an m by k matrix,
      // op( B ) a k by n matrix and C an m by n matrix.
      char transa = 'N';
      char transb;
      int m = (int)ncols;
      int blas_n = Y->size(n);
      int k;
      int lda = (int)ncols;
      int ldb = strideU;
      int ldc = (int)ncols;
      double alpha = 1;
      double beta = 0;
      if(Utransp) {
        transb = 'N';
        k = Unrows;
      } else {
        transb = 'T';
        k = Uncols;
      }
      dgemm_(&transa, &transb, &m, &blas_n, &k, &alpha,
          X->data()+i*k*m, &lda, Uptr, &ldb, &beta,
          Y->data()+i*m*blas_n, &ldc);
    }
  }
}


MetricData* computeSliceMetrics(const Tensor* Y, const int mode, const int metrics)
{
  // If there are no slices, calling this function was a bad idea
  int numSlices = Y->size(mode);
  if(numSlices <= 0) {
    std::ostringstream oss;
    oss << "Tucker::computeSliceMetrics(const Tensor* Y, const int mode, const int metrics): "
        << "numSlices = " << numSlices << " <= 0";
    throw std::runtime_error(oss.str());
  }

  // Allocate memory for the result
  MetricData* result = MemoryManager::safe_new<MetricData>(metrics, numSlices);

  // Initialize the result
  double* delta;
  int* nArray;
  if((metrics & MEAN) || (metrics & VARIANCE)) {
    delta = MemoryManager::safe_new_array<double>(numSlices);
    nArray = MemoryManager::safe_new_array<int>(numSlices);
  }
  for(int i=0; i<numSlices; i++) {
    if(metrics & MIN) {
      result->getMinData()[i] = std::numeric_limits<double>::max();
    }
    if(metrics & MAX) {
      result->getMaxData()[i] = std::numeric_limits<double>::lowest();
    }
    if(metrics & SUM) {
      result->getSumData()[i] = 0;
    }
    if((metrics & MEAN) || (metrics & VARIANCE)) {
      result->getMeanData()[i] = 0;
      nArray[i] = 0;
    }
    if(metrics & VARIANCE) {
      result->getVarianceData()[i] = 0;
    }
  } // end for(int i=0; i<numSlices; i++)

  if(Y->getNumElements() == 0) {
    if((metrics & MEAN) || (metrics & VARIANCE)) {
      MemoryManager::safe_delete_array<double>(delta,numSlices);
      MemoryManager::safe_delete_array<int>(nArray,numSlices);
    }
    return result;
  }

  // Compute the result
  int ndims = Y->N();
  size_t numContig = Y->size().prod(0,mode-1,1); // Number of contiguous elements in a slice
  size_t numSetsContig = Y->size().prod(mode+1,ndims-1,1); // Number of sets of contiguous elements per slice
  size_t distBetweenSets = Y->size().prod(0,mode); // Distance between sets of contiguous elements

  const double* dataPtr;
  int slice;
  size_t i, c;
  #pragma omp parallel for default(shared) private(slice,i,c,dataPtr)
  for(slice=0; slice<numSlices; slice++)
  {
    dataPtr = Y->data() + slice*numContig;
    for(c=0; c<numSetsContig; c++)
    {
      for(i=0; i<numContig; i++)
      {
        if(metrics & MIN) {
          result->getMinData()[slice] = std::min(result->getMinData()[slice],dataPtr[i]);
        }
        if(metrics & MAX) {
          result->getMaxData()[slice] = std::max(result->getMaxData()[slice],dataPtr[i]);
        }
        if(metrics & SUM) {
          result->getSumData()[slice] += dataPtr[i];
        }
        if((metrics & MEAN) || (metrics & VARIANCE)) {
          delta[slice] = dataPtr[i] - result->getMeanData()[slice];
          nArray[slice]++;
          result->getMeanData()[slice] += (delta[slice]/nArray[slice]);
        }
        if(metrics & VARIANCE) {
          result->getVarianceData()[slice] +=
              (delta[slice]*(dataPtr[i]-result->getMeanData()[slice]));
        }
      } // end for(i=0; i<numContig; i++)
      dataPtr += distBetweenSets;
    } // end for(c=0; c<numSetsContig; c++)
  } // end for(slice=0; slice<numSlices; slice++)

  if((metrics & MEAN) || (metrics & VARIANCE)) {
    MemoryManager::safe_delete_array<double>(delta,numSlices);
    MemoryManager::safe_delete_array<int>(nArray,numSlices);
  }
  if(metrics & VARIANCE) {
    size_t sizeOfSlice = numContig*numSetsContig;
    for(int i=0; i<numSlices; i++) {
      result->getVarianceData()[i] /= (double)sizeOfSlice;
    }
  }

  return result;
}


// Shift is applied before scale
// We divide by scaleVals, not multiply
void transformSlices(Tensor* Y, int mode, const double* scales, const double* shifts)
{
  // If the tensor has no entries, no transformation is necessary
  size_t numEntries = Y->getNumElements();
  if(numEntries == 0)
    return;

  // Compute the result
  int ndims = Y->N();
  int numSlices = Y->size(mode);
  size_t numContig = Y->size().prod(0,mode-1,1); // Number of contiguous elements in a slice
  size_t numSetsContig = Y->size().prod(mode+1,ndims-1,1); // Number of sets of contiguous elements per slice
  size_t distBetweenSets = Y->size().prod(0,mode); // Distance between sets of contiguous elements

  double* dataPtr;
  int slice;
  size_t i, c;
  #pragma omp parallel for default(shared) private(slice,i,c,dataPtr)
  for(slice=0; slice<numSlices; slice++)
  {
    dataPtr = Y->data() + slice*numContig;
    for(c=0; c<numSetsContig; c++)
    {
      for(i=0; i<numContig; i++)
        dataPtr[i] = (dataPtr[i] + shifts[slice]) / scales[slice];
      dataPtr += distBetweenSets;
    }
  }
}


void normalizeTensorMinMax(Tensor* Y, int mode, const char* scale_file)
{
  MetricData* metrics = computeSliceMetrics(Y, mode, MAX+MIN);
  int sizeOfModeDim = Y->size(mode);
  double* scales = MemoryManager::safe_new_array<double>(sizeOfModeDim);
  double* shifts = MemoryManager::safe_new_array<double>(sizeOfModeDim);
  for(int i=0; i<sizeOfModeDim; i++) {
    scales[i] = metrics->getMaxData()[i] - metrics->getMinData()[i];
    shifts[i] = -metrics->getMinData()[i];
  }
  transformSlices(Y,mode,scales,shifts);
  if(scale_file) writeScaleShift(mode,sizeOfModeDim,scales,shifts,scale_file);
  MemoryManager::safe_delete_array<double>(scales,sizeOfModeDim);
  MemoryManager::safe_delete_array<double>(shifts,sizeOfModeDim);
  MemoryManager::safe_delete<MetricData>(metrics);
}

// \todo This function is not being tested at all
void normalizeTensorMax(Tensor* Y, int mode, const char* scale_file)
{
  Tucker::MetricData* metrics = computeSliceMetrics(Y, mode,
      Tucker::MIN + Tucker::MAX);
  int sizeOfModeDim = Y->size(mode);
  double* scales = Tucker::MemoryManager::safe_new_array<double>(sizeOfModeDim);
  double* shifts = Tucker::MemoryManager::safe_new_array<double>(sizeOfModeDim);
  for(int i=0; i<sizeOfModeDim; i++) {
    double scaleval = std::max(std::abs(metrics->getMinData()[i]),
        std::abs(metrics->getMaxData()[i]));
    scales[i] = scaleval;
    shifts[i] = 0;
  }
  transformSlices(Y,mode,scales,shifts);
  if(scale_file) writeScaleShift(mode,sizeOfModeDim,scales,shifts,scale_file);
  Tucker::MemoryManager::safe_delete_array<double>(scales,sizeOfModeDim);
  Tucker::MemoryManager::safe_delete_array<double>(shifts,sizeOfModeDim);
  Tucker::MemoryManager::safe_delete<Tucker::MetricData>(metrics);
}

void normalizeTensorStandardCentering(Tensor* Y, int mode, double stdThresh, const char* scale_file)
{
  MetricData* metrics = computeSliceMetrics(Y, mode, MEAN+VARIANCE);
  int sizeOfModeDim = Y->size(mode);
  double* scales = MemoryManager::safe_new_array<double>(sizeOfModeDim);
  double* shifts = MemoryManager::safe_new_array<double>(sizeOfModeDim);
  for(int i=0; i<sizeOfModeDim; i++) {
    scales[i] = sqrt(metrics->getVarianceData()[i]);
    shifts[i] = -metrics->getMeanData()[i];

    if(scales[i] < stdThresh) {
      scales[i] = 1;
    }
  }
  transformSlices(Y,mode,scales,shifts);
  if(scale_file) writeScaleShift(mode,sizeOfModeDim,scales,shifts,scale_file);
  MemoryManager::safe_delete_array<double>(scales,sizeOfModeDim);
  MemoryManager::safe_delete_array<double>(shifts,sizeOfModeDim);
  MemoryManager::safe_delete<MetricData>(metrics);
}

// \todo This function is not being tested
void writeScaleShift(const int mode, const int sizeOfModeDim, const double* scales,
    const double* shifts, const char* scale_file)
{
  std::ofstream outStream(scale_file);

  outStream << mode << std::endl;

  outStream << std::fixed << std::setprecision(16);
  for(int i=0; i<sizeOfModeDim; i++)
  {
    outStream << scales[i] << " " << shifts[i] << std::endl;
  }

  outStream.close();
}

// \todo This function is not being tested
void readTensorBinary(Tensor* Y, const char* filename)
{
  // Count the number of filenames
  std::ifstream inStream(filename);

  std::string temp;
  int nfiles = 0;
  while(inStream >> temp) {
    nfiles++;
  }

  inStream.close();

  if(nfiles == 1) {
    importTensorBinary(Y,temp.c_str());
  }
  else {
    int ndims = Y->N();
    if(nfiles != Y->size(ndims-1)) {
      std::ostringstream oss;
      oss << "Tucker::readTensorBinary(Tensor* Y, const char* filename: "
          << "The number of filenames you provided is "
          << nfiles << ", but the dimension of the tensor's last mode is "
          << Y->size(ndims-1);

      throw std::runtime_error(oss.str());
    }
    importTimeSeries(Y,filename);
  }
}


Tensor* importTensor(const char* filename)
{
  // Open file
  std::ifstream ifs;
  ifs.open(filename);
  assert(ifs.is_open());

  // Read the type of object
  // If the type is not "tensor", that's bad
  std::string tensorStr;
  ifs >> tensorStr;
  assert(tensorStr == "tensor" || tensorStr == "matrix");

  // Read the number of dimensions
  int ndims;
  ifs >> ndims;

  // Create a SizeArray of that length
  SizeArray sz(ndims);

  // Read the dimensions
  for(int i=0; i<ndims; i++) {
    ifs >> sz[i];
  }

  // Create a tensor using that SizeArray
  Tensor* t = MemoryManager::safe_new<Tensor>(sz);

  // Read the entries of the tensor
  size_t numEntries = sz.prod();
  double * data = t->data();
  for(size_t i=0; i<numEntries; i++) {
    ifs >> data[i];
  }

  // Close the file
  ifs.close();

  // Return the tensor
  return t;
}

void importTensorBinary(Tensor* t, const char* filename)
{
  // Get the maximum file size we can read
  const std::streamoff MAX_OFFSET =
      std::numeric_limits<std::streamoff>::max();
//  std::cout << "The maximum file size is " << MAX_OFFSET << " bytes\n";

  // Open file
  std::ifstream ifs;
  ifs.open(filename, std::ios::in | std::ios::binary);
  assert(ifs.is_open());

  // Get the size of the file
  std::streampos begin, end, size;
  begin = ifs.tellg();
  ifs.seekg(0, std::ios::end);
  end = ifs.tellg();
  size = end - begin;
//  std::cout << "Reading " << size << " bytes...\n";

  // Assert that this size is consistent with the number of tensor entries
  size_t numEntries = t->getNumElements();
  assert(size == numEntries*sizeof(double));

  // Read the file
  double* data = t->data();
  ifs.seekg(0, std::ios::beg);
  ifs.read((char*)data,size);

  // Close the file
  ifs.close();
}

// \todo This function never gets tested
void importTimeSeries(Tensor* Y, const char* filename)
{
   // Open the file
   std::ifstream ifs;
   ifs.open(filename);

   // Define data layout parameters
   int ndims = Y->N();

   int nsteps = Y->size(ndims-1);
   double* dataPtr = Y->data();
   size_t count = Y->size().prod(0,ndims-2);
   assert(count <= std::numeric_limits<int>::max());

   for(int step=0; step<nsteps; step++) {
     std::string stepFilename;
     ifs >> stepFilename;
     std::cout << "Reading file " << stepFilename << std::endl;

     std::ifstream bifs;
     bifs.open(stepFilename.c_str(), std::ios::in | std::ios::binary);
     assert(bifs.is_open());

     // Get the size of the file
     std::streampos begin, end, size;
     begin = bifs.tellg();
     bifs.seekg(0, std::ios::end);
     end = bifs.tellg();
     size = end - begin;

     // Assert that this size is consistent with the number of tensor entries
     size_t numEntries = Y->getNumElements();
     assert(size == numEntries*sizeof(double));

     // Read the file
     bifs.seekg(0, std::ios::beg);
     bifs.read((char*)dataPtr,size);

     bifs.close();

     // Increment the pointer
     dataPtr += count;
   }

   ifs.close();
}


Matrix* importMatrix(const char* filename)
{
  // Open file
  std::ifstream ifs;
  ifs.open(filename);
  assert(ifs.is_open());

  // Read the type of object
  // If the type is not "tensor", that's bad
  std::string tensorStr;
  ifs >> tensorStr;
  assert(tensorStr == "tensor" || tensorStr == "matrix");

  // Read the number of dimensions
  int ndims;
  ifs >> ndims;
  assert(ndims == 2);

  // Read the dimensions
  int nrows, ncols;
    ifs >> nrows >> ncols;

  // Create a matrix
  Matrix* m = MemoryManager::safe_new<Matrix>(nrows,ncols);

  // Read the entries of the tensor
  size_t numEntries = nrows*ncols;
  double * data = m->data();
  for(size_t i=0; i<numEntries; i++) {
    ifs >> data[i];
  }

  // Close the file
  ifs.close();

  // Return the tensor
  return m;
}


SparseMatrix* importSparseMatrix(const char* filename)
{
  // Open file
  std::ifstream ifs;
  ifs.open(filename);
  assert(ifs.is_open());

  // Read the type of object
  // If the type is not "sptensor", that's bad
  std::string tensorStr;
  ifs >> tensorStr;
  assert(tensorStr == "sptensor");

  // Read the number of dimensions
  int ndims;
  ifs >> ndims;
  assert(ndims == 2);

  // Read the dimensions
  int nrows, ncols, nnz;
    ifs >> nrows >> ncols >> nnz;

  // Create a matrix
  SparseMatrix* m = MemoryManager::safe_new<SparseMatrix>(nrows,ncols,nnz);

  // Read the entries of the tensor
  int* rows = m->rows();
  int* cols = m->cols();
  double* vals = m->vals();
  for(size_t i=0; i<nnz; i++) {
    ifs >> rows[i] >> cols[i] >> vals[i];
    rows[i]--;
    cols[i]--;
  }

  // Close the file
  ifs.close();

  // Return the sparse matrix
  return m;
}

// \todo This function never gets tested
void writeTensorBinary(const Tensor* Y, const char* filename)
{
  // Count the number of filenames
   std::ifstream inStream(filename);

   std::string temp;
   int nfiles = 0;
   while(inStream >> temp) {
     nfiles++;
   }

   inStream.close();

   if(nfiles == 1) {
     exportTensorBinary(Y,temp.c_str());
   }
   else {
     int ndims = Y->N();
     if(nfiles != Y->size(ndims-1)) {
       std::ostringstream oss;
       oss << "Tucker::writeTensorBinary(const Tensor* Y, const char* filename: "
           << "The number of filenames you provided is "
           << nfiles << ", but the dimension of the tensor's last mode is "
           << Y->size(ndims-1);

       throw std::runtime_error(oss.str());
     }
     exportTimeSeries(Y,filename);
   }
}

// \todo This function never gets tested
void exportTensor(const Tensor* Y, const char* filename)
{
  // Open the file
  std::ofstream ofs;
  ofs.open(filename);

  // Write the type of object
  ofs << "tensor\n";

  // Write the number of dimensions of the tensor
  int ndims = Y->size().size();
  ofs << ndims << std::endl;

  // Write the size of each dimension
  for(int i=0; i<ndims; i++) {
    ofs << Y->size(i) << " ";
  }
  ofs << std::endl;

  // Write the elements of the tensor
  size_t numEntries = Y->size().prod();
  const double* data = Y->data();
  for(size_t i=0; i<numEntries; i++) {
    ofs << data[i] << std::endl;
  }

  // Close the file
  ofs.close();
}

void exportTensorBinary(const Tensor* Y, const char* filename)
{
  // Get the maximum file size we can write
  const std::streamoff MAX_OFFSET =
      std::numeric_limits<std::streamoff>::max();
//  std::cout << "The maximum file size is " << MAX_OFFSET << " bytes\n";

  // Determine how many bytes we are writing
  size_t numEntries = Y->getNumElements();
//  std::cout << "Writing " << numEntries*sizeof(double) << " bytes...\n";

  // Open file
  std::ofstream ofs;
  ofs.open(filename, std::ios::out | std::ios::binary);
  assert(ofs.is_open());

  // Write the file
  const double* data = Y->data();
  ofs.write((char*)data,numEntries*sizeof(double));

  // Close the file
  ofs.close();
}

// \todo This function never gets tested
void exportTimeSeries(const Tensor* Y, const char* filename)
{
  // Open the file
  std::ifstream ifs;
  ifs.open(filename);

  // Determine how many bytes we are writing per file
  int N = Y->N();
  size_t numEntriesPerTimestep = Y->size().prod(0,N-2);

  int nsteps = Y->size(N-1);
  size_t offset = 0;
  for(int step=0; step<nsteps; step++) {
    std::string stepFilename;
    ifs >> stepFilename;
    std::cout << "Writing file " << stepFilename << std::endl;

    std::ofstream ofs;
    ofs.open(stepFilename.c_str(), std::ios::out | std::ios::binary);
    assert(ofs.is_open());

    const double* data = Y->data() + offset;
    ofs.write((char*)data,numEntriesPerTimestep*sizeof(double));
    ofs.close();

    offset += numEntriesPerTimestep;
  }
}


void premultByDiag(const Vector* diag, Matrix* mat)
{
  double* mydata = mat->data();
  int myrows = mat->nrows();
  int mycols = mat->ncols();

  assert(myrows == diag->nrows());

  for(int r=0; r<myrows; r++) {
    for(int c=0; c<mycols; c++) {
      mydata[r+c*myrows] *= (*diag)[r];
    }
  }
}

} // end namespace Tucker
