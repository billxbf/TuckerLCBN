/*
 * SparseMatrix.cpp
 *
 *  Created on: Mar 14, 2017
 *      Author: amklinv
 */


#include "Tucker_SparseMatrix.hpp"
#include "Tucker_Util.hpp"
#include <cassert>

namespace Tucker {

SparseMatrix::SparseMatrix(const int nrows, const int ncols, const int nnz)
{
  nrows_ = nrows;
  ncols_ = ncols;
  nnz_ = nnz;

  rows_ = Tucker::MemoryManager::safe_new_array<int>(nnz);
  cols_ = Tucker::MemoryManager::safe_new_array<int>(nnz);
  vals_ = Tucker::MemoryManager::safe_new_array<double>(nnz);
}

SparseMatrix::~SparseMatrix() {
  Tucker::MemoryManager::safe_delete_array<int>(rows_,nnz_);
  Tucker::MemoryManager::safe_delete_array<int>(cols_,nnz_);
  Tucker::MemoryManager::safe_delete_array<double>(vals_,nnz_);
}

Matrix* SparseMatrix::multiply(const Matrix* m, bool transp)
{
  assert(m != NULL);

  // Ensure that the dimensions are consistent
  if(transp) {
    assert(m->nrows() == nrows_);
  }
  else {
    assert(m->nrows() == ncols_);
  }

  // Allocate memory for result
  Matrix* result;
  if(transp) {
    result = Tucker::MemoryManager::safe_new<Matrix>(ncols_,m->ncols());
  }
  else {
    result = Tucker::MemoryManager::safe_new<Matrix>(nrows_,m->ncols());
  }

  // Initialize result to 0
  result->initialize();

  double* resultData = result->data();
  const double* mData = m->data();
  for(int i=0; i<nnz_; i++) {
    for(int c=0; c<m->ncols(); c++) {
      if(transp) {
        // result(col(i),c) += val(i) m(row(i),c)
        resultData[cols_[i]+ncols_*c] += vals_[i]*mData[rows_[i]+m->nrows()*c];
      }
      else {
        // result(row(i),c) += val(i) m(col(i),c)
        resultData[rows_[i]+nrows_*c] += vals_[i]*mData[cols_[i]+m->nrows()*c];
      }
    } // end for(int c=0; c<m->ncols(); c++)
  } // end for(int i=0; i<nnz_; i++)

  return result;
}

Vector* SparseMatrix::multiply(const Vector* v, bool transp)
{
  assert(v != NULL);

  // Ensure that the dimensions are consistent
  if(transp) {
    assert(v->nrows() == nrows_);
  }
  else {
    assert(v->nrows() == ncols_);
  }

  // Allocate memory for result
  Vector* result;
  if(transp) {
    result = Tucker::MemoryManager::safe_new<Vector>(ncols_);
  }
  else {
    result = Tucker::MemoryManager::safe_new<Vector>(nrows_);
  }

  // Initialize result to 0
  result->initialize();

  double* resultData = result->data();
  const double* vData = v->data();
  for(int i=0; i<nnz_; i++) {
    if(transp) {
      // result(col(i)) += val(i) m(row(i))
      resultData[cols_[i]] += vals_[i]*vData[rows_[i]];
    }
    else {
      // result(row(i),c) += val(i) m(col(i),c)
      resultData[rows_[i]] += vals_[i]*vData[cols_[i]];
    }
  } // end for(int i=0; i<nnz_; i++)

  return result;
}

} /* namespace Tucker */
