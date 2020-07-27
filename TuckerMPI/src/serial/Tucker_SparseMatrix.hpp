/*
 * SparseMatrix.hpp
 *
 *  Created on: Mar 14, 2017
 *      Author: amklinv
 */

#ifndef SERIAL_SPARSEMATRIX_HPP_
#define SERIAL_SPARSEMATRIX_HPP_

#include "Tucker_Matrix.hpp"
#include "Tucker_Vector.hpp"

namespace Tucker {

class SparseMatrix {
public:
  SparseMatrix(const int nrows, const int ncols, const int nnz);
  ~SparseMatrix();

  Matrix* multiply(const Matrix* m, bool transp=false);
  Vector* multiply(const Vector* v, bool transp=false);

  int* rows() { return rows_; }
  int* cols() { return cols_; }
  double* vals() { return vals_; }
  int nrows() const { return nrows_; }
  int ncols() const { return ncols_; }

private:
  /// @cond EXCLUDE
  SparseMatrix(const SparseMatrix&);
  /// @endcond

  int nrows_;
  int ncols_;
  int nnz_;
  int* rows_;
  int* cols_;
  double* vals_;
};

} /* namespace Tucker */

#endif /* SERIAL_TUCKER_SPARSEMATRIX_HPP_ */
