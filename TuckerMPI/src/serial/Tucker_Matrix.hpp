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
 * \brief The current storage for dense matrices/factors
 *
 * @author Alicia Klinvex
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cassert>
#include <iostream>
#include "Tucker_Tensor.hpp"

/// @cond EXCLUDE
// Copy from one array to another
extern "C" void dcopy_(const int*, const double*, const int*,
    double*, const int*);
/// @endcond

namespace Tucker {

/** \brief A sequential dense matrix
 *
 * Since a matrix is a 2-dimensional tensor, this class
 * inherits from the tensor class.
 */
class Matrix : public Tensor {
public:
  /** \brief Constructor
   * \param[in] nrows Number of rows
   * \param[in] ncols Number of columns
   */
  Matrix(const int nrows, const int ncols) :
    Tensor(nrows,ncols)
  {

  }

  /// Returns the number of rows
  int nrows() const
  {
    return I_[0];
  }

  /// Returns the number of columns
  int ncols() const
  {
    return I_[1];
  }

  Matrix* getSubmatrix(const int rbegin, const int rend) const
  {
    const int ONE = 1;
    int new_nrows = rend-rbegin+1;
    int old_nrows = nrows();
    int myncols = ncols();
    Matrix* newMat = MemoryManager::safe_new<Matrix>(new_nrows,myncols);

    for(int c=0; c<myncols; c++)
    {
      dcopy_(&new_nrows, data()+c*old_nrows+rbegin, &ONE,
          newMat->data()+c*new_nrows, &ONE);
    }

    return newMat;
  }

protected:
  Matrix(const int nrows) :
    Tensor(nrows)
  {

  }

private:
  /// @cond EXCLUDE
  // Disables the copy constructor
  Matrix(const Matrix& m);
  /// @endcond
};

} // end of namespace Tucker

#endif /* MATRIX_HPP_ */
