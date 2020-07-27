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
 * \brief Contains the parallel Matrix class
 *
 * @author: Alicia Klinvex
 */

#ifndef MATRIX_MPI_HPP_
#define MATRIX_MPI_HPP_

#include "TuckerMPI_Map.hpp"
#include "Tucker_Matrix.hpp"
#include "mpi.h"

namespace TuckerMPI {
/** \brief Parallel matrix class
 *
 * Used by the new Gram routine
 *
 * \todo Block row distribution is never tested
 */
class Matrix {
public:
  /** \brief Constructor
   *
   * \param nrows Number of global rows
   * \param ncols Number of global columns
   * \param comm MPI communicator to distribute this matrix over
   * \param isBlockRow Describes whether this matrix has a block row
   * or block column distribution.
   */
  Matrix(int nrows, int ncols, const MPI_Comm& comm, bool isBlockRow = true);

  //! Destructor
  virtual ~Matrix();

  /** \brief Returns the locally owned portion of the matrix as a sequential
   * Tucker::Matrix pointer
   */
  Tucker::Matrix* getLocalMatrix();

  /** \brief Returns the locally owned portion of the matrix as a sequential
   * Tucker::Matrix pointer
   */
  const Tucker::Matrix* getLocalMatrix() const;

  //! Returns the number of entries owned by the calling process
  size_t getLocalNumEntries() const;

  //! Returns the global number of rows
  int getGlobalNumRows() const;

  //! Returns the number of rows owned by the calling process
  int getLocalNumRows() const;

  //! Returns the global number of columns
  int getGlobalNumCols() const;

  //! Returns the number of columns owned by the calling process
  int getLocalNumCols() const;

  /** \brief Returns the map describing the data distribution
   *
   * If the matrix has a block row distribution, this is the row map.
   * Otherwise, it is the column map.
   */
  const Map* getMap() const;

  //! Prints the entries of the matrix as a single dimensional array
  void print() const;
private:
  /// @cond EXCLUDE
  Matrix(const Matrix& m);
  /// @endcond

  //! The locally owned portion of the matrix
  Tucker::Matrix* localMatrix_;
  //! The map describing the distribution
  Map* map_;
  //! The global number of rows
  int globalRows_;
  //! The global number of columns
  int globalCols_;
  //! The MPI communicator
  const MPI_Comm& comm_;
};

} /* namespace TuckerMPI */
#endif /* MATRIX_MPI_HPP_ */
