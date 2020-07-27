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
 * \author Alicia Klinvex
 */

#include "TuckerMPI_Matrix.hpp"

namespace TuckerMPI {

// TODO: isBlockRow currently does NOTHING
Matrix::Matrix(int nrows, int ncols, const MPI_Comm& comm, bool isBlockRow) :
  globalRows_(nrows),
  globalCols_(ncols),
  comm_(comm)
{
  // Get the communicator information
  int rank, nprocs;
  MPI_Comm_rank(comm_,&rank);
  MPI_Comm_size(comm_,&nprocs);

  // Create a map
  if(isBlockRow) {
    map_ = Tucker::MemoryManager::safe_new<Map>(nrows,comm);
  }
  else {
    map_ = Tucker::MemoryManager::safe_new<Map>(ncols,comm);
  }

  // Get the local number of rows and columns
  int localRows, localCols;
  if(isBlockRow) {
    localRows = map_->getLocalNumEntries();
    localCols = ncols;
  }
  else {
    localRows = nrows;
    localCols = map_->getLocalNumEntries();
  }

  // Create the local portion of the matrix
  localMatrix_ = Tucker::MemoryManager::safe_new<Tucker::Matrix>(localRows,localCols);
}


Matrix::~Matrix()
{
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(localMatrix_);
  Tucker::MemoryManager::safe_delete<Map>(map_);
}


Tucker::Matrix* Matrix::getLocalMatrix()
{
  return localMatrix_;
}


const Tucker::Matrix* Matrix::getLocalMatrix() const
{
  return localMatrix_;
}


size_t Matrix::getLocalNumEntries() const
{
  return localMatrix_->getNumElements();
}


int Matrix::getGlobalNumRows() const
{
  return globalRows_;
}


int Matrix::getLocalNumRows() const
{
  return localMatrix_->nrows();
}


int Matrix::getGlobalNumCols() const
{
  return globalCols_;
}


int Matrix::getLocalNumCols() const
{
  return localMatrix_->ncols();
}


const Map* Matrix::getMap() const
{
  return map_;
}


void Matrix::print() const
{
  localMatrix_->print();
}

} /* namespace TuckerMPI */
