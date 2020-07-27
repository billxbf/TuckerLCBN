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
 * \brief Defines a parallel distribution of MPI processes on a Cartesian grid
 *
 * @author Alicia Klinvex
 */

#include "TuckerMPI_ProcessorGrid.hpp"
#include "Tucker_Util.hpp"

namespace TuckerMPI {

ProcessorGrid::ProcessorGrid(const Tucker::SizeArray& sz,
    const MPI_Comm& comm) :
        size_(sz.size()),
        squeezed_(false),
        cartComm_squeezed_(MPI_COMM_NULL),
        rowcomms_squeezed_(0),
        colcomms_squeezed_(0)

{
  int ndims = sz.size();

  for(int i=0; i<ndims; i++) {
    size_[i] = sz[i];
  }

  // Get the number of MPI processes
  int nprocs;
  MPI_Comm_size(comm, &nprocs);

  // Assert that the number of processes is consistent
  int nprocsRequested = 1;
  for(int i=0; i<ndims; i++) {
    nprocsRequested *= sz[i];
  }
  if(nprocsRequested != nprocs) {
    std::cerr << "ERROR in ProcessorGrid constructor: the processor grid "
        << "supplied is inconsistent with the total number of processes\n";
  }

  // Create a virtual topology MPI communicator
  int* periods = Tucker::MemoryManager::safe_new_array<int>(ndims);
  for(int i=0; i<ndims; i++) periods[i] = 1;
  int reorder = 0;
  MPI_Cart_create(comm, ndims, (int*)sz.data(), periods,
      reorder, &cartComm_);
  Tucker::MemoryManager::safe_delete_array<int>(periods,ndims);

  // Allocate memory for subcommunicators
  rowcomms_ = Tucker::MemoryManager::safe_new_array<MPI_Comm>(ndims);
  colcomms_ = Tucker::MemoryManager::safe_new_array<MPI_Comm>(ndims);

  // Get the subcommunicators
  int* remainDims = Tucker::MemoryManager::safe_new_array<int>(ndims);
  for(int i=0; i<ndims; i++) remainDims[i] = 0;
  for(int i=0; i<ndims; i++)
  {
    remainDims[i] = 1;
    MPI_Cart_sub(cartComm_, remainDims, &(colcomms_[i]));
    remainDims[i] = 0;
  }

  for(int i=0; i<ndims; i++) remainDims[i] = 1;
  for(int i=0; i<ndims; i++)
  {
    remainDims[i] = 0;
    MPI_Cart_sub(cartComm_, remainDims, &(rowcomms_[i]));
    remainDims[i] = 1;
  }
  Tucker::MemoryManager::safe_delete_array<int>(remainDims,ndims);
}


ProcessorGrid::~ProcessorGrid()
{
  int ndims = size_.size();
  int finalized;
  MPI_Finalized(&finalized);

  if(!finalized) {
    for(int i=0; i<ndims; i++) {
      MPI_Comm_free(rowcomms_+i);
      MPI_Comm_free(colcomms_+i);
    }
  }
  Tucker::MemoryManager::safe_delete_array<MPI_Comm>(rowcomms_,ndims);
  Tucker::MemoryManager::safe_delete_array<MPI_Comm>(colcomms_,ndims);
  if(squeezed_) {
    if(!finalized) {
      MPI_Comm_free(&cartComm_squeezed_);
      for(int i=0; i<ndims; i++) {
        MPI_Comm_free(rowcomms_squeezed_+i);
        MPI_Comm_free(colcomms_squeezed_+i);
      }
    }
    Tucker::MemoryManager::safe_delete_array<MPI_Comm>(rowcomms_squeezed_,ndims);
    Tucker::MemoryManager::safe_delete_array<MPI_Comm>(colcomms_squeezed_,ndims);
  }
}

const MPI_Comm& ProcessorGrid::getComm(bool squeezed) const
{
  if(squeezed && squeezed_) {
    return cartComm_squeezed_;
  }
  return cartComm_;
}


void ProcessorGrid::getCoordinates(int* coords) const
{
  int globalRank;
  MPI_Comm_rank(cartComm_, &globalRank);
  getCoordinates(coords, globalRank);
}

void ProcessorGrid::getCoordinates(int* coords, int globalRank) const
{
  int ndims = size_.size();
  MPI_Cart_coords(cartComm_, globalRank, ndims, coords);
}


const MPI_Comm& ProcessorGrid::getRowComm(const int d, bool squeezed) const
{
  if(squeezed && squeezed_) {
    return rowcomms_squeezed_[d];
  }
  return rowcomms_[d];
}


const MPI_Comm& ProcessorGrid::getColComm(const int d, bool squeezed) const
{
  if(squeezed && squeezed_) {
    return colcomms_squeezed_[d];
  }
  return colcomms_[d];
}


int ProcessorGrid::getRank(const int* coords) const
{
  int rank;
  MPI_Cart_rank(cartComm_,(int*)coords,&rank);
  return rank;
}


int ProcessorGrid::getNumProcs(int d, bool squeezed) const
{
  int nprocs;
  if(squeezed && squeezed_) {
    MPI_Comm_size(colcomms_squeezed_[d],&nprocs);
  }
  else {
    MPI_Comm_size(colcomms_[d],&nprocs);
  }
  return nprocs;
}

void ProcessorGrid::squeeze(const Tucker::SizeArray& sz, const MPI_Comm& comm)
{
  squeezed_ = true;
  int ndims = size_.size();

  // Get the number of MPI processes
  int nprocs;
  MPI_Comm_size(comm, &nprocs);

  // Assert that the number of processes is consistent
  int nprocsRequested = 1;
  for(int i=0; i<ndims; i++) {
    nprocsRequested *= sz[i];
  }
  if(nprocsRequested != nprocs) {
    std::cerr << "ERROR in ProcessorGrid::squeeze: the processor grid "
        << "supplied is inconsistent with the total number of processes\n";
  }

  // Create a virtual topology MPI communicator
  int* periods = Tucker::MemoryManager::safe_new_array<int>(ndims);
  for(int i=0; i<ndims; i++) periods[i] = 1;
  int reorder = 0;
  MPI_Cart_create(comm, ndims, (int*)sz.data(), periods,
      reorder, &cartComm_squeezed_);
  Tucker::MemoryManager::safe_delete_array<int>(periods,ndims);

  // Allocate memory for subcommunicators
  rowcomms_squeezed_ = Tucker::MemoryManager::safe_new_array<MPI_Comm>(ndims);
  colcomms_squeezed_ = Tucker::MemoryManager::safe_new_array<MPI_Comm>(ndims);

  // Get the subcommunicators
  int* remainDims = Tucker::MemoryManager::safe_new_array<int>(ndims);
  for(int i=0; i<ndims; i++) remainDims[i] = 0;
  for(int i=0; i<ndims; i++)
  {
    remainDims[i] = 1;
    MPI_Cart_sub(cartComm_squeezed_, remainDims, &(colcomms_squeezed_[i]));
    remainDims[i] = 0;
  }

  for(int i=0; i<ndims; i++) remainDims[i] = 1;
  for(int i=0; i<ndims; i++)
  {
    remainDims[i] = 0;
    MPI_Cart_sub(cartComm_squeezed_, remainDims, &(rowcomms_squeezed_[i]));
    remainDims[i] = 1;
  }
  Tucker::MemoryManager::safe_delete_array<int>(remainDims,ndims);
}

} /* namespace Tucker */
