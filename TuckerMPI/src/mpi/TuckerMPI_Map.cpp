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
 * \brief Describes a 1D distribution over a set of MPI processes
 *
 * @author Alicia Klinvex
 */

#include "TuckerMPI_Map.hpp"
#include <cassert>
#include <limits>

namespace TuckerMPI {

Map::Map(int globalNumEntries, const MPI_Comm& comm) :
  comm_(comm),
  globalNumEntries_(globalNumEntries),
  removedEmptyProcs_(false)
{
  // Get the number of MPI processes
  int myRank, nprocs;
  MPI_Comm_rank(comm,&myRank);
  MPI_Comm_size(comm,&nprocs);

  // Assert that the global number of entries is bigger
  // than the number of MPI processes
//  assert(globalNumEntries > nprocs);

  // Determine the number of entries owned by each process
  numElementsPerProc_ = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nprocs);

  for(int rank=0; rank<nprocs; rank++) {
    (*numElementsPerProc_)[rank] = globalNumEntries/nprocs;
    if(rank < globalNumEntries%nprocs)
      (*numElementsPerProc_)[rank]++;
  }

  // Determine the row offsets for each process
  offsets_ = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nprocs+1);

  (*offsets_)[0] = 0;
  for(int rank=1; rank<=nprocs; rank++) {
    (*offsets_)[rank] = (*offsets_)[rank-1] + (*numElementsPerProc_)[rank-1];
  }

  // Determine the starting and ending indices for THIS process
  indexBegin_ = (*offsets_)[myRank];
  indexEnd_ = (*offsets_)[myRank+1]-1;

  localNumEntries_ = 1 + indexEnd_ - indexBegin_;
}

Map::~Map()
{
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(numElementsPerProc_);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(offsets_);

  int finalized;
  MPI_Finalized(&finalized);
  if(removedEmptyProcs_ && getLocalNumEntries() > 0 && !finalized) {
    MPI_Comm_free(&comm_);
  }
}

int Map::getLocalIndex(int globalIndex) const
{
  if(globalIndex > indexEnd_ || globalIndex < indexBegin_) {
    return -1;
  }
  return globalIndex - indexBegin_;
}

int Map::getGlobalIndex(int localIndex) const
{
  if(localIndex < 0 || localIndex >= localNumEntries_) {
    return -1;
  }
  return indexBegin_+localIndex;
}

int Map::getLocalNumEntries() const
{
  return localNumEntries_;
}

int Map::getGlobalNumEntries() const
{
  return globalNumEntries_;
}


int Map::getMaxNumEntries() const
{
  int maxNumEntries = 0;
  for(int i=0; i<numElementsPerProc_->size(); i++) {
    if((*numElementsPerProc_)[i] > maxNumEntries) {
      maxNumEntries = (*numElementsPerProc_)[i];
    }
  }
  return maxNumEntries;
}


int Map::getNumEntries(int rank) const
{
  return (*numElementsPerProc_)[rank];
}

int Map::getOffset(int rank) const
{
  return (*offsets_)[rank];
}

const MPI_Comm& Map::getComm() const
{
  return comm_;
}

// \todo This function doesn't seem to get used, and it's definitely not tested
const Tucker::SizeArray* Map::getNumElementsPerProc() const
{
  return numElementsPerProc_;
}

// \todo This function doesn't seem to get used, and it's definitely not tested
const Tucker::SizeArray* Map::getOffsets() const
{
  return offsets_;
}

void Map::removeEmptyProcs()
{
  // Determine which processes are empty
  std::vector<int> emptyProcs;
  for(int rank=0; rank<numElementsPerProc_->size(); rank++) {
    if((*numElementsPerProc_)[rank] == 0) {
      emptyProcs.push_back(rank);
    }
  }

  // If none are empty, there's nothing to be done
  if(emptyProcs.size() == 0) {
    return;
  }

  // Remove those from numElementsPerProc
  size_t newNumProcs = numElementsPerProc_->size() - emptyProcs.size();
  size_t i=0;
  int src=0;
  assert(newNumProcs <= std::numeric_limits<int>::max());
  Tucker::SizeArray* newSize = Tucker::MemoryManager::safe_new<Tucker::SizeArray>((int)newNumProcs);
  for(int dest=0; dest<(int)newNumProcs; dest++) {
    while(i < emptyProcs.size() && src == emptyProcs[i]) {
      src++;
      i++;
    }
    (*newSize)[dest] = (*numElementsPerProc_)[src];
    src++;
  }
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(numElementsPerProc_);
  numElementsPerProc_ = newSize;

  // Remove them from offsets too
  Tucker::SizeArray* newOffsets = Tucker::MemoryManager::safe_new<Tucker::SizeArray>((int)newNumProcs);
  i=0;
  src=0;
  for(int dest=0; dest<(int)newNumProcs; dest++) {
    while(i < emptyProcs.size() && src == emptyProcs[i]) {
      src++;
      i++;
    }
    (*newOffsets)[dest] = (*offsets_)[src];
    src++;
  }
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(offsets_);
  offsets_ = newOffsets;

  assert(emptyProcs.size() <= std::numeric_limits<int>::max());

  // Remove them from the communicator too
  MPI_Group old_group, new_group;
  MPI_Comm new_comm;
  MPI_Comm_group(comm_, &old_group);
  MPI_Group_excl (old_group, (int)emptyProcs.size(),
      emptyProcs.data(), &new_group);
  MPI_Comm_create (comm_, new_group, &new_comm);
  comm_ = new_comm;
  removedEmptyProcs_ = true;
}

} /* namespace TuckerMPI */
