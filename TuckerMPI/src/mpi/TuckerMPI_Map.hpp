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

#ifndef MAP_HPP_
#define MAP_HPP_

#include "mpi.h"
#include "Tucker_SizeArray.hpp"

namespace TuckerMPI {

/** \brief Describes a 1D distribution over a set of MPI processes.
 *
 *  Data is divided into sequential blocks and assigned to
 *  the MPI processes accordingly.  For instance, the set of
 *  integers 0..7 would be divided amongst 3 processes as
 *  follows:
 *  - process 0 owns {0,1,2}
 *  - process 1 owns {3,4,5}
 *  - process 2 owns {6,7}
 *
 *  In the event that the size of the dataset (n) cannot be
 *  divided evenly by np (number of processors), the first
 *  processors will receive ceil(n/np) objects, and the last
 *  ones will receive floor(n/np) objects.
 */
class Map {
public:
  /** \brief Constructor
   *
   * \param globalNumEntries the size of the dataset
   * \param comm a communicator defining the set of MPI processes
   * the data will be divided among
   */
  Map(int globalNumEntries, const MPI_Comm& comm);

  //! Destructor
  virtual ~Map();

  /** \brief Given a local index, returns the corresponding global index
   *
   * If the local index is invalid, returns -1
   */
  int getLocalIndex(int globalIndex) const;

  /** \brief Given a global index, returns the corresponding local index
   *
   * If the global index is invalid, returns -1
   */
  int getGlobalIndex(int localIndex) const;

  //! Returns the number of entries owned by the calling process
  int getLocalNumEntries() const;

  //! Returns the total number of entries owned by ALL processes
  int getGlobalNumEntries() const;

  /** \brief Returns the maximum number of entries owned by any process
   *
   * This function is called by oldGram, which uses it to allocate
   * a single buffer for the transferred data which can be reused
   * again and again, rather than allocating a new buffer each time
   * data needs to be received.
   */
  int getMaxNumEntries() const;

  /** \brief Returns the number of entries owned by a particular rank
   *
   * This function is useful for MPI communication, to determine how
   * much data will be received from another process.
   */
  int getNumEntries(int rank) const;

  /** \brief Returns the offset corresponding to a particular rank
   *
   * This function is useful for MPI collective communication calls
   * where each MPI process owns a different amount of data,
   * such as MPI_Alltoallv.  For those functions, you must pass
   * an integer array of length np containing the amount of data each
   * MPI process will be sending, along with the displacements (or
   * offsets) for placing the data in the receive buffer.
   */
  int getOffset(int rank) const;

  /** \brief Returns an array of the number of entries owned by each process
   *
   * Used by driver.cpp
   */
  const Tucker::SizeArray* getNumElementsPerProc() const;

  /** \brief Returns the offset array
   *
   * Used by driver.cpp
   */
  const Tucker::SizeArray* getOffsets() const;

  //! Returns the MPI communicator
  const MPI_Comm& getComm() const;

  /** \brief Removes empty MPI processes from the communicator
   *
   * Also removes them from the arrays that store the number of
   * elements owned by each process and the offsets.
   */
  void removeEmptyProcs();
private:
  /// @cond EXCLUDE
  Map(const Map& m);
  /// @endcond

  //! MPI communicator
  MPI_Comm comm_;
  //! Number of elements owned by each process
  Tucker::SizeArray* numElementsPerProc_;
  //! Offset/displacement array
  Tucker::SizeArray* offsets_;
  //! First index owned by this MPI process
  int indexBegin_;
  //! Last index owned by this MPI process
  int indexEnd_;
  //! Number of entries owned by this MPI process
  int localNumEntries_;
  //! Total number of entries
  int globalNumEntries_;

  bool removedEmptyProcs_;
};

} /* namespace TuckerMPI */
#endif /* MAP_HPP_ */
