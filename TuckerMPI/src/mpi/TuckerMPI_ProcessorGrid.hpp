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

#ifndef PROCESSORGRID_HPP_
#define PROCESSORGRID_HPP_

#include "mpi.h"
#include "Tucker_SizeArray.hpp"

namespace TuckerMPI {

/** \brief Defines a parallel distribution of MPI processes on a Cartesian grid
 *
 * \todo This and much of the other MPI code assumes that the global
 * communicator is MPI_COMM_WORLD.  Is that a safe assumption to make?
 */
class ProcessorGrid {
public:
  /** \brief Constructor
   *
   * Creates a set of MPI communicators using the cool virtual topology
   * routines in the MPI standard
   * (https://computing.llnl.gov/tutorials/mpi/#Virtual_Topologies).
   * \param[in] sz Defines the number of processes in each direction
   * \note This class stores a COPY of \a sz
   * \param[in] comm MPI communicator
   */
  ProcessorGrid(const Tucker::SizeArray& sz, const MPI_Comm& comm);

  //! Destructor
  ~ProcessorGrid();

  /**
   * \brief Returns the MPI communicator
   *
   * \param[in] squeezed If false, return the entire communicator.
   * If true, return only a subset of MPI processes: the ones that
   * were not eliminated by the squeeze function.
   */
  const MPI_Comm& getComm(bool squeezed) const;

  /** \brief Returns the row communicator for dimension d
   *
   * Let the calling process be mapped to location
   * (s_0 ... s_{d-1} s_d s_{d+1} ... s_{N-1})
   * The row communicator contains the processes mapped to
   * (* ... * s_d * ... *)
   *
   * \param[in] d Dimension
   * \param[in] squeezed If false, return the entire communicator.
   * If true, return only a subset of MPI processes: the ones that
   * were not eliminated by the squeeze function.
   */
  const MPI_Comm& getRowComm(const int d, bool squeezed) const;

  /** \brief Returns the row communicator for dimension d
   *
   * Let the calling process be mapped to location
   * (s_0 ... s_{d-1} s_d s_{d+1} ... s_{N-1})
   * The row communicator contains the processes mapped to
   * (s_0 ... s_{d-1} * s_{d+1} ... s_{N-1})
   *
   * \param[in] d Dimension
   * \param[in] squeezed If false, return the entire communicator.
   * If true, return only a subset of MPI processes: the ones that
   * were not eliminated by the squeeze function.
   */
  const MPI_Comm& getColComm(const int d, bool squeezed) const;

  /** \brief Returns the cartesian coordinates of the calling process
   * in the grid
   *
   * This function is not required by any of the Tucker code
   * and is merely provided for the sake of debugging.
   */
  void getCoordinates(int* coords) const;

  void getCoordinates(int* coords, int globalRank) const;

  /** \brief Returns the rank of the MPI process at a given coordinate
   *
   * This function is not required by any of the Tucker code
   * and is merely provided for the sake of debugging.
   */
  int getRank(const int* coords) const;

  //! Gets the number of MPI processors in a given dimension
  int getNumProcs(int d, bool squeezed) const;

  void squeeze(const Tucker::SizeArray& sz, const MPI_Comm& comm);

  const Tucker::SizeArray& size() const
  {
    return size_;
  }

private:
  /// @cond EXCLUDE
  ProcessorGrid(const ProcessorGrid& g);
  /// @endcond

  Tucker::SizeArray size_;

  //! MPI communicator storing the Cartesian grid information
  MPI_Comm cartComm_;

  //! Array of row communicators
  MPI_Comm* rowcomms_;

  //! Array of column communicators
  MPI_Comm* colcomms_;

  bool squeezed_;

  //! MPI communicator storing the Cartesian grid information
  MPI_Comm cartComm_squeezed_;

  //! Array of row communicators
  MPI_Comm* rowcomms_squeezed_;

  //! Array of column communicators
  MPI_Comm* colcomms_squeezed_;
};

} /* namespace TuckerMPI */
#endif /* PROCESSORGRID_HPP_ */
