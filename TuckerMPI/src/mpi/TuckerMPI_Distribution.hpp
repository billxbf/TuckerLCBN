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
 * \brief Describes the distribution of a tensor over a set of MPI processes
 *
 * @author: Alicia Klinvex
 */

#ifndef DISTRIBUTION_HPP_
#define DISTRIBUTION_HPP_

#include <memory>
#include "TuckerMPI_Map.hpp"
#include "TuckerMPI_ProcessorGrid.hpp"

namespace TuckerMPI {

//! Describes the distribution of a Cartesian grid over a set of MPI processes
class Distribution {
public:
  /** \brief Constructor
   *
   * The constructor will create a set of N maps, where N is the number
   * of dimensions, each of which describes how the data will be mapped
   * to the supplied processor grid.
   *
   * \param dims an array whose i-th entry is the number of points in the
   * Cartesian grid in dimension i.
   * \param procs an array whose i-th entry is the number of MPI processes
   * in the i-th dimension.
   */
  Distribution(const Tucker::SizeArray& dims,
      const Tucker::SizeArray& procs);

  //! Destructor
  ~Distribution();

  //! Returns the dimensions of the locally owned portion of the N-dimensional grid
  const Tucker::SizeArray& getLocalDims() const;

  //! Returns the dimensions of the N-dimensional grid
  const Tucker::SizeArray& getGlobalDims() const;

  //! Returns the processor grid
  const ProcessorGrid* getProcessorGrid() const;

  //! Returns the map of a given dimension
  const Map* getMap(int dimension, bool squeezed) const;

  const MPI_Comm& getComm(bool squeezed) const;

  bool ownNothing() const;

private:
  /// @cond EXCLUDE
  Distribution(const Distribution& d);
  /// @endcond

  void createMaps();

  /** \brief Finds and eliminates processes that don't have work
   *
   * Sets nonEmptyComm_.
   * Returns true if it found any such processes.
   */
  void findAndEliminateEmptyProcs(MPI_Comm& comm);

  //! Creates new processor grid without the processes that don't have work
  void updateProcessorGrid(const MPI_Comm& comm);

  /** \brief Size of the local grid; number of entries owned in each
   * dimension.
   *
   * This is different from what the map stores.  With the map, the
   * dimension is constant; with this array, the MPI process rank is constant.
   */
  Tucker::SizeArray localDims_;

  //! The global Cartesian grid size
  Tucker::SizeArray globalDims_;

  //! Maps MPI processes to a grid
  ProcessorGrid* grid_;

  /** \brief The maps describing the parallel distribution in each
   * dimension
   *
   * It is an array of Map*s rather than an array of Maps because Map
   * has no default constructor
   */
  Map** maps_;

  Map** maps_squeezed_;

  bool ownNothing_;

  bool squeezed_;
};

} /* namespace TuckerMPI */
#endif /* DISTRIBUTION_HPP_ */
