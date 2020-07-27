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
 * \brief Contains a class for storing a distributed tensor
 *
 * @author Alicia Klinvex
 */

#ifndef TENSOR_MPI_HPP_
#define TENSOR_MPI_HPP_

#include "TuckerMPI_Distribution.hpp"
#include "Tucker_Tensor.hpp"

namespace TuckerMPI {

//! Distributed tensor
class Tensor
{
public:
  /** \brief Constructor
   *
   * \param dist Distribution object containing the processor grid
   * and size of the tensor
   *
   * \note \a dist is owned by the tensor and will be freed by it
   */
  Tensor(const Distribution* dist);

  //! Destructor
  ~Tensor();

  //! Returns the sequential locally owned portion of the tensor
  Tucker::Tensor* getLocalTensor();

  //! Returns the sequential locally owned portion of the tensor
  const Tucker::Tensor* getLocalTensor() const;

  //! Returns the global size
  const Tucker::SizeArray& getGlobalSize() const;

  //! Returns the size of the locally owned tensor
  const Tucker::SizeArray& getLocalSize() const;

  //! Returns the global size of dimension n
  int getGlobalSize(int n) const;

  //! Returns the local size of dimension n
  int getLocalSize(int n) const;

  //! Returns the number of dimensions
  int getNumDimensions() const;

  //! Returns the distribution object
  const Distribution* getDistribution() const;

  //! Returns the local number of entries
  size_t getLocalNumEntries() const;

  size_t getGlobalNumEntries() const;

  //! Compute the norm squared
  double norm2() const;

  //! Print the entries of the tensor as a one-dimensional array
  void print() const;

  void rand();

  //! Returns \a this - \a t
  Tensor* subtract(const Tensor* t) const;

  //! Returns the maximum entry of the tensor
  double maxEntry() const;

  //! Returns the minimum entry of the tensor
  double minEntry() const;

private:
  /// @cond EXCLUDE
  Tensor(const Tensor& t);
  /// @endcond

  //! Sequential Tensor of locally owned data
  Tucker::Tensor* localTensor_;
  //! Distribution object
  const Distribution* dist_;
};

bool isApproxEqual(const Tensor* t1, const Tensor* t2,
    double tol);

} /* namespace TuckerMPI */
#endif /* TENSOR_MPI_HPP_ */
