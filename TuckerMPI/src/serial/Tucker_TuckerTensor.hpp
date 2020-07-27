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
 * \brief Stores a %Tucker decomposition
 *
 * @author Alicia Klinvex
 */

#ifndef TUCKERTENSOR_HPP_
#define TUCKERTENSOR_HPP_

#include "Tucker_Tensor.hpp"
#include "Tucker_Matrix.hpp"
#include "Tucker_Timer.hpp"

namespace Tucker {

/** \brief A structure for storing a %Tucker decomposition
 *
 * Allows users to access the data directly
 */
class TuckerTensor {
public:
  /** \brief Constructor
   *
   * Allocates memory for the array of factors, the array of eigenvalues,
   * and all timers.
   * \param[in] ndims Number of dimensions
   *
   * \note This function does not allocate memory for the core tensor (since
   * the size may not be known in advance), nor does it allocate memory
   * for the individual factors (for the same reason).  It also does not
   * allocate memory for each individual array of eigenvalues.
   *
   * \exception std::runtime_error \a ndims <= 0
   */
  TuckerTensor(const int ndims)
  {
    if(ndims <= 0) {
      std::ostringstream oss;
      oss << "Tucker::TuckerTensor::TuckerTensor(const int ndims): ndims = "
          << ndims << " <= 0";
      throw std::runtime_error(oss.str());
    }

    N = ndims;
    U = MemoryManager::safe_new_array<Matrix*>(N);
    eigenvalues = MemoryManager::safe_new_array<double*>(N);
    G = 0;

    for(int i=0; i<N; i++) {
      eigenvalues[i] = 0;
      U[i] = 0;
    }

    gram_timer_ = MemoryManager::safe_new_array<Timer>(ndims);
    eigen_timer_ = MemoryManager::safe_new_array<Timer>(ndims);
    ttm_timer_ = MemoryManager::safe_new_array<Timer>(ndims);
  }

  /** \brief Destructor
   *
   * Deallocates memory for the core tensor, the factor matrices,
   * the eigenvalues, and the timers.
   */
  ~TuckerTensor()
  {
    if(G) MemoryManager::safe_delete<Tensor>(G);
    for(int i=0; i<N; i++) {
      if(eigenvalues[i]) MemoryManager::safe_delete_array<double>(eigenvalues[i],U[i]->nrows());
      if(U[i]) MemoryManager::safe_delete<Matrix>(U[i]);
    }
    MemoryManager::safe_delete_array<Matrix*>(U,N);
    MemoryManager::safe_delete_array<double*>(eigenvalues,N);

    MemoryManager::safe_delete_array<Timer>(gram_timer_,N);
    MemoryManager::safe_delete_array<Timer>(eigen_timer_,N);
    MemoryManager::safe_delete_array<Timer>(ttm_timer_,N);
  }

  /** \brief Prints some runtime information
   *
   */
  void printTimers() const
  {
    for(int i=0; i<N; i++) {
      std::cout << "Gram(" << i << ")      : " << std::scientific
          << gram_timer_[i].duration() << std::endl;

      std::cout << "Eigensolve(" << i << "): " << std::scientific
          << eigen_timer_[i].duration() << std::endl;

      std::cout << "TTM(" << i << ")       : " << std::scientific
          << ttm_timer_[i].duration() << std::endl;
    }

    std::cout << "Total        : " << std::scientific
        << total_timer_.duration() << std::endl;
  }

  //! The core tensor
  Tensor* G;

  /** \brief Factors
   *
   * #U[n] is a pointer to the n-th factor
   *
   * \note These are stored as pointers because Matrix has no
   * default constructor.
   */
  Matrix** U;

  //! The number of factors
  int N;

  /** \brief Eigenvalues
   *
   * #eigenvalues[n] is an array which holds the eigenvalues of
   * the n-th Gram matrix
   */
  double** eigenvalues;

  /** \note STHOSVD has been declared as a friend function of
   * TuckerTensor so that the timers can remain private
   */
  friend const struct TuckerTensor* STHOSVD(const Tensor* const X,
      const double epsilon, bool flipSign);

  /** \note STHOSVD has been declared as a friend function of
   * TuckerTensor so that the timers can remain private
   */
  friend const struct TuckerTensor* STHOSVD(const Tensor* const X,
      const SizeArray* const reducedI, bool flipSign);

private:
  /// @cond EXCLUDE
  TuckerTensor(const TuckerTensor& tt);
  /// @endcond

  /// \brief Array of timers for Gram matrix computation
  Timer* gram_timer_;

  /// \brief Array of timers for eigensolver computation
  Timer* eigen_timer_;

  /// \brief Array of timers for TTM computation
  Timer* ttm_timer_;

  /// \brief Total ST-HOSVD runtime
  Timer total_timer_;
};

} // end namespace Tucker

#endif /* TUCKERTENSOR_HPP_ */
