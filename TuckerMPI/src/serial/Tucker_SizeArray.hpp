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
 *
 * \brief A simple class for holding N-way size data.
 *
 * We use int rather than size_t for ease of compatibility with MPI
 * Likewise, we use an int* rather than an STL vector class for
 * guaranteed compatibility with MPI.
 *
 * @author Alicia Klinvex
 */

#ifndef SIZEARRAY_HPP_
#define SIZEARRAY_HPP_

#include<cassert>
#include<string.h>
#include<iostream>
#include "Tucker_Util.hpp"

namespace Tucker {

/** \brief A simple class for holding N-way size data.
 *
 * We use int rather than size_t for ease of compatibility with MPI.
 */
class SizeArray
{
public:
  /** \brief Constructor
   *
   * \param[in] n The number of entries in the array
   *
   * \exception std::invalid_argument n was not positive
   */
  SizeArray(const int n) :
    nsz_(n)
  {
    if(n <= 0)
      throw std::invalid_argument("n must be positive");
    sz_ = MemoryManager::safe_new_array<int>(nsz_);
  }

  //! Destructor
  ~SizeArray() {
    MemoryManager::safe_delete_array<int>(sz_,nsz_);
  }

  /** \brief Length/Number of entries in the array
   *
   * Named size to be consistent with C++ vector
   */
  int size() const {
    return nsz_;
  }

  /// Data pointer
  int* data() {
    return sz_;
  }

  /// Const data pointer
  const int* data() const {
    return sz_;
  }

  /** \brief Element access
   *
   * \param[in] i Index of entry to be returned
   *
   * \exception std::out_of_range \a i is not in the range [0, nsz_)
   */
  int& operator[](const int i) {
    if(i < 0 || i >= nsz_)
      throw std::out_of_range("invalid index");
    return sz_[i];
  }

  /** \brief Const element access
   *
   * \param[in] i Index of entry to be returned
   *
   * \exception std::out_of_range \a i is not in the range [0, nsz_)
   */
  const int& operator[](const int i) const {
    if(i < 0 || i >= nsz_)
      throw std::out_of_range("invalid index");
    return sz_[i];
  }

  //! Returns the product of all entries in the SizeArray
  size_t prod() const {
    return prod(0, nsz_-1);
  }
    
  /** \brief Inclusive product
   *
   * Returns the product of entries SizeArray[\a low] to SizeArray[\a high].
   * If the range is invalid (meaning \a low > \a high, \a low < 0,
   * or \a high >= #nsz_), \a defaultReturnVal is returned instead.
   *
   * \param[in] low Subscript of beginning of range
   * \param[in] high Subscript of end of range
   * \param[in] defaultReturnVal Value to be returned if input is invalid.
   * -1 by default.
   *
   * \note We chose to make defaultReturnVal a parameter because there are
   * some cases where you want it to be 1 by default and others where you would
   * want it to be 0 by default.
   */
  size_t prod(const int low, const int high, const int defaultReturnVal = -1) const {
    if(low < 0 || high >= nsz_) {
      std::cerr << "ERROR: prod(" << low << "," << high
          << ") is invalid because indices must be in the range [0,"
          << nsz_ << ").  Returning " << defaultReturnVal << std::endl;
        return defaultReturnVal;
    }
    if(low > high) {
      return defaultReturnVal;
    }
    size_t result = 1;
    for(int j = low; j <= high; j++)
      result *= sz_[j];
    return result;
  }

  /** \brief Returns whether two SizeArrays are exactly equal
   *
   * They are exactly equal if they have the same dimension, and every
   * entry is the same.
   *
   * \param[in] sz1 First SizeArray
   * \param[in] sz2 Second SizeArray
   */
  friend bool operator==(const SizeArray& sz1, const SizeArray& sz2)
  {
    if(sz1.size() != sz2.size())
      return false;

    for(int i=0; i<sz1.size(); i++) {
      if(sz1[i] != sz2[i] ){
        return false;
      }
    }
    return true;
  }

  /** \brief Returns whether two SizeArrays are not equal
   *
   * They are exactly equal if they have the same dimension, and every
   * entry is the same.  Otherwise, they are not.
   *
   * \param[in] sz1 First SizeArray
   * \param[in] sz2 Second SizeArray
   */
  friend bool operator!=(const SizeArray& sz1, const SizeArray& sz2)
  {
    return !(sz1 == sz2);
  }

  /** \brief Allows the user to print the SizeArray
   *
   * Usage: \n
   * <tt>
   * %SizeArray sz(4); \n
   * sz[0] = 3; sz[1] = 5; sz[2] = 7; sz[3] = 11; \n
   * std::cout << sz << std::endl;
   * </tt>
   *
   * produces the output \n
   * <tt>
   * 3 5 7 11
   * </tt>
   *
   * \param[in,out] os Output stream
   * \param[in] sz SizeArray to be printed
   */
  friend std::ostream& operator<<(std::ostream& os,
      const SizeArray& sz)
  {
    for(int i=0; i<sz.size(); i++) {
      os << sz[i] << " ";
    }
    return os;
  }

private:
  /// @cond EXCLUDE
  // Copy constructor
  // We define a private one with no implementation to disable it
  SizeArray(const SizeArray& sa);
  /// @endcond

  /** \brief Length of size array
   *
   * Size of the array cannot change after construction.
   */
  const int nsz_;

  /** \brief Size array entries
   *
   * Since the size of the array is not allowed to change after
   * construction, this can be allocated in the constructor and
   * destroyed in the destructor.
   */
  int* sz_;
};

} // end of namespace

#endif /* SIZEARRAY_HPP_ */
