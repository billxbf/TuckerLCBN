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
 * \brief Defines a timer class
 *
 * @author Alicia Klinvex
 */

#ifndef TUCKER_TIMER_HPP_
#define TUCKER_TIMER_HPP_

#include<chrono>

namespace Tucker {

/// Represents an interval of time
class Timer {
public:
  /** \brief Constructor
   *
   * Initializes the elapsed time to 0
   */
  Timer();

  /** \brief Destructor
   *
   * Does nothing because this class does no dynamic memory allocation.
   */
  ~Timer();

  /// Starts the timer
  void start();

  /** \brief Stops the timer and computes the elapsed time
   *
   * Uses std::chrono::steady_clock
   */
  void stop();

  /// Elapsed time
  double duration() const;

private:
  /// @cond EXCLUDE
  Timer(const Timer&);
  /// @endcond

  /// \brief Whether start has been called
  bool start_called_;

  /// \brief Time when start was called
  std::chrono::steady_clock::time_point start_;

  /// \brief Time when stop was called
  std::chrono::steady_clock::time_point end_;

  /// \brief Number of elapsed seconds
  double duration_;
};

} /* namespace Tucker */
#endif /* TUCKER_TIMER_HPP_ */
