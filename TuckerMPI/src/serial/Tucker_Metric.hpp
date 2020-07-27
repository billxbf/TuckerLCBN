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
 * \brief Contains an enumerated metric type and a struct that stores
 * all possible metrics
 *
 * @author Alicia Klinvex
 */

#ifndef METRIC_HPP_
#define METRIC_HPP_

#include<exception>

namespace Tucker {

/** \brief Statistic to be computed
 *
 * All values are powers of 2, so they can be summed together
 * to produce a unique combination of metrics.  This type is
 * used as an input to Tucker::computeSliceMetrics.
 *
 * \warning Actual values may change at any time.  Please do not
 * assume that MIN will always be 1, etc.
 */
enum Metric {
  MIN = 1, ///< Minimum value
  MAX = 2, ///< Maximum value
  SUM = 4, ///< Sum of all values
  NORM1 = 8, ///< 1-norm
  NORM2 = 16, ///< 2-norm
  MEAN = 32, ///< Mean
  VARIANCE = 64 ///< Variance
};

/** \brief Class that stores all metrics
 *
 * Returned by Tucker::computeSliceMetrics.
 */
class MetricData
{
public:
  /** \brief Constructor
   *
   * Allocates memory for an array of metrics.
   * Only allocates memory for the requested metrics.
   * For instance, if \a metrics = MIN + MAX + NORM1,
   * minData, maxData, and normData will all be allocated.
   * All other metrics will be null pointers.
   *
   * \param[in] metrics A sum of Metric values
   * \param[in] dimension Size of array
   *
   * \todo NORM1 and NORM2 are never tested
   */
  MetricData(const int metrics, const int dimension) :
    dimension_(dimension)
  {
    assert(dimension > 0);

    if(metrics & MIN)
      minData_ = MemoryManager::safe_new_array<double>(dimension);
    else
      minData_ = 0;

    if(metrics & MAX)
      maxData_ = MemoryManager::safe_new_array<double>(dimension);
    else
      maxData_ = 0;

    if(metrics & SUM)
      sumData_ = MemoryManager::safe_new_array<double>(dimension);
    else
      sumData_ = 0;

    if(metrics & NORM1)
      norm1Data_ = MemoryManager::safe_new_array<double>(dimension);
    else
      norm1Data_ = 0;

    if(metrics & NORM2)
      norm2Data_ = MemoryManager::safe_new_array<double>(dimension);
    else
      norm2Data_ = 0;

    if(metrics & MEAN)
      meanData_ = MemoryManager::safe_new_array<double>(dimension);
    else
      meanData_ = 0;

    if(metrics & VARIANCE)
      varianceData_ = MemoryManager::safe_new_array<double>(dimension);
    else
      varianceData_ = 0;
  }

  /** \brief Destructor
   *
   * Deallocates all allocated arrays
   */
  ~MetricData()
  {
    if(minData_) MemoryManager::safe_delete_array<double>(minData_,dimension_);
    if(maxData_) MemoryManager::safe_delete_array<double>(maxData_,dimension_);
    if(sumData_) MemoryManager::safe_delete_array<double>(sumData_,dimension_);
    if(norm1Data_) MemoryManager::safe_delete_array<double>(norm1Data_,dimension_);
    if(norm2Data_) MemoryManager::safe_delete_array<double>(norm2Data_,dimension_);
    if(meanData_) MemoryManager::safe_delete_array<double>(meanData_,dimension_);
    if(varianceData_) MemoryManager::safe_delete_array<double>(varianceData_,dimension_);
  }

  /** \brief Returns the min-data array
   *
   * \exception std::runtime_error min-data was never allocated
   */
  double* getMinData()
  {
    if(minData_) return minData_;
    throw std::runtime_error("min data was never allocated!");
  }

  /** \brief Returns the max-data array
   *
   * \exception std::runtime_error max-data was never allocated
   */
  double* getMaxData()
  {
    if(maxData_) return maxData_;
    throw std::runtime_error("max data was never allocated!");
  }

  /** \brief Returns the sum-data array
   *
   * \exception std::runtime_error sum-data was never allocated
   */
  double* getSumData()
  {
    if(sumData_) return sumData_;
    throw std::runtime_error("sum data was never allocated!");
  }

  /** \brief Returns the norm1-data array
   *
   * \exception std::runtime_error norm1-data was never allocated
   */
  double* getNorm1Data()
  {
    if(norm1Data_) return norm1Data_;
    throw std::runtime_error("norm1 data was never allocated!");
  }

  /** \brief Returns the norm2-data array
   *
   * \exception std::runtime_error norm2-data was never allocated
   */
  double* getNorm2Data()
  {
    if(norm2Data_) return norm2Data_;
    throw std::runtime_error("norm2 data was never allocated!");
  }

  /** \brief Returns the mean-data array
   *
   * \exception std::runtime_error mean-data was never allocated
   */
  double* getMeanData()
  {
    if(meanData_) return meanData_;
    throw std::runtime_error("mean data was never allocated!");
  }

  /** \brief Returns the variance array
   *
   * \exception std::runtime_error variance-data was never allocated
   */
  double* getVarianceData()
  {
    if(varianceData_) return varianceData_;
    throw std::runtime_error("variance data was never allocated!");
  }

private:
  /// @cond EXCLUDE
  // Disable copy constructor
  MetricData(const MetricData& m);
  /// @endcond

  double* minData_; ///< Minimum value array
  double* maxData_; ///< Maximum value array
  double* sumData_; ///< Sum of all values array
  double* norm1Data_; ///< 1-norm array
  double* norm2Data_; ///< 2-norm array
  double* meanData_; ///< Mean array
  double* varianceData_; ///< Variance array
  int dimension_;
};

}

#endif /* METRIC_HPP_ */
