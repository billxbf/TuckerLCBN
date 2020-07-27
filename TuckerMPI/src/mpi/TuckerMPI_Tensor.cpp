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

#include <limits>
#include "TuckerMPI_Tensor.hpp"

namespace TuckerMPI {

Tensor::Tensor(const Distribution* dist) :
    dist_(dist)
{
  localTensor_ = Tucker::MemoryManager::safe_new<Tucker::Tensor>(dist->getLocalDims());
}


Tensor::~Tensor()
{
  Tucker::MemoryManager::safe_delete<const Distribution>(dist_);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(localTensor_);
}


Tucker::Tensor* Tensor::getLocalTensor()
{
  return localTensor_;
}


const Tucker::Tensor* Tensor::getLocalTensor() const
{
  return localTensor_;
}


const Tucker::SizeArray& Tensor::getGlobalSize() const
{
  return dist_->getGlobalDims();
}


const Tucker::SizeArray& Tensor::getLocalSize() const
{
  return dist_->getLocalDims();
}


int Tensor::getGlobalSize(int n) const
{
  return dist_->getGlobalDims()[n];
}


int Tensor::getLocalSize(int n) const
{
  return dist_->getLocalDims()[n];
}


int Tensor::getNumDimensions() const
{
  return dist_->getGlobalDims().size();
}


// This function is necessary for TTM
// TTM needs to know how many entries EACH processor has,
// not just THIS processor.
// We don't technically need to return the entire distribution;
// that part is flexible.
const Distribution* Tensor::getDistribution() const
{
  return dist_;
}


size_t Tensor::getLocalNumEntries() const
{
  return localTensor_->getNumElements();
}


size_t Tensor::getGlobalNumEntries() const
{
  return dist_->getGlobalDims().prod();
}


// Compute the norm squared
double Tensor::norm2() const
{
  // Compute the local portion
  double localNorm2 = localTensor_->norm2();

  // Perform a reduction
  double globalNorm2;
  MPI_Allreduce(&localNorm2, &globalNorm2, 1,
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return globalNorm2;
}


void Tensor::print() const
{
  localTensor_->print();
}

void Tensor::rand()
{
  localTensor_->rand();
}

// \todo This function is never tested
Tensor* Tensor::subtract(const Tensor* t) const
{
  Tensor* sub = Tucker::MemoryManager::safe_new<Tensor>(dist_);

  size_t nnz = getLocalNumEntries();
  if(nnz > 0) {
    double* subdata = sub->localTensor_->data();
    double* thisdata = localTensor_->data();
    double* tdata = t->localTensor_->data();

    for(size_t i=0; i<nnz; i++) {
      subdata[i] = thisdata[i] - tdata[i];
    }
  }
  return sub;
}

// \todo This function is never tested
double Tensor::maxEntry() const
{
  double localMax = std::numeric_limits<double>::lowest();
  size_t nnz = getLocalNumEntries();
  if(nnz > 0) {
    double* data = localTensor_->data();
    for(size_t i=0; i<nnz; i++) {
      localMax = std::max(localMax,data[i]);
    }
  }

  double globalMax;
  MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE,
            MPI_MAX, MPI_COMM_WORLD);
  return globalMax;
}

// \todo This function is never tested
double Tensor::minEntry() const
{
  double localMin = std::numeric_limits<double>::max();
  size_t nnz = getLocalNumEntries();
  if(nnz > 0) {
    double* data = localTensor_->data();
    for(size_t i=0; i<nnz; i++) {
      localMin = std::min(localMin,data[i]);
    }
  }

  double globalMin;
  MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE,
            MPI_MIN, MPI_COMM_WORLD);
  return globalMin;
}

bool isApproxEqual(const Tensor* t1, const Tensor* t2,
    double tol)
{
  if(t1->getGlobalSize() != t2->getGlobalSize()) {
    std::cerr << "t1 and t2 have different global sizes\n";
    std::cerr << "t1: " << t1->getGlobalSize() << std::endl;
    std::cerr << "t2: " << t2->getGlobalSize() << std::endl;
    return false;
  }
  return (isApproxEqual(t1->getLocalTensor(),
      t2->getLocalTensor(), tol));
}

} /* namespace TuckerMPI */
