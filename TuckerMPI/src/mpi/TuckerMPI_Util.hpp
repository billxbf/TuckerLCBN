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
 * \brief Stores utility routines
 *
 * @author Alicia Klinvex
 */

#ifndef UTIL_MPI_HPP_
#define UTIL_MPI_HPP_

#include "Tucker_Timer.hpp"
#include "TuckerMPI_Tensor.hpp"
#include "TuckerMPI_Matrix.hpp"

namespace TuckerMPI {

/// @cond EXCLUDE
extern "C" void dcopy_(const int*, const double*, const int*,
    double*, const int*);

extern "C" void dsyrk_(const char*, const char*, const int*, const int*,
    const double*, const double*, const int*, const double*, double*,
    const int*);

extern "C" void dgemm_(const char*, const char*, const int*,
    const int*, const int*, const double*, const double*,
    const int*, const double*, const int*, const double*,
    double*, const int*);
/// @endcond

/** \brief Determines whether packing is necessary for computing the Gram matrix
 */
bool isPackForGramNecessary(int n, const Map* origMap, const Map* redistMap);

/** \brief Packs the data for computing the Gram matrix
 */
const double* packForGram(const Tensor* Y, int n, const Map* redistMat);

/** \brief Determines whether unpacking is necessary for computing the Gram matrix
 */
bool isUnpackForGramNecessary(int n, int ndims, const Map* origMap,
    const Map* redistMap);

/** \brief Unpacks the data for computing the Gram matrix
 */
void unpackForGram(int n, int ndims, Matrix* redistMat,
    const double* dataToUnpack, const Map* origMap);

/** \brief Redistributes the tensor for computing the Gram matrix using the
 * new algorithm
 */
const Matrix* redistributeTensorForGram(const Tensor* Y, int n,
    Tucker::Timer* pack_timer=0, Tucker::Timer* alltoall_timer=0,
    Tucker::Timer* unpack_timer=0);

/** \brief Local rank-k update for computing the Gram matrix
 */
const Tucker::Matrix* localRankKForGram(const Matrix* Y, int n, int ndims);

/** \brief Local matrix-matrix multiply for computing the Gram matrix
 */
void localGEMMForGram(const double* Y1, int nrowsY1, int n,
    const Tensor* Y2, double* result);

/** \brief Perform a reduction for the Gram matrix computation
 *
 */
Tucker::Matrix* reduceForGram(const Tucker::Matrix* U);

/** \brief Packs the tensor for the TTM
 *
 */
void packForTTM(Tucker::Tensor* Y, int n, const Map* map);

} // end namespace TuckerMPI

#endif /* UTIL_MPI_HPP_ */
