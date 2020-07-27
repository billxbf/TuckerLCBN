/*
 * TuckerMPI_ttm.hpp
 *
 *  Created on: Nov 2, 2016
 *      Author: amklinv
 */

#ifndef MPI_TUCKERMPI_TTM_HPP_
#define MPI_TUCKERMPI_TTM_HPP_

#include<cmath>
#include "Tucker.hpp"
#include "TuckerMPI_Util.hpp"

namespace TuckerMPI
{
/** \brief Parallel tensor times matrix computation
 *
 * \param X A parallel tensor
 * \param n The dimension for the tensor unfolding
 * \param U A sequential matrix
 * \param Utransp Whether to compute X * U^T or X * U
 * \param mult_timer Timer for the local multiplication
 * \param pack_timer Timer for packing the data
 * \param reduce_scatter_timer Timer for the reduce-scatter
 * \param reduce_timer Timer for the reduction
 */
Tensor* ttm(const Tensor* X, const int n,
    const Tucker::Matrix* const U, bool Utransp=false,
    Tucker::Timer* mult_timer=0, Tucker::Timer* pack_timer=0,
    Tucker::Timer* reduce_scatter_timer=0,
    Tucker::Timer* reduce_timer=0, size_t nnz_limit=0);

} // end namespace TuckerMPI

#endif /* MPI_TUCKERMPI_TTM_HPP_ */
