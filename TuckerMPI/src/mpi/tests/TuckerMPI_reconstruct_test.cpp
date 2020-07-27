/*
 * reconstruct_test.cpp
 *
 *  Created on: Jul 25, 2016
 *      Author: amklinv
 */

#include<cstdlib>
#include "TuckerMPI.hpp"

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Create the SizeArray
  Tucker::SizeArray* size =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(3);
  (*size)[0] = 2;
  (*size)[1] = 3;
  (*size)[2] = 5;

  // Create the MPI processor grid
  int ndims = 3;
  Tucker::SizeArray* nprocsPerDim =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(ndims);
  (*nprocsPerDim)[0] = 1; (*nprocsPerDim)[1] = 2; (*nprocsPerDim)[2] = 3;

  // Create the distribution object
  TuckerMPI::Distribution* dist =
        Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*size,*nprocsPerDim);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(size);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(nprocsPerDim);

  // Create the tensor
  TuckerMPI::Tensor* tensor =
      Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(dist);

  // Read the entries from a file
  std::string filename = "input_files/tensor24.mpi";
  TuckerMPI::importTensorBinary(filename.c_str(),tensor);

  const struct TuckerMPI::TuckerTensor* factorization =
      TuckerMPI::STHOSVD(tensor,1e-6);

  // Reconstruct the original tensor
  TuckerMPI::Tensor* temp = TuckerMPI::ttm(factorization->G,0,
      factorization->U[0]);
  TuckerMPI::Tensor* temp2 = TuckerMPI::ttm(temp,1,
        factorization->U[1]);
  TuckerMPI::Tensor* temp3 = TuckerMPI::ttm(temp2,2,
        factorization->U[2]);

  bool eq = isApproxEqual(temp3,tensor,1e-10);
  if(!eq) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // Free some memory
  Tucker::MemoryManager::safe_delete<const struct TuckerMPI::TuckerTensor>(factorization);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(temp);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(temp2);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(temp3);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(tensor);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}



