/*
 * matrix_read_test.cpp
 *
 *  Created on: Jul 12, 2016
 *      Author: amklinv
 */

#include<cstdlib>
#include "TuckerMPI.hpp"

bool checkUTEqual(const double* arr1, const double* arr2, int numRows);

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Get rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // Create a distribution object
  int ndims = 3;
  Tucker::SizeArray* sz =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(ndims);
  (*sz)[0] = 2; (*sz)[1] = 6; (*sz)[2] = 3;
  Tucker::SizeArray* nprocsPerDim =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(ndims);
  (*nprocsPerDim)[0] = 1; (*nprocsPerDim)[1] = 3; (*nprocsPerDim)[2] = 1;
  TuckerMPI::Distribution* dist =
        Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*sz,*nprocsPerDim);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(sz);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(nprocsPerDim);

  // Create a tensor
  TuckerMPI::Tensor* tensor =
      Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(dist);

  // Read the entries from a file
  std::string filename = "input_files/tensor36.mpi";
  TuckerMPI::importTensorBinary(filename.c_str(),tensor);

  // Compute the gram matrix in dimension 1
  const Tucker::Matrix* mat = TuckerMPI::newGram(tensor,1);

  double trueData[36] = {1515, 1665, 1815, 1965, 2115, 2265,
                         1665, 1839, 2013, 2187, 2361, 2535,
                         1815, 2013, 2211, 2409, 2607, 2805,
                         1965, 2187, 2409, 2631, 2853, 3075,
                         2115, 2361, 2607, 2853, 3099, 3345,
                         2265, 2535, 2805, 3075, 3345, 3615};

  bool equal = checkUTEqual(mat->data(), trueData, 6);
  if(!equal) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<const Tucker::Matrix>(mat);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

bool checkUTEqual(const double* arr1, const double* arr2, int numRows)
{
  for(int r=0; r<numRows; r++) {
    for(int c=r; c<numRows; c++) {
      int ind = r+c*numRows;
      if(arr1[ind] != arr2[ind]) {
        return false;
      }
    }
  }
  return true;
}
