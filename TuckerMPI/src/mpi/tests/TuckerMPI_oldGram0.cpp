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
  (*sz)[0] = 4; (*sz)[1] = 4; (*sz)[2] = 4;
  Tucker::SizeArray* nprocsPerDim =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(ndims);
  (*nprocsPerDim)[0] = 2; (*nprocsPerDim)[1] = 2; (*nprocsPerDim)[2] = 2;
  TuckerMPI::Distribution* dist =
        Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*sz,*nprocsPerDim);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(sz);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(nprocsPerDim);

  // Create a tensor
  TuckerMPI::Tensor* tensor =
      Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(dist);

  // Read the entries from a file
  std::string filename = "input_files/tensor64.mpi";
  TuckerMPI::importTensorBinary(filename.c_str(),tensor);

  // Compute the gram matrix in dimension 0
  const Tucker::Matrix* mat = TuckerMPI::oldGram(tensor,0);

  double trueData0[36] = {19840, 20320, 20800, 21280,
                         20320, 20816, 21312, 21808,
                         20800, 21312, 21824, 22336,
                         21280, 21808, 22336, 22864};

  bool equal = checkUTEqual(mat->data(), trueData0, 4);
  if(!equal) {
    std::cerr << "Mode 0 failure\n";
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<const Tucker::Matrix>(mat);

  // Compute the gram matrix in dimension 1
  mat = TuckerMPI::oldGram(tensor,1);

  double trueData1[36] = {15544, 17176, 18808, 20440,
                          17176, 19064, 20952, 22840,
                          18808, 20952, 23096, 25240,
                          20440, 22840, 25240, 27640};

  equal = checkUTEqual(mat->data(), trueData1, 4);
  if(!equal) {
    std::cerr << "Mode 1 failure\n";
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<const Tucker::Matrix>(mat);

  // Compute the gram matrix in dimension 2
  mat = TuckerMPI::oldGram(tensor,2);

  double trueData2[36] = {1240, 3160, 5080, 7000,
                          3160, 9176, 15192, 21208,
                          5080, 15192, 25304, 35416,
                          7000, 21208, 35416, 49624};

  equal = checkUTEqual(mat->data(), trueData2, 4);
  if(!equal) {
    std::cerr << "Mode 2 failure\n";
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
//        std::cerr << "arr[" << ind << "] (" << r << "," << c
//            << ") are not equal(" << arr1[ind] << " "
//            << arr2[ind] << ")\n";
        return false;
      }
    }
  }
  return true;
}
