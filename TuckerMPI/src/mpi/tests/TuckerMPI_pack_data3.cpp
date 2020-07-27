/*
 * pack_data.cpp
 *
 * This tests the mode-2 packData function.
 * Since this example only has three modes, it should
 * do nothing.
 *
 *  Created on: Jul 12, 2016
 *      Author: amklinv
 */

#include <cstdlib>
#include "mpi.h"
#include "TuckerMPI.hpp"

bool checkArrayEqual(double* arr1, double* arr2, int numEl);

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Create a 3x4x3 tensor
  Tucker::SizeArray* sa =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(3);
  (*sa)[0] = 3; (*sa)[1] = 4; (*sa)[2] = 3;
  Tucker::Tensor* tensor =
      Tucker::MemoryManager::safe_new<Tucker::Tensor>(*sa);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(sa);

  // Fill it with the entries 0:35
  double* data = tensor->data();
  for(int i=0; i<36; i++)
    data[i] = i;

  // Create a map describing the distribution
  TuckerMPI::Map* map =
      Tucker::MemoryManager::safe_new<TuckerMPI::Map>(4,MPI_COMM_WORLD);

  // Pack the tensor
  packForTTM(tensor, 2, map);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Map>(map);

  tensor->print();

  // Check the result
  double trueResult[36];
  for(int i=0; i<36; i++) trueResult[i] = i;
  bool equal = checkArrayEqual(tensor->data(),trueResult,36);
  if(!equal) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Call MPI_Finalize
  MPI_Finalize();
  return EXIT_SUCCESS;
}

bool checkArrayEqual(double* arr1, double* arr2, int numEl)
{
  for(int i=0; i<numEl; i++) {
    if(arr1[i] != arr2[i]) {
      return false;
    }
  }
  return true;
}
