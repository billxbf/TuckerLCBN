/*
 * TuckerMPI_converter.cpp
 *
 *  Created on: Sep 1, 2016
 *      Author: amklinv
 */

#include "TuckerMPI.hpp"
#include "Tucker.hpp"

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  Tucker::Tensor* t = Tucker::importTensor(argv[1]);
  TuckerMPI::exportTensorBinary(argv[2],t);

  MPI_Finalize();
}


