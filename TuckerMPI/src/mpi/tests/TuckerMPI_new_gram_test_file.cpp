/*
 * norm_test.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: amklinv
 */

#include <cmath>
#include "TuckerMPI.hpp"

bool checkUTEqual(const double* arr1, const double* arr2, int numRows);
bool runSim(Tucker::SizeArray& procs);

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  bool success;

  // Create a processor grid with the number of MPI processes in each dimension
  Tucker::SizeArray* procs =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(4);
  (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
  if(nprocs == 2) {
    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 3) {
    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 4) {
    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 5) {
    (*procs)[0] = 5; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 5; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 5; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 6) {
    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 7) {
    (*procs)[0] = 7; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 7; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 7; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 7;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 8) {
    (*procs)[0] = 8; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 8; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 8; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 9) {
    (*procs)[0] = 9; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 9; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 9; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 9;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 10) {
    (*procs)[0] = 10; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 10; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 10; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 10;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 5; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 5; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 5; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 5; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 5; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 5; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 11) {
    (*procs)[0] = 11; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 11; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 11; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 11;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 12) {
    (*procs)[0] = 12; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 12; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 12; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 12;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 13) {
    (*procs)[0] = 13; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 13; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 13; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 13;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 14) {
    (*procs)[0] = 14; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 14; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 14; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 14;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 7; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 7; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 7; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 7; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 7; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 7; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 7; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 7; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 7;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 7; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 7;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 7;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 15) {
    (*procs)[0] = 15; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 15; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 15; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 15;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 5; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 5; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 5; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 5; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 5; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 5; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 5; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 5;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 16) {
    (*procs)[0] = 16; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 16; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 16; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 16;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 8; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 8; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 8; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 8; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 8; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 8; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 4; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 4; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else if(nprocs == 24) {
    (*procs)[0] = 24; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 24; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 24; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 24;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 12; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 12; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 12; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 12; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 12; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 12; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 8; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 8; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 8; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 8; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 6; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 6; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 4; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 4; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 8; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 8; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 8; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 8;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 4; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 4; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 3; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 12; (*procs)[2] = 1; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 12; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 1; (*procs)[3] = 12;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 12; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 12;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 12;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 6; (*procs)[2] = 2; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 6; (*procs)[2] = 1; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 6; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 6; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 4; (*procs)[2] = 3; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 4; (*procs)[2] = 1; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 4; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 4; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 4; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 1; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 3; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 4;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 3; (*procs)[2] = 2; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 6; (*procs)[3] = 1;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 1; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 1; (*procs)[2] = 2; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 1; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 6;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 3; (*procs)[3] = 2;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;

    (*procs)[0] = 2; (*procs)[1] = 2; (*procs)[2] = 2; (*procs)[3] = 3;
    success = runSim(*procs); if(!success) return EXIT_FAILURE;
  }
  else {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(procs);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}


bool checkUTEqual(const double* arr1, const double* arr2, int numRows)
{
  const double TOL = 1e-10;

  for(int r=0; r<numRows; r++) {
    for(int c=r; c<numRows; c++) {
      int ind = r+c*numRows;
      if(std::abs(arr1[ind]-arr2[ind]) > TOL) {
        std::cerr << "ERROR: The true solution is " << arr2[ind]
                  << ", but the computed solution was " << arr1[ind]
                  << ", a difference of " << std::abs(arr1[ind]-arr2[ind])
                  << std::endl;
        return false;
      }
    }
  }
  return true;
}

bool runSim(Tucker::SizeArray& procs)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  Tucker::Matrix* matrix;
  bool matchesTrueSol;

  if(rank == 0) std::cout << procs << std::endl;

  // Create a SizeArray with the dimensions of the tensor
  Tucker::SizeArray dims(4);
  dims[0] = 3; dims[1] = 5; dims[2] = 7; dims[3] = 11;

  // Create a distribution object
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);

  // Create a tensor with that distribution
  TuckerMPI::Tensor tensor(dist);

  // Read the tensor from a binary file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11.mpi",&tensor);

  // Compute the Gram matrix
  matrix = TuckerMPI::newGram(&tensor,0);

  const double TRUE_SOLUTION_0[3*3] =
    {30.443306716415002, 2.839326384970902, 3.302846757455287,
     2.839326384970902, 30.254535550071385, -1.525701935166856,
     3.302846757455286,  -1.525701935166856,  33.203090378426815};

  matchesTrueSol = checkUTEqual(matrix->data(), TRUE_SOLUTION_0, 3);
  if(!matchesTrueSol)
    return false;

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(matrix);

  // Compute the Gram matrix
  matrix = TuckerMPI::newGram(&tensor,1);

  const double TRUE_SOLUTION_1[5*5] =
    {19.818638147579669, -1.262410163567289, -0.702854185504209, -1.530487842705398, 0.076595773936855,
     -1.262410163567289, 19.055467182423385, -0.066700982927269,  2.017798574327055, -2.178581571284073,
     -0.702854185504209, -0.066700982927269, 18.565387462362573, -1.776252595103098, 1.511348003921888,
     -1.530487842705398,  2.017798574327055, -1.776252595103098, 19.948082945419085, -1.701919532297699,
      0.076595773936855, -2.178581571284071,  1.511348003921888, -1.701919532297699, 16.513356907128461};

  matchesTrueSol = checkUTEqual(matrix->data(), TRUE_SOLUTION_1, 5);

  if(!matchesTrueSol)
    return false;

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(matrix);

  // Compute the Gram matrix
  matrix = TuckerMPI::newGram(&tensor,2);

  const double TRUE_SOLUTION_2[7*7] =
    {15.008467435803363, 0.759679780649146, 0.036286707752565, -0.289429119623268, 0.297697996043826, 1.775651754316252, 0.290583279984489,
     0.759679780649146, 13.534258753086892, 0.216613144804639, -0.909444710671252, -0.440040933466265, -1.316163513868651, 0.331266618564718,
     0.036286707752565, 0.216613144804639, 13.107978854355352, 0.701142578929455, 1.325848667245948, -0.355644067560675, 0.307781943716007,
     -0.289429119623268, -0.909444710671252, 0.701142578929455, 13.223174991133867, 1.085334988154543, 1.098880462506622, 2.998592432064371,
     0.297697996043825, -0.440040933466265, 1.325848667245947, 1.085334988154542, 13.872205816335024, -0.682375996350724, 1.263779626695686,
     1.775651754316251, -1.316163513868651, -0.355644067560674, 1.098880462506623, -0.682375996350724, 12.100552292592960, 0.747070331494893,
     0.290583279984489, 0.331266618564718, 0.307781943716007, 2.998592432064369, 1.263779626695687, 0.747070331494893, 13.054294501605717};

  matchesTrueSol = checkUTEqual(matrix->data(), TRUE_SOLUTION_2, 7);

  if(!matchesTrueSol)
    return false;

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(matrix);

  // Compute the Gram matrix
  matrix = TuckerMPI::newGram(&tensor,3);

  const double TRUE_SOLUTION_3[11*11] =
    {9.161015785664386, 0.727680925086397, -0.027808461977047, -1.858926623169289, 2.161910031054029, 0.538498274148853, -1.513458574761186, 1.560326634563823, 0.092898181083783, 0.193490784986641, -1.989938728094985,
        0.727680925086397, 8.579166359218990, 0.141184182846715, 0.895099070075687, -0.217980995935905, -0.257633332960250, 0.874112582420890, -1.046151056779062, 0.285393081978074, -0.864025655243411, 0.241418995156843,
        -0.027808461977047, 0.141184182846715, 8.487659393716891, -0.196516290753068, 0.855490962332431, 0.483991332254019, -0.096130754960697, 1.054249161028251, 0.718345150564935, -0.962619165522184, -1.619740748464319,
        -1.858926623169289, 0.895099070075687, -0.196516290753068, 8.222091258079500, -0.657900464867635, -0.649610924287497, 1.661744576911002, -0.193946562374636, 0.072191263966888, 0.586444766094151, -0.340892667710611,
        2.161910031054029, -0.217980995935905, 0.855490962332431, -0.657900464867635, 6.873280265688026, 1.816757520854934, -1.262051827038299, 0.678477183630724, 0.135986851830552, -0.228764676015497, -1.666495122122524,
        0.538498274148853, -0.257633332960250, 0.483991332254019, -0.649610924287497, 1.816757520854934, 8.227968322612568, 0.861291545524726, -0.353270918783090, -0.570666332802020, -1.789631071892280, 0.442503856536748,
        -1.513458574761186, 0.874112582420890, -0.096130754960697, 1.661744576911002, -1.262051827038299, 0.861291545524726, 9.269694148981491, -1.574380517040051, 1.099770875691088, 0.069588373131300, -0.717801546637968,
        1.560326634563823, -1.046151056779062, 1.054249161028251, -0.193946562374636, 0.678477183630724, -0.353270918783090, -1.574380517040051, 9.070501064224349, -0.973741303421328, -0.752651299409004, 0.151673397796920,
        0.092898181083782, 0.285393081978074, 0.718345150564935, 0.072191263966888, 0.135986851830551, -0.570666332802020, 1.099770875691088, -0.973741303421328, 7.925894762690504, 0.447258387553726, -0.591921295268904,
        0.193490784986642, -0.864025655243411, -0.962619165522184, 0.586444766094151, -0.228764676015497, -1.789631071892280, 0.069588373131299, -0.752651299409004, 0.447258387553726, 8.517370143817459, -0.523119424574941,
        -1.989938728094985, 0.241418995156843, -1.619740748464319, -0.340892667710611, -1.666495122122525, 0.442503856536748, -0.717801546637968, 0.151673397796920, -0.591921295268904, -0.523119424574941, 9.566291140219018};

  matchesTrueSol = checkUTEqual(matrix->data(), TRUE_SOLUTION_3, 11);

  if(!matchesTrueSol)
    return false;

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(matrix);

  return true;
}
