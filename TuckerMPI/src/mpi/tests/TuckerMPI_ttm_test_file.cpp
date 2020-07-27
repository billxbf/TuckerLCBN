/*
 * ttm_test_file.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: amklinv
 */

#include "Tucker.hpp"
#include "TuckerMPI.hpp"

bool runSim(Tucker::SizeArray& procs);

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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


bool runSim(Tucker::SizeArray& procs)
{
  Tucker::Matrix* mat;
  TuckerMPI::Tensor* mySol;
  TuckerMPI::Distribution* trueDist;
  TuckerMPI::Tensor* trueSol;
  bool approxEqual;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank == 0)
    std::cout << procs << std::endl;

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

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/3x2.txt");

  // Read the true solution
  dims[0] = 2; dims[1] = 5; dims[2] = 7; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/3x2_mult_transp.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 0, mat, true);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/4x3.txt");

  // Read the true solution
  dims[0] = 4; dims[1] = 5; dims[2] = 7; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/4x3_mult.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 0, mat);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/5x8.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 8; dims[2] = 7; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/5x8_mult_transp.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 1, mat, true);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/2x5.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 2; dims[2] = 7; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/2x5_mult.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 1, mat);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/7x1.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 5; dims[2] = 1; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/7x1_mult_transp.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 2, mat, true);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/1x7.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 5; dims[2] = 1; dims[3] = 11;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/1x7_mult.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 2, mat);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/11x25.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 5; dims[2] = 7; dims[3] = 25;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/11x25_mult_transp.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 3, mat, true);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/17x11.txt");

  // Read the true solution
  dims[0] = 3; dims[1] = 5; dims[2] = 7; dims[3] = 17;
  trueDist = Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  trueSol = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(trueDist);
  TuckerMPI::importTensorBinary("input_files/17x11_mult.mpi",trueSol);

  // Compute the TTM
  mySol = TuckerMPI::ttm(&tensor, 3, mat);

  // Compare the computed solution to the true solution
  approxEqual = isApproxEqual(trueSol, mySol, 1e-10);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(mySol);
  if(!approxEqual)
    return false;

  return true;
}


