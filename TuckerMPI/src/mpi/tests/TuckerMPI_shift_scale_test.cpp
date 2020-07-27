/*
 * shift_scale_test.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: amklinv
 */

#include "TuckerMPI.hpp"

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


bool runSim(Tucker::SizeArray& procs)
{
  const TuckerMPI::Map* map;
  int index;

  // Create a SizeArray with the dimensions of the tensor
  Tucker::SizeArray dims(4);
  dims[0] = 3; dims[1] = 5; dims[2] = 7; dims[3] = 11;

  std::cout << procs << std::endl;

  // Create a distribution object
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);
  TuckerMPI::Distribution* dist2 =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(dims,procs);

  // Create a tensor with that distribution
  TuckerMPI::Tensor tensor(dist);
  TuckerMPI::Tensor true_sol(dist2);

  // Read the tensor from a binary file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11.mpi",&tensor);

  double shifts0[3] = {-0.450368098480610,
      -0.408899824586946, 0.094037031444121};

  double scales0[3] = {-0.258915944830978,
      0.341369101972144, 0.357212764090606};

  map = tensor.getDistribution()->getMap(0,false);
  index = map->getGlobalIndex(0);

  // Read true solution from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11_ss0.mpi",&true_sol);

  // Call shift-scale
  TuckerMPI::transformSlices(&tensor,0,scales0+index,shifts0+index);

  if(!isApproxEqual(&tensor,&true_sol,1e-10)) {
    return false;
  }

  // Read tensor from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11.mpi",&tensor);

  double shifts1[5] = {0.463612200951355, -0.011100213839396,
      -0.279689899431367, -0.273791359158714, 0.036787804512826};

  double scales1[5] = {0.262109709211147, -0.152432849551241,
      -0.038768240608501, 0.139323762199356, 0.417336040866845};

  map = tensor.getDistribution()->getMap(1,false);
  index = map->getGlobalIndex(0);

  // Read true solution from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11_ss1.mpi",&true_sol);

  // Call shift-scale
  TuckerMPI::transformSlices(&tensor,1,scales1+index,shifts1+index);

  if(!isApproxEqual(&tensor,&true_sol,1e-10)) {
    return false;
  }

  // Read tensor from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11.mpi",&tensor);

  double shifts2[7] = {-0.338427426109669, 0.215635404167474,
      0.077738876192409, -0.066701022790881, 0.384242782631094,
     -0.106948244623087, -0.321024847372268};

  double scales2[7] = {0.133333580320122, 0.124000554312344,
     -0.172058403026751, 0.302965315958237, 0.499477858635892,
      0.480978160932146, -0.372963057805143};

  map = tensor.getDistribution()->getMap(2,false);
  index = map->getGlobalIndex(0);

  // Read true solution from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11_ss2.mpi",&true_sol);

  // Call shift-scale
  TuckerMPI::transformSlices(&tensor,2,scales2+index,shifts2+index);

  if(!isApproxEqual(&tensor,&true_sol,1e-10)) {
    return false;
  }

  // Read tensor from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11.mpi",&tensor);

  double shifts3[11] = {-0.267759854038207, -0.476367533341775,
       0.107432610401855, -0.389190678712850, -0.092540492121601,
       0.384076806661962, 0.048132777476588, -0.130996923288383,
      -0.291654017186653, -0.059056723475676, 0.456196152175878};

  double scales3[11] = {-0.375974083351957, -0.029236754133043,
       0.356896327782193, -0.456609528433047, 0.191625145201306,
       0.478985466675038, -0.216732101507863, -0.366219500005577,
       0.185279684412687, 0.409454555749395, 0.110868982383243};

  map = tensor.getDistribution()->getMap(3,false);
  index = map->getGlobalIndex(0);

  // Read true solution from file
  TuckerMPI::importTensorBinary("input_files/3x5x7x11_ss3.mpi",&true_sol);

  // Call shift-scale
  TuckerMPI::transformSlices(&tensor,3,scales3+index,shifts3+index);

  if(!isApproxEqual(&tensor,&true_sol,1e-10)) {
    return false;
  }

  return true;
}


