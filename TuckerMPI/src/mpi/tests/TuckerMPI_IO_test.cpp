#include "TuckerMPI.hpp"
#include<limits>

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
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank == 0) std::cout << procs << std::endl;

  Tucker::SizeArray sz(4);
  sz[0] = 3; sz[1] = 5; sz[2] = 7; sz[3] = 11;
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(sz, procs);
  TuckerMPI::Tensor Y(dist);
  for(size_t i=0; i<Y.getLocalNumEntries(); i++)
    Y.getLocalTensor()->data()[i] = std::numeric_limits<double>::signaling_NaN();
  std::string filenames = "input_files/input_test_filenames.txt";
  TuckerMPI::readTensorBinary(filenames,Y);

  std::string filename = "input_files/input_test_filename.txt";
  TuckerMPI::Distribution* dist2 =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(sz, procs);
  TuckerMPI::Tensor Y2(dist2);
  for(size_t i=0; i<Y2.getLocalNumEntries(); i++)
    Y2.getLocalTensor()->data()[i] = std::numeric_limits<double>::signaling_NaN();
  TuckerMPI::readTensorBinary(filename,Y2);

  filename = "input_files/output_test_filename.txt";
  TuckerMPI::writeTensorBinary(filename,Y);

  filenames = "input_files/output_test_filenames.txt";
  TuckerMPI::writeTensorBinary(filenames,Y);

  TuckerMPI::Distribution* dist3 =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(sz, procs);
  TuckerMPI::Tensor Y3(dist3);
  for(size_t i=0; i<Y3.getLocalNumEntries(); i++)
    Y3.getLocalTensor()->data()[i] = std::numeric_limits<double>::signaling_NaN();
  TuckerMPI::readTensorBinary(filename,Y3);

  TuckerMPI::Distribution* dist4 =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(sz, procs);
  TuckerMPI::Tensor Y4(dist4);
  for(size_t i=0; i<Y4.getLocalNumEntries(); i++)
    Y4.getLocalTensor()->data()[i] = std::numeric_limits<double>::signaling_NaN();
  TuckerMPI::readTensorBinary(filenames,Y4);

  if(!isApproxEqual(&Y, &Y2, 1e-10)) {
    std::cout << "Y and Y2 are not equal\n";
    return false;
  }
  if(!isApproxEqual(&Y, &Y3, 1e-10)) {
    std::cout << "Y and Y3 are not equal\n";
    return false;
  }
  if(!isApproxEqual(&Y, &Y4, 1e-10)) {
    std::cout << "Y and Y4 are not equal\n";
    return false;
  }

  return true;
}
