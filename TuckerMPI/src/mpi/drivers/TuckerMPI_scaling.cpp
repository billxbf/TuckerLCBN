/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "TuckerMPI.hpp"
#include "Tucker.hpp"
#include "Tucker_IO_Util.hpp"
#include "TuckerMPI_IO_Util.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "assert.h"

int main(int argc, char* argv[])
{
  //
  // Initialize MPI
  //
  MPI_Init(&argc, &argv);

  //
  // Get the rank of this MPI process
  // Only rank 0 will print to stdout
  //
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  //
  // Get the name of the input file
  //
  std::string paramfn = Tucker::parseString(argc, (const char**)argv,
      "--parameter-file", "paramfile.txt");

  //
  // Parse parameter file
  // Put's each line as a string into a vector ignoring empty lines
  // and comments
  //
  std::vector<std::string> fileAsString = Tucker::getFileAsStrings(paramfn);
  bool boolPrintOptions                 = Tucker::stringParse<bool>(fileAsString, "Print options", false);
  bool boolUseOldGram                   = Tucker::stringParse<bool>(fileAsString, "Use old Gram", false);

  Tucker::SizeArray* I_dims             = Tucker::stringParseSizeArray(fileAsString, "Global dims");
  Tucker::SizeArray* R_dims             = Tucker::stringParseSizeArray(fileAsString, "Ranks");
  Tucker::SizeArray* proc_grid_dims     = Tucker::stringParseSizeArray(fileAsString, "Grid dims");

  std::string timing_file               = Tucker::stringParse<std::string>(fileAsString, "Timing file", "runtime.csv");

  int nd = I_dims->size();

  //
  // Assert that the user has supplied their own ranks
  //
  if(!R_dims) {
    std::cerr << "ERROR: Please provide the desired core tensor size via the Ranks parameter\n";
    return EXIT_FAILURE;
  }

  //
  // Print options
  //
  if (rank == 0 && boolPrintOptions) {
    std::cout << "The global dimensions of the tensor to be scaled or compressed\n";
    std::cout << "- Global dims = " << *I_dims << std::endl << std::endl;

    std::cout << "The global dimensions of the processor grid\n";
    std::cout << "- Grid dims = " << *proc_grid_dims << std::endl << std::endl;

    std::cout << "Global dimensions of the desired core tensor\n";
    std::cout << "- Ranks = " << *R_dims << std::endl << std::endl;

    std::cout << "If true, use the old Gram algorithm; otherwise use the new one\n";
    std::cout << "- Use old Gram = " << (boolUseOldGram ? "true" : "false") << std::endl << std::endl;

    std::cout << "Name of the CSV file holding the timing results\n";
    std::cout << "- Timing file = " << timing_file << std::endl << std::endl;

    std::cout << "If true, print the parameters\n";
    std::cout << "- Print options = " << (boolPrintOptions ? "true" : "false") << std::endl << std::endl;

    std::cout << std::endl;
  }

  assert(R_dims->size() == nd);

  ///////////////////////
  // Check array sizes //
  ///////////////////////

  // Does |grid| == nprocs?
  if ((int)proc_grid_dims->prod() != nprocs){
    if (rank==0) {
      std::cerr << "Processor grid dimensions do not multiply to nprocs" << std::endl;
      std::cout << "Processor grid dimensions: " << *proc_grid_dims << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (nd != proc_grid_dims->size()) {
    if (rank == 0) {
      std::cerr << "Error: The size of global dimension array (" << nd;
      std::cerr << ") must be equal to the size of the processor grid ("
          << proc_grid_dims->size() << ")" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (R_dims->size() != 0 && R_dims->size() != nd) {
    if (rank == 0) {
      std::cerr << "Error: The size of the ranks array (" << R_dims->size();
      std::cerr << ") must be 0 or equal to the size of the processor grid (" << nd << ")" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  ///////////////////////////
  // Set up processor grid //
  ///////////////////////////
  if (rank == 0) {
    std::cout << "Creating process grid" << std::endl;
  }
  ////////////////////////////////
  // Set up distribution object //
  ////////////////////////////////
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*I_dims, *proc_grid_dims);

  /////////////////////////////////
  // Generate random tensor data //
  /////////////////////////////////
  TuckerMPI::Tensor X(dist);
  X.rand();

  if(rank == 0) {
    size_t local_nnz = X.getLocalNumEntries();
    size_t global_nnz = X.getGlobalNumEntries();
    std::cout << "Local input tensor size: " << X.getLocalSize() << ", or ";
    Tucker::printBytes(local_nnz*sizeof(double));
    std::cout << "Global input tensor size: " << X.getGlobalSize() << ", or ";
    Tucker::printBytes(global_nnz*sizeof(double));
  }

  /////////////////////
  // Perform STHOSVD //
  /////////////////////
  const TuckerMPI::TuckerTensor* solution = TuckerMPI::STHOSVD(&X, R_dims, boolUseOldGram);

  // Send the timing information to a CSV
  solution->printTimers(timing_file);

  double xnorm = std::sqrt(X.norm2());
  double gnorm = std::sqrt(solution->G->norm2());
  if(rank == 0) {
    std::cout << "Norm of input tensor: " << xnorm << std::endl;
    std::cout << "Norm of core tensor: " << gnorm << std::endl;

    // Compute the error bound based on the eigenvalues
    double eb =0;
    for(int i=0; i<nd; i++) {
      for(int j=solution->G->getGlobalSize(i); j<X.getGlobalSize(i); j++) {
        eb += solution->eigenvalues[i][j];
      }
    }
    std::cout << "Error bound: " << std::sqrt(eb)/xnorm << std::endl;
  }

  //
  // Free memory
  //
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  if(R_dims) Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(R_dims);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(proc_grid_dims);

  if(rank == 0) {
    Tucker::MemoryManager::printMaxMemUsage();
  }

  // Finalize MPI
  MPI_Finalize();
}
