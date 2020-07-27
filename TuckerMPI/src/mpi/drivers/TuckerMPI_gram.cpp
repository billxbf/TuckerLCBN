/*
 * TuckerMPI_gram.cpp
 *
 *  Created on: May 8, 2017
 *      Author: amklinv
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
  int mode                              = Tucker::stringParse<int>(fileAsString, "Mode", 0);
  std::string timing_file               = Tucker::stringParse<std::string>(fileAsString, "Timing file", "runtime.csv");
  Tucker::SizeArray* I_dims             = Tucker::stringParseSizeArray(fileAsString, "Global dims");
  Tucker::SizeArray* proc_grid_dims     = Tucker::stringParseSizeArray(fileAsString, "Grid dims");
  std::string in_fns_file               = Tucker::stringParse<std::string>(fileAsString, "Input file list", "raw.txt");

  int nd = I_dims->size();

  //
  // Print options
  //
  if (rank == 0 && boolPrintOptions) {
    std::cout << "The global dimensions of the tensor to be scaled or compressed\n";
    std::cout << "- Global dims = " << *I_dims << std::endl << std::endl;

    std::cout << "The global dimensions of the processor grid\n";
    std::cout << "- Grid dims = " << *proc_grid_dims << std::endl << std::endl;

    std::cout << "If true, use the old Gram algorithm; otherwise use the new one\n";
    std::cout << "- Use old Gram = " << (boolUseOldGram ? "true" : "false") << std::endl << std::endl;

    std::cout << "Mode for Gram computation\n";
    std::cout << "- Mode = " << mode << std::endl << std::endl;

    std::cout << "Name of the CSV file holding the timing results\n";
    std::cout << "- Timing file = " << timing_file << std::endl << std::endl;

    std::cout << "List of filenames of raw data to be read\n";
    std::cout << "- Input file list = " << in_fns_file << std::endl << std::endl;

    std::cout << "If true, print the parameters\n";
    std::cout << "- Print options = " << (boolPrintOptions ? "true" : "false") << std::endl << std::endl;

    std::cout << std::endl;
  }

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

  ////////////////////////////////
  // Set up distribution object //
  ////////////////////////////////
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*I_dims, *proc_grid_dims);

  ///////////////////////////
  // Read full tensor data //
  ///////////////////////////
  Tucker::Timer readTimer;
  MPI_Barrier(MPI_COMM_WORLD);
  readTimer.start();
  TuckerMPI::Tensor X(dist);
  TuckerMPI::readTensorBinary(in_fns_file,X);
  MPI_Barrier(MPI_COMM_WORLD);
  readTimer.stop();

  if(rank == 0) {
    std::cout << "Time to read tensor: " << readTimer.duration() << std::endl;
  }

  ///////////////////////////////////////
  // Compute the specified Gram matrix //
  ///////////////////////////////////////
  Tucker::Timer mult_timer, shift_timer, allreduce_timer, allgather_timer,
                  pack_timer, alltoall_timer, unpack_timer, total_timer;

  Tucker::Matrix* temp;

  MPI_Barrier(MPI_COMM_WORLD);
  total_timer.start();
  if(boolUseOldGram) {
    temp = TuckerMPI::oldGram(&X, mode, &mult_timer, &shift_timer,
        &allreduce_timer, & allgather_timer);
  }
  else {
    temp = TuckerMPI::newGram(&X, mode, &mult_timer, &pack_timer,
        &alltoall_timer, &unpack_timer, &allreduce_timer);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  total_timer.stop();

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(temp);

  if(rank == 0) {
    std::cout << "Gram runtime: " << total_timer.duration() << std::endl;
  }

  ////////////////////////////////////////
  // Determine the maximum memory usage //
  ////////////////////////////////////////
  size_t max_mem = Tucker::MemoryManager::maxMemUsage;

  if(rank == 0) {
    Tucker::MemoryManager::printMaxMemUsage();
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Gather the resource usage information to process 0 and write it to a file //
  ///////////////////////////////////////////////////////////////////////////////
  double mult_time = mult_timer.duration();
  double shift_time = shift_timer.duration();
  double allreduce_time = allreduce_timer.duration();
  double allgather_time = allgather_timer.duration();
  double pack_time = pack_timer.duration();
  double alltoall_time = alltoall_timer.duration();
  double unpack_time = unpack_timer.duration();
  double total_time = total_timer.duration();

  double *mult_times, *shift_times, *allreduce_times, *allgather_times,
      *pack_times, *alltoall_times, *unpack_times, *total_times;
  size_t* max_mems;

  if(rank == 0) {
    mult_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    shift_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    allreduce_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    allgather_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    pack_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    alltoall_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    unpack_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    total_times = Tucker::MemoryManager::safe_new_array<double>(nprocs);
    max_mems = Tucker::MemoryManager::safe_new_array<size_t>(nprocs);
  }

  MPI_Gather(&mult_time, 1, MPI_DOUBLE, mult_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&shift_time, 1, MPI_DOUBLE, shift_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&allreduce_time, 1, MPI_DOUBLE, allreduce_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&allgather_time, 1, MPI_DOUBLE, allgather_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&pack_time, 1, MPI_DOUBLE, pack_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&alltoall_time, 1, MPI_DOUBLE, alltoall_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&unpack_time, 1, MPI_DOUBLE, unpack_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&total_time, 1, MPI_DOUBLE, total_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&max_mem, sizeof(size_t), MPI_CHAR, max_mems, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);

  if(rank == 0) {
    // Send the data to a file
    std::ofstream os(timing_file);

    // Create the header row
    os << "mult,shift,allreduce,allgather,pack,alltoall,unpack,total,mem\n";

    // For each MPI process
    for(int r=0; r<nprocs; r++) {
      os << mult_times[r] << "," << shift_times[r] << "," << allreduce_times[r] << ","
         << allgather_times[r] << "," << pack_times[r] << "," << alltoall_times[r] << ","
         << unpack_times[r] << "," << total_times[r] << "," << max_mems[r] << std::endl;
    }

    os.close();

    Tucker::MemoryManager::safe_delete_array<double>(mult_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(shift_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(allreduce_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(allgather_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(pack_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(alltoall_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(unpack_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<double>(total_times,nprocs);
    Tucker::MemoryManager::safe_delete_array<size_t>(max_mems,nprocs);
  }

  /////////////////
  // Free memory //
  /////////////////
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(proc_grid_dims);

  //////////////////
  // Finalize MPI //
  //////////////////
  MPI_Finalize();
}
