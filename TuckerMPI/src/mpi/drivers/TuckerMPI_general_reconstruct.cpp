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

  Tucker::SizeArray* proc_grid_dims     = Tucker::stringParseSizeArray(fileAsString, "Grid dims");
  Tucker::SizeArray* rec_order          = Tucker::stringParseSizeArray(fileAsString, "Reconstruction order");

  std::string sthosvd_dir               = Tucker::stringParse<std::string>(fileAsString, "STHOSVD directory", "compressed");
  std::string sthosvd_fn                = Tucker::stringParse<std::string>(fileAsString, "STHOSVD file prefix", "sthosvd");
  std::string out_fns_file              = Tucker::stringParse<std::string>(fileAsString, "Output file list", "rec.txt");

  /////////////////////////////////////////////////
  // Assert that none of the SizeArrays are null //
  /////////////////////////////////////////////////
  if(proc_grid_dims == NULL) {
    if(rank == 0)
      std::cerr << "Error: Grid dims is a required parameter\n";
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  ///////////////////
  // Print options //
  ///////////////////
  if (rank == 0 && boolPrintOptions) {
    std::cout << "Global dimensions of the processor grid\n";
    std::cout << "- Grid dims = " << *proc_grid_dims << std::endl << std::endl;

    std::cout << "Directory location of ST-HOSVD output files\n";
    std::cout << "- STHOSVD directory = " << sthosvd_dir << std::endl << std::endl;

    std::cout << "Base name of the ST-HOSVD output files\n";
    std::cout << "- STHOSVD file prefix = " << sthosvd_fn << std::endl << std::endl;

    std::cout << "File containing a list of filenames to output the reconstructed data into\n";
    std::cout << "- Output file list = " << out_fns_file << std::endl << std::endl;

    if(rec_order != NULL) {
      std::cout << "Mode order for reconstruction\n";
      std::cout << "NOTE: if left unspecified, the memory-optimal one will be automatically selected\n";
      std::cout << "- Reconstruction order = " << *rec_order << std::endl << std::endl;
    }

    std::cout << "If true, print the parameters\n";
    std::cout << "- Print options = " << (boolPrintOptions ? "true" : "false") << std::endl << std::endl;

    std::cout << std::endl;
  }

  ///////////////////////
  // Check array sizes //
  ///////////////////////
  int nd = proc_grid_dims->size();

  // Does |grid| == nprocs?
  if ((int)proc_grid_dims->prod() != nprocs){
    if (rank==0) {
      std::cerr << "Processor grid dimensions do not multiply to nprocs" << std::endl;
      std::cout << "Processor grid dimensions: " << *proc_grid_dims << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (rec_order != NULL && nd != rec_order->size()) {
    if (rank == 0) {
      std::cerr << "Error: The size of the rec_order array (" << rec_order->size();
      std::cerr << ") must be equal to the size of the processor grid ("
          << nd << ")" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  ////////////////////////////////////
  // Read the core size from a file //
  ////////////////////////////////////
  Tucker::SizeArray* coreSize = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nd);
  if(rank == 0)
  {
    std::string dimFilename = sthosvd_dir + "/" + sthosvd_fn +
        "_ranks.txt";
    std::ifstream ifs(dimFilename);

    if(!ifs.is_open()) {
      if(rank == 0) {
        std::cerr << "Failed to open core size file: " << dimFilename
            << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for(int mode=0; mode<nd; mode++) {
      ifs >> (*coreSize)[mode];
    }
    ifs.close();
  }
  MPI_Bcast(coreSize->data(),nd,MPI_INT,0,MPI_COMM_WORLD);

  //////////////////////////////////////
  // Read the global size from a file //
  //////////////////////////////////////
  Tucker::SizeArray* I_dims = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nd);
  if(rank == 0) {
    std::string sizeFilename = sthosvd_dir + "/" + sthosvd_fn +
        "_size.txt";
    std::ifstream ifs(sizeFilename);

    if(!ifs.is_open()) {
      std::cerr << "Failed to open global size file: " << sizeFilename
                << std::endl;

      return EXIT_FAILURE;
    }

    for(int mode=0; mode<nd; mode++) {
      ifs >> (*I_dims)[mode];
    }
    ifs.close();
  }
  MPI_Bcast(I_dims->data(),nd,MPI_INT,0,MPI_COMM_WORLD);

  //////////////////////////////////////////////
  // Make sure the core size data makes sense //
  //////////////////////////////////////////////
  for(int i=0; i<nd; i++) {
    if((*coreSize)[i] <= 0) {
      if(rank == 0) {
        std::cerr << "coreSize[" << i << "] = " << (*coreSize)[i]
                  << " <= 0\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if((*coreSize)[i] > (*I_dims)[i]) {
      if(rank == 0) {
        std::cerr << "coreSize[" << i << "] = " << (*coreSize)[i]
                  << " > I_dims[" << (*I_dims)[i] << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  /////////////////////////////////
  // Set up factorization object //
  /////////////////////////////////
  TuckerMPI::TuckerTensor* fact =
      Tucker::MemoryManager::safe_new<TuckerMPI::TuckerTensor>(nd);

  /////////////////////////////////////////////
  // Set up distribution object for the core //
  /////////////////////////////////////////////
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*coreSize, *proc_grid_dims);

  ///////////////////////////
  // Read core tensor data //
  ///////////////////////////
  std::string coreFilename = sthosvd_dir + "/" + sthosvd_fn +
            "_core.mpi";
  fact->G = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(dist);
  TuckerMPI::importTensorBinary(coreFilename.c_str(),fact->G);
  if(rank == 0) {
    size_t local_nnz = fact->G->getLocalNumEntries();
    size_t global_nnz = fact->G->getGlobalNumEntries();
    std::cout << "Local core tensor size: " << fact->G->getLocalSize() << ", or ";
    Tucker::printBytes(local_nnz*sizeof(double));
    std::cout << "Global core tensor size: " << fact->G->getGlobalSize() << ", or ";
    Tucker::printBytes(global_nnz*sizeof(double));
  }

  //////////////////////////
  // Read factor matrices //
  //////////////////////////
  for(int mode=0; mode<nd; mode++)
  {
    std::ostringstream ss;
    ss << sthosvd_dir << "/" << sthosvd_fn << "_mat_" << mode << ".mpi";

    fact->U[mode] = Tucker::MemoryManager::safe_new<Tucker::Matrix>((*I_dims)[mode],(*coreSize)[mode]);
    TuckerMPI::importTensorBinary(ss.str().c_str(), fact->U[mode]);
  }

  ///////////////////////////////////////////////////////
  // Scale the appropriate factor matrix               //
  // This step only happens if the scaling file exists //
  ///////////////////////////////////////////////////////
  int scale_mode, scale_size=0;
  Tucker::Vector* scales = NULL;
  Tucker::Vector* shifts = NULL;
  if(rank == 0) {
    std::string scaleFilename = sthosvd_dir + "/" + sthosvd_fn +
        "_scale.txt";
    std::ifstream ifs(scaleFilename);

    // If the scaling file exists...
    if(ifs.is_open()) {
      // Read the mode being scaled
      ifs >> scale_mode;
      std::cout << "Scaling mode " << scale_mode << std::endl;

      // Read the scales and shifts
      scale_size = fact->U[scale_mode]->nrows();
      scales = Tucker::MemoryManager::safe_new<Tucker::Vector>(scale_size);
      shifts = Tucker::MemoryManager::safe_new<Tucker::Vector>(scale_size);
      for(int i=0; i<scale_size; i++) {
        ifs >> (*scales)[i] >> (*shifts)[i];
      }
      ifs.close();
    }
    else {
      std::cerr << "Failed to open scaling and shifting file: " << scaleFilename
          << "\nAssuming no scaling and shifting was performed\n";
      scale_mode = nd;
    }
  }

  MPI_Bcast(&scale_mode,1,MPI_INT,0,MPI_COMM_WORLD);

  // If we are actually scaling
  if(scale_mode < nd) {
    scale_size = fact->U[scale_mode]->nrows();
    if(rank != 0) {
      scales = Tucker::MemoryManager::safe_new<Tucker::Vector>(scale_size);
      shifts = Tucker::MemoryManager::safe_new<Tucker::Vector>(scale_size);
    }

    MPI_Bcast(scales->data(),scale_size,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(shifts->data(),scale_size,MPI_INT,0,MPI_COMM_WORLD);

    // Scale the appropriate factor matrix
    Tucker::premultByDiag(scales, fact->U[scale_mode]);
  }


  ////////////////////////////////////
  // Read the C matrices from files //
  ////////////////////////////////////
  Tucker::SparseMatrix** Cs = Tucker::MemoryManager::safe_new_array<Tucker::SparseMatrix*>(nd);
  for(int mode=0; mode<nd; mode++) {
    std::ostringstream ss;
    ss << "C" << mode << ".txt";

    // If the file exists, read it and replace U with C^T U
    // Otherwise, leave U alone
    std::ifstream ifs(ss.str());
    if(ifs.is_open()) {
      ifs.close();

      // Read C from file
      Cs[mode] = Tucker::importSparseMatrix(ss.str().c_str());

      // Compute C^T U
      Tucker::Matrix* temp = Cs[mode]->multiply(fact->U[mode],true);

      // Replace U with temp
      Tucker::MemoryManager::safe_delete<Tucker::Matrix>(fact->U[mode]);
      fact->U[mode] = temp;
    }
    else {
      Cs[mode] = NULL;
    }
  }

  ////////////////////////////////////////////////////////////
  // Create the optimal reconstruction order if unspecified //
  ////////////////////////////////////////////////////////////
  if(rec_order == NULL) {
    // Create the SizeArray
    rec_order = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nd);
    for(int i=0; i<nd; i++) {
      (*rec_order)[i] = i;
    }

    // Compute the ratios of reconstructed size to core size
    int* rec_size = Tucker::MemoryManager::safe_new_array<int>(nd);
    double* ratios = Tucker::MemoryManager::safe_new_array<double>(nd);
    for(int i=0; i<nd; i++) {
      rec_size[i] = fact->U[i]->nrows();
      ratios[i] = (double)rec_size[i] / (*coreSize)[i];
    }

    // Sort the ratios
    for(int i=1; i<nd; i++) {
      for(int j=0; j<nd-i; j++) {
        if(ratios[j] > ratios[j+1]) {
          std::swap(ratios[j],ratios[j+1]);
          std::swap((*rec_order)[j],(*rec_order)[j+1]);
        }
      }
    }
    if(rank == 0) std::cout << "Reconstruction order: " << *rec_order << std::endl;

    // Free the memory
    Tucker::MemoryManager::safe_delete_array<int>(rec_size,nd);
    Tucker::MemoryManager::safe_delete_array<double>(ratios,nd);
  }

  //////////////////////////////////////////////////////////
  // Make sure the reconstruction order array makes sense //
  //////////////////////////////////////////////////////////
  for(int i=0; i<nd; i++) {
    if((*rec_order)[i] < 0) {
      if(rank == 0) {
        std::cerr << "Error: rec_order[" << i << "] = "
            << (*rec_order)[i] << " < 0\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if((*rec_order)[i] >= nd) {
      if(rank == 0) {
        std::cerr << "Error: rec_order[" << i << "] = "
            << (*rec_order)[i] << " >= nd = " << nd << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for(int j=i+1; j<nd; j++) {
      if((*rec_order)[i] == (*rec_order)[j]) {
        if(rank == 0) {
          std::cerr << "Error: rec_order[" << i << "] == rec_order["
              << j << "] = " << (*rec_order)[i] << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }

  ////////////////////////////////////////////////////
  // Reconstruct the requested pieces of the tensor //
  ////////////////////////////////////////////////////
  TuckerMPI::Tensor* result = fact->G;

  // Compute the nnz of the largest tensor piece being stored by any process
  size_t max_lcl_nnz = 1;
  for(int i=0; i<nd; i++) {
    max_lcl_nnz *= result->getDistribution()->getMap(i,false)->getMaxNumEntries();
  }

  for(int i=0; i<nd; i++)
  {
    int mode = (*rec_order)[i];

    // Perform the TTM
    TuckerMPI::Tensor* temp = TuckerMPI::ttm(result,mode,fact->U[mode],false,0,0,0,0,max_lcl_nnz);

    if(result != fact->G) {
      Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(result);
    }
    result = temp;

    double nrm = result->norm2();
    if(rank == 0) {
      size_t local_nnz = result->getLocalNumEntries();
      size_t global_nnz = result->getGlobalNumEntries();
      std::cout << "Local tensor size after reconstruction iteration "
          << i << ": " << result->getLocalSize() << ", or ";
      Tucker::printBytes(local_nnz*sizeof(double));
      std::cout << "Global tensor size after reconstruction iteration "
          << i << ": " << result->getGlobalSize() << ", or ";
      Tucker::printBytes(global_nnz*sizeof(double));
      std::cout << "Norm: " << nrm << std::endl;
    }
  }

  /////////////////////////////////////
  // Shift, if scaling was performed //
  /////////////////////////////////////
//  if(shifts != NULL) {
//    for(int i=0; i<scale_size; i++) {
//      (*scales)[i] = 1;
//      (*shifts)[i] = -(*shifts)[i];
//    }
//    if(Cs[scale_mode] != NULL) {
//      Tucker::Vector* temp = Cs[scale_mode]->multiply(shifts,true);
//      Tucker::MemoryManager::safe_delete<Tucker::Vector>(shifts);
//      shifts = temp;
//    }
//    int row_begin = result->getDistribution()->getMap(scale_mode,false)->getGlobalIndex(0);
//    if(row_begin >= 0) TuckerMPI::transformSlices(result, scale_mode, scales->data()+row_begin, shifts->data()+row_begin);
//  }

  ////////////////////////////////////////////
  // Write the reconstructed tensor to disk //
  ////////////////////////////////////////////
  TuckerMPI::writeTensorBinary(out_fns_file, *result);

  /////////////////
  // Free memory //
  /////////////////
  Tucker::MemoryManager::safe_delete<TuckerMPI::TuckerTensor>(fact);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(coreSize);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(proc_grid_dims);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(rec_order);
  Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(result);
  if(scales != NULL) Tucker::MemoryManager::safe_delete<Tucker::Vector>(scales);
  if(shifts != NULL) Tucker::MemoryManager::safe_delete<Tucker::Vector>(shifts);
  for(int i=0; i<nd; i++) {
    if(Cs[i] != NULL) {
      Tucker::MemoryManager::safe_delete<Tucker::SparseMatrix>(Cs[i]);
    }
  }
  Tucker::MemoryManager::safe_delete_array<Tucker::SparseMatrix*>(Cs,nd);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  if(rank == 0) {
    Tucker::MemoryManager::printMaxMemUsage();
  }

  //////////////////
  // Finalize MPI //
  //////////////////
  MPI_Finalize();
}
