/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include "Tucker_Vector.hpp"
#include "Tucker_IO_Util.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "assert.h"

int main(int argc, char* argv[])
{
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

  Tucker::SizeArray* rec_order          = Tucker::stringParseSizeArray(fileAsString, "Reconstruction order");

  std::string sthosvd_dir               = Tucker::stringParse<std::string>(fileAsString, "STHOSVD directory", "compressed");
  std::string sthosvd_fn                = Tucker::stringParse<std::string>(fileAsString, "STHOSVD file prefix", "sthosvd");
  std::string out_fns_file              = Tucker::stringParse<std::string>(fileAsString, "Output file list", "rec.txt");

  ///////////////////
  // Print options //
  ///////////////////
  if (boolPrintOptions) {

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

  ///////////////////////////////////
  // Determine the number of modes //
  ///////////////////////////////////
  std::string dimFilename = sthosvd_dir + "/" + sthosvd_fn +
      "_ranks.txt";
  std::ifstream ifs(dimFilename);

  if(!ifs.is_open()) {
    std::cerr << "Failed to open core size file: " << dimFilename
              << std::endl;

    return EXIT_FAILURE;
  }

  int nd=-1;
  do {
    int temp;
    ifs >> temp;
    nd++;
  } while(!ifs.eof());
  ifs.close();


  ////////////////////////////////////
  // Read the core size from a file //
  ////////////////////////////////////
  Tucker::SizeArray* coreSize = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nd);
  ifs.open(dimFilename);

  for(int mode=0; mode<nd; mode++) {
    ifs >> (*coreSize)[mode];
  }
  ifs.close();

  //////////////////////////////////////
  // Read the global size from a file //
  //////////////////////////////////////
  Tucker::SizeArray* I_dims = Tucker::MemoryManager::safe_new<Tucker::SizeArray>(nd);
  std::string sizeFilename = sthosvd_dir + "/" + sthosvd_fn +
      "_size.txt";
  ifs.open(sizeFilename);

  if(!ifs.is_open()) {
    std::cerr << "Failed to open global size file: " << sizeFilename
              << std::endl;

    return EXIT_FAILURE;
  }

  for(int mode=0; mode<nd; mode++) {
    ifs >> (*I_dims)[mode];
  }
  ifs.close();

  //////////////////////////////////////////////
  // Make sure the core size data makes sense //
  //////////////////////////////////////////////
  for(int i=0; i<nd; i++) {
    if((*coreSize)[i] <= 0) {
      std::cerr << "coreSize[" << i << "] = " << (*coreSize)[i]
                << " <= 0\n";

      return EXIT_FAILURE;
    }

    if((*coreSize)[i] > (*I_dims)[i]) {
      std::cerr << "coreSize[" << i << "] = " << (*coreSize)[i]
                << " > I_dims[" << (*I_dims)[i] << std::endl;

      return EXIT_FAILURE;
    }
  }

  /////////////////////////////////
  // Set up factorization object //
  /////////////////////////////////
  Tucker::TuckerTensor* fact =
      Tucker::MemoryManager::safe_new<Tucker::TuckerTensor>(nd);

  ///////////////////////////
  // Read core tensor data //
  ///////////////////////////
  Tucker::Timer* readTimer = Tucker::MemoryManager::safe_new<Tucker::Timer>();
  readTimer->start();
  std::string coreFilename = sthosvd_dir + "/" + sthosvd_fn +
      "_core.mpi";
  fact->G = Tucker::MemoryManager::safe_new<Tucker::Tensor>(*coreSize);
  Tucker::importTensorBinary(fact->G, coreFilename.c_str());

  size_t nnz = fact->G->getNumElements();
  std::cout << "Core tensor size: " << fact->G->size() << ", or ";
  Tucker::printBytes(nnz*sizeof(double));

  //////////////////////////
  // Read factor matrices //
  //////////////////////////
  for(int mode=0; mode<nd; mode++)
  {
    std::ostringstream ss;
    ss << sthosvd_dir << "/" << sthosvd_fn << "_mat_" << mode << ".mpi";

    fact->U[mode] = Tucker::MemoryManager::safe_new<Tucker::Matrix>((*I_dims)[mode],(*coreSize)[mode]);
    Tucker::importTensorBinary(fact->U[mode], ss.str().c_str());
  }

  ///////////////////////////////////////////////////////
  // Scale the appropriate factor matrix               //
  // This step only happens if the scaling file exists //
  ///////////////////////////////////////////////////////
  Tucker::Timer* scaleTimer = Tucker::MemoryManager::safe_new<Tucker::Timer>();
  scaleTimer->start();
  std::string scaleFilename = sthosvd_dir + "/" + sthosvd_fn +
      "_scale.txt";
  ifs.open(scaleFilename);

  Tucker::Vector* scales = NULL;
  Tucker::Vector* shifts = NULL;
  int scale_size = 0;
  int scale_mode;

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

    // Scale the appropriate factor matrix
    Tucker::premultByDiag(scales, fact->U[scale_mode]);
  }
  else {
    std::cerr << "Failed to open scaling and shifting file: " << scaleFilename
        << "\nAssuming no scaling and shifting was performed\n";
  }
  scaleTimer->stop();
  std::cout << "Time spent scaling: " << scaleTimer->duration() << "s\n";
  Tucker::MemoryManager::safe_delete<Tucker::Timer>(scaleTimer);


  ////////////////////////////////////
  // Read the C matrices from files //
  ////////////////////////////////////
  Tucker::SparseMatrix** Cs = Tucker::MemoryManager::safe_new_array<Tucker::SparseMatrix*>(nd);
  for(int mode=0; mode<nd; mode++) {
    std::ostringstream ss;
    ss << "C" << mode << ".txt";

    // If the file exists, read it and replace U with C^T U
    // Otherwise, leave U alone
    ifs.open(ss.str());
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
  readTimer->stop();
  std::cout << "Time spent reading: " << readTimer->duration() << "s\n";
  Tucker::MemoryManager::safe_delete<Tucker::Timer>(readTimer);

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
    for(int mode=0; mode<nd; mode++) {
      rec_size[mode] = fact->U[mode]->nrows();
      ratios[mode] = (double)rec_size[mode] / (*coreSize)[mode];
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

    std::cout << "Reconstruction order: " << *rec_order << std::endl;

    // Free the memory
    Tucker::MemoryManager::safe_delete_array<int>(rec_size,nd);
    Tucker::MemoryManager::safe_delete_array<double>(ratios,nd);
  }

  //////////////////////////////////////////////////////////
  // Make sure the reconstruction order array makes sense //
  //////////////////////////////////////////////////////////
  for(int i=0; i<nd; i++) {
    if((*rec_order)[i] < 0) {
      std::cerr << "Error: rec_order[" << i << "] = "
          << (*rec_order)[i] << " < 0\n";

      return EXIT_FAILURE;
    }

    if((*rec_order)[i] >= nd) {
        std::cerr << "Error: rec_order[" << i << "] = "
            << (*rec_order)[i] << " >= nd = " << nd << std::endl;

      return EXIT_FAILURE;
    }

    for(int j=i+1; j<nd; j++) {
      if((*rec_order)[i] == (*rec_order)[j]) {
        std::cerr << "Error: rec_order[" << i << "] == rec_order["
            << j << "] = " << (*rec_order)[i] << std::endl;

        return EXIT_FAILURE;
      }
    }
  }

  ////////////////////////////////////////////////////
  // Reconstruct the requested pieces of the tensor //
  ////////////////////////////////////////////////////
  Tucker::Timer* reconstructTimer = Tucker::MemoryManager::safe_new<Tucker::Timer>();
  reconstructTimer->start();
  Tucker::Tensor* result = fact->G;
  for(int i=0; i<nd; i++)
  {
    int mode = (*rec_order)[i];

    // Perform the TTM
    Tucker::Tensor* temp = Tucker::ttm(result,mode,fact->U[mode]);

    if(result != fact->G) {
      Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);
    }
    result = temp;

    size_t nnz = result->getNumElements();
    std::cout << "Tensor size after reconstruction iteration "
        << i << ": " << result->size() << ", or ";
    Tucker::printBytes(nnz*sizeof(double));
  }
  reconstructTimer->stop();
  std::cout << "Time spent reconstructing: " << reconstructTimer->duration() << "s\n";
  Tucker::MemoryManager::safe_delete<Tucker::Timer>(reconstructTimer);

  /////////////////////////////////////
  // Shift, if scaling was performed //
  /////////////////////////////////////
  if(shifts != NULL) {
    for(int i=0; i<scale_size; i++) {
      (*scales)[i] = 1;
      (*shifts)[i] = -(*shifts)[i];
    }
    if(Cs[scale_mode] != NULL) {
      Tucker::Vector* temp = Cs[scale_mode]->multiply(shifts,true);
      Tucker::MemoryManager::safe_delete<Tucker::Vector>(shifts);
      shifts = temp;
    }
    Tucker::transformSlices(result, scale_mode, scales->data(), shifts->data());
  }

  ////////////////////////////////////////////
  // Write the reconstructed tensor to disk //
  ////////////////////////////////////////////
  Tucker::Timer* writeTimer = Tucker::MemoryManager::safe_new<Tucker::Timer>();
  writeTimer->start();
  Tucker::writeTensorBinary(result, out_fns_file.c_str());
  writeTimer->stop();
  std::cout << "Time spent writing: " << writeTimer->duration() << "s\n";
  Tucker::MemoryManager::safe_delete<Tucker::Timer>(writeTimer);

  /////////////////
  // Free memory //
  /////////////////
  Tucker::MemoryManager::safe_delete<Tucker::TuckerTensor>(fact);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(coreSize);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(rec_order);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);
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
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::printMaxMemUsage();

  return EXIT_SUCCESS;
}
