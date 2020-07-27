/*
 * Tucker_generate.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: amklinv
 */

#include<random>
#include<chrono>
#include "Tucker.hpp"
#include "Tucker_IO_Util.hpp"

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

  Tucker::SizeArray* I_dims             = Tucker::stringParseSizeArray(fileAsString, "Global dims");
  Tucker::SizeArray* R_dims             = Tucker::stringParseSizeArray(fileAsString, "Ranks");
  std::string out_fns_file              = Tucker::stringParse<std::string>(fileAsString, "Output file list", "rec.txt");
  unsigned seed                         = Tucker::stringParse<unsigned>(fileAsString, "RNG seed",
                                          std::chrono::system_clock::now().time_since_epoch().count());
  double eps                            = Tucker::stringParse<double>(fileAsString, "Noise", 1e-8);

  /////////////////////////////////////////////////
  // Assert that none of the SizeArrays are null //
  /////////////////////////////////////////////////
  if(I_dims == NULL) {
    std::cerr << "Error: Global dims is a required parameter\n";
    return EXIT_FAILURE;
  }
  if(R_dims == NULL) {
    std::cerr << "Error: Ranks is a required parameter\n";
    return EXIT_FAILURE;
  }

  ///////////////////
  // Print options //
  ///////////////////
  if (boolPrintOptions) {
    std::cout << "Global dimensions of the original tensor\n";
    std::cout << "- Global dims = " << *I_dims << std::endl << std::endl;

    std::cout << "Global dimensions of the desired core tensor\n";
    std::cout << "- Ranks = " << *R_dims << std::endl << std::endl;

    std::cout << "File containing a list of filenames to output the constructed data into\n";
    std::cout << "- Output file list = " << out_fns_file << std::endl << std::endl;

    std::cout << "Seed for RNG\n";
    std::cout << "- RNG seed = " << seed << std::endl << std::endl;

    std::cout << "Amount of noise to be added\n";
    std::cout << "- Noise = " << eps << std::endl << std::endl;

    std::cout << "If true, print the parameters\n";
    std::cout << "- Print options = " << (boolPrintOptions ? "true" : "false") << std::endl << std::endl;

    std::cout << std::endl;
  }

  ///////////////////////
  // Check array sizes //
  ///////////////////////
  int nd = I_dims->size();

  if (nd != R_dims->size()) {
    std::cerr << "Error: The size of the core tensor dimension array (" << R_dims->size();
    std::cerr << ") must be equal to the size of the global dimension array ("
        << nd << ")" << std::endl;
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////
  // Create the normal distribution generator //
  // This is a cool C++11 feature             //
  //////////////////////////////////////////////
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution;

  ///////////////////////////////////
  // Generate a random core tensor //
  ///////////////////////////////////
  std::cout << "Generating a random core tensor...\n";
  Tucker::Timer coreTimer;
  coreTimer.start();
  Tucker::TuckerTensor fact(nd);
  fact.G = Tucker::MemoryManager::safe_new<Tucker::Tensor>(*R_dims);
  size_t nnz = R_dims->prod();
  double* dataptr = fact.G->data();
  for(size_t i=0; i<nnz; i++) {
    dataptr[i] = distribution(generator);
  }
  coreTimer.stop();
  std::cout << "Time spent generating core tensor: " << coreTimer.duration() << "s\n";

  ////////////////////////////////
  // Create the factor matrices //
  ////////////////////////////////
  Tucker::Timer factorTimer;
  factorTimer.start();
  for(int d=0; d<nd; d++) {
    std::cout << "Generating factor matrix " << d << "...\n";
    int nrows = (*I_dims)[d];
    int ncols = (*R_dims)[d];
    fact.U[d] = Tucker::MemoryManager::safe_new<Tucker::Matrix>(nrows,ncols);
    nnz = nrows*ncols;
    dataptr = fact.U[d]->data();
    for(size_t i=0; i<nnz; i++) {
      dataptr[i] = distribution(generator);
    }
  }
  factorTimer.stop();
  std::cout << "Time spent generating factor matrices: " << factorTimer.duration() << "s\n";

  ////////////////////////////////////////////////////////
  // Construct the global tensor using a series of TTMs //
  ////////////////////////////////////////////////////////
  Tucker::Tensor* product = fact.G;
  for(int d=0; d<nd; d++) {
    Tucker::Timer ttmTimer;
    ttmTimer.start();
    std::cout << "Performing mode " << d << " TTM...\n";
    Tucker::Tensor* temp = Tucker::ttm(product, d, fact.U[d]);

    if(product != fact.G) {
      Tucker::MemoryManager::safe_delete<Tucker::Tensor>(product);
    }
    product = temp;
    ttmTimer.stop();
    std::cout << "Time spent performing mode " << d << " TTM: " << ttmTimer.duration() << "s\n";
  }

  /////////////////////////////////////////////////////////////////////
  // Compute the norm of the global tensor                           //
  // \todo This could be more efficient; see Bader/Kolda for details //
  /////////////////////////////////////////////////////////////////////
  std::cout << "Computing the global tensor norm...\n";
  Tucker::Timer normTimer;
  normTimer.start();
  double normM = std::sqrt(product->norm2());
  normTimer.stop();
  std::cout << "Time spent computing the tensor norm: " << normTimer.duration() << "s\n";

  ///////////////////////////////////////////////////////////////////
  // Compute the estimated norm of the noise matrix                //
  // The average of each element squared is the standard deviation //
  // squared, so this quantity should be sqrt(nnz * stdev^2)       //
  ///////////////////////////////////////////////////////////////////
  nnz = I_dims->prod();
  double normN = std::sqrt(nnz);

  ///////////////////
  // Compute alpha //
  ///////////////////
  double alpha = eps*normM/normN;

  //////////////////////////////////////////////////////////////////////
  // For each entry of the global tensor, add alpha*randn             //
  // Note that this is space-efficient, as we do not store the entire //
  // noise tensor                                                     //
  //////////////////////////////////////////////////////////////////////
  std::cout << "Adding noise...\n";
  Tucker::Timer noiseTimer;
  noiseTimer.start();
  dataptr = product->data();
  for(size_t i=0; i<nnz; i++) {
    dataptr[i] += alpha*distribution(generator);
  }
  noiseTimer.stop();
  std::cout << "Time spent adding noise: " << noiseTimer.duration() << "s\n";

  //////////////////////////////
  // Write the tensor to disk //
  //////////////////////////////
  std::cout << "Writing tensor to disk...\n";
  Tucker::Timer writeTimer;
  writeTimer.start();
  Tucker::writeTensorBinary(product, out_fns_file.c_str());
  writeTimer.stop();
  std::cout << "Time spent writing to disk: " << writeTimer.duration() << "s\n";

  std::cout << "Done!\n";
  return EXIT_SUCCESS;
}
