/*
 * TuckerMPI_generate.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: amklinv
 */

#include<random>
#include<chrono>
#include "TuckerMPI.hpp"
#include "Tucker_IO_Util.hpp"

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
  Tucker::SizeArray* I_dims             = Tucker::stringParseSizeArray(fileAsString, "Global dims");
  Tucker::SizeArray* R_dims             = Tucker::stringParseSizeArray(fileAsString, "Ranks");
  std::string out_fns_file              = Tucker::stringParse<std::string>(fileAsString, "Output file list", "rec.txt");
  unsigned int seed                     = Tucker::stringParse<unsigned int>(fileAsString, "RNG seed",
                                          std::chrono::system_clock::now().time_since_epoch().count());
  double eps                            = Tucker::stringParse<double>(fileAsString, "Noise", 1e-8);

  /////////////////////////////////////////////////
  // Assert that none of the SizeArrays are null //
  /////////////////////////////////////////////////
  if(I_dims == NULL) {
    if(rank == 0)
      std::cerr << "Error: Global dims is a required parameter\n";
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(R_dims == NULL) {
    if(rank == 0)
      std::cerr << "Error: Ranks is a required parameter\n";
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
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
    std::cout << "Global dimensions of the original tensor\n";
    std::cout << "- Global dims = " << *I_dims << std::endl << std::endl;

    std::cout << "Global dimensions of the processor grid\n";
    std::cout << "- Grid dims = " << *proc_grid_dims << std::endl << std::endl;

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
  if(nd != proc_grid_dims->size()) {
    if (rank==0) {
      std::cerr << "Processor grid dimensions do not multiply to nprocs" << std::endl;
      std::cout << "Processor grid dimensions: " << *proc_grid_dims << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  ///////////////////////////////////////////////////////////////
  // Generate the seeds for each MPI process and scatter them  //
  ///////////////////////////////////////////////////////////////
  int myseed;
  if(rank == 0) {
    unsigned* seeds = Tucker::MemoryManager::safe_new_array<unsigned>(nprocs);

    srand(seed);
    for(int i=0; i<nprocs; i++) {
      seeds[i] = rand();
    }

    MPI_Scatter(seeds,1,MPI_INT,&myseed,1,MPI_INT,0,MPI_COMM_WORLD);

    Tucker::MemoryManager::safe_delete_array<unsigned>(seeds,nprocs);
  }
  else {
    MPI_Scatter(NULL,1,MPI_INT,&myseed,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  //////////////////////////////////////////////
  // Create the normal distribution generator //
  //////////////////////////////////////////////
  std::default_random_engine generator(myseed);
  std::normal_distribution<double> distribution;

  /////////////////////////////////////////////
  // Set up distribution object for the core //
  /////////////////////////////////////////////
  TuckerMPI::Distribution* dist =
      Tucker::MemoryManager::safe_new<TuckerMPI::Distribution>(*R_dims, *proc_grid_dims);

  ///////////////////////////////////
  // Generate a random core tensor //
  ///////////////////////////////////
  if(rank == 0) std::cout << "Generating a random core tensor...\n";
  Tucker::Timer coreTimer;
  coreTimer.start();
  TuckerMPI::TuckerTensor fact(nd);
  fact.G = Tucker::MemoryManager::safe_new<TuckerMPI::Tensor>(dist);
  size_t nnz = dist->getLocalDims().prod();
  double* dataptr = fact.G->getLocalTensor()->data();
  for(size_t i=0; i<nnz; i++) {
    dataptr[i] = distribution(generator);
  }
  coreTimer.stop();
  if(rank == 0) std::cout << "Time spent generating core tensor: " << coreTimer.duration() << "s\n";

  //////////////////////////////////////////////////////////////////
  // Create the factor matrices on process 0, then broadcast them //
  // This could tecnically be a single communication              //
  //////////////////////////////////////////////////////////////////
  Tucker::Timer factorTimer;
  factorTimer.start();
  for(int d=0; d<nd; d++) {
    if(rank == 0) std::cout << "Generating factor matrix " << d << "...\n";
    int nrows = (*I_dims)[d];
    int ncols = (*R_dims)[d];
    fact.U[d] = Tucker::MemoryManager::safe_new<Tucker::Matrix>(nrows,ncols);
    nnz = nrows*ncols;
    dataptr = fact.U[d]->data();
    if(rank == 0) {
      for(size_t i=0; i<nnz; i++) {
        dataptr[i] = distribution(generator);
      }
    }

    MPI_Bcast(dataptr,nnz,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  factorTimer.stop();
  if(rank == 0) std::cout << "Time spent generating factor matrices: " << factorTimer.duration() << "s\n";


  ////////////////////////////////////////////////////////
  // Construct the global tensor using a series of TTMs //
  ////////////////////////////////////////////////////////
  TuckerMPI::Tensor* product = fact.G;
  for(int d=0; d<nd; d++) {
    Tucker::Timer ttmTimer;
    ttmTimer.start();
    if(rank == 0) std::cout << "Performing mode " << d << " TTM...\n";
    TuckerMPI::Tensor* temp = TuckerMPI::ttm(product, d, fact.U[d]);

    if(product != fact.G) {
      Tucker::MemoryManager::safe_delete<TuckerMPI::Tensor>(product);
    }
    product = temp;
    ttmTimer.stop();
    if(rank == 0) std::cout << "Time spent performing mode " << d << " TTM: " << ttmTimer.duration() << "s\n";
  }

  /////////////////////////////////////////////////////////////////////
  // Compute the norm of the global tensor                           //
  // \todo This could be more efficient; see Bader/Kolda for details //
  /////////////////////////////////////////////////////////////////////
  if(rank == 0) std::cout << "Computing the global tensor norm...\n";
  Tucker::Timer normTimer;
  normTimer.start();
  double normM = std::sqrt(product->norm2());
  normTimer.stop();
  if(rank == 0) std::cout << "Time spent computing the tensor norm: " << normTimer.duration() << "s\n";

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
  if(rank == 0) std::cout << "Adding noise...\n";
  Tucker::Timer noiseTimer;
  noiseTimer.start();
  dataptr = product->getLocalTensor()->data();
  nnz = dist->getLocalDims().prod();
  for(size_t i=0; i<nnz; i++) {
    dataptr[i] += alpha*distribution(generator);
  }
  noiseTimer.stop();
  if(rank == 0) std::cout << "Time spent adding noise: " << noiseTimer.duration() << "s\n";

  //////////////////////////////
  // Write the tensor to disk //
  //////////////////////////////
  if(rank == 0) std::cout << "Writing tensor to disk...\n";
  Tucker::Timer writeTimer;
  writeTimer.start();
  TuckerMPI::writeTensorBinary(out_fns_file, *product);
  writeTimer.stop();
  if(rank == 0) std::cout << "Time spent writing to disk: " << writeTimer.duration() << "s\n";

  MPI_Finalize();

  if(rank == 0) std::cout << "Done!\n";
  return EXIT_SUCCESS;
}
