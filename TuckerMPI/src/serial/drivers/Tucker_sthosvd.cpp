/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include "Tucker_IO_Util.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "assert.h"

int main(int argc, char* argv[])
{
  Tucker::Timer totalTimer;
  totalTimer.start();

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
  bool boolAuto                         = Tucker::stringParse<bool>(fileAsString, "Automatic rank determination", false);
  bool boolSTHOSVD                      = Tucker::stringParse<bool>(fileAsString, "Perform STHOSVD", false);
  bool boolWriteSTHOSVD                 = Tucker::stringParse<bool>(fileAsString, "Write STHOSVD result", false);
  bool boolPrintOptions                 = Tucker::stringParse<bool>(fileAsString, "Print options", false);
  bool boolWritePreprocessed            = Tucker::stringParse<bool>(fileAsString, "Write preprocessed data", false);

  double tol                            = Tucker::stringParse<double>(fileAsString, "SV Threshold", 1e-6);
  double stdThresh                      = Tucker::stringParse<double>(fileAsString, "STD Threshold", 1e-9);

  Tucker::SizeArray* I_dims             = Tucker::stringParseSizeArray(fileAsString, "Global dims");
  Tucker::SizeArray* R_dims = 0;
  if(!boolAuto)  R_dims                 = Tucker::stringParseSizeArray(fileAsString, "Ranks");

  std::string scaling_type              = Tucker::stringParse<std::string>(fileAsString, "Scaling type", "None");
  std::string sthosvd_dir               = Tucker::stringParse<std::string>(fileAsString, "STHOSVD directory", "compressed");
  std::string sthosvd_fn                = Tucker::stringParse<std::string>(fileAsString, "STHOSVD file prefix", "sthosvd");
  std::string sv_dir                    = Tucker::stringParse<std::string>(fileAsString, "SV directory", ".");
  std::string sv_fn                     = Tucker::stringParse<std::string>(fileAsString, "SV file prefix", "sv");
  std::string in_fns_file               = Tucker::stringParse<std::string>(fileAsString, "Input file list", "raw.txt");
  std::string pre_fns_file              = Tucker::stringParse<std::string>(fileAsString, "Preprocessed output file list", "pre.txt");
  std::string stats_file                = Tucker::stringParse<std::string>(fileAsString, "Stats file", "stats.txt");

  int nd = I_dims->size();
  int scale_mode                        = Tucker::stringParse<int>(fileAsString, "Scale mode", nd-1);


  //
  // Assert that we either have automatic rank determination or the user
  // has supplied their own ranks
  //
  if(!boolAuto && !R_dims) {
    std::cerr << "ERROR: Please either enable Automatic rank determination, "
              << "or provide the desired core tensor size via the Ranks parameter\n";
    return EXIT_FAILURE;
  }

  //
  // Print options
  //
  if (boolPrintOptions) {
    std::cout << "The global dimensions of the tensor to be scaled or compressed\n";
    std::cout << "- Global dims = " << *I_dims << std::endl << std::endl;

    std::cout << "If true, automatically determine rank; otherwise, use the user-defined ranks\n";
    std::cout << "- Automatic rank determination = " << (boolAuto ? "true" : "false") << std::endl << std::endl;

    std::cout << "Used for automatic rank determination; the desired error rate\n";
    std::cout << "- SV Threshold = " << tol << std::endl << std::endl;

    if(!boolAuto) {
      std::cout << "Global dimensions of the desired core tensor\n";
      std::cout << "Not used if \"Automatic rank determination\" is enabled\n";
      std::cout << "- Ranks = " << *R_dims << std::endl << std::endl;
    }

    std::cout << "List of filenames of raw data to be read\n";
    std::cout << "- Input file list = " << in_fns_file << std::endl << std::endl;

    std::cout << "How to scale the tensor\n";
    std::cout << "- Scaling type = " << scaling_type << std::endl << std::endl;

    std::cout << "Which mode's hyperslices will be scaled\n";
    std::cout << "- Scale mode = " << scale_mode << std::endl << std::endl;

    std::cout << "Threshold for standard deviation before we simply set it to 1\n";
    std::cout << "Used in StandardCentering scaling\n";
    std::cout << "- STD Threshold = " << stdThresh << std::endl << std::endl;

    std::cout << "If true, perform ST-HOSVD\n";
    std::cout << "- Perform STHOSVD = " << (boolSTHOSVD ? "true" : "false") << std::endl << std::endl;

    std::cout << "Location of statistics file containing min, max, mean, and std of each hyperslice\n";
    std::cout << "- Stats file = " << stats_file << std::endl << std::endl;

    std::cout << "If true, write the preprocessed data to a file\n";
    std::cout << "- Write preprocessed data = " << (boolWritePreprocessed ? "true" : "false") << std::endl << std::endl;

    std::cout << "File containing a list of filenames to output the scaled data into\n";
    std::cout << "- Preprocessed output file list = " << pre_fns_file << std::endl << std::endl;

    std::cout << "If true, record the result of ST-HOSVD (the core tensor and all factors)\n";
    std::cout << "- Write STHOSVD result = " << (boolWriteSTHOSVD ? "true" : "false") << std::endl << std::endl;

    std::cout << "Directory location of ST-HOSVD output files\n";
    if(boolWriteSTHOSVD) std::cout << "NOTE: Please ensure that this directory actually exists!\n";
    std::cout << "- STHOSVD directory = " << sthosvd_dir << std::endl << std::endl;

    std::cout << "Base name of ST-HOSVD output files\n";
    std::cout << "- STHOSVD file prefix = " << sthosvd_fn << std::endl << std::endl;

    std::cout << "Directory to place singular value files into\n";
    if(boolWriteSTHOSVD) std::cout << "NOTE: Please ensure that this directory actually exists!\n";
    std::cout << "- SV directory = " << sv_dir << std::endl << std::endl;

    std::cout << "Base name for writing the singular value files\n";
    std::cout << "- SV file prefix = " << sv_fn << std::endl << std::endl;

    std::cout << "If true, print the parameters\n";
    std::cout << "- Print options = " << (boolPrintOptions ? "true" : "false") << std::endl << std::endl;

    std::cout << std::endl;
  }

  assert(boolAuto || R_dims->size() == nd);

  ///////////////////////
  // Check array sizes //
  ///////////////////////
  if (!boolAuto && R_dims->size() != 0 && R_dims->size() != nd) {
    std::cerr << "Error: The size of the ranks array (" << R_dims->size();
    std::cerr << ") must be 0 or equal to the size of the global dimensions (" << nd << ")" << std::endl;

    return EXIT_FAILURE;
  }

  ///////////////////////////
  // Read full tensor data //
  ///////////////////////////
  Tucker::Timer readTimer;
  readTimer.start();
  Tucker::Tensor* X = Tucker::MemoryManager::safe_new<Tucker::Tensor>(*I_dims);
  Tucker::readTensorBinary(X,in_fns_file.c_str());
  readTimer.stop();

  size_t nnz = X->getNumElements();
  std::cout << "Input tensor size: " << X->size() << ", or ";
  Tucker::printBytes(nnz*sizeof(double));

  ////////////////////////
  // Compute statistics //
  ////////////////////////
  Tucker::MetricData* metrics = Tucker::computeSliceMetrics(X,
      scale_mode,
      Tucker::MIN+Tucker::MAX+Tucker::MEAN+Tucker::VARIANCE);

  std::cout << "Writing file " << stats_file << std::endl;

  std::ofstream statStream(stats_file);
  statStream << std::setw(5) << "Mode"
      << std::setw(13) << "Mean"
      << std::setw(13) << "Stdev"
      << std::setw(13) << "Min"
      << std::setw(13) << "Max"
      << std::endl;

  for(int i=0; i<X->size(scale_mode); i++) {
    double stdev = sqrt(metrics->getVarianceData()[i]);

    if(stdev < stdThresh) {
      std::cout << "Slice " << i
          << " is below the cutoff. True value is: "
          << stdev << std::endl;
      stdev = 1;
    }

    statStream << std::setw(5) << i
        << std::setw(13) << metrics->getMeanData()[i]
        << std::setw(13) << stdev
        << std::setw(13) << metrics->getMinData()[i]
        << std::setw(13) << metrics->getMaxData()[i] << std::endl;
  }

  statStream.close();

  ///////////////////////////
  // Perform preprocessing //
  ///////////////////////////
  Tucker::Timer preprocessTimer;
  preprocessTimer.start();
  std::string scale_file = sthosvd_dir + "/" + sthosvd_fn +
      "_scale.txt";
  if(scaling_type == "Max") {
    std::cout << "Normalizing the tensor by maximum entry - mode "
              << scale_mode << std::endl;
    normalizeTensorMax(X, scale_mode, scale_file.c_str());
  }
  else if(scaling_type == "MinMax") {
    std::cout << "Normalizing the tensor using minmax scaling - mode "
              << scale_mode << std::endl;
    normalizeTensorMinMax(X, scale_mode, scale_file.c_str());
  }
  else if(scaling_type == "StandardCentering") {
    std::cout << "Normalizing the tensor using standard centering - mode "
              << scale_mode << std::endl;
    normalizeTensorStandardCentering(X, scale_mode, stdThresh, scale_file.c_str());
  }
  else if(scaling_type == "None") {
    std::cout << "Not scaling the tensor\n";
  }
  else {
    std::cerr << "Error: invalid scaling type: " << scaling_type << std::endl;
  }
  preprocessTimer.stop();

  if(boolWritePreprocessed) {
    Tucker::writeTensorBinary(X,pre_fns_file.c_str());
  }

  /////////////////////
  // Perform STHOSVD //
  /////////////////////
  Tucker::Timer sthosvdTimer, writeTimer;
  if(boolSTHOSVD) {
    const Tucker::TuckerTensor* solution;

    sthosvdTimer.start();
    if(boolAuto) {
      solution = Tucker::STHOSVD(X, tol);
    }
    else {
      solution = Tucker::STHOSVD(X, R_dims);
    }
    sthosvdTimer.stop();

    // Write the eigenvalues to files
    std::string filePrefix = sv_dir + "/" + sv_fn + "_mode_";
    Tucker::printEigenvalues(solution, filePrefix);

    double xnorm = std::sqrt(X->norm2());
    double gnorm = std::sqrt(solution->G->norm2());
    std::cout << "Norm of input tensor: " << xnorm << std::endl;
    std::cout << "Norm of core tensor: " << gnorm << std::endl;

    // Compute the error bound based on the eigenvalues
    double eb =0;
    for(int i=0; i<nd; i++) {
      for(int j=solution->G->size(i); j<X->size(i); j++) {
        eb += solution->eigenvalues[i][j];
      }
    }
    std::cout << "Error bound: " << std::sqrt(eb)/xnorm << std::endl;

    writeTimer.start();
    if(boolWriteSTHOSVD) {
      // Write dimension of core tensor
      std::string dimFilename = sthosvd_dir + "/" + sthosvd_fn +
          "_ranks.txt";
      std::cout << "Writing core tensor dimensions to " << dimFilename << std::endl;
      std::ofstream of(dimFilename);
      assert(of.is_open());
      for(int mode=0; mode<nd; mode++) {
        of << solution->G->size(mode) << std::endl;
      }
      of.close();

      // Write dimension of global tensor
      std::string sizeFilename = sthosvd_dir + "/" + sthosvd_fn +
          "_size.txt";
      std::cout << "Writing global tensor dimensions to " << sizeFilename << std::endl;
      of.open(sizeFilename);
      assert(of.is_open());
      for(int mode=0; mode<nd; mode++) {
        of << (*I_dims)[mode] << std::endl;
      }
      of.close();

      // Write core tensor
      std::string coreFilename = sthosvd_dir + "/" + sthosvd_fn +
          "_core.mpi";
      std::cout << "Writing core tensor to " << coreFilename << std::endl;
      Tucker::exportTensorBinary(solution->G, coreFilename.c_str());

      // Write each factor
      for(int mode=0; mode<nd; mode++) {
        // Create the filename by appending the mode #
        std::ostringstream ss;
        ss << sthosvd_dir << "/" << sthosvd_fn << "_mat_" << mode
            << ".mpi";       // Open the file
        std::cout << "Writing factor " << mode << " to " << ss.str() << std::endl;
        Tucker::exportTensorBinary(solution->U[mode], ss.str().c_str());
      }
    }
    writeTimer.stop();
  }

  //
  // Free memory
  //
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(I_dims);
  if(R_dims) Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(R_dims);

  Tucker::MemoryManager::printMaxMemUsage();

  totalTimer.stop();
  std::cout << "Read time: " << readTimer.duration() << std::endl;
  std::cout << "Preprocessing time: " << preprocessTimer.duration() << std::endl;
  std::cout << "STHOSVD time: " << sthosvdTimer.duration() << std::endl;
  std::cout << "Write time: " << writeTimer.duration() << std::endl;
  std::cout << "Total time: " << totalTimer.duration() << std::endl;

  return EXIT_SUCCESS;
}
