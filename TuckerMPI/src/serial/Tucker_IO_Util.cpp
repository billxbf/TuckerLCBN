/** \copyright
 * Copyright (2016) Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 * \n\n
 * BSD 2-Clause License
 * \n\n
 * Copyright (c) 2016, Sandia Corporation
 * All rights reserved.
 * \n\n
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * \n\n
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * \n\n
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * .
 * \n
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @file
 * \brief Contains I/O functions for %Tucker code
 *
 * \authors Woody Austin
 * \authors Grey Ballard
 * \authors Alicia Klinvex
 */

#include "Tucker_IO_Util.hpp"
#include<fstream>
#include<sstream>

namespace Tucker
{

std::string parseString(const int argc, const char* argv[],
    const std::string& cl_arg, const std::string& default_value)
{
  using std::string;

  int arg=1;
  string tmp;
  while (arg < argc) {
    if (cl_arg == string(argv[arg])) {
      // get next cl_arg
      arg++;
      if (arg >= argc)
        return "";
      // convert to string
      tmp = string(argv[arg]);
      // return tkr_real if everything is OK
      return tmp;
    }
    arg++;
  }
  // return default value if not specified on command line
  return default_value;
}

std::vector<std::string> getFileAsStrings(const std::string& paramfn) {
  std::string line;
  std::vector<std::string> fileLines;
  std::ifstream myFile(paramfn);

  if (!myFile) { std::cerr << "Error opening parameter file: " << paramfn << std::endl; }

  while (std::getline(myFile, line)) {
    // Lines starting with # are comments. The second part tests to see if the line is blank
    if (!(line.find_first_of("#") == 0) && line.find_first_not_of(" \t") != std::string::npos) {
      fileLines.push_back(line);
    }
  }

  return fileLines;
}

SizeArray* stringParseSizeArray(const std::vector<std::string>& lines,
    const std::string& keyword) {
  std::vector<int> tmp; // Placeholder for values in IndxArray
  int value;

  for (auto line : lines) {
    // If the keyword is in the string then use that value
    if (line.find(keyword) != std::string::npos) {
      // Find the equal sign
      std::size_t equalPos = line.find("=");
      std::stringstream valueStream(line.substr(equalPos+1));

      // The value should be one "word", extract it from the string
      while(valueStream >> value) {
        tmp.push_back(value);
      }
    }
  }

  assert(tmp.size() <= std::numeric_limits<int>::max());

  if(tmp.size() == 0)
    return NULL;

  // Copy tmp vector into the IndxArray
  SizeArray* arr = MemoryManager::safe_new<SizeArray>((int)tmp.size());
  for (int i = 0; i < (int)tmp.size(); i++) {
    (*arr)[i] = tmp[i];
  }

  return arr; // Returns empty array if nothing is ever pushed onto tmp vector
}

void printEigenvalues(const TuckerTensor* factorization,
    const std::string& filePrefix)
{
  // For each mode...
  int nmodes = factorization->N;
  for(int mode=0; mode<nmodes; mode++) {
    // Create the filename by appending the mode #
    std::ostringstream ss;
    ss << filePrefix << mode << ".txt";

    // Open the file
    std::ofstream ofs(ss.str());

    // Determine the number of eigenvalues for this mode
    int nevals = factorization->U[mode]->nrows();
    for(int i=0; i<nevals; i++) {
      ofs << factorization->eigenvalues[mode][i] << std::endl;
    }

    ofs.close();
  }
}

}
