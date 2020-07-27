/*
 * Tucker_Util.cpp
 *
 *  Created on: Dec 8, 2016
 *      Author: amklinv
 */

#include "Tucker_Util.hpp"
#include<iomanip>

namespace Tucker {

void printBytes(size_t bytes)
{
  const size_t KB = 1e3;
  const size_t MB = 1e6;
  const size_t GB = 1e9;
  const size_t TB = 1e12;

  if(bytes > TB) {
    std::cout << std::setprecision(5) << bytes / (double)TB << " TB\n";
  }
  else if(bytes > GB) {
    std::cout << std::setprecision(5) << bytes / (double)GB << " GB\n";
  }
  else if(bytes > MB) {
    std::cout << std::setprecision(5) << bytes / (double)MB << " MB\n";
  }
  else if(bytes > KB) {
    std::cout << std::setprecision(5) << bytes / (double)KB << " KB\n";
  }
  else {
    std::cout << bytes << " bytes\n";
  }
}

bool MemoryManager::verbose = false;
size_t MemoryManager::curMemUsage = 0;
size_t MemoryManager::maxMemUsage = 0;

} // end namespace Tucker
