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
 * \brief Contains the utility functions supporting the tucker
 * decomposition code
 *
 * @author Alicia Klinvex
 */

#ifndef TUCKER_UTIL_HPP_
#define TUCKER_UTIL_HPP_

#include<string>
#include<iostream>
#include<vector>
#include<sstream>
#include<stdexcept>

namespace Tucker {

void printBytes(size_t bytes);

class MemoryManager {

public:

template<typename T, typename... Args>
static T* safe_new(Args&&... args)
{
  T* allocatedPtr;
  try {
    allocatedPtr = new T(std::forward<Args>(args)...);
  }
  catch(std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
  }
  if(verbose) {
    std::cout << "Allocated ";
    printBytes(sizeof(T));
  }

  curMemUsage += sizeof(T);
  maxMemUsage = std::max(maxMemUsage,curMemUsage);
  return allocatedPtr;
}

/** \brief Allocates memory within a try/catch block
 *
 * \param[in] numToAllocate Number of entries to allocate
 * \exception std::runtime_error \a numToAllocate <= 0
 */
template<typename T>
static T* safe_new_array(const size_t numToAllocate)
{
  if(numToAllocate <= 0) {
    std::ostringstream oss;
    oss << "Tucker::safe_new(const size_t numToAllocate): numToAllocate = "
        << numToAllocate << " <= 0";
    throw std::runtime_error(oss.str());
  }

  T* allocatedPtr;
  try {
    allocatedPtr = new T[numToAllocate];
  }
  catch(std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
  }
  if(verbose) {
    std::cout << "Allocated ";
    printBytes(sizeof(T)*numToAllocate);
  }
  curMemUsage += (sizeof(T)*numToAllocate);
  maxMemUsage = std::max(maxMemUsage,curMemUsage);
  return allocatedPtr;
}

template<typename T>
static void safe_delete(T*& t)
{
  try {
    delete t;
    t=0;
  }
  catch(std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
  }
  if(verbose) {
    std::cout << "Deallocated ";
    printBytes(sizeof(T));
  }
  curMemUsage -= sizeof(T);
}

template<typename T>
static void safe_delete_array(T*& t, const size_t numToDealloc)
{
  try {
    delete[] t;
    t=0;
  }
  catch(std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
  }
  if(verbose) {
    std::cout << "Deallocated ";
    printBytes(numToDealloc*sizeof(T));
  }
  curMemUsage -= (numToDealloc*sizeof(T));
}

static void printCurrentMemUsage()
{
  std::cout << "Current memory usage: ";
  printBytes(curMemUsage);
}

static void printMaxMemUsage()
{
  std::cout << "Maximum local memory usage: ";
  printBytes(maxMemUsage);
}

static size_t curMemUsage;
static size_t maxMemUsage;
static bool verbose;

private:
MemoryManager();

}; // end class MemoryManager

} // end namespace Tucker

#endif /* TUCKER_UTIL_HPP_ */
