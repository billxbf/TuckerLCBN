/**
 * @file
 * \example
 *
 * \author Alicia Klinvex
 */

#include <cmath>
#include "Tucker.hpp"

int main()
{
  const double TRUE_SOLUTION = 9.690249359274157;

  // Read the matrix from a binary file
  Tucker::Tensor* tensor =
      Tucker::importTensor("input_files/3x5x7x11.txt");

  // Compute its norm
  double norm = sqrt(tensor->norm2());

  // Free memory
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);

  // Compare computed solution to true solution
  double diff = std::abs(norm-TRUE_SOLUTION);
  if(diff < 1e-10)
    return EXIT_SUCCESS;
  else {
    std::cerr << "ERROR: The true solution is " << TRUE_SOLUTION
              << ", but the computed solution was " << norm
              << ", a difference of " << diff << std::endl;
  }

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_FAILURE;
}

