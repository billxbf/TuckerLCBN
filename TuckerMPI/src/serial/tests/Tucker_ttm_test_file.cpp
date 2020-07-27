/*
 * ttm_test_file.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: amklinv
 */

#include "Tucker.hpp"

int main()
{
  // Read the tensor from a file
  Tucker::Tensor* tensor =
      Tucker::importTensor("input_files/3x5x7x11.txt");

  // Read a matrix to multiply
  Tucker::Matrix* mat =
      Tucker::importMatrix("input_files/3x2.txt");

  // Read the true solution
  Tucker::Tensor* trueSol =
      Tucker::importTensor("input_files/3x2_mult_transp.txt");

  // Compute the TTM
  Tucker::Tensor* mySol = Tucker::ttm(tensor, 0, mat, true);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/4x3.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/4x3_mult.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 0, mat);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/5x8.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/5x8_mult_transp.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 1, mat, true);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/2x5.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/2x5_mult.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 1, mat);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/7x1.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/7x1_mult_transp.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 2, mat, true);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/1x7.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/1x7_mult.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 2, mat);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/11x25.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/11x25_mult_transp.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 3, mat, true);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);

  // Read a matrix to multiply
  mat = Tucker::importMatrix("input_files/17x11.txt");

  // Read the true solution
  trueSol = Tucker::importTensor("input_files/17x11_mult.txt");

  // Compute the TTM
  mySol = Tucker::ttm(tensor, 3, mat);

  // Compare the computed solution to the true solution
  if(!isApproxEqual(trueSol, mySol, 1e-10))
  {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(mySol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(trueSol);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


