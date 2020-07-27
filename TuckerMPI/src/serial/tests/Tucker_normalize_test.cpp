/*
 * normalize_test.cpp
 *
 *  Created on: Sep 1, 2016
 *      Author: amklinv
 */

#include "Tucker.hpp"

int main()
{
  Tucker::Tensor* tensor;
  Tucker::Tensor* true_sol;

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorMinMax(tensor,0);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_mm0.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorMinMax(tensor,1);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_mm1.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorMinMax(tensor,2);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_mm2.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorMinMax(tensor,3);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_mm3.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorStandardCentering(tensor,0,1e-10);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_sc0.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorStandardCentering(tensor,1,1e-10);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_sc1.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorStandardCentering(tensor,2,1e-10);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_sc2.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  // Read tensor from file
  tensor = Tucker::importTensor("input_files/3x5x7x11.txt");

  // Normalize the tensor
  Tucker::normalizeTensorStandardCentering(tensor,3,1e-10);

  // Read true solution from file
  true_sol = Tucker::importTensor("input_files/3x5x7x11_sc3.txt");

  // Compare the computed solution to the true solution
  if(!isApproxEqual(tensor,true_sol,1e-10)) {
    return EXIT_FAILURE;
  }

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(tensor);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(true_sol);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
