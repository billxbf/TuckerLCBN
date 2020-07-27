/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include <iostream>
#include <cstdlib>

bool checkUTEqual(const double* arr1, const double* arr2, int numRows);

int main()
{
  Tucker::SizeArray* size =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(4);
  for(int i=0; i<4; i++)
    (*size)[i] = 2;
  Tucker::Tensor* t =
      Tucker::MemoryManager::safe_new<Tucker::Tensor>(*size);
  double* data = t->data();
  for(int i=0; i<16; i++)
    data[i] = i+1;

  Tucker::Matrix* Gram;
  Gram = Tucker::computeGram(t,0);
  double trueSol0[4];
  trueSol0[0] = 680; trueSol0[1] = 744;
  trueSol0[2] = 744; trueSol0[3] = 816;
  bool isEqual = checkUTEqual(Gram->data(), trueSol0, 2);
  if(!isEqual) return EXIT_FAILURE;
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(Gram);

  Gram = Tucker::computeGram(t,1);
  double trueSol1[4];
  trueSol1[0] = 612; trueSol1[1] = 732;
  trueSol1[2] = 732; trueSol1[3] = 884;
  isEqual = checkUTEqual(Gram->data(), trueSol1, 2);
  if(!isEqual) return EXIT_FAILURE;
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(Gram);

  Gram = Tucker::computeGram(t,2);
  double trueSol2[4];
  trueSol2[0] = 476; trueSol2[1] = 684;
  trueSol2[2] = 684; trueSol2[3] = 1020;
  isEqual = checkUTEqual(Gram->data(), trueSol2, 2);
  if(!isEqual) return EXIT_FAILURE;
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(Gram);

  Gram = Tucker::computeGram(t,3);
  double trueSol3[4];
  trueSol3[0] = 204; trueSol3[1] = 492;
  trueSol3[2] = 492; trueSol3[3] = 1292;
  isEqual = checkUTEqual(Gram->data(), trueSol3, 2);
  if(!isEqual) return EXIT_FAILURE;
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(Gram);

  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(t);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(size);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

bool checkUTEqual(const double* arr1, const double* arr2, int numRows)
{
  for(int r=0; r<numRows; r++) {
    for(int c=r; c<numRows; c++) {
      int ind = r+c*numRows;
      if(arr1[ind] != arr2[ind]) {
        return false;
      }
    }
  }
  return true;
}
