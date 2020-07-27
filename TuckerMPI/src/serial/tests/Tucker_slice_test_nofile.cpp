/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include <iostream>
#include <cstdlib>

int main()
{
  Tucker::SizeArray* size =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(4);
  (*size)[0] = 2;
  (*size)[1] = 3;
  (*size)[2] = 5;
  (*size)[3] = 7;
  Tucker::Tensor* t =
      Tucker::MemoryManager::safe_new<Tucker::Tensor>(*size);

  size_t nnz = size->prod();
  double* data = t->data();
  for(size_t i=0; i<nnz; i++) {
    data[i] = (double)i+1;
  }

  double trueData[7][4][3];
  trueData[0][0][0] = 209;
  trueData[0][0][1] = 1;
  trueData[0][0][2] = 11025;
  trueData[1][0][0] = 210;
  trueData[1][0][1] = 2;
  trueData[1][0][2] = 11130;
  trueData[0][1][0] = 206;
  trueData[0][1][1] = 1;
  trueData[0][1][2] = 7245;
  trueData[1][1][0] = 208;
  trueData[1][1][1] = 3;
  trueData[1][1][2] = 7385;
  trueData[2][1][0] = 210;
  trueData[2][1][1] = 5;
  trueData[2][1][2] = 7525;
  trueData[0][2][0] = 186;
  trueData[0][2][1] = 1;
  trueData[0][2][2] = 3927;
  trueData[1][2][0] = 192;
  trueData[1][2][1] = 7;
  trueData[1][2][2] = 4179;
  trueData[2][2][0] = 198;
  trueData[2][2][1] = 13;
  trueData[2][2][2] = 4431;
  trueData[3][2][0] = 204;
  trueData[3][2][1] = 19;
  trueData[3][2][2] = 4683;
  trueData[4][2][0] = 210;
  trueData[4][2][1] = 25;
  trueData[4][2][2] = 4935;
  trueData[0][3][0] = 30;
  trueData[0][3][1] = 1;
  trueData[0][3][2] = 465;
  trueData[1][3][0] = 60;
  trueData[1][3][1] = 31;
  trueData[1][3][2] = 1365;
  trueData[2][3][0] = 90;
  trueData[2][3][1] = 61;
  trueData[2][3][2] = 2265;
  trueData[3][3][0] = 120;
  trueData[3][3][1] = 91;
  trueData[3][3][2] = 3165;
  trueData[4][3][0] = 150;
  trueData[4][3][1] = 121;
  trueData[4][3][2] = 4065;
  trueData[5][3][0] = 180;
  trueData[5][3][1] = 151;
  trueData[5][3][2] = 4965;
  trueData[6][3][0] = 210;
  trueData[6][3][1] = 181;
  trueData[6][3][2] = 5865;

  for(int i=0; i<size->size(); i++) {
    Tucker::MetricData* mets = Tucker::computeSliceMetrics(t,i,Tucker::MIN + Tucker::MAX + Tucker::SUM);
    for(int j=0; j<(*size)[i]; j++) {
      std::cout << "The maximum of slice " << j << " of mode "
          << i << " is " << mets->getMaxData()[j] << std::endl;
      std::cout << "The minimum of slice " << j << " of mode "
          << i << " is " << mets->getMinData()[j] << std::endl;
      std::cout << "The sum of slice " << j << " of mode "
          << i << " is " << mets->getSumData()[j] << std::endl;

      if(mets->getMaxData()[j] != trueData[j][i][0]) {
        std::cout << mets->getMaxData()[j] << " != " << trueData[j][i][0]
                  << "; the difference is " << mets->getMaxData()[j] - trueData[j][i][0]
                  << std::endl;
        return EXIT_FAILURE;
      }
      if(mets->getMinData()[j] != trueData[j][i][1]) {
        std::cout << mets->getMinData()[j] << " != " << trueData[j][i][1]
                  << "; the difference is " << mets->getMinData()[j] - trueData[j][i][1]
                  << std::endl;
        return EXIT_FAILURE;
      }
      if(mets->getSumData()[j] != trueData[j][i][2]) {
        std::cout << mets->getSumData()[j] << " != " << trueData[j][i][2]
                  << "; the difference is " << mets->getSumData()[j] - trueData[j][i][2]
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    Tucker::MemoryManager::safe_delete<Tucker::MetricData>(mets);
  }

  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(size);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(t);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

