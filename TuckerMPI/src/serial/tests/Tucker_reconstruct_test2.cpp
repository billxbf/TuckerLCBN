/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>

bool approxEqual(double a, double b, double tol)
{
  if(std::abs(a-b) > tol)
    return false;
  return true;
}

int main()
{
  Tucker::SizeArray* size =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(3);
  (*size)[0] = 2;
  (*size)[1] = 3;
  (*size)[2] = 5;
  Tucker::Tensor* t =
      Tucker::MemoryManager::safe_new<Tucker::Tensor>(*size);
  double* data = t->data();
  for(int i=0; i<30; i++)
    data[i] = i+1;

  Tucker::SizeArray* newSize =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(3);
  (*newSize)[0] = 2;
  (*newSize)[1] = 1;
  (*newSize)[2] = 3;
  const struct Tucker::TuckerTensor* factorization = Tucker::STHOSVD(t,newSize);

  // Reconstruct the original tensor
  Tucker::Tensor* temp = Tucker::ttm(factorization->G,0,
      factorization->U[0]);
  Tucker::Tensor* temp2 = Tucker::ttm(temp,1,
        factorization->U[1]);
  Tucker::Tensor* temp3 = Tucker::ttm(temp2,2,
        factorization->U[2]);

  data = temp3->data();

  double trueData[30];
  trueData[0] = 2.802710268427637; trueData[1] = 3.697422679857209; trueData[2] = 3.112029908952288;
  trueData[3] = 4.105486783765727; trueData[4] = 3.421349549476938; trueData[5] = 4.513550887674245;
  trueData[6] = 8.170984737005284; trueData[7] = 9.065697148434888; trueData[8] = 9.072771157833156;
  trueData[9] = 10.066228032646631; trueData[10] = 9.974557578661029; trueData[11] = 11.066758916858374;
  trueData[12] = 13.539259205582903; trueData[13] = 14.433971617012585; trueData[14] = 15.033512406713992;
  trueData[15] = 16.026969281527553; trueData[16] = 16.527765607845087; trueData[17] = 17.619966946042524;
  trueData[18] = 18.907533674160593; trueData[19] = 19.802246085590220; trueData[20] = 20.994253655594914;
  trueData[21] = 21.987710530408407; trueData[22] = 23.080973637029231; trueData[23] = 24.173174975226598;
  trueData[24] = 24.275808142738295; trueData[25] = 25.170520554167854; trueData[26] = 26.954994904475836;
  trueData[27] = 27.948451779289265; trueData[28] = 29.634181666213383; trueData[29] = 30.726383004410678;

  for(int i=0; i<30; i++) {
    if(!approxEqual(data[i],trueData[i],1e-10)) {
      std::cerr << "data[" << i << "] (" << data[i] << ") != trueData["
          << i << "] (" << trueData[i] << "); the difference is "
          << data[i] - trueData[i] << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Free some memory
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(t);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(size);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(newSize);
  Tucker::MemoryManager::safe_delete<const struct Tucker::TuckerTensor>(factorization);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(temp);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(temp2);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(temp3);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
