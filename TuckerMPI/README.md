[![build status](https://gitlab.com/tensors/TuckerMPI/badges/master/build.svg)](https://gitlab.com/tensors/TuckerMPI/commits/master)
[![codecov](https://codecov.io/gl/tensors/tuckermpi/branch/master/graph/badge.svg?token=FO5DT6r7wc)](https://codecov.io/gl/tensors/tuckermpi)
<a href="https://scan.coverity.com/projects/tuckermpi">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/10762/badge.svg"/>
</a>

This is the GIT repo for the work on building a parallel Tucker for combustion data.                                                   

For more information:  
[Grey Ballard](mailto:ballard@wfu.edu)  
[Tamara Kolda](mailto:tgkolda@sandia.gov)  
[Hemanth Kolla](mailto:hnkolla@sandia.gov)  

We would live to also acknowledge important contributions to this code by the following persons:
* Woody Austin - Author of the original code, built using an older tensor codebase from Sandia
* Alicia Klinvex - Author of second version of code, which was rewritten from scratch on a new codebase
* Casey Battaglino - Important bug fix for computing TTMs in the case where the result is bigger than the input tensor

WARNING
-------
This code is still in development, but we welcome evaluation by friendly expert users.  Please contact us if you have any questions, or submit an issue if you find a bug or wish to request a new feature.

Requirements
------------
MPI implementation (We use openMPI MPICH2, and MVAPICH2)  
BLAS implementation  
LAPACK implementation  
C++11 or greater  

Documentation
-------------
Please see https://tensors.gitlab.io/TuckerMPI

Papers
------
Parallel Tensor Compression for Large-Scale Scientific Data  
Woody Austin, Grey Ballard, and Tamara G. Kolda  
IPDPS'16, [doi:10.1109/IPDPS.2016.67](https://doi.org/10.1109/IPDPS.2016.67)

[Find it on arXiv!](https://arxiv.org/abs/1510.06689)  
[Cite this!](latex/citations.txt)

Funding Statement
-----------------
The development of this software was supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Applied Mathematics program and a Sandia Truman Postdoctoral Fellowship (LDRD funding). Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation, a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE–AC04–94AL85000.