INCLUDE_DIRECTORIES(.. ../../serial)

ADD_EXECUTABLE(reconstruct TuckerMPI_reconstruct.cpp)
TARGET_LINK_LIBRARIES(reconstruct mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(general_reconstruct TuckerMPI_general_reconstruct.cpp)
TARGET_LINK_LIBRARIES(general_reconstruct mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(sthosvd TuckerMPI_sthosvd.cpp)
TARGET_LINK_LIBRARIES(sthosvd mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(scaling TuckerMPI_scaling.cpp)
TARGET_LINK_LIBRARIES(scaling mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(generate TuckerMPI_generate.cpp)
TARGET_LINK_LIBRARIES(generate mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(gram TuckerMPI_gram.cpp)
TARGET_LINK_LIBRARIES(gram mpi_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})