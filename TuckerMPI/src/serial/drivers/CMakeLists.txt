INCLUDE_DIRECTORIES(..)

ADD_EXECUTABLE(Tucker_reconstruct Tucker_reconstruct.cpp)
TARGET_LINK_LIBRARIES(Tucker_reconstruct serial_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(Tucker_general_reconstruct Tucker_general_reconstruct.cpp)
TARGET_LINK_LIBRARIES(Tucker_general_reconstruct serial_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(Tucker_sthosvd Tucker_sthosvd.cpp)
TARGET_LINK_LIBRARIES(Tucker_sthosvd serial_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(LCBN_sthosvd LCBN_sthosvd.cpp)
TARGET_LINK_LIBRARIES(LCBN_sthosvd serial_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(Tucker_generate Tucker_generate.cpp)
TARGET_LINK_LIBRARIES(Tucker_generate serial_tucker ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
