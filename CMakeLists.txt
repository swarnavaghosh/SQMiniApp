cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(miniapp LANGUAGES CXX CUDA)

find_package(MPI REQUIRED) 
add_library(lsdft STATIC
  src/lsdft_initObjs.cc
  src/lsdft_matvec.cu
  src/lsdft_scf.cc
  src/lsdft_readfiles.cc
  src/lsdft_tools.cc
  src/lsdft_nl.cc
  )


include_directories(SYSTEM ${MPI_INCLUDE_PATH}) 
include_directories(include ${PETSC_DIR}/lib/petsc/conf/variables) 
include_directories(include ${PETSC_DIR}/lib/petsc/conf/rules) 
include_directories(include ${PETSC_DIR}/lib/petsc/conf/test) 

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to 
# particles will also build with -std=c++11
target_compile_features(lsdft PUBLIC cxx_std_11)
 
# We need to explicitly state that we need all CUDA files in the 
# particle library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( lsdft
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
 
target_link_libraries(lsdft petsc mpi_ibm)
add_executable(lsdft_test src/main.cc)
 
set_property(TARGET lsdft_test 
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(lsdft_test PRIVATE lsdft)
