set(CMAKE_MODULE_PATH "/nfs/sw/modules/cmake/3.26.0")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(cryptoTools_DIR "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/cryptoTools/cmake")
set(NTL_DIR "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/NTL")
set(NTL_INCLUDE_DIR "${NTL_DIR}/include")
set(NTL_LIBRARY_DIR "${NTL_DIR}/lib")
set(NTL_LIBRARY "${NTL_LIBRARY_DIR}/libntl.a")
set(GMP_DIR "/nfs/sw/gmp/gmp-6.1.2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
cmake_minimum_required(VERSION 3.18)

project(securesort)
set(CMAKE_CXX_COMPILER "/nfs/sw/gcc/gcc-11.2.0/bin/g++")

add_executable(securesort sharing.cpp mpc.cpp utils.cpp)

find_package(cryptoTools REQUIRED)
find_package(GMP REQUIRED PATHS "/nfs/sw/gmp/gmp-6.1.2")
link_directories(${NTL_LIBRARY_DIR})
message(${NTL_LIBRARY_DIR})


target_include_directories(securesort PRIVATE ${GMP_DIR}/include ${NTL_INCLUDE_DIR})
target_link_libraries(securesort oc::cryptoTools ${GMP_LIBRARIES} ${NTL_LIBRARY})

target_compile_features(securesort PRIVATE cxx_std_14)

# find_package(NTL REQUIRED PATHS "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/modules")
# find_library(NTL_LIBRARY NAMES libntl.a PATHS ${NTL_DIR}/lib)
# find_library(GMP_LIBRARY NAMES libgmp.a PATHS ${GMP_DIR}/lib)
# find_library(GMP_LIBRARY NAMES libgmp.a PATHS ${GMP_DIR}/gmp)
# target_link_libraries(crypto ${NTL_LIBRARY})
# add_executable(mpc sharing.cpp)
# target_link_libraries(mpc oc::cryptoTools)
