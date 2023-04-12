set(CMAKE_MODULE_PATH "/nfs/sw/modules/cmake/3.26.0")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(cryptoTools_DIR "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/cryptoTools/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
cmake_minimum_required(VERSION 3.18)

project(securesort)
set(CMAKE_CXX_COMPILER "/nfs/sw/gcc/gcc-11.2.0/bin/g++")

add_executable(securesort sharing.cpp)

find_package(cryptoTools REQUIRED)

target_link_libraries(securesort oc::cryptoTools -lstdc++)
# target_link_libraries(securesort -lcryptotools -lgmp -lgmpxx)
target_compile_features(securesort PRIVATE cxx_std_14)
# add_executable(mpc sharing.cpp)
# target_link_libraries(mpc oc::cryptoTools)
# target_link_libraries(securesort mpc.cpp)