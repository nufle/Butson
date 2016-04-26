cmake_minimum_required (VERSION 2.6)

set(CMAKE_BUILD_TYPE Release)
project (BUTSON)
set(SOURCES main.cpp matrix.h ObjectPool.h)
add_executable(butson ${SOURCES})

include_directories("${PROJECT_BINARY_DIR}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -std=c++11")

