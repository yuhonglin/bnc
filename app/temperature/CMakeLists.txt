cmake_minimum_required(VERSION 3.8)
set (CMAKE_CXX_STANDARD 11)
add_executable(sample ./temperature.cpp ../../src/dist/brmath/build/libbrmath.a ../../src/rng/BasicRNG.cc ../../src/util/logger.cc)
target_include_directories(sample PUBLIC ../../src/)
