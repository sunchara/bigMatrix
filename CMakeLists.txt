cmake_minimum_required(VERSION 2.8) 
project(BigMatrix)
set(CMAKE_CXX_STANDARD 20)
set(TEST_SOURSES tests.cpp BigMatrix.hpp)

add_executable(tests ${TEST_SOURSES})