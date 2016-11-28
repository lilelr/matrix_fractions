cmake_minimum_required(VERSION 3.6)
project(Matrix_Fractorization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(Matrix_Fractorization ${SOURCE_FILES})