cmake_minimum_required(VERSION 3.5.2)
project(imm_cpp)

set(CMAKE_CXX_STANDARD 11)

add_executable(mpi_parallel main.cpp test.cpp
        test.h task.cpp task.h scheme.cpp scheme.h solver.h solver.cpp)