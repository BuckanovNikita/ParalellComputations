cmake_minimum_required(VERSION 3.5.2)
project(imm_cpp)

set(CMAKE_CXX_STANDARD 11)

#add_executable(imm_cpp_sequential sequential_main.cpp test.cpp
#        test.h task.cpp task.h scheme.cpp scheme.h solver.cpp solver.h)

add_executable(mpi_parallel main.cpp test.cpp
        test.h task.cpp task.h scheme.cpp scheme.h solver.cpp solver.h)
		
add_executable(mpi_solver parallel_main.cpp test.cpp
        test.h task.cpp task.h scheme.cpp scheme.h solver.cpp solver.h)
		
add_executable(openmp_solver parallel_main.cpp test.cpp
        test.h task.cpp task.h scheme.cpp scheme.h solver.cpp solver.h)

#add_executable(imm_cpp pseudo_parallel_thomas.cpp test.cpp
#        test.h task.cpp task.h scheme.cpp scheme.h solver.cpp solver.h)