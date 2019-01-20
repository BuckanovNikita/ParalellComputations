#ifndef IMM_CPP_SOLVER_H
#define IMM_CPP_SOLVER_H

#include <vector>
#include <functional>
#include <iostream>
#include "test.h"
#include "mpi.h"
#include "omp.h"

using namespace std;

vector<double> SequentialThomasSolver(const size_t &N,
                                      const function<double(size_t)> &A,
                                      const function<double(size_t)> &B,
                                      const function<double(size_t)> &C,
                                      const function<double(size_t)> &F);

vector<double>
PseudoParallelThomasSolver(const size_t &N, const function<double(size_t)> &A, const function<double(size_t)> &B,
                           const function<double(size_t)> &C, const function<double(size_t)> &F, const size_t &mp,
                           vector<double> B_V_, vector<double> F_V_, vector<double> L_V_,
                           vector<double> R_V_, vector<double> BUFFER_A_, vector<double> BUFFER_B_,
                           vector<double> BUFFER_C_, vector<double> BUFFER_F_, size_t np_);
										  
vector<double> MPISolver(const size_t &N, 
						 const function<double(size_t)> &A,
                         const function<double(size_t)> &B,
                         const function<double(size_t)> &C,
                         const function<double(size_t)> &F,
						 const size_t &np,
						 const size_t &mp);
						 
						 
vector<double> MPI_OMP_Solver(const size_t &N,
						 const function<double(size_t)> &A,
                         const function<double(size_t)> &B,
                         const function<double(size_t)> &C,
                         const function<double(size_t)> &F,
						 const size_t &np,
						 const size_t &mp,
						 const size_t &mt);

#endif //IMM_CPP_SOLVER_H