//
// Created by nikita on 16.01.2019.
//

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

vector<double> PseudoParallelThomasSolver(const size_t &N,
                                          const function<double(size_t)> &A,
                                          const function<double(size_t)> &B,
                                          const function<double(size_t)> &C,
                                          const function<double(size_t)> &F, const size_t &mp,
										  const vector<double> B_V_,
										  const vector<double> C_V_,
										  const vector<double> F_V_,
										  const vector<double> L_V_,
										  const vector<double> R_V_,
										  const vector<double> BUFFER_A_,
										  const vector<double> BUFFER_B_,
										  const vector<double> BUFFER_C_,
										  const vector<double> BUFFER_F_, size_t np_);
										  

						 
						 
vector<double> MPI_OMP_Solver(const size_t &N,
						 const function<double(size_t)> &A,
                         const function<double(size_t)> &B,
                         const function<double(size_t)> &C,
                         const function<double(size_t)> &F,
						 const size_t &np,
						 const size_t &mp,
						 const size_t &mt);
						 
vector<double> MPISolver(const size_t &N,
						 const function<double(size_t)> &A,
                         const function<double(size_t)> &B,
                         const function<double(size_t)> &C,
                         const function<double(size_t)> &F,
						 const size_t &np,
						 const size_t &mp);

#endif //IMM_CPP_SOLVER_H