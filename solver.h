//
// Created by nikita on 16.01.2019.
//

#ifndef IMM_CPP_SOLVER_H
#define IMM_CPP_SOLVER_H

#include <vector>
#include <functional>
using namespace std;

vector<double> SequentialThomasSolver(const size_t &N,
        const function<double(size_t)> &A,
        const function<double(size_t)> &B,
        const function<double(size_t)> &C,
        const function<double(size_t)> &F);

#endif //IMM_CPP_SOLVER_H
