//
// Created by nikita on 16.01.2019.
//


#ifndef IMM_CPP_TEST_H
#define IMM_CPP_TEST_H

#include "test.h"
#include "task.h"
#include "scheme.h"
#include <cassert>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <functional>

using namespace std;

void ThomasSolutionTest(vector<double> &U,
                        const size_t &N,
                        const function<double(size_t)> &A,
                        const function<double(size_t)> &B,
                        const function<double(size_t)> &C,
                        const function<double(size_t)> &F);

void DiagonalDominance(const size_t &N,
                      const function<double(size_t)> &A,
                      const function<double(size_t)> &B,
                      const function<double(size_t)> &C,
                      const function<double(size_t)> &F);

void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h);

void PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double a, double h, size_t N);

#endif //IMM_CPP_TEST_H
