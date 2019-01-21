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

void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double h);

double PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double h, size_t N);

void TopTriangleCheck(const size_t &l, const size_t &r, const size_t &np,
                      const vector<double> &L_V,
                      const vector<double> &B_V,
                      const vector<double> &C_V,
                      const vector<double> &F_V,
                      const vector<double> &correct);

void MatrixCheck(const size_t &l, const size_t &r, const size_t &np, const size_t &mp, const size_t &N,
                 const vector<double> &L_V,
                 const vector<double> &B_V,
                 const vector<double> &C_V,
                 const vector<double> &F_V,
                 const vector<double> &R_V,
                 const vector<double> &correct);

#endif //IMM_CPP_TEST_H
