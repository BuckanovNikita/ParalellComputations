//
// Created by nikita on 16.01.2019.
//

#ifndef IMM_CPP_SCHEME_H
#define IMM_CPP_SCHEME_H

#include "task.h"
#include <vector>

using namespace std;

double A_1(const size_t &i, const size_t &N, const double &h, const double &tau);

double B_1(const size_t &i, const size_t &N, const double &h, const double &tau);

double C_1(const size_t &i, const size_t &N, const double &h, const double &tau);

double F_1(const vector<vector<double>> &U, const size_t &i, const size_t &j, const size_t &N,
           const double &h,
           const double &tau);

double A_2(const size_t &i, const size_t &N, const double &h, const double &tau);

double B_2(const size_t &j, const size_t &N, const double &h, const double &tau);

double C_2(const size_t &j, const size_t &N, const double &h, const double &tau);

double F_2(const vector<vector<double>> &U, const size_t &j, const size_t &i, const size_t &N,
           const double &h, const double &tau);

#endif //IMM_CPP_SCHEME_H
