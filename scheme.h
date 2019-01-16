//
// Created by nikita on 16.01.2019.
//

#ifndef IMM_CPP_SCHEME_H
#define IMM_CPP_SCHEME_H

#include "task.h"
#include <stddef.h>
#include <vector>

using namespace std;

double A_(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau);

double B_(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau);

double C_(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau);

double F_(const vector<vector<double>> &U, const size_t &i, const size_t &j, const size_t &N, const double &a, const double &h,
          const double &tau);

double A1(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau);

double B1(const size_t &j, const size_t &N, const double &a, const double &h, const double &tau);

double C1(const size_t &j, const size_t &N, const double &a, const double &h, const double &tau);

double F1(vector<vector<double>> &U, const size_t &i, const size_t &j, const size_t &N, const double &a,
          const double &h, const double &tau);

#endif //IMM_CPP_SCHEME_H
