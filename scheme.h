//
// Created by nikita on 16.01.2019.
//

#ifndef IMM_CPP_SCHEME_H
#define IMM_CPP_SCHEME_H

#include "task.h"
#include <stddef.h>
#include <vector>

using namespace std;

double A_(size_t i, size_t N, double a, double h, double tau);

double B_(size_t i, size_t N, double a, double h, double tau);

double C_(size_t i, size_t N, double a, double h, double tau);

double F_(vector<vector<double>> &U, size_t i, size_t j, size_t N, double a, double h, double tau);

double A1(size_t i, size_t N, double a, double h, double tau);

double B1(size_t j, size_t N, double a, double h, double tau);

double C1(size_t j, size_t N, double a, double h, double tau);

double F1(vector<vector<double>> &U, size_t i, size_t j, size_t N, double a, double h, double tau);

#endif //IMM_CPP_SCHEME_H
