//
// Created by nikita on 16.01.2019.
//


#ifndef IMM_CPP_TEST_H
#define IMM_CPP_TEST_H

#include "test.h"
#include "task.h"
#include "scheme.h"
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

using namespace std;

void DiagonalMajority(size_t N, double a, double h, double tau, char *msg, ofstream &output);
void SolutionTest_I(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h, double tau,
                    char *msg, ofstream &output);
void SolutionTest_J(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h, double tau,
                    char *msg, ofstream &output);
void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h, char *msg,
                        ofstream &output);
void PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double a, double h, size_t it, size_t N,
                   char *fmt1, char *fmt2, ofstream &output);

#endif //IMM_CPP_TEST_H
