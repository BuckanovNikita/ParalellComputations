#include <iostream>
#include <vector>
#include <fstream>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "task.h"
#include "scheme.h"
#include "test.h"

#define PRECISION 1e-8
#define INFO
#define TEST

using namespace std;
int main(int argc, char **argv) {

    if(argc < 3) {
        cout << "Using N, tau, max_iteration filename" <<endl;
        return -1;
    }

    char* err_marker = nullptr;
    size_t N = strtoul(argv[1], &err_marker, 10);
    if(*err_marker != '\0') {
        cout << "Invalid N parameter" << endl;
        return -1;
    }

    const double tau = strtod(argv[2], &err_marker);
    if(*err_marker != '\0') {
        cout << "Invalid tau parameter" << endl;
        return -1;
    }

    size_t max_iter = strtoul(argv[3], &err_marker, 10);
    if(*err_marker != '\0') {
        cout << "Invalid max_iteration parameter" <<endl;
        return -1;
    }

    auto &logfile = (ofstream&)cout;
    size_t it = 0;
    const double a = 0;
    const double h = 1.0 / N;

    if (argc > 4) {
        logfile = ofstream();
        logfile.open(argv[4],ios::in);
    }

    cout << "Grid size: " << N << endl;

    vector<vector<double>> U, U_1;
    for (size_t i = 0; i <= N; i++) {
        U.emplace_back(vector<double>(N));
        U_1.emplace_back(vector<double>(N));
    }

    for(size_t i = 0; i<=N; i++) {
        U[i][0] = U_1[i][0] = u_x_1_0(a + h*i);
        U[i][N] = U_1[i][N] = u_x_1_1(a + h*i);
    }

    vector<double> s(N+1);
    vector<double> t(N+1);

    for (it = 0; it < max_iter; it++) {

#ifdef TEST
        DiagonalMajority(N, a, h, tau, (char*)"Matrix is diagonal majority\n", logfile);
        BorderEqualityTest(U, U_1, N, a, h, (char*)"Border is correct\n", logfile);
#endif

        for (size_t j = 1; j < N; j++) {

            s[0] = C_(0, N, a, h, tau) / B_(0, N, a, h, tau);
            t[0] = -F_(U, 0, j, N, a, h, tau) / B_(0, N, a, h, tau);
            for (size_t i = 1; i <= N; i++) {
                s[i] = C_(i, N, a, h, tau) / (B_(i, N, a, h, tau) - A_(i, N, a, h, tau) * s[i - 1]);
                t[i] = (A_(i, N, a, h, tau) * t[i - 1] - F_(U, i, j, N, a, h, tau)) /
                       (B_(i, N, a, h, tau) - A_(i, N, a, h, tau) * s[i - 1]);
            }
            U_1[N][j] = t[N];
            for (size_t i = N - 1; i > 0; i--) {
                U_1[i][j] = U_1[i + 1][j] * s[i] + t[i];
            }
            U_1[0][j] = U_1[1][j] * s[0] + t[0];
        }

#ifdef TEST
        SolutionTest_I(U, U_1, N, a, h, tau, (char*)"Progonka by columns is correct", logfile);
#endif
        swap(U, U_1);

#ifdef TEST
        DiagonalMajority(N, a, h, tau, (char*)"Matrix is diagonal majority", logfile);
        BorderEqualityTest(U, U_1, N, a, h,(char*)"Border is correct", logfile);
#endif

#ifdef INFO
        PrecisionInfo(U, U_1, a, h, it, N, (char*)"\nIteration: I", (char*)"max error \nmax_i  max_j ", logfile);
#endif

        for (size_t i = 1; i < N; i++) {
            s[0] = C1(0, N, a, h, tau) / B1(0, N, a, h, tau);
            t[0] = -F1(U, i, 0, N, a, h, tau) / B1(0, N, a, h, tau);
            for (size_t j = 1; j <= N; j++) {
                s[j] = C1(j, N, a, h, tau) / (B1(j, N, a, h, tau) - A1(j, N, a, h, tau) * s[j - 1]);
                t[j] = (A1(j, N, a, h, tau) * t[j - 1] - F1(U, i, j, N, a, h, tau)) /
                       (B1(j, N, a, h, tau) - A1(j, N, a, h, tau) * s[j - 1]);
            }

            U_1[i][N] = t[N];
            for (size_t j = N - 1; j > 0; j--) {
                U_1[i][j] = U_1[i][j + 1] * s[j] + t[j];
            }
            U_1[i][0] =  U_1[i][1] * s[0] + t[0];
        }
#ifdef TEST
        SolutionTest_J(U, U_1, N, a, h, tau, (char*)"Progonka by rows is correct", logfile);
#endif

        swap(U, U_1);

#ifdef INFO
        PrecisionInfo(U, U_1, a, h, it, N, (char*)"Iteration: J", (char*)"max error \nmax_i %d max_j %d", logfile);
#endif
    }
    return 0;
}
