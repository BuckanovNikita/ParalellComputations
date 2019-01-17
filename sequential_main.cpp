#include <iostream>
#include <vector>
#include <fstream>


#include "task.h"
#include "scheme.h"
#include "test.h"
#include "solver.h"

#define PRECISION 1e-8
#define INFO
#define TEST

using namespace std;

int main(int argc, char **argv) {

    if (argc < 3) {
        cout << "Using N, tau, max_iteration filename" << endl;
        return -1;
    }

    char *err_marker = nullptr;
    size_t N = strtoul(argv[1], &err_marker, 10);
    if (*err_marker != '\0') {
        cout << "Invalid N parameter" << endl;
        return -1;
    }

    const double tau = strtod(argv[2], &err_marker);
    if (*err_marker != '\0') {
        cout << "Invalid tau parameter" << endl;
        return -1;
    }

    size_t max_iter = strtoul(argv[3], &err_marker, 10);
    if (*err_marker != '\0') {
        cout << "Invalid max_iteration parameter" << endl;
        return -1;
    }

    size_t it = 0;
    const double a = 0;
    const double h = 1.0 / N;


    cout << "Grid size: " << N << endl;

    vector<vector<double>> U, U_1;
    for (size_t i = 0; i <= N; i++) {
        U.emplace_back(vector<double>(N));
        U_1.emplace_back(vector<double>(N));
    }

    for (size_t i = 0; i <= N; i++) {
        U[i][0] = U_1[i][0] = u_x_1_0(a + h * i);
        U[i][N] = U_1[i][N] = u_x_1_1(a + h * i);
    }

    for (it = 0; it < max_iter; it++) {

#ifdef TEST
        BorderEqualityTest(U, U_1, N, a, h);
#endif
        for (size_t j = 1; j < N; j++) {
#ifdef TEST
            DiagonalDominance(N,
                    [=](size_t i) { return A_1(i, N, a, h, tau); },
                    [=](size_t i) { return B_1(i, N, a, h, tau); },
                    [=](size_t i) { return C_1(i, N, a, h, tau); },
                    [&U, j, N, a, h, tau](size_t i) {
                        return F_1(U, i, j, N, a, h, tau);
                    });
#endif
            vector<double> tmp = SequentialThomasSolver(N,
                                                        [=](size_t i) { return A_1(i, N, a, h, tau); },
                                                        [=](size_t i) { return B_1(i, N, a, h, tau); },
                                                        [=](size_t i) { return C_1(i, N, a, h, tau); },
                                                        [&U, j, N, a, h, tau](size_t i) {
                                                            return F_1(U, i, j, N, a, h, tau);
                                                        });

            for (int i = 0; i <= N; i++) {
                U_1[i][j] = tmp[i];
            }
#ifdef TEST
            ThomasSolutionTest(tmp, N,
                               [=](size_t i) { return A_1(i, N, a, h, tau); },
                               [=](size_t i) { return B_1(i, N, a, h, tau); },
                               [=](size_t i) { return C_1(i, N, a, h, tau); },
                               [&U, j, N, a, h, tau](size_t i) { return F_1(U, i, j, N, a, h, tau); });
#endif
        }

        swap(U, U_1);

#ifdef TEST
        BorderEqualityTest(U, U_1, N, a, h);
#endif

#ifdef INFO
        PrecisionInfo(U, U_1, a, h, N);
#endif

        for (size_t j = 1; j < N; j++) {
#ifdef TEST
            DiagonalDominance(N,
                             [=](size_t i) { return A_2(i, N, a, h, tau); },
                             [=](size_t i) { return B_2(i, N, a, h, tau); },
                             [=](size_t i) { return C_2(i, N, a, h, tau); },
                             [&U, j, N, a, h, tau](size_t i) {
                                 return F_2(U, j, i, N, a, h, tau);
                             });
#endif
            vector<double> tmp = SequentialThomasSolver(N,
                                                        [=](size_t i) { return A_2(i, N, a, h, tau); },
                                                        [=](size_t i) { return B_2(i, N, a, h, tau); },
                                                        [=](size_t i) { return C_2(i, N, a, h, tau); },
                                                        [&U, j, N, a, h, tau](size_t i) {
                                                            return F_2(U, j, i, N, a, h, tau);
                                                        });

            for (int i = 0; i <= N; i++)
                U_1[j][i] = tmp[i];
#ifdef TEST
            ThomasSolutionTest(tmp, N,
                               [=](size_t i) { return A_2(i, N, a, h, tau); },
                               [=](size_t i) { return B_2(i, N, a, h, tau); },
                               [=](size_t i) { return C_2(i, N, a, h, tau); },
                               [&U, j, N, a, h, tau](size_t i) { return F_2(U, j, i, N, a, h, tau); });
#endif
        }
        swap(U, U_1);
#ifdef INFO
        PrecisionInfo(U, U_1, a, h, N);
#endif
    }
    return 0;
}