#include <fstream>
#include "test.h"

#define PRECISION 1e-8

void DiagonalDominance(const size_t &N,
                       const function<double(size_t)> &A,
                       const function<double(size_t)> &B,
                       const function<double(size_t)> &C,
                       const function<double(size_t)> &F) {
    for(size_t i = 0; i <= N; i++)
        assert(fabs(A(i) + C(i)) <= fabs(B(i)));
    cout<<"Diagonal dominance: OK"<<endl;
}

void ThomasSolutionTest(vector<double> &U,
                        const size_t &N,
                        const function<double(size_t)> &A,
                        const function<double(size_t)> &B,
                        const function<double(size_t)> &C,
                        const function<double(size_t)> &F) {
    assert(fabs(-B(0) * U[0] + C(0) * U[1] - F(0)) < PRECISION);

    for (size_t i = 1; i < N; i++)
    {
        double res = fabs(A(i) * U[i - 1] - B(i) * U[i] + C(i) * U[i + 1] - F(i));
        assert(fabs(A(i) * U[i - 1] - B(i) * U[i] + C(i) * U[i + 1] - F(i)) < PRECISION);
    }
    assert(fabs(A(N) * U[N - 1] - B(N) * U[N] + -F(N)) < PRECISION);
}

void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h) {
    for (size_t i = 0; i <= N; i++) {
        assert(U[i][0] == u_x_1_0(a + h * i));
        assert(U_1[i][0] == u_x_1_0(a + h * i));
        assert(U[i][N] == u_x_1_1(a + h * i));
        assert(U_1[i][N] == u_x_1_1(a + h * i));
    }
}

void PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double a, double h, size_t N) {
    size_t max_i = 0;
    size_t max_j = 0;
    double diff = -2;

    for (size_t i = 0; i <= N; i++)
        for (size_t j = 0; j <= N; j++)
            if (fabs(U[i][j] - solution(a + h * i, a + h * j)) > diff) {
                diff = fabs(U[i][j] - solution(a + h * i, a + h * j));
                max_i = i;
                max_j = j;
            }
    cout << "Max error: " << diff << endl;
    cout << "I:" << max_i << " J:" << max_j << endl;
}
