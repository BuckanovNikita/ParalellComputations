#include <fstream>
#include "test.h"

#define PRECISION 1e-8

void DiagonalDominance(const size_t &N,
                       const function<double(size_t)> &A,
                       const function<double(size_t)> &B,
                       const function<double(size_t)> &C,
                       const function<double(size_t)> &F) {
    for (size_t i = 0; i <= N; i++)
        assert(fabs(A(i) + C(i)) <= fabs(B(i)));
}

void ThomasSolutionTest(vector<double> &U,
                        const size_t &N,
                        const function<double(size_t)> &A,
                        const function<double(size_t)> &B,
                        const function<double(size_t)> &C,
                        const function<double(size_t)> &F) {
    assert(fabs(-B(0) * U[0] + C(0) * U[1] - F(0)) < PRECISION);

    for (size_t i = 1; i < N; i++) {
        assert(fabs(A(i) * U[i - 1] - B(i) * U[i] + C(i) * U[i + 1] - F(i)) < PRECISION);
    }
    assert(fabs(A(N) * U[N - 1] - B(N) * U[N] + -F(N)) < PRECISION);
}

void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double h) {
    for (size_t i = 0; i <= N; i++) {
        assert(U[i][0] == u_x_1_0(h * i));
        assert(U_1[i][0] == u_x_1_0(h * i));
        assert(U[i][N] == u_x_1_1(h * i));
        assert(U_1[i][N] == u_x_1_1(h * i));
    }
}

double PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double h, size_t N) {
    double diff = -2;

    for (size_t i = 0; i <= N; i++)
        for (size_t j = 0; j <= N; j++)
            if (fabs(U[i][j] - solution(h * i, h * j)) > diff) {
                diff = fabs(U[i][j] - solution(h * i, h * j));
            }
	return diff;
}

void TopTriangleCheck(const size_t &l, const size_t &r, const size_t &np,
                      const vector<double> &L_V,
                      const vector<double> &B_V,
                      const vector<double> &C_V,
                      const vector<double> &F_V,
                      const vector<double> &correct) {
    if (np == 0)
        for (size_t i = l; i <= r; i++)
            assert(abs(B_V[i] * correct[i]
                       + C_V[i] * correct[i + 1] - F_V[i]) < PRECISION);
    else
        for (size_t i = l; i <= r; i++)
            assert(abs(L_V[i] * correct[l - 1] + B_V[i] * correct[i]
                       + C_V[i] * correct[i + 1] - F_V[i]) < PRECISION);
}

void MatrixCheck(const size_t &l, const size_t &r, const size_t &np, const size_t &mp, const size_t &N,
                 const vector<double> &L_V,
                 const vector<double> &B_V,
                 const vector<double> &C_V,
                 const vector<double> &F_V,
                 const vector<double> &R_V,
                 const vector<double> &correct) {
    if (np == 0) {
        for (size_t i = l; i < r; i++)
            assert(abs(B_V[i] * correct[i]
                       + R_V[i] * correct[r] - F_V[i]) < PRECISION);
        //assert(abs(B_V[r] * correct[r]
         //          + R_V[r] * correct[min(r + (N + 1) / mp, N)] - F_V[r]) < PRECISION);
    } else {
        for (size_t i = l; i < r; i++)
            assert(abs(L_V[i] * correct[l - 1] + B_V[i] * correct[i]
                       + R_V[i] * correct[r] - F_V[i]) < PRECISION);

        /*if (np != mp - 1) {
            size_t tmp_r = r + (N + 1) / mp;
            if(np == mp-2)
                tmp_r = N;
        assert(abs(L_V[r] * correct[l - 1] + B_V[r] * correct[r]
                   + R_V[r] * correct[tmp_r] - F_V[r]) < PRECISION);
        }
        else
            assert(abs(L_V[r] * correct[l - 1] + B_V[r] * correct[r] - F_V[r]) < PRECISION);*/
    }
}