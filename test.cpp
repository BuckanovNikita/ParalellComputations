//
// Created by nikita on 15.01.2019.
//

#include <fstream>
#include "test.h"
#define PRECISION 1e-8

void DiagonalMajority(size_t N, double a, double h, double tau)
{
    size_t i;
    for(i=0; i<=N; i++)
        assert(fabs(A_(i, N, a, h, tau) + C_(i, N, a, h, tau)) <= fabs(B_(i, N, a, h, tau)));
}

void SolutionTest_I(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h, double tau)
{
    size_t j, i;
    for(j=1; j<N; j++){
        assert(fabs(U_1[0][j])<PRECISION);
        assert(fabs(U_1[N][j])<PRECISION);
        for(i=1; i<N; i++)
        {
            if((fabs(A_(i, N, a, h, tau) * U_1[i - 1][j] - B_(i, N, a, h, tau) * U_1[i][j] +
                     C_(i, N, a, h, tau) * U_1[i + 1][j] -
                     F_(U, i, j, N, a, h, tau)) > PRECISION))
                fprintf(stdout,"%lf \n", fabs(
                        A_(i, N, a, h, tau) * U_1[i - 1][j] - B_(i, N, a, h, tau) * U_1[i][j] +
                        C_(i, N, a, h, tau) * U_1[i + 1][j] -
                        F_(U, i, j, N, a, h, tau)));
            assert(fabs(
                    A_(i, N, a, h, tau) * U_1[i - 1][j] - B_(i, N, a, h, tau) * U_1[i][j] +
                    C_(i, N, a, h, tau) * U_1[i + 1][j] -
                    F_(U, i, j, N, a, h, tau)) < PRECISION);
        }
    }
}

void SolutionTest_J(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h, double tau)
{
    size_t i, j;
    for (i=1; i<N;i++)
    {
        assert(fabs(U_1[i][0] - u_x_1_0(a + h * i)) < PRECISION);
        assert(fabs(U_1[i][N] - u_x_1_1(a + h * i)) < PRECISION);
        for (j = 1; j < N; j++) {
            if (fabs(A1(j, N, a, h, tau) * U_1[i][j - 1] - B1(j, N, a, h, tau) * U_1[i][j] +
                     C1(j, N, a, h, tau) * U_1[i][j + 1] - F1(U, i, j, N, a, h, tau)) > PRECISION)
                fprintf(stdout, "%lf \n", fabs(A1(j, N, a, h, tau) * U_1[i][j - 1] - B1(j, N, a, h, tau) * U_1[i][j] +
                                               C1(j, N, a, h, tau) * U_1[i][j + 1] - F1(U, i, j, N, a, h, tau)));
            assert(fabs(A1(j, N, a, h, tau) * U_1[i][j - 1] - B1(j, N, a, h, tau) * U_1[i][j] +
                        C1(j, N, a, h, tau) * U_1[i][j + 1] - F1(U, i, j, N, a, h, tau)) < PRECISION);
        }

    }
}

void BorderEqualityTest(vector<vector<double>> &U, vector<vector<double>> &U_1, size_t N, double a, double h) {
    size_t i;
    for(i=0; i<=N; i++)
    {
        assert(U[i][0] == u_x_1_0(a + h*i));
        assert(U_1[i][0] == u_x_1_0(a + h*i));
        assert(U[i][N] == u_x_1_1(a + h*i));
        assert(U_1[i][N] == u_x_1_1(a + h*i));
    }
}

void PrecisionInfo(vector<vector<double>> &U, vector<vector<double>> &U_1, double a, double h, size_t it, size_t N)
{
    int i = 0;
    int j = 0;
    int max_i = 0;
    int max_j = 0;
    double diff = -2;

    for(i=0; i<N; i++)
    {
        for(j=0;j<N; j++)
        {
            if (fabs(U[i][j] - f(a + h * i, a + h * j)) > diff)
            {
                diff = fabs(U[i][j] - solution(a + h * i, a + h * j));
                max_i = i;
                max_j = j;
            }
        }
    }
}
