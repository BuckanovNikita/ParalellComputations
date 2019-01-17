//
// Created by nikita on 15.01.2019.
//

#include "scheme.h"

double A_1(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau) {
    if (i == 0 || i == N)
        return 0;

    return k_1(a + h * (i - 0.5)) * tau / (2 * h * h);
}

double B_1(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau) {
    if (i == 0 || i == N)
        return -1;

    return 1 + (k_1(a + h * (i + 0.5)) + k_1(a + h * (i - 0.5))) * tau / (2 * h * h);
}

double C_1(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau) {
    if (i == 0 || i == N)
        return 0;

    return k_1(a + h * (i + 0.5)) * tau / (2 * h * h);
}

double
F_1(const vector<vector<double>> &U, const size_t &i, const size_t &j, const size_t &N, const double &a,
    const double &h,
    const double &tau) {
    if (i == 0 || i == N)
        return 0;

    if (j == 0)
        return u_x_1_0(a + h * i);

    if (j == N)
        return u_x_1_1(a + h * i);

    return -U[i][j] - tau * (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (2 * h * h) - tau / 2 * f(a + h * i, a + h * j);
}

double A_2(const size_t &i, const size_t &N, const double &a, const double &h, const double &tau) {
    if (i == 0 || i == N)
        return 0;

    return tau / (2 * h * h);
}

double B_2(const size_t &j, const size_t &N, const double &a, const double &h, const double &tau) {
    if (j == 0 || j == N)
        return -1;

    return 1 + tau / (h * h);
}

double C_2(const size_t &j, const size_t &N, const double &a, const double &h, const double &tau) {
    if (j == 0 || j == N)
        return 0;

    return tau / (2 * h * h);
}

double F_2(const vector<vector<double>> &U, const size_t &j, const size_t &i, const size_t &N, const double &a,
           const double &h, const double &tau) {
    if (i == 0)
        return u_x_1_0(a + h * j);

    if (i == N)
        return u_x_1_1(a + h * j);

    return -U[j][i] - (k_1(a + h * (j + 0.5)) * (U[j + 1][i] - U[j][i])
                       - k_1(a + h * (j - 0.5)) * (U[j][i] - U[j - 1][i])) * tau / (2 * h * h) -
           tau * f(a + h * j, a + h * i) / 2;
}