//
// Created by nikita on 16.01.2019.
//

#include "solver.h"

vector<double> SequentialThomasSolver(const size_t &N,
                                      const function<double(size_t)> &A,
                                      const function<double(size_t)> &B,
                                      const function<double(size_t)> &C,
                                      const function<double(size_t)> &F) {
    vector<double> s(N + 1);
    vector<double> t(N + 1);
    vector<double> result(N + 1);
    s[0] = C(0) / B(0);
    t[0] = -F(0) / B(0);
    for (size_t i = 1; i <= N; i++) {
        s[i] = C(i) / (B(i) - A(i) * s[i - 1]);
        t[i] = (A(i) * t[i - 1] - F(i)) / (B(i) - A(i) * s[i - 1]);
    }
    result[N] = t[N];
    for (size_t i = N - 1; i > 0; i--) {
        result[i] = result[i + 1] * s[i] + t[i];
    }
    result[0] = result[1] * s[0] + t[0];
    return result;
}
