//
// Created by nikita on 16.01.2019.
//

#include "solver.h"

#define TEST
#define PRECISION 1e-8

vector<double> SequentialThomasSolver(const size_t &N,
                                      const function<double(size_t)> &A,
                                      const function<double(size_t)> &B,
                                      const function<double(size_t)> &C,
                                      const function<double(size_t)> &F) {
#ifdef TEST
    DiagonalDominance(N, A, B, C, F);
#endif

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
#ifdef TEST
    ThomasSolutionTest(result, N, A, B, C, F);
#endif
    return result;
}

vector<double> PseudoParallelThomasSolver(const size_t &N,
                                          const function<double(size_t)> &A,
                                          const function<double(size_t)> &B,
                                          const function<double(size_t)> &C,
                                          const function<double(size_t)> &F, const size_t &mp) {
#ifdef TEST
    vector<double> correct = SequentialThomasSolver(N, A, B, C, F);
    ThomasSolutionTest(correct, N, A, B, C, F);
#endif
    vector<double> BUFFER_A(mp), BUFFER_B(mp), BUFFER_C(mp), BUFFER_F(mp);
    vector<double> ANSWER(N + 1);
    vector<double> B_V(N + 1), C_V(N + 1), F_V(N + 1), L_V(N + 1), R_V(N + 1);

    for (size_t i = 0; i < N + 1; i++) {
        B_V[i] = -B(i);
        C_V[i] = C(i);
        F_V[i] = F(i);
    }

    for (size_t np = 0; np < mp; np++) {
        size_t l = np * (N + 1) / mp;
        size_t r = (np + 1) * (N + 1) / mp - 1;
        if (np == mp) {
            r = N + 1;
        }
        L_V[l] = A(l);
        for (size_t i = l + 1; i <= r; i++) {
            double tmp = A(i) / B_V[i - 1];
            B_V[i] -= tmp * C_V[i - 1];
            F_V[i] -= tmp * F_V[i - 1];
            if (np != 0 && i != l)
                L_V[i] -= tmp * L_V[i - 1];
        }
#ifdef TEST
        TopTriagleCheck(l, r, np, L_V, B_V, C_V, F_V, correct);
#endif
        R_V[r - 1] = C_V[r - 1];

        for (size_t i = r - 2; i >= l; i--) {
            if (i == 0 && np == 0)
                break;
            double tmp = C_V[i] / B_V[i + 1];
            C_V[i] -= tmp * B_V[i + 1];
            R_V[i] -= tmp * R_V[i + 1];
            L_V[i] -= tmp * L_V[i + 1];
            F_V[i] -= tmp * F_V[i + 1];
        }
        if (np != 0) {
            size_t i = l - 1;
            double tmp = C_V[i] / B_V[i + 1];
            C_V[i] -= tmp * B_V[i + 1];
            R_V[i] -= tmp * R_V[i + 1];
            B_V[i] -= tmp * L_V[i + 1];
            F_V[i] -= tmp * F_V[i + 1];
        }
    }

    for (size_t np = 0; np < mp; np++) {
        size_t l = np * (N + 1) / mp;
        size_t r = (np + 1) * (N + 1) / mp - 1;
        if (np == mp)
            r = N + 1;
        BUFFER_A[np] = L_V[r];
        BUFFER_B[np] = B_V[r];
        BUFFER_C[np] = R_V[r];
        BUFFER_F[np] = F_V[r];

#ifdef TEST
        MatrixCheck(l, r, np, mp, N, L_V, B_V, C_V, F_V, R_V, correct);
#endif
    }

    vector<double> small_solution = SequentialThomasSolver(mp - 1,
                                                           [=](size_t i) { return BUFFER_A[i]; },
                                                           [=](size_t i) { return -BUFFER_B[i]; },
                                                           [=](size_t i) { return BUFFER_C[i]; },
                                                           [=](size_t i) { return BUFFER_F[i]; });

    for (size_t np = 0; np < mp; np++) {
        size_t l = np * (N + 1) / mp;
        size_t r = (np + 1) * (N + 1) / mp - 1;
        if (np == mp)
            r = N + 1;
        ANSWER[r] = small_solution[np];
        for (ssize_t i = r - 1; i >= (ssize_t) l; i--) {
            if (np == 0) {
                ANSWER[i] = (F_V[i] - ANSWER[r] * R_V[i]) / B_V[i];
            } else {
                ANSWER[i] = (F_V[i] - ANSWER[l - 1] * L_V[i] - ANSWER[r] * R_V[i]) / B_V[i];
            }
        }
    }
#ifdef TEST
    {
        double max_d = -1;
        for (size_t i = 0; i <= N; i++) {
            max_d = max(max_d, abs(ANSWER[i] - correct[i]));
            assert(abs(ANSWER[i] - correct[i]) < PRECISION);
        }
        ThomasSolutionTest(ANSWER, N, A, B, C, F);
    };
#endif
    return ANSWER;

};