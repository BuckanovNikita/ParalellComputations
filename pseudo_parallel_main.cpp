#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>


#include "task.h"
#include "scheme.h"
#include "test.h"
#include "solver.h"

#define PARALLEL_INFO
//#define INFO
#define PRECISION 1e-8
//#define MATRIX_INFO
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
    size_t mp = 3;
    size_t it = 0;
    const double a = 0;
    const double h = 1.0 / N;


    cout << "Grid size: " << N << endl;

    vector<vector<double>> U, U_1;
    for (size_t i = 0; i <= N; i++) {
        U.emplace_back(vector<double>(N+1));
        U_1.emplace_back(vector<double>(N+1));
    }

    for (size_t i = 0; i <= N; i++) {
        U[i][0] = U_1[i][0] = u_x_1_0(a + h * i);
        U[i][N] = U_1[i][N] = u_x_1_1(a + h * i);
    }

    vector<double> BUFFER_A(mp), BUFFER_B(mp), BUFFER_C(mp), BUFFER_F(mp);

    vector<double> BUFFER_U(N+1);

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
            //Параллельная прогонка по строкам
            auto A = [=](size_t i) { return A_1(i, N, a, h, tau); };
            auto B = [=](size_t i) { return B_1(i, N, a, h, tau); };
            auto C = [=](size_t i) { return C_1(i, N, a, h, tau); };
            auto F = [&U, j, N, a, h, tau](size_t i) { return F_1(U, i, j, N, a, h, tau); };

#ifdef TEST
            vector<double> correct = SequentialThomasSolver(N,A,B,C,F);
            ThomasSolutionTest(correct, N, A, B, C, F);
            cout<<"Sequential solution: OK"<<endl;
#endif
            vector<double> A_V(N+1), B_V(N+1), C_V(N+1), F_V(N+1), L_V(N+1), R_V(N+1);
            for(size_t i=0; i < N+1; i++)
            {
                A_V[i] = A(i);
                B_V[i] = -B(i);
                C_V[i] = C(i);
                F_V[i] = F(i);
            }

            for(size_t np=0; np < mp; np++)
            {
                size_t l = np * (N + 1) / mp;
                size_t r = (np + 1) * (N + 1) / mp - 1;
                if (np == mp) {
                    r = N + 1;
                }
                cout << "Block: " << l << " " << r << endl;
                L_V[l] = A(l);
                for(size_t i=l+1; i<=r; i++)
                {
                    double tmp = A_V[i]/B_V[i-1];
                    A_V[i] -= tmp * B_V[i-1];
                    B_V[i] -= tmp * C_V[i-1];
                    F_V[i] -= tmp * F_V[i-1];
                    if(np!=0 && i!=l)
                        L_V[i] -= tmp * L_V[i-1];
                }
#ifdef TEST
                if(np == 0)
                    for(size_t i=l; i<=r; i++)
                        assert(abs(B_V[i]*correct[i]
                        + C_V[i]*correct[i+1] - F_V[i])<PRECISION);
                else
                    for(size_t i=l; i<=r; i++)
                        assert(abs(L_V[i] * correct[l-1] + B_V[i]*correct[i]
                        + C_V[i]*correct[i+1] - F_V[i])<PRECISION);
                 cout<<"Top triangle shape: OK"<<endl;
#endif
                 R_V[r] = B_V[r];
                 R_V[r-1]= C_V[r-1];

                 for(size_t i=r-2; i>=l; i--)
                 {
                     if(i == 0 && np == 0)
                         break;
                     double tmp = C_V[i]/B_V[i+1];
                     C_V[i] -= tmp * B_V[i+1];
                     R_V[i] -= tmp * R_V[i+1];
                     L_V[i] -= tmp * L_V[i+1];
                     F_V[i] -= tmp * F_V[i+1];
                 }
#ifdef TEST
                if(np == 0)
                    for(size_t i=l; i<r; i++)
                        assert(abs(B_V[i]*correct[i]
                                   + R_V[i]*correct[r] - F_V[i])<PRECISION);
                else
                    for(size_t i=l; i<r; i++)
                        assert(abs(L_V[i] * correct[l-1] + B_V[i]*correct[i]
                                   + R_V[i]*correct[r] - F_V[i])<PRECISION);
                cout<<"Parallel shape: OK"<<endl;
#endif
            }
        }
    }
    return 0;
}