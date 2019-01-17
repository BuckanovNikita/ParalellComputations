#include <iostream>
#include <vector>
#include <fstream>


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

    vector<double> B_V[mp], C_V[mp], F_V[mp], L[mp], R[mp];

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
            //TODO Параллельная прогонка решение в переменной tmp
            auto A = [=](size_t i) { return A_1(i, N, a, h, tau); };
            auto B = [=](size_t i) { return -B_1(i, N, a, h, tau); };
            auto B_ = [=](size_t i) { return B_1(i, N, a, h, tau); };
            auto C = [=](size_t i) { return C_1(i, N, a, h, tau); };
            auto F = [&U, j, N, a, h, tau](size_t i) { return F_1(U, i, j, N, a, h, tau); };

#ifdef TEST
            vector<double> correct = SequentialThomasSolver(N,A,B_,C,F);
            ThomasSolutionTest(correct, N, A, B_, C, F);
            cout<<"Sequential solution: OK"<<endl;
#endif

#ifdef MATRIX_INFO
            for (size_t i = 0; i <= N; i++)
                cout << A(i) << " " << B(i) << " " << C(i) << " " << F(i) << endl;
#endif

            for (size_t np = 0; np < mp; np++) {
                size_t l = np * (N + 1) / mp;
                size_t r = (np + 1) * (N + 1) / mp - 1;
                if (np == mp) {
                    r = N + 1;
                }
                B_V[np] = vector<double>(r - l + 1);
                C_V[np] = vector<double>(r - l + 1);
                F_V[np] = vector<double>(r - l + 1);
                L[np] = vector<double>(r - l + 1);
                R[np] = vector<double>(r - l + 2);
            }

            for (size_t np = 0; np < mp; np++) {
                size_t l = np * (N + 1) / mp;
                size_t r = (np + 1) * (N + 1) / mp - 1;
                if (np == mp) {
                    r = N + 1;
                }
                cout << "Block: " << l << " " << r << endl;
                for (size_t k = 0; k <= r - l; k++) {
                    B_V[np][k] = B(l + k);
                    C_V[np][k] = C(l + k);
                    F_V[np][k] = F(l + k);
                }

                L[np][0] = A(l);
                for (size_t k = 1; k <= r - l; k++) {
                    B_V[np][k] -= C_V[np][k - 1] * A(l + k) / B_V[np][k - 1];
                    F_V[np][k] -= F_V[np][k - 1] * A(l + k) / B_V[np][k - 1];
                    if (np != 0) {
                        L[np][k] = -L[np][k - 1] * A(l + k) / B_V[np][k - 1];
                     }
                }
#ifdef TEST
                if( np == 0 )
                    for(size_t i = 0; i <= r-l; i++)
                    {
                        double res = fabs(B_V[np][i]*correct[i] + C_V[np][i]*correct[i+1] - F_V[np][i]);
                        assert(res < PRECISION);
                    }
                else {

                    for (size_t i = 0; i <= r - l; i++) {
                        double res = fabs(L[np][i] * correct[l - 1] + B_V[np][i] * correct[l + i] +
                                          C_V[np][i] * correct[l + i + 1] - F_V[np][i]);
                        assert(res < PRECISION);
                    }
                }
                cout<<"Upper triangle form: OK"<<endl;
#endif

                R[np][r-l+1] = B_V[np][r-l];
                R[np][r-l] = C_V[np][r-l-1];
                for(size_t k= r-l-1; k>0; k--)
                {
                    R[np][k] = -C_V[np][k-1]/B_V[np][k]*C_V[np][k];
                    L[np][k-1] = L[np][k-1]-C_V[np][k-1]/B_V[np][k]*L[np][k];
                    F_V[np][k-1] =  F_V[np][k-1]-C_V[np][k-1]/B_V[np][k]*F_V[np][k];
                }

                //TODO Обмен коэффиецентами

                if(np!=0)
                {
                    B_V[np-1].back() -= C(l-1)/B(l)*L[np][0];
                    R[np][0] = -C(l-1)/B(l)*R[np][1];
                }
            }

            for (size_t np = 0; np < mp; np++) {
                size_t l = np * (N + 1) / mp;
                size_t r = (np + 1) * (N + 1) / mp - 1;
                if (np == mp) {
                    r = N + 1;
                }
                cout << "Block: " << l << " " << r << endl;
                if(np!=0)
                    BUFFER_A[np] = L[np][r-l];

                BUFFER_B[np] = B_V[np][r-l];
                BUFFER_F[np] = F_V[np][r-l];

                if(np!=mp-1)
                    BUFFER_C[np] = R[np+1][0];

            }
            vector<double> tmp = SequentialThomasSolver(mp,
                                                        [=](size_t i) {return BUFFER_A[i];},
                                                        [=](size_t i) {return BUFFER_B[i];},
                                                        [=](size_t i) {return BUFFER_C[i];},
                                                        [=](size_t i) {return BUFFER_F[i];});
#ifdef TEST
            ThomasSolutionTest(tmp, mp,[=](size_t i) {return BUFFER_A[i];},
                               [=](size_t i) {return BUFFER_B[i];},
                               [=](size_t i) {return BUFFER_C[i];},
                               [=](size_t i) {return BUFFER_F[i];});
            cout<<"Solution of small system: OK"<<endl;
#endif
            for (size_t np = 0; np < mp; np++) {
                size_t l = np * (N + 1) / mp;
                size_t r = (np + 1) * (N + 1) / mp - 1;
                if (np == mp) {
                    r = N + 1;
                }
                cout << "Block: " << l << " " << r << endl;
                BUFFER_U[r] = tmp[np];
                for(ssize_t i=r-l-1; i>=0; i--) {
                    BUFFER_U[l+i] = F_V[np][i];
                    if(np!=0)
                        BUFFER_U[l+i] -= L[np][i]*tmp[np-1];
                    if(np!=mp-1)
                        BUFFER_U[l+i] -= R[np][i+1]*tmp[np];
                    BUFFER_U[l+i] /= B_V[np][i];
                }
            }
#ifdef TEST
            ThomasSolutionTest(BUFFER_U, N,
                               [=](size_t i) { return A_1(i, N, a, h, tau); },
                               [=](size_t i) { return B_1(i, N, a, h, tau); },
                               [=](size_t i) { return C_1(i, N, a, h, tau); },
                               [&U, j, N, a, h, tau](size_t i) { return F_1(U, i, j, N, a, h, tau); });
#endif
            }

        }
        cout<<BUFFER_U[0]<<endl;



        /*swap(U, U_1);

#ifdef TEST
        BorderEqualityTest(U, U_1, N, a, h);
#endif

#ifdef INFO
        PrecisionInfo(U, U_1, a, h, N);
#endif

        for (size_t j = 1; j < N; j++) {
#ifdef TEST
            DiagonalDominance(N,
                              [=](size_t i) { return A1(i, N, a, h, tau); },
                              [=](size_t i) { return B1(i, N, a, h, tau); },
                              [=](size_t i) { return C1(i, N, a, h, tau); },
                              [&U, j, N, a, h, tau](size_t i) {
                                  return F1(U, j, i, N, a, h, tau);
                              });
#endif
            //TODO Параллельная прогонка решение в переменной tmp
#ifdef TEST
            ThomasSolutionTest(tmp, N,
                               [=](size_t i) { return A1(i, N, a, h, tau); },
                               [=](size_t i) { return B1(i, N, a, h, tau); },
                               [=](size_t i) { return C1(i, N, a, h, tau); },
                               [&U, j, N, a, h, tau](size_t i) { return F1(U, j, i, N, a, h, tau); });
#endif
        }
        swap(U, U_1);
#ifdef INFO
        PrecisionInfo(U, U_1, a, h, N);
#endif*/
    return 0;
}