#include <iostream>
#include <vector>
#include <unistd.h>

#include "task.h"
#include "scheme.h"
#include "test.h"
#include "solver.h"
#include "mpi.h"

#define PRECISION 1e-8
#define TEST
#define MY_TAG 777

using namespace std;

int MyNetInit(int *argc, char ***argv, int *np, int *mp,
              char *processor_name, double *tick) {

    int i = MPI_Init(argc, argv);
    int n1;
    if (i != 0) {
        cerr << "MPI initialization error" << endl;
        exit(i);
    }

    MPI_Comm_size(MPI_COMM_WORLD, mp);
    MPI_Comm_rank(MPI_COMM_WORLD, np);
    MPI_Get_processor_name(processor_name, &n1);

    *tick = MPI_Wtick();
    sleep(1);
    return 0;
}

int main(int argc, char **argv) {

    double tick;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int np_, mp_;

    MyNetInit(&argc, &argv, &np_, &mp_, processor_name, &tick);

    auto np = (size_t) np_, mp = (size_t) mp_;
    cout << "Process ID: " << np << " Max ID: " << mp << " Name: " << processor_name << endl;
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
    MPI_Status status;


    cout << "Grid size: " << N << endl;

    vector<vector<double>> U, U_1;
    for (size_t i = 0; i <= N; i++) {
        U.emplace_back(vector<double>(N + 1));
        U_1.emplace_back(vector<double>(N + 1));
    }

    for (size_t i = 0; i <= N; i++) {
        U[i][0] = U_1[i][0] = u_x_1_0(a + h * i);
        U[i][N] = U_1[i][N] = u_x_1_1(a + h * i);
    }

    for (it = 0; it < max_iter; it++) {
#ifdef TEST
        if (np == 0)
            BorderEqualityTest(U, U_1, N, a, h);
#endif
        for (size_t j = 1; j < N; j++) {
            const size_t l = np * (N + 1) / mp;
            const size_t r = np == mp - 1 ? N : (np + 1) * (N + 1) / mp - 1;

            auto A = [=](size_t i) { return A_1(i, N, a, h, tau); };
            auto B = [=](size_t i) { return B_1(i, N, a, h, tau); };
            auto C = [=](size_t i) { return C_1(i, N, a, h, tau); };
            auto F = [&U, j, N, a, h, tau](size_t i) { return F_1(U, i, j, N, a, h, tau); };
#ifdef TEST
            vector<double> correct;
            correct = SequentialThomasSolver(N, A, B, C, F);
#endif
            vector<double> BUFFER_A(mp), BUFFER_B(mp), BUFFER_C(mp), BUFFER_F(mp);
            vector<double> ANSWER(N + 1);
            vector<double> B_V(N + 1), C_V(N + 1), F_V(N + 1), L_V(N + 1), R_V(N + 1);

            for (size_t i = 0; i < N + 1; i++) {
                B_V[i] = -B(i);
                C_V[i] = C(i);
                F_V[i] = F(i);
            }

            L_V[l] = A(l);
            for (size_t i = l + 1; i <= r; i++) {
                double tmp_ = A(i) / B_V[i - 1];
                B_V[i] -= tmp_ * C_V[i - 1];
                F_V[i] -= tmp_ * F_V[i - 1];
                if (np != 0 && i != l)
                    L_V[i] -= tmp_ * L_V[i - 1];
            }
#ifdef TEST
            TopTriangleCheck(l, r, np, L_V, B_V, C_V, F_V, correct);
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
                MPI_Send(&L_V[l], 1, MPI_DOUBLE, (int) np - 1, MY_TAG, MPI_COMM_WORLD);
                MPI_Send(&B_V[l], 1, MPI_DOUBLE, (int) np - 1, MY_TAG, MPI_COMM_WORLD);
                MPI_Send(&R_V[l], 1, MPI_DOUBLE, (int) np - 1, MY_TAG, MPI_COMM_WORLD);
                MPI_Send(&F_V[l], 1, MPI_DOUBLE, (int) np - 1, MY_TAG, MPI_COMM_WORLD);
            }
            if (np != mp - 1) {
                MPI_Recv(L_V.data() + r + 1, 1, MPI_DOUBLE, (int) np + 1, MY_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(B_V.data() + r + 1, 1, MPI_DOUBLE, (int) np + 1, MY_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(C_V.data() + r + 1, 1, MPI_DOUBLE, (int) np + 1, MY_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(F_V.data() + r + 1, 1, MPI_DOUBLE, (int) np + 1, MY_TAG, MPI_COMM_WORLD, &status);
                double tmp = C_V[r] / B_V[r + 1];
                C_V[r] -= tmp * B_V[r + 1];
                R_V[r] -= tmp * R_V[r + 1];
                B_V[r] -= tmp * L_V[r + 1];
                F_V[r] -= tmp * F_V[r + 1];
            }
#ifdef TEST
            MatrixCheck(l, r, np, mp, N, L_V, B_V, C_V, F_V, R_V, correct);
#endif
            for (int i = 0; i < mp; i++) {
                if (i != np) {
                    MPI_Send(&L_V[r], 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD);
                    MPI_Send(&B_V[r], 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD);
                    MPI_Send(&R_V[r], 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD);
                    MPI_Send(&F_V[r], 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD);
                } else {
                    BUFFER_A[np] = L_V[r];
                    BUFFER_B[np] = B_V[r];
                    BUFFER_C[np] = R_V[r];
                    BUFFER_F[np] = F_V[r];
                }
            }

            for (int i = 0; i < mp; i++) {
                if (i != np) {
                    MPI_Recv(BUFFER_A.data() + i, 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD, &status);
                    MPI_Recv(BUFFER_B.data() + i, 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD, &status);
                    MPI_Recv(BUFFER_C.data() + i, 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD, &status);
                    MPI_Recv(BUFFER_F.data() + i, 1, MPI_DOUBLE, i, MY_TAG, MPI_COMM_WORLD, &status);
                }
            }
            vector<double> small_solution = SequentialThomasSolver(mp - 1,
                                                                   [=](size_t i) { return BUFFER_A[i]; },
                                                                   [=](size_t i) { return -BUFFER_B[i]; },
                                                                   [=](size_t i) { return BUFFER_C[i]; },
                                                                   [=](size_t i) { return BUFFER_F[i]; });

            for (int i = 0; i < mp; i++) {
                size_t r_ = i == mp - 1 ? N : (i + 1) * (N + 1) / mp - 1;
                ANSWER[r_] = small_solution[i];
            }
            for (ssize_t i = r - 1; i >= (ssize_t) l; i--) {
                if (np == 0) {
                    ANSWER[i] = (F_V[i] - ANSWER[r] * R_V[i]) / B_V[i];
                } else {
                    ANSWER[i] = (F_V[i] - ANSWER[l - 1] * L_V[i] - ANSWER[r] * R_V[i]) / B_V[i];
                }
            }
#ifdef TEST
            double max_d = -1;
            for (size_t i = l; i <= r; i++) {
                max_d = max(max_d, abs(ANSWER[i] - correct[i]));
                assert(abs(ANSWER[i] - correct[i]) < PRECISION);
            }
#endif
            for (int i = 0; i <= N; i++) {
                U_1[i][j] = ANSWER[i];
            }
        }
    }
    MPI_Finalize();
    return 0;
}