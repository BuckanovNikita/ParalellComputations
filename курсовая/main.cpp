#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>

#include "task.h"
#include "scheme.h"
#include "test.h"
#include "solver.h"

#define PRECISION 1e-8

using namespace std;

int MyNetInit(int *argc, char ***argv, int *np, int *mp,
              char *processor_name, double *tick) {

    int i = MPI_Init(argc, argv);
    int n1;
    if (i != 0) {
        cout << "MPI initialization error" << endl;
        exit(i);
    }

    MPI_Comm_size(MPI_COMM_WORLD, mp);
    MPI_Comm_rank(MPI_COMM_WORLD, np);
    MPI_Get_processor_name(processor_name, &n1);

    *tick = MPI_Wtime();
    sleep(1);
    return 0;
}

int main(int argc, char **argv) {

    double tick;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int np_, mp_;

    MyNetInit(&argc, &argv, &np_, &mp_, processor_name, &tick);

    auto np = (size_t) np_, mp = (size_t) mp_;

    cout << processor_name << endl;
    if (argc < 3) {
        cout << "Using N, max_iteration, num_threads" << endl;
        return -1;
    }

    char *err_marker = nullptr;
    size_t N = strtoul(argv[1], &err_marker, 10);
    if (*err_marker != '\0') {
        cout << "Invalid N parameter" << endl;
        return -1;
    }

    size_t max_iter = strtoul(argv[2], &err_marker, 10);
    if (*err_marker != '\0') {
        cout << "Invalid max_iteration parameter" << endl;
        return -1;
    }

    size_t mt = strtoul(argv[3], &err_marker, 10);
    if (*err_marker != '\0') {
        cout << "Invalid Thread number parameter" << endl;
        return -1;
    }

    const double h = 1.0 / N;
	double d1 = 16/1.01/(h*h)*sin(M_PI*h/2)*sin(M_PI*h/2);
	double D1 = 1600/(h*h)*cos(M_PI*h/2)*cos(M_PI*h/2);
	double d2 = 4/(h*h)*sin(M_PI*h/2)*sin(M_PI*h/2);
	double D2 = 4/(h*h)*cos(M_PI*h/2)*cos(M_PI*h/2);
    const double tau = 2/sqrt(min(d1,d2)*max(D1,D2));

    cout << "Grid size: " << N << endl;

    vector<vector<double>> U, U_1;
    for (size_t i = 0; i <= N; i++) {
        U.emplace_back(vector<double>(N + 1));
        U_1.emplace_back(vector<double>(N + 1));
    }

    for (size_t i = 0; i <= N; i++) {
        U[i][0] = U_1[i][0] = u_x_1_0(h * i);
        U[i][N] = U_1[i][N] = u_x_1_1(h * i);
    }

    for (size_t it = 0; it < max_iter; it++) {
        if (np == 0) {
            cout << "Iteration: " << it << endl;
            cout << "Max error: " << PrecisionInfo(U, U_1, h, N) << endl;
            cout << "Work time: " << MPI_Wtime() - tick << endl;
            cout << endl;
        }
#ifdef TEST
        BorderEqualityTest(U, U_1, N, h);
#endif
        for (size_t j = 1; j < N; j++) {
            vector<double> tmp;
            if (mp == 1 && mt == 1)
                tmp = SequentialThomasSolver(N,
                                             [=](size_t i) { return A_1(i, N, h, tau); },
                                             [=](size_t i) { return B_1(i, N, h, tau); },
                                             [=](size_t i) { return C_1(i, N, h, tau); },
                                             [&U, j, N, h, tau](size_t i) { return F_1(U, i, j, N, h, tau); });
            else if (mt == 1)
                tmp = MPISolver(N,
                                [=](size_t i) { return A_1(i, N, h, tau); },
                                [=](size_t i) { return B_1(i, N, h, tau); },
                                [=](size_t i) { return C_1(i, N, h, tau); },
                                [&U, j, N, h, tau](size_t i) { return F_1(U, i, j, N, h, tau); },
                                np, mp);
            else
                tmp = MPI_OMP_Solver(N,
                                     [=](size_t i) { return A_1(i, N, h, tau); },
                                     [=](size_t i) { return B_1(i, N, h, tau); },
                                     [=](size_t i) { return C_1(i, N, h, tau); },
                                     [&U, j, N, h, tau](size_t i) { return F_1(U, i, j, N, h, tau); },
                                     np, mp, mt);

            for (size_t i = 0; i <= N; i++) {
                U_1[i][j] = tmp[i];
            }
#ifdef TEST
            BorderEqualityTest(U, U_1, N, h);
#endif
        }

        swap(U, U_1);

        for (size_t j = 1; j < N; j++) {

            vector<double> tmp;
            if (mp == 1 && mt == 1)
                tmp = SequentialThomasSolver(N,
                                             [=](size_t i) { return A_2(i, N, h, tau); },
                                             [=](size_t i) { return B_2(i, N, h, tau); },
                                             [=](size_t i) { return C_2(i, N, h, tau); },
                                             [&U, j, N, h, tau](size_t i) {
                                                 return F_2(U, j, i, N, h, tau);
                                             });
            else if (mt == 1)
                tmp = MPISolver(N,
                                [=](size_t i) { return A_2(i, N, h, tau); },
                                [=](size_t i) { return B_2(i, N, h, tau); },
                                [=](size_t i) { return C_2(i, N, h, tau); },
                                [&U, j, N, h, tau](size_t i) { return F_2(U, j, i, N, h, tau); },
                                np, mp);
            else
                tmp = MPI_OMP_Solver(N,
                                     [=](size_t i) { return A_2(i, N, h, tau); },
                                     [=](size_t i) { return B_2(i, N, h, tau); },
                                     [=](size_t i) { return C_2(i, N, h, tau); },
                                     [&U, j, N, h, tau](size_t i) { return F_2(U, j, i, N, h, tau); },
                                     np, mp, mt);
            for (size_t i = 0; i <= N; i++)
                U_1[j][i] = tmp[i];
#ifdef TEST
            BorderEqualityTest(U, U_1, N, h);
#endif
        }
        swap(U, U_1);

    }
    if (np == 0) {

        cout << "Max error: " << PrecisionInfo(U, U_1, h, N) << endl;
        cout << "Work time: " << MPI_Wtime() - tick << endl;
        cout << endl;
    }
    MPI_Finalize();
    return 0;
}