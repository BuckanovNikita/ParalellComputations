#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

#include "task.h"
#include "scheme.h"
#include "test.h"
#include "solver.h"
#include "mpi.h"

#define PRECISION 1e-8
#define INFO

#define TEST
//#define SEQUENTIAL

#ifndef SEQUENTIAL
#ifndef  PSEUDO_PARALLEL
#ifndef MPI_PARALLEL
#define MPI_OMP_PARALLEL
#endif
#endif
#endif

using namespace std;

int MyNetInit(int* argc, char*** argv, int* np, int* mp,
        char* processor_name, double* tick)
{

    int i = MPI_Init(argc, argv);
    int n1;
    if (i != 0){
        cerr << "MPI initialization error"<<endl;
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

    auto np = (size_t)np_, mp = (size_t)mp_;

    cout<<processor_name<<endl;
    if (argc < 4) {
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
	
	size_t mt = strtoul(argv[4], &err_marker, 10);
	if (*err_marker != '\0') {
        cout << "Invalid Thread number parameter" << endl;
        return -1;
    }

    size_t it = 0;
    const double a = 0;
    const double h = 1.0 / N;


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
        BorderEqualityTest(U, U_1, N, a, h);
#endif
        for (size_t j = 1; j < N; j++) {

#ifdef SEQUENTIAL
            vector<double> tmp = SequentialThomasSolver(N,
                                                        [=](size_t i) { return A_1(i, N, a, h, tau); },
                                                        [=](size_t i) { return B_1(i, N, a, h, tau); },
                                                        [=](size_t i) { return C_1(i, N, a, h, tau); },
                                                        [&U, j, N, a, h, tau](size_t i) {
                                                            return F_1(U, i, j, N, a, h, tau);
                                                        });
#endif

#ifdef PSEUDO_PARALLEL
            vector<double> tmp = PseudoParallelThomasSolver(N,
                                                            [=](size_t i) { return A_1(i, N, a, h, tau); },
                                                            [=](size_t i) { return B_1(i, N, a, h, tau); },
                                                            [=](size_t i) { return C_1(i, N, a, h, tau); },
                                                            [&U, j, N, a, h, tau](size_t i) {
                                                                return F_1(U, i, j, N, a, h, tau);
                                                            });
#endif

#ifdef MPI_PARALLEL
            vector<double> tmp = MPISolver(N,
                                        [=](size_t i) { return A_1(i, N, a, h, tau); },
                                        [=](size_t i) { return B_1(i, N, a, h, tau); },
                                        [=](size_t i) { return C_1(i, N, a, h, tau); },
                                        [&U, j, N, a, h, tau](size_t i) {return F_1(U, i, j, N, a, h, tau);},
										np, mp);
#endif

#ifdef MPI_OMP_PARALLEL
            vector<double> tmp = MPI_OMP_Solver(N,
                                        [=](size_t i) { return A_1(i, N, a, h, tau); },
                                        [=](size_t i) { return B_1(i, N, a, h, tau); },
                                        [=](size_t i) { return C_1(i, N, a, h, tau); },
                                        [&U, j, N, a, h, tau](size_t i) {return F_1(U, i, j, N, a, h, tau);},
										np, mp, mt);
#endif

            for (int i = 0; i <= N; i++) {
                U_1[i][j] = tmp[i];
            }
#ifdef TEST
            BorderEqualityTest(U, U_1, N, a, h);
#endif
        }

        swap(U, U_1);

#ifdef INFO
if(np == 0)
	PrecisionInfo(U, U_1, a, h, N);
#endif

        for (size_t j = 1; j < N; j++) {

#ifdef SEQUENTIAL
            vector<double> tmp = SequentialThomasSolver(N,
                                                        [=](size_t i) { return A_2(i, N, a, h, tau); },
                                                        [=](size_t i) { return B_2(i, N, a, h, tau); },
                                                        [=](size_t i) { return C_2(i, N, a, h, tau); },
                                                        [&U, j, N, a, h, tau](size_t i) {
                                                            return F_2(U, j, i, N, a, h, tau);
                                                        });
#endif

#ifdef PSEUDO_PARALLEL
            vector<double> tmp = PseudoParallelThomasSolver(N,
                                                            [=](size_t i) { return A_2(i, N, a, h, tau); },
                                                            [=](size_t i) { return B_2(i, N, a, h, tau); },
                                                            [=](size_t i) { return C_2(i, N, a, h, tau); },
                                                            [&U, j, N, a, h, tau](size_t i) {
                                                                return F_2(U, j, i, N, a, h, tau);
                                                            });
#endif

#ifdef MPI_PARALLEL
            vector<double> tmp = MPISolver(N,
                                        [=](size_t i) { return A_2(i, N, a, h, tau); },
                                        [=](size_t i) { return B_2(i, N, a, h, tau); },
                                        [=](size_t i) { return C_2(i, N, a, h, tau); },
                                        [&U, j, N, a, h, tau](size_t i) {return F_2(U, j, i, N, a, h, tau);},
										np, mp, mt);
#endif

#ifdef MPI_OMP_PARALLEL
            vector<double> tmp = MPI_OMP_Solver(N,
                                        [=](size_t i) { return A_2(i, N, a, h, tau); },
                                        [=](size_t i) { return B_2(i, N, a, h, tau); },
                                        [=](size_t i) { return C_2(i, N, a, h, tau); },
                                        [&U, j, N, a, h, tau](size_t i) {return F_2(U, j, i, N, a, h, tau);},
										np, mp, mt);
#endif
            for (int i = 0; i <= N; i++)
                U_1[j][i] = tmp[i];
#ifdef TEST
            BorderEqualityTest(U, U_1, N, a, h);
#endif
        }
        swap(U, U_1);
#ifdef INFO
if(np == 0)
{
	PrecisionInfo(U, U_1, a, h, N);
	cout << endl;
}
#endif
    }
	MPI_Finalize();
    return 0;
}