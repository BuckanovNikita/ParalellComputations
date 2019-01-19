//
// Created by nikita on 15.01.2019.
//

#include "task.h"


double k_1(double x) {
    return 4 / (0.01 + sin(M_PI * x) * sin(M_PI * x));
}

double f(double x_1, double x_2) {
    return 8 * M_PI * M_PI * sin(2 * M_PI * x_1) * sinh(M_PI * (x_2 - 0.5)) *
           (-0.125 + 2 / (0.01 + sin(M_PI * x_1) * sin(M_PI * x_1)) +
            cos(2 * M_PI * x_1) /
            ((0.01 + sin(M_PI * x_1) * sin(M_PI * x_1)) * (0.01 + sin(M_PI * x_1) * sin(M_PI * x_1))));
}

double u_x_1_0(double x_1) {
    return sin(2 * M_PI * x_1) * sinh(-M_PI * 0.5);
}

double u_x_1_1(double x_1) {
    return sin(2 * M_PI * x_1) * sinh(M_PI * 0.5);
}

double solution(double x_1, double x_2) {
    return sin(2 * M_PI * x_1) * sinh(M_PI * (x_2 - 0.5));
}