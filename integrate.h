#pragma once
#include <iostream>

double integrate_1_dim(double (*f)(double), double start, double end, int n);
double integrate_multi_dim(double (*f)(double[]), double *start, double *end, int dim, int n,
                           double *farg=nullptr, double *w=nullptr, double *x=nullptr);

void get_weights_points(double *weights, double *points, int n, double const a = -1., double const b = 1.); 