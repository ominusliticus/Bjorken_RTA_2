
#ifndef DERIVATIVE_H_
#define DERIVATIVE_H_


void compute_first_derivative(double *fprime, double *f, int Nx, double dx);
void compute_second_derivative(double *f2prime, double *f, int Nx, double dx);
void compute_third_derivative(double *f3prime, double *f, int Nx, double dx);
void compute_fourth_derivative(double *f4prime, double *f, int Nx, double dx);

#endif