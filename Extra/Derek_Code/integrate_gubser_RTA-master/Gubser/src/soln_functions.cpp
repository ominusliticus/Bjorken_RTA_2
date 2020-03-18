#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "coord_functions.cpp"
#include "Parameter.h"

float f_eq(float p_omega_hat, float p_eta_hat, float rho, float T_hat)
{
  float p_rho_hat = sqrt( (p_omega_hat*p_omega_hat / (cosh(rho)*cosh(rho)) ) + p_eta_hat*p_eta_hat );
  return exp(-p_rho_hat / T_hat);
}

float DampingFunction(float rho2, float rho1, gsl_spline *T_spline, gsl_interp_accel *acc, parameters &params)
{
  //define an integration step size
  float drho = params.delta_rho / params.resolution;
  float c = params.c;
  float result = 0.;
  float rho = rho1;
  while (rho < rho2)
  {
    result += gsl_spline_eval(T_spline, rho, acc);
    rho += drho;
  }
  result *= drho;
  return exp(-1.0 * result / c);
}

//NOTE THIS SOLUTION ASSUMES THE INITIAL DISTRIBUTION IS EQUIL. f0 = feq(tau_0)
float f_solution(float rho, float p_omega_hat, float p_eta_hat, gsl_spline *T_spline, gsl_interp_accel *acc, parameters &params)
{
  float drho = params.delta_rho / params.resolution;
  float rho0 = params.rho0;
  float T_hat0 = params.T_hat0;

  //the first term
  float I1 = DampingFunction(rho, rho0, T_spline, acc, params) * f_eq(p_omega_hat, p_eta_hat, rho0, T_hat0);
  //the second term
  float I2 = 0.;
  float rho_p = rho0;
  while (rho_p < rho)
  {
    float T_hat_rho_p = gsl_spline_eval(T_spline, rho_p, acc);
    I2 += DampingFunction(rho, rho_p, T_spline, acc, params) * T_hat_rho_p * f_eq(p_omega_hat, p_eta_hat, rho_p, T_hat_rho_p);
    rho_p += drho;
  }
  I2 *= drho;
  return I1 + I2;
}
