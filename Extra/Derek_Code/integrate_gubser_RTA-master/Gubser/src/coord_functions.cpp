#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

float rho_func(float tau, float r, float q)
{
    float num = 1. - (q*q) * (tau*tau) + (q*q) * (r*r);
    float den = 2. * q * tau + 1e-10;
    return -asinh(num / den);
}

float theta_func(float tau, float r, float q)
{
  float num = 2. * q * r;
  float den = 1. + (q*q) * (tau*tau) - (q*q) * (r*r) + 1e-10;
  return atan(num / den);
}

float kappa_func(float tau, float r, float q)
{
  float num = 2. * q*q * tau * r;
  float den = 1. + q*q * ( tau*tau + r*r) + 1e-10;
  return atanh(num / den);
}

float p_theta_hat_func(float p_tau, float p_r, float q)
{
  float num = 2. * q * p_r;
  float den = 1. + (q*q) * (p_tau*p_tau) - (q*q) * (p_r*p_r) + 1e-10;
  return atan(num / den);
}

float tau_func(float rho, float theta, float q)
{
  float num = 1. / cosh(rho);
  float den = cos(theta) - tanh(rho) + 1e-10;
  return num / den / q;
}

float r_func(float rho, float theta, float q)
{
  float num = sin(theta);
  float den = cos(theta) - tanh(rho) + 1e-10;
  return num / den / q;
}

float p_omega_hat_func(float p_theta_hat, float p_phi_hat, float theta)
{
  float c = (p_phi_hat / (sin(theta) + 1e-10) );
  return p_theta_hat*p_theta_hat + c*c;
}

float p_theta_func(float rho, float theta, float q, float p_tau, float p_r)
{
  float del_tau_del_theta_num = - sin(theta) / cosh(rho);
  float c = (cos(theta) - tanh(rho));
  float den = c*c + 1e-10;
  float del_tau_del_theta = del_tau_del_theta_num / den / q;

  float del_r_del_theta_num = cos(theta)*(cos(theta) - tanh(rho)) + sin(theta)*sin(theta);
  float del_r_del_theta = del_r_del_theta_num / den / q;

  return del_tau_del_theta * p_tau + del_r_del_theta * p_r;
}
