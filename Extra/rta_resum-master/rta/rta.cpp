// RTA_BLBE - A solver for the Relaxation Time Approximation (RTA) in the Bjorken Limit (BL) (i.e. using the Bjorken symmetry) of the Boltzmann Equation (BE)
// Author: Gojko Vujanovic 
// Licence: GNU General Public Licence version 3, see details at http://www.gnu.org/licenses/gpl-3.0.html
// Copyright (c) 2018 
// This code was used to obtain the solution in Physcal Review C 97, no.6, 064909 (2018)
//
// Last edited: 1/8/20 Mike McNelis
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include "include/hypergeometric.h"
#include "include/derivative.h"
#include "include/anisobjorken.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
using namespace std;

// constants
const double prefactor = 3. / (M_PI * M_PI);              // energy density prefactor
const double hbarc = 0.197326938;                         // hbarc [GeV.fm]


// free parameters
//::::::::::::::::::::::::::::::::::::::::
const double T0_GeV = 0.6;                // initial temperature (GeV)
const int N_iterations = 10;              // number of Landau matching iterations
const bool tau_R_constant = false;        // constant relaxation time switch (true => constant, false => conformal)
const double tau_R = 0.5;                 // constant relaxation time value (the particles are still massless; no particle mass parameter built in)

const double etas = 3. / (4. * M_PI);     // shear viscosity
const double xi0 = 0.;                    // anisotropy parameter xi = (-1, infty)

const int N_tau = 4001;                   // longitudinal proper time points (uniform grid)
const double tau_min = 0.25;
const double tau_max = 20.25;

const int gauss_pts = 48;                 // gauss-legendre integration points for the second term (32, 48 or 100)
//::::::::::::::::::::::::::::::::::::::::


// cubic spline interpolation
double cubic_spline(double x, gsl_spline *f_spline)
{
  gsl_interp_accel * accel = gsl_interp_accel_alloc();

  double interpolation = gsl_spline_eval(f_spline, x, accel);

  gsl_interp_accel_free(accel);

  return interpolation;
}


// relaxation time: tau_R
double tau_R_function(double T)
{
  if(tau_R_constant) return tau_R;
  else return (5. * etas) / T;
}


// inverse relaxation time: 1/tau_R
double tau_R_inverse(double T)
{
  if(tau_R_constant) return 1./tau_R;
  else return T / (5. * etas);
}

// inverse shear relaxation time: 1/tau_pi
double tau_pi_inverse(double t, gsl_spline *T_spline)
{
  double T = cubic_spline(t, T_spline);
  return T / (5. * etas);
}


// z(t2,t1)
double z_function(double t2, double t1, int pts, double *root, double *weight, gsl_spline *T_spline)
{
  if(tau_R_constant)                                  // constant tau_R case
  {
    return (t2 - t1) / tau_R;
  }
  else if(fabs(t2 - t1) <= 1.e-16)                    // time-dependent tau_R case
  {
    return 0.;
  }
  else
  {
    double z = 0.;

    double average = (t1 + t2)/2.;
    double delta = (t2 - t1)/2.;

    for(int i = 0; i < pts; i++)                      // gauss-legendre integrate z
    {
      double t = average  +  delta * root[i];         // gauss-legendre transform

      z += weight[i] * tau_pi_inverse(t, T_spline);
    }
      z *= delta;                                     // prefactor from substitution

      return z;
  }
}


// dampening function: D(t2,t1) = exp[-z(t2,t1)]
double D_function(double t2, double t1, bool tau_R_constant, double tau_R, int pts, double *root, double *weight, gsl_spline *T_spline)
{
  if(tau_R_constant)                                  // constant tau_R case
  {
    return exp(- (t2 - t1) / tau_R);
  }
  else if(fabs(t2 - t1) <= 1.e-16)                    // time-dependent tau_R case
  {
    return 1.;
  }
  else
  {
    double z = 0.;

    double average = (t1 + t2)/2.;
    double delta = (t2 - t1)/2.;

    for(int i = 0; i < pts; i++)                      // gauss-legendre integrate z
    {
      double t = average  +  delta * root[i];         // gauss-legendre transform

      z += weight[i] * tau_pi_inverse(t, T_spline);
    }
      z *= delta;                                     // prefactor from substitution

      return exp(-z);
  }
}



int main()
{
  // load gauss data
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  double root[gauss_pts];
  double weight[gauss_pts];

  char file[255] = "";
  sprintf(file, "tables/gauss_legendre_%dpts.dat", gauss_pts);
  FILE * gauss_file = fopen(file, "r");
  if(gauss_file == NULL) printf("Error: couldn't open gauss legendre file\n");

  for(int i = 0; i < gauss_pts; i++)
  {
    fscanf(gauss_file, "%lf\t%lf", &root[i], &weight[i]);
  }
  fclose(gauss_file);
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  // uniform longitudinal proper time grid
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  double tau[N_tau];

  const double dtau = (tau_max - tau_min)  /  ((double)(N_tau - 1));

  for(int i = 0; i < N_tau; i++)
  {
    tau[i] = tau_min  +  dtau * ((double)i);
  }
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  // starting approximation for temperature profile
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  double T0 = T0_GeV / hbarc;                               // initial temperature [fm^-1]
  double Temp[N_tau];                                       // temperature grid [fm^-1]
  Temp[0] = T0;

  double L0 = T0 / pow(HE_function(1./sqrt(1.+xi0)), 0.25); // initial effective temperature [fm^-1]

  if(!tau_R_constant)                                       // aniso hydro approximation
  {
    // initial longitudinal pressure
    double pl0 = 1.5 * pow(L0, 4) / (M_PI * M_PI) * HL_function(1./sqrt(1. + xi0));

    run_aniso_bjorken(Temp, pl0, tau_min, dtau, N_tau, prefactor, etas);

    FILE *fp_temp0;
    fp_temp0 = fopen("results/T_aniso.dat", "w");
    for(int i = 0; i < N_tau; i++)
    {
      fprintf(fp_temp0, "%.5f %.8f\n", tau[i], Temp[i]);

    }
    fclose(fp_temp0);
  }
  else                                                      // ideal hydro approximation
  {
    for(int i = 1; i < N_tau; i++)
    {
      Temp[i] = T0 * pow(tau_min / tau[i], 1./3.);
    }
  }

  FILE *fp_temp_ideal;
  fp_temp_ideal = fopen("results/T_ideal.dat", "w");
  for(int i = 0; i < N_tau; i++) fprintf(fp_temp_ideal, "%.5f %.8f\n", tau[i], T0 * pow(tau_min / tau[i], 1./3.));
  fclose(fp_temp_ideal);

  double Energy[N_tau];                                     // energy density grid [fm^-4]
  Energy[0] = prefactor * pow(T0, 4);
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  // evolve Bjorken energy density and update temperature
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for(int n = 0; n < N_iterations; n++)                     // loop over iterations
  {
    gsl_spline * T_spline;                                  // construct cubic spline interpolation
    T_spline = gsl_spline_alloc(gsl_interp_cspline, N_tau);
    gsl_spline_init(T_spline, tau, Temp, N_tau);

    for(int itau = 1; itau < N_tau; itau++)                 // loop over proper time points
    {
      double t0 = tau_min;
      double t = tau[itau];

      double s_min = t0;                                    // s = tprime [t0, t]
      double s_max = t;

      double s_avg = (s_max + s_min)/2.;
      double s_delta = (s_max - s_min)/2.;

      double integral = 0.;

      for(int is = 0; is < gauss_pts; is++)                 // gauss integration s = tprime from [t0, t]
      {
        double s = s_avg  +  s_delta * root[is];            // gauss-legendre transform
        double T_s = cubic_spline(s, T_spline);             // T(s)

        double energy_s = prefactor * T_s * T_s * T_s * T_s;
        double D_s = D_function(t, s, tau_R_constant, tau_R, gauss_pts, root, weight, T_spline);
        double HE = HE_function(s / t);

        integral += weight[is] * D_s * tau_R_inverse(T_s) * energy_s * HE;
      } // tprime loop

      integral *= s_delta;                                  // prefactor from substitution

      double D = D_function(t, t0, tau_R_constant, tau_R, gauss_pts, root, weight, T_spline);
      double HE = HE_function(t0 / t / sqrt(1.+xi0));
      double HE0 = HE_function(1. / sqrt(1.+xi0));

      Energy[itau] = D * Energy[0] * HE / HE0 +  integral;

    } // t loop

    for(int itau = 0; itau < N_tau; itau++)                 // update temperature
    {
      Temp[itau] = pow(Energy[itau] / prefactor, 0.25);
    }

    gsl_spline_free(T_spline);

    printf("Finished iteration %d\n", n + 1);

  } // iteration loop
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  gsl_spline * T_spline;                                    // construct cubic spline interpolation
  T_spline = gsl_spline_alloc(gsl_interp_cspline, N_tau);   // of final iteration
  gsl_spline_init(T_spline, tau, Temp, N_tau);


  // compute the normalized shear stress pibar = pi / Peq
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  double pibar[N_tau];
  pibar[0] = pow(L0 / T0, 4) *  (HT_function(1./sqrt(1.+xi0))/2. - HL_function(1./sqrt(1. + xi0)));


  for(int itau = 1; itau < N_tau; itau++)
  {
    double t0 = tau_min;
    double t = tau[itau];
    double T = Temp[itau];

    double integral = 0.;

    double s_min = t0;                                      // s = tprime [t0, t]
    double s_max = t;

    double s_avg = (s_max + s_min)/2.;
    double s_delta = (s_max - s_min)/2.;

    for(int is = 0; is < gauss_pts; is++)                   // gauss integration s = tprime from [t0, t]
    {
      double s = s_avg  +  s_delta * root[is];              // gauss-legendre transform
      double T_s = cubic_spline(s, T_spline);               // T(s)
      double T4_s = T_s * T_s * T_s * T_s;

      double D_s = D_function(t, s, tau_R_constant, tau_R, gauss_pts, root, weight, T_spline);
      double HT = HT_function(s / t);
      double HL = HL_function(s / t);

      integral += weight[is] * D_s * tau_R_inverse(T_s) * T4_s * (HT/2. - HL);
    } // tprime loop

    integral *= (s_delta / pow(T,4));

    double D = D_function(t, t0, tau_R_constant, tau_R, gauss_pts, root, weight, T_spline);
    double HT = HT_function(t0 / t / sqrt(1.+xi0));
    double HL = HL_function(t0 / t / sqrt(1.+xi0));

    pibar[itau] = D * pow(L0 / T, 4) *  (HT/2. - HL)  +  integral;
  }
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  // compute derivatives of pibar
  double pibar_derivative[N_tau];
  double pibar_second_derivative[N_tau];
  double pibar_third_derivative[N_tau];
  double pibar_fourth_derivative[N_tau];

  compute_first_derivative(pibar_derivative, pibar, N_tau, dtau);
  compute_second_derivative(pibar_second_derivative, pibar, N_tau, dtau);
  compute_third_derivative(pibar_third_derivative, pibar, N_tau, dtau);
  compute_fourth_derivative(pibar_fourth_derivative, pibar, N_tau, dtau);

  // write results to file
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  FILE *fp_temp, *fp_tau_R, *fp_z;
  FILE *fp_pibar, *fp_pibar_derivative, *fp_pibar_2nd_derivative, *fp_pibar_3rd_derivative, *fp_pibar_4th_derivative;
  fp_temp = fopen("results/T_exact.dat", "w");
  fp_tau_R = fopen("results/tau_r_exact.dat", "w");
  fp_z = fopen("results/z_exact.dat", "w");

  fp_pibar = fopen("results/pibar_exact.dat", "w");
  fp_pibar_derivative = fopen("results/pibar_derivative_exact.dat", "w");
  fp_pibar_2nd_derivative = fopen("results/pibar_2nd_derivative_exact.dat", "w");
  fp_pibar_3rd_derivative = fopen("results/pibar_3rd_derivative_exact.dat", "w");
  fp_pibar_4th_derivative = fopen("results/pibar_4th_derivative_exact.dat", "w");

  for(int i = 0; i < N_tau; i++)
  {
    double t = tau[i];
    double T = Temp[i];
    double tau_R = tau_R_function(T);
    double z = z_function(t, tau_min, gauss_pts, root, weight, T_spline);

    fprintf(fp_temp, "%.5f %.8f\n", t, T);
    fprintf(fp_tau_R, "%.5f %.16f\n", t, tau_R);
    fprintf(fp_z, "%.5f %.16f\n", t, z);
    fprintf(fp_pibar, "%.5f %.16f\n", t, pibar[i]);
    fprintf(fp_pibar_derivative, "%.5f %.16f\n", t, pibar_derivative[i]);
    fprintf(fp_pibar_2nd_derivative, "%.5f %.16f\n", t, pibar_second_derivative[i]);
    fprintf(fp_pibar_3rd_derivative, "%.5f %.16f\n", t, pibar_third_derivative[i]);
    fprintf(fp_pibar_4th_derivative, "%.5f %.16f\n", t, pibar_fourth_derivative[i]);
  }

  fclose(fp_temp);
  fclose(fp_tau_R);
  fclose(fp_z);
  fclose(fp_pibar);
  fclose(fp_pibar_derivative);
  fclose(fp_pibar_2nd_derivative);
  fclose(fp_pibar_3rd_derivative);
  fclose(fp_pibar_4th_derivative);

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  gsl_spline_free(T_spline);
  printf("\n");
  return 0;
}








