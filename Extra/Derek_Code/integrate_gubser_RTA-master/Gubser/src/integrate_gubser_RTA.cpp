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
#include <chrono>
#include <thread>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "coord_functions.cpp"
#include "soln_functions.cpp"
#include "Parameter.h"
#include "FileIO.cpp"
#include "Memoryf.cpp"

int main(void)
{
  std::cout << "####################### "<< "\n";
  printf("This code solves the distribution function f(tau, x, y; |p|, phi_p, v_z)... \n");
  printf("using the solution to the RTA Boltzmann EQN for Gubser symmetry. \n");

  //read in the parameters
  parameters params;
  std::cout << "reading in parameters" << "\n";
  readInParameters(params);

  int n_rho = params.n_rho; //grid size for rho
  //float c = params.c; // 5 eta / s constant
  //float rho0 = params.rho0;
  //float T_hat0 = params.T_hat0;
  int n_grid_rho = params.n_grid_rho;
  float rho_min = params.rho_min;
  float delta_rho = params.delta_rho;
  int n_grid_p_omega_hat = params.n_grid_p_omega_hat;
  float p_omega_hat_min = params.p_omega_hat_min;
  float delta_p_omega_hat = params.delta_p_omega_hat;
  int n_grid_p_eta_hat = params.n_grid_p_eta_hat;
  float p_eta_hat_min = params.p_eta_hat_min;
  float delta_p_eta_hat = params.delta_p_eta_hat;

  std::cout << "allocating space "<< "\n";
  double rho[n_rho]; //gubser time rho
  double e_rho[n_rho]; //real part of hatted energy density
  double T_rho[n_rho];
  float dummy;

  std::cout << "Reading ^e profile from input/energy.dat "<< "\n";
  std::ostringstream edat_stream;
  edat_stream << "input/energy.dat";
  std::ifstream edat(edat_stream.str().c_str());
  for (int i = 0; i < n_rho; i++)
  {
    // file format is :  rho Re{e} Im{e}
    edat >> rho[i] >> e_rho[i] >> dummy;
    //std::cout << rho[i] << "\t" << e_rho[i] << "\n";
  }
  //now find the hatted Temperature
  for (int i = 0; i < n_rho; i++) T_rho[i] = pow(M_PI * M_PI * e_rho[i] / 3., 0.25);

  //now create a gsl interp of T_hat(rho)
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *T_spline = gsl_spline_alloc(gsl_interp_cspline, n_rho);
  gsl_spline_init(T_spline, rho, T_rho, n_rho);

  //calculate the distribution function for various values of rho, p_omega_hat and p_eta_hat
  std::cout << "generating f(rho, p_omega_hat, p_eta_hat) "<< "\n";
  float f_soln_gubser[n_grid_rho][n_grid_p_omega_hat][n_grid_p_eta_hat];
  auto t1 = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for collapse(3)
  for (int irho = 0; irho < n_grid_rho; irho++)
  {
    float rho = rho_min + (float)(irho) * delta_rho;
    //std::cout << "rho = " << rho << "\n";
    for (int ip_omega_hat = 0; ip_omega_hat < n_grid_p_omega_hat; ip_omega_hat++)
    {
      float p_omega_hat = p_omega_hat_min + (float)(ip_omega_hat) * delta_p_omega_hat;
      for (int ip_eta_hat = 0; ip_eta_hat < n_grid_p_eta_hat; ip_eta_hat++)
      {
        float p_eta_hat = p_eta_hat_min + (float)(ip_eta_hat) * delta_p_eta_hat;
        f_soln_gubser[irho][ip_omega_hat][ip_eta_hat] = f_solution(rho, p_omega_hat, p_eta_hat, T_spline, acc, params);
      } // for (int ip_eta_hat = 0; ...
    } // for (int ip_omega_hat = 0; ...
  } // for (int irho = 0; ...

  //now dump the solution to file
  std::ofstream myfile;
  myfile.open("output/f_rho_pomegahat_petahat.dat");
  for (int irho = 0; irho < n_grid_rho; irho++)
  {
    float rho = rho_min + (float)(irho) * delta_rho;
    for (int ip_omega_hat = 0; ip_omega_hat < n_grid_p_omega_hat; ip_omega_hat++)
    {
      float p_omega_hat = p_omega_hat_min + (float)(ip_omega_hat) * delta_p_omega_hat;
      for (int ip_eta_hat = 0; ip_eta_hat < n_grid_p_eta_hat; ip_eta_hat++)
      {
        float p_eta_hat = p_eta_hat_min + (float)(ip_eta_hat) * delta_p_eta_hat;
        myfile << rho << " " << p_omega_hat << " " << p_eta_hat << " " << f_soln_gubser[irho][ip_omega_hat][ip_eta_hat] << "\n";
      } // for (int ip_eta_hat = 0; ...
    } // for (int ip_omega_hat = 0; ...
  } // for (int irho = 0; ...
  myfile.close();
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Done generating f(rho, p_omega_hat, p_eta_hat) "<< "\n";
  std::cout << "took " << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds\n";

  int n_r = params.n_r;
  float delta_r = params.delta_r;
  int n_p_mag = params.n_p_mag;
  float delta_p_mag = params.delta_p_mag;
  int n_phip = params.n_phip;
  float delta_phip = 2. * M_PI / (n_phip);
  int n_vz = params.n_vz;
  float delta_vz = 1.0 / n_vz;
  float tau0 = params.tau0;
  float q0 = params.q0;


  //Now generate the distribution function in milne coordinates at a fixed proper time
  // f(tau_0, x, y; |p|, phi_p, v_z)
  //because of the rotational symmetry, we can compute the distribution function at some fixed coordinate
  //transverse angle phi, and then compute the distribution function as a function of momentum phi_p where phi_p is
  //defined w.r.t. phi
  //then we only need to compute it as a function of radius, at a fixed angle....
  /*
  float f_soln_milne[n_r][n_p_mag][n_phip][n_vz];
  t1 = std::chrono::high_resolution_clock::now();
  std::cout << "generating f(tau_0, r; |p|, phi_p, v_z) "<< "\n";
  #pragma omp parallel for collapse(4)
  for (int ir = 0; ir < n_r; ir++)
  {
    //spacetime coordinates in milne
    float r = ir * delta_r;
    //phi = 0.
    //sin_phi = np.sin(phi)
    //cos_phi = np.cos(phi)
    float sin_phi = 0.;
    float cos_phi = 1.;
    //spacetime coordinates in Gubser
    float rho = rho_func(tau0, r, q0);
    float theta = theta_func(tau0, r, q0);
    for (int ip_mag = 0; ip_mag < n_p_mag; ip_mag++)
    {
      float p_mag = ip_mag * delta_p_mag;
      for (int ivz = 0; ivz < n_vz; ivz++)
      {
        float vz = ivz * delta_vz;
        float sin_theta_p = sqrt(1. - vz*vz);
        float cos_theta_p = vz;
        for (int iphip = 0; iphip < n_phip; iphip++)
        {
          float phip = iphip * delta_phip;
          float vx = cos_theta_p * cos(phip);
          float vy = cos_theta_p * sin(phip);
          //use mass shell condition p_mu p^mu = 0
          float p_T = p_mag * sin_theta_p; // p_T = |p| sin \theta_p
          float p_eta = p_mag * tau0 * vz;
          float p_tau = sqrt( p_T*p_T + ( (p_eta*p_eta) / (tau0*tau0) ) );
          float p_phi = p_mag * (-sin_phi * vx + cos_phi * vy);
          float p_theta_hat = p_theta_func(rho, theta, q0, p_tau, p_T);
          float p_omega_hat = p_omega_hat_func(p_theta_hat, p_phi, theta);
          //check that arguments lie inside bounds
          //p_omega_hat = np.minimum(p_omega_hat, p_omega_hat_grid[-1])
          f_soln_milne[ir][ip_mag][iphip][ivz] = f_solution(rho, p_omega_hat, p_eta, T_spline, acc, params);
        } // for (iphip = 0; iphip < n_phip; iphip++)
      } //for (ivz = 0; ivz < n_vz; ivz++)
    }
  } //for (ir = 0; ir < n_r; ir++)
  t2 = std::chrono::high_resolution_clock::now();

  //now write the solution to a file
  myfile.open("output/f_r_p_phip_vz.dat");
  for (int ir = 0; ir < n_r; ir++)
  {
    float r = ir * delta_r;
    for (int ip_mag = 0; ip_mag < n_p_mag; ip_mag++)
    {
      float p_mag = ip_mag * delta_p_mag;
      for (int ivz = 0; ivz < n_vz; ivz++)
      {
        float vz = ivz * delta_vz;
        for (int iphip = 0; iphip < n_phip; iphip++)
        {
          float phip = iphip * delta_phip;
          myfile << r << " " << p_mag << " " << phip << " " << vz << " " << f_soln_milne[ir][ip_mag][iphip][ivz] << "\n";
        } // for (iphip = 0; iphip < n_phip; iphip++)
      } //for (ivz = 0; ivz < n_vz; ivz++)
    }
  } //for (ir = 0; ir < n_r; ir++)
  myfile.close();
  std::cout << "Done generating f(tau_0, r; |p|, phi_p, v_z) "<< "\n";
  std::cout << "took " << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds\n";
  */

  //Now generate the solution on a regular grid in x,y in the positive quadrant of the xy plane
  float *****f_soln_milne_2 = NULL;
  f_soln_milne_2 = calloc5dArrayf(f_soln_milne_2, n_r, n_r, n_p_mag, n_phip, n_vz);

  t1 = std::chrono::high_resolution_clock::now();
  std::cout << "generating f(tau_0, x, y; |p|, phi_p, v_z) "<< "\n";
  #pragma omp parallel for collapse(5)
  for (int ix = 0; ix < n_r; ix++)
  {
    float x = ix * delta_r;
    for (int iy = 0; iy < n_r; iy++)
    {
      float y = iy * delta_r;
      //spacetime coordinates in milne
      float r = sqrt(x*x + y*y);
      float phi = atan2(y, x);
      float sin_phi = sin(phi);
      float cos_phi = cos(phi);
      //spacetime coordinates in Gubser
      float rho = rho_func(tau0, r, q0);
      float theta = theta_func(tau0, r, q0);
      for (int ip_mag = 0; ip_mag < n_p_mag; ip_mag++)
      {
        float p_mag = ip_mag * delta_p_mag;
        for (int ivz = 0; ivz < n_vz; ivz++)
        {
          float vz = ivz * delta_vz;
          float sin_theta_p = sqrt(1. - vz*vz);
          float cos_theta_p = vz;
          for (int iphip = 0; iphip < n_phip; iphip++)
          {
            float phip = iphip * delta_phip;
            float vx = cos_theta_p * cos(phip);
            float vy = cos_theta_p * sin(phip);
            //use mass shell condition p_mu p^mu = 0
            float p_T = p_mag * sin_theta_p; // p_T = |p| sin \theta_p
            float p_eta = p_mag * tau0 * vz;
            float p_tau = sqrt( p_T*p_T + ( (p_eta*p_eta) / (tau0*tau0) ) );
            float p_phi = p_mag * (-sin_phi * vx + cos_phi * vy);
            float p_theta_hat = p_theta_func(rho, theta, q0, p_tau, p_T);
            float p_omega_hat = p_omega_hat_func(p_theta_hat, p_phi, theta);
            //check that arguments lie inside bounds
            //p_omega_hat = np.minimum(p_omega_hat, p_omega_hat_grid[-1])
            f_soln_milne_2[ix][iy][ip_mag][iphip][ivz] = f_solution(rho, p_omega_hat, p_eta, T_spline, acc, params);
          } // for (iphip = 0; iphip < n_phip; iphip++)
        } //for (ivz = 0; ivz < n_vz; ivz++)
      }
    }
  } //for (ir = 0; ir < n_r; ir++)
  t2 = std::chrono::high_resolution_clock::now();

  //now write the solution to a file
  myfile.open("output/f_x_y_p_phip_vz.dat");
  for (int ix = 0; ix < n_r; ix++)
  {
    float x = ix * delta_r;
    for (int iy = 0; iy < n_r; iy++)
    {
      float y = iy * delta_r;
      for (int ip_mag = 0; ip_mag < n_p_mag; ip_mag++)
      {
        float p_mag = ip_mag * delta_p_mag;
        for (int ivz = 0; ivz < n_vz; ivz++)
        {
          float vz = ivz * delta_vz;
          for (int iphip = 0; iphip < n_phip; iphip++)
          {
            float phip = iphip * delta_phip;
            myfile << x << " " << y << " " << p_mag << " " << phip << " " << vz << " " << f_soln_milne_2[ix][iy][ip_mag][iphip][ivz] << "\n";
          } // for (iphip = 0; iphip < n_phip; iphip++)
        } //for (ivz = 0; ivz < n_vz; ivz++)
      }
    }
  } //for (ir = 0; ir < n_r; ir++)
  myfile.close();

  std::cout << "Done generating f(tau_0, x, y; |p|, phi_p, v_z) "<< "\n";
  std::cout << "took " << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds\n";


  std::cout << "Finished. Goodbye! "<< "\n";
  gsl_spline_free(T_spline);
  gsl_interp_accel_free(acc);

  free5dArrayf(f_soln_milne_2, n_r, n_r, n_p_mag, n_phip);

}
