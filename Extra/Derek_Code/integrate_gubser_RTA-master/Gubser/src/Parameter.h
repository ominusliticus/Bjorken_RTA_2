#pragma once
struct parameters
{
  //parameters for constructing solution in gubser coordinates
  int n_rho; //grid size for rho
  float c; // 5 eta / s constant
  float rho0; //initial gubse time
  float T_hat0; //initial ^T
  int n_grid_rho; //number of points to compute in rho
  float rho_min; //min value of rho
  float delta_rho; //rho spacing
  int n_grid_p_omega_hat; // number of points in p_omega_hat
  float p_omega_hat_min; //min value of p_omega_hat
  float delta_p_omega_hat; //p_omega_hat spacing
  int n_grid_p_eta_hat;
  float p_eta_hat_min;
  float delta_p_eta_hat;
  //parameters for constructing solution in milne coords
  int n_r;
  float delta_r;
  int n_p_mag;
  float delta_p_mag;
  int n_phip;
  int n_vz;
  float tau0;
  float q0;
  float resolution;
};
