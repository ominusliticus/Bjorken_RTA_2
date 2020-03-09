/** 
 * Hydrodynamic equations for Bjorken in RTA
 * The equation is taken from PRC 100, 034901 (2019) Eq. 28
 * The various coefficients determine the approximation 
 * being used.
 * 
 * Theory 	a 		lambda		gamma		f
 * ============================================
 * MIS		4/15	0			4/3			0
 * DNMR		4/15	10/21		4/3			0
 * aHydro 	5/12	4/3			4/3			3/4 
 *
 * The equations for aHydro were taken from PRD 94, 125003 (2016)
 * Author: Kevin Ingles
 * -------------------------------------------------
 */
 
#ifndef __hydroTheories_h__
#define __hydroTheories_h__

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

typedef struct 
{
	double a;
	double lambda;
	double gamma;
	double f;
} hydro_coeff;

// Equations of evolution for IS and DMNR
double evolve_pi_bar(double init_tau,double tau_end,double dtau,double pi_bar0,double e0,hydro_coeff coeffs);
double calc_dpi_bar_dtau(double tau,double pi_bar,double e,hydro_coeff coeffs);
double calc_de_dtau(double tau,double pi_bar,double e); 
double cal_F(double z); // Parameterization for this function was done by Chandrodoy
					// and will be published in a letter coming soon.
 #endif 