/** 
 * Implementation file for hydroTheories.h
 *
 * Author: Kevin Ingles
 * ----------------------------------------
 */
 
#include "hydroTheories.h"

extern double eb, T0, g, pi, hbarc;

// RK4 solution to Eq. 11, returns the value of pi_bar at tau_end
double evolve_pi_bar(double init_tau,double tau_end,double dtau,double pi_bar0,double e0,hydro_coeff coeffs)
{
	double currentStep = init_tau;
	double pi_bar = pi_bar0; ///< Correct
	double e = e0;
	double e1, e2, e3, e4;
	double pibar1, pibar2, pibar3, pibar4;
	
	std::fstream fout;
	if (coeffs.lambda == 0)
		fout.open("output/hydroTheories/MIS_RK4.dat",std::fstream::out);
	else if (coeffs.lambda == 10./21.)
		fout.open("output/hydroTheories/DNMR_RK4.dat",std::fstream::out);
	else 
		fout.open("output/hydroTheories/aHydro_RK4.dat",std::fstream::out);
	
	while (currentStep < tau_end)
	{
		e1     = dtau * calc_de_dtau     (currentStep           ,pi_bar             ,e);
		pibar1 = dtau * calc_dpi_bar_dtau(currentStep           ,pi_bar             ,e         ,coeffs);
		
		e2     = dtau * calc_de_dtau     (currentStep + 0.5*dtau,pi_bar + 0.5*pibar1,e + 0.5*e1);
		pibar2 = dtau * calc_dpi_bar_dtau(currentStep + 0.5*dtau,pi_bar + 0.5*pibar1,e + 0.5*e1,coeffs);
		
		e3     = dtau * calc_de_dtau     (currentStep + 0.5*dtau,pi_bar + 0.5*pibar2,e + 0.5*e2);
		pibar3 = dtau * calc_dpi_bar_dtau(currentStep + 0.5*dtau,pi_bar + 0.5*pibar2,e + 0.5*e2,coeffs);
		
		e4     = dtau * calc_de_dtau     (currentStep + dtau    ,pi_bar + pibar3    ,e +     e3);
		pibar4 = dtau * calc_dpi_bar_dtau(currentStep + dtau    ,pi_bar + pibar3    ,e +     e3,coeffs);
		
		e      += (1.0/6.0)*(e1     + 2.0*e2     + 2.0*e3      + e4);
		pi_bar += (1.0/6.0)*(pibar1 + 2.0*pibar2 + 2.0*pibar3 + pibar4);
		fout << currentStep << "\t" << e << "\t" << pi_bar << std::endl;
		currentStep += dtau;
	}
	
	return pi_bar;
}

// Calculate value for Eq. 28 RHS
double calc_dpi_bar_dtau(double tau,double pi_bar,double e,hydro_coeff coeffs)
{
	double a      = coeffs.a;
	double lambda = coeffs.lambda;
	double gamma  = coeffs.gamma;
	double f      = coeffs.f;
	double pi_bar2 = pi_bar*pi_bar;
	
	double temp = pow(pi*pi*e/6,0.25);
	double tau_pi = 5.0*eb/temp;
	
	double term1 = -pi_bar/tau_pi + (1.0/tau)*(a - lambda*pi_bar - gamma*pi_bar2 - f*cal_F(pi_bar));
	
	return term1;
}

// Calculate value for Eq. 27 RHS
double calc_de_dtau(double tau,double pi_bar,double e)
{
	double result = (4./3.)*(e/tau)*(pi_bar - 1.0);
	return result;
}

double cal_F(double z)
{
	double z2 = z*z;
	double z3 = z*z*z;
	// Values received from Chandrodoy
	double a0 = 0.2;
	double a1 = -1.18526;
	double a2 = 1.30385;
	double a3 = 0.948719;
	
	// Vlaue where to keep function constant
	double z_flat = 0.333286;
	if (z <= -0.5)
		return 1.0;
	else if (z >= z_flat)
		return -0.0150767;
	else 
		return a0 + a1*z + a2*z2 + a3*z3;
}

/*
	So the plan after debugging this:
		1. Need to adjust exact code to calculate the pressure pi eta eta component.
		2. Determine how I want to infer the value of eta bar from the calculated quanities
		
*/