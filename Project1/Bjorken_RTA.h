#pragma once

#include "Algo.h"

#include <string>
#include <vector>

static double hbarc = 0.1973269804;  ///< GeV fm

typedef struct
{
	double a;
	double lambda;
	double gamma;
	double f;
} hydro_coeff;

class Bjorken_RTA
{
public:
	// Constructor and destructor:
	// 
	// param_file_name - contains the input parameters to Strickland code
	// data_files_name - contains ouput from Strickland code
	// momentum_;imits - limits of integration
	Bjorken_RTA(std::string param_file_name, std::string data_file_name, double momentum_limits[2]); 
	~Bjorken_RTA();

	// Functions relevant to calculating quantities related to SPD:
	//
	// tau_r - the relaxation time function
	// Decay factor - exp( integrate( 1.0 / tau_r, tau2, tau1 ) )
	// R - thermodynamic integral from Strickland paper
	// feq - equilibrium Boltzmann SPD
	// finit - Romatschke-Strickland parameterization
	// SPD - the solution to the BE in RTA
	double tau_r(double ttau);
	double Decay_factor(double tau1, double tau2);
	double R(double z);
	double feq(double ttau, double ww, double pT);
	double finit(double ww, double pT);
	double SPD(double ttau, double ww, double pT);

	// Moments of the SPD:
	// We only care about time-time-component, of the second moment - stress-energy tensor
	double moment(double ttau, int* index_array, int deg); // deg - which moment of the SPD

	// For parameter processing
	void SetParameter(const char* key, const char* val);

	// Add sampling methods...

	// Hydrodynamic differential equations solvers
	// Equations of evolution for IS and DMNR
	double NavierStokes_Energy_Evolution(double init_tau, double tau_end, double dtau, double e0);
	double evolve_pi_bar(double init_tau, double tau_end, double dtau, double pi_bar0, double e0, hydro_coeff coeffs);
	double calc_dpi_bar_dtau(double tau, double pi_bar, double e, hydro_coeff coeffs);
	double calc_de_dtau(double tau, double pi_bar, double e);
	double cal_F(double z); // Parameterization for this function was done by Chandrodoy
							// and will be published in a letter coming soon.

	// Vector containing data from Strickland code
	std::vector<double> tau;
	std::vector<double> T_eff;

	inline int GetCount(void) const
	{
		return _count;
	}

	inline gsl_spline* getSplinePtr(void) const
	{
		return _temp_spline;
	}

private:
	// These are the variables that get used by the code from Strickland
	// They are included here to make it possible to caluclate the single
	// particle distribution (SPD) function using the same input
	int _grid_points;
	double _T0;
	double _teq;
	double _fpieb;
	double _t0;
	double _tf;
	double _xi0;
	double _rta_const;

	// Variable containing the data ouput from Strickland code
	int _count;

	// These are the parameters that are need to evaluate the integrals for
	// the moments of the SPD
	double _mom_range_pT;
	double _mom_range_w;

	gsl_spline* _temp_spline;
};

