#include "Bjorken_RTA.h"
#include "ErrorMessages.h"

#include <omp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;

// An idea is speed things up is to make the 
// temperature a parameter you pass. rather than
// calculate from interpolation


void Bjorken_RTA::SetParameter(const char* key, const char* val)
{
	if (strcmp(key, "grid_points") == 0)
		_grid_points = atoi(val);
	if (strcmp(key, "T0") == 0)
		_T0 = atof(val);
	if (strcmp(key, "teq") == 0)
		_teq = atof(val);
	if (strcmp(key, "fpieb") == 0)
		_fpieb = atof(val);
	if (strcmp(key, "t0") == 0)
		_t0 = atof(val);
	if (strcmp(key, "tf") == 0)
		_tf = atof(val);
	if (strcmp(key, "xi0") == 0)
		_xi0 = atof(val);
	return;
}

Bjorken_RTA::Bjorken_RTA(string param_file_name, string data_file_name, double momentum_limits[2])
{
	// Set all private variables to zero
	_grid_points = 0;
	_T0 = 0;
	_teq = 0;
	_fpieb = 0;
	_t0 = 0;
	_tf = 0;
	_xi0 = 0;
	_rta_const = 0;
	_count = 0;
	_mom_range_pT = 0;
	_mom_range_w = 0;


	// Open and store parameters: 
	cout << "Reading in and printing out simulation parameters." << endl;
	string commentmarker = "//";
	char space = ' ';
	char tab = '\t';

	const int maxline = 128; 
	char buffer[maxline];
	ifstream paramFile(param_file_name.c_str());
	if (!paramFile.is_open())
		fatalError("Failed to open " + param_file_name);

	while (!paramFile.eof()) {
		paramFile.getline(buffer, maxline, '\n');
		string line = buffer; 
		int length = strlen(buffer);
		if (line.substr(0, commentmarker.length()) != commentmarker && line.length() > 0) 
		{
			char key[64] = "", value[64] = ""; 
			int founddelim = 0;
			for (int i = 0; i < length; i++) 
			{
				if (buffer[i] == space || buffer[i] == tab) 
					founddelim = 1;
				else 
				{
					if (founddelim == 0) 
						key[strlen(key)] = buffer[i];
					else 
						value[strlen(value)] = buffer[i];
				}
			}
			if (strlen(key) > 0 && strlen(value) > 0) 
			{
				SetParameter(key, value);
				cout << key << " = " << value << endl;
			}
		}
	}
	cout << endl;
	paramFile.close();
	_rta_const = 5.0 * _fpieb / (16.0 * atan(1.0));

	// Open and store output from Strickland code
	cout << "Reading in and storing temperature file." << endl;
	ifstream dataFile(data_file_name.c_str());
	if (!dataFile.is_open())
		fatalError("Failed to open file " + data_file_name);

	do
	{
		int c;
		double a, b;
		dataFile >> c >> a >> b;
		// cout << a << "\t" << b << endl;
		tau.push_back(a);
		T_eff.push_back(b);
	} while (!dataFile.eof());
	dataFile.close();
	_count = (int)tau.size(); 

	cout << tau.size() << endl;

	// Set limits of integration
	cout << "Setting integration limits." << endl;
	_mom_range_pT = momentum_limits[1];
	_mom_range_w = momentum_limits[0];

	// initialize cubic spline
	double* _tau = new double[_count];
	double* _T_eff = new double[_count];
	_tau = tau.data();
	_T_eff = T_eff.data();
	
	_temp_spline = gsl_spline_alloc(gsl_interp_cspline, _count);
	gsl_spline_init(_temp_spline, _tau, _T_eff, _count);

	delete[] _tau;
	delete[] _T_eff;
}

Bjorken_RTA::~Bjorken_RTA()
{
	// free gsl_spline
	gsl_spline_free(_temp_spline);
}

double Bjorken_RTA::tau_r(double ttau)
{
	// Note that the units are fm/c
	double temp = LagrangeInterp(ttau, _temp_spline);
	return (hbarc * _rta_const / temp);
}

// 24-point gauss quadrature points and weights
const int NSUM48 = 24;
double x48[NSUM48] = { 0.0323801709628694,
					0.0970046992094627,
					0.1612223560688917,
					0.2247637903946891,
					0.2873624873554556,
					0.3487558862921608,
					0.4086864819907167,
					0.4669029047509584,
					0.5231609747222330,
					0.5772247260839727,
					0.6288673967765136,
					0.6778723796326639,
					0.7240341309238146,
					0.7671590325157404,
					0.8070662040294426,
					0.8435882616243935,
					0.8765720202742479,
					0.9058791367155696,
					0.9313866907065543,
					0.9529877031604309,
					0.9705915925462473,
					0.9841245837228269,
					0.9935301722663508,
					0.9987710072524261 };
double w48[NSUM48] = { 0.0647376968126839,
					0.0644661644359501,
					0.0639242385846482,
					0.0631141922862540,
					0.0620394231598927,
					0.0607044391658939,
					0.0591148396983956,
					0.0572772921004032,
					0.0551995036999842,
					0.0528901894851937,
					0.0503590355538545,
					0.0476166584924905,
					0.0446745608566943,
					0.0415450829434647,
					0.0382413510658307,
					0.0347772225647704,
					0.0311672278327981,
					0.0274265097083569,
					0.0235707608393244,
					0.0196161604573555,
					0.0155793157229438,
					0.0114772345792345,
					0.0073275539012763,
					0.0031533460523058 };

const int NSUM = 12;
double x24[NSUM] = { 0.064056892862605626085,
					0.191118867473616309159,
					0.315042679696163374387,
					0.433793507626045138487,
					0.545421471388839535658,
					0.648093651936975569252,
					0.740124191578554364244,
					0.820001985973902921954,
					0.886415527004401034213,
					0.938274552002732758524,
					0.974728555971309498198,
					0.995187219997021360180 };
double w24[NSUM] = { 0.127938195346752156974,
					0.125837456346828296121,
					0.121670472927803391204,
					0.115505668053725601353,
					0.107444270115965634783,
					0.097618652104113888270,
					0.086190161531953275917,
					0.073346481411080305734,
					0.059298584915436780746,
					0.044277438817419806169,
					0.028531388628933663181,
					0.012341229799987199547 };

double Bjorken_RTA::Decay_factor(double tau1, double tau2)
{
	double result = 0;

	///////////////////////////////////////////////////////////////////////////
	//  This it of code was introduced to see if the convergence at large    //
	//  tau could be improved. The idea is to partition the integral based   //
	//  on the multiple of the floor ratio of the upper limit to lower       //
	//  limit minus one. The rest is just machinery to ensure the integral   //
	//  is evaluated correctly.                                              //
	///////////////////////////////////////////////////////////////////////////
	/***************************Begin Code Snippet****************************/
	int its = 1; // (int)(tau2 / tau1);

	double ttaup = (tau2 - tau1) / (double)(its);

	// #pragma omp parallel for reduction(+: result)
	for (int k = 0; k < its; k++)
	{
		double ttau1 = tau1 + (double)k * ttaup;
		double ttau2 = tau1 + (double)(k + 1) * ttaup;

		double resultt = 0.0;
	/***************************End Code Snippet*******************************/
		//#pragma omp parallel for reduction(+: result)
		for (int n = 0; n < NSUM48; n++)
		{
			double tneg = ((ttau2 - ttau1) * (-x48[n]) + (ttau2 + ttau1)) / 2.0;
			double tpos = ((ttau2 - ttau1) * (x48[n]) + (ttau2 + ttau1)) / 2.0;

			resultt += w48[n] * (1.0 / tau_r(tpos) + 1.0 / tau_r(tneg));
		}
		result += resultt*(ttau2 - ttau1) / 2.0;
	}
	return exp(-result);
}

double Bjorken_RTA::R(double z)
{
	double result;
	if (z == 0.0)
		result = 2.0;
	else
		result = 1.0 / (1.0 + z) + atan(sqrt(z)) / sqrt(z);

	return result * 0.5;
}

double Bjorken_RTA::finit(double ww, double pT)
{
	double pi = 4.0 * atan(1.0);
	double temp = T_eff[0];
	double Lambda = temp / pow(R(_xi0), 0.25);
	double result = 2.0 / pow(2.0 * pi, 3.0);
	result *= exp(-sqrt((1.0 + _xi0) * ww * ww + pT * pT * tau[0] * tau[0]) / Lambda / tau[0]);
	return result;
}

double Bjorken_RTA::feq(double ttau, double ww, double pT)
{
	double pi = 4.0 * atan(1.0);
	double temp = LagrangeInterp(ttau, _temp_spline);
	if (ttau == tau[0])
		temp = T_eff[0];

	double result = 2.0 / pow(2.0 * pi, 3.0);
	result *= exp(-sqrt(ww * ww + pT * pT * ttau * ttau) / temp / ttau);
	return result;
}

double Bjorken_RTA::SPD(double ttau, double ww, double pT)
{
	if (ttau == tau[0])
		return finit(ww, pT);

	double result = 0.;
	double xneg, xpos;
	double tau0 = tau[0];

	///////////////////////////////////////////////////////////////////////////
	//  This it of code was introduced to see if the convergence at large    //
	//  tau could be improved. The idea is to partition the integral based   //
	//  on the multiple of the floor ratio of the upper limit to lower       //
	//  limit minus one. The rest is just machinery to ensure the integral   //
	//  is evaluated correctly.                                              //
	///////////////////////////////////////////////////////////////////////////
	/***************************Begin Code Snippet****************************/
	int its = 1; // (int)floor(ttau / tau[0]);
	 
	 double ttaup = (ttau - tau[0]) / (double)(its);
	 
	 // #pragma omp parallel for reduction(+: result)
	 for (int k = 0; k < its; k++)
	 {
	 	double ttau1 = tau0 + (double)k * ttaup;
		double ttau2 = tau0 + (double)(k + 1) * ttaup;
	 
	 	double resultt = 0.0;
	/***************************End Code Snippet*******************************/
		for (int n = 0; n < NSUM48; n++)
		{
			xneg = ((ttau2 - ttau1) * (-x48[n]) + (ttau2 + ttau1)) / 2.0;
			xpos = ((ttau2 - ttau1) * (x48[n]) + (ttau2 + ttau1)) / 2.0;
			// xneg = ((ttau - tau0) * (-x48[n]) + (ttau + tau0)) / 2.0;
			// xpos = ((ttau - tau0) * (x48[n]) + (ttau + tau0)) / 2.0;

			double trpos = 1.0 / tau_r(xpos);
			double trneg = 1.0 / tau_r(xneg);

			double dfpos = Decay_factor(xpos, ttau);
			double dfneg = Decay_factor(xneg, ttau);

			double fqpos = feq(xpos, ww, pT);
			double fqneg = feq(xneg, ww, pT);

			resultt += w48[n] * ((trpos * dfpos * fqpos) + (trneg * dfneg * fqneg));
			// result += w48[n] * ((trpos * dfpos * fqpos) + (trneg * dfneg * fqneg));
		}
		result += resultt * (ttau2 - ttau1) / 2.0;
		// result *= (ttau - tau0) / 2.0;
		// cout << resultt << "\t" << result << endl;
	}
	

	double res = Decay_factor(tau0, ttau) * finit(ww, pT) + result;
	// cout << "res = " << res << endl;
	
	return res;
}

double Bjorken_RTA::moment(double ttau, int* index_array, int deg)
{
	// if deg = 2 and index_arry = {0,0}, then we are integrating 
	// over all momenta with the zero component of the momentum 4-vectors
	// against the SPD, i.e the energy density of the system.

	double hbarc3 = pow(hbarc, 3.0);
	double PI = 4.0 * atan(1.0);

	double result = 0.0;
	double wpos, wneg;
	double pTpos, pTneg;
	double temperature = LagrangeInterp(ttau, _temp_spline);
	// Add code to change to Gaussian Quadrature

	// private(temp_resultneg, temp_resultpos)
	// #pragma omp parallel for reduction(+: result) 
	// int repeats = 1;
	// double split_pT = _mom_range_pT / double(repeats);
	// double split_w = _mom_range_w / double(repeats);
	// for (int n = 0; n < repeats; n++)
	// {
	// 	int np1 = n + 1;
	// 	double w1 = double(n) * split_w;
	// 	double w2 = double(np1) * split_w;
	// 	double resultt = 0.0;
	// 	for (int m = 0; m < repeats; m++)
	// 	{
	// 		int mp1 = m + 1;
	// 		double pT1 = double(m) * split_pT;
	// 		double pT2 = double(mp1) * split_pT;
	// 		double resulttt = 0.0;
	// 		//for (int ix = 0; ix < NSUM; ix++)
			for (int ix = 0; ix < NSUM48; ix++)
			{
				wneg = (_mom_range_w) * (-x48[ix] + 1.0) / 2.0;
				wpos = (_mom_range_w) * (x48[ix] + 1.0) / 2.0;
				// wneg = ((w2 - w1) * (-x48[ix]) + (w2 + w1)) / 2.0;
				// wpos = ((w2 - w1) * (x48[ix]) + (w2 + w1)) / 2.0;

				double temp_resultneg = 0;
				double temp_resultpos = 0;
				// for (int iy = 0; iy < NSUM; iy++)
				for (int iy = 0; iy < NSUM48; iy++)
				{
					pTneg = (_mom_range_pT) * (-x48[iy] + 1.0) / 2.0;
					pTpos = (_mom_range_pT) * (x48[iy] + 1.0) / 2.0;
					// pTneg = ((pT2 - pT1) * (-x48[iy]) + (pT2 + pT1)) / 2.0;
					// pTpos = ((pT2 - pT2) * (x48[iy]) + (pT2 + pT1)) / 2.0;

					double vpospos = sqrt(wpos * wpos + pTpos * pTpos * ttau * ttau);
					double vposneg = sqrt(wpos * wpos + pTneg * pTneg * ttau * ttau);
					double vnegpos = sqrt(wneg * wneg + pTpos * pTpos * ttau * ttau);
					double vnegneg = sqrt(wneg * wneg + pTneg * pTneg * ttau * ttau);

					double argnegneg = 1.0;
					double argnegpos = 1.0;
					double argposneg = 1.0;
					double argpospos = 1.0;

					for (int n = 0; n < deg; n++)
					{
						if (index_array[n] == 0)
						{
							argnegneg *= vnegneg / ttau;
							argnegpos *= vnegpos / ttau;
							argposneg *= vposneg / ttau;
							argpospos *= vpospos / ttau;
						}
						else if (index_array[n] == 1 || index_array[n] == 2)
						{
							argnegneg *= pTneg;
							argnegpos *= pTpos;
							argposneg *= pTneg;
							argpospos *= pTpos;
						}
						else
						{
							argnegneg *= wneg / ttau / ttau;
							argnegpos *= wneg / ttau / ttau;
							argposneg *= wpos / ttau / ttau;
							argpospos *= wpos / ttau / ttau;
						} // end else
					} // end loop over deg

					// double t = temperature;
					// double xpo = 4.0;
					double t = 1;
					double xpo = 0.0;
					double Inegneg = argnegneg * pTneg * SPD(ttau, wneg * t, pTneg * t) / vnegneg;
					double Inegpos = argnegpos * pTpos * SPD(ttau, wneg * t, pTpos * t) / vnegpos;
					double Iposneg = argposneg * pTneg * SPD(ttau, wpos * t, pTneg * t) / vposneg;
					double Ipospos = argpospos * pTpos * SPD(ttau, wpos * t, pTpos * t) / vpospos;

					// temp_resultneg += 4.0 * PI * w24[iy] * pow(temperature, xpo) * (Inegneg + Inegpos) / hbarc3;
					// temp_resultpos += 4.0 * PI * w24[iy] * pow(temperature, xpo) * (Iposneg + Ipospos) / hbarc3;

					temp_resultneg += w48[iy] * pow(temperature, xpo) * (Inegneg + Inegpos) / hbarc3;
					temp_resultpos += w48[iy] * pow(temperature, xpo) * (Iposneg + Ipospos) / hbarc3;
				} // end loop over iy 
				result += w48[ix] * (temp_resultpos + temp_resultneg) * _mom_range_pT * _mom_range_w / 4.0;
				// resulttt += w48[ix] * (temp_resultpos + temp_resultneg);
			} // end loop ix
	// 		resultt += resulttt * split_pT / 2.0;
	// 	} // end loop m
	// 	result += resultt * split_w / 2.0;
	// } // end loop n
	return 4.0 * PI * result;
}

///////////////////////////////////////////////////////////////////////
///																	///	
///		Numerical Hydrodynamic Equation Solves for DNMR and IS      ///
///																	///
///////////////////////////////////////////////////////////////////////

double Bjorken_RTA::NavierStokes_Energy_Evolution(double init_tau, double tau_end, double dtau, double e0)
{
	double e1, e2, e3, e4;
	double e = e0;
	double tau = init_tau;
	double teq = _rta_const / 5.0;

	while (tau < tau_end)
	{
		cout << tau << "\t" << e << endl;
		e1 = -dtau * (4.0 / 3.0 / tau) * e * (1.0 - teq * (4.0 / 3.0) / tau);
		e2 = -dtau * (4.0 / 3.0 / (tau + 0.5 * dtau)) * (e + 0.5 * e1) * (1.0 - teq * (4.0 / 3.0) / (tau + 0.5 * dtau));
		e3 = -dtau * (4.0 / 3.0 / (tau + 0.5 * dtau)) * (e + 0.5 * e2) * (1.0 - teq * (4.0 / 3.0) / (tau + 0.5 * dtau));
		e4 = -dtau * (4.0 / 3.0 / (tau + dtau)) * (e + e3) * (1.0 - teq * (4.0 / 3.0) / (tau + dtau));

		e += (1.0) / (6.0) * (e1 + 2.0 * e2 + 2.0 * e3 + e4);
		tau += dtau;
	}
	return 0.0;
}

// RK4 solution to Eq. 11, returns the value of pi_bar at tau_end
double Bjorken_RTA::evolve_pi_bar(double init_tau, double tau_end, double dtau, double pi_bar0, double e0, hydro_coeff coeffs)
{
	double currentStep = init_tau;
	double pi_bar = pi_bar0; ///< Correct
	double e = e0;
	double e1, e2, e3, e4;
	double pibar1, pibar2, pibar3, pibar4;

	std::fstream fout;
	if (coeffs.lambda == 0)
		fout.open("output/hydroTheories/MIS_RK4.dat", std::fstream::out);
	else if (coeffs.lambda == 10. / 21.)
		fout.open("output/hydroTheories/DNMR_RK4.dat", std::fstream::out);
	else
		fout.open("output/hydroTheories/aHydro_RK4.dat", std::fstream::out);

	while (currentStep < tau_end)
	{
		e1     = dtau * calc_de_dtau(currentStep, pi_bar, e);
		pibar1 = dtau * calc_dpi_bar_dtau(currentStep, pi_bar, e, coeffs);

		e2     = dtau * calc_de_dtau(currentStep + 0.5 * dtau, pi_bar + 0.5 * pibar1, e + 0.5 * e1);
		pibar2 = dtau * calc_dpi_bar_dtau(currentStep + 0.5 * dtau, pi_bar + 0.5 * pibar1, e + 0.5 * e1, coeffs);

		e3     = dtau * calc_de_dtau(currentStep + 0.5 * dtau, pi_bar + 0.5 * pibar2, e + 0.5 * e2);
		pibar3 = dtau * calc_dpi_bar_dtau(currentStep + 0.5 * dtau, pi_bar + 0.5 * pibar2, e + 0.5 * e2, coeffs);

		e4     = dtau * calc_de_dtau(currentStep + dtau, pi_bar + pibar3, e + e3);
		pibar4 = dtau * calc_dpi_bar_dtau(currentStep + dtau, pi_bar + pibar3, e + e3, coeffs);

		e      += (1.0 / 6.0) * (e1 + 2.0 * e2 + 2.0 * e3 + e4);
		pi_bar += (1.0 / 6.0) * (pibar1 + 2.0 * pibar2 + 2.0 * pibar3 + pibar4);
		fout << currentStep << "\t" << e << "\t" << pi_bar << std::endl;
		currentStep += dtau;
	}

	return pi_bar;
}

// Calculate value for Eq. 28 RHS
double Bjorken_RTA::calc_dpi_bar_dtau(double tau, double pi_bar, double e, hydro_coeff coeffs)
{
	double pi = 4.0 * atan(1.0);
	double eb = _fpieb / (4.0 * pi);

	double a = coeffs.a;
	double lambda = coeffs.lambda;
	double gamma = coeffs.gamma;
	double f = coeffs.f;
	double pi_bar2 = pi_bar * pi_bar;

	double temp = pow(pi * pi * e / 6, 0.25);
	double tau_pi = 5.0 * eb / temp;

	double term1 = -pi_bar / tau_pi + (1.0 / tau) * (a - lambda * pi_bar - gamma * pi_bar2 - f * cal_F(pi_bar));

	return term1;
}

// Calculate value for Eq. 27 RHS
double Bjorken_RTA::calc_de_dtau(double tau, double pi_bar, double e)
{
	double result = (4. / 3.) * (e / tau) * (pi_bar - 1.0);
	return result;
}

double Bjorken_RTA::cal_F(double z)
{
	double z2 = z * z;
	double z3 = z * z * z;
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
		return a0 + a1 * z + a2 * z2 + a3 * z3;
}