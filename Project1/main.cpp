#include "ErrorMessages.h"
#include "Bjorken_RTA.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
	cout << "Running program: Bjorken_RTA Exact Solution" << endl;
	cout << "Author: Kevin Ingles" << endl;
	cout << "Contribution from the Ulrich Heinz Group and Michael Strickland" << endl;
	cout << "All rights reserved." << endl << endl;

	int a[2] = { 0, 0 };
	int d = 2;

	double mom_limits[2] = { 10.0, 10.0 }; ///< GeV
	string paramFile = "params.txt";
	string dataFile = "output/T_500.dat";
	Bjorken_RTA *myClass = new Bjorken_RTA(paramFile, dataFile, mom_limits);

	cout << "Class initialized!" << endl;

	double t0 = 1.0;
	double t = 4.90;
	double out0 = LagrangeInterp(t, myClass->getSplinePtr());
	double out1 = myClass->tau_r(t); 			
	double out2 = myClass->Decay_factor(t0, t); 
	double out3 = myClass->SPD(t, 2.0, 2.0);	
	double out4 = myClass->feq(t, 2.0, 2.0);	
	double out5 = myClass->finit(2.0, 2.0);		
	// double out4 = myClass->moment(10.4, a, d);

	cout << "temp:         " << out0 << endl;
	cout << "tau_r:        " << out1 << endl;
	cout << "Decay_factor: " << out2 << endl;
	cout << "SPD:          " << out3 << endl;
	cout << "feq:          " << out4 << endl;
	cout << "finit:        " << out5 << endl;

	/*double pTmax1 = 10;
	double pTmax2 = 50;
	ofstream pTdist1, pTdist2;
	pTdist1.open("output/pTdist_10GeV.dat");
	pTdist2.open("output/pTdist_50GeV.dat");
	int steps = 500;
	for (int i = 0; i < steps; i++)
	{
		double pt1 = double(i) * pTmax1 / double(steps);
		double pt2 = double(i) * pTmax2 / double(steps);

		pTdist1 << pt1 << "\t" << myClass->SPD(4.9, 0, pt1) << endl;
		pTdist2 << pt2 << "\t" << myClass->SPD(4.9, 0, pt2) << endl;
	}
	pTdist1.close();
	pTdist2.close();*/

	const int NSUM = 24;
	double x24[NSUM] = { 0.0323801709628694,
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
	double w24[NSUM] = { 0.0647376968126839,
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

	double pi = 4.0 * atan(1.0);
	double tau0 = myClass->tau[0];
	double tauf = 50;  ///< fm/c
	double dtau = 0.025;
	double e0 = 6.0 * pow(myClass->T_eff[0], 4.0) / pi / pi / pow(hbarc, 3.0);

	// myClass->NavierStokes_Energy_Evolution(tau0, tauf, dtau, e0);

	int a2[1] = { 1 };
	int a3[1] = { 3 };
	
	const int count = myClass->GetCount();
	double *tau = new double[count];
	double* T_eff = new double[count];
	tau = myClass->tau.data();
	T_eff = myClass->T_eff.data();

	// ofstream t_evolution, spd_surface;
	// t_evolution.open("output/t_evolution.dat");
	// spd_surface.open("output/spd_surface.dat");
	// double p = 2.0;
	// double pmax = 10;
	// for (int i = 0; i < count; i++)
	// {
	// 	double t = tau[i];
	// 	double feq = myClass->feq(t, p, p);
	// 	double teq = myClass->tau_r(t);
	// 	double decay = myClass->Decay_factor(tau0, t);
	// 	t_evolution << t << "\t" << feq << "\t" << teq << "\t" << decay << "\t" << endl;
	// 	for (int j = 0; j < count; j++)
	// 	{
	// 		double pp = double(j) * pmax / double(count);
	// 		double spd = myClass->SPD(t, 0, pp);
	// 		spd_surface << t << "\t" << pp << "\t" << spd << endl;
	// 	}
	// }
	// t_evolution.close();
	// spd_surface.close();

	ofstream fout("output/shear_Ingles_alt.dat");
	ifstream fin_PL("output/pl.dat"), fin_PT("output/pt.dat");
	if (!fin_PL.is_open() || !fin_PT.is_open())
		fatalError("Failed to open pt.dat or pl.dat");
	 
	cout << "tau\t e \t u \t pt \t PT \t pl \t pL" << endl;
	for (int i = 100; i < count; i++)
	{
		double t = tau[i];
		double T = T_eff[i];
		double u = 6.0 * pow(T, 4.0) / pi / pi / pow(hbarc, 3.0);
		double e = myClass->moment(t, a, d);
		// double pt = myClass->moment(t, a2, 1);
		double pl = myClass->moment(t, a3, 1);
		double pt2 = 0.5 * (e - pl);
	
		double PT, PL, temp;
		fin_PL >> temp >> PL;
		fin_PT >> temp >> PT;
	
		PT *= 2.0 * pow(hbarc, -3.0);
		PL *= 2.0 * pow(hbarc, -3.0);
		// cout << tau << "\t" << e << "\t" << u << endl;
		cout << t << "\t" << e << "\t" << u << "\t" /*<< pt << "\t"*/ << pt2 << "\t" << PT << "\t" << pl << "\t" << PL << endl;
		fout << t << "\t" << e << "\t" << u << "\t" /*<< pt << "\t"*/ << pt2 << "\t" << PT << "\t" << pl << "\t" << PL << endl;
	
		// double shear = (2.0 / 3.0) * (pl - pt);
		// double shear2 = (2.0 / 3.0) * (pl - pt2);
		// double SHEAR = (2.0 / 3.0) * (PL - PT);
		// cout << t << "\t" << shear << "\t" << shear2 << "\t" << SHEAR << endl;
		// fout << t << "\t" << shear << "\t" << shear2 << "\t" << SHEAR << endl;
	}
	
	// myClass->~Bjorken_RTA();
	return 0;
}