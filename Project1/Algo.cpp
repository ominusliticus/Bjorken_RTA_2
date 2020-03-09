#include "Algo.h"

#include <iostream>
#include <algorithm>



using namespace std;

int FindPoint(double x, vector<double> xIn)
{
	// find appropriate pointd for interpolation.
	// to speed things up, we will not loop over the entire input,
	// but only the five nearest points
	bool found = false;
	int i = 0;
	int NPTS = xIn.size();
	while (true)
	{
		if (x <= xIn[0] || x <= xIn[1])
		{
			i = 2;
			return i;
		}
		else if (x <= xIn[i] && i < NPTS-2)
		{
			return i;
		}
		else if (x >= xIn[NPTS - 2] || x >= xIn[NPTS - 1] || x >= xIn[NPTS - 3])
		{
			return NPTS - 2;;
		}
		else if (i > NPTS)
		{
			cout << "Unknown error. . .Terminating\n";
			exit(-1);
		}
		else
		{
			i++;
		}
	}
}

double LagrangeInterp(double x, vector<double> xIn, vector<double> yIn)
{
	int i = FindPoint(x, xIn);
	double sum;
	if (x == xIn[i])
	{
		return yIn[i];
	}
	else
	{
		// need to update to spline
		double a = xIn[i + 1];
		double b = xIn[i - 1];
		sum = yIn[i - 1] * (x - a) / (b - a) + yIn[i + 1] * (x - b) / (a - b);
	}
	return sum;
}

double LagrangeInterp(double x, gsl_spline* spline)
{
	gsl_interp_accel* a = gsl_interp_accel_alloc();
	double interp = gsl_spline_eval(spline, x, a);
	gsl_interp_accel_free(a);
	return interp;
}

struct Compare_Ascend2
{
	vector<double> fData;
	Compare_Ascend2(vector<double> d) : fData(d) {}

	bool operator()(int i1, int i2)
	{
		return fData[i1] < fData[i2];
	}
};

struct Compare_Descend2
{
	vector<double> fData;
	Compare_Descend2(vector<double> d) : fData(d) {}

	bool operator()(int i1, int i2)
	{
		return fData[i1] > fData[i2];
	}
};

void SortBins(int bins, vector<double> data, int* index, bool down)
{
	for (int i = 0; i < bins; i++)
		index[i] = i;

	if (down)
		sort(index, index + bins, Compare_Descend2(data));
	else
		sort(index, index + bins, Compare_Ascend2(data));
}