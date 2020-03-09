#pragma once

#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

// Function that finds between which two points the value x lies
// in the array xIn
int FindPoint(double x, std::vector<double> xIn);

// Function that provides an interpolation polynomial for a given 
// array xIn, and evaluates is x
double LagrangeInterp(double x, std::vector<double> xIn, std::vector<double> yIn);
double LagrangeInterp(double x, gsl_spline* spline);

// Function that sorts an array array of doubles, by providing 
// the correct sequence of indices in an indexing array
void  SortBins(int bins, std::vector<double> data, int* index, bool down);


