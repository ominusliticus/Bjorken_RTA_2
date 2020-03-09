/*
 
 outputroutines.h
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
*/

#ifndef __outputroutines_h__
#define __outputroutines_h__

const int dwidth = 24;

void outputMeasurements(const int iter);
void outputTemperatureSnapshot(double *data, int iter, string label);
void outputArray(double *data, string label);
void outputArray(double *data, string label, int step);
void outputScaledArray(double *data, double *temperature, double fac, string label);
void outputScaledArray(double *data, double *temperature, double fac, string label, int step);
void outputDistribution(double*** array);
void print_line();

#endif /* __outputroutines_h__ */
