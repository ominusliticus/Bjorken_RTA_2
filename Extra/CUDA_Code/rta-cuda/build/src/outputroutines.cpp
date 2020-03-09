/*
 
 outputroutines.cpp
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

#include "rta.h"
#include "outputroutines.h"

void outputMeasurements(const int iter) {
    cout.width(dwidth); cout << iter;
    cout << setprecision(16);
    cout.width(dwidth); cout << pow(t4[0],0.25);
    cout.width(dwidth); cout << pow(t4[num/2],0.25);
    cout.width(dwidth); cout << pow(t4[num-1],0.25);
    cout << endl;
}

void outputTemperatureSnapshot(double *lt4, int iter, string label) {
    fstream out;
    char fname[255];
    sprintf(fname,"output/%s_%d.dat",label.c_str(),iter);
    out.open(fname, ios::out);
    out << setprecision(16);
    for (int i=0; i<num; i++) {
        out << t[i] << "\t" << pow(lt4[i],0.25) << endl;
    }
    out.close();
    return;
}

void outputArray(double *data, string label) {
    outputArray(data,label,1);
}

void outputArray(double *data, string label, int step) {
    fstream out;
    char fname[255];
    sprintf(fname,"output/%s.dat",label.c_str());
    out.open(fname, ios::out);
    out << setprecision(16);
    for (int i=0; i<num; i+=step) {
        out << t[i] << "\t" << data[i] << endl;
    }
    out.close();
    return;
}

void outputScaledArray(double *data, double *lt4, double fac, string label) {
    outputScaledArray(data,lt4,fac,label,1);
}
void outputScaledArray(double *data, double *lt4, double fac, string label, int step) {
    fstream out;
    char fname[255];
    sprintf(fname,"output/%s.dat",label.c_str());
    out.open(fname, ios::out);
    out << setprecision(16);
    for (int i=0; i<num; i+=step) {
        out << t[i]*pow(lt4[i],0.25)/HBARC/fac << "\t" << data[i] << endl;
    }
    out.close();
    return;
}

void outputDistribution(double*** array) {
    // output f in binary format
    ofstream myFile ("output/f.bin", ios::out | ios::binary);
    for (int k=0; k<num/fStep; k++) {
        for (int i=0; i<numPZ; i++) {
            for (int j=0; j<numPT; j++) {
                myFile.write((char*) &array[k][i][j], sizeof(double));
            }
        }
    }
    myFile.close();
}

void print_line() {
    for (int i=0;i<5*dwidth;i++) cout << "-"; cout << endl;
    return;
}

