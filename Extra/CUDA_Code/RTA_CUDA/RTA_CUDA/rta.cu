/*
 
 rta.cu
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <helper_cuda.h>
#include <helper_functions.h>

using namespace std;

#include "rta.h"
#include "outputroutines.h"
#include "paramreader.h"
#include "memory.h"

// defines
#define TIDX(i,j) (j + i*(i + 1)/2)
#define BLOCKSIZE1  256
#define BLOCKSIZE2  128

// constants that are shared with the GPU
__constant__ int NUM;
__constant__ double DTAU,A_0,T_0,EB;
__constant__ int N,M;
__constant__ double PZ,PT;
__constant__ double M_PI;
__constant__ double hbarc;

// these global vars are initialized from parameters file
// defaults set here are overridden by that file
int    num = 100, maxiters = 10, update = 10, snapupdate = 20;

double fpieb = 1; // 4 Pi eta / S
double t0 = 0.25; // initial time in fm/c
double tf = 20; // final time in fm/c
double T0 = 0.6; // initial temperature in GeV
double a0 = 1; // initial anisotropy a0 = 1/sqrt(1+xi0)


// time step
double dt;

// this holds the current values of T^4
double *t4;

// this holds the updated values of T^4
double *T4;

// this holds the integration abscissae (timeGrid)
double *t;

// parameters for moment computation
int computeMoments=0, maxN=4, maxM=4;

// this will hold the final solution for the distribution function f for a fix w and pt
double *f;
int computeDist=0,numPZ=40,numPT=40, fStep=1;;
double maxPT=2, maxPZ=2;

// these hold the values of hnm and the initial value array for the general moment equation
double *hnm,*hnm0;

// these are pointers for the device memory
double *dev_t4, *dev_T4, *dev_time, *dev_d, *dev_h, *dev_hnm, *dev_hnm0, *dev_m, *dev_f;

/*----------------------------------------------------------------------------------------------------*/
// Special functions
/*----------------------------------------------------------------------------------------------------*/

__device__ double H(double y) {
    if (y==1) return 2;
    if (fabs(y)<1) return y*(fabs(y) + asin(sqrt(1-y*y))/sqrt(1-y*y));
    if (fabs(y)>1) return y*(fabs(y) + asinh(sqrt(y*y-1))/sqrt(y*y-1));
    return 0;
}

double hostH(double y) {
    if (y==1) return 2;
    if (fabs(y)<1) return y*(fabs(y) + asin(sqrt(1-y*y))/sqrt(1-y*y));
    if (fabs(y)>1) return y*(fabs(y) + asinh(sqrt(y*y-1))/sqrt(y*y-1));
    return 0;
}

double my2F1(double a, double b, double c, double z)
{
	if (fabs(z)<=1) return gsl_sf_hyperg_2F1(a,b,c,z);
	if (z<-1) return pow(1-z,-a)*gsl_sf_hyperg_2F1(a,c-b,c,z/(z-1));
	else { cout << "mu2F1 err" << endl; exit(-1); }
}

double H(int n, int m, double y) {
    if (n==1) return 2*pow(y,2*m+1)/(2*m+1);
    if (y==0) return 0;
    if (y==1) return 2./(2*m+1);
    return 2*pow(y,2*m+1)*my2F1(0.5+m, 0.5*(1-n), 1.5+m, 1-y*y)/(2*m+1);
}

/*----------------------------------------------------------------------------------------------------*/
// Damping function
/*----------------------------------------------------------------------------------------------------*/

__device__ double D(int i2, int i1, double *lt4, double *lt) {
    if (i1==i2) return 1;
    double res = 0, w = 1;
    for (int j = i1; j <= i2; j++) {
        if (j==i1 || j==i2) w = 0.5;
        else w = 1.0;
        res += w*pow(lt4[j],0.25)*lt[j];
    }

    res *= DTAU/hbarc/EB/5.;
    return exp(-res);
}

/*----------------------------------------------------------------------------------------------------*/
// Device routines for T^4 iterative computation
/*----------------------------------------------------------------------------------------------------*/

// right hand side for t4 update
__device__ double rhs(int i, double *lt4, double *lt, double *ld, double *lh) {
    double res = 0;
    double w = 1;
    // second term
    if (i>0) {
        for (int ip = 0; ip <= i; ip++) {
            if (ip==0 || ip==i) w = 0.5;
            else w = 1.0;
            res += w*ld[TIDX(i,ip)]*lh[TIDX(i,ip)]*pow(lt4[ip],1.25)*lt[ip];
        }
        res *= DTAU/hbarc/EB/10.;
    }
    // first term
    res += ld[TIDX(i,0)]*pow(T_0,4.)*H(A_0*lt[0]/lt[i])/H(A_0);
    // return result
    return res;
}

// makes one iteration
__global__ void makeIteration(double *lt4, double *lT4, double *lt, double *ld, double *lh) {
    //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid<NUM) {
        lT4[tid] = rhs(tid,lt4,lt,ld,lh);
        tid += blockDim.x * gridDim.x;
    }
}

// load damping function
__global__ void loadDampingFunction(double *lt4, double *lt, double *ld) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < NUM*(NUM+1)/2) {
        int row = floor(-0.5 + sqrt(0.25 + 2 * tid));
        int triangularNumber = row * (row + 1) / 2;
        int column = tid - triangularNumber;
        ld[tid] = D(row,column,lt4,lt);
        tid += blockDim.x * gridDim.x;
    }
}

// load H function
__global__ void loadHFunction(double *lt, double *lh) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < NUM*(NUM+1)/2) {
        int row = floor(-0.5 + sqrt(0.25 + 2 * tid));
        int triangularNumber = row * (row + 1) / 2;
        int column = tid - triangularNumber;
        lh[tid] = H(lt[column]/lt[row]);
        tid += blockDim.x * gridDim.x;
    }
}

/*----------------------------------------------------------------------------------------------------*/
// Device routines for general moment computation
/*----------------------------------------------------------------------------------------------------*/

// right hand side for mnm update
__device__ double rhsMNM(int i, double *lt4, double *lt, double *ld, double *lhnm, double *lhnm0, double lH0) {
    double res = 0;
    double w = 1;
    int r = N + 2*M + 2;
    // second term
    if (i>0) {
        for (int ip = 0; ip <= i; ip++) {
            if (ip==0 || ip==i) w = 0.5;
            else w = 1.0;
            res += w*ld[TIDX(i,ip)]*lhnm[TIDX(i,ip)]*pow(lt4[ip],0.25*(1+r))*lt[ip];
        }
        res *= DTAU/hbarc/EB/5.;
    }
    // first term
    res += pow(2.,0.25*r)*ld[TIDX(i,0)]*pow(T_0,r)*lhnm0[i]/pow(lH0,0.25*r);
    // return result
    return tgamma((double)r)*res/2/2/M_PI/M_PI;
}

// makes one iteration; this is a "kernel"
__global__ void computeMNM(double *lm, double *lt4, double *lt, double *ld, double *lhnm, double *lh, double *lhnm0, double H0) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid<NUM) {
        lm[tid] = rhsMNM(tid,lt4,lt,ld,lhnm,lhnm0,H0);
        tid += blockDim.x * gridDim.x;
    }
}

/*----------------------------------------------------------------------------------------------------*/
// Host routines for calculating f
/*----------------------------------------------------------------------------------------------------*/

// right hand side for f update
__device__ double rhsF(int i, double *lt4, double *lt, double *ld) {
    double res = 0, feq=0, T=1;
    double w = 1;
    // second term
    if (i>0) {
        for (int ip = 0; ip <= i; ip++) {
            if (ip==0 || ip==i) w = 0.5;
            else w = 1.0;
            T = pow(lt4[ip],0.25);
            feq = exp(-sqrt(PZ*PZ+PT*PT));
            res += w*ld[TIDX(i,ip)]*feq*T*lt[ip];
        }
        res *= DTAU/hbarc/EB/5.;
    }
    // first term
    T = pow(lt4[i],0.25);
    double l0 = pow(2./H(A_0),0.25)*T_0;
    double f0 = exp(-sqrt(pow(PZ*lt[i]/(A_0*lt[0]),2) + PT*PT)/(l0/T));
    res += ld[TIDX(i,0)]*f0;
    // return result
    return res;
}

// makes one iteration; this is a "kernel"
__global__ void computeF(double *lf, double *lt4, double *lt, double *ld) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid<NUM) {
        lf[tid] = rhsF(tid,lt4,lt,ld);
        tid += blockDim.x * gridDim.x;
    }
}

/*----------------------------------------------------------------------------------------------------*/
// Host routines for T^4 iterations
/*----------------------------------------------------------------------------------------------------*/

void makeIterations(double *lt4, double *lT4, double *lt, double *ld, double *lh) {

    outputMeasurements(0);
    outputTemperatureSnapshot(t4,0,"T");

    // load H function
    loadHFunction<<<num*(num+1)/2/BLOCKSIZE1,BLOCKSIZE1>>>(lt,lh);
    checkCudaErrors(cudaDeviceSynchronize());

    for (int i=1; i<=maxiters;i++) {

        // load D function
        loadDampingFunction<<<num*(num+1)/2/BLOCKSIZE1,BLOCKSIZE1>>>(lt4,lt,ld);
        cudaDeviceSynchronize();
        
        // make an iteration
        makeIteration<<<num/BLOCKSIZE2,BLOCKSIZE2>>>(lt4,lT4,lt,ld,lh);
        cudaDeviceSynchronize();
        
        // swap pointers to make old <-> new
        swapPointers(&lt4,&lT4);
        
        // output some stuff if appropriate
        if (i%update==0) {
            cudaMemcpy(t4, lt4, num*sizeof(double), cudaMemcpyDeviceToHost);
            outputMeasurements(i);
        }
        if (i%snapupdate==0) {
            cudaMemcpy(t4, lt4, num*sizeof(double), cudaMemcpyDeviceToHost);
            outputTemperatureSnapshot(t4,i,"T");
        }
    }
    // load the device d function based on final result and copy t4 back to host for subsequent use
    loadDampingFunction<<<num*(num+1)/2/BLOCKSIZE1,BLOCKSIZE1>>>(lt4,lt,ld);
    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpy(t4, lt4, num*sizeof(double), cudaMemcpyDeviceToHost));
}

// loads integration abscissae
void loadTimeGrid() {
    cout << "==> Loading time grid" << endl;
    double ltf = log(tf);
    double lt0 = log(t0);
    dt = (ltf-lt0)/(num-1);
    for (int i = 0; i < num; i++) t[i] = exp(lt0 + i*dt);
}

// initializes t4 array
void initializeT4() {
    cout << "==> Initializing T^4 array" << endl;
    t4[0] = T0*T0*T0*T0;
    for (int i=1; i < num; i++) {
        t4[i] = T0*T0*T0*T0*pow(t0/t[i],4./3.);
    }
}

/*----------------------------------------------------------------------------------------------------*/
// Host routines for general moment computation
/*----------------------------------------------------------------------------------------------------*/

// initializes hnm array
void setupHNM(int n, int m) {
    for (int idx=0; idx < num*(num+1)/2; idx++) {
        int row = floor(-0.5 + sqrt(0.25 + 2 * idx));
        int triangularNumber = row * (row + 1) / 2;
        int column = idx - triangularNumber;
        hnm[idx] = H(n,m,t[column]/t[row]);
    }
    for (int idx=0; idx < num; idx++)
        hnm0[idx] = H(n,m,t[0]*a0/t[idx]);
}

// computes a general moment based on the current iterations results for t4
double* computeMoment(int n, int m) {
    cout << "==> Computing M(" << n << "," << m << ")" << endl;

    cudaMemcpyToSymbol(&N, &n, sizeof(int));
    cudaMemcpyToSymbol(&M, &m, sizeof(int));
    
    setupHNM(n,m);
    cudaMemcpy(dev_hnm, hnm, sizeof(double)*num*(num+1)/2, cudaMemcpyHostToDevice); // transfer to device
    cudaMemcpy(dev_hnm0, hnm0, sizeof(double)*num, cudaMemcpyHostToDevice); // transfer to device
    
    double *lm;
    lm = allocate1DArray();
    computeMNM<<<num/BLOCKSIZE2,BLOCKSIZE2>>>(dev_m, dev_t4, dev_time, dev_d, dev_hnm, dev_h, dev_hnm0, hostH(a0));
    cudaMemcpy(lm, dev_m, num*sizeof(double), cudaMemcpyDeviceToHost);
    
    return lm;
}

// computes a general moment based on an equilbrium form with t4
inline double computeEQMoment(int n, int m, int i) {
    int r = n+2*m+2;
    return gsl_sf_gamma(r)*pow(t4[i],0.25*r)*2/(2*m+1)/2/2/M_PI/M_PI;
}

/*----------------------------------------------------------------------------------------------------*/
// Host routines for f computation
/*----------------------------------------------------------------------------------------------------*/

// computes f based on the current iterations results for t4
double* computeDistributionFunction(double pz, double pt) {
    //cout << "==> Computing f(" << pz << "," << pt << ")" << endl;
    
    cudaMemcpyToSymbol(&PZ, &pz, sizeof(double));
    cudaMemcpyToSymbol(&PT, &pt, sizeof(double));

    double *lf;
    lf = allocate1DArray();
    computeF<<<num/BLOCKSIZE2,BLOCKSIZE2>>>(dev_f, dev_t4, dev_time, dev_d);
    cudaMemcpy(lf, dev_f, num*sizeof(double), cudaMemcpyDeviceToHost);
    
    return lf;
}

/*----------------------------------------------------------------------------------------------------*/
// Main routine
/*----------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
    const double m_pi = 4.0 * atan(1.0);

    char fname[20]; // for later use
    
    print_line();
    
    // read parameters from file and command line
    readParametersFromFile("params.txt",1);
    if (argc>1) {
        print_line();
        cout << "Parameters from commandline" << endl;
        print_line();
        readParametersFromCommandLine(argc,argv,1);
    }
    // perform any processing of parameters necessary
    processParameters();
    
    print_line();
    print_line();
    print_line();
    
    //setup
    allocateMemory();
    loadTimeGrid();
    initializeT4();
        
    print_line();
    // copy grid and initial conditions to device
    checkCudaErrors(cudaMemcpy(dev_t4, t4, sizeof(double)*num, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_T4, T4, sizeof(double)*num, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_time, t, sizeof(double)*num, cudaMemcpyHostToDevice));

    print_line();
    // copy parameters to device memory
    checkCudaErrors(cudaMemcpyToSymbol(NUM, &num, sizeof(int)));
    checkCudaErrors(cudaMemcpyToSymbol(A_0, &a0, sizeof(double)));
    checkCudaErrors(cudaMemcpyToSymbol(T_0, &T0, sizeof(double)));
    double eb = fpieb/m_pi/4.;
    checkCudaErrors(cudaMemcpyToSymbol(EB, &eb, sizeof(double)));
    checkCudaErrors(cudaMemcpyToSymbol(DTAU, &dt, sizeof(double)));

    print_line();
    // Copy to grid global constants
    checkCudaErrors(cudaMemcpyToSymbol(M_PI, &m_pi, sizeof(double)));
    checkCudaErrors(cudaMemcpyToSymbol(hbarc, &HBARC, sizeof(double)));
    
    // print some stuff
    print_line();
    cout.width(dwidth); cout << "iteration";
    cout.width(dwidth); cout << "T[0]";
    cout.width(dwidth); cout << "T[num/2]";
    cout.width(dwidth); cout << "T[num-1]";
    cout << endl;
    print_line();
    
    /*----------------------------------------------------------------------------------------------------*/
    // Iterations
    /*----------------------------------------------------------------------------------------------------*/

    cout << "4 pi eta / S: " << fpieb << endl;
    makeIterations(dev_t4, dev_T4, dev_time, dev_d, dev_h);
    
    /*----------------------------------------------------------------------------------------------------*/
    // Compute some things with the solution
    /*----------------------------------------------------------------------------------------------------*/

    print_line();
    double *ed, *pl,*pt,*plopt;
    ed = allocate1DArray();
    pl = computeMoment(0,1);
    pt = allocate1DArray();
    plopt = allocate1DArray();
    for (int i=0; i<num; i++) {
        ed[i] = 3*t4[i]/m_pi/m_pi;
        pt[i] = 0.5*(ed[i] - pl[i]);
        plopt[i] = pl[i]/pt[i];
        cout << ed[i] << "\t" << pl[i] << "\t" << pt[i] << endl;
    }

    outputArray(ed,"ed");
    outputArray(pl,"pl");
    outputArray(pt,"pt");
    outputArray(plopt,"pratio");


    // compute distribution function
    if (computeDist==1) {
        print_line();
        cout << "==> Computing f ";
        double ***f3DArray;
        f3DArray = allocate3DArray(num/fStep,numPZ,numPT);
        double dpz = maxPZ/(numPZ-1);
        double dpt = maxPT/(numPT-1);
        for (int i=0; i<numPZ; i++) {
            for (int j=0; j<numPT; j++) {
                double *f;
                f = computeDistributionFunction(i*dpz,j*dpt);
                for (int k=0; k<num/fStep; k++) f3DArray[k][i][j] = f[k*fStep]; // load into f array for later binary output
                free1DArray(f);
            }
            cout << "." << std::flush;
        }
        cout << endl;
        outputDistribution(f3DArray); // output f in binary format
        free3DArray(f3DArray,num/fStep,numPZ,numPT);
    }
    
    if (computeMoments==1) {
        // loop over moments
        print_line();
        for (int n=0; n<=maxN; n++) {
            for (int m=0; m<=maxM; m++) {
                double *mom;
                mom = computeMoment(n,m);
                sprintf(fname,"moms/m-%d-%d",n,m);
                outputArray(mom,fname);
                for (int i=0; i<num; i++) mom[i] /= computeEQMoment(n,m,i);
                sprintf(fname,"moms/m-%d-%d-scaled",n,m);
                outputScaledArray(mom,t4,5*eb,fname);
                free1DArray(mom);
            }
        }
    }
    
    free1DArray(ed);
    free1DArray(pl);
    free1DArray(pt);
    free1DArray(plopt);

    /*----------------------------------------------------------------------------------------------------*/

    // print some more stuff
    print_line();
    cout << "Done.\n";
    print_line();
    
    // free memory
    freeMemory();
    
    return 0;
}
