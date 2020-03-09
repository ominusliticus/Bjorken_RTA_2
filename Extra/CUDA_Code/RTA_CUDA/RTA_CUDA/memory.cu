/*
 
 memory.cpp
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <cuda_runtime.h>

#include <helper_cuda.h>
#include <helper_functions.h>

using namespace std;

#include "rta.h"

double* allocate1DArray() {
    double *tmp;
    tmp = new double[num];
    return tmp;
}

// allocate flattened triangle array
double* allocateFTArray() {
    double *tmp;
    tmp = new double[num*(num+1)];
    return tmp;
}

double*** allocate3DArray(int n1, int n2, int n3) {
    double ***tmp;
    tmp = new double**[n1];
    for (int sx=0;sx<n1;sx++) tmp[sx] = new double*[n2];
    for (int sx=0;sx<n1;sx++) for (int sy=0;sy<n2;sy++) tmp[sx][sy] = new double[n3];
    return tmp;
}

void free1DArray(double *array) {
    delete[] array;
}

void free3DArray(double*** array, int n1, int n2, int n3) {
    for (int sx=0;sx<n1;sx++) for (int sy=0;sy<n2;sy++) free(array[sx][sy]);
    for (int sx=0;sx<n1;sx++) free(array[sx]);
    free(array);
    return;
}

void allocateMemory()
{
    cout << "==> Allocating memory\n";
    
    // allocate host memory
    t4 = allocate1DArray();
    T4 = allocate1DArray();
    t = allocate1DArray();
    hnm = allocateFTArray();
    hnm0 = allocate1DArray();
    f = allocate1DArray();
    
    // allocate device memory (GPU)
    checkCudaErrors(cudaMalloc((void**)&dev_t4, sizeof(double)*num));
    checkCudaErrors(cudaMalloc((void**)&dev_T4, sizeof(double)*num));
    checkCudaErrors(cudaMalloc((void**)&dev_time, sizeof(double)*num));
    checkCudaErrors(cudaMalloc((void**)&dev_d, sizeof(double)*num*(num+1)/2));
    checkCudaErrors(cudaMalloc((void**)&dev_h, sizeof(double)*num*(num+1)/2));
    checkCudaErrors(cudaMalloc((void**)&dev_hnm, sizeof(double)*num*(num+1)/2));
    checkCudaErrors(cudaMalloc((void**)&dev_hnm0, sizeof(double)*num));
    checkCudaErrors(cudaMalloc((void**)&dev_m, sizeof(double)*num*(num+1)/2));
    checkCudaErrors(cudaMalloc((void**)&dev_f, sizeof(double)*num));
}

void freeMemory()
{
    // free host memory
    free1DArray(t4);
    free1DArray(T4);
    free1DArray(t);
    free1DArray(hnm);
    free1DArray(hnm0);
    free1DArray(f);
    
    // free device memory
    cudaFree( dev_t4 );
    cudaFree( dev_T4 );
    cudaFree( dev_time );
    cudaFree( dev_d );
    cudaFree( dev_h );
    cudaFree( dev_hnm );
    cudaFree( dev_hnm0 );
    cudaFree( dev_m );
    cudaFree( dev_f );
}

void swapPointers(double **pt1, double **pt2) {
    double *tmp = *pt1;
    *pt1 = *pt2;
    *pt2 = tmp;
}
