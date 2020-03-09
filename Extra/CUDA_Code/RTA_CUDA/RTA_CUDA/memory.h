/*
 
 memory.h
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
*/

#ifndef __memory_h__
#define __memory_h__

void swapPointers(double **pt1, double **pt2);
void allocateMemory();
void freeMemory();
double* allocate1DArray();
double*** allocate3DArray(int n1, int n2, int n3);
void free1DArray(double *m);
void free3DArray(double ***m, int n1, int n2, int n3);

#endif /* __memory_h__ */
