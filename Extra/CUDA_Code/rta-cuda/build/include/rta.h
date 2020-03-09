/*
 
 dse.h
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
*/

#ifndef __dse_h__
#define __dse_h__

extern int num, maxiters, update, snapupdate, computeMoments, computeDist, maxN, maxM;
extern double fpieb, t0, tf, T0, a0;

extern double *t4, *T4, *t, *hnm, *hnm0, *f;
extern double *dev_t4, *dev_T4, *dev_time, *dev_d, *dev_h, *dev_hnm, *dev_hnm0, *dev_m, *dev_f;

extern int numPZ,numPT,fStep;
extern double maxPZ, maxPT;

const double HBARC = 0.197326938; // GeV fm

#endif /* __dse_h__ */
