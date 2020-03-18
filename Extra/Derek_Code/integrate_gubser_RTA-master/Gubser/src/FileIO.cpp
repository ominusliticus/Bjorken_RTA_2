#pragma once
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameter.h"

void readInParameters(struct parameters &params)
{
  char dummyChar[255];
  int dummyInt;
  float dummyFloat;

  FILE *fileIn;
  std::stringstream paramsStream;
  paramsStream << "parameters.dat";
  fileIn = fopen(paramsStream.str().c_str(),"r");

  if (fileIn == NULL)
  {
    printf("Couldn't open parameters.dat\n");
  }

  else
  {
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_rho = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.c = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.rho0 = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.T_hat0 = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_grid_rho = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.rho_min = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.delta_rho = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_grid_p_omega_hat = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.p_omega_hat_min = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.delta_p_omega_hat = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_grid_p_eta_hat = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.p_eta_hat_min = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.delta_p_eta_hat = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_r = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.delta_r = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_p_mag = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.delta_p_mag = dummyFloat;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_phip = dummyInt;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.n_vz = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.tau0 = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.q0 = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.resolution = dummyFloat;

    fclose(fileIn);
  }
}
