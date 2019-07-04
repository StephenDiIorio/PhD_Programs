/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  Function list + headers
  --------------------------------------------------
*/

#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "constants.h"

int initpos(double **pos);
int initmom(double **mom);
int initEfield(double ****Efield);
int initEfieldFromFile(double ****Efield, const char *xFile, const char *yFile, const char *zFile);
int initBfield(double ****Bfield);

int WeightF(double ****EVec, double **parVec, double **FVec);

int pushmomentum(double **FVec, double **momVec, const double dt_m);
int pushxpos(double **momVec, double **parVec, const double dt_m);
int Lorentz(double **momVec, double **BVec, const double dt_m);
int singlepushxpos(double **momVec, double **parVec, const double dist);
