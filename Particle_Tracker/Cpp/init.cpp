/*
   3D particle tracking code
  AGRT 2010

  Initialization functions for PARTRAC
*/

#ifndef PARTRAC_INIT
#define PARTRAC_INIT
#include "constants.h"

//Initialize particle positions
int initpos(double **pos) {
  unsigned long i, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three position components
    for (i = 0; i < Npar; ++i) {
      pos[x1][i] = 0.0;
    }
  }
  // pos[0][0] = 2.0;
  // pos[1][0] = 2.0;
  // pos[2][0] = 2.0;

  return 0;
}

//Initialize particle momenta
int initmom(double **mom) {
  unsigned long i, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three position components
    for (i = 0; i < Npar; ++i) {
      mom[x1][i] = 0.0;
    }
    // mom[0][0] = 1.0;
  }

  return 0;
}

//Initialize Efield
int initEfield(double ****Efield) {
  unsigned long i, j, k, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three position components
    for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) {
        for (k = 0; k < Nz; ++k) {
          Efield[x1][i][j][k] = (1.0 * k - 2.0) * 0.0;
        }
      }
    }
  }

  return 0;
}

//Initialize Bfield
int initBfield(double ****Bfield) {
  unsigned long i, j, k, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three position components
    for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) {
        for (k = 0; k < Nz; ++k) {
          Bfield[x1][i][j][k] = 0.0;
          Bfield[2][i][j][k] = 1.0;
        }
      }
    }
  }

  return 0;
}

#endif
