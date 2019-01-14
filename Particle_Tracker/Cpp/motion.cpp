/*
  --------------------------------------------------
 3D particle tracking code
  AGRT 2010

  Equations of motion

  p^n+1/2 = p^n-1/2 + F^n*dt
  x^n+1 = x^n + p^n+1/2 *dt/sqrt(1+p^2)
  --------------------------------------------------
*/

#ifndef PARTRAC_MOTION
#define PARTRAC_MOTION
#include "constants.h"
#include "iostream"
// dt_m allows for 1/2 step

//Electric push
int pushmomentum(double **FVec, double **momVec, const double dt_m) {
  unsigned long i, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three momentum components
    for (i = 0; i < Npar; ++i) {
      momVec[x1][i] += FVec[x1][i] * dt_m;
    }
  }

  return 0;
}

int pushxpos(double **momVec, double **parVec, const double dt_m) {
  double psquared, igamma;
  unsigned long i, x1;

  for (x1 = 0; x1 < 3; ++x1) { // three position components
    for (i = 0; i < Npar; ++i) {
      psquared = momVec[0][i] * momVec[0][i] + momVec[1][i] * momVec[1][i] + momVec[2][i] * momVec[2][i];
      igamma = 1.0 / sqrt(1.0 + psquared);
      parVec[x1][i] += momVec[x1][i] * dt_m * igamma;
    }
  }

  return 0;
}

//Magnetic rotation
int Lorentz(double *momVec[], double *BVec[], const double dt_m) {
  long unsigned parnum, x1;
  double igamma, psquared, Bsquared;

  // t vector in Boris method
  double tt[3], ttsquared;
  // s vector in Boris method
  double ss[3];
  // vstar vector in Boris method
  double vstar[3];
  // perpendicular component of v in Boris method
  double vperp[3];

  for (parnum = 0; parnum < Npar; ++parnum) {
    Bsquared = (BVec[0][parnum] * BVec[0][parnum] + BVec[1][parnum] * BVec[1][parnum] + BVec[2][parnum] * BVec[2][parnum]);

    if (Bsquared) {
      psquared = momVec[0][parnum] * momVec[0][parnum] + momVec[1][parnum] * momVec[1][parnum] + momVec[2][parnum] * momVec[2][parnum];
      igamma = 1.0 / sqrt(1.0 + psquared);
      ttsquared = 0.0;
      for (x1 = 0; x1 < 3; ++x1) {
        tt[x1] = BVec[x1][parnum] * dt_m * 0.5;
        ttsquared += tt[x1] * tt[x1];
      }

      for (x1 = 0; x1 < 3; ++x1) {
        ss[x1] = 2.0 * tt[x1] / (1.0 + ttsquared);
      }

      // calculated vperp
      for (x1 = 0; x1 < 3; ++x1) {
        vperp[x1] = momVec[x1][parnum] * (1.0 - BVec[x1][parnum] / sqrt(Bsquared));
      }

      //calculate vstar
      // for efficiency, component by component
      vstar[0] = vperp[0] + (vperp[1] * tt[2] - vperp[2] * tt[1]) * igamma;
      vstar[1] = vperp[1] + (vperp[2] * tt[0] - vperp[0] * tt[2]) * igamma;
      vstar[2] = vperp[2] + (vperp[0] * tt[1] - vperp[1] * tt[0]) * igamma;
      //Finally update momentum
      momVec[0][parnum] += (vstar[1] * ss[2] - vstar[2] * ss[1]) * igamma;
      momVec[1][parnum] += (vstar[2] * ss[0] - vstar[0] * ss[2]) * igamma;
      momVec[2][parnum] += (vstar[0] * ss[1] - vstar[1] * ss[0]) * igamma;
    }
  }

  return 0;
}

#endif
