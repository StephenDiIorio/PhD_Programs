/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  Equations of motion

  p^n+1/2 = p^n-1/2 + F^n*dt
  x^n+1 = x^n + p^n+1/2 *dt/sqrt(1+p^2)
  --------------------------------------------------
*/

#ifndef PARTRAC_MOTION
#define PARTRAC_MOTION
#include "constants.h"

/**
 * @brief Electric push
 *
 * @param FVec
 * @param momVec
 * @param dt_m
 */
int pushmomentum(double **FVec, double **momVec, const double dt_m) {
    unsigned long i, x1;

    for (x1 = 0; x1 < ndims; ++x1) { // three momentum components
        for (i = 0; i < Npar; ++i) {
            momVec[x1][i] += FVec[x1][i] * dt_m;
        }
    }

    return 0;
}

/**
 * @brief
 *
 * @param momVec
 * @param parVec
 * @param dt_m
 */
int pushxpos(double **momVec, double **parVec, const double dt_m) {
    double psquared, igamma;
    unsigned long i, x1;

    for (x1 = 0; x1 < ndims; ++x1) { // three position components
        for (i = 0; i < Npar; ++i) {
            psquared = momVec[xIndx][i] * momVec[xIndx][i] + momVec[yIndx][i] * momVec[yIndx][i] + momVec[zIndx][i] * momVec[zIndx][i];
            igamma = 1.0 / sqrt(1.0 + psquared);
            parVec[x1][i] += momVec[x1][i] * dt_m * igamma;
        }
    }

    return 0;
}

/**
 * @brief Magnetic rotation
 *
 * @param momVec
 * @param BVec
 * @param dt_m
 */
int Lorentz(double *momVec[], double *BVec[], const double dt_m) {
    long unsigned parnum, x1;
    double igamma, psquared, Bsquared;

    // t vector in Boris method
    double tt[ndims], ttsquared;
    // s vector in Boris method
    double ss[ndims];
    // vstar vector in Boris method
    double vstar[ndims];
    // perpendicular component of v in Boris method
    double vperp[ndims];

    for (parnum = 0; parnum < Npar; ++parnum) {
        Bsquared = (BVec[xIndx][parnum] * BVec[xIndx][parnum] + BVec[yIndx][parnum] * BVec[yIndx][parnum] + BVec[zIndx][parnum] * BVec[zIndx][parnum]);

        if (Bsquared) {
            psquared = momVec[xIndx][parnum] * momVec[xIndx][parnum] + momVec[yIndx][parnum] * momVec[yIndx][parnum] + momVec[zIndx][parnum] * momVec[zIndx][parnum];
            igamma = 1.0 / sqrt(1.0 + psquared);
            ttsquared = 0.0;
            for (x1 = 0; x1 < ndims; ++x1) {
                tt[x1] = BVec[x1][parnum] * dt_m * 0.5;
                ttsquared += tt[x1] * tt[x1];
            }

            for (x1 = 0; x1 < ndims; ++x1) {
                ss[x1] = 2.0 * tt[x1] / (1.0 + ttsquared);
            }

            // calculated vperp
            for (x1 = 0; x1 < ndims; ++x1) {
                vperp[x1] = momVec[x1][parnum] * (1.0 - BVec[x1][parnum] / sqrt(Bsquared));
            }

            //calculate vstar
            // for efficiency, component by component
            vstar[xIndx] = vperp[xIndx] + (vperp[yIndx] * tt[zIndx] - vperp[zIndx] * tt[yIndx]) * igamma;
            vstar[yIndx] = vperp[yIndx] + (vperp[zIndx] * tt[xIndx] - vperp[xIndx] * tt[zIndx]) * igamma;
            vstar[zIndx] = vperp[zIndx] + (vperp[xIndx] * tt[yIndx] - vperp[yIndx] * tt[xIndx]) * igamma;
            //Finally update momentum
            momVec[xIndx][parnum] += (vstar[yIndx] * ss[zIndx] - vstar[zIndx] * ss[yIndx]) * igamma;
            momVec[yIndx][parnum] += (vstar[zIndx] * ss[xIndx] - vstar[xIndx] * ss[zIndx]) * igamma;
            momVec[zIndx][parnum] += (vstar[xIndx] * ss[yIndx] - vstar[yIndx] * ss[xIndx]) * igamma;
    }
  }

  return 0;
}

int singlepushxpos(double **momVec, double **parVec, const double dist) {
    unsigned long i;
    double t;

    for (i = 0; i < Npar; ++i) {
        t = (dist - parVec[xIndx][i]) / momVec[xIndx][i];

        parVec[xIndx][i] += momVec[xIndx][i] * t;
        parVec[yIndx][i] += momVec[yIndx][i] * t;
        parVec[zIndx][i] += momVec[zIndx][i] * t;
    }

    return 0;
}

#endif
