/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  Weighting routines
  conditions
  --------------------------------------------------
*/

#ifndef PARTRAC_WEIGHTING
#define PARTRAC_WEIGHTING
#include "constants.h"

/**
 * @brief Function for weighting electric or magnetic forces to particle position
 *
 * @param EVec Field matrix (E or B) to interpolate
 * @param parVec matrix containing the position of particles
 * @param FVec Interpolated force matrix to fill
 * @return int Whether or not all of the particles are outside the electric or magnetic field array
 */
int WeightF(double ****EVec, double **parVec, double **FVec) {
    bool outsideDomain = true;
    unsigned long parnum, x1, x2, l, m, n, ii[3];
    double weight[3];
    double volumes[2][2][2]; // for volume weighting

    for (x1 = 0; x1 < ndims; ++x1) {
        for (parnum = 0; parnum < Npar; ++parnum) {
            FVec[x1][parnum] = 0.0; //reset F

            // First get nearest gridpoints
            for (x2 = 0; x2 < ndims; ++x2) {
                weight[x2] = parVec[x2][parnum];

                weight[x2] *= idx[x2];
                ii[x2] = (int)(weight[x2]); // cast float -> int

                weight[x2] -= (double)ii[x2]; // now weight is difference between xi and xj
            }

            // Now get the volumes
            for (l = 0; l < 2; ++l) {
                for (m = 0; m < 2; ++m) {
                    for (n = 0; n < 2; ++n) {
                        volumes[l][m][n] = fabs(((1.0 - l) - weight[0]) * ((1.0 - m) - weight[1]) * ((1.0 - n) - weight[2]));
                    }
                }
            }

            for (l = 0; l < 2; ++l) {
                for (m = 0; m < 2; ++m) {
                    for (n = 0; n < 2; ++n) {
                        if (((ii[0] + l) < Nx) && ((ii[1] + m) < Ny) && ((ii[2] + n) < Nz) && ((ii[0] + l) >= 0) && ((ii[1] + m) >= 0) && ((ii[2] + n) >= 0)) {
                            FVec[x1][parnum] += volumes[l][m][n] * EVec[x1][ii[0] + l][ii[1] + m][ii[2] + n];
                            outsideDomain = false;
                        }
                    }
                }
            }

            FVec[x1][parnum] *= qoverm;
        }
    }

    return outsideDomain;
}

#endif
