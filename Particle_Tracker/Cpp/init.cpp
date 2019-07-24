/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  Initialization functions for PARTRAC
  --------------------------------------------------
*/

#ifndef PARTRAC_INIT
#define PARTRAC_INIT
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "constants.h"

/**
 * @brief Initialize particle positions
 *
 * @param pos
 * @return int
 */
int initpos(double **pos) {
    unsigned long i, x1;

    const double centerX = 0.0;
    const double centerY = Ny * dx[yIndx] / 2.;
    const double centerZ = Nz * dx[zIndx] / 2.;

    // Spread particles in 40e-6 range from 0 to cover domain
    const double spreadY = 31.81359646532949;
    const double spreadZ = 0.0;

    for (i = 0; i < Npar; ++i) {
        pos[xIndx][i] = centerX;
        pos[yIndx][i] = (2 * centerY) *
            ((double)rand() / (RAND_MAX)) + (centerY - spreadY);
        pos[zIndx][i] = centerZ;//(2 * centerZ) * ((double)rand() / (RAND_MAX)) + (centerZ - spreadZ);
    }

    return 0;
}

/**
 * @brief Initialize particle momenta
 *
 * @param mom
 * @return int
 */
int initmom(double **mom) {
    unsigned long i, x1;

    const double lowerMomTot = 0.6410483775684676;
    const double upperMomTot = 0.8035966691084001;

    const double lowerPy = -0.00242294 * 2.; // init to 5 eV
    const double upperPy = 0.00242294 * 2.;
    const double lowerPz = 0.0; // -0.00242294;
    const double upperPz = 0.0; // 0.00242294;

    double px, py, pz, momenta;

    for (i = 0; i < Npar; ++i) {
        momenta = (upperMomTot - lowerMomTot) *
            ((double)rand() / (RAND_MAX)) +
            lowerMomTot;
        py = (upperPy - lowerPy) *
            ((double)rand() / (RAND_MAX)) + lowerPy;
        pz = (upperPz - lowerPz) *
            ((double)rand() / (RAND_MAX)) + lowerPz;

        px = sqrt(pow(momenta, 2.) - pow(py, 2.) - pow(pz, 2.));

        mom[xIndx][i] = px;
        mom[yIndx][i] = py;
        mom[zIndx][i] = pz;
    }

    return 0;
}

/**
 * @brief Initialize Efield
 *
 * @param Efield
 * @return int
 */
int initEfield(double ****Efield) {
    unsigned long i, j, k, x1;

    for (x1 = 0; x1 < ndims; ++x1) { // three position components
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

/**
 * @brief Initialize Efield from a text file
 *
 * @param Efield
 * @return int
 */
int initEfieldFromFile(double ****Efield, const char *xFile, const char *yFile, const char *zFile) {
    unsigned long i, j, k, x1;

    std::ifstream xinputFile(xFile);//"xfield.dat"); // Input file stream object
    std::ifstream yinputFile(yFile);//"yfield.dat");
    std::ifstream zinputFile(zFile);//"zfield.dat");

    // Check if exists and then open the file.
    if (xinputFile.good() && yinputFile.good() && zinputFile.good()) {
        double xVal = 0., yVal = 0., zVal = 0.;
        double n0_const = sqrt(n0 * 1.e-6);
        double ef_n0_coeff = 9.613e-10;

        for (i = 0; i < Nx; ++i) {
            for (j = 0; j < Ny; ++j) {
                xinputFile >> xVal;
                yinputFile >> yVal;
                // zinputFile >> zVal;
                for (k = 0; k < Nz; ++k) {
                    Efield[xIndx][i][j][k] = xVal / 1.e11 /
                        ef_n0_coeff / n0_const; // convert to osiris units

                    Efield[yIndx][i][j][k] = yVal / 1.e11 /
                        ef_n0_coeff / n0_const; // convert to osiris units

                    Efield[zIndx][i][j][k] = zVal / 1.e11 /
                        ef_n0_coeff / n0_const; // convert to osiris units
                }
            }
        }

        // Close the file.
        xinputFile.close();
        yinputFile.close();
        zinputFile.close();

        // Display the numbers read:
        // printf("The numbers are:\n");
        // for (i = 0; i < Nx; ++i) {
        //     for (j = 0; j < Ny; ++j) {
        //         for (k = 0; k < Nz; ++k) {
        //             printf("%f ", Efield[xIndx][i][j][k]);
        //         }
        //     }
        //     printf("\n");
        // }
        // printf("\n");
    } else {
        printf("Error reading file!\n");
        exit(0);
    }

    return 0;
}

/**
 * @brief Initialize Bfield
 *
 * @param Bfield
 * @return int
 */
int initBfield(double ****Bfield) {
    unsigned long i, j, k, x1;

    for (x1 = 0; x1 < ndims; ++x1) { // three position components
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j < Ny; ++j) {
                for (k = 0; k < Nz; ++k) {
                    Bfield[x1][i][j][k] = 0.0;
                    // Bfield[2][i][j][k] = 1.0;
                }
            }
        }
    }

  return 0;
}

#endif
