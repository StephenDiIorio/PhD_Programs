/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  ********* PARTRAC *****************
  This is a  relativisitic  particle tracking code
  for charged particles in the presence of gridded
  fields including a second order weighting

  The finite difference equations use Boris scheme

  The fields are on a large 3D array in a Cartesian
  geometry
  --------------------------------------------------
*/

// User included
#include "partrac.h"

int main(int argnum, char **argstr) {
    unsigned long i, j, x1, f;
    const double tmax = 125; //sqrt(2.0) * 2.0 * PI;
    double t = 0;
    const double dt = 0.01;
    const double dist = 143002.;
    bool outsideDomain;

    double start_time = clock();

    //----------------------------------------
    // Front material
    printf("\n---- Partrac v");
    printf(version);
    printf(" ----\n\n");

    //----------------------------------------
    // For File Output
    FILE *trajOutFile = fopen(trajFilename, "w");
    FILE *sliceOutFile = fopen(sliceFilename, "w");

    //----------------------------------------
    // Arrays and counters for data dumps
    unsigned long t_count;
    static double *t_array;
    t_array = new double[(unsigned long)(tmax / dt)];

    static double ***pos_data;
    pos_data = new double **[(unsigned long)(tmax / dt)];
    for (i = 0; i < (unsigned long)(tmax / dt); ++i) {
        pos_data[i] = new double *[Npar];
        for (j = 0; j < Npar; ++j) {
            pos_data[i][j] = new double [ndims];
        }
    }

    static double **finalYPos;
    finalYPos = new double *[numFiles];
    for (f = 0; f < numFiles; ++f) {
        finalYPos[f] = new double [Npar];
    }

    //----------------------------------------
    // Particle positions and momenta
    static double **pos; //position
    static double **mom; //momentum
    static double **fce; //electric force at particle
    static double **bpa; //magnetic force at particle
    pos = new double *[ndims];
    mom = new double *[ndims];
    fce = new double *[ndims];
    bpa = new double *[ndims];

    for (x1 = 0; x1 < ndims; ++x1) {
        pos[x1] = new double [Npar];
        mom[x1] = new double [Npar];
        fce[x1] = new double [Npar];
        bpa[x1] = new double [Npar];
    }

    // Field cubes - 3 dimensional (the '3' here refers to Ex, Ey, Ez at every coordinate space)
    static double ****Efield;
    static double ****Bfield;
    Efield = new double ***[ndims];
    Bfield = new double ***[ndims];
    for (x1 = 0; x1 < ndims; ++x1) {
        Efield[x1] = new double **[Nx];
        Bfield[x1] = new double **[Nx];
        for (i = 0; i < Nx; ++i) {
            Efield[x1][i] = new double *[Ny];
            Bfield[x1][i] = new double *[Ny];
            for (j = 0; j < Ny; ++j) {
                Efield[x1][i][j] = new double [Nz];
                Bfield[x1][i][j] = new double [Nz];
            }
        }
    }

    // NOTE: Even if input files are same, will get different results
    // becuase of rand() in pos initialization. Mitigate this by calling
    // srand(1) in initpos.
    for (f = 0; f < numFiles; ++f) {
        sprintf(xInput, "x%lu_field.dat", f);
        sprintf(yInput, "y%lu_field.dat", f);
        sprintf(zInput, "z%lu_field.dat", f);

        // initialize particles
        initpos(pos);
        initmom(mom);

        //initialize grids
        initEfieldFromFile(Efield, xInput, yInput, zInput);
        initBfield(Bfield);
        outsideDomain = WeightF(Efield, pos, fce);
        outsideDomain = WeightF(Bfield, pos, bpa);

        t = 0;
        t_count = 0;

        //----------------------------------------
        // Main loop
        while (t<=tmax) {
            //Weight Electric fields to particles
            outsideDomain = WeightF(Efield, pos, fce);
            outsideDomain = WeightF(Bfield, pos, bpa);

            if (outsideDomain) {
                printf("Breaking time loop. All particles are out of box. t=%f\n", t);
                break;
            }

            // Equations of motion
            pushmomentum(fce, mom, dt * 0.5);
            Lorentz(mom, bpa, dt);
            pushmomentum(fce, mom, dt * 0.5);
            pushxpos(mom, pos, dt);

            // Store data
            t_array[t_count] = t;
            for (i = 0; i < Npar; ++i) {
                for (x1 = 0; x1 < ndims; ++x1) {
                    pos_data[t_count][i][x1] = pos[x1][i];
                }
            }
            ++t_count;

            t += dt;
        }

        printf("Were all the particles outside the domain?: %i (0=No, 1=Yes)\n", outsideDomain);

        // Final push to target distance and data dump
        singlepushxpos(mom, pos, dist);
        t_array[t_count] = t;
        for (i = 0; i < Npar; ++i) {
            for (x1 = 0; x1 < ndims; ++x1) {
                pos_data[t_count][i][x1] = pos[x1][i];
            }
        }
        ++t_count;

        for (i = 0; i < Npar; ++i) {
            finalYPos[f][i] = pos[yIndx][i];
        }
    }

    //----------------------------------------
    // Write data to file
    for (i = 0; i < t_count; ++i) {
        fprintf(trajOutFile, "%e", t_array[i]);
        for (j = 0; j < Npar; ++j) {
            fprintf(trajOutFile, "\t%e\t%e\t%e", pos_data[i][j][xIndx], pos_data[i][j][yIndx], pos_data[i][j][zIndx]);
        }
        fprintf(trajOutFile, "\n");
    }

    for (f = 0; f < numFiles; ++f) {
        for (i = 0; i < Npar; ++i) {
            fprintf(sliceOutFile, "%f\t", finalYPos[f][i]);
        }
        fprintf(sliceOutFile, "\n");
    }

    //----------------------------------------
    // Finish
    fclose(trajOutFile);
    fclose(sliceOutFile);
    double end_time = clock();
    printf("finished\n");
    printf("\nEND OF CALCULATION\n");
    printf("Calculation took %f\n", (end_time - start_time) / CLOCKS_PER_SEC);

    //----------------------------------------
    // Clean up
    for (i = 0; i < (uint) (tmax / dt); ++i) {
        for (j = 0; j < Npar; ++j) {
            delete[] pos_data[i][j];
        }
        delete[] pos_data[i];
    }
    delete[] pos_data;
    delete[] t_array;

    for (f = 0; f < numFiles; ++f) {
        delete[] finalYPos[f];
    }
    delete[] finalYPos;

    for (x1 = 0; x1 < ndims; ++x1) {
        for (i = 0; i < Nx; ++i) {
            for (j = 0; j < Ny; ++j) {
                delete[] Efield[x1][i][j];
                delete[] Bfield[x1][i][j];
            }
            delete[] Efield[x1][i];
            delete[] Bfield[x1][i];
        }
        delete[] Efield[x1];
        delete[] Bfield[x1];

        delete[] pos[x1];
        delete[] mom[x1];
        delete[] fce[x1];
        delete[] bpa[x1];
    }
    delete[] Efield;
    delete[] Bfield;

    delete[] pos;
    delete[] mom;
    delete[] fce;
    delete[] bpa;

    return 0;
}
