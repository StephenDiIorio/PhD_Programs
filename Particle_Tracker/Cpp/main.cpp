/*
  ********* PARTRAC *****************
  This is a  relativisitic  particle tracking code
  for charged particles in the presence of gridded
  fields including a second order weighting

  The finite difference equations use Boris scheme

  The fields are on a large 3D array in a Cartesian
  geometry
*/

// User included
#include "partrac.h"

int main(int argnum, char **argstr) {
  unsigned long i, j, x1;
  double start_time = clock();

  //----------------------------------------
  // Front material
  printf("\n---- Partrac v");
  printf(version);
  printf(" ----\n\n");

  // for file output
  FILE *outFile=fopen(filename, "w");

  //----------------------------------------
  // Particle positions and momenta
  static double **pos; //position
  static double **mom; //momentum
  static double **fce; //electric force at particle
  static double **bpa; //magnetic force at particle
  pos = new double *[3];
  mom = new double *[3];
  fce = new double *[3];
  bpa = new double *[3];

  for (x1 = 0; x1 < 3; ++x1) {
    pos[x1] = new double [Npar];
    mom[x1] = new double [Npar];
    fce[x1] = new double [Npar];
    bpa[x1] = new double [Npar];
  }

  // Field cubes - 3 dimensional (the '3' here refers to Ex, Ey, Ez at every coordinate space)
  double ****Efield;
  double ****Bfield;
  Efield = new double ***[3];
  Bfield = new double ***[3];
  for (x1 = 0; x1 < 3; ++x1) {
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

  // initialize particles
  initpos(pos);
  initmom(mom);

  //initialize grids
  initEfield(Efield);
  initBfield(Bfield);
  WeightF(Efield, pos, fce);
  WeightF(Bfield, pos, bpa);

  const double tmax = sqrt(2.0) * 2.0 * PI;
  double t = 0;
  const double dt = 0.1;

  // Main loop
  while (t<=tmax) {
    //Weight Electric fields to particles
    WeightF(Efield, pos, fce);
    WeightF(Bfield, pos, bpa);

    // Equations of motion
    pushmomentum(fce, mom, dt * 0.5);
    Lorentz(mom, bpa, dt);
    pushmomentum(fce, mom, dt * 0.5);
    pushxpos(mom, pos, dt);

    fprintf(outFile, "%e\t%e\t%e\t%e\n", t, mom[0][0], mom[1][0], fce[0][0]);

    t += dt;
  }

  /*// three position components
  for (i = 0; i < Nx; ++i) {
    for (j = 0; j < Ny; ++j) {
      for (k = 0; k < Nz; ++k) {
        for (x1 = 0; x1 < 3; ++x1) {
          printf("component=%i, pos=%i,%i,%i E=%f\n", (int)x1, (int)i, (int)j, (int)k, Bfield[x1][i][j][k]);
        }
      }
    }
  }

  long unsigned parnum;
  for (x1 = 0; x1 < 3; ++x1) {
    for (parnum = 0; parnum < Npar; ++parnum) {
      printf("component=%i, number=%i, pos=%f, mom=%f, fce=%f\n", (int)x1, (int)parnum, pos[x1][parnum], mom[x1][parnum], fce[x1][parnum]);
    }
  }*/

  //----------------------------------------
  // Finish
  double end_time = clock();
  printf("finished\n");
  printf("\nEND OF CALCULATION\n");
  printf("Calculation took %f\n", (end_time - start_time) / CLOCKS_PER_SEC);

  return 0;
}
