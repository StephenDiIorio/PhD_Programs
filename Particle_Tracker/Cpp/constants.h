/*
  --------------------------------------------------
  PARTRAC v1.0: 3D particle tracking code
  AGRT 2010
  Modernized by Stephen DiIorio 2019

  Global Constants
  --------------------------------------------------
*/

#include <math.h>

static const unsigned int ndims = 3;
static const unsigned int xIndx = 0;
static const unsigned int yIndx = 1;
static const unsigned int zIndx = 2;

//Constants
static const double PI = acos(-1.0);
static const double e = 1.602176e-19;
static const double me = 9.10938e-31;
static const double mu0 = 4.0*PI*1e-7;
static const double c = 299792458.;

// particle definitions

// charge to mass ratio - normalized to |e|/m_e
static const double qoverm = -1.0;

// number of particles
static const unsigned long Npar = 5000;

// Field grid size
static const unsigned long Nx = 250;
static const unsigned long Ny = 250;
static const unsigned long Nz = 3; // need to provide some buffer cells so
                                   // so weighting stays within grid

// Domain of -40e-6 to 40e-6 with 250 grid cells in Epoch
// normalized Osiris units uisng n0 = 0.17863390738e26
static const double n0 = 0.17863390738e26; // in m^-3
static const double dx[ndims] = {0.2545087717226359, 0.2545087717226359, 0.2545087717226359};
static const double idx[ndims] = {1.0 / dx[xIndx], 1.0 / dx[yIndx], 1.0 / dx[zIndx]}; //inverse of dx

//------------------------------------------------------
//Messages, input, and output
static const int Ndumps = 1000;
static const int Nmessages = 10;

static const int negTimeOffset = 40;

static char xInput[15];
static char yInput[15];
static char zInput[15];
static const int numFiles = 2;//401; //+1 to account for t=0

static const char trajFilename[9] = "traj.txt";
static const char sliceFilename[10] = "slice.txt";
static const char version[4] = "1.0";
