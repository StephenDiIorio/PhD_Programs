// GLOBAL CONSTANTS
//*****************************************************************
#include <math.h>

//Constants
const double PI = acos(-1.0);
const double e = 1.602176e-19;
const double me = 9.10938e-31;
const double mu0 = 4.0*PI*1e-7;
const double c = 299792458;

// particle definitions

// charge to mass ratio - normalized to |e|/m_e
const double qoverm = -1.0;

// number of particles
const unsigned long Npar = 1;

// Field grid size
const unsigned long Nx = 100;
const unsigned long Ny = 100;
const unsigned long Nz = 100;

const double dx[3] = {1.0, 1.0, 1.0}; //the step size
const double idx[3] = {1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2]}; //inverse of dx

//------------------------------------------------------
//Messages and output
const int Ndumps = 1000;
const int Nmessages = 10;
const char filename[] = "data.txt";
const char version[5] = "1.0";
