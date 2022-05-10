#include <stdlib.h>
#include <iostream>

#include "SmileiCrossSection.h"

#include "SmileiIonizationTables.h"
#include "MyIonizationTables.h"

#include <cmath>


using namespace std;

// Coefficients used for energy interpolation
// The list of energies has to be in logarithmic scale,
//  with Emin=1eV, Emax=10MeV and npoints=100.
const int    SmileiCrossSection::npoints = 100;
const double SmileiCrossSection::npointsm1 = ( double )( npoints-1 );
const double SmileiCrossSection::a1 = 510998.9 ; // = me*c^2/Emin
const double SmileiCrossSection::a2 = 6.142165 ; // = (npoints-1) / ln( Emax/Emin )

// Constructor
SmileiCrossSection::SmileiCrossSection(int Z, double reference_angular_frequency_SI )
{
    atomic_number = Z;
    if( Z>0 ) {
        dataBaseIndex = createDatabase( reference_angular_frequency_SI );
        assignDatabase( dataBaseIndex );
    }
}

// Static members
vector<int> SmileiCrossSection::DB_Z;
vector<vector<vector<double> > > SmileiCrossSection::DB_crossSection;
vector<vector<vector<double> > > SmileiCrossSection::DB_transferredEnergy;
vector<vector<vector<double> > > SmileiCrossSection::DB_lostEnergy;

// Initializes the databases (by patch master only)
unsigned int SmileiCrossSection::createDatabase( double reference_angular_frequency_SI )
{
    // Leave if the database already exists with same atomic number
    for( unsigned int i=0; i<DB_Z.size(); i++ ) {
        if( atomic_number == DB_Z[i] ) {
            return i;
        }
    }

    // Otherwise, create the arrays:
    // For each ionization state, calculate the tables of integrated cross-sections
    // Pérez et al., Phys. Plasmas 19, 083104 (2012)
    energies = vector<double>(npoints, 0);
    vector<vector<double> > cs; // cross section
    vector<vector<double> > te; // transferred energy
    vector<vector<double> > le; // lost energy
    cs.resize( atomic_number );
    te.resize( atomic_number );
    le.resize( atomic_number );
    double e, ep, bp, up, ep2, betae2, betab2, betau2, s0, A1, A2, A3, sk, wk, ek;
    int N; // occupation number
    // double normalization = 2.81794e-15 * reference_angular_frequency_SI / ( 2.*299792458. ); // r_e omega / 2c
    double normalization = (2. * 3.1415) * (2.81794e-15 * 2.81794e-15) * (reference_angular_frequency_SI / 299792458.) * (reference_angular_frequency_SI / 299792458.);
    for( int Zstar=0; Zstar<atomic_number; Zstar++ ) { // For each ionization state
        cs[Zstar].resize( npoints, 0. );
        te[Zstar].resize( npoints, 0. );
        le[Zstar].resize( npoints, 0. );
        for( int i=0; i<npoints; i++ ) { // For each incident electron energy
            ep = exp( double( i )/a2 ) / a1; // = incident electron energy
            energies[i] = ep;
            N = 1;
            for( int k=0; k<atomic_number-Zstar; k++ ) { // For each orbital
                bp = SmileiIonizationTables::binding_energy( atomic_number, Zstar, k );
                // If next orbital is on same level, then continue directly to next
                std::cout << "binding energy: " << bp << std::endl;
                if( k<atomic_number-Zstar-1 ) {
                    if( bp == SmileiIonizationTables::binding_energy( atomic_number, Zstar, k+1 ) ) {
                        N++;
                        continue;
                    }
                }
                // If electron energy below the ionization energy, then skip to next level
                e = ep/bp;
                if( e>1. ) {
                    up = bp; // we assume up=bp because we don't have exact tables
                    betae2 = 1. - 1./( ( 1.+ep )*( 1.+ep ) );
                    betab2 = 1. - 1./( ( 1.+bp )*( 1.+bp ) );
                    betau2 = 1. - 1./( ( 1.+up )*( 1.+up ) );
                    s0 = normalization * N /( bp * ( betae2 + betab2 + betau2 ) );
                    ep2 = 1./( 1.+ep*0.5 );
                    ep2 *= ep2;
                    A1 = ( 1.+2.*ep )/( 1.+e )*ep2;
                    A2 = ( e-1. )*bp*bp*0.5*ep2;
                    A3 = log( betae2/( 1.-betae2 ) ) - betae2 - log( 2.*bp );
                    sk = s0*( 0.5*A3*( 1.-1./( e*e ) ) + 1. - 1./e + A2 - A1*log( e ) );
                    wk = s0 * ( 0.5*A3*( e-1. )*( e-1. )/e/( e+1. )  + 2.*log( 0.5*( e+1. ) ) - log( e )
                                + 0.25*A2*( e-1. ) - A1*( e*log( e )-( e+1. )*log( 0.5*( e+1. ) ) ) );
                    ek = wk + sk;
                    // Sum these data to the total ones
                    cs[Zstar][i] += sk;
                    // if (Zstar == 0)
                    // {
                    //     cs[Zstar][i] = 1.62957969e-6 * 1e4;
                    // }
                    // else if (Zstar == 1)
                    // {
                    //     cs[Zstar][i] = 1.41176233e-6 * 1e4;
                    // }
                    te[Zstar][i] += wk * bp;
                    le[Zstar][i] += ek * bp;
                }
                // Reset occupation number for next level
                N = 1;
            }
            // std::cout << Zstar << std::endl;
            std::cout << cs[Zstar][i] << std::endl;

            // The transferred and lost energies are averages over the orbitals
            if( cs[Zstar][i]>0. ) {
                te[Zstar][i] /= cs[Zstar][i];
                le[Zstar][i] /= cs[Zstar][i];
            }
        }
    }

    // Add the new arrays to the static database
    DB_Z                .push_back( atomic_number );
    DB_crossSection     .push_back( cs );
    DB_transferredEnergy.push_back( te );
    DB_lostEnergy       .push_back( le );

    return DB_Z.size()-1;
}


// Assign the correct databases
void SmileiCrossSection::assignDatabase( unsigned int index )
{

    crossSection      = &( DB_crossSection     [index] );
    transferredEnergy = &( DB_transferredEnergy[index] );
    lostEnergy        = &( DB_lostEnergy       [index] );

}
