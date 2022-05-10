#include "MyIonizationTables.h"

// Gets the ionization energy of a given ion
double MyIonizationTables::ionization_energy( int atomic_number, int Zstar )
{
    return ionizationEnergy[( atomic_number*( atomic_number-1 ) )/2 + Zstar];
}

// Gets the azimuthal atomic number of a given ion
int MyIonizationTables::azimuthal_atomic_number( int atomic_number, int Zstar )
{
    return ( int )azimuthalQuantumNumber[( atomic_number*( atomic_number-1 ) )/2 + Zstar];
}

// Gets the k-th binding energy in any neutral or ionized atom with atomic number Z and charge Zstar
// We use the formula by Carlson et al., At. Data Nucl. Data Tables 2, 63 (1970)
double MyIonizationTables::binding_energy( int atomic_number, int Zstar, int k )
{
    // int offset = ( atomic_number*( atomic_number-1 ) )/2;
    // return ( (ionizationEnergy[offset + Zstar                ] /510998.9)
    //          - bindingEnergy   [offset + atomic_number-Zstar-1]
    //          + bindingEnergy   [offset + k                    ]
    //        );///510998.9; // converted to mc^2
    return bindingEnergy[(atomic_number * (atomic_number - 1)) / 2 + Zstar];
}

double MyIonizationTables::screening_constant( int atomic_number, int Zstar )
{
	return screeningConstants[( atomic_number*( atomic_number-1 ) )/2 + Zstar];
}
