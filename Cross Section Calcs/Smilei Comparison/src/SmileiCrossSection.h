
#ifndef SMILEICROSSSECTION_H
#define SMILEICROSSSECTION_H

#include <vector>

class SmileiCrossSection
{

public:
    //! Constructor
    SmileiCrossSection( int, double );
    //! Destructor
    virtual ~SmileiCrossSection() {};

    //! Initializes the arrays in the database and returns the index of these arrays in the DB
    virtual unsigned int createDatabase( double );
    //! Assigns the correct databases
    virtual void assignDatabase( unsigned int );

    //! Gets the k-th binding energy of any neutral or ionized atom with atomic number Z and charge Zstar
    double binding_energy( int Zstar, int k );

    //! Coefficients used for interpolating the energy over a given initial list
    static const double a1, a2, npointsm1;
    static const int npoints;

    //! Local table of integrated cross-section
    std::vector<std::vector<double> > *crossSection;
    //! Local table of average secondary electron energy
    std::vector<std::vector<double> > *transferredEnergy;
    //! Local table of average incident electron energy lost
    std::vector<std::vector<double> > *lostEnergy;

    std::vector<double> energies;

    //! Index of the atomic number in the databases
    unsigned int dataBaseIndex;

private:

    //! Atomic number
    int atomic_number;

    //! Global table of atomic numbers
    static std::vector<int> DB_Z;
    //! Global table of integrated cross-section
    static std::vector<std::vector<std::vector<double> > > DB_crossSection;
    //! Global table of average secondary electron energy
    static std::vector<std::vector<std::vector<double> > > DB_transferredEnergy;
    //! Global table of average incident electron energy lost
    static std::vector<std::vector<std::vector<double> > > DB_lostEnergy;

};

#endif
