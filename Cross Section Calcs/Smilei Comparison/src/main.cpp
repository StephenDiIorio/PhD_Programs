#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "SmileiCrossSection.h"

int main()
{
    int z = 28;
    double wp = 2.62191378283848e16;

    // FILE *outFile = fopen("cross.txt", "w");

    std::ofstream outFile;

    outFile.open("cross.txt");

    SmileiCrossSection cs(z, wp);

    for (std::size_t i = 0; i < cs.energies.size(); ++i)
    {
        outFile << cs.energies[i] << '\t';
    }
    outFile << '\n';

    for (std::size_t i = 0; i < cs.crossSection->size(); ++i)
    {
        for (std::size_t j = 0; j < cs.crossSection->at(i).size(); ++j)
        {
            outFile << cs.crossSection->at(i).at(j) << '\t';
        }
        outFile << '\n';
    }

    outFile.close();

    return 0;
}
