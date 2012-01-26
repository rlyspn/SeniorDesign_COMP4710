#include <cuda.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "metroCudaUtil.cuh"
#include "metroParallelUtil.h"

/**
Monte Carlo molecular dynamic simulation of krypton based on code from
http://siremol.org/largefiles/metropolis.cpp provided by Dr. Orlando Acevedo

Date: 1/26/2012
Author(s): Riley Spahn, Seth Wooten, Alexander Luchs
*/

// number of atoms in the simulation
const int numberOfAtoms = 25;
// number of moves to run the simulation
const int numberOfMoves = 10000;
// maximum translation of an atom per move in any direction
const double maxTranslation = 0.5;
// the tempeature of the simulation (Kelvin)
const double temperature = 298.15;
// the sigma value of krypton used in the LJ simulation
const double kryptonSigma = 3.624;
// the epsilon value of krypton used in the LJ simulation
const double kryptonEpsilon = 0.317;
// boltzman constant
const double kBoltz = 1.987206504191549E-003;
const double kT = kBoltz * temperature;
// lengths of the box in each dimension
const double xLength = 10.0;
const double yLength = xLength;
const double zLength = xLength;


double randomFloat(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}


int main(int argc, char ** argv){
    //initialize random number generator;
    srand(time(NULL));

    //Generate Environment
    Environment enviro = createEnvironment(xLength, yLength, zLength,
            maxTranslation, temperature, numberOfAtoms);
    //Generate atoms with random positions
    Atom atoms[numberOfAtoms];
    generatePoints(atoms, &enviro);

    printf("Generated Atoms: \n");
    printAtoms(atoms, numberOfAtoms);

    int accepted; // number of accepted moves
    int rejected; // number of rejected moves
   
    atoms[0].sigma = kryptonSigma;
    atoms[0].epsilon = kryptonEpsilon;

    for(int move = 0; move < numberOfMoves; move++){
        double oldEnergy = calcEnergyWrapper(atoms, enviro);

        int atomIndex = randomFloat(0, numberOfAtoms);
        Atom oldAtom = atoms[atomIndex];
       
        //From here ========== to 
        const double deltaX = randomFloat(-maxTranslation, maxTranslation);
        const double deltaY = randomFloat(-maxTranslation, maxTranslation);
        const double deltaZ = randomFloat(-maxTranslation, maxTranslation);

        atoms[atomIndex] = createAtom((unsigned long) atomIndex, oldAtom.x +
        deltaX, oldAtom.y + deltaY, oldAtom.z + deltaZ, kryptonSigma, kryptonEpsilon);
        //here ===== could be its own function

        double newEnergy = calcEnergyWrapper(atoms, enviro);

        bool accept = false;

        if(newEnergy < oldEnergy){
            accept = true;
        }
        else{
            double x = exp(-(newEnergy - oldEnergy) / kT);

            if(x >= randomFloat(0.0, 1.0)){
                accept = true;
            }
            else{
                accept = false;
            }
        }

        if(accept){
            accepted++;
        }
        else{
            rejected++;
            //restore previous configuration
            atoms[atomIndex] = oldAtom;
        }

        printf("accepted: %d\nrejected: %d\n\n", accepted, rejected);

    }
    
}
