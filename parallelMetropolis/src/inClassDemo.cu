#include <cuda.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "metroCudaUtil.cuh"
#include "../../Utilities/src/metroUtil.h"

/**
Monte Carlo molecular dynamic simulation of krypton based on code from
http://siremol.org/largefiles/metropolis.cpp provided by Dr. Orlando Acevedo

Date: 1/26/2012
Author(s): Riley Spahn, Seth Wooten, Alexander Luchs
*/

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

void runParallel(Atom *atoms, Environment *enviro, int numberOfSteps){
    int accepted = 0; // number of accepted moves
    int rejected = 0; // number of rejected moves
    int numberOfAtoms = enviro->numOfAtoms;

    for(int move = 0; move < numberOfSteps; move++){
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


    }

}

int main(int argc, char ** argv){
    int atomSizes[] = {500};
    int numberOfMoves = 1000;

    //initialize random number generator;
    srand(time(NULL));

    for(int i = 0; i < 1; i++){
       int numberOfAtoms = atomSizes[i];

        //Generate Environment
        Environment enviro = createEnvironment(xLength, yLength, zLength,
                maxTranslation, temperature, numberOfAtoms);
        //Generate atoms with random positions
        Atom atoms[numberOfAtoms];
        printf("Generating %d atoms.\n", numberOfAtoms);
        generatePoints(atoms, &enviro);
        atoms[0].sigma = kryptonSigma;
        atoms[0].epsilon = kryptonEpsilon;

        //Run the simulation parallel

        clock_t startTime, endTime;
        printf("Simulating %d atoms for %d steps\n", numberOfAtoms, numberOfMoves);
        
        startTime = clock();
        runParallel(atoms, &enviro, numberOfMoves);
        endTime = clock();

        double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;
        printf("Time = %f seconds.\n\n", diffTime);


    
    }



   


    
}
