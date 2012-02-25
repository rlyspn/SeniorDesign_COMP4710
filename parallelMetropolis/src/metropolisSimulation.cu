#include <cuda.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Zmatrix_Scan.h"

#include "metroCudaUtil.cuh"
/**
Will run simulations of any single atom system.  Can run from either z matrix
file or from a state file.

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

void runParallel(Molecule *molecules, Environment *enviro, int numberOfSteps){
    int accepted = 0; // number of accepted moves
    int rejected = 0; // number of rejected moves
    int numberOfAtoms = enviro->numOfAtoms;

    Atom *atoms = (Atom *)malloc(sizeof(Atom) * numberOfAtoms);
    //create array of atoms from arrays in the molecules
    int atomIndex = 0;
    for(int i = 0; i < enviro->numOfMolecules; i++){
        for(int j = 0; j < molecules[i].numOfAtoms; j++){
            atoms[atomIndex] = molecules[i].atoms[j];
            atomIndex++;
        }
    }

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

        double newEnergy = calcEnergyWrapper(molecules, enviro);

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

/**
  ./metroSim flag path/to/config/file
    flags:
        -z run from z matrix file spcecified in the configuration file
        -s run from state input file specified in the configuration file
*/
int main(int argc, char ** argv){
    if(argc != 3){
        printf("Error.  Expected metro flag path.\n");
        exit(1);
    }
    
    // z matrix or state flage
    string flag = argv[1];
    //path to the configuration file
    string configPath = argv[2];
    
    //Configuration file scanner
    Config_Scan configScan(configPath);
    configScan.readInConfig();

    //Environment for the simulation
    Environment enviro;
    long simulationSteps = configScan.getSteps();
    Molecule *molecules;

    //Simulation will run based on the zMatrix and configuration Files
    if(flag.compare("-z") == 0){
        printf("Running simulation based on zMatrixFile\n");
        //get environment from the config file
        enviro = configScan.getEnviro();
        //set up Opls scan and zMatrixScan
        Opls_Scan oplsScan (configScan.getOplsusaparPath());
        Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan);
        zMatrixScan.scanInZmatrix();
        //Convert molecule vectors into an array
        molecules = (Molecule *)malloc(sizeof(Molecule) * enviro.numOfMolecules);
        int moleculeIndex = 0;
        while(moleculeIndex < enviro.numOfMolecules){
            vector<Molecule> molecVec = zMatrixScan.buildMolecule(moleculeIndex);
            for(int j = 0; j < molecVec.size(); j++){
                molecules[moleculeIndex] = molecVec[j];
                moleculeIndex++;
            }
        }
    }       
    //Simulation will run based on the state file
    else if(flag.compare("-s") == 0){
        printf("Running simulation based on state file\n");
        //get path for the state file
        string statePath = configScan.getStatePath();
        //get environment from the state file
        enviro = readInEnvironment(statePath);
        //get vector of molecules from the state file
        vector<Molecule> molecVec = readInMolecules(statePath);
        enviro.numOfMolecules = molecVec.size();
        //convert vector of molecules to array
        molecules = (Molecule *)malloc(sizeof(Molecule) * molecVec.size());
        for(int i = 0; i < molecVec.size(); i++){
            molecules[i] = molecVec[i];
        }
        

    }
    else{
        printf("Error, Unknown flag.\n");
        exit(1);
    }
    runParallel(molecules, &enviro, simulationSteps); 
}