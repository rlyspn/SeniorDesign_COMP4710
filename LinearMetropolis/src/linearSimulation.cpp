#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"
#include "../src/metroLinearUtil.h"
   using namespace std;

// boltzman constant
   const double kBoltz = 1.987206504191549E-003;
	
	const double maxRotation = 10.0; // degrees

    double randomFloat(const double start, const double end)
   {
      return (end-start) * (double(rand()) / RAND_MAX) + start;
   }
   
    void runLinear(Molecule *molecules, Environment *enviro, int numberOfSteps, string stateFile, string pdbFile){
      int accepted = 0; // number of accepted moves
      int rejected = 0; // number of rejected moves
      int numberOfAtoms = enviro->numOfAtoms;
      double maxTranslation = enviro->maxTranslation;
      double temperature = enviro->temperature;
      double kT = kBoltz * temperature;
      double finalEnergy = 0.0;
      
   	generatePoints(molecules, enviro);
   	cout << "points generated" << endl;
   
      Atom *atoms;
      atoms = (Atom *)malloc(sizeof(Atom) * numberOfAtoms);
    //create array of atoms from arrays in the molecules
      cout << "Allocated array" << endl;
      int atomIndex = 0;
      for(int i = 0; i < enviro->numOfMolecules; i++){
         for(int j = 0; j < molecules[i].numOfAtoms; j++){
            atoms[atomIndex] = molecules[i].atoms[j];
            atomIndex++;
            
         }
      }
   
      cout << "array assigned " << endl;
		int atomTotal = 0;
      int aIndex = 0;
      int mIndex = 0;
      
      /** Not Needed only?? 
   	//copies back atoms into molecules
      int atomTotal = 0;
      int aIndex = 0;
      int mIndex = 0;
      while(atomTotal < numberOfAtoms){
         molecules[mIndex].atoms[aIndex] = atoms[atomTotal];
         atomTotal++;
         aIndex++;
         if(aIndex == molecules[mIndex].numOfAtoms){
            aIndex = 0;
            mIndex++;
         }
      }*/
    
      printState(enviro, molecules, enviro->numOfMolecules, "initialState");
      
      for(int move = 0; move < numberOfSteps; move++){
         double oldEnergy = calcEnergyWrapper(atoms, enviro);

         //Pick a molecule to move
         int moleculeIndex = randomFloat(0, enviro->numOfMolecules);
         Molecule toMove = molecules[moleculeIndex];
         Molecule oldToMove;
         copyMolecule(&oldToMove, &toMove);

         //Pick an atom in the molecule about which to rotate
         int atomIndex = randomFloat(0, molecules[moleculeIndex].numOfAtoms);
         Atom vertex = molecules[moleculeIndex].atoms[atomIndex];
      
         const double deltaX = randomFloat(-maxTranslation, maxTranslation);
         const double deltaY = randomFloat(-maxTranslation, maxTranslation);
         const double deltaZ = randomFloat(-maxTranslation, maxTranslation);
      
         const double degreesX = randomFloat(-maxRotation, maxRotation);
         const double degreesY = randomFloat(-maxRotation, maxRotation);
         const double degreesZ = randomFloat(-maxRotation, maxRotation); 
        
         toMove = moveMolecule(toMove, vertex, deltaX, deltaY, deltaZ,
                degreesX, degreesY, degreesZ);
      
         keepMoleculeInBox(&toMove, enviro);

         molecules[moleculeIndex] = toMove;
         
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
            if (move == numberOfSteps - 1){
                finalEnergy = newEnergy;
            }
         }
         else{
            rejected++;
            //restore previous configuration
            //atoms[atomIndex] = oldAtom;
            molecules[moleculeIndex] = oldToMove;
            if (move == numberOfSteps - 1){
                finalEnergy = oldEnergy;
            }
         }
      
        
        //Print the state every 100 moves.
         if(move % 100 == 0){
            atomTotal = 0;
            aIndex = 0;
            mIndex = 0;
            while(atomTotal < numberOfAtoms){
               molecules[mIndex].atoms[aIndex] = atoms[atomTotal];
               atomTotal++;
               aIndex++;
               if(aIndex == molecules[mIndex].numOfAtoms){
                  aIndex = 0;
                  mIndex++;
               }
            }
         
            printState(enviro, molecules, enviro->numOfMolecules, stateFile);
         }
      //        cout << "Accepted: " << accepted << endl;
      }
      atomTotal = 0;
      aIndex = 0;
      mIndex = 0;
      while(atomTotal < numberOfAtoms){
         molecules[mIndex].atoms[aIndex] = atoms[atomTotal];
         atomTotal++;
         aIndex++;
         if(aIndex == molecules[mIndex].numOfAtoms){
            aIndex = 0;
            mIndex++;
         }
      }
      cout << "Final Energy = " << finalEnergy << endl;
   }

	


    int main(int argc, char ** argv){
    
    
    /***===================
      TEMPORARILY WITHOUT COMMANDLINE ARGUMENTS
    if(argc != 3){
        printf("Error.  Expected metro flag path.\n");
        exit(1);
    }
   
    // z matrix or state flage
    string flag = argv[1];
    //path to the configuration file
    string configPath = argv[2];
    ====================*/
   
      string flag = "-z";
      string configPath = "bin/demoConfiguration.txt";
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
         cout << "Configuration Path = " << configScan.getConfigPath() << endl;
        //set up Opls scan and zMatrixScan
         cout << "OPLS File Path = " << configScan.getOplsusaparPath() << endl;
         string oplsPath = configScan.getOplsusaparPath();
         Opls_Scan oplsScan (oplsPath);
         oplsScan.scanInOpls(oplsPath);
         cout << "Created oplsScan" << endl;
        
         Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan);
         if (zMatrixScan.scanInZmatrix() == -1){
            cerr << "Error, Could not open: " << configScan.getZmatrixPath() << endl;
            exit(1);
         }
         cout << "Opened zMatrix file" << endl;
        //Convert molecule vectors into an array
         molecules = (Molecule *)malloc(sizeof(Molecule) * enviro.numOfMolecules);
         int moleculeIndex = 0;
         int atomCount = 0;
         while(moleculeIndex < enviro.numOfMolecules){
            vector<Molecule> molecVec = zMatrixScan.buildMolecule(atomCount);
            //cout << "Vector size = " << molecVec.size() << endl;
            //cycle through the number of molecules from the zMatrix
            for(int j = 0; j < molecVec.size(); j++){
                //Copy data from vector to molecule
               Molecule molec1 = molecVec[j];
               molecules[moleculeIndex].id = molec1.id;
               molecules[moleculeIndex].atoms = (Atom *)malloc(sizeof(Atom) * molec1.numOfAtoms);
               molecules[moleculeIndex].bonds = (Bond *)malloc(sizeof(Bond) * molec1.numOfBonds);
               molecules[moleculeIndex].angles = (Angle *)malloc(sizeof(Angle) * molec1.numOfAngles);
               molecules[moleculeIndex].dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * molec1.numOfDihedrals);
               molecules[moleculeIndex].hops = (Hop *)malloc(sizeof(Hop) * molec1.numOfHops);
               molecules[moleculeIndex].numOfAtoms = molec1.numOfAtoms;                
               molecules[moleculeIndex].numOfBonds = molec1.numOfBonds;
               molecules[moleculeIndex].numOfAngles = molec1.numOfAngles;
               molecules[moleculeIndex].numOfDihedrals = molec1.numOfDihedrals;
               molecules[moleculeIndex].numOfHops = molec1.numOfHops;
            
                //get the atoms from the vector molecule
               for(int k = 0; k < molec1.numOfAtoms; k++){
                  molecules[moleculeIndex].atoms[k] = molec1.atoms[k];
               }
            	 //get the bonds from the vector molecule
               for(int k = 0; k < molec1.numOfBonds; k++){
                  molecules[moleculeIndex].bonds[k] = molec1.bonds[k];
               }
            	  //get the angles from the vector molecule
               for(int k = 0; k < molec1.numOfAngles; k++){
                  molecules[moleculeIndex].angles[k] = molec1.angles[k];
               }
            	  //get the dihedrals from the vector molecule
               for(int k = 0; k < molec1.numOfDihedrals; k++){
                  molecules[moleculeIndex].dihedrals[k] = molec1.dihedrals[k];
               }
            	  //get the hops from the vector molecule
               for(int k = 0; k < molec1.numOfHops; k++){
                  molecules[moleculeIndex].hops[k] = molec1.hops[k];
               }               
                
                //cout << "AtomIndex ID: " << molecules[moleculeIndex].atoms[0].id << endl;
               atomCount += molecules[moleculeIndex].numOfAtoms;
                //cout << "MolecIndex " << moleculeIndex << endl;
               
               moleculeIndex++;
            }
         }
         enviro.numOfAtoms = atomCount;
         cout << "Created Molecule Array" << endl;
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
      cout << "Beginning simulation with: " << endl;
      printf("%d atoms\n%d molecules\n%d steps\n", enviro.numOfAtoms,
            enviro.numOfMolecules, simulationSteps);
    runLinear(molecules, &enviro, simulationSteps, configScan.getStateOutputPath(),
           configScan.getPdbOutputPath());
    
   }
