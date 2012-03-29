#ifndef METROLINEARUTIL_H
#define METROLINEARUTIL_H

#include <cstdlib>
#include <math.h>
#include "../../Utilities/src/metroUtil.h"
#include "../test/baseTests.h"
#include <curand.h>
#define THREADS_PER_BLOCK 128
#define PI 3.14159265358979323

/**
  @param idx - the index in the 1 dimensional array of energies
  @return - the id of the X atom
*/
int getXFromIndex(int idx);

/**
  @param x - the id of the x atom
  @param idx - the index in the 1 dimensional array of energies
*/
int getYFromIndex(int x, int idx);

/**
  @param x - the variable to make periodic
  @param box - the size of the period
  @return - x after being made periodic
*/
double makePeriodic(double x, const double box);

/**
  @param x - the value to continue on the other side of the box
  @param box - the length of one side of the box (cube)
*/
double wrapBox(double x, double box);

/**
    Keeps a molecule intact within the box
    @param molecule - molecule to keep in box
    @param enviro - environment structure to keep bounds
*/
void keepMoleculeInBox(Molecule *molecule, Environment *enviro);

/**
  Calculates the energy between 2 atoms
  @param atom1 - the first atom in the pair
  @param atom2 - the second atom in the pair
  @param enviro - the environmental variables
*/
double calc_lj(Atom atom1, Atom atom2, Environment enviro); 

/**
    Global function for GPU to assign random doubles to atom positions.
    @param dev_doubles - array of randomly generated doubles 0.0 to 1.0
    @param atoms - array of atoms to be positioned
    @param enviro - Environment of the simulation
*/
void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro);

/**
  Generate random positions for atoms in the box
  nVidia CURAND reference: http://developer.download.nvidia.com/compute/cuda/5_0/toolkit/docs/CURAND_Library.pdf
  @param molec - array of molecules to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatePoints(Molecule *molec, Environment *enviro);

/**
  This is a wrapper function for the calcEnergy kernel.
  @param molec - the array of molecules
  @param enviro - the environmental variables
  @return - the total energy of the system.
*/
double calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules=NULL);

/**
  Wrapper function for the calcEnergy kernel.
  @param molecules - array of molecules in the system
  @param enviro - the environment of the system.
*/
double calcEnergyWrapper(Molecule *molecules, Environment *enviro);

/**
  This calculates energy between two atoms on host for fringe
  disagreeable atoms.
  @param atom1 - first atom
  @param atom2 - second atom
  @parm enviro - the environment for the system
*/
//double calcEnergyOnHost(Atom atom1, Atom atom2, Environment *enviro);


/**
  Calculates the energy between n atoms where n is the
  the number of threads in the block. The block's sum is then stored 
  in the energySum array at its block index position.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @param *energySum - the array of block sums
*/
void calcEnergy(Atom *atoms, Environment *enviro, double *energySum);

/**
  Calculates the charge portion of the force field energy calculation between two atoms.
  @param atom1 - the first atom in the calculation
  @param atom2 - the second atom in the calculation
  @return - the charge portion of the force field.
  
  Assigned to Alex
*/
double calcCharge(Atom atom1, Atom atom2, Environment *enviro);

/**
  Returns the "fudge factor" to be used in force field calculation.
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @return - 1.0 if the atoms are in seperate molecules
            .5 if the bond traversal distance is greater or equal to 4
            0.0 if the bond traversal is less than 4

  Assigned to TBD
*/
double getFValue(Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro);

/**
  Return if the two atom ids are have a hop value >=3
  returns 1 if true and 0 if false
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param molecule - the molecule that contains atom1 and atom 2
*/
int hopGE3(int atom1, int atom2, Molecule *molecule);

/**
  Returns sqrt(d1 * d2)
  @param d1 - the first double in the calculation
  @param d2 - the second double in the calculation
  @return - sqrt(d1 * d2)
*/
double calcBlending(double d1, double d2);

/**
  Returns the molecule id from the atomid
  @param atom - the atom from which to find the molecule
  @param molecules - the list of molecules to be searched
  @param return - returns a pointer to the molecule
*/
Molecule* getMoleculeFromAtomID(Atom *a1, Molecule *molecules, Environment *enviro);

/**
  Rotates a molecule about a given atom a random amount
  @param molecule - the molecule to be rotated
  @param pivotAtom - the atom that the molecule is rotated about
  @param maxRotation - the maximum number of degrees for each axis of rotation

*/
void rotateMolecule(Molecule molecule, Atom pivotAtom, double maxRotation);

/****************************
  Begin Stubs for outputs
****************************/

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double solventAccessibleSurfaceArea();

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double soluteSolventDistributionFunction();

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double atomAtomDistributionFunction();

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double solventSolventTotalEnergy();

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double soluteSolventTotalEnergy();

#endif
