#ifndef BASETEST_H
#define BASETEST_H

#include <math.h>
#include <sys/time.h>
#include "../../Utilities/src/metroUtil.h"

double make_periodic(double x, double box);

double wrap_into_box(double x, double box);

double calculate_energy(double **coords, int n_atoms, double *box_size, double sigma, double epsilon);

/**
  Calculates the energy assuming that all atoms are the same element.
  @param atoms - array of atoms that will be used to calculate energy
  @param atomNumber - the number of atoms in the system.
  @param boxSize- the length of one side of the box
*/

double calculate_energy(Atom *atoms, Environment *enviro, Molecule *molecules=NULL);

//retrieves the difference in time
long timevaldiff(struct timeval *starttime, struct timeval *finishtime);

double calc_r_value(Atom a1, Atom a2, Environment enviro);

double calc_charge(Atom a1, Atom a2, Environment enviro);

/**
  Returns the molecule id from the atomid (on host)
  @param atom - the atom from which to find the molecule
  @param molecules - the list of molecules to be searched
  @param return - returns the id of the molecule
*/
int getMoleculeFromIDLinear(Atom a1, Molecule *molecules, Environment enviro);

/**
  Returns the "fudge factor" to be used in force field calculation. (on host)
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @return - 1.0 if the atoms are in seperate molecules
            .5 if the bond traversal distance is greater or equal to 4
            0.0 if the bond traversal is less than 4
*/
double getFValueLinear(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);

/**
  Return if the two atom ids are have a hop value >=3
  returns 1 if true and 0 if false (on host)
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param molecule - the molecule that contains atom1 and atom 2
*/
int hopGE3Linear(int atom1, int atom2, Molecule molecule);

#endif //BASETEST_H
