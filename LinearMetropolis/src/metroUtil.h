#ifndef METROUTIL_H
#define METROUTIL_H
#include <stdio.h>
#include<math.h>

/**
  Structure repesenting an atom.
*/

struct Atom{
    double x; // x coord of the atom
    double y; // y coord of the atom
    double z; // z coord of the atom

    unsigned long id; // unique id of the atom

    double sigma; // sigma value for the atom for the LJ calculations
    double epsilon; // epsilon value for the atom for the LJ calculation
};

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon);
Atom createAtom(unsigned long id, double x, double y, double z);

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
  Calculates the energy assuming that all atoms are the same element.
  @param atoms - array of atoms that will be used to calculate energy
  @param atomNumber - the number of atoms in the system.
  @param boxSize- the length of one side of the box
*/
double calculateEnergy(Atom *atoms, int atomNumber, double boxSize);


#endif //METROUTIL_H
