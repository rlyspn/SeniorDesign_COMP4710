#ifndef METROPARALLELUTIL_H
#define METROPARALLELUTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

struct Atom{
    double x; // x coord of the atom
    double y; // y coord of the atom
    double z; // z coord of the atom

    unsigned long id; // unique id of the atom
    string name; // name of the atom or it's symbol or some other id.

    double sigma; // sigma value for the atom for the LJ calculations
    double epsilon; // epsilon value for the atom for the LJ calculation
    double charge; // epsilon value for the atom for the LJ calculation
};

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge);
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon);
Atom createAtom(string newName, unsigned long id, double x, double y, double z);
Atom createAtom(unsigned long id, double x, double y, double z);

struct Environment{
    double x; // length of the box in the x direction
    double y; // length of the box in the y direction
    double z; // length of the box in the z direction

    double maxTranslation; // the maximum distance that an atom can move
    double temperature; // the temperature of the box

    int numOfAtoms; // the number of atoms in the environment
};

Environment createEnvironment(double x, double y, double z, double maxTrans, double temp, int numOfAtoms);


/**
  Prints the position of all the atoms in an array.
  @param atoms - array of atoms
  @param count - number of atoms
*/
void printAtoms(Atom *atoms, int count);

/**
  Writes the list of atoms to the file named filename
  @param atoms - list of atoms to be written
  @param enviro - environmental information
  @param filename - name of the file to be written
*/
void writeOutAtoms(Atom *atoms, Environment *enviro, string filename);

#endif //METROPARALLELUTIL_H
