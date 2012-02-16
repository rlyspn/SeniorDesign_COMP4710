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
Atom createAtom(string newName, unsigned long id, double x, double y, double z, double sigma, double epsilon);
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

// Structure to respresent bonds between atoms in a 
struct Bond{
    int atom1;
    int atom2;

    double distance;
    bool variable;
};

/**
  @param atom1 - the id of the first atom in the bond
  @param atom2 - the id of the second atom in the bond
  @param distance - the distance between the the two atoms
  @param variable - boolean if distance is variable
*/
Bond createBond(int atom1, int atom2, double distance, bool variable);

struct Angle{
    int atom1; // the first atom in the angle
    int atom2; // the second atom in the angle
    double value; // the angle between the atoms
    bool variable; // if the angle is variable
};

/**
  @param atom1 - the first atom in the angle
  @param atom2 - the second atom in the angle
  @param value - the value of the angle in degrees
  @param variable - if the angle between the atoms can change
*/
Angle createAngle(int atom1, int atom2, double value, bool variable);

struct Dihedral{
    int atom1; // the first atom in the dihedral
    int atom2; // the second atom in the dihedral
    double distance; // the distance between the atoms
    bool variable; // if the distance between atoms is variable
};

/**
  @param atom1 - the first atom in the dihedral
  @param atom2 - the second atom in the dihedral
  @param distance - the distance between the atoms
  @param variable - if the dihedral is variable
*/
Dihedral createDihedral(int atom1, int atom2, double distance, bool variable);


// Structure to represent 
struct Molecule{
    int id; // the name of the molecule
    string name;
    
    Atom *atoms; // array of atoms in the molecule
    Bond *bonds; // array of bonds of the atoms in the molecule.
    Angle *angles; // angles in the molecule between atoms
    Dihedral *dihedrals; // array of dihedrals in the molecule

    int atomCount; // the number of atoms in the molecule
    int bondCount; // the number of bonds in the molecule
    int angleCount; // the number of angles in the molecule
    int dihedralCount; // the number of dihedrals in the atom
};

/**
    @param id - the integer id of the molecule
    @param name - the name of the molecule
    @param atoms - an array of the atoms in the molecule
    @param bonds - an array of the bonds in the atom
    @pararm dihedrals - array of dihedrals in the atom
    @param atomCount - the number of atoms in the molecule
    @param bondCount - the number of bonds in the atom
    @param dihedralCount - the number of dihedrals in the molecule
*/
Molecule createMolecule(int id, string name,
                        Atom *atoms, Angle *angles, Dihedral *dihedrals,
                        int atomCount, int angleCount, int dihedralCount );
/**
    @param id - the integer id of the molecule
    @param name - the name of the molecule
    @param atoms - an array of the atoms in the molecule
    @param atomCount - the number of atoms in the molecule
*/
Molecule createMolecule(int id, string name,
                        Atom *atoms,
                        int atomCount);

#endif //METROPARALLELUTIL_H
