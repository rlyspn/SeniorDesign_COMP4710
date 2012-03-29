#ifndef METROUTIL_H
#define METROUTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

struct Atom{
    double x; // x coord of the atom
    double y; // y coord of the atom
    double z; // z coord of the atom

    unsigned long id; // unique id of the atom

    double sigma; // sigma value for the atom for the LJ calculations
    double epsilon; // epsilon value for the atom for the LJ calculation
    double charge; // epsilon value for the atom for the LJ calculation
};

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge);
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon);
Atom createAtom(unsigned long id, double x, double y, double z);

struct Environment{
    double x; // length of the box in the x direction
    double y; // length of the box in the y direction
    double z; // length of the box in the z direction

    double maxTranslation; // the maximum distance that an atom can move
    double temperature; // the temperature of the box

    int numOfAtoms; // the number of atoms in the environment
    int numOfMolecules; // the number of molecues in the environment
};

Environment createEnvironment(double x, double y, double z, double maxTrans, double temp, int numOfAtoms);
/**
  Writes the list of atoms to the file named filename
  @param atoms - list of atoms to be written
  @param enviro - environmental information
  @param filename - name of the file to be written
  @pararm accepts - the number of moves accepted by the system.
  @param rejects - the number of moves rejected by the system.
*/
void writeOutAtoms(Atom *atoms, Environment *enviro, string filename, int accepts, int rejects, double totalEnergy);

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
    double value; // the distance between the atoms
    bool variable; // if the distance between atoms is variable
};

/**
  @param atom1 - the first atom in the dihedral
  @param atom2 - the second atom in the dihedral
  @param value - the distance between the atoms
  @param variable - if the dihedral is variable
*/
Dihedral createDihedral(int atom1, int atom2, double value, bool variable);

/**
  Atom pairs and their node distance(hops) away from each other
  used in the fudge factor, for total energy calculations
  @param atom1 - the starting atom
  @param atom2 - the ending atom
  @param hop  - the number of nodes between start and finish
*/
struct Hop{
    int atom1;
	 int atom2;
	 int hop;
};

Hop createHop(int atom1, int atom2, int hops);

// Structure to represent 
struct Molecule{
    int id; // the name of the molecule
    
    Atom *atoms; // array of atoms in the molecule
    Bond *bonds; // array of bonds of the atoms in the molecule.
    Angle *angles; // angles in the molecule between atoms
    Dihedral *dihedrals; // array of dihedrals in the molecule
	 Hop *hops; //array containing a list of atoms that are less than 4 nodes away

    int numOfAtoms; // the number of atoms in the molecule
    int numOfBonds; // the number of bonds in the molecule
    int numOfAngles; // the number of angles in the molecule
    int numOfDihedrals; // the number of dihedrals in the atom
	 int numOfHops; // the number of Hops or pairs of atoms that are less than 4 nodes away
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
Molecule createMolecule(int id,
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals,
                        int atomCount, int angleCount, int bondCount, int dihedralCount );

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
Molecule createMolecule(int id,
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals, Hop *hops,
                        int atomCount, int angleCount, int bondCount, int dihedralCount, int hopsCount );

								
/**
    @param id - the integer id of the molecule
    @param name - the name of the molecule
    @param atoms - an array of the atoms in the molecule
    @param atomCount - the number of atoms in the molecule
*/
Molecule createMolecule(int id, 
                        Atom *atoms,
                        int atomCount);

/**
    Copies by value molec2 into molec1
    @param molec1 - the destination molecule
    @param molec2 - the source molecule
 */
void copyMolecule(Molecule *molec1, Molecule *molec2);

/**
  Prints the position of all the atoms in an array.
  @param atoms - array of atoms
  @param count - number of atoms
*/
void printAtoms(Atom *atoms, int count);

/**
  @param atoms - array of atoms to be written to the file
  @param enviro - structure holding environmental data.
  @param filename - the name of the file that is to be written
*/
void writePDB(Atom *atoms, Environment enviro, string filename);

#define DEFAULT 0
#define START 1
#define OPLS 2
#define Z_MATRIX 3
/**
  Logs output to the OutputLog file
  @param text - the text to be written to the output log
*/

void writeToLog(string text, int stamp=0 );

void printMolecule(Molecule *molec);

#endif //METROUTIL_H
