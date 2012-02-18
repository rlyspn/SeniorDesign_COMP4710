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


// Structure to represent 
struct Molecule{
    int id; // the name of the molecule
    
    Atom *atoms; // array of atoms in the molecule
    Bond *bonds; // array of bonds of the atoms in the molecule.
    Angle *angles; // angles in the molecule between atoms
    Dihedral *dihedrals; // array of dihedrals in the molecule

    int numOfAtoms; // the number of atoms in the molecule
    int numOfBonds; // the number of bonds in the molecule
    int numOfAngles; // the number of angles in the molecule
    int numOfDihedrals; // the number of dihedrals in the atom
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
    @param atomCount - the number of atoms in the molecule
*/
Molecule createMolecule(int id, 
                        Atom *atoms,
                        int atomCount);

/**
  Prints the position of all the atoms in an array.
  @param atoms - array of atoms
  @param count - number of atoms
*/
void printAtoms(Atom *atoms, int count);

/**
  @param enviro - the environment state
  @param molecules - array of molecules to be printed out
  @param numOfMolecules - the number of molecules to be written out
  @param fileName - the name of the file to be written
*/
void printState(Environment *enviro, Molecule *molecules, int numOfMolecules, string filename);

/**
  @param filename - the name of the state file
  @return - the environment recorded in the state file.
*/
Environment readInEnvironment(string filename);

/**
  @pararm filename - the name of the state file
  @return - an array of molecules
*/
vector<Molecule> readInMolecules(string filename);

/**
  input line is of the format:
  "id x y z sigma epsilon"
  @param line - the line from the state file that contains atom information.
  @return - atom containing information from the line
*/
Atom getAtomFromLine(string line);

/**
  input line is of the format:
  "x y z numOfAtoms"
  @param line - the line to be parsed
  @return - Environment read from the line
*/
Environment getEnvironmentFromLine(string line);

/**
  expected input line:
  "atom1 atom2 distance [0|1]"
  0 represents not variable
  1 represents a variable bond
  @param line - the line to be read.
  @return - bond representing the information read from the line.
*/
Bond getBondFromLine(string line);

/**
  expected input line:
  "atom1 atom2 value [0|1]"
  0 represents not variable
  1 represents variable angle

  @param line - line containing information about the angle
  @return - returns a bond 
*/
Angle getAngleFromLine(string line);

/**
  expected input line:
  "atom1 atom2 value [0|1]"
  0 represents not variable
  1 represents variable dihedral 

  @param line - line containing information about the dihedral 
  @return - returns the dihedral represented on the line
*/
Dihedral getDihedralFromLine(string line);

#endif //METROUTIL_H
