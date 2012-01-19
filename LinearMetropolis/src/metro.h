//Header file for the metropolis simulation
#ifndef METRO_H
#define METRO_H

#define ATOMS 25 // Number of atoms in the box
#define MOVES 10 // Number of moves to run
#define MAX_MOVE (double) 0.5 // Angstroms

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>


//Struct to represent a specific atom in a Simulation
struct Atom{
    //sigma and epsilon are specific leonard jones parameters and will vary
    //by element.
    double sigma;
    double epsilon;

    //Coordinates in a 3d plane of the specific atom.
    double x;
    double y;
    double z;

};

//Struct to represent the "box of atoms" in a monte carlo simulation.
struct Box{
   
    //The dimensions of the box.
    double height; // y direction
    double width; // x direction
    double depth; // z direction

    //The maximum translation of any atom in the box
    double maxTranslation;

    //The temperature of the box (kelvin)
    double temperature;
    
};

/*
@param newSig - the sigma parameter for this element
@param newEp - the epsilon parameter for this element
@param newX - the x coordinate for this element in the box
@param newY - the y coordinate for this element in the box
@param newX - the z coordinate for this element in the box
@return - returns an instance of an Atom Struct
*/
Atom createAtom(double newSig, double newEp, double newX, double newY, double newZ);

/*
* @param h - the height of the box
* @param w - the width of the box
* @param d - the height of the box
* @return - a Box struct
*/
Box createBox(double h, double w, double d, double trans, double newTemp);

/*
* @param atoms - a pointer to an array of atoms used to calculate the energy
* @param box - the box in which the atoms reside
* @return - total energy of the atoms in the box
*/
double calcTotalEnergy(Atom *atoms, Box box);

/*
* @param delta - the distance the molecule moved in a direction
* @param length - the length of the box in that direction
*/
double periodic(double delta, double length);

/**
* Generates a random point [0, sideLength]
* @param sideLength - upper bound of the range
* @return - a random number [0, sideLength]
*/
double randomPoint(double sideLength);

/**
* @param atom - the atom to be translated
* @param box - the box where the atom resides
* @return - the atom translated.
*/
Atom translateAtome(Atom atom, Box box);

/**
* The box is treated as a torus
* @param atom - the atom to be wrapped around the torus
* @param box - the box that will that contains the atom
* @return - the atom translated by a random amount
*/
Atom wrap(Atom atom, Box box);

void printAtom(Atom atom, int i);
void printAtom(Atom atom);

void printAtoms(Atom *atoms);


int main(int argc, const char **argv);
#endif //METRO_H
