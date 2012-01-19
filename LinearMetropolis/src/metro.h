//Header file for the metropolis simulation
#ifndef METRO_H
#define METRO_H

#define ATOMS 25 // Number of atoms in the box
#define MOVES 10000000 // Number of moves to run
#define MAX_MOVE (double) 0.5 // Angstroms

#include <stdio.h>
#include <math.h>
#include <assert.h>


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


int main(int argc, const char **argv);
#endif //METRO_H
