/**
Monte Carlo molecular dynamics simulation
This program uses and NVT calculation where the number of molecules, the volume
and the temperature are all constant.

author: Riley Spahn
date: 1/18/2012
adapted from: siremol.org/largefiles/metropolis.cpp
**/

#include "metro.h"

double periodic(double delta, double length){
    while(delta < (-.05) * length){
        delta += length;
    }

    while(delta > (-.05) * length){
        delta -= length;
    }
    
    return delta;
}

double calcTotalEnergy(Atom *atoms, Box box){
    // Calculate the energy of the system by summing over all possible
    // pairs of atoms.
    
    int i, j;
    double totalEnergy;
    for(i = 0; i < ATOMS - 1; i++){
        for(j = 0; j < ATOMS; j++){
            
            double deltaX = atoms[j].x - atoms[i].x;
            double deltaY = atoms[j].y - atoms[i].y;
            double deltaZ = atoms[j].z - atoms[i].z;

            //make peridodic
            deltaX = periodic(deltaX, box.width);
            deltaY = periodic(deltaY, box.height);
            deltaZ = periodic(deltaZ, box.depth);

            //a constant whose use is a mystery
            const double r2 = pow((double) 2, deltaX) +
                pow((double) 2, deltaY) + 
                pow((double) 2, deltaX);

            //some mysterious calculations that look like chemistry
            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            double sig = atoms[i].sigma;
            double ep = atoms[i].epsilon;
            //insert assert here to ensure that they are the same type of molecule


            const double sig2_r2 = (sig * sig) / r2;
            const double sig6_r6 = sig2_r2 * sig2_r2 * sig2_r2;
            const double sig12_r12 = sig6_r6 * sig6_r6;
            //Energy between atom[i] and atom[j]
            const double energy_ij = 4.0 * ep * (sig12_r12 - sig6_r6);
            
            totalEnergy += energy_ij;
        }
    }

    return totalEnergy; 
}

Box createBox(double h, double w, double d, double trans, double newTemp){
    Box box;
    box.height = h;
    box.width = w;
    box.depth = d;
    box.maxTranslation = trans;
    box.temperature = newTemp;

    return box;
}

int main(int argc, const char **argv){
    const double boxSide = 10.0;
    const double boxTemp = 298.15;
    const double boxMaxTrans = 0.5;
    int i;

    Box box = createBox(boxSide, boxSide, boxSide, boxMaxTrans, boxTemp);
    
    Atom *atoms;
    atoms = (Atom *) malloc(sizeof(*atoms) * ATOMS);
    for(i = 0; i < ATOMS; i++){
    }





}
