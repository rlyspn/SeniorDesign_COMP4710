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

Atom createAtom(double newSig, double newEp, double newX, double newY, double newZ){
    Atom atom;

    atom.sigma = newSig;
    atom.epsilon = newEp;
    atom.x = newX;
    atom.y = newY;
    atom.z = newZ;

    return atom;
}

double randomPoint(double sideLength){
    return sideLength * (double(rand()) / RAND_MAX);
}

double wrap(double x, double max){
    if(x > max){
        x -= max;
    }
    else if(x < 0){
        x += max;
    }
    return x;
}

Atom translateAtom(Atom atom, Box box){
    double deltaX = randomPoint(box.maxTranslation * 2.0) - box.maxTranslation;
    double deltaY = randomPoint(box.maxTranslation * 2.0) - box.maxTranslation;
    double deltaZ = randomPoint(box.maxTranslation * 2.0) - box.maxTranslation;

    double newX = wrap(deltaX + atom.x, box.width);
    double newY = wrap(deltaY + atom.y, box.height);
    double newZ = wrap(deltaZ + atom.z, box.depth);
    
    printf("Moved %f, %f, %f\n", deltaX, deltaY, deltaZ);

    atom.x = newX;
    atom.y = newY;
    atom.z = newZ;

    return atom;

}

void printAtom(Atom atom, int i){
    printf("%d: (%f, %f, %f)\n", i, atom.x, atom.y, atom.z);
}

void printAtom(Atom atom){
    printf("(%f, %f, %f)\n", atom.x, atom.y, atom.z);
}


void printAtoms(Atom *atoms){
    for(int i = 0; i < ATOMS; i++){
       printAtom(atoms[i], i); 
    }
}


int main(int argc, const char **argv){
    srand(time(NULL));

    const double boxSide = 10.0;
    const double boxTemp = 298.15;
    const double boxMaxTrans = 0.5;

    const double kryptonSigma = 3.624;
    const double kryptonEpsilon = 0.317;

    const double kBoltz = 1.987206504191549E-003;
    const double kT = kBoltz * boxTemp;
    int i;

    Box box = createBox(boxSide, boxSide, boxSide, boxMaxTrans, boxTemp);
    
    Atom *atoms;
    atoms = (Atom *) malloc(sizeof(*atoms) * ATOMS);
   
    //this is potentially parallelizable
    //Create the atoms that are going to be used for this simulation.
    for(i = 0; i < ATOMS; i++){
        double x = randomPoint(box.width);
        double y = randomPoint(box.height);
        double z = randomPoint(box.depth);

        atoms[i] = createAtom(kryptonSigma, kryptonEpsilon, x, y, z);
    }

    //printAtoms(atoms);

    int acceptedMoves = 0;
    int rejectedMoves = 0;

    for(i = 0; i < MOVES; i++){
        const double oldEnergy = calcTotalEnergy(atoms, box);
        
        //index of the atom that is randomly chosen
        int atomIndex = randomPoint(ATOMS);
        Atom oldAtom = atoms[atomIndex];
        
        atoms[atomIndex] = translateAtom(atoms[atomIndex], box);
       // printAtom(oldAtom);
        //printAtom(atoms[atomIndex]);

        const double newEnergy = calcTotalEnergy(atoms, box);
        

        bool accept = false;

        if(newEnergy <= oldEnergy){
            accept = true;
        }
        else{
            //Monte Carlo it
            double mcValue = exp( -(newEnergy - oldEnergy) / kT);
            double randPoint = randomPoint(1.0);
            if(mcValue >= randPoint){
                accept = true;
            }
            else
                accept = false;

        }

        //restore old move if not accepted
        if(accept){
            acceptedMoves++;
        }
        else{
            rejectedMoves++;
            atoms[atomIndex] = oldAtom;
        }
        
    }

   // printAtoms(atoms);
    printf("accepted: %d\nrejected: %d\n", acceptedMoves, rejectedMoves);

}
