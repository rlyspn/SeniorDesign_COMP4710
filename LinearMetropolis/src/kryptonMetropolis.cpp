#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "metroUtil.h"


const int atomNumber = 25;

const int numMoves = 10;

const double boxSize = 10.0;

const double maxTranslation = 0.5;

const double temperature = 298.15;



double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

void printAtoms(Atom *atoms, int count){
    for(int i = 0; i < count; i++){
        printf("%f, %f, %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

int main(int argc, char **argv){
    double sigma = 3.624;

    double epsilon = 0.317;
    
    srand(time(NULL)); 
    Atom *atoms = new Atom[atomNumber];

    //atoms = (Atom *)malloc(sizeof(*atoms) * atomNumber);

    for(int i = 0; i < atomNumber; i++){
        double x = rand(0, boxSize);
        double y = rand(0, boxSize);
        double z = rand(0, boxSize);

      //  Atom newAtom = createAtom((unsigned long) i, x, y, z, sigma, epsilon);
       // atoms[i] = newAtom;
    }

    printAtoms(atoms, atomNumber);

    int accepted = 0;
    int rejected = 0;

    int move;
    for(move = 0; move < numMoves; move++){
        const double oldEnergy = calculateEnergy(atoms, atomNumber, boxSize); 

        //save old atom
        int atomIndex = int(rand(0, atomNumber));
        Atom oldAtom = atoms[atomIndex];

        const double deltaX = rand(-maxTranslation, maxTranslation);
        const double deltaY = rand(-maxTranslation, maxTranslation);
        const double deltaZ = rand(-maxTranslation, maxTranslation);

        atoms[atomIndex] = createAtom((unsigned long) atomIndex, oldAtom.x +
                deltaX, oldAtom.y + deltaY, oldAtom.z + deltaZ, sigma, epsilon);

        double newEnergy = calculateEnergy(atoms, atomNumber, boxSize);


    }
}
