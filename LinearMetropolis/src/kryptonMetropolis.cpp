#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "metroUtil.h"


const int atomNumber = 3;

const int numMoves = 1000;

const double boxSize = 10.0;

const double maxTranslation = 0.5;

const double temperature = 298.15;



double randomFloat(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

int main(int argc, char **argv){
    double sigma = 3.624;

    double epsilon = 0.317;

    double kBoltz = 1.987206504191549E-003;

    double kT = kBoltz * temperature;

    int atomSizes[] = {10, 50, 100, 500, 1000, 5000, 10000};
    
    clock_t startTime, endTime;
    for(int index = 0; index < 7; index++){
        int atomNumber = atomSizes[index];
        printf("Simulating %d atoms for %d steps\n", atomNumber, numMoves);
        startTime = clock();

        srand(time(NULL)); 
        Atom *atoms = new Atom[atomNumber];

        //atoms = (Atom *)malloc(sizeof(*atoms) * atomNumber);

        for(int i = 0; i < atomNumber; i++){
            double x = randomFloat(0, boxSize);
            double y = randomFloat(0, boxSize);
            double z = randomFloat(0, boxSize);

            Atom newAtom = createAtom((unsigned long) i, x, y, z, sigma, epsilon);
            atoms[i] = newAtom;
        }


        int accepted = 0;
        int rejected = 0;

        int move;
        for(move = 0; move < numMoves; move++){
            //printAtoms(atoms, atomNumber);
            const double oldEnergy = calculateEnergy(atoms, atomNumber, boxSize); 

            //save old atom
            int atomIndex = int(randomFloat(0, atomNumber));
            Atom oldAtom = atoms[atomIndex];

            const double deltaX = randomFloat(-maxTranslation, maxTranslation);
            const double deltaY = randomFloat(-maxTranslation, maxTranslation);
            const double deltaZ = randomFloat(-maxTranslation, maxTranslation);

            atoms[atomIndex] = createAtom((unsigned long) atomIndex, oldAtom.x +
                    deltaX, oldAtom.y + deltaY, oldAtom.z + deltaZ, sigma, epsilon);

            double newEnergy = calculateEnergy(atoms, atomNumber, boxSize);
            //printf("newEnergy: %f\n", newEnergy);
            //printf("oldEnergy: %f\n", oldEnergy);
            bool accept = false;
            
            if(newEnergy < oldEnergy){
                accept = true;
            }
            else{
                double x = exp(-(newEnergy - oldEnergy) / kT);
                
                if(x >= randomFloat(0.0, 1.0)){
                    accept = true;
                }
                else{
                    accept = false;
                }
            }

            if(accept){
                accepted++;
            }
            else{
                rejected++;
                atoms[atomIndex] = oldAtom;
            }
            //printf("old energy: %f\nnew energy: %f\n", oldEnergy, newEnergy);
            //printf("accepted: %d\nrejected: %d\n\n", accepted, rejected);
        }

        endTime = clock();

        double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;
        printf("Time = %f seconds.\n\n", diffTime);

    }
}
