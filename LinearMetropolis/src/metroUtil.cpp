#include "metroUtil.h"

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = sigma;
    atom.epsilon = epsilon;

    return atom;
}
Atom createAtom(unsigned long id, double x, double y, double z){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = 0.0;
    atom.epsilon = 0.0;

    return atom;
}

double makePeriodic(double x, const double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

double calculateEnergy(Atom *atoms, int atomNumber, double boxSize){
    double sigma = atoms[0].sigma;
    double epsilon = atoms[0].epsilon;

    double totalEnergy = 0;

    int i;
    for(i = 0; i < atomNumber - 1; i++){
        int j;
        for(j = i + 1; j < atomNumber; j++){
            double deltaX = atoms[i].x - atoms[j].x;
            double deltaY = atoms[i].y - atoms[j].y;
            double deltaZ = atoms[i].z - atoms[j].z;

  //          //printf("preperiodic: %f\n", deltaX);
            deltaX = makePeriodic(deltaX, boxSize);
            deltaY = makePeriodic(deltaY, boxSize);
            deltaZ = makePeriodic(deltaZ, boxSize);
//            //printf("postperiodic: %f\n\n", deltaX);

            const double r2 = (deltaX * deltaX) +
                              (deltaY * deltaY) + 
                              (deltaZ * deltaZ);

            const double sig2OverR2 = pow(sigma, 2) / r2;
            const double sig6OverR6 = pow(sig2OverR2, 3);
            const double sig12OverR12 = pow(sig6OverR6, 2);
            //printf("%f\n", sig2OverR2); 
            //printf("%f\n", sig6OverR6);
            //printf("%f\n", sig12OverR12);
            const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
            //printf("%f\n", energy);
            totalEnergy += energy;
            //printf("%f\n\n", totalEnergy); 
        }
    }
    return totalEnergy;

}
