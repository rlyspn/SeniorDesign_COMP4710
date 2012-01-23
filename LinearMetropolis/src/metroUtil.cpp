#include "metroUtil.h"

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = sigma;
    atom.epsilon = epsilon;
}
Atom createAtom(unsigned long id, double x, double y, double z){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = 0.0;
    atom.epsilon = 0.0;
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
        for(j = ++i; j < atomNumber; j++){
            double deltaX = atoms[i].x - atoms[j].x;
            double deltaY = atoms[i].y - atoms[j].y;
            double deltaZ = atoms[i].z - atoms[j].z;

            deltaX = makePeriodic(deltaX, boxSize);
            deltaY = makePeriodic(deltaY, boxSize);
            deltaZ = makePeriodic(deltaZ, boxSize);

            const double r2 = (deltaX * deltaX) +
                              (deltaY * deltaY) + 
                              (deltaZ * deltaZ);

            const double sig2OverR2 = pow(sigma, 2) / r2;
            const double sig6OverR6 = pow(sig2OverR2, 3);
            const double sig12OverR12 = pow(sig6OverR6, 2);

            const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);

            totalEnergy += energy;

            
        }
    }
    return totalEnergy;

}
