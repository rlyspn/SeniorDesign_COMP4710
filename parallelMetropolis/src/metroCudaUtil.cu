#include "metroCudaUtil.cuh"


__device__ int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv;
}

__device__ int getYFromIndex(int x, int idx){
    return idx - x;
}

__device__ double makePeriodic(double x, const double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

__device__ double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    double sigma = atom1.sigma;
    double epsilon = atom2.epsilon;
    
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    deltaX = makePeriodic(deltaX, enviro.x);
    deltaY = makePeriodic(deltaY, enviro.y);
    deltaZ = makePeriodic(deltaZ, enviro.z);

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

    return energy;
}
