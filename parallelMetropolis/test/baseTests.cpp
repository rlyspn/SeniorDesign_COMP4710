#include "baseTests.h"

double calculate_energy(double **coords,  int n_atoms,  double *box_size,
                         double sigma,  double epsilon){

    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double total_energy = 0;

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {
            double delta_x = coords[j][0] - coords[i][0];
            double delta_y = coords[j][1] - coords[i][1];
            double delta_z = coords[j][2] - coords[i][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

             double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);

            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
             double sig2_over_r2 = (sigma*sigma) / r2;
             double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
             double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

             double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

            total_energy = total_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return total_energy;
}

//Same as above but uses the new enviorment struct
double calculate_energy(Atom *atoms, Environment *enviro){
    int atomNumber = enviro->numOfAtoms;
    double sigma = atoms[0].sigma;
    double epsilon = atoms[0].epsilon;

    double totalEnergy = 0;
    
    int i;
    for(i = 0; i < atomNumber - 1; i++){
        int j;
        for(j = i + 1; j < atomNumber; j++){
            double deltaX = atoms[j].x - atoms[i].x;
            double deltaY = atoms[j].y - atoms[i].y;
            double deltaZ = atoms[j].z - atoms[i].z;
          
           // printf("J: %d\nI: %d\n", i, j);
            //printf("atoms[j].x: %f\n", atoms[j].x);
            //printf("atoms[i].x: %f\n", atoms[i].x);
            //printf("deltaX: %f\n", deltaX);
            //printf("deltaY: %f\n", deltaY);
            //printf("deltaZ: %f\n", deltaZ);


  //          //printf("preperiodic: %f\n", deltaX);
            deltaX = make_periodic(deltaX, enviro->x);
            deltaY = make_periodic(deltaY, enviro->y);
            deltaZ = make_periodic(deltaZ, enviro->z);
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
            printf("energy[%d]: %f\n",j, energy);
            totalEnergy += energy;
            //printf("%f\n\n", totalEnergy); 
            //printf("AtomNumber: %d\n", atomNumber);
        }
    }
    printf("number of atoms: %d", atomNumber);
    //printf("Sigma: %f\n", sigma);
    //printf("Epsilon: %f\n", epsilon);
    return totalEnergy;

}


// Subroutine to apply periodic boundaries
double make_periodic(double x,  double box)
{
    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }

    return x;
}

// Subroutine to wrap the coordinates into a box
double wrap_into_box(double x, double box)
{
    while (x > box)
    {
        x = x - box;
    }

    while (x < 0)
    {
        x = x + box;
    }

    return x;
}
long timevaldiff(struct timeval *starttime, struct timeval *finishtime)
{
    long msec;
    msec = (finishtime->tv_sec - starttime->tv_sec)*1000;
    msec += (finishtime->tv_usec - starttime->tv_usec)/1000;
    return msec;
}
