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
double calculate_energy(Atom *atoms, Environment *enviro, Molecule *molecules){
    int atomNumber = enviro->numOfAtoms;
    double sigma = atoms[0].sigma;
    double epsilon = atoms[0].epsilon;
    
    const double e = 1.602176565 * pow(10.f,-19.f);

    double totalEnergy = 0;
    
    int i;
    for(i = 0; i < atomNumber - 1; i++){
        int j;
        for(j = i + 1; j < atomNumber; j++){
            double deltaX = atoms[j].x - atoms[i].x;
            double deltaY = atoms[j].y - atoms[i].y;
            double deltaZ = atoms[j].z - atoms[i].z;
          
            deltaX = make_periodic(deltaX, enviro->x);
            deltaY = make_periodic(deltaY, enviro->y);
            deltaZ = make_periodic(deltaZ, enviro->z);

            const double r2 = (deltaX * deltaX) +
                              (deltaY * deltaY) + 
                              (deltaZ * deltaZ);

            const double r = sqrt(r2);

            const double sig2OverR2 = pow(sigma, 2) / r2;
            const double sig6OverR6 = pow(sig2OverR2, 3);
            const double sig12OverR12 = pow(sig6OverR6, 2);
            const double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);

            const double charge_energy = (atoms[i].charge * atoms[j].charge * pow(e,2) / r);
            double fValue = 1.0;
            if (molecules != NULL)
            {
                fValue = getFValueLinear(atoms[i], atoms[j], molecules, enviro);
            }
            
            int N = (j * j - j ) / 2 + i;

            //printf("energy[%d]: %f | r2 = %f ", N, energy, r2);
            //printf(" | atomX = {id:%d, x:%f, y:%f, z:%f} ", atoms[j].id, atoms[j].x, atoms[j].y, atoms[j].z);
            //printf(" | atomY = {id:%d, x:%f, y:%f, z:%f}\n", atoms[i].id, atoms[i].x, atoms[i].y, atoms[i].z);
            totalEnergy += fValue * (lj_energy + charge_energy);
        }
    }
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

double calc_r_value(Atom a1, Atom a2, Environment enviro){
    double dx = make_periodic(a1.x - a2.x, enviro.x);
    double dy = make_periodic(a1.y - a2.y, enviro.y);
    double dz = make_periodic(a1.z - a2.z, enviro.z);

    return sqrt(dx * dx + dy * dy + dz * dz);

}

double calc_charge(Atom a1, Atom a2, Environment enviro){
    double e = 1.602176565 * pow(10.f,-19.f);  

    return (a1.charge * a2.charge * pow(e, 2)) / calc_r_value(a1, a2, enviro);
}

//returns the molecule that contains a given atom
int getMoleculeFromIDLinear(Atom a1, Molecule *molecules, Environment enviro){
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }
    return molecId;

}

double getFValueLinear(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro){
    int m1 = getMoleculeFromIDLinear(atom1, molecules, *enviro);
    int m2 = getMoleculeFromIDLinear(atom2, molecules, *enviro);
    Molecule molec = molecules[0];
    for(int i = 0; i < enviro->numOfMolecules; i++){
        if(molecules[i].id == m1){
            molec = molecules[i];
            break;
        }
    }

    if(m1 != m2)
        return 1.0;
	 else if( hopGE3Linear(atom1.id, atom2.id, molecules[m1]) )     
		  return 0.5;
	 else
		  return 0.0;
}

int hopGE3Linear(int atom1, int atom2, Molecule molecule){
    for(int x=0; x< molecule.numOfHops; x++){
		      Hop myHop = molecule.hops[x];
				if(myHop.atom1==atom1 && myHop.atom2==atom2)
				    return 1;
	 }
	 return 0;
}
