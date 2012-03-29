#include "metroLinearUtil.h"

//calculates X (larger indexed atom) for energy calculation based on index in atom array
int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

//calculates Y (smaller indexed atom) for energy calculation based on index in atom array
int getYFromIndex(int x, int idx){
    return idx - (x * x - x) / 2;
}

//apply periodic boundaries
double makePeriodic(double x, double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

//keep coordinates with box
double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

void keepMoleculeInBox(Molecule *molecule, Environment *enviro){

    double maxX = 0.0;
    double maxY = 0.0;
    double maxZ = 0.0;

    double minX = 0.0;
    double minY = 0.0;
    double minZ = 0.0;

    //determine extreme boundaries for molecule
    for (int i = 0; i < molecule->numOfAtoms; i++){
        double currentX = molecule->atoms[i].x;
        double currentY = molecule->atoms[i].y;
        double currentZ = molecule->atoms[i].z;

        if (currentX > maxX)
           maxX = currentX;
        else if (currentX < minX)
           minX = currentX;

        if (currentY > maxY)
            maxY = currentY;
        else if (currentY < minY)
            minY = currentY;

        if (currentZ > maxZ)
            maxZ = currentZ;
        else if (currentZ < minZ)
            minZ = currentZ;
    
    }

    bool isFullyOutX = (minX > enviro->x || maxX < 0);
    bool isFullyOutY = (minY > enviro->y || maxY < 0);
    bool isFullyOutZ = (minZ > enviro->z || maxZ < 0);

    //for each axis, determine if the molecule escapes the environment 
    //and wrap it around into the environment
    for (int i = 0; i < molecule->numOfAtoms; i++){
        double* currentX = &(molecule->atoms[i].x);
        double* currentY = &(molecule->atoms[i].y);
        double* currentZ = &(molecule->atoms[i].z);
        if (maxX > enviro->x){
            if (!isFullyOutX)
                *currentX += (enviro->x - minX);
            *currentX = wrapBox(*currentX, enviro->x);
        }
        else if (minX < 0){
            if (!isFullyOutX)
                *currentX -= maxX;
            *currentX = wrapBox(*currentX, enviro->x);
        }

        if (maxY > enviro->y){
            if (!isFullyOutY)
                *currentY += (enviro->y - minY);
            *currentY = wrapBox(*currentY, enviro->y);
        }
        else if (minY < 0){
            if (!isFullyOutY)
                *currentY -= maxY;
            *currentY = wrapBox(*currentY, enviro->y);
        }

        if (maxZ > enviro->z){
            if (!isFullyOutZ)
                *currentZ += (enviro->z - minZ);
            *currentZ = wrapBox(*currentZ, enviro->z);
        }
        else if (minZ < 0){
            if (!isFullyOutZ)
                *currentZ -= maxZ;
            *currentZ = wrapBox(*currentZ, enviro->z);
        }
    }
}

//calculate Lennard-Jones energy between two atoms
double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    //store LJ constants locally
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro.x);
    deltaY = makePeriodic(deltaY, enviro.y);
    deltaZ = makePeriodic(deltaZ, enviro.z);

    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

    //calculate terms
    const double sig2OverR2 = pow(sigma, 2) / r2;
    const double sig6OverR6 = pow(sig2OverR2, 3);
    const double sig12OverR12 = pow(sig6OverR6, 2);
    const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
    
    if (r2 == 0.0){
        return 0.0;
    }
    else{
        return energy;
    }
}

void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro){
    //Translates each Molecule a random X,Y,and Z direction
	 //By translating every atom in that molecule by that translation

    //for each Molecule...
	 for(int i=0; i<enviro->numOfMolecules; i++){
	     for(int a=0; a<molec[i].numOfAtoms;a++){
		      Atom myAtom  =  molec[i].atoms[a];
		      myAtom.x =  dev_doublesX[i] * enviro->x + myAtom.x;
				myAtom.y =  dev_doublesY[i] * enviro->y + myAtom.y;
				myAtom.z =  dev_doublesZ[i] * enviro->z + myAtom.z;
		  }
		   keepMoleculeInBox(&molec[i],enviro);
    }
}

//generate coordinate data for the atoms
void generatePoints(Molecule *molec, Environment *enviro){
    double *devXDoubles = (double*)malloc(enviro->numOfMolecules * sizeof(double));
    double *devYDoubles = (double*)malloc(enviro->numOfMolecules * sizeof(double));
    double *devZDoubles = (double*)malloc(enviro->numOfMolecules * sizeof(double));
   
    srand((unsigned int) time(NULL));
	 //for each Molecule assaign a new XYZ
    for (int i = 0; i < enviro->numOfMolecules; i = i++){
        double newDouble = ((double) rand()) / ((double) (RAND_MAX));
        devXDoubles[i] = newDouble;
        newDouble = ((double) rand()) / ((double) (RAND_MAX));
        devYDoubles[i] = newDouble;
        newDouble = ((double) rand()) / ((double) (RAND_MAX));
        devZDoubles[i] = newDouble;
    }

    //assign the doubles to the coordinates
    assignAtomPositions(devXDoubles, devYDoubles, devZDoubles, molec, enviro);
}

//Calculates the energy of system using molecules
double calcEnergyWrapper(Molecule *molecules, Environment *enviro){
    
    Atom *atoms = (Atom *) malloc(sizeof(Atom) * enviro->numOfAtoms);
    int atomIndex = 0;
    for(int i = 0; i < enviro->numOfMolecules; i++){
        Molecule currentMolecule = molecules[i];
        for(int j = 0; j < currentMolecule.numOfAtoms; j++){
            atoms[atomIndex] = currentMolecule.atoms[j];
            atomIndex++;
        }
    }

    return calcEnergyWrapper(atoms, enviro, molecules);
}

double calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules){
    //setup CUDA storage
    double totalEnergy = 0.0;
    double *energySum_device;

    //calculate CUDA thread mgmt
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;	 
    size_t energySumSize = N * sizeof(double); 
	 
	double* energySum = (double*) malloc(energySumSize);

    calcEnergy(atoms, enviro, energySum);
    
    for(int i = 0; i < N; i++){

        int c = -2 * i;
        int discriminant = 1 - 4 * c;
        int qv = (-1 + sqrtf(discriminant)) / 2;
        int atomXid = qv + 1;
        
        int atomYid =  i - (atomXid * atomXid - atomXid) / 2;

        double xx = atoms[atomXid].x;
        double xy = atoms[atomXid].y;
        double xz = atoms[atomXid].z;

        double yx = atoms[atomYid].x;
        double yy = atoms[atomYid].y;
        double yz = atoms[atomYid].z;
        if (molecules != NULL){
            energySum[i] = energySum[i] * getFValue(atoms[atomXid], atoms[atomYid], molecules, enviro); 
        }
        totalEnergy += energySum[i];

    }
    
    return totalEnergy;
}

void calcEnergy(Atom *atoms, Environment *enviro, double *energySum){
    double lj_energy,charge_energy, fValue;

    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;
    for(int idx=0; idx<N; idx++){
    //calculate the x and y positions in the Atom array
        int xAtom_pos, yAtom_pos;
        xAtom_pos = getXFromIndex(idx);
        yAtom_pos = getYFromIndex(xAtom_pos, idx);

        Atom xAtom, yAtom;
        xAtom = atoms[xAtom_pos];
        yAtom = atoms[yAtom_pos];
        if (xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0){
            energySum[idx] = 0.0;
        }
        else{
            lj_energy = calc_lj(xAtom,yAtom,*enviro);
            charge_energy = calcCharge(xAtom, yAtom, enviro);
            fValue = 1.0;
            
            energySum[idx] = fValue * (lj_energy + charge_energy);
        }
	 }
}

double calcCharge(Atom atom1, Atom atom2, Environment *enviro){
    const double e = 1.602176565 * pow(10.f,-19.f);
 
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro->x);
    deltaY = makePeriodic(deltaY, enviro->y);
    deltaZ = makePeriodic(deltaZ, enviro->z);

    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);
    
    const double r = sqrt(r2);

    if (r == 0.0){
        return 0.0;
    }
    else{
        return (atom1.charge * atom2.charge * pow(e,2) / r);
    }
}

double calcBlending(double d1, double d2){
    return sqrt(d1 * d2);
}

//returns the molecule that contains a given atom
int getMoleculeFromAtomID(Atom a1, Molecule *molecules, Environment enviro){
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }
    return molecId;

}

double getFValue(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro){
    int m1 = getMoleculeFromAtomID(atom1, molecules, *enviro);
    int m2 = getMoleculeFromAtomID(atom2, molecules, *enviro);

    if(m1 != m2)
        return 1.0;
	 else if( hopGE3(atom1.id, atom2.id,molecules[m1]) )     
		  return 0.5;
	 else
		  return 0.0;
}

int hopGE3(int atom1, int atom2, Molecule molecule){
    for(int x=0; x< molecule.numOfHops; x++){
		      Hop myHop = molecule.hops[x];
				if(myHop.atom1==atom1 && myHop.atom2==atom2)
				    return 1;
	 }
	 return 0;
}

void rotateMolecule(Molecule molecule, Atom pivotAtom, double maxRotation){
    //save pivot atom coordinates because they will change
    double pivotAtomX = pivotAtom.x;
    double pivotAtomY = pivotAtom.y;
    double pivotAtomZ = pivotAtom.z;

    //translate entire molecule to place pivotAtom at origin
    for (int i = 0; i < molecule.numOfAtoms; i++){
        molecule.atoms[i].x -= pivotAtomX;
        molecule.atoms[i].y -= pivotAtomY;
        molecule.atoms[i].z -= pivotAtomZ;
    }

    srand(time(NULL));
    double dtr = PI / 180.0;

    //rotate molecule about origin
    for (int axis = 0; axis < 3; axis++){
        double rotation = ((double) rand() / (double) RAND_MAX) * maxRotation * dtr;
        double sinrot = sin(rotation);
        double cosrot = cos(rotation);
        if (axis == 0){ //rotate about x-axis
            for (int i = 0; i < molecule.numOfAtoms; i++){
                Atom *thisAtom = &(molecule.atoms[i]);
                double oldY = thisAtom->y;
                double oldZ = thisAtom->z;
                thisAtom->y = cosrot * oldY + sinrot * oldZ;
                thisAtom->z = cosrot * oldZ - sinrot * oldY;
            }
        }
        else if (axis == 1){ //rotate about y-axis
            for (int i = 0; i < molecule.numOfAtoms; i++){
                Atom *thisAtom = &(molecule.atoms[i]);
                double oldX = thisAtom->x;
                double oldZ = thisAtom->z;
                thisAtom->x = cosrot * oldX - sinrot * oldZ;
                thisAtom->z = cosrot * oldZ + sinrot * oldX;
            }
        }
        if (axis == 2){ //rotate about z-axis
            for (int i = 0; i < molecule.numOfAtoms; i++){
                Atom *thisAtom = &(molecule.atoms[i]);
                double oldX = thisAtom->x;
                double oldY = thisAtom->y;
                thisAtom->x = cosrot * oldX + sinrot * oldY;
                thisAtom->y = cosrot * oldY - sinrot * oldX;
            }
        }
    }

    //translate entire molecule back based on original pivot point
    for (int i = 0; i < molecule.numOfAtoms; i++){
        molecule.atoms[i].x += pivotAtomX;
        molecule.atoms[i].y += pivotAtomY;
        molecule.atoms[i].z += pivotAtomZ;
    }
}

/**
  This  is currently a stub pending information from Dr. Acevedo
*/
double solventAccessibleSurfaceArea(){
    return -1.f;
}

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double soluteSolventDistributionFunction(){
    return -1.f;
}

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double atomAtomDistributionFunction(){
    return -1.f;
}

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double solventSolventTotalEnergy(){
    return -1.f;
}

/**
  This is currently a stub pending information from Dr. Acevedo
*/
double soluteSolventTotalEnergy(){
    return -1.f;
}

