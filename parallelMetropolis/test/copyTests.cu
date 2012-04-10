#include "copyTests.cuh"

void testCopyMolecules(){
    Molecule *molecs;
    int numOfMolec = 10;
    size_t molecSize = sizeof(Molecule) * numOfMolec;
    molecs = (Molecule *)malloc(molecSize);

    for(int i = 0; i < numOfMolec; i++){
        Molecule m = molecs[i];
        
        int atomCount = 3;
        size_t atomSize = sizeof(Atom) * atomCount;
        m.atoms = (Atom *)malloc(atomSize);
        m.atoms[0] = createAtom(1, 1, 1, 1);
        m.atoms[1] = createAtom(2, 1, 1, 1);
        m.atoms[2] = createAtom(3, 1, 2, 3);

        int bondCount = 2;
        size_t bondSize = sizeof(Bond) * bondCount;
        m.bonds = (Bond *)malloc(bondSize);
        m.bonds[0] = createBond(1, 2, 1.2, false);
        m.bonds[1] = createBond(2, 3, 3.1, true);

        int angleCount = 2;
        size_t angleSize = sizeof(Angle) * angleCount;
        m.angles = (Angle *)malloc(angleSize);
        m.angles[0] = createAngle(1, 2, 86, false);
        m.angles[1] = createAngle(1, 3, 180, true);

        int dihedralCount = 2;
        size_t dihedralSize = sizeof(Dihedral) * dihedralCount;
        m.dihedrals = (Dihedral *)malloc(dihedralSize);
        m.dihedrals[0] = createDihedral(1, 2, 65, true);
        m.dihedrals[1] = createDihedral(1, 3, 43, false);

        int hopCount = 2;
        size_t hopSize = sizeof(Hop) * hopCount;
        m.hops = (Hop *)malloc(hopSize);
        m.hops[0] = createHop(1,2,1);
        m.hops[1] = createHop(2,3,1);
        
        m.id = i;
        m.numOfAtoms = atomCount;
        m.numOfBonds = bondCount;
        m.numOfAngles = angleCount;
        m.numOfDihedrals = dihedralCount;
        m.numOfHops = hopCount;

        molecs[i] = m;
    }

    //start cuda-ing
    Molecule *molec_d;
    Molecule molec2;
    printf("Allocating on the device.\n");
    allocateOnDevice(molec_d, molecs, numOfMolec);

    printf("Copying to the device\n");
    //moleculeDeepCopyToDevice(molec_d, molecs, numOfMolec);
    //cudaMalloc((void**)&molec_d, molecSize);
    cudaMemcpy(molec_d, molecs, molecSize, cudaMemcpyHostToDevice);
    /******
      Tests and assert statements.
    ******/
    Molecule *copiedMolecs;
    copiedMolecs = (Molecule *)malloc(molecSize);
    cudaMemcpy(copiedMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost); 
    molec2 = copiedMolecs[0];
    //moleculeDeepCopyToHost(&molec2, molec_d);

    printf("molec.id = %d, copiedMolecs[0].id = %d\n", molecs[0].id, copiedMolecs[0].id);
    printf("molec.numOfAtoms = %d, copiedMolecs[0].numOfAtoms = %d\n", molecs[0].numOfAtoms, copiedMolecs[0].numOfAtoms);
    printf("molec.numOfBonds = %d, copiedMolecs[0].numOfBonds = %d\n", molecs[0].numOfBonds, copiedMolecs[0].numOfBonds);
    printf("molec.numOfAngles = %d, copiedMolecs[0].numOfAngles = %d\n", molecs[0].numOfAngles, copiedMolecs[0].numOfAngles);
    printf("molec.numOfDihedrals = %d, copiedMolecs[0].numOfDihedrals = %d\n", molecs[0].numOfDihedrals, copiedMolecs[0].numOfDihedrals);
    printf("molec.numOfHops = %d, copiedMolecs[0].numOfHops = %d\n", molecs[0].numOfHops, copiedMolecs[0].numOfHops);
    
    
    /*assert(molec.id == molec2.id);
    assert(molec.numOfAtoms == molec2.numOfAtoms);
    assert(molec.numOfBonds == molec2.numOfBonds);
    assert(molec.numOfAngles == molec2.numOfAngles);
    assert(molec.numOfDihedrals == molec2.numOfDihedrals);
    assert(molec.numOfHops == molec2.numOfHops);*/
}

void testAllocateMemory(){
    int numOfMolecules = 3;

    Molecule *molec;
    Molecule *molec_d;
    Molecule *copiedMolecs;

    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molec = (Molecule *)malloc(molecSize);
    copiedMolecs = (Molecule *)malloc(molecSize);
    cudaMalloc((void **) &molec_d, molecSize);
    
    for(int i = 0; i < numOfMolecules; i++){
        printf("Creating atom %d\n.", i);

        Molecule m = molec[i];
        
        int atomCount = 3;
        size_t atomSize = sizeof(Atom) * atomCount;
        m.atoms = (Atom *)malloc(atomSize);
        m.atoms[0] = createAtom(1, 1, 1, 1);
        m.atoms[1] = createAtom(2, 1, 1, 1);
        m.atoms[2] = createAtom(3, 1, 2, 3);

        int bondCount = 2;
        size_t bondSize = sizeof(Bond) * bondCount;
        m.bonds = (Bond *)malloc(bondSize);
        m.bonds[0] = createBond(1, 2, 1.2, false);
        m.bonds[1] = createBond(2, 3, 3.1, true);

        int angleCount = 2;
        size_t angleSize = sizeof(Angle) * angleCount;
        m.angles = (Angle *)malloc(angleSize);
        m.angles[0] = createAngle(1, 2, 86, false);
        m.angles[1] = createAngle(1, 3, 180, true);

        int dihedralCount = 2;
        size_t dihedralSize = sizeof(Dihedral) * dihedralCount;
        m.dihedrals = (Dihedral *)malloc(dihedralSize);
        m.dihedrals[0] = createDihedral(1, 2, 65, true);
        m.dihedrals[1] = createDihedral(1, 3, 43, false);

        int hopCount = 2;
        size_t hopSize = sizeof(Hop) * hopCount;
        m.hops = (Hop *)malloc(hopSize);
        m.hops[0] = createHop(1,2,1);
        m.hops[1] = createHop(2,3,1);
       
        m.id = i + 1;
        m.numOfAtoms = atomCount;
        m.numOfBonds = bondCount;
        m.numOfAngles = angleCount;
        m.numOfDihedrals = dihedralCount;
        m.numOfHops = hopCount;
    
        molec[i] = m;
    }
    printf("molecSize = %d\n", molecSize);    
    
    //TODO
    //allocateOnDevice(molec_d, molec, numOfMolecules);
    //Seg fault here
    cudaMemcpy(molec_d, molec, molecSize, cudaMemcpyHostToDevice);
    
    cudaMemcpy(copiedMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost);
    printf("molec.id = %d, copiedMolecs[0].id = %d\n", molec[0].id, copiedMolecs[0].id);
    printf("molec.numOfAtoms = %d, copiedMolecs[0].numOfAtoms = %d\n", molec[0].numOfAtoms, copiedMolecs[0].numOfAtoms);
    printf("molec.numOfBonds = %d, copiedMolecs[0].numOfBonds = %d\n", molec[0].numOfBonds, copiedMolecs[0].numOfBonds);
    printf("molec.numOfAngles = %d, copiedMolecs[0].numOfAngles = %d\n", molec[0].numOfAngles, copiedMolecs[0].numOfAngles);
    printf("molec.numOfDihedrals = %d, copiedMolecs[0].numOfDihedrals = %d\n", molec[0].numOfDihedrals, copiedMolecs[0].numOfDihedrals);
    printf("molec.numOfHops = %d, copiedMolecs[0].numOfHops = %d\n", molec[0].numOfHops, copiedMolecs[0].numOfHops);
    
    printf("allocateOnArray finished.\n");
}

void testFreeMemory(){
    //TODO
}
