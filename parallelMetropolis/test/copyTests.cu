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
    moleculeDeepCopyToDevice(molec_d, molecs, numOfMolec);

    /******
      Tests and assert statements.
    ******/

    printf("molec2.id = %d,  before copy\n", molec2.id);
    cudaMemcpy(&molec2, molec_d, sizeof(Molecule), cudaMemcpyDeviceToHost); 
    //moleculeDeepCopyToHost(&molec2, molec_d);

    printf("molec.id = %d, molec2.id = %d\n", molecs[0].id, molec2.id);
    printf("molec.numOfAtoms = %d, molec2.numOfAtoms = %d\n", molecs[0].numOfAtoms, molec2.numOfAtoms);
    printf("molec.numOfBonds = %d, molec2.numOfBonds = %d\n", molecs[0].numOfBonds, molec2.numOfBonds);
    printf("molec.numOfAngles = %d, molec2.numOfAngles = %d\n", molecs[0].numOfAngles, molec2.numOfAngles);
    printf("molec.numOfDihedrals = %d, molec2.numOfDihedrals = %d\n", molecs[0].numOfDihedrals, molec2.numOfDihedrals);
    printf("molec.numOfHops = %d, molec2.numOfHops = %d\n", molecs[0].numOfHops, molec2.numOfHops);
    
    
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
    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molec = (Molecule *)malloc(molecSize);

    for(int i = 0; i < numOfMolecules; i++){
        printf("Creating %dth array.\n", i);
        Atom *atoms;
        int atomCount = 3;
        size_t atomSize = sizeof(Atom) * atomCount;
        atoms = (Atom *)malloc(atomSize);
        atoms[0] = createAtom(1, 1, 1, 1);
        atoms[1] = createAtom(2, 1, 1, 1);
        atoms[2] = createAtom(3, 1, 2, 3);

        Bond *bonds;
        int bondCount = 2;
        size_t bondSize = sizeof(Bond) * bondCount;
        bonds = (Bond *)malloc(bondSize);
        bonds[0] = createBond(1, 2, 1.2, false);
        bonds[1] = createBond(2, 3, 3.1, true);

        Angle *angles;
        int angleCount = 2;
        size_t angleSize = sizeof(Angle) * angleCount;
        angles = (Angle *)malloc(angleSize);
        angles[0] = createAngle(1, 2, 86, false);
        angles[1] = createAngle(1, 3, 180, true);


        Dihedral *dihedrals;
        int dihedralCount = 2;
        size_t dihedralSize = sizeof(Dihedral) * dihedralCount;
        dihedrals = (Dihedral *)malloc(dihedralSize);
        dihedrals[0] = createDihedral(1, 2, 65, true);
        dihedrals[1] = createDihedral(1, 3, 43, false);

        Hop *hops;
        int hopCount = 2;
        size_t hopSize = sizeof(Hop) * hopCount;
        hops = (Hop *)malloc(hopSize);
        hops[0] = createHop(1,2,1);
        hops[1] = createHop(2,3,1);
        
        molec[i] = createMolecule(i + 1,
                atoms, angles, bonds, dihedrals, hops,
                atomCount, angleCount, bondCount, dihedralCount, hopCount);

        free(atoms);
        free(angles);
        free(bonds);
        free(dihedrals);
        free(hops);
    }
    
    allocateOnDevice(molec_d, molec, numOfMolecules);
    printf("Allocated on Device\n");
}

void testFreeMemory(){
    //TODO
}
