#include "copyTests.cuh"

void testCopyMolecules(){
    printf("testCopyMolecules()\n");
    //Molecules on host
    Molecule *molecs;
    //molecules on device
    Molecule *molec_d;
    //Used to test after deepCopyToDevice
    Molecule *copiedMolecs;
    //molecules that have been deep copied back to the host
    Molecule *deepMolecs;
    
    int numOfMolecules = 3;
    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molecs = (Molecule *)malloc(molecSize);
    copiedMolecs = (Molecule *)malloc(molecSize);
    deepMolecs   = (Molecule *)malloc(molecSize);

    for(int i = 0; i < numOfMolecules; i++){
        printf("Creating atom %d\n.", i);
        
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
        
        m.id = i + 1;
        m.numOfAtoms = atomCount;
        m.numOfBonds = bondCount;
        m.numOfAngles = angleCount;
        m.numOfDihedrals = dihedralCount;
        m.numOfHops = hopCount;

        molecs[i] = m;
    }

    //start cuda-ing
    Molecule molec2;
    printf("Allocating on the device.\n");

    molec_d = allocateOnDevice(molec_d, molecs, numOfMolecules);
    //cudaMemcpy(copiedMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost);
    Molecule *m = (Molecule *)malloc(molecSize);
    cudaMemcpy(m, molec_d, molecSize, cudaMemcpyDeviceToHost);
    printf("atoms = %d\n", m[0].atoms);
    printf("bonds = %d\n", m[0].bonds);
    printf("hops = %d\n", m[0].hops);
    printf("Copying to the device\n");
    
    moleculeDeepCopyToDevice(molec_d, molecs, numOfMolecules);
    
    /******
      Tests and assert statements.
    ******/
    printf("Copying %d bytes from device.\n", molecSize); 
    cudaMemcpy(deepMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost);
   
    printf("atoms = %d\n", copiedMolecs[0].atoms);
    printf("bonds = %d\n", copiedMolecs[0].bonds);
    printf("hops = %d\n", copiedMolecs[0].hops);

    printf("After deep copy to device.\n");
    printf("molec.id = %d, copiedMolecs[0].id = %d\n", molecs[0].id, copiedMolecs[0].id);
    printf("molec.numOfAtoms = %d, copiedMolecs[0].numOfAtoms = %d\n", molecs[0].numOfAtoms, copiedMolecs[0].numOfAtoms);
    printf("molec.numOfBonds = %d, copiedMolecs[0].numOfBonds = %d\n", molecs[0].numOfBonds, copiedMolecs[0].numOfBonds);
    printf("molec.numOfAngles = %d, copiedMolecs[0].numOfAngles = %d\n", molecs[0].numOfAngles, copiedMolecs[0].numOfAngles);
    printf("molec.numOfDihedrals = %d, copiedMolecs[0].numOfDihedrals = %d\n", molecs[0].numOfDihedrals, copiedMolecs[0].numOfDihedrals);
    printf("molec.numOfHops = %d, copiedMolecs[0].numOfHops = %d\n", molecs[0].numOfHops, copiedMolecs[0].numOfHops);
    
    for(int i = 0; i < numOfMolecules; i++){
        printf("Testing %dth molecule after copy to device.\n", i);

        assert(molecs[i].id == copiedMolecs[i].id);
        assert(molecs[i].numOfAtoms == copiedMolecs[i].numOfAtoms);
        assert(molecs[i].numOfBonds == copiedMolecs[i].numOfBonds);
        assert(molecs[i].numOfAngles == copiedMolecs[i].numOfAngles);
        assert(molecs[i].numOfDihedrals == copiedMolecs[i].numOfDihedrals);
        assert(molecs[i].numOfHops == copiedMolecs[i].numOfHops);
    }

    //moleculeDeepCopyToHost(deepMolecs, molec_d, numOfMolecules);

    printf("After deep copy to host\n");
    printf("molec.id = %d, deepMolecs[0].id = %d\n", molecs[0].id, deepMolecs[0].id);
    printf("molec.numOfAtoms = %d, deepMolecs[0].numOfAtoms = %d\n", molecs[0].numOfAtoms, deepMolecs[0].numOfAtoms);
    printf("molec.numOfBonds = %d, deepMolecs[0].numOfBonds = %d\n", molecs[0].numOfBonds, deepMolecs[0].numOfBonds);
    printf("molec.numOfAngles = %d, deepMolecs[0].numOfAngles = %d\n", molecs[0].numOfAngles, deepMolecs[0].numOfAngles);
    printf("molec.numOfDihedrals = %d, deepMolecs[0].numOfDihedrals = %d\n", molecs[0].numOfDihedrals, deepMolecs[0].numOfDihedrals);
    printf("molec.numOfHops = %d, deepMolecs[0].numOfHops = %d\n", molecs[0].numOfHops, deepMolecs[0].numOfHops);
    
    for(int i = 0; i < numOfMolecules; i++){
        printf("Testing %dth molecule.\n", i);

        assert(molecs[i].id == deepMolecs[i].id);
        assert(molecs[i].numOfAtoms == deepMolecs[i].numOfAtoms);
        assert(molecs[i].numOfBonds == deepMolecs[i].numOfBonds);
        assert(molecs[i].numOfAngles == deepMolecs[i].numOfAngles);
        assert(molecs[i].numOfDihedrals == deepMolecs[i].numOfDihedrals);
        assert(molecs[i].numOfHops == deepMolecs[i].numOfHops);
    }


}

void testAllocateMemory(){
    int numOfMolecules = 3;

    Molecule *molec;
    Molecule *molec_d;
    Molecule *copiedMolecs;

    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molec = (Molecule *)malloc(molecSize);
    copiedMolecs = (Molecule *)malloc(molecSize);
    
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
    printf("molec_d = %d before.\n", molec_d); 
    molec_d = allocateOnDevice(molec_d, molec, numOfMolecules);
    
    printf("molec_d = %d before.\n", molec_d); 
    
    cudaMemcpy(copiedMolecs, molec_d, molecSize, cudaMemcpyDeviceToHost);

    printf("atomAddr = %d\n", copiedMolecs[0].atoms);
    printf("bondAddr = %d\n", copiedMolecs[0].angles);
    printf("angleAddr = %d\n", copiedMolecs[0].angles);

    printf("molec.id = %d, copiedMolecs[0].id = %d\n", molec[0].id, copiedMolecs[0].id);
    printf("molec.numOfAtoms = %d, copiedMolecs[0].numOfAtoms = %d\n", molec[0].numOfAtoms, copiedMolecs[0].numOfAtoms);
    printf("molec.numOfBonds = %d, copiedMolecs[0].numOfBonds = %d\n", molec[0].numOfBonds, copiedMolecs[0].numOfBonds);
    printf("molec.numOfAngles = %d, copiedMolecs[0].numOfAngles = %d\n", molec[0].numOfAngles, copiedMolecs[0].numOfAngles);
    printf("molec.numOfDihedrals = %d, copiedMolecs[0].numOfDihedrals = %d\n", molec[0].numOfDihedrals, copiedMolecs[0].numOfDihedrals);
    printf("molec.numOfHops = %d, copiedMolecs[0].numOfHops = %d\n", molec[0].numOfHops, copiedMolecs[0].numOfHops);
    
    for(int i = 0; i < numOfMolecules; i++){
        printf("Testing %dth molecule.\n", i);
        assert(molec[i].id == copiedMolecs[i].id);
        assert(molec[i].numOfAtoms == copiedMolecs[i].numOfAtoms);
        assert(molec[i].numOfBonds == copiedMolecs[i].numOfBonds);
        assert(molec[i].numOfAngles == copiedMolecs[i].numOfAngles);
        assert(molec[i].numOfDihedrals == copiedMolecs[i].numOfDihedrals);
        assert(molec[i].numOfHops == copiedMolecs[i].numOfHops);
    
        assert(molec[i].atoms != 0);
        assert(molec[i].bonds != 0);
        assert(molec[i].angles != 0);
        assert(molec[i].dihedrals != 0);
        assert(molec[i].hops != 0);
    }


    printf("allocateOnArray finished.\n");
}

void testFreeMemory(){
    //TODO
}
