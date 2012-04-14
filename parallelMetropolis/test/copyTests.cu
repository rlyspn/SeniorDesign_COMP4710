#include "copyTests.cuh"

void testCopyMolecules(){
    printf("testCopyMolecules()\n");
    //Molecules on host
    Molecule *molecs;
    Molecule *copiedMolecs;
    
    int numOfMolecules = 3;
    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molecs = (Molecule *)malloc(molecSize);
    copiedMolecs = (Molecule *)malloc(molecSize);

    int angleCount = 2;
    int dihedralCount = 2;
    int bondCount = 2;
    int atomCount = 3;
    int hopCount = 2;

    for(int i = 0; i < numOfMolecules; i++){
        printf("Creating molecule %d\n.", i);
        
        Molecule m = molecs[i];
        
        size_t atomSize = sizeof(Atom) * atomCount;
        copiedMolecs[i].atoms = (Atom *)malloc(atomSize);

        m.atoms = (Atom *)malloc(atomSize);
        m.atoms[0] = createAtom(1, 1, 1, 1);
        m.atoms[1] = createAtom(2, 2, 2, 2);
        m.atoms[2] = createAtom(3, 3, 3, 3);

        size_t bondSize = sizeof(Bond) * bondCount;
        copiedMolecs[i].bonds = (Bond *)malloc(bondSize);

        m.bonds = (Bond *)malloc(bondSize);
        m.bonds[0] = createBond(1, 2, 1.2, false);
        m.bonds[1] = createBond(2, 3, 3.1, true);

        size_t angleSize = sizeof(Angle) * angleCount;
        copiedMolecs[i].angles = (Angle *)malloc(angleSize);

        m.angles = (Angle *)malloc(angleSize);
        m.angles[0] = createAngle(1, 2, 86, false);
        m.angles[1] = createAngle(1, 3, 180, true);

        size_t dihedralSize = sizeof(Dihedral) * dihedralCount;
        copiedMolecs[i].dihedrals = (Dihedral *)malloc(dihedralSize);

        m.dihedrals = (Dihedral *)malloc(dihedralSize);
        m.dihedrals[0] = createDihedral(1, 2, 65, true);
        m.dihedrals[1] = createDihedral(1, 3, 43, false);

        size_t hopSize = sizeof(Hop) * hopCount;
        copiedMolecs[i].hops = (Hop *)malloc(hopSize);
        m.hops = (Hop *)malloc(hopSize);
        
        m.hops[0] = createHop(1,3,2);
        m.hops[1] = createHop(2,3,1);
        
        m.id = i * atomCount + 57;
        m.numOfAtoms = atomCount;
        m.numOfBonds = bondCount;
        m.numOfAngles = angleCount;
        m.numOfDihedrals = dihedralCount;
        m.numOfHops = hopCount;

        molecs[i] = m;
        printMolecule(&molecs[i]);
    }

    printf("Testing deep copy to device.\n");
    size_t atomSize = numOfMolecules * atomCount * sizeof(Atom);
    size_t bondSize = numOfMolecules * bondCount * sizeof(Bond);
    size_t angleSize = numOfMolecules * angleCount * sizeof(Angle);
    size_t dihedralSize = numOfMolecules * dihedralCount * sizeof(Dihedral);
    size_t hopSize = numOfMolecules * hopCount * sizeof(Hop);
    size_t deviceMolecSize = numOfMolecules * sizeof(DeviceMolecule);

    Atom *atoms_d;
    Bond *bonds_d;
    Angle *angles_d;
    Dihedral *dihedrals_d;
    Hop *hops_d;
    DeviceMolecule *molec_d;
/*
    allocateOnDevice(molecs, molec_d, numOfMolecules, atoms_d, bonds_d, 
           angles_d, dihedrals_d, hops_d);
*/
    printf("molecs[0].atoms[0].x = %f\n", molecs[0].atoms[0].x);
 
    cudaMalloc((void **) &molec_d, deviceMolecSize);
    cudaMalloc((void **) &atoms_d, atomSize);
    cudaMalloc((void **) &bonds_d, bondSize);
    cudaMalloc((void **) &angles_d, angleSize);
    cudaMalloc((void **) &dihedrals_d, dihedralSize);
    cudaMalloc((void **) &hops_d, hopSize);
    
    moleculeDeepCopyToDevice(molec_d, molecs, numOfMolecules, atoms_d,
            bonds_d, angles_d, dihedrals_d, hops_d);

    moleculeDeepCopyToHost(copiedMolecs, molec_d, numOfMolecules, atoms_d, bonds_d,
            angles_d, dihedrals_d, hops_d);

    printf("Testing molecules.\n");
    for(int i = 0; i < numOfMolecules; i++){
        Molecule m = molecs[i];
        Molecule dm = copiedMolecs[i];
        
        printf("id = %d, %d\n", dm.id, m.id);
        printf("numOfAtoms = %d, %d\n", dm.numOfAtoms, m.numOfAtoms);
        printf("numOfBonds = %d, %d\n", dm.numOfBonds, m.numOfBonds);
        printf("numOfAngles = %d, %d\n", dm.numOfAngles, m.numOfAngles);
        printf("numOfDihedrals = %d, %d\n", dm.numOfDihedrals, m.numOfDihedrals);
        printf("numOfHops = %d, %d\n", dm.numOfHops, m.numOfHops);
        assert(dm.id == m.id);
        assert(dm.numOfAtoms == m.numOfAtoms);
        assert(dm.numOfBonds == m.numOfBonds);
        assert(dm.numOfAngles == m.numOfAngles);
        assert(dm.numOfDihedrals == m.numOfDihedrals);
        assert(dm.numOfHops == m.numOfHops);
        
        printf("Atoms: \n");
        for(int j = 0; j < copiedMolecs[i].numOfAtoms; j++){
            Atom a1 = copiedMolecs[i].atoms[j];
            Atom a2 = molecs[i].atoms[j];
            printf("id = %d, %d\n", a1.id, a2.id);
            printf("x = %f, %f\n", a1.x, a2.x);
            printf("y = %f, %f\n", a1.y, a2.y);
            printf("z = %f, %f\n", a1.z, a2.z);

        }
        printf("Bonds: \n");
        for(int j = 0; j < copiedMolecs[i].numOfBonds; j++){
            Bond a1 = copiedMolecs[i].bonds[j];
            Bond a2 = molecs[i].bonds[j];
            printf("atom1 = %d, %d\n", a1.atom1, a2.atom1);
            printf("atom2 = %d, %d\n", a1.atom2, a2.atom2);
            printf("distance = %f, %f\n", a1.distance, a2.distance);

        }
        printf("Angles: \n");
        for(int j = 0; j < copiedMolecs[i].numOfAngles; j++){
            Angle a1 = copiedMolecs[i].angles[j];
            Angle a2 = molecs[i].angles[j];
            printf("atom1 = %d, %d\n", a1.atom1, a2.atom1);
            printf("atom2 = %d, %d\n", a1.atom2, a2.atom2);
            printf("value = %f, %f\n", a1.value, a2.value);

        }
        printf("Dihedrals: \n");
        for(int j = 0; j < copiedMolecs[i].numOfDihedrals; j++){
            Dihedral a1 = copiedMolecs[i].dihedrals[j];
            Dihedral a2 = molecs[i].dihedrals[j];
            printf("atom1 = %d, %d\n", a1.atom1, a2.atom1);
            printf("atom2 = %d, %d\n", a1.atom2, a2.atom2);
            printf("value = %f, %f\n", a1.value, a2.value);

        }
        printf("Hops: \n");
        for(int j = 0; j < copiedMolecs[i].numOfBonds; j++){
            Hop a1 = copiedMolecs[i].hops[j];
            Hop a2 = molecs[i].hops[j];
            printf("atom1 = %d, %d\n", a1.atom1, a2.atom1);
            printf("atom2 = %d, %d\n", a1.atom2, a2.atom2);
            printf("hop = %d, %d\n", a1.hop, a2.hop);

        }

    }
    printf("testCopyMolecules completed successfully.\n");
}

void testAllocateMemory(){
    int numOfMolecules = 3;

    Molecule *molec;
    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molec = (Molecule *)malloc(molecSize);
   

    for(int i = 0; i < numOfMolecules; i++){
        printf("Creating atom %d\n", i);

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
    
    Atom *atoms_d;
    Bond *bonds_d;
    Angle *angles_d;
    Dihedral *dihedrals_d;
    Hop *hops_d;
    DeviceMolecule *molec_d;
    
    allocateOnDevice(molec, molec_d, numOfMolecules, atoms_d, bonds_d, 
           angles_d, dihedrals_d, hops_d);
    
    printf("molec_d = %p\n", molec_d);
    printf("atoms_d = %p\n", atoms_d);
    printf("bonds_d = %p\n", bonds_d);
    printf("angles_d = %p\n", angles_d);
    printf("dihedrals_d = %p\n", dihedrals_d);
    printf("hops_d = %p\n", hops_d);

    assert(molec_d != 0);
    assert(atoms_d != 0);
    assert(bonds_d != 0);
    assert(angles_d != 0);
    assert(dihedrals_d != 0);
    assert(hops_d != 0);
    
    printf("allocateOnDevice passed tests.\n");
}

void testFreeMemory(){
    //TODO
}
