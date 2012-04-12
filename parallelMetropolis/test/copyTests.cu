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
        printf("Creating atom %d\n.", i);
        
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

    allocateOnDevice(molecs, molec_d, numOfMolecules, atoms_d, bonds_d, 
           angles_d, dihedrals_d, hops_d);

    moleculeDeepCopyToDevice(molec_d, molecs, numOfMolecules, atoms_d,
            bonds_d, angles_d, dihedrals_d, hops_d);
    moleculeDeepCopyToHost(copiedMolecs, molec_d, numOfMolecules, atoms_d, bonds_d,
            angles_d, dihedrals_d, hops_d);

    for(int i = 0; i < 1; i++){
        Molecule m = copiedMolecs[i];
        
        printf("id = %d", m.id);
        printf("numOfAtoms = %d\n", m.numOfAtoms);
        printf("numOfBonds = %d\n", m.numOfBonds);
        printf("numOfAngles =%d\n", m.numOfAngles);
        printf("numOfDihedrals = %d\n", m.numOfDihedrals);
        printf("numOfHops = %d\n", m.numOfHops);
    
    
    }
    /**
    DeviceMolecule *copiedDM = (DeviceMolecule *)malloc(deviceMolecSize);
    Atom *copiedAtoms = (Atom *)malloc(atomSize);
    Bond *copiedBonds = (Bond *)malloc(bondSize);
    Angle *copiedAngles = (Angle *)malloc(angleSize);
    Dihedral *copiedDihedrals = (Dihedral *)malloc(dihedralSize);
    Hop *copiedHops = (Hop *)malloc(hopSize);
    
    
    cudaMemcpy(copiedDM, molec_d, deviceMolecSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedAtoms, atoms_d, atomSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedBonds, bonds_d, bondSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedAngles, angles_d, angleSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedDihedrals, dihedrals_d, dihedralSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(copiedHops, hops_d, hopSize, cudaMemcpyDeviceToHost);

    for(int i = 0; i < numOfMolecules; i++){
        Molecule m = molecs[i];
        DeviceMolecule dm = copiedDM[i];
        
        printf("id = %d, %d\n", dm.id, m.id);
        printf("numOfAtoms = %d, %d\n", dm.numOfAtoms, m.numOfAtoms);
        printf("numOfBonds = %d, %d\n", dm.numOfBonds, m.numOfBonds);
        printf("numOfAngles = %d, %d\n", dm.numOfAngles, m.numOfAngles);
        printf("numOfDihedrals = %d, %d\n", dm.numOfDihedrals, m.numOfDihedrals);
        printf("numOfHops = %d, %d\n", dm.numOfHops, m.numOfHops);
        //assert(dm.id == m.id);
        assert(dm.numOfAtoms == m.numOfAtoms);
        assert(dm.numOfBonds == m.numOfBonds);
        assert(dm.numOfAngles == m.numOfAngles);
        assert(dm.numOfDihedrals == m.numOfDihedrals);
        assert(dm.numOfHops == m.numOfHops);
    }

    printf("Testing atoms.\n");

    int atomIndex = 0;
    for(int i = 0; i < numOfMolecules; i++){
        Molecule m = molecs[i];
        for(int j = 0; j < m.numOfAtoms; j++){
            Atom a1 = copiedAtoms[atomIndex];
            Atom a2 = molecs[i].atoms[j];
        printf("id = %d, %d\n", a1.id, a2.id);
        printf("x = %f, %f\n", a1.x, a2.x);
        printf("y = %f, %f\n", a1.y, a2.y);
        printf("z = %f, %f\n", a1.z, a2.z);
        atomIndex++;
        }
    }
    int bondIndex = 0;
    for(int i = 0; i < numOfMolecules; i++){
        Molecule m = molecs[i];
        for(int j = 0; j < m.numOfBonds; j++){
            Bond a1 = copiedBonds[bondIndex];
            Bond a2 = molecs[i].bonds[j];
            printf("atom1 = %d, %d\n", a1.atom1, a2.atom1);
            printf("atom2 = %d, %d\n", a1.atom2, a2.atom2);
            printf("dist  = %f, %f\n", a1.distance, a2.distance);
            
            bondIndex++;
        }
    }*/

}

void testAllocateMemory(){
    int numOfMolecules = 3;

    Molecule *molec;
    size_t molecSize = sizeof(Molecule) * numOfMolecules;
    molec = (Molecule *)malloc(molecSize);
   

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
    
    Atom *atoms_d;
    Bond *bonds_d;
    Angle *angles_d;
    Dihedral *dihedrals_d;
    Hop *hops_d;
    DeviceMolecule *molec_d;
    
    allocateOnDevice(molec, molec_d, numOfMolecules, atoms_d, bonds_d, 
           angles_d, dihedrals_d, hops_d);
    
    printf("molec_d = %d\n", molec_d);
    printf("atoms_d = %d\n", atoms_d);
    printf("bonds_d = %d\n", bonds_d);
    printf("angles_d = %d\n", angles_d);
    printf("dihedrals_d = %d\n", dihedrals_d);
    printf("hops_d = %d\n", hops_d);

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
