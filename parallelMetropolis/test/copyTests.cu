#include "copyTests.cuh"

void testCopyMolecules(){
    Molecule molec;
    size_t molecSize = sizeof(Molecule);

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
    
    molec = createMolecule(1,
            atoms, angles, bonds, dihedrals, hops,
            atomCount, angleCount, bondCount, dihedralCount, hopCount);

    //start cuda-ing
    Molecule *molec_d;
    printf("Allocating on the device.\n");
    allocateOnDevice(molec_d, &molec);

    printf("Copying to the device\n");
    moleculeDeepCopyToDevice(molec_d, &molec);

    printf("Copying from device\n");
    Molecule molec2;
    /*molec2.numOfAtoms = atomCount;
    molec2.numOfAngles = angleCount;
    molec2.numOfBonds = bondCount;
    molec2.numOfDihedrals = dihedralCount;
    molec2.numOfHops = hopCount;*/

    molec2.atoms = (Atom *)malloc(atomSize);
    molec2.angles = (Angle *)malloc(angleSize);
    molec2.bonds = (Bond *)malloc(bondSize);
    molec2.dihedrals = (Dihedral *)malloc(dihedralSize);
    molec2.hops = (Hop *)malloc(hopSize);
    moleculeDeepCopyToHost(&molec2, molec_d);

    printf("molec.id = %d, molec2.id = %d\n", molec.id, molec2.id);
    assert(molec.id == molec2.id);
    assert(molec.numOfAtoms == molec2.numOfAtoms);
    assert(molec.numOfBonds == molec2.numOfBonds);
    assert(molec.numOfAngles == molec2.numOfAngles);
    assert(molec.numOfDihedrals == molec2.numOfDihedrals);
    assert(molec.numOfHops == molec2.numOfHops);
}

void testAllocateMemory(){
    //TODO
}

void testFreeMemory(){
    //TODO
}
