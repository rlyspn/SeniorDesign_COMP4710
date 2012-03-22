#include "metroUtil.h"


//create an instance of an Atom struct
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = sigma;
    atom.epsilon = epsilon;
    atom.charge = charge;

    return atom;
}
//create an instance of an Atom struct
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = sigma;
    atom.epsilon = epsilon;

    return atom;
}
//create an instance of an Atom struct w/o sigma/epsilon constants
Atom createAtom(unsigned long id, double x, double y, double z){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = 0.0;
    atom.epsilon = 0.0;
    atom.charge = 0.0;

    return atom;
}


//create an instance of the Environment struct
Environment createEnvironment(double x, double y, double z, double maxTrans, double temp, int numOfAtoms){
    Environment enviro;
    enviro.x = x;
    enviro.y = y;
    enviro.z = z;
    enviro.maxTranslation = maxTrans;
    enviro.temperature = temp;
    enviro.numOfAtoms = numOfAtoms;

    return enviro;
}

//returns an instance of a molecule object
Molecule createMolecule(int id, 
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals, 
                        int atomCount, int angleCount, int bondCount, int dihedralCount){
    Molecule molecule;
    molecule.id = id;

    molecule.atoms = atoms;
    molecule.angles = angles;
    molecule.bonds = bonds;
    molecule.dihedrals = dihedrals;

    molecule.numOfAtoms = atomCount;
    molecule.numOfAngles = angleCount;
    molecule.numOfBonds = bondCount;
    molecule.numOfDihedrals = dihedralCount;

    return molecule;
}

//returns an instance of a molecule object
Molecule createMolecule(int id, 
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals, Hop *hops, 
                        int atomCount, int angleCount, int bondCount, int dihedralCount, int hopCount){
    Molecule molecule;
    molecule.id = id;

    molecule.atoms = atoms;
    molecule.angles = angles;
    molecule.bonds = bonds;
    molecule.dihedrals = dihedrals;
	 molecule.hops = hops;

    molecule.numOfAtoms = atomCount;
    molecule.numOfAngles = angleCount;
    molecule.numOfBonds = bondCount;
    molecule.numOfDihedrals = dihedralCount;
	 molecule.numOfHops = hopCount; 	

    return molecule;
}


// returns an instance of the molecule struct
Molecule createMolecule(int id, Atom *atoms, int atomCount){
    Molecule molecule;
    molecule.id = id;
    molecule.atoms = atoms;
    molecule.numOfAtoms = atomCount;

    return molecule;
}

//returns an instance of a bond
Bond createBond(int atom1, int atom2, double distance, bool variable){
    Bond bond;
    bond.atom1 = atom1;
    bond.atom2 = atom2;
    bond.distance = distance;
    bond.variable = variable;

    return bond;
}

//returns an instance of a dihedral
Dihedral createDihedral(int atom1, int atom2, double value, bool variable){
    Dihedral dihedral;

    dihedral.atom1 = atom1;
    dihedral.atom2 = atom2;
    dihedral.value = value;
    dihedral.variable = variable;

    return dihedral;
}

Angle createAngle(int atom1, int atom2, double value, bool variable){
    Angle angle;
    angle.atom1 = atom1;
    angle.atom2 = atom2;
    angle.value = value;
    angle.variable = variable;

    return angle;
}

Hop createHop(int atom1, int atom2, int hop){
    Hop hops;
    hops.atom1 = atom1;
    hops.atom2 = atom2;
    hops.hop = hop;

    return hops;
}


//utility to print off Atom coordinate data
void printAtoms(Atom *atoms, int count){
    for(int i = 0; i < count; i++){
        printf("%f, %f, %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

//writes to the listed filename in the protein databank format
void writePDB(Atom *atoms, Environment enviro, string filename){
    ofstream outputFile;
    outputFile.open(filename.c_str());
    for(int i = 0; i < enviro.numOfAtoms; i++){
        Atom currentAtom = atoms[i];
        //ATOM number name residueName residueNumber chain x y z occupancy temp
        outputFile << "ATOM " << currentAtom.id << " NAME" << " residueName residueNumber chain "
            << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z << endl;
    }
    outputFile.close();
}

// writes atom information to a file.  Based on the siremol files.
void writeOutAtoms(Atom *atoms, Environment *enviro, string filename, int accepts, int rejects, double totalEnergy){
   ofstream outputFile;
   outputFile.open(filename.c_str());
   outputFile << "Total Energy: " << totalEnergy << endl;
   outputFile << "Acceptance Rate: " << (double)((double) accepts / (double) rejects) << endl;
   for(int i = 0; i < enviro->numOfAtoms; i++){
       Atom currentAtom = atoms[i];
       outputFile <<  currentAtom.id << " " << currentAtom.x << " " << currentAtom. y
           << " " << currentAtom.z << " " << endl;
   }
   outputFile.close();
}
/**
  Copies by value the values in molec2 into molec1
*/
void copyMolecule(Molecule *molec1, Molecule *molec2){
    molec1->atoms = (Atom *)malloc(sizeof(Atom) * molec2->numOfAtoms);
    molec1->bonds = (Bond *)malloc(sizeof(Bond) * molec2->numOfBonds);
    molec1->angles = (Angle *)malloc(sizeof(Angle) * molec2->numOfAngles);
    molec1->dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * molec2->numOfDihedrals);
	 molec1->hops = (Hop *)malloc(sizeof(Hop) * molec2->numOfHops);

    molec1->numOfAtoms = molec2->numOfAtoms;
    molec1->numOfBonds = molec2->numOfBonds;
    molec1->numOfAngles = molec2->numOfAngles;
    molec1->numOfDihedrals = molec2->numOfDihedrals;
	 molec1->numOfHops =  molec2->numOfHops;
    molec1->id = molec2->id;

    for(int i = 0; i < molec1->numOfAtoms; i++){
        molec1->atoms[i] = molec2->atoms[i];
    }
    for(int i = 0; i < molec1->numOfBonds; i++){
        molec1->bonds[i] = molec2->bonds[i];
    }
    for(int i = 0; i < molec1->numOfAngles; i++){
        molec1->angles[i] = molec2->angles[i];
    }
    for(int i = 0; i < molec1->numOfDihedrals; i++){
        molec1->dihedrals[i] = molec2->dihedrals[i];
    }
	 for(int i = 0; i < molec1->numOfHops; i++){
        molec1->hops[i] = molec2->hops[i];
    }



}


