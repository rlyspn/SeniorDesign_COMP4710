#include "metroParallelUtil.h"


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

    return atom;
}

Atom createAtom(string name, unsigned long id, double x, double y, double z){
    Atom atom;
    atom.id = id;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.sigma = 0.0;
    atom.epsilon = 0.0;
    atom.name = name;

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

//returns a instance of a bond structure
Bond createBond(Atom *atom1, Atom *atom2, int bondCount){
    Bond bond;
    bond.bondCount = bondCount;
    bond.atom1 = atom1;
    bond.atom2 = atom2;
    return bond;
}

//returns an instance of a molecule object
Molecule createMolecule(string name, Atom *atoms, Bond *bonds, int atomCount, int bondCount){
    Molecule molecule;
    molecule.name = name;
    molecule.atoms = atoms;
    molecule.bonds = bonds;
    molecule.atomCount = atomCount;
    molecule.bondCount = bondCount;

    return molecule;
}


//utility to print off Atom coordinate data
void printAtoms(Atom *atoms, int count){
    for(int i = 0; i < count; i++){
        printf("%f, %f, %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

// writes atom information to a file.  Based on the siremol files.
void writeOutAtoms(Atom *atoms, Environment *enviro, string filename){
    FILE *atomOutput = fopen(filename.c_str(), "w");

    for(int i = 0; i < enviro->numOfAtoms; i++){
       fprintf(atomOutput, "ATOM  %5d  %s   %s    1    %8.3f%8.3f%8.3f  1.00  0.00          %s\n",
               i+1, atoms[i].name.c_str(), atoms[i].name.c_str(), atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].name.c_str());
    }
    fclose(atomOutput);
}

