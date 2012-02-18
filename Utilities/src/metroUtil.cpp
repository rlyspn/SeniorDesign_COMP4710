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

//utility to print off Atom coordinate data
void printAtoms(Atom *atoms, int count){
    for(int i = 0; i < count; i++){
        printf("%f, %f, %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
    }
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

void printState(Environment *enviro, Molecule *molecules, int numOfMolecules, string filename){
    ofstream outFile;
    outFile.open(filename.c_str());
    //print the environment
    outFile << enviro->x << " " << enviro->y << " " << enviro->z << " " << enviro->numOfAtoms << endl;
    outFile << endl; // blank line
    for(int i = 0; i < numOfMolecules; i++){
        Molecule currentMol = molecules[i];
        outFile << currentMol.id << endl;
        outFile << "=" << endl; // delimiter
    
        //write atoms
        for(int j = 0; j < currentMol.numOfAtoms; j++){
            Atom currentAtom = currentMol.atoms[j];
            outFile << currentAtom.id << " "
                << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z
                << currentAtom.sigma << " " << currentAtom.epsilon << endl;
        }
        outFile << "=" << endl; // delimiter
        
        //write bonds
        for(int j = 0; j < currentMol.numOfBonds; j++){
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " " << currentBond.atom2 << " "
                << currentBond.distance << " ";
            if(currentBond.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;

        }

        for(int j = 0; j < currentMol.numOfDihedrals; j++){
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " " << currentDi.atom2 << " "
                << currentDi.value << " ";
            if(currentDi.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }

        //print angles
        for(int j = 0; j < currentMol.numOfAngles; j++){
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " " << currentAngle.atom2 << " "
                << currentAngle.value << " ";
            if(currentAngle.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }

        //write a == line
        outFile << "==" << endl;
    }
    outFile.close();
}

// Returns an Environment read from the line of the state file.
Environment getEnvironmentFromLine(string line){
    Environment enviro;
    //tokenize input line
    char *tokens; 
    char *charLine;
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    //extract data from line
    int tokenNumber = 0;
    double x, y, z;
    int numOfAtoms;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0:
                enviro.x = atof(tokens);
                break;
            case 1:
                enviro.y = atof(tokens);
                break;
            case 2:
                enviro.z = atof(tokens);
                break;
            case 3:
                enviro.numOfAtoms = atoi(tokens);
                break;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }

    return enviro;
}

// Returns an Atom read from a line of the state file
Atom getAtomFromLine(string line){
    Atom returnAtom = createAtom(-1, -1, -1, -1, -1, -1);
    char *tokens;
    char *charLine;
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    
    int tokenNumber = 0;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0: // id
                returnAtom.id = atoi(tokens);
                break;
            case 1: // x
                returnAtom.x = atof(tokens);
                break;
            case 2: // y
                returnAtom.y = atof(tokens);
                break;
            case 3: // z
                returnAtom.z = atof(tokens);
                break;
            case 4: // sigma
                returnAtom.sigma = atof(tokens);
                break;
            case 5: // epsilon
                returnAtom.epsilon = atof(tokens);
                break;
    
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return returnAtom;
}

// reads information about a bond from the line and returns a bond
Bond getBondFromLine(string line){
    Bond bond = createBond( -1, -1, -1, false);

    char *tokens;
    char *charLine;
    strcpy(charLine, line.c_str());
    int tokenNumber = 0;
    
    while(tokens != NULL){
        switch(tokenNumber){
            case 0: // atom 1
                bond.atom1 = atoi(tokens);
                break;
            case 1: // atom 2
                bond.atom2 = atoi(tokens);
                break;
            case 2: // value
                bond.distance = atof(tokens);
                break;
            case 3: // variable
                if(atoi(tokens) == 1)
                    bond.variable = true;
                else
                    bond.variable = false;
                break;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return bond;
}

Angle getAngleFromLine(string line){
    Angle angle = createAngle(-1, -1, -1, false);

    char *tokens;
    char *charLine;
    strcpy(charLine, line.c_str());
    int tokenNumber = 0;

    while(tokens != NULL){
       switch(tokenNumber){
           case 0:
               angle.atom1 = atoi(tokens);
               break;
            case 1:
               angle.atom2 = atoi(tokens);
               break;
            case 2:
               angle.value = atof(tokens);
               break;
            case 3:
               if(atoi(tokens) == 1)
                   angle.variable = true;
               else
                   angle.variable = false;
       }
       tokens = strtok(NULL, " ");
       tokenNumber++;
    }
    return angle;
}

Dihedral getDihedralFromLine(string line){
    Dihedral dihedral = createDihedral(-1, -1, -1, false);

    char *tokens;
    char *charLine;
    strcpy(charLine, line.c_str());
    int tokenNumber = 0;

    while(tokens != NULL){
        switch(tokenNumber){
            case 0:
                dihedral.atom1 = atoi(tokens);
                break;
            case 1:
                dihedral.atom2 = atoi(tokens);
                break;
            case 2:
                dihedral.value = atof(tokens);
                break;
            case3:
                if(atoi(tokens) == 1)
                    dihedral.variable = true;
                else
                    dihedral.variable = false;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return dihedral;

}

//Read in the environment from the provided file
Environment readInEnvironment(string filename){
    ifstream inFile (filename.c_str());
    string line;
    Environment enviro;

    if(inFile.is_open()){
        getline(inFile, line);
        enviro = getEnvironmentFromLine(line);    
    }
    else
        return enviro;

    return enviro;
        
}

Molecule readMoleculeSection(ifstream inFile){
    vector<Bond> bonds;
    vector<Angle> angles;
    vector<Atom> atoms;
    vector<Dihedral> dihedrals;
    int lines = 0;
    string line;

    while(inFile.good()){
        
    }

    Bond *bondArray = (Bond *) malloc(sizeof(Bond) * bonds.size());
    Angle *angleArray = (Angle *) malloc(sizeof(Angle) * angles.size());
    Atom *atomArray = (Atom *) malloc(sizeof(Atom) * atoms.size());
    Dihedral *dihedralArray = (Dihedral *) malloc(sizeof(Dihedral) * dihedrals.size());

    for(int i = 0; i < bonds.size(); i++)
        bondArray[i] = bonds[i];
    for(int i = 0; i < angles.size(); i++)
        angleArray[i] = angles[i];
    for(int i = 0; i < atoms.size(); i++)
        atomArray[i] = atoms[i];
    for(int i = 0; i < dihedrals.size(); i++)
        dihedralArray[i] = dihedrals[i];

    return createMolecule(0, atomArray, angleArray, bondArray, dihedralArray,
            atoms.size(), angles.size(), bonds.size(), dihedrals.size());
}

Molecule* readInMolecules(string filename){
    vector<Molecule> molecules;
    ifstream inFile(filename.c_str());
    string line;
    
    if(inFile.is_open()){
        //the third line starts the first molecule
        getline(inFile, line);
        getline(inFile, line);
        
        
        
    }
    else{
        Molecule *molec;
        return molec;
    }

   Molecule *moleculeArray;
   moleculeArray = (Molecule *)malloc(sizeof(Molecule) * molecules.size());
   for(int i = 0; i < molecules.size(); i++){
        moleculeArray[i] = molecules[i];
   }

   return moleculeArray;

}


