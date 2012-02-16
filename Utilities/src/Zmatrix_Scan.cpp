#include "Zmatrix_Scan.h"

//constructor
Zmatrix_Scan::Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef){
   fileName = filename;
    oplsScanner = oplsScannerRef;
}

Zmatrix_Scan::~Zmatrix_Scan(){
}

//scans in the zmatrix file and adds entries to the table
int Zmatrix_Scan::scanInZmatrix(string filename){
    int numOfLines=0;
   ifstream zmatrixScanner(filename.c_str());
   if( !zmatrixScanner.is_open() )
      return -1;
   else {
      string line; 
      while( zmatrixScanner.good() )
      {
        numOfLines++;
        getline(zmatrixScanner,line);

        //check if it is a commented line,
		  //or if it is a title line
        if(line.at(0) != '#' && numOfLines > 1)
            parseLine(line,numOfLines);
      }
        Atom* atomArray;
        Bond* bondArray;
        Angle* angleArray;
        Dihedral* dihedralArray;

        atomArray = (Atom*) malloc(sizeof(Atom) * atomVector.size());
        bondArray = (Bond*) malloc(sizeof(Bond) * bondVector.size());
        angleArray = (Angle*) malloc(sizeof(Angle) * angleVector.size());
        dihedralArray = (Dihedral*) malloc(sizeof(Dihedral) * dihedralVector.size());

        for (int i = 0; i < atomVector.size(); i++){
            atomArray[i] = atomVector[i];
        }
        for (int i = 0; i < bondVector.size(); i++){
            bondArray[i] = bondVector[i];
        }
        for (int i = 0; i < angleVector.size(); i++){
            angleArray[i] = angleVector[i];
        }
        for (int i = 0; i < dihedralVector.size(); i++){
            dihedralArray[i] = dihedralVector[i];
        }

        moleculePattern = createMolecule(-1, atomArray, angleArray, bondArray, dihedralArray, 
                            atomVector.size(), angleVector.size(), bondVector.size(), dihedralVector.size());

      zmatrixScanner.close();
   }
}

// adds a string line into the table
void Zmatrix_Scan::parseLine(string line, int numOfLines){

    string atomID, atomType, oplsA, oplsB, bondWith, bondDistance, angleWith, angleMeasure, dihedralWith, dihedralMeasure;
	 
	 //check if line contains correct format
	 if(checkFormat(line)){
	     stringstream ss ; 
        ss << line;    	
        ss >> atomID >> atomType >> oplsA >> oplsB >> bondWith >> bondDistance >> angleWith >> angleMeasure >> dihedralWith >> dihedralMeasure;
	 
        Atom lineAtom;
        Bond lineBond;
        Angle lineAngle;
        Dihedral lineDihedral;

        if (oplsA.compare("-1") != 0)
        {
	     lineAtom = oplsScanner->getAtom(oplsA);
         lineAtom.id = atoi(atomID.c_str());
         atomVector.push_back(lineAtom);

         idMap[lineAtom.id] = idMap.size();

         if (bondWith.compare("0") != 0){
            lineBond.atom1 = lineAtom.id;
            lineBond.atom2 = atoi(bondWith.c_str());
            lineBond.distance = atof(bondDistance.c_str());
            bondVector.push_back(lineBond);
         }

         if (angleWith.compare("0") != 0){
            lineAngle.atom1 = lineAtom.id;
            lineAngle.atom2 = atoi(angleWith.c_str());
            lineAngle.value = atof(angleMeasure.c_str());
            angleVector.push_back(lineAngle);
         }

         if (dihedralWith.compare("0") != 0){
            lineDihedral.atom1 = lineAtom.id;
            lineDihedral.atom2 = atoi(dihedralWith.c_str());
            lineDihedral.value = atof(dihedralMeasure.c_str());
            dihedralVector.push_back(lineDihedral);
         }
        }
    }
}

// check if line contains the right format...
bool Zmatrix_Scan::checkFormat(string line){ 
    stringstream iss(line);
	 int count =0;
	 string word;
	 while (iss >> word)
	     count ++;
	 return (count >= 11);
}

//return molecule
Molecule Zmatrix_Scan::buildMolecule(int startingID){
    Molecule newMolecule = moleculePattern;
    newMolecule.numOfAtoms = idMap.size();

    for(int i = 0; i < atomVector.size(); i++){
        int atomID = newMolecule.atoms[i].id;
        newMolecule.atoms[i].id = idMap[atomID] + startingID;
    }
    for (int i = 0; i < bondVector.size(); i++){
        int atom1ID = newMolecule.bonds[i].atom1;
        int atom2ID = newMolecule.bonds[i].atom2;
        
        newMolecule.bonds[i].atom1 = atom1ID + startingID;
        newMolecule.bonds[i].atom2 = atom2ID + startingID;
    }
    for (int i = 0; i < bondVector.size(); i++){
        int atom1ID = newMolecule.angles[i].atom1;
        int atom2ID = newMolecule.angles[i].atom2;
        
        newMolecule.angles[i].atom1 = atom1ID + startingID;
        newMolecule.angles[i].atom2 = atom2ID + startingID;
    }
    for (int i = 0; i < bondVector.size(); i++){
        int atom1ID = newMolecule.dihedrals[i].atom1;
        int atom2ID = newMolecule.dihedrals[i].atom2;
        
        newMolecule.dihedrals[i].atom1 = atom1ID + startingID;
        newMolecule.dihedrals[i].atom2 = atom2ID + startingID;
    }

    return newMolecule;
}
