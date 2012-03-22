#include "Zmatrix_Scan.h"

//constructor
Zmatrix_Scan::Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef){
   fileName = filename;
    oplsScanner = oplsScannerRef;
    startNewMolecule = false;
}

Zmatrix_Scan::~Zmatrix_Scan(){
}

//scans in the zmatrix file and adds entries to the table
int Zmatrix_Scan::scanInZmatrix(){
    int numOfLines=0;
   ifstream zmatrixScanner(fileName.c_str());
   if( !zmatrixScanner.is_open() )
      return -1;
   else {
      string line; 
      while( zmatrixScanner.good() )
      {
        numOfLines++;
        getline(zmatrixScanner,line);

        Molecule workingMolecule;

        //check if it is a commented line,
		  //or if it is a title line
        try{
            if(line.at(0) != '#' && numOfLines > 1)
                parseLine(line,numOfLines);
        }
        catch(std::out_of_range& e){}

        if (startNewMolecule){
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

            moleculePattern.push_back(createMolecule(-1, atomArray, angleArray, bondArray, dihedralArray, 
                                 atomVector.size(), angleVector.size(), bondVector.size(), dihedralVector.size()));

            atomVector.clear();
            bondVector.clear();
            angleVector.clear();
            dihedralVector.clear();

            startNewMolecule = false;
        } 
      }

      zmatrixScanner.close();
   }
}

// adds a string line into the table
void Zmatrix_Scan::parseLine(string line, int numOfLines){

    string atomID, atomType, oplsA, oplsB, bondWith, bondDistance, angleWith, angleMeasure, dihedralWith, dihedralMeasure;
	
    stringstream ss;
     
	 //check if line contains correct format
	 int format = checkFormat(line);
	 //cout <<"-- " <<line<<endl; //DEBUG
	 //cout <<"Format: "<< format << endl; //DEBUG
	 
	 if(format == 1){
	    //read in strings in columns and store the data in temporary variables
        ss << line;    	
        ss >> atomID >> atomType >> oplsA >> oplsB >> bondWith >> bondDistance >> angleWith >> angleMeasure >> dihedralWith >> dihedralMeasure;
	     
		  //setup structures for permanent encapsulation
        Atom lineAtom;
        Bond lineBond;
        Angle lineAngle;
        Dihedral lineDihedral;

        bool hasBond = false;
        bool hasAngle = false;
        bool hasDihedral = false;

        if (oplsA.compare("-1") != 0)
        {
            lineAtom = oplsScanner->getAtom(oplsA);
            lineAtom.id = atoi(atomID.c_str());
        //    atomVector.push_back(lineAtom);
        }
        else//dummy atom
        {
            lineAtom = createAtom(atoi(atomID.c_str()), -1, -1, -1);
        //    atomVector.push_back(lineAtom); 
        }

        if (bondWith.compare("0") != 0){
            hasBond = true;

            lineBond.atom1 = lineAtom.id;
            lineBond.atom2 = atoi(bondWith.c_str());
            lineBond.distance = atof(bondDistance.c_str());
				lineBond.variable = false;
            bondVector.push_back(lineBond);
        }

        if (angleWith.compare("0") != 0){
            hasAngle = true;

            lineAngle.atom1 = lineAtom.id;
            lineAngle.atom2 = atoi(angleWith.c_str());
            lineAngle.value = atof(angleMeasure.c_str());
				lineAngle.variable = false;
            angleVector.push_back(lineAngle);
        }

        if (dihedralWith.compare("0") != 0){
            hasDihedral = true;

            lineDihedral.atom1 = lineAtom.id;
            lineDihedral.atom2 = atoi(dihedralWith.c_str());
            lineDihedral.value = atof(dihedralMeasure.c_str());
				lineDihedral.variable = false;
            dihedralVector.push_back(lineDihedral);
        }
       
        /******************************************
            BUILDING MOLECULE WITH CORRECT POSITIONS
         ******************************************/        
         // must not be a dummy atom
         if(lineAtom.z != -1 && lineAtom.y != -1 && lineAtom.x != -1){
            if(hasBond){

                // Get other atom in bond
                unsigned long otherID = getOppositeAtom(lineBond, lineAtom.id);
                Atom otherAtom = getAtom(atomVector, otherID);
                if(otherAtom.id == -1 && otherAtom.x == -1 && otherAtom.y == -1
                        && otherAtom.z == -1){
                    
                    // this should be an error but I don't know what kind.
                    cout << "Other atom not found. Error?" << endl;
                }
                
                // Move newAtom bond distance away from other atom in y direction.
                lineAtom.x = otherAtom.x;
                lineAtom.y = otherAtom.y + lineBond.distance;
                lineAtom.z = otherAtom.z;
            }
            if(hasAngle){
                // Get other atom in angle
                Atom otherAtom = createAtom(-1, -1, -1, -1);
                unsigned long otherID = getOppositeAtom(lineAngle, lineAtom.id);
                otherAtom = getAtom(atomVector, otherID);
                if(otherAtom.id == -1 && otherAtom.x == -1 && otherAtom.y == -1
                        && otherAtom.z == -1){
                    
                    // this should be an error but I don't know what kind.
                    cout << "Other atom not found. Error?" << endl;
                }
                
                // Get common atom that lineAtom and otherAtom are bonded to
                unsigned long commonID = getCommonAtom(bondVector, lineAtom.id,
                       otherID);
                Atom commonAtom = getAtom(atomVector, commonID);

                double currentAngle = getAngle(lineAtom, commonAtom, otherAtom); 
                double angleChange = lineAngle.value - currentAngle;

                lineAtom = rotateAtom(lineAtom, commonAtom, otherAtom, angleChange);
            }
            if(hasDihedral){
                //get other atom in the dihedral
                unsigned long otherID = getOppositeAtom(lineDihedral, lineAtom.id);
                Atom otherAtom = getAtom(atomVector, otherID);

                /*********
                  There are guranteed to be 4 atoms involved in the dihedral
                  because it takes atleast 4 atoms to define two non equal
                  planes.
                **********/
                
                //get all of the atoms bonded to lineAtom
                vector<unsigned long> bondedToLineAtom = getAllBonds(bondVector, lineAtom.id);
                //get all of the atoms bonded to  otherAtom
                vector<unsigned long> bondedToOtherAtom = getAllBonds(bondVector, otherAtom.id);
                
                /*find the intersection of the two previous vectors.
                because of how the zMatrix file is set up the intersection
                should be of size 1*/
                vector<unsigned long> intersection = getIntersection(bondedToLineAtom, bondedToOtherAtom);
                                 
                //find bond that bonds together two of the atoms in the intersection
                Bond linkingBond = createBond(-1, -1, -1, false);
                for(int i = 0; i < intersection.size() - 1; i++){
                    for(int j = i + 1; i < intersection.size(); i++){
                        if(getOppositeAtom(bondVector, intersection[i]) == intersection[j]){
                            //linkingBond = getBond(bondVector, intersection[i], intersection[j]);
                        }
                    }
                } 
         
                /**
                plane 1 is lineAtom and atoms in linking bond
                plane 2 is otherAtom and atoms in linking bond
                the bond creates the vector about which lineAtom will be rotated.
                */    
                
                //find the angle between the planes.
                //find the angle needed to rotate.
                //rotate lineAtom needed degrees about linkbond.


        }

         atomVector.push_back(lineAtom);

    }
    else if(format == 2)
        startNewMolecule = true;
	 else if(format == 3){
	     startNewMolecule = true;
	 }
	 if (previousFormat >= 3 && format == -1)
	     handleZAdditions(line, previousFormat);
	 		  
	     
	 previousFormat = format;

    }
}
// check if line contains the right format...
int Zmatrix_Scan::checkFormat(string line){
    int format =-1; 
    stringstream iss(line);
	 stringstream iss2(line);
	 string atomType, someLine;
	 int atomID, oplsA, oplsB, bondWith, angleWith,dihedralWith,extra;
	 double bondDistance, angleMeasure, dihedralMeasure;	 
	 
	 // check if it is the normal 11 line format
	  if( iss >> atomID >> atomType >> 
	      oplsA >> oplsB >> 
			bondWith >> bondDistance >> 
			angleWith >> angleMeasure >> 
			dihedralWith >> dihedralMeasure >> extra)
			format = 1;
	 else{
	     someLine = line;
	     if(someLine.find("TERZ")!=string::npos)
		      format = 2;
		  else if(someLine.find("Geometry Variations follow")!=string::npos)
		      format = 3;
		  else if(someLine.find("Variable Bonds follow")!=string::npos)
		      format = 4;
		  else if(someLine.find("Additional Bonds follow")!=string::npos)
		      format = 5;
		  else if(someLine.find("Harmonic Constraints follow")!=string::npos)
		      format = 6;
		  else if(someLine.find("Variable Bond Angles follow")!=string::npos)
		      format = 7;
		  else if(someLine.find("Additional Bond Angles follow")!=string::npos)
		      format = 8;
		  else if(someLine.find("Variable Dihedrals follow")!=string::npos)
		      format = 9;
		  else if(someLine.find("Additional Dihedrals follow")!=string::npos)
		      format = 10;
		 else if(someLine.find("Domain Definitions follow")!=string::npos)
		      format = 11;
	 	 else if(someLine.find("Final blank line")!=string::npos)
		      format = -2;  
	 }	 
	 return format;
}

void Zmatrix_Scan::handleZAdditions(string line, int cmdFormat){
    vector<int> atomIds;
	 int id;
	 stringstream tss(line.substr(0,15) );
	 if(line.find("AUTO")!=string::npos){
	 }
	 else{
        while(tss >> id){
            atomIds.push_back(id);
            if(tss.peek()=='-'||tss.peek()==','||tss.peek()==' ')
                tss.ignore();
        }
		  int start, end=0;
		  if( atomIds.size()== 1){
		     start = atomIds[0];
			  end = atomIds[0]; 
			  }
		  else if(atomIds.size() == 2){
		     start = atomIds[0]; end = atomIds[1];
			 }
		  //cout << "start: " << start<< " end: " <<end <<endl; //DEBUG
            switch(cmdFormat){
	             case 3:
    				// Geometry Variations follow 
		              break;
    		      case 4:
    				// Variable Bonds follow
					    for(int i=0; i< moleculePattern[0].numOfBonds; i++){
					        if(  moleculePattern[0].bonds[i].atom1 >= start &&  moleculePattern[0].bonds[i].atom1 <= end){
							       //cout << "Bond Atom1: "<<  moleculePattern[0].bonds[i].atom1 << " : " <<  moleculePattern[0].bonds[i].variable<<endl;//DEBUG
			                   moleculePattern[0].bonds[i].variable = true;
									 }
						 }
    		          break;
    		      case 5:
    				//  Additional Bonds follow 
    		          break;
    		      case 6:
    				// Harmonic Constraints follow 
    		          break;
    		      case 7:
    				//  Variable Bond Angles follow
					    for(int i=0; i<  moleculePattern[0].numOfAngles; i++){
					        if(  moleculePattern[0].angles[i].atom1 >= start && moleculePattern[0].angles[i].atom1 <= end){
							      //cout << "Angle Atom1: "<<  moleculePattern[0].angles[i].atom1 << " : " << moleculePattern[0].angles[i].variable << endl;//DEBUG
					            moleculePattern[0].angles[i].variable = true;
									}
						 }
    		          break;
    		      case 8:
    				// Additional Bond Angles follow
    		          break;
    		      case 9:
    				// Variable Dihedrals follow
					    for(int i=0; i< moleculePattern[0].numOfDihedrals; i++){
					        if(  moleculePattern[0].dihedrals[i].atom1 >= start &&  moleculePattern[0].dihedrals[i].atom1 <= end ) {
							      //cout << "Dihedral Atom1: "<<  moleculePattern[0].dihedrals[i].atom1 << " : " <<   moleculePattern[0].dihedrals[i].variable << endl;//DEBUG
					            moleculePattern[0].dihedrals[i].variable = true;
									}
						 }
    		          break;
    		      case 10:
    				//  Domain Definitions follow
    		          break;
					default:
					//Do nothing
					    break;
		  }
	 }
}

//returns a vector containing all the node to node hops that 
//are more than 3.
vector<Hop> Zmatrix_Scan::calculateHops(Molecule molec){
    vector<Hop> newHops;
    int **graph;
	 int size = molec.numOfAtoms;
    
	 
	 //cout << "ZMATRIX-- Creating Graph "<<endl;
	 buildAdjacencyMatrix(graph,molec);
	 //cout << "ZMATRIX-- Creating Graph Sucessful " <<endl;
	 
    /* cout << "  ";
    for(int r=0;r<size; r++)
       cout<< r+1 << " ";
    cout << endl;
      
   	
    for(int r=0;r<size; r++){
       cout << r+1 << " ";
       for(int c=0;c<size; c++)
          cout << graph[r][c]<<" ";
       cout << endl;
    }	// DEBUG */
	 
	 for(int atom1=0; atom1<size; atom1++){
	     for(int atom2=atom1+1; atom2<size; atom2++){
		     //cout << "ZMATRIX-- Finding Hops "<< atom1+1<<" - "<<atom2+1<<endl;
		     int distance = findHopDistance(atom1,atom2,size,graph);
			  //cout << "ZMATRIX-- Hop Distance: "<< distance << endl;
			  if(distance >3){
			      Hop tempHop = createHop(atom1+1,atom2+1,distance); //+1 because atoms start at 1
				   newHops.push_back(tempHop);					
			      //cout << "ZMATRIX-- Creating and Adding new Hop atom1: "<< atom1+1<<" atom2: "<< atom2+1<<" \n--Distance: " << distance<<endl;
			  }  		      
		  }
	 }
	 return newHops; 
}

//checks to see if int item is pressent in the vector
//liear search
bool Zmatrix_Scan::contains(vector<int> &vect, int item){
     for(int i=0; i<vect.size(); i++){
        if(vect[i]==item)
           return true;
     }
     return false;
}


//returns the node hop distance between two atoms
int Zmatrix_Scan::findHopDistance(int atom1,int atom2,int size, int **graph){
    map<int,int> distance;
    queue<int> Queue;
    vector<int> checked;
    vector<int> bonds;
   
      
    Queue.push(atom1);
    checked.push_back(atom1);
    distance.insert( pair<int,int>(atom1,0) );	
   	
    while(!Queue.empty()){
       int target = Queue.front();
       Queue.pop();
       if(target == atom2)
          return distance[target];
       //if(distance[target]>4)
       //   return -1;
      	
    	//get/push all bonds that are conected to target
       bonds.clear();
       for(int col=0;col<size;col++){
          if( graph[target][col]==1 )
             bonds.push_back(col);
       }
         
       for(int x=0; x<bonds.size();x++){
          int currentBond = bonds[x];
          if(!contains(checked,currentBond) ){
             checked.push_back(currentBond);
             int newDistance = distance[target]+1;
             distance.insert(pair<int,int>(currentBond, newDistance));
             Queue.push(currentBond);
          }
       }
    }
}

//Creates a graph or Adjacency matrix so the createHops function knows
//which nodes/atoms are neighbors/linked
void Zmatrix_Scan::buildAdjacencyMatrix(int **&graph, Molecule molec){
    int size = molec.numOfAtoms;	
	 graph =  new int*[size]; //create colums
	 for(int i=0; i<size; i++) //create rows
	      graph[i]=new int[size];	
	 
	 //fill with zero
    for(int c=0; c<size; c++)
        for(int r=0; r<size; r++)
            graph[c][r]=0;
	 
	 //cout << "ZMATRIX-- Number of Bonds: "<< molec.numOfBonds <<endl;
	 //fill with adjacent array with bonds
	 for(int x=0; x<molec.numOfBonds; x++){
	     Bond bond = molec.bonds[x];
		  //cout << "ZMATRIX-- Bonds: "<< x<< " atom1: "<< molec.bonds[x].atom1<< " atom2: "<< molec.bonds[x].atom2 <<endl;
	     graph[bond.atom1-1][bond.atom2-1]=1;
		  graph[bond.atom2-1][bond.atom1-1]=1;
	 }
}

//return a vector of all the molecules scanned in.
//resed the id's of all molecules and atoms to be the location in Atom array
//based on the startingID. Atoms adn Molecules should be stored one behind the other.
vector<Molecule> Zmatrix_Scan::buildMolecule(int startingID){
    vector<Molecule> newMolecules;
	 //need a deep copy of molecule pattern incase it is modified.
	  for (int i = 0; i < moleculePattern.size(); i++){
	      Atom *atomCopy = new Atom[ moleculePattern[i].numOfAtoms] ;
			for(int a=0; a <  moleculePattern[i].numOfAtoms ; a++)
			    atomCopy[a]=  moleculePattern[i].atoms[a];
	      
			Bond *bondCopy = new Bond[ moleculePattern[i].numOfBonds] ;
			for(int a=0; a <  moleculePattern[i].numOfBonds ; a++)
			    bondCopy[a]=  moleculePattern[i].bonds[a];

			Angle *angleCopy = new Angle[ moleculePattern[i].numOfAngles] ;
			for(int a=0; a <  moleculePattern[i].numOfAngles ; a++)
			    angleCopy[a]=  moleculePattern[i].angles[a];

	      Dihedral *dihedCopy = new Dihedral[ moleculePattern[i].numOfDihedrals];
			for(int a=0; a <  moleculePattern[i].numOfDihedrals ; a++)
			    dihedCopy[a]=  moleculePattern[i].dihedrals[a];
			
			//calculate and add array of Hops to the molecule
			vector<Hop> calculatedHops;
			//cout << "ZMATRIX-- calculating hops " << startingID <<endl; 
			calculatedHops = calculateHops(moleculePattern[i]);
			int numOfHops = calculatedHops.size();
			Hop *hopCopy = new Hop[numOfHops];
			for(int a=0; a < numOfHops; a++)
			    hopCopy[a] = calculatedHops[a];

			
			Molecule molecCopy = createMolecule(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
			                                     moleculePattern[i].numOfAtoms, 
															 moleculePattern[i].numOfAngles,
															 moleculePattern[i].numOfBonds,
															 moleculePattern[i].numOfDihedrals,
															 numOfHops);
															 
		   newMolecules.push_back(molecCopy);	      
	  }
	 
	 // cout << "ZMATRIX-- startingID: " << startingID <<endl;    
	 //cout << "ZMATRIX-- newMolecules size: " <<  moleculePattern.size() <<endl;
	
	 //cout << "ZMATRIX-- newMolecules[i] ID: " <<  newMolecules[0].id <<endl;
	 for (int i = 0; i < newMolecules.size(); i++)
    {
	     if(i == 0){
		      newMolecules[i].id = startingID;
		  }
		  else
		      newMolecules[i].id = newMolecules[i-1].id + newMolecules[i-1].numOfAtoms; 
    }
	 
    for (int j = 0; j < newMolecules.size(); j++)
    {
        Molecule newMolecule = newMolecules[j];
        //map unique IDs to atoms within structs based on startingID
        for(int i = 0; i < newMolecules[j].numOfAtoms; i++){
		      //cout << "ZMATRIX-- newMolecules[j].atoms[i] ID: " <<  moleculePattern[j].atoms[i].id <<endl;
            int atomID = newMolecule.atoms[i].id - 1;
				//cout << "ZMATRIX-- atomID: " <<  atomID <<endl;
				//cout << "ZMATRIX-- atomID memory Loc: " << &newMolecules[j].atoms[i]  <<endl;
            newMolecule.atoms[i].id = atomID + newMolecule.id;
				//cout << "ZMATRIX-- newMolecule.atoms[i] ID + Molec ID: " <<  newMolecule.atoms[i].id <<endl;
        }
        for (int i = 0; i < newMolecule.numOfBonds; i++){
            int atom1ID = newMolecule.bonds[i].atom1 - 1;
            int atom2ID = newMolecule.bonds[i].atom2 - 1;
            
            newMolecule.bonds[i].atom1 = atom1ID + newMolecule.id;
            newMolecule.bonds[i].atom2 = atom2ID + newMolecule.id;
        }
        for (int i = 0; i < newMolecule.numOfAngles; i++){
            int atom1ID = newMolecule.angles[i].atom1 - 1;
            int atom2ID = newMolecule.angles[i].atom2 - 1;
            
            newMolecule.angles[i].atom1 = atom1ID + newMolecule.id;
            newMolecule.angles[i].atom2 = atom2ID + newMolecule.id;
        }
        for (int i = 0; i < newMolecule.numOfDihedrals; i++){
            int atom1ID = newMolecule.dihedrals[i].atom1 - 1;
            int atom2ID = newMolecule.dihedrals[i].atom2 - 1;
            
            newMolecule.dihedrals[i].atom1 = atom1ID + newMolecule.id;
            newMolecule.dihedrals[i].atom2 = atom2ID + newMolecule.id;
        }
		   for (int i = 0; i < newMolecule.numOfHops; i++){
            int atom1ID = newMolecule.hops[i].atom1 - 1;
            int atom2ID = newMolecule.hops[i].atom2 - 1;
            
            newMolecule.hops[i].atom1 = atom1ID + newMolecule.id;
            newMolecule.hops[i].atom2 = atom2ID + newMolecule.id;
        }

    }

    return newMolecules;
}
