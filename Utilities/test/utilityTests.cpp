#include "../src/Opls_Scan.h"
#include "../src/Zmatrix_Scan.h"
#include "../src/State_Scan.h"
#include "stateTest.h"
#include "configurationTest.h"
#include "geometricTest.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;
string oplsPath = "oplsua.par";

string atomNumber1 = "418";
string atomNumber2 = "1";
string atomNumber3 = "53";

double atom1Charge =  0.365;
double atom1Sigma = 3.850;
double atom1Epsilon = 0.080;

double atom2Charge = 0.5;
double atom2Sigma = 3.75;
double atom2Epsilon = 0.105;
double atom2V[4] = {};	

double atom3Charge = .640;
double atom3Sigma = 2.250;
double atom3Epsilon = 0.05;

void testGetAtom(Opls_Scan scan){

    cout << "Testing Opls_Scan.getAtom\n" << endl;
    Atom atom1 = scan.getAtom(atomNumber1);
    Atom atom2 = scan.getAtom(atomNumber2);
    Atom atom3 = scan.getAtom(atomNumber3);
    
    //printf("%f, %f, %f\n", atom1.charge, atom1.sigma, atom1.epsilon);
    assert(atom1.charge == atom1Charge);
    assert(atom1.sigma == atom1Sigma);
    assert(atom1.epsilon == atom1Epsilon);
    
    assert(atom2.charge == atom2Charge);
    assert(atom2.sigma == atom2Sigma);
    assert(atom2.epsilon == atom2Epsilon);
    
    assert(atom3.charge == atom3Charge);
    assert(atom3.sigma == atom3Sigma);
    assert(atom3.epsilon == atom3Epsilon);
    cout << "Testing Opls_Scan.getAtom Completed\n" << endl;
}

void testGetSigma(Opls_Scan scan){
    cout << "Testing Opls_Scan.getSigma\n" << endl;
    assert(scan.getSigma(atomNumber1) == atom1Sigma); 
    assert(scan.getSigma(atomNumber2) == atom2Sigma); 
    assert(scan.getSigma(atomNumber3) == atom3Sigma); 
    cout << "Testing Opls_Scan.getSigma Completed\n" << endl;
}

void testGetEpsilon(Opls_Scan scan){
    cout << "Testing Opls_Scan.getEpsilon\n" << endl;
    assert(scan.getEpsilon(atomNumber1) == atom1Epsilon);
    assert(scan.getEpsilon(atomNumber2) == atom2Epsilon);
    assert(scan.getEpsilon(atomNumber3) == atom3Epsilon);
    cout << "Testing Opls_Scan.getEpsilon Completed\n" << endl;
}

void testGetCharge(Opls_Scan scan){
    cout << "Testing Opls_Scan.getCharge\n" << endl;
    assert(scan.getCharge(atomNumber1) == atom1Charge);
    assert(scan.getCharge(atomNumber2) == atom2Charge);
    assert(scan.getCharge(atomNumber3) == atom3Charge);
    cout << "Testing Opls_Scan.getCharge Completed\n" << endl;
}

string atomNumber4 ="415";
double atom4Fourier[4] = {0.0 , -2.500, 1.250, 3.100 };

string atomNumber5 ="002";
double atom5Fourier[4] = { 1.363, 0.343, -0.436, -1.121};
	
string atomNumber6 ="073";
double atom6Fourier[4] = { 0.0, -4.0396, 1.2261, 3.5637 };

void testGetFourier(Opls_Scan scan){
    cout << "Testing Opls_Scan.getVvalues\n" << endl;
	 Fourier f = scan.getFourier(atomNumber4);
	 for(int i=0; i<4; i++){
	     assert(f.vValues[i]==atom4Fourier[i]);	 
	 }
	 
	 f = scan.getFourier(atomNumber5);
	 for(int i=0; i<4; i++){
	     assert(f.vValues[i]==atom5Fourier[i]);	 
	 }
	 
	 f = scan.getFourier(atomNumber6);
	 for(int i=0; i<4; i++){
	     assert(f.vValues[i]==atom6Fourier[i]);	 
	 }
    cout << "Testing Opls_Scan.getVvalues Completed\n" << endl;
}

void testPDBoutput(){
    Atom* pdbAtoms;
    Environment pdbEnviro;

    pdbAtoms = (Atom *) malloc(sizeof(Atom)*2);
    pdbEnviro.numOfAtoms = 2;

    for (int i = 0; i < 2; i++){
        pdbAtoms[i].x = 1.1 * (i + 1);
        pdbAtoms[i].y = 2.2 * (i + 1);
        pdbAtoms[i].z = 3.3 * (i + 1);

        pdbAtoms[i].id = i;
    }

    string fileName = "test.pdb";

    writePDB(pdbAtoms, pdbEnviro, fileName);
    
    string lines[2];
    ifstream readPDB;
    readPDB.open(fileName.c_str());

    for (int i = 0; i < 2; i++){
        stringstream ss;
        string tokens[9];

        getline(readPDB, lines[i]);
        ss << lines[i];
        ss >> tokens[0] >> tokens[1] >> tokens[2] >> tokens[3] >> tokens[4] >> tokens[5] >> tokens[6] >> tokens[7] >> tokens[8];

        assert(tokens[0].compare("ATOM") == 0);
        assert(atoi(tokens[1].c_str()) == pdbAtoms[i].id);
        assert(atof(tokens[6].c_str()) == pdbAtoms[i].x);
        assert(atof(tokens[7].c_str()) == pdbAtoms[i].y);
        assert(atof(tokens[8].c_str()) == pdbAtoms[i].z);
    }

    readPDB.close();

    cout << "testPDBoutput successful." << endl;
}

bool asserTwoBool(bool b1, bool b2){
    if(b1 && b2)
        return true;
    else if(!b1 && !b2)
        return true;
    else
        return false;
}

bool percentDifference(double d1, double d2){
    double difference = d2-d1;
    double average = (d2+d1)/d2;
    double percentDiff = (difference/average)*100;
    //cout <<"Percent Diff: " << percentDiff << endl;;
    return percentDiff < 3;
}

//Test Z-matrix files using mesh.z
Molecule createMeshZMolecules(Opls_Scan scanner){

    //1 S    200    0    0    0.000000   0    0.000000   0    0.000000        0  
    Atom atom1=scanner.getAtom("200");
    atom1.id=1;

    // 2 DUM   -1    0    1    0.500000   0    0.000000   0    0.000000        0
    Atom atom2=createAtom(2,-1,-1,-1,-1,-1,-1);
	 Bond bond2=createBond(2,1,0.5,false);
	     
    //3 DUM   -1    0    2    0.500000   1   90.000000   0    0.000000        0 
    Atom atom3=createAtom(3,-1,-1,-1,-1,-1,-1);
	 Bond bond3=createBond(3,2,0.5,false);
	 Angle angle3=createAngle(3,1,90,false);
	 
	 //4 hH    204    0    1    1.336532   2   90.000000   3  180.000000        0    
    Atom atom4=scanner.getAtom("204");
        atom4.id=4;
	 Bond bond4=createBond(4,1,1.336532,true);
	 Angle angle4=createAngle(4,2,90,false);
	 Dihedral dihed4=createDihedral(4,3,180,false);
	 
    //5 C    217    0    1    1.811119   4   96.401770   2  180.000000        0
    Atom atom5=scanner.getAtom("217");
        atom5.id=5;
	 Bond bond5=createBond(5,1,1.8119,true);
	 Angle angle5=createAngle(5,4,96.401770,true);
	 Dihedral dihed5=createDihedral(5,2,180,false);
	 
    //6 HC   140    0    5    1.090187   1  110.255589   4  179.999947        0
	 Atom atom6=scanner.getAtom("140");
        atom6.id=6;
	 Bond bond6=createBond(6,5,1.090187,true);
	 Angle angle6=createAngle(6,1,110.255589,true);
	 Dihedral dihed6=createDihedral(6,4,179.999947,true);    

    //7 HC   140    0    5    1.090135   6  108.527646   1  121.053891        0
	 Atom atom7=scanner.getAtom("140");
        atom7.id=7;
	 Bond bond7=createBond(7,5,1.090135,true);
	 Angle angle7=createAngle(7,6,108.527646,true);
	 Dihedral dihed7=createDihedral(7,1,121.053891,true);    

    //8 HC   140    0    5    1.090135   6  108.527646   1  238.946114        0
	 Atom atom8=scanner.getAtom("140");
        atom8.id=8;
	 Bond bond8=createBond(8,5,1.090135,true);
	 Angle angle8=createAngle(8,6,108.527646,true);
	 Dihedral dihed8=createDihedral(8,1,238.946114,true);
	 
	 /* HOPS and BONDS Diagram of Mesh.z
	 //      1--2--3
	 //      |\
    //      | \
    //      4  5
	 //        /|\
    //       / | \
    //      6  7  8
	 */
	 //All hops that have a hop distance > 4
	 Hop hop1= createHop(3,6,4);
	 Hop hop2= createHop(3,7,4);
	 Hop hop3= createHop(3,8,4);
	 Hop hop4= createHop(3,5,3);
	 Hop hop5= createHop(3,4,3);
	 Hop hop6= createHop(2,6,3);
	 Hop hop7= createHop(2,7,3);
	 Hop hop8= createHop(2,8,3);
	 Hop hop9= createHop(4,6,3);
	 Hop hop10= createHop(4,7,3);
	 Hop hop11= createHop(4,8,3);
	 
	 Atom *atomPtr = new Atom[8];
    Bond *bondPtr = new Bond[7];
    Angle *anglePtr = new Angle[6];
    Dihedral *dihedPtr = new Dihedral[5];
	 Hop *hopPtr = new Hop[11];
	 

    atomPtr[0]=atom1;  atomPtr[1]=atom2; atomPtr[2]=atom3; atomPtr[3]=atom4; atomPtr[4]=atom5; atomPtr[5]=atom6; atomPtr[6]=atom7; atomPtr[7]=atom8;
    bondPtr[0]=bond2; bondPtr[1]=bond3; bondPtr[2]=bond4; bondPtr[3]=bond5; bondPtr[4]=bond6; bondPtr[5]=bond7; bondPtr[6]=bond8;
    anglePtr[0]=angle3; anglePtr[1]=angle4; anglePtr[2]=angle5; anglePtr[3]=angle6; anglePtr[4]=angle7; anglePtr[5]=angle8;
    dihedPtr[0]=dihed4; dihedPtr[1]=dihed5; dihedPtr[2]=dihed6; dihedPtr[3]=dihed7; dihedPtr[4]=dihed8;
	 hopPtr[0]=hop1; hopPtr[1]=hop2; hopPtr[2]=hop3; hopPtr[3]=hop4; hopPtr[4]=hop5; hopPtr[5]=hop6;  hopPtr[6]=hop7; hopPtr[7]=hop8; hopPtr[8]=hop9; hopPtr[9]=hop10; hopPtr[10]=hop11;
 
 
	 return createMolecule(1,atomPtr,anglePtr,bondPtr,dihedPtr,hopPtr,8,6,7,5,11);
}

//checks to see if two molecule objects are equal
//uses assert
void compareTestMolecules(Molecule molec1, Molecule molec2){
    // check if id's are equal
    //cout << "-- Moleculde1 id: " << molec1.id << " Moleculde2 id: " << molec2.id << endl; 
    assert(molec1.id == molec2.id);	 

	 //check atoms
        //cout << endl<< "Testing Atoms Array" << endl;
        //cout << "-- Moleculde1 # atoms: " << molec1.numOfAtoms << " Moleculde2 # atoms: " << molec2.numOfAtoms << endl;
	 assert(molec1.numOfAtoms == molec2.numOfAtoms);
	 for(int i=0; i< molec1.numOfAtoms; i++){
            //cout << "-- Atom1#: " << i <<endl;
            //cout << "-- Atom1 Id: " << molec1.atoms[i].id << " Atom2 Id: " << molec2.atoms[i].id <<endl;
            //cout << "-- Atom1 charge: " << molec1.atoms[i].charge << " Atom2 charge: " << molec2.atoms[i].charge <<endl;
            //cout << "-- Atom1 sigma: " << molec1.atoms[i].sigma << " Atom2 sigma: " << molec1.atoms[i].sigma <<endl;
            //cout << "-- Atom1 epsilon: " << molec1.atoms[i].epsilon << " Atom2 epsilon: " << molec1.atoms[i].epsilon<<endl;
	     
            assert(molec1.atoms[i].id==molec2.atoms[i].id);
	     assert(percentDifference(molec1.atoms[i].charge,molec2.atoms[i].charge));
	     assert(percentDifference(molec1.atoms[i].sigma,molec2.atoms[i].sigma));
	     assert(percentDifference(molec1.atoms[i].epsilon,molec2.atoms[i].epsilon));		
	 }
	 //check bond
        //cout << endl << "Testing Bonds Array" << endl;
        //cout << "-- Moleculde1 # bonds: " << molec1.numOfBonds << " Moleculde2 # bonds: " << molec2.numOfBonds << endl;
	 assert(molec1.numOfBonds == molec2.numOfBonds);
	 for(int i=0; i< molec1.numOfBonds; i++){
            //cout << "-- Atom1#: " << i <<endl;
            //cout << "-- Atom1 atom1: " << molec1.bonds[i].atom1 << " Atom2 atom1: " << molec2.bonds[i].atom1 <<endl;
            //cout << "-- Atom1 atom2: " << molec1.bonds[i].atom2 << " Atom2 atom2: " << molec2.bonds[i].atom2 <<endl;
            //cout << "-- Atom1 distance: " << molec1.bonds[i].distance << " Atom2 distance: " << molec2.bonds[i].distance <<endl;
            //cout << "-- Atom1 variable: " << molec1.bonds[i].variable << " Atom2 variable: " << molec2.bonds[i].variable <<endl;
                         
	     assert(molec1.bonds[i].atom1 == molec2.bonds[i].atom1);
	     assert(molec1.bonds[i].atom2 == molec2.bonds[i].atom2);
	     assert(percentDifference(molec1.bonds[i].distance, molec2.bonds[i].distance));
	     assert(asserTwoBool( molec1.bonds[i].variable,molec2.bonds[i].variable));
	 }
	 //check angles
        //cout << endl << "Testing Angles Array" << endl;
        //cout << "-- Moleculde1 # angles: " << molec1.numOfBonds << " Moleculde2 # angles: " << molec2.numOfBonds << endl;
	 assert(molec1.numOfAngles == molec2.numOfAngles);
	 for(int i=0; i< molec1.numOfAngles; i++){
            //cout << "-- Atom1#: " << i <<endl;
            //cout << "-- Atom1 atom1: " << molec1.angles[i].atom1 << " Atom2 atom1: " << molec2.angles[i].atom1 <<endl;
            //cout << "-- Atom1 atom2: " << molec1.angles[i].atom2 << " Atom2 atom2: " << molec2.angles[i].atom2 <<endl;
            //cout << "-- Atom1 distance: " << molec1.angles[i].value << " Atom2 distance: " << molec2.angles[i].value <<endl;
            //cout << "-- Atom1 variable: " << molec1.angles[i].variable << " Atom2 variable: " << molec2.angles[i].variable <<endl;
	     
            assert(molec1.angles[i].atom1 == molec2.angles[i].atom1);
	     assert(molec1.angles[i].atom2 == molec2.angles[i].atom2);
	     assert(percentDifference(molec1.angles[i].value,molec2.angles[i].value));
	     assert(asserTwoBool(molec1.angles[i].variable, molec2.angles[i].variable));
	 }
	  //check dihederals
        //cout << endl << "Testing Dihederals Array" << endl;
        //cout << "-- Moleculde1 # dihedrals: " << molec1.numOfBonds << " Moleculde2 # dihedrals: " << molec2.numOfBonds << endl;
	 assert(molec1.numOfDihedrals == molec2.numOfDihedrals);
	 for(int i=0; i< molec1.numOfDihedrals; i++){
            //cout << "-- Atom1#: " << i <<endl;
            //cout << "-- Atom1 atom1: " << molec1.dihedrals[i].atom1 << " Atom2 atom1: " << molec2.dihedrals[i].atom1 <<endl;
            //cout << "-- Atom1 atom2: " << molec1.dihedrals[i].atom2 << " Atom2 atom2: " << molec2.dihedrals[i].atom2 <<endl;
            //cout << "-- Atom1 distance: " << molec1.dihedrals[i].value << " Atom2 distance: " << molec2.dihedrals[i].value <<endl;
            //cout << "-- Atom1 variable: " << molec1.dihedrals[i].variable << " Atom2 variable: " << molec2.dihedrals[i].variable <<endl;

	     assert(molec1.dihedrals[i].atom1 == molec2.dihedrals[i].atom1);
	     assert(molec1.dihedrals[i].atom2 == molec2.dihedrals[i].atom2);
	     assert(percentDifference(molec1.dihedrals[i].value,molec2.dihedrals[i].value));
	     assert(asserTwoBool(molec1.dihedrals[i].variable,molec2.dihedrals[i].variable));
	 }
	 //check hops
        //cout << endl << "Testing Hops Array" << endl;
        //cout << "-- Moleculde1 # hops: " << molec1.numOfHops << " Moleculde2 # hops: " << molec2.numOfHops << endl;
	 assert(molec1.numOfHops == molec2.numOfHops);
	 for(int i=0; i< molec1.numOfHops; i++){
            //cout << "-- Atom1#: " << i <<endl;
            //cout << "-- Atom1 atom1: " << molec1.hops[i].atom1 << " Atom2 atom1: " << molec2.hops[i].atom1 <<endl;
            //cout << "-- Atom1 atom2: " << molec1.hops[i].atom2 << " Atom2 atom2: " << molec2.hops[i].atom2 <<endl;
            //cout << "-- Atom1 distance: " << molec1.hops[i].value << " Atom2 distance: " << molec2.hops[i].value <<endl;
            //cout << "-- Atom1 variable: " << molec1.hops[i].variable << " Atom2 variable: " << molec2.hops[i].variable <<endl;

	     assert(molec1.hops[i].atom1 == molec2.hops[i].atom1);
	     assert(molec1.hops[i].atom2 == molec2.hops[i].atom2);
	     assert(molec1.hops[i].hop == molec2.hops[i].hop);
	 }

}

void testZmatrixScanner(Opls_Scan opls){
    cout << endl <<"Testing Z-matrix scanner\n"<< endl;
    string zMatrixFile1 = "../Utilities/bossFiles/mesh.z";
    Molecule meshZ= createMeshZMolecules(opls);
    Zmatrix_Scan zScan (zMatrixFile1,&opls);
    vector<Molecule> scannedInMolecules;
	 int open = zScan.scanInZmatrix();
	 if(open == -1)
	     cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
    else{
	     //for(int x=0; x<5; x++)
	         scannedInMolecules = zScan.buildMolecule(1);
		  }
    //cout << "-- # Scanned in molecules:" << scannedInMolecules.size() <<endl;
    assert(scannedInMolecules.size()==1);
    compareTestMolecules(scannedInMolecules[0],meshZ);  
    cout << "Testing Z-matrix scanner Completed\n"<< endl;  
}

void testZmatrixScanner_multpleSingle(Opls_Scan opls){
    
	 cout << "Testing Z-matrix scanner with reuse 500 molec : 1 atom each\n"<< endl;
    string zMatrixFile1 = "../Utilities/bossFiles/testZ.z";
	 Zmatrix_Scan zScan (zMatrixFile1,&opls);
	 vector<Molecule> scannedInMolecules;
	 int open = zScan.scanInZmatrix();
	 if(open == -1){
	     cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
		  assert(open >= 0 );
	 }
	 Molecule molec;
	 for(int i=0; i<500; i++){
	     scannedInMolecules = zScan.buildMolecule(i);
		   //cout << "-- # Scanned in molecules:" << scannedInMolecules.size() <<endl;
         assert(scannedInMolecules.size()==1);
			molec = scannedInMolecules[0];
			//cout << "\nIndex: " << i << endl;
			//cout << "Molecule Id: "<< molec.id << endl;
			assert(molec.id == i);
			assert(molec.numOfAtoms==1);
			assert(molec.atoms[0].id == i);
			// cout << "# of Atoms: " << molec.numOfAtoms <<endl;
// 			for(int j=0; j< molec.numOfAtoms; j++){
// 			    cout<< "Atom Id: " << molec.atoms[j].id << endl;
// 			}	     
	 } 
	 cout << "Testing Z-matrix scanner with reuse 500 molec : 1 atom each Complete\n"<< endl;	 
}

void testZmatrixScanner_multpleAmount(Opls_Scan opls){
    
	 cout << "Testing Z-matrix scanner with reuse 500 molec : 8 atoms each\n"<< endl;
    string zMatrixFile1 = "../Utilities/bossFiles/mesh.z";
	 Zmatrix_Scan zScan (zMatrixFile1,&opls);
	 vector<Molecule> scannedInMolecules;
	 int open = zScan.scanInZmatrix();
	 if(open == -1){
	     cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
		  assert(open >= 0 );
	 }
	 int runningNumOfAtoms =0;	 
	 Molecule* myMolecArray = new Molecule[500];
	 Molecule molec;
	 for(int i=0; i<500; i++){	
		   scannedInMolecules = zScan.buildMolecule(runningNumOfAtoms);
		   //cout << "-- # Scanned in molecules:" << scannedInMolecules.size() <<endl;
         assert(scannedInMolecules.size()==1);
			myMolecArray[i] = scannedInMolecules[0];
			molec = 	myMolecArray[i] ;
			//cout << "\nIndex: " << i << endl;
			//cout << "Molecule Id: "<< molec.id << endl;		
			//cout << "# of Atoms: " << molec.numOfAtoms <<endl;
			assert(molec.id == i*8);
			assert(molec.numOfAtoms==8);
			for(int j=0; j< molec.numOfAtoms; j++){
			    //cout<< "Atom Id: " << molec.atoms[j].id << endl;
				 assert( molec.atoms[j].id == molec.id+j); 
			}
			runningNumOfAtoms += molec.numOfAtoms;	     
	 } 
	 cout << "Testing Z-matrix scanner with reuse 500 molec : 8 atoms each Complete\n"<< endl;	 
}




int main(){
    runStateTests();
    testConfigScan();    
    
    Opls_Scan scanner(oplsPath);
    int returnInt = scanner.scanInOpls(oplsPath);
    cout << "Attempting to open " << oplsPath << endl;
    if(returnInt == -1){
        cout << "Failed to open Opls file." << endl;
        exit(0);
    }
    testGetAtom(scanner);
    testGetSigma(scanner);
    testGetEpsilon(scanner);
    testGetCharge(scanner);
    testGetFourier(scanner);
    testPDBoutput();
    testZmatrixScanner(scanner);
	testZmatrixScanner_multpleSingle(scanner);
	testZmatrixScanner_multpleAmount(scanner);
    testTranslateAtom();
    testRotateAboutX();
    testRotateAboutY();
    testRotateAboutZ();

}
