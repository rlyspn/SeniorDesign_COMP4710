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
string oplsPath = "bin/oplsua.par";

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
    cout << "Testing Opls_Scan.getAtom" << endl;
    Atom atom1 = scan.getAtom(atomNumber1);
    Atom atom2 = scan.getAtom(atomNumber2);
    Atom atom3 = scan.getAtom(atomNumber3);
    
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
    cout << "Testing Opls_Scan.getSigma" << endl;
    assert(scan.getSigma(atomNumber1) == atom1Sigma); 
    assert(scan.getSigma(atomNumber2) == atom2Sigma); 
    assert(scan.getSigma(atomNumber3) == atom3Sigma); 
    cout << "Testing Opls_Scan.getSigma Completed\n" << endl;
}

void testGetEpsilon(Opls_Scan scan){
    cout << "Testing Opls_Scan.getEpsilon" << endl;
    assert(scan.getEpsilon(atomNumber1) == atom1Epsilon);
    assert(scan.getEpsilon(atomNumber2) == atom2Epsilon);
    assert(scan.getEpsilon(atomNumber3) == atom3Epsilon);
    cout << "Testing Opls_Scan.getEpsilon Completed\n" << endl;
}

void testGetCharge(Opls_Scan scan){
    cout << "Testing Opls_Scan.getCharge" << endl;
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
    cout << "Testing Opls_Scan.getVvalues" << endl;
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
    cout << "Testing PDBoutput" << endl;
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

    cout << "Testing PDBoutput Completed\n" << endl;
}

void testLogOutput(){
    cout << "Testing OutputLog writer" <<endl;
    
	 system("find ../ -name OutputLog | xargs rm");
    string line1 = "This is line1 text";
	 string line2 = "This is line2 text";
	 string line3 = "This is line3 text";
	 
	 writeToLog(line1);
	 writeToLog(line2,Z_MATRIX);
	 writeToLog(line3,OPLS);


    string outPutFile= "OutputLog";
	 ifstream fileReader;
	 fileReader.open(outPutFile.c_str());
	 assert( fileReader.good());
	 
	 string readInLine;
	 getline(fileReader,readInLine);
	 assert(readInLine.compare(line1)==0);
	 
	 getline(fileReader,readInLine);
	 string temp= "--Z_Matrix: ";
	 temp+=line2;
	 assert(readInLine.compare(temp)==0);
	 //getline(fileReader,readInLine);
	 //assert(readInLine.compare(line2)==0);
	 
	 getline(fileReader,readInLine);
	 temp= "--OPLS: ";
	 temp+=line3;
	 assert(readInLine.compare(temp)==0);
	 //getline(fileReader,readInLine);
	 //assert(readInLine.compare(line3)==0);
	 
	 cout << "Testing OutputLog writer Completed\n" <<endl;	 
}


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
    Bond bond5=createBond(5,1,1.811119,true);
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

    Hop hop1= createHop(2,6,3);
    Hop hop2= createHop(2,7,3);
    Hop hop3= createHop(2,8,3);
    Hop hop4= createHop(3,4,3);
    Hop hop5= createHop(3,5,3);
    Hop hop6= createHop(3,6,4);
    Hop hop7= createHop(3,7,4);
    Hop hop8= createHop(3,8,4);
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
    hopPtr[0]=hop1; hopPtr[1]=hop2; hopPtr[2]=hop3; hopPtr[3]=hop4; hopPtr[4]=hop5; hopPtr[5]=hop6;  hopPtr[6]=hop7; hopPtr[7]=hop8; hopPtr[8]=hop9; 
    hopPtr[9]=hop10; hopPtr[10]=hop11;


    return createMolecule(1,atomPtr,anglePtr,bondPtr,dihedPtr,hopPtr,8,6,7,5,11);
}

void addPositionsToMeshZ(Molecule *meshMolec){
    //the hardcoded positions(x y z) of the mezh.z molecule
    // Atom#   X    Y    Z
	 // 0       0    0    0
	 meshMolec->atoms[0].x = 0;
	 meshMolec->atoms[0].y = 0;
    meshMolec->atoms[0].z = 0;
	 // 1       0    0.5    0
	 meshMolec->atoms[1].x = 0;
	 meshMolec->atoms[1].y = 0.5;
    meshMolec->atoms[1].z = 0;
	 // 2     -0.5    0.5   -0.5
	 meshMolec->atoms[2].x = -0.5;
	 meshMolec->atoms[2].y = 0.5;
    meshMolec->atoms[2].z = -0.5;
	 // 3    1.33653   2.39894e-9   1.33653
	 meshMolec->atoms[3].x = 1.33653;
	 meshMolec->atoms[3].y = 0.00000000239894;
    meshMolec->atoms[3].z = 1.33653;
	 // 4   -0.857136    -1.36313   -0.823996
	 meshMolec->atoms[4].x = -0.857136 ;
	 meshMolec->atoms[4].y = -1.36313;
    meshMolec->atoms[4].z = -0.823996;
	 // 5   -1.92671    -1.51275   -0.405374
	 meshMolec->atoms[5].x = -1.92671;
	 meshMolec->atoms[5].y = -1.51275;
    meshMolec->atoms[5].z = -0.405374;
	 // 6   -0.878138    -1.16034   -1.94557
	 meshMolec->atoms[6].x = -0.878138;
	 meshMolec->atoms[6].y = -1.16034;
    meshMolec->atoms[6].z = -1.94557;
	 // 7    -0.244983    -2.314   -0.710251
	 meshMolec->atoms[7].x = -0.244983;
	 meshMolec->atoms[7].y = -2.314;
    meshMolec->atoms[7].z = -0.710251;
}

void compareTestMolecules(Molecule molec1, Molecule molec2){
    // check if id's are equal
    assert(molec1.id == molec2.id);	 

    //check atoms
    assert(molec1.numOfAtoms == molec2.numOfAtoms);
    for(int i=0; i< molec1.numOfAtoms; i++){
        assert(molec1.atoms[i].id==molec2.atoms[i].id);
        assert(percentDifference(molec1.atoms[i].charge, molec2.atoms[i].charge));
        assert(percentDifference(molec1.atoms[i].sigma, molec2.atoms[i].sigma));
        assert(percentDifference(molec1.atoms[i].epsilon, molec2.atoms[i].epsilon));
        assert(percentDifference(molec1.atoms[i].x, molec2.atoms[i].x));
        assert(percentDifference(molec1.atoms[i].y, molec2.atoms[i].y));
        assert(percentDifference(molec1.atoms[i].z, molec2.atoms[i].z));			
    }

    //check bond
    assert(molec1.numOfBonds == molec2.numOfBonds);
    for(int i=0; i< molec1.numOfBonds; i++){
        assert(molec1.bonds[i].atom1 == molec2.bonds[i].atom1);
        assert(molec1.bonds[i].atom2 == molec2.bonds[i].atom2);
        assert(percentDifference(molec1.bonds[i].distance, molec2.bonds[i].distance));
        assert(asserTwoBool( molec1.bonds[i].variable,molec2.bonds[i].variable));
    }

    //check angles
    assert(molec1.numOfAngles == molec2.numOfAngles);
    for(int i=0; i< molec1.numOfAngles; i++){
        assert(molec1.angles[i].atom1 == molec2.angles[i].atom1);
        assert(molec1.angles[i].atom2 == molec2.angles[i].atom2);
        assert(percentDifference(molec1.angles[i].value,molec2.angles[i].value));
        assert(asserTwoBool(molec1.angles[i].variable, molec2.angles[i].variable));
    }

    //check dihederals
    assert(molec1.numOfDihedrals == molec2.numOfDihedrals);
    for(int i=0; i< molec1.numOfDihedrals; i++){
        assert(molec1.dihedrals[i].atom1 == molec2.dihedrals[i].atom1);
        assert(molec1.dihedrals[i].atom2 == molec2.dihedrals[i].atom2);
        assert(percentDifference(molec1.dihedrals[i].value,molec2.dihedrals[i].value));
        assert(asserTwoBool(molec1.dihedrals[i].variable,molec2.dihedrals[i].variable));
    }

    //check hops
    assert(molec1.numOfHops == molec2.numOfHops);
    for(int i=0; i< molec1.numOfHops; i++){
        assert(molec1.hops[i].atom1 == molec2.hops[i].atom1);
        assert(molec1.hops[i].atom2 == molec2.hops[i].atom2);
        assert(molec1.hops[i].hop == molec2.hops[i].hop);
    }

}

void testZmatrixScanner(Opls_Scan opls){
    cout << "Testing Z-matrix scanner"<< endl;
    string zMatrixFile1 = "Utilities/bossFiles/mesh.z";
    Molecule meshZ= createMeshZMolecules(opls);	 
	 addPositionsToMeshZ(&meshZ);
	 
    Zmatrix_Scan zScan (zMatrixFile1,&opls);
    vector<Molecule> scannedInMolecules;

    int open = zScan.scanInZmatrix();
    if(open == -1)
        cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
    else{
        scannedInMolecules = zScan.buildMolecule(1);
    }

    assert(scannedInMolecules.size()==1);
    compareTestMolecules(scannedInMolecules[0],meshZ);  
    cout << "Testing Z-matrix scanner Completed\n"<< endl;  
}

void testZmatrixScanner_multpleSingle(Opls_Scan opls){
    cout << "Testing Z-matrix scanner with reuse 500 molec : 1 atom each"<< endl;
    string zMatrixFile1 = "Utilities/bossFiles/testZ.z";
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
        assert(scannedInMolecules.size()==1);
        molec = scannedInMolecules[0];
        assert(molec.id == i);
        assert(molec.numOfAtoms==1);
        assert(molec.atoms[0].id == i);
    }
     
    cout << "Testing Z-matrix scanner with reuse 500 molec : 1 atom each Complete\n"<< endl;	 
}

void testZmatrixScanner_multpleAmount(Opls_Scan opls){
    cout << "Testing Z-matrix scanner with reuse 500 molec : 8 atoms each"<< endl;
    string zMatrixFile1 =  "Utilities/bossFiles/mesh.z";
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
        assert(scannedInMolecules.size()==1);

        myMolecArray[i] = scannedInMolecules[0];
        molec = myMolecArray[i] ;
        assert(molec.id == i*8);
        assert(molec.numOfAtoms==8);

        for(int j=0; j< molec.numOfAtoms; j++){
            assert( molec.atoms[j].id == molec.id+j); 
        }
        runningNumOfAtoms += molec.numOfAtoms;	     
    }
    cout << "Testing Z-matrix scanner with reuse 500 molec : 8 atoms each Complete\n"<< endl;	 
}



int main(){
    cout<< "----- Starting Utility Test ----\n" << endl;
    runStateTests();
    testConfigScan();    
    
    Opls_Scan scanner(oplsPath);
    int returnInt = scanner.scanInOpls(oplsPath);
    cout << "--Attempting to open " << oplsPath << endl;
    if(returnInt == -1){
        cout << "Failed to open Opls file." << endl;
        exit(0);
    }

    //test OPLS
    testGetAtom(scanner);
    testGetSigma(scanner);
    testGetEpsilon(scanner);
    testGetCharge(scanner);
    testGetFourier(scanner);

    testPDBoutput();
	 testLogOutput();
	 
	 //Test Zmatrix
    testZmatrixScanner(scanner); 
	 testZmatrixScanner_multpleSingle(scanner);
    testZmatrixScanner_multpleAmount(scanner);

    testGeometric();
}
