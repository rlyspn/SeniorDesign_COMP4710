#include "../src/Opls_Scan.h"
#include <assert.h>
#include <iostream>

using namespace std;

string oplsPath = "../Utilities/bossFiles/oplsua.par";

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
    
    printf("%f, %f, %f\n", atom1.charge, atom1.sigma, atom1.epsilon);
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

int main(){
    Opls_Scan scanner(oplsPath);
    cout << scanner.scanInOpls(oplsPath) << endl;
    cout << "Reading file: " << oplsPath << endl;
    
    testGetAtom(scanner);
    testGetSigma(scanner);
    testGetEpsilon(scanner);
    testGetCharge(scanner);
    testGetFourier(scanner);
}
