#include "../src/Opls_scan.h"
#include <assert.h>
#include <iostream>

using namespace std;
int atomNumber1 = 420;
int atomNumber2 = 1;
int atomNumber3 = 53;

double atom1Charge = -2.0;
double atom1Sigma = 0.5;
double atom1Epsilon = 3.25;

double atom2Charge = 0.5;
double atom2Sigma = 3.75;
double atom2Epsilon = 0.105;

double atom3Charge = .640;
double atom3Sigma = 2.250;
double atom3Epsilon = 0.50;

void testGetAtom(Opls_scan scan){
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
}

void testGetSigma(Opls_scan scan){
    assert(scan.getSigma(atomNumber1) == atom1Sigma); 
    assert(scan.getSigma(atomNumber2) == atom2Sigma); 
    assert(scan.getSigma(atomNumber3) == atom3Sigma); 
}

void testGetEpsilon(Opls_scan scan){
    assert(scan.getEpsilon(atomNumber1) == atom1Epsilon);
    assert(scan.getEpsilon(atomNumber2) == atom2Epsilon);
    assert(scan.getEpsilon(atomNumber3) == atom3Epsilon);
}

void testGetCharge(Opls_scan scan){
    assert(scan.getCharge(atomNumber1) == atom1Charge);
    assert(scan.getCharge(atomNumber2) == atom2Charge);
    assert(scan.getCharge(atomNumber3) == atom3Charge);
}

int main(){
    string fileName = "../bossFiles/oplsua.par";
    Opls_Scan scanner(fileName);

    testGetAtom(scanner);
    testGetSigma(scanner);
    testGetEpsilon(scanner);
    testGetCharge(scanner);
}
