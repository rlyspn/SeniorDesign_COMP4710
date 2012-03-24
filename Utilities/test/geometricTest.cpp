#include "geometricTest.h"


void testGetNormal(){
    Atom atom1, atom2, atom3;
    atom1.x = 3.2;
    atom1.y = 4.6;
    atom1.z = 0.2;
    atom2.x = 7.1;
    atom2.y = 1.4;
    atom2.z = 10.0;
    atom3.x = 0.0;
    atom3.y = 0.0;
    atom3.z = 0.0;

    Plane testPlane = createPlane(atom1, atom2, atom3);

    Point expected = createPoint(45.72, -30.58, -28.18);

    Point test = getNormal(testPlane);
    
    test.x = ((double) ((int) (test.x * 100))) / 100.0;
    test.y = ((double) ((int) (test.y * 100))) / 100.0;
    test.z = ((double) ((int) (test.z * 100))) / 100.0;
    assert(test.x == expected.x);
    assert(test.y == expected.y);
    assert(test.z == expected.z);

    printf("testGetNormal completed successfully.\n");
}

void testGetAngleBetweenPlanes(){
    //TODO after getAngle test is written
}

void testGetBond(){
    vector<Bond> bonds;

    Bond bond1 = createBond(1, 2, 3.0, false);
    Bond bond2 = createBond(3, 4, 3.0, false);
    Bond bond3 = createBond(2, 3, 3.0, false);

    bonds.push_back(bond1);
    bonds.push_back(bond2);
    bonds.push_back(bond3);

    Bond testBond1 = getBond(bonds, 1, 2);
    Bond testBond2 = getBond(bonds, 3, 4);
    Bond testBond3 = getBond(bonds, 2, 3);
    Bond testBond4 = getBond(bonds, 1, 4);

    assert(bond1.atom1 == testBond1.atom1 && bond1.atom2 == testBond1.atom2 && bond1.distance == testBond1.distance);
    assert(bond2.atom1 == testBond2.atom1 && bond2.atom2 == testBond2.atom2 && bond2.distance == testBond2.distance);
    assert(bond3.atom1 == testBond3.atom1 && bond3.atom2 == testBond3.atom2 && bond3.distance == testBond3.distance);
    assert(testBond4.atom1 == -1 && testBond4.atom2 == -1 && testBond4.distance == -1.0);

    printf("testGetBond completed successfully.\n");
}

void testGetAllBonds(){

}

void testGetIntersection(){

}

void testIsMember(){

}

void testTranslateAtom(){
    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++){
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i){
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double oldX = testAtom.x;
    double oldY = testAtom.y;
    double oldZ = testAtom.z;

    double translateX = ((double) rand() / RAND_MAX) * 15; 
    double translateY = ((double) rand() / RAND_MAX) * 15; 
    double translateZ = ((double) rand() / RAND_MAX) * 15; 
    
    testAtom = translateAtom(testAtom, translateX, translateY, translateZ);

    assert(testAtom.x == oldX + translateX);
    assert(testAtom.y == oldY + translateY);
    assert(testAtom.z == oldZ + translateZ);

    cout << "testTranslateAtom completed successfully" << endl;    
}

void testRotateAboutX(){
    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++){
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i){
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldX = testAtom.x;
    double oldY = testAtom.y;
    double oldZ = testAtom.z;
 
    double theta = ((double) rand() / RAND_MAX) * 180;
    
    testAtom = rotateAboutX(testAtom, theta);
    
    assert(testAtom.y == cos(theta*dtr) * oldY + sin(theta*dtr) * oldZ);
    assert(testAtom.z == cos(theta*dtr) * oldZ - sin(theta*dtr) * oldY);


    cout << "testRotateAboutX completed successfully" << endl;    
}

void testRotateAboutY(){

    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++){
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i){
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldX = testAtom.x;
    double oldY = testAtom.y;
    double oldZ = testAtom.z;
 
    double theta = ((double) rand() / RAND_MAX) * 180;
    
    testAtom = rotateAboutY(testAtom, theta);
    
    assert(testAtom.x == cos(theta*dtr) * oldX - sin(theta*dtr) * oldZ);
    assert(testAtom.z == cos(theta*dtr) * oldZ + sin(theta*dtr) * oldX);


    cout << "testRotateAboutY completed successfully" << endl;    
}

void testRotateAboutZ(){

    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++){
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i){
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldX = testAtom.x;
    double oldY = testAtom.y;
    double oldZ = testAtom.z;
 
    double theta = ((double) rand() / RAND_MAX) * 180;
    
    testAtom = rotateAboutZ(testAtom, theta);
    
    assert(testAtom.x == cos(theta*dtr) * oldX + sin(theta*dtr) * oldY);
    assert(testAtom.y == cos(theta*dtr) * oldY - sin(theta*dtr) * oldX);


    cout << "testRotateAboutZ completed successfully" << endl;    
}

void testRotateAboutVector(){
    double rotation = 90;
    Atom vertex = createAtom(-1, 0, 0, 0);
    Atom head = createAtom(-1, 1, 0, 0);
    Atom toRotate = createAtom(-1, 0, 1, 0);

    printAtoms(&toRotate, 1);
    Atom rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - 0.0) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - 1.0) < .01);

    rotation = 45;
    
    rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - sqrt(.5)) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - sqrt(.5)) < .01);


    printf("rotateAtomAboutVectorPassed()\n");
}

void testRotateInPlane(){
     double rotation = 90;
    Atom vertex = createAtom(-1, 0, 0, 0);
    Atom head = createAtom(-1, 0, 0, -1);
    Atom toRotate = createAtom(-1, 0, 1, 0);

    printAtoms(&toRotate, 1);
    Atom rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - 0.0) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - 1.0) < .01);
 
    rotation = 45;
    
    rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - sqrt(.5)) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - sqrt(.5)) < .01);
  
    printf("rotateAtomInPlanePassed()\n");
}
