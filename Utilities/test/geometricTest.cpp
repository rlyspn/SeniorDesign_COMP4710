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
    //TODO by alex
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
    vector<Bond> bonds;

    Bond bond1 = createBond(1, 2, 3.0, false);
    Bond bond2 = createBond(3, 4, 3.0, false);
    Bond bond3 = createBond(2, 3, 3.0, false);

    bonds.push_back(bond1);
    bonds.push_back(bond2);
    bonds.push_back(bond3);

    vector<unsigned long> testBonds1 = getAllBonds(bonds, 1);
    vector<unsigned long> testBonds2 = getAllBonds(bonds, 2);
    vector<unsigned long> testBonds3 = getAllBonds(bonds, 3);
    vector<unsigned long> testBonds4 = getAllBonds(bonds, 4);
    vector<unsigned long> testBonds5 = getAllBonds(bonds, 5);

    assert(testBonds1[0] == 2 && testBonds1.size() == 1);
    assert(testBonds2[0] == 1 && testBonds2[1] == 3 && testBonds2.size() == 2);
    assert(testBonds3[0] == 4 && testBonds3[1] == 2 && testBonds3.size() == 2);
    assert(testBonds4[0] == 3 && testBonds4.size() == 1);
    assert(testBonds5.size() == 0);

    printf("testGetAllBonds completed successfully.\n");
}

void testGetIntersection(){
    vector<unsigned long> section1, section2;

    section1.push_back(1);
    section1.push_back(2);
    section1.push_back(3);
    section1.push_back(4);

    section2.push_back(2);
    section2.push_back(4);
    section2.push_back(6);
    section2.push_back(8);

    vector<unsigned long> intersection = getIntersection(section1, section2);

    assert(intersection[0] == 2 && intersection[1] == 4 && intersection.size() == 2);

    printf("testGetIntersection completed successfully.\n");

}

void testIsMember(){
    vector<unsigned long> section;

    section.push_back(1);
    section.push_back(2);
    section.push_back(3);
    section.push_back(4);

    assert(isMember(section, 1));
    assert(isMember(section, 2));
    assert(isMember(section, 3));
    assert(isMember(section, 4));
    assert(!isMember(section, 5));

    printf("testIsMember completed successfully.\n");
}

bool percentDifferenceG(double d1, double d2){
    double difference = d2-d1;
    double average = (d2+d1)/d2;
    double percentDiff = (difference/average)*100;
    //cout <<"Percent Diff: " << percentDiff << endl;;
    return percentDiff < 3;
}


void testD2RandR2D(){
    //test converting Degrees to Radians
    assert(percentDifferenceG(degreesToRadians(1),0.0174532925) );
	 assert(percentDifferenceG(degreesToRadians(45),0.785398163) );
	 assert(percentDifferenceG(degreesToRadians(254),4.4331363) );
	 assert(percentDifferenceG(degreesToRadians(360),6.283185307) );
	 assert(percentDifferenceG(degreesToRadians(15),0.261799388) );
	 
	 //test converting Radians to Degrees
	 assert(percentDifferenceG( radiansToDegrees(3.14),179.908747671) );
	 assert(percentDifferenceG( radiansToDegrees(.1234567),7.073547863) );
	 assert(percentDifferenceG( radiansToDegrees(0.174532925),10) );
	 assert(percentDifferenceG( radiansToDegrees(1.745329252),100) );
	 assert(percentDifferenceG( radiansToDegrees(6.195918845),355) );
	 
	 cout << "test degrees2radins and radian2degrees complete" << endl;

}

void testGetOppositeAtom(){
    Bond testBond = createBond(1,3,5.5,true);
	 Angle testAngle = createAngle(2,6,30,false);
	 Dihedral testDihed = createDihedral(3,7,180,true);
	 
	 //test with bonds
	 assert(getOppositeAtom(testBond,1)==3);
	 assert(getOppositeAtom(testBond,3)==1);
	 assert(getOppositeAtom(testBond,4)==-1);
	 
	 //test with bonds
	 assert(getOppositeAtom(testAngle,2)==6);
	 assert(getOppositeAtom(testAngle,6)==2);
	 assert(getOppositeAtom(testAngle,5)==-1);
	 
	 //test with bonds
	 assert(getOppositeAtom(testDihed,3)==7);
	 assert(getOppositeAtom(testDihed,7)==3);
	 assert(getOppositeAtom(testDihed,6)==-1);
	 
	 testBond = createBond(11,12,3.5,true);
	 testAngle = createAngle(1,22,30,false);
	 testDihed = createDihedral(5,9,180,true);
	 
	 //test with bonds
	 assert(getOppositeAtom(testBond,11)==12);
	 assert(getOppositeAtom(testBond,12)==11);
	 assert(getOppositeAtom(testBond,13)==-1);
	 
	 //test with bonds
	 assert(getOppositeAtom(testAngle,1)==22);
	 assert(getOppositeAtom(testAngle,22)==1);
	 assert(getOppositeAtom(testAngle,15)==-1);
	 
	 //test with bonds
	 assert(getOppositeAtom(testDihed,5)==9);
	 assert(getOppositeAtom(testDihed,9)==5);
	 assert(getOppositeAtom(testDihed,62)==-1);
	 
	 cout << "test getOppositeAtom complete" << endl;

}

void testGetCommonAtom(){
    vector<Bond> bondVect(7); 
    bondVect[0] =createBond(2,1,0.5,false);
	 bondVect[1] =createBond(3,2,0.5,false);
	 bondVect[2] =createBond(4,1,1.336532,true);
	 bondVect[3] =createBond(5,1,1.8119,true);
	 bondVect[4] =createBond(6,5,1.090187,true);
	 bondVect[5] =createBond(7,5,1.090135,true);
	 bondVect[6] =createBond(8,5,1.090135,true);
	 
	 assert(getCommonAtom(bondVect,4,5)==1 );
	 assert(getCommonAtom(bondVect,1,3)==2 );
	 assert(getCommonAtom(bondVect,6,8)==5 );
	 assert(getCommonAtom(bondVect,2,5)==1 );
	 assert(getCommonAtom(bondVect,3,5)==-1 );
	 assert(getCommonAtom(bondVect,8,2)==-1 );
	 
	 cout << "test getCommonAtom complete" << endl;
}

void testGetDistance(){
    Atom atom1= createAtom(1,5,6,7);
	 Atom atom2= createAtom(2,10,11,12);
	 assert(percentDifferenceG(getDistance(atom1,atom2),8.6) );
	 
	 atom1= createAtom(1,8,12,21);
	 atom2= createAtom(2,4,5,10);
	 assert(percentDifferenceG(getDistance(atom1,atom2),13.638181) );
	 
	 atom1= createAtom(1,45,2,22);
	 atom2= createAtom(2,37,22,18);
	 assert(percentDifferenceG(getDistance(atom1,atom2),21.9089023002) );
	 
	 cout << "test getDistance complete" << endl;
}

void testGetAngle(){
    Atom atom1= createAtom(1,5,6,7);
	 Atom atom2= createAtom(2,10,11,12);
	 Atom atom3= createAtom(3,14,22,9);
	 assert(percentDifferenceG(getAngle(atom1,atom2,atom3),124.986));
	 
	 atom1= createAtom(1,15,23,8);
	 atom2= createAtom(2,5,3,12);
	 atom3= createAtom(3,9,18,7);
	 assert(percentDifferenceG(getAngle(atom1,atom2,atom3),13.6609));
	 
	 cout << "test getAngle complete" << endl;
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
   
    double expectedY = cos(theta * dtr) * oldY + sin(theta * dtr) * oldZ;
    double expectedZ = cos(theta * dtr) * oldZ - sin(theta * dtr) * oldY;

    printf("expectedY = %f\nactual Y = %f\n", expectedY, testAtom.y);
    printf("expectedZ = %f\nactual Z = %f\n", expectedZ, testAtom.z);

    assert((testAtom.y - expectedY) / expectedY < PRECISION);
    assert((testAtom.z - expectedZ) / expectedZ < PRECISION);

    /*assert(testAtom.y == cos(theta*dtr) * oldY + sin(theta*dtr) * oldZ);
    assert(testAtom.z == cos(theta*dtr) * oldZ - sin(theta*dtr) * oldY);*/


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
    
    double expectedX = cos(theta * dtr) * oldX - sin(theta * dtr) * oldZ;
    double expectedZ = cos(theta * dtr) * oldZ + sin(theta * dtr) * oldX;
    
    assert((testAtom.x - expectedX) / expectedX < PRECISION);
    assert((testAtom.z - expectedZ) / expectedZ < PRECISION);
    
    /*assert(testAtom.x == cos(theta*dtr) * oldX - sin(theta*dtr) * oldZ);
    assert(testAtom.z == cos(theta*dtr) * oldZ + sin(theta*dtr) * oldX);*/


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
    
    double expectedX = cos(theta * dtr) * oldX + sin(theta * dtr) * oldY;
    double expectedY = cos(theta * dtr) * oldY - sin(theta * dtr) * oldX;
    
    assert((testAtom.x - expectedX) / expectedX < PRECISION);
    assert((testAtom.y - expectedY) / expectedY < PRECISION);
    
    /*assert(testAtom.x == cos(theta*dtr) * oldX + sin(theta*dtr) * oldY);
    assert(testAtom.y == cos(theta*dtr) * oldY - sin(theta*dtr) * oldX);*/


    cout << "testRotateAboutZ completed successfully" << endl;    
}

void testRotateAboutVector(){
    double rotation = 90;
    Atom vertex = createAtom(-1, 0, 0, 0);
    Atom head = createAtom(-1, 1, 0, 0);
    Atom toRotate = createAtom(-1, 0, 1, 0);
    
    printf("Testing At 90 degrees.\n");
    printAtoms(&toRotate, 1);
    Atom rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - 0.0) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - 1.0) < .01);

    rotation = 45;
    
    printf("Testing At 45 degrees.\n");
    rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - sqrt(.5)) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - sqrt(.5)) < .01);

    rotation = 90;
    //test rotating about atom with 0 intial angle
    printf("Before rotation\n");
    toRotate = createAtom(-1,2,0,0);
    printAtoms(&toRotate, 1);
    rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    printf("After rotation\n");
    printAtoms(&rotated, 1);
    assert(fabs(rotated.x - 2) < .01);
    assert(fabs(rotated.y - 0) < .01);
    assert(fabs(rotated.z - 0) < .01);


    printf("rotateAtomAboutVectorPassed()\n");
}

void testRotateInPlane(){
    printf("\nTesting Rotation in Plane.\n");
     double rotation = 90;
    Atom vertex = createAtom(-1, 0, .5, 0);
    Atom head = createAtom(-1, 0, 0, 0);
    Atom toRotate = createAtom(-1, 0, 1, 0);
    
    printAtoms(&toRotate, 1);
    printAtoms(&vertex, 1);
    printAtoms(&head, 1);


    Atom rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);
    
    printAtoms(&rotated, 1);
    printAtoms(&vertex, 1);
    printAtoms(&head, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - .5) < .01);
 
    rotation = 45;
    
    printf("rotating 45 degrees\n");

    vertex = createAtom(-1, 0, 0, 0);
    head = createAtom(-1, 0, 0, -1);
    toRotate = createAtom(-1, 0, 1, 0);

    rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);
    printAtoms(&rotated, 1);
    printf("Angle = %f\n", getAngle(rotated, vertex, toRotate));
    assert(fabs(rotated.y - sqrt(.5)) < .01);
    assert(fabs(rotated.x - 0.0) < .01);
    assert(fabs(rotated.z - sqrt(.5)) < .01);
  
    printf("rotateAtomInPlanePassed()\n");
}
