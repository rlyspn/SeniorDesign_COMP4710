#include "geometricTest.h"

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

