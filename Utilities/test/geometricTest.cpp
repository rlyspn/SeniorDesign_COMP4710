#include "geometricTest.h"

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

void testRotateAtom(){
    //TODO
}

