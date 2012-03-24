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
