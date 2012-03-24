#ifndef GEOMETRICTEST_H
#define GEOMETRICTEST_H

#include "../src/geometricUtil.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

/**
 test degrees to radians and radians to degrees
*/
void testD2RandR2D();

/**
 test getOpositeAtom
*/
void testGetOppositeAtom();

/**
 test getOpositeAtom
*/
void testGetCommonAtom();

/**
 test getDistance
*/
void testGetDistance();

/**
 test getAngle
*/
void testGetAngle();

/**
    tests the translateAtom function
*/
void testTranslateAtom();

/**
    tests the rotateAboutX function
*/
void testRotateAboutX();

/**
    tests the rotateAboutY function
*/
void testRotateAboutY();

/**
    tests the rotateAboutZ function
*/
void testRotateAboutZ();

 

#endif //GEOMETRICTEST_H
