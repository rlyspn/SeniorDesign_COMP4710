#ifndef GEOMETRICTEST_H
#define GEOMETRICTEST_H

#include "../src/geometricUtil.h"
#include "../src/metroUtil.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

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


/**
  tests rotating an atom about a vector.
*/
void testRotateAboutVector();

/**
  tests rotating within a plane
*/
void testRotateInPlane();


#endif //GEOMETRICTEST_H
