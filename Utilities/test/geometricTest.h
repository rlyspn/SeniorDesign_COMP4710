#ifndef GEOMETRICTEST_H
#define GEOMETRICTEST_H

#include "../src/geometricUtil.h"
#include "../src/metroUtil.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

/*
    tests the getNormal function
*/
void testGetNormal();

/*
    tests the getAngle(Plane, Plane) function
*/
void testGetAngleBetweenPlanes();

/*
    tests the getBond function
*/
void testGetBond();

/*
    tests the getAllBonds function
*/
void testGetAllBonds();

/*
    tests the getIntersection function
*/
void testGetIntersection();

/*
    tests the isMember function
*/
void testIsMember();


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
