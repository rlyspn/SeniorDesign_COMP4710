#ifndef STATETEST_H
#define STATETEST_H

#include "../src/metroUtil.h"
#include <assert.h>
#include <iostream>

/**
  tests the getDihedralFromLine function
*/
void testGetDihedralFromLine();

/**
  tests the getAngleFromLine function
*/
void testGetAngleFromLine();

/**
  tests the getBondFromLine function
*/
void testGetBondFromLine();

/**
  tests the getEnvironmentFromLine function
*/
void testGetEnvironmentFromLine();

/**
  tests the getAtomFromLine function
*/
void testGetAtomFromLine();

/**
  tests the getEnvironmentFromLine function
*/
void testGetEnvironmentFromLine();

/**
  runs the tests associated with the state reading/writing
*/
void runStateRead_WriteTests();
#endif //STATETEST_H