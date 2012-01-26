#ifndef PARALLELTEST_H
#define PARALLELTEST_H

#include <assert.h>
#include <cuda.h>
#include "baseTests.h"
#include "../src/metroCudaUtil.cuh"

/**
Tests the getXFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetXFromIndex();


/**
Tests the getYFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetYFromIndex();

/**
Tests the makePeriodic function. The global Kernel is used to setup the tests.
*/
void setupMakePeriodic();

/**
Tests the wrapBox Function. The global Kernel is used to setup the tests.
*/
void setupWrapBox();

/**
Tests the calc_lj Function. The global Kernel is used to setup the tests.
*/
void setupCalc_lj();

void testGeneratePoints();

void testCalcEnergy();


#endif //PARALLELTEST_H
