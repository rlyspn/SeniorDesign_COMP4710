#ifndef PARALLELTEST_H
#define PARALLELTEST_H

#include <assert.h>
#include <cuda.h>
#include "baseTests.h"
#include "metroCudaUtil.cuh"

/**
Tests the getXFromIndex function. The global Kernel is used to setup the tests.
*/
void testGetXFromIndex();
__global__ void testGetXKernel();

/**
Tests the getYFromIndex function. The global Kernel is used to setup the tests.
*/
void testGetYFromIndex();
__global__ void testGetXKernel();

/**
Tests the makePeriodic function. The global Kernel is used to setup the tests.
*/
void testMakePeriodic();
__global__ void testGetYKernel();

/**
Tests the wrapBox Function. The global Kernel is used to setup the tests.
*/
void testWrapBox();
__global__ void testMakePeriodicKernel();

/**
Tests the calc_lj Function. The global Kernel is used to setup the tests.
*/
void testCalc_lj();
__global__ void calcLJKernel();

void testGeneratePoints();

void testCalcEnergy();


#endif //PARALLELTEST_H
