#ifndef PARALLELTEST_H
#define PARALLELTEST_H

#include <assert.h>
#include <cuda.h>
#include "baseTests.h"
#include "../src/metroCudaUtil.cuh"
#include <sys/time.h>

/**
Tests the getXFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetXFromIndex();
__global__ void testGetXKernel(int *xValues);

/**
Tests the getYFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetYFromIndex();
__global__ void testGetYKernel();

/**
Tests the makePeriodic function. The global Kernel is used to setup the tests.
*/
void setupMakePeriodic();
__global__ void testMakePeriodicKernel();

/**
Tests the wrapBox Function. The global Kernel is used to setup the tests.
*/
void setupWrapBox();
__global__ void testWrapBoxKernel();

/**
Tests the calc_lj Function. The global Kernel is used to setup the tests.
*/
void setupCalc_lj();
__global__ void testCalcLJKernel();

void testGeneratePoints();

void testCalcEnergy();


#endif //PARALLELTEST_H
