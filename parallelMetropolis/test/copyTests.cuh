#ifndef COPYTEST_H
#define COPYTEST_H

#include <assert.h>
#include "../src/metroCudaUtil.cuh"

/**
  Tests the moleculeDeepCopyToDevice() function
  Tests the moleculeDeepCopyToHost() function.
*/
void testCopyMolecules();

/**
  Tests the allocateOnDeviceFunction
*/
void testAllocateMemory();

/**
  tests the freeMOleculeOnDevice() function.
*/
void testFreeMemory();


#endif // COPYTEST_H
