#ifndef METROCUDAUTIL_CUH
#define METROCUDAUTIL_CUH

#define DEBUG

#include <cuda.h>
#include <curand_kernel.h>
#include <math.h>
#include "metroParallelUtil.h"
#include <curand.h>

/**
  @param idx - the index in the 1 dimensional array of energies
  @return - the id of the X atom
*/
__device__ int getXFromIndex(int idx);

/**
  @param x - the id of the x atom
  @param idx - the index in the 1 dimensional array of energies
*/
__device__ int getYFromIndex(int x, int idx);

/**
  @param x - the variable to make periodic
  @param box - the size of the period
  @return - x after being made periodic
*/
__device__ double makePeriodic(double x, const double box);

/**
  @param x - the value to continue on the other side of the box
  @param box - the length of one side of the box (cube)
*/
__device__ double wrapBox(double x, double box);

/**
  Calculates the energy between 2 atoms
  @param atom1 - the first atom in the pair
  @param atom2 - the second atom in the pair
  @param enviro - the environmental variables
*/
__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro); 

/**
    Global function for GPU to assign random doubles to atom positions.
    @param dev_doubles - array of randomly generated doubles 0.0 to 1.0
    @param atoms - array of atoms to be positioned
    @param enviro - Environment of the simulation
*/
__global__ void assignAtomPositions(double *dev_doubles, Atom *atoms, Environment *enviro);

/**
  Generate random positions for atoms in the box
  nVidia CURAND reference: http://developer.download.nvidia.com/compute/cuda/5_0/toolkit/docs/CURAND_Library.pdf
  @param atoms - array of atoms to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatePoints(Atom *atoms, Environment *enviro);

/**
  This is a wrapper function for the calcEnergy kernel.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @return - the total energy of the system.
*/
double calcEnergyWrapper(Atom *atoms, Environment *enviro);

/**
  Calculates the energy between n atoms where n is the
  the number of threads in the block. The block's sum is then stored 
  in the energySum array at its block index position.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @param *energySum - the array of block sums
*/
__global__ void calcEnergy(Atom *atoms, Environment *enviro, double *energySum);

#ifdef DEBUG
/**
    Kernel call that will be used to test the getXFromIndexFunction
    @apram xValues - a series of xValues used to test
*/
__global__ void testGetXKernel(int *xValues, int totalTests);

/**
  Kernel function that will be used to test getYFromIndex device function
*/
__global__ void testGetYKernel(int *xValues, int *yValues, int numberOfTests);

/**
  Kernel function used to facilitate tests of the makePeriodic device function
*/
__global__ void testMakePeriodicKernel(double *x, double *box, int numberOfTests);

/**
  Kernel function used to facilitate tess of the wrap box device function
*/
__global__ void testWrapBoxKernel(double *x, double *box, int numberOfTests);

#endif //DEBUG

#endif //METROCUDAUTIL_H
