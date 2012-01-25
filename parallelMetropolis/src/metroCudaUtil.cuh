#ifndef METROCUDAUTIL_CUH
#define METROCUDAUTIL_CUH

#include <cuda.h>
#include <curand_kernel.h>
#include <math.h>
#include "metroParallelUtil.h"

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
  Initializes CURAND random number generator
  Borrowed from: http://aresio.blogspot.com/2011/05/cuda-random-numbers-inside-kernels.html
  @param globalState - state of the generator
  @param seed - seed for the generator
*/
__global__ void setup_generator(curandState *globalState, unsigned long seed);

/**
  Generate random positions for atoms in the box
  Some code borrowed from: http://aresio.blogspot.com/2011/05/cuda-random-numbers-inside-kernels.html
  nVidia CURAND reference: http://developer.download.nvidia.com/compute/cuda/5_0/toolkit/docs/CURAND_Library.pdf
  @param globalState - the state of the generator
  @param atoms - array of atoms to generate positions
  @param enviro - enviroment structure defining the box
*/
__global__ void generatePoints(curandState *globalState, Atom *atom, Environment *enviro);

/**
  This is a wrapper function for the calcEnergy kernel.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @return - the total energy of the system.
*/
double calcEnergyWrapper(Atom *atoms, Environment enviro);

/**
  Calculates the energy between n atoms where n is the
  the number of threads in the block. The block's sum is then stored 
  in the energySum array at its block index position.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @param *energySum - the array of block sums
*/
__global__ void calcEnergy(Atom *atoms, Environment enviro, double *energySum);


#endif //METROCUDAUTIL_H
