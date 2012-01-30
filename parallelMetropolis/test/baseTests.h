#ifndef BASETEST_H
#define BASETEST_H

#include <math.h>
#include <sys/time.h>
#include "../src/metroParallelUtil.h"

double make_periodic(double x, double box);

double wrap_into_box(double x, double box);

double calculate_energy(double **coords, int n_atoms, double *box_size, double sigma, double epsilon);

/**
  Calculates the energy assuming that all atoms are the same element.
  @param atoms - array of atoms that will be used to calculate energy
  @param atomNumber - the number of atoms in the system.
  @param boxSize- the length of one side of the box
*/

double calculate_energy(Atom *atoms, Environment enviro);

//retrieves the difference in time
long timevaldiff(struct timeval *starttime, struct timeval *finishtime);

#endif //BASETEST_H
