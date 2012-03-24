#ifndef BASETEST_H
#define BASETEST_H

#include <math.h>

double make_periodic(double x, double box);

double wrap_into_box(double x, double box);

double calculate_energy(double **coords, int n_atoms, double *box_size, double sigma, double epsilon);

#endif //BASETEST_H
