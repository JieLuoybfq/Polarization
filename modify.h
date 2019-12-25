#pragma once
#include "shape.h"

/* 
Function distr(double, Modificator const *) and 
struct Modificator describe the details of particle density distribution. 
They must be defined together. */

typedef struct Modificator {
	double A;
	double B;
	double f_m;  // remaining_dipoles / all_dipoles
	double a;    // biggest dipole distance from the mass center
	double r_c;  // radius of the homogenous core
	double C;
	int n;
	double * initial_shape_volumes;
} Modificator;

typedef double (*Density_func)(double, Modificator const *);

double density(double, Modificator const *);

void modificator_set(Modificator *, Density_func, double, double, double, Shape const *, double, int, double);

void modificator_delete(Modificator *);

void modificator_print_parameters(Modificator const *);

Shape get_modified_shape(Shape const *, Density_func, Modificator const *);

void check_distribution(Shape const *, int, Modificator const * m);
