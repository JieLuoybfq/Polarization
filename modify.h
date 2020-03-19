#pragma once
#include "shape.h"

/* 
Function distr(double, Modificator const *) and 
struct Modificator describe the details of particle density distribution. 
They must be defined together. */

// Current distribution:
// homogeneous core with raduis r_c
// for r > C ro(r) = A / r + B
// A, B, f_m are set, r_c is calculated
typedef struct Modificator {
	double A;
	double B;
	double f_m;                      // remaining_dipoles / all_dipoles
	double a;                        // biggest dipole distance from the mass center
	double r_c;                      // radius of the homogenous core
	double C;                        // core density
	int n;                           // accuracy (number of layers)
	double * initial_shape_volumes;  // initial volumes of the layers (not necessarily 
	                                 // spherical because the particle can be asymmetrtic)
									 // but the initial particle is assumed homogeneous
} Modificator;

typedef double (*Density_func)(double, Modificator const *);

double density(double, Modificator const *);

void modificator_set(Modificator *, Density_func, double, double, double, Shape const *, double, int, double);

void modificator_delete(Modificator *);

void modificator_print_parameters(Modificator const *);

Shape get_modified_shape(Shape const *, Density_func, Modificator const *);

void check_distribution(Shape const *, int, Modificator const * m);
