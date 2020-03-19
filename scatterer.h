#pragma once

#include "shape.h"

typedef struct Scatterer {
	
	Shape shape;
	double r_eq;
	
} Scatterer;

// constructor
// shape is set from the adda file (.geom), size(r_eq) is set manually
void scat_set(Scatterer *, char const *, double);

// constructor
// set from the init file
void scat_set_from_file(Scatterer *, char const *);

// destructor
void scat_delete(Scatterer *);

// print parameters
void scat_print(Scatterer const *);
