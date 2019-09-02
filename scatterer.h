#pragma once

#include "shape.h"

typedef struct Scatterer {
	
	Shape shape;
	double m_re;
	double m_im;
	double r_eq;
	
} Scatterer;

void scat_set(Scatterer *, char const *, double, double, double);

void scat_set_from_file(Scatterer *, char const *);

void scat_delete(Scatterer *);

void scat_print(Scatterer const *);
