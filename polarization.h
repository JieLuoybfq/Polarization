#pragma once

#include "adda.h"
#include "spectrum.h"

double polarization_cross_section(Scatterer const *, double, double, 
int, int, Adda *, Spectrum_point const *, char const *);
double polarization_cross_section_dir(Scatterer const *, double, double, 
int, int, char const *, int, int, const char *);

