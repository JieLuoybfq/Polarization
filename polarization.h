#pragma once

#include "adda.h"
double polarization_cross_section(Scatterer const *, double, double, int, int, Adda *);
double polarization_cross_section_dir(Scatterer const *, double, double, int, int, char const *, int *);

