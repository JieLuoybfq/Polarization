#pragma once

#include "scatterer.h"
#include "orientation.h"
#include "vector.h"

typedef struct Adda {
	double lambda;
	char * run_path;
	char * dir;
	EulerOrientation euler;
	Vector prop;
} Adda;

void adda_set_euler(Adda *, EulerOrientation const *);

void adda_set_prop(Adda *, Vector const *);

void adda_set(Adda *, double, char const *);

void adda_set_from_file(Adda *, char const *);

void adda_delete(Adda *);

void adda_print_parameters(Adda const *);

char * adda_run(Adda const *, Scatterer const *);
