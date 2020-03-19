#pragma once

#include "scatterer.h"
#include "orientation.h"
#include "vector.h"
#include "spectrum.h"

typedef struct Adda {
	char * run_path;
	char * dir;
	EulerOrientation euler;
	Vector prop;
} Adda;

// default constructor
void adda_set(Adda *);

// set specific fields
void adda_set_euler(Adda *, EulerOrientation const *);
void adda_set_prop(Adda *, Vector const *);
void adda_set_dir(Adda *, char const *);

// constructor from the init file
void adda_set_from_file(Adda *, char const *);

// destructor
void adda_delete(Adda *);

// print metainformation int the standard output
void adda_print_parameters(Adda const *);

// run ADDA with parameters set into the adda structure on the specific scatterer. 
// Returns the name of a directory to which ADDA saved the calculated files
char * adda_run(Adda const *, Scatterer const *, Spectrum_point const *);
