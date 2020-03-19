#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shape.h"
#include "constants.h"
#include "scatterer.h"

// constructor
// shape is set from the adda file (.geom), size(r_eq) is set manually
void scat_set(Scatterer * sc, char const * shape_file, double r_eq) {
	shape_set(&sc->shape, shape_file);
	sc->r_eq = r_eq;
}

// constructor
// set from the init file
void scat_set_from_file(Scatterer * sc , char const * filename) {
	FILE * file = fopen(filename, "r");
    char line[STR_SIZE];
	printf("reading parameters of the scatterer from file %s\n", filename);
	while (fscanf(file, "%s", line) != EOF) {
		if (strcmp(line, "r_eq") == 0) {
			fscanf(file, "%s %lf", line, &sc->r_eq);
		} else if (strcmp(line, "source_file") == 0) {
			fscanf(file, "%s", line);
			fscanf(file, "%s", line);
			shape_set(&sc->shape, line);
		} 
	}
	fclose(file);
}

// destructor
void scat_delete(Scatterer * sc) {
	shape_delete(&sc->shape);
}

// print parameters
void scat_print(Scatterer const * sc) {
	printf("Scatterer:\n");
	printf("r_eq = %lf\n", sc->r_eq);
	printf("Used shape:\n");
	shape_print_parameters(&sc->shape);
}
