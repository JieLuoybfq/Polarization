#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shape.h"
#include "constants.h"
#include "scatterer.h"


void scat_set(Scatterer * sc, char const * shape_file, double m_re, double m_im, double r_eq) {
	shape_set(&sc->shape, shape_file);
	sc->m_re = m_re;
	sc->m_im = m_im;
	sc->r_eq = r_eq;
}


void scat_set_from_file(Scatterer * sc , char const * filename) {
	FILE * file = fopen(filename, "r");
    char line[STR_SIZE];
	int x_set = 0;
	double lambda = 0.0;
	printf("reading scat parameters from file %s\n", filename);
	while (fscanf(file, "%s", &line[0]) != EOF) {
		if(strcmp(line, "m_re") == 0) {
			fscanf(file, "%s %lf", &line[0], &sc->m_re);
		} else if (strcmp(line, "m_im") == 0) {
			fscanf(file, "%s %lf", &line[0], &sc->m_im);
		} else if (strcmp(line, "x") == 0) {
			fscanf(file, "%s %lf", &line[0], &sc->r_eq);
			x_set = 1;
		} else if (strcmp(line, "lambda") == 0) {
			fscanf(file, "%s %lf", &line[0], &lambda);
		} else if (strcmp(line, "r_eq") == 0) {
			fscanf(file, "%s %lf", &line[0], &sc->r_eq);
		} else if (strcmp(line, "source_file") == 0) {
			fscanf(file, "%s", &line[0]);
			fscanf(file, "%s", &line[0]);
			shape_set(&sc->shape, line);
		} 
	}
	if(x_set == 1) {
		if(lambda  < EPS) {
			printf("Lambda is not provided. It is set to 2 * PI.\n");
		} else {
			sc->r_eq *= 2.0 * PI / lambda;
		}
	}
	fclose(file);
}

void scat_delete(Scatterer * sc) {
//	printf("deleting\n");
	shape_delete(&sc->shape);
}

void scat_print(Scatterer const * sc) {
	printf("Scatterer:\n");
	printf("m = %lf + %lf * i\n", sc->m_re, sc->m_im);
	printf("r_eq = %lf\n", sc->r_eq);
	printf("Used shape:\n");
	shape_print_parameters(&sc->shape);
}
