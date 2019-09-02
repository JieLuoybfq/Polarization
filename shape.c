#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <math.h>

#include "shape.h"
#include "geometry.h"
#include "constants.h"

void shape_set_mass_center(Shape * sh) {
	Vector center = vector_average(sh->dipoles, sh->dipoles + sh->number);
	vector_copy(&sh->mass_center, &center);
}

double shape_momentum_z(Shape const * sh, EulerOrientation const * f, Vector const * offset) {
	double im = 0.0;
	for(Vector * v = sh->dipoles; v != sh->dipoles + sh->number; ++v) {
		Vector u = vector_difference(v, offset);
		Vector w = rotate_euler(&u, f);
		im += (w.x * w.x + w.y * w.y);
	}
	return im;
}

void shape_set_max_momentum_axis(Shape * sh, size_t nb) {
	EulerOrientation cur = {0.0, 0.0, 0.0};
	euler_set(&sh->max_momentum_axis, 0.0, 0.0, 0.0);
	double imax = shape_momentum_z(sh, &cur, &sh->mass_center);
	double db = PI / nb;
	for(size_t i = 1; i < nb; ++i) {
		cur.beta = db * i;
		double da = db * sin(cur.beta);
		int na = (int)(2.0 * PI / da);
		for(size_t j = 0; j < na; ++j) {
			cur.alpha = j * da;
			double im = shape_momentum_z(sh, &cur, &sh->mass_center);
			if(im > imax) {
				sh->max_momentum_axis = cur;
				imax = im;
			}
		}
	}
}

void shape_set(Shape * sh, char const * fname) {
	sh->source_file = (char *)malloc(STR_SIZE * sizeof(char));
	strcpy(sh->source_file, fname);
	sh->number = 0;
	sh->size = MIN_DIPOLE_NUMBER;
	sh->dipoles = (Vector *)malloc((sh->size) * sizeof(Vector));
	FILE * file = fopen(fname, "r");
    char line[STR_SIZE];
	while (fgets(line, sizeof(line), file)) {
		if(strlen(line) == 0 || line[0] == '#')
			continue;
		if(sh->size == sh->number) {
			sh->size *= 2;
			sh->dipoles = (Vector *)realloc(sh->dipoles, (sh->size) * sizeof(Vector));
		}
		sscanf(line, "%lf %lf %lf", &sh->dipoles[sh->number].x, 
		&sh->dipoles[sh->number].y, &sh->dipoles[sh->number].z);
		++(sh->number);
	}
	fclose(file);
	shape_set_mass_center(sh);
	shape_set_max_momentum_axis(sh, DIVISION_NUMBER);
	printf("A shape with %ld dipoles was created from %s file.\n", sh->number, sh->source_file);
}

Shape shape_copy(Shape const * sh) {
	Shape new = {NULL, sh->number, sh->size, sh->mass_center, sh->max_momentum_axis, NULL};
	new.source_file = (char*)malloc(STR_SIZE * sizeof(char));
	strcpy(new.source_file, sh->source_file);
	new.dipoles = (Vector *)malloc(new.size * sizeof(Vector));
	for(size_t i = 0; i < new.number; ++i)
		new.dipoles[i] = sh->dipoles[i];
	return new;
}

void shape_delete(Shape * sh) {
	free(sh->source_file);
	free(sh->dipoles);
	sh->number = 0;
	sh->size = MIN_DIPOLE_NUMBER;
	vector_set(&sh->mass_center, 0.0, 0.0, 0.0);
	euler_set(&sh->max_momentum_axis, 0.0, 0.0, 0.0);
}

void shape_print_parameters(Shape const * sh) {
	printf("Parameters of shape read from file %s.\n", sh->source_file);
	printf("Mass center: x = %lf y = %lf z = %lf\n", sh->mass_center.x, sh->mass_center.y, sh->mass_center.z);
	printf("Max momentum axis orientation(degrees): alpha = %lf beta = %lf gamma = %lf\n", 
	sh->max_momentum_axis.alpha * 180 / PI, 
	sh->max_momentum_axis.beta * 180 / PI,
	sh->max_momentum_axis.gamma * 180 / PI);
}

void shape_print_dipoles(Shape const * sh) {
	printf("%ld dipoles of shape read from file %s.\n", sh->number, sh->source_file);
	for(size_t i = 0; i < sh->number; ++i) {
		printf("%lf %lf %lf\n", sh->dipoles[i].x, sh->dipoles[i].y, sh->dipoles[i].z);
	}
}

void shape_move(Shape * sh, Vector const * r) {
	for(size_t i = 0; i < sh->number; ++i) {
		vector_set(sh->dipoles + i, sh->dipoles[i].x + r->x, 
		sh->dipoles[i].y + r->y, sh->dipoles[i].z + r->z);
	}
}

double shape_max_center_dist(Shape const * sh) {
	double len2 = 0.0;
	for(int i = 0; i < sh->number; ++i) {
		Vector dif = vector_difference(sh->dipoles + i, &sh->mass_center);
		double cur = vector_length_sqr(&dif);
		if(len2 < cur)
			len2 = cur;
	}
	return sqrt(len2);
}
