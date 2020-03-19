#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <math.h>

#include "modify.h"
#include "shape.h"
#include "vector.h"
#include "geometry.h"
#include "constants.h"

// (supposed) volume density in dipoles / volume (adda units)
double density(double r, Modificator const * m) {
	if(r < m->r_c)
		return m->C;
	if(r > m->a)
		return 0.0;
	return m->A / r + m->B;
}

// number of remaining dipoles in layer from r1 to r2
double fraction(Density_func df, double r1, double r2, Modificator const * m, Shape const * sh, int k) {
// volume of the layer inside a particle * required density of particles = 
// supposed number of dipoles to remain
	return m->initial_shape_volumes[k] * df((r1 + r2) / 2.0, m); 
}

void modificator_set(
		Modificator * m, 
		Density_func df, 
		double A, double B, 
		double f_m, 
		Shape const * sh, 
		double x,                 // core radius / full radius (maximum axis if not spherical)
		int n, 
		double dipole_volume) {
// set what we know
	m->A = A;
	m->B = B;
	m->f_m = f_m;
	m->a = shape_max_center_dist(sh);
	m->r_c = x * m->a;
	m->n = n;
	
	m->C = 0.0;
	m->initial_shape_volumes = (double*)malloc((m->n) * sizeof(double));
	int remain = sh->number * f_m;
	int shell = 0;
	int core = 0;
	int d[n + 1];
	for(int i = 0; i <= n; ++i)
		d[i] = 0;
	double da = m->a / n;
// count dipoles by layers
	for(int i = 0; i < sh->number; ++i) {
		double len = vector_distance(sh->dipoles + i, &sh->mass_center);
		++d[(int)(len / m->a * m->n)];
	}
	d[n - 1] += d[n];
// here d[i] = number of initial dipoles in the ith layer
// set initial_shape_volumes
	for(int i = 0; i < n; ++i)
		m->initial_shape_volumes[i] = dipole_volume * d[i];
// calculate shell = number of dipoles supposed to remain outside of the core
	int i = n - 1;
	for(i = n - 1; i > 0 && i * da >= m->r_c; --i) {
		int need = fraction(df, i * da, (i + 1) * da, m, sh, i);
		if(d[i] < need) {
			printf("we need %d dipoles in %d layer, but there are only %d\n", need, i, d[i]);
		}
		shell += (int)fmin(d[i], need);
	}
// if shell > remain we can't make the shape with required density distribution and 
// f_m out of the initial shape just by removing dipoles
	if(remain <= shell) {
		printf("Too much materialin the shell. The core is empty.\n");
		printf("Can't achieve the required factor %lf. Least possible number of dipoles is %d.\n \
Least possible factor is %lf.\n", m->f_m,shell, (double)shell/sh->number);
		shell = remain;
	}
	for(int j = i; j >= 0; --j)
		core += d[j];
// everything in remain = f_m * initial_dipoles - shell must be in the core
// if initially there were less material than remain - shell,
// that we can't just remove dipoles
	if(core < remain - shell) {
		printf("Too much material for the core.\n");
	} else {
		core = remain - shell;
	}
// set the core density
	m->C = (double)core / volume_spherical_layer(0, m->r_c);
}

// destructor
void modificator_delete(Modificator * m) {
	free(m->initial_shape_volumes);
}

// print
void modificator_print_parameters(Modificator const * m) {
	printf("Modificator:\n");
	printf("A = %lf\n", m->A);
	printf("B = %lf\n", m->B);
	printf("a = %lf\n", m->a);
	printf("r_c = %lf\n", m->r_c);
	printf("f_m = %lf\n", m->f_m);
	printf("C = %lf\n", m->C);
}

// remove some dipoles if the shape's so that the remaining dipoles are distributed according 
// to the df with m denoting parameters of the distribution
Shape get_modified_shape(Shape const * base_shape, Density_func df, Modificator const * m) {

	Shape new_shape = shape_copy(base_shape);
	Shape * sh = &new_shape;

// sort the dipoles by their distance to the center
	Vector rc = vector_opposite(&sh->mass_center);
	shape_move(sh, &rc);
	qsort(sh->dipoles, sh->number, sizeof(Vector), vector_compare_length);
	double dr = vector_length(sh->dipoles + sh->number - 1) / m->n;
	int j = 0, k = 0;
	srand (time(NULL));
	for(size_t i = 0; i < m->n; ++i) {
		
//		interval [j..k) of sh->dipoles contains all dipoles 
//		with r from [i * dr, (i + 1) * dr) (one spherical layer)
		j = k;
		for(; k < sh->number && 
		vector_length_sqr(sh->dipoles + k) < (i + 1) * (i + 1) * dr * dr; ++k);
		if(i == m->n - 1)
			k = sh->number;
		int remain = (int)(fraction(df, i * dr, (i + 1) * dr, m, sh, i));
// to remove
		int dp = (int)fmax((double)(k - j) - remain, 0.0);
		printf("layer = %ld %d here %d remains %d to delete\n", i, k - j, remain, dp);
		for(int p = 0; p < dp; ) {
			int q = rand() % (k - j);
			if(sh->dipoles[j + q].x != MAX_COORD) {
				sh->dipoles[j + q].x = sh->dipoles[j + q].y 
				= sh->dipoles[j + q].z = MAX_COORD;
				++p;
			}
		}
	}
	qsort(sh->dipoles, sh->number, sizeof(Vector), vector_compare_length);
	for(; sh->number > 0 && sh->dipoles[sh->number - 1].x == MAX_COORD; --sh->number);
	if(sh->number < 4) {
		printf("Warning! The remaining dipole set is too small. Only %ld dipoles.\n", sh->number);
	}
	shape_move(sh, &sh->mass_center);
	qsort(sh->dipoles, sh->number, sizeof(Vector), vector_compare_zyx);

//  Modify name of the resulting file.
	char * place = sh->source_file + strlen(sh->source_file);
	for(; place != sh->source_file && *place != '.'; --place);
	strcpy(place, "_mod.geom");

	FILE * base_file = fopen(base_shape->source_file, "r");
	FILE * new_file = fopen(sh->source_file, "w");
	
//  Copy the ADDA information from the original file.
	char line[STR_SIZE];
	while (fgets(line, sizeof(line), base_file) && strlen(line) > 0 && line[0] == '#') {
		fputs(line, new_file);
	}	
//  Print dipoles' coordinates.
	for(size_t i = 0; i < sh->number; ++i) {
		fprintf(new_file, "%d %d %d\n", (int)sh->dipoles[i].x, (int)sh->dipoles[i].y, (int)sh->dipoles[i].z);
	}
	
	fclose(base_file);
	fclose(new_file);
	
	return new_shape;
}

// prints the distribution of the shape's dipoles over r
void check_distribution(Shape const * sh, int n, Modificator const * m) {
	int * d = (int *)malloc((n + 1) * sizeof(int));
	double da = m->a / n;
// count the dipoles by layers
	for(int i = 0; i <= n; ++i)
		d[i] = 0;
	for(Vector * p = sh->dipoles; p != sh->dipoles + sh->number; ++p) {
		++d[(int)floor(vector_distance(p, &sh->mass_center) / m->a * n)];
	}
	d[n - 1] += d[n];

	printf("Shape contains %ld dipoles. Distribution:\n", sh->number);
	for(int i = 0; i < n - 1; ++i) {
		printf("r : [%f, %f) %d dipoles or %f of total, density = %f, expected %f\n", 
		i * da, da * (i + 1), d[i], (double)d[i] / sh->number, 
		(double)d[i] / m->initial_shape_volumes[i], density(da * (i + 0.5), m));
	}
	printf("r : [%f, %f] %d dipoles or %f of total, density = %f, expected %f\n", 
	da * (n - 1), m->a, d[n - 1], (double)d[n - 1] / sh->number, 
	(double)d[n - 1] / m->initial_shape_volumes[n - 1], density(da * (n - 0.5), m));
	free(d);
}
