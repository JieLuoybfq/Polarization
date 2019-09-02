#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "vector.h"

// Vector creation
void vector_set(Vector * v, double x, double y, double z) {
	v->x = x;
	v->y = y;
	v->z = z;
}

Vector vector_make(double x, double y, double z) {
	return (Vector){x, y, z};
}

Vector vector_make_from_spherical(double theta, double lambda) {
	return (Vector){cos(lambda) * sin(theta), sin(lambda) * sin(theta), cos(theta)};
}

void vector_copy(Vector * v, Vector const * u) {
	v->x = u->x;
	v->y = u->y;
	v->z = u->z;
}

// Vector functions

// Vector * scalar
void vector_multiply(Vector * v, double k) {
	v->x *= k;
	v->y *= k;
	v->z *= k;
}

// Vector * Vector
double vector_scalar_product(Vector const * v, Vector const * u) {
	return v->x * u->x + v->y * u->y + v->z * u->z;
}

double vector_length_sqr(Vector const * v) {
	return vector_scalar_product(v, v);
}

double vector_length(Vector const * v) {
	return sqrt(vector_length_sqr(v));
}

// +/-
Vector vector_opposite(Vector const * v) {
	return (Vector){-v->x, -v->y, -v->z};
}

void vector_add(Vector * v, Vector const * u) {
	v->x += u->x;
	v->y += u->y;
	v->z += u->z;
}

void vector_substract(Vector * v, Vector const * u) {
	v->x -= u->x;
	v->y -= u->y;
	v->z -= u->z;
}

Vector vector_sum(Vector const * v, Vector const * u) {
	return (Vector){v->x + u->x, v->y + u->y, v->z + u->z};
}

Vector vector_difference(Vector const * v, Vector const * u) {
	return (Vector){v->x - u->x, v->y - u->y, v->z - u->z};
}

double vector_distance(Vector const * v, Vector const * u) {
	Vector s = vector_difference(v, u);
	return vector_length(&s);
}

// Geometric center
Vector vector_average(Vector const * begin, Vector const * end) {
	ptrdiff_t n = end - begin;
	if(n < 1) {
		printf("No vectors to average!\n");
	}
	Vector c = {0.0, 0.0, 0.0};
	for(Vector const * v = begin; v != end; ++v) {
		vector_add(&c, v);
	}
	if(n < 1) {
		printf("No vectors to average!\n");
	} else {
		vector_multiply(&c, 1.0 / n);
	}
	return c;
}

// i/o
void vector_print(Vector const * v) {
	printf("Vector: x = %lf y = %lf z = %lf\n", v->x, v->y, v->z);
}

// Comporators

// Comparator for Vector structure to sort by their length. 
int vector_compare_length(void const * v, void const * u) {
	double rv = vector_length_sqr((Vector *)v);
	double ru = vector_length_sqr((Vector *)u);
	if(rv < ru) return -1;
	if(rv > ru) return 1;
	return 0;
}

// Comparator for Vector structure to sort by z, then y, then x coordinates - the default order for adda files. 
int vector_compare_zyx(void const * v, void const * u)  {
	Vector * fv = ((Vector *)v);
	Vector * fu = ((Vector *)u);
	if(fv->z < fu->z) return -1;
	if(fv->z > fu->z) return 1;
	if(fv->y < fu->y) return -1;
	if(fv->y > fu->y) return 1;
	if(fv->x < fu->x) return -1;
	if(fv->x > fu->x) return 1;
	return 0;
}
