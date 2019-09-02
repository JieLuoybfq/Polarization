#pragma once

typedef struct Vector {
	double x;
	double y;
	double z;
} Vector;

// Vector creation
void vector_set(Vector * v, double x, double y, double z);
Vector vector_make(double x, double y, double z);
Vector vector_make_from_spherical(double theta, double lambda);
void vector_copy(Vector * v, Vector const * u);

// Vector functions

// Vector * scalar
void vector_multiply(Vector * v, double k);

// Vector * Vector
double vector_scalar_product(Vector const * v, Vector const * u);
double vector_length_sqr(Vector const * v);
double vector_length(Vector const * v);

// +/-
Vector vector_opposite(Vector const * v);
void vector_add(Vector * v, Vector const * u);
void vector_substract(Vector * v, Vector const * u);
Vector vector_sum(Vector const * v, Vector const * u);
Vector vector_difference(Vector const * v, Vector const * u);
double vector_distance(Vector const * v, Vector const * u);

// Geometric center
Vector vector_average(Vector const * begin, Vector const * end);

// i/o
void vector_print(Vector const * v);

// Comporators
typedef int (*Compare_func)(void const *, void const *);

// Comparator for Vector structure to sort by their length. 
int vector_compare_length(void const * v, void const * u);

// Comparator for Vector structure to sort by z, then y, then x coordinates - the default order for adda files. 
int vector_compare_zyx(void const * v, void const * u);
