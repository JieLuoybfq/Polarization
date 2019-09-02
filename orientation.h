#pragma once

typedef struct EulerOrientation {
	double alpha;
	double beta;
	double gamma;
} EulerOrientation;

typedef struct Quaternion {
	double w;
	double x;
	double y;
	double z;
} Quaternion;


// EulerOrientation creation 
void euler_set(EulerOrientation *, double, double, double);

EulerOrientation euler_make(double, double, double);

EulerOrientation euler_make_from_quater(Quaternion const *);

// EulerOrientation functions
void euler_print(EulerOrientation const *);

// Quaternion creation
void quater_set(Quaternion *, double, double, double, double);

Quaternion quater_make(double, double, double, double);
#include "vector.h"
Quaternion quater_make_from_pair(double, Vector const *);

Quaternion quater_make_from_euler(EulerOrientation const *);

// Quaternion functions
Quaternion quater_composition(Quaternion const *, Quaternion const *);

void quater_print(Quaternion const *);

