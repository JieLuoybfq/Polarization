#pragma once

#include "vector.h"
#include "orientation.h"

// Rotations
Vector vector_rotate(Vector const * v, double const m[3][3]);

Vector rotate_euler(Vector const * v, EulerOrientation const * f);
void rotate_set_euler(Vector * v, EulerOrientation const * f);

Vector rotate_quater(Vector const * v, Quaternion const * q);
void rotate_set_quater(Vector * v, Quaternion const * q);

// Volumes
double volume_spherical_layer(double r1, double r2);
