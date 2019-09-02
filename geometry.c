#include <math.h>

#include "geometry.h"
#include "vector.h"
#include "orientation.h"
#include "constants.h"

// Rotations
Vector rotate(Vector const * v, double const m[3][3]) {
	double base[3] = {v->x, v->y, v->z};
	double res[3] = {0.0, 0.0, 0.0};
	for(size_t i = 0; i < 3; ++i)
		for(size_t j = 0; j < 3; ++j)
			res[i] += m[i][j] * base[j];
	return (Vector){res[0], res[1], res[2]};
}

Vector rotate_euler(Vector const * v, EulerOrientation const * f) {
	double r[3][3] = {
					  {cos(f->alpha) * cos(f->beta) * cos(f->gamma) - sin(f->alpha) * sin(f->gamma), 
					  -cos(f->alpha) * cos(f->beta) * sin(f->gamma) - sin(f->alpha) * cos(f->gamma), 
					   cos(f->alpha) * sin(f->beta)}, 
					  {sin(f->alpha) * cos(f->beta) * cos(f->gamma) + cos(f->alpha) * sin(f->gamma), 
					  -sin(f->alpha) * cos(f->beta) * sin(f->gamma) + cos(f->alpha) * cos(f->gamma), 
					   sin(f->alpha) * sin(f->beta)}, 
					  {-sin(f->beta) * cos(f->gamma),  
					   sin(f->beta) * sin(f->gamma), 
					   cos(f->beta)}
					 };
	return rotate(v, r);
}

void rotate_set_euler(Vector * v, EulerOrientation const * f) {
	Vector u = rotate_euler(v, f);
	vector_copy(v, &u);
}

Vector rotate_quater(Vector const * v, Quaternion const * q){
	double r[3][3] = {
					  {q->w * q->w + q->x * q->x - q->y * q->y - q->z * q->z, 
					   2.0 * (q->x * q->y - q->w * q->z), 
					   2.0 * (q->w * q->y + q->x * q->z)}, 
					  {2.0 * (q->x * q->y + q->w * q->z),  
					   q->w * q->w - q->x * q->x + q->y * q->y - q->z * q->z,  
					   2.0 * (q->y * q->z - q->w * q->x)}, 
					  {2.0 * (-q->w * q->y + q->x * q->z),  
					   2.0 * (q->y * q->z + q->w * q->x), 
					   q->w * q->w - q->x * q->x - q->y * q->y + q->z * q->z}
					 };
	return rotate(v, r);
}

void rotate_set_quater(Vector * v, Quaternion const * q) {
	Vector u = rotate_quater(v, q);
	vector_copy(v, &u);
}

// Volumes

double volume_spherical_layer(double r1, double r2) {
	return 4.0 * PI / 3.0 * fabs(r2 * r2 * r2 - r1 * r1 * r1);
}
