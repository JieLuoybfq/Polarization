#include <stdio.h>
#include <math.h>

#include "vector.h"
#include "constants.h"
#include "orientation.h"

// EulerOrientation creation 
void euler_set(EulerOrientation * eo, double a, double b, double c) {
	eo->alpha = a;
	eo->beta = b;
	eo->gamma = c;
}

EulerOrientation euler_make(double a, double b, double c) {
	return (EulerOrientation){a, b, c};
}

EulerOrientation euler_make_from_quater(Quaternion const * q) {
	double cosa = 1.0;
	double cosb = (q->w * q->w + q->z * q->z - q->x * q->x - q->y * q->y) / (q->w * q->w + q->z * q->z + q->x * q->x + q->y * q->y);
	if(cosb < -1.0) cosb = -1.0;
	if(cosb > 1.0) cosb = 1.0;
	double beta = acos(cosb);
	double cosg = 1.0;
	if(fabs(sin(beta)) > EPS) {
		cosg = -2.0 * (q->z * q->x - q->y * q->w) / (q->w * q->w + q->z * q->z + q->x * q->x + q->y * q->y) / sin(beta);
		cosa = 2.0 * (q->y * q->w + q->z * q->x) / (q->w * q->w + q->z * q->z + q->x * q->x + q->y * q->y) / sin(beta);
	} else {
		cosa = (q->w * q->w - q->x * q->x + q->y * q->y - q->z * q->z) / (q->w * q->w + q->z * q->z + q->x * q->x + q->y * q->y);
	}
	if(cosa < -1.0) cosa = -1.0;
	if(cosg < -1.0) cosg = -1.0;
	if(cosa > 1.0) cosa = 1.0;
	if(cosg > 1.0) cosg = 1.0;
	return euler_make(acos(cosa), acos(cosb), acos(cosg));
}

// EulerOrientation functions
void euler_print(EulerOrientation const * eo) {
	printf("Euler angels: alpha = %lf beta = %lf gamma = %lf\n", 
	eo->alpha, eo->beta, eo->gamma);
}

// Quaternion creation
void quater_set(Quaternion * q, double w, double x, double y, double z) {
	q->w = w;
	q->x = x;
	q->y = y;
	q->z = z;
}

Quaternion quater_make(double w, double x, double y, double z) {
	return (Quaternion) {w, x, y, z};
}

Quaternion quater_make_from_pair(double theta, Vector const * v) {
	return (Quaternion) {cos(theta / 2.0), 
	                     sin(theta / 2.0) * v->x, 
	                     sin(theta / 2.0) * v->y, 
						 sin(theta / 2.0) * v->z};
}

Quaternion quater_make_from_euler(EulerOrientation const * eo) {
	Vector pz = {0.0, 0.0, 1.0};
	Vector py = {0.0, 1.0, 0.0};
	Quaternion qa = quater_make_from_pair(eo->alpha, &pz);
	Quaternion qb = quater_make_from_pair(eo->beta, &py);
	Quaternion qg = quater_make_from_pair(eo->gamma, &pz);
	Quaternion q = quater_composition(&qa, &qb);
	return quater_composition(&q, &qg);
}

// Quaternion functions
Quaternion quater_composition(Quaternion const * q1, Quaternion const * q2) {
	return quater_make(q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z,
	                   q1->w * q2->x + q2->w * q1->x + q1->y * q2->z - q2->y * q1->z,
	                   q1->w * q2->y + q2->w * q1->y + q1->z * q2->x - q2->z * q1->x, 
	                   q1->w * q2->z + q2->w * q1->z + q1->x * q2->y - q1->y * q2->x);
}

void quater_print(Quaternion const * q) {
	printf("Quaternion: w = %lf x = %lf y = %lf z = %lf\n", q->w, q->x, q->y, q->z);
}
