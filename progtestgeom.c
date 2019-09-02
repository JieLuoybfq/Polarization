#include <stdio.h>

#include "vector.h"
#include "orientation.h"
#include "geometry.h"
#include "constants.h"

int main() {
	Vector p = {1.0, 2.0, 3.0};
/*	Point rp = point_opposite(&p);
	point_print(&rp);
	point_assign(&p, 2.0, -3.5, 11.0);
	point_print(&p);
	Point prp = point_move(&p, &rp);
	point_print(&prp);
*/	
	EulerOrientation eo = euler_make(PI / 23.0, PI / 13.0, PI / 16.0);
	euler_print(&eo);
	Quaternion q = quater_make_from_euler(&eo);
	quater_print(&q);
	EulerOrientation eo2 = euler_make_from_quater(&q);
	euler_print(&eo2);
/*	euler_assign(&eo2, -PI / 16.0, -PI / 13.0, -PI / 23.0);
	Quaternion q2 = quater_make_from_euler(&eo2);
	quater_print(&q2);
	Quaternion q3 = quater_composition(&q, &q2);
	quater_print(&q3);
	EulerOrientation eo3 = euler_make_from_quater(&q3);
	euler_print(&eo3);
*/	Vector p2 = rotate_euler(&p, &eo);
	vector_print(&p2);
	Vector p3 = rotate_quater(&p, &q);
	vector_print(&p3);
	return 0;
}