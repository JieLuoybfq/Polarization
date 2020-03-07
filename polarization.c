#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "scatterer.h"
#include "adda.h"
#include "orientation.h"
#include "constants.h"

/* Following functions are used to calculate C_pol. */

double get_Cabs(char const * filename) {
	FILE * file = fopen(filename, "r");
	char curstr[STR_SIZE];
	double Cabs = 0.0;
	while(strcmp(curstr, "Cabs") != 0) {
		fscanf(file, "%s", &curstr[0]);
	}
	fscanf(file, "%s %le", &curstr[0], &Cabs);
	fclose(file);
	return Cabs;
}

double differential_polarization(Adda const * ad, Scatterer const * sc) {
	
	char * dirname = adda_run(ad, sc);
	double Cabs_par = 0.0;
	double Cabs_per = 0.0;
	char file_x[STR_SIZE];
	char file_y[STR_SIZE];
	strcpy(file_x, dirname);
	strcpy(file_y, dirname);
	strcat(file_x, "/CrossSec-X");
	strcat(file_y, "/CrossSec-Y");
	Cabs_par = get_Cabs(file_x);
	Cabs_per = get_Cabs(file_y);
	free(dirname);
	return (Cabs_par - Cabs_per);
}

/* Здесь нужно одно через другое выразить*/
/* Directory version*/
double differential_polarization_dir(char const * dir, Scatterer const * sc) {
	
	char const * dirname = dir;
//	printf("dirname = %s\n", dirname);
	double Cabs_par = 0.0;
	double Cabs_per = 0.0;
	char file_x[STR_SIZE];
	char file_y[STR_SIZE];
//	printf("here\n");
	strcpy(file_x, dirname);
	strcpy(file_y, dirname);
	strcat(file_x, "/CrossSec-X");
	strcat(file_y, "/CrossSec-Y");
//	printf("here\n");
//	printf("file_x = %s\n", file_x);
//	printf("file_y = %s\n", file_y);
	if( access( file_x, R_OK ) != -1 ) {
		Cabs_par = get_Cabs(file_x);
		Cabs_per = get_Cabs(file_y);
	} else {
		return 0;
	}
/*	free(dirname);

*/
//	printf("here\n");
	return (Cabs_par - Cabs_per);
}

// Integration

typedef struct ScatPosition {
	Scatterer const * scat; 
	Adda * adda; 
	EulerOrientation * system;
	EulerOrientation * magnet;
} ScatPosition;

void sp_set(ScatPosition * sp, Scatterer * sc, Adda * ad, EulerOrientation * s, EulerOrientation * m) {
	sp->scat = sc;
	sp->adda = ad;
	sp->system = s;
	sp->magnet = m;
}

typedef double (*INTEGR_FUNC)(double, double, ScatPosition * sp);


double pol_integr_func(double theta, double omega, ScatPosition * sp) {
	sp->magnet->alpha = omega;
	Quaternion qmag = quater_make_from_euler(sp->magnet);
	sp->system->alpha = theta;
	Quaternion qsys = quater_make_from_euler(sp->system);
	Quaternion qtot = quater_composition(&qsys, &qmag);
	EulerOrientation total = euler_make_from_quater(&qtot);
	adda_set_euler(sp->adda, &total);
	return differential_polarization(sp->adda, sp->scat);
}

// Integration

double integr_rect_2d(double a, double b, size_t n, double c, double d, size_t m, INTEGR_FUNC f, ScatPosition * sp) {
	double res = 0.0;
	double dx = (b - a) / n;
	double dy = (d - c) / m;
	for(size_t i = 0; i < n; ++i) {
		for(size_t j = 0; j < m; ++j) {
//			printf("i = %d j = %d\n", i, j);
			res += f(a + i * dx, c + j * dy, sp); 
		}
	}
	return dx * dy * res;
}

double integr_trap_2d(double a, double b, size_t n, double c, double d, size_t m, INTEGR_FUNC f, ScatPosition * sp) {
	double res = 0.0;
	double dx = (b - a) / n;
	double dy = (d - c) / m;
	for(size_t i = 1; i < n; ++i) {
		for(size_t j = 1; j < m; ++j) {
			res += f(a + i * dx, c + j * dy, sp); 
		}
	}
	for(int i = 1; i < n; ++i) {
		res += 0.5 * (f(i * dx, c, sp) + f(i * dx, d, sp));
	}
	for(int j = 1; j < m; ++j) {
		res += 0.5 * (f(a, j * dy, sp) + f(b, j * dy, sp));
	}
	res += 0.25 * (f(a, c, sp) + f(a, d, sp) + f(b, c, sp) + f(b, d, sp));
	return res * dx * dy;
}

double polarization_cross_section(Scatterer const * sc, double beta, double gamma, int ntheta, int nomega, Adda * ad) {
//	printf("Calculating polarization cross-section.\n");
	Vector prop = vector_make_from_spherical(gamma, 0.0);
	adda_set_prop(ad, &prop);
//		printf("here\n");
//	adda_print_parameters(ad);
	EulerOrientation system = euler_make(0.0, PI - sc->shape.max_momentum_axis.beta, 2.0 * PI - sc->shape.max_momentum_axis.alpha);
	EulerOrientation magnet = euler_make(0.0, beta, 0.0);
/*	ScatPosition sp = {sc, ad, &system, &magnet};
	double result = integr_trap_2d(0.0, 2.0 * PI, ntheta, 0.0, 2.0 * PI, nomega, pol_integr_func, &sp);
*/	double dtheta = 2.0 * PI / ntheta;
	double domega = 2.0 * PI / nomega;
	double result  = 0.0;
	for(int i = 0; i < nomega; ++i) {
		magnet.alpha = i * domega;
		Quaternion qmag = quater_make_from_euler(&magnet);
		for(int j = 0; j < ntheta; ++j) {
//			printf("i = %d j = %d\n", i, j);
			system.alpha = j * dtheta;
			Quaternion qsys = quater_make_from_euler(&system);
			Quaternion qtot = quater_composition(&qsys, &qmag);
			EulerOrientation total = euler_make_from_quater(&qtot);
			adda_set_euler(ad, &total);
			result += dtheta * domega * differential_polarization(ad, sc);
		}
	}
	
	return result / (4.0 * PI * PI);
}

/*Directory version*/
double polarization_cross_section_dir(Scatterer const * sc, double beta, double gamma, int ntheta, int nomega, char const * dir, int * pcnt) {
//	printf("Calculating polarization cross-section.\n");
			char dirname[STR_SIZE];
//	Vector prop = vector_make_from_spherical(gamma, 0.0);
//	adda_set_prop(ad, &prop);
//		printf("here\n");
//	adda_print_parameters(ad);
//	EulerOrientation system = euler_make(0.0, PI - sc->shape.max_momentum_axis.beta, 2.0 * PI - sc->shape.max_momentum_axis.alpha);
//	EulerOrientation magnet = euler_make(0.0, beta, 0.0);
/*	ScatPosition sp = {sc, ad, &system, &magnet};
	double result = integr_trap_2d(0.0, 2.0 * PI, ntheta, 0.0, 2.0 * PI, nomega, pol_integr_func, &sp);
*/	double dtheta = 2.0 * PI / ntheta;
	double domega = 2.0 * PI / nomega;
	double result  = 0.0;
	for(int i = 0; i < nomega; ++i) {
//		magnet.alpha = i * domega;
//		Quaternion qmag = quater_make_from_euler(&magnet);
		for(int j = 0; j < ntheta; ++j) {
//			printf("i = %d j = %d\n", i, j);
//			system.alpha = j * dtheta;
//			Quaternion qsys = quater_make_from_euler(&system);
//			Quaternion qtot = quater_composition(&qsys, &qmag);
//			EulerOrientation total = euler_make_from_quater(&qtot);
//			adda_set_euler(ad, &total);
			sprintf(dirname, "%s/adda_run_%d", dir, *pcnt);
//			printf("%s\n", dirname);
			++(*pcnt);
			result += dtheta * domega * differential_polarization_dir(dirname, sc);
		}
	}
	
	return result / (4.0 * PI * PI);
}
