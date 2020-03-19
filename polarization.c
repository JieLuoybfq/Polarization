#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "scatterer.h"
#include "adda.h"
#include "orientation.h"
#include "constants.h"
#include "spectrum.h"

/* Following functions are used to calculate C_pol. */

// extracts absorption cross_section from the file written in the output ADDA format
double get_Cabs(char const * filename) {
	FILE * file = fopen(filename, "r");
	char curstr[STR_SIZE];
	double Cabs = 0.0;
	while(strcmp(curstr, "Cabs") != 0) {
		fscanf(file, "%s", curstr);
	}
	fscanf(file, "%s %le", curstr, &Cabs);
	fclose(file);
	return Cabs;
}

// get differential polarization cross-section from a certain directory
double differential_polarization_dir(char const * dir) {
	
	double Cabs_par = 0.0;
	double Cabs_per = 0.0;
	char file_x[STR_SIZE];
	char file_y[STR_SIZE];
	strcpy(file_x, dir);
	strcpy(file_y, dir);
	strcat(file_x, "/CrossSec-X");
	strcat(file_y, "/CrossSec-Y");

// if CrossSec-X file was not created the scatterer is symmetric from this direction
// and the differential cross-section equals 0
	if( access( file_x, R_OK ) == -1 ) {
		return 0;
	}
	
	Cabs_par = get_Cabs(file_x);
	Cabs_per = get_Cabs(file_y);
	return (Cabs_par - Cabs_per);
}

// get differential polarization cross-section of a certain scatterer using ADDA
double differential_polarization(Adda const * ad, Scatterer const * sc, Spectrum_point const * sp) {
	
	char * dirname = adda_run(ad, sc, sp);
	double dCpol = differential_polarization_dir(dirname);
	free(dirname);
	return dCpol;
}

// roots and weights for gauss integration
const double xg[6][6] = {
			{0.0          ,  0.0         ,  0.0         , 0.0         , 0.0         , 0.0         },
			{-0.5773502692,  0.5773502692,  0.0         , 0.0         , 0.0         , 0.0         },
			{-0.7745966692,  0.0         ,  0.7745966692, 0.0         , 0.0         , 0.0         },
			{-0.8611363115, -0.3399810436,  0.3399810436, 0.8611363115, 0.0         , 0.0         },
			{-0.9061798459, -0.5384693101,  0.0         , 0.5384693101, 0.9061798459, 0.0         }, 
			{-0.9324695142, -0.6612093864, -0.2386191861, 0.2386191861, 0.6612093864, 0.9324695142}};
			
const double wg[6][6] = {
			{ 2.0         ,  0.0         ,  0.0         , 0.0         , 0.0         , 0.0         },
			{ 1.0         ,  1.0         ,  0.0         , 0.0         , 0.0         , 0.0         },
			{ 0.5555555556,  0.8888888888,  0.5555555556, 0.0         , 0.0         , 0.0         },
			{ 0.3478548451,  0.6521451549,  0.6521451549, 0.3478548451, 0.0         , 0.0         },
			{ 0.2369268851,  0.4786286705,  0.5688888888, 0.4786286705, 0.2369268851, 0.0         },
			{ 0.1713244924,  0.3607615730,  0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924}};
/*
double testf(double x, double y) {
	return pow(x, 9.0) + 3.0 * pow(y, 8.0) - 2.0 * pow(x, 3.0) * pow(y, 2.0) - y + 1.5;
}
*/
double testf(double x, double y) {
	return x * sin(y / 2.0);
}
// calculate polarization_cross_section of a scatterer, orientated by beta, gamma,
// in a spectrum_point sp, with an integration method intmet (gauss 1-6 "g[1-6]", rectangles otherwise)
// with precision: the intervals are divided into ntheta, nomega parts
double polarization_cross_section(
			Scatterer const * sc, 
			double beta, double gamma, 
			int ntheta, int nomega, 
			Adda * ad,
			Spectrum_point const * sp,
			char const * intmet) {
// beta = angle(B, max_momentum_axis)
// gamma = angle(B, hnu)
// ADDA works in the particle system
	Vector prop = vector_make_from_spherical(gamma, 0.0);
// angle(hnu, z) = gamma
	adda_set_prop(ad, &prop);
// after rotating by "system" the z axis = max_momentum_axis
	EulerOrientation system = euler_make(0.0, PI - sc->shape.max_momentum_axis.beta, 2.0 * PI - sc->shape.max_momentum_axis.alpha);
// after rotating by magnet after the system z axis = B
	EulerOrientation magnet = euler_make(0.0, beta, 0.0);

	double dtheta = 2.0 * PI / ntheta;
	double domega = 2.0 * PI / nomega;
	double result  = 0.0;

// integration over omega - around B and theta - around max_momentum_axis
	if(strlen(intmet) == 2 && intmet[0] == 'g' && intmet[1] > '0' && intmet[1] < '7') {
// gauss
		int order = intmet[1] - '0';
		--order;		
		for(int i = 0; i < nomega; ++i) {
			for(int k = 0; k <= order; ++k) {
				magnet.alpha = domega / 2.0 * (xg[order][k] + 2.0 * i + 1.0);
				Quaternion qmag = quater_make_from_euler(&magnet);
				for(int j = 0; j < ntheta; ++j) {
					for(int m = 0; m <= order; ++m) {
						system.alpha = dtheta / 2.0 * (xg[order][m] + 2.0 * j + 1.0);
						Quaternion qsys = quater_make_from_euler(&system);
						Quaternion qtot = quater_composition(&qsys, &qmag);
						EulerOrientation total = euler_make_from_quater(&qtot);
						adda_set_euler(ad, &total);
						double f = differential_polarization(ad, sc, sp);
//@						double f = testf(magnet.alpha, system.alpha);
//						printf("f = %lf\n", f);
						result += f * wg[order][k] * wg[order][m];
					}
				}
			}
		}
		return result * domega * dtheta / 16.0 / PI / PI;
	} else {
// rectangles
		for(int i = 0; i < nomega; ++i) {
			magnet.alpha = (i + 0.5) * domega;
			Quaternion qmag = quater_make_from_euler(&magnet);
			for(int j = 0; j < ntheta; ++j) {
				system.alpha = (j + 0.5) * dtheta;
				Quaternion qsys = quater_make_from_euler(&system);
				Quaternion qtot = quater_composition(&qsys, &qmag);
				EulerOrientation total = euler_make_from_quater(&qtot);
				adda_set_euler(ad, &total);
				double f = differential_polarization(ad, sc, sp);
//@				double f = testf(magnet.alpha, system.alpha);
				result += dtheta * domega * f;
			}
		}	
		return result / (4.0 * PI * PI);
	}
}

// calculate from a directory full of subdirectories adda_run_n
// from number "begin" to number "end" (included) with integration method intmet
double polarization_cross_section_dir(
		Scatterer const * sc, 
		double beta, double gamma, 
		int ntheta, int nomega, 
		char const * dir, 
		int begin, int end,
		char const * intmet) {

	char dirname[STR_SIZE];
	double result  = 0.0;
	double dtheta = 2.0 * PI / ntheta;
	double domega = 2.0 * PI / nomega;
	int cnt = begin;

	if(strlen(intmet) == 2 && intmet[0] == 'g' && intmet[1] > '0' && intmet[1] < '7') {
// gauss
		int order = intmet[1] - '0';
		if(end - begin + 1 != ntheta * nomega * order) {
			printf("Inconsistent directory numbers. %d directories in the range, %d expected for \
			integration.", end - begin + 1, ntheta * nomega * order);
			return 0.0;
		}
		--order;		
		for(int i = 0; i < nomega; ++i) {
			for(int k = 0; k <= order; ++k) {
				for(int j = 0; j < ntheta; ++j) {
					for(int m = 0; m <= order; ++m) {
						sprintf(dirname, "%s/adda_run_%d", dir, cnt);
						++cnt;
						result += differential_polarization_dir(dirname) * wg[order][k] * wg[order][m];
					}
				}
			}
		}
		return result * domega * dtheta / 16.0 / PI / PI;
	} else {
// rectangles
		if(end - begin + 1 != ntheta * nomega) {
			printf("Inconsistent directory numbers. %d directories in the range, %d expected for \
			integration.", end - begin + 1, ntheta * nomega);
			return 0.0;
		}
		for(int i = 0; i < nomega; ++i) {
			for(int j = 0; j < ntheta; ++j) {
				sprintf(dirname, "%s/adda_run_%d", dir, cnt);
				++cnt;
				result += dtheta * domega * differential_polarization_dir(dirname);
			}
		}	
		return result / (4.0 * PI * PI);
	}
}
