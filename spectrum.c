#include <stdlib.h>
#include <stdio.h>

#include "spectrum.h"
#include "constants.h"

// constructor
void set(Spectrum_point * sp, double lam, double m_re, double m_im) {
	sp->lambda = lam;
	sp->m_re = m_re;
	sp->m_im = m_im;
}

// Read a list of n of (lambda, m) pairs from a file.
// if there are less than n pairs available (accounting for the fact we save every kth)
// the rest is filled with (2*PI, 1 + 0*i)
Spectrum_point * spec_read_from_file(
		char const * filename, 
		double const lam0,                    // wavelength range 
		int n, int k) {                       // maximum number of pairs, save every kth of them
	Spectrum_point * sp = (Spectrum_point *)malloc(n * sizeof(Spectrum_point));
	for(int i = 0; i < n; ++i) {
			sp[i] = (Spectrum_point){2.0 * PI, 1.0, 0.0};
	}
	FILE * file = fopen(filename, "r");
	int j = 0, p = k;
	char line[STR_SIZE];
	fgets(line, sizeof(line), file);
	int begin = 0;
	for(int i = 0; j < n &&
	fscanf(file, "%le %le %le", &sp[j].lambda, &sp[j].m_re, &sp[j].m_im) != EOF; ++i) {
		if(!begin && sp[j].lambda >= lam0) {
			begin = 1;
		}
		if(begin) {
			if(p == k) {
				++j;
				p = 1;
			} else {
				++p;
			}
		}
	}
	if(j < n) {
		printf("Incomplete wavelength range. Only %d points read from file %s.\n", j, filename);
	} 
	fclose(file);
	return sp;
}
