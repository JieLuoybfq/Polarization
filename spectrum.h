#pragma once

// The pair (lambda, m) for which the calculation is performed
typedef struct Spectrum_point {
	double lambda;
	double m_re;
	double m_im;
} Spectrum_point;


// constructor
void set(Spectrum_point *, double, double, double);

// Read a list of (lambda, m) pairs from a file.  
Spectrum_point * spec_read_from_file(char const *, double, int, int);