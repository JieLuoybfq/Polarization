#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "polarization.h"
#include "constants.h"
#include "adda.h"
#include "scatterer.h"
#include "spectrum.h"


/* This program calculates the polarization cross_section of a given scatterer over an interval of 
light propagation angles [bet_left, bet_right] x [gam_left, gam_right] 
and wavelength [lam_left, lam_right] taken from file recieved through command line option -init
and saves the result into a file taken from command line option -result.
If the result file is not chosen it is set to default "result.dat" */

typedef struct Boundary {
	double left;
	double right; 
	int n;
} Boundary;

void calculate_polarization(Scatterer const * sc, Adda * ad,
		Boundary const * bet, 
		Boundary const * gam, 
		Spectrum_point const * sp, 
		int nlam,
		int tet_p,
		int omg_p,
		char const * resultfile,
		char const * integration_method) {
	
	double dbet = (bet->n == 1 ? 0.0 : (bet->right - bet->left) / (bet->n - 1));
	double dgam = (gam->n == 1 ? 0.0 : (gam->right - gam->left) / (gam->n - 1));

	FILE * file = fopen(resultfile, "w");	
	fprintf(file, "beta  %9s","");
	for(int i = 0; i < bet->n; ++i)
		for(int j = 0; j < gam->n; ++j)
			fprintf(file, "%+9.7le ", (bet->left + i * dbet) * 180.0 / PI);
	fprintf(file, "\n");

	fprintf(file, "gamma %9s","");
	for(int i = 0; i < bet->n; ++i)
		for(int j = 0; j < gam->n; ++j)
			fprintf(file, "%+9.7le ", (gam->left + j * dgam) * 180.0 / PI);
	fprintf(file, "\n");
	
	fprintf(file, "lambda\n");
	printf("Calculation starts.\n");
	for(int k = 0; k < nlam; ++k) {
		fprintf(file, "%9.7le  ", sp[k].lambda);
		for(int i = 0; i < bet->n; ++i) {
			for(int j = 0; j < gam->n; ++j) {
				printf("Calculating lambda = %lf bet = %lf gam = %lf ", 
				sp[k].lambda, bet->left + dbet * i, gam->left + dgam * j);
				double Cpol = polarization_cross_section(
					sc, 
					bet->left + dbet * i, 
					gam->left + dgam * j, 
					tet_p, 
					omg_p, 
					ad, 
					sp + k,
					integration_method);
				fprintf(file, "%+9.7le ", Cpol);
				printf(" Cpol = %+15.13le \n", Cpol);
			}
		}
		fprintf(file, "\n");
	}
	fclose(file);
	printf("Calculation ends.\n");
}

void print_help() {
	printf("Available options:\n");
	printf("  -init initfile,\n");
	printf("  -result resultfile,\n");
	printf("  -int integration_method (r for rectangles, gx for gauss, x from 1 to 6 - degree. (default r),\n");
	printf("  -sp specfile (with lambda m_re m-im),\n");
	printf("  -lam0 lam0 - starting wavelength if specfile is set (default lam0 = 0),\n");
	printf("  -n n - number of wavelength points if specfile is set (default n = 1),\n");
	printf("  -k k - use every kth point starting from lam0 (default k = 1),\n");
	printf("  -lam lam - wavelength for calculation (if specfile is not set),\n");
	printf("  -m m_re m_im - refraction index for calculation if specfile is not set,\n");
	printf("  -h - view list of options.\n");
}

int main(int argc, char *argv[]) {
	if(argc == 2 && strcmp(argv[1], "-h") == 0) {
		print_help();
		return 0;
	}
	
	char initfile[STR_SIZE];      // file with the ADDA information
	strcpy(initfile, "");
	char resultfile[STR_SIZE];    // file with the resulting table
	strcpy(resultfile, "result.dat");
	char intmet[STR_SIZE];        // integration method flag. -r for rectangles -gx for gauss, x = degree
	strcpy(intmet, "r");
	Spectrum_point singlesp = (Spectrum_point){-1.0, -1.0, -1.0}; // (lambda, m), used if fully set and spfile isn't set
	char specfile[STR_SIZE];
	strcpy(specfile, "");         // file with (lambda, m), used if set
	char ** endptr = 0;
	double lam0 = 0.0;
	int nlam = 1, klam = 1;
	for(int i = 1; i < argc - 1; ) {
		if (strcmp(argv[i], "-init") == 0) {
			strcpy(initfile, argv[i + 1]);
		} else if (strcmp(argv[i], "-result") == 0) {
			strcpy(resultfile, argv[i + 1]);
		} else if (strcmp(argv[i], "-int") == 0) {
			strcpy(intmet, argv[i + 1]);
		} else if (strcmp(argv[i], "-lam") == 0) {
			singlesp.lambda = strtod(argv[i + 1], endptr);
		} else if (strcmp(argv[i], "-lam0") == 0) {
			lam0 = strtod(argv[i + 1], endptr);
		} else if (strcmp(argv[i], "-n") == 0) {
			nlam = strtol(argv[i + 1], endptr, 10);
		} else if (strcmp(argv[i], "-k") == 0) {
			klam = strtol(argv[i + 1], endptr, 10);
		} else if (i < argc - 2 && strcmp(argv[i], "-m") == 0) {
			singlesp.m_re = strtod(argv[i + 1], endptr);
			singlesp.m_im = strtod(argv[i + 2], endptr);
			i += 3;
			continue;
		} else if (strcmp(argv[i], "-sp") == 0) {
			strcpy(specfile, argv[i + 1]);
		} else {
			printf("%s\n", argv[i]);
			printf("Unrecognized option. Please, use -h to view the list of options.\n");
			return 0;
		} 
		i += 2;
	}
	if (strcmp(initfile, "") == 0) {
		printf("Please, set the init file:");
		scanf("%s", initfile);
	}
	printf("initfile = %s\n", initfile);
	
// setting the scatterer and adda parameters from initfile	
	Scatterer my_scat;
	scat_set_from_file(&my_scat, initfile);
	scat_print(&my_scat);
	
	Adda my_adda;
	adda_set_from_file(&my_adda, initfile);
	adda_print_parameters(&my_adda);

	if((singlesp.lambda < 0 || singlesp.m_re < 0 || singlesp.m_im < 0) && strcmp(specfile, "") == 0) {
		printf("Please, set the spectral file:");
		scanf("%s", specfile);
	}
	if(strcmp(specfile, "") != 0) {
		printf("%d points will be read from the file %s starting from lambda = %lf. \
		k = %d.\n", nlam, specfile, lam0, klam);
	} else {
		printf("Single point calculation. Lambda = %lf, m = %lf + %lf * i\n", singlesp.lambda, 
		singlesp.m_re, singlesp.m_im);
		nlam = 1;
	}
	Spectrum_point * sp = &singlesp;
	if(strcmp(specfile, "") != 0) {
		sp = spec_read_from_file(specfile, lam0, nlam, klam);
	}
// setting the calculation boundaries and precision from initfile
	Boundary bet = {0.0, PI / 2, 2};
	Boundary gam = {0.0, PI / 2, 2};
	int tet_p = 10;
	int omg_p = 10;
	
	FILE * file = fopen(initfile, "r");
    char line[STR_SIZE];
	while (fscanf(file, "%s", line) != EOF) {
		if (strcmp(line, "beta") == 0) {
			fscanf(file, "%s %lf %lf %d", line, &bet.left, &bet.right, &bet.n);
			bet.left *= PI / 180.0;
			bet.right *= PI / 180.0;
		} else if (strcmp(line, "gamma") == 0) {
			fscanf(file, "%s %lf %lf %d", line, &gam.left, &gam.right, &gam.n);
			gam.left *= PI / 180.0;
			gam.right *= PI / 180.0;
		} if (strcmp(line, "precision") == 0) {
			fscanf(file, "%s %d %d", line, &tet_p, &omg_p);
		}
	}
	fclose(file);

// performing the calculation
	calculate_polarization(&my_scat, &my_adda, &bet, &gam, sp, nlam, tet_p, omg_p, resultfile, intmet);
	
	scat_delete(&my_scat);
	adda_delete(&my_adda);
	
	return 0;
}
