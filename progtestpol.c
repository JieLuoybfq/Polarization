#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "polarization.h"
#include "constants.h"
#include "adda.h"
#include "scatterer.h"


//Program for comparison with formulas for small prolate spheroid from prolate_ellipsoidz.geom

int main(int argc, char *argv[]) {
	char initfile[STR_SIZE];
	for(int i = 1; i < argc - 1; i += 2) {
		if (strcmp(argv[i], "-init") == 0) {
			strcpy(initfile, argv[i + 1]);
		} else {
			printf("%s\n", argv[i]);
			printf("Unrecognized option. Please, use -adda, -dir or -file.\n");
			return 0;
		} 
	}
	printf("initfile = %s\n", initfile);
	
	Scatterer my_scat;
	scat_set_from_file(&my_scat, initfile);
	scat_print(&my_scat);
	printf("creating adda\n");
	Adda my_adda;
	adda_set_from_file(&my_adda, initfile);
	adda_print_parameters(&my_adda);
	int n = 3;
	int m = 3;
	double dbeta = PI / 2.0 / n;
	double dgamma = PI / 2.0 / m;
//	taken from Crosssec-X and Crosssec-Y after the run 
//	adda -shape read shapes/prolate_ellipsoidz.geom -m 1.7 0.1 -lambda 100 -eq_rad 0.1 -prop 0 1 0 
//	double Cx = 5.082981788e-05;
//	double Cy = 2.879236556e-05;
//	taken from Crosssec-X and Crosssec-Y after the run 
//	adda -shape ellipsoid 1.0 4.0 -m 1.7 0.1 -lambda 100 -eq_rad 0.1 -prop 1 0 0 -save_geom test_prolate.geom
	double Cx = 6.873160629e-05;
	double Cy = 2.592325308e-05;
//	printf("we are here\n");
	printf("beta \\ gamma ");
	for(int j = 0; j <= m; ++j)
		printf(" %31f ", j * dgamma);
	printf("\n");
//	printf("here\n");
	for(int i = 0; i <= n; ++i) {
//		printf("here\n");
		printf("beta = %6f ", i * dbeta);
//	printf("here\n");
		for(int j = 0; j <= m; ++j) {
			double Cpol = polarization_cross_section(&my_scat, dbeta * i, dgamma * j, 10, 10, &my_adda);
			printf(" % 8.7e / % 8.7e ", Cpol, (Cy - Cx) / 2.0 * sin(dgamma * j) * sin(dgamma * j) 
			* 1.5 * (cos(dbeta * i) * cos(dbeta * i) - 1.0 / 3.0));
		}
		printf("\n");
	}
	scat_delete(&my_scat);
	adda_delete(&my_adda);
	return 0;
}
