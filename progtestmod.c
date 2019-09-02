#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shape.h"
#include "modify.h"
#include "constants.h"

// Program for generation of a modified shape with f(r) = mod
int main(int argc, char *argv[]) {
	char filename[STR_SIZE];
	for(int i = 1; i < argc - 1; i += 2) {
		if (strcmp(argv[i], "-file") == 0) {
			strcpy(filename, argv[i + 1]);
		} else {
			printf("%s\n", argv[i]);
			printf("Unrecognized option. Please, use -adda, -dir or -file.\n");
			return 0;
		} 
	}
	printf("filename = %s\n", filename);
	
	Shape my_shape;
	shape_set(&my_shape, filename);
//	printf("file000 = %s\n", my_shape.source_file);
	shape_print_parameters(&my_shape);
//	shape_print_dipoles(&my_shape);
	Modificator m;
//	printf("file000 = %s\n", my_shape.source_file);
//	modificator_set(&m, density, 0.1, 0.002, 0.2, &my_shape, 0.2, 10);
	modificator_set(&m, density, 0.06, 0.002, 0.1, &my_shape, 0.2, 10);
	modificator_print_parameters(&m);
	check_distribution(&my_shape, 10, &m);
//	printf("file000 = %s\n", my_shape.source_file);
	Shape mod_shape = get_modified_shape(&my_shape, density, &m);
	printf("done\n");
	shape_print_parameters(&mod_shape);
	check_distribution(&mod_shape, 10, &m);	
	shape_delete(&my_shape);
	shape_delete(&mod_shape);
	return 0;
}
