#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "adda.h"
#include "scatterer.h"
#include "constants.h"

void adda_set_euler(Adda * ad, EulerOrientation const * orient) {
	ad->euler.alpha = orient->alpha;
	ad->euler.beta = orient->beta;
	ad->euler.gamma = orient->gamma;
}

void adda_set_prop(Adda * ad, Vector const * p) {
	ad->prop.x = p->x;
	ad->prop.y = p->y;
	ad->prop.z = p->z;
}

void adda_set(Adda * ad, double lam, char const * path) {		
	ad->lambda = lam;
	ad->run_path = (char *)malloc(STR_SIZE * sizeof(char));
	ad->dir = (char *)malloc(STR_SIZE * sizeof(char));
	strcpy(ad->run_path, path);
	strcpy(ad->dir, "");
	EulerOrientation or = {0.0, 0.0, 0.0};
	Vector prop = {0.0, 0.0, 1.0};
	adda_set_euler(ad, &or);
	adda_set_prop(ad, &prop);
}

void adda_set_dir(Adda * ad, char const * dirname) {
	int result = mkdir(dirname, 0777);
	if(result == 0) {
		char filename[128] = "";
		strcpy(filename, dirname);
		strcat(filename, "/adda_counter");
		FILE * fout = fopen(filename, "w");
		fprintf(fout, "0\n");
		fclose(fout);
	}
	strcpy(ad->dir, dirname);
}

void adda_set_from_file(Adda * ad, char const * filename){
	FILE * file = fopen(filename, "r");
	adda_set(ad, 2.0 * PI, "");
//	ad->run_path = (char *)malloc(STR_SIZE * sizeof(char));
//	ad->dir = (char *)malloc(STR_SIZE * sizeof(char));
	printf("basic done\n");
    char line[STR_SIZE], line2[STR_SIZE];
	while (fscanf(file, "%s", &line) != EOF) {
//		printf("s = %s\n", line);
		if (strcmp(line, "lambda") == 0) {
			fscanf(file, "%s %lf", &line, &ad->lambda);
		} else if (strcmp(line, "alpha") == 0) {
			fscanf(file, "%s %lf", &line, &ad->euler.alpha);
		} else if (strcmp(line, "beta") == 0) {
			fscanf(file, "%s %lf", &line, &ad->euler.beta);
		} else if (strcmp(line, "gamma") == 0) {
			fscanf(file, "%s %lf", &line, &ad->euler.gamma);
		} else if (strcmp(line, "run_path") == 0) {
			fscanf(file, "%s %s", &line, &line2);
			strcpy(ad->run_path, line2);
//			printf("rp = %s\n", ad->run_path);
		} else if (strcmp(line, "dir") == 0) {
//			printf("reading dir\n");
			fscanf(file, "%s %s", &line, &line2);
			strcpy(ad->dir, line2);
//			printf("dir = %s\n", ad->dir);
			adda_set_dir(ad, ad->dir);
		}  
	}
	fclose(file);
}

void adda_delete(Adda * ad) {
	free(ad->run_path);
	free(ad->dir);
}

void adda_print_parameters(Adda const * ad) {
	printf("Parameters ADDA:\n");
	printf("lambda = %lf\n", ad->lambda);
	printf("path to run ADDA = %s\n", ad->run_path);
	if(strcmp(ad->dir, ".") == 0)
		printf("result directory = ADDA default\n");
	else
		printf("result directory = %s\n", ad->dir);
	printf("Orientation of the scatterer:\n");
	printf("alpha = %lf\n", ad->euler.alpha);
	printf("beta = %lf\n", ad->euler.beta);
	printf("gamma = %lf\n", ad->euler.gamma);
	printf("Propagation direction of the beam:\n");
	printf("x = %lf\n", ad->prop.x);
	printf("y = %lf\n", ad->prop.y);
	printf("z = %lf\n", ad->prop.z);
	printf("\n");
}

char * adda_run(Adda const * ad, Scatterer const * sc) {
	char * result = (char *)malloc(STR_SIZE * sizeof(char));
	
	char m[64] = "";
	char eu[64] = "";
	char lam[64] = "";
	char r[64] = "";
	char prop[64] = "";
	sprintf(m, "%lf %lf", sc->m_re, sc->m_im);
	sprintf(lam, "%lf", ad->lambda);
	sprintf(eu, "%lf %lf %lf", 
	(ad->euler.alpha) * 180 / PI, 
	(ad->euler.beta) * 180 / PI, 
	(ad->euler.gamma) * 180 / PI);
	sprintf(r, "%lf", sc->r_eq);
	sprintf(prop, "%lf %lf %lf", 
	ad->prop.x, ad->prop.y, ad->prop.z); 
	
	char command[1024] = "./";
	strcat(command, ad->run_path);
	
	strcat(command, " -shape read ");
	strcat(command, sc->shape.source_file);
	
	strcat(command, " -m ");
	strcat(command, m);
	strcat(command, " -lambda ");
	strcat(command, lam);
	strcat(command, " -eq_rad ");
	strcat(command, r);
	
	strcat(command, " -orient ");
	strcat(command, eu);
	strcat(command, " -prop ");
	strcat(command, prop);
	
	if(strcmp(ad->dir, "") != 0) {
		char dir[64] = "";
		char filename[128] = "";
		int c = 0;
		strcpy(filename, ad->dir);
		strcat(filename, "/adda_counter");
		FILE * fin = fopen(filename, "r");
		fscanf(fin, "%d", &c);
		fclose(fin);
		char number[8] = "";
		sprintf(number, "%d", c);
		strcat(dir, ad->dir);
		strcat(dir, "/adda_run_");
		strcat(dir, number);
		FILE * fout = fopen(filename, "w");
		fprintf(fout, "%d\n", c + 1);
		fclose(fout);
		
		strcat(command, " -dir ");
		strcat(command, dir);
		strcpy(result, dir);
	}
	
	strcat(command, " > addalog\n");
	
//	printf("command = %s", command);
	system(command);
	if(strcmp(ad->dir, "") == 0) {
		FILE * addalog = fopen("addalog", "r");
		char dir[64] = "";
		for(; dir[0] != '\''; )
			fscanf(addalog, "%s", &dir);
		for(int i = 1; i < strlen(dir); ++i)
			dir[i - 1] = dir[i];
		dir[strlen(dir) - 2] = 0;
		fclose(addalog);
		strcpy(result, dir);
	}
	return result;	
}