zip : mainprog.c shape.c geometry.c constants.c adda.c modify.c shape.h geometry.h constants.h adda.h modify.h makefile init
		zip program.zip $^

constants.o : constants.c 
		gcc -c $<

vector.o : vector.c
		gcc -c $<
orientation.o : orientation.c orientation.h constants.h
		gcc -c $<
geometry.o : geometry.c vector.h orientation.h constants.h geometry.h
		gcc -c $<
progtestgeom.o : progtestgeom.c vector.h orientation.h geometry.h constants.h 
		gcc -c $<
testgeom : progtestgeom.o vector.o orientation.o geometry.o constants.o
		gcc $^ -lm -o testgeom

shape.o : shape.c shape.h constants.h vector.h orientation.h geometry.h
		gcc -c $<
scatterer.o : scatterer.c scatterer.h shape.h constants.h
		gcc -c $<
adda.o : adda.c adda.h scatterer.h constants.h
		gcc -c $<
polarization.o : polarization.c adda.h scatterer.h orientation.h constants.h
		gcc -c $<
progtestpol.o : progtestpol.c polarization.h constants.h adda.h scatterer.h
		gcc -c $<
testpol : progtestpol.o polarization.o geometry.o orientation.o vector.o constants.o shape.o adda.o scatterer.o
		gcc -std=c99 $^ -lm -o testpol
		
modify.o : modify.c modify.h constants.h shape.h geometry.h vector.h
		gcc -c $<
progtestmod.o : progtestmod.c shape.h modify.h constants.h
		gcc -c $<
testmod : progtestmod.o modify.o geometry.o orientation.o vector.o constants.o shape.o
		gcc -std=c99 $^ -lm -o testmod