pol : progpol.o spectrum.o polarization.o geometry.o orientation.o vector.o constants.o shape.o adda.o scatterer.o
		gcc -std=c99 $^ -lm -o pol
testpol : progtestpol.o polarization.o geometry.o orientation.o vector.o constants.o shape.o adda.o scatterer.o
		gcc -std=c99 $^ -lm -o testpol
testpoldir : progtestpoldir.o polarization.o geometry.o orientation.o vector.o constants.o shape.o adda.o scatterer.o
		gcc -std=c99 $^ -lm -o testpoldir

constants.o : constants.c 
		gcc -std=c99 -c $<

vector.o : vector.c
		gcc -std=c99 -c $<
orientation.o : orientation.c orientation.h constants.h
		gcc -std=c99 -c $<
geometry.o : geometry.c vector.h orientation.h constants.h geometry.h
		gcc -std=c99 -c $<
progtestgeom.o : progtestgeom.c vector.h orientation.h geometry.h constants.h 
		gcc -std=c99 -c $<
testgeom : progtestgeom.o vector.o orientation.o geometry.o constants.o
		gcc -std=c99 $^ -lm -o testgeom

shape.o : shape.c shape.h constants.h vector.h orientation.h geometry.h
		gcc -std=c99 -c $<
scatterer.o : scatterer.c scatterer.h shape.h constants.h
		gcc -std=c99 -c $<
adda.o : adda.c adda.h scatterer.h constants.h spectrum.h
		gcc -std=c99 -c $<
polarization.o : polarization.c adda.h scatterer.h orientation.h constants.h
		gcc -std=c99 -c $<
spectrum.o : spectrum.c spectrum.h 
		gcc -std=c99 -c $<
progpol.o : progpol.c polarization.h constants.h adda.h scatterer.h
		gcc -std=c99 -c $<
progtestpol.o : progtestpol.c polarization.h constants.h adda.h scatterer.h
		gcc -std=c99 -c $<
progtestpoldir.o : progtestpoldir.c polarization.h constants.h adda.h scatterer.h
		gcc -std=c99 -c $<
		
modify.o : modify.c modify.h constants.h shape.h geometry.h vector.h
		gcc -std=c99 -c $<
progtestmod.o : progtestmod.c shape.h modify.h constants.h
		gcc -std=c99 -c $<
testmod : progtestmod.o modify.o geometry.o orientation.o vector.o constants.o shape.o
		gcc -std=c99 $^ -lm -o testmod