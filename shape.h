#pragma once

#include "vector.h"
#include "orientation.h"

// Represents a set of dipoles that constitute a particle.
typedef struct Shape {
	
	Vector * dipoles;                  
	size_t number;
	size_t size;                     // size of the dipoles array. can be bigger than number. In that case
	                                 // the remaining points are junk.
	Vector mass_center;
	EulerOrientation max_momentum_axis;   // Orientation of the axis relative to which the momentum of inertia is maximum.
	char * source_file;              
	
} Shape;

// constructor from .geom file
void shape_set(Shape *, char const *);

// constructor from .geom file
Shape shape_copy(Shape const *);

// destructor
void shape_delete(Shape *);

// Printing 
// metainformation
void shape_print_parameters(Shape const *);
// dipoles coordinates (like the source .geom without meta)
void shape_print_dipoles(Shape const *); 

// move shape by vector
void shape_move(Shape *, Vector const *);

// calculate the spacial parameter of the shape (in ADDA units)
double shape_max_center_dist(Shape const *);