#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"


// INPUT OUTER PARAMETERS HERE
void initialize_outer_parameters( struct outer_field *f ){


  // Grid point dimensions of the outer flow field
  f->dims[0] 				= 1;
  f->dims[1] 				= 10;
  f->dims[2] 				= 10;
  f->dims[3] 				= 10;

  // Length dimensions of the outer flow field (in scale units)
  f->l[0]    				= 1;
  f->l[1]    				= 1;
  f->l[2]    				= 10;

  // Duration of simulation
  f->duration  				= 30.0;

  // Numer of internal steps per outer step
  f->internal_steps_per_outer_steps 	= 100;

  f->mu					= 0;

  // Important scale units
  f->length_scale	= 1.0;
  f->time_scale		= 1.0;
  f->mol_scale		= 1.0;

  return;

}


// INPUT INNER PARAMETERS HERE
void initialize_bacteria_parameters( struct bacteria_field *b){

  // Grid point dimensions of the bacteria field
  b->dims[0] 				= 1;
  b->dims[1] 				= 5;
  b->dims[2] 				= 30;

  // Length dimensions of the bacteria (in scale units)
  b->l[0]		= 0.25;
  b->l[1]		= 1.50;


  // Important scale units
  b->length_scale	= 1e-6;
  b->time_scale		= 1e-4;
  b->mol_scale		= 1e-15;

  return;

}


// INITIALIIZE DIFFUSION CONSTANTS
void initialize_outer_diffusion_constants( struct outer_field *f){

  // Outer field diffusion constant in terms of outer field scale units
  f->D[0] = 0.88e-9;

  return;
}


void initialize_bacteria_diffusion_constants( struct bacteria_field *b){

  // Bacteria field diffusion constants in terms of bacteria field scale units
  b->D[0] = 0.022; // 0.22*10^-9 m2/s

  return;
}




// INITIALIZE OUTER FIELD FLOW
void initialize_outer_flow( struct outer_field *f ){

    

  return;
}


// INITIALIZE OUTER CONCENTRATION FIELD
void initialize_outer_conc( struct outer_field *f ){


  return;
}




// INITIALIZE OUTER CONCENTRATION FIELD
void initialize_outer_pop( struct outer_field *f ){


  return;
}






// INITIALIZE BACTERIA CONCENTRATION FIELD
void initialize_bacteria_conc( struct bacteria_field *b ){

  int i,j;

  for( i = 0 ; i < b->dims[1] ; i++) {
  for( j = 0 ; j < b->dims[2] ; j++) {
    b->conc[0][i][j] = 1.0;
   
  }
  }

  return;
}





