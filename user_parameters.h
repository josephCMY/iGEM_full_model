#ifndef USER_PARAMETERS_HEADER
#define USER_PARAMETERS_HEADER

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"variables.h"


  // GENERIC INITIALIZATION FUNCTIONS
  void initialize_outer_parameters( struct outer_field *f );
  void initialize_bacteria_parameters( struct bacteria_field *b );


  // INITIALIZE OUTER DIFFUSION-FLOW FIELD QUANTITIES
  void initialize_outer_flow( struct outer_field *f );
  void initialize_outer_conc( struct outer_field *f );
  void initialize_outer_pop(  struct outer_field *f );

  // INITIALIZE DIFFUSION CONSTANTS
  void initialize_outer_diffusion_constants( struct outer_field *f);
  void initialize_bacteria_diffusion_constants( struct bacteria_field *b);


  void initialize_bacteria_conc( struct bacteria_field *b );
  void user_bacteria_interactions( struct bacteria_field *b);

#endif
