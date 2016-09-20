#ifndef INITIALIZATION_HEADER
#define INITIALIZATION_HEADER

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"variables.h"


  // GENERIC INITIALIZATION FUNCTIONS
  void initialize_all( struct outer_field *f );
  void allocate_4d( double *****f, int d0, int d1, int d2, int d3 );
  void allocate_3d( double  ****f, int d0, int d1, int d2         );
  void allocate_2d( double   ***f, int d0, int d1                 );
  void allocate_1d( double    **f, int d0                         );

  void initializing_bacteria_field( struct bacteria_field *b );
  void initializing_outer_field( struct outer_field *f );

#endif
