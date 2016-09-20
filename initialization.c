#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"
#include "user_parameters.c"
#include "initialization.h"

void initialize_all( struct outer_field *f){

  int i, j, k;

  // Reading in the user specified outer field parameters
  initialize_outer_parameters( f );

  /* OUTER FIELD INITIALIZATION */
  initializing_outer_field( f );


  /* HANDLING BACTERIA FIELD MEMORY ALLOCATION */
  for ( i = 0 ; i < f->dims[1] ; i++){
  for ( j = 0 ; j < f->dims[2] ; j++){
  for ( k = 0 ; k < f->dims[3] ; k++){

    // Reading in user specified bacteria field parameters
    initialize_bacteria_parameters( &(f->b[i][j][k]) );
    initializing_bacteria_field( &(f->b[i][j][k]) );


  }
  }
  }


  // Tests
  printf("Random concentration point: %f\n",f->b[1][1][1].conc[0][2][2]);

  
  return;
}




void initializing_bacteria_field( struct bacteria_field *b ){

    // Initializing the concentration field 
    allocate_3d( &(b->conc) , b->dims[0], b->dims[1], b->dims[2]);
    initialize_bacteria_conc( b );

    // Reading in user-specified bacteria field diffusion constants
    allocate_1d( &(b->D), b->dims[0]);
    initialize_bacteria_diffusion_constants( b );


  return;
}



void initializing_outer_field( struct outer_field *f ){

  int i, j;

  // Allocating memory for the tracer concentration field
  allocate_4d( &(f->conc), f->dims[0], f->dims[1], f->dims[2], f->dims[3]);

  // Allocating memory for the population field
  allocate_3d( &(f->pop),  f->dims[1], f->dims[2], f->dims[3]);


  // Allocating memory for the flow field
  allocate_3d( &(f->u),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->v),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->w),  f->dims[1], f->dims[2], f->dims[3]);


  // Reading in user specified outer field diffusion constants
  f->D = (double*) malloc( sizeof(double) * f->dims[0]);
  initialize_outer_diffusion_constants( f );

  











  // Allocating the bacteria field structures
  f->b = (struct bacteria_field ***) malloc( sizeof(struct bacteria_field **) * f->dims[1]);

  for( i = 0 ; i < f->dims[1] ; i++ ){
    f->b[i] = (struct bacteria_field **) malloc( sizeof(struct bacteria_field *) * f->dims[2]);
    for ( j = 0 ; j < f->dims[2] ; j++ ){
      f->b[i][j] = (struct bacteria_field *) malloc( sizeof(struct bacteria_field ) * f->dims[3]);
    }
  }


  return;
}





void allocate_1d( double **f, int d0 ){

  int i;


  f[0] = (double*) malloc( sizeof(double) * d0);

  return;
}



void allocate_2d( double ***f, int d0, int d1 ){

  int i, j;

  f[0] = (double**) malloc( sizeof(double*) * d0);

  for ( i = 0 ; i < d0; i++){
    f[0][i] = (double*) malloc( sizeof(double) * d1);
  }

  return;
}





void allocate_3d( double ****f, int d0, int d1, int d2 ){

  int i, j;

  f[0] = (double***) malloc( sizeof(double**) * d0);

  for ( i = 0 ; i < d0; i++){
    f[0][i] = (double**) malloc( sizeof(double*) * d1);
    for ( j = 0 ; j < d1; j++){
      f[0][i][j] = (double*) malloc( sizeof(double) * d2);
    }
  }

  return;
}



void allocate_4d( double *****f, int d0, int d1, int d2, int d3 ){

  int i, j, k;

  f[0] = (double****) malloc( sizeof(double***) * d0);

  for ( i = 0 ; i < d0; i++){
    f[0][i] = (double***) malloc( sizeof(double**) * d1);
    for ( j = 0 ; j < d1; j++){
      f[0][i][j] = (double**) malloc( sizeof(double*) * d2);
      for ( k = 0 ; k < d2; k++){
        f[0][i][j][k] = (double*) malloc( sizeof(double) * d3);
      }
    }
  }

  return;
}


