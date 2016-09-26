#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"
#include "user_parameters.c"
#include "initialization.h"
#define PI 3.14159265358979323846

void initialize_all( struct outer_field *f){

  int i, j, k, n;

  // Reading in the user specified outer field parameters
  printf("Initializing outer parameters\n");
  initialize_outer_parameters( f );


  /* OUTER FIELD INITIALIZATION */
  initializing_outer_field( f );
  f->bac_D = 4e-10 * f->time_scale / pow(f->length_scale,2) ; 
  // E. Coli in Motion by Howard C Berg, page 56, chapter 6: physical constraints
  printf("Initialized outer parameters\n\n");

  /* HANDLING BACTERIA FIELD MEMORY ALLOCATION */
  printf("Initializing bacteria parameters\n");
  for ( i = 0 ; i < f->dims[1] ; i++){
  for ( j = 0 ; j < f->dims[2] ; j++){
  for ( k = 0 ; k < f->dims[3] ; k++){

    // Reading in user specified bacteria scalar parameters
    initialize_bacteria_parameters( &(f->b[i][j][k]) );

    // Read in user specified information about initial conditions and allocate
    // memory to all fields
    initializing_bacteria_field( &(f->b[i][j][k]) );

    // Linking the environmental concentration to the bacteria
    for ( n = 0 ; n < f->dims[0] ; n++){
      f->b[i][j][k].conc_env[n] = &(f->conc[n][i][j][k]);
    }

    f->b[i][j][k].out_vol = f->vol;
    f->b[i][j][k].outM = f->mol_scale;
    f->b[i][j][k].outL = f->length_scale;



  }
  }
  }
  printf("Initialized bacteria parameters\n\n");

  // Tests
//  printf("Random bacteria concentration point: %f\n",f->b[1][1][1].conc[0][2][2]);
//  printf("Random outer    concentration point: %f\n",f->conc[0][2][2][1]);


  
  return;
}




void initializing_bacteria_field( struct bacteria_field *b ){

    int i;

    // Initializing the concentration field 
    allocate_3d( &(b->conc) , b->dims[0], b->dims[1], b->dims[2]);

    // Reading in user-specified bacteria field diffusion constants
    allocate_1d( &(b->D), b->dims[0]);
    initialize_bacteria_diffusion_constants( b );

    allocate_3d( &(b->i_conc_dt), b->dims[0], b->dims[1], b->dims[2]);
    allocate_3d( &(b->d_conc_dt), b->dims[0], b->dims[1], b->dims[2]);

    // Allocating memory for the environmental concentrations
    b->conc_env = (double**) malloc( sizeof(double*) * b->dims[0]);

    // Grid intervals
    for ( i = 0 ; i < 2 ; i++) b->h[i] = b->l[i] / b->dims[i+1] ;    

    // Ring coeff
    allocate_1d( &(b->ringCoeff), b->dims[1] - 1);
    for ( i = 0 ; i < b->dims[1]-1 ; i++){
      b->ringCoeff[i] = 2 * PI * ( (i+1)*b->h[0] ) * b->h[1];
    }

    // Base coeff
    allocate_1d( &(b->baseCoeff), b->dims[1]);
    for ( i = 0 ; i < b->dims[1] ; i++){
      b->baseCoeff[i] = PI * pow( (i+1)*b->h[0], 2 ) - PI * pow( (i)*b->h[0], 2 );
    }


    // Volume of each ring
    allocate_1d( &(b->vol), b->dims[1]);
    for ( i = 0 ; i < b->dims[1] ; i++){
      b->vol[i] = PI * pow( (i+1)*b->h[0], 2 ) - PI * pow( (i)*b->h[0], 2 );
      b->vol[i] *= b->h[1];
    }

    initialize_bacteria_conc( b );

  return;
}



void initializing_outer_field( struct outer_field *f ){

  int i, j;

  // Allocating memory for the tracer concentration field
  allocate_4d( &(f->conc), f->dims[0], f->dims[1], f->dims[2], f->dims[3]);
  allocate_4d( &(f->d_conc_dt), f->dims[0], f->dims[1], f->dims[2], f->dims[3]);
  allocate_4d( &(f->a_conc_dt), f->dims[0], f->dims[1], f->dims[2], f->dims[3]);

  // Allocating memory for the population field
  allocate_3d( &(f->pop),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->d_pop_dt),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->a_pop_dt),  f->dims[1], f->dims[2], f->dims[3]);

  // Allocating memory for the flow field
  allocate_3d( &(f->u),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->v),  f->dims[1], f->dims[2], f->dims[3]);
  allocate_3d( &(f->w),  f->dims[1], f->dims[2], f->dims[3]);


  // Reading in user specified outer field diffusion constants
  f->D = (double*) malloc( sizeof(double) * f->dims[0]);
  initialize_outer_diffusion_constants( f );


  // Compute grid spacing
  for ( i = 0 ; i < 3 ; i++) f->h[i] = f->l[i] / f->dims[i+1] ;


  // Employing user specified outer field initializations
  initialize_outer_flow( f );
  initialize_outer_conc( f );
  initialize_outer_pop(  f );

  
  // Allocating the bacteria field structures
  f->b = (struct bacteria_field ***) malloc( sizeof(struct bacteria_field **) * f->dims[1]);

  for( i = 0 ; i < f->dims[1] ; i++ ){
    f->b[i] = (struct bacteria_field **) malloc( sizeof(struct bacteria_field *) * f->dims[2]);
    for ( j = 0 ; j < f->dims[2] ; j++ ){
      f->b[i][j] = (struct bacteria_field *) malloc( sizeof(struct bacteria_field ) * f->dims[3]);
    }
  }


  // Handling the 3 types of surfaces for each cuboid volume element
  f->xyEle = f->h[0]*f->h[1];
  f->yzEle = f->h[1]*f->h[2];
  f->xzEle = f->h[0]*f->h[2];

  // Volume element
  f->vol = f->h[0]*f->h[1]*f->h[2];

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


