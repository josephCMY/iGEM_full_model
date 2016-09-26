#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"
#define PI 3.14159265358979323846

// INPUT OUTER PARAMETERS HERE
void initialize_outer_parameters( struct outer_field *f ){


  // Grid point dimensions of the outer flow field
  f->dims[0] 				= 1;
  f->dims[1] 				= 1;
  f->dims[2] 				= 10;
  f->dims[3] 				= 10;

  // Length dimensions of the outer flow field (in scale units)
  f->l[0]    				= 1;
  f->l[1]    				= 1;
  f->l[2]    				= 1;

  // Duration of simulation
  f->duration  				= 10;
  // Output interval
  f->outputInterval 			= 0.1;

  // Numer of internal steps per outer step
  f->internal_steps_per_outer_steps 	= 1;

  f->mu					= 0;

  // Important scale units
  f->length_scale	= 1.0;
  f->time_scale		= 1e-4;
  f->mol_scale		= 1.0;




  return;

}


// INPUT INNER PARAMETERS HERE
void initialize_bacteria_parameters( struct bacteria_field *b){

  // Grid point dimensions of the bacteria field
  b->dims[0] 				= 1;
  b->dims[1] 				= 5;
  b->dims[2] 				= 20;

  // Length dimensions of the bacteria (in scale units)
  b->l[0]		= 0.25;
  b->l[1]		= 1.5;


  // Time step size (Note that f->duration must be a multiple of b->dt)
  b->dt			= 0.01;

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

  // Time invariant flow field
  int i, j, k;
  for ( i = 0; i < f->dims[1]; i++){
  for ( j = 0; j < f->dims[2]; j++){
  for ( k = 0; k < f->dims[3]; k++){

    f->u[i][j][k] = 0.0;
    f->v[i][j][k] = 0.0;
    f->w[i][j][k] = 0.0;

  }
  }
  }

  return;
}


// INITIALIZE OUTER CONCENTRATION FIELD
void initialize_outer_conc( struct outer_field *f ){

  // Initial outer field concentrations
  int i, j, k,n;
  for ( n = 0; n < f->dims[0]; n++){
  for ( i = 0; i < f->dims[1]; i++){
  for ( j = 0; j < f->dims[2]; j++){
  for ( k = 0; k < f->dims[3]; k++){

    f->conc[n][i][j][k] = 1e-4;

  }
  }
  }
  }




  return;
}




// INITIALIZE OUTER CONCENTRATION FIELD
void initialize_outer_pop( struct outer_field *f ){

  int i, j, k;

  for ( i = 0 ; i < f->dims[1] ; i++){
  for ( j = 0 ; j < f->dims[2] ; j++){
  for ( k = 0 ; k < f->dims[3] ; k++){

    f->pop[i][j][k] = 0;

  }
  }
  }

  return;
}






// INITIALIZE BACTERIA CONCENTRATION FIELD
void initialize_bacteria_conc( struct bacteria_field *b ){

  int i,j;
  double exponentPow, t0=1e-1, localConc;

  for( i = 0 ; i < b->dims[1] ; i++) {
  for( j = 0 ; j < b->dims[2] ; j++) {

    exponentPow =  pow( (i+0.5)*b->h[0] , 2  ) + pow( (j-b->dims[2]/2)*b->h[1] , 2  );
    exponentPow /= - 4 * b->D[0] * t0;
    localConc = 1*exp(  exponentPow  );

    if (localConc > 1e-1){
      b->conc[0][i][j] = localConc ;
    }
    else b->conc[0][i][j] = 0;

   
  }
  }


  return;
}


// USER SPECIFIED BACTERIA INTERACTION (RATE OF CHANGE OF CONCENTRATION)
void user_bacteria_interactions( struct bacteria_field *b){

  int n, i, j;
  for (n = 0 ; n < b->dims[0] ; n++){
  for (i = 0 ; i < b->dims[1] ; i++){
  for (j = 0 ; j < b->dims[2] ; j++){

    b->i_conc_dt[n][i][j] = 0.0;

  }
  }
  }


  return;
}



void user_bacteria_wall_flux ( struct bacteria_field *b){

  int n, i, j;
  double flux;
  double fluxCoeff1[1] = {10};
  double fluxCoeff2[1] = {3};
  double outerCUnit = b->outM / pow( b->outL, 3);
  double innerCUnit = b->mol_scale/pow( b->length_scale,3);
  double outerConc, innerConc;
  double outerConcDecrease;

  for ( n = 0 ; n < b->dims[0]; n++ ){

    outerConc = *( b->conc_env[n]);
    outerConc *= outerCUnit;

    for ( i = 0 ; i < b->dims[1]; i++ ){

      innerConc = b->conc[n][i][0];
      innerConc *= innerCUnit;

      flux = (fluxCoeff1[n]*outerConc - fluxCoeff2[n]*innerConc)* b->baseCoeff[i];
      b->i_conc_dt[n][i][0] +=  (flux/innerCUnit)/b->vol[i];
      outerConcDecrease -= b->dt*(flux/outerCUnit)/b->out_vol;


      innerConc = b->conc[n][i][b->dims[2]-1];
      innerConc *= innerCUnit;

      flux = (fluxCoeff1[n]*outerConc - fluxCoeff2[n]*innerConc)* b->baseCoeff[i];
      b->i_conc_dt[n][i][b->dims[2]-1] +=  (flux/innerCUnit)/b->vol[i];
      outerConcDecrease -= b->dt*(flux/outerCUnit)/b->out_vol;
      
    }


    for ( j = 0 ; j < b->dims[2]; j++ ){

      innerConc = b->conc[n][b->dims[1]-1][j];
      innerConc *= innerCUnit;

      flux = (fluxCoeff1[n]*outerConc - fluxCoeff2[n]*innerConc)* b->ringCoeff[j];
      b->i_conc_dt[n][b->dims[1]-1][j] +=  (flux/innerCUnit)/b->vol[b->dims[1]-1];
      outerConcDecrease -= b->dt*(flux/outerCUnit)/b->out_vol;
      
    }

  b->conc_env[n][0] += outerConcDecrease;

  }





  return;
}


