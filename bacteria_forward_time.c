#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"
#include "bacteria_forward_time.h"


void bacteria_diffusion_reaction_euler( struct bacteria_field *b){

  int n, i , j;

  // Computing the time differentials.
  bacteria_diffusion( b );
  user_bacteria_interactions( b );
  user_bacteria_wall_flux( b );
  

  for ( n = 0 ; n < b->dims[0] ; n++){
  for ( i = 0 ; i < b->dims[1] ; i++){
  for ( j = 0 ; j < b->dims[2] ; j++){
    b->conc[n][i][j] += b->dt * b->d_conc_dt[n][i][j];
    b->conc[n][i][j] += b->dt * b->i_conc_dt[n][i][j];
  }
  }
  }


  return;
}



void bacteria_diffusion( struct bacteria_field *b ){

  int n, i, j;
  double fluxDensity, flux;

  // Zeroing the diffusion dt array
  for (n = 0 ; n < b->dims[0] ; n++){
  for (i = 0 ; i < b->dims[1] ; i++){
  for (j = 0 ; j < b->dims[2] ; j++){

    b->d_conc_dt[n][i][j] = 0.0;

  }
  }
  }
 


  // Determining the mass flux through surface of each volume element
  for (n = 0 ; n < b->dims[0] ; n++){

    // Flow through the curved surfaces
    for (i = 0 ; i < b->dims[1] - 1 ; i++){
    for (j = 0 ; j < b->dims[2]     ; j++){
      fluxDensity = (b->conc[n][i+1][j] - b->conc[n][i][j]) / b->h[0] ; 
      flux        = b->D[n] * fluxDensity * b->ringCoeff[i];
      b->d_conc_dt[n][i+1][j] -= flux;
      b->d_conc_dt[n][i  ][j] += flux;
    }
    }

    // Flow through flat surfaces of volume element
    for ( i = 0 ; i < b->dims[1]     ; i++){
    for ( j = 0 ; j < b->dims[2] - 1 ; j++){

      fluxDensity = (b->conc[n][i][j+1] - b->conc[n][i][j]) / b->h[1];
      flux        = b->D[n] * fluxDensity * b->baseCoeff[i];
      b->d_conc_dt[n][i][j+1] -= flux;
      b->d_conc_dt[n][i][j  ] += flux;
    }
    }


    // Dividing by volume element
    for ( i = 0 ; i < b->dims[1] ; i++){
    for ( j = 0 ; j < b->dims[2] ; j++){

      b->d_conc_dt[n][i][j] /= b->vol[i];

    }
    }   

  }


  return;
}

/*
  // Handling the off-boundary parts
  for ( n = 0 ; n < b->dims[0]   ; n++){
  for ( i = 1 ; i < b->dims[1]-1 ; i++){
  for ( j = 1 ; j < b->dims[2]-1 ; j++){

    b->d_conc_dt[n][i][j] = ( b->conc[n][i][j] - b->conc[n][i-1][j] - (b->conc[n][i+1][j] - b->conc[n][i][j]) ) /  pow(b->h[0],2);

    b->d_conc_dt[n][i][j] += ( b->conc[n][i][j] - b->conc[n][i][j-1] - (b->conc[n][i][j+1] - b->conc[n][i][j]) ) /  pow(b->h[1],2);

    b->d_conc_dt[n][i][j] *= b->D[n];

  }
  }
  }


  // Handling the off-corner parts

  for ( n = 0 ; n < b->dims[0] ; n++){
  for ( i = 1 ; i < b->dims[1] -1; i++){

    j = 0;
    b->d_conc_dt[n][i][0] = ( b->conc[n][i][j] - b->conc[n][i-1][j] - (b->conc[n][i+1][j] - b->conc[n][i][j]) ) /  pow(b->h[0],2);
    b->d_conc_dt[n][i][0] += ( - (b->conc[n][i][j+1] - b->conc[n][i][j]) ) /  pow(b->h[1],2);
    b->d_conc_dt[n][i][0] *= b->D[n];

    j = b->dims[2]-1;
    b->d_conc_dt[n][i][j] = ( b->conc[n][i][j] - b->conc[n][i-1][j] - (b->conc[n][i+1][j] - b->conc[n][i][j]) ) /  pow(b->h[0],2);
    b->d_conc_dt[n][i][j] += ( b->conc[n][i][j] - b->conc[n][i][j-1] ) /  pow(b->h[1],2);
    b->d_conc_dt[n][i][j] *= b->D[n];

  }
  }


  for ( n = 0 ; n < b->dims[0] ; n++){
  for ( j = 1 ; j < b->dims[2]-1 ; j++){

    i = 0;
    b->d_conc_dt[n][i][j] = (  - (b->conc[n][i+1][j] - b->conc[n][i][j]) ) /  pow(b->h[0],2);

    b->d_conc_dt[n][i][j] += ( b->conc[n][i][j] - b->conc[n][i][j-1] - (b->conc[n][i][j+1] - b->conc[n][i][j]) ) /  pow(b->h[1],2);

    b->d_conc_dt[n][i][j] *= b->D[n];

    i = b->dims[1]-1;

    b->d_conc_dt[n][i][j] = ( b->conc[n][i][j] - b->conc[n][i-1][j] ) /  pow(b->h[0],2);

    b->d_conc_dt[n][i][j] += ( b->conc[n][i][j] - b->conc[n][i][j-1] - (b->conc[n][i][j+1] - b->conc[n][i][j]) ) /  pow(b->h[1],2);

    b->d_conc_dt[n][i][j] *= b->D[n];

  }
  }



  // Handling corner points
  for ( n = 0 ; n < b->dims[0] ; n++) {
    b->d_conc_dt[n][0][0] = (  - (b->conc[n][1][0] - b->conc[n][0][0]) ) /  pow(b->h[0],2);
    b->d_conc_dt[n][0][0] += ( - (b->conc[n][0][1] - b->conc[n][0][0]) ) /  pow(b->h[1],2);

    b->d_conc_dt[n][b->dims[1]-1][b->dims[2]-1] =  (  b->conc[n][b->dims[1]-1][b->dims[2]-1] - b->conc[n][b->dims[1]-2][b->dims[2]-1] ) /  pow(b->h[0],2);
    b->d_conc_dt[n][b->dims[1]-1][b->dims[2]-1] += (  b->conc[n][b->dims[1]-1][b->dims[2]-1] - b->conc[n][b->dims[1]-1][b->dims[2]-2] ) /  pow(b->h[1],2);

    b->d_conc_dt[n][0][b->dims[2]-1] =  ( - ( b->conc[n][1][b->dims[2]-1] - b->conc[n][0][b->dims[2]-1] ) ) /  pow(b->h[0],2);
    b->d_conc_dt[n][0][b->dims[2]-1] += (     b->conc[n][0][b->dims[2]-1] - b->conc[n][0][b->dims[2]-2] ) /  pow(b->h[1],2);

    b->d_conc_dt[n][b->dims[1]-1][0] =  (     b->conc[n][b->dims[1]-1][0] - b->conc[n][b->dims[1]-2][0] ) /  pow(b->h[0],2);
    b->d_conc_dt[n][b->dims[1]-1][0] += ( - ( b->conc[n][b->dims[1]-1][1] - b->conc[n][b->dims[1]-1][0] ) ) /  pow(b->h[1],2);

  }

*/
