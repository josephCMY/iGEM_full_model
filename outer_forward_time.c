#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "user_parameters.h"
#include "outer_forward_time.h"

void outer_euler_solver(struct outer_field *f ){


  int n, i, j, k;

  // Computing rate of change due to advection and diffusion;
  outer_conc_diffusion_advection( f );
  outer_population_diffusion_advection( f );


  // Population growth rate (FUTURE)

  // Euler solver for tracer field
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->conc[n][i][j][k] += (f->a_conc_dt[n][i][j][k] + f->d_conc_dt[n][i][j][k])* f->dt;
  }
  }
  }
  }

  // Euler solver for bacteria field
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->pop[i][j][k] += (f->a_pop_dt[i][j][k] + f->d_pop_dt[i][j][k])* f->dt;
  }
  }
  }




  return;
}


void outer_conc_diffusion_advection(struct outer_field *f){

  int n, i, j, k;
  double fluxDensity, flux;


  // Zeroing rate of change array
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->d_conc_dt[n][i][j][k] = 0.0;
    f->a_conc_dt[n][i][j][k] = 0.0;
  }
  }
  }
  }


  // Diffusion and advection through yz element
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]-1 ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){

    // Diffusion
    fluxDensity = (f->conc[n][i+1][j][k] - f->conc[n][i][j][k]) / f->h[0];
    flux = fluxDensity * f->D[n] * f->yzEle;
    f->d_conc_dt[n][i+1][j][k] -= flux;
    f->d_conc_dt[n][i  ][j][k] += flux;

    // Advection
    flux = 0.5*(f->conc[n][i+1][j][k] + f->conc[n][i][j][k]) * f->yzEle;
    flux *= (f->u[i+1][j][k] + f->u[i][j][k])*0.5;
    f->a_conc_dt[n][i+1][j][k] += flux;
    f->a_conc_dt[n][i  ][j][k] -= flux;


  }
  }
  }
  }


  // Diffusion and advection through xz element
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]-1 ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    // Diffusion
    fluxDensity = (f->conc[n][i][j+1][k] - f->conc[n][i][j][k]) / f->h[1];
    flux = fluxDensity * f->D[n] * f->xzEle;
    f->d_conc_dt[n][i][j+1][k] -= flux;
    f->d_conc_dt[n][i][j  ][k] += flux;

    // Advection
    flux = 0.5*(f->conc[n][i][j+1][k] + f->conc[n][i][j][k]) * f->xzEle;
    flux *= (f->v[i][j+1][k] + f->v[i][j][k])*0.5;
    f->a_conc_dt[n][i][j+1][k] += flux;
    f->a_conc_dt[n][i][j  ][k] -= flux;

  }
  }
  }
  }


  // Diffusion through xy element
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]-1 ; k++){

    // Diffusion
    fluxDensity = (f->conc[n][i][j][k+1] - f->conc[n][i][j][k]) / f->h[2];
    flux = fluxDensity * f->D[n] * f->xyEle;
    f->d_conc_dt[n][i][j][k+1] -= flux;
    f->d_conc_dt[n][i][j][k  ] += flux;

    // Advection
    flux = 0.5*(f->conc[n][i][j][k+1] + f->conc[n][i][j][k]) * f->xyEle;
    flux *= (f->w[i][j][k+1] + f->w[i][j][k])*0.5;
    f->a_conc_dt[n][i][j][k+1] += flux;
    f->a_conc_dt[n][i][j][k  ] -= flux;


  }
  }
  }
  }


  // Dividing by volume element
  for ( n = 0 ; n < f->dims[0]   ; n++){
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->d_conc_dt[n][i][j][k] /= f->vol;
    f->a_conc_dt[n][i][j][k] /= f->vol;
  }
  }
  }
  }
  

  return;
}





void outer_population_diffusion_advection(struct outer_field *f){

  int i, j, k;
  double fluxDensity, flux;


  // Zeroing rate of change array
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->d_pop_dt[i][j][k] = 0.0;
    f->a_pop_dt[i][j][k] = 0.0;
  }
  }
  }



  // Diffusion and advection through yz element
  for ( i = 0 ; i < f->dims[1]-1 ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){

    // Diffusion
    fluxDensity = (f->pop[i+1][j][k] - f->pop[i][j][k]) / f->h[0];
    flux = fluxDensity * f->bac_D * f->yzEle;
    f->d_pop_dt[i+1][j][k] -= flux;
    f->d_pop_dt[i  ][j][k] += flux;

    // Advection
    flux = 0.5*(f->pop[i+1][j][k] + f->pop[i][j][k]) * f->yzEle;
    flux *= (f->u[i+1][j][k] + f->u[i][j][k])*0.5;
    f->a_pop_dt[i+1][j][k] += flux;
    f->a_pop_dt[i  ][j][k] -= flux;


  }
  }
  }



  // Diffusion and advection through xz element
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]-1 ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    // Diffusion
    fluxDensity = (f->pop[i][j+1][k] - f->pop[i][j][k]) / f->h[1];
    flux = fluxDensity * f->bac_D * f->xzEle;
    f->d_pop_dt[i][j+1][k] -= flux;
    f->d_pop_dt[i][j  ][k] += flux;

    // Advection
    flux = 0.5*(f->pop[i][j+1][k] + f->pop[i][j][k]) * f->xzEle;
    flux *= (f->v[i][j+1][k] + f->v[i][j][k])*0.5;
    f->a_pop_dt[i][j+1][k] += flux;
    f->a_pop_dt[i][j  ][k] -= flux;

  }
  }
  }



  // Diffusion through xy element
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]-1 ; k++){

    // Diffusion
    fluxDensity = (f->pop[i][j][k+1] - f->pop[i][j][k]) / f->h[2];
    flux = fluxDensity * f->bac_D * f->xyEle;
    f->d_pop_dt[i][j][k+1] -= flux;
    f->d_pop_dt[i][j][k  ] += flux;

    // Advection
    flux = 0.5*(f->pop[i][j][k+1] + f->pop[i][j][k]) * f->xyEle;
    flux *= (f->u[i][j][k+1] + f->u[i][j][k])*0.5;
    f->a_pop_dt[i][j][k+1] += flux;
    f->a_pop_dt[i][j][k  ] -= flux;


  }
  }
  }



  // Dividing by volume element
  for ( i = 0 ; i < f->dims[1]   ; i++){
  for ( j = 0 ; j < f->dims[2]   ; j++){
  for ( k = 0 ; k < f->dims[3]   ; k++){
    f->d_pop_dt[i][j][k] /= f->vol;
    f->a_pop_dt[i][j][k] /= f->vol;
  }
  }
  }

  

  return;
}

