#ifndef VARIABLE_HEADER
#define VARIABLE_HEADER

#include<stdio.h>
#include<stdlib.h>


struct outer_field{

  // Dimensions of the concentration field
  double l[3], h[3];
  int dims[4];
  // dims is number of grid points per dimension. 
  // l is the length of the entire bacteria
  // h is the inter-point distances

  // Bacteria population
  double ***pop;

  // Cartesian flow vectors
  double ***u;
  double ***v;
  double ***w;

  // Tracer concentration
  double ****conc;

  // Newtonian viscosity
  double mu;

  // Diffusion constants
  double *D;

  // Bacteria fields
  struct bacteria_field ***b;
  // Dimension (x,y,z)

  // Time steps
  int internal_steps_per_outer_steps;

  double duration;


  // Scale units (in meters, seconds and mol)
  double length_scale, time_scale, mol_scale;

};


struct bacteria_field{

  // Dimensions of the concentration field
  double l[2], h[2];
  int dims[3];
  // dims is number of grid points per dimension. 
  // l is the length of the entire bacteria
  // h is the inter-point distances

  // Size of each time step
  double dt;

  // Concentration field
  double ***conc; 
  // Dimensions ( species, r, z )

  // Rate of change due to diffusion
  double ***d_conc_dt;
  // Dimensions ( species, r, z )

  // Rate of change due to interactions
  double ***i_conc_dt;
  // Dimensions ( species, r, z )

  // Diffusion constants
  double *D;

  // Arrays of premultiplication coefficients (diffusion)
  /* The cylindrical diffusion equation, under azimuthal symmetry has three terms:
     first order differential in radius, second order differential in radius,
     second order differential in elevation*/

  double *d1r_left, *d1r_right;
  // Finite centered difference 1st order r differential discretization has 2 terms

  double *d2r_left, *d2r_center, *d2r_right;
  // Finite centered difference 2nd order r differential discretization has 3 terms

  double *d2z_bot, *d2z_center, *d2z_top;
  // Finite centered difference 2nd order z differential discretization has 3 terms


  // Commonly used indices
  int i0, i1, i2, i3, i4, i5;

  // Some double variables (no need to define more memory)
  double d0, d1, d2, d3, d4, d5;

  // Scale units (in meters, seconds and mol)
  double length_scale, time_scale, mol_scale;

};

#endif



