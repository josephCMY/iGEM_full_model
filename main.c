#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"
#include "initialization.h"
#include "initialization.c"
#include "bacteria_forward_time.h"
#include "bacteria_forward_time.c"
#include "outer_forward_time.h"
#include "outer_forward_time.c"
#define PI 3.14159265358979323846


void save_data( struct outer_field f, int tt);



int main(){

  struct outer_field f;
  initialize_all( &f );

  // Computing inner and outer time parameter loops
  int t1, t2;

  int nT = ( (f.duration * f.time_scale) / (f.b[0][0][0].dt * f.b[0][0][0].time_scale) )+0.5;



  int tRatio = f.internal_steps_per_outer_steps;
  int nTOuter = nT/tRatio;


  f.dt = tRatio * f.b[0][0][0].dt * f.b[0][0][0].time_scale / f.time_scale;

  printf("nT : %d,%f\n\n", nT, (f.duration * f.time_scale) / (f.b[0][0][0].dt * f.b[0][0][0].time_scale) );

  int i,j,k;

  double outInterval1 = (f.outputInterval/ f.duration) * nTOuter;
  printf("outInterval : %f\n\n",outInterval1);
  int outInterval = outInterval1;

  printf("Output outer loop intervals: %d\n\n", outInterval);
  // Looping over the systems now


  for ( t1 = 0; t1 < nTOuter ; t1++){

    // Saving data
    if ( t1 % outInterval == 0 ) {
      save_data( f, t1*tRatio);

      printf("Saved data at outer step: %06d\n",t1);

      printf("Real time: %e seconds\n\n",t1*tRatio*f.b[0][0][0].dt*f.b[0][0][0].time_scale);

    }

    for (t2 = 0; t2 < tRatio ; t2++){
     
      for ( i = 0; i < f.dims[1] ; i++){
      for ( j = 0; j < f.dims[2] ; j++){
      for ( k = 0; k < f.dims[3] ; k++){

        bacteria_diffusion_reaction_euler( &(f.b[i][j][k]) );

      }
      }
      }
    }

    outer_euler_solver( &f );

  }
  




  save_data( f, nT);
  printf("Saved data at outer step: %06d\n",t1);
  printf("Real time: %e seconds\n\n",nT*f.b[0][0][0].dt*f.b[0][0][0].time_scale);

  return;
}
















void save_data( struct outer_field f, int tt){


  // Output file names will be output.txt

  char outname[100]; 
  int n, i0, j0, k0, i, j;

  sprintf( outname, "output_%06d.csv",tt);

  FILE *outfile = fopen( outname , "w");

  fprintf(outfile, "CONCENTRATION FIELDS IN BACTERIA\n\n");
  for ( i0= 0 ; i0 < f.dims[1] ; i0++){
  for ( j0= 0 ; j0 < f.dims[2] ; j0++){
  for ( k0= 0 ; k0 < f.dims[3] ; k0++){

    fprintf(outfile, "Bacteria at: %d, %d, %d\n", i0, j0, k0);

    for ( n = 0 ; n < f.b[0][0][0].dims[0] ; n++ ){
      fprintf(outfile, "Species index: %d\n", n);
  
      for ( i = 0 ; i < f.b[0][0][0].dims[1] ; i++ ){
        for ( j = 0 ; j < f.b[0][0][0].dims[2] ; j++ ){

          fprintf(outfile, "%5e", f.b[i0][j0][k0].conc[n][i][j]);
          if (j < f.b[0][0][0].dims[2]-1) fprintf(outfile,", ");
        }
        fprintf(outfile, "\n");
      }
    }
  
  }
  }
  }


  return;
}


