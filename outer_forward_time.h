#ifndef OUTER_FORWARD_HEADER
#define OUTER_FORWARD_HEADER

#include<stdio.h>
#include<stdlib.h>
#include "variables.h"


void outer_euler_solver(struct outer_field *f );
void outer_conc_diffusion_advection(struct outer_field *f);
void outer_population_diffusion_advection(struct outer_field *f);

#endif


