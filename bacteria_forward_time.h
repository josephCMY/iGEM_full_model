#ifndef BACTERIA_FORWARD_HEADER
#define BACTERIA_FORWARD_HEADER

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "variables.h"


void bacteria_diffusion( struct bacteria_field  *b );
void bacteria_diffusion_reaction_euler( struct bacteria_field *b );

#endif
