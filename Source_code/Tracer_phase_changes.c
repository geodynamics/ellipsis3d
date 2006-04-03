/*

Copyright (C) 1995 The GeoFramework Consortium

This file is part of Ellipsis3D.

Ellipsis3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, version 2,
as published by the Free Software Foundation.

Ellipsis3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Author:
  Louis Moresi <louis.moresi@sci.monash.edu>

*/


#include <math.h>
/* #include <malloc.h> */
/* #include <rpc/rpc.h> */

#if (defined __sunos__)
#include <string.h>
#else
#if (!defined __GNUC__)
#include <strings.h> 
#endif
#endif
/*
#if (! defined __GNUC__)
#include <rpc/xdr.h> 
#endif
*/
#include "element_definitions.h"
#include "global_defs.h"

int get_eq_phase(
		 struct All_variables *E,
		 int tr
		 )
{
  int phase_no,found,phase_b;
  standard_precision Depth,Temp,DPhase;

  if(E->tracer.Phases[E->tracer.property_group[tr]] == 1) {
    E->tracer.phase_function[tr] = 1.0;
    return(0);
  }

  if(E->data.grav_acc == 0.0) {
    report(E,"Phase changes cannot be calculated because gravitational acceleration is zero");
    if(E->control.verbose)
      fprintf(stderr,"Phase changes cannot be calculated because gravitational acceleration is zero");
    return(0);
  }

  /* Search in phase below each phase boundary, keep track of
     phase if found and phase boundary number. If not found,
     must be in the highest pressure phase above the highest
     phase boundary */
  phase_no=0;
  phase_b=0;
  found=-1;

  do {
    Depth = E->tracer.PZ0[E->tracer.property_group[tr]*MAX_MATERIAL_PHASES+phase_no] + 
      (E->tracer.T[tr] -  E->tracer.PT0[E->tracer.property_group[tr]*MAX_MATERIAL_PHASES+phase_no] ) *
      E->tracer.Clapeyron[E->tracer.property_group[tr]*MAX_MATERIAL_PHASES+phase_no] / (E->data.grav_acc * 0.5 *
      (E->tracer.Density[E->tracer.property_group[tr]*MAX_MATERIAL_PHASES+phase_no]+
      E->tracer.Density[E->tracer.property_group[tr]*MAX_MATERIAL_PHASES+phase_no+1])); 

    if(E->tracer.tz[tr] < Depth) {
      found=phase_no;
      phase_b=phase_no;
      break;
    }
  } while (++phase_no < E->tracer.Phases[E->tracer.property_group[tr]]-1);

  if(found==-1) {
    phase_no=E->tracer.Phases[E->tracer.property_group[tr]]-1;  /* it must be in the highest pressure phase */
    phase_b=E->tracer.Phases[E->tracer.property_group[tr]]-2;
  }
  else {
    phase_no=found;
    phase_b=found;
  }

   E->tracer.phase_function[tr] = 0.5 * (1.0 + tanh(fabs(E->tracer.tz[tr] - Depth) / 250.0));
    
  if(E->tracer.tz[tr] < Depth) 
    E->tracer.phase_function[tr] *= -1.0;

  return(phase_no);
}
