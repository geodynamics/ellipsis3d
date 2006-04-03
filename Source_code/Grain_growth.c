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
#include "element_definitions.h"
#include "global_defs.h"
#include "function_prototypes.h"
#include <stdlib.h>

/* Functions to compute grain growth to use in rheological law.
   Applies if viscosity.GRAINSIZE != 0
   Setup routines called from viscosity routines.
 */

void grain_growth_initialize(
     struct All_variables *E
)
{
  int i,j,el,level,rheo,material;
  char default_str[40];

  /* Which formulation is to be used .... growth/nucleation or empirical */

  input_int("Grain_size_model",&(E->tracer.grain_size_model),"1");

  for(material=0;material<=E->tracer.NUM_MATERIALS;material++) {

    if(E->tracer.grain_size_model==1) {
      sprintf(default_str,"Material_%d_grsz_epsT",material);
      input_std_precision(default_str,&(E->tracer.grain[material].grsz_epsT),"0.0");
      sprintf(default_str,"Material_%d_grsz_B",material);
      input_std_precision(default_str,&(E->tracer.grain[material].grsz_B),"1.0");
      sprintf(default_str,"Material_%d_grsz_stsexp",material);
      input_std_precision(default_str,&(E->tracer.grain[material].grsz_stsexp),"1.0");

    }

    else if (E->tracer.grain_size_model==2) {

      sprintf(default_str,"Material_%d_reduction_factor_e",material);
      input_std_precision(default_str,&(E->tracer.grain[material].reduction_factor_equil),"1.0");
      sprintf(default_str,"Material_%d_reduction_factor_m",material);
      input_std_precision(default_str,&(E->tracer.grain[material].reduction_factor_metas),"1.0");
      
      sprintf(default_str,"Material_%d_grain_T_dep",material);
      input_int_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].GRAIN_T_dep,1);
      sprintf(default_str,"Material_%d_grain_gr_a",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].ggrw_a,0.0);
      sprintf(default_str,"Material_%d_grain_gr_m",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].ggrw_m,0.0);
      sprintf(default_str,"Material_%d_grain_gr_Q",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].ggrw_Q,0.0);
      sprintf(default_str,"Material_%d_grain_gr_T",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].ggrw_T,0.0);
      sprintf(default_str,"Material_%d_grain_gr_T0",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].ggrw_T0,0.0);
      sprintf(default_str,"Material_%d_grain_nu_a",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].gnuc_a,0.0);
      sprintf(default_str,"Material_%d_grain_nu_m",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].gnuc_m,0.0);
      sprintf(default_str,"Material_%d_grain_nu_meps",material);
      input_std_precision_vector(default_str,E->tracer.Phases[material],E->tracer.grain[material].gnuc_meps,0.0);
    }
  }
  return;
}


void grow_grains(
     struct All_variables *E
)

{
  int m,material,phase;
  standard_precision exp_term;
  standard_precision delta_grain_size;
  standard_precision grain_size_half;
  standard_precision d_inf;

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(E->tracer.property_group[m] < 0)
      continue;

    material = E->tracer.property_group[m]; 
    phase = E->tracer.Current_phase[m];
    
    /* Exponential component of grain growth may be
       linearized to match equivalent term in viscosity
       law */

    if(E->tracer.grain_size_model == 1) {

      if(E->tracer.grain[material].grsz_epsT <= 0.0)  /* fast way to switch off this effect */
	continue;

      d_inf = E->tracer.grain[material].grsz_B * 
	pow(E->tracer.strd[m], -E->tracer.grain[material].grsz_stsexp);

      /* fprintf(stderr,"%d: d_inf = %g (stress = %g) \n",m, d_inf, E->tracer.strd[m] ); */

      delta_grain_size = E->advection.timestep * 
	E->tracer.edot[m] / E->tracer.grain[material].grsz_epsT * (d_inf - E->tracer.grain_size[m]); 
    }
    
    else if (E->tracer.grain_size_model == 2) {

    /* 2nd Order method ... */

    switch (E->tracer.grain[material].GRAIN_T_dep[phase]) {
    case 1:  /* Linearized exponential */ 
      exp_term = exp(E->tracer.grain[material].ggrw_T[phase] * E->tracer.T[m]);
      break;
    case 2:
      exp_term = exp(-E->tracer.grain[material].ggrw_Q[phase]/
		     (E->tracer.grain[material].ggrw_T[phase] * 
		      (E->tracer.T[m] + E->tracer.grain[material].ggrw_T0[phase])));
      break;
    default:
      fprintf(stderr,"Grain growth exponential term undefined\n");
      exit(-1);

    }
 
    delta_grain_size = E->advection.timestep * 
      ( /* growth terms */
       E->tracer.grain[material].ggrw_a[phase] * 
       pow(E->tracer.grain_size[m],E->tracer.grain[material].ggrw_m[phase]) * 
       exp_term -

       /* reduction terms */
       E->tracer.grain[material].gnuc_a[phase] *  
       pow(E->tracer.grain_size[m],E->tracer.grain[material].gnuc_m[phase]) * 
       pow(E->tracer.edot[m],E->tracer.grain[material].gnuc_meps[phase])
       );
    }
  
    else
      delta_grain_size = 0.0;
 
    E->tracer.grain_size[m] +=  delta_grain_size;
  }
  return;
}


void phase_boundary_grain_size(
  struct All_variables *E,
  int tr,
  int meta
)
     
{
  int material;
  
  /* in simplified form, only change grain
     size for transformation from metastable state  */

  material = E->tracer.property_group[tr]; 

  /* E->tracer.grain_size[1] = 0.1; */

  /* 1. Check to see if a phase change has recently occurred */

  if(E->tracer.time_since_phase_change[tr] == 0.0) {
    if(meta)
      E->tracer.grain_size[tr] = E->tracer.grain[material].reduction_factor_metas;  /* *= */
    else
      E->tracer.grain_size[tr] = E->tracer.grain[material].reduction_factor_equil;  /* *= */
  }
  return;
}
