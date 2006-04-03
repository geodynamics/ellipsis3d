/*

Copyright (C) 2003 The GeoFramework Consortium

This file is part of Ellipsis3D.

Ellipsis3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, version 2,
as published by the Free Software Foundation.

Ellipsis3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Authors:
  Louis Moresi <louis.moresi@sci.monash.edu>
  Richard Albert <richard.albert@exxonmobil.com>

*/




/*  Here are the routines which process the results of each buoyancy solution, and call
    any relevant output routines. Much of the information has probably been output along
    with the velocity field. (So the velocity vectors and other data are fully in sync).
    However, heat fluxes and temperature averages are calculated here (even when they
    get output the next time around the velocity solver);
    */

#include <math.h>
/* #include <stdlib.h> */ /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_new_buoyancy_field(
 struct All_variables *E,
    int ii
)
{ 
    void heat_flux();

  if(E->control.verbose)
    fprintf(stderr,"Process temperature info\n");


    if(E->advection.ADVECTION) {
	heat_flux(E);
    }

    return_horiz_ave(E,E->T,E->Have.T);
    if(E->Have.T[E->mesh.noz] != 0.0)
      E->monitor.Nusselt /= E->Have.T[E->mesh.noz];

    if(E->control.verbose)
      fprintf(stderr,"Process temperature info ... done\n");

    return;
}

/* ===================
    Surface heat flux  
   =================== */

void heat_flux(
    struct All_variables *E
)
{
    standard_precision *dTdz,*dTdzb;
/*    standard_precision *uT,*uTb;  RAA, 26/3/01 where are these variables used? */
    int i,j,node,lnode;
    standard_precision return_layer_value();
    standard_precision adv_hfl,basal_adv_hfl;
    standard_precision z1,z2,F1,F2;

    static int been_here=0;
    
    dTdz  = (standard_precision *) Malloc0((sizeof(standard_precision)) * (1 + E->mesh.nox*(1+E->mesh.noy)));
    dTdzb = (standard_precision *) Malloc0((sizeof(standard_precision)) * (1 + E->mesh.nox*(1+E->mesh.noy)));
 
    for(j=1;j<=E->mesh.noy;j++) 
	for(i=1;i<=E->mesh.nox;i++)   {
	    lnode = i + (j-1)*E->mesh.nox;
	    node = 1 + (i-1)*E->mesh.noz + (j-1)*E->mesh.noz*E->mesh.nox;

	    dTdz[lnode]  = (E->T[node] - E->T[node+1]) / (E->sx[2][node] - E->sx[2][node+1]); 
	    E->slice.shflux[lnode] = dTdz[lnode] - 
	      0.25 * (E->V[2][node]+E->V[2][node+1]) * (E->T[node] + E->T[node+1]);  
	    
	    E->slice.shflux[lnode] = E->Tdot[node]*(E->x[2][node+1]-E->x[2][node])*
	     0.5/E->advection.temp_iterations;

	    node += E->mesh.noz-1;

	    dTdzb[lnode] = (E->T[node] - E->T[node-1]) / (E->sx[2][node] - E->sx[2][node-1]);
	    E->slice.bhflux[lnode] = dTdzb[lnode] - 
	      0.25 * (E->V[2][node]+E->V[2][node-1]) * (E->T[node] + E->T[node-1]);
	
	    E->slice.bhflux[lnode] = -E->Tdot[node]*(E->x[2][node]-E->x[2][node-1])*0.5/E->advection.temp_iterations;

	}
     
    if(been_here++==0)
      E->monitor.Nusselt = E->monitor.F_surface = E->monitor.F_base = 0.0; /* not solved T yet ! */
    else {
      E->monitor.Nusselt = 0.0;
      E->monitor.F_surface = return_layer_value(E,E->slice.shflux,1,1);
      E->monitor.F_base = return_layer_value(E,E->slice.bhflux,1,1);

     }
  free((void *) dTdz);
  free((void *) dTdzb);


  return;
  }
