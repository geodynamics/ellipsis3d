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



#include <math.h> 


#if (defined __sunos__)
#include <string.h>
#else
#if (!defined __GNUC__)
#include <strings.h> 
#endif
#endif

/* #if (! defined __GNUC__)
   #include <rpc/xdr.h> 
   #endif
*/

#include "element_definitions.h"
#include "global_defs.h"

void tracer_time_dep_terms(
			   struct All_variables *E
			   )
{
  int i,j,k,m,el,node,n;
  standard_precision *divUe,*divUt,*divUn;
  standard_precision eta1,eta2,eta3,dOmega;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision timestep,tp;
  const int level = E->mesh.levmax;  
  standard_precision V1,V2,V3;

  standard_precision *dUn1,*dUn2,*dUn3;
  standard_precision *Ut1,*Ut2,*Ut3;

  const int dims = E->mesh.nsd;

  higher_precision *node_R;
 
  static int been_here = 0;
  divUt = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision));
  divUn = (standard_precision *) Malloc0((E->mesh.nno + 1) * sizeof(standard_precision));
  divUe = (standard_precision *) Malloc0((E->mesh.nel + 1) * sizeof(standard_precision));

  if(3==dims) {
    dUn3 = (standard_precision *) Malloc0((E->mesh.nno + 1) * sizeof(standard_precision));
    Ut3 = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision));
  }

  /* Find average divergence of velocity over last time increment */

  for(el=1;el<=E->mesh.nel;el++) {
    get_global_v_x_shape_fn(E,el,lNx,&dOmega,0.0,0.0,0.0,E->mesh.levmax);
    divUe[el] = 0.0;
    for(i=1;i<=enodes[E->mesh.nsd];i++) {
      node = E->IEN[level][el].node[i];
  
      j = E->IEN[level][el].node[i];
      if(E->advection.timesteps == 1) {
        if(2==dims) /*RAA: 12/12/02, added this distinction*/
	  divUe[el] += lNx[1][i] * E->V[1][j]  + lNx[2][i] * E->V[2][j];
	else if(3==dims)  /*RAA: 12/12/02, added this part*/
	  divUe[el] += lNx[1][i] * E->V[1][j]  + lNx[2][i] * E->V[2][j] + lNx[3][i] * E->V[3][j];
      }
      else {
        if(2==dims)  /*RAA: 12/12/02, added this distinction*/
	   divUe[el] += 0.5 * (lNx[1][i] * (E->V1[1][j]+E->V[1][j])+
	        	       lNx[2][i] * (E->V1[2][j]+E->V[2][j]));
	else if(3==dims)  /*RAA: 12/12/02, added this part, make sure it is 1/3 & not 0.5*/
	   divUe[el] += (1.0/3.0) * (lNx[1][i] * (E->V1[1][j]+E->V[1][j])+
	        	             lNx[2][i] * (E->V1[2][j]+E->V[2][j])+
	        	             lNx[3][i] * (E->V1[3][j]+E->V[3][j]));
      }
    }
  /*raa: fprintf(E->fp1,"el: %d; timesteps: %d   DivUe = %g\n",el,E->advection.timesteps,divUe[el]);*/
  }

  sp_to_nodes(E,divUe,divUn,E->mesh.levmax);
  nodes_to_tracers(E,divUn,divUt,E->mesh.levmax);

    /* Elastic bulk modulus for numerical efficiency */

    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      if( E->tracer.visc[E->tracer.property_group[m]].Num_Elas_Bulk_mod > 0.0) {
	E->tracer.DQ[m] += E->tracer.visc[E->tracer.property_group[m]].Num_Elas_Bulk_mod * divUt[m] * E->advection.timestep;
	E->tracer.Density_change[m] *= 1 - divUt[m] * E->advection.timestep ;
      }
    }

  free((void *) divUt);
  free((void *) divUe);
  free((void *) divUn);
  /*RAA: 01/07/02 - added these lines to avoid memory leakage */
  if(3==dims) {
    free((void *) dUn3);
    free((void *) Ut3);
  }

  return;
}

void tracer_inertial_terms(
			   struct All_variables *E,
			   standard_precision *uu1,
			   standard_precision *uu2,
			   standard_precision *uu3,
			   int kk,
			   int level
)
{
  int i,j,k,m,el,node,n;
  standard_precision eta1,eta2,eta3,dOmega;
  standard_precision lN[ELNMAX+1];
  standard_precision timestep,tp;
  standard_precision T1,T2;

  standard_precision *Ut1,*Ut2,*Ut3;

  standard_precision Rey;

  standard_precision mag;

  const int dims = E->mesh.nsd;

  higher_precision *node_R;
 
  static int been_here = 0;
   
  Ut1 = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision));
  Ut2 = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision));
 
  if(3==dims) {
     Ut3 = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision));
  }
  /*  fprintf(stderr,"Inertial terms (K=%d\n",kk); */

  /* Acceleration of tracers ... */
 
  nodes_to_tracers(E,uu1,Ut1,level);
  nodes_to_tracers(E,uu2,Ut2,level);
  if(3==dims)  /*RAA: 5/4/01, added this call for 3D */ 
    nodes_to_tracers(E,uu3,Ut3,level);

  Rey = 0;

  mag = 0.0;
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {/* That'll become a Reynolds number or inverse of Prandtl number or something */
    if(E->advection.timesteps == 0) {
      E->tracer.A1[m] = 0.0;
      E->tracer.A2[m] = 0.0;
      if(3==dims)
	E->tracer.A3[m] = 0.0;
    }
    else if(E->advection.timesteps == 1) { 
      E->tracer.A1[m] =  Rey  *  (Ut1[m]-E->tracer.U1[m]) / E->advection.timestep;  
      E->tracer.A2[m] =  Rey  *  (Ut2[m]-E->tracer.U2[m]) / E->advection.timestep;

      mag +=  E->tracer.A1[m] *  E->tracer.A1[m] + E->tracer.A2[m] *  E->tracer.A2[m];

      if(3==dims) {
	E->tracer.A3[m] = 1.0e2 * (Ut3[m]-E->tracer.U3[m]) / E->advection.previous_timestep;
      }
    }
     else if(kk==1) { 
      E->tracer.A1[m] = 0.5 * E->tracer.A1[m] + 0.5 * Rey  *  (Ut1[m]-E->tracer.U1[m]) / E->advection.timestep;  
      E->tracer.A2[m] = 0.5 * E->tracer.A2[m] + 0.5 * Rey  *  (Ut2[m]-E->tracer.U2[m]) / E->advection.timestep;

      mag += E->tracer.A1[m] * E->tracer.A1[m] + E->tracer.A2[m] * E->tracer.A2[m];

      if(3==dims) {
	E->tracer.A3[m] = 1.0e2 * (Ut3[m]-E->tracer.U3[m]) / E->advection.previous_timestep;
      }
    }

    else if(E->advection.timesteps >= 1) {

      T1 =  E->advection.timestep;
      T2 =  E->advection.timestep + E->advection.previous_timestep;
 
      E->tracer.A1[m] = 0.0 * E->tracer.A1[m] - 1.0 * Rey *
	(-Ut1[m] * (T2*T2 - T1*T1) + T2*T2*E->tracer.U1[m] - T1*T1*E->tracer.UU1[m]) / ( T1 * T2 * (T2 - T1)) ;

      E->tracer.A2[m] = 0.0 * E->tracer.A2[m] - 1.0 * Rey *
	(-Ut2[m] * (T2*T2 - T1*T1) + T2*T2*E->tracer.U2[m] - T1*T1*E->tracer.UU2[m]) / ( T1 * T2 * (T2 - T1)) ;

      if(0 && m > 1000 && m < 1005) {
	fprintf(stderr,"%d: Accn term X ... V = %g  %g  %g  -> %g \n",m,Ut1[m],E->tracer.U1[m],E->tracer.UU1[m],E->tracer.A1[m] );
  	fprintf(stderr,"%d: Accn term Z ... V = %g  %g  %g  -> %g \n",m,Ut2[m],E->tracer.U2[m],E->tracer.UU2[m],E->tracer.A2[m] );
      }

      mag +=  E->tracer.A1[m] *  E->tracer.A1[m] + E->tracer.A2[m] *  E->tracer.A2[m];
      
      if(3==dims) {
	E->tracer.A3[m] = 1.0e2 * (Ut3[m]-E->tracer.U3[m]) / E->advection.previous_timestep;
      }
    }
  } /* next tracer */
  
  mag = sqrt(mag / E->tracer.NUM_TRACERS);

  free((void *) Ut1);
  free((void *) Ut2);

  if(3==dims) {
    free((void *) Ut3);
  }
  return;
}

