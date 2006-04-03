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

/* Functions to handle Hans' chemical transport model based on
   mobile + immobile phase pseudo-porous flow */



/* 1. get stress normal to layering at tracers */

void update_sigma_n( 
		    struct All_variables *E
		    )

{
  int m,i,j,k;
  standard_precision tau[7][7];
  standard_precision D[7][7];
  
  standard_precision pressure;
  standard_precision Dkk,*divU;
  standard_precision Hvisc;

  divU = (standard_precision *) Malloc0((E->mesh.nel+1) * sizeof(standard_precision));

  assemble_div_us3(E,E->V[1],E->V[2],E->V[3],divU,E->mesh.levmax);

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio == 1.0) {
      E->tracer.sigma_n[m] = 0.0;
      continue;
    }
 

    /* If elastic, then this is already available, essentially, in S11 etc */

    tracer_deviatoric_stress_and_pressure(E,tau,D,&pressure,&Dkk,&Hvisc,m,E->mesh.levmax);

    E->tracer.sigma_n[m] = 
      E->tracer.n1[m] *  E->tracer.n1[m] * tau[1][1] + 
      E->tracer.n2[m] *  E->tracer.n2[m] * tau[2][2] + 
      2.0 * E->tracer.n2[m] *  E->tracer.n1[m] * tau[1][2]  - 
      (pressure - E->tracer.visc[E->tracer.property_group[m]].Pen_bulk*E->tracer.Visc[m] * divU[E->tracer.tracer_elt[E->mesh.levmax][m]]) ;


    if(fabs(E->tracer.sigma_n[m]) > 1.0e10)
      fprintf(stderr,"Tracer %d, sigma_n = %g  (%g,%g,%g,%g,%g)\n",
	      m,E->tracer.sigma_n[m],E->tracer.n1[m],E->tracer.n2[m],
	      tau[1][1],tau[2][2],tau[1][2]); 

  }
  free((void *) divU) ;
  return;
}


void update_porosity(
		     struct All_variables *E
		     )
{

  standard_precision eltk[9][9];
  standard_precision S[4][4];
  standard_precision dOmega;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision time1;
  standard_precision eta1,eta2,eta3;

  standard_precision *dphi_n;
  standard_precision *dphi_t;
  standard_precision *sigman_n;
  standard_precision dphi_e[9];

  int el,a,b,i,j,m;

  const int level = E->mesh.levmax;
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];


  dphi_n = (standard_precision *)  Malloc0((E->mesh.nno + 1) * sizeof(standard_precision)); 
  dphi_t = (standard_precision *)  Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(standard_precision)); 
  sigman_n = (standard_precision *)  Malloc0((E->mesh.nno + 1) * sizeof(standard_precision)); 

  gs_tracers_to_nodes(E,sigman_n,NULL,NULL,NULL,E->tracer.sigma_n,level,0);

  for(i=1;i<=E->mesh.nno;i++) {
    fprintf(stderr,"Sigman[%d] = %g\n",i,sigman_n[i]);
    dphi_n[i] = 0.0;
  }

  for(el=1;el<=E->mesh.nel;el++) {

    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
 
      get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,NULL,level);

      dOmega  = 
	E->tracer.tracer_weight_fn[level][m] 
	* E->tracer.tracer_jacobian[level][m] * 
	E->tracer.permeability_factor[E->tracer.property_group[m]];

   
      S[1][1] = 1.0 - E->tracer.n1[m] * E->tracer.n1[m]; 
      S[2][2] = 1.0 - E->tracer.n2[m] * E->tracer.n2[m]; 
      S[1][2] =  - E->tracer.n2[m] * E->tracer.n1[m]; 

      for(a=1;a<=ends;a++) {
	dphi_e[a] = 0.0;
	for(b=1;b<=ends;b++) 
	  eltk[a][b] = 0.0;
      }
      for(a=1;a<=ends;a++) {
	for(b=1;b<=ends;b++) {
	  eltk[a][b] += dOmega * ( lNx[1][a] * S[1][1] * lNx[1][b] +
				   lNx[2][a] * S[1][2] * lNx[1][b] +
				   lNx[1][a] * S[1][2] * lNx[2][b] +
				   lNx[2][a] * S[2][2] * lNx[2][b] );
	
	}
      }
    }


    for(a=1;a<=ends;a++) 
      for(b=1;b<=ends;b++) 
	dphi_e[a] += eltk[a][b] * sigman_n[E->ien[el].node[b]];
    
    for(a=1;a<=ends;a++) {
      dphi_n[E->ien[el].node[a]] += dphi_e[a];
    }

  } /* Next element */


  time1 = E->advection.timestep;
  
  for(i=1;i<=E->mesh.nno;i++) {
    dphi_n[i] = - time1 * dphi_n[i] / E->mass[i];
  }


    /* Tracer value of delta phi */

    nodes_to_tracers(E,dphi_n,dphi_t,level);
    
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      E->tracer.volfraction[m] += dphi_t[m];

      if(E->tracer.volfraction[m] > 1.0)
	E->tracer.volfraction[m] = 1.0;
      if(E->tracer.volfraction[m] < 0.0)
	E->tracer.volfraction[m] = 0.0;
    }


  free((void *) dphi_n);
  free((void *) dphi_t);
  free((void *) sigman_n);

  return;
  }
