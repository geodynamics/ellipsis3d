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


/*  */ 
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


/* Copy one tracer (at location tr1) to another tracer (at location tr2) */


void tracer_copy(
		 struct All_variables *E,
		 int tr1,
		 int tr2
		 )
{
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  int level,q;

  E->tracer.tx[tr2] = E->tracer.tx[tr1] ;
  E->tracer.tz[tr2] = E->tracer.tz[tr1] ;
  if(3==dims)
     E->tracer.ty[tr2] = E->tracer.ty[tr1];

  E->tracer.tx1[tr2] = E->tracer.tx1[tr1] ;
  E->tracer.tz1[tr2] = E->tracer.tz1[tr1] ;
  if(3==dims)
     E->tracer.ty1[tr2] = E->tracer.ty1[tr1];

  E->tracer.T[tr2] = E->tracer.T[tr1];
  E->tracer.Q[tr2] = E->tracer.Q[tr1];
  E->tracer.Hvisc[tr2] = E->tracer.Hvisc[tr1];
  E->tracer.weight[tr2] = E->tracer.weight[tr1];
  E->tracer.Density_change[tr2] = E->tracer.Density_change[tr1];
  E->tracer.edot[tr2] = E->tracer.edot[tr1];
  E->tracer.strs[tr2] = E->tracer.strs[tr1];
  E->tracer.strd[tr2] = E->tracer.strd[tr1];
  E->tracer.depl[tr2] = E->tracer.depl[tr1]; /*RAA: 24/09/02, O'Neill - melting stuff */
  E->tracer.dFdot[tr2] = E->tracer.dFdot[tr1]; /*RAA: 10/07/03, O'Neill - melting stuff */
  E->tracer.edot_integrated[tr2] = E->tracer.edot_integrated[tr1];
  E->tracer.edotp_integrated[tr2] = E->tracer.edotp_integrated[tr1];
  E->tracer.time_since_phase_change[tr2] = E->tracer.time_since_phase_change[tr1];
  E->tracer.phase_function[tr2] = E->tracer.phase_function[tr1];
  E->tracer.grain_size[tr2] = E->tracer.grain_size[tr1];

  E->tracer.yielded[tr2] = E->tracer.yielded[tr1];
  E->tracer.property_group[tr2] = E->tracer.property_group[tr1];

  E->tracer.tracer_weight_fn[E->mesh.levmax][tr2] = 
    E->tracer.tracer_weight_fn[E->mesh.levmax][tr1];
   E->tracer.tracer_std_weighting[E->mesh.levmax][tr2] = 
    E->tracer.tracer_std_weighting[E->mesh.levmax][tr1];
  E->tracer.tracer_jacobian[E->mesh.levmax][tr2] = 
    E->tracer.tracer_jacobian[E->mesh.levmax][tr1];

  E->tracer.tr_dim[tr2] =  E->tracer.tr_dim[tr1];
  E->tracer.DQ[tr2] =  E->tracer.DQ[tr1];
  E->tracer.DQ1[tr2] =  E->tracer.DQ1[tr1];

  if(E->control.ORTHOTROPY) {
    E->tracer.n1[tr2] = E->tracer.n1[tr1];
    E->tracer.n2[tr2] = E->tracer.n2[tr1];
    if(3==dims) {
      E->tracer.n3[tr2] = E->tracer.n3[tr1];
    }
  }

#if defined (CHEM_TRANS_MODULE_INSTALLED)
  if(E->control.CHEM_TRANS) {
    E->tracer.volfraction[tr2] = E->tracer.volfraction[tr1];
    E->tracer.sigma_n[tr2] = E->tracer.sigma_n[tr1];
  }
  
#endif

  if(E->control.ELASTICITY) {
    E->tracer.Pt[tr2] =  E->tracer.Pt[tr1];
    E->tracer.S11[tr2] =  E->tracer.S11[tr1];
    E->tracer.S12[tr2] =  E->tracer.S12[tr1];
    E->tracer.S22[tr2] =  E->tracer.S22[tr1];
    

    if(3==dims & 3==dofs) {
      E->tracer.S13[tr2] =  E->tracer.S13[tr1];
      E->tracer.S23[tr2] =  E->tracer.S23[tr1];
      E->tracer.S33[tr2] =  E->tracer.S33[tr1];
    }
    else if(2==dims & 3==dofs) {
      E->tracer.S21[tr2] =  E->tracer.S21[tr1];
      E->tracer.M31[tr2] =  E->tracer.M31[tr1];
      E->tracer.M32[tr2] =  E->tracer.M32[tr1];
    }
    else if(6==dofs) {
      E->tracer.S13[tr2] =  E->tracer.S13[tr1];
      E->tracer.S31[tr2] =  E->tracer.S31[tr1];
      E->tracer.S23[tr2] =  E->tracer.S23[tr1];
      E->tracer.S32[tr2] =  E->tracer.S32[tr1];
      E->tracer.S33[tr2] =  E->tracer.S33[tr1];
      E->tracer.S21[tr2] =  E->tracer.S21[tr1];
      E->tracer.M31[tr2] =  E->tracer.M31[tr1];
      E->tracer.M32[tr2] =  E->tracer.M32[tr1];
      E->tracer.M12[tr2] =  E->tracer.M12[tr1];
      E->tracer.M13[tr2] =  E->tracer.M13[tr1];
      E->tracer.M23[tr2] =  E->tracer.M23[tr1];
      E->tracer.M21[tr2] =  E->tracer.M21[tr1];
    }
  }


  E->tracer.A1[tr2] =  E->tracer.A1[tr1];
  E->tracer.A2[tr2] =  E->tracer.A2[tr1];
  if(3==dims) 
    E->tracer.A3[tr2] =  E->tracer.A3[tr1];

  E->tracer.U1[tr2] =  E->tracer.U1[tr1];
  E->tracer.U2[tr2] =  E->tracer.U2[tr1];
  if(3==dims) 
    E->tracer.U3[tr2] =  E->tracer.U3[tr1];

  E->tracer.UU1[tr2] =  E->tracer.UU1[tr1];
  E->tracer.UU2[tr2] =  E->tracer.UU2[tr1];
  if(3==dims) 
    E->tracer.UU3[tr2] =  E->tracer.UU3[tr1];

  E->tracer.dX11[tr2] = E->tracer.dX11[tr1];
  E->tracer.dX12[tr2] = E->tracer.dX12[tr1];
  E->tracer.dX21[tr2] = E->tracer.dX21[tr1];
  E->tracer.dX22[tr2] = E->tracer.dX22[tr1];

  if(3==dims) {
    E->tracer.dX13[tr2] = E->tracer.dX13[tr1];
    E->tracer.dX23[tr2] = E->tracer.dX23[tr1];
    E->tracer.dX31[tr2] = E->tracer.dX31[tr1];
    E->tracer.dX32[tr2] = E->tracer.dX32[tr1];
    E->tracer.dX33[tr2] = E->tracer.dX33[tr1];
  }
 
  E->tracer.Current_phase[tr2] = E->tracer.Current_phase[tr1];
  E->tracer.Equilibrium_phase[tr2] = E->tracer.Equilibrium_phase[tr1];

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
    E->tracer.tracer_elt[level][tr2] = E->tracer.tracer_elt[level][tr1];
 //printf("In tracer copy, no free statements here \n");

  return;
}

/* Merge two tracers tr1, tr2 to location tr1 (provided they are of the same property group).
   Then fill in the hole left by tr2 with the tracer at the top of the stack */

void tracer_merge(
		  struct All_variables *E,
		  int tr1,
		  int tr2
		  )
{
  void tracer_copy();

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;

  if(E->tracer.property_group[tr1] != E->tracer.property_group[tr2]) {
    fprintf(stderr,"ERROR trying to merge tracers of different kinds %d & %d \n",tr1,tr2);
    return;
  }

  /* Volume (mass) is conserved */
  E->tracer.tracer_weight_fn[E->mesh.levmax][tr1] = 
    (E->tracer.tracer_weight_fn[E->mesh.levmax][tr1] + E->tracer.tracer_weight_fn[E->mesh.levmax][tr2]);
  E->tracer.tracer_std_weighting[E->mesh.levmax][tr1] = 
    E->tracer.tracer_std_weighting[E->mesh.levmax][tr1]+E->tracer.tracer_std_weighting[E->mesh.levmax][tr2];
  E->tracer.tracer_jacobian[E->mesh.levmax][tr1] = 
    (E->tracer.tracer_jacobian[E->mesh.levmax][tr1] + E->tracer.tracer_weight_fn[E->mesh.levmax][tr2]);

  
  /* These properties are averaged */

  E->tracer.yielded[tr1] = 0.5 * (E->tracer.yielded[tr1] + E->tracer.yielded[tr2]);

  E->tracer.tx[tr1] = 0.5 * (E->tracer.tx[tr1]+E->tracer.tx[tr2]);
  E->tracer.tz[tr1] = 0.5 * (E->tracer.tz[tr1]+E->tracer.tz[tr2]);
  if(3==dims)
    E->tracer.ty[tr1] = 0.5 * (E->tracer.ty[tr1]+E->tracer.ty[tr2]);

  E->tracer.tx1[tr1] = 0.5 * (E->tracer.tx1[tr1]+E->tracer.tx1[tr2]);
  E->tracer.tz1[tr1] = 0.5 * (E->tracer.tz1[tr1]+E->tracer.tz1[tr2]);
  if(3==dims)
    E->tracer.ty1[tr1] = 0.5 * (E->tracer.ty1[tr1]+E->tracer.ty1[tr2]);


  if(E->control.ELASTICITY) {

    E->tracer.S11[tr1] = 0.5 * (E->tracer.S11[tr1] + E->tracer.S11[tr2]);
    E->tracer.S12[tr1] = 0.5 * (E->tracer.S12[tr1] + E->tracer.S12[tr2]);
    E->tracer.S22[tr1] = 0.5 * (E->tracer.S22[tr1] + E->tracer.S22[tr2]);
    if(3==dims & 3==dofs) {
      E->tracer.S13[tr1] = 0.5 * (E->tracer.S13[tr1] + E->tracer.S13[tr2]);
      E->tracer.S23[tr1] = 0.5 * (E->tracer.S23[tr1] + E->tracer.S23[tr2]);
      E->tracer.S33[tr1] = 0.5 * (E->tracer.S33[tr1] + E->tracer.S33[tr2]) ;
    }
    else if(2==dims & 3==dofs) {
      E->tracer.S21[tr1] = 0.5 * (E->tracer.S21[tr1] + E->tracer.S21[tr2]);
      E->tracer.M31[tr1] = 0.5 * (E->tracer.M31[tr1] + E->tracer.M31[tr2]);
      E->tracer.M32[tr1] = 0.5 * (E->tracer.M32[tr1] + E->tracer.M32[tr2]);
    }
    else if(6==dofs) {
      E->tracer.S13[tr1] = 0.5 * (E->tracer.S13[tr1] + E->tracer.S13[tr2]);
      E->tracer.S31[tr1] = 0.5 * (E->tracer.S31[tr1] + E->tracer.S31[tr2]);
      E->tracer.S23[tr1] = 0.5 * (E->tracer.S23[tr1] + E->tracer.S23[tr2]);
      E->tracer.S32[tr1] = 0.5 * (E->tracer.S32[tr1] + E->tracer.S32[tr2]);
      E->tracer.S33[tr1] = 0.5 * (E->tracer.S33[tr1] + E->tracer.S33[tr2]);
      E->tracer.S21[tr1] = 0.5 * (E->tracer.S21[tr1] + E->tracer.S21[tr2]);
      E->tracer.M31[tr1] = 0.5 * (E->tracer.M31[tr1] + E->tracer.M31[tr2]);
      E->tracer.M32[tr1] = 0.5 * (E->tracer.M32[tr1] + E->tracer.M32[tr2]);
      E->tracer.M12[tr1] = 0.5 * (E->tracer.M12[tr1] + E->tracer.M12[tr2]);
      E->tracer.M13[tr1] = 0.5 * (E->tracer.M13[tr1] + E->tracer.M13[tr2]);
      E->tracer.M23[tr1] = 0.5 * (E->tracer.M23[tr1] + E->tracer.M23[tr2]);
      E->tracer.M21[tr1] = 0.5 * (E->tracer.M21[tr1] + E->tracer.M21[tr2]);
    }

    E->tracer.Pt[tr1] =  0.5 * (E->tracer.Pt[tr1]+E->tracer.Pt[tr2]);
  }

  if(E->control.ORTHOTROPY) {
     E->tracer.n1[tr1] =  0.5 * (E->tracer.n1[tr1]+E->tracer.n1[tr2]);
     E->tracer.n2[tr1] =  0.5 * (E->tracer.n2[tr1]+E->tracer.n2[tr2]);
     if(3==dims) 
       E->tracer.n3[tr1] =  0.5 * (E->tracer.n3[tr1]+E->tracer.n3[tr2]);
  }


#if defined (CHEM_TRANS_MODULE_INSTALLED)
  if(E->control.CHEM_TRANS) {
    E->tracer.volfraction[tr1] =  0.5 * (E->tracer.volfraction[tr1]+E->tracer.volfraction[tr2]);
    E->tracer.sigma_n[tr1] =  0.5 * (E->tracer.sigma_n[tr1]+E->tracer.sigma_n[tr2]);
  }
#endif

  E->tracer.T[tr1] = 0.5 * (E->tracer.T[tr1]+E->tracer.T[tr2]);
  E->tracer.Q[tr1] = 0.5 * (E->tracer.Q[tr1]+E->tracer.Q[tr2]);
  E->tracer.Hvisc[tr1] = 0.5 * (E->tracer.Hvisc[tr1]+E->tracer.Hvisc[tr2]);
  E->tracer.weight[tr1] = 0.5 * (E->tracer.weight[tr1]+E->tracer.weight[tr2]);
  E->tracer.Density_change[tr1] = 0.5 * (E->tracer.Density_change[tr1]+E->tracer.Density_change[tr2]);
  E->tracer.edot[tr1] = 0.5 * (E->tracer.edot[tr1]+E->tracer.edot[tr2]);
  E->tracer.strs[tr1] = 0.5 * (E->tracer.strs[tr1]+E->tracer.strs[tr2]);
  E->tracer.strd[tr1] = 0.5 * (E->tracer.strd[tr1]+E->tracer.strd[tr2]);
  E->tracer.depl[tr1] = 0.5 * (E->tracer.depl[tr1]+E->tracer.depl[tr2]); /*RAA: 24/09/02, C. O'Neill - melting stuff */
  E->tracer.dFdot[tr1] = 0.5 * (E->tracer.dFdot[tr1]+E->tracer.dFdot[tr2]);/*RAA: 10/07/03, O'Neill - melting stuff */
  E->tracer.edot_integrated[tr1] = 0.5 * (E->tracer.edot_integrated[tr1]+E->tracer.edot_integrated[tr2]);
  E->tracer.edotp_integrated[tr1] = 0.5 * (E->tracer.edotp_integrated[tr1]+E->tracer.edotp_integrated[tr2]);

  E->tracer.tr_dim[tr2] =  0.5 * (E->tracer.tr_dim[tr1]+E->tracer.tr_dim[tr2]);
  E->tracer.DQ[tr2] =  0.5 * (E->tracer.DQ[tr1]+E->tracer.DQ[tr2]);
  E->tracer.DQ1[tr2] =  0.5 * (E->tracer.DQ1[tr1]+E->tracer.DQ1[tr2]);
  
  E->tracer.A1[tr2] =  0.5 * (E->tracer.A1[tr1]+E->tracer.A1[tr2]);
  E->tracer.A2[tr2] =  0.5 * (E->tracer.A2[tr1]+E->tracer.A2[tr2]);
  if(3==dims) 
    E->tracer.A3[tr2] =  0.5 * (E->tracer.A3[tr1]+E->tracer.A3[tr2]);

  E->tracer.U1[tr2] =  0.5 * (E->tracer.U1[tr1]+E->tracer.U1[tr2]);
  E->tracer.U2[tr2] =  0.5 * (E->tracer.U2[tr1]+E->tracer.U2[tr2]);
  if(3==dims) 
    E->tracer.U3[tr2] =  0.5 * (E->tracer.U3[tr1]+E->tracer.U3[tr2]);

  E->tracer.UU1[tr2] =  0.5 * (E->tracer.UU1[tr1]+E->tracer.UU1[tr2]);
  E->tracer.UU2[tr2] =  0.5 * (E->tracer.UU2[tr1]+E->tracer.UU2[tr2]);
  if(3==dims) 
    E->tracer.UU3[tr2] =  0.5 * (E->tracer.UU3[tr1]+E->tracer.UU3[tr2]);
  
  E->tracer.dX11[tr1] = 0.5 * (E->tracer.dX11[tr1]+E->tracer.dX11[tr2]);
  E->tracer.dX12[tr1] = 0.5 * (E->tracer.dX12[tr1]+E->tracer.dX12[tr2]);
  E->tracer.dX21[tr1] = 0.5 * (E->tracer.dX21[tr1]+E->tracer.dX21[tr2]);
  E->tracer.dX22[tr1] = 0.5 * (E->tracer.dX22[tr1]+E->tracer.dX22[tr2]);

  if(3==dims) {
    E->tracer.dX13[tr1] = 0.5 * (E->tracer.dX13[tr1]+E->tracer.dX13[tr2]);
    E->tracer.dX23[tr1] = 0.5 * (E->tracer.dX23[tr1]+E->tracer.dX23[tr2]);
    E->tracer.dX31[tr1] = 0.5 * (E->tracer.dX31[tr1]+E->tracer.dX31[tr2]);
    E->tracer.dX32[tr1] = 0.5 * (E->tracer.dX32[tr1]+E->tracer.dX32[tr2]);
    E->tracer.dX33[tr1] = 0.5 * (E->tracer.dX33[tr1]+E->tracer.dX33[tr2]);
  }

  /* Well, this might be ok, might need some refining */
  E->tracer.Current_phase[tr1] = (int) (0.5 * (E->tracer.Current_phase[tr1]+E->tracer.Current_phase[tr2]));
  E->tracer.Equilibrium_phase[tr1] = (int) 
    (0.5 * (E->tracer.Equilibrium_phase[tr1]+E->tracer.Equilibrium_phase[tr2]));
  E->tracer.time_since_phase_change[tr1] = 0.5 * 
    (E->tracer.time_since_phase_change[tr1] + E->tracer.time_since_phase_change[tr2]);
  E->tracer.phase_function[tr1] = 0.5 * 
    (E->tracer.phase_function[tr1] + E->tracer.phase_function[tr2]);
  E->tracer.grain_size[tr1] = 0.5 * 
    (E->tracer.grain_size[tr1] + E->tracer.grain_size[tr2]);


  /* Eta, shape functions etc are not merged - these must be
     recalculated !! */

  tracer_copy(E,E->tracer.NUM_TRACERS,tr2);
  E->tracer.NUM_TRACERS--;
  printf("In tracer merge, no frees \n");
}


/* Routine to compute tracers splitting and 
   eating each other.

   To minimize the searches, we assume that 
   tracers get hungry when they are about to
   split into two. This also keeps the number of
   tracers close to being constant

   Like black widow spiders, the bloated
   tracer goes in search
   of a mate and then eats him. 
*/


/* ?? CAN THIS BE DONE INDEPENDENTLY
   FOR THE TRACERS IN EACH ELEMENT */

void tracer_black_widow 
(
 struct All_variables *E,
 int m 
 )
{
  standard_precision distortion1;
  standard_precision distortion2;
  standard_precision distortion3;
  standard_precision distance;

  standard_precision delX,delZ,delY;
  standard_precision element_size2;

  standard_precision eta1,eta2,eta3;

  int horribly_distorted;

  const int dims = E->mesh.nsd;

  int el,i,j,tr;

  /* First obtain a measure of distortion along
     each of the original principal directions */

  el = E->tracer.tracer_elt[E->mesh.levmax][m];


  if(3==dims) {
    distortion1 = 
      E->tracer.dX11[m]*E->tracer.dX11[m]+
      E->tracer.dX12[m]*E->tracer.dX12[m]+
      E->tracer.dX13[m]*E->tracer.dX13[m];
    distortion2 = 
      E->tracer.dX21[m]*E->tracer.dX21[m]+
      E->tracer.dX22[m]*E->tracer.dX22[m]+
      E->tracer.dX23[m]*E->tracer.dX23[m];
    distortion3 = 
      E->tracer.dX31[m]*E->tracer.dX31[m]+
      E->tracer.dX32[m]*E->tracer.dX32[m]+
      E->tracer.dX33[m]*E->tracer.dX33[m];
    element_size2 = 0.3333 * (
      E->eco[el].size[1] * E->eco[el].size[1] +
      E->eco[el].size[2] * E->eco[el].size[2] +
      E->eco[el].size[3] * E->eco[el].size[3]);
  }
  else {
    distortion1 = 
      E->tracer.dX11[m]*E->tracer.dX11[m]+
      E->tracer.dX12[m]*E->tracer.dX12[m];
    distortion2 = 
      E->tracer.dX21[m]*E->tracer.dX21[m]+
      E->tracer.dX22[m]*E->tracer.dX22[m];
    distortion3 = 0.0;
    element_size2 = 0.1 * (
      E->eco[el].size[1] * E->eco[el].size[1] +
      E->eco[el].size[2] * E->eco[el].size[2]) ; 
  }

   element_size2 =  4.0 * E->tracer.tr_dim[m] * E->tracer.tr_dim[m]; 

  /* Otherwise, the particle is distorted and it
     is HUNGRY. If there is another particle
     nearby it will eat it.
  */

#if 1
  for(i=0;i<E->tracer.tr_in_element_number[E->mesh.levmax][el];i++) {
      tr = E->tracer.tr_in_element[E->mesh.levmax][i+E->tracer.tr_in_element_offset[E->mesh.levmax][el]];
      if(tr != m && tr != -1 && tr <= E->tracer.NUM_TRACERS) { 
	/*  m: Obviously it's bad to kill and eat yourself.
	    -1: Ghost- has been eaten already */

	distance = 
	  (E->tracer.tx[m]-E->tracer.tx[tr]) * (E->tracer.tx[m]-E->tracer.tx[tr]) +
	  (E->tracer.tz[m]-E->tracer.tz[tr]) * (E->tracer.tz[m]-E->tracer.tz[tr]) +
	  ((2==dims) ? 0.0 :  (E->tracer.ty[m]-E->tracer.ty[tr]) * (E->tracer.ty[m]-E->tracer.ty[tr]));

	if((E->tracer.property_group[m] == E->tracer.property_group[tr]) &&
	   (E->tracer.Current_phase[m] == E->tracer.Current_phase[tr] ) && 
	   (distance < 
	    (E->tracer.tr_dim[m] + E->tracer.tr_dim[tr]) *
	    (E->tracer.tr_dim[m] + E->tracer.tr_dim[tr]) * 
	    E->tracer.particle_appetite )) {


	  tracer_merge(E,m,tr);
	  	  
	  E->tracer.tr_in_element[E->mesh.levmax][i+E->tracer.tr_in_element_offset[E->mesh.levmax][el]] = -1;
	  break; /* don't be greedy  */
	}
      }
  } 
#endif

  if ((distortion1 < element_size2) &&
      (distortion2 < element_size2) &&
      (distortion3 < element_size2) ) {
    
    /* Fat, round, happy  particle */
	// printf("Black widow routine \n"); 
    return;
  }


  /* fprintf(stderr,"Breaking apart particle %d\n",m); */


  /* Now split particle according to its distortion. We work on
     the basis that the individual directions are unlikely to have
     all distorted by the same amount and let only one of the directions
     be triggered. This shouldn't matter because if any two are equal they
     will end up being aligned 
  */

  if(3 == dims) {
    if(distortion1 > element_size2) {
      delX = (E->tracer.dX11[m] *= 0.5);
      delZ = (E->tracer.dX12[m] *= 0.5);
      delY = (E->tracer.dX13[m] *= 0.5);
    }
    else if (distortion2 > element_size2) {
      delX = (E->tracer.dX21[m] *= 0.5);
      delZ = (E->tracer.dX22[m] *= 0.5);
      delY = (E->tracer.dX23[m] *= 0.5);
    }
    else /* (distortion3 > element_size2) */ {
      delX = (E->tracer.dX31[m] *= 0.5);
      delZ = (E->tracer.dX32[m] *= 0.5);
      delY = (E->tracer.dX33[m] *= 0.5);
    }

    E->tracer.tracer_weight_fn[E->mesh.levmax][m] *= 0.5;
    E->tracer.tracer_std_weighting[E->mesh.levmax][m] *= 0.5;
    E->tracer.NUM_TRACERS++;
    tracer_allocate_memory(E);
    tracer_copy(E,m,E->tracer.NUM_TRACERS);  


    /* Shift new and old tracers apart */

    E->tracer.tx[E->tracer.NUM_TRACERS]  += delX;
    E->tracer.tx1[E->tracer.NUM_TRACERS] += delX;
    E->tracer.tz[E->tracer.NUM_TRACERS]  += delZ;
    E->tracer.tz1[E->tracer.NUM_TRACERS] += delZ;
    E->tracer.ty[E->tracer.NUM_TRACERS]  += delY;
    E->tracer.ty1[E->tracer.NUM_TRACERS] += delY;

    E->tracer.tx[m]  -= delX;
    E->tracer.tx1[m] -= delX;
    E->tracer.tz[m]  -= delZ; 
    E->tracer.tz1[m] -= delZ; 
    E->tracer.ty[m]  -= delY;
    E->tracer.ty1[m] -= delY; 
    
    /* The nature of this procedure is that it 
       can violate the otherwise perfect obedience to
       boundary conditions 
       (also the particle-track crossing I suppose)
       so this must be checked.
    */

    el = get_tracers_element(E,m,&eta1,&eta2,&eta3,E->mesh.levmax);
    if(el==-1) {
      E->tracer.tx[m] += delX;
      E->tracer.tz[m] += delZ; 
      E->tracer.ty[m] += delY; 
      E->tracer.tx1[m] += delX;
      E->tracer.tz1[m] += delZ; 
      E->tracer.ty1[m] += delY; 
      get_tracers_element(E,m,&eta1,&eta2,&eta3,E->mesh.levmax);
    }

    el = get_tracers_element(E,E->tracer.NUM_TRACERS,&eta1,&eta2,&eta3,E->mesh.levmax);
    if(el==-1) {
      E->tracer.tx[E->tracer.NUM_TRACERS] -= delX;
      E->tracer.tz[E->tracer.NUM_TRACERS] -= delZ; 
      E->tracer.ty[E->tracer.NUM_TRACERS] -= delY; 
      E->tracer.tx1[E->tracer.NUM_TRACERS] -= delX;
      E->tracer.tz1[E->tracer.NUM_TRACERS] -= delZ; 
      E->tracer.ty1[E->tracer.NUM_TRACERS] -= delY; 
      get_tracers_element(E,E->tracer.NUM_TRACERS,&eta1,&eta2,&eta3,E->mesh.levmax);
    }
  }
  else /* (2 == dims) */ {                         /* TWO DIMENSIONAL */
    if(distortion1 > element_size2) {
      delX = 1.0 * (E->tracer.dX11[m]);
      delZ = 1.0 * (E->tracer.dX12[m]);
      E->tracer.dX11[m] *= 0.5;
      E->tracer.dX12[m] *= 0.5;
    }
    else /* (distortion2 > element_size2) */ {
      delX = 1.0 * (E->tracer.dX21[m]);
      delZ = 1.0 * (E->tracer.dX22[m]);
      E->tracer.dX21[m] *= 0.5;
      E->tracer.dX22[m] *= 0.5;
   }

    E->tracer.tracer_weight_fn[E->mesh.levmax][m] *= 0.5;
    E->tracer.tracer_std_weighting[E->mesh.levmax][m] *= 0.5;
    E->tracer.NUM_TRACERS++;
    tracer_allocate_memory(E);
    tracer_copy(E,m,E->tracer.NUM_TRACERS);  

    /* fprintf(stderr,"%d: Copy from %d: edotp = %g, (%g,%g) yielded=%g\n",
	      E->tracer.NUM_TRACERS,m,
	      E->tracer.edotp_integrated[m],
	      E->tracer.tx[m],
	      E->tracer.tz[m],
	      E->tracer.yielded[m]); */

    
    /* Shift new and old tracers apart */
  
    E->tracer.tx[E->tracer.NUM_TRACERS]  += delX;
    E->tracer.tx1[E->tracer.NUM_TRACERS] += delX;
    E->tracer.tz[E->tracer.NUM_TRACERS]  += delZ;
    E->tracer.tz1[E->tracer.NUM_TRACERS] += delZ;
 
    E->tracer.tx[m]  -= delX;
    E->tracer.tx1[m] -= delX;
    E->tracer.tz[m]  -= delZ; 
    E->tracer.tz1[m] -= delZ; 
     
    /* The nature of this procedure is that it 
       can violate the otherwise perfect obedience to
       boundary conditions 
       (also the particle-track crossing I suppose)
       so this must be checked.
    */

   

    el = get_tracers_element(E,m,&eta1,&eta2,&eta3,E->mesh.levmax);
    if(el==-1) {
      E->tracer.tx[m] += delX;
      E->tracer.tz[m] += delZ; 
      E->tracer.tx1[m] += delX;
      E->tracer.tz1[m] += delZ; 
      get_tracers_element(E,m,&eta1,&eta2,&eta3,E->mesh.levmax);
    }

    el = get_tracers_element(E,E->tracer.NUM_TRACERS,&eta1,&eta2,&eta3,E->mesh.levmax);
    if(el==-1) {
       E->tracer.tx[E->tracer.NUM_TRACERS] -= delX;
      E->tracer.tz[E->tracer.NUM_TRACERS] -= delZ; 
      E->tracer.tx1[E->tracer.NUM_TRACERS] -= delX;
      E->tracer.tz1[E->tracer.NUM_TRACERS] -= delZ; 
      get_tracers_element(E,E->tracer.NUM_TRACERS,&eta1,&eta2,&eta3,E->mesh.levmax);
    }
  }

 
  
  return;
}

/* Function to parcel up the element to give a good idea of
   the representative volume occupied by each tracer. 

   This version computes a massive number of distances to see
   which particles are close to which sub-regions of the element. 
   Although this seems wasteful it does give a result which exactly
   fills the element without scaling (useful ?). It is an approximation
   to computing the Voronoi diagram and doesn't worry about the edges.
   It's also not as slow as you might think at first !

*/


void tracer_representative_volumes (
				    struct All_variables *E,
				    standard_precision *ETA1,
				    standard_precision *ETA2,
				    standard_precision *ETA3 
				    )
{
  int el;
  int tr,m,n,i,j,k,q;
  int level;
  int subel,near_t,trac,new;
  int subelt_tracer2[25][25];
  int subelt_tracer3[25][25][25];
  
  standard_precision xx,zz,yy,rad;
  standard_precision subelt_dist2[25][25];
  standard_precision subelt_dist3[25][25][25];
  standard_precision eta1,eta2,eta3;
  standard_precision radx,radz;
 
  const int SUBD = 9;
  const standard_precision subd2_1 = 4.0 / (SUBD * SUBD);
  const standard_precision subd3_1 = 8.0 / (SUBD * SUBD * SUBD);
  const standard_precision subd1_1 =  1.0 / SUBD;


  /* What works best ?   At the moment it seems best to
     use the initial particle volumes but adjust them (elsewhere)
     at each timestep to fit the constraints.

     The alternative - to recompute the particle volumes each time -
     works about as well, but takes considerable time. However,
     this may be different in special circumstances so keep both
     options available 
  */

#if 0

  if(3==E->mesh.nsd) {
      for(el=1;el<=E->mesh.nel;el++) {
      for(i=1;i<=SUBD;i++)
	for(j=1;j<=SUBD;j++) 
	  for(k=1;k<=SUBD;k++) {
	  subelt_dist3[i][j][k] = 1.0e32;
	  subelt_tracer3[i][j][k]=0;
	}
    
      for(tr=0;tr<E->tracer.tr_in_element_number[E->mesh.levmax][el];tr++) {
      
	q = E->tracer.tr_in_element_offset[E->mesh.levmax][el] + tr;
	m = E->tracer.tr_in_element[E->mesh.levmax][q]; 
    
	eta1 = ETA1[m];
	eta2 = ETA2[m];
	eta3 = ETA3[m];
      
	E->tracer.tracer_weight_fn[E->mesh.levmax][m] = 0.0;

	for(i=1;i<=SUBD;i++) {
	  xx = (2.0 * i - 1.0) * subd1_1 - 1.0;   /* sub-element centre */
	  radx = (xx-eta1)*(xx-eta1);

	  for(j=1;j<=SUBD;j++) {
	    zz = (2.0 * j - 1.0) * subd1_1 - 1.0; 
	    radz = (zz-eta2)*(zz-eta2);

	    for(k=1;k<=SUBD;k++)   {
	      yy = (2.0 * k - 1.0) * subd1_1 - 1.0; 
	      
	      rad = radx + radz + (yy-eta3)*(yy-eta3);
	      
	      if(rad < subelt_dist3[i][j][k]) {
		subelt_dist3[i][j][k] = rad;
		subelt_tracer3[i][j][k] = m;
	      }
	    }
	  }
	}
      }
    
      for(i=1;i<=SUBD;i++)
	for(j=1;j<=SUBD;j++) {
	  for(k=1;k<=SUBD;k++)
	    if(subelt_tracer3[i][j][k]!=0)
	      E->tracer.tracer_weight_fn[E->mesh.levmax][subelt_tracer3[i][j][k]] += 
		subd3_1;
	}
    }
  }
  else {                                 /* TWO DIMENSIONS */
    for(el=1;el<=E->mesh.nel;el++) {
      for(i=1;i<=SUBD;i++)
	for(j=1;j<=SUBD;j++) {
	  subelt_dist2[i][j] = 1.0e32;
	  subelt_tracer2[i][j]=0;
	}
    
      for(tr=0;tr<E->tracer.tr_in_element_number[E->mesh.levmax][el];tr++) {
      
	q = E->tracer.tr_in_element_offset[E->mesh.levmax][el] + tr;
	m = E->tracer.tr_in_element[E->mesh.levmax][q]; 
 
	eta1 = ETA1[m];
	eta2 = ETA2[m];
      
	E->tracer.tracer_weight_fn[E->mesh.levmax][m] = 0.0;

	for(i=1;i<=SUBD;i++) {
	  xx = (2.0 * i - 1.0) * subd1_1 - 1.0;   /* sub-element centre */
	  radx = (xx-eta1)*(xx-eta1); 

	  for(j=1;j<=SUBD;j++) {
	    zz = (2.0 * j - 1.0) * subd1_1 - 1.0; 
	    
	    rad = radx + (zz-eta2)*(zz-eta2);
	    
	    if(rad < subelt_dist2[i][j]) {
	      subelt_dist2[i][j] = rad;
	      subelt_tracer2[i][j] = m;
	    }
	  }
	}
      }
      
      for(i=1;i<=SUBD;i++)
	for(j=1;j<=SUBD;j++) {
	  /* fprintf(stderr,"el %d, sub %d,%d -> %g to tracer %d\n",
	     el,i,j,4.0 * subd2_1,subelt_tracer2[i][j]); */
	  if(subelt_tracer2[i][j]!=0)
	    E->tracer.tracer_weight_fn[E->mesh.levmax][subelt_tracer2[i][j]] += 
	      subd2_1;
	}
    } /* Next element */    
  }
 
#else

 for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
   E->tracer.tracer_weight_fn[E->mesh.levmax][m] = E->tracer.tracer_std_weighting[E->mesh.levmax][m] ;
	/* fprintf(stderr,"Tracer %d; weight %g\n",m,E->tracer.tracer_weight_fn[E->mesh.levmax][m] ); */
}
#endif

 /* Relative weights cascade to all levels, but will need rescaling
     for the coarser meshes since the number of tracers increases
     while the size of the master domain is the same */

 if(3==E->mesh.nsd)
    for(level=E->mesh.levmax;level>E->mesh.levmin;level--) 
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	E->tracer.tracer_weight_fn[level-1][m] = 0.125 * E->tracer.tracer_weight_fn[level][m];
      }
  else
   for(level=E->mesh.levmax;level>E->mesh.levmin;level--) 
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	E->tracer.tracer_weight_fn[level-1][m] =  0.25 * E->tracer.tracer_weight_fn[level][m];
      }
  return;
}

void tracer_constrain_weightings (
				  struct All_variables *E,
				  standard_precision **ETA1,
				  standard_precision **ETA2,
				  standard_precision **ETA3
				  )
{
  int i,j,k;
  int ii;
  int level;
  int el;
  int m,tr,q;
  
  standard_precision t1,t2,t3,t4,t5,t6,t7;
  standard_precision tt1,tt2,tt3,tt4,tt5,tt6,tt7;
  standard_precision W1,W2,W3,W4,W5,W6,W7;
  standard_precision master_volume;
  standard_precision eta1,eta2,eta3;

  const int dims = E->mesh.nsd;

 if(3==dims) 
    master_volume = 8.0;
  else
    master_volume = 4.0;

  /* First, and for all levels, we ensure that the 
     weights combine to give the correct element volume */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
    for(el=1;el<=E->mesh.NEL[level];el++) {
      t1=0.0;
      for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	  q = E->tracer.tr_in_element_offset[level][el] + tr;
	  m = E->tracer.tr_in_element[level][q]; 
   	  t1 += E->tracer.tracer_weight_fn[level][m];
      }
      if(t1 != 0.0) {
	for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	  q = E->tracer.tr_in_element_offset[level][el] + tr;
	  m = E->tracer.tr_in_element[level][q];  
	  E->tracer.tracer_weight_fn[level][m] *= master_volume / t1;
	}
      }
    }
  }
  
  /* Now apply the various constraints (as best 
     we can - it is not possible exactly) one at a time
     and iteratively. Ideally this should be done at 
     all levels, but if done only at the highest, it
     saves an enormous amount of time (we can use
     the existing ETA array).
  */    



#if 1
  if(2==dims) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	
      for(el=1;el<=E->mesh.NEL[level];el++) {
	for(ii=1;ii<=100;ii++) {
	  t1 = t2 = t3 = 0.0;
	  W1 = W2 = W3 = 0.0;
	
	  for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	    q = E->tracer.tr_in_element_offset[level][el] + tr;
	    m = E->tracer.tr_in_element[level][q]; 
	  
	    eta1 = ETA1[level][m];
	    eta2 = ETA2[level][m];
	  
	    W1 +=  E->tracer.tracer_weight_fn[level][m] * eta1;
	    W2 +=  E->tracer.tracer_weight_fn[level][m] * eta2;
	    W3 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2;

	    t1 += eta1 * eta1;
	    t2 += eta2 * eta2;
	    t3 += eta1 * eta1 * eta2 * eta2;
	  }
	
	  if(t1 != 0.0)
	    tt1 = - W1 / t1;
	  if(t2 != 0.0)
	    tt2 = - W2 / t2;
	  if(t3 != 0.0)
	    tt3 = - W3 / t3;

	  t1 = 0.0;
	  W1 = W2 = W3 = 0.0;

	  for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	    q = E->tracer.tr_in_element_offset[level][el] + tr;
	    m = E->tracer.tr_in_element[level][q]; 
     
	    eta1 = ETA1[level][m];
	    eta2 = ETA2[level][m];
	
	    E->tracer.tracer_weight_fn[level][m] += tt1 * eta1 + tt2 * eta2 + tt3 * eta1 * eta2;
	  
	    if(E->tracer.tracer_weight_fn[level][m] < 0.0)
	      E->tracer.tracer_weight_fn[level][m] = 0.0;

	    W1 +=  E->tracer.tracer_weight_fn[level][m] * eta1;
	    W2 +=  E->tracer.tracer_weight_fn[level][m] * eta2;
	    W3 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2;

	    t1 += E->tracer.tracer_weight_fn[level][m];
	  }
    
	  /* And rescale the weights to maintain element volume */

	  if(t1 != 0.0)
	    for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	      q = E->tracer.tr_in_element_offset[level][el] + tr;
	      m = E->tracer.tr_in_element[level][q];
      
	      E->tracer.tracer_weight_fn[level][m] *= 4.0 / t1;
	   
	   	 /*  fprintf(stderr,"Tracer %d/%d reweighed -> %g (%g,%g)\n",m,
	   	  	level,E->tracer.tracer_weight_fn[level][m],eta1,eta2); */

	    }
      
	  /* The total of W1,W2,W3 should be an indicator of the 
	     residual from this process ... */

	  if((fabs(W1)+fabs(W2)+fabs(W3)) < 0.1 * E->control.accuracy * t1)
	    break; 
	}
      }
    }
  }
  else /* (3==dims) */ {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	
      for(el=1;el<=E->mesh.NEL[level];el++) {	
	for(ii=1;ii<=100;ii++) {							
	  
	  t1 = t2 = t3 = t4 = t5 = t6 = t7 = 0.0;
	  W1 = W2 = W3 = W4 = W5 = W6 = W7 = 0.0;
	
	  for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {     
	    q = E->tracer.tr_in_element_offset[level][el] + tr;
	    m = E->tracer.tr_in_element[level][q]; 
	  
	    eta1 = ETA1[level][m];
	    eta2 = ETA2[level][m];
	    eta3 = ETA3[level][m];
	  
	    W1 +=  E->tracer.tracer_weight_fn[level][m] * eta1;
	    W2 +=  E->tracer.tracer_weight_fn[level][m] * eta2;
	    W3 +=  E->tracer.tracer_weight_fn[level][m] * eta3;
	    W4 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2;
	    W5 +=  E->tracer.tracer_weight_fn[level][m] * eta2 * eta3;
	    W6 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta3;
	    W7 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2 * eta3;
	 
	    t1 += eta1 * eta1;
	    t2 += eta2 * eta2;
	    t3 += eta3 * eta3;
	    t4 += eta1 * eta2 * eta1 * eta2;
	    t5 += eta2 * eta3 * eta2 * eta3;
	    t6 += eta1 * eta3 * eta1 * eta3;
	    t7 += eta1 * eta1 * eta2 * eta2 * eta3 * eta3;
	  }

	  if(t1 != 0.0)
	    tt1 = - W1 / t1;
	  if(t2 != 0.0)
	    tt2 = - W2 / t2;
	  if(t3 != 0.0)
	    tt3 = - W3 / t3;
	  if(t4 != 0.0)
	    tt4 = - W4 / t4;
	  if(t5 != 0.0)
	    tt5 = - W5 / t5;
	  if(t6 != 0.0)
	    tt6 = - W6 / t6;
	  if(t7 != 0.0)
	    tt7 = - W7 / t7;

	  t1 = t2 = t3 = t4 = t5 = t6 = t7 = 0.0;
	  W1 = W2 = W3 = W4 = W5 = W6 = W7 = 0.0;

	  for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {
	    q = E->tracer.tr_in_element_offset[level][el] + tr;
	    m = E->tracer.tr_in_element[level][q]; 
     
	    eta1 = ETA1[level][m];
	    eta2 = ETA2[level][m];
	    eta3 = ETA3[level][m];

	    E->tracer.tracer_weight_fn[level][m] += 
	      tt1 * eta1 +
	      tt2 * eta2 +
	      tt3 * eta3 +
	      tt4 * eta1 * eta2 +
	      tt5 * eta2 * eta3 +
	      tt6 * eta1 * eta3 +
	      tt7 * eta1 * eta2 * eta3;
	 
	    if(E->tracer.tracer_weight_fn[level][m] < 0.0)
	      E->tracer.tracer_weight_fn[level][m] = 0.0;

	    W1 +=  E->tracer.tracer_weight_fn[level][m] * eta1;
	    W2 +=  E->tracer.tracer_weight_fn[level][m] * eta2;
	    W3 +=  E->tracer.tracer_weight_fn[level][m] * eta3;
	    W4 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2;
	    W5 +=  E->tracer.tracer_weight_fn[level][m] * eta2 * eta3;
	    W6 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta3;
	    W7 +=  E->tracer.tracer_weight_fn[level][m] * eta1 * eta2 * eta3;
	
	    t1 += E->tracer.tracer_weight_fn[level][m];
	  }

	  /* And rescale the weights once more to maintain element volume */

	  if(t1 != 0.0)
	    for(tr=0;tr<E->tracer.tr_in_element_number[level][el];tr++) {   
	      q = E->tracer.tr_in_element_offset[level][el] + tr;
	      m = E->tracer.tr_in_element[level][q];
      
	      E->tracer.tracer_weight_fn[level][m] *= master_volume / t1;
	    }
      
	  /* The total of W1,W2 ... W7 should be an indicator of the 
	     residual from this process ... */

	  if((fabs(W1)+fabs(W2)+fabs(W3)+fabs(W4)+fabs(W5)+fabs(W6)+fabs(W7)) < 
	     0.2 * E->control.accuracy * t1)
	    break; 
	}
      }
    }
  }
#endif
}
