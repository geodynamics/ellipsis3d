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
#include "element_definitions.h"
#include "global_defs.h"
#include <sys/time.h>
#include <sys/resource.h>


/* =================================================================
   Assemble a grad_P vector element by element from particle locations
   =================================================================  */

void assemble_grad_qst3(
			struct All_variables *E,
			standard_precision *Q,
			standard_precision *gradQ1,
			standard_precision *gradQ2,
			standard_precision *gradQ3,
			int level
			)
{
  int e,i,j,p,a,node,el,m;
  void s_strip_bcs_from_residual();
 
  const int nel=E->mesh.NEL[level];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int trcs=E->tracer.NUM_TRACERS;


  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight;
  standard_precision g1,g2,g3;

  higher_precision *node_R;

/*  if(3==dims)*/
    for(i=1;i<=E->mesh.NNO[level];i++)
      gradQ1[i] = gradQ2[i] = gradQ3[i] = 0.0;
/*  else
    for(i=1;i<=E->mesh.NNO[level];i++)
      gradQ1[i] = gradQ2[i] = 0.0;*/
  
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(E->tracer.property_group[m] < 0)
      continue;

    el = E->tracer.tracer_elt[level][m];

    eta1=eta2=eta3=0.0;
    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);  
    dOmega = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];

    if(dims==2) {
      for(a=1;a<=ENODES2D;a++) {
	node = E->IEN[level][el].node[a];
	gradQ1[node] += lNx[1][a] * Q[m] * dOmega;
	gradQ2[node] += lNx[2][a] * Q[m] * dOmega;
      }
    }
    else {
      for(a=1;a<=ENODES3D;a++) {
	node = E->IEN[level][el].node[a];
	gradQ1[node] += lNx[1][a] * Q[m] * dOmega;
	gradQ2[node] += lNx[2][a] * Q[m] * dOmega;
	gradQ3[node] += lNx[3][a] * Q[m] * dOmega;
      }
    }
  }
  /* Now may need to rotate out of cartesian domain 
     to accomodate any skewed bc's */

  if(E->control.HAVE_SKEWBCS) {
    for(i=1;i<=E->mesh.NNO[level];i++)
      if(E->NODE[level][i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[level][i];
	if(dims==2) {
	  g1 = node_R[0 + dims * 0] * gradQ1[i] +  node_R[0 + dims * 1] * gradQ2[i];
	  g2 = node_R[1 + dims * 0] * gradQ1[i] +  node_R[1 + dims * 1] * gradQ2[i];
	  gradQ1[i] = g1;
	  gradQ2[i] = g2;
	}
	else {
	  g1 = node_R[0 + dims * 0] * gradQ1[i] + node_R[0 + dims * 1] * gradQ2[i] + node_R[0 + dims * 2] * gradQ3[i];
	  g2 = node_R[1 + dims * 0] * gradQ1[i] + node_R[1 + dims * 1] * gradQ2[i] + node_R[1 + dims * 2] * gradQ3[i];
	  g3 = node_R[2 + dims * 0] * gradQ1[i] + node_R[2 + dims * 1] * gradQ2[i] + node_R[2 + dims * 2] * gradQ3[i];
	  gradQ1[i] = g1;
	  gradQ2[i] = g2;
	  gradQ3[i] = g3;
	}
      }
  }
  strip_bcs_from_residual_6(E,gradQ1,gradQ2,gradQ3,NULL,NULL,NULL,level);
  return;
}

/* =================================================================
   Assemble a grad_P vector element by element from particle locations
   =================================================================  */

void assemble_Mq_t(
		   struct All_variables *E,
		   standard_precision *Q,
		   standard_precision *MQ,
		   int level
		   )
{
  int e,i,j,p,a,el,m;

  const int nel=E->mesh.NEL[level];
  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;

  standard_precision dOmega,bulk_visc;

  /* This should be pressure shape function * P / bulk_visc but
     simplifies for constant pressure elements  */

  for(i=1;i<=nel;i++)
    MQ[i] = 0.0;

#if 0
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    /* Allow for perfectly incompressible material */
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio < 1.0)
      continue;
    
    el = E->tracer.tracer_elt[level][m]; 
    /* constant pressure ... no shape function calculation needed 
       for this special case !! */
    
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0) {  /* Otherwise assume incompressible */
      if(dims==2 && dofs==2)
	bulk_visc = E->tracer.BulkVisc[m] - E->tracer.Visc[m];
      if(dims==3 && dofs==3)
	bulk_visc = E->tracer.BulkVisc[m] - 0.6666666666 * E->tracer.Visc[m] ;
    }
    else
      bulk_visc = 0.0;

    dOmega = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];

    if(bulk_visc > 0.0)
      MQ[el] +=  Q[el] * dOmega / bulk_visc;
  }
#else
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    /* Allow for perfectly incompressible material */
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio < (2.0/dims))
      continue;
    
    el = E->tracer.tracer_elt[level][m]; 
    /* constant pressure ... no shape function calculation needed 
       for this special case !! */
    bulk_visc = E->tracer.BulkVisc[m] - (2.0/dims)*E->tracer.Visc[m];
    dOmega = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];
    MQ[el] +=  Q[el] * dOmega / bulk_visc;
  }
#endif
  return;
}

void add_tracers_to_global_f(
			     struct All_variables *E,
			     standard_precision *glob_F,
			     int k,
			     int level
			     )
{
  int a,b,p,m,el,qd,q,node,i,tr;
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,volume;

  standard_precision *ftrmag;
  standard_precision *dF1,dF2[4];
  standard_precision magnitude;

  higher_precision *node_R;
 
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const int nno=E->mesh.NNO[level];
  const int nel=E->mesh.NEL[level];
  const int neq=E->mesh.NEQ[level];
  const int trcs = E->tracer.NUM_TRACERS;

  standard_precision nodelx,nodelz,nodely;
  standard_precision trdelx,trdelz,trdely,trdelr;
  standard_precision expans;
  standard_precision sintheta,costheta;

  struct TRACER_ELT_WEIGHT *lN;
  struct IEN *IEN1;
  struct ID *ID1;
  int *Phases;
  int *Current;
  int *Equil;
  int *Propgrp;
  int *offset;
  int *in_elt;
  standard_precision *Phase_fn;
  standard_precision *Density;
  standard_precision *Temp;
  standard_precision *Therm_exp;
  standard_precision *Weight;
  standard_precision *Depl;    /*RAA: 24/09/02, C. O'Neill melting stuff*/  
  standard_precision *Depl_exp;  /*RAA: 24/09/02, C. O'Neill melting stuff*/  

  void Vis_Elas_RHS() ;
  void Vis_Elas_LD_RHS() ;
  /* Pointer aliases to allow optimizer, vectorizer to 
     see the locally unchanging, unaliased (ironically)
     nature of the various tracer arrays */

  Phases = E->tracer.Phases;
  Current = E->tracer.Current_phase;
  Equil = E->tracer.Equilibrium_phase;
  Phase_fn = E->tracer.phase_function;
  Temp = E->tracer.T;
  Depl = E->tracer.depl;  /*RAA: 24/09/02, C. O'Neill melting stuff*/
  Depl_exp = E->tracer.Depl_exp;    /*RAA: 24/09/02, C. O'Neill melting stuff*/
  Therm_exp = E->tracer.Therm_exp;
  Density = E->tracer.Density;
  Propgrp = E->tracer.property_group;
  Weight = E->tracer.weight;
  IEN1 = E->IEN[level];
  ID1 = E->ID[level];
  lN = E->tracer.sfn_values[level];
  offset = E->tracer.tr_in_element_offset[level];
  in_elt = E->tracer.tr_in_element[level];
 
  dF1 = (standard_precision *) Malloc0((E->mesh.NEQ[level] + 1) * sizeof(standard_precision));  
  ftrmag = (standard_precision *) Malloc0((trcs + 1) * sizeof(standard_precision));  

  if(E->control.verbose)
    fprintf(stderr,"Compute Global F from tracers\n");

  /* tracer_inertial_terms(E,u1,u2,u3,k,level); */

  for(p=0;p<neq;p++) 
    dF1[p] = 0.0;

  /* 1. Materials where there is no phase change possible */
  
  for(m=1;m<=trcs;m++) { 
    if(Phases[Propgrp[m]] == 1) {
      /* expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m]);   */
      expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m] - Depl_exp[Propgrp[m]]*Depl[m]); /*O'Neill: melting stuff*/ 
      Weight[m] = Density[Propgrp[m]*MAX_MATERIAL_PHASES+0] * expans ;
    }
  }
  /* 2. Materials which are in a metastable state (rare case) */

  for(m=1;m<=trcs;m++) { 
    if(Phases[Propgrp[m]] != 1 && 
       (Current[m]  !=  Equil[m])) {
      /* expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m]);   */
      expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m] - Depl_exp[Propgrp[m]]*Depl[m]);   /*O'Neill: melting stuff*/ 
      Weight[m] = Density[Propgrp[m]*MAX_MATERIAL_PHASES+Current[m]] * expans ; 
    }
  }

  /* 3. Materials in equilibrium phase, phase function
     is less than zero - need to smooth to phase boundary */

 for(m=1;m<=trcs;m++) { 
    if(Phases[Propgrp[m]] != 1 && 
       (Current[m]  ==  Equil[m]) &&
       (Phase_fn[m] > 0.0)) {
      /* expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m]);   */
      expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m] - Depl_exp[Propgrp[m]]*Depl[m]);   /*O'Neill: melting stuff*/ 
      Weight[m] = 
	(Phase_fn[m] * 
	 Density[Propgrp[m]*MAX_MATERIAL_PHASES+Current[m]] +
	 (1.0 - Phase_fn[m]) * 
	 Density[Propgrp[m]*MAX_MATERIAL_PHASES+Current[m]-1]) * expans ; 
    }
 }

 /* 4. As 3, but other side of phase boundary. */

 for(m=1;m<=trcs;m++) { 
    if(Phases[Propgrp[m]] != 1 && 
	 (Current[m]  ==  Equil[m]) &&
	 (Phase_fn[m] <= 0.0)) { 
      /* expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m]);   */
      expans = (1.0 - Therm_exp[Propgrp[m]] * Temp[m] - Depl_exp[Propgrp[m]]*Depl[m]);  /*O'Neill: melting stuff*/ 
      Weight[m] = 
	(-Phase_fn[m] * 
	 Density[Propgrp[m]*MAX_MATERIAL_PHASES+Current[m]] +
	 (1.0+Phase_fn[m]) * Density[Propgrp[m]*MAX_MATERIAL_PHASES+Current[m]+1])
	    * expans ;
    }
 }

 sintheta = sin(E->data.grav_theta * M_PI); 
 costheta = cos(E->data.grav_theta * M_PI);

 for(m=1;m<=trcs;m++) {
   ftrmag[m] = 
     E->tracer.tracer_weight_fn[level][m] *
     E->tracer.tracer_jacobian[level][m] * 
     E->data.grav_acc * E->tracer.weight[m] * E->tracer.Density_change[m] ;
     /* fprintf(stderr,"Tracer %d: weight = %g, jacobian = %g\n",m,
     	E->tracer.tracer_weight_fn[level][m],E->tracer.tracer_jacobian[level][m]); */
 }

 if(!E->control.AXI) {
   for(el=1;el<=nel;el++){
     for(a=1;a<=ends;a++) {
       node=IEN1[el].node[a];    

       for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
	 m = in_elt[i+offset[el]];
	 
	 /* translation degrees of freedom */
	 
	 /* fprintf(stderr,"Tracer %d contributes %g (%g,%g,%g) to F1[%d] \n",
	 m,lN[m].node[a] * ftrmag[m] * sintheta,lN[m].node[a],ftrmag[m] , sintheta,
		ID1[node].doff[1]); */
	
	 
	 dF1[ID1[node].doff[1]] +=  lN[m].node[a] * ftrmag[m] * sintheta;
	 dF1[ID1[node].doff[2]] +=  lN[m].node[a] * ftrmag[m] * costheta;
	 /* if Non-gravitational force terms, or gravity pointing at another
	    funny angle, then need to add a third component if 3==dims */
	  if(3==dims) /*RAA */
            dF1[ID1[node].doff[3]] +=  0.0;
       }
     }
   }
 } 
 else { /*AXI*/
   for(el=1;el<=nel;el++) {
     for(a=1;a<=ends;a++) {
       node=E->IEN[level][el].node[a];
       
       for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
	 m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
	
	 magnitude = E->tracer.sfn_values[level][m].node[a] * 2 * M_PI *  E->tracer.tx[m] *
	   E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];

           dF1[E->ID[level][node].doff[1]] += 
	  magnitude * E->data.grav_acc * E->tracer.weight[m] * ( sin(E->data.grav_theta * M_PI));
	   dF1[E->ID[level][node].doff[2]] += 
	  magnitude * E->data.grav_acc * E->tracer.weight[m] * ( cos(E->data.grav_theta * M_PI));

       }
     }
   }
 }


 if(E->control.ELASTICITY) {
  Vis_Elas_RHS(E,dF1,level) ;
  /*  if(E->control.deformation)
     Vis_Elas_LD_RHS(E,dF1,level) ;*/ /**/
 }


  /* Transfer to global array, change coordinates for
     skewed nodes */

  if(E->control.HAVE_SKEWBCS) 
    for(p=1;p<=E->mesh.NNO[level];p++) {
      if(E->NODE[level][p] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[level][p];
	dF2[1] = 
	  node_R[0 + dims * 0] * dF1[E->ID[level][p].doff[1]] +
	  node_R[0 + dims * 1] * dF1[E->ID[level][p].doff[2]];  
	dF2[2] = 
	  node_R[1 + dims * 0] * dF1[E->ID[level][p].doff[1]] + 
	  node_R[1 + dims * 1] * dF1[E->ID[level][p].doff[2]];  
	dF1[E->ID[level][p].doff[1]] = dF2[1]; 
	dF1[E->ID[level][p].doff[2]] = dF2[2];
      }
    }
  
#pragma loop novrec glob_F,dF1
  for(p=0;p<E->mesh.NEQ[level];p++) {
    glob_F[p] += dF1[p];
    /* fprintf(stderr,"%d (%dx%d): F[%d] = %g (%g)\n",level,
       E->mesh.ELX[level],E->mesh.ELZ[level],p,glob_F[p],dF1[p]); */
  } 
  
  if(E->control.verbose)
    fprintf(stderr,"Computed Global F from tracers\n");
  
  free((void *) dF1);
  free((void *) ftrmag);
  return;
}

void Vis_Elas_RHS(
		  struct All_variables *E,
		  standard_precision *F,
		  int level
		  )
{
  standard_precision eta1,eta2,eta3,dOmega ;
  standard_precision Shearscale,Bulkscale;
  int m,node,a,el,i,jj ;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision elastime1;

  const int dims=E->mesh.nsd ;
  const int dofs=E->mesh.dof ;
  const int nel=E->mesh.NEL[level] ;
  const int ends=enodes[dims];

  if(!E->control.ELASTICITY)
    return;

  jj=0;
  for(el=1;el<=nel;el++){
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
  
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
      if(3==dims)
	eta3 = E->tracer.eta3[level][m];
      
      get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);
      v_shape_fn(E,el,lN,eta1,eta2,eta3,level);
      dOmega  = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];
      
      Shearscale = Bulkscale = 0.0 ;       
      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	elastime1 =  E->advection.elastic_timestep * E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
      	Shearscale = E->tracer.Visc[m]  / elastime1; 
      }
     
      if(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod > 0.0) {
	elastime1 =  E->advection.elastic_timestep * E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod;
	Bulkscale = E->tracer.BulkVisc[m]  / elastime1; /* Need to be fixed */
      }
       	

      for(a=1;a<=ends;a++) {
	node=E->IEN[level][el].node[a];
	
	if(2==dofs && 2==dims) {
	  F[E->ID[level][node].doff[1]] -= (Shearscale * (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m]) - 
					    Bulkscale * lNx[1][a]*E->tracer.Pt[m] ) * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (Shearscale * (lNx[1][a]*E->tracer.S12[m] + lNx[2][a]*E->tracer.S22[m])- 
					    Bulkscale * lNx[2][a]*E->tracer.Pt[m] )  * dOmega ;
 	}
	else if(3==dofs && 2==dims) { /* ???????? */
	  F[E->ID[level][node].doff[1]] -= (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m]) * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (lNx[2][a]*E->tracer.S22[m] + lNx[1][a]*E->tracer.S21[m]) * dOmega ;
	  F[E->ID[level][node].doff[3]] -= (lN[a]*E->tracer.S12[m] - lN[a]*E->tracer.S21[m] 
	    + lNx[1][a]*E->tracer.M31[m] + lNx[2][a]*E->tracer.M32[m]) * dOmega ;
	}
	else if(3==dims && 3==dofs) {
	  F[E->ID[level][node].doff[1]] -= (Shearscale * (lNx[1][a]*E->tracer.S11[m] +
							  lNx[2][a]*E->tracer.S12[m] + 
							  lNx[3][a]*E->tracer.S13[m]) - 
					    Bulkscale * lNx[1][a]*E->tracer.Pt[m] ) * dOmega ;

	  F[E->ID[level][node].doff[2]] -= (Shearscale * (lNx[1][a]*E->tracer.S21[m] +
							  lNx[2][a]*E->tracer.S22[m] + 
							  lNx[3][a]*E->tracer.S23[m]) - 
					    Bulkscale * lNx[2][a]*E->tracer.Pt[m] ) * dOmega ;

	  F[E->ID[level][node].doff[3]] -= (Shearscale * (lNx[1][a]*E->tracer.S31[m] +
							  lNx[2][a]*E->tracer.S32[m] + 
							  lNx[3][a]*E->tracer.S33[m]) - 
					    Bulkscale * lNx[3][a]*E->tracer.Pt[m] ) * dOmega ;
	}
	else /*if(6==dofs)*/ { /* ?????? */
	  F[E->ID[level][node].doff[1]] -= (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m] + lNx[3][a]*E->tracer.S13[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (lNx[2][a]*E->tracer.S22[m] + lNx[1][a]*E->tracer.S21[m] + lNx[3][a]*E->tracer.S23[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[3]] -= (lNx[3][a]*E->tracer.S33[m] + lNx[1][a]*E->tracer.S31[m] + lNx[2][a]*E->tracer.S32[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[4]] -= (lN[a]*(E->tracer.S23[m]-E->tracer.S32[m]) + 
	    lNx[2][a]*E->tracer.M12[m] + lNx[3][a]*E->tracer.M13[m]) * dOmega ;
	  F[E->ID[level][node].doff[5]] -= (lN[a]*(E->tracer.S13[m]-E->tracer.S31[m]) + 
	    lNx[1][a]*E->tracer.M21[m] + lNx[3][a]*E->tracer.M23[m]) * dOmega ;
	  F[E->ID[level][node].doff[6]] -= (lN[a]*(E->tracer.S12[m]-E->tracer.S21[m]) + 
	    lNx[1][a]*E->tracer.M31[m] + lNx[2][a]*E->tracer.M32[m]) * dOmega ;
	}
      }
    }
  }
return ;
}

void Vis_Elas_LD_RHS(
		     struct All_variables *E,
		     standard_precision *F,
		     int level
		     )
{
  standard_precision eta1,eta2,eta3,dOmega,omega,omegamax,V12,V21 ;
  int m,node,a,el,i,j ;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];

  const int dims=E->mesh.nsd ;
  const int dofs=E->mesh.dof ;
  const int nel=E->mesh.NEL[level] ;
  const int ends=enodes[dims];    
  
  if(!E->control.ELASTICITY)
    return;
  
  omegamax = 0.0 ;

  if(E->control.verbose)
    fprintf(stderr,"Updating stresses in Large deformation \n") ;



  /* Shouldn't it respect rotated BC's - presumably
     the rotation rate doesn't change though ?? */


  for(el=1;el<=nel;el++){      
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
      if(3==dims)
	eta3 = E->tracer.eta3[level][m];

      eta1=eta2=0.0; /**/
      if(3==dims) /*RAA - addeded these 2 lines, but should verify*/
	eta3 = 0.0;

      get_global_v_x_shape_fn(E,E->tracer.tracer_elt[level][m],lNx,
			      &dOmega,eta1,eta2,eta3,level);
     
      V12 = V21 =0.0 ;
      for(j=1;j<=ends;j++) {
	V12 +=  E->VV[level][1][E->IEN[level][el].node[j]] * lNx[2][j]; 
	V21 +=  E->VV[level][2][E->IEN[level][el].node[j]] * lNx[1][j]; 
	}
      
      omega = 0.5 * (V12 - V21);

      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	dOmega  = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * 
	  E->tracer.Visc[m] / E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;	
      }
      else
	continue;

      for(a=1;a<=ends;a++) {
	node=E->IEN[level][el].node[a];

	if(2==dofs && 2==dims) {
	  F[E->ID[level][node].doff[1]] -= 
	    ( 2.0 * lNx[1][a]*E->tracer.S12[m] + lNx[2][a]*(E->tracer.S22[m]-E->tracer.S11[m])) * omega * dOmega ;
	  F[E->ID[level][node].doff[2]] -= 
	    (-2.0 * lNx[2][a]*E->tracer.S12[m] + lNx[1][a]*(E->tracer.S22[m]-E->tracer.S11[m])) * omega * dOmega ;
 	}
	else if(3==dofs && 2==dims) {
	  F[E->ID[level][node].doff[1]] -= (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m]) * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (lNx[2][a]*E->tracer.S22[m] + lNx[1][a]*E->tracer.S21[m]) * dOmega ;
	  F[E->ID[level][node].doff[3]] -= (lN[a]*E->tracer.S12[m] - lN[a]*E->tracer.S21[m] 
	    + lNx[1][a]*E->tracer.M31[m] + lNx[2][a]*E->tracer.M32[m]) * dOmega ;
	}
	else if(3==dims && 3==dofs) {  /* ????? */
	  F[E->ID[level][node].doff[1]] -= (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m] + lNx[3][a]*E->tracer.S13[m])
	    * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (lNx[2][a]*E->tracer.S22[m] + lNx[3][a]*E->tracer.S23[m] + lNx[1][a]*E->tracer.S12[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[3]] -= (lNx[3][a]*E->tracer.S33[m] + lNx[2][a]*E->tracer.S23[m] + lNx[1][a]*E->tracer.S13[m]) 
	    * dOmega ;
	}
	else /*if(6==dofs)*/ {         /* NOT FINISHED */
	  F[E->ID[level][node].doff[1]] -= (lNx[1][a]*E->tracer.S11[m] + lNx[2][a]*E->tracer.S12[m] + lNx[3][a]*E->tracer.S13[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[2]] -= (lNx[2][a]*E->tracer.S22[m] + lNx[1][a]*E->tracer.S21[m] + lNx[3][a]*E->tracer.S23[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[3]] -= (lNx[3][a]*E->tracer.S33[m] + lNx[1][a]*E->tracer.S31[m] + lNx[2][a]*E->tracer.S32[m]) 
	    * dOmega ;
	  F[E->ID[level][node].doff[4]] -= (lN[a]*(E->tracer.S23[m]-E->tracer.S32[m]) + 
	    lNx[2][a]*E->tracer.M12[m] + lNx[3][a]*E->tracer.M13[m]) * dOmega ;
	  F[E->ID[level][node].doff[5]] -= (lN[a]*(E->tracer.S13[m]-E->tracer.S31[m]) + 
	    lNx[1][a]*E->tracer.M21[m] + lNx[3][a]*E->tracer.M23[m]) * dOmega ;
	  F[E->ID[level][node].doff[6]] -= (lN[a]*(E->tracer.S12[m]-E->tracer.S21[m]) + 
	    lNx[1][a]*E->tracer.M31[m] + lNx[2][a]*E->tracer.M32[m]) * dOmega ;
	}
      }
    }
  }
  if(E->control.verbose)
    fprintf(stderr,"Stresses in Large deformation updated \n") ;
return ;
}
