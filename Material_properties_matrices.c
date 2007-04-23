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

void construct_Bi(
  struct All_variables *E ,
  standard_precision lNx[4][ELNMAX+1] ,
  int i ,
  standard_precision BiV[16][7] ,
  standard_precision one_over_x ,
  standard_precision lN[ELNMAX+1]
)
{
  const int dims=E->mesh.nsd ;
  int k,j ;

  for(k=1;k<16;k++)
    for(j=1;j<7;j++)
      BiV[k][j] = 0.0 ;

  if(2==dims) { 
    BiV[1][1] = lNx[1][i] ; 
    BiV[2][2] = lNx[2][i] ;
    if(E->control.model==1) {
      BiV[3][1] = lNx[2][i] 	; BiV[3][2] = lNx[1][i] ;
    }
    if(E->control.model==2) {
      BiV[3][1] = lNx[2][i] ;                         BiV[3][3] = lN[i] ;
      BiV[4][2] = lNx[1][i] ; BiV[4][3] = -lN[i] ;
      BiV[5][3] = lNx[1][i] ;
      BiV[6][3] = lNx[2][i] ;
    }
    if(E->control.AXI) {
      BiV[4][1] =lN[i]*one_over_x ;
    }
  }
  else {
    BiV[1][1] = lNx[1][i] ;
    BiV[2][2] = lNx[2][i] ;
    BiV[3][3] = lNx[3][i] ;
    if(E->control.model==1) {
      BiV[4][2] = lNx[3][i] ; BiV[4][3] = lNx[2][i] ;
      BiV[5][1] = lNx[3][i] ;                         BiV[5][3] = lNx[1][i] ;
      BiV[6][1] = lNx[2][i] ; BiV[6][2] = lNx[1][i] ;
    }
    else {
      BiV[4][1] = lNx[2][i] ; BiV[4][6] = lN[i] ;
      BiV[5][2] = lNx[1][i] ; BiV[5][6] = -lN[i] ;
      BiV[6][1] = lNx[3][i] ; BiV[6][5] = lN[i] ;
      BiV[7][3] = lNx[1][i] ; BiV[7][5] = -lN[i] ;
      BiV[8][2] = lNx[3][i] ; BiV[8][4] = lN[i] ;
      BiV[9][3] = lNx[2][i] ; BiV[9][4] = -lN[i] ;
      BiV[10][4] = lNx[2][i] ;
      BiV[11][5] = lNx[1][i] ;
      BiV[12][4] = lNx[3][i] ;
      BiV[13][6] = lNx[1][i] ;
      BiV[14][5] = lNx[3][i] ;
      BiV[15][6] = lNx[2][i] ;
    }
  }
return ;
}

void construct_D_Vel(
		     struct All_variables *E,
		     standard_precision DV[16][16],
		     standard_precision dOmega,
		     standard_precision penalty,
		     int m,
		     int level
		     )
{
  int i,j ;
  const int dims=E->mesh.nsd ;

  for(i=1;i<16;i++)
    for(j=1;j<16;j++)
      DV[i][j] = 0.0 ;

  if(E->control.model==2) {
    DV[1][1] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A11 - dOmega * penalty ;
    DV[2][2] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A22 - dOmega * penalty ;
    DV[3][3] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G11 ;
    DV[4][4] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G22 ;
    DV[3][4] = DV[4][3] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G12 ;
    DV[5][5] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B31 ;
    DV[6][6] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B32 ;
    if(3==dims) { /* must be checked */
      DV[1][1] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].A11 - dOmega * penalty ;
      DV[2][2] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].A22 - dOmega * penalty ;
      DV[3][3] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].A33 - dOmega * penalty ;
      DV[3][3] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].G11 ;
      DV[4][4] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].G22 ;
      DV[3][4] = DV[4][3] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].G12 ;
      DV[5][5] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].B31 ;
      DV[6][6] = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	E->tracer.coss[E->tracer.property_group[m]].B32 ;
    }
  }
  else {
    DV[1][1] = 2.0*dOmega ;
    DV[2][2] = 2.0*dOmega ;
    DV[3][3] = dOmega ;
    if(E->control.AXI)
      DV[4][4] = 2.0*dOmega ;
    if(3==dims) {
    DV[3][3] = 2.0*dOmega ;
    DV[4][4] = dOmega ;
    DV[5][5] = dOmega ;
    DV[6][6] = dOmega ;
    }
  }
  return;
}

void construct_D_Pres(
		      struct All_variables *E,
		      standard_precision DP[16][16],
		      standard_precision dOmegap,
		      standard_precision penalty
		      )
{
  int i,j ;
  const int dims=E->mesh.nsd ;
  standard_precision prod;

  for(i=1;i<16;i++)
    for(j=1;j<16;j++)
      DP[i][j] = 0.0 ;


  prod = dOmegap*penalty ;

  if(2==dims){
  
    DP[1][1] = prod ; DP[1][2] = prod ;
    DP[2][1] = prod ; DP[2][2] = prod ;
    
    if(E->control.AXI) {
      prod = dOmegap*penalty ;
      DP[1][1] = prod ; DP[1][2] = prod ;              DP[1][4] = prod ;
      DP[2][1] = prod ; DP[2][2] = prod ;	       DP[2][4] = prod ;
      DP[4][1] = prod ; DP[4][2] = prod ;              DP[4][4] = prod ;
    }
  }
  else {
    DP[1][1] = prod ; DP[1][2] = prod ; DP[1][3] = prod ;
    DP[2][1] = prod ; DP[2][2] = prod ; DP[2][3] = prod ;
    DP[3][1] = prod ; DP[3][2] = prod ; DP[3][3] = prod ;
  }
  return;
}





/* All this really has to do is copy the current
   stress calculated with the stress history, into the
   stress history. The routine calculating the stress
   takes care of large deformation measure or not */

void stress_update
(
 struct All_variables *E
 )
{
  standard_precision tau[7][7];
  standard_precision D[7][7];
  standard_precision pressure,Dkk;
  standard_precision smoothing_ratio;
  standard_precision Hvisc;

  int m,i;

  const int dofs = E->mesh.dof ;
  const int dims = E->mesh.nsd ;
  const int ends = enodes[dims] ;
  struct IEN *IEN = E->ien;

  void tracer_deviatoric_stress_and_pressure();

 


  /* The smoothing ratio is not defined for the first timestep.

     In the first pass, the timestep is either fixed, or defaults to
     1.0e-16 (supposedly a small number). The default value must also
     be used for advection.

     Beware of any elastic case where the default initial
     timestep is too large !

  */

  if(E->control.ELASTICITY) {
    smoothing_ratio = E->advection.timestep / E->advection.elastic_timestep ;     
    fprintf(E->fp,"Timestep = %g v %g, smoothing = %g\n",E->advection.timestep,E->advection.elastic_timestep,
	  smoothing_ratio);
  }

  if(2==dims && 2==dofs) {
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) { 
      
      tracer_deviatoric_stress_and_pressure(E,tau,D,&pressure,&Dkk,&Hvisc,m,E->mesh.levmax);
      E->tracer.Hvisc[m] = Hvisc;
      
      if(E->control.ELASTICITY && E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	E->tracer.S11[m] = smoothing_ratio * tau[1][1] + (1.0 -  smoothing_ratio) * E->tracer.S11[m];
	E->tracer.S22[m] = smoothing_ratio * tau[2][2] + (1.0 -  smoothing_ratio) * E->tracer.S22[m];
	E->tracer.S12[m] = smoothing_ratio * tau[1][2] + (1.0 -  smoothing_ratio) * E->tracer.S12[m];
	
	/*E->tracer.Pt[m] = pressure - 
	  E->tracer.visc[E->tracer.property_group[m]].Pen_bulk*E->tracer.Visc[m] * Dkk ;*/
	E->tracer.Pt[m] = smoothing_ratio * pressure + (1.0 -  smoothing_ratio) * E->tracer.Pt[m] ;
	/* C.O'N: FD's version has this line here instead of above: E->tracer.Pt[m] = smoothing_ratio * pressure + (1.0 -  smoothing_ratio) * E->tracer.Pt[m] ;
	 * 
 */
 	if (m == 436) {
	    fprintf(stderr,"HERE IS TRACER PRESSURE (3) %g tracer.Pt %g smoothing_ratio %g \n",pressure,E->tracer.Pt[m],smoothing_ratio);
    }

	
      }
    }
  }
  else if(2==dims && 3==dofs) {
    /* ???? should be similar, right ???? */
  }
  else if(3==dims && 3==dofs) {
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    

        /*fprintf(stderr,"Testing pressures going into this: %g \n",pressure);*/
	tracer_deviatoric_stress_and_pressure(E,tau,D,&pressure,&Dkk,&Hvisc,m,E->mesh.levmax);
/*	fprintf(stderr,"Testing pressures coming out of this: %g \n",pressure);*/
	
	/* C.O'N: this was an experiment - putting FD's version back in*/
	  if(E->control.ELASTICITY && E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0){
	  E->tracer.S11[m] = smoothing_ratio * tau[1][1] + (1.0 -  smoothing_ratio) * E->tracer.S11[m];
	  E->tracer.S22[m] = smoothing_ratio * tau[1][2] + (1.0 -  smoothing_ratio) * E->tracer.S22[m];
	  E->tracer.S33[m] = smoothing_ratio * tau[3][3] + (1.0 -  smoothing_ratio) * E->tracer.S33[m];
	  E->tracer.S12[m] = smoothing_ratio * tau[1][2] + (1.0 -  smoothing_ratio) * E->tracer.S12[m];
	  E->tracer.S13[m] = smoothing_ratio * tau[1][3] + (1.0 -  smoothing_ratio) * E->tracer.S13[m];
	  E->tracer.S23[m] = smoothing_ratio * tau[2][3] + (1.0 -  smoothing_ratio) * E->tracer.S23[m];

	  E->tracer.Pt[m] = smoothing_ratio * pressure + (1.0 -  smoothing_ratio) * E->tracer.Pt[m] ; 
	  /* C.O'N: eliminating these terms get rid of previous instabilities, and is consistent with FDs version (well, not know), 
	     but benchmark is not right. - 
	    E->tracer.visc[E->tracer.property_group[m]].Pen_bulk*E->tracer.Visc[m] * Dkk; */
	/*if(E->control.ELASTICITY && E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0){
		E->tracer.S11[m] = tau[1][1];
		E->tracer.S22[m] = tau[2][2];
		E->tracer.S33[m] = tau[3][3];
		E->tracer.S12[m] = tau[1][2];
		E->tracer.S13[m] = tau[1][3];
		E->tracer.S23[m] = tau[2][3];

		E->tracer.Pt[m] = pressure ;*/

	
	}
      }
    
  }
  else /*if(6==dofs)*/{
    /* ???? should be similar, right ???? */
  }
  fprintf(stderr,"HERE IS TRACER.PT (B)   %g \n",E->tracer.Pt[436]) ;	
  return ;
}

