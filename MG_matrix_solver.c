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


#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>



#define gauss_seidel6 gauss_seidel6_nbn


/* ==============================================================================
   Multigrid solver - takes RHS (R1...), assumes NO initial guess, and returns
   the solution in U1,U2,U3 after one MG cycle. Residual is returned in R
   ============================================================================== */


int multigrid_solve_del2_u(
			   struct All_variables *E,
			   standard_precision *U1,
			   standard_precision *U2,
			   standard_precision *U3,
			   standard_precision *U4,
			   standard_precision *U5,
			   standard_precision *U6,
			   standard_precision *RR1,
			   standard_precision *RR2,
			   standard_precision *RR3,
			   standard_precision *RR4,
			   standard_precision *RR5,
			   standard_precision *RR6,
			   standard_precision acc,
			   int LEVEL,
			   int relax
			   )
{
  void projection_integral();
  void interpolation_6();
  void interpolation_4pt_6();
  void n_assemble_del2_u_6();
  void direct_low_level();
  void gauss_seidel6();

  standard_precision residual_del2_u6();
  standard_precision return_bulk_value_l();
  standard_precision vsselfdot();

  standard_precision residual,r0,prior_residual,prior_residual2,delta_res,delta_res2;
  standard_precision AudotAu,AudotR,alpha;
  standard_precision c1,c2,aveVx;

  int counter[MAX_LEVELS],completed_cycle;
  int count,cycles,convergent,loop;
  int i,j,k,l,m;
  int level,current_level;
  
  standard_precision CPU_time(),initial_time,time,omeg;

  static standard_precision *V1[MAX_LEVELS],*V2[MAX_LEVELS],*V3[MAX_LEVELS];
  static standard_precision *V4[MAX_LEVELS],*V5[MAX_LEVELS],*V6[MAX_LEVELS];
  static standard_precision *F1[MAX_LEVELS],*F2[MAX_LEVELS],*F3[MAX_LEVELS];
  static standard_precision *F4[MAX_LEVELS],*F5[MAX_LEVELS],*F6[MAX_LEVELS];
  static standard_precision *R1[MAX_LEVELS],*R2[MAX_LEVELS],*R3[MAX_LEVELS];
  static standard_precision *R4[MAX_LEVELS],*R5[MAX_LEVELS],*R6[MAX_LEVELS];
  static standard_precision *v1,*v2,*v3,*v4,*v5,*v6;
  static standard_precision *Au1,*Au2,*Au3,*Au4,*Au5,*Au6;
  static int been_here = 0;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int max_level=LEVEL;

  int count1,count2;

  standard_precision Vmean1,Vmean2,ave ;

  Vmean1 = Vmean2 = 0.0 ;
  
  if(been_here++==0) {
    E->control.total_iteration_cycles = 0;
    E->control.total_v_solver_calls = 0;

    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
      switch(dofs) {
      case 2:
	V1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	break ;
      case 3:
	V1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	break ;
      case 6:
	V1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V4[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V5[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	V6[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R4[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R5[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	R6[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F1[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F2[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F3[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F4[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F5[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	F6[i] = (standard_precision *)Malloc0((E->mesh.NNO[i] + 20)*sizeof(standard_precision));
	break ;
      }
    }
    Au1 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    Au2 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    Au3 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    Au4 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    Au5 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    Au6 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v1 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v2 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v3 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v4 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v5 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
    v6 = (standard_precision *)Malloc0((E->mesh.nno+2) * sizeof(standard_precision));
  }

  /* 1. convert format */
  r0=0.0;
  switch(dofs) {
  case 2:
    for(i=1;i<=E->mesh.NNO[max_level];i++) {
      U1[i] = U2[i] = 
	V1[max_level][i] = V2[max_level][i] = 0.0;
      F1[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR1[i];
      F2[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR2[i];
      
      r0 += RR1[i] *  RR1[i] +  RR2[i] * RR2[i];
    }
    break ;
  case 3:
    for(i=1;i<=E->mesh.NNO[max_level];i++) {
      U1[i] = U2[i] = U3[i] = 
	V1[max_level][i] = V2[max_level][i] = V3[max_level][i] = 0.0;
      
      F1[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR1[i];
      F2[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR2[i];
      F3[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR3[i];
      r0 += RR1[i] *  RR1[i] +  RR2[i] * RR2[i] +  RR3[i] * RR3[i];
    }
    break ;
  case 6:
    for(i=1;i<=E->mesh.NNO[max_level];i++) {
      U1[i] = U2[i] = U3[i] =  U4[i] = U5[i] = U6[i] = 
	V1[max_level][i] = V2[max_level][i] = V3[max_level][i] =
	V4[max_level][i] = V5[max_level][i] = V6[max_level][i] = 0.0;
      
      F1[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR1[i];
      F2[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR2[i];
      F3[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR3[i];
      F4[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR4[i];
      F5[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR5[i];
      F6[max_level][i] = ((E->NODE[max_level][i] & ( OFFSIDE)) == 0) * RR6[i];

      r0 += RR1[i] *  RR1[i] +  RR2[i] * RR2[i] +  RR3[i] * RR3[i] 
	+ RR4[i] *  RR4[i] +  RR5[i] * RR5[i] +  RR6[i] * RR6[i];
    }
    break ;
  }

  residual=r0=sqrt(r0/E->mesh.NNO[max_level]);

  if(0.0 == r0) {
    fprintf(stderr,"No significant residual, abandoning MG solve !\n");
    return(0);
  }

  prior_residual2 = prior_residual = 2 * residual; 
  count = 0;  
  initial_time=CPU_time();

  /* If the accuracy is larger than the initial residual, then
     the iteration may make the residual diverge ! */

  acc = min(acc,0.9999*r0);
  
  /* Execute the V, W, V3 cycle until residual is small enough */

  if(E->mesh.levmax==E->mesh.levmin) {
    direct_low_level(E,V1[E->mesh.levmin],V2[E->mesh.levmin],V3[E->mesh.levmin],
		     V4[E->mesh.levmin],V5[E->mesh.levmin],V6[E->mesh.levmin],
		     F1[E->mesh.levmin],F2[E->mesh.levmin],F3[E->mesh.levmin],
		     F4[E->mesh.levmin],F5[E->mesh.levmin],F6[E->mesh.levmin]); 
  }
  else
  do {
    /* Set up cycle counters (for V, W etc) */
    for(i=E->mesh.levmin;i<=max_level;i++)
      counter[i]=0;
    
    level=max_level;
    completed_cycle=0;
 
    /*  presmoothing */
    gauss_seidel6(E,V1[max_level],V2[max_level],V3[max_level],V4[max_level],V5[max_level],V6[max_level],
		  F1[max_level],F2[max_level],F3[max_level],F4[max_level],F5[max_level],F6[max_level],
		  Au1,Au2,Au3,Au4,Au5,Au6,acc,1.0,
		  E->control.vsolver_relaxations,max_level,relax);
    residual = residual_del2_u6(E,
				F1[max_level],F2[max_level],F3[max_level],F4[max_level],F5[max_level],F6[max_level],
				V1[max_level],V2[max_level],V3[max_level],V4[max_level],V5[max_level],V6[max_level],
				R1[max_level],R2[max_level],R3[max_level],R4[max_level],R5[max_level],R6[max_level],
				Au1,Au2,Au3,Au4,Au5,Au6,
				max_level,1);
    /* Main multilevel loop */
  
    while(! completed_cycle) {
      if((counter[level] < E->control.mg_cycle && level > E->mesh.levmin)) {	
	projection_integral(E,R1[level],R2[level],R3[level],R4[level],R5[level],R6[level],
			    F1[level-1],F2[level-1],F3[level-1],F4[level-1],F5[level-1],F6[level-1],
			    level); 

	/* We have projected only the residual force term so that the velocity at this level 
	   is entirely unknown. Our best guess is therefore Zero. */
	 
	for(i=1;i<=E->mesh.NNO[level-1];i++) {
	  V1[level-1][i] = V2[level-1][i] = 0.0;
	  if(3==dofs) 
	    V3[level-1][i] = 0.0;
	  else if(6==dofs) 
	    V3[level-1][i] = V4[level-1][i] = V5[level-1][i] = V6[level-1][i] = 0.0;
	}

	if(level==E->mesh.levmin+1  /* && !(E->mesh.periodic_x || E->mesh.periodic_y) */ ) {
	  /* Problems here with periodic bcs */
	  direct_low_level(E,V1[E->mesh.levmin],V2[E->mesh.levmin],V3[E->mesh.levmin],
			   V4[E->mesh.levmin],V5[E->mesh.levmin],V6[E->mesh.levmin],
			   F1[E->mesh.levmin],F2[E->mesh.levmin],F3[E->mesh.levmin],
			   F4[E->mesh.levmin],F5[E->mesh.levmin],F6[E->mesh.levmin]); 
	  if(E->control.vector_optimization)
	    e_assemble_del2_u_6(E,
				V1[E->mesh.levmin],V2[E->mesh.levmin],V3[E->mesh.levmin],
				V4[E->mesh.levmin],V5[E->mesh.levmin],V6[E->mesh.levmin],
				Au1,Au2,Au3,Au4,Au5,Au6,
				E->mesh.levmin,1);
	  else
	    n_assemble_del2_u_6(E,
				V1[E->mesh.levmin],V2[E->mesh.levmin],V3[E->mesh.levmin],
				V4[E->mesh.levmin],V5[E->mesh.levmin],V6[E->mesh.levmin],
				Au1,Au2,Au3,Au4,Au5,Au6,
				E->mesh.levmin,1);
	}
	else if (level==E->mesh.levmin+1) {/* Periodic "exact solution" */ 
	/*RAA: 30/4/01, no difference between this 'if' and the one above it when periodic part is commented out*/
	  gauss_seidel6(E,V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
			F1[level-1],F2[level-1],F3[level-1],F4[level-1],F5[level-1],F6[level-1],
			Au1,Au2,Au3,Au4,Au5,Au6,
			acc*0.1,0.1,
			E->control.vsolver_relaxations,level-1,relax);
	}
	else {
	  gauss_seidel6(E,V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
			F1[level-1],F2[level-1],F3[level-1],F4[level-1],F5[level-1],F6[level-1],
			Au1,Au2,Au3,Au4,Au5,Au6,
			acc,1.0,
			E->control.vsolver_relaxations,level-1,relax);
	}

	residual= residual_del2_u6(E,
				   F1[level-1],F2[level-1],F3[level-1],F4[level-1],F5[level-1],F6[level-1],
				   V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
				   R1[level-1],R2[level-1],R3[level-1],R4[level-1],R5[level-1],R6[level-1],
				   Au1,Au2,Au3,Au4,Au5,Au6,
				   level-1,0);
	counter[level]++; 
	level--;
      }
      else if (level != max_level) {     /* go up one level and smooth */  
	counter[level]=0;
	level++;
	if(level==max_level) 
	  completed_cycle=1;
	if((E->mesh.NOX[level-1] <=4) || (E->mesh.NOZ[level-1] <=4) || (3==dims && (E->mesh.NOY[level-1] <=4))) {

	  interpolation_6(E,v1,v2,v3,v4,v5,v6,
			  V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
			  level);
	}
	else {
#if 0
	  interpolation_6(E,v1,v2,v3,v4,v5,v6,
			  V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
			  level);
#else
	  interpolation_4pt_6(E,v1,v2,v3,v4,v5,v6,
			      V1[level-1],V2[level-1],V3[level-1],V4[level-1],V5[level-1],V6[level-1],
			      level);
#endif
	}
	/* Smooth the resulting interpolation and then calculate if it will be of
	   any benefit to use it (the alpha factor should be close to unity if
	   the solution is good). If it is rubbish, and we are at a lower level,
	   then iterate some more to save time at higher levels later. We are most
	   anxious to avoid cases where small alpha gets all the MG work thrown away.
	   This explains the choice of parameters */
	
	alpha=-1.0;
	loop=0;
	
	while(loop++ < 1 * (E->mesh.levmax-level+1) &&  (alpha < 0.0 || alpha > 5.0))  {
	  gauss_seidel6(E,v1,v2,v3,v4,v5,v6,
			R1[level],R2[level],R3[level],R4[level],R5[level],R6[level],
			Au1,Au2,Au3,Au4,Au5,Au6,
			acc,0.95,
			E->control.vsolver_relaxations,level,relax); 	   
	  AudotAu= AudotR=0.0;
	  if(2==dofs)
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      AudotR += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * (R1[level][i] * Au1[i] + R2[level][i] * Au2[i]);
	      AudotAu += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * (Au1[i] * Au1[i] + Au2[i] * Au2[i]) ;
	    }
	  else if(dofs==3)
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      AudotR += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * 
		(R1[level][i] * Au1[i] + R2[level][i] * Au2[i] + R3[level][i] * Au3[i]) ;
	      AudotAu += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * 
		(Au1[i] * Au1[i] + Au2[i] * Au2[i] + Au3[i] * Au3[i]) ;
	    }
	  else /* if(dofs==6) */
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      AudotR += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * 
		(R1[level][i] * Au1[i] + R2[level][i] * Au2[i] + R3[level][i] * Au3[i] 
		 + R4[level][i] * Au4[i] + R5[level][i] * Au5[i] + R6[level][i] * Au6[i]) ;
	      AudotAu += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * 
		(Au1[i] * Au1[i] + Au2[i] * Au2[i] + Au3[i] * Au3[i] 
		 + Au4[i] * Au4[i] + Au5[i] * Au5[i] + Au6[i] * Au6[i]) ;
	    }
	  if(fabs(AudotAu) > 1.0e-32 && fabs(AudotR / AudotAu) > 1.0e-32)
	    alpha = AudotR / AudotAu;
	  else {
	    alpha = 0.0;
	    break;
	  }
	}
	if(alpha > 0.0 && alpha < 100.0) {
	  if(2==dofs)
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      V1[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v1[i];
	      V2[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v2[i];
	    }
	  else if(dofs==3)
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      V1[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v1[i];
	      V2[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v2[i];
	      V3[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v3[i];
	    }
	  else /*if(dofs==6)*/
	    for(i=1;i<=E->mesh.NNO[level];i++) {
	      V1[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v1[i];
	      V2[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v2[i];
	      V3[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v3[i];
	      V4[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v4[i];
	      V5[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v5[i];
	      V6[level][i] += ((E->NODE[level][i] & ( OFFSIDE)) == 0) * alpha * v6[i];
	    }
	}
      }
      else {
	completed_cycle=1;
      }
    }

    if(E->control.vector_optimization) 
      e_assemble_del2_u_6(E,V1[max_level],V2[max_level],V3[max_level],V4[max_level],V5[max_level],V6[max_level],
			  Au1,Au2,Au3,Au4,Au5,Au6,max_level,1);
    else
      n_assemble_del2_u_6(E,V1[max_level],V2[max_level],V3[max_level],V4[max_level],V5[max_level],V6[max_level],
			  Au1,Au2,Au3,Au4,Au5,Au6,max_level,1);

    residual = residual_del2_u6(E,
				F1[max_level],F2[max_level],F3[max_level],F4[max_level],F5[max_level],F6[max_level],
				V1[max_level],V2[max_level],V3[max_level],V4[max_level],V5[max_level],V6[max_level],
				R1[max_level],R2[max_level],R3[max_level],R4[max_level],R5[max_level],R6[max_level],
				Au1,Au2,Au3,Au4,Au5,Au6,
				max_level,0); 

    delta_res = fabs(prior_residual - residual);
    delta_res2 = fabs(prior_residual2 -residual);

    /* Residual not declining very fast */
    if((count > 10) && ((residual > r0*10.0) || (delta_res < acc*0.01 && delta_res2 < acc*0.01))) {
      convergent=0;      
    }
    else {
      convergent=1;
      prior_residual2 = prior_residual;
      prior_residual=residual;
    }
    
    if(count > 25 || E->control.print_convergence)
      fprintf(stderr,"%s %02d residual (%03d) = %.3e from %.3e to %.3e (%.3e) in %5.2f secs - l=%d \n",
	      (convergent ? " * ":"!!!"),max_level,count,residual,r0,acc,delta_res,
	      CPU_time()-initial_time,E->control.vsolver_relaxations);

    count++;

    if ((count%100)==0)
      fprintf(E->fp,"Warning, V solver's done %d V cycles at level %d, residual reduced from %g to %g\n",
	      count,max_level,r0,residual);




  } while ((residual > acc) && (count < E->control.max_vel_iterations) && convergent );

  if(count >= E->control.max_vel_iterations) {
    fprintf(stderr,"Exceeded %d V cycles @ high level %d (%d)\n",E->control.max_vel_iterations,max_level,count);
    convergent=0;
  }
  
  /* 3. convert format */
 
  if(E->mesh.periodic_x && E->mesh.periodic_rm_vx) {
    aveVx = return_bulk_value_l(E,V1[max_level],1,max_level);
    for(i=1;i<=E->mesh.NNO[max_level];i++)
      V1[max_level][i] -= aveVx;

    if(E->control.verbose)
      fprintf(stderr,"Removing %g from average x velocity \n",aveVx);

    fprintf(E->fp,"[%d] Removing %g from average x velocity \n",max_level,aveVx);
    /*RAA: 5/4/02, why are the lines below repeated? - seems incorrect, so
      I will comment them out*/ 
/*
    aveVx = return_bulk_value_l(E,V1[max_level],1,max_level);
    for(i=1;i<=E->mesh.NNO[max_level];i++)
      V1[max_level][i] -= aveVx;

    if(E->control.verbose)
      fprintf(stderr,"Removing %g from average x velocity \n",aveVx);

    fprintf(E->fp,"[%d] Removing %g from average x velocity \n",max_level,aveVx);
 */
  }

  /*RAA: 02/04/02, here is the periodic_y addition, but use aveVx to save space*/
  if(E->mesh.periodic_y && E->mesh.periodic_rm_vy) {
    aveVx = return_bulk_value_l(E,V3[max_level],1,max_level);
    for(i=1;i<=E->mesh.NNO[max_level];i++)
      V3[max_level][i] -= aveVx;

    if(E->control.verbose)
     fprintf(stderr,"Removing %g from average y velocity \n",aveVx);

    fprintf(E->fp,"[%d] Removing %g from average y velocity \n",max_level,aveVx);
  }


  
  /* 0 if couple-stress applied */
  /*  if(E->control.outputcase==4) {
      omeg = V3[max_level][1] ;
      }
      else*/
  omeg = 0.0 ;

  for(i=1;i<=E->mesh.NNO[max_level];i++) {
    U1[i] = V1[max_level][i];
    U2[i] = V2[max_level][i];
    if(dofs==3 && dims==2) {
      U3[i] = V3[max_level][i] - omeg;
    }
    else if(dofs==3 && dims==3)
      U3[i] = V3[max_level][i];  /*RAA: 9/5/01, this part missing for 3D */
    else if(dofs==6) {
      U3[i] = V3[max_level][i];
      U4[i] = V4[max_level][i];
      U5[i] = V5[max_level][i];
      U6[i] = V6[max_level][i];
    }
  }
  E->control.total_iteration_cycles += count;
  E->control.total_v_solver_calls +=1;

#if 0
  /* Double check ... does claimed residual value add up ? */

  n_assemble_del2_u_6(E,U1,U2,U3,U4,U5,U6,Au1,Au2,Au3,Au4,Au5,Au6,max_level,1);
  c1=c2=0.0;

  for(i=1;i<=E->mesh.NNO[max_level];i++) {
    if(E->NODE[max_level][i] & ( OFFSIDE  ))
      continue;
    c1 += (F1[max_level][i] - Au1[i]) * (F1[max_level][i] - Au1[i]) + 
      (F2[max_level][i] - Au2[i]) * (F2[max_level][i] - Au2[i]);
    c2 += F1[max_level][i] * F1[max_level][i] + F2[max_level][i] * F2[max_level][i];
  }
  
  fprintf(stderr,"MG test 1: residual = %g (relative %g)\n",
	  sqrt(c1/E->mesh.NNO[max_level]),
	  sqrt(c1/c2));
#endif
  return(convergent);
}

/* Direct solver for low level MG iteration */

void direct_low_level(
  struct All_variables *E,
  standard_precision *V1,
  standard_precision *V2,
  standard_precision *V3,
  standard_precision *V4,
  standard_precision *V5,
  standard_precision *V6,
  standard_precision *R1,
  standard_precision *R2,
  standard_precision *R3,
  standard_precision *R4,
  standard_precision *R5,
  standard_precision *R6
)
{
  static standard_precision *A;
  static standard_precision *F;
  static standard_precision *U;
  
  static int been_here=0;

  int i,j,k;
  int eq1,eq2,eq3,eq4,eq5,eq6;
  int Eq1,Eq2,Eq3,Eq4,Eq5,Eq6;
  int node1;
  
  double t1,t2;

  const int neq = E->mesh.neqd;
  const int nno = E->mesh.NNO[E->mesh.levmin];
  const int dofs = E->mesh.dof;
  const int dims = E->mesh.nsd;
  const int level = E->mesh.levmin;
  const int max_node = max_node_interaction[dims];


  if(!been_here++) {
    A = (standard_precision *) Malloc0((neq*neq+1) * sizeof (standard_precision));
    U = (standard_precision *) Malloc0((neq+1) * sizeof (standard_precision));
    F = (standard_precision *) Malloc0((neq+1) * sizeof (standard_precision));
  }

  /* need to assemble the matrix A, and the vector F */

  for(i=0;i<neq;i++) {
    F[i] = 0.0;
    for(j=0;j<neq;j++)
      A[i+j*neq] = 0.0;
  }

  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  ))
      continue;
 
    eq1=E->idd[i].doff[1];
    eq2=E->idd[i].doff[2];
    eq3=E->idd[i].doff[3];
    eq4=E->idd[i].doff[4];
    eq5=E->idd[i].doff[5];
    eq6=E->idd[i].doff[6];
    for(j=0;j<max_node;j++) {
      node1=E->Node_map_3[level][i*max_node+j];
      if(node1==0)
	continue;

      if(E->NODE[level][node1] & ( OFFSIDE  ))  /* shouldn't happen though */
	continue;

      switch(dofs) {
      case 2:
	Eq1=E->idd[node1].doff[1];
	Eq2=E->idd[node1].doff[2];
	if(eq1 != -1 && Eq1 != -1)
	  A[eq1+neq*Eq1] += E->Node_k11[level][i*max_node+j];
	if(eq2 != -1 && Eq1 != -1)  
	  A[eq2+neq*Eq1] += E->Node_k21[level][i*max_node+j];
	if(eq1 != -1 && Eq2 != -1)   
	  A[eq1+neq*Eq2] += E->Node_k12[level][i*max_node+j];
	if(eq2 != -1 && Eq2 != -1)  
	  A[eq2+neq*Eq2] += E->Node_k22[level][i*max_node+j];
	break ;
      case 3:
	Eq1=E->idd[node1].doff[1];
	Eq2=E->idd[node1].doff[2];
	Eq3=E->idd[node1].doff[3];
	if(eq1 != -1 && Eq1 != -1)
	  A[eq1+neq*Eq1] += E->Node_k11[level][i*max_node+j];
	if(eq2 != -1 && Eq1 != -1)  
	  A[eq2+neq*Eq1] += E->Node_k21[level][i*max_node+j];
	if(eq1 != -1 && Eq2 != -1)   
	  A[eq1+neq*Eq2] += E->Node_k12[level][i*max_node+j];
	if(eq2 != -1 && Eq2 != -1)  
	  A[eq2+neq*Eq2] += E->Node_k22[level][i*max_node+j];
	eq3=E->idd[i].doff[3];
	
	if(eq1 != -1 && Eq3 != -1)
	  A[eq1+neq*Eq3] += E->Node_k13[level][i*max_node+j];
	if(eq2 != -1 && Eq3 != -1)
	  A[eq2+neq*Eq3] += E->Node_k23[level][i*max_node+j];
	if(eq3 != -1 && Eq1 != -1)
	  A[eq3+neq*Eq1] += E->Node_k31[level][i*max_node+j];
	if(eq3 != -1 && Eq2 != -1)
	  A[eq3+neq*Eq2] += E->Node_k32[level][i*max_node+j];
	if(eq3 != -1 && Eq3 != -1)
	  A[eq3+neq*Eq3] += E->Node_k33[level][i*max_node+j];
	break ;
      case 6:
	Eq1=E->idd[node1].doff[1];
	Eq2=E->idd[node1].doff[2];
	Eq3=E->idd[node1].doff[3];
      	Eq4=E->idd[node1].doff[4];
	Eq5=E->idd[node1].doff[5];
	Eq6=E->idd[node1].doff[6];
      
	if(eq1 != -1 && Eq1 != -1)
	  A[eq1+neq*Eq1] += E->Node_k11[level][i*max_node+j];
	if(eq2 != -1 && Eq1 != -1)  
	  A[eq2+neq*Eq1] += E->Node_k21[level][i*max_node+j];
	if(eq3 != -1 && Eq1 != -1)
	  A[eq3+neq*Eq1] += E->Node_k31[level][i*max_node+j];
	if(eq4 != -1 && Eq1 != -1)
	  A[eq4+neq*Eq1] += E->Node_k41[level][i*max_node+j];
	if(eq5 != -1 && Eq1 != -1)  
	  A[eq5+neq*Eq1] += E->Node_k51[level][i*max_node+j];
	if(eq6 != -1 && Eq1 != -1)
	  A[eq6+neq*Eq1] += E->Node_k61[level][i*max_node+j];

	if(eq1 != -1 && Eq2 != -1)
	  A[eq1+neq*Eq2] += E->Node_k12[level][i*max_node+j];
	if(eq2 != -1 && Eq2 != -1)  
	  A[eq2+neq*Eq2] += E->Node_k22[level][i*max_node+j];
	if(eq3 != -1 && Eq2 != -1)
	  A[eq3+neq*Eq2] += E->Node_k32[level][i*max_node+j];
	if(eq4 != -1 && Eq2 != -1)
	  A[eq4+neq*Eq2] += E->Node_k42[level][i*max_node+j];
	if(eq5 != -1 && Eq2 != -1)  
	  A[eq5+neq*Eq2] += E->Node_k52[level][i*max_node+j];
	if(eq6 != -1 && Eq2 != -1)
	  A[eq6+neq*Eq2] += E->Node_k62[level][i*max_node+j];
	
      if(eq1 != -1 && Eq3 != -1)
	A[eq1+neq*Eq3] += E->Node_k13[level][i*max_node+j];
      if(eq2 != -1 && Eq3 != -1)  
	A[eq2+neq*Eq3] += E->Node_k23[level][i*max_node+j];
      if(eq3 != -1 && Eq3 != -1)
        A[eq3+neq*Eq3] += E->Node_k33[level][i*max_node+j];
      if(eq4 != -1 && Eq3 != -1)
	A[eq4+neq*Eq3] += E->Node_k43[level][i*max_node+j];
      if(eq5 != -1 && Eq3 != -1)  
	A[eq5+neq*Eq3] += E->Node_k53[level][i*max_node+j];
      if(eq6 != -1 && Eq3 != -1)
	A[eq6+neq*Eq3] += E->Node_k63[level][i*max_node+j];

      if(eq1 != -1 && Eq4 != -1)
	A[eq1+neq*Eq4] += E->Node_k14[level][i*max_node+j];
      if(eq2 != -1 && Eq4 != -1)  
	A[eq2+neq*Eq4] += E->Node_k24[level][i*max_node+j];
      if(eq3 != -1 && Eq4 != -1)
        A[eq3+neq*Eq4] += E->Node_k34[level][i*max_node+j];
      if(eq4 != -1 && Eq4 != -1)
	A[eq4+neq*Eq4] += E->Node_k44[level][i*max_node+j];
      if(eq5 != -1 && Eq4 != -1)  
	A[eq5+neq*Eq4] += E->Node_k54[level][i*max_node+j];
      if(eq6 != -1 && Eq4 != -1)
	A[eq6+neq*Eq4] += E->Node_k64[level][i*max_node+j];

      if(eq1 != -1 && Eq5 != -1)
	A[eq1+neq*Eq5] += E->Node_k15[level][i*max_node+j];
      if(eq2 != -1 && Eq5 != -1)  
	A[eq2+neq*Eq5] += E->Node_k25[level][i*max_node+j];
      if(eq3 != -1 && Eq5 != -1)
        A[eq3+neq*Eq5] += E->Node_k35[level][i*max_node+j];
      if(eq4 != -1 && Eq5 != -1)
	A[eq4+neq*Eq5] += E->Node_k45[level][i*max_node+j];
      if(eq5 != -1 && Eq5 != -1)  
	A[eq5+neq*Eq5] += E->Node_k55[level][i*max_node+j];
      if(eq6 != -1 && Eq5 != -1)
	A[eq6+neq*Eq5] += E->Node_k65[level][i*max_node+j];

      if(eq1 != -1 && Eq6 != -1)
	A[eq1+neq*Eq6] += E->Node_k16[level][i*max_node+j];
      if(eq2 != -1 && Eq6 != -1)  
	A[eq2+neq*Eq6] += E->Node_k26[level][i*max_node+j];
      if(eq3 != -1 && Eq6 != -1)
        A[eq3+neq*Eq6] += E->Node_k36[level][i*max_node+j];
      if(eq4 != -1 && Eq6 != -1)
	A[eq4+neq*Eq6] += E->Node_k46[level][i*max_node+j];
      if(eq5 != -1 && Eq6 != -1)  
	A[eq5+neq*Eq6] += E->Node_k56[level][i*max_node+j];
      if(eq6 != -1 && Eq6 != -1)
	A[eq6+neq*Eq6] += E->Node_k66[level][i*max_node+j];

	break ;
      }
    }
  } /* next node */

 
  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  ))
      continue;

    switch(dofs) {
    case 2:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      if(eq1 != -1)
	F[eq1] = R1[i];
      if(eq2 != -1)
	F[eq2] = R2[i];
      break ;
    case 3:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      eq3=E->idd[i].doff[3];
      if(eq1 != -1)
	F[eq1] = R1[i];
      if(eq2 != -1)
	F[eq2] = R2[i];
      if(eq3 != -1)
	F[eq3] = R3[i];
      break ;
    case 6:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      eq3=E->idd[i].doff[3];
      eq4=E->idd[i].doff[4];
      eq5=E->idd[i].doff[5];
      eq6=E->idd[i].doff[6];
      if(eq1 != -1)
	F[eq1] = R1[i];
      if(eq2 != -1)
	F[eq2] = R2[i];
      if(eq3 != -1)
	F[eq3] = R3[i];
      if(eq4 != -1)
	F[eq4] = R4[i];
      if(eq5 != -1)
	F[eq5] = R5[i];
      if(eq6 != -1)
	F[eq6] = R6[i];
      break ;
    }
  }

  /*RAA: check*/
  /*   for(i=1;i<=neq*neq;i++) 
      if(3==dims)
      fprintf(E->fp1,"w/i direct_low': Ai  %d  %g\n",i,A[i]); 
      */

  /* have built the matrix problem, now it can be solved (by elimination, for example) */

  /* Factorization loops */
  
  for(j=0;j<neq;j++)  {
    for(i=1;i<j;i++) {
      t1=0.0;
      for(k=0;k<i;k++)
	t1 += A[k+i*neq] * A[k+j*neq];
      A[i+j*neq] -= t1;
    }
 
    for(i=0;i<j;i++) {
      if(A[i+i*neq] == 0.0)
	continue;
      t1=A[i+j*neq];
      A[i+j*neq] /= A[i+i*neq];
      A[j+j*neq] -= t1 * A[i+j*neq];
    }
  }

  /* Forward reduction */
  
  for(j=1;j<neq;j++) {
    t1 = 0.0;
    for(i=0;i<j;i++)
      t1 += A[i+j*neq] * F[i];
    F[j] -= t1;
  }
  
  /* Diagonal Scaling */

  for(j=0;j<neq;j++) {
    if((t1=A[j+j*neq]) != 0.0)
      F[j] /= t1;
  }

  /* Back substitution */
  for(j=neq-1;j>0;j--) {
    t1=F[j];
    for(i=0;i<j;i++)
      F[i] -= A[i+j*neq] * t1;
  }
  
  /* This is the solution, now extract into the return format */


  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  )) 
      continue;
      
    switch(dofs) {
    case 2:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      if(eq1 != -1)
	V1[i] = F[eq1];
      else
	V1[i] = 0.0;
      if(eq2 != -1)
	V2[i] = F[eq2];
      else
	V2[i] = 0.0;
      break ;
    case 3:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      eq3=E->idd[i].doff[3];
      if(eq1 != -1)
	V1[i] = F[eq1];
      else
	V1[i] = 0.0;
      if(eq2 != -1)
	V2[i] = F[eq2];
      else
	V2[i] = 0.0;
      if(eq3 != -1)
	V3[i] = F[eq3];
      else 
	V3[i] = 0.0;
      break ;
    case 6:
      eq1=E->idd[i].doff[1];
      eq2=E->idd[i].doff[2];
      eq3=E->idd[i].doff[3];
      eq4=E->idd[i].doff[4];
      eq5=E->idd[i].doff[5];
      eq6=E->idd[i].doff[6];
      if(eq1 != -1)
	V1[i] = F[eq1];
      else
	V1[i] = 0.0;
      if(eq2 != -1)
	V2[i] = F[eq2];
      else
	V2[i] = 0.0;
      if(eq3 != -1)
	V3[i] = F[eq3];
      else 
	V3[i] = 0.0;
      if(eq4 != -1)
	V4[i] = F[eq4];
      else
	V4[i] = 0.0;
      if(eq5 != -1)
	V5[i] = F[eq5];
      else
	V5[i] = 0.0;
      if(eq6 != -1)
	V6[i] = F[eq6];
      else 
	V6[i] = 0.0;
      break ;
    }
    }

return;
}
