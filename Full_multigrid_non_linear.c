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

void solve_constrained_flow_iterative(
				      struct All_variables *E
				      )
{
  standard_precision *U;
  higher_precision residual_ddash;
  standard_precision vsselfdot(),lpselfdot();

  static int been_here = 0;
   
  standard_precision solve_Ahat_q_fhat();

  void s_strip_bcs_from_residual();
  void sp_to_nodes();
  void remove_surf_horiz_press_ave();

  int cycles;
  int i,j,m,vel_cycles_previous,vel_calls_previous;
  
  standard_precision time,CPU_time();
  standard_precision work1,work2,work3;
  higher_precision *node_R;
  
  const int npno = E->mesh.npno;
  const int nno = E->mesh.nno;
  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int neq = E->mesh.neq;
 

  if(E->control.verbose)
    fprintf(stderr,"Incompressibility via CG Uzawa iteration\n");
  

  E->monitor.elapsed_time_vsoln1 =  E->monitor.elapsed_time_vsoln;
  E->monitor.elapsed_time_vsoln = E->monitor.elapsed_time;

  /* First copy existing fine-grid velocity field for later use */

  if(E->advection.timestep > 0)
      for(i=1;i<=nno;i++)
	for(j=1;j<=dofs;j++)
	  E->V1[j][i] =  E->V[j][i] ;

/* For changing BC */

#if 0
  if(E->advection.timesteps==30) {
    for(i=1;i<=nno;i++) {
      if(!(E->NODE[E->mesh.levmax][i] & OFFSIDE)) {
	if (E->NODE[E->mesh.levmax][i] & BC1)
	  E->Vb[1][E->mesh.levmax][i] = 0.0; 
	if (E->NODE[E->mesh.levmax][i] & BC2)
	  E->Vb[2][E->mesh.levmax][i] = 0.0; 
      }
    }
    for(i=E->mesh.levmax;i>E->mesh.levmin;i--) {
      inject_node_values(E,i,E->Vb[1][i],E->Vb[1][i-1]);
      inject_node_values(E,i,E->Vb[2][i],E->Vb[2][i-1]);
    }
  }
#endif

   /* Need to convert the standard cartesian velocity field into
     the mixed cartesian/curvilinear in which boundary conditions 
     are applied */

  if(3==dims && 3==dofs)  { /* Classical 3D */
   /* fprintf(stderr,"IN IF DIMS=3 STATEMENT  \n") ;	*/
    for(i=1;i<=nno;i++) {
      if(E->node[i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[E->mesh.levmax][i];
	work1 = node_R[0*dims+0] * E->V[1][i] + node_R[1*dims+0] * E->V[2][i] + node_R[2*dims+0] * E->V[3][i];
	work2 = node_R[0*dims+1] * E->V[1][i] + node_R[1*dims+1] * E->V[2][i] + node_R[2*dims+1] * E->V[3][i];
	work3 = node_R[0*dims+2] * E->V[1][i] + node_R[1*dims+2] * E->V[2][i] + node_R[2*dims+2] * E->V[3][i];
	E->V[1][i] = work1 ;
	E->V[2][i] = work2 ;
	E->V[3][i] = work3 ;
      }
    }
  }
  else if(3==dims && 6==dofs) /* Cosserat 3D */
    for(i=1;i<=nno;i++) {
      if(E->node[i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[E->mesh.levmax][i];
	work1 = node_R[0*dims+0] * E->V[1][i] + node_R[1*dims+0] * E->V[2][i] + node_R[2*dims+0] * E->V[3][i];
	work2 = node_R[0*dims+1] * E->V[1][i] + node_R[1*dims+1] * E->V[2][i] + node_R[2*dims+1] * E->V[3][i];
	work3 = node_R[0*dims+2] * E->V[1][i] + node_R[1*dims+2] * E->V[2][i] + node_R[2*dims+2] * E->V[3][i];
	E->V[1][i] = work1 ;
	E->V[2][i] = work2 ;
	E->V[3][i] = work3 ;
	work1 = node_R[0*dims+0] * E->V[4][i] + node_R[1*dims+0] * E->V[5][i] + node_R[2*dims+0] * E->V[6][i];
	work2 = node_R[0*dims+1] * E->V[4][i] + node_R[1*dims+1] * E->V[5][i] + node_R[2*dims+1] * E->V[6][i];
	work3 = node_R[0*dims+2] * E->V[4][i] + node_R[1*dims+2] * E->V[5][i] + node_R[2*dims+2] * E->V[6][i];
	E->V[4][i] = work1 ;
	E->V[5][i] = work2 ;
	E->V[6][i] = work3 ;
      }
    }
  else if(2==dims)   /* Classical 2D (and Cosserat: Omega3 unchanged by rotation about y axis) */
    for(i=1;i<=E->mesh.nno;i++) {
      if(E->node[i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[E->mesh.levmax][i];
	work1 = node_R[0*dims+0] * E->V[1][i] + node_R[1*dims+0] * E->V[2][i] ;
	work2 = node_R[0*dims+1] * E->V[1][i] + node_R[1*dims+1] * E->V[2][i] ;
	E->V[1][i] = work1 ;
	E->V[2][i] = work2 ;
      }
    }
  
   /* and signal this to the rest of the program */
   E->control.CURRENT_SKEWED_Vs = 1;
   /*fprintf(stderr,"GOING INTO BCS HERE  \n") ;	 */
   velocities_conform_bcs_6(E,E->V[1],E->V[2],E->V[3],E->V[4],E->V[5],E->V[6],E->mesh.levmax);
   fprintf(stderr,"OUT OF BCS HERE  \n") ;	
   time=CPU_time();
   vel_cycles_previous=E->control.total_iteration_cycles;
   vel_calls_previous=E->control.total_v_solver_calls;

   cycles=E->control.p_iterations;
 /* fprintf(stderr,"HERE IS TRACER.PT (0.5)   %g \n",E->tracer.Pt[436]) ;	 */
  if(E->control.verbose)
    fprintf(stderr,"Solve for q,v\n");

  /*fprintf(stderr,"Full_multigrid, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);*/

  residual_ddash=solve_Ahat_q_fhat(E,1.0);

  /* And convert any transformed boundary nodes back to 
     regular cartesian domain */
 
  if(3==dims) 
    for(i=1;i<=nno;i++) {
      if(E->node[i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[E->mesh.levmax][i];
	work1 = node_R[0*dims+0] * E->V[1][i] + node_R[0*dims+1] * E->V[2][i] + node_R[0*dims+2] * E->V[3][i];
	work2 = node_R[1*dims+0] * E->V[1][i] + node_R[1*dims+1] * E->V[2][i] + node_R[1*dims+2] * E->V[3][i];
	work3 = node_R[2*dims+0] * E->V[1][i] + node_R[2*dims+1] * E->V[2][i] + node_R[2*dims+2] * E->V[3][i];
	E->V[1][i] = work1 ;
	E->V[2][i] = work2 ;
	E->V[3][i] = work3 ;
	if(6==dofs) {
	  work1 = node_R[0*dims+0] * E->V[4][i] + node_R[0*dims+1] * E->V[5][i] + node_R[0*dims+2] * E->V[6][i];
	  work2 = node_R[1*dims+0] * E->V[4][i] + node_R[1*dims+1] * E->V[5][i] + node_R[1*dims+2] * E->V[6][i];
	  work3 = node_R[2*dims+0] * E->V[4][i] + node_R[2*dims+1] * E->V[5][i] + node_R[2*dims+2] * E->V[6][i];
	  E->V[4][i] = work1 ;
	  E->V[5][i] = work2 ;
	  E->V[6][i] = work3 ;
	}
      }
    }
  else if(2==dims) 
    for(i=1;i<=E->mesh.nno;i++) {
      if(E->node[i] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[E->mesh.levmax][i];
	work1 = node_R[0*dims+0] * E->V[1][i] + node_R[0*dims+1] * E->V[2][i] ;
	work2 = node_R[1*dims+0] * E->V[1][i] + node_R[1*dims+1] * E->V[2][i] ;
     	 E->V[1][i] = work1 ;
	 E->V[2][i] = work2 ;
      }
    }
  
  /* ... and signal this conversion to the rest of the program */
   E->control.CURRENT_SKEWED_Vs = 0;
  fprintf(stderr,"FM a1:  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	 
  if(E->control.verbose)
    fprintf(stderr,"Solve for q,v ...  done\n");

  sp_to_nodes(E,E->Q,E->nQ,E->mesh.levmax);
  fprintf(stderr,"FM a2:  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	
     remove_surf_horiz_press_ave(E,E->nQ,E->mesh.levmax);
  /*   * C.O'N: this call completely stuffs NQ and the elastic benchmark
   *   - don't want to fix surface pressure if no grav
   *   Take this call out if doing elastic benchmark*/
  fprintf(stderr,"FM a3:  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	
  been_here=1;
  return; 
}

standard_precision solve_Ahat_q_fhat (
				      struct All_variables *E,
				      standard_precision imp
				      )
{ 
  int i,j,count,convergent,valid,problems;
  int level;

  standard_precision *q0,*q1;

  standard_precision *UU1,*UU2,*UU3,*UU4,*UU5,*UU6 ;
  standard_precision *UUU1,*UUU2,*UUU3,*UUU4,*UUU5,*UUU6;

  standard_precision *Ah1,*Ah2,*Ah3,*Au1,*Au2,*Au3,*Au4,*Au5,*Au6,*R1,*R2,*R3,*DU1,*DU2,*DU3,*r1;
  standard_precision *AAu1,*AAu2,*AAu3,*AAu4,*AAu5,*AAu6;
  standard_precision *Ahv1,*Ahv2,*Ahv3;
  higher_precision r0dotr0;
  higher_precision residual, initial_residual,v_res;

  standard_precision vsselfdot(),lpselfdot(),vsselfdot6();

  static int been_here=0,convergence_problem=0;
  static int iteration_limit[MAX_LEVELS];

  standard_precision time,CPU_time();

  void interpolation_4pt_6();
  void interpolate_q();
  
  int solve_eta_v_q_3();

  const int dofs=E->mesh.dof ;
  const int dims=E->mesh.nsd ;

 
  q0 = (standard_precision *)Malloc0((E->mesh.npno+1)*sizeof(standard_precision));
  q1 = (standard_precision *)Malloc0((E->mesh.npno+1)*sizeof(standard_precision));
 
  r1  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
 
  Ah1  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Ah2  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Ah3  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision)); 
  
  Ahv1  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Ahv2  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Ahv3  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision)); 

  Au1  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Au2  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Au3  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Au4  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Au5  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  Au6  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));

  UU1  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UU2  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UU3  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UU4  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UU5  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UU6  = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));

  UUU1 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UUU2 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UUU3 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UUU4 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UUU5 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  UUU6 = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  
  /* Initialize */

  time=CPU_time();

 if(been_here==0) {
    switch(dofs) {
    case 2:
      if(E->mesh.periodic_x) { 
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
      }
      break;
    case 3:
      if(E->mesh.periodic_x || E->mesh.periodic_y) {
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
      }
      break;
    case 6:
      if(E->mesh.periodic_x) {
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[4],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[5],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[6],E->mesh.levmax);
      }
    }

    for(level=E->mesh.levmax;level>=E->mesh.q_levmin;level--) {
      inject06(E,level,
	       E->VV[level][1],E->VV[level][2],E->VV[level][3],  /* Uses only as many as needed */
	       E->VV[level][4],E->VV[level][5],E->VV[level][6],
	       E->VV[level-1][1],E->VV[level-1][2],E->VV[level-1][3], 
	       E->VV[level-1][4],E->VV[level-1][5],E->VV[level-1][6]);

    }
 }

#if 0
  /* Ignore previous timestep !! */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
    for(i=1;i<=E->mesh.NNO[level];i++) {
      UU1[i] = E->VV[level][1][i]=0.0;
      UU2[i] = E->VV[level][2][i]=0.0;
      /* UU3[i] = E->VV[level][3][i]=0.0; /* if 3D */
    }

    for(i=1;i<=E->mesh.NPNO[level];i++) 
      q0[i] = E->QQ[level][i] = 0.0;
  }
#endif

  /* Store old value of velocity/pressure for comparison with new */

  for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
    UU1[i] = -E->VV[E->mesh.q_levmin][1][i];
    UU2[i] = -E->VV[E->mesh.q_levmin][2][i];
  }
  if(dofs==3) {
    for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
      UU3[i] = -E->VV[E->mesh.q_levmin][3][i];
    }
  }
  else if(dofs==6) {
    for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
      UU3[i] = -E->VV[E->mesh.q_levmin][3][i];
      UU4[i] = -E->VV[E->mesh.q_levmin][4][i];
      UU5[i] = -E->VV[E->mesh.q_levmin][5][i];
      UU6[i] = -E->VV[E->mesh.q_levmin][6][i];
    }
  }

  for(i=1;i<=E->mesh.NPNO[E->mesh.q_levmin];i++) {
    q0[i] = -E->QQ[E->mesh.q_levmin][i];
  }

  if(E->control.verbose)
     fprintf(stderr,"Solve Visc/Vel/Q\n");

  convergent=solve_eta_v_q_3(E,E->mesh.q_levmin);
 
  if(E->control.verbose)
     fprintf(stderr,"Solve Visc/Vel/Q ... done one cycle\n");


  /* Currently we assume that non-convergence is either not a problem
     (fixed in the subsequent timesteps), or too big a problem to fix
     by automated tinkering ... i.e. restarting with different parameters
     would be the best solution */

  for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
    UU1[i] += E->VV[E->mesh.q_levmin][1][i];
    UU2[i] += E->VV[E->mesh.q_levmin][2][i];
  }
    if(dofs==3) {
      for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
	UU3[i] += E->VV[E->mesh.q_levmin][3][i];
      }
    }
    else if(dofs==6) {
      for(i=1;i<=E->mesh.NNO[E->mesh.q_levmin];i++) {
	UU3[i] += E->VV[E->mesh.q_levmin][3][i];
	UU4[i] += E->VV[E->mesh.q_levmin][4][i];
	UU5[i] += E->VV[E->mesh.q_levmin][5][i];
	UU6[i] += E->VV[E->mesh.q_levmin][6][i];
      }
    }

    for(i=1;i<=E->mesh.NPNO[E->mesh.q_levmin];i++) {
      q0[i] += E->QQ[E->mesh.q_levmin][i];
    }
  

  if(E->control.verbose) 
     fprintf(stderr,"RAA: NEXT COMMENT, HAS THE MAIN LOOP BEEN ENTERED at level %d -> %d\n",level,level+1);

  /* Now main loop */
  for(level=E->mesh.q_levmin;level<E->mesh.levmax;level++) {
   
    if(E->control.verbose)
      fprintf(stderr,"RAA: WELCOME TO THE MAIN LOOP at level %d -> %d\n",level,level+1);

    if(E->control.verbose) 
      fprintf(stderr,"Solving at level %d -> %d\n",level,level+1);

    
    interpolate_q(E,q1,q0,level+1);

#if 1
     interpolation_6(E,UUU1,UUU2,UUU3,UUU4,UUU5,UUU6,
		     UU1,UU2,UU3,UU4,UU5,UU6,level+1);
#else
     interpolation_4pt_6(E,UUU1,UUU2,UUU3,UUU4,UUU5,UUU6,
			 UU1,UU2,UU3,UU4,UU5,UU6,level+1);
#endif

     if(E->control.verbose) 
      fprintf(stderr,"Interpolated %d -> %d\n",level,level+1);

     
     if(dofs==2) {
       for(i=1;i<=E->mesh.NNO[level+1];i++) {
	 UU1[i] = -E->VV[level+1][1][i]; 
	 UU2[i] = -E->VV[level+1][2][i]; 
	 E->VV[level+1][1][i] += UUU1[i];
	 E->VV[level+1][2][i] += UUU2[i];
       }
     }
     else if(dofs==3) {
       for(i=1;i<=E->mesh.NNO[level+1];i++) {
	 UU1[i] = -E->VV[level+1][1][i]; 
	 UU2[i] = -E->VV[level+1][2][i]; 
	 UU3[i] = -E->VV[level+1][3][i]; 
	 E->VV[level+1][1][i] += UUU1[i];
	 E->VV[level+1][2][i] += UUU2[i];
   	 E->VV[level+1][3][i] += UUU3[i];
       }
     }
     else if(dofs==6) {
       for(i=1;i<=E->mesh.NNO[level+1];i++) {
    	 UU1[i] = -E->VV[level+1][1][i]; 
	 UU2[i] = -E->VV[level+1][2][i]; 
	 UU3[i] = -E->VV[level+1][3][i]; 
	 UU4[i] = -E->VV[level+1][4][i]; 
	 UU5[i] = -E->VV[level+1][5][i]; 
	 UU6[i] = -E->VV[level+1][6][i]; 
	 E->VV[level+1][1][i] += UUU1[i];
	 E->VV[level+1][2][i] += UUU2[i];
   	 E->VV[level+1][3][i] += UUU3[i];
 	 E->VV[level+1][4][i] += UUU4[i];
	 E->VV[level+1][5][i] += UUU5[i];
   	 E->VV[level+1][6][i] += UUU6[i];
       }
     }

     velocities_conform_bcs_6(E,
			     E->VV[level+1][1],E->VV[level+1][2],  
			     E->VV[level+1][3],E->VV[level+1][4], 
			     E->VV[level+1][5],E->VV[level+1][6], 
			     level+1); 

     for(i=1;i<=E->mesh.NPNO[level+1];i++) {
       q0[i] = -E->QQ[level+1][i];
       E->QQ[level+1][i] += q1[i];
     }

    convergent=solve_eta_v_q_3(E,level+1);   

    /* Get total change in velocity/pressure relative
       to the previous estimate (includes interpolation
       from lower levels, and correction at this level */
 
   switch(dofs) {
    case 2:
      for(i=1;i<=E->mesh.NNO[level+1];i++) {
	UU1[i] += E->VV[level+1][1][i];
	UU2[i] += E->VV[level+1][2][i];
      }
      break ;
    case 3:
      for(i=1;i<=E->mesh.NNO[level+1];i++) {
  	UU1[i] += E->VV[level+1][1][i];
	UU2[i] += E->VV[level+1][2][i];
 	UU3[i] += E->VV[level+1][3][i];
      }
      break ;
    case 6:
      for(i=1;i<=E->mesh.NNO[level+1];i++) {
 	UU1[i] += E->VV[level+1][1][i];
	UU2[i] += E->VV[level+1][2][i];
 	UU3[i] += E->VV[level+1][3][i];
 	UU4[i] += E->VV[level+1][4][i];
	UU5[i] += E->VV[level+1][5][i];
 	UU6[i] += E->VV[level+1][6][i];
      }
      break;
    }

    for(i=1;i<=E->mesh.NPNO[level+1];i++) {
      q0[i] += E->QQ[level+1][i];
    }
  }

  switch(dofs) {
    case 2:
      if(E->mesh.periodic_x) {
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
      }
      break ;
    case 3:
      if(E->mesh.periodic_x || E->mesh.periodic_y) {
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
      }
      break ;
    case 6:
      if(E->mesh.periodic_x) {
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[4],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[5],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[6],E->mesh.levmax);
      }
      break ;
    }

  /* Get formal residual */

  if(E->control.vector_optimization)
    e_assemble_del2_u_6(E,E->V[1],E->V[2],E->V[3],E->V[4],E->V[5],E->V[6],
			Au1,Au2,Au3,Au4,Au5,Au6,E->mesh.levmax,1); 
  else
    n_assemble_del2_u_6(E,E->V[1],E->V[2],E->V[3],E->V[4],E->V[5],E->V[6],
			Au1,Au2,Au3,Au4,Au5,Au6,E->mesh.levmax,1);


  assemble_grad_qs3(E,E->QQ[E->mesh.levmax],Ah1,Ah2,Ah3,E->mesh.levmax);
  assemble_grad_qst3(E,E->tracer.DQ,Ahv1,Ahv2,Ahv3,E->mesh.levmax);

  if(3==dofs && 3==dims) 
    for(i=1;i<=E->mesh.nno;i++) {
      if(E->node[i] & ( OFFSIDE  )) /*RAA: 30/5/01, these two lines added for 3D periodic*/
          continue;
      UU1[i] = E->F[E->id[i].doff[1]] - Au1[i] - Ahv1[i] - Ah1[i];
      UU2[i] = E->F[E->id[i].doff[2]] - Au2[i] - Ahv2[i] - Ah2[i];
      UU3[i] = E->F[E->id[i].doff[3]] - Au3[i] - Ahv3[i] - Ah3[i];
    }
  else if(3==dofs && 2==dims) {
    for(i=1;i<=E->mesh.nno;i++) {
      UU1[i] = E->F[E->id[i].doff[1]] - Au1[i] - Ahv1[i] - Ah1[i];
      UU2[i] = E->F[E->id[i].doff[2]] - Au2[i] - Ahv2[i] - Ah2[i];
      UU3[i] = E->F[E->id[i].doff[3]] - Au3[i] ;
    }
  }
  else if(2==dofs)
    for(i=1;i<=E->mesh.nno;i++) { 
      if(E->node[i] & ( OFFSIDE  ))
	continue;	
      UU1[i] = E->F[E->id[i].doff[1]] - Au1[i] - Ahv1[i] - Ah1[i];
      UU2[i] = E->F[E->id[i].doff[2]] - Au2[i] - Ahv2[i] - Ah2[i];
    }
  else /*if (6==dofs)*/ {
    for(i=1;i<=E->mesh.nno;i++) {
      UU1[i] = E->F[E->id[i].doff[1]] - Au1[i] - Ahv1[i] - Ah1[i];
      UU2[i] = E->F[E->id[i].doff[2]] - Au2[i] - Ahv2[i] - Ah2[i];
      UU3[i] = E->F[E->id[i].doff[3]] - Au3[i] - Ahv3[i] - Ah3[i];
      UU4[i] = E->F[E->id[i].doff[4]] - Au4[i];
      UU5[i] = E->F[E->id[i].doff[5]] - Au5[i];
      UU6[i] = E->F[E->id[i].doff[6]] - Au6[i];
    }
  }

  /* assemble_div_us3(E,u1[E->mesh.levmax],u2[E->mesh.levmax],u3[E->mesh.levmax],r1,E->mesh.levmax);

  E->monitor.vdotv = vsselfdot6(E,u1[E->mesh.levmax],u2[E->mesh.levmax],u3[E->mesh.levmax],
				u4[E->mesh.levmax],u5[E->mesh.levmax],u6[E->mesh.levmax],E->mesh.levmax);

  residual = sqrt(lpselfdot(E,r1,E->mesh.levmax)*E->mesh.nno/(E->monitor.vdotv*E->mesh.npno));
  E->monitor.incompressibility = residual; */

  if (E->control.print_convergence)
    fprintf(stderr,"Velocity residual %.6e - after  %gs (%gs total) \n",
	    sqrt(vsselfdot6(E,UU1,UU2,UU3,UU4,UU5,UU6,E->mesh.levmax)/
		 max(vsselfdot(E,E->F,E->mesh.levmax),vsselfdot6(E,Au1,Au2,Au3,Au4,Au5,Au6,E->mesh.levmax))),
	    CPU_time()-time,CPU_time());
 
  fprintf(E->fp,"Vel residual %.6e - after  %gs (%gs total) \n",
	  sqrt(vsselfdot6(E,UU1,UU2,UU3,UU4,UU5,UU6,E->mesh.levmax)/
	       max(vsselfdot(E,E->F,E->mesh.levmax),vsselfdot6(E,Au1,Au2,Au3,Au4,Au5,Au6,E->mesh.levmax))),
	  CPU_time()-time,CPU_time());


 
  free((void *) q0);
  free((void *) q1);
  free((void *) r1);
  free((void *) Ah1);
  free((void *) Ah2);
  free((void *) Ah3);
  free((void *) Ahv1);
  free((void *) Ahv2);
  free((void *) Ahv3);
  free((void *) Au1);
  free((void *) Au2);
  free((void *) Au3);
  free((void *) Au4);
  free((void *) Au5);
  free((void *) Au6);
  free((void *) UU1);
  free((void *) UU2);
  free((void *) UU3);
  free((void *) UU4);
  free((void *) UU5);
  free((void *) UU6);
  free((void *) UUU1);
  free((void *) UUU2);
  free((void *) UUU3);
  free((void *) UUU4);
  free((void *) UUU5);
  free((void *) UUU6);

  been_here++;
 
  /* After one cycle, q_levmin -> levmax */

    E->mesh.q_levmin = E->mesh.levmax; /* */

  return;
}

int solve_eta_v_q_3 (
		     struct All_variables *E,
		     int level
		     )
{
  int i,j,k,m,convergent,max_iterations;
  int continuity_satisfied;
  standard_precision *R1,*R2,*R3,*R4,*R5,*R6,*Eta;
  standard_precision *u1,*u2,*u3,*u4,*u5,*u6,*uu1,*uu2,*uu3,*uu4,*uu5,*uu6,*sigmat1;
  standard_precision v_res,delta_u,delta_sigma;

  standard_precision under;
  standard_precision t1,ave;

  int maxnode;

  standard_precision vsselfdot(),vsselfdot6(),lpselfdot();
  standard_precision resid_level_corr;
 
  int solve_q_velocity_3();

  void s_strip_bcs_from_residual();
  void sp_to_nodes();
  void remove_surf_horiz_press_ave();
  void remove_horiz_press_ave();

  static int been_here=0;
  static int been_here_level[MAX_LEVELS];

  const int npno=E->mesh.NPNO[level];
  const int nno=E->mesh.NNO[level];
  const int dofs = E->mesh.dof;
  const int vpts = vpoints[E->mesh.nsd];
 
  R1   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R2   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R3   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision)); 
  R4   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R5   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R6   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision)); 
  u1   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u2   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u3   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u4   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u5   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u6   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu1   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu2   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu3   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu4   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu5   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  uu6   = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));

  sigmat1 = (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));


  resid_level_corr = pow(E->control.delta_accuracy_factor,(double)(E->mesh.levmax-level));
  convergent=1;
  max_iterations = 50*(E->mesh.levmax-level+1);
 fprintf(stderr,"In solve_eta_v_q_3\n") ;	
 /*fprintf(stderr,"0.  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	*/
  if(dofs==2)
    for(i=1;i<=E->mesh.NNO[level];i++) {  
	u1[i] = 0.0;
	u2[i] = 0.0;
    }
  else if(dofs==3)
    for(i=1;i<=E->mesh.NNO[level];i++) {
	u1[i] = 0.0;
	u2[i] = 0.0;
	u3[i] = 0.0;
    }
  else if(dofs==6)
    for(i=1;i<=E->mesh.NNO[level];i++) {
	u1[i] = 0.0;
	u2[i] = 0.0;
	u3[i] = 0.0;
	u4[i] = 0.0;
	u5[i] = 0.0;
	u6[i] = 0.0;
    }

  
  for(i=1;i<= E->tracer.NUM_TRACERS;i++)
    sigmat1[i] = 0.0;

  /* If non-linear rheology, boundary conditions or anything like that,
     add it in to make sure the iterations are done */

  if(/*E->advection.timesteps > 1 ||*/ E->viscosity.SDEPV || E->viscosity.YIELD)
    max_iterations = 50*(E->mesh.levmax-level+1);
  else
    max_iterations = 10*(E->mesh.levmax-level+1)/**/; 

  continuity_satisfied = 0;
  


  for(k=1;k<=max_iterations;k++) {
    delta_u = 0.0;

    if(1 || k==1 || E->viscosity.SDEPV || E->viscosity.YIELD) {
      
      fprintf(E->fp1,"############RAA: iter: %d  ###############################\n",k);

      get_viscosity_all_elements(E,level); 

      /* Change in viscosity means that any velocity bc's must be 
	 recalculated and therefore E->FF might be out of date */

      assemble_forces(E,k,level);      
      s_strip_bcs_from_residual(E,E->FF[level],level);
    }

    if(dofs==2) {
        for(i=1;i<=E->mesh.NNO[level];i++) {
	  R1[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[1]];
	  R2[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[2]];
	  uu1[i] = E->VV[level][1][i];
	  uu2[i] = E->VV[level][2][i];
      }
    }
    else if(dofs==3)
      for(i=1;i<=E->mesh.NNO[level];i++) {
	  R1[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[1]];
	  R2[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[2]];
	  R3[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[3]];
	  uu1[i] = E->VV[level][1][i];
	  uu2[i] = E->VV[level][2][i];
  	  uu3[i] = E->VV[level][3][i];
      }
    else if(dofs==6)
      for(i=1;i<=E->mesh.NNO[level];i++) {
	  R1[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[1]];
	  R2[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[2]];
	  R3[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[3]];
	  R4[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[4]];
	  R5[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[5]];
	  R6[i] = ((E->NODE[level][i] & ( OFFSIDE)) == 0) * E->FF[level][E->ID[level][i].doff[6]];
   	  uu1[i] = E->VV[level][1][i];
	  uu2[i] = E->VV[level][2][i];
  	  uu3[i] = E->VV[level][3][i];
 	  uu4[i] = E->VV[level][4][i];
	  uu5[i] = E->VV[level][5][i];
  	  uu6[i] = E->VV[level][6][i];
   }
    
   /* fprintf(stderr,"HERE IS TRACER.PT (1)   %g  NQ %g \n",E->tracer.Pt[436],E->NQ[1][31]) ;	*/
    if(E->control.verbose)
      fprintf(stderr,"Solve Q & V (1) at level %d \n",level);
  
 /*RAA: check out the velocities */
 /*    for(i=1;i<=E->mesh.nno;i++) 
    fprintf(E->fp1,"-->VV velocities (before solve_q_v_3): node#, VV1, VV2, VV3 &iter#: %d %g %g %g %d\n",i,E->VV[level][1][i],E->VV[level][2][i],E->VV[level][3][i],k); */


 /*RAA: 11/5/01, check out the pressures */
 /*    for(i=1;i<=E->mesh.NPNO[level];i++) 
    fprintf(E->fp1,"CHECK OUT PRESSURES before 'solve_q_velocity_3': level npno QQ %d %d %g\n",level,i,E->QQ[level][i]);   */
    
    /*fprintf(stderr,"1.  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	*/
    
    continuity_satisfied = solve_q_velocity_3(E,E->VV[level][1],E->VV[level][2],
					      E->VV[level][3],E->VV[level][4],
					      E->VV[level][5],E->VV[level][6],
					      E->QQ[level],
					      R1,R2,R3,R4,R5,R6,1.0,k,level);
   

    /* Pressure to nodes */

    /*fprintf(stderr,"2.  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	*/
    sp_to_nodes(E,E->QQ[level],E->NQ[level],level);	
    remove_surf_horiz_press_ave(E,E->NQ[level],level);
    nodes_to_tracers(E,E->NQ[level],E->tracer.DQ1,level); /* This'll do ! */
  
   /* fprintf(stderr,"3.  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	
    fprintf(stderr,"HERE IS TRACER.PT (2)   %g NQ %g \n",E->tracer.Pt[436],E->NQ[1][31]) ;	*/
    if(E->control.verbose)
      fprintf(stderr,"Solved Q & V - one cycle\n");

    if(k>5) {
      under = 0.95;
    } 
    else {
      under = 1.0;
    }
  
    /* Shouldn't Q be underrelaxed too ? */

    delta_u=0.0; 
    if(3==dofs) {
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(!(E->NODE[level][i] & OFFSIDE)) {
	  delta_u += 
	    (uu1[i] - E->VV[level][1][i]) * (uu1[i] - E->VV[level][1][i]) + 
	    (uu2[i] - E->VV[level][2][i]) * (uu2[i] - E->VV[level][2][i]) + 
	    (uu3[i] - E->VV[level][3][i]) * (uu3[i] - E->VV[level][3][i]) ;
	  uu1[i] = E->VV[level][1][i] = under * E->VV[level][1][i] + (1.0-under) * uu1[i];
	  uu2[i] = E->VV[level][2][i] = under * E->VV[level][2][i] + (1.0-under) * uu2[i];
	  uu3[i] = E->VV[level][3][i] = under * E->VV[level][3][i] + (1.0-under) * uu3[i];
	}
      }
    }
    else if (2==dofs) {
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(!(E->NODE[level][i] & OFFSIDE)) {
	   delta_u += 
	    (uu1[i] - E->VV[level][1][i]) * (uu1[i] - E->VV[level][1][i]) + 
	    (uu2[i] - E->VV[level][2][i]) * (uu2[i] - E->VV[level][2][i]) ;
	   uu1[i] = E->VV[level][1][i] = under * E->VV[level][1][i] + (1.0-under) * uu1[i];
	   uu2[i] = E->VV[level][2][i] = under * E->VV[level][2][i] + (1.0-under) * uu2[i];
	}
      }
    }
    else /* if (6==dofs) */ {
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(!(E->NODE[level][i] & OFFSIDE)) {
	  delta_u +=  
	    (uu1[i] - E->VV[level][1][i]) * (uu1[i] - E->VV[level][1][i]) + 
	    (uu2[i] - E->VV[level][2][i]) * (uu2[i] - E->VV[level][2][i]) + 
	    (uu3[i] - E->VV[level][3][i]) * (uu3[i] - E->VV[level][3][i]) +
	    (uu4[i] - E->VV[level][4][i]) * (uu4[i] - E->VV[level][4][i]) + 
	    (uu5[i] - E->VV[level][5][i]) * (uu5[i] - E->VV[level][5][i]) + 
	    (uu6[i] - E->VV[level][6][i]) * (uu6[i] - E->VV[level][6][i]) ;
	  uu1[i] = E->VV[level][1][i];
	  uu2[i] = E->VV[level][2][i];
	  uu3[i] = E->VV[level][3][i];
	  uu4[i] = E->VV[level][4][i];
	  uu5[i] = E->VV[level][5][i];
	  uu6[i] = E->VV[level][6][i];
	}
      }
    }

    t1 = vsselfdot6(E,E->VV[level][1],E->VV[level][2],E->VV[level][3],
                    E->VV[level][4],E->VV[level][5],E->VV[level][6],level); 
  /*RAA: 9/4/01, Guard against division by zero, same as C. Wijns version*/
    if(t1==0.0) {
      /*RAA: output this problem (& correction) to screen */
      fprintf(stderr,"Denominator t1 is zero in solve_eta_v_p_3: delta_u corrected. \n");
      if(delta_u!=0.0) 
         delta_u = 1.0e32;
    }
    else 
      delta_u = sqrt(delta_u/t1);
 
/*    delta_u = sqrt(delta_u/vsselfdot6(E,
				      E->VV[level][1],E->VV[level][2],E->VV[level][3],
				      E->VV[level][4],E->VV[level][5],E->VV[level][6],level)); */


    if(E->control.print_convergence) {
      fprintf(stderr,"Level %d ---> delta u(%d) = %g Continuity satisfied - %s\n",
	      level,k,delta_u,(continuity_satisfied ? "Yes" : "no")); 
    }

    if(0 && level==E->mesh.levmax) {
      tracer_corrector(E,
		       E->VV[level][1],E->VV[level][2],E->VV[level][3],
		       E->VV[level][4],E->VV[level][5],E->VV[level][6],
		       E->advection.timestep);  
      get_tracer_elts_and_weights(E,1,E->tracer.NUM_TRACERS); 
    }

    /* TEST for completion */

     if((k != 1) && /* i.e. didn't do anything yet ! */
	((delta_u <= E->control.accuracy * resid_level_corr)) &&
	/* (delta_visc <= E->control.accuracy * resid_level_corr) && */
	continuity_satisfied ) {
       convergent = 1;
       break; 
    }

     fprintf(E->fp,"Level %d ---> delta u(%d) = %g Continuity satisfied - %s\n",
	     level,k-1,delta_u,(continuity_satisfied ? "Yes" : "no"));

  }

  if(max_iterations != 1 && k>=max_iterations) {
    /* issue a warning */ 
    convergent=0;
    fprintf(E->fp,"Level %d ---> delta u(%d) = %g \n",level,k,delta_u); 
    fprintf(E->fp,"V,Q,eta - Non convergence at level %d (%d iterations)\n",level,max_iterations);
    fflush(E->fp); 
  }

  fprintf(E->fp,"Level %d ---> delta u(%d) = %g Continuity satisfied - %s\n",
	  level,k-1,delta_u,(continuity_satisfied ? "Yes" : "no"));
  fprintf(stderr,"FM 4:  QQ  %g NQ %g \n",E->QQ[1][31],E->NQ[1][31]) ;	
 
  free((void *) u1);
  free((void *) u2);
  free((void *) u3);
  free((void *) u4);
  free((void *) u5);
  free((void *) u6);
  free((void *) uu1);
  free((void *) uu2);
  free((void *) uu3);
  free((void *) uu4);
  free((void *) uu5);
  free((void *) uu6);
  free((void *) R1);
  free((void *) R2);
  free((void *) R3);
  free((void *) R4);
  free((void *) R5);
  free((void *) R6);
  /*RAA: 01/07/02, added 'free' below for sigmat1 to avoid a memory leak */
  free((void *) sigmat1);

  return(convergent);
}
