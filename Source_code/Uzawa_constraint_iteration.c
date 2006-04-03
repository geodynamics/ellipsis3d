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


/*   Functions which solve for the velocity and
     pressure fields using Uzawa-type iteration loop.  */


#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

/* =========================================================================
   Solves the Au + Bp = f equation of motion assuming no initial guess for p 
   but that V contains helpful best-guess info of velocity 
   ========================================================================= */

int solve_q_velocity_3(
			      struct All_variables *E,
			      standard_precision *V1,
			      standard_precision *V2,
			      standard_precision *V3,
			      standard_precision *V4,
			      standard_precision *V5,
			      standard_precision *V6,
			      standard_precision *Q,
			      standard_precision *F1,
			      standard_precision *F2,
			      standard_precision *F3,
			      standard_precision *F4,
			      standard_precision *F5,
			      standard_precision *F6,
			      standard_precision v_res,
			      int V_GUESS,
			      int level
			      )
{
  int i,j,m,count,loops;
  int steps_max,relax;
  int continuity_satisfied;

  standard_precision *u1,*u2,*u3,*u4,*u5,*u6;
  standard_precision *Ah1,*Ah2,*Ah3,*R1,*R2,*R3,*R4,*R5,*R6,*Au1,*Au2,*Au3;
  
  standard_precision *Ahv1,*Ahv2,*Ahv3,*QQ;
  standard_precision alpha,delta,dotAhat;
 
  standard_precision imp_fact,this_res,v_mag,q_res;
  standard_precision time;
  standard_precision this_imp;

  standard_precision continuity;
  
  static int been_here=0;
  static int solved_this_level[MAX_LEVELS];

  void assemble_div_us3();
  void s_strip_bcs_from_residual();
  /* void remove_horiz_press_ave();
  void remove_nodal_horiz_press_ave();
  void remove_surf_horiz_press_ave(); */  /*RAA: 17/4/01, these 3 functions not called*/
  
  standard_precision return_bulk_value_l();
  standard_precision CPU_time();
  standard_precision lpdot(),vsdot(),vsselfdot(),vsselfdot6(),lpselfdot(),return_bulk_value_l();
  /* int solve_del2_u(); */ /*RAA: old function, not called*/
  int apply_continuity();

  const int npno=E->mesh.NPNO[level];
  const int nno=E->mesh.NNO[level];
  const int neq=E->mesh.NEQ[level];
  const int dofs = E->mesh.dof;
  const int dims = E->mesh.nsd;

  int count1,count2;

  standard_precision aveVx ;

  count1 = count2 = 0 ;


  if(been_here++ == 0) {
    for(i=E->mesh.q_levmin;i<=E->mesh.levmax;i++) 
      solved_this_level[i]=0;
  }

  QQ =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));

  u1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u4  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u5  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u6  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  
  R1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R4  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R5  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R6  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));

  Au1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Au2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Au3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));

  Ah1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ah2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ah3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));

  Ahv1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ahv2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ahv3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  
  steps_max=E->control.p_iterations;
  this_imp = pow(E->control.delta_accuracy_factor,(double)(E->mesh.levmax-level));
  if(this_imp < 0.01)
    this_imp = 0.01;
  if(this_imp > 10.0)
    this_imp=10.0;
  
  time=CPU_time();

  if(E->control.vector_optimization)
    e_assemble_del2_u_6(E,V1,V2,V3,V4,V5,V6,u1,u2,u3,u4,u5,u6,level,1);
  else 
    n_assemble_del2_u_6(E,V1,V2,V3,V4,V5,V6,u1,u2,u3,u4,u5,u6,level,1);


/*RAA: 28/11/01, check out the load vector*/
  /* for(i=1;i<=nno;i++) {
       fprintf(stderr,"-->FF's at node#, F1,F2,F3: %d %g %g %g\n",i,F1[i],F2[i],F3[i]); 
     }
  */
  
  assemble_grad_qst3(E,E->tracer.DQ,Ah1,Ah2,Ah3,level);
  assemble_grad_qs3(E,Q,Ahv1,Ahv2,Ahv3,level);

  if(2==dofs) {
    for(i=1;i<=nno;i++) {
      R1[i] = (F1[i] - u1[i] - Ah1[i] - Ahv1[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R2[i] = (F2[i] - u2[i] - Ah2[i] - Ahv2[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
    }
  }
  else if(3==dofs && 3==dims) {
    for(i=1;i<=nno;i++) {
      R1[i] = (F1[i] - u1[i] - Ah1[i] - Ahv1[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R2[i] = (F2[i] - u2[i] - Ah2[i] - Ahv2[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R3[i] = (F3[i] - u3[i] - Ah3[i] - Ahv3[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
    }    
  }
  else if(3==dofs && 2==dims) {
    for(i=1;i<=nno;i++) {
      R1[i] = (F1[i] - u1[i] - Ah1[i] - Ahv1[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R2[i] = (F2[i] - u2[i] - Ah2[i] - Ahv2[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R3[i] = (F3[i] - u3[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0) ;
    }
  }
  else {
    for(i=1;i<=nno;i++) {
      R1[i] = (F1[i] - u1[i] - Ah1[i] - Ahv1[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R2[i] = (F2[i] - u2[i] - Ah2[i] - Ahv2[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R3[i] = (F3[i] - u3[i] - Ah3[i] - Ahv3[i]) * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R4[i] = F4[i] - u4[i] * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R5[i] = F5[i] - u5[i] * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
      R6[i] = F6[i] - u6[i] * ((E->NODE[level][i] & ( OFFSIDE)) == 0);
    }    
  }
  
  strip_bcs_from_residual_6(E,R1,R2,R3,R4,R5,R6,level);
  v_res=sqrt(vsselfdot6(E,F1,F2,F3,F4,F5,F6,level)/nno);

/*
  for(i=1;i<=E->mesh.nno;i++) {
     if(3==dims) {
        fprintf(E->fp1,"b4 multisolve:  node  V1[n] V2[n] V3[n]  %d  %g  %g  %g \n",i,V1[i],V2[i],V3[i]); 
        fprintf(E->fp1,"b4 multisolve:  node  u1[n] u2[n] u3[n]  %d  %g  %g  %g \n",i,u1[i],u2[i],u3[i]);
        fprintf(E->fp1,"b4 multisolve:  node  F1[n] F2[n] F3[n]  %d  %g  %g  %g \n",i,F1[i],F2[i],F3[i]); 
     }
    else
        fprintf(E->fp1,"b4 multisolve:  node  F1[n] F2[n] %d   %g  %g \n",i,F1[i],F2[i]);
    } */

/* fprintf(stderr,"Vres = %g (%g)\n",v_res,vsselfdot6(E,F1,F2,F3,F4,F5,F6,level)); */

  if(1 || sqrt(vsselfdot6(E,R1,R2,R3,R4,R5,R6,level)/nno)  > 0.1 * E->control.accuracy * v_res * this_imp) {
    /* 1. Initial improvement of velocity residual */
    multigrid_solve_del2_u(E,u1,u2,u3,u4,u5,u6,R1,R2,R3,R4,R5,R6,0.1*E->control.accuracy*v_res* this_imp,level,1);
    
    v_mag = sqrt(vsselfdot6(E,V1,V2,V3,V4,V5,V6,level)/nno);
    if(v_mag == 0.0) {
      v_mag = sqrt(vsselfdot6(E,u1,u2,u3,u4,u5,u6,level)/nno);
    }

    switch(dofs) {
    case 2:
      for(i=1;i<=nno;i++) {
	if(E->NODE[level][i] & ( OFFSIDE  ))
	  continue;
	V1[i] +=  u1[i];
	V2[i] +=  u2[i];
      }
      break ;
    case 3:
      for(i=1;i<=nno;i++) {
	if(E->NODE[level][i] & ( OFFSIDE  ))
	 continue;
       V1[i] += u1[i];
       V2[i] += u2[i];
       V3[i] += u3[i];
     }    
     break ;
   case 6:
     for(i=1;i<=nno;i++) {
       if(E->NODE[level][i] & ( OFFSIDE  ))
	 continue;
       V1[i] += u1[i];
       V2[i] += u2[i];
       V3[i] += u3[i];
       V4[i] += u4[i];
       V5[i] += u5[i];
       V6[i] += u6[i];
     }
     break ;
   }
   velocity_conform_bcs_6(E,V1,V2,V3,V4,V5,V6,level);
 }
 v_mag = sqrt(vsselfdot6(E,V1,V2,V3,V4,V5,V6,level)/nno);

/* fprintf(stderr,"Vmag = %g, Vres = %g\n",v_mag,v_res); */

 /* 2 magnitude of pressure equation terms 
    (approx - just one V cycle to get an estimate. Otherwise 
    this can be really hard to solve at certain times) */

 multigrid_solve_del2_u(E,u1,u2,u3,u4,u5,u6,F1,F2,F3,F4,F5,F6,1000.0*E->control.accuracy*v_res,level,1);
 assemble_div_us3(E,u1,u2,u3,QQ,level);

 q_res = 0.0;
 for(i=1;i<=E->mesh.NPNO[level];i++) {
   q_res += QQ[i] * QQ[i];
 }
 q_res = sqrt(q_res / E->mesh.NPNO[level]);

 /* fprintf(stderr,"%d: |Fhat| = %g\n",level,p_res); */
 
 /* Solve CG form of Uzawa's iteration to impose continuity */
 continuity_satisfied = apply_continuity(E,V1,V2,V3,Q,v_mag,v_res,q_res,level);
 velocity_conform_bcs_6(E,V1,V2,V3,NULL,NULL,NULL,level);

/* continuity_satisfied=1 ;/**/
 for(m=1;m<=E->tracer.NUM_TRACERS;m++)
   E->tracer.Q[m] = E->tracer.DQ1[m] /* + E->tracer.DQ[m]*/;

 free((void *) u1);
 free((void *) u2);
 free((void *) u3);
 free((void *) u4);
 free((void *) u5);
 free((void *) u6);
 free((void *) R1); 
 free((void *) R2); 
 free((void *) R3); 
 free((void *) R4); 
 free((void *) R5); 
 free((void *) R6); 
 free((void *) Ah1);
 free((void *) Ahv1);
 free((void *) Au1);
 free((void *) Ah2);
 free((void *) Ahv2);
 free((void *) Au2);
 free((void *) Ah3);
 free((void *) Ahv3);
 free((void *) Au3);
 
 free((void *) QQ);

 printf("In Uzawa, first free finished\n");

 return(continuity_satisfied);
}

/* =========================================================================
   Solves the Au + Bp = f equation of motion assuming no initial guess for p 
   but that V contains helpful best-guess info of velocity 
   ========================================================================= */

int apply_continuity(
		     struct All_variables *E,
		     standard_precision *V1,
		     standard_precision *V2,
		     standard_precision *V3,
		     standard_precision *Q,
		     standard_precision v_mag,
		     standard_precision v_res,
		     standard_precision q_res,
		     int level
		     )
{
  int i,j,m,count,loops;
  int steps_max,relax;

  standard_precision *r0,*r1,*r2;
  standard_precision *u1,*u2,*u3;

  standard_precision *Ah1,*Ah2,*Ah3,*R1,*R2,*R3;
  standard_precision *QQ,*ZZ1,*ZZ2,*QQ1,*w1;

  standard_precision alpha,delta,dotAhat,r0dotr0,r1dotr1;
  standard_precision imp_fact,this_res,v_res2;
  standard_precision time;
  standard_precision this_imp;
  standard_precision continuity;
 
  static int been_here=0;

  void assemble_div_us3();
  standard_precision lpdot(),vsdot(),vsselfdot(),vsselfdot3(),vsselfdot6(),lpselfdot();
  standard_precision CPU_time();

  /* int solve_del2_u(); */ /*RAA: old function, not called*/
   
  const int npno=E->mesh.NPNO[level];
  const int nno=E->mesh.NNO[level];
  const int neq=E->mesh.NEQ[level];
  const int dofs = E->mesh.dof;
  const int dims = E->mesh.nsd;

 
  this_imp = pow(E->control.delta_accuracy_factor,(double)(E->mesh.levmax-level));
  if(this_imp < 0.01)
    this_imp = 0.01;
  if(this_imp > 10.0)
    this_imp=10.0;


/*  QQ =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  ZZ1 =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  QQ1 =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  r1 =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  w1 =  (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));*/
 
  QQ =  (standard_precision *)Malloc0((npno+1)*sizeof(standard_precision));
  ZZ1 =  (standard_precision *)Malloc0((npno+1)*sizeof(standard_precision));
  QQ1 =  (standard_precision *)Malloc0((npno+1)*sizeof(standard_precision));
  r1 =  (standard_precision *)Malloc0((npno+1)*sizeof(standard_precision));
  w1 =  (standard_precision *)Malloc0((npno+1)*sizeof(standard_precision));
 
  u1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  u3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));  

  R1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  R3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));

  Ah1  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ah2  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  Ah3  = (standard_precision *)Malloc0((nno+1)*sizeof(standard_precision));
  
  steps_max=E->control.p_iterations;
  time=CPU_time();

  /* Preconditioner */

  assemble_div_us3(E,V1,V2,V3,r1,level); /* Compute Bt.U */
  assemble_Mq_t(E,Q,QQ,level); /* Compute QQ = M.q */

  for(i=1;i<=E->mesh.NPNO[level];i++) {
    r1[i] -= QQ[i]; /* First residue r = Bt.U + M.q */
    ZZ1[i] = r1[i] * E->BPI[level][i];  /* preconditionner M.z = r with M=(Ahat-1) */
  }

  count=0;

  do {
    if(count==0) {
      r0dotr0 = 0.0;
      for(i=1;i<=E->mesh.NPNO[level];i++) {
 	QQ1[i] = ZZ1[i];
	r0dotr0 += r1[i] * ZZ1[i];
      }
    }
    else {
      r1dotr1=0.0;
      for(i=1;i<=E->mesh.NPNO[level];i++) 
	r1dotr1 += r1[i] * ZZ1[i];

      delta = r1dotr1 / r0dotr0; 
      r0dotr0 = r1dotr1;

      for(i=1;i<=E->mesh.NPNO[level];i++) 
	QQ1[i] = ZZ1[i] + delta * QQ1[i];
 
    }

    /* Velocity search direction */
    assemble_grad_qs3(E,QQ1,Ah1,Ah2,Ah3,level);
    v_res=sqrt(vsselfdot3(E,Ah1,Ah2,Ah3,level)/nno);
    multigrid_solve_del2_u(E,u1,u2,u3,0,0,0,Ah1,Ah2,Ah3,0,0,0,0.5 * E->control.accuracy*v_res*this_imp,level,2);
    strip_bcs_from_residual_6(E,u1,u2,u3,0,0,0,level);

    /* Gives pressure search direction */

    assemble_div_us3(E,u1,u2,u3,w1,level); /* Compute Bt.u or Bt.A-1.B.p */
    assemble_Mq_t(E,QQ1,QQ,level); /* Compute PP = M.QQ1 */
    
    dotAhat=0.0;
    for(i=1;i<=E->mesh.NPNO[level];i++) {
      w1[i] += QQ[i]; /* Bt.A-1.B.p - M.p or Ahat.p*/
      dotAhat += QQ1[i] * w1[i];
    }
    
    if(dotAhat == 0.0)  
      alpha = 0.0;
     else
      alpha = r0dotr0/dotAhat; /* rT.z / p.Ahat.p */
  
    for(i=1;i<=E->mesh.NPNO[level];i++) {
      Q[i]  += alpha * QQ1[i];
      r1[i] -= alpha * w1[i]; /* r = r - alpha.Ahat.p */
      ZZ1[i] = r1[i] * E->BPI[level][i];
    }
    if(3==dims)
      for(i=1;i<=nno;i++) {
	if(E->NODE[level][i] & OFFSIDE  ) /*RAA: 30/5/01, these 2 lines added for 3D periodic */
	  continue;
	V1[i] -= alpha * u1[i];
   	V2[i] -= alpha * u2[i];
   	V3[i] -= alpha * u3[i];
      }
    else
      for(i=1;i<=nno;i++) {
	if(E->NODE[level][i] & OFFSIDE  )
	  continue;
	V1[i] -= alpha *  u1[i];
   	V2[i] -= alpha *  u2[i];
      }

    v_res2 = 0.0;
    for(i=1;i<=E->mesh.NPNO[level];i++) {
      v_res2 += r1[i] * r1[i];
    }

    /* This definition of the continuity residual is prefered since it is stable
     and simplifies to divU/U for incompressible version of the algorithm */

    continuity = sqrt(v_res2 / E->mesh.NPNO[level]);

    if(E->control.print_convergence)
      fprintf(stderr,"%d: Continuity Residual = %g (relative) = %g\n",count,continuity, continuity / 
	      q_res);

   count++;
  } while ((count <E->control.p_iterations) && (/* (pressure_change > E->control.accuracy) ||*/ 
			     /*(continuity > 0.001 * E->control.accuracy + 0.0 * continuity_residual)*/
			     (continuity > 0.5 * q_res * E->control.accuracy * this_imp) ));

  fprintf(E->fp,"%d: Continuity Residual = %g (relative) = %g\n",count,continuity, continuity / 
	      q_res);
  
  free((void *) u1);
  free((void *) Ah1);
  free((void *) R1); 
  free((void *) u2);
  free((void *) Ah2);
  free((void *) R2); 
  free((void *) u3);
  free((void *) Ah3);

  free((void *) R3); 
  free((void *) QQ); 
  free((void *) ZZ1); 
  free((void *) QQ1); 
  free((void *) r1); 
  free((void *) w1); 

  printf("In Uzawa, second free done \n");
 
   return((continuity / q_res) < E->control.accuracy * 
			      pow(E->control.delta_accuracy_factor,(double)(-E->mesh.levmax+level)) );
}
