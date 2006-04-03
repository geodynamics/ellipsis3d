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
/* #include <malloc.h> */
/* #include <rpc/rpc.h> */

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


/* =================================================
   tracer_advection: Midpoint method (2nd order RK)
   =================================================  */

void passive_tracer_advection(
			      struct All_variables *E,
			      standard_precision timestep
			      )
{
  /* Use the general tracer advection routine to

     1) Apply a corrector to the previous particle locations
     if this is not the first timestep 

     2) Shoot the current particle positions forward using the
     latest velocity solution (predictor)

  */
  int n;
  printf("passive: 1 doing general_tracer_advection \n");
  void general_tracer_advection();
  printf("passive: 2 doing standard_tracer_advection \n");
  void standard_tracer_advection();

  /* Corrector step */

  if(0 && E->advection.timesteps > 1)  { 
    /* This is called after the first timestep value has been calculated */
    general_tracer_advection(E,E->advection.previous_timestep,
			     E->tracer.sample_x1,E->tracer.sample_z1,E->tracer.sample_y1,
			     E->tracer.sample_x, E->tracer.sample_z, E->tracer.sample_y ,
			     E->V1[1], E->V1[2], E->V1[3],
			     E->V[1], E->V[2], E->V[3],
			     0,E->tracer.SAMPLE_PTS,
			     E->tracer.sample_in_element,E->tracer.sample_lagrangian);
  }
  /* Updated tracer positions should now become the stored positions
     for the next correction */
  
  for(n=0;n<E->tracer.SAMPLE_PTS;n++) {
    E->tracer.sample_x1[n] = E->tracer.sample_x[n];
    E->tracer.sample_z1[n] = E->tracer.sample_z[n];
    if(E->mesh.nsd==3)
      E->tracer.sample_y1[n] = E->tracer.sample_y[n];
  }

  /* Predictor step */
  printf("passive: 3 doing general_tracer_advection \n");
  general_tracer_advection(E,E->advection.timestep,
			   E->tracer.sample_x1,E->tracer.sample_z1,E->tracer.sample_y1,
			   E->tracer.sample_x, E->tracer.sample_z, E->tracer.sample_y ,
			   E->V[1], E->V[2], E->V[3], 
			   E->V[1], E->V[2], E->V[3],
			   0,E->tracer.SAMPLE_PTS,
			   E->tracer.sample_in_element,E->tracer.sample_lagrangian);

  printf("passive: Done!!! \n");

  return;
}

/* =================================================
   tracer_advection: Midpoint method (2nd order RK)
   =================================================  */

void tracer_advection(
		      struct All_variables *E,
		      standard_precision timestep
		      )
{
  /* Use the general tracer advection routine to
     
     1) Apply a corrector to the previous particle locations
     if this is not the first timestep 
     
     2) Shoot the current particle positions forward using the
     latest velocity solution (predictor)
  */
  
  standard_precision *XX,*ZZ,*XX1,*ZZ1;
  standard_precision *dX1,*dX2,*dX3;

  int *el_list;
  const int dims = E->mesh.nsd;
  int n,m,tr;

  static int here = 0;

  void standard_tracer_advection();
  void general_tracer_advection();

  dX1 = (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  dX2 = (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  dX3 = (standard_precision *)Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(standard_precision));
  
  el_list = (int *) Malloc0((E->tracer.NUM_TRACERS+1)*sizeof(int));
  printf("tracer_adv: 1 just alloced memory \n");

 

 /* Store current particle positions */

   for(n=1;n<=E->tracer.NUM_TRACERS;n++) {
     E->tracer.tx1[n] = E->tracer.tx[n];
     E->tracer.tz1[n] = E->tracer.tz[n];
     if(E->mesh.nsd==3)
       E->tracer.ty1[n] = E->tracer.ty[n];
   }


  if(E->control.verbose)
    fprintf(stderr,"Particle advection\n");


  /* Predictor step */


 if(E->control.verbose)
    fprintf(stderr,"Particle advection ... 1\n");

  printf("tracer_adv: 2 doing standard_tracer_advection \n");

  standard_tracer_advection(E,E->advection.timestep,
			   E->tracer.tx1,E->tracer.tz1,E->tracer.ty1,
			   E->tracer.tx, E->tracer.tz, E->tracer.ty ,
			   E->V[1], E->V[2], E->V[3],   
			   E->V[1], E->V[2], E->V[3],
			   1,E->tracer.NUM_TRACERS,
			   E->tracer.tracer_elt[E->mesh.levmax],NULL);

  /* Local strain information */

  /* X dirn */

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    dX1[m] = E->tracer.tx1[m] + E->tracer.dX11[m];
    dX2[m] = E->tracer.tz1[m] + E->tracer.dX12[m];
    if(3==dims)
      dX3[m] = E->tracer.ty1[m] + E->tracer.dX13[m];

    el_list[m] = E->tracer.tracer_elt[E->mesh.levmax][m];
    general_tracer_within_boundaries(E,dX1,dX2,dX3,m);
  }
  printf("tracer_adv: 3 doing general_tracer_advection \n");

  general_tracer_advection(E,E->advection.timestep,
			    dX1,dX2,dX3,
			    dX1,dX2,dX3,
			    E->V[1], E->V[2], E->V[3],   
			    E->V[1], E->V[2], E->V[3],
			    1,E->tracer.NUM_TRACERS,
			    el_list,NULL);

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(el_list[m] == -1) {
      fprintf(stderr,"Tracer %d's X shadow is a little lost %g,%g\n",m,E->tracer.dX11[m],E->tracer.dX12[m]); 
      E->tracer.dX11[m] *= 0.9;
      E->tracer.dX12[m] *= 0.9;
      if(3==dims)
	E->tracer.dX13[m]  *= 0.9; 
    }
    else { 
      E->tracer.dX11[m] = dX1[m] - E->tracer.tx[m];
      E->tracer.dX12[m] = dX2[m] - E->tracer.tz[m];
      if(3==dims)
	E->tracer.dX13[m] = dX3[m] - E->tracer.ty[m];
    }
  }
  /* Z dirn */


   for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    dX1[m] = E->tracer.tx1[m] + E->tracer.dX21[m];
    dX2[m] = E->tracer.tz1[m] + E->tracer.dX22[m];
    if(3==dims)
      dX3[m] = E->tracer.ty1[m] + E->tracer.dX23[m];
    el_list[m] = E->tracer.tracer_elt[E->mesh.levmax][m];
    general_tracer_within_boundaries(E,dX1,dX2,dX3,m);
  }
   printf("tracer_adv: 4 doing general_tracer_advection \n");

   general_tracer_advection(E,E->advection.timestep,
			     dX1,dX2,dX3,
			     dX1,dX2,dX3,
			     E->V[1], E->V[2], E->V[3],   
			     E->V[1], E->V[2], E->V[3],
			     1,E->tracer.NUM_TRACERS,
			     el_list,NULL);
 
 
   for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
     if(el_list[m] == -1) {
       /* fprintf(stderr,"Tracer %d's Z shadow is a little lost\n",m); */
       E->tracer.dX21[m] *= 0.9;
       E->tracer.dX22[m] *= 0.9;
       if(3==dims)
	 E->tracer.dX23[m] *= 0.9; 
     }
     else {
       E->tracer.dX21[m] = dX1[m] - E->tracer.tx[m];
       E->tracer.dX22[m] = dX2[m] - E->tracer.tz[m];
       if(3==dims)
	 E->tracer.dX23[m] = dX3[m] - E->tracer.ty[m];
     }
  }
 
   /* A difficulty occurs when periodic bc's are in force.
      The tracker particles can be wrapped around away from 
      their host particles (or vice versa). There's not really 
      an elegant solution to this, given that the local velocity
      gradient tracker is pretty awful */

   if(E->mesh.periodic_x) {
     for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
       if(E->tracer.dX21[m] > (E->x[1][E->mesh.nno]-E->x[1][1]) * 0.9)
	 E->tracer.dX21[m] -= (E->x[1][E->mesh.nno]-E->x[1][1]);
       else if(E->tracer.dX21[m] < (E->x[1][E->mesh.nno]-E->x[1][1]) * -0.9)
	 E->tracer.dX21[m] += (E->x[1][E->mesh.nno]-E->x[1][1]);

      if(E->tracer.dX11[m] > (E->x[1][E->mesh.nno]-E->x[1][1]) * 0.9)
	 E->tracer.dX11[m] -= (E->x[1][E->mesh.nno]-E->x[1][1]);
       else if(E->tracer.dX11[m] < (E->x[1][E->mesh.nno]-E->x[1][1]) * -0.9)
	 E->tracer.dX11[m] += (E->x[1][E->mesh.nno]-E->x[1][1]);
     }

   }
printf("tracer_adv: freeing memory \n");

free((void *) dX1); 
free((void *) dX2); 
free((void *) dX3); 
free((void *) el_list);

   if(E->control.verbose)
    fprintf(stderr,"Particle advection ... 2\n");

printf("tracer_adv: Done!!! \n");

  return;
}

/* =================================================
   tracer_advection: Midpoint method (2nd order RK)
   =================================================  */

void tracer_corrector(
		      struct All_variables *E,
		      standard_precision *U1,
		      standard_precision *U2,
		      standard_precision *U3,
		      standard_precision timestep
		      )
{
  /* Use the general tracer advection routine to

     1) Apply a corrector to the previous particle locations
     if this is not the first timestep 

     2) Shoot the current particle positions forward using the
     latest velocity solution (predictor)
  */

  int n;

  void general_tracer_advection();
  void standard_tracer_advection();

  /* Corrector step */
  if(E->advection.timesteps > 1)  { /* This is called after the first timestep value has been calculated */
    
    standard_tracer_advection(E,timestep,
			      E->tracer.tx1,E->tracer.tz1,E->tracer.ty1,
			      E->tracer.tx, E->tracer.tz, E->tracer.ty ,
			      E->V1[1], E->V1[2], E->V1[3],
			      U1,U2,U3,
			      1,E->tracer.NUM_TRACERS,
			      E->tracer.tracer_elt[E->mesh.levmax],
			      NULL
			      );
  }
  
  return;
}

/* =================================================
   Generalized tracer advection routine ... to 
   take existing particle positions and update them
   based on velocities supplied for beginning and
   end of the timestep. This routine can be 
   used for both the predictor and the corrector
   phases of the update if the appropriate velocity
   fields are supplied.
   =================================================  */

void general_tracer_advection(
			      struct All_variables *E,
			      standard_precision timestep,
			      standard_precision *X1,
			      standard_precision *Z1,
			      standard_precision *Y1,   /* previous location */
			      standard_precision *X,
			      standard_precision *Z,
			      standard_precision *Y,      /* new location  */
			      standard_precision *U1,
			      standard_precision *W1,
			      standard_precision *V1,   /* previous velocity */
			      standard_precision *U,
			      standard_precision *W,
			      standard_precision *V,      /* new velocity */
			      int Nstart,
			      int NQ,
			      int *in_element,
			      int *lagrangian
			      )
{
  int n;
  int i,j,k,m;
  int n_x,n_y,n_z;
  int node1;
  int iteration;
  int level;
  int tr,el,el1;
     
  standard_precision eta1,eta2,eta3;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision dOmega;
  standard_precision kx1,kx2,kx3,kx4;
  standard_precision kz1,kz2,kz3,kz4;
  standard_precision ky1,ky2,ky3,ky4;
  standard_precision vx,vz,vy,vzz;
  standard_precision x0,z0,y0;

  standard_precision half_delta_t;
  standard_precision CPU_time(),time;

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int ends = enodes[dims];

  struct IEN *IEN = E->ien;

  half_delta_t = timestep * 0.5;
  time=CPU_time();

  /* Loop over each tracer in order */

  for(n=Nstart;n<Nstart+NQ;n++) {
   if((lagrangian != NULL) && (!lagrangian[n]))
      continue;
  
    /* Initial location */

    x0 = X1[n] ;
    z0 = Z1[n] ;
    if(3==dims)
      y0 = Y1[n] ;
    else
      y0 = 0.0 ;

    /* Velocity at initial location */

    /* printf("In tracer advection calling tracers_element 1st \n"); */

    el1 = in_element[n];
    el1 = general_tracers_element(E,el1,x0,z0,y0,&eta1,&eta2,&eta3,E->mesh.levmax);

    if(0 && el1 != in_element[n])
      fprintf(stderr,"Tracer %d is different element from expectations ! (%g,%g -> %d v %d) \n",n,
	     x0,z0,el1,in_element[n] );

    if(el1 == -1) {  /* Asked to advect something not in the mesh ! */
      /* fprintf(stderr,"Warning, tracer %d at %g,%g appears to be lost\n",n,x0,z0); */
      in_element[n] = -1;
      continue;
    }
    else
      in_element[n] = el1;

    v_shape_fn(E,el1,lN,eta1,eta2,eta3,E->mesh.levmax);
 
    vx=vz=vy=0.0;
    for(m=1;m<=ends;m++) {
      node1 = E->ien[el1].node[m];

      vx += lN[m] * U1[node1];
      vz += lN[m] * W1[node1];
      if(dims==3)
	vy += lN[m] * V1[node1];
    }

    /* Determine midpoint values */

    kx1 = timestep * vx;
    kz1 = timestep * vz;
    X[n] = x0 + kx1 * 0.5;
    Z[n] = z0 + kz1 * 0.5;
    if(dims==3) {
      ky1 = timestep * vy;
      Y[n] = y0 + ky1 * 0.5;
    }

    if(E->mesh.periodic_x) {
      if(X[n]  > E->x[1][E->mesh.nno]) {
	X[n]  -= (E->x[1][E->mesh.nno]-E->x[1][1]);
      }
      if(X[n]  < E->x[1][1]) {
	X[n]  += (E->x[1][E->mesh.nno]-E->x[1][1]);
      }
    }
 
  /*RAA: 18/10/01, need to add these ~10 lines for periodic_y */
    if(3==dims && E->mesh.periodic_y) {  
      if(Y[n]  > E->x[3][E->mesh.nno]) {
	Y[n]  -= (E->x[3][E->mesh.nno]-E->x[3][1]);
      }
      if(Y[n]  < E->x[3][1]) {
	Y[n]  += (E->x[3][E->mesh.nno]-E->x[3][1]);
      }
    }
 
    /* Velocity at x + k1 /2 */

    /* printf("Hello all. In trace_adv still, calling 2nd time \n"); */

    el1 = general_tracers_element(E,el1,X[n],Z[n],((E->mesh.nsd==3) ? Y[n] : 0.0),
				  &eta1,&eta2,&eta3,E->mesh.levmax);

    if(el1 == -1)  {
      el1 = in_element[n];
      /* has moved outside grid ... but maybe not accurate as this
	 is the mid point guess of the real step, 
	 so use orginal element as best guess 
	 until correction is made */
      fprintf(stderr,"Tracer %d is being reset to old element %d\n",n,el1);

      get_element_coords(E,el1,n,X,Z,Y,&eta1,&eta2,&eta3,E->mesh.levmax); /* Should be OK */
    }
    else
      in_element[n] = el1;

    v_shape_fn(E,el1,lN,eta1,eta2,eta3,E->mesh.levmax);

    vx=vz=vy=0.0;
    for(m=1;m<=ends;m++) {
      node1 = E->ien[el1].node[m];   

      vx += lN[m] * 0.5 * (U1[node1] + U[node1]);
      vz += lN[m] * 0.5 * (W1[node1] + W[node1]);
      if(3==dims)
	vy += lN[m] * 0.5 * (V1[node1] + V[node1]);
    }

    /* Final location */
    
    kx2 = timestep * vx;
    kz2 = timestep * vz;
    X[n] = x0 + kx2;
    Z[n] = z0 + kz2;
    if(3==dims) {
      ky2 = timestep * vy;
      Y[n] = y0 + ky2;
    }
 
    /* if periodic, then tracers can scoot around out of the
       box and need to be scooped up again */

    if(E->mesh.periodic_x) {
      if(X[n]  > E->x[1][E->mesh.nno]) {
	X[n]  -= (E->x[1][E->mesh.nno]-E->x[1][1]);
      }
      if(X[n]  < E->x[1][1]) {
	X[n]  += (E->x[1][E->mesh.nno]-E->x[1][1]);
      }
    }
  /*RAA: 18/10/01, need to add these ~10 lines for periodic_y */
    if(3==dims && E->mesh.periodic_y) {  
      if(Y[n]  > E->x[3][E->mesh.nno]) {
	Y[n]  -= (E->x[3][E->mesh.nno]-E->x[3][1]);
      }
      if(Y[n]  < E->x[3][1]) {
	Y[n]  += (E->x[3][E->mesh.nno]-E->x[3][1]);
      }
    }
 
  }
  return;
}


/* =================================================
   Generalized tracer advection routine ... to 
   take existing particle positions and update them
   based on velocities supplied for beginning and
   end of the timestep. This routine can be 
   used for both the predictor and the corrector
   phases of the update if the appropriate velocity
   fields are supplied.

   This is the vector version - move tracers 
   incrementally all at the same time.

   MOD - need to supply and return shape function
   information as well as element distributions.

   =================================================  */

void standard_tracer_advection(
			       struct All_variables *E,
			       standard_precision timestep,
			       standard_precision *X1,
			       standard_precision *Z1,
			       standard_precision *Y1,   /* previous location */
			       standard_precision *X,
			       standard_precision *Z,
			       standard_precision *Y,      /* new location  */
			       standard_precision *U1,
			       standard_precision *W1,
			       standard_precision *V1,   /* previous velocity */
			       standard_precision *U,
			       standard_precision *W,
			       standard_precision *V,      /* new velocity */
			       int N1,
			       int N2,
			       int *el_list,
			       int *lagrangian
			       )
{
  int n;
  int el;
  int i,j,k,m;

  struct TRACER_ELT_WEIGHT *lN;

  standard_precision *VX,*VZ,*VY;
  standard_precision *KX,*KZ,*KY;
  standard_precision *XX0,*ZZ0,*YY0;
  standard_precision *eta1;
  standard_precision *eta2;
  standard_precision *eta3;
  standard_precision perx0,perx1;
  standard_precision pery0,pery1; /*RAA: 18/10/01, add this line*/
  standard_precision half_delta_t;
  standard_precision CPU_time(),time;

  int *old_el_list;
  struct IEN *IEN;

  int all_tracers_elts_and_sfns();

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int ends = enodes[dims];


  if(N2-N1 == 0 || N2+1 == 0)   /* Nothing wrong with this situation, but nothing
		      to do if this is the case */
    return;

  lN = (struct TRACER_ELT_WEIGHT *) Malloc0((N2+1) * sizeof(struct TRACER_ELT_WEIGHT));
  VX = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  VZ = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  VY = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  KX = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  KZ = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  KY = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  XX0 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  ZZ0 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  YY0 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  eta1 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  eta2 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  eta3 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));

  old_el_list = (int *) Malloc0((N2+1) * sizeof(int));

  /* Shorthand for some tracers to help vectorizer */

  IEN = E->ien;

  if(E->mesh.periodic_x) {   /* How should this be generalized ? */
    perx0 = E->x[1][1];
    perx1 = E->x[1][E->mesh.nno];
  }
  /*RAA: 18/10/01, add this for periodic_y*/
  if(3==dims && E->mesh.periodic_y) {   /* How should this be generalized ? */
    pery0 = E->x[3][1];
    pery1 = E->x[3][E->mesh.nno];
  }


  time=CPU_time();

  /* use a copy of the element list in case we need to backtrack */

/* #pragma loop novrec old_el_list,el_list */
  for(m=N1;m<=N2;m++) {
    old_el_list[m] = el_list[m];
    if(el_list[m] == -1)
      fprintf(stderr,"tracer %d/%d already on its own (%d)  %g,%g\n",m,N2,old_el_list[m],X[m],Z[m]);
  }

  /* Get shape functions in order to obtain velocity field.
     NOTE - this information is probably already available
     and could be passed into this routine. */

  /* tr_local_coords(E,el_list,lN,X1,Z1,Y1,eta1,eta2,eta3,
     N1,N2,E->mesh.levmax); */

  k = all_tracers_elts_and_sfns(E,el_list,lN,X1,Z1,Y1,
				eta1,eta2,eta3,N1,N2,E->mesh.levmax);

  for(m=N1;m<=N2;m++) {
    if(fabs(eta1[m]) > 1.0 || fabs(eta2[m]) > 1.0)
      fprintf(stderr,"tracer %d/%d not found (%d)  %g,%g\n",m,N2,old_el_list[m],X[m],Z[m]); 
  }

  /* Determine velocity ... */
  if(2==dims) {
/* #pragma loop novrec VX,VZ */
    for(m=N1;m<=N2;m++) {
      VX[m] = 0.0;
      VZ[m] = 0.0;
    }

    for(k=1;k<=ends;k++) {
/* #pragma loop novrec VX,VZ,U1,W1,lN,IEN,el_list */
	  for(m=N1;m<=N2;m++) {
	    el = el_list[m];
	    VX[m] += U1[IEN[el].node[k]] * lN[m].node[k];
	    VZ[m] += W1[IEN[el].node[k]] * lN[m].node[k];
	  }
	}
  }
  else {
/* #pragma loop novrec VX,VZ,VY */
    for(m=N1;m<=N2;m++) {
      VX[m] = 0.0;
      VZ[m] = 0.0;
      VY[m] = 0.0;
    }
    
    for(k=1;k<=ends;k++) {
/* #pragma loop novrec VX,VZ,VY,U1,V1,W1,lN,IEN,el_list */
	  for(m=N1;m<=N2;m++) {
	    el = el_list[m];
	    VX[m] += U1[IEN[el].node[k]] * lN[m].node[k];
	    VZ[m] += W1[IEN[el].node[k]] * lN[m].node[k];
	    VY[m] += V1[IEN[el].node[k]] * lN[m].node[k];
	  }
	}
  }

  /* Determine midpoint positions */

  if(2==dims)  
/* #pragma loop novrec XX0,ZZ0,X1,Z1,VX,VZ */
    for(m=N1;m<=N2;m++) {
      XX0[m] = X1[m] + 0.5 * timestep * VX[m];
      ZZ0[m] = Z1[m] + 0.5 * timestep * VZ[m];
    }
  else
/* #pragma loop novrec XX0,ZZ0,YY0,X1,Z1,Y1,VX,VZ,VY */
    for(m=N1;m<=N2;m++) {
      XX0[m] = X1[m] + 0.5 * timestep * VX[m];
      ZZ0[m] = Z1[m] + 0.5 * timestep * VZ[m];
      YY0[m] = Y1[m] + 0.5 * timestep * VY[m];
    }
  
  /* Catch wanderers if periodic bc's */

  if(E->mesh.periodic_x) {   /* general case - this would call a function
				to wraparound tracers - periodic bc's are only
				meaningful when some high-symmetry geometry exists */
    for(m=N1;m<=N2;m++) {
/* #pragma loop novrec XX0 */
      XX0[m] = fmod(perx1+XX0[m],perx1);  /* should be (XX0[m] - perx0) */
    }
  }

/*RAA: 18/10/01, add this for periodic_y*/
  if(3==dims && E->mesh.periodic_y) {
    for(m=N1;m<=N2;m++) {
/* #pragma loop novrec YY0 */
      YY0[m] = fmod(pery1+YY0[m],pery1);  /* should be (YY0[m] - pery0), RAA: verify*/
    }
  }

  k = all_tracers_elts_and_sfns(E,el_list,lN,XX0,ZZ0,YY0,
			    eta1,eta2,eta3,N1,N2,E->mesh.levmax);

  /* Some may have wandered out of the mesh with this method ... */

  if(k)
    for(m=N1;m<=N2;m++) {  /* to vectorize we would need a list-based version of 
			      the element coords and v_shape_fns */
      if(el_list[m] == -1) {
	fprintf(stderr,"tracer %d wandered off on its own (%d) %g,%g <- %g,%g (%g,%g)\n",
		m,old_el_list[m],XX0[m],ZZ0[m],X[m],Z[m],VX[m],VZ[m]); 
	/* Restore it to old element */
	el_list[m] = old_el_list[m];
	get_element_coords(E,el_list[m],m,XX0,ZZ0,YY0,
			   &(eta1[m]),&(eta2[m]),&(eta3[m]),E->mesh.levmax); /* Should be OK */
	v_shape_fn(E,el_list[m],&(lN[m]),&(eta1[m]),&(eta2[m]),&(eta3[m]),E->mesh.levmax);
      }
    }

 /* Velocity at x + k1 /2 */

  if(2==dims) {
/* #pragma loop novrec VX,VZ */
    for(m=N1;m<=N2;m++) {
      VX[m] = 0.0;
      VZ[m] = 0.0;
    }

    for(k=1;k<=ends;k++) {
/* #pragma loop novrec VX,VZ,U1,W1,U,W,lN,IEN */
      for(m=N1;m<=N2;m++) {
	el = el_list[m];
	VX[m] += 0.5 * (U1[IEN[el].node[k]]+U[IEN[el].node[k]]) * lN[m].node[k];
	VZ[m] += 0.5 * (W1[IEN[el].node[k]]+W[IEN[el].node[k]]) * lN[m].node[k];
      }
    }
  }
  else {
/* #pragma loop novrec VX,VZ,VY */
    for(m=N1;m<=N2;m++) {
      VX[m] = 0.0;
      VZ[m] = 0.0;
      VY[m] = 0.0;
    }
    
    for(k=1;k<=ends;k++) {
/* #pragma loop novrec VX,VZ,VY,U1,W1,V1,U,W,V,lN,IEN */
      for(m=N1;m<=N2;m++) {
	el = el_list[m];
	VX[m] += 0.5 * (U1[IEN[el].node[k]]+U[IEN[el].node[k]]) * lN[m].node[k];
	VZ[m] += 0.5 * (W1[IEN[el].node[k]]+W[IEN[el].node[k]]) * lN[m].node[k];
	VY[m] += 0.5 * (V1[IEN[el].node[k]]+V[IEN[el].node[k]]) * lN[m].node[k];
      }
    }
  }

  /* And then use this velocity to update the particle over
     the whole timestep */

  if(lagrangian == NULL) {
    if(2==dims)
/* #pragma loop novrec X,Z,X1,Z1,VX,VZ */
      for(m=N1;m<=N2;m++) {
	X[m] = X1[m] + timestep * VX[m];
	Z[m] = Z1[m] + timestep * VZ[m];
    }
    else
/* #pragma loop novrec X,Z,Y,X1,Z1,Y1,VX,VZ,VY */
      for(m=N1;m<=N2;m++) {
	X[m] = X1[m] + timestep * VX[m];
	Z[m] = Z1[m] + timestep * VZ[m];
	Y[m] = Y1[m] + timestep * VY[m];
      }
  }
  else {  /* The possibility that some tracers are
	     not to be updated */
     if(2==dims)  
/* #pragma loop novrec X,Z,X1,Z1,VX,VZ */
       for(m=N1;m<=N2;m++) {
	 if(lagrangian[m]) {
	   X[m] = X1[m] + timestep * VX[m];
	   Z[m] = Z1[m] + timestep * VZ[m];
	 }
       }
    else
/* #pragma loop novrec X,Z,Y,X1,Z1,Y1,VX,VZ,VY */
      for(m=N1;m<=N2;m++) {
	 if(lagrangian[m]) {
	   X[m] = X1[m] + timestep * VX[m];
	   Z[m] = Z1[m] + timestep * VZ[m];
	   Y[m] = Y1[m] + timestep * VY[m];

	 }
      }
  }

/* Again, catch wanderers if periodic bc's */
    
  if(E->mesh.periodic_x) {   
/* #pragma loop novrec X */
    for(m=N1;m<=N2;m++) {
      X[m] = fmod(perx1+X[m],perx1);  /* should be (XX0[m] - perx0) */
    }
  }
  /*RAA: 18/10/01, add this for periodic_y*/
  if(3==dims && E->mesh.periodic_y) {
/* #pragma loop novrec Y */
    for(m=N1;m<=N2;m++) {
      Y[m] = fmod(pery1+Y[m],pery1);  /* should be (YY0[m] - perY0), RAA: verify this*/
    }
  }

  free((void *) lN);
  free((void *) VX);
  free((void *) VZ);
  free((void *) VY);
  free((void *) KX);
  free((void *) KZ);
  free((void *) KY);
  free((void *) XX0);
  free((void *) ZZ0);
  free((void *) YY0);
  free((void *) eta1);
  free((void *) eta2);
  free((void *) eta3);
  free((void *) old_el_list);

  return;
}



/* Rotation of the director for Orthotropic materials */

void rotate_director (  
		      struct All_variables *E
		      )
{
  
  higher_precision *node_R;

  standard_precision V11,V22,V33,V12,V21,V31,V13,V23,V32;

  standard_precision vvx,vvz;

  standard_precision eta1,eta2,eta3,dOmega;
  standard_precision lNx[4][ELNMAX+1],lN[ELNMAX+1];
  standard_precision U1,U2,U3;
  standard_precision n1dot,n2dot,n3dot;
  standard_precision x1,z1,y1;
  standard_precision mag;

  const int dofs = E->mesh.dof ;
  const int dims = E->mesh.nsd ;
  const int ends = enodes[dims] ;
  const int level = E->mesh.levmax;

  int node,m,i,el;

  struct IEN *IEN = E->ien;
  
  /* 2D Classical */

  if(2==dims && 2==dofs) {

    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      if((E->control.CHEM_TRANS && E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio == 1.0) ||
	 E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio == 1.0) 
	continue;


      /* 1. get deformation rate tensor */

      V11 = V12 = V21 = V22 = 0.0;
      vvx=vvz=0.0;

      el = E->tracer.tracer_elt[E->mesh.levmax][m];
      x1 = 0.5 * (E->tracer.tx[m] + E->tracer.tx1[m]);
      z1 = 0.5 * (E->tracer.tz[m] + E->tracer.tz1[m]);

      /* printf("3rd \n"); */

      general_tracers_element(E,el,x1,z1,NULL,&eta1,&eta2,NULL,E->mesh.levmax); /* */

      /*
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
      /* */
 
      get_global_v_x_shape_fn(E,E->tracer.tracer_elt[level][m],lNx,&dOmega,eta1,eta2,0.0,level);

      for(i=1;i<=ends;i++) {
	node = IEN[E->tracer.tracer_elt[level][m]].node[i];

	if(!(E->control.CURRENT_SKEWED_Vs && E->NODE[level][node] & SKEWBC)) {
	  U1 = E->VV[level][1][node];
	  U2 = E->VV[level][2][node]; 
	}
	else {
	  node_R = E->curvilinear.NODE_R[level][node];
	  U1 = node_R[0*dims+0] * E->VV[level][1][node] + node_R[0*dims+1] * E->VV[level][2][node];
	  U2 = node_R[1*dims+0] * E->VV[level][1][node] + node_R[1*dims+1] * E->VV[level][2][node];
	} 

	vvx+=U1*E->tracer.sfn_values[level][m].node[i];
	vvz+=U2*E->tracer.sfn_values[level][m].node[i];

	V11 += U1 * lNx[1][i];
	V21 += U2 * lNx[1][i];
	V12 += U1 * lNx[2][i];
	V22 += U2 * lNx[2][i];
      }

      /* n1dot = V11 *  E->tracer.n1[m] * (E->tracer.n1[m] *  E->tracer.n1[m] -  1.0) 
	- V21 * E->tracer.n2[m] + (V21 + V12) * E->tracer.n1[m] *  E->tracer.n1[m] * E->tracer.n2[m]
	+ V22 * E->tracer.n1[m] *  E->tracer.n2[m] * E->tracer.n2[m];
     
      n2dot = V22 *  E->tracer.n2[m] * (E->tracer.n2[m] *  E->tracer.n2[m] -  1.0) 
	- V12 * E->tracer.n1[m] + (V21 + V12) * E->tracer.n1[m] *  E->tracer.n2[m] * E->tracer.n2[m]
	+ V11 * E->tracer.n1[m] *  E->tracer.n1[m] * E->tracer.n2[m]; */
     
      n1dot = -V11 *  E->tracer.n1[m] - V21 * E->tracer.n2[m];
      n2dot = -V12 *  E->tracer.n1[m] - V22 * E->tracer.n2[m];

 
      E->tracer.n1[m] += n1dot * E->advection.timestep;
      E->tracer.n2[m] += n2dot * E->advection.timestep;
    
      mag = sqrt(E->tracer.n1[m]*E->tracer.n1[m]+E->tracer.n2[m]*E->tracer.n2[m]);

      if(mag != 0.0)
	mag = 1.0 / mag;
      else
	mag = 1.0;

      E->tracer.n1[m] *= mag;
      E->tracer.n2[m] *= mag;


      /*  if( E->tracer.n1[m] * E->tracer.n1[m] +  E->tracer.n2[m] * E->tracer.n2[m] > 1.05)
	fprintf(stderr,"Tracer %d, n = (%g,%g), |n| = %g\n", m,
		E->tracer.n1[m],E->tracer.n2[m],
		sqrt(E->tracer.n1[m] * E->tracer.n1[m] +  E->tracer.n2[m] * E->tracer.n2[m])); */

    }
  }


    return;

}






