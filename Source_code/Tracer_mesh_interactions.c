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
#include <stdlib.h>

#if (defined __sunos__) || defined(__GNUC__)
#include <string.h>
#else
#include <strings.h> 
#endif

/*
#if (! defined __GNUC__)
#include <rpc/xdr.h> 
#endif
*/

#include "element_definitions.h"
#include "global_defs.h"/******************************************************************/



void tracer_within_boundaries(
     struct All_variables *E,
     int tr
)

{
  /* if periodic, then tracers can scoot around out of the
     box and need to be scooped up again */

  if(E->mesh.periodic_x) {
    if(E->tracer.tx[tr]  > E->x[1][E->mesh.nno])
      E->tracer.tx[tr]  -= (E->x[1][E->mesh.nno]-E->x[1][1]);

    if(E->tracer.tx[tr]  < E->x[1][1])
      E->tracer.tx[tr]  += (E->x[1][E->mesh.nno]-E->x[1][1]);
  }

   if(E->mesh.periodic_x) {
    if(E->tracer.tx1[tr]  > E->x[1][E->mesh.nno])
      E->tracer.tx1[tr]  -= (E->x[1][E->mesh.nno]-E->x[1][1]);

    if(E->tracer.tx1[tr]  < E->x[1][1])
      E->tracer.tx1[tr]  += (E->x[1][E->mesh.nno]-E->x[1][1]);
  }
  return;
}

void general_tracer_within_boundaries(
     struct All_variables *E,
     standard_precision *tx,
     standard_precision *tz,
     standard_precision *ty,
     int tr
)

{
  /* if periodic, then tracers can scoot around out of the
     box and need to be scooped up again */

  if(E->mesh.periodic_x) {
    if(tx[tr]  > E->x[1][E->mesh.nno])
      tx[tr]  -= (E->x[1][E->mesh.nno]-E->x[1][1]);

    if(tx[tr]  < E->x[1][1])
      tx[tr]  += (E->x[1][E->mesh.nno]-E->x[1][1]);
  }

  /*RAA: 10/04/02, this function is now used in Linux version, 
   * so add per_y stuff*/
  if(E->mesh.periodic_y) {
    if(ty[tr] > E->x[3][E->mesh.nno])
       ty[tr] -= (E->x[3][E->mesh.nno]-E->x[3][1]);

    if(ty[tr] < E->x[3][1])
       ty[tr]  += (E->x[3][E->mesh.nno]-E->x[3][1]);
  }

  return;
}




void gs_tracers_to_nodes (		 
 struct All_variables *E,
 standard_precision *field,
 standard_precision *BCval,
 int *BCinfo,
 int BCmask,
 standard_precision *property,
 int level,
 int guess
 )
{
  int p,q,m,n,i,j,k,l,node,el,node1;
  int M,N;

  standard_precision *rk0,*rk1;
  standard_precision *xk0,*pk0;
  standard_precision *Ap,*zk0,*zk1;
  standard_precision *field1;
  standard_precision rdot1,rdot2,pdotAp,alpha,beta;
  standard_precision initial_res,res,res_mag;
  standard_precision time,CPU_time();

  standard_precision extra_acc;

  rk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  rk1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  xk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  pk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  zk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  zk1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  Ap = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  field1 = (standard_precision *)Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));

  time=CPU_time();
  M=1;
  N=E->tracer.NUM_TRACERS;
  extra_acc=1.0;

  /* Least squares improvement of nodal value fit. (Using Gauss-Seidel) */
 
  for(n=1;n<=E->mesh.NNO[level];n++) 
    field[n] = field1[n]=0.0;

  for(el = 1; el <= E->mesh.NEL[level]; el++) {
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      if(E->tracer.property_group[m] < 0)
	continue;

      if(m > E->tracer.NUM_TRACERS)
	fprintf(stderr,"Tracer %d/%d of element %d is global %d/%d, %g / %g (%g,%g)\n",
		i,E->tracer.tr_in_element_number[level][el],el,m,E->tracer.NUM_TRACERS,
		property[m],E->tracer.tracer_weight_fn[E->mesh.levmax][m],
		E->tracer.T[m],E->tracer.weight[m]);
     
      for(n=1;n<=enodes[E->mesh.nsd];n++) {

    field[E->IEN[level][el].node[n]] += property[m] * E->tracer.sfn_values[level][m].node[n] *
	  E->tracer.tracer_weight_fn[E->mesh.levmax][m];
    field1[E->IEN[level][el].node[n]] += E->tracer.sfn_values[level][m].node[n] *
	  E->tracer.tracer_weight_fn[E->mesh.levmax][m];
      }
    }
  }

  for(n=1;n<=E->mesh.NNO[level];n++) 
    if(field1[n] != 0.0)
      field[n] /= field1[n];
  
  /* If boundary conditions are supplied then they must be
     applied to the input field first */

  if(BCinfo != NULL && BCval != NULL)
    for(i=1;i<=E->mesh.NNO[level];i++)
      if(BCinfo[i] & BCmask)
	field[i] = BCval[i];



#if 0

 /* Initialize */

  for(i=1;i<=E->mesh.NNO[level];i++) {
    rk0[i] = 0.0;
    xk0[i] = field[i];
    Ap[i] = 0.0;
  }

  for(m=M;m<=N;m++) {
    el = E->tracer.tracer_elt[level][m];
    for(j=1;j<=enodes[E->mesh.nsd];j++) {
      node = E->IEN[level][el].node[j];
      if((BCinfo != NULL && BCval != NULL) && (BCinfo[node] & BCmask)) 
	xk0[i] = 0.0;
      else 
	/* rk0 */ xk0[node] += E->tracer.sfn_values[level][m].node[j] * property[m];
    }
  }



  /* Typical magnitude of RHS is stored in Ap */

  rdot1=0.0;
  for(i=1;i<=E->mesh.NNO[level];i++) 
    rdot1 += rk0[i] * rk0[i];

  res_mag = initial_res = 1.0e-16 + sqrt(rdot1)/E->mesh.NNO[level];

   /* GS */

  k=0;
    do {
      k++;
       
  /* get A*pk0 term */ 
    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;

    for(m=M;m<=N;m++) {
      el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	if((BCinfo != NULL && BCval != NULL) &&  (BCinfo[node] & BCmask))
	  Ap[node] = 0.0;
	else
	  for(p=1;p<=enodes[E->mesh.nsd];p++) {
	    node1 = E->IEN[level][el].node[p]; 
	    Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	      E->tracer.sfn_values[level][m].node[p] * xk0[node1];
	  }
      }
    }

   /* In 2D do a four-colour ordering - four passes on nodes
       which do not have any interaction with each other. This eliminates
       searching for tracers in elements as we just compute the 
       residual 4 times over (yes, this is quite a bit quicker !) */

    /* PT 1 */
      
    for(i=1;i<=E->mesh.NOX[level];i+=2)
      for(j=1;j<=E->mesh.NOZ[level];j+=2) {
	node = j + (i-1) * E->mesh.NOZ[level];
	xk0[node] += (rk0[node] - Ap[node]) * E->tracer.ls_diag[level][node];
      }

    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;

    for(m=M;m<=N;m++) {
      el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	for(p=1;p<=enodes[E->mesh.nsd];p++) {
	    node1 = E->IEN[level][el].node[p]; 
	    Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	      E->tracer.sfn_values[level][m].node[p] * xk0[node1];
	  }
      }
    }

    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++)
	if(BCinfo[i] & BCmask)
	  Ap[i] = 0.0;


    
    /* PT 2 */
      
    for(i=2;i<=E->mesh.NOX[level];i+=2)
      for(j=1;j<=E->mesh.NOZ[level];j+=2) {
	node = j + (i-1) * E->mesh.NOZ[level];
	xk0[node] += (rk0[node] - Ap[node]) * E->tracer.ls_diag[level][node];
      }
    
    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;
    
    for(m=M;m<=N;m++) {
      el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	for(p=1;p<=enodes[E->mesh.nsd];p++) {
	  node1 = E->IEN[level][el].node[p]; 
	  Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	    E->tracer.sfn_values[level][m].node[p] * xk0[node1];
	  }
      }
    }
    
    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++)
	if(BCinfo[i] & BCmask)
	  Ap[i] = 0.0;
    
    
    /* PT 3 */
      
     for(i=1;i<=E->mesh.NOX[level];i+=2)
      for(j=2;j<=E->mesh.NOZ[level];j+=2) {
	node = j + (i-1) * E->mesh.NOZ[level];
	xk0[node] += (rk0[node] - Ap[node]) * E->tracer.ls_diag[level][node];
      }
    
    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;
    
    for(m=M;m<=N;m++) {
      el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	for(p=1;p<=enodes[E->mesh.nsd];p++) {
	  node1 = E->IEN[level][el].node[p]; 
	  Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	    E->tracer.sfn_values[level][m].node[p] * xk0[node1];
	  }
      }
    }
    
    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++)
	if(BCinfo[i] & BCmask)
	  Ap[i] = 0.0;
    
  /* PT 4 */
      
    for(i=2;i<=E->mesh.NOX[level];i+=2)
      for(j=2;j<=E->mesh.NOZ[level];j+=2) {
	node = j + (i-1) * E->mesh.NOZ[level];
	xk0[node] += (rk0[node] - Ap[node]) * E->tracer.ls_diag[level][node];
      }
    
    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;
    
    for(m=M;m<=N;m++) {
      el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	for(p=1;p<=enodes[E->mesh.nsd];p++) {
	  node1 = E->IEN[level][el].node[p]; 
	  Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	    E->tracer.sfn_values[level][m].node[p] * xk0[node1];
	  }
      }
    }
    
    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++)
	if(BCinfo[i] & BCmask)
	  Ap[i] = 0.0;
    

    /* Compute residual */
  
    rdot1=0.0;
    for(i=1;i<=E->mesh.NNO[level];i++) {
      rdot1 += (rk0[i]-Ap[i]) * (rk0[i]-Ap[i]);
    }

    res = sqrt(rdot1)/E->mesh.NNO[level];

     fprintf(stderr,"%d: LS residual [%d] = %g (%g) --- %gs\n",
	  level,k,res,res/initial_res,CPU_time()-time );

  }   while (k <= 500 && res >  extra_acc * E->control.accuracy * res_mag)  ;


	       
    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(BCinfo[i] & BCmask)
	  field[i] = BCval[i];
	else {
	  if(i<=40) {
	    fprintf(stderr,"Nodal value %g v %g\n",field[i],xk0[i]);
	      }
	  field[i] = xk0[i];
	}
      }

#endif

  free((void *) rk0);
  free((void *) rk1);
  free((void *) xk0);
  free((void *) pk0);
  free((void *) zk0);
  free((void *) zk1);
  free((void *) Ap);
  free((void *) field1);	       
  return;
}

void tracers_to_nodes(
  struct All_variables *E,
  standard_precision *field,
  standard_precision *BCval,
  int *BCinfo,
  int BCmask,
  standard_precision *property,
  int level,
  int guess
)
{
  int p,q,m,i,j,k,l,node,el,node1;

  standard_precision *rk0,*rk1;
  standard_precision *xk0,*pk0;
  standard_precision *Ap,*zk0,*zk1;
  standard_precision rdot1,rdot2,pdotAp,alpha,beta;
  standard_precision initial_res,res,res_mag;

  rk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  rk1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  xk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  pk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  zk0 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  zk1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  Ap = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));


#if 0
  if( !guess ) {
  
    for(p=1;p<=E->mesh.NNO[level];p++) {
      field[p] = 0.0;
    }
    
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
       if(E->tracer.property_group[m] < 0)
	continue;

      el = E->tracer.tracer_elt[level][m];
      for(q=1;q<=enodes[E->mesh.nsd];q++) {
	p = E->IEN[level][el].node[q];
	field[p] += E->tracer.tracer_int_weights[level][m].node[q] * property[m];
      }
      
    }
  

    if(BCinfo != NULL && BCval != NULL)
      for(i=1;i<=E->mesh.NNO[level];i++)
	if(BCinfo[i] & BCmask)
	  field[i] = BCval[i];
  }
#endif

  /* Least squares improvement of nodal value fit. (Try CG method) */

  /* Initialize */
 
  for(i=1;i<=E->mesh.NNO[level];i++) {
    rk0[i] = 0.0;
    xk0[i] = field[i];
    Ap[i] = 0.0;
  }

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
     if(E->tracer.property_group[m] < 0)
	continue;


    el = E->tracer.tracer_elt[level][m];
    for(j=1;j<=enodes[E->mesh.nsd];j++) {
      node = E->IEN[level][el].node[j];
      if((BCinfo != NULL && BCval != NULL) &&  (BCinfo[node] & BCmask)) 
	rk0[i] = 0.0;

      else {
	rk0[node] += E->tracer.sfn_values[level][m].node[j] * property[m];
      	Ap[node] += E->tracer.sfn_values[level][m].node[j] * property[m];

	for(p=1;p<=enodes[E->mesh.nsd];p++) {
	  node1 = E->IEN[level][el].node[p]; 
	  rk0[node] -= E->tracer.sfn_values[level][m].node[j] * 
	    E->tracer.sfn_values[level][m].node[p] * field[node1];
	} 
      }
    }
  }
  
  rdot1=0.0;
  res_mag=0.0;
  for(i=1;i<=E->mesh.NNO[level];i++) {
    zk0[i] = rk0[i] * E->tracer.ls_diag[level][i];
    rdot1 += rk0[i] * zk0[i];
    res_mag += Ap[i] * Ap[i];
  }

  initial_res = sqrt(rdot1)/E->mesh.NNO[level];
  res_mag = sqrt(res_mag)/E->mesh.NNO[level];

   /* CG loop */

  k=0;
  do {
    k++;
    for(i=1;i<=E->mesh.NNO[level];i++)
      zk0[i] = rk0[i] * E->tracer.ls_diag[level][i];


      /* Search directions */

      if(k==1) {
	for(i=1;i<=E->mesh.NNO[level];i++)
	 pk0[i] = rk0[i];
      }
      else {

	if(rdot2 != 0.0)
	  beta = rdot1 / rdot2;
	else {
	  beta = 1.0;
	}
	  
	for(i=1;i<=E->mesh.NNO[level];i++) {
	  if(0 && (E->NODE[level][i] & TBD))
	    pk0[i] = 0.0;
	  else
	    pk0[i] = zk1[i] + beta * pk0[i]; 
	}
      }

    /* get A*pk0 term */ 

    for(i=1;i<=E->mesh.NNO[level];i++)
      Ap[i] = 0.0;

    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
       el = E->tracer.tracer_elt[level][m];
      
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	node = E->IEN[level][el].node[j]; 
	if( (BCinfo != NULL && BCval != NULL) &&  (BCinfo[node] & BCmask))
	  Ap[node] = 0.0;
	else
	  for(p=1;p<=enodes[E->mesh.nsd];p++) {
	    node1 = E->IEN[level][el].node[p]; 
	    Ap[node] += E->tracer.sfn_values[level][m].node[j] *
	      E->tracer.sfn_values[level][m].node[p] * pk0[node1];
	  }
      }
    }

    /* update solution */

    pdotAp = 0.0;
    for(i=1;i<=E->mesh.NNO[level];i++) {
      pdotAp += pk0[i] * Ap[i];
    }
   
    if( pdotAp != 0.0)
      alpha = rdot1 / pdotAp;
    else
      alpha = 1.0;

    for(i=1;i<=E->mesh.NNO[level];i++) {
      xk0[i] += alpha * pk0[i];
      rk1[i] = rk0[i];
      zk1[i] = zk0[i];
      rk0[i] -= alpha * Ap[i];
    }
    
    /* test residual */

    rdot2=rdot1;
    rdot1=0.0;
    for(i=1;i<=E->mesh.NNO[level];i++)
      rdot1 += rk0[i] * zk0[i];

    res = sqrt(rdot1)/E->mesh.NNO[level];

  
  } while (k <= 100 && res > 0.01 * initial_res) ;

  
  fprintf(E->fp,"%d: Tracer 2 node residual [%d] = %g (%g) in %g\n",level,k,res,res/initial_res,res_mag );
    

  for(i=1;i<=E->mesh.NNO[level];i++)
    field[i] = xk0[i];

  free((void *) rk0);
  free((void *) rk1);
  free((void *) xk0);
  free((void *) pk0);
  free((void *) zk0);
  free((void *) zk1);
  free((void *) Ap);

  return;
}


void nodes_to_tracers (
 struct All_variables *E,
 standard_precision *field,
 standard_precision *property,
 int level
 )
{
  int i,el,j,node;
 
  standard_precision lN[ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega;
 
  const int dims = E->mesh.nsd;

  for(i=1;i<=E->tracer.NUM_TRACERS;i++) {

    el = E->tracer.tracer_elt[level][i];
    /* get_element_coords(E,el,i,E->tracer.tx,E->tracer.tz,
       E->tracer.ty,&eta1,&eta2,&eta3,level); */
    
    eta1 = E->tracer.eta1[level][i];
    eta2 = E->tracer.eta2[level][i];
    if(3==dims) {
      eta3 = E->tracer.eta3[level][i];
    }

    /* v_shape_fn(E,el,lN,eta1,eta2,eta3,level); */

    property[i] = 0.0;
    for(j=1;j<=enodes[E->mesh.nsd];j++) {
      node = E->IEN[level][el].node[j];
      property[i] += field[node] * E->tracer.sfn_values[level][i].node[j];
    }
  }
  return;
}





/*
  Note, that since this routine assumes tracer.S11 etc and tracer.Pt
  it should be called with temporary variables and the results copied
  to S11 ... , Pt later. Later we also compute the real pressure as we need
  divU on element.
*/


void tracer_deviatoric_stress_and_pressure (
 struct All_variables *E,
 standard_precision tau[7][7],
 standard_precision D[7][7],
 standard_precision *pressure,
 standard_precision *Dkk,
 standard_precision *Hvisc,
 int tracer,
 int level
)
{

  higher_precision *node_R;

  standard_precision K12,K21,K13,K31,K23,K32 ;
  standard_precision V12,V21,V31,V13,V23,V32;
  standard_precision Q,traceD;
  standard_precision eta1,eta2,eta3,dOmega,elastime,twovisc,elas1,velastime;
  standard_precision omega12,omega13,omega23;
  standard_precision lNx[4][ELNMAX+1],lN[ELNMAX+1];
  standard_precision U1,U2,U3;
  standard_precision t11,t12,t13,t22,t23,t33;
  standard_precision delta,a1,a0;
  standard_precision H,phi;
  standard_precision h1,h2,h3,h4,h5,h6;
  standard_precision visc1;
  standard_precision gamma1;
  standard_precision n1,n2,n3;

  const int dofs = E->mesh.dof ;
  const int dims = E->mesh.nsd ;
  const int ends = enodes[dims] ;

  int node,m,i;
  struct IEN *IEN = E->IEN[level];

  m = tracer;
  
  /* 2D Classical */

  /*fprintf(stderr,"C.O'N: In tracer_deviatoric_stress_and_pressure \n");
  fprintf(stderr,"C.O'N:  TESTING ELAS_SHEAR_MOD %g dofs %d dims %d \n",E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod,dofs,dims); */
  if(2==dims && 2==dofs) {
    
    /* 1. get deformation rate tensor, and split into deviatoric/volumetric */

    D[1][1] = D[2][2] = D[1][2] = 0.0 ;
    V12 = V21 = 0.0;
    Q = 0.0;
 
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
 
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

      /*  fprintf(stderr,"level %d: tracer %d: velocity: %g,%g, P %g,  node %d/%d (elt %d/%d,node %d)\n",
	  level,m,U1,U2,E->NQ[level][node],node,E->mesh.NNO[level],E->tracer.tracer_elt[level][m],E->mesh.NEL[level],i); */

      Q += E->tracer.sfn_values[level][m].node[i] * E->NQ[level][node];
      if (m == 436) {
      fprintf(stderr,"WHAT IS Q: sfn_value %g NQ %g Q %g level %d node %d m %d i %d \n",E->tracer.sfn_values[level][m].node[i],E->NQ[level][node],Q,level,node,m,i);
      }
      
      D[1][1] += lNx[1][i] * U1 ;
      D[2][2] += lNx[2][i] * U2 ;
      D[1][2] += 0.5 * (
	lNx[1][i] * U2 +
	lNx[2][i] * U1 );

      V12 += lNx[2][i] * U1 ;
      V21 += lNx[1][i] * U2 ;
    }
	
    omega12 = 0.5 * (V12 - V21);
    
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0)  /* ????? */
      traceD = (D[1][1] + D[2][2]);
    else
      traceD = 0.0;

    D[1][1] -= 0.5 * traceD;
    D[2][2] -= 0.5 * traceD; 

    *Dkk = traceD; 


    /* 2. Now the instantaneous stress from this */

    tau[1][1] = 2.0 * E->tracer.Visc[m] * D[1][1];
    tau[2][2] = 2.0 * E->tracer.Visc[m] * D[2][2];
    tau[1][2] = 2.0 * E->tracer.Visc[m] * D[1][2];


    /* NOTE: assumptions are that the elastic orthotropic case has the 
       same relaxation time for normal and shear moduli. 
       i.e. \mu_{ijkl} = \eta_{ijkl} / \alpha.

       This also means that the stored stress term is unchanged, but it is 
       a special case. If the elastic modulus is assumed constant, for example,
       then it is necessary to consider two relaxation times in this part, and 
       also in the stored stress update */

   if(E->control.CHEM_TRANS  && 
       E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio != 1.0 ) {

     H = E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio;
     phi = E->tracer.volfraction[m];
     
     delta = E->tracer.Visc[m] * (1.0 - H / ((phi * H + 1.0 - phi) * (phi + H - phi * H))); 
     
     a0 = 4.0 * delta *
       E->tracer.n1[m] *  E->tracer.n1[m] *  E->tracer.n2[m] *  E->tracer.n2[m];
     a1 = 2.0 * delta * E->tracer.n1[m] *  E->tracer.n2[m] *
       ( E->tracer.n2[m] *  E->tracer.n2[m] - E->tracer.n1[m] *  E->tracer.n1[m] );

      
      tau[1][1] += -a0 *  D[1][1] + a0 * D[2][2] - 2.0 * a1 * D[1][2];
      tau[2][2] +=  a0 *  D[1][1] - a0 * D[2][2] + 2.0 * a1 * D[1][2];
      tau[1][2] += -a1 *  D[1][1] + a1 * D[2][2] - 2.0 * (a0-delta) * D[1][2];

   }
   else if(E->control.ORTHOTROPY && 
	   E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio < 1.0) {
     
     delta = E->tracer.Visc[m] * (1.0 - E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio);
     
     a0 = 4.0 * delta *
       E->tracer.n1[m] *  E->tracer.n1[m] *  E->tracer.n2[m] *  E->tracer.n2[m];
     a1 = 2.0 * delta * E->tracer.n1[m] *  E->tracer.n2[m] *
       ( E->tracer.n2[m] *  E->tracer.n2[m] - E->tracer.n1[m] *  E->tracer.n1[m] );

     
     tau[1][1] += -a0 *  D[1][1] + a0 * D[2][2] - 2.0 * a1 * D[1][2];
     tau[2][2] +=  a0 *  D[1][1] - a0 * D[2][2] + 2.0 * a1 * D[1][2];
     tau[1][2] += -a1 *  D[1][1] + a1 * D[2][2] + 2.0 * (a0-delta) * D[1][2];
     
     H=phi=-1.0; 
   }

   if(fabs(tau[1][1]) > 1.0e10) {
     fprintf(stderr,"Tracer %d: Visc %g, Tau11=%g, n = %g,%g,  d = %g, a0 = %g, a1 =%g, H = %g, phi = %g\n",
	     m, E->tracer.Visc[m], tau[1][1],E->tracer.n1[m],E->tracer.n2[m],delta,a0,a1,H,phi);
   }
    if (m == 436)
	    fprintf(stderr,"HERE IS TRACER PRESSURE %g \n",*pressure);
    
    if(2==E->mesh.nsd) { /*RAA 9/7/02, will commenting out Q match up pressures with cIItcom?*/
	*pressure = Q  - (1.0 + E->tracer.visc[E->tracer.property_group[m]].Pen_bulk) * E->tracer.Visc[m] * traceD ;
	/* Following line from RAA's version, inserted FD's above instead). 
      *pressure = Q  - E->tracer.Visc[m] * traceD ; */
    }
    else if(3==E->mesh.nsd) /*RAA: this can never be accessed */
      *pressure = Q  - 0.666666666667 * E->tracer.Visc[m] * traceD;


    if (m == 436) {
	    fprintf(stderr,"HERE IS TRACER PRESSURE (2) %g Q %g Penalty %g Visc %g traceD %g\n",*pressure,Q,E->tracer.visc[E->tracer.property_group[m]].Pen_bulk,E->tracer.Visc[m],traceD);
    }

    if(E->control.ELASTICITY) {
      /* 3. and add the stress history to these values */
      
      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	elastime = 1.0/(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod * E->advection.elastic_timestep);
	
	if(E->control.deformation)
	  elas1 = 1.0 / (E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod);
	else
	  elas1 = 0.0;

	tau[1][1] += E->tracer.Visc[m] * (E->tracer.S11[m] * elastime + 2.0 * omega12 * E->tracer.S12[m] * elas1 );
        tau[2][2] += E->tracer.Visc[m] * (E->tracer.S22[m] * elastime - 2.0 * omega12 * E->tracer.S12[m] * elas1 );
        tau[1][2] += E->tracer.Visc[m] * (E->tracer.S12[m] * elastime +
						omega12 * (E->tracer.S22[m]-E->tracer.S11[m]) * elas1);
	/* C.O'N: taking out the t11 = tau[1][1] += */
      }
    }
    
    /* Probably should check that the Bulkvisc is > 0.0 otherwise a 
       mistake in parameter values could lead to catastrophic errors */

    if(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod > 0.0) {    
      *pressure += E->tracer.BulkVisc[m] / E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod *
	  (E->tracer.Pt[m] / E->advection.elastic_timestep );
	
      /* C.O'N: shouldn't make difference, but put in FD's above
       * *pressure += E->tracer.BulkVisc[m] * E->tracer.Pt[m] / 
	(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod * E->advection.elastic_timestep); */
    }
    if (m == 436) {
	    fprintf(stderr,"HERE IS TRACER PRESSURE (3) %g Bulk Visc %g Bulk Mod %g tracer.Pt %g elastic timestep %g \n",*pressure,E->tracer.BulkVisc[m],E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod,E->tracer.Pt[m],E->advection.elastic_timestep);
    }

  }

  /* ?? Cosserat ?? */


  /* 3D, classical */
  
  if(3==dims && 3==dofs) {

   /* fprintf(stderr,"Made it into 3D %g \n",E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod);	*/
    /* 1. get deformation rate tensor, and split into deviatoric/volumetric */

    D[1][1] = D[2][2] = D[3][3] = D[1][2] = D[1][3] = D[2][3] = 0.0 ;
    V12 = V21 = V13 = V23 = V32 = V31 = 0.0;
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    eta3 = E->tracer.eta3[level][m];

    /*RAA: 20/03/02, this Q=0 initialization was missing for 3D, but may
     *     not matter depending on the proper definition of *pressure */
    Q = 0.0; 

    get_global_v_x_shape_fn(E,E->tracer.tracer_elt[level][m],lNx,&dOmega,eta1,eta2,eta3,level);
	
    for(i=1;i<=ends;i++) {
      node = IEN[E->tracer.tracer_elt[level][m]].node[i];


      if(!(E->control.CURRENT_SKEWED_Vs && E->NODE[level][node] & SKEWBC)) {
	U1 = E->VV[level][1][node];
	U2 = E->VV[level][2][node]; 
 	U3 = E->VV[level][3][node]; 
      }
      else {
	node_R = E->curvilinear.NODE_R[level][node];
	U1 = node_R[0*dims+0]*E->VV[level][1][node] + node_R[0*dims+1]*E->VV[level][2][node] + node_R[0*dims+2] * E->VV[level][3][node];
	U2 = node_R[1*dims+0]*E->VV[level][1][node] + node_R[1*dims+1]*E->VV[level][2][node] + node_R[1*dims+2] * E->VV[level][3][node];
	U3 = node_R[2*dims+0]*E->VV[level][1][node] + node_R[2*dims+1]*E->VV[level][2][node] + node_R[2*dims+2] * E->VV[level][3][node];
      } 

      /* ?? */

      Q += E->tracer.sfn_values[level][m].node[i] * E->NQ[level][node];

      if (m == 1739) {
      fprintf(stderr,"WHAT IS Q: sfn_value %g NQ %g Q %g level %d node %d m %d i %d \n",E->tracer.sfn_values[level][m].node[i],E->NQ[level][node],Q,level,node,m,i);
      }

      
      D[1][1] += lNx[1][i] * U1;
      D[2][2] += lNx[2][i] * U2;
      D[3][3] += lNx[3][i] * U3;
      D[1][2] += 0.5 * (lNx[1][i] * U2 + lNx[2][i] * U1);
      D[1][3] += 0.5 * (lNx[1][i] * U3 + lNx[3][i] * U1);
      D[2][3] += 0.5 * (lNx[2][i] * U3 + lNx[3][i] * U2);

      V12 += lNx[2][i] * U1;
      V21 += lNx[1][i] * U2;
      V13 += lNx[3][i] * U1;
      V31 += lNx[1][i] * U3;
      V32 += lNx[2][i] * U3;
      V23 += lNx[3][i] * U3;
   }
		
    /*RAA: 20/03/02, why were there no 0.5 factors here?*/
    /* C.O'N: took out 0.5's to conform with FDs stuff - not sure who's right...
     * This is what's published for material spin tensor, though: */
    omega12 = 0.5*(V12 - V21);
    omega13 = 0.5*(V13 - V31);
    omega23 = 0.5*(V23 - V32);
    
    /*RAA: 20/03/02, the following lines present for 2D but were missing for 3D 
    if(E->tracer.BulkVisc[m] > 0.0) 
      traceD = (D[1][1] + D[2][2] + D[3][3]);
    else
      traceD = 0.0;
	*/
    traceD = (D[1][1] + D[2][2] + D[3][3]);
          
    /*  traceD = (D[1][1] + D[2][2] + D[3][3]); */ /*RAA - comment this out*/
    
    D[1][1] -= 0.333333333333 * traceD;
    D[2][2] -= 0.333333333333 * traceD;
    D[3][3] -= 0.333333333333 * traceD;
 
    *Dkk = traceD; 

    /* 2. Now the instantaneous stress from this lot */

    tau[1][1] = 2.0 * E->tracer.Visc[m] * D[1][1];
    tau[2][2] = 2.0 * E->tracer.Visc[m] * D[2][2];
    tau[3][3] = 2.0 * E->tracer.Visc[m] * D[3][3];
    tau[1][2] = 2.0 * E->tracer.Visc[m] * D[1][2];
    tau[1][3] = 2.0 * E->tracer.Visc[m] * D[1][3];
    tau[2][3] = 2.0 * E->tracer.Visc[m] * D[2][3];

    /*RAA: 15/03/02, added comments around Q to conform with 2D case*/
    /* C.O'N: took out comments to conform with FDs stuff  - put back in now... 
     * Also changed BulkVisc to plain Visc below
     * Also now put in the 2d line ie. 0.666666 -> 1+Pen_bulk */
    if (m == 1739)
	    fprintf(stderr,"HERE IS TRACER PRESSURE %g \n",*pressure);

    *pressure = Q -  0.666666666666666666666666 * E->tracer.Visc[m] * traceD;

    if (m == 1739)
	    fprintf(stderr,"HERE IS TRACER PRESSURE %g Q %g Visc %g and traceD %g\n",*pressure,Q,E->tracer.Visc[m],traceD);

    
    /*fprintf(stderr,"Bout to head into Elasticity %g \n",E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod);	*/
    if(E->control.ELASTICITY) {
      /* 3. and add the stress history to these values */

      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod != 0.0) {
	elastime = 1.0 / (E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod * E->advection.elastic_timestep);

	/* C.O'N: changed this: */
	if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod != 0.0)
	  elas1 = 1.0 / (E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod);
	else
	  elas1 = 0.0;


/*	fprintf(stderr,"What is Elas_shear_mod? %g \n",E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod);	
 *	These are equation 13 in Moresi et al. (ACES_WS2 extended abstract) */
	tau[1][1] += E->tracer.Visc[m] * (E->tracer.S11[m] * elastime + 
					  (omega12 * E->tracer.S12[m] + omega13*E->tracer.S13[m]) * elas1);
      
	tau[2][2] += E->tracer.Visc[m] * (E->tracer.S22[m] * elastime + 
					  (omega23 * E->tracer.S23[m] - omega12*E->tracer.S12[m]) * elas1);
	
	tau[3][3] += E->tracer.Visc[m] * (E->tracer.S33[m] * elastime + 
					  (omega23 * E->tracer.S23[m] + omega13*E->tracer.S13[m]) * elas1);
      
	tau[1][2] += E->tracer.Visc[m] * (E->tracer.S12[m] * elastime + 
					  (omega12 * (E->tracer.S22[m]-E->tracer.S11[m])  +
					   omega13 * E->tracer.S23[m] + 
					   omega23 * E->tracer.S13[m] ) * elas1 );
	
	tau[1][3] += E->tracer.Visc[m] * (E->tracer.S13[m] * elastime + 
					  (omega13 * (E->tracer.S33[m]-E->tracer.S11[m])  +
					   omega12 * E->tracer.S23[m] - 
					   omega23 * E->tracer.S12[m] ) * elas1 );
	
	tau[2][3] += E->tracer.Visc[m] * (E->tracer.S23[m] * elastime + 
					  (omega23 * (E->tracer.S33[m]-E->tracer.S22[m])  +
					   omega13 * E->tracer.S12[m] - 
					   omega12 * E->tracer.S13[m] ) * elas1 );
      }
    }
    /* This is equation 14 in Moresi et al. */
    if(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod > 0.0) {    
      *pressure += E->tracer.BulkVisc[m] * E->tracer.Pt[m] / 
	(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod * E->advection.elastic_timestep);
    }
  
    if (m == 1739)
	    fprintf(stderr,"HERE IS TRACER PRESSURE %g \n",*pressure);

    /* Need to do the time-average stabilization for 3D as well */

  }

  /* Viscous shear heating - in VE case, requires an estimate of the
     viscous strain rate (as opposed to viscoelastic total strain rate) - this comes
     about via the stress and inversion of the original viscosity tensor */

  if(Hvisc != NULL) {
    if (E->control.SHEAR_HEATING) { 

    /*RAA: - skip the 3D version of viscous shear heating - probably won't be needed*/
    if(3==dims) {
     fprintf(stderr,"3D shear heating not properly coded. Adios, muchacho.\n");
     exit(1);
    }

      /* OK, the 2D case */

      visc1 = 0.5 / E->tracer.Visc0[m];
    
      h1 = visc1 * tau[1][1] * tau[1][1];
      h2 = visc1 * tau[2][2] * tau[2][2];
      h3 = visc1 * tau[1][2] * tau[1][2];
    
      if (E->control.ORTHOTROPY && 
	  E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio < 1.0) {
     
	gamma1 = (1-E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio) /
	  E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio;
    
	n1 = E->tracer.n1[m];
	n2 = E->tracer.n2[m];  
    
	h1 += visc1 * tau[1][1] * gamma1 * 2.0 * ( n1*n1*(1.0-n1*n1) * tau[1][1] - 
						   n1*n1*n2*n2 * tau[2][2] + 
						   n1*n2*(1-2.0*n1*n1) * tau[1][2] );
      
	h2 += visc1 * tau[2][2] * gamma1 * 2.0 * ( n2*n2*(1.0-n2*n2) * tau[2][2] - 
						   n1*n1*n2*n2 * tau[1][1] + 
						   n1*n2*(1-2.0*n2*n2) * tau[1][2] );
 
	h3 += visc1 * tau[1][2] * gamma1 * ( n1*n2*(1.0-2.0*n1*n1) * tau[1][1] +
					     n1*n2*(1.0-2.0*n2*n2) * tau[2][2] +
					     (n1*n1+n2*n2-4.0*n1*n1*n2*n2) * tau[1][2] );

      }

      if (h1 < 0.0 || h2 < 0.0 || h3 < 0.0)
	fprintf(stderr,"Tracer %d: visc heating term in anisotropic case is -ve: %g,%g,%g\n",m,h1,h2,h3);
  
      *Hvisc = E->tracer.Hv0[E->tracer.property_group[m]] * (h1 + h2 + h3 + h3);
    }
    else {  /* No shear heating, set the value to zero anyway, for safety */
      *Hvisc = 0.0;
    }

  }

  /* fprintf(stderr,"HERE IS TRACER.PT (C)   %g \n",E->tracer.Pt[436]) ;	*/
 
  return;
  
}


void get_D (
	    struct All_variables *E,
	    standard_precision tau[7][7],
	    standard_precision D[7][7],
	    standard_precision *Dkk,
	    int tracer,
	    int level
	    )
{
  higher_precision *node_R;
  
  standard_precision V12,V21,V31,V13,V23,V32;
  standard_precision traceD;
  standard_precision eta1,eta2,eta3,dOmega,elastime;
  standard_precision omega12,omega13,omega23;
  standard_precision lNx[4][ELNMAX+1];
  standard_precision U1,U2,U3;

  const int dims = E->mesh.nsd ;
  const int ends = enodes[dims] ;

  int node,m,i;
  struct IEN *IEN = E->IEN[level];

  m = tracer;
  
  /* 2D Classical */
  
  if(2==dims) {
    /* 1. get deformation rate tensor, and split into deviatoric/volumetric */
    
    D[1][1] = D[2][2] = D[1][2] = 0.0 ;
    V12 = V21 = 0.0 ;
    
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    
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
      
      D[1][1] += lNx[1][i] * U1 ;
      D[2][2] += lNx[2][i] * U2 ;
      D[1][2] += 0.5 * (lNx[1][i] * U2 + lNx[2][i] * U1 );
      
      V12 += lNx[2][i] * U1 ;
      V21 += lNx[1][i] * U2 ;
    }
    
    omega12 = 0.5 * (V12 - V21);
    
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0)
      traceD = (D[1][1] + D[2][2]);
    else
      traceD = 0.0;
    
    D[1][1] -= 0.5 * traceD;
    D[2][2] -= 0.5 * traceD;

    *Dkk = traceD;
    /* 2. Now the instantaneous stress from this */

    tau[1][1] = D[1][1];
    tau[2][2] = D[2][2];
    tau[1][2] = D[1][2];

    if(E->control.ELASTICITY) {
      /* 3. and add the stress history to these values */
      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	elastime = 1.0/(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod * E->advection.elastic_timestep);
	
	tau[1][1] += E->tracer.S11[m] * elastime + 2.0 * omega12 * E->tracer.S12[m]/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
        tau[2][2] += E->tracer.S22[m] * elastime - 2.0 * omega12 * E->tracer.S12[m]/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod ;
        tau[1][2] += E->tracer.S12[m] * elastime + omega12 * (E->tracer.S22[m]-E->tracer.S11[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
      }
    }
  }
  else if(3==dims) {

    /* 1. get deformation rate tensor, and split into deviatoric/volumetric */

    D[1][1] = D[2][2] = D[3][3] = D[1][2] = D[1][3] = D[2][3] = 0.0 ;
    V12 = V21 = V13 = V23 = V32 = V31 = 0.0;
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    eta3 = E->tracer.eta3[level][m];

    get_global_v_x_shape_fn(E,E->tracer.tracer_elt[level][m],lNx,&dOmega,eta1,eta2,eta3,level);
	
    for(i=1;i<=ends;i++) {
      node = IEN[E->tracer.tracer_elt[level][m]].node[i];

      if(!(E->control.CURRENT_SKEWED_Vs && E->NODE[level][node] & SKEWBC)) {
	U1 = E->VV[level][1][node];
	U2 = E->VV[level][2][node]; 
 	U3 = E->VV[level][3][node]; 
      }
      else {
	node_R = E->curvilinear.NODE_R[level][node];
	U1 = node_R[0*dims+0]*E->VV[level][1][node] + node_R[0*dims+1]*E->VV[level][2][node] + node_R[0*dims+2] * E->VV[level][3][node];
	U2 = node_R[1*dims+0]*E->VV[level][1][node] + node_R[1*dims+1]*E->VV[level][2][node] + node_R[1*dims+2] * E->VV[level][3][node];
	U3 = node_R[2*dims+0]*E->VV[level][1][node] + node_R[2*dims+1]*E->VV[level][2][node] + node_R[2*dims+2] * E->VV[level][3][node];
      } 

      D[1][1] += lNx[1][i] * U1;
      D[2][2] += lNx[2][i] * U2;
      D[3][3] += lNx[3][i] * U3;
      D[1][2] += 0.5 * (lNx[1][i] * U2 + lNx[2][i] * U1);
      D[1][3] += 0.5 * (lNx[1][i] * U3 + lNx[3][i] * U1);
      D[2][3] += 0.5 * (lNx[2][i] * U3 + lNx[3][i] * U2);

      V12 += lNx[2][i] * U1;
      V21 += lNx[1][i] * U2;
      V13 += lNx[3][i] * U1;
      V31 += lNx[1][i] * U3;
      V32 += lNx[2][i] * U3;
      V23 += lNx[3][i] * U3;
   }
		
    omega12 = (V12 - V21);
    omega13 = (V13 - V31);
    omega23 = (V23 - V32);
    
    traceD = (D[1][1] + D[2][2] + D[3][3]);
    D[1][1] -= 0.333333333333 * traceD;
    D[2][2] -= 0.333333333333 * traceD;
    D[3][3] -= 0.333333333333 * traceD;
 
    *Dkk = traceD; 

    /* 2. Now the instantaneous stress from this lot */

    tau[1][1] = D[1][1];
    tau[2][2] = D[2][2];
    tau[3][3] = D[3][3];
    tau[1][2] = D[1][2];
    tau[1][3] = D[1][3];
    tau[2][3] = D[2][3];

    if(E->control.ELASTICITY) {
      /* 3. and add the stress history to these values */

      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod != 0.0) {
	elastime = 1.0 / (E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod * E->advection.elastic_timestep);

	tau[1][1] += E->tracer.S11[m] * elastime + (omega12 * E->tracer.S12[m] + omega13*E->tracer.S13[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	tau[2][2] += E->tracer.S22[m] * elastime + (omega23 * E->tracer.S23[m] - omega12*E->tracer.S12[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	tau[3][3] += E->tracer.S33[m] * elastime + (omega23 * E->tracer.S23[m] + omega13*E->tracer.S13[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	tau[1][2] += E->tracer.S12[m] * elastime + 
	  (omega12 * (E->tracer.S22[m]-E->tracer.S11[m]) + omega13 * E->tracer.S23[m] + omega23 * E->tracer.S13[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	tau[1][3] += E->tracer.S13[m] * elastime + 
	  (omega13 * (E->tracer.S33[m]-E->tracer.S11[m]) + omega12 * E->tracer.S23[m] - omega23 * E->tracer.S12[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	tau[2][3] += E->tracer.S23[m] * elastime + 
	  (omega23 * (E->tracer.S33[m]-E->tracer.S22[m]) + omega13 * E->tracer.S12[m] - omega12 * E->tracer.S13[m])/
	  E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
      }
    }
  }
  return;
}
