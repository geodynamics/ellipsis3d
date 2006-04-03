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


/*========================================================
  Function to make the IEN array for a mesh of given 
  dimension. IEN is an externally defined structure array

  NOTE: this is not really general enough for new elements:
  it should be done through a pre-calculated lookup table.
  ======================================================== */

void construct_ien(
     struct All_variables *E
)
{	
  int lev,p,q,r,rr,e1,e2,i,a,node,node2,e;
  int element,start,nel,nno;
  int pmax,qmax,rmax;
  int pnmax,qnmax,rnmax;
  
  const int  dims=E->mesh.nsd;
  const int  dofs=E->mesh.dof;
  const int ends=enodes[dims];

  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {

    if(E->control.verbose) {
      fprintf(stderr,"Constructing IEN arrays for level %d\n",lev);
    }

      pmax = E->mesh.ELZ[lev];
      qmax = E->mesh.ELX[lev];
      rmax = E->mesh.ELY[lev];
      pnmax = E->mesh.NOZ[lev];
      qnmax = E->mesh.NOX[lev];
      rnmax = E->mesh.NOY[lev];
      nel=E->mesh.NEL[lev];
      nno=E->mesh.NNO[lev];

      for(p=1;p<=pmax;p++)
	for(q=1;q<=qmax;q++)
	      for(r=1;r<=rmax;r++) {
		  element = (r-1)*pmax*qmax + (q-1)*pmax  + p;	
		  start = (r-1)*pnmax*qnmax + (q-1)*pnmax + p; /*assumes 2 nodes per edge*/
		  for(rr=1;rr<=ends;rr++)
		      E->IEN[lev][element].node[rr]= start 
			  + offset[rr].vector[1]
			  + offset[rr].vector[0]*pnmax
			  + offset[rr].vector[2]*pnmax*qnmax;
	      }

      
      for(element=1;element<=nel;element++)
	for(node=1;node<=ends;node++) {
	  E->IENP[lev][element].node[node] =  E->IEN[lev][element].node[node];
	}
      

      /*RAA: 26/10/01, added !per_y to line directly below to help in organizing*/
      if(E->mesh.periodic_x && !E->mesh.periodic_y)  /* make end nodes vanish */
	  for(p=1;p<=pmax;p++)
	      for(r=1;r<=rmax;r++)  {
		  e1=p+(r-1)*E->mesh.ELX[lev]*E->mesh.ELZ[lev];
		  e2=e1+(E->mesh.ELX[lev]-1)*E->mesh.ELZ[lev];

		  E->NODE[lev][E->IEN[lev][e1].node[1]] = E->NODE[lev][E->IEN[lev][e1].node[1]] | PER_OFFSIDE;
		  E->NODE[lev][E->IEN[lev][e1].node[2]] = E->NODE[lev][E->IEN[lev][e1].node[2]] | PER_OFFSIDE;
		  
		  E->IEN[lev][e2].node[4] = E->IEN[lev][e1].node[1];
		  E->IEN[lev][e2].node[3] = E->IEN[lev][e1].node[2];
		  
		  if(3==E->mesh.nsd) { /*RAA: 31/5/01,26/10/01 next 4 lines added for proper PER_OFFSIDE-ing*/
		      E->IEN[lev][e2].node[4] = E->IEN[lev][e1].node[1];
		      E->IEN[lev][e2].node[3] = E->IEN[lev][e1].node[2];
		      E->NODE[lev][E->IEN[lev][e1].node[5]] = E->NODE[lev][E->IEN[lev][e1].node[5]] | PER_OFFSIDE;
		      E->NODE[lev][E->IEN[lev][e1].node[6]] = E->NODE[lev][E->IEN[lev][e1].node[6]] | PER_OFFSIDE;
		      E->IEN[lev][e2].node[8] = E->IEN[lev][e1].node[5];
		      E->IEN[lev][e2].node[7] = E->IEN[lev][e1].node[6];
		  }
	      }
	  
      if(3==dims) { /*RAA: 26/10/01, added !per_x to line directly below to help in organizing*/
	  if(E->mesh.periodic_y && !E->mesh.periodic_x) /* end nodes vanish */
	      for(p=1;p<=pmax;p++)
		  for(q=1;q<=qmax;q++) {
		      e1=p+(q-1)*E->mesh.ELZ[lev];
		      e2=e1+(E->mesh.ELY[lev]-1)*E->mesh.ELZ[lev]*E->mesh.ELX[lev];

		      E->NODE[lev][E->IEN[lev][e1].node[4]] = E->NODE[lev][E->IEN[lev][e1].node[4]] | PER_OFFSIDE;
		      E->NODE[lev][E->IEN[lev][e1].node[3]] = E->NODE[lev][E->IEN[lev][e1].node[3]] | PER_OFFSIDE;
		      E->NODE[lev][E->IEN[lev][e1].node[2]] = E->NODE[lev][E->IEN[lev][e1].node[2]] | PER_OFFSIDE;
		      E->NODE[lev][E->IEN[lev][e1].node[1]] = E->NODE[lev][E->IEN[lev][e1].node[1]] | PER_OFFSIDE;

		      E->IEN[lev][e2].node[5] = E->IEN[lev][e1].node[1];
		      E->IEN[lev][e2].node[6] = E->IEN[lev][e1].node[2];
		      E->IEN[lev][e2].node[7] = E->IEN[lev][e1].node[3];
		      E->IEN[lev][e2].node[8] = E->IEN[lev][e1].node[4];
		}    
      }
    
      /* & if both ... ? */
      /*RAA: 26/10/01, let's try to get this right, in spite of the ambiguity*/
     if(3==dims && E->mesh.periodic_y && E->mesh.periodic_x) { 
          /*periodic_x part*/
	  for(p=1;p<=pmax;p++)
	      for(r=1;r<=rmax;r++)  {
		  e1=p+(r-1)*E->mesh.ELX[lev]*E->mesh.ELZ[lev];
		  e2=e1+(E->mesh.ELX[lev]-1)*E->mesh.ELZ[lev];

		  E->NODE[lev][E->IEN[lev][e1].node[1]] = E->NODE[lev][E->IEN[lev][e1].node[1]] | PER_OFFSIDE;
		  E->NODE[lev][E->IEN[lev][e1].node[2]] = E->NODE[lev][E->IEN[lev][e1].node[2]] | PER_OFFSIDE;
		  E->NODE[lev][E->IEN[lev][e1].node[5]] = E->NODE[lev][E->IEN[lev][e1].node[5]] | PER_OFFSIDE;
		  E->NODE[lev][E->IEN[lev][e1].node[6]] = E->NODE[lev][E->IEN[lev][e1].node[6]] | PER_OFFSIDE;
		  
		  E->IEN[lev][e2].node[4] = E->IEN[lev][e1].node[1];
		  E->IEN[lev][e2].node[3] = E->IEN[lev][e1].node[2];
		  E->IEN[lev][e2].node[8] = E->IEN[lev][e1].node[5];
		  E->IEN[lev][e2].node[7] = E->IEN[lev][e1].node[6];
	      }
          /*periodic_y part*/
	  for(p=1;p<=pmax;p++)
              for(q=1;q<=qmax;q++) {
	          e1=p+(q-1)*E->mesh.ELZ[lev];
	          e2=e1+(E->mesh.ELY[lev]-1)*E->mesh.ELZ[lev]*E->mesh.ELX[lev];

	          E->NODE[lev][E->IEN[lev][e1].node[4]] = E->NODE[lev][E->IEN[lev][e1].node[4]] | PER_OFFSIDE;
	          E->NODE[lev][E->IEN[lev][e1].node[3]] = E->NODE[lev][E->IEN[lev][e1].node[3]] | PER_OFFSIDE;
	          E->NODE[lev][E->IEN[lev][e1].node[2]] = E->NODE[lev][E->IEN[lev][e1].node[2]] | PER_OFFSIDE;
	          E->NODE[lev][E->IEN[lev][e1].node[1]] = E->NODE[lev][E->IEN[lev][e1].node[1]] | PER_OFFSIDE;

	          E->IEN[lev][e2].node[5] = E->IEN[lev][e1].node[1];
	          E->IEN[lev][e2].node[6] = E->IEN[lev][e1].node[2];
	          E->IEN[lev][e2].node[7] = E->IEN[lev][e1].node[3];
	          E->IEN[lev][e2].node[8] = E->IEN[lev][e1].node[4];
	      }    
        /*RAA: Now at this point the back right edge nodes are not PER_OFFSIDE, 
          and to make them so has caused the program to continously loop.
          So here I equate the front right edge nodes with those at the 
          back right edge without making those nodes PER_OFFSIDE. Is this ok?*/
  	  for(p=1;p<=pmax;p++) {
	      e1=p+(E->mesh.ELX[lev]-1)*E->mesh.ELZ[lev];
	      e2=e1+(E->mesh.ELY[lev]-1)*E->mesh.ELZ[lev]*E->mesh.ELX[lev];

	  /*    E->NODE[lev][E->IEN[lev][e1].node[3]] = E->NODE[lev][E->IEN[lev][e1].node[3]] | PER_OFFSIDE;
	      E->NODE[lev][E->IEN[lev][e1].node[4]] = E->NODE[lev][E->IEN[lev][e1].node[4]] | PER_OFFSIDE;
           */ 
	      E->IEN[lev][e2].node[7] = E->IEN[lev][e1].node[3] + (E->mesh.NOX[lev]-1)*(E->mesh.NOZ[lev]);
	      E->IEN[lev][e2].node[8] = E->IEN[lev][e1].node[4] + (E->mesh.NOX[lev]-1)*(E->mesh.NOZ[lev]);
          }
     } /*end of per_x and per_y*/

/*RAA: 31/5/01, check PER_OFFSIDE stuff*/  
/*    if(E->control.verbose) {
      for(p=1;p<=E->mesh.nno;p++)
        if(E->node[p] & (PER_OFFSIDE))
          fprintf(stderr,"(construct_mesh_arrays.c)!! A PER_OFFSIDE NODE!!: %d\n",p);
    }
*/
  /*RAA: 31/5/01, check PER_OFFSIDE stuff*/  
/*    if(E->control.verbose) {
      for(p=1;p<=E->mesh.nel;p++)
          fprintf(stderr,"(construct_mesh_arrays.c)!! element #:  %d, gn for ln7: %d , gn for ln8: %d\n",p,E->IEN[0][p].node[7],E->IEN[0][p].node[8]);
    }
*/
  
      for(i=1;i<=nno;i++)
	  E->NEI[lev].nels[i] = 0; 

      for(e=1;e<=nel;e++)
	  for(a=1;a<=ends;a++) {
	      node=E->IEN[lev][e].node[a];
	      E->NEI[lev].nels[node]++;
	      E->NEI[lev].element[(node-1)*ends+E->NEI[lev].nels[node]-1] = e;
	      E->NEI[lev].lnode[(node-1)*ends+E->NEI[lev].nels[node]-1] = a;
	  }   
  
      if(E->mesh.periodic_x) {  /* this should be redundant, the offside nodes are not used */
	  for(p=1;p<=pnmax;p++)
	      for(r=1;r<=rnmax;r++)  {
		  node=p+(r-1)*pnmax*qnmax;
		  node2=node+(qnmax-1)*pnmax;
		  E->NEI[lev].nels[node2]=E->NEI[lev].nels[node];
		  
		  for(a=1;a<=E->NEI[lev].nels[node2];a++) {
		      E->NEI[lev].element[(node2-1)*ends+a] = E->NEI[lev].element[(node-1)*ends+a];
		      E->NEI[lev].lnode[(node2-1)*ends+a] = E->NEI[lev].lnode[(node-1)*ends+a];
		  }
	      }	      
      }
      if(E->control.verbose) {
      fprintf(stderr,"Contructed IEN arrays for level %d\n",lev);
    }
  } /* next level */
  return;
}

/*============================================
  Function to make the ID array for above case
  ============================================ */

void construct_id(
		  struct All_variables *E
		  )
{ 
    int i,j,k;	
    int eqn_count,node_count,eqn_countd;
    unsigned int type,doff;
    int lev;
    int nox,noy,noz;

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;

    for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
      eqn_count = 0;
      nox=E->mesh.NOX[lev];
      noz=E->mesh.NOZ[lev];
      noy=E->mesh.NOY[lev];
      for(k=1;k<=noy;k++)
	for(i=1;i<=nox;i++)
	  for(j=1;j<=noz;j++)
	    for(doff=1;doff<=dofs;doff++) {
	      if(E->mesh.periodic_x && i==nox)  
		continue;
	      if(E->mesh.periodic_y && k==noy) /*RAA: 8/10/01, added these 2 lines */  
		continue;
	      /*RAA: needs correction for both periodic_x and periodic_y, too?, should be ok as is */
	      node_count = j+(i-1)*noz+(k-1)*nox*noz;
	      E->ID[lev][node_count].doff[doff] = eqn_count++;
	    }
      if(E->mesh.periodic_x)
	for(k=1;k<=noy;k++)
	  for(j=noz;j>=1;j--) {
	    E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[1] = E->ID[lev][j+(k-1)*nox*noz].doff[1];
	    E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[2] = E->ID[lev][j+(k-1)*nox*noz].doff[2];
	    if(3==dofs)
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[3] = E->ID[lev][j+(k-1)*nox*noz].doff[3];
	    else if(6==dofs) {
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[3] = E->ID[lev][j+(k-1)*nox*noz].doff[3];
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[4] = E->ID[lev][j+(k-1)*nox*noz].doff[4];
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[5] = E->ID[lev][j+(k-1)*nox*noz].doff[5];
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[6] = E->ID[lev][j+(k-1)*nox*noz].doff[6];
	    }
	  }
      if(E->mesh.periodic_y && 3==dofs)  /*RAA: 8/10/01, added this periodic_y stuff */  
	for(i=1;i<=nox;i++)  
	  for(j=noz;j>=1;j--) {
	    E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[1] = E->ID[lev][j+(i-1)*noz].doff[1];
	    E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[2] = E->ID[lev][j+(i-1)*noz].doff[2];
	    E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[3] = E->ID[lev][j+(i-1)*noz].doff[3];
	    if(6==dofs) {
	      E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[4] = E->ID[lev][j+(i-1)*noz].doff[4];
	      E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[5] = E->ID[lev][j+(i-1)*noz].doff[5];
	      E->ID[lev][j+(i-1)*noz+(noy-1)*nox*noz].doff[6] = E->ID[lev][j+(i-1)*noz].doff[6];
	    }
	  }
      /*RAA: 22/10/01, probably needs a correction for both periodic_x and periodic_y, 
         so now we will RE-DO THE PERIODIC_X STEP, eg., see flogical_mesh_to_real( ) */
      if(E->mesh.periodic_x && E->mesh.periodic_y && 3==dofs)  
	for(k=1;k<=noy;k++)
	  for(j=noz;j>=1;j--) {
	    E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[1] = E->ID[lev][j+(k-1)*nox*noz].doff[1];
	    E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[2] = E->ID[lev][j+(k-1)*nox*noz].doff[2];
	    E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[3] = E->ID[lev][j+(k-1)*nox*noz].doff[3];
	    if(6==dofs) {
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[4] = E->ID[lev][j+(k-1)*nox*noz].doff[4];
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[5] = E->ID[lev][j+(k-1)*nox*noz].doff[5];
	      E->ID[lev][j+(nox-1)*noz+(k-1)*nox*noz].doff[6] = E->ID[lev][j+(k-1)*nox*noz].doff[6];
	    }
	  }

      E->mesh.NEQ[lev] = eqn_count;
    }
    
    E->mesh.neq = E->mesh.NEQ[E->mesh.levmax];  /*  Total NUMBER of independent variables  */
 
    /* Now do the low level direct solver stuff */

    eqn_countd = 0;

     for(j=1;j<=E->mesh.NOZ[E->mesh.levmin];j++)
      for(k=1;k<=E->mesh.NOY[E->mesh.levmin];k++)
	for(i=1;i<=E->mesh.NOX[E->mesh.levmin];i++)	{ 
	  if(E->mesh.periodic_x && i==E->mesh.NOX[E->mesh.levmin])  
	      continue;
	  if(E->mesh.periodic_y && k==E->mesh.NOY[E->mesh.levmin]) /*RAA: 8/10/01, added these 2 lines */   
	      continue;
	  /*RAA: needs correction for both periodic_x and periodic_y, too?  should be ok as is.*/
	  node_count = j + (i-1)*E->mesh.NOZ[E->mesh.levmin]+
	    (k-1)*E->mesh.NOX[E->mesh.levmin] * E->mesh.NOZ[E->mesh.levmin];
	  
	  if((E->NODE[E->mesh.levmin][node_count] & BC1) == 0)
	    E->idd[node_count].doff[1] = eqn_countd++;
	  else
	    E->idd[node_count].doff[1] = -1;
	
	  if((E->NODE[E->mesh.levmin][node_count] & BC2) == 0)
	    E->idd[node_count].doff[2] = eqn_countd++;
	  else
	    E->idd[node_count].doff[2] = -1;

	  if(3==E->mesh.dof)  {
	    if((E->NODE[E->mesh.levmin][node_count] & BC3) == 0)
	      E->idd[node_count].doff[3] = eqn_countd++;
	    else
	      E->idd[node_count].doff[3] = -1;
	  }

	}
    E->mesh.neqd = eqn_countd; 
    return;
 }

/*==========================================================
  Function to construct  the LM array from the ID and IEN arrays 
  ========================================================== */

void construct_lm(
     struct All_variables *E
)
{	
  return;	
}


/* ============================================
   Function to set up the boundary condition
   masks and other indicators.
   ============================================  */

void construct_masks(		/* Add lid/edge masks/nodal weightings */
     struct All_variables *E
)
{	
  int i,j,k,l,node,elt,n1,n2;
  int lev,elx,elz,ely;
  
  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--){
      elx = E->mesh.ELX[lev];
      elz = E->mesh.ELZ[lev];
      ely = E->mesh.ELY[lev];
 
      if(E->mesh.periodic_x)
	  for(i=1;i<=E->mesh.NOZ[lev];i++)
	      for(j=1;j<=E->mesh.NOY[lev];j++) {
		  n1=i+(j-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
		  n2=n1+(E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
  
		  E->NODE[lev][n2] = E->NODE[lev][n2] | OFFSIDE;
	      } 

    /*RAA: 5/10/01, added periodic_y stuff below to create proper OFFSIDE nodes on front face*/
      if(E->mesh.periodic_y)
	  for(i=1;i<=E->mesh.NOX[lev];i++)
	      for(j=1;j<=E->mesh.NOZ[lev];j++) {
		  n1=j+(i-1)*E->mesh.NOZ[lev];
		  n2=n1+(E->mesh.NOY[lev]-1)*E->mesh.NOZ[lev]*E->mesh.NOX[lev];
  
		  E->NODE[lev][n2] = E->NODE[lev][n2] | OFFSIDE;
	      } 
    /*RAA: end of periodic_y insert*/

      for(i=1;i<=E->mesh.NNO[lev];i++)
	  E->TW[lev][i] = 0.0;

      for(i=1;i<=elz;i++)
	  for(j=1;j<=elx;j++)
	      for(k=1;k<=ely;k++)
		  for(l=1;l<=enodes[E->mesh.nsd];l++) {
		      elt =  i + (j-1) * elz + (k-1) * elz * elx;
		      node = E->IEN[lev][elt].node[l];
		      E->TW[lev][node] += 1.0;
		  }

      for(i=1;i<=E->mesh.NNO[lev];i++) {
	  if(E->NODE[lev][i] & ( OFFSIDE  ))
	      continue;

	  if( E->TW[lev][i] == 0.0 )
	      fprintf(stderr,"Weightings broken at level %d, node %d\n",lev,i);
	  E->TW[lev][i] = 1.0/(E->TW[lev][i]); 
      }      
  }
 				/* Edge masks  */
  
  for(i=1;i<=E->mesh.nox;i++)	/* Horizontal  */
    { for(j=1;j<=E->mesh.noy;j++)	
	{ node = 1+(i-1)*E->mesh.noz+(j-1)*E->mesh.noz*E->mesh.nox;
	  E->node[node] = E->node[node] | TZEDGE;
	  E->node[node] = E->node[node] | VZEDGE;
	  node += E->mesh.noz-1;;
	  E->node[node] = E->node[node] | VZEDGE;
	  E->node[node] = E->node[node] | TZEDGE;  
	}
    }

  if (E->mesh.dof == 3) /* not appropriate otherwise */
    for(i=1;i<=E->mesh.noz;i++)	/* vertical edge, x normal */
      { for(j=1;j<=E->mesh.noy;j++)	
	  { node = i + (j-1) * E->mesh.nox * E->mesh.noz;
	    E->node[node] = E->node[node] | TXEDGE;
	    E->node[node] = E->node[node] | VXEDGE;
	    node = i+(E->mesh.nox-1)*E->mesh.noz + (j-1) * E->mesh.nox * E->mesh.noz;
	    E->node[node] = E->node[node] | TXEDGE; 
	    E->node[node] = E->node[node] | VXEDGE; } }

  for(i=1;i<=E->mesh.noz;i++)	/* vertical edge, y normal */
    { for(j=1;j<=E->mesh.nox;j++)	
	{ node = i + (j-1) * E->mesh.noz;
	  E->node[node] = E->node[node] | TYEDGE;
	  E->node[node] = E->node[node] | VYEDGE;
	  node = i+(E->mesh.noy-1)*E->mesh.noz*E->mesh.nox + (j-1) * E->mesh.noz;
	  E->node[node] = E->node[node] | TYEDGE; 
	  E->node[node] = E->node[node] | VYEDGE; } }
  
  return;  }

/*   ==========================================
     build the sub-element reference matrices
     ==========================================   */


void construct_sub_element(
     struct All_variables *E
)

{    int i,j,k,l;
     int lev,elx,elz,ely,elzu,elxu,elt,eltu;

     for(lev=E->mesh.levmax-1;lev>=E->mesh.levmin;lev--) { 

        if(E->control.verbose) {
	  fprintf(stderr,"Constructing Sub-element arrays for level %d\n",lev);
	}

	 elx = E->mesh.ELX[lev];
	 elz = E->mesh.ELZ[lev];
	 ely = E->mesh.ELY[lev];
	 elzu = 2 * elz;
	 elxu = 2 * elx;

	  for(i=1;i<=elx;i++)
	    for(j=1;j<=elz;j++)
	      for(k=1;k<=ely;k++)
		{ elt = j + (i-1)*elz +(k-1)*elz*elx;
		  eltu = (j*2-1) + elzu *2*(i-1) + elxu*elzu*2*(k-1);

		  for(l=1;l<=enodes[E->mesh.nsd];l++) {
		    E->EL[lev][elt].sub[l] = eltu
		      + offset[l].vector[1] 
		      + offset[l].vector[0] * elzu
		      + offset[l].vector[2] * elzu * elxu; 
		  }
		}  
       
	  if(E->control.verbose) {
	    fprintf(stderr,"Contructed Sub-element arrays for level %d\n",lev);
	  }
     }
     return; 
}

void construct_elem(
     struct All_variables *E
)
{	
  int lev;
  int nel;
  int i;
  
  const int dims = E->mesh.nsd;

/*--------------------------------------------------*/
/*RAA: 22/3/01, correction here for 3D with name_v, etc 
  N.B. a more elegant solution gave trouble with redefining name_v, p */
  int name_v2=FOUR_NODES_QUAD;
  int name_p2=ONE_NODE_QUAD;
  int name_v3=EIGHT_NODES_CUBIC;
  int name_p3=ONE_NODE_CUBIC;

  if(E->control.verbose) {
    if(2==dims) 
      fprintf(stderr,"***** ELEM type is --------> %d node 2D, hybrid\n",name_v2);
    else
      fprintf(stderr,"***** ELEM type is --------> %d node 3D, hybrid\n",name_v3);
  }
			  
/*--------------------------------------------------*/

  if(2==dims) {
    E->control.ELEMENT_TYPE = FOUR_NODES_QUAD;
    E->control.ELEMENT_TYPE_P = ONE_NODE_QUAD;
  }
  else {
    E->control.ELEMENT_TYPE = EIGHT_NODES_CUBIC;
    E->control.ELEMENT_TYPE_P = ONE_NODE_CUBIC;
  }

  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {

    if(E->control.verbose) {
      fprintf(stderr,"Constructing ELEM arrays for level %d\n",lev);
    }

      nel=E->mesh.NEL[lev];
      for(i=1;i<=nel;i++){
        if(2==dims) {
           E->ELEM[lev][i].type_v = name_v2;
           E->ELEM[lev][i].type_p = name_p2;
	}
	else if(3==dims) {
           E->ELEM[lev][i].type_v = name_v3;
           E->ELEM[lev][i].type_p = name_p3;
	}
      }
  }

return ;
}
