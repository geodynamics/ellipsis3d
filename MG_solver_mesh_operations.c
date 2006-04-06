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

void inject_node_values(
			struct All_variables *E,
			int start_lev,
			standard_precision *AU1,
			standard_precision *AD1  /* data on upper/lower mesh  */
			)
{
  int i,j,k;
  int ii,jj,kk;
  int n,nn;
  standard_precision temp;
  
  for(k=kk=1;k<=E->mesh.NOY[start_lev];k+=2,kk++)
    for(j=jj=1;j<=E->mesh.NOZ[start_lev];j+=2,jj++)
      for(i=ii=1;i<=E->mesh.NOX[start_lev];i+=2,ii++) {
	n=j+(i-1)*E->mesh.NOZ[start_lev]+(k-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	nn=jj+(ii-1)*E->mesh.NOZ[start_lev-1]+(kk-1)*E->mesh.NOZ[start_lev-1]*E->mesh.NOX[start_lev-1];
	AD1[nn] = AU1[n];
      }
  return;
}

void inject_node_int_values(
			    struct All_variables *E,
			    int start_lev,
			    int *AU1,
			    int *AD1  /* data on upper/lower mesh  */
			    )
{
  int i,j,k;
  int ii,jj,kk;
  int n,nn;
  int temp;
  
  for(k=kk=1;k<=E->mesh.NOY[start_lev];k+=2,kk++)
    for(j=jj=1;j<=E->mesh.NOZ[start_lev];j+=2,jj++)
      for(i=ii=1;i<=E->mesh.NOX[start_lev];i+=2,ii++) {
	n=j+(i-1)*E->mesh.NOZ[start_lev]+(k-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	nn=jj+(ii-1)*E->mesh.NOZ[start_lev-1]+(kk-1)*E->mesh.NOZ[start_lev-1]*E->mesh.NOX[start_lev-1];
	AD1[nn] = AU1[n];
      }
  return;
}



void inject06(
     struct All_variables *E,
     int start_lev,
     standard_precision *AU1,
     standard_precision *AU2,
     standard_precision *AU3,
     standard_precision *AU4,
     standard_precision *AU5,
     standard_precision *AU6,

     standard_precision *AD1,
     standard_precision *AD2,
     standard_precision *AD3,
     standard_precision *AD4,
     standard_precision *AD5,
     standard_precision *AD6  /* data on upper/lower mesh  */
)

{
  int i,j,k;
  int ii,jj,kk;
  int n,nn;
  const int dofs=E->mesh.dof ;

  if(start_lev==E->mesh.levmin)
    return;

  if(dofs==2) {
    for(k=kk=1;k<=E->mesh.NOY[start_lev];k+=2,kk++)
      for(j=jj=1;j<=E->mesh.NOZ[start_lev];j+=2,jj++)
	for(i=ii=1;i<=E->mesh.NOX[start_lev];i+=2,ii++) {
	  n=j+(i-1)*E->mesh.NOZ[start_lev]+(k-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	  nn=jj+(ii-1)*E->mesh.NOZ[start_lev-1]+(kk-1)*E->mesh.NOZ[start_lev-1]*E->mesh.NOX[start_lev-1];
	  AD1[nn] = AU1[n];
	  AD2[nn] = AU2[n];
	}
  }
  else if(3==dofs) {
    for(k=kk=1;k<=E->mesh.NOY[start_lev];k+=2,kk++)
      for(j=jj=1;j<=E->mesh.NOZ[start_lev];j+=2,jj++)
	for(i=ii=1;i<=E->mesh.NOX[start_lev];i+=2,ii++) {
	  n=j+(i-1)*E->mesh.NOZ[start_lev]+(k-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	  nn=jj+(ii-1)*E->mesh.NOZ[start_lev-1]+(kk-1)*E->mesh.NOZ[start_lev-1]*E->mesh.NOX[start_lev-1];
	  AD1[nn] = AU1[n];
	  AD2[nn] = AU2[n];
	  AD3[nn] = AU3[n];
	}
  }
  else {
    for(k=kk=1;k<=E->mesh.NOY[start_lev];k+=2,kk++)
      for(j=jj=1;j<=E->mesh.NOZ[start_lev];j+=2,jj++)
	for(i=ii=1;i<=E->mesh.NOX[start_lev];i+=2,ii++) {
	  n=j+(i-1)*E->mesh.NOZ[start_lev]+(k-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	  nn=jj+(ii-1)*E->mesh.NOZ[start_lev-1]+(kk-1)*E->mesh.NOZ[start_lev-1]*E->mesh.NOX[start_lev-1];
	  AD1[nn] = AU1[n];
	  AD2[nn] = AU2[n];
	  AD3[nn] = AU3[n];
	  AD4[nn] = AU4[n];
	  AD5[nn] = AU5[n];
	  AD6[nn] = AU6[n];
	}
  }
  return;
}

#if 0
void project01(
     struct All_variables *E,
     int start_lev,
     standard_precision *AU1,
     standard_precision *AD1  /* data on upper/lower mesh  */
)

{
    int i,j;
    int el,node,node1;
    int eqn1,eqn_minus1;
    int eqn2,eqn_minus2;
    int eqn3,eqn_minus3;
    double average1,average2,average3,w;
    standard_precision CPU_time(),time;
 
    const int sl_minus = start_lev-1;
    const int neq_minus=E->mesh.NEQ[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight= 1.0; /* (double) 1.0/ends; */

    /* on the lower level the average value of data in upper level
       ELEMENTS are slid across to the nodes.  */
 
    if(3==dims) {
      for(i=1;i<=E->mesh.NNO[sl_minus];i++)
	AD1[i] = 0.0;

	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES3D;i++) {
		average1=average2=average3=0.0;
		node1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES3D;j++) {
		    node=E->IEN[start_lev][node1].node[j];
		    average1 += AU1[node];
		}     
	       
		node= E->IEN[sl_minus][el].node[i];
		w=weight * E->TW[sl_minus][node]; 

		AD1[node] += w * average1; 
	    }
    }
   else {

     for(i=1;i<=E->mesh.NNO[sl_minus];i++)
       AD1[i]  = 0.0;

	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES2D;i++)	{ 
		average1=average2=0.0;
		node1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES2D;j++)    {
		    node =E->IEN[start_lev][node1].node[j];
		    average1 += AU1[node];
		}    	
		node=E->IEN[sl_minus][el].node[i];
		w = weight * E->TW[sl_minus][node]; 
		AD1[node] += w * average1; 
	    }
   }
/* Thats's all */

return;  }
#endif
#if 0

void project03(
	       struct All_variables *E,
	       int start_lev,
	       standard_precision *AU1,
	       standard_precision *AU2,
	       standard_precision *AU3,
	       standard_precision *AD1,
	       standard_precision *AD2,
	       standard_precision *AD3  /* data on upper/lower mesh  */
	       )

{
  int i,j;
  int el,node,node1;
  int eqn1,eqn_minus1;
  int eqn2,eqn_minus2;
  int eqn3,eqn_minus3;
  double average1,average2,average3,w;
  standard_precision CPU_time(),time;
 
  const int sl_minus = start_lev-1;
  const int neq_minus=E->mesh.NEQ[start_lev-1];
  const int nels_minus=E->mesh.NEL[start_lev-1];
  const int  dims=E->mesh.nsd;
  const int ends=enodes[E->mesh.nsd];
  const double weight=1.0 /* (double) 1.0/ends */ ;

  /* on the lower level the average value of data in upper level
     ELEMENTS are slid across to the nodes.  */
 
  if(3==dims) {
    for(i=1;i<=E->mesh.NNO[sl_minus];i++)
      AD1[i] = AD2[i] = AD3[i] = 0.0;

    for(el=1;el<=nels_minus;el++)
      for(i=1;i<=ENODES3D;i++) {
	average1=average2=average3=0.0;
	node1 = E->EL[sl_minus][el].sub[i];
	for(j=1;j<=ENODES3D;j++) {
	  node=E->IEN[start_lev][node1].node[j];
	  average1 += AU1[node];
	  average2 += AU2[node];
	  average3 += AU3[node];
	}     
	       
	node= E->IEN[sl_minus][el].node[i];
	w=weight * E->TW[sl_minus][node]; 

	AD1[node] += w * average1; 
	AD2[node] += w * average2; 
	AD3[node] += w * average3; 
      }
  }
  else {

    for(i=1;i<=E->mesh.NNO[sl_minus];i++)
      AD1[i]  = AD2[i] = 0.0;

    for(el=1;el<=nels_minus;el++)
      for(i=1;i<=ENODES2D;i++)	{ 
	average1=average2=0.0;
	node1 = E->EL[sl_minus][el].sub[i];
	for(j=1;j<=ENODES2D;j++)    {
	  node =E->IEN[start_lev][node1].node[j];
	  average1 += AU1[node];
	  average2 += AU2[node];
	}    	
	node=E->IEN[sl_minus][el].node[i];
	w = weight * E->TW[sl_minus][node]; 
	AD1[node] += w * average1; 
	AD2[node] += w * average2; 
      }
  }
  /* Thats's all */

  return;  }

#endif




/* =====================================================================================================
   Projection operation which knows about node arrangement & spacing. This is less `pure' finite element
   approach than using only element information, but it's so much more efficient that I couldn't stop
   myself. My appologies to FE purists everywhere 
   =====================================================================================================*/

#if 0
void project_visc(
     struct All_variables *E,
     int start_lev,
     standard_precision *VVU,
     standard_precision *VD  /* data on upper/lower mesh  */
)
{

  void inject_node_values();
   
    int i,j,k,el;
    standard_precision x1,x2;
    standard_precision n1,n2;
    standard_precision temp;
    int node0,node1,node2;
    int eqn0,eqn1,eqn2;

    const double alpha=0.5;
    const double one_minus_alpha=(1.0-alpha)*0.5;

    const int level = start_lev;
    const int dims =E->mesh.nsd;
    const int ends= enodes[dims];
    
    const int nox = E->mesh.NOX[level];
    const int noz = E->mesh.NOZ[level];
    const int noy = E->mesh.NOY[level];
    const int nno = E->mesh.NNO[level];

    standard_precision *VU;

    VU=(standard_precision *)Malloc0((1+nno)*sizeof(standard_precision));
   
    for(i=1;i<=nno;i++)
	VU[i] = VVU[i];

    /* distribute data from surplus nodes in x direction */
    for(k=1;k<=noy;k++)
	for(j=1;j<=noz;j++)
	    for(i=2;i<nox;i+=2) {
		node0 =  j + (i-1)*noz + (k-1)*nox*noz; /* this node */
		node1 = node0 - noz;
		node2 = node0 + noz;
		VU[node1] = alpha*VU[node1] + one_minus_alpha*VU[node0];
		VU[node2] = alpha*VU[node2] + one_minus_alpha*VU[node0];
	    }

    /* distribute data from surplus nodes in z direction */
    for(k=1;k<=noy;k++)
	for(j=2;j<noz;j+=2)
	    for(i=1;i<=nox;i+=2) {
		node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
		node1 = node0 - 1;
		node2 = node0 + 1;
		VU[node1] = alpha*VU[node1] + one_minus_alpha*VU[node0];
		VU[node2] = alpha*VU[node2] + one_minus_alpha*VU[node0];
	    }
	     
    /* distribute data from surplus nodes in y direction */
    if(3==dims){
	for(k=2;k<noy;k+=2)
	    for(j=1;j<=noz;j+=2)
		for(i=1;i<=nox;i+=2) {
		node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
		node1 = node0 - nox*noz;
		node2 = node0 + nox*noz;
		VU[node1] = alpha*VU[node1] + one_minus_alpha*VU[node0];
		VU[node2] = alpha*VU[node2] + one_minus_alpha*VU[node0];	
		}
    }

    /*
    for(el=1;el<=E->mesh.NEL[start_lev-1];el++)
	for(i=1;i<=ends;i++) {
	    node1 = E->IEN[start_lev-1][el].node[i];
	    node2=E->IEN[start_lev][E->EL[start_lev-1][el].sub[i]].node[i];
	    VD[node1] = VU[node2];
	}
	*/

    inject_node_values(E,start_lev,VU,VD);

	

    free((void *)VU);

    return;
}
#endif
#if 0
void un_inject(
  struct All_variables *E,
  int start_lev,
  standard_precision *AU,
  standard_precision *AD,  /* data on upper/lower mesh  */
  standard_precision scale
)
{
  int i,j,k;
  int ii,jj,kk;
  int n,nn;
  int e,ee;


  if(3==E->mesh.dof) {
    for(k=kk=1;k<=E->mesh.NOY[start_lev+1];k+=2,kk++)
      for(j=jj=1;j<=E->mesh.NOZ[start_lev+1];j+=2,jj++)
	for(i=ii=1;i<=E->mesh.NOX[start_lev+1];i+=2,ii++) {
	  n=j+(i-1)*E->mesh.NOZ[start_lev+1]+(k-1)*E->mesh.NOZ[start_lev+1]*E->mesh.NOX[start_lev+1];
	  nn=jj+(ii-1)*E->mesh.NOZ[start_lev]+(kk-1)*E->mesh.NOZ[start_lev]*E->mesh.NOX[start_lev];
	  e = E->ID[start_lev+1][n].doff[1];
	  ee = E->ID[start_lev][nn].doff[1];
	  AU[e] = AD[ee];
	  e = E->ID[start_lev+1][n].doff[2];
	  ee = E->ID[start_lev][nn].doff[2];
	  AU[e] = AD[ee];
	  e = E->ID[start_lev+1][n].doff[3];
	  ee = E->ID[start_lev][nn].doff[3];
	  AU[e] = AD[ee];
	 
	}
  }
  else {
    for(j=jj=1;j<=E->mesh.NOZ[start_lev+1];j+=2,jj++)
	for(i=ii=1;i<=E->mesh.NOX[start_lev+1];i+=2,ii++) {
	  n=j+(i-1)*E->mesh.NOZ[start_lev+1];
	  nn=jj+(ii-1)*E->mesh.NOZ[start_lev];
	  e = E->ID[start_lev+1][n].doff[1];
	  ee = E->ID[start_lev][nn].doff[1];
	  AU[e] = AD[ee];
	  e = E->ID[start_lev+1][n].doff[2];
	  ee = E->ID[start_lev][nn].doff[2];
	  AU[e] = AD[ee]; 
	}
  }

  return;
}
#endif

/*  ==============================================
    function to project viscosity down to all the 
    levels in the problem. (no gaps for vbcs)
    ==============================================  */

void project_q(
  struct All_variables *E,
  standard_precision *Pfine,
  standard_precision *Pcoarse,
  int level
)
{ 
    int lv,i,j,k,el;
    standard_precision *P1,*P2;
 
    const int vpts = vpoints[E->mesh.nsd];
  
    P1 = (standard_precision *)Malloc0((1+ E->mesh.NNO[level]) * sizeof(standard_precision));
    P2 = (standard_precision *)Malloc0((1+ E->mesh.NNO[level]) * sizeof(standard_precision));
    
    sp_to_nodes(E,Pfine,P1,level);
    inject_node_values(E,level,P1,P2);
    sp_to_centres(E,P2,Pcoarse,level-1); 
  
    free((void *) P1);
    free((void *) P2);

    return;  
}

/*  ==============================================
    function to project viscosity down to all the 
    levels in the problem. (no gaps for vbcs)
    ==============================================  */

void interpolate_q(
			  struct All_variables *E,
			  standard_precision *Pfine,
			  standard_precision *Pcoarse,
			  int level
			  )
{ 
  int i,j,k,elf,elc,ii,jj,kk;
 
  const int dims = E->mesh.nsd;
   
    /* D1 - D4 are dummy variables but should be set to some value to prevent
       unitialized variable type errors */

   for(k=kk=1;k<=E->mesh.ELY[level];kk++,k+=2)
     for(i=ii=1;i<=E->mesh.ELX[level];ii++,i+=2)
       for(j=jj=1;j<=E->mesh.ELZ[level];jj++,j+=2) {
	 elf=j+(i-1)*E->mesh.ELZ[level]+(k-1)*E->mesh.ELZ[level]*E->mesh.ELX[level];
	 elc=jj+(ii-1)*E->mesh.ELZ[level-1]+(kk-1)*E->mesh.ELZ[level-1]*E->mesh.ELX[level-1];
	 
	 Pfine[elf] = Pfine[elf+1] = 
	   Pfine[elf+E->mesh.ELZ[level]] = Pfine[elf+E->mesh.ELZ[level]+1] = Pcoarse[elc];
 
	 if(3==dims)
	   Pfine[elf+E->mesh.ELZ[level]*E->mesh.ELX[level]] =
	     Pfine[elf+E->mesh.ELZ[level]*E->mesh.ELX[level]+1] = 
	     Pfine[elf+E->mesh.ELZ[level]+E->mesh.ELZ[level]*E->mesh.ELX[level]] = 
	     Pfine[elf+E->mesh.ELZ[level]+E->mesh.ELZ[level]*E->mesh.ELX[level]+1] = Pcoarse[elc];
       }

   /*
   sp_to_nodes(E,Pcoarse,P1,level-1);
   interpolation_4pt(E,P2,D1,D2,P1,D3,D4,level);
   sp_to_centres(E,P2,P1,level);
   
   for(i=1;i<=E->mesh.NPNO[level];i++)
     Pfine[i] = 0.5 * Pfine[i] + 0.5 * P1[i];
     */
    return;  
}

/* Given an ordered set of points XX and values there (YY),
   determine the interpolation at point x using four local values */
#if 0
standard_precision interpolation_1d_4pt(
  struct All_variables *E,
  standard_precision *XX,
  standard_precision *YY,
  int points,
  standard_precision x
)
{
  int locn;
  int n,node,node0,node1,node2,node3;
  standard_precision x0,x1,x2,x3;
  standard_precision xx0,xx1,xx2,xx3;
  standard_precision interpolant;

  /* First search for the nearest neighbour nodes */

  locn=0;
  for(n=1;n<=points;n++)
    if (XX[n] > x ) {
      locn = n;
      break;
    }
  
  if(locn==1 || locn==2) {  /* It's a point before all the others or right up at the beginning */
    node0=1;
    node1=2;
    node2=3;
    node3=4;
  }
  else if(locn==points || locn==0) {  /* Last point or beyond it */
      node0=points;
      node1=points-1;
      node2=points-2;
      node3=points-3;
    }

  else {
    node0=locn-2;
    node1=locn-1;
    node2=locn;
    node3=locn+1;
  }

   x0 = XX[node0];
   x1 = XX[node1];
   x2 = XX[node2];
   x3 = XX[node3];

   /* if any are equal there will be trouble ! */

   xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
   xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
   xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
   xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

   interpolant = xx0 * YY[node0] +  xx1 * YY[node1] +  xx2 * YY[node2] +  xx3 * YY[node3];

   if(isnan(interpolant))
    fprintf(stderr,"locn %d/%d P(%g) = %g, based on P(%g) = %g, P(%g) = %g, P(%g) = %g, P(%g) = %g [%d/%d/%d/%d]\n",
	   locn,points,x,interpolant,x0,YY[node0],x1,YY[node1],x2,YY[node2],x3,YY[node3],
	    node0,node1,node2,node3); 

   return(interpolant);
}

#endif


/* ========================================================================
   Projection operation
   ======================================================================== */

void projection_integral(
  struct All_variables *E,
  standard_precision *V1,
  standard_precision *V2,
  standard_precision *V3,
  standard_precision *V4,
  standard_precision *V5,
  standard_precision *V6,
  standard_precision *U1,
  standard_precision *U2,
  standard_precision *U3,
  standard_precision *U4,
  standard_precision *U5,
  standard_precision *U6,
  int level
)
{
  standard_precision *ev1,*ev2,*ev3,*ev4,*ev5,*ev6,*ea;
  standard_precision *work1,*work2,*work3,*work4,*work5,*work6;
  standard_precision area;
  
  higher_precision *node_R;
  
  int e,n,nn,i,j,k,ii,jj,kk;

  struct Shape_function GN;
  struct Shape_function_dx GNx;
  struct Shape_function_dA dOmega;

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];

 
  ev1 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ev2 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ev3 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ev4 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ev5 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ev6 = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  ea  = (standard_precision *)Malloc0((E->mesh.NEL[level] +1) * sizeof(standard_precision));
  work1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work2 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work3 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work4 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work5 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work6 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));


  /* 1 get element average values */
  if(2==dofs && 2==dims) { /* 2D classical */

    if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,V1,level);
       flogical_mesh_to_real(E,V2,level);
     }

    for(i=1;i<=E->mesh.NNO[level];i++) {
      work1[i] = V1[i];
      work2[i] = V2[i];
    }

   if(E->control.HAVE_SKEWBCS)
     for(i=1;i<=E->mesh.NNO[level];i++) {  /* rotate to cartesian for interpolation */
       if(E->NODE[level][i] & SKEWBC) {
	 node_R = E->curvilinear.NODE_R[level][i];
	 work1[i] = node_R[0*dims+0] * V1[i] + node_R[0*dims+1] * V2[i] ;
	 work2[i] = node_R[1*dims+0] * V1[i] + node_R[1*dims+1] * V2[i] ;
       }
     }

    for(e=1;e<=E->mesh.NEL[level];e++) {
      /* get_global_shape_fn(E,e,&GN,&GNx,&dOmega,0,level); */
      area = E->ECO[level][e].area;
      ev1[e] = ev2[e] = ea[e] = 0.0;
      ea[e] = area;
      for(j=1;j<=ENODES2D;j++) {
	i = E->IEN[level][e].node[j];
	for(k=1;k<=VPOINTS2D;k++) {
	  ev1[e] += work1[i] * E->N.vpt[GNVINDEX(j,k)] * area /* dOmega.vpt[k]*/;
	  ev2[e] += work2[i] * E->N.vpt[GNVINDEX(j,k)] * area /* dOmega.vpt[k]*/;
	  /*ea[e]  +=         E->N.vpt[GNVINDEX(j,k)] * dOmega.vpt[k];*/
	}
      }
    }
    
    for(k=kk=1;k<=noy;kk++,k+=2)
      for(i=ii=1;i<=nox;ii++,i+=2) 
	for(j=jj=1;j<=noz;jj++,j+=2) {
	  n= j + (i-1) * noz + (k-1) * noz * nox;
	  nn = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1] ;
	
	  U1[nn] = U2[nn] = 0.0;

	  if(E->NODE[level-1][nn] & ( OFFSIDE ))
	    continue;

	  area = 0.0;
	  for(e=1;e<=E->NEI[level].nels[n];e++) {
	    U1[nn] += ev1[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	    U2[nn] += ev2[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	    area +=   ea[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	  }
	  U1[nn] /= area;
	  U2[nn] /= area; 
      }
  
    for(i=1;i<=E->mesh.NNO[level-1];i++) {
      if(E->NODE[level-1][i] & ( OFFSIDE ))
	continue;
 
      if(E->control.HAVE_SKEWBCS) /* and rotate back to skewed coord system */ {
	if(E->NODE[level-1][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level-1][i];
	  work1[i] = node_R[0*dims+0] * U1[i] + node_R[1*dims+0] * U2[i] ;
	  work2[i] = node_R[0*dims+1] * U1[i] + node_R[1*dims+1] * U2[i] ;
	  U1[i] = work1[i] ;
	  U2[i] = work2[i] ;
	}
      }

      if((E->NODE[level-1][i] & BC1) != 0)
	U1[i] = 0.0;
      if((E->NODE[level-1][i] & BC2) != 0)
	U2[i] = 0.0;
    }
  }

  else if (dims==3 && dofs==3) {  /* 3D CLASSICAL */

    /*RAA: 1/6/01, this periodic part missing until now*/
    if(E->mesh.periodic_x || E->mesh.periodic_y) {
       flogical_mesh_to_real(E,V1,level);
       flogical_mesh_to_real(E,V2,level);
       flogical_mesh_to_real(E,V3,level);
    }

	  
    for(i=1;i<=E->mesh.NNO[level];i++) {
      work1[i] = V1[i];
      work2[i] = V2[i];
      work3[i] = V3[i];
    }


    if(E->control.HAVE_SKEWBCS) {
      for(i=1;i<=E->mesh.NNO[level];i++) {  /* rotate to cartesian for interpolation */
	if(E->NODE[level][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][i];
	  work1[i] = node_R[0*dims+0] * V1[i] + node_R[0*dims+1] * V2[i] + node_R[0*dims+2] * V3[i];
	  work2[i] = node_R[1*dims+0] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[1*dims+2] * V3[i];
	  work3[i] = node_R[2*dims+0] * V1[i] + node_R[2*dims+1] * V2[i] + node_R[2*dims+2] * V3[i];
	}
      }
    }


    for(e=1;e<=E->mesh.NEL[level];e++) {
      /*get_global_shape_fn(E,e,&GN,&GNx,&dOmega,0,level); */
      area = E->ECO[level][e].area;
      ev1[e] = ev2[e] = ev3[e] = ea[e] = 0.0;
      ea[e] = area;
	for(j=1;j<=ENODES3D;j++) {
	  i = E->IEN[level][e].node[j];
	  for(k=1;k<=VPOINTS3D;k++) {
	    ev1[e] += work1[i] * E->N.vpt[GNVINDEX(j,k)] * area /* dOmega.vpt[k] */;
	    ev2[e] += work2[i] * E->N.vpt[GNVINDEX(j,k)] * area /* dOmega.vpt[k] */;
	    ev3[e] += work3[i] * E->N.vpt[GNVINDEX(j,k)] * area /* dOmega.vpt[k] */;
	    /*ea[e]  +=         E->N.vpt[GNVINDEX(j,k)] * dOmega.vpt[k];*/
	  }
	}
    }
    
    for(k=kk=1;k<=noy;kk++,k+=2)
      for(i=ii=1;i<=nox;ii++,i+=2) 
	for(j=jj=1;j<=noz;jj++,j+=2) {
	  n= j + (i-1) * noz + (k-1) * noz * nox;
	  nn = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1] ;
	  
	  U1[nn] = U2[nn] = U3[nn] = 0.0;

	  if(E->NODE[level-1][nn] & ( OFFSIDE ))
	    continue;

	  area = 0.0;
	  for(e=1;e<=E->NEI[level].nels[n];e++) {
	    U1[nn] += ev1[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U2[nn] += ev2[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U3[nn] += ev3[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    area +=   ea[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	  }
	  U1[nn] /= area;
	  U2[nn] /= area;
  	  U3[nn] /= area;
      }
  
    for(i=1;i<=E->mesh.NNO[level-1];i++) {
      if(E->NODE[level-1][i] & ( OFFSIDE ))
	continue;

      if(E->control.HAVE_SKEWBCS) /* and rotate back to skewed coord system */ {
	if(E->NODE[level-1][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level-1][i];
	  work1[i] = node_R[0*dims+0] * U1[i] + node_R[1*dims+0] * U2[i] + node_R[2*dims+0] * U3[i];
	  work2[i] = node_R[0*dims+1] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[2*dims+1] * U3[i];
	  work3[i] = node_R[0*dims+2] * U1[i] + node_R[1*dims+2] * U2[i] + node_R[2*dims+2] * U3[i];
	  U1[i] = work1[i] ;
	  U2[i] = work2[i] ;
	  U3[i] = work3[i] ;
	}
      }
    
      if((E->NODE[level-1][i] & BC1) != 0)
	U1[i] = 0.0;
      if((E->NODE[level-1][i] & BC2) != 0)
	U2[i] = 0.0;
      if((E->NODE[level-1][i] & BC3) != 0)
	U3[i] = 0.0;
    }
  }
  else if (dims==2 && dofs==3) {  /* 2D COSSERAT */

    for(i=1;i<=E->mesh.NNO[level];i++) {
      work1[i] = V1[i];
      work2[i] = V2[i];
    }

    if(E->control.HAVE_SKEWBCS) {
      for(i=1;i<=E->mesh.NNO[level];i++) {  /* rotate to cartesian for interpolation */
	if(E->NODE[level][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][i];
	  work1[i] = node_R[0*dims+0] * V1[i] + node_R[0*dims+1] * V2[i] ;
	  work2[i] = node_R[1*dims+0] * V1[i] + node_R[1*dims+1] * V2[i] ;
	}
      }
    }

    for(e=1;e<=E->mesh.NEL[level];e++) {
      area = E->ECO[level][e].area;
      ev1[e] = ev2[e] = ev3[e] = 0.0;
      ea[e] = area;
      for(j=1;j<=ENODES2D;j++) {
	i = E->IEN[level][e].node[j];
	for(k=1;k<=VPOINTS2D;k++) {
	  ev1[e] += work1[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	  ev2[e] += work2[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	  ev3[e] += V3[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	}
      }
    }
    
    for(k=kk=1;k<=noy;kk++,k+=2)
      for(i=ii=1;i<=nox;ii++,i+=2) 
	for(j=jj=1;j<=noz;jj++,j+=2) {
	  n= j + (i-1) * noz + (k-1) * noz * nox;
	  nn = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1] ;
	  
	  U1[nn] = U2[nn] = U3[nn] = 0.0;

	  if(E->NODE[level-1][nn] & ( OFFSIDE ))
	    continue;
	  
	  area = 0.0;
	  for(e=1;e<=E->NEI[level].nels[n];e++) {
	    U1[nn] += ev1[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	    U2[nn] += ev2[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	    U3[nn] += ev3[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	    area +=   ea[E->NEI[level].element[(n-1)*ENODES2D+e-1]];
	  }
	  U1[nn] /= area;
	  U2[nn] /= area;
	  U3[nn] /= area;
	}
    
    for(i=1;i<=E->mesh.NNO[level-1];i++) {
      if(E->NODE[level-1][i] & ( OFFSIDE ))
	continue;
      
      if(E->control.HAVE_SKEWBCS) /* and rotate back to skewed coord system */ {
	if(E->NODE[level-1][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level-1][i];
	  work1[i] = node_R[0*dims+0] * U1[i] + node_R[1*dims+0] * U2[i] ;
	  work2[i] = node_R[0*dims+1] * U1[i] + node_R[1*dims+1] * U2[i] ;
	  U1[i] = work1[i] ;
	  U2[i] = work2[i] ;
	}
      }
    
      if((E->NODE[level-1][i] & BC1) != 0)
	U1[i] = 0.0;
      if((E->NODE[level-1][i] & BC2) != 0)
	U2[i] = 0.0;
      if((E->NODE[level-1][i] & BC3) != 0)
	U3[i] = 0.0;
    }
  }
    else if (dims==3 && dofs==6) {  /* 3D COSSERAT */

    for(i=1;i<=E->mesh.NNO[level];i++) {
      work1[i] = V1[i];
      work2[i] = V2[i];
      work3[i] = V3[i];
      work4[i] = V4[i];
      work5[i] = V5[i];
      work6[i] = V6[i];
    }

    if(E->control.HAVE_SKEWBCS) {
      for(i=1;i<=E->mesh.NNO[level];i++) {  /* rotate to cartesian for interpolation */
	if(E->NODE[level][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][i];
	  work1[i] = node_R[0*dims+0] * V1[i] + node_R[0*dims+1] * V2[i] + node_R[0*dims+2] * V3[i];
	  work2[i] = node_R[1*dims+0] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[1*dims+2] * V3[i];
	  work3[i] = node_R[2*dims+0] * V1[i] + node_R[2*dims+1] * V2[i] + node_R[2*dims+2] * V3[i];
	  work4[i] = node_R[0*dims+0] * V4[i] + node_R[0*dims+1] * V5[i] + node_R[0*dims+2] * V6[i];
	  work5[i] = node_R[1*dims+0] * V4[i] + node_R[1*dims+1] * V5[i] + node_R[1*dims+2] * V6[i];
	  work6[i] = node_R[2*dims+0] * V4[i] + node_R[2*dims+1] * V5[i] + node_R[2*dims+2] * V6[i];
	}
      }
    }

    for(e=1;e<=E->mesh.NEL[level];e++) {
      area = E->ECO[level][e].area;
      ev1[e] = ev2[e] = ev3[e] = ev4[e] = ev5[e] = ev6[e] = 0.0;
      ea[e] = area;
	for(j=1;j<=ENODES3D;j++) {
	  i = E->IEN[level][e].node[j];
	  for(k=1;k<=VPOINTS3D;k++) {
	    ev1[e] += work1[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	    ev2[e] += work2[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	    ev3[e] += work3[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	    ev4[e] += work4[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	    ev5[e] += work5[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	    ev6[e] += work6[i] * E->N.vpt[GNVINDEX(j,k)] * area ;
	  }
	}
    }
    
    for(k=kk=1;k<=noy;kk++,k+=2)
      for(i=ii=1;i<=nox;ii++,i+=2) 
	for(j=jj=1;j<=noz;jj++,j+=2) {
	  n= j + (i-1) * noz + (k-1) * noz * nox;
	  nn = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1] ;

	  U1[nn] = U2[nn] = U3[nn] = U4[nn] = U5[nn] = U6[nn] = 0.0;
	  
	  if(E->NODE[level-1][nn] & ( OFFSIDE ))
	    continue;

	  area = 0.0;
	  for(e=1;e<=E->NEI[level].nels[n];e++) {
	    U1[nn] += ev1[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U2[nn] += ev2[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U3[nn] += ev3[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U4[nn] += ev4[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U5[nn] += ev5[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    U6[nn] += ev6[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	    area +=   ea[E->NEI[level].element[(n-1)*ENODES3D+e-1]];
	  }
	  U1[nn] /= area;
	  U2[nn] /= area;
  	  U3[nn] /= area;
	  U4[nn] /= area;
	  U5[nn] /= area;
  	  U6[nn] /= area;
      }
  
    for(i=1;i<=E->mesh.NNO[level-1];i++) {
      if(E->NODE[level-1][i] & ( OFFSIDE ))
	continue;

      if(E->control.HAVE_SKEWBCS) /* and rotate back to skewed coord system */ {
	if(E->NODE[level-1][i] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level-1][i];
	  work1[i] = node_R[0*dims+0] * U1[i] + node_R[1*dims+0] * U2[i] + node_R[2*dims+0] * U3[i];
	  work2[i] = node_R[0*dims+1] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[2*dims+1] * U3[i];
	  work3[i] = node_R[0*dims+2] * U1[i] + node_R[1*dims+2] * U2[i] + node_R[2*dims+2] * U3[i];
	  work4[i] = node_R[0*dims+0] * U4[i] + node_R[1*dims+0] * U5[i] + node_R[2*dims+0] * U6[i];
	  work5[i] = node_R[0*dims+1] * U4[i] + node_R[1*dims+1] * U5[i] + node_R[2*dims+1] * U6[i];
	  work6[i] = node_R[0*dims+2] * U4[i] + node_R[1*dims+2] * U5[i] + node_R[2*dims+2] * U6[i];
	  U1[i] = work1[i] ;
	  U2[i] = work2[i] ;
	  U3[i] = work3[i] ;
	  U4[i] = work4[i] ;
	  U5[i] = work5[i] ;
	  U6[i] = work6[i] ;
	}
      }
    
      if((E->NODE[level-1][i] & BC1) != 0)
	U1[i] = 0.0;
      if((E->NODE[level-1][i] & BC2) != 0)
	U2[i] = 0.0;
      if((E->NODE[level-1][i] & BC3) != 0)
	U3[i] = 0.0;
      if((E->NODE[level-1][i] & BC4) != 0)
	U4[i] = 0.0;
      if((E->NODE[level-1][i] & BC5) != 0)
	U5[i] = 0.0;
      if((E->NODE[level-1][i] & BC6) != 0)
	U6[i] = 0.0;
    }
  }
  
  free((void *) work1);
  free((void *) work2);
  free((void *) work3);
  free((void *) work4);
  free((void *) work5);
  free((void *) work6);
  free((void *) ev1);
  free((void *) ev2);
  free((void *) ev3);
  free((void *) ev4);
  free((void *) ev5);
  free((void *) ev6);
  free((void *) ea);
}


/* ========================================================================
   Interpolation operation
   ======================================================================== */

void interpolation_6(
  struct All_variables *E,
  standard_precision *V1,
  standard_precision *V2,
  standard_precision *V3,
  standard_precision *V4,
  standard_precision *V5,
  standard_precision *V6,
  standard_precision *U1,
  standard_precision *U2,
  standard_precision *U3,
  standard_precision *U4,
  standard_precision *U5,
  standard_precision *U6,
  int level
)
{
  standard_precision *work1,*work2,*work3,*work4,*work5,*work6;
  standard_precision x1,x2,x;
  standard_precision xx,xx1,xx2;
  int i,j,k;
  int ii,jj,kk;
  int n,node,node1,node2;

  higher_precision *node_R;

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];

  work1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work2 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work3 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work4 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work5 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work6 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));

  /* The coarse grid is uninjected to the fine and then the result is interpolated */

   if(2==dims && 2==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
     }

     /* if cylindrical, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;	  
	 }
       }

     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {
	 node = j + (i-1) * noz;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1];
	
	 work1[node] = V1[node] = U1[node1];
	 work2[node] = V2[node] = U2[node1];
	 
       }
 
     for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2) {
 	 node=j + (i-1) * noz;
	 node1 = node-noz;
	 node2 = node+noz;
	 x = E->SX[level][1][node]; 
	 x1 = E->SX[level][1][node1];
	 x2 = E->SX[level][1][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
 
     for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
	 node=j + (i-1) * noz;
	 node1 = node-1;
	 node2 = node+1;
	 x = E->SX[level][2][node]; 
	 x1 = E->SX[level][2][node1];
	 x2 = E->SX[level][2][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
    
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i];
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i];
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;
	
	 }
       }
   }

   if(2==dims && 3==dofs) { /* 2D COSSERAT */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     /* if cylindrical, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;	  
	 }
       }

     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {
	 node = j + (i-1) * noz;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1];
	 work1[node] = V1[node] = U1[node1];
	 work2[node] = V2[node] = U2[node1];
	 V3[node] = U3[node1];
	 
       }
 
     for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2) {
 	 node=j + (i-1) * noz;
	 node1 = node-noz;
	 node2 = node+noz;
	 x = E->SX[level][1][node]; 
	 x1 = E->SX[level][1][node1];
	 x2 = E->SX[level][1][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V3[node] = (V3[node1] * xx2 + V3[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
 
     for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
	 node=j + (i-1) * noz;
	 node1 = node-1;
	 node2 = node+1;
	 x = E->SX[level][2][node]; 
	 x1 = E->SX[level][2][node1];
	 x2 = E->SX[level][2][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V3[node] = (V3[node1] * xx2 + V3[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
    
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i];
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i];
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;	
	 }
       }
   }
   
   else if(3==dims && 3==dofs) { /* 3D CLASSICAL */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	 }
       }
     
     for(k=kk=1; k<=noy; kk++,k+=2)
       for(j=jj=1; j<=noz; jj++,j+=2)
	 for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	}

     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   x = E->SX[level][1][node]; 
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	 }

     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-1;
	   node2 = node+1;
	   x = E->SX[level][2][node]; 
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	 }

   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  x = E->SX[level][3][node]; 
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  xx = 1.0/(x2-x1);
	  xx1= x-x1;
	  xx2= x2-x;
	  work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	  work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	  work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	}

      /* if sphere, rotate boundary nodes back to rtf */
     
   if(E->control.HAVE_SKEWBCS) 
     for(i=1;i<=E->mesh.NNO[level];i++) {
       if(E->NODE[level][i] & SKEWBC) {
	 node_R = E->curvilinear.NODE_R[level][i];
	 work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	 work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	 work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	 V1[i] = work1[i] ;
	 V2[i] = work2[i] ;
	 V3[i] = work3[i] ;
       }
     }
   }
    else if(3==dims && 6==dofs) { /* 3D COSSERAT */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
       flogical_mesh_to_real(E,U4,level-1);
       flogical_mesh_to_real(E,U5,level-1);
       flogical_mesh_to_real(E,U6,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   work1[i] = node_R[0*dims+0] * U4[i] + node_R[0*dims+1] * U5[i] + node_R[0*dims+2] * U6[i];
	   work2[i] = node_R[1*dims+0] * U4[i] + node_R[1*dims+1] * U5[i] + node_R[1*dims+2] * U6[i];
	   work3[i] = node_R[2*dims+0] * U4[i] + node_R[2*dims+1] * U5[i] + node_R[2*dims+2] * U6[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	   U4[i] = work4[i] ;
	   U5[i] = work5[i] ;
	   U6[i] = work6[i] ;
	 }
       }
     
     for(k=kk=1; k<=noy; kk++,k+=2)
       for(j=jj=1; j<=noz; jj++,j+=2)
	 for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	  work4[node] = V4[node] = U4[node1];
	  work5[node] = V5[node] = U5[node1];
	  work6[node] = V6[node] = U6[node1];
	}

     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   x = E->SX[level][1][node]; 
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	   work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	   work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
	   V6[node] = work6[node];
	 }

     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-1;
	   node2 = node+1;
	   x = E->SX[level][2][node]; 
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	   work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	   work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
	   V6[node] = work6[node];
	 }

   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  x = E->SX[level][3][node]; 
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  xx = 1.0/(x2-x1);
	  xx1= x-x1;
	  xx2= x2-x;
	  work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	  work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	  work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	  work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	  work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	  work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	  V4[node] = work4[node];
	  V5[node] = work5[node];
	  V6[node] = work6[node];
	}

      /* if sphere, rotate boundary nodes back to rtf */
     
   if(E->control.HAVE_SKEWBCS) 
     for(i=1;i<=E->mesh.NNO[level];i++) {
       if(E->NODE[level][i] & SKEWBC) {
	 node_R = E->curvilinear.NODE_R[level][i];
	 work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	 work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	 work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	 work4[i] = node_R[0*dims+0] * V4[i] + node_R[1*dims+0] * V5[i] + node_R[2*dims+0] * V6[i];
	 work5[i] = node_R[0*dims+1] * V4[i] + node_R[1*dims+1] * V5[i] + node_R[2*dims+1] * V6[i];
	 work6[i] = node_R[0*dims+2] * V4[i] + node_R[1*dims+2] * V5[i] + node_R[2*dims+2] * V6[i];
	 V1[i] = work1[i] ;
	 V2[i] = work2[i] ;
	 V3[i] = work3[i] ;
	 V4[i] = work4[i] ;
	 V5[i] = work5[i] ;
	 V6[i] = work6[i] ;
       }
     }
   }
   free((void *) work1);
   free((void *) work2);
   free((void *) work3);
   free((void *) work4);
   free((void *) work5);
   free((void *) work6);
}

#if 0
/* ========================================================================
   Interpolation operation
   ======================================================================== */

void test_interpolation_6(
  struct All_variables *E,
  standard_precision *V1,
  standard_precision *V2,
  standard_precision *V3,
  standard_precision *V4,
  standard_precision *V5,
  standard_precision *V6,
  standard_precision *U1,
  standard_precision *U2,
  standard_precision *U3,
  standard_precision *U4,
  standard_precision *U5,
  standard_precision *U6,
  int level
)

{
  standard_precision *work1,*work2,*work3,*work4,*work5,*work6;
  standard_precision x1,x2,x;
  standard_precision xx,xx1,xx2;
  int i,j,k;
  int ii,jj,kk;
  int n,node,node1,node2;

  higher_precision *node_R;

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];

  work1 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work2 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work3 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work4 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work5 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));
  work6 = (standard_precision *) Malloc0((E->mesh.NNO[level]+1) * sizeof(standard_precision));

  /* The coarse grid is uninjected to the fine and then the result is interpolated */

   if(2==dims && 2==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
     }


     for(i=1;i<=E->mesh.NNO[level];i++) {
       work1[i] = work2[i] = V1[i] = V2[i] = 0.0;

     }

     /* if cylindrical, rotate boundary nodes to xyz 
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;	  
	 }
       }
     */
       

     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {  /* This could be redone in a non-mesh-dependent way */
	 node = j + (i-1) * noz;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1];
	
	 work1[node] = 1.0;
	 V1[node] = U1[node1]; /* First test this part ! */
	 work2[node] = 1.0;
	 V2[node] = U2[node1];
	 
       }

     /* Use node-k relationships to distribute information to surrounding nodes */

     
     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {  /* This could be redone in a non-mesh-dependent way */
	 node = j + (i-1) * noz;
	 
	 ii = node * MAX_NODE_INT_2D;

	 for(k=0;k<MAX_NODE_INT_2D;k++) {
	   node1 = E->Node_map_3[level][ii+k];
	   /* fprintf(stderr,"Interpolating from node %d to %d\n",node,node1); */
	   if(node != node1 && node1 != 0) {
	     V1[node1] += V1[node]  * E->Node_k11[level][ii+k] + V2[node] * E->Node_k21[level][ii+k];
	     V2[node1] += V1[node]  * E->Node_k12[level][ii+k] + V2[node] * E->Node_k22[level][ii+k];
	     work1[node1] +=  0.0 + 1.0 * (E->Node_k11[level][ii+k] + E->Node_k21[level][ii+k]);
	     work2[node1] +=  0.0 + 1.0 * (E->Node_k12[level][ii+k] + E->Node_k22[level][ii+k]);
	   }
	 }
       }

     for(i=1;i<=E->mesh.NNO[level];i++) {
       if(work1[i] != 0.0)
	 V1[i]  /= work1[i] ;
       if(work2[i] != 0.0)
	 V2[i]  /= work2[i] ;
       
       /*  fprintf(stderr,"Node %d - V1 = %g,%g  weight %g,%g \n",i,V1[i],V2[i],work1[i],work2[i]); */
     }


   }

   if(2==dims && 3==dofs) { /* 2D COSSERAT */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     /* if cylindrical, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;	  
	 }
       }

     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {
	 node = j + (i-1) * noz;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1];
	 work1[node] = V1[node] = U1[node1];
	 work2[node] = V2[node] = U2[node1];
	 V3[node] = U3[node1];
	 
       }
 
     for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2) {
 	 node=j + (i-1) * noz;
	 node1 = node-noz;
	 node2 = node+noz;
	 x = E->SX[level][1][node]; 
	 x1 = E->SX[level][1][node1];
	 x2 = E->SX[level][1][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V3[node] = (V3[node1] * xx2 + V3[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
 
     for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
	 node=j + (i-1) * noz;
	 node1 = node-1;
	 node2 = node+1;
	 x = E->SX[level][2][node]; 
	 x1 = E->SX[level][2][node1];
	 x2 = E->SX[level][2][node2];
	 xx = 1.0/(x2-x1);
	 xx1= x-x1;
	 xx2= x2-x;
	 work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	 work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	 V3[node] = (V3[node1] * xx2 + V3[node2] * xx1) * xx;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
    
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i];
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i];
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;	
	 }
       }
   }
   
   else if(3==dims && 3==dofs) { /* 3D CLASSICAL */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	 }
       }
     
     for(k=kk=1; k<=noy; kk++,k+=2)
       for(j=jj=1; j<=noz; jj++,j+=2)
	 for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	}

     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   x = E->SX[level][1][node]; 
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	 }

     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-1;
	   node2 = node+1;
	   x = E->SX[level][2][node]; 
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	 }

   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  x = E->SX[level][3][node]; 
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  xx = 1.0/(x2-x1);
	  xx1= x-x1;
	  xx2= x2-x;
	  work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	  work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	  work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	}

      /* if sphere, rotate boundary nodes back to rtf */
     
   if(E->control.HAVE_SKEWBCS) 
     for(i=1;i<=E->mesh.NNO[level];i++) {
       if(E->NODE[level][i] & SKEWBC) {
	 node_R = E->curvilinear.NODE_R[level][i];
	 work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	 work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	 work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	 V1[i] = work1[i] ;
	 V2[i] = work2[i] ;
	 V3[i] = work3[i] ;
       }
     }
   }
    else if(3==dims && 6==dofs) { /* 3D COSSERAT */
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
       flogical_mesh_to_real(E,U4,level-1);
       flogical_mesh_to_real(E,U5,level-1);
       flogical_mesh_to_real(E,U6,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   work1[i] = node_R[0*dims+0] * U4[i] + node_R[0*dims+1] * U5[i] + node_R[0*dims+2] * U6[i];
	   work2[i] = node_R[1*dims+0] * U4[i] + node_R[1*dims+1] * U5[i] + node_R[1*dims+2] * U6[i];
	   work3[i] = node_R[2*dims+0] * U4[i] + node_R[2*dims+1] * U5[i] + node_R[2*dims+2] * U6[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	   U4[i] = work4[i] ;
	   U5[i] = work5[i] ;
	   U6[i] = work6[i] ;
	 }
       }
     
     for(k=kk=1; k<=noy; kk++,k+=2)
       for(j=jj=1; j<=noz; jj++,j+=2)
	 for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	  work4[node] = V4[node] = U4[node1];
	  work5[node] = V5[node] = U5[node1];
	  work6[node] = V6[node] = U6[node1];
	}

     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   x = E->SX[level][1][node]; 
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	   work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	   work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
	   V6[node] = work6[node];
	 }

     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-1;
	   node2 = node+1;
	   x = E->SX[level][2][node]; 
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   xx = 1.0/(x2-x1);
	   xx1= x-x1;
	   xx2= x2-x;
	   work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	   work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	   work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	   work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	   work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	   work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
	   V6[node] = work6[node];
	 }

   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  x = E->SX[level][3][node]; 
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  xx = 1.0/(x2-x1);
	  xx1= x-x1;
	  xx2= x2-x;
	  work1[node] = (work1[node1] * xx2 + work1[node2] * xx1) * xx;
	  work2[node] = (work2[node1] * xx2 + work2[node2] * xx1) * xx;
	  work3[node] = (work3[node1] * xx2 + work3[node2] * xx1) * xx;
	  work4[node] = (work4[node1] * xx2 + work4[node2] * xx1) * xx;
	  work5[node] = (work5[node1] * xx2 + work5[node2] * xx1) * xx;
	  work6[node] = (work6[node1] * xx2 + work6[node2] * xx1) * xx;
	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	  V4[node] = work4[node];
	  V5[node] = work5[node];
	  V6[node] = work6[node];
	}

      /* if sphere, rotate boundary nodes back to rtf */
     
   if(E->control.HAVE_SKEWBCS) 
     for(i=1;i<=E->mesh.NNO[level];i++) {
       if(E->NODE[level][i] & SKEWBC) {
	 node_R = E->curvilinear.NODE_R[level][i];
	 work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	 work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	 work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	 work4[i] = node_R[0*dims+0] * V4[i] + node_R[1*dims+0] * V5[i] + node_R[2*dims+0] * V6[i];
	 work5[i] = node_R[0*dims+1] * V4[i] + node_R[1*dims+1] * V5[i] + node_R[2*dims+1] * V6[i];
	 work6[i] = node_R[0*dims+2] * V4[i] + node_R[1*dims+2] * V5[i] + node_R[2*dims+2] * V6[i];
	 V1[i] = work1[i] ;
	 V2[i] = work2[i] ;
	 V3[i] = work3[i] ;
	 V4[i] = work4[i] ;
	 V5[i] = work5[i] ;
	 V6[i] = work6[i] ;
       }
     }
   }
   free((void *) work1);
   free((void *) work2);
   free((void *) work3);
   free((void *) work4);
   free((void *) work5);
   free((void *) work6);
}
#endif

/* ========================================================================
   Interpolation operation  from level-1 -> level  (V fine, U coarse)
   ======================================================================== */

void interpolation_4pt_6(
  struct All_variables *E,
  standard_precision *V1, /* output */
  standard_precision *V2,
  standard_precision *V3,
  standard_precision *V4,
  standard_precision *V5,
  standard_precision *V6,
  standard_precision *U1, /* input */
  standard_precision *U2,
  standard_precision *U3,
  standard_precision *U4,
  standard_precision *U5,
  standard_precision *U6,
  int level
)
{
  standard_precision *work1,*work2,*work3,*work4,*work5,*work6;
  standard_precision x0,x1,x2,x3,x;
  standard_precision xx,xx0,xx1,xx2,xx3;
  int i,j,k;
  int ii,jj,kk;
  int n,node,node0,node1,node2,node3;

  higher_precision *node_R;

  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];
  const int nno = E->mesh.NNO[level] ;

  work1 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));
  work2 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));
  work3 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));
  work4 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));
  work5 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));
  work6 = (standard_precision *) Malloc0((nno+1) * sizeof(standard_precision));

  /* The coarse grid is uninjected to the fine and then the result is interpolated */

   if(2==dims && 2==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
     }

     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	 }
       }
     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {
	 node = j + (i-1) * noz;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1];
	 work1[node] = V1[node] = U1[node1];
	 work2[node] = V2[node] = U2[node1];	 
       }
 
     /* Fill in the X direction */
     for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2) {
 	 node=j + (i-1) * noz;
	 node1 = node-noz;
	 node2 = node+noz;	
	 node0 = (i != 2) ?  node-3*noz : node+5*noz; /* Use a different point if at edge of grid */
	 node3 = (i != (nox-1)) ? node+3*noz : node-5*noz ;

	 x = E->SX[level][1][node]; 
	 x0 = E->SX[level][1][node0];
	 x1 = E->SX[level][1][node1];
	 x2 = E->SX[level][1][node2];
	 x3 = E->SX[level][1][node3];

	 xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	 xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	 xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	 xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	 work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	 work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }
  
     /* Fill in the Z direction */
     for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
	 node=j + (i-1) * noz;
	 node1 = node-1;
	 node2 = node+1;
	 node0 = (j != 2) ?  node-3 : node+5; /* Use a different point if at edge of grid */
	 node3 = (j != noz-1) ? node+3 : node-5;

	 x = E->SX[level][2][node]; 
	 x0 = E->SX[level][2][node0];
	 x1 = E->SX[level][2][node1];
	 x2 = E->SX[level][2][node2];
	 x3 = E->SX[level][2][node3];

	 xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	 xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	 xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	 xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	 work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	 work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
       }

     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] ;
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] ;
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;
	 }
       }
   }
   else if(2==dims && 3==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     if(E->control.HAVE_SKEWBCS)
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	 }
       }
     /* Transfer common points */

     for(j=jj=1; j<=noz; jj++,j+=2)
       for(i=ii=1; i<=nox; ii++,i+=2) {
	 node = j + (i-1) * noz ;
	 node1 = jj + (ii-1) * E->mesh.NOZ[level-1] ;
	   work1[node] = V1[node] = U1[node1];
	   work2[node] = V2[node] = U2[node1];
	   work3[node] = V3[node] = U3[node1];
     }

     /* Fill in the X direction */
     for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2) {
 	 node=j + (i-1) * noz;
	 node1 = node-noz;
	 node2 = node+noz;	
	 node0 = (i != 2) ?  node-3*noz : node+5*noz; /* Use a different point if at edge of grid */
	 node3 = (i != (nox-1)) ? node+3*noz : node-5*noz ;

	 x = E->SX[level][1][node]; 
	 x0 = E->SX[level][1][node0];
	 x1 = E->SX[level][1][node1];
	 x2 = E->SX[level][1][node2];
	 x3 = E->SX[level][1][node3];

	 xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	 xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	 xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	 xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	 work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	 work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	 work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
	 V3[node] = work3[node];
       }

     /* Fill in the Z direction */
     for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
	 node=j + (i-1) * noz;
	 node1 = node-1;
	 node2 = node+1;
	 node0 = (j != 2) ?  node-3 : node+5; /* Use a different point if at edge of grid */
	 node3 = (j != noz-1) ? node+3 : node-5;

	 x = E->SX[level][2][node]; 
	 x0 = E->SX[level][2][node0];
	 x1 = E->SX[level][2][node1];
	 x2 = E->SX[level][2][node2];
	 x3 = E->SX[level][2][node3];

	 xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	 xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	 xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	 xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	 work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	 work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	 work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	 V1[node] = work1[node];
	 V2[node] = work2[node];
	 V3[node] = work3[node];
       }

     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] ;
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] ;
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;
	 }
       }
    }
   else  if(3==dims && 3==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	 }
       }

     for(k=kk=1; k<=noy; kk++,k+=2)  /* Initial un-injection */
      for(j=jj=1; j<=noz; jj++,j+=2)
	for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	}

     /* Fill in the X direction */
     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   node1 = node-noz;
	   node2 = node+noz;	
	   node0 = (i != 2) ?  node-3*noz : node+5*noz; /* Use a different point if at edge of grid */
	   node3 = (i != (nox-1)) ? node+3*noz : node-5*noz ;

	   x = E->SX[level][1][node]; 
	   x0 = E->SX[level][1][node0];
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   x3 = E->SX[level][1][node3];

	   xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	   xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	   xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	   xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	   work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	   work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	   work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
  	   V3[node] = work3[node];
      }

     /* Fill in the Z direction */
     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;

	   node1 = node-1;
	   node2 = node+1;
	   node0 = (j != 2) ?  node-3 : node+5; /* Use a different point if at edge of grid */
	   node3 = (j != noz-1) ? node+3 : node-5;

	   x = E->SX[level][2][node]; 
	   x0 = E->SX[level][2][node0];
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   x3 = E->SX[level][2][node3];

	   xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	   xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	   xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	   xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	   work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	   work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	   work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;

	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
      }

     /* Fill in the Y direction */
   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  node0 = (k != 2) ?  node-3* nox * noz : node+5* nox * noz; /* Use a different point if at edge of grid */
	  node3 = (k != noy-1) ? node+3* nox * noz : node-5* nox * noz;

	  x = E->SX[level][3][node]; 
	  x0 = E->SX[level][3][node0];
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  x3 = E->SX[level][3][node3];

	  xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	  xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	  xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	  xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	  work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	  work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	  work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;

	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	}
  
    /* if sphere, rotate boundary nodes back to rtf */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	   work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;
	   V3[i] = work3[i] ;
	 }
       }
   }


   else  if(3==dims && 6==dofs) {
     if(E->mesh.periodic_x || E->mesh.periodic_y) { 
       flogical_mesh_to_real(E,U1,level-1);
       flogical_mesh_to_real(E,U2,level-1);
       flogical_mesh_to_real(E,U3,level-1);
       flogical_mesh_to_real(E,U4,level-1);
       flogical_mesh_to_real(E,U5,level-1);
       flogical_mesh_to_real(E,U6,level-1);
     }

     /* if sphere, rotate boundary nodes to xyz */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level-1];i++) {
	 if(E->NODE[level-1][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level-1][i];
	   work1[i] = node_R[0*dims+0] * U1[i] + node_R[0*dims+1] * U2[i] + node_R[0*dims+2] * U3[i];
	   work2[i] = node_R[1*dims+0] * U1[i] + node_R[1*dims+1] * U2[i] + node_R[1*dims+2] * U3[i];
	   work3[i] = node_R[2*dims+0] * U1[i] + node_R[2*dims+1] * U2[i] + node_R[2*dims+2] * U3[i];
	   work4[i] = node_R[0*dims+0] * U4[i] + node_R[0*dims+1] * U5[i] + node_R[0*dims+2] * U6[i];
	   work5[i] = node_R[1*dims+0] * U4[i] + node_R[1*dims+1] * U5[i] + node_R[1*dims+2] * U6[i];
	   work6[i] = node_R[2*dims+0] * U4[i] + node_R[2*dims+1] * U5[i] + node_R[2*dims+2] * U6[i];
	   U1[i] = work1[i] ;
	   U2[i] = work2[i] ;
	   U3[i] = work3[i] ;
	   U4[i] = work4[i] ;
	   U5[i] = work5[i] ;
	   U6[i] = work6[i] ;
	 }
       }

     for(k=kk=1; k<=noy; kk++,k+=2)  /* Initial un-injection */
      for(j=jj=1; j<=noz; jj++,j+=2)
	for(i=ii=1; i<=nox; ii++,i+=2) {
	  node = j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = jj + (ii-1) * E->mesh.NOZ[level-1] + (kk-1) * E->mesh.NOZ[level-1] * E->mesh.NOX[level-1];
	  work1[node] = V1[node] = U1[node1];
	  work2[node] = V2[node] = U2[node1];
	  work3[node] = V3[node] = U3[node1];
	  work4[node] = V4[node] = U4[node1];
	  work5[node] = V5[node] = U5[node1];
	  work6[node] = V6[node] = U6[node1];
	}

     /* Fill in the X direction */
     for(k=1;k<=noy;k+=2)
       for(j=1;j<=noz;j+=2)
	 for(i=2;i<nox;i+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;
	   node1 = node-noz;
	   node2 = node+noz;
	   node1 = node-noz;
	   node2 = node+noz;	
	   node0 = (i != 2) ?  node-3*noz : node+5*noz; /* Use a different point if at edge of grid */
	   node3 = (i != (nox-1)) ? node+3*noz : node-5*noz ;

	   x = E->SX[level][1][node]; 
	   x0 = E->SX[level][1][node0];
	   x1 = E->SX[level][1][node1];
	   x2 = E->SX[level][1][node2];
	   x3 = E->SX[level][1][node3];

	   xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	   xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	   xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	   xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	   work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	   work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	   work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	   work4[node] = work4[node0]*xx0 + work4[node1]*xx1 + work4[node2]*xx2 + work4[node3]*xx3;
	   work5[node] = work5[node0]*xx0 + work5[node1]*xx1 + work5[node2]*xx2 + work5[node3]*xx3;
	   work6[node] = work6[node0]*xx0 + work6[node1]*xx1 + work6[node2]*xx2 + work6[node3]*xx3;
	   V1[node] = work1[node];
	   V2[node] = work2[node];
  	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
  	   V6[node] = work6[node];
      }

     /* Fill in the Z direction */
     for(k=1;k<=noy;k+=2)
       for(i=1;i<=nox;i++)
	 for(j=2;j<noz;j+=2) {
	   node=j + (i-1) * noz + (k-1) * nox * noz;

	   node1 = node-1;
	   node2 = node+1;
	   node0 = (j != 2) ?  node-3 : node+5; /* Use a different point if at edge of grid */
	   node3 = (j != noz-1) ? node+3 : node-5;

	   x = E->SX[level][2][node]; 
	   x0 = E->SX[level][2][node0];
	   x1 = E->SX[level][2][node1];
	   x2 = E->SX[level][2][node2];
	   x3 = E->SX[level][2][node3];

	   xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	   xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	   xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	   xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	   work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	   work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	   work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	   work4[node] = work4[node0]*xx0 + work4[node1]*xx1 + work4[node2]*xx2 + work4[node3]*xx3;
	   work5[node] = work5[node0]*xx0 + work5[node1]*xx1 + work5[node2]*xx2 + work5[node3]*xx3;
	   work6[node] = work6[node0]*xx0 + work6[node1]*xx1 + work6[node2]*xx2 + work6[node3]*xx3;

	   V1[node] = work1[node];
	   V2[node] = work2[node];
	   V3[node] = work3[node];
	   V4[node] = work4[node];
	   V5[node] = work5[node];
	   V6[node] = work6[node];
      }

     /* Fill in the Y direction */
   for(i=1;i<=nox;i++)
      for(j=1;j<=noz;j++)  
	for(k=2;k<noy;k+=2) {
	  node=j + (i-1) * noz + (k-1) * nox * noz;
	  node1 = node-nox*noz;
	  node2 = node+nox*noz;
	  node0 = (k != 2) ?  node-3* nox * noz : node+5* nox * noz; /* Use a different point if at edge of grid */
	  node3 = (k != noy-1) ? node+3* nox * noz : node-5* nox * noz;

	  x = E->SX[level][3][node]; 
	  x0 = E->SX[level][3][node0];
	  x1 = E->SX[level][3][node1];
	  x2 = E->SX[level][3][node2];
	  x3 = E->SX[level][3][node3];

	  xx0= (x-x1)*(x-x2)*(x-x3) / ((x0-x1)*(x0-x2)*(x0-x3));
	  xx1= (x-x0)*(x-x2)*(x-x3) / ((x1-x0)*(x1-x2)*(x1-x3));
	  xx2= (x-x0)*(x-x1)*(x-x3) / ((x2-x0)*(x2-x1)*(x2-x3));
	  xx3= (x-x0)*(x-x1)*(x-x2) / ((x3-x0)*(x3-x1)*(x3-x2));	 

	  work1[node] = work1[node0]*xx0 + work1[node1]*xx1 + work1[node2]*xx2 + work1[node3]*xx3;
	  work2[node] = work2[node0]*xx0 + work2[node1]*xx1 + work2[node2]*xx2 + work2[node3]*xx3;
	  work3[node] = work3[node0]*xx0 + work3[node1]*xx1 + work3[node2]*xx2 + work3[node3]*xx3;
	  work4[node] = work4[node0]*xx0 + work4[node1]*xx1 + work4[node2]*xx2 + work4[node3]*xx3;
	  work5[node] = work5[node0]*xx0 + work5[node1]*xx1 + work5[node2]*xx2 + work5[node3]*xx3;
	  work6[node] = work6[node0]*xx0 + work6[node1]*xx1 + work6[node2]*xx2 + work6[node3]*xx3;

	  V1[node] = work1[node];
	  V2[node] = work2[node];
	  V3[node] = work3[node];
	  V4[node] = work4[node];
	  V5[node] = work5[node];
	  V6[node] = work6[node];
	}
  
    /* if sphere, rotate boundary nodes back to rtf */
     
     if(E->control.HAVE_SKEWBCS) 
       for(i=1;i<=E->mesh.NNO[level];i++) {
	 if(E->NODE[level][i] & SKEWBC) {
	   node_R = E->curvilinear.NODE_R[level][i];
	   work1[i] = node_R[0*dims+0] * V1[i] + node_R[1*dims+0] * V2[i] + node_R[2*dims+0] * V3[i];
	   work2[i] = node_R[0*dims+1] * V1[i] + node_R[1*dims+1] * V2[i] + node_R[2*dims+1] * V3[i];
	   work3[i] = node_R[0*dims+2] * V1[i] + node_R[1*dims+2] * V2[i] + node_R[2*dims+2] * V3[i];
	   work4[i] = node_R[0*dims+0] * V4[i] + node_R[1*dims+0] * V5[i] + node_R[2*dims+0] * V6[i];
	   work5[i] = node_R[0*dims+1] * V4[i] + node_R[1*dims+1] * V5[i] + node_R[2*dims+1] * V6[i];
	   work6[i] = node_R[0*dims+2] * V4[i] + node_R[1*dims+2] * V5[i] + node_R[2*dims+2] * V6[i];
	   V1[i] = work1[i] ;
	   V2[i] = work2[i] ;
	   V3[i] = work3[i] ;
	   V4[i] = work4[i] ;
	   V5[i] = work5[i] ;
	   V6[i] = work6[i] ;
	 }
       }
   }

   free((void *) work1);
   free((void *) work2);
   free((void *) work3);
   free((void *) work4);
   free((void *) work5);
   free((void *) work6);
return;
}

standard_precision residual_del2_u6(
  struct All_variables *E,
  standard_precision *F1,
  standard_precision *F2,
  standard_precision *F3,
  standard_precision *F4,
  standard_precision *F5,
  standard_precision *F6,
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
  standard_precision *R6,
  standard_precision *Au1,
  standard_precision *Au2,
  standard_precision *Au3,
  standard_precision *Au4,
  standard_precision *Au5,
  standard_precision *Au6,
  int level,
  int sd
)
{
  standard_precision res,res2,x1,x2,x3,au1,au2,au3;
  standard_precision alpha,AudotAu,AudotR;
  int i;
  const int dofs=E->mesh.dof ;
  
  AudotR = 0.0;
  AudotAu = 0.0;

  if(sd) {
    switch(dofs) {
    case 2:
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(E->NODE[level][i] & ( OFFSIDE ))
	  continue;
	AudotR += Au1[i] * F1[i] + Au2[i] * F2[i];
	AudotAu += Au1[i] * Au1[i] + Au2[i] * Au2[i];
      }
      break ;
    case 3:
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(E->NODE[level][i] & ( OFFSIDE ))
	  continue;
	AudotR += Au1[i] * F1[i] + Au2[i] * F2[i]+ Au3[i] * F3[i];
	AudotAu += Au1[i] * Au1[i] + Au2[i] * Au2[i] + Au3[i] * Au3[i];
      }
      break ;
    case 6:
      for(i=1;i<=E->mesh.NNO[level];i++) {
	if(E->NODE[level][i] & ( OFFSIDE ))
	  continue;
	AudotR += Au1[i]*F1[i] + Au2[i]*F2[i]+ Au3[i]*F3[i] + Au4[i]*F4[i] + Au5[i]*F5[i]+ Au6[i]*F6[i];
	AudotAu += Au1[i]*Au1[i] + Au2[i]*Au2[i] + Au3[i]*Au3[i] + Au4[i]*Au4[i] + Au5[i]*Au5[i] + Au6[i]*Au6[i];
      }
      break ;
    }

    if(AudotAu != 0.0)
      alpha = AudotR / AudotAu   ;
    else 
      alpha = 1.0;
  }
  else
    alpha = 1.0;

  res=0.0;

  switch(dofs) {
  case 2:
    for(i=1;i<=E->mesh.NNO[level];i++) {
      if(E->NODE[level][i] & ( OFFSIDE )) {
	  R1[i] = R2[i] =0.0;
	continue;
      }
    
      V1[i] *= alpha;
      V2[i] *= alpha;

      Au1[i] *= alpha;
      Au2[i] *= alpha; 

      R1[i] = F1[i] - Au1[i];
      R2[i] = F2[i] - Au2[i];

      /*  if(level==E->mesh.levmin)
	  fprintf(stderr,"Node %d: F = %g, Res = %g,%g \n",i,F1[i], R1[i], R2[i]); */


      res += R1[i]*R1[i]+R2[i]*R2[i];
    }
    break ;
  case 3:
    for(i=1;i<=E->mesh.NNO[level];i++) {
      if(E->NODE[level][i] & ( OFFSIDE )) {
	  R1[i] = R2[i] = R3[i] = 0.0;
	continue;
      }
       V1[i] *= alpha;
       V2[i] *= alpha;
       V3[i] *= alpha;

       Au1[i] *= alpha;
       Au2[i] *= alpha; 
       Au3[i] *= alpha; 

       R1[i] = F1[i] - Au1[i];
       R2[i] = F2[i] - Au2[i];
       R3[i] = F3[i] - Au3[i];
       res += R1[i]*R1[i]+R2[i]*R2[i] + R3[i]*R3[i];
    }
    break ;
  case 6:
    for(i=1;i<=E->mesh.NNO[level];i++) {
      if(E->NODE[level][i] & ( OFFSIDE )) {
	  R1[i] = R2[i] = R3[i] = R4[i] = R5[i] = R6[i] = 0.0;
	continue;
      }
       V1[i] *= alpha;
       V2[i] *= alpha;
       V3[i] *= alpha;
       V4[i] *= alpha;
       V5[i] *= alpha;
       V6[i] *= alpha;

       Au1[i] *= alpha;
       Au2[i] *= alpha; 
       Au3[i] *= alpha; 
       Au4[i] *= alpha;
       Au5[i] *= alpha; 
       Au6[i] *= alpha; 

       R1[i] = F1[i] - Au1[i];
       R2[i] = F2[i] - Au2[i];
       R3[i] = F3[i] - Au3[i];
       R4[i] = F4[i] - Au4[i];
       R5[i] = F5[i] - Au5[i];
       R6[i] = F6[i] - Au6[i];

   res += R1[i]*R1[i] + R2[i]*R2[i] + R3[i]*R3[i] + R4[i]*R4[i] + R5[i]*R5[i] + R6[i]*R6[i];
    }
    break ;
  }

  res = sqrt(res/E->mesh.NNO[level]);
  return(res);
}
