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


 

/*  Functions which construct the shape function values at all of the gauss
    points in the element (including the reduced quadrature points). The element in question is
    biquadratic in the velocities and therefore bilinear in the pressures. 
    
    To change elements it is necessary to change this file: Shape_functions.c,
    and the element-data header file : element_definitions.h  but it should not be
    necessary to change the main calculation/setup/solving machinery.		 */

#include <math.h>
#include "element_definitions.h"				
#include "global_defs.h"
 
/*  =======================================================
    Function creating shape_fn data in form of a structure
    =======================================================*/

void construct_shape_functions(
     struct All_variables *E
)
{	
  higher_precision lpoly(),lpolydash();
  int i,j,k,d,dd;
  int remapj,remapk;

  /* first zero ALL entries, even those not used in 2d. */
  /* GNVI = 8*8 because 8 is the maximum number of node in the element 
     and 8 is also the number of velocity variable in the element */

  for(i=0;i<GNVI;i++)
    { E->N.vpt[i] = 0.0; 
      E->Nx.vpt[i] = 0.0;
      E->Nx.vpt[GNVI+i] = 0.0;
      E->Nx.vpt[2*GNVI+i] = 0.0; 
    }

  /* GNPI = 8*1 because 8 is the maximum number of node in the element
     and 1 is the number of pressure variable in the element */
  
   for(i=0;i<GNPI;i++)
    { E->N.ppt[i] = 0.0; 
      E->Nx.ppt[i] = 0.0;
      E->Nx.ppt[GNPI+i] = 0.0;
      E->Nx.ppt[2*GNPI+i] = 0.0; 
    }
  /* GN1VI = 4*4 */

  for(i=0;i<GN1VI;i++)
    { E->M.vpt[i] = 0.0; 
      E->Mx.vpt[i] = 0.0;
      E->Mx.vpt[GN1VI+i] = 0.0;
    }

  /* GN1PI = 4*1 */
  
   for(i=0;i<GN1PI;i++)
    { E->M.ppt[i] = 0.0; 
      E->Mx.ppt[i] = 0.0;
      E->Mx.ppt[GN1PI+i] = 0.0;
    }
  
  for(i=1;i<=enodes[E->mesh.nsd];i++) {
      for(j=1;j<=vpoints[E->mesh.nsd];j++)  {
	  /* for each integration point vpoints={0,2,4,8} */

	  E->N.vpt[GNVINDEX(i,j)] = 1.0;
	  for(d=1;d<=E->mesh.nsd;d++)  {
	      E->N.vpt[GNVINDEX(i,j)] *=  
		  lpoly(bb[d-1][i],g_point[j].x[d-1]);
	  }
	  for(dd=1;dd<=E->mesh.nsd;dd++) {
	      E->Nx.vpt[GNVXINDEX(dd-1,i,j)] = lpolydash(bb[dd-1][i],g_point[j].x[dd-1]);
	      for(d=1;d<=E->mesh.nsd;d++)
		  if (d != dd)
		      E->Nx.vpt[GNVXINDEX(dd-1,i,j)] *= lpoly(bb[d-1][i],g_point[j].x[d-1]);
	  }
      }
      /* next velocity point */
 
      for(j=1;j<=ppoints[E->mesh.nsd];j++)  {/*2*/
	  /* for each p-integration point ppoints={0,1,1,1}  */
	/* never used for the time being */

	  E->N.ppt[GNPINDEX(i,j)] = 1.0;
	  for(d=1;d<=E->mesh.nsd;d++){
	      E->N.ppt[GNPINDEX(i,j)] *=  
		  lpoly(bb[d-1][i],p_point[j].x[d-1]);
	  }
	  for(dd=1;dd<=E->mesh.nsd;dd++) {
	      E->Nx.ppt[GNPXINDEX(dd-1,i,j)] = lpolydash(bb[dd-1][i],p_point[j].x[dd-1]);
	      for(d=1;d<=E->mesh.nsd;d++)
		  if (d != dd) {
		      E->Nx.ppt[GNPXINDEX(dd-1,i,j)] *= lpoly(bb[d-1][i],p_point[j].x[d-1]); 
	  }
	}
    }/* next pressure node */
  }/* next node */
      
  /*	1d cases ... set up M (Mx is not used in FE formulation but is used for dxdxsi &c)  */

  /* note: naughty boy wedgie, this works best for spatial ordering (z,x,y)
     so we actually need to read the z part if it's one-d and the x,z part if
     it's two. ' Comes about from the node numbering really and so we need the
     remapping of j in the 1d case. */

  for(j=1;j<=onedvpoints[E->mesh.nsd];j++) {
      remapj = ccc[E->mesh.nsd-2][j];
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++) {
	  remapk = ccc[E->mesh.nsd-2][k];
	  E->M.vpt[GMVINDEX(j,k)] = 1.0;
	  for(d=1;d<=E->mesh.nsd-1;d++)
	      E->M.vpt[GMVINDEX(j,k)] *= lpoly(bb[d-1][remapj],g_1d[remapk].x[d-1]);
	  for(dd=1;dd<=E->mesh.nsd-1;dd++) {
	      E->Mx.vpt[GMVXINDEX(dd-1,j,k)] = lpolydash(bb[dd-1][remapj],g_1d[remapk].x[dd-1]); 
	      for(d=1;d<=E->mesh.nsd-1;d++)
		if (d != dd)
		  E->Mx.vpt[GMVXINDEX(dd-1,j,k)] *= lpoly(bb[d-1][remapj],g_1d[remapk].x[d-1]); 
	      }
	  }
    }
  return; }
		
higher_precision lpoly(
     int p,	   /*   selects lagrange polynomial , 1d: node p */
     higher_precision y  /*   coordinate in given direction to evaluate poly */
)
{
  higher_precision value;
  
  switch (p)
    {
    case 1:
      value = 0.5 * (1-y) ;
      break;
    case 2:
      value = 0.5 * (1+y) ;
      break;
    default:
      value = 0.0;
    }
  return(value);
}
	
higher_precision lpolydash(
     int p,
     higher_precision y
)
{	
  higher_precision value;
  switch (p)
    {
    case 1:
      value = -0.5 ;
      break;
    case 2:
      value =  0.5 ;
      break;
    default:
      value = 0.0;
    }
  return(value);
}

void sp_to_nodes(
  struct All_variables *E,
  standard_precision *P,
  standard_precision *PN,
  int lev
)
{
  int p,element,node,node1,i,j,k;

  const int nox = E->mesh.NOX[lev];
  const int noz = E->mesh.NOZ[lev];
  const int noy = E->mesh.NOY[lev];

  standard_precision *WN;

  WN = (standard_precision *) Malloc0((E->mesh.NNO[lev] + 1) * sizeof(standard_precision));

  for(p=1;p<=E->mesh.NNO[lev];p++) PN[p] = WN[p] = 0.0;
  
  for(p=1;p<=E->mesh.NEL[lev];p++){
    if( E->tracer.tr_in_element_number[lev][p] == 0) continue;

    element = p;
    for(j=1;j<=enodes[E->mesh.nsd];j++) {
      node = E->IEN[lev][element].node[j];
      PN[node] += P[element];  
      WN[node] += 1.0;
    }
  }

  for(p=1;p<=E->mesh.NNO[lev];p++) {
    if( WN[p] != 0.0) 
      PN[p] /= WN[p];
  }

  free((void *) WN);
  return;

  /* Edge correction */

   /* X direction */

  for(k=1;k<=noy;k++)
    for(j=1;j<=noz;j++) {
      node  = j + (k-1) * noz * nox ;
      node1 = j + (nox-1) * noz + (k-1) * noz * nox ;
      PN[node]  = 2 * PN[node]  - PN[node+noz];
      PN[node1] = 2 * PN[node1] - PN[node1-noz];
    }

  /* Z direction */

  for(k=1;k<=noy;k++)
    for(i=1;i<=nox;i++) { 
      node = 1 + (i-1) * noz + (k-1) * noz * nox ;
      node1 = node + noz - 1;
      PN[node]  = 2 * PN[node]  - PN[node+1];
      PN[node1] = 2 * PN[node1] - PN[node1-1];
    }

   /* Y direction */

  if(3==E->mesh.nsd)
    for(i=1;i<=nox;i++) 
      for(j=1;j<=noz;j++) {
	node  = j + (i-1) * noz;
	node1 = j + (i-1) * noz + (noy-1) * noz * nox ;
      PN[node]  = 2 * PN[node]  - PN[node+noz*nox];
      PN[node1] = 2 * PN[node1] - PN[node1-noz*nox];
    }

  /* Corners are still wrong ! */

  /* TWO DIM correction !! */
  /* top left */
  PN[1] = PN[2] - 0.5 * PN[3] + PN[1+noz] - 0.5 * PN[1+2*noz]; 
  /* bottom left */
  PN[noz] = PN[noz-1] - 0.5 * PN[noz-2] + PN[noz+noz] - 0.5 * PN[noz+2*noz];
  /* top right */
  PN[1+(nox-1)*noz] = PN[1+(nox-1)*noz+1] - 0.5 * PN[1+(nox-1)*noz+2] + 
    PN[1+(nox-2)*noz] - 0.5 * PN[1+(nox-3)*noz]; 
  /* bottom right */
  PN[noz*nox] =  PN[noz*nox -1] - 0.5 * PN[noz*nox -2] + PN[noz*nox-noz] - 0.5 * PN[noz*nox-2*noz];

  /* NOTE this is still not correct for three dimensional case !!! */
     return;
 }

void sp_to_centres(
     struct All_variables *E,
     standard_precision *PN,
     standard_precision *P,
     int lev
)
{  int p,element,node,j;
   higher_precision weight;

   for(p=1;p<=E->mesh.NEL[lev];p++)
     P[p] = 0.0;

   weight=1.0/((higher_precision)enodes[E->mesh.nsd]) ;
   
   for(p=1;p<=E->mesh.NEL[lev];p++) {
     for(j=1;j<=enodes[E->mesh.nsd];j++) {
        P[p] +=  PN[E->IEN[lev][p].node[j]] * weight;
     }
   }
   return;  
}

void v_shape_fn(
		struct All_variables *E,
		int el,
		standard_precision *lN,
		standard_precision eta1,
		standard_precision eta2,
		standard_precision eta3,
		int level
		)
      /* For further details on this implementation see BATHE's book on page 200 */
{
  const int dims = E->mesh.nsd;
  int i;

  if(2==dims) {
    for(i=1;i<=9;i++)
      lN[i] = 0.0 ;

    switch(E->ELEM[level][el].type_v) {
    case NINE_NODES_QUAD:
      lN[9] = (1-eta1*eta1) * (1-eta2*eta2) ;
    case EIGHT_NODES_QUAD:
      lN[8] = 0.5 * (1-eta1*eta1) * (1-eta2) - 0.5*lN[9] ;
    case SEVEN_NODES_QUAD:
      lN[6] = 0.5 * (1-eta1*eta1) * (1+eta2) - 0.5*lN[9] ;
    case SIX_NODES_QUAD: /* The sixth node is on the opposite side than the fifth,
			    that's why we compute the node-7 shape function  */
      lN[7] = 0.5 * (1+eta1) * (1-eta2*eta2) - 0.5*lN[9] ;

    case FIVE_NODES_QUAD:
      lN[5] = 0.5 * (1-eta1) * (1-eta2*eta2) - 0.5*lN[9] ;

/*    case SIX_NODES_TRIA : */
    case FOUR_NODES_QUAD:
      lN[4] = 0.25 * (1.0+eta1) * (1.0-eta2) - 0.5*lN[6] - 0.5*lN[8] - 0.25*lN[9] ;
      lN[3] = 0.25 * (1.0+eta1) * (1.0+eta2) - 0.5*lN[6] - 0.5*lN[7] - 0.25*lN[9] ;
      lN[2] = 0.25 * (1.0-eta1) * (1.0+eta2) - 0.5*lN[5] - 0.5*lN[7] - 0.25*lN[9] ;
      lN[1] = 0.25 * (1.0-eta1) * (1.0-eta2) - 0.5*lN[5] - 0.5*lN[8] - 0.25*lN[9] ;
    }
  }
  else {
    for(i=1;i<=27;i++) lN[i] = 0.0 ;
    switch(E->ELEM[level][el].type_v) {
    case TWENTY_SEVEN_NODES_CUBIC:
     lN[27] = (1-eta1*eta1) * (1-eta2*eta2) * (1-eta3*eta3) ;
    case TWENTY_SIX_NODES_CUBIC:
     lN[26] = 0.5 * (1-eta1*eta1) * (1-eta2*eta2) * (1-eta3) - 0.5*lN[27] ;
     lN[25] = 0.5 * (1-eta1*eta1) * (1-eta2*eta2) * (1+eta3) - 0.5*lN[27] ;
     lN[24] = 0.5 * (1+eta1) * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lN[27] ;
     lN[23] = 0.5 * (1-eta1*eta1) * (1+eta2) * (1-eta3*eta3) - 0.5*lN[27] ;
     lN[22] = 0.5 * (1-eta1) * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lN[27] ;
     lN[21] = 0.5 * (1-eta1*eta1) * (1-eta2) * (1-eta3*eta3) - 0.5*lN[27] ;
    case TWENTY_NODES_CUBIC:
     lN[20] = 0.25 * (1+eta1) * (1-eta2) * (1-eta3*eta3) - 0.5*(lN[21]+lN[24]) - 0.25*lN[27] ;
     lN[19] = 0.25 * (1+eta1) * (1+eta2) * (1-eta3*eta3) - 0.5*(lN[23]+lN[24]) - 0.25*lN[27] ;
     lN[18] = 0.25 * (1-eta1) * (1+eta2) * (1-eta3*eta3) - 0.5*(lN[22]+lN[23]) - 0.25*lN[27] ;
     lN[17] = 0.25 * (1-eta1) * (1-eta2) * (1-eta3*eta3) - 0.5*(lN[21]+lN[22]) - 0.25*lN[27] ;
     lN[16] = 0.25 * (1-eta1*eta1) * (1-eta2) * (1+eta3) - 0.5*(lN[21]+lN[25]) - 0.25*lN[27] ;
     lN[15] = 0.25 * (1+eta1) * (1-eta2*eta2) * (1+eta3) - 0.5*(lN[24]+lN[25]) - 0.25*lN[27] ;
     lN[14] = 0.25 * (1-eta1*eta1) * (1+eta2) * (1+eta3) - 0.5*(lN[23]+lN[25]) - 0.25*lN[27] ;
     lN[13] = 0.25 * (1-eta1) * (1-eta2*eta2) * (1+eta3) - 0.5*(lN[22]+lN[25]) - 0.25*lN[27] ;
     lN[12] = 0.25 * (1-eta1*eta1) * (1-eta2) * (1-eta3) - 0.5*(lN[21]+lN[26]) - 0.25*lN[27] ;
     lN[11] = 0.25 * (1+eta1) * (1-eta2*eta2) * (1-eta3) - 0.5*(lN[24]+lN[26]) - 0.25*lN[27] ;
     lN[10] = 0.25 * (1-eta1*eta1) * (1+eta2) * (1-eta3) - 0.5*(lN[23]+lN[26]) - 0.25*lN[27] ;
     lN[9] = 0.25 * (1-eta1) * (1-eta2*eta2) * (1-eta3) - 0.5*(lN[22]+lN[26]) - 0.25*lN[27] ;
    case EIGHT_NODES_CUBIC: /*RAA, comment out higher order node contributions as they are unnecessary*/
     lN[8] = 0.125 * (1.0+eta1) * (1.0-eta2) * (1.0+eta3); /* - 0.5*(lN[15]+lN[16]+lN[20]) - 0.25*(lN[21]+lN[24]+lN[25]) - 0.125*lN[27] ; */
     lN[7] = 0.125 * (1.0+eta1) * (1.0+eta2) * (1.0+eta3); /* - 0.5*(lN[14]+lN[15]+lN[19]) - 0.25*(lN[23]+lN[24]+lN[25]) - 0.125*lN[27] ; */
     lN[6] = 0.125 * (1.0-eta1) * (1.0+eta2) * (1.0+eta3); /* - 0.5*(lN[13]+lN[14]+lN[18]) - 0.25*(lN[22]+lN[23]+lN[25]) - 0.125*lN[27] ; */
     lN[5] = 0.125 * (1.0-eta1) * (1.0-eta2) * (1.0+eta3); /* - 0.5*(lN[13]+lN[16]+lN[17]) - 0.25*(lN[21]+lN[22]+lN[25]) - 0.125*lN[27] ; */
     lN[4] = 0.125 * (1.0+eta1) * (1.0-eta2) * (1.0-eta3); /* - 0.5*(lN[11]+lN[12]+lN[20]) - 0.25*(lN[21]+lN[24]+lN[26]) - 0.125*lN[27] ; */
     lN[3] = 0.125 * (1.0+eta1) * (1.0+eta2) * (1.0-eta3); /* - 0.5*(lN[10]+lN[11]+lN[19]) - 0.25*(lN[23]+lN[24]+lN[26]) - 0.125*lN[27] ; */
     lN[2] = 0.125 * (1.0-eta1) * (1.0+eta2) * (1.0-eta3); /* - 0.5*(lN[9]+lN[10]+lN[18]) - 0.25*(lN[22]+lN[23]+lN[26]) - 0.125*lN[27] ; */
     lN[1] = 0.125 * (1.0-eta1) * (1.0-eta2) * (1.0-eta3); /* - 0.5*(lN[9]+lN[12]+lN[17]) - 0.25*(lN[21]+lN[22]+lN[26]) - 0.125*lN[27] ; */
  }
  }
  return;
}

/* Vector version of shape functions */


void all_v_shape_fn(
  struct All_variables *E,
  struct TRACER_ELT_WEIGHT *lN,
  standard_precision *eta1,
  standard_precision *eta2,
  standard_precision *eta3,
  int N1,
  int N2,
  int level
)
      /* For further details on this implementation see BATHE's book on page 200 */
{
  const int dims = E->mesh.nsd;
  int i,m;

  if(2==dims) {
    for(i=1;i<=9;i++)
      for(m=N1;m<=N2;m++) {
	lN[m].node[i] = 0.0;
	

      }

    switch(E->control.ELEMENT_TYPE) {
    case NINE_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) 
	lN[m].node[9] = (1-eta1[m]*eta1[m]) * (1-eta2[m]*eta2[m]) ;
    case EIGHT_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) 
	lN[m].node[8] = 0.5 * (1-eta1[m]*eta1[m]) * (1-eta2[m]) - 0.5 * lN[m].node[9] ;
    case SEVEN_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) 
	lN[m].node[6] = 0.5 * (1-eta1[m]*eta1[m]) * (1+eta2[m]) - 0.5*lN[m].node[9] ;
    case SIX_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) 
	lN[m].node[7] = 0.5 * (1+eta1[m]) * (1-eta2[m]*eta2[m]) - 0.5*lN[m].node[9] ;
    case FIVE_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) 
	lN[m].node[5] = 0.5 * (1-eta1[m]) * (1-eta2[m]*eta2[m]) - 0.5*lN[m].node[9] ;
/*  case SIX_NODES_TRIA : */
    case FOUR_NODES_QUAD:
#pragma loop novrec lN,eta1,eta2
      for(m=N1;m<=N2;m++) {
       lN[m].node[4] = 0.25 * (1.0+eta1[m]) * (1.0-eta2[m]) - 
	 0.5*lN[m].node[6] - 0.5*lN[m].node[8] - 0.25*lN[m].node[9] ;
       lN[m].node[3] = 0.25 * (1.0+eta1[m]) * (1.0+eta2[m]) - 
	 0.5*lN[m].node[6] - 0.5*lN[m].node[7] - 0.25*lN[m].node[9] ;
       lN[m].node[2] = 0.25 * (1.0-eta1[m]) * (1.0+eta2[m]) - 
	 0.5*lN[m].node[5] - 0.5*lN[m].node[7] - 0.25*lN[m].node[9] ;
       lN[m].node[1] = 0.25 * (1.0-eta1[m]) * (1.0-eta2[m]) - 
	 0.5*lN[m].node[5] - 0.5*lN[m].node[8] - 0.25*lN[m].node[9] ;
     }
    }
  }

/*#if 0 */  /*RAA: 23/3/01, fixed from #if 0 and vectorized as with 2D case */
  else if(3==dims) {
    for(i=1;i<=27;i++)
       for(m=N1;m<=N2;m++) 
          lN[m].node[i] = 0.0;

    switch(E->control.ELEMENT_TYPE) {
    case TWENTY_SEVEN_NODES_CUBIC:
#pragma loop novrec lN,eta1,eta2,eta3
      for(m=N1;m<=N2;m++) {
         lN[m].node[27] = (1-eta1[m]*eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]*eta3[m]) ;
      }
    case TWENTY_SIX_NODES_CUBIC:
#pragma loop novrec lN,eta1,eta2,eta3
      for(m=N1;m<=N2;m++) {
	 lN[m].node[26] = 0.5 * (1-eta1[m]*eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]) - 0.5*lN[m].node[27] ;
	 lN[m].node[25] = 0.5 * (1-eta1[m]*eta1[m]) * (1-eta2[m]*eta2[m]) * (1+eta3[m]) - 0.5*lN[m].node[27] ;
	 lN[m].node[24] = 0.5 * (1+eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*lN[m].node[27] ;
	 lN[m].node[23] = 0.5 * (1-eta1[m]*eta1[m]) * (1+eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*lN[m].node[27] ;
	 lN[m].node[22] = 0.5 * (1-eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*lN[m].node[27] ;
	 lN[m].node[21] = 0.5 * (1-eta1[m]*eta1[m]) * (1-eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*lN[m].node[27] ;
      }
    case TWENTY_NODES_CUBIC:
#pragma loop novrec lN,eta1,eta2,eta3
      for(m=N1;m<=N2;m++) {
	 lN[m].node[20] = 0.25 * (1+eta1[m]) * (1-eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*(lN[m].node[21]+lN[m].node[24]) - 0.25*lN[m].node[27] ;
	 lN[m].node[19] = 0.25 * (1+eta1[m]) * (1+eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*(lN[m].node[23]+lN[m].node[24]) - 0.25*lN[m].node[27] ;
	 lN[m].node[18] = 0.25 * (1-eta1[m]) * (1+eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*(lN[m].node[22]+lN[m].node[23]) - 0.25*lN[m].node[27] ;
	 lN[m].node[17] = 0.25 * (1-eta1[m]) * (1-eta2[m]) * (1-eta3[m]*eta3[m]) - 0.5*(lN[m].node[21]+lN[m].node[22]) - 0.25*lN[m].node[27] ;
	 lN[m].node[16] = 0.25 * (1-eta1[m]*eta1[m]) * (1-eta2[m]) * (1+eta3[m]) - 0.5*(lN[m].node[21]+lN[m].node[25]) - 0.25*lN[m].node[27] ;
	 lN[m].node[15] = 0.25 * (1+eta1[m]) * (1-eta2[m]*eta2[m]) * (1+eta3[m]) - 0.5*(lN[m].node[24]+lN[m].node[25]) - 0.25*lN[m].node[27] ;
	 lN[m].node[14] = 0.25 * (1-eta1[m]*eta1[m]) * (1+eta2[m]) * (1+eta3[m]) - 0.5*(lN[m].node[23]+lN[m].node[25]) - 0.25*lN[m].node[27] ;
	 lN[m].node[13] = 0.25 * (1-eta1[m]) * (1-eta2[m]*eta2[m]) * (1+eta3[m]) - 0.5*(lN[m].node[22]+lN[m].node[25]) - 0.25*lN[m].node[27] ;
	 lN[m].node[12] = 0.25 * (1-eta1[m]*eta1[m]) * (1-eta2[m]) * (1-eta3[m]) - 0.5*(lN[m].node[21]+lN[m].node[26]) - 0.25*lN[m].node[27] ;
	 lN[m].node[11] = 0.25 * (1+eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]) - 0.5*(lN[m].node[24]+lN[m].node[26]) - 0.25*lN[m].node[27] ;
	 lN[m].node[10] = 0.25 * (1-eta1[m]*eta1[m]) * (1+eta2[m]) * (1-eta3[m]) - 0.5*(lN[m].node[23]+lN[m].node[26]) - 0.25*lN[m].node[27] ;
	 lN[m].node[9] = 0.25 * (1-eta1[m]) * (1-eta2[m]*eta2[m]) * (1-eta3[m]) - 0.5*(lN[m].node[22]+lN[m].node[26]) - 0.25*lN[m].node[27] ;
      }
    case EIGHT_NODES_CUBIC:
#pragma loop novrec lN,eta1,eta2,eta3
      for(m=N1;m<=N2;m++) {
	 lN[m].node[8] = 0.125 * (1.0+eta1[m]) * (1.0-eta2[m]) * (1.0+eta3[m]);/* - 0.5*(lN[m].node[15]+lN[m].node[16]+lN[m].node[20]) - 0.25*(lN[m].node[21]+lN[m].node[24]+lN[m].node[25]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[7] = 0.125 * (1.0+eta1[m]) * (1.0+eta2[m]) * (1.0+eta3[m]); /*- 0.5*(lN[m].node[14]+lN[m].node[15]+lN[m].node[19]) - 0.25*(lN[m].node[23]+lN[m].node[24]+lN[m].node[25]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[6] = 0.125 * (1.0-eta1[m]) * (1.0+eta2[m]) * (1.0+eta3[m]);/* - 0.5*(lN[m].node[13]+lN[m].node[14]+lN[m].node[18]) - 0.25*(lN[m].node[22]+lN[m].node[23]+lN[m].node[25]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[5] = 0.125 * (1.0-eta1[m]) * (1.0-eta2[m]) * (1.0+eta3[m]);/* - 0.5*(lN[m].node[13]+lN[m].node[16]+lN[m].node[17]) - 0.25*(lN[m].node[21]+lN[m].node[22]+lN[m].node[25]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[4] = 0.125 * (1.0+eta1[m]) * (1.0-eta2[m]) * (1.0-eta3[m]);/* - 0.5*(lN[m].node[11]+lN[m].node[12]+lN[m].node[20]) - 0.25*(lN[m].node[21]+lN[m].node[24]+lN[m].node[26]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[3] = 0.125 * (1.0+eta1[m]) * (1.0+eta2[m]) * (1.0-eta3[m]);/* - 0.5*(lN[m].node[10]+lN[m].node[11]+lN[m].node[19]) - 0.25*(lN[m].node[23]+lN[m].node[24]+lN[m].node[26]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[2] = 0.125 * (1.0-eta1[m]) * (1.0+eta2[m]) * (1.0-eta3[m]);/* - 0.5*(lN[m].node[9]+lN[m].node[10]+lN[m].node[18]) - 0.25*(lN[m].node[22]+lN[m].node[23]+lN[m].node[26]) - 0.125*lN[m].node[27] ; */
	 lN[m].node[1] = 0.125 * (1.0-eta1[m]) * (1.0-eta2[m]) * (1.0-eta3[m]);/* - 0.5*(lN[m].node[9]+lN[m].node[12]+lN[m].node[17]) - 0.25*(lN[m].node[21]+lN[m].node[22]+lN[m].node[26]) - 0.125*lN[m].node[27] ; */
      }
  }
  }
/* #endif */

  return;
}

void v_x_shape_fn(
  struct All_variables *E,
  int el,
  standard_precision lNx[4][ELNMAX+1],
  standard_precision eta1,
  standard_precision eta2,
  standard_precision eta3,
  int level
)
      /* For further details on this implementation see BATHE's book on page 200 */
{
  const int dims = E->mesh.nsd;
  int i,j;

  if(2==dims) {

    for(j=1;j<=2;j++) for(i=1;i<=9;i++) lNx[j][i] = 0.0 ;
    switch(E->ELEM[level][el].type_v) {
    case NINE_NODES_QUAD:
      lNx[1][9] = -2.*eta1 * (1-eta2*eta2) ;
      lNx[2][9] = (1-eta1*eta1) * (-2.*eta2) ;
    case EIGHT_NODES_QUAD:
      lNx[1][8] = -eta1 * (1-eta2) - 0.5*lNx[1][9] ;
      lNx[2][8] = -0.5 * (1-eta1*eta1) - 0.5*lNx[2][9] ;
    case SEVEN_NODES_QUAD:
      lNx[1][6] = -eta1 * (1+eta2) - 0.5*lNx[1][9] ;
      lNx[2][6] = 0.5 * (1-eta1*eta1) - 0.5*lNx[2][9] ;
    case SIX_NODES_QUAD:
      lNx[1][7] = 0.5 * (1-eta2*eta2) - 0.5*lNx[1][9] ;
      lNx[2][7] = (1+eta1) * (-eta2) - 0.5*lNx[2][9] ;
    case FIVE_NODES_QUAD:
      lNx[1][5] = -0.5 * (1-eta2*eta2) - 0.5*lNx[1][9] ;
      lNx[2][5] = (1-eta1) * (-eta2) - 0.5*lNx[2][9] ;

/*    case SIX_NODES_TRIA : */
    case FOUR_NODES_QUAD:
      lNx[1][4] = 0.25 * (1.0-eta2) - 0.5*lNx[1][6] - 0.5*lNx[1][8] - 0.25*lNx[1][9] ;
      lNx[1][3] = 0.25 * (1.0+eta2) - 0.5*lNx[1][6] - 0.5*lNx[1][7] - 0.25*lNx[1][9] ;
      lNx[1][2] = -0.25 * (1.0+eta2) - 0.5*lNx[1][5] - 0.5*lNx[1][7] - 0.25*lNx[1][9] ;
      lNx[1][1] = -0.25 * (1.0-eta2) - 0.5*lNx[1][5] - 0.5*lNx[1][8] - 0.25*lNx[1][9] ;
      lNx[2][4] = -0.25 * (1.0+eta1) - 0.5*lNx[2][6] - 0.5*lNx[2][8] - 0.25*lNx[2][9] ;
      lNx[2][3] = 0.25 * (1.0+eta1) - 0.5*lNx[2][6] - 0.5*lNx[2][7] - 0.25*lNx[2][9] ;
      lNx[2][2] = 0.25 * (1.0-eta1) - 0.5*lNx[2][5] - 0.5*lNx[2][7] - 0.25*lNx[2][9] ;
      lNx[2][1] = -0.25 * (1.0-eta1) - 0.5*lNx[2][5] - 0.5*lNx[2][8] - 0.25*lNx[2][9] ;
  }
  }
  else {
    for(j=1;j<=3;j++) for(i=1;i<=27;i++) lNx[j][i] = 0.0 ;
    switch(E->ELEM[level][el].type_v) {
    case TWENTY_SEVEN_NODES_CUBIC:
     lNx[1][27] = -2.*eta1 * (1-eta2*eta2) * (1-eta3*eta3) ;
     lNx[2][27] = (1-eta1*eta1) * (-2.*eta2) * (1-eta3*eta3) ;
     lNx[3][27] = (1-eta1*eta1) * (1-eta2*eta2) * (-2.*eta3) ;
    case TWENTY_SIX_NODES_CUBIC:
     lNx[1][26] = -eta1 * (1-eta2*eta2) * (1-eta3) - 0.5*lNx[1][27] ;
     lNx[1][25] = -eta1 * (1-eta2*eta2) * (1+eta3) - 0.5*lNx[1][27] ;
     lNx[1][24] = 0.5 * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][23] = -eta1 * (1+eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][22] = -0.5 * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][21] = -eta1 * (1-eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;

     lNx[2][26] = (1-eta1*eta1) * (-eta2) * (1-eta3) - 0.5*lNx[2][27] ;
     lNx[2][25] = (1-eta1*eta1) * (-eta2) * (1+eta3) - 0.5*lNx[2][27] ;
     lNx[2][24] = (1+eta1) * (-eta2) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][23] = 0.5 * (1-eta1*eta1) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][22] = (1-eta1) * (-eta2) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][21] = -0.5 * (1-eta1*eta1) * (1-eta3*eta3) - 0.5*lNx[2][27] ;

     lNx[3][26] = -0.5 * (1-eta1*eta1) * (1-eta2*eta2) - 0.5*lNx[3][27] ;
     lNx[3][25] = 0.5 * (1-eta1*eta1) * (1-eta2*eta2) - 0.5*lNx[3][27] ;
     lNx[3][24] = (1+eta1) * (1-eta2*eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][23] = (1-eta1*eta1) * (1+eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][22] = (1-eta1) * (1-eta2*eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][21] = (1-eta1*eta1) * (1-eta2) * (-eta3) - 0.5*lNx[3][27] ;
    case TWENTY_NODES_CUBIC:
     lNx[1][20] = 0.25*(1-eta2)*(1-eta3*eta3)-0.5*(lNx[1][21]+lNx[1][24]) - 0.25*lNx[1][27] ;
     lNx[1][19] = 0.25*(1+eta2)*(1-eta3*eta3)- 0.5*(lNx[1][23]+lNx[1][24]) - 0.25*lNx[1][27] ;
     lNx[1][18] = -0.25*(1+eta2)*(1-eta3*eta3)-0.5*(lNx[1][22]+lNx[1][23]) - 0.25*lNx[1][27] ;
     lNx[1][17] = -0.25*(1-eta2)*(1-eta3*eta3)-0.5*(lNx[1][21]+lNx[1][22]) - 0.25*lNx[1][27] ;
     lNx[1][16] = 0.5*(-eta1)*(1-eta2)*(1+eta3)-0.5*(lNx[1][21]+lNx[1][25])-0.25*lNx[1][27] ;
     lNx[1][15] = 0.25*(1-eta2*eta2)*(1+eta3)-0.5*(lNx[1][24]+lNx[1][25]) - 0.25*lNx[1][27] ;
     lNx[1][14] = 0.5*(-eta1)*(1+eta2)*(1+eta3)-0.5*(lNx[1][23]+lNx[1][25])- 0.25*lNx[1][27] ;
     lNx[1][13] = -0.25*(1-eta2*eta2)*(1+eta3)-0.5*(lNx[1][22]+lNx[1][25]) - 0.25*lNx[1][27] ;
     lNx[1][12] = 0.5*(-eta1)*(1-eta2)*(1-eta3)-0.5*(lNx[1][21]+lNx[1][26])- 0.25*lNx[1][27] ;
     lNx[1][11] = 0.25*(1-eta2*eta2)*(1-eta3)-0.5*(lNx[1][24]+lNx[1][26]) - 0.25*lNx[1][27] ;
     lNx[1][10] = 0.5*(-eta1)*(1+eta2)*(1-eta3)-0.5*(lNx[1][23]+lNx[1][26])- 0.25*lNx[1][27] ;
     lNx[1][9] = -0.25*(1-eta2*eta2)*(1-eta3)-0.5*(lNx[1][22]+lNx[1][26]) - 0.25*lNx[1][27] ;

     lNx[2][20] = -0.25*(1+eta1)*(1-eta3*eta3)-0.5*(lNx[2][21]+lNx[2][24])-0.25*lNx[2][27] ;
     lNx[2][19] = 0.25*(1+eta1)*(1-eta3*eta3)-0.5*(lNx[2][23]+lNx[2][24])-0.25*lNx[2][27] ;
     lNx[2][18] = 0.25*(1-eta1)*(1-eta3*eta3)-0.5*(lNx[2][22]+lNx[2][23])-0.25*lNx[2][27] ;
     lNx[2][17] = -0.25*(1-eta1)*(1-eta3*eta3)-0.5*(lNx[2][21]+lNx[2][22])-0.25*lNx[2][27] ;
     lNx[2][16] = -0.25*(1-eta1*eta1)*(1+eta3)-0.5*(lNx[2][21]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][15] = 0.5*(1+eta1)*(-eta2)*(1+eta3)-0.5*(lNx[2][24]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][14] = 0.25*(1-eta1*eta1)*(1+eta3)-0.5*(lNx[2][23]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][13] = 0.5*(1-eta1)*(-eta2)*(1+eta3)-0.5*(lNx[2][22]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][12] = -0.25*(1-eta1*eta1)*(1-eta3)-0.5*(lNx[2][21]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][11] = 0.5*(1+eta1)*(-eta2)*(1-eta3)-0.5*(lNx[2][24]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][10] = 0.25*(1-eta1*eta1)*(1-eta3)-0.5*(lNx[2][23]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][9] = 0.5*(1-eta1)*(-eta2)*(1-eta3)-0.5*(lNx[2][22]+lNx[2][26])-0.25*lNx[2][27] ;

     lNx[3][20] = 0.5*(1+eta1)*(1-eta2)*(-eta3)-0.5*(lNx[3][21]+lNx[3][24])-0.25*lNx[3][27] ;
     lNx[3][19] = 0.5*(1+eta1)*(1+eta2)*(-eta3)-0.5*(lNx[3][23]+lNx[3][24])-0.25*lNx[3][27] ;
     lNx[3][18] = 0.5*(1-eta1)*(1+eta2)*(-eta3)-0.5*(lNx[3][22]+lNx[3][23])-0.25*lNx[3][27] ;
     lNx[3][17] = 0.5*(1-eta1)*(1-eta2)*(-eta3)-0.5*(lNx[3][21]+lNx[3][22])-0.25*lNx[3][27] ;
     lNx[3][16] = 0.25*(1-eta1*eta1)*(1-eta2)-0.5*(lNx[3][21]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][15] = 0.25*(1+eta1)*(1-eta2*eta2)-0.5*(lNx[3][24]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][14] = 0.25*(1-eta1*eta1)*(1+eta2)-0.5*(lNx[3][23]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][13] = 0.25*(1-eta1)*(1-eta2*eta2)-0.5*(lNx[3][22]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][12] = -0.25*(1-eta1*eta1)*(1-eta2)-0.5*(lNx[3][21]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][11] = -0.25*(1+eta1)*(1-eta2*eta2)-0.5*(lNx[3][24]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][10] = -0.25*(1-eta1*eta1)*(1+eta2)-0.5*(lNx[3][23]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][9] = -0.25*(1-eta1)*(1-eta2*eta2)-0.5*(lNx[3][22]+lNx[3][26])-0.25*lNx[3][27] ;

    case EIGHT_NODES_CUBIC:  /*RAA, comment out higher order node contributions as they are unnecessary*/
     lNx[1][8] = 0.125 * (1.0-eta2) * (1.0+eta3); /* - 0.5*(lNx[1][15]+lNx[1][16]+lNx[1][20]) - 0.25*(lNx[1][21]+lNx[1][24]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][7] = 0.125  * (1.0+eta2) * (1.0+eta3); /* - 0.5*(lNx[1][14]+lNx[1][15]+lNx[1][19]) - 0.25*(lNx[1][23]+lNx[1][24]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][6] = -0.125 * (1.0+eta2) * (1.0+eta3); /* - 0.5*(lNx[1][13]+lNx[1][14]+lNx[1][18]) - 0.25*(lNx[1][22]+lNx[1][23]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][5] = -0.125 * (1.0-eta2) * (1.0+eta3); /* - 0.5*(lNx[1][13]+lNx[1][16]+lNx[1][17]) - 0.25*(lNx[1][21]+lNx[1][22]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][4] = 0.125 * (1.0-eta2) * (1.0-eta3); /* - 0.5*(lNx[1][11]+lNx[1][12]+lNx[1][20]) - 0.25*(lNx[1][21]+lNx[1][24]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][3] = 0.125 * (1.0+eta2) * (1.0-eta3); /* - 0.5*(lNx[1][10]+lNx[1][11]+lNx[1][19]) - 0.25*(lNx[1][23]+lNx[1][24]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][2] = -0.125 * (1.0+eta2) * (1.0-eta3); /* - 0.5*(lNx[1][9]+lNx[1][10]+lNx[1][18]) - 0.25*(lNx[1][22]+lNx[1][23]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][1] = -0.125 * (1.0-eta2) * (1.0-eta3); /* - 0.5*(lNx[1][9]+lNx[1][12]+lNx[1][17]) - 0.25*(lNx[1][21]+lNx[1][22]+lNx[1][26]) - 0.125*lNx[1][27] ; */

     lNx[2][8] = -0.125 * (1.0+eta1) * (1.0+eta3); /* - 0.5*(lNx[2][15]+lNx[2][16]+lNx[2][20]) - 0.25*(lNx[2][21]+lNx[2][24]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][7] = 0.125 * (1.0+eta1) * (1.0+eta3); /* - 0.5*(lNx[2][14]+lNx[2][15]+lNx[2][19]) - 0.25*(lNx[2][23]+lNx[2][24]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][6] = 0.125 * (1.0-eta1) * (1.0+eta3); /* - 0.5*(lNx[2][13]+lNx[2][14]+lNx[2][18]) - 0.25*(lNx[2][22]+lNx[2][23]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][5] = -0.125 * (1.0-eta1) * (1.0+eta3); /* - 0.5*(lNx[2][13]+lNx[2][16]+lNx[2][17]) - 0.25*(lNx[2][21]+lNx[2][22]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][4] = -0.125 * (1.0+eta1) * (1.0-eta3); /* - 0.5*(lNx[2][11]+lNx[2][12]+lNx[2][20]) - 0.25*(lNx[2][21]+lNx[2][24]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][3] = 0.125 * (1.0+eta1) * (1.0-eta3); /* - 0.5*(lNx[2][10]+lNx[2][11]+lNx[2][19]) - 0.25*(lNx[2][23]+lNx[2][24]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][2] = 0.125 * (1.0-eta1) * (1.0-eta3); /* - 0.5*(lNx[2][9]+lNx[2][10]+lNx[2][18]) - 0.25*(lNx[2][22]+lNx[2][23]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][1] = -0.125 * (1.0-eta1) * (1.0-eta3); /* - 0.5*(lNx[2][9]+lNx[2][12]+lNx[2][17]) - 0.25*(lNx[2][21]+lNx[2][22]+lNx[2][26]) - 0.125*lNx[2][27] ; */

     lNx[3][8] = 0.125 * (1.0+eta1) * (1.0-eta2); /*- 0.5*(lNx[3][15]+lNx[3][16]+lNx[3][20]) - 0.25*(lNx[3][21]+lNx[3][24]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][7] = 0.125 * (1.0+eta1) * (1.0+eta2); /*- 0.5*(lNx[3][14]+lNx[3][15]+lNx[3][19]) - 0.25*(lNx[3][23]+lNx[3][24]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][6] = 0.125 * (1.0-eta1) * (1.0+eta2); /*- 0.5*(lNx[3][13]+lNx[3][14]+lNx[3][18]) - 0.25*(lNx[3][22]+lNx[3][23]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][5] = 0.125 * (1.0-eta1) * (1.0-eta2); /*- 0.5*(lNx[3][13]+lNx[3][16]+lNx[3][17]) - 0.25*(lNx[3][21]+lNx[3][22]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][4] = -0.125 * (1.0+eta1) * (1.0-eta2); /*- 0.5*(lNx[3][11]+lNx[3][12]+lNx[3][20]) - 0.25*(lNx[3][21]+lNx[3][24]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][3] = -0.125 * (1.0+eta1) * (1.0+eta2); /*- 0.5*(lNx[3][10]+lNx[3][11]+lNx[3][19]) - 0.25*(lNx[3][23]+lNx[3][24]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][2] = -0.125 * (1.0-eta1) * (1.0+eta2); /*- 0.5*(lNx[3][9]+lNx[3][10]+lNx[3][18]) - 0.25*(lNx[3][22]+lNx[3][23]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][1] = -0.125 * (1.0-eta1) * (1.0-eta2); /*- 0.5*(lNx[3][9]+lNx[3][12]+lNx[3][17]) - 0.25*(lNx[3][21]+lNx[3][22]+lNx[3][26]) - 0.125*lNx[3][27] ; */
  }
  }
  return;
}

void p_shape_fn(
  struct All_variables *E,
  int el,
  standard_precision *lN,
  standard_precision eta1,
  standard_precision eta2,
  standard_precision eta3,
  int level
)
      /* For further details on this implementation see BATHE's book on page 200 */
{
  const int dims = E->mesh.nsd;
  int i;

  if(2==dims) {

    for(i=1;i<=9;i++) lN[i] = 0.0 ;
    switch(E->ELEM[level][el].type_p) {
    case FOUR_NODES_QUAD:
    case ONE_NODE_QUAD:
     lN[1] = 1.0 ;
  }
  }
  else {
    for(i=1;i<=27;i++) lN[i] = 0.0 ;
    switch(E->ELEM[level][el].type_p) {
    case ONE_NODE_QUAD:
     lN[1] = 1.0 ;
  }
  }
  return;
}

void p_x_shape_fn(
  struct All_variables *E,
  int el,
  standard_precision lNx[4][ELNMAX+1],
  standard_precision eta1,
  standard_precision eta2,
  standard_precision eta3,
  int level
)
      /* For further details on this implementation see BATHE's book on page 200 */
{
  const int dims = E->mesh.nsd;
  int i,j;

  if(2==dims) {

    for(j=1;j<=2;j++) for(i=1;i<=9;i++) lNx[j][i] = 0.0 ;
    switch(E->ELEM[level][el].type_p) {
    case NINE_NODES_QUAD:
      lNx[1][9] = -2.*eta1 * (1-eta2*eta2) ;
      lNx[2][9] = (1-eta1*eta1) * (-2.*eta2) ;
    case EIGHT_NODES_QUAD:
      lNx[1][8] = -eta1 * (1-eta2) - 0.5*lNx[1][9] ;
      lNx[2][8] = -0.5 * (1-eta1*eta1) - 0.5*lNx[2][9] ;
    case SEVEN_NODES_QUAD:
      lNx[1][6] = -eta1 * (1+eta2) - 0.5*lNx[1][9] ;
      lNx[2][6] = 0.5 * (1-eta1*eta1) - 0.5*lNx[2][9] ;
    case SIX_NODES_QUAD:
      lNx[1][7] = 0.5 * (1-eta2*eta2) - 0.5*lNx[1][9] ;
      lNx[2][7] = (1+eta1) * (-eta2) - 0.5*lNx[2][9] ;
    case FIVE_NODES_QUAD:
      lNx[1][5] = -0.5 * (1-eta2*eta2) - 0.5*lNx[1][9] ;
      lNx[2][5] = (1-eta1) * (-eta2) - 0.5*lNx[2][9] ;

/*    case SIX_NODES_TRIA : */
    case FOUR_NODES_QUAD:
      lNx[1][4] = 0.25 * (1.0-eta2) - 0.5*lNx[1][6] - 0.5*lNx[1][8] - 0.25*lNx[1][9] ;
      lNx[1][3] = 0.25 * (1.0+eta2) - 0.5*lNx[1][6] - 0.5*lNx[1][7] - 0.25*lNx[1][9] ;
      lNx[1][2] = -0.25 * (1.0+eta2) - 0.5*lNx[1][5] - 0.5*lNx[1][7] - 0.25*lNx[1][9] ;
      lNx[1][1] = -0.25 * (1.0-eta2) - 0.5*lNx[1][5] - 0.5*lNx[1][8] - 0.25*lNx[1][9] ;
      lNx[2][4] = -0.25 * (1.0+eta1) - 0.5*lNx[2][6] - 0.5*lNx[2][8] - 0.25*lNx[2][9] ;
      lNx[2][3] = 0.25 * (1.0+eta1) - 0.5*lNx[2][6] - 0.5*lNx[2][7] - 0.25*lNx[2][9] ;
      lNx[2][2] = 0.25 * (1.0-eta1) - 0.5*lNx[2][5] - 0.5*lNx[2][7] - 0.25*lNx[2][9] ;
      lNx[2][1] = -0.25 * (1.0-eta1) - 0.5*lNx[2][5] - 0.5*lNx[2][8] - 0.25*lNx[2][9] ;
  }
  }
  else {
    for(j=1;j<=3;j++) for(i=1;i<=27;i++) lNx[j][i] = 0.0 ;
    switch(E->ELEM[level][el].type_p) {
    case TWENTY_SEVEN_NODES_CUBIC:
     lNx[1][27] = -2.*eta1 * (1-eta2*eta2) * (1-eta3*eta3) ;
     lNx[2][27] = (1-eta1*eta1) * (-2.*eta2) * (1-eta3*eta3) ;
     lNx[3][27] = (1-eta1*eta1) * (1-eta2*eta2) * (-2.*eta3) ;
    case TWENTY_SIX_NODES_CUBIC:
     lNx[1][26] = -eta1 * (1-eta2*eta2) * (1-eta3) - 0.5*lNx[1][27] ;
     lNx[1][25] = -eta1 * (1-eta2*eta2) * (1+eta3) - 0.5*lNx[1][27] ;
     lNx[1][24] = 0.5 * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][23] = -eta1 * (1+eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][22] = -0.5 * (1-eta2*eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;
     lNx[1][21] = -eta1 * (1-eta2) * (1-eta3*eta3) - 0.5*lNx[1][27] ;

     lNx[2][26] = (1-eta1*eta1) * (-eta2) * (1-eta3) - 0.5*lNx[2][27] ;
     lNx[2][25] = (1-eta1*eta1) * (-eta2) * (1+eta3) - 0.5*lNx[2][27] ;
     lNx[2][24] = (1+eta1) * (-eta2) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][23] = 0.5 * (1-eta1*eta1) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][22] = (1-eta1) * (-eta2) * (1-eta3*eta3) - 0.5*lNx[2][27] ;
     lNx[2][21] = -0.5 * (1-eta1*eta1) * (1-eta3*eta3) - 0.5*lNx[2][27] ;

     lNx[3][26] = -0.5 * (1-eta1*eta1) * (1-eta2*eta2) - 0.5*lNx[3][27] ;
     lNx[3][25] = 0.5 * (1-eta1*eta1) * (1-eta2*eta2) - 0.5*lNx[3][27] ;
     lNx[3][24] = (1+eta1) * (1-eta2*eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][23] = (1-eta1*eta1) * (1+eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][22] = (1-eta1) * (1-eta2*eta2) * (-eta3) - 0.5*lNx[3][27] ;
     lNx[3][21] = (1-eta1*eta1) * (1-eta2) * (-eta3) - 0.5*lNx[3][27] ;
    case TWENTY_NODES_CUBIC:
     lNx[1][20] = 0.25*(1-eta2)*(1-eta3*eta3)-0.5*(lNx[1][21]+lNx[1][24]) - 0.25*lNx[1][27] ;
     lNx[1][19] = 0.25*(1+eta2)*(1-eta3*eta3)- 0.5*(lNx[1][23]+lNx[1][24]) - 0.25*lNx[1][27] ;
     lNx[1][18] = -0.25*(1+eta2)*(1-eta3*eta3)-0.5*(lNx[1][22]+lNx[1][23]) - 0.25*lNx[1][27] ;
     lNx[1][17] = -0.25*(1-eta2)*(1-eta3*eta3)-0.5*(lNx[1][21]+lNx[1][22]) - 0.25*lNx[1][27] ;
     lNx[1][16] = 0.5*(-eta1)*(1-eta2)*(1+eta3)-0.5*(lNx[1][21]+lNx[1][25])-0.25*lNx[1][27] ;
     lNx[1][15] = 0.25*(1-eta2*eta2)*(1+eta3)-0.5*(lNx[1][24]+lNx[1][25]) - 0.25*lNx[1][27] ;
     lNx[1][14] = 0.5*(-eta1)*(1+eta2)*(1+eta3)-0.5*(lNx[1][23]+lNx[1][25])- 0.25*lNx[1][27] ;
     lNx[1][13] = -0.25*(1-eta2*eta2)*(1+eta3)-0.5*(lNx[1][22]+lNx[1][25]) - 0.25*lNx[1][27] ;
     lNx[1][12] = 0.5*(-eta1)*(1-eta2)*(1-eta3)-0.5*(lNx[1][21]+lNx[1][26])- 0.25*lNx[1][27] ;
     lNx[1][11] = 0.25*(1-eta2*eta2)*(1-eta3)-0.5*(lNx[1][24]+lNx[1][26]) - 0.25*lNx[1][27] ;
     lNx[1][10] = 0.5*(-eta1)*(1+eta2)*(1-eta3)-0.5*(lNx[1][23]+lNx[1][26])- 0.25*lNx[1][27] ;
     lNx[1][9] = -0.25*(1-eta2*eta2)*(1-eta3)-0.5*(lNx[1][22]+lNx[1][26]) - 0.25*lNx[1][27] ;

     lNx[2][20] = -0.25*(1+eta1)*(1-eta3*eta3)-0.5*(lNx[2][21]+lNx[2][24])-0.25*lNx[2][27] ;
     lNx[2][19] = 0.25*(1+eta1)*(1-eta3*eta3)-0.5*(lNx[2][23]+lNx[2][24])-0.25*lNx[2][27] ;
     lNx[2][18] = 0.25*(1-eta1)*(1-eta3*eta3)-0.5*(lNx[2][22]+lNx[2][23])-0.25*lNx[2][27] ;
     lNx[2][17] = -0.25*(1-eta1)*(1-eta3*eta3)-0.5*(lNx[2][21]+lNx[2][22])-0.25*lNx[2][27] ;
     lNx[2][16] = -0.25*(1-eta1*eta1)*(1+eta3)-0.5*(lNx[2][21]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][15] = 0.5*(1+eta1)*(-eta2)*(1+eta3)-0.5*(lNx[2][24]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][14] = 0.25*(1-eta1*eta1)*(1+eta3)-0.5*(lNx[2][23]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][13] = 0.5*(1-eta1)*(-eta2)*(1+eta3)-0.5*(lNx[2][22]+lNx[2][25])-0.25*lNx[2][27] ;
     lNx[2][12] = -0.25*(1-eta1*eta1)*(1-eta3)-0.5*(lNx[2][21]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][11] = 0.5*(1+eta1)*(-eta2)*(1-eta3)-0.5*(lNx[2][24]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][10] = 0.25*(1-eta1*eta1)*(1-eta3)-0.5*(lNx[2][23]+lNx[2][26])-0.25*lNx[2][27] ;
     lNx[2][9] = 0.5*(1-eta1)*(-eta2)*(1-eta3)-0.5*(lNx[2][22]+lNx[2][26])-0.25*lNx[2][27] ;

     lNx[3][20] = 0.5*(1+eta1)*(1-eta2)*(-eta3)-0.5*(lNx[3][21]+lNx[3][24])-0.25*lNx[3][27] ;
     lNx[3][19] = 0.5*(1+eta1)*(1+eta2)*(-eta3)-0.5*(lNx[3][23]+lNx[3][24])-0.25*lNx[3][27] ;
     lNx[3][18] = 0.5*(1-eta1)*(1+eta2)*(-eta3)-0.5*(lNx[3][22]+lNx[3][23])-0.25*lNx[3][27] ;
     lNx[3][17] = 0.5*(1-eta1)*(1-eta2)*(-eta3)-0.5*(lNx[3][21]+lNx[3][22])-0.25*lNx[3][27] ;
     lNx[3][16] = 0.25*(1-eta1*eta1)*(1-eta2)-0.5*(lNx[3][21]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][15] = 0.25*(1+eta1)*(1-eta2*eta2)-0.5*(lNx[3][24]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][14] = 0.25*(1-eta1*eta1)*(1+eta2)-0.5*(lNx[3][23]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][13] = 0.25*(1-eta1)*(1-eta2*eta2)-0.5*(lNx[3][22]+lNx[3][25])-0.25*lNx[3][27] ;
     lNx[3][12] = -0.25*(1-eta1*eta1)*(1-eta2)-0.5*(lNx[3][21]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][11] = -0.25*(1+eta1)*(1-eta2*eta2)-0.5*(lNx[3][24]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][10] = -0.25*(1-eta1*eta1)*(1+eta2)-0.5*(lNx[3][23]+lNx[3][26])-0.25*lNx[3][27] ;
     lNx[3][9] = -0.25*(1-eta1)*(1-eta2*eta2)-0.5*(lNx[3][22]+lNx[3][26])-0.25*lNx[3][27] ;

    case EIGHT_NODES_CUBIC:
     lNx[1][8] = 0.125 * (1.0-eta2) * (1.0+eta3); /*- 0.5*(lNx[1][15]+lNx[1][16]+lNx[1][20]) - 0.25*(lNx[1][21]+lNx[1][24]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][7] = 0.125  * (1.0+eta2) * (1.0+eta3); /*- 0.5*(lNx[1][14]+lNx[1][15]+lNx[1][19]) - 0.25*(lNx[1][23]+lNx[1][24]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][6] = -0.125 * (1.0+eta2) * (1.0+eta3); /*- 0.5*(lNx[1][13]+lNx[1][14]+lNx[1][18]) - 0.25*(lNx[1][22]+lNx[1][23]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][5] = -0.125 * (1.0-eta2) * (1.0+eta3); /*- 0.5*(lNx[1][13]+lNx[1][16]+lNx[1][17]) - 0.25*(lNx[1][21]+lNx[1][22]+lNx[1][25]) - 0.125*lNx[1][27] ; */
     lNx[1][4] = 0.125 * (1.0-eta2) * (1.0-eta3); /*- 0.5*(lNx[1][11]+lNx[1][12]+lNx[1][20]) - 0.25*(lNx[1][21]+lNx[1][24]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][3] = 0.125 * (1.0+eta2) * (1.0-eta3); /*- 0.5*(lNx[1][10]+lNx[1][11]+lNx[1][19]) - 0.25*(lNx[1][23]+lNx[1][24]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][2] = -0.125 * (1.0+eta2) * (1.0-eta3); /*- 0.5*(lNx[1][9]+lNx[1][10]+lNx[1][18]) - 0.25*(lNx[1][22]+lNx[1][23]+lNx[1][26]) - 0.125*lNx[1][27] ; */
     lNx[1][1] = -0.125 * (1.0-eta2) * (1.0-eta3); /*- 0.5*(lNx[1][9]+lNx[1][12]+lNx[1][17]) - 0.25*(lNx[1][21]+lNx[1][22]+lNx[1][26]) - 0.125*lNx[1][27] ; */

     lNx[2][8] = -0.125 * (1.0+eta1) * (1.0+eta3); /*- 0.5*(lNx[2][15]+lNx[2][16]+lNx[2][20]) - 0.25*(lNx[2][21]+lNx[2][24]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][7] = 0.125 * (1.0+eta1) * (1.0+eta3); /*- 0.5*(lNx[2][14]+lNx[2][15]+lNx[2][19]) - 0.25*(lNx[2][23]+lNx[2][24]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][6] = 0.125 * (1.0-eta1) * (1.0+eta3); /*- 0.5*(lNx[2][13]+lNx[2][14]+lNx[2][18]) - 0.25*(lNx[2][22]+lNx[2][23]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][5] = -0.125 * (1.0-eta1) * (1.0+eta3); /*- 0.5*(lNx[2][13]+lNx[2][16]+lNx[2][17]) - 0.25*(lNx[2][21]+lNx[2][22]+lNx[2][25]) - 0.125*lNx[2][27] ; */
     lNx[2][4] = -0.125 * (1.0+eta1) * (1.0-eta3); /*- 0.5*(lNx[2][11]+lNx[2][12]+lNx[2][20]) - 0.25*(lNx[2][21]+lNx[2][24]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][3] = 0.125 * (1.0+eta1) * (1.0-eta3); /*- 0.5*(lNx[2][10]+lNx[2][11]+lNx[2][19]) - 0.25*(lNx[2][23]+lNx[2][24]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][2] = 0.125 * (1.0-eta1) * (1.0-eta3); /*- 0.5*(lNx[2][9]+lNx[2][10]+lNx[2][18]) - 0.25*(lNx[2][22]+lNx[2][23]+lNx[2][26]) - 0.125*lNx[2][27] ; */
     lNx[2][1] = -0.125 * (1.0-eta1) * (1.0-eta3); /*- 0.5*(lNx[2][9]+lNx[2][12]+lNx[2][17]) - 0.25*(lNx[2][21]+lNx[2][22]+lNx[2][26]) - 0.125*lNx[2][27] ; */

     lNx[3][8] = 0.125 * (1.0+eta1) * (1.0-eta2); /*- 0.5*(lNx[3][15]+lNx[3][16]+lNx[3][20]) - 0.25*(lNx[3][21]+lNx[3][24]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][7] = 0.125 * (1.0+eta1) * (1.0+eta2); /*- 0.5*(lNx[3][14]+lNx[3][15]+lNx[3][19]) - 0.25*(lNx[3][23]+lNx[3][24]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][6] = 0.125 * (1.0-eta1) * (1.0+eta2); /*- 0.5*(lNx[3][13]+lNx[3][14]+lNx[3][18]) - 0.25*(lNx[3][22]+lNx[3][23]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][5] = 0.125 * (1.0-eta1) * (1.0-eta2); /*- 0.5*(lNx[3][13]+lNx[3][16]+lNx[3][17]) - 0.25*(lNx[3][21]+lNx[3][22]+lNx[3][25]) - 0.125*lNx[3][27] ; */
     lNx[3][4] = -0.125 * (1.0+eta1) * (1.0-eta2); /*- 0.5*(lNx[3][11]+lNx[3][12]+lNx[3][20]) - 0.25*(lNx[3][21]+lNx[3][24]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][3] = -0.125 * (1.0+eta1) * (1.0+eta2); /*- 0.5*(lNx[3][10]+lNx[3][11]+lNx[3][19]) - 0.25*(lNx[3][23]+lNx[3][24]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][2] = -0.125 * (1.0-eta1) * (1.0+eta2); /*- 0.5*(lNx[3][9]+lNx[3][10]+lNx[3][18]) - 0.25*(lNx[3][22]+lNx[3][23]+lNx[3][26]) - 0.125*lNx[3][27] ; */
     lNx[3][1] = -0.125 * (1.0-eta1) * (1.0-eta2); /*- 0.5*(lNx[3][9]+lNx[3][12]+lNx[3][17]) - 0.25*(lNx[3][21]+lNx[3][22]+lNx[3][26]) - 0.125*lNx[3][27] ; */
  }
  }
  return;
}

