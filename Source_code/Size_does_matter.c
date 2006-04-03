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



 /*   This is where the scaling functions and grid related things are kept.
	Louis Moresi aka LUIGI   6.xii.1989                */

#include <math.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"

/*===========================================================================================
  Function to give the global shape function from the local: Accelerated for  ORTHOGONAL MESH 
  ===========================================================================================      */
	
void get_global_shape_fn(
     struct All_variables *E,
     int el,
     struct Shape_function *GN,
     struct Shape_function_dx *GNx,
     struct Shape_function_dA *dOmega,
     int pressure,
     int level
)
{
  int i,j,k,d,e;
  higher_precision scale1,scale2,scale3;
  higher_precision area;
  higher_precision jacobian;
  higher_precision determinant();
  higher_precision cofactor();
 
  higher_precision dxda[4][4],cof[4][4];
 
  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const int ppts=ppoints[dims];
  
  if  (E->control.AXI || E->control.ORTHO)  /*accelerated for simple elements */  {
      scale1 = E->ECO[level][el].recip_size[1];
      scale2 = E->ECO[level][el].recip_size[2];
      if(3==dims)
	scale3 = E->ECO[level][el].recip_size[3];
        area = E->ECO[level][el].area;
     
      if(pressure < 2) /* velocity point info is required */	{
	for(i=0;i<GNVI;i++)    {
	  GNx->vpt[GNVXSHORT(0,i)] = E->Nx.vpt[GNVXSHORT(0,i)] * scale1;
	  GNx->vpt[GNVXSHORT(1,i)] = E->Nx.vpt[GNVXSHORT(1,i)] * scale2;
	  if(3==dims)
	    GNx->vpt[GNVXSHORT(2,i)] = E->Nx.vpt[GNVXSHORT(2,i)] * scale3;
	}
	  
	for(i=1;i<=vpts;i++)
	  dOmega->vpt[i] = area;
	  
	if (E->control.AXI) {
	  for(i=1;i<=vpts;i++)
	    dOmega->vpt[i] *= 2.0*M_PI*(g_point[i].x[0]*E->ECO[level][el].size[1]+E->ECO[level][el].centre[1]);
	}
      }
      
      if(pressure > 0) /* pressure point info is required */ {
	for(i=0;i<GNPI;i++) {
	  GNx->ppt[GNPXSHORT(0,i)] = E->Nx.ppt[GNPXSHORT(0,i)] * scale1;
	  GNx->ppt[GNPXSHORT(1,i)] = E->Nx.ppt[GNPXSHORT(1,i)] * scale2;
	  if(3==dims)
	    GNx->ppt[GNPXSHORT(2,i)] = E->Nx.ppt[GNPXSHORT(2,i)] * scale3;
	}
	  
	for(i=1;i<=ppts;i++)
	  dOmega->ppt[i] = area; 
	  
	if (E->control.AXI)
	  for(i=1;i<=ppts;i++)
	    dOmega->ppt[i] *= 2 * M_PI * (p_point[i].x[0] * E->ECO[level][el].size[1] + E->ECO[level][el].centre[1]); 
      }
  }
  
  else /* not axisymmetry/orthogonal */  { 

    if(pressure < 2) {
      for(k=1;k<=vpts;k++) /* all of the vpoints */	    {
 
	for(d=1;d<=dims;d++)
	  for(e=1;e<=dims;e++)
	    dxda[d][e]=0.0;
	   
	for(i=1;i<=ends;i++)
	  for(d=1;d<=dims;d++)
	    for(e=1;e<=dims;e++)
	      dxda[d][e] += E->X[level][e][E->IEN[level][el].node[i]] 
		* E->Nx.vpt[GNVXINDEX(d-1,i,k)];
/* NOTE:This is correct, Hughes book is in  error and one needs to exchange (d<->e) wrt to his notation. */
	jacobian = determinant(dxda,dims);  
	dOmega->vpt[k] = jacobian;

	for(d=1;d<=dims;d++)
	  for(e=1;e<=dims;e++)
	    cof[d][e]=cofactor(dxda,d,e,dims);
 
	for(j=1;j<=ends;j++)
	  for(d=1;d<=dims;d++) {
	    GNx->vpt[GNVXINDEX(d-1,j,k)] = 0.0;
	    for(e=1;e<=dims;e++)
	      GNx->vpt[GNVXINDEX(d-1,j,k)] += 
		E->Nx.vpt[GNVXINDEX(e-1,j,k)] *cof[e][d];

	    GNx->vpt[GNVXINDEX(d-1,j,k)] /= jacobian;
	  }
      }
    }

    if(pressure > 0) {
      for(k=1;k<=ppts;k++) /* all of the ppoints */  {
	for(d=1;d<=dims;d++)
	  for(e=1;e<=dims;e++)
	    dxda[d][e]=0.0;
	      
	for(i=1;i<=ends;i++)
	  for(d=1;d<=dims;d++)
	    for(e=1;e<=dims;e++)
	      dxda[d][e] += E->X[level][e][E->IEN[level][el].node[i]] * E->Nx.ppt[GNPXINDEX(d-1,i,k)];

	jacobian = determinant(dxda,E->mesh.nsd);     
	dOmega->ppt[k] = jacobian;
	      
	for(d=1;d<=dims;d++)
	  for(e=1;e<=dims;e++)
	    cof[d][e]=cofactor(dxda,d,e,E->mesh.nsd); 
	      
	for(j=1;j<=ends;j++)
	  for(d=1;d<=dims;d++)   {
	    GNx->ppt[GNPXINDEX(d-1,j,k)]=0.0;
	    for(e=1;e<=dims;e++)
	      GNx->ppt[GNPXINDEX(d-1,j,k)] += E->Nx.ppt[GNPXINDEX(e-1,j,k)]*cof[e][d]; 
	    GNx->ppt[GNPXINDEX(d-1,j,k)] /= jacobian;
	  }
      }
    } 
  }
  return;
}

/* Tracer jacobian */

standard_precision get_tracer_jacobian(
				       struct All_variables *E,
				       int m,
				       standard_precision eta1,
				       standard_precision eta2,
				       standard_precision eta3,
				       int level
				       )
{
  int d,e,el,i;

  higher_precision jacobian;
  higher_precision determinant();
  higher_precision cofactor();

  standard_precision lNx[4][ELNMAX+1];
  higher_precision dxda[4][4],cof[4][4];

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

  el =  E->tracer.tracer_elt[level][m];

  if ((E->control.AXI || E->control.ORTHO)) {
    jacobian = E->ECO[level][el].area;
  }
  else {
    v_x_shape_fn(E,el,lNx,eta1,eta2,eta3,level);
    for(d=1;d<=dims;d++)
      for(e=1;e<=dims;e++)
	dxda[d][e]=0.0;
	   
    for(i=1;i<=ends;i++)
      for(d=1;d<=dims;d++)
	for(e=1;e<=dims;e++)
	  dxda[d][e] += E->X[level][e][E->IEN[level][el].node[i]] 
	    * lNx[d][i];                

    /* NOTE:  This is correct, Hughes book p147 is in
       error and one needs to exchange (d<->e)   wrt to his notation. */

    jacobian = determinant(dxda,dims);  
  }
  return(jacobian);
}

void get_global_v_x_shape_fn(
			     struct All_variables *E,
			     int el,
			     standard_precision lNx[4][ELNMAX+1],
			     standard_precision *dOmega,
			     standard_precision eta1,
			     standard_precision eta2,
			     standard_precision eta3,
			     int level
			     )
{
  int i,j,k,d,e;
  higher_precision scale1,scale2,scale3;
  higher_precision area;
  higher_precision jacobian;
  higher_precision determinant();
  higher_precision cofactor();
 
  higher_precision dxda[4][4],cof[4][4];

  standard_precision llNx[4][ELNMAX+1];

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const int ppts=ppoints[dims];

  scale1 = E->ECO[level][el].recip_size[1];
  scale2 = E->ECO[level][el].recip_size[2];
  if(3==dims)
    scale3 = E->ECO[level][el].recip_size[3];
  
  if(2==dims) {
    *dOmega = 4.0 *  E->ECO[level][el].area;

    lNx[1][1] = (llNx[1][1] = -0.25 * (1.0-eta2)) * scale1;
    lNx[1][2] = (llNx[1][2] = -0.25 * (1.0+eta2)) * scale1;
    lNx[1][3] = (llNx[1][3] =  0.25 * (1.0+eta2)) * scale1;
    lNx[1][4] = (llNx[1][4] =  0.25 * (1.0-eta2)) * scale1;
    lNx[2][1] = (llNx[2][1] = -0.25 * (1.0-eta1)) * scale2;
    lNx[2][2] = (llNx[2][2] =  0.25 * (1.0-eta1)) * scale2;
    lNx[2][3] = (llNx[2][3] =  0.25 * (1.0+eta1)) * scale2;
    lNx[2][4] = (llNx[2][4] = -0.25 * (1.0+eta1)) * scale2;
  }
  else {
    *dOmega = 8.0 *  E->ECO[level][el].area;

    lNx[1][1] = (llNx[1][1] = -0.125 * (1.0-eta2) * (1.0-eta3)) * scale1;
    lNx[1][2] = (llNx[1][2] = -0.125 * (1.0+eta2) * (1.0-eta3)) * scale1;
    lNx[1][3] = (llNx[1][3] =  0.125 * (1.0+eta2) * (1.0-eta3)) * scale1;
    lNx[1][4] = (llNx[1][4] =  0.125 * (1.0-eta2) * (1.0-eta3)) * scale1;
    lNx[1][5] = (llNx[1][5] = -0.125 * (1.0-eta2) * (1.0+eta3)) * scale1;
    lNx[1][6] = (llNx[1][6] = -0.125 * (1.0+eta2) * (1.0+eta3)) * scale1;
    lNx[1][7] = (llNx[1][7] =  0.125 * (1.0+eta2) * (1.0+eta3)) * scale1;
    lNx[1][8] = (llNx[1][8] =  0.125 * (1.0-eta2) * (1.0+eta3)) * scale1;
    
    lNx[2][1] = (llNx[2][1] = -0.125 * (1.0-eta1) * (1.0-eta3)) * scale2;
    lNx[2][2] = (llNx[2][2] =  0.125 * (1.0-eta1) * (1.0-eta3)) * scale2;
    lNx[2][3] = (llNx[2][3] =  0.125 * (1.0+eta1) * (1.0-eta3)) * scale2;
    lNx[2][4] = (llNx[2][4] = -0.125 * (1.0+eta1) * (1.0-eta3)) * scale2;
    lNx[2][5] = (llNx[2][5] = -0.125 * (1.0-eta1) * (1.0+eta3)) * scale2;
    lNx[2][6] = (llNx[2][6] =  0.125 * (1.0-eta1) * (1.0+eta3)) * scale2;
    lNx[2][7] = (llNx[2][7] =  0.125 * (1.0+eta1) * (1.0+eta3)) * scale2;
    lNx[2][8] = (llNx[2][8] = -0.125 * (1.0+eta1) * (1.0+eta3)) * scale2;
  
    lNx[3][1] = (llNx[3][1] = -0.125 * (1.0-eta1) * (1.0-eta2)) * scale3;
    lNx[3][2] = (llNx[3][2] = -0.125 * (1.0-eta1) * (1.0+eta2)) * scale3;
    lNx[3][3] = (llNx[3][3] = -0.125 * (1.0+eta1) * (1.0+eta2)) * scale3;
    lNx[3][4] = (llNx[3][4] = -0.125 * (1.0+eta1) * (1.0-eta2)) * scale3;
    lNx[3][5] = (llNx[3][5] =  0.125 * (1.0-eta1) * (1.0-eta2)) * scale3;
    lNx[3][6] = (llNx[3][6] =  0.125 * (1.0-eta1) * (1.0+eta2)) * scale3;
    lNx[3][7] = (llNx[3][7] =  0.125 * (1.0+eta1) * (1.0+eta2)) * scale3;
    lNx[3][8] = (llNx[3][8] =  0.125 * (1.0+eta1) * (1.0-eta2)) * scale3;
    }
  
  
  /*RAA: check*/
/* if(3==dims) { 
      fprintf(E->fp1,"el:  %d  lNx1s:  %g   %g  %g  %g   %g   %g   %g  %g\n",el,lNx[1][1],lNx[1][2],lNx[1][3],lNx[1][4],lNx[1][5],lNx[1][6],lNx[1][7],lNx[1][8]);
      fprintf(E->fp1,"el:  %d  lNx2s:  %g   %g  %g  %g   %g   %g   %g  %g\n",el,lNx[2][1],lNx[2][2],lNx[2][3],lNx[2][4],lNx[2][5],lNx[2][6],lNx[2][7],lNx[2][8]);
      fprintf(E->fp1,"el:  %d  lNx3s:  %g   %g  %g  %g   %g   %g   %g  %g\n",el,lNx[3][1],lNx[3][2],lNx[3][3],lNx[3][4],lNx[3][5],lNx[3][6],lNx[3][7],lNx[3][8]); 
      fprintf(E->fp1,"el:  %d  scales:  %g   %g  %g\n",el,scale1,scale2,scale3);
   } 
 */ 

    /* If not distorted element, then the lNx values are correct.
       Otherwise, use the llNx copy of the local shape function derivatives
       to get lNx, the global value. */
  
    if  (E->control.AXI  || E->control.ORTHO) {
      return;
    }
 
    /* Distorted element */
    
    for(d=1;d<=dims;d++)
      for(e=1;e<=dims;e++)
	dxda[d][e]=0.0;
	   
    for(i=1;i<=ends;i++)
      for(d=1;d<=dims;d++)
	for(e=1;e<=dims;e++)
	  dxda[d][e] += E->X[level][e][E->IEN[level][el].node[i]] 
	    *  llNx[d][i];                      
    /* NOTE:  This is correct, Hughes book is in
       error and one needs to exchange (d<->e)
       wrt to his notation. */

    jacobian = determinant(dxda,dims);   /* may have this already ... check */  

    for(d=1;d<=dims;d++)
      for(e=1;e<=dims;e++)
	cof[d][e]=cofactor(dxda,d,e,dims); 
    
    for(j=1;j<=ends;j++)
      for(d=1;d<=dims;d++) {
	lNx[d][j] = 0.0;
	for(e=1;e<=dims;e++)
	  lNx[d][j] += 
	    llNx[e][j] * cof[e][d];

	lNx[d][j] /= jacobian;

	/* fprintf(stderr,"lNx[%d][%d] = %g v %g\n",d,j,lNx[d][j],llNx[d][j] * scale1); */
      }    
  return;
}
/*   ======================================================================
     Function to produce the appropriate one-dimensional global shape 
     function for a particular element edge. Referenced by the element/local
     node number.
     ======================================================================  */
	
void get_global_1d_shape_fn(
  struct All_variables *E,
  int el,
  struct Shape_function1 *GM,
  struct Shape_function1_dA *dGammax,
  int level
)
{ 
  int i,k,d,e;
  int dirn,locn,node[5];
  int collapsed_dirn[2];
  higher_precision scale[4];

  higher_precision jacobian;
  higher_precision determinant();
  higher_precision cofactor();
 
  void get_neighbour_nodes();

  static higher_precision dxda[4][4],cof[4][4];
  static int been_here = 0;
  
  if(E->control.AXI || E->control.ORTHO)  /* faster version */   {
    for(d=1;d<=E->mesh.nsd;d++)
	 scale[d] = E->ECO[level][el].recip_size[d];
      
	   for(i=1;i<=onedvpoints[E->mesh.nsd];i++) {
	       for(d=1;d<=E->mesh.nsd;d++)
		   dGammax->vpt[GMVGAMMA(d-1,i)] = E->ECO[level][el].area * scale[d];
	       for(d=1;d<=E->mesh.nsd;d++)
		   dGammax->vpt[GMVGAMMA(d-1+E->mesh.nsd,i)] = dGammax->vpt[GMVGAMMA(d-1,i)];
	     
	       if ( E->control.AXI) { /* change up and down dGamma for axi case [1, (nsd+1)] */
		   dGammax->vpt[GMVGAMMA(1,i)] *= 
		     2*M_PI*(E->ECO[level][el].size[1]*g_1d[i].x[0]+E->ECO[level][el].centre[1]);
		   dGammax->vpt[GMVGAMMA(E->mesh.nsd+1,i)] *= 
		     2*M_PI*(E->ECO[level][el].size[1]*g_1d[i].x[0]+E->ECO[level][el].centre[1]); 
	       }
	   }
     }

  else  /* Non orthogonal meshes:
	   The trick is to make sure the axis ordering is correct. The pairs should be
	   Z-X, X-Y, Y-Z for Y,Z,X normals 
	 */    {
      for(locn=0;locn<=1;locn++) /* top/bottom, front/back, left/right */ 
	  for(dirn=1;dirn<=E->mesh.nsd;dirn++) {
	      get_neighbour_nodes(node,dirn,locn);

	      if(3==E->mesh.nsd)
		  switch (dirn) {
		  case 1:  /* Y-Z */
		      collapsed_dirn[0]=3;
		      collapsed_dirn[1]=2;
		      break;
		  case 2:  /* X-Y */
		      collapsed_dirn[0]=1;
		      collapsed_dirn[1]=3;
		      break;
		  case 3:  /* Z-X */
		      collapsed_dirn[0]=2;
		      collapsed_dirn[1]=1;
		      break;
		  }
	      else
		  switch (dirn) {
		  case 1:  /* Z integral */
		      collapsed_dirn[0]=2;
		      break;
		  case 2:  /* X integral */
		      collapsed_dirn[0]=1;
		      break;
		  } 
    
	      for(k=1;k<=onedvpoints[E->mesh.nsd];k++) /* all of the vpoints */ {
		  for(d=1;d<=E->mesh.nsd-1;d++)
		      for(e=1;e<=E->mesh.nsd-1;e++)
			  dxda[d][e]=0.0;
	  
		for(i=1;i<=onedvpoints[E->mesh.nsd];i++) /* nodes */  {
		    for(d=1;d<=E->mesh.nsd-1;d++)
			for(e=1;e<=E->mesh.nsd-1;e++)
			    dxda[d][e] += E->X[level][collapsed_dirn[e-1]][E->IEN[level][el].node[node[i]]]*
				E->Mx.vpt[GMVXINDEX(d-1,i,k)];
		 }

		jacobian = determinant(dxda,E->mesh.nsd-1); 
		dGammax->vpt[GMVGAMMA(dirn-1+E->mesh.nsd*locn,k)] = jacobian;
	      }
	  }
    }
  
  return;
}

void get_neighbour_nodes(
     int node[5],
     int dirn,
     int locn  /* dirn is normal to the surface, loc is front/back/top/bottom/left/right*/
)
{ int a;  /* reference node, then use same (proven) scheme as 1d integration */

  switch(dirn)
    { case 1:    /* x vector normal */
	a = loc[1].node_nebrs[0][locn]; 
	node[1] = loc[loc[a].node_nebrs[2][0]].node_nebrs[1][0];
	node[2] = loc[loc[a].node_nebrs[2][0]].node_nebrs[1][1];
	node[4] = loc[loc[a].node_nebrs[2][1]].node_nebrs[1][0];
	node[3] = loc[loc[a].node_nebrs[2][1]].node_nebrs[1][1];
	break;
  
      case 2:    /* z vector normal */
	a = loc[1].node_nebrs[1][locn]; 
	node[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
	node[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
	node[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
	node[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];
	break;
	
      case 3:    /* y vector normal */
	a = loc[1].node_nebrs[2][locn]; 
	node[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[1][0];
	node[2] = loc[loc[a].node_nebrs[0][0]].node_nebrs[1][1];
	node[4] = loc[loc[a].node_nebrs[0][1]].node_nebrs[1][0];
	node[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[1][1];
	break;

      }

return;
}

/* Find the local coordinates for element `el', and tracer `num' */


void get_element_coords(
			struct All_variables *E,
			int el,
			int num,
			standard_precision *x,
			standard_precision *z,
			standard_precision *y,
			standard_precision *eta1,
			standard_precision *eta2,
			standard_precision *eta3,
			int level
			)
{
  int k,kk;
  int node;
  int lnode[28]; /* what's the #defined variable for the max nodes/element ? */
 
  standard_precision xx1,xx2,xx3;
  standard_precision x1,x2,x3;
  standard_precision etadash1,etadash2,etadash3;
  standard_precision distance;
  standard_precision dirn[5][4],mag;
  standard_precision lN[ELNMAX+1];
  
  standard_precision area_1;

  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];


  /* initial guess */

  *eta1 = *eta2 = *eta3 = 0.0;
  kk = 0;

  /* If periodic, we want actual (not wrapped around) coordinates 
     so we build this first */


/*RAA: 3/9/01 fixes below for x3, and if(2==dims), (3==dims), etc. 
      N.B., for 3D ..area is 1/8th of the element volume. */
  
  for(k=1;k<=ends;k++) {
    node = E->IEN[level][el].node[k];
    if((E->NODE[level][node] & (PER_OFFSIDE/* | OFFSIDE*/))) { /*RAA: 1/11/01, added (..| OFFSIDE) for perx and y */
      /* This node has ambiguous coordinates ! */
      
      x1 = E->X[level][1][node] - E->ECO[level][el].centre[1];
      x2 = E->X[level][2][node] - E->ECO[level][el].centre[2];
      if(3==dims) 
	x3 = E->X[level][3][node] - E->ECO[level][el].centre[3];
      
      if(2==dims) { /*RAA: added this distinction*/
	if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area) 
	    node += E->mesh.NOZ[level] * (E->mesh.NOX[level]-1);
      }
/*RAA, 3D part below is not finished yet! current fix is for aspect ratios of < 3.0 in all directions!*/
/*  and works by checking 2 dimensional squared distance versus 2D area, not 3D volume*/
      else if(3==dims && E->mesh.periodic_x && !E->mesh.periodic_y) {
        if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) {
	  node += E->mesh.NOZ[level] * (E->mesh.NOX[level]-1);
  
          /*RAA: 5/6/01, check out the periodic node numbers*/
/*        if(E->control.verbose)
            fprintf(stderr,"'if' condition met!!: Element %d, area: %g  local node, node num %d %d , x1, x2, x3: %g %g %g , center1, center2, center3:  %g %g %g\n",el,E->ECO[level][el].area,k,node,x1,x2,x3,E->ECO[level][el].centre[1],E->ECO[level][el].centre[2],E->ECO[level][el].centre[3]); */     

        }
      }
  
      else if(3==dims && E->mesh.periodic_y && !E->mesh.periodic_x)  {
        if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) {
	  node += (E->mesh.NOY[level]-1)* E->mesh.NOX[level] * E->mesh.NOZ[level];
        }
      }
  
      /*RAA: 22/10/01, added this part for both per_x and _y. The point is to figure out
         which PER_OFFSIDE face you are dealing with, and then make the appropriate node 
         # adjustment. The local node #s are set for either the front face or right face,
         since these are the faces that have the prospect of satisfying the 'if' condition
         for 4*area, while the back & left shouldn't. This fix assumes that the right side
         PER_OFFSIDE has a x coord of 0.0, and the front side PER_OFFSIDE has a y coord of 
         0.0.  N.B., perx and _y needed OFFSIDE added to line above.*/
      else if(3==dims && E->mesh.periodic_x && E->mesh.periodic_y)  {
        if ((k==5 || k==6 || k==7 || k==8) && (E->X[level][3][node]==0.0 && E->X[level][1][node]!=0.0)) {
          if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) {
	     node += (E->mesh.NOY[level]-1)* E->mesh.NOX[level] * E->mesh.NOZ[level];
          }
        }
        if ((k==3 || k==4 || k==7 || k==8) && (E->X[level][1][node]==0.0 && E->X[level][3][node]!=0.0)) {
          if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) {
	     node += E->mesh.NOZ[level] * (E->mesh.NOX[level]-1);
          }
        }
        if (E->X[level][1][node]==0.0 && E->X[level][3][node]==0.0)  {  
           if (k==3 || k==4 || k==7 || k==8) {   
              if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) 
	         node += E->mesh.NOZ[level] * (E->mesh.NOX[level]-1);
           }
           else if (k==5 || k==6)  {  /*this part gets the left front edge*/
              if( x1*x1 + x2*x2 > 4.0 * E->ECO[level][el].area/fabs(x3)) 
	         node += (E->mesh.NOY[level]-1)* E->mesh.NOX[level] * E->mesh.NOZ[level];
           }
        }  


      }  /*end of 'if' per_x and per_y */
    }  /*end of 'if' PER_OFFSIDE*/
    lnode[k] = node;
  }

  if(2==dims) {
    do {
      v_shape_fn(E,el,lN,*eta1,*eta2,*eta3,level);

      /* fprintf(stderr,"Shape function for tracer %d at %g, %g in element %d = %g,%g,%g,%g\n",
			num,*eta1,*eta2,el,lN[1],lN[2],lN[3],lN[4]); */

      xx1=xx2=0.0;
      for(k=1;k<=ends;k++) {
	node = lnode[k];  
	xx1 += E->X[level][1][node] * lN[k];
	xx2 += E->X[level][2][node] * lN[k];
      }

      x1 = x[num] - xx1;
      x2 = z[num] - xx2;

      distance = (x1*x1+x2*x2);
      
      /* fprintf(stderr,"tracer %d, element %d, mismatch %g,%g (%g-%g,%g-%g)\n",
      		num,el,x1,x2,x[num],xx1,z[num],xx2); */

      etadash1 = ( x1 * E->ECO[level][el].ntl_dirns[1][1] + 
		   x2 * E->ECO[level][el].ntl_dirns[1][2] ) * E->ECO[level][el].ntl_recip_size[1];

      etadash2 = ( x1 * E->ECO[level][el].ntl_dirns[2][1] + 
		   x2 * E->ECO[level][el].ntl_dirns[2][2] ) * E->ECO[level][el].ntl_recip_size[2];

      if(kk != 0) { /* Damping */
	*eta1 += 0.8 * etadash1;
	*eta2 += 0.8 * etadash2;
      }
      else {
	*eta1 += etadash1;
	*eta2 += etadash2;
      }
    
      if(/* (level==E->mesh.levmax && num==2064) || */++kk > 99)
	fprintf(stderr,"%d ... Tracer %d/%d in element %d ... eta (%g,%g v %g,%g) -> distance %g (%g,%g)\n",kk,
		num,level,el,*eta1,*eta2,x[num],z[num],distance,xx1,xx2); 

      /* Only need to iterate if this is marginal. If eta > distortion of 
	 an individual element then almost certainly the tracer is in a different element ...
	 or the mesh is terrible !  */

    } while((distance > E->ECO[level][el].area * E->control.accuracy * E->control.accuracy) && 
	    (fabs(*eta1) < 5.0) && (fabs(*eta2) < 5.0) &&
	    (kk < 100));
  }
  else /* 3==dims */ {

    do {     
      v_shape_fn(E,el,lN,*eta1,*eta2,*eta3,level);

/*------------------------------------------------------*/    
/* fprintf(E->fp1,"%d ...G'day!  Tracer %d/%d in element %d ... eta (%g,%g,%g v %g,%g,%g) -> distance %g (%g,%g,%g)\n",kk,num,level,el,*eta1,*eta2,*eta3,x[num],z[num],y[num],distance,xx1,xx2,xx3); 

    fprintf(E->fp1,"Tracer %d in el %d ... eta (%g,%g,%g)\n",num,el,*eta1,*eta2,*eta3); 
*/
/*------------------------------------------------------*/    

      /* If periodic, we want actual (not wrapped around) coordinates
	 NB - this currently assumes only periodic in x direction - 
	 and will need to be extended to y direction
      */

      xx1=xx2=xx3=0.0;
      for(k=1;k<=ends;k++) {
	node = lnode[k];
	
	xx1 += E->X[level][1][node] * lN[k];
	xx2 += E->X[level][2][node] * lN[k];
	xx3 += E->X[level][3][node] * lN[k];
      }
      
      x1 = x[num] - xx1;
      x2 = z[num] - xx2;
      x3 = y[num] - xx3;

      distance = (x1*x1+x2*x2+x3*x3);

      etadash1 = ( x1 * E->ECO[level][el].ntl_dirns[1][1] + 
		   x2 * E->ECO[level][el].ntl_dirns[1][2] +
		   x3 * E->ECO[level][el].ntl_dirns[1][3] ) * E->ECO[level][el].ntl_recip_size[1];

      etadash2 = ( x1 * E->ECO[level][el].ntl_dirns[2][1] + 
		   x2 * E->ECO[level][el].ntl_dirns[2][2] +
		   x3 * E->ECO[level][el].ntl_dirns[2][3] ) * E->ECO[level][el].ntl_recip_size[2];
   
      etadash3 = ( x1 * E->ECO[level][el].ntl_dirns[3][1] + 
		   x2 * E->ECO[level][el].ntl_dirns[3][2] +
		   x3 * E->ECO[level][el].ntl_dirns[3][3] ) * E->ECO[level][el].ntl_recip_size[3];

      if(kk == 0) {
	*eta1 += etadash1;
	*eta2 += etadash2;
	*eta3 += etadash3;
      }
      else /* Damping */{
	*eta1 += 0.8 * etadash1;
	*eta2 += 0.8 * etadash2;
	*eta3 += 0.8 * etadash3;
      }
    
      if(++kk > 10)
	fprintf(stderr,"%d ... Tracer %d/%d in element %d ... eta (%g,%g,%g v %g,%g,%g) -> distance %g (%g,%g,%g)\n",kk,
		num,level,el,*eta1,*eta2,*eta3,x[num],z[num],y[num],distance,xx1,xx2,xx3); 

      /* Only need to iterate if this is marginal. If eta > distortion of 
	 an individual element then almost certainly the tracer is in a different element ...
	 or the mesh is terrible !  */

    } while((distance > E->ECO[level][el].area * E->control.accuracy * E->control.accuracy) &&	
	    (fabs(*eta1) < 5.0) && (fabs(*eta2) < 5.0) && (fabs(*eta3) < 5.0) &&
	    (kk < 100)); /*RAA: changed 1.5 to 5.0 in this line, to correspond with 2D case */

  }

  return;
}

/* Find the local coordinates for all tracers for element list specified */

#if 1
void tr_local_coords(
   struct All_variables *E,
   int *el_list,
   struct TRACER_ELT_WEIGHT *lN,
   standard_precision *x,
   standard_precision *z,
   standard_precision *y,
   standard_precision *eta1,
   standard_precision *eta2,
   standard_precision *eta3,
   int N1,  /* expect 0 or 1 ... rather than some huge number */
   int N2,
   int level
)
{
  int k,kk,i,j,el,m;
  int node;
 
  void get_element_coords();
  standard_precision xx1,xx2,xx3;
  standard_precision x1,x2,x3;
  standard_precision distance;
  standard_precision dirn[5][4],mag;
  standard_precision damping;

  standard_precision maxD;

  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int elts = E->mesh.NEL[level];
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];


  standard_precision *ENX1,*ENX2,*ENX3;
  standard_precision *XX1,*XX2,*XX3;
  standard_precision *D;
  standard_precision *Etadash1,*Etadash2,*Etadash3;
  
  int *loop_again,loop_check;
  int *per_offset;


  struct COORD *ECO1;
  struct IEN *IEN1;
  unsigned int *NODE1;
  standard_precision *EX1,*EX2,*EX3;


  /* Vector arrays */
  

  ENX1 = (standard_precision *) Malloc0((elts*ends+1) * sizeof(standard_precision));
  ENX2 = (standard_precision *) Malloc0((elts*ends+1) * sizeof(standard_precision));
  ENX3 = (standard_precision *) Malloc0((elts*ends+1) * sizeof(standard_precision));

  XX1 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  XX2 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));
  XX3 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));

  D = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision));

  Etadash1 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision)); 
  Etadash2 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision)); 
  Etadash3 = (standard_precision *) Malloc0((N2+1) * sizeof(standard_precision)); 

  loop_again = (int *) Malloc0((N2+1) * sizeof(int)); 
  per_offset = (int *) Malloc0(((elts*ends+1)) * sizeof(int)); 


  /* untangle pointers a little */

  ECO1 = E->ECO[level];
  IEN1 = E->IENP[level];
  EX1 = E->X[level][1];
  EX2 = E->X[level][2];
  EX3 = E->X[level][3];
  NODE1 = E->NODE[level];

 
  
  if(2==dims)
/* #pragma loop novrec ENX1,ENX2,EX1,EX2,IEN1 */
    for(i=0;i<ends*elts;i++) {
      el   = i/ends+1;
      node = i%ends+1;  /* index as (el-1)*ends+node-1] */
      
      ENX1[i] = EX1[IEN1[el].node[node]];
      ENX2[i] = EX2[IEN1[el].node[node]]; 
    }
  else
/* #pragma loop novrec ENX1,ENX2,ENX3,EX1,EX2,EX3,IEN1 */
    for(i=0;i<ends*elts;i++) {
      el   = i/ends+1;
      node = i%ends+1;  /* index as (el-1)*ends+node-1] */
      ENX1[i] = EX1[IEN1[el].node[node]];
      ENX2[i] = EX2[IEN1[el].node[node]]; 
      ENX3[i] = EX3[IEN1[el].node[node]]; 
    }

 


  /* initial guess is that we know nothing about eta1,2,3 */

  if(2==dims)
/* #pragma loop novrec eta1,eta2 */
    for(i=N1;i<=N2;i++) {
      eta1[i] = eta2[i] = 0.0;
    }
  else
/* #pragma loop novrec eta1,eta2,eta3 */
    for(i=N1;i<=N2;i++) {
      eta1[i] = eta2[i] = eta3[i] = 0.0;
    }


  kk = 0;
  damping = 1.0;
  do {		

    all_v_shape_fn(E,lN,eta1,eta2,eta3,N1,N2,level);
   
      if(2==dims) {
/* #pragma loop novrec XX1,XX2 */
	for(m=N1;m<=N2;m++) {
	  XX1[m] = 0.0;
	  XX2[m] = 0.0;
	}
	for(k=1;k<=ends;k++) {
/* #pragma loop novrec XX1,XX2,ENX1,ENX2,lN */
	  for(m=N1;m<=N2;m++) {
	    el = el_list[m];
	    XX1[m] += ENX1[(el-1)*ends+k-1] * lN[m].node[k];
	    XX2[m] += ENX2[(el-1)*ends+k-1] * lN[m].node[k];
	  }
	}

/* #pragma loop novrec XX1,XX2,x,z */
	for(m=N1;m<=N2;m++) {
	  XX1[m] = x[m] - XX1[m];
	  XX2[m] = z[m] - XX2[m];
	}

/* #pragma loop novrec D,XX1,XX2,ECO1 */
	for(m=N1;m<=N2;m++) {
	  el = el_list[m];
	  D[m] = /* Distance from actual position */
	    (XX1[m]*XX1[m] + XX2[m]*XX2[m]);
	}

/* #pragma loop novrec Etadash1,Etadash2,ECO1,XX1,XX2 */
	for(m=N1;m<=N2;m++) {
	  el = el_list[m];
	  Etadash1[m] = ( XX1[m] * ECO1[el].ntl_dirns[1][1] + 
			  XX2[m] * ECO1[el].ntl_dirns[1][2] ) * ECO1[el].ntl_recip_size[1];
	  Etadash2[m] = ( XX1[m] * ECO1[el].ntl_dirns[2][1] + 
			  XX2[m] * ECO1[el].ntl_dirns[2][2] ) * ECO1[el].ntl_recip_size[2];
	}
/* #pragma loop novrec eta1,eta2,Etadash1,Etadash2 */
	for(m=N1;m<=N2;m++) {
	  eta1[m] += damping * Etadash1[m];
	  eta2[m] += damping * Etadash2[m];
	}

      }
      else {
/*#pragma loop novrec XX1,XX2,XX3 */
	for(m=N1;m<=N2;m++) {
	  XX1[m] = 0.0;
	  XX2[m] = 0.0;
	  XX3[m] = 0.0;
	}
	for(k=1;k<=ends;k++) {
/* #pragma loop novrec XX1,XX2,XX3,ENX1,ENX2,ENX3,lN */
	  for(m=N1;m<=N2;m++) {
	    el = el_list[m];
	    XX1[m] += ENX1[(el-1)*ends+k-1] * lN[m].node[k];
	    XX2[m] += ENX2[(el-1)*ends+k-1] * lN[m].node[k];
	    XX3[m] += ENX3[(el-1)*ends+k-1] * lN[m].node[k];
	  }
	}
/* #pragma loop novrec XX1,XX2,XX3,x,z,y */
	for(m=N1;m<=N2;m++) {
	  XX1[m] = x[m] - XX1[m];
	  XX2[m] = z[m] - XX2[m];
	  XX3[m] = y[m] - XX3[m];
	}

/* #pragma loop novrec D,XX1,XX2,XX3,ECO1 */
	for(m=N1;m<=N2;m++) {
	  el = el_list[m];
	  D[m] = /* Distance from actual position */
	    (XX1[m]*XX1[m] + XX2[m]*XX2[m] + XX3[m]*XX3[m]);
	}

/* #pragma loop novrec Etadash1,Etadash2,Etadash3,XX1,XX2,XX3,ECO1 */
	for(m=N1;m<=N2;m++) {
	  el = el_list[m];
	  Etadash1[m] = ( XX1[m] * ECO1[el].ntl_dirns[1][1] + 
			  XX2[m] * ECO1[el].ntl_dirns[1][2] +
			  XX3[m] * ECO1[el].ntl_dirns[1][3] ) * ECO1[el].ntl_recip_size[1];
	  Etadash2[m] = ( XX1[m] * ECO1[el].ntl_dirns[2][1] + 
			  XX2[m] * ECO1[el].ntl_dirns[2][2] +
			  XX3[m] * ECO1[el].ntl_dirns[2][3] ) * ECO1[el].ntl_recip_size[2];
	  Etadash3[m] = ( XX1[m] * ECO1[el].ntl_dirns[3][1] + 
			  XX2[m] * ECO1[el].ntl_dirns[3][2] +
			  XX3[m] * ECO1[el].ntl_dirns[3][3] ) * ECO1[el].ntl_recip_size[3];
	}

/*#pragma loop novrec eta1,eta2,eta3,Etadash1,Etadash2,Etadash3 */	
	for(m=N1;m<=N2;m++) {
	  eta1[m] += damping * Etadash1[m];
	  eta2[m] += damping * Etadash2[m];
	  eta3[m] += damping * Etadash3[m];
	}
      }

      /* Now test to see if D[m]'s are OK */


      mag = E->control.accuracy * E->control.accuracy;

/* #pragma loop novrec loop_again,D,ECO1,eta1,eta2,el_list     */
      for(m=N1;m<=N2;m++) {
	 el = el_list[m];
	 loop_again[m] = 
	   (D[m] > mag * ECO1[el].area) *
	   (fabs(eta1[m]) < 1.5) * 
	   (fabs(eta2[m]) < 1.5);
      }

      loop_check = 0;

      for(m=N1;m<=N2;m++) {
	loop_check += loop_again[m];
      }


      damping = 0.8;

      /* Only need to iterate if this is marginal. If eta > distortion of 
	 an individual element then almost certainly the tracer is in a different element ...
	 or the mesh is terrible !  */

  } while(loop_check &&
	  (kk < 100));


  free((void *) ENX1); 
  free((void *) ENX2); 
  free((void *) ENX3); 
 
  free((void *) XX1); 
  free((void *) XX2); 
  free((void *) XX3); 
 
  free((void *) D); 
 
  free((void *) Etadash1); 
  free((void *) Etadash2); 
  free((void *) Etadash3); 

  free((void *) loop_again);
  free((void *) per_offset);

  return;
}

#endif




/* The input is the natural particle coordinates(x,y,z,), the output is the local one (eta).
   This version is only approximate if the element is non-orthogonal. */ 


void general_element_coords(
   struct All_variables *E,
   int el,
   standard_precision x,
   standard_precision z,
   standard_precision y,
   standard_precision *eta1,
   standard_precision *eta2,
   standard_precision *eta3,
   int level
)
{
  const int dims = E->mesh.nsd;
  standard_precision xx1,xx2,xx3;

  if(3==dims) {
    xx1 = x - E->ECO[level][el].centre[1];
    xx2 = z - E->ECO[level][el].centre[2];
    xx3 = y - E->ECO[level][el].centre[3];
  
    *eta1 = ( xx1 * E->ECO[level][el].ntl_dirns[1][1] + 
	      xx2 * E->ECO[level][el].ntl_dirns[1][2] + 
	      xx3 * E->ECO[level][el].ntl_dirns[1][3] ) * E->ECO[level][el].ntl_recip_size[1];

    *eta2 = ( xx1 * E->ECO[level][el].ntl_dirns[2][1] + 
	      xx2 * E->ECO[level][el].ntl_dirns[2][2] + 
	      xx3 * E->ECO[level][el].ntl_dirns[2][3] ) * E->ECO[level][el].ntl_recip_size[2];

    *eta3 = ( xx1 * E->ECO[level][el].ntl_dirns[3][1] + 
	      xx2 * E->ECO[level][el].ntl_dirns[3][2] + 
	      xx3 * E->ECO[level][el].ntl_dirns[3][3] ) * E->ECO[level][el].ntl_recip_size[3];
  }
  else {
    xx1 = x - E->ECO[level][el].centre[1];
    xx2 = z - E->ECO[level][el].centre[2];

    *eta1 = ( xx1 * E->ECO[level][el].ntl_dirns[1][1] + 
	      xx2 * E->ECO[level][el].ntl_dirns[1][2] ) * E->ECO[level][el].ntl_recip_size[1];

    *eta2 = ( xx1 * E->ECO[level][el].ntl_dirns[2][1] + 
	      xx2 * E->ECO[level][el].ntl_dirns[2][2] ) * E->ECO[level][el].ntl_recip_size[2];

  }
  return;
}
