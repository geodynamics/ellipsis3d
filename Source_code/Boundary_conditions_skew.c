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


#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>


/* Rotation matrix to rotate boundary node velocities
   to a more natural reference frame where degrees of
   freedom can be eliminated */

void compute_node_R(
  struct All_variables *E
)
{
  int i,j,k,ii,jj,kk;
  int lev;
  
  higher_precision *node_R;
  higher_precision cos_theta,sin_theta,cos_phi,sin_phi;

  static int been_here=0;

  const int dofs = E->mesh.dof;
  const int dims = E->mesh.nsd;

  /* At present, no skew bc's are used in
     purely cartesian problems */

  if(E->control.CART2D || E->control.CART3D)
    return;

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {

    for(kk=1;kk<=E->mesh.NOY[lev];kk++)
      for(ii=1;ii<=E->mesh.NOX[lev];ii++)
	for(jj=1;jj<=E->mesh.NOZ[lev];jj++)   {

	  i = jj + (ii-1) * E->mesh.NOZ[lev] + (kk-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];

	  if(E->NODE[lev][i] & SKEWBC) { /* in principle, rotation matrix required */

	      /* Allocate memory */
	    if(! been_here)
	      E->curvilinear.NODE_R[lev][i] = (higher_precision *) Malloc0(dofs*dofs*sizeof(higher_precision));
	    node_R = E->curvilinear.NODE_R[lev][i];

	    /* fill matrix */
	
	    if(E->control.SPHERE) {
	     
	      cos_phi=E->curvilinear.cosph[lev][ii];
	      sin_phi=E->curvilinear.sinph[lev][ii];
	      cos_theta=E->curvilinear.cost[lev][kk];
	      sin_theta=E->curvilinear.sint[lev][kk];
 
	      node_R[0*dims+0] =  cos_phi;
	      node_R[0*dims+1] =  cos_theta*sin_phi;
	      node_R[0*dims+2] = -sin_theta*sin_phi;
	      node_R[1*dims+0] = -sin_phi;
	      node_R[1*dims+1] =  cos_theta*cos_phi;
	      node_R[1*dims+2] = -sin_theta*cos_phi;
	      node_R[2*dims+0] =  0.0;
	      node_R[2*dims+1] =  sin_theta;
	      node_R[2*dims+2] =  cos_theta;
  	    }
      
	    if(E->control.CYLINDER) {
	      cos_phi=E->curvilinear.cosph[lev][ii];
	      sin_phi=E->curvilinear.sinph[lev][ii];
	      
	      node_R[0*dims+0] =  cos_phi;
	      node_R[0*dims+1] =  sin_phi;
	      node_R[1*dims+0] = -sin_phi;
	      node_R[1*dims+1] =  cos_phi;
	    }
	  }
	}
  }
  been_here=1;
  return;
}

/* Rotation matrix to rotate boundary node velocities
   to a more natural reference frame where degrees of
   freedom can be eliminated */

void store_node_R(
  struct All_variables *E,
  int node,
  standard_precision phi,
  standard_precision theta,
  int level
)
{
  int i,j,k,ii,jj,kk;
 
  higher_precision *node_R;
  higher_precision cos_theta,sin_theta,cos_phi,sin_phi;
  
  static int been_here=0;
  
  const int dims = E->mesh.nsd;

  if(!(E->NODE[level][node] & SKEWBC))  /* rotation matrix not required */
    return;
  
  node_R = E->curvilinear.NODE_R[level][node];

  cos_phi = cos(phi);
  sin_phi = sin(phi);

  /* fill matrix */
	
  if(3==E->mesh.nsd) {

    cos_theta = cos(theta);
    sin_theta = sin(theta);

    node_R[0*dims+0] =  cos_phi;
    node_R[0*dims+1] =  cos_theta*sin_phi;
    node_R[0*dims+2] = -sin_theta*sin_phi;
    node_R[1*dims+0] = -sin_phi;
    node_R[1*dims+1] =  cos_theta*cos_phi;
    node_R[1*dims+2] = -sin_theta*cos_phi;
    node_R[2*dims+0] =  0.0;
    node_R[2*dims+1] =  sin_theta;
    node_R[2*dims+2] =  cos_theta;
  }
  else {
    node_R[0*dims+0] =  cos_phi;
    node_R[0*dims+1] =  sin_phi;
    node_R[1*dims+0] = -sin_phi;
    node_R[1*dims+1] =  cos_phi;
  }
  return;
}
