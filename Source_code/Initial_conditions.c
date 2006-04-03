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


#include <signal.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

/*
Initialize values which are not problem dependent.
NOTE: viscosity may be a function of all previous input fields (temperature,
pressure, velocity, chemistry) and so is always to be done last.
*/
void common_initial_fields(
    struct All_variables *E
)
{
    void initial_qterm();
    void read_viscosity_option();
    void initial_velocity();

    report(E,"Initialize velocity");
    initial_velocity(E);
    report(E,"Initialize pressure field");
    fprintf(stderr,"Initial cons, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);
    initial_qterm(E);
    fprintf(stderr,"Initial cons, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);

    
    return;
}

void initial_qterm(
     struct All_variables *E
)
{
    int i;
    standard_precision *Q;
    FILE *fp;
    void p_to_centres();

    Q =  (standard_precision *)Malloc0((E->mesh.nno+2)*sizeof(standard_precision)); 
   
    if(read_previous_field(E,Q,"pressure","Pres"))
	sp_to_centres(E,Q,E->Q,E->mesh.levmax);
    
    else
	for(i=1;i<=E->mesh.npno;i++)
	    E->Q[i]=0.0;

    for(i=E->mesh.levmax;i>E->mesh.levmin;i--) {
      project_q(E,E->QQ[i],E->QQ[i-1],i);
    }
    free((void *)Q);
    return; 
}

void initial_velocity(
     struct All_variables *E
)
{
  int i,node,ii,level;
  FILE *fp;
  
  if(!read_previous_field(E,E->V[1],"velocity","Velx"))
    for(i=1;i<=E->mesh.nno;i++)
      E->V1[1][i]=E->V[1][i]=0.0;
  if(!read_previous_field(E,E->V[2],"velocity","Velz"))
    for(i=1;i<=E->mesh.nno;i++)
      E->V1[2][i]=E->V[2][i]=0.0;

  switch(E->control.model) {
  case 1:
    if(3==E->mesh.dof && !read_previous_field(E,E->V[3],"velocity","Vely"))
      for(i=1;i<=E->mesh.nno;i++)
	E->V1[3][i]=E->V[3][i]=0.0;  
    break;
  case 2:
    if(3==E->mesh.dof && 2==E->mesh.nsd && !read_previous_field(E,E->V[3],"rotation","Roty"))
      for(i=1;i<=E->mesh.nno;i++)
	E->V1[3][i]=E->V[3][i]=0.0;
    if(6==E->mesh.dof) {
      if(!read_previous_field(E,E->V[3],"velocity","Vely"))
	for(i=1;i<=E->mesh.nno;i++)
	  E->V1[3][i]=E->V[3][i]=0.0;  
      if(!read_previous_field(E,E->V[4],"rotation","Rotx"))
	for(i=1;i<=E->mesh.nno;i++)
	  E->V1[4][i]=E->V[4][i]=0.0;  
      if(!read_previous_field(E,E->V[5],"rotation","Rotz"))
	for(i=1;i<=E->mesh.nno;i++)
	  E->V1[5][i]=E->V[5][i]=0.0;
      if(!read_previous_field(E,E->V[6],"rotation","Roty"))
	for(i=1;i<=E->mesh.nno;i++)
	  E->V1[6][i]=E->V[6][i]=0.0;
    }
    break ;
  }

  /* Now inject high level velocity field
     to all levels in the mesh (V -> VV) */

  for(level=E->mesh.levmax;level>E->mesh.q_levmin;level--) {
    inject06(E,level,E->VV[level][1],E->VV[level][2],E->VV[level][3],
	     E->VV[level][4],E->VV[level][5],E->VV[level][6],
	     E->VV[level-1][1],E->VV[level-1][2],E->VV[level-1][3],
	     E->VV[level-1][4],E->VV[level-1][5],E->VV[level-1][6]
	     );
  }
  return; 
}
