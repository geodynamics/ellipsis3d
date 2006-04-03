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


/*  */ 
/*  Here are the routines which process the results of each velocity solution, and call
    the relevant output routines. At this point, the velocity and pressure fields have
    been calculated and stored at the nodes. The only properties of the velocity field
    which are already known are those required to check convergence of the iterative
    scheme and so on. */

#include <math.h>
/* #include <stdlib.h> */ /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_new_velocity(
    struct All_variables *E,
    int ii
)
{ 
    int i,j,k,el,node;

    void surface_observables();
    void calculate_stream_function();
    void velocity_averages();

    higher_precision cos_theta,sin_theta,cos_phi,sin_phi;

    standard_precision CPU_time();
    
/*    standard_precision *P1,*P2,*P3;*/

    static int been_here=0;
    const int vpts = vpoints[E->mesh.nsd];
    const int dims = E->mesh.nsd ;

    /* if spherical, convert VX to spherical coordinates */

#if 1
    if(E->control.SPHERE ) {
      for(i=1;i<=E->mesh.nox;i++) {
	
	cos_phi=cos(E->sx[1][1+(i-1)*E->mesh.noz]);
	sin_phi=sin(E->sx[1][1+(i-1)*E->mesh.noz]);
	
	for(k=1;k<=E->mesh.noy;k++) {
	
	  cos_theta=cos(E->sx[3][1+(k-1)*E->mesh.noz*E->mesh.nox]);
	  sin_theta=sin(E->sx[3][1+(k-1)*E->mesh.noz*E->mesh.nox]);

	  for(j=1;j<=E->mesh.noz;j++) {
	  
	    node = j + (i-1) * E->mesh.noz + (k-1) * E->mesh.noz * E->mesh.nox;
	    if(0 && E->node[node] & (BC1 | BC3 | BC2)) /* rtf already (temporarily) */ {
	      E->sv[1][node] = E->V[1][node];
	      E->sv[2][node] = E->V[2][node];
	      E->sv[3][node] = E->V[3][node];
	    }
	    else {
	      E->sv[1][node] = 
		 cos_phi * E->V[1][node] +
		-sin_phi * E->V[2][node];
	      E->sv[2][node] =
		 cos_theta*sin_phi * E->V[1][node] +
	         cos_theta*cos_phi * E->V[2][node] +
		       sin_theta * E->V[3][node];
	      E->sv[3][node] =
		-sin_theta*sin_phi * E->V[1][node] +
	        -sin_theta*cos_phi * E->V[2][node] +
		       cos_theta * E->V[3][node];
	    }
	  }
	}
      }
      velocities_conform_bcs_6(E,E->sv[1],E->sv[2],E->sv[3],0,0,0,E->mesh.levmax);   
    }

    if(E->control.CYLINDER ) {
      for(i=1;i<=E->mesh.nox;i++) {
	
	cos_phi=cos(E->sx[1][1+(i-1)*E->mesh.noz]);
	sin_phi=sin(E->sx[1][1+(i-1)*E->mesh.noz]);
		
	for(j=1;j<=E->mesh.noz;j++) {
	  
	  node = j + (i-1) * E->mesh.noz;
	  if(0 && E->node[node] & (BC1 | BC2)) /* rtf already (temporarily) */ {
	    E->sv[1][node] = E->V[1][node];
	    E->sv[2][node] = E->V[2][node];
	  }
	  else {
	    E->sv[1][node] = 
	      cos_phi * E->V[1][node] +
	      -sin_phi * E->V[2][node];
	    E->sv[2][node] =
	      sin_phi * E->V[1][node] +
	      cos_phi * E->V[2][node];
	  }
	}
      }
      velocities_conform_bcs_6(E,E->sv[1],E->sv[2],E->sv[3],0,0,0,E->mesh.levmax);
    }

 
/*    P1=(standard_precision *)Malloc0((1+E->mesh.nel)*vpoints[E->mesh.nsd]*sizeof(standard_precision));
    P2=(standard_precision *)Malloc0((1+E->mesh.nel)*vpoints[E->mesh.nsd]*sizeof(standard_precision));
    P3=(standard_precision *)Malloc0((1+E->mesh.nel)*vpoints[E->mesh.nsd]*sizeof(standard_precision));
 */
    for(k=1;k<=E->mesh.noy;k++)
      for(i=1;i<=E->mesh.nox;i++) {
	E->slice.vxsurf[1][i+(k-1)*E->mesh.noy] = E->sv[1][1+(i-1)*E->mesh.noz+(k-1)*E->mesh.nox*E->mesh.noz];
	E->slice.vxsurf[2][i+(k-1)*E->mesh.noy] = E->sv[2][1+(i-1)*E->mesh.noz+(k-1)*E->mesh.nox*E->mesh.noz];
	if(3==dims)
		E->slice.vxsurf[3][i+(k-1)*E->mesh.noy] = E->sv[3][1+(i-1)*E->mesh.noz+(k-1)*E->mesh.nox*E->mesh.noz];
      }

/*    free((void *)P1);
    free((void *)P2);
    free((void *)P3);*/

#endif

  /*RAA: 29/3/01, surface observables needs to be corrected for 3D - do this later  (some things done on 6/11/02)*/
   if(2==dims)  {
      surface_observables(E,ii); 
   }
   else {
      report(E,"WARNING: surface observables for 3D are likely not completely updated/fixed");
      surface_observables(E,ii); 
      report(E,"Done surface observables\n");
   }

    /* calculate new stream function if 2d */
 
    if(!E->control.CYLINDER && E->mesh.nsd==2) 
      calculate_stream_function(E,E->Psi);

    report(E,"Into velocity averages\n");
    velocity_averages(E);
    /*fprintf(stderr,"HERE IS TRACER.PT (3)   %g \n",E->tracer.Pt[436]) ;	*/
    report(E,"Calculated velocity averages\n");

    return;
}

/* Obtain the various average properties of the velocity field (and
   other fields based on this information */


void velocity_averages(
		       struct All_variables *E
		       )
{ 
  int i,m,el,j,k,tr,sample,element;
    
  standard_precision return_bulk_value();
  standard_precision *v,vmax,vx,vz,vy;
  standard_precision bulk_val;
  standard_precision average;
  standard_precision eta1,eta2,eta3;
  standard_precision dudx[4][4];
  standard_precision lN[ELNMAX+1];

  standard_precision xx,zz/*,z1,z2,x1,x2*/;  /*RAA: 2/1/02, N.B., x1,x2,z1,z2 not used*/
  standard_precision yy; /*RAA: 2/1/02, added this for sampling*/
  standard_precision phi;
  standard_precision *sample_variable;
  standard_precision max,min;
  standard_precision maxsto,minsto;
 

  standard_precision dXdash11,dXdash12,dXdash21,dXdash22;
  standard_precision dXdash13,dXdash31,dXdash23,dXdash32,dXdash33; /*RAA: 7/12/01*/

  standard_precision *Node,*Elt,*Elt1;

  static int been_here = 0;
  static int sample_element[MAX_SAMPLE_PTS];

  const int dims = E->mesh.nsd;
 /* fprintf(stderr,"HERE IS TRACER.PT (4)   %g \n",E->tracer.Pt[436]) ;	*/
  if(E->control.verbose)
      fprintf(stderr,"Velocity processing\n");

  if(been_here++ == 0) {
    for(sample=0;sample<E->tracer.SAMPLE_PTS;sample++) {
      sample_element[sample] = 1;
    }
  }
 
  v = (standard_precision *)Malloc0((E->mesh.nno+2)*sizeof(standard_precision));
  Node = (standard_precision *)Malloc0((E->mesh.nno+2)*sizeof(standard_precision));

  vmax=0.0;

  for(i=1;i<=E->mesh.nno;i++) {
    v[i] = (E->V[1][i]*E->V[1][i] + E->V[2][i]*E->V[2][i] + ((3==dims) ?  E->V[3][i]*E->V[3][i] : 0.0));
    if(vmax < v[i]) 
      vmax=v[i];
  }
  return_horiz_ave(E,v,E->Have.vrms);
  for(i=1;i<=E->mesh.noz;i++)
    E->Have.vrms[i] = sqrt(E->Have.vrms[i]);
  
  E->monitor.Vmax=sqrt(vmax);
  E->monitor.Vsrms=E->Have.vrms[1];
  E->monitor.Vrms=sqrt(return_bulk_value(E,v,1));
  E->monitor.Vrms_surface=E->Have.vrms[1];
  E->monitor.Vrms_base=E->Have.vrms[E->mesh.noz];
   
  if(3==dims) {
    for(i=1;i<=E->mesh.nno;i++)
      v[i] = E->V[3][i]*E->V[3][i];
    E->monitor.Vyrms=sqrt(return_bulk_value(E,v,1));
  }
      
  return_horiz_rms(E,E->V[1],E->Have.V[1]);
  return_horiz_rms(E,E->V[2],E->Have.V[2]);
  if(3==dims)
    return_horiz_rms(E,E->V[3],E->Have.V[3]);


  /* Tracer strain measures */

  /* I. Scalar, damage-like quantity */

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	
    E->tracer.edot_integrated[m] += E->tracer.edot[m] * E->advection.timestep;
    
    if(E->viscosity.YIELD) {
      /* Increasing plastic strain measure */
      if (E->tracer.yielded[m] != 0.0) 
	E->tracer.edotp_integrated[m] += E->tracer.edot[m] * E->advection.timestep;
      
      /* Decay of strain by healing over time */
      if(E->tracer.visc[E->tracer.property_group[m]].yield_stress_E0dt != 0.0)
	E->tracer.edotp_integrated[m] *= 
	  max(0.0,1.0 -  E->advection.timestep / E->tracer.visc[E->tracer.property_group[m]].yield_stress_E0dt);
      
      /* Decay of strain by healing due to heating */
      if(E->tracer.T[m] > E->tracer.visc[E->tracer.property_group[m]].yield_stress_ET)
	E->tracer.edotp_integrated[m] = 0.0;  
    }
  }
  
 

 
  /* Obtain a continuous and smooth representation
     of edot for plotting purposes */
  

 

#if 1

  /* Nodal representation of strain-rate and deviatoric stress invariants
      (RAA...and depletion)*/
gs_tracers_to_nodes(E,E->edot,NULL,NULL,NULL,E->tracer.edot,E->mesh.levmax,0); 
gs_tracers_to_nodes(E,E->strs,NULL,NULL,NULL,E->tracer.strs,E->mesh.levmax,0); 
gs_tracers_to_nodes(E,E->strd,NULL,NULL,NULL,E->tracer.strd,E->mesh.levmax,0); 
gs_tracers_to_nodes(E,E->strd1,NULL,NULL,NULL,E->tracer.strd1,E->mesh.levmax,0); 
gs_tracers_to_nodes(E,E->depl,NULL,NULL,NULL,E->tracer.depl,E->mesh.levmax,0); 
/*RAA: 24/09/02, C. O'Neill - melting stuff */



#endif  

/* Store data for profiles ... */
  
for(sample=0;sample<E->tracer.SAMPLE_PTS;sample++) {
  /* Which  variable are we going to plot ? */
  
  switch(E->tracer.sample_type[sample]) {
    case 1:
      sample_variable=E->T;
      break;
    case 2:
      sample_variable=E->V[1];
      break;
    case 3:
      sample_variable=E->V[2];
      break;
    case 4:  
      sample_variable=E->nQ;
      break;
    case 5:  
      sample_variable=E->edot;
      break;
  case 6:  
    sample_variable=E->strd;
    break;
  case 7:  
    sample_variable=E->strs;
    break;
    
    
  case 8:
    sample_variable=E->V[3];
    break;
  case 9:
      sample_variable=E->V[4];
      break;
  case 10:
    sample_variable=E->V[5];
    break;
  case 11:
      sample_variable=E->V[6];
      break;
  case 12:  
    gs_tracers_to_nodes(E,Node,NULL,NULL,NULL,E->tracer.grain_size,E->mesh.levmax,0); 
    sample_variable=Node;
      break;
  case 13:
    gs_tracers_to_nodes(E,Node,NULL,NULL,NULL,E->tracer.sigma_n,E->mesh.levmax,0); 
    sample_variable=Node;
    break;
  case 14:  /*RAA: 24/09/02, C. O'Neill - melting stuff */
    sample_variable=E->depl;
    break; 
  case 15:
    sample_variable=E->strd1;
    break;
  }
    
    E->tracer.sampled_data[0][101] = E->tracer.sample_x[sample];
    E->tracer.sampled_data[1][101] = E->tracer.sample_z[sample];
    if(3==dims) /*RAA: 2/1/02, added this part for y */
      E->tracer.sampled_data[2][101] = E->tracer.sample_y[sample];
        
    if(E->tracer.sample_direction[sample] == 1) {
      zz = E->tracer.sample_z[sample]; 
      if(3==dims) /*RAA: 2/1/02, added this part for y, should be ok*/
	yy = E->tracer.sample_y[sample];
      max = -1.0e32;
      min =  1.0e32;
      
      for(i=0;i<=101;i++) {
	/*  fprintf(stderr,"Sample locn %d\n",i); */
	
	xx = E->tracer.sampled_data[0][i];
	/*printf("Processing velocity averages, calling general_tracers_element \n");*/

	if(2==dims) /*RAA: 2/1/02, added this distinction.*/ 
	   element = sample_element[sample] =
	       general_tracers_element(E,sample_element[sample],xx,zz,0.0,&eta1,&eta2,&eta3,E->mesh.levmax);
	else if(3==dims) /*RAA: 2/1/02, added this distinction.*/ 
	   element = sample_element[sample] =
	       general_tracers_element(E,sample_element[sample],xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);
	
	if(element == -1) 
	  continue;
	
	v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);
	
	phi=0.0;
	for(j=1;j<=enodes[E->mesh.nsd];j++) {
	  phi += sample_variable[E->ien[element].node[j]] * lN[j];
	}
	
	if(max < phi) 
	  max = phi;
	if(min > phi)
	  min = phi;
	
	E->tracer.sampled_data[3+sample][i] = phi;

      }
    }
    else if(E->tracer.sample_direction[sample] == 2) { /*RAA: 2/1/02, added this distinction for z-direction*/
      xx = E->tracer.sample_x[sample]; 
      if(3==dims) /*RAA: 2/1/02, added this part for y.*/
	 yy = E->tracer.sample_y[sample];
      max = -1.0e32;
      min =  1.0e32;
      printf("2 processing velcoties: calling general_tracers_element \n");
     /* fprintf(stderr,"HERE IS TRACER.PT (5)   %g \n",E->tracer.Pt[436]) ;	*/
      for(i=0;i<=101;i++) {
	zz = E->tracer.sampled_data[1][i];
	if(2==dims) /*RAA: 2/1/02, added this distinction.*/ 
	   element = sample_element[sample] =
	       general_tracers_element(E,sample_element[sample],xx,zz,0.0,&eta1,&eta2,&eta3,E->mesh.levmax);
	else if(3==dims) /*RAA: 2/1/02, added this distinction. */
	   element = sample_element[sample] =
	       general_tracers_element(E,sample_element[sample],xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);

	if(element == -1) 
	   continue;
	
	v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);
	
	phi=0.0;
	for(j=1;j<=enodes[E->mesh.nsd];j++)
	  phi += sample_variable[E->ien[element].node[j]] * lN[j];

	if(max < phi) 
	  max = phi;
	if(min > phi)
	  min = phi;
	
	E->tracer.sampled_data[3+sample][i] = phi;
      }
    }
    else if(3==dims && E->tracer.sample_direction[sample] == 3) { /*RAA: 2/1/02, added this section for y-direction*/
      xx = E->tracer.sample_x[sample]; 
      zz = E->tracer.sample_z[sample];
      max = -1.0e32;
      min =  1.0e32;

      printf("2 processing velcoties: calling general_tracers_element \n");
   /*   fprintf(stderr,"HERE IS TRACER.PT (6)   %g \n",E->tracer.Pt[436]) ;	*/
      for(i=0;i<=101;i++) {
        yy = E->tracer.sampled_data[2][i];
        element = sample_element[sample] =
        general_tracers_element(E,sample_element[sample],xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);
        if(element == -1) 
          continue;

        v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);

        phi=0.0;
        for(j=1;j<=enodes[E->mesh.nsd];j++)
          phi += sample_variable[E->ien[element].node[j]] * lN[j];

        if(max < phi) 
          max = phi;
        if(min > phi)
          min = phi;

        E->tracer.sampled_data[3+sample][i] = phi;
      }
    } /*RAA: end of else if for y-direction sampling*/

    if(E->tracer.sample_normalize[sample]) {
      E->tracer.sample_plotmax[sample] = max;
      E->tracer.sample_plotmin[sample] = min;
    }

    E->tracer.sample_value[sample] = E->tracer.sampled_data[3+sample][101];
  } /*RAA: end of 'for' w/ sample*/

free((void *) v);
free((void *) Node);


/*  fprintf(stderr,"HERE IS TRACER.PT (7) \n") ;	*/
  if(E->control.verbose)
      fprintf(stderr,"Velocity processing ... done\n");

  printf("Velocity processing ... done\n");

  
  return;
}

/*  PSI: only acceptable if 2d solutions */
void calculate_stream_function(
    struct All_variables *E,
    standard_precision *Psi
)

{
    int i,j,that_node,this_node;
    standard_precision x1,x2;
    if (E->mesh.nsd != 2) 
	return;

    Psi[E->mesh.noz] = 0.0 ;	  /*  Define it so */

    for(i=2;i<=E->mesh.nox;i++)   /* work along lowest row */ {
	that_node = (i-1) * E->mesh.noz;
	this_node = i * E->mesh.noz;
	x1 = x2 = 1.0;
	    if (E->control.AXI) {
		x1 = E->x[1][this_node];
		x2 = E->x[1][that_node]; 
	    }    
	Psi[this_node] = Psi[that_node]
	  + (E->x[1][this_node] - E->x[1][that_node] ) 
	    * (x1 * E->V[2][this_node] + x2 * E->V[2][that_node]) /2.0;
    }

    for(i=1;i<=E->mesh.nox;i++)
      for(j=E->mesh.noz;j>1;j--){
	  that_node = (i-1) * E->mesh.noz + j;
	  this_node = that_node - 1;
	  x1 = x2 = 1.0;
	  if (E->control.AXI)  {
	      x1 = E->x[1][this_node];
	      x2 = E->x[1][that_node];
	  }
	  Psi[this_node] = Psi[that_node]-		
	      (E->x[2][this_node] - E->x[2][that_node])*
	      (x1 * E->V[1][this_node]+ x2 * E->V[1][that_node])/2.0; 
      }
    return;
}
