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




/* Functions relating to the building and use of mesh locations ... */

#include "config.h"

#include <math.h>

#if HAVE_STRING_H
#include <string.h>
#endif

#if HAVE_STRINGS_H
#include <strings.h>
#endif

#if HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"

int M_node_cols;
int *M_node_col_num[MAX_LEVELS];
int **M_node_col_list[MAX_LEVELS];


extern int Emergency_stop;

/* =================================================
   Standard node positions including mesh refinement 
   =================================================  */

void node_locations(
		    struct All_variables *E,
		    int step
		    )
{ 
  int i,j,k,ijk[4],d,node,lv;
 
  standard_precision *XX[4];
  standard_precision smallest,this_space;
  void set_spacing_functions();
  void inject_node_values();

  const int dims = E->mesh.nsd ;
  const int dofs = E->mesh.dof ;

  set_spacing_functions(E);
 

  for(d=1;d<=dims;d++) {
    if(E->control.verbose)
      fprintf(stderr,"Computing nodes locations for dirn %d\n",d);
    XX[d] = (standard_precision *)Malloc0((2+E->mesh.nnx[d])*sizeof(standard_precision)); 
  }
  /* First build nsd arrays along each edge & then convolve them (orthogonal) 
     Initially the coordinates are determined on the unit interval and then scaled/offset */

   for(d=1;d<=dims;d++) {
       XX[d][1] = 0.0;
       for(i=2;i<=E->mesh.nnx[d];i++) {
	   XX[d][i] = XX[d][i-1] +
	       (E->node_space_function[d-1])(E,(i-1.5)/(E->mesh.nnx[d]-1),d,i); }
   }

  /* Correct for scaling errors. This is simpler than obtaining the integral of the
     scaling function. */
 
  for(d=1;d<=dims;d++)    {
      for(i=1;i<=E->mesh.nnx[d];i++) {
	  XX[d][i] = E->mesh.layer0[d] + XX[d][i] * (E->mesh.layer1[d]-E->mesh.layer0[d])/(XX[d][E->mesh.nnx[d]]); 
	}
  }

  /* find smallest spacing in each direction */ 

  E->mesh.equiv_even_nodes[3] = 1; /* useful to define this when doing 2d */
  for(d=1;d<=dims;d++) {
    smallest=E->mesh.layer1[d];
    for(i=1;i<E->mesh.nnx[d];i++) {
      this_space = XX[d][i+1]-XX[d][i];
      if(this_space < smallest)
	smallest = this_space;
    }
     E->mesh.equiv_even_nodes[d] = (int)((E->mesh.layer1[d]-E->mesh.layer0[d])/smallest) + 1;
     E->mesh.fine_spacing[d] = E->mesh.layer1[d] / (E->mesh.equiv_even_nodes[d]-1);
  }

  E->mesh.equiv_even_nsf = E->mesh.equiv_even_nodes[1] *  E->mesh.equiv_even_nodes[3];
  
  /* 1.  If the problem is Cartesian, then this set of coordinates is
     used as is.  Define Cartesian mesh from point distributions on 3 orthogonal axes */

  if (E->control.CART2D || E->control.CART2pt5D || E->control.CART3D || E->control.AXI) {
  
    for(ijk[1]=1;ijk[1]<=E->mesh.nox;ijk[1]++)
      for(ijk[2]=1;ijk[2]<=E->mesh.noz;ijk[2]++)
	for(ijk[3]=1;ijk[3]<=E->mesh.noy;ijk[3]++)	{ 
	  node = ijk[2]+(ijk[1]-1)*E->mesh.noz+(ijk[3]-1)*E->mesh.noz*E->mesh.nox;
	  for(d=1;d<=dims;d++) {
	    E->x[d][node] = XX[d][ijk[d]]; 
	  }
	}

    /* AD HOC adjustment for funnel problem */

#if 0
    E->control.ORTHO = 0;

    for(i=1;i<=E->mesh.nno;i++) {
      if(E->x[2][i] > 2.76 && E->x[2][i] <= 3.375 ) 
      E->x[1][i] += 1.0 * (E->x[2][i]-2.75) * (1.0 - E->x[1][i]);
      if(E->x[2][i] > 3.375)
	E->x[1][i] += 0.625 * (1.0 - E->x[1][i]);
    }
#endif    

    /* AD HOC adjustment for interface problem */

#if 0
    E->control.ORTHO = 0;
    for(i=1;i<=E->mesh.nno;i++) {
      if( fabs(E->x[2][i] - (0.5 + 0.01 * cos(2.*M_PI * E->x[1][i] / 1.284)) ) < 0.4999*fabs(E->x[2][1] - E->x[2][2])) {
        E->x[2][i] = 0.5 + 0.01 * cos(2.*M_PI * E->x[1][i] / 1.284) ;
        fprintf(stderr,"Node N0 %d has been moved to z = %.10e \n",i,E->x[2][i]) ;
      }
    }
#endif    

    
    /* AD HOC adjustment for Rayleigh-Taylor instability problem */

#if 0
    E->control.ORTHO = 0;
    for(i=1;i<=E->mesh.nno;i++) {
      if(E->x[2][i] == 0.5) {
	for(j=i-0.5*(sqrt(E->mesh.nno)-1);j<=i;j++)
	  E->x[2][j] += 0.01 * 0.5 * cos(2.*M_PI * E->x[1][i] / 1.284) * (E->x[2][j]/0.5) ;
	for(j=i+1;j<=i+0.5*(sqrt(E->mesh.nno)-1);j++)
	  E->x[2][j] += 0.01 * 0.5 * cos(2.*M_PI * E->x[1][i] / 1.284) * ((1-E->x[2][j])/0.5) ;
      	fprintf(stderr,"Node N0 %d has been moved to z = %.10e \n",i,E->x[2][i]) ;
      }
    }
#endif    

    /* AD HOC adjustment for Buckling / Folding instability problem */

#if 0
    /* Here you can put ad hoc modifications for specific problems !!! */
    input_std_precision_vector("test_geometry_parameters",5,geom);  /* What if fewer entries ?? */
    E->control.secret_case_fl[1] = geom[0] ;
    E->control.secret_case_fl[2] = geom[1] ;
    E->control.secret_case_fl[3] = geom[2] ;
    fprintf(stderr,"For buckling problem (nodes) w0 = %g k = %g and h = %g \n",geom[0],geom[1],geom[2]) ;
    E->control.ORTHO = 0;
#if 1
/* All nodes are moved with a linear interpolation between the edge and the centered line */
    for(i=1;i<=E->mesh.NOX[E->mesh.levmax];i++) {
      for(j=(i-1)*E->mesh.NOZ[E->mesh.levmax]+1;j<=i*E->mesh.NOZ[E->mesh.levmax];j++) {
	if(E->x[2][j]<=E->x[2][E->mesh.nno]*0.5) {
	  E->x[2][j] += (geom[0] * cos(geom[1]*M_PI/3.0*E->x[1][j])) * E->x[2][j] / (E->x[2][E->mesh.nno]*0.5) ;
	}
	else {
	  E->x[2][j] += (geom[0] * cos(geom[1]*M_PI/3.0*E->x[1][j])) * (E->x[2][E->mesh.nno] - E->x[2][j]) / (E->x[2][E->mesh.nno]*0.5) ;
	}
      }
    }
#else
/* Only the node on the middle line are slightly moved */
    for(i=1;i<=E->mesh.nno;i++) {
      if(E->x[2][i] <= (0.5*E->x[2][E->mesh.nno]-geom[2]/2.+0.001) && E->x[2][i] >= (0.5*E->x[2][E->mesh.nno]-geom[2]/2.-0.001)) {
	E->x[2][i] += geom[0] * cos(geom[1]*M_PI/3.0*E->x[1][i]) ;
      }
      else if(E->x[2][i] <= (0.5*E->x[2][E->mesh.nno]+geom[2]/2.+0.001) && E->x[2][i] >= (0.5*E->x[2][E->mesh.nno]+geom[2]/2.-0.001)) {
	E->x[2][i] += geom[0] * cos(geom[1]*M_PI/3.0*E->x[1][i]) ;
      }
    }
#endif
#endif  


    /* Now inject to all levels of the multigrid formulation */
    
    for(i=E->mesh.levmax;i>E->mesh.levmin;i--) {
      inject_node_values(E,i,E->X[i][1],E->X[i-1][1]); 
      inject_node_values(E,i,E->X[i][2],E->X[i-1][2]);
      if(dims==3)
	inject_node_values(E,i,E->X[i][3],E->X[i-1][3]);
    }

    /* and alias the curvilinear mesh coords to the cartesian ones (for bc's etc) */
    
    E->sx[1] =  E->x[1]; 
    E->sx[2] =  E->x[2]; 
    E->sx[3] =  E->x[3]; 

    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
      E->SX[i][1] =  E->X[i][1];
      E->SX[i][2] =  E->X[i][2];
      E->SX[i][3] =  E->X[i][3];
    }
  }

  /* 2. If mesh is spherical or cylindrical then there are two sets of coordinates: the
     curvilinear ones and a second representation in Cartesian geometry to use for
     element computations */
  
  else {   /* Logical mesh coordinates are in phi,r,theta space (E->sx) */
  
    for(ijk[1]=1;ijk[1]<=E->mesh.nox;ijk[1]++)
      for(ijk[2]=1;ijk[2]<=E->mesh.noz;ijk[2]++)
	for(ijk[3]=1;ijk[3]<=E->mesh.noy;ijk[3]++) { 
	  node = ijk[2]+(ijk[1]-1)*E->mesh.noz+(ijk[3]-1)*E->mesh.noz*E->mesh.nox;
	  for(d=1;d<dims;d++) {
	    E->sx[d][node] = XX[d][ijk[d]]; 
	  }
	}

    /* inject to all levels */

    for(i=E->mesh.levmax;i>E->mesh.levmin;i--) {
      inject_node_values(E,i,E->SX[i][1],E->SX[i-1][1]); 
      inject_node_values(E,i,E->SX[i][2],E->SX[i-1][2]);
      if(dims==3)
	inject_node_values(E,i,E->SX[i][3],E->SX[i-1][3]);
    }
    
    /* Obtain cartesian coordinates from here */
    /* SPHERICAL */
    if(E->control.SPHERE) {
      for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
	for(i=1;i<=E->mesh.NOX[lv];i++) {
	  node = 1 + (i-1) * E->mesh.NOZ[lv];
	  E->curvilinear.cosph[lv][i] = cos(E->SX[lv][1][node]);
	  E->curvilinear.sinph[lv][i] = sin(E->SX[lv][1][node]);
	}
	  
	for(k=1;k<=E->mesh.NOY[lv];k++) {
	  node = 1 + (k-1) * E->mesh.NOZ[lv]*E->mesh.NOX[lv];
	  E->curvilinear.cost[lv][k] = cos(E->SX[lv][3][node]);
	  E->curvilinear.sint[lv][k] = sin(E->SX[lv][3][node]);
	}
      
	      
      for(i=1;i<=E->mesh.NOX[lv];i++) 
	for(k=1;k<=E->mesh.NOY[lv];k++)
	  for(j=1;j<=E->mesh.NOZ[lv];j++)  {
	    node = j + (i-1) * E->mesh.NOZ[lv]  + (k-1) * E->mesh.NOZ[lv] * E->mesh.NOX[lv] ;
	    
	    E->X[lv][1][node] =  E->SX[lv][2][node] * E->curvilinear.cost[lv][k] * E->curvilinear.sinph[lv][i];  /* X */ 
	    E->X[lv][2][node] =  E->SX[lv][2][node] * E->curvilinear.cost[lv][k] * E->curvilinear.cosph[lv][i];   /* Z */
	    E->X[lv][3][node] =  E->SX[lv][2][node] * E->curvilinear.sint[lv][k];                         /* Y */
	  }
      }
    }

    if(E->control.CYLINDER) {

      for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
	for(i=1;i<=E->mesh.NOX[lv];i++) {
	  node = 1 + (i-1) * E->mesh.NOZ[lv];
	  E->curvilinear.cosph[lv][i] = cos(E->SX[lv][1][node]);
	  E->curvilinear.sinph[lv][i] = sin(E->SX[lv][1][node]);
	}
	
	for(i=1;i<=E->mesh.NOX[lv];i++) 
	  for(j=1;j<=E->mesh.NOZ[lv];j++) {
	    node = j + (i-1) * E->mesh.NOZ[lv];
	    
	    E->X[lv][1][node] =  E->SX[lv][2][node] * E->curvilinear.sinph[lv][i];   /* X */ 
	    E->X[lv][2][node] =  E->SX[lv][2][node] * E->curvilinear.cosph[lv][i];   /* Z */

	    /* fprintf(stderr,"%d, Coords: %d/%d (%g,%g) -> (%g,%g)\n",lv,node,
		    E->mesh.NNO[lv],E->SX[lv][2][node],E->SX[lv][2][node],
		    E->X[lv][2][node],E->X[lv][2][node]); */
	  }
      }
    }
  }

  if(E->control.verbose)
      fprintf(stderr,"Computing all node locations\n");

  for(d=1;d<=dims;d++)
    free((void *)XX[d]);
 
return;  }


void set_spacing_functions(
     struct All_variables *E
)
{ 
  standard_precision regular_spacing();
  standard_precision bound_lyr_spacing();
  standard_precision region_spacing();
  standard_precision file_spacing();

  int i;
  char parsing_string[100];
  char direction[3][2];

  const int dims = E->mesh.nsd ;

  strcpy(direction[0],"x");
  strcpy(direction[1],"z");
  strcpy(direction[2],"y");

  for(i=0;i<dims;i++){
    switch(E->control.GRID_TYPE[i]) {
    case 1:
      E->node_space_function[i]=region_spacing;
      sprintf(parsing_string,"regions_%s",direction[i]);
      input_int(parsing_string,&(E->mesh.spacing[i].numb),"0,0,39");

      fprintf(E->fp,"Regional (%d) spacing of nodes in %s direction\n",E->mesh.spacing[i].numb,direction[i]);
      fflush(E->fp);

      sprintf(parsing_string,"region_offset_%s",direction[i]);
      input_std_precision_vector(parsing_string,E->mesh.spacing[i].numb,E->mesh.spacing[i].offset);
      sprintf(parsing_string,"region_compression_%s",direction[i]);
      input_std_precision_vector(parsing_string,E->mesh.spacing[i].numb,E->mesh.spacing[i].compression);
      sprintf(parsing_string,"region_width_%s",direction[i]);
      input_std_precision_vector(parsing_string,E->mesh.spacing[i].numb,E->mesh.spacing[i].width);
      sprintf(parsing_string,"region_edge_width_%s",direction[i]);
      input_std_precision_vector(parsing_string,E->mesh.spacing[i].numb,E->mesh.spacing[i].edge_width);

      break;
    case 2:
      E->node_space_function[i]=file_spacing;
      fprintf(E->fp,"Spacing of nodes from files for %s direction\n",direction[i]);fflush(E->fp);
      sprintf(parsing_string,"node_locations_file_%s",direction[i]);
      input_string(parsing_string,E->mesh.gridfile[i],"");

      break;
    }
  }
  return;
}

/* location is a coordinate in direction d (1,2,3) at edge node i */


standard_precision region_spacing(
     struct All_variables *E,
     standard_precision locn,
     int d,
     int i
)
{
  standard_precision location,delta,scale,radius,spacing[40],output,recip_width,width_scale;
  int ii;
  
  /* The commented expression causes the code to misinterpret grid element compression instructions. */
  width_scale = E->mesh.layer1[d] ;   /* /E->mesh.layer1[2]; */
  recip_width = 1.0/width_scale;
  locn *= width_scale;
  
  scale = 1.0/(E->mesh.nnx[d]-1.0); 
  
  for(ii=0;ii<E->mesh.spacing[d-1].numb;ii++) {
    radius = (locn-E->mesh.spacing[d-1].offset[ii])*(locn-E->mesh.spacing[d-1].offset[ii]);
    if(radius < E->mesh.spacing[d-1].width[ii]*E->mesh.spacing[d-1].width[ii]) {
      spacing[ii]=scale*E->mesh.spacing[d-1].compression[ii];
    }
    else {
      radius -= E->mesh.spacing[d-1].width[ii]*E->mesh.spacing[d-1].width[ii];
      spacing[ii]= scale*(1.0 - (1.0-E->mesh.spacing[d-1].compression[ii]) *
			  exp(-radius/(E->mesh.spacing[d-1].edge_width[ii]*E->mesh.spacing[d-1].edge_width[ii])) );
    }
  }
  output=scale;
  for(ii=0;ii<E->mesh.spacing[d-1].numb;ii++) {
    output=min(output,spacing[ii]);
  }
  return(output);
}

standard_precision file_spacing(
     struct All_variables *E,
     standard_precision locn,
     int d,
     int i
)
{ /* note, the format here is a little bit odd because the function should return the current spacing not
     the absolute position, hmmm. Also, the layer dimensions ought to be reset to reflect the values in the
     files and not any previously set values */

  static int been_here[4];
  const int dims=E->mesh.nsd; 
  const int dofs=E->mesh.dof;
  standard_precision *location[4];
  float temporary;
  int dd,ii,err;
  FILE *fp;

  if(0==been_here[d]) {
	location[d]=(standard_precision *)Malloc0((E->mesh.nox+1)*sizeof(standard_precision));
		
	/* open files */
	err=0.0;
	if((fp=fopen(E->mesh.gridfile[d-1],"r")) == NULL)   { 
	    fprintf(E->fp,"Error reading the grid information file: %s\n",E->mesh.gridfile[d-1]); 
	    fflush(E->fp);
	    exit(1);
	}
	else {
	  for(ii=1;ii<=E->mesh.nnx[d];ii++) {
		fscanf(fp,"%f",&(temporary));
		location[d][ii] = (standard_precision) temporary;
	  }
	    fclose(fp);
	}
  
	E->mesh.layer1[d] = location[d][E->mesh.nnx[d]];
	
	been_here[d]++;
  }
  return(location[d][i]-location[d][i-1]);
}


void dlogical_mesh_to_real(
     struct All_variables *E,
     higher_precision *data,
     int level
)
{ int i,j,n1,n2;

  if(E->mesh.periodic_x)
    for(i=1;i<=E->mesh.NOZ[level];i++)
      for(j=1;j<=E->mesh.NOY[level];j++)
	{ n1 = i + (j-1)* E->mesh.NOX[level]*E->mesh.NOZ[level];
	  n2 = n1 + (E->mesh.NOX[level]-1)*E->mesh.NOZ[level];
	  data[n2] = data[n1];
	}

  if(E->mesh.periodic_y)
    for(i=1;i<=E->mesh.NOZ[level];i++)
      for(j=1;j<=E->mesh.NOX[level];j++)
	{ n1 = i + (j-1)* E->mesh.NOZ[level];
	  n2 = n1 + (E->mesh.NOY[level]-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
	  data[n2] = data[n1];
	}

   if(E->mesh.periodic_y && E->mesh.periodic_x)    /* then need to do the 1st one again */ 
     for(i=1;i<=E->mesh.NOZ[level];i++)
       for(j=1;j<=E->mesh.NOY[level];j++)
	 { n1 = i + (j-1)* E->mesh.NOX[level]*E->mesh.NOZ[level];
	   n2 = n1 + (E->mesh.NOX[level]-1)*E->mesh.NOZ[level];
	   data[n2] = data[n1];
	 }
  return;
}

void flogical_mesh_to_real(
			   struct All_variables *E,
			   standard_precision *data,
			   int level
			   )
{ 
  int i,j,n1,n2;

  if(E->mesh.periodic_x)
    for(i=1;i<=E->mesh.NOZ[level];i++)
      for(j=1;j<=E->mesh.NOY[level];j++) {
	n1 = i + (j-1)* E->mesh.NOX[level]*E->mesh.NOZ[level];
	n2 = n1 + (E->mesh.NOX[level]-1)*E->mesh.NOZ[level];
	data[n2] = data[n1];
      }

  if(E->mesh.periodic_y)
    for(i=1;i<=E->mesh.NOZ[level];i++)
      for(j=1;j<=E->mesh.NOX[level];j++)
	{ n1 = i + (j-1)* E->mesh.NOZ[level];
	  n2 = n1 + (E->mesh.NOY[level]-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
	  data[n2] = data[n1];
	}

   if(E->mesh.periodic_y && E->mesh.periodic_x)    /* then need to do the 1st one again */ 
     for(i=1;i<=E->mesh.NOZ[level];i++)
       for(j=1;j<=E->mesh.NOY[level];j++)
	 { n1 = i + (j-1)* E->mesh.NOX[level]*E->mesh.NOZ[level];
	   n2 = n1 + (E->mesh.NOX[level]-1)*E->mesh.NOZ[level];

	   data[n2] = data[n1];
	 }
  return;
}

standard_precision nodal_interpolated_value(
     struct All_variables *E,
     standard_precision *field,
     standard_precision x,
     standard_precision z,
     standard_precision y
)
{
 int found,element,count,level,i,j;
 standard_precision eta1,eta2,eta3;
 standard_precision lN[ELNMAX+1];
 standard_precision dOmega;
 standard_precision interpolant;
 
 const int dims = E->mesh.nsd;

 level=E->mesh.levmax;
 element=1;
 count=0;

 do {
   count++;
   found=1;
   general_element_coords(E,element,x,z,y,&eta1,&eta2,&eta3,level);
     
   if((eta1-1.0) > 1.0e-12) {  /* find in X direction */
      element += E->mesh.ELZ[level];
      found = 0;
   }
   else if((eta1+1.0) < -1.0e-12) {
     element -= E->mesh.ELZ[level];
     found = 0;
   }
	
   if((eta2-1.0) > 1.0e-12) {  /* find in Z direction */
     element += 1;
     found = 0;
   }
   else if((eta2+1.0) < -1.0e-12) {
     element -= 1;
     found = 0;
   }
   
   if(3==dims && (eta3-1.0) > 1.0e-12) {  /* find in Y direction */
     element += E->mesh.ELX[level] * E->mesh.ELZ[level];
     found = 0;
   }
   else if(3==dims && (eta3+1.0) < -1.0e-12) { /*RAA: 23/01/03, was missing '-' here*/
     element -= E->mesh.ELX[level] * E->mesh.ELZ[level];
     found = 0;
   }
   
   if(element < 1) 
     element = 1;
   else if(element > E->mesh.NEL[level])
     element = E->mesh.NEL[level];
   
 } while (found == 0 && count <= E->mesh.NEL[level]);
 
 if(!found) {
   fprintf(stderr,"%g,%g,%g is not in the domain !\n",x,z,y);
   exit(-1);
 }
   
 v_shape_fn(E,element,lN,eta1,eta2,eta3,level);
 interpolant=0.0;
 for(j=1;j<=enodes[dims];j++) {
  interpolant += field[E->IEN[level][element].node[j]] * lN[j];
 }
 
 /*
 fprintf(stderr,"%g %g %g -> %g\n",x,z,y,interpolant);
 */

 return(interpolant);

}


/* ===========================================================================
   Function to create the element locations from  the node positions. 
   ===========================================================================	*/
	
void element_locations_from_nodes(
     struct All_variables *E
)
{	
  int i,j,k,d,element,node;
  int elx,elz,ely,lev;
  int xhighnum,xlownum;
  int n1,n2,n3,n4,n5,n6,n7,n8;

  standard_precision xlowmean,xhighmean,xmean;
  standard_precision zlowmean,zhighmean;
  standard_precision ylowmean,yhighmean;
  standard_precision area;

  void inject_node_values();

  const int dims = E->mesh.nsd;

  for(i=1;i<=E->mesh.elx;i++)
    for(j=1;j<=E->mesh.elz;j++)
      for(k=1;k<=E->mesh.ely;k++)     {
	element =  j + (i-1)*E->mesh.elz + (k-1)*E->mesh.elz*E->mesh.elx;
	area = 1.0;
	for(d=1;d<=dims;d++) {

	  /* Get the mean of the coordinates in direction d, and
	     the mean of the differences from the centroid which give us
	     a measure of the size of the element */
	  
	  xmean = 0.0;
	  for(node=1;node<=enodes[dims];node++)  {
	    n1=E->ien[element].node[node];
	    /*RAA: 23/10/01, added !per_y part to line directly below*/
            if(E->mesh.periodic_x && !E->mesh.periodic_y && i == E->mesh.elx && loc[node].plus[0]==0)
	      n1 += (E->mesh.nox-1)*E->mesh.noz;
            /*RAA: 9/10/01, 24/10/01, added rest of statements below */
	    else if(E->mesh.periodic_y && !E->mesh.periodic_x && k == E->mesh.ely && (node==5 || node==6 || node==7 || node==8)) 
	      n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
	    else if(3==dims && E->mesh.periodic_x && E->mesh.periodic_y)  {
	       if(i == E->mesh.elx && k != E->mesh.ely && loc[node].plus[0]==0) 
	         n1 += (E->mesh.nox-1)*E->mesh.noz;
	       else if(k == E->mesh.ely && i != E->mesh.elx && (node==5 || node==6 || node==7 || node==8)) 
	         n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
	       else if(k == E->mesh.ely && i == E->mesh.elx) {
                 if(d==1 || d==2) {
                    if(node==3 || node==4 || node == 7 || node == 8)  
	               n1 += (E->mesh.nox-1)*E->mesh.noz;
                    else if(node==5 || node==6)  
	               n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
                 }
                 else if(d==3) {
                    if(node==3 || node==4)  
	               n1 += (E->mesh.nox-1)*E->mesh.noz;
                    else if(node==5 || node==6 || node == 7 || node == 8)  
	               n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
                 }
               }
            } /*end else if per_x and per_y*/
	    
	    xmean += E->x[d][n1];
	  }
	  
	  xmean /= (standard_precision) enodes[dims];
	  
	  xlowmean = xhighmean = 0.0;
	  xlownum = xhighnum = 0;
	  
	  for(node=1;node<=enodes[dims];node++)  {
	    n1=E->ien[element].node[node];
            /*RAA: 23/10/01, added !per_y part to line directly below*/
	    if(E->mesh.periodic_x && !E->mesh.periodic_y && i == E->mesh.elx && loc[node].plus[0]==0) 
	      n1 += (E->mesh.nox-1)*E->mesh.noz;
            /*RAA: 9/10/01, 24/10/01, added rest of statements below */
	    else if(E->mesh.periodic_y && !E->mesh.periodic_x && k == E->mesh.ely && (node==5 || node==6 || node==7 || node==8)) 
	      n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
	    else if(3==dims && E->mesh.periodic_x && E->mesh.periodic_y)  {
	       if(i == E->mesh.elx && k != E->mesh.ely && loc[node].plus[0]==0) 
	         n1 += (E->mesh.nox-1)*E->mesh.noz;
	       else if(k == E->mesh.ely && i != E->mesh.elx && (node==5 || node==6 || node==7 || node==8)) 
	         n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
	       else if(k == E->mesh.ely && i == E->mesh.elx) {
                 if(d==1 || d==2) {
                    if(node==3 || node==4 || node == 7 || node == 8)  
	               n1 += (E->mesh.nox-1)*E->mesh.noz;
                    else if(node==5 || node==6)  
	               n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
                 }
                 else if(d==3) {
                    if(node==3 || node==4)  
	               n1 += (E->mesh.nox-1)*E->mesh.noz;
                    else if(node==5 || node==6 || node == 7 || node == 8)  
	               n1 += (E->mesh.noy-1)*E->mesh.nox*E->mesh.noz;
                 }
               }
            } /*end else if per_x and per_y*/
	    
	    if (E->x[d][n1] <= xmean) {
	      xlownum++;
	      xlowmean += E->x[d][n1];
	    }
	    else		  {
	      xhighnum++;
	      xhighmean += E->x[d][n1];
	    }
	  }
	  
	  if( xlownum > 0)
	    xlowmean /= (standard_precision) xlownum;
	  else
	    xlowmean = xmean;
	  
	  if( xhighnum > 0)
	    xhighmean /= (standard_precision) xhighnum;
	  else
	    xhighmean = xmean;	      
	  E->eco[element].size[d] = 0.5 * (xhighmean-xlowmean);
	  E->eco[element].recip_size[d] = 1.0/E->eco[element].size[d];
	  E->eco[element].centre[d] = xmean;
	  area *= E->eco[element].size[d];
	}
	E->eco[element].area = area; 
	
	/* if(area < 0.0) */ 
	/* fprintf(stderr,"element %d - area = %g (%g x %g) (%g , %g)\n",element,area,
		E->eco[element].size[1],E->eco[element].size[2],
		E->eco[element].centre[1],E->eco[element].centre[2]	); */
      }/*end of the big loop*/

  /* now project them to lower levels where appropriate */
 	
  for(lev=E->mesh.levmax-1;lev>=E->mesh.levmin;lev--)    { 
    elx = E->mesh.ELX[lev];
    elz = E->mesh.ELZ[lev];
    ely = E->mesh.ELY[lev];
    
    for(i=1;i<=elx;i++)
      for(j=1;j<=elz;j++)
	for(k=1;k<=ely;k++) {
	  element = (k-1)*elx*elz + (i-1)*elz + j;
	  area = 1.0;
	  for(d=1;d<=dims;d++)	{ 
	     
	    /* Get The Mean Of The Coordinates In Direction D, And
	       The Mean Of The Differences From The Centroid Which Give Us
	       A Measure Of The Size Of The Element */
	    
	    xmean = 0.0;
	    for(node=1;node<=enodes[dims];node++)  {
	      n1=E->IEN[lev][element].node[node];
              /*RAA: 23/10/01, added !per_y part to line directly below*/
	      if(E->mesh.periodic_x && !E->mesh.periodic_y && i == E->mesh.ELX[lev] && loc[node].plus[0]==0) 
		n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
            /*RAA: 9/10/01, 24/10/01, added the rest of the statements below*/
	      else if(E->mesh.periodic_y && !E->mesh.periodic_x && k == E->mesh.ELY[lev] && (node==5 || node==6 || node==7 || node==8)) 
		n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
	      else if(3==dims && E->mesh.periodic_x && E->mesh.periodic_y)  {
	         if(i == E->mesh.ELX[lev] && k != E->mesh.ELY[lev] && loc[node].plus[0]==0) 
		    n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
	         else if(k == E->mesh.ELY[lev] && i != E->mesh.ELX[lev] && (node==5 || node==6 || node==7 || node==8)) 
		    n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
	         else if(k == E->mesh.ELY[lev] && i == E->mesh.ELX[lev]) {
                    if(d==1 || d==2) {
                       if(node==3 || node==4 || node==7 || node==8)  
		          n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
                       else if(node==5 || node==6)  
		          n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
                    }
                    else if(d==3) {
                       if(node==3 || node==4)  
		          n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
                       else if(node==5 || node==6 || node==7 || node==8)  
		          n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
                    }
                 }
              } /*end else if per_x and per_y*/
	      
	      xmean += E->X[lev][d][n1];
	    }
	    xmean /= (standard_precision) enodes[dims] ;

	    xlowmean = xhighmean = 0.0;
	    xlownum = xhighnum = 0;
	    
	    for(node=1;node<=enodes[dims];node++)  {
	      n1=E->IEN[lev][element].node[node];
              /*RAA: 23/10/01, added !per_y part to line directly below*/
	      if(E->mesh.periodic_x && !E->mesh.periodic_y && i == E->mesh.ELX[lev] && loc[node].plus[0]==0) 
		n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
            /*RAA: 9/10/01, 24/10/01, added the rest of the statements below*/
	      else if(E->mesh.periodic_y && !E->mesh.periodic_x && k == E->mesh.ELY[lev] && (node==5 || node==6 || node==7 || node==8)) 
		n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
	      else if(3==dims && E->mesh.periodic_x && E->mesh.periodic_y)  {
	         if(i == E->mesh.ELX[lev] && k != E->mesh.ELY[lev] && loc[node].plus[0]==0) 
		    n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
	         else if(k == E->mesh.ELY[lev] && i != E->mesh.ELX[lev] && (node==5 || node==6 || node==7 || node==8)) 
		    n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
	         else if(k == E->mesh.ELY[lev] && i == E->mesh.ELX[lev]) {
                    if(d==1 || d==2) {
                       if(node==3 || node==4 || node==7 || node==8)  
		          n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
                       else if(node==5 || node==6)  
		          n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
                    }
                    else if(d==3) {
                       if(node==3 || node==4)  
		          n1 += (E->mesh.NOX[lev]-1)*E->mesh.NOZ[lev];
                       else if(node==5 || node==6 || node==7 || node==8)  
		          n1 += (E->mesh.NOY[lev]-1)*E->mesh.NOX[lev]*E->mesh.NOZ[lev];
                    }
                 }
              } /*end else if per_x and per_y*/
	      
	      if (E->X[lev][d][n1] <= xmean) {
		xlownum++;
		xlowmean += E->X[lev][d][n1];
	      }
	      else		  {
		xhighnum++;
		xhighmean += E->X[lev][d][n1];
	      }
	    }
	      
	    if( xlownum > 0)
	      xlowmean /= (standard_precision) xlownum;
	    else
	      xlowmean = xmean;
	    
	      if( xhighnum > 0)
		xhighmean /= (standard_precision) xhighnum;
	      else
		xhighmean = xmean;	      
	    
	      E->ECO[lev][element].size[d] = 0.5 * (xhighmean-xlowmean);
	      E->ECO[lev][element].recip_size[d] = 1.0/E->ECO[lev][element].size[d];
	      E->ECO[lev][element].centre[d] = xmean;
	      area *= E->ECO[lev][element].size[d]; 
	      }

	      E->ECO[lev][element].area = area; 
	    }
       /* onto next level */
    } 
  /* Generate data for the advection scheme: element characteristic sizes and directions ... */
  
  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
    elx = E->mesh.ELX[lev];
    elz = E->mesh.ELZ[lev];
    ely = E->mesh.ELY[lev];
    
    for(i=1;i<=elx;i++)
      for(j=1;j<=elz;j++)
	for(k=1;k<=ely;k++) {
	  element = (k-1)*elx*elz + (i-1)*elz + j;

	  n1=E->IEN[lev][element].node[1];
	  n2=E->IEN[lev][element].node[2];
	  n3=E->IEN[lev][element].node[3];
	  n4=E->IEN[lev][element].node[4];

	  if(3==dims){
	    n5=E->IEN[lev][element].node[5];
	    n6=E->IEN[lev][element].node[6];
	    n7=E->IEN[lev][element].node[7];
	    n8=E->IEN[lev][element].node[8];
	  }

	  /*RAA: 24/10/01,  added !per_y to statement directly below */
          if(i==elx && E->mesh.periodic_x && !E->mesh.periodic_y) {
	    /* fprintf(stderr,"%d: el %d, node 3 %d ->%d [%d]\n",lev,element,
		    n3,n3+(E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev],E->mesh.NNO[lev]); */
	    n3 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	    n4 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];

	    if(3==dims) {
	      n7 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	      n8 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	    }
	  }
          /*RAA: 9/10/01, 24/10/01, added statements below */
	  else if(3==dims && E->mesh.periodic_y && !E->mesh.periodic_x && k==ely) {     
	    n5 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	    n6 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	    n7 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	    n8 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	  }
          /*RAA: 24/10/01, separate out this section for clarity */
	  else if(3==dims && E->mesh.periodic_y && E->mesh.periodic_x) {     
	     if(k!=ely && i==elx)  {    
	        n3 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n4 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n7 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n8 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
             }
	     else if(k==ely && i!=elx)  {    
	        n5 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
                n6 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	        n7 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	        n8 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
             }
	     else if(k==ely && i==elx)  {    
	        n3 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n4 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n5 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
                n6 += (E->mesh.NOY[lev]-1) * E->mesh.NOX[lev] * E->mesh.NOZ[lev];
	     /*   n7 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
	        n8 += (E->mesh.NOX[lev]-1) * E->mesh.NOZ[lev];
             */
             }
          } /*RAA: end of per_x and per_y addition*/
	  
	  /* 1: x direction */

	  if(2==dims) {
	    xlowmean =  (E->X[lev][1][n1] + 
			 E->X[lev][1][n2] )  * 0.5;
	    xhighmean = (E->X[lev][1][n3] + 
			 E->X[lev][1][n4] )  * 0.5;
	    zlowmean =  (E->X[lev][2][n1] + 
			 E->X[lev][2][n2] )  * 0.5;
	    zhighmean = (E->X[lev][2][n3] + 
			 E->X[lev][2][n4] )  * 0.5;


	    E->ECO[lev][element].ntl_dirns[1][1] = xhighmean - xlowmean;
	    E->ECO[lev][element].ntl_dirns[1][2] = zhighmean - zlowmean;

	    E->ECO[lev][element].ntl_size[1] =  sqrt(E->ECO[lev][element].ntl_dirns[1][1] *
						     E->ECO[lev][element].ntl_dirns[1][1] + 
						     E->ECO[lev][element].ntl_dirns[1][2] *
						     E->ECO[lev][element].ntl_dirns[1][2] ); 
						     
	    E->ECO[lev][element].ntl_recip_size[1] = 1.0 / E->ECO[lev][element].ntl_size[1] ;
	   
	    E->ECO[lev][element].ntl_dirns[1][1] *= E->ECO[lev][element].ntl_recip_size[1];
	    E->ECO[lev][element].ntl_dirns[1][2] *= E->ECO[lev][element].ntl_recip_size[1];

	    E->ECO[lev][element].ntl_size[1] *= 0.5;   /* consistent with other element definitions */
	    E->ECO[lev][element].ntl_recip_size[1] *= 2.0;

	  }
	  else {
	    xlowmean =  (E->X[lev][1][n1] + 
			 E->X[lev][1][n2] +
			 E->X[lev][1][n5] +
			 E->X[lev][1][n6] )  * 0.25;
	    xhighmean = (E->X[lev][1][n3] + 
			 E->X[lev][1][n4] +
			 E->X[lev][1][n7] +
			 E->X[lev][1][n8] )  * 0.25;
	    zlowmean =  (E->X[lev][2][n1] + 
			 E->X[lev][2][n2] +
			 E->X[lev][2][n5] +
			 E->X[lev][2][n6] )  * 0.25;
	    zhighmean = (E->X[lev][2][n3] + 
			 E->X[lev][2][n4] +
			 E->X[lev][2][n7] +
			 E->X[lev][2][n8] )  * 0.25;
	    ylowmean =  (E->X[lev][3][n1] + 
			 E->X[lev][3][n2] +
			 E->X[lev][3][n5] +
			 E->X[lev][3][n6] )  * 0.25;
	    yhighmean = (E->X[lev][3][n3] + 
			 E->X[lev][3][n4] +
			 E->X[lev][3][n7] +
			 E->X[lev][3][n8] )  * 0.25;

	    E->ECO[lev][element].ntl_dirns[1][1] = xhighmean - xlowmean;
	    E->ECO[lev][element].ntl_dirns[1][2] = zhighmean - zlowmean;
	    E->ECO[lev][element].ntl_dirns[1][3] = yhighmean - ylowmean;

	    E->ECO[lev][element].ntl_size[1] =  sqrt(E->ECO[lev][element].ntl_dirns[1][1] *
							   E->ECO[lev][element].ntl_dirns[1][1] + 
							   E->ECO[lev][element].ntl_dirns[1][2] *
							   E->ECO[lev][element].ntl_dirns[1][2] + 
							   E->ECO[lev][element].ntl_dirns[1][3] *
							   E->ECO[lev][element].ntl_dirns[1][3] ); 
						     
	    E->ECO[lev][element].ntl_recip_size[1] = 1.0 / E->ECO[lev][element].ntl_size[1] ;
	   
	    E->ECO[lev][element].ntl_dirns[1][1] *= E->ECO[lev][element].ntl_recip_size[1];
	    E->ECO[lev][element].ntl_dirns[1][2] *= E->ECO[lev][element].ntl_recip_size[1];
	    E->ECO[lev][element].ntl_dirns[1][3] *= E->ECO[lev][element].ntl_recip_size[1];

	    E->ECO[lev][element].ntl_size[1] *= 0.5;
	    E->ECO[lev][element].ntl_recip_size[1] *= 2.0;
	  }
	  /* 2: z direction */

	  if(2==dims) {
	    xlowmean =  (E->X[lev][1][n1] + 
			 E->X[lev][1][n4] )  * 0.5;
	    xhighmean = (E->X[lev][1][n2] + 
			 E->X[lev][1][n3] )  * 0.5;
	    zlowmean =  (E->X[lev][2][n1] + 
			 E->X[lev][2][n4] )  * 0.5;
	    zhighmean = (E->X[lev][2][n2] + 
			 E->X[lev][2][n3] )  * 0.5;

	    E->ECO[lev][element].ntl_dirns[2][1] = xhighmean - xlowmean;
	    E->ECO[lev][element].ntl_dirns[2][2] = zhighmean - zlowmean;

	    E->ECO[lev][element].ntl_size[2] =  sqrt(E->ECO[lev][element].ntl_dirns[2][1] *
						     E->ECO[lev][element].ntl_dirns[2][1] + 
						     E->ECO[lev][element].ntl_dirns[2][2] *
						     E->ECO[lev][element].ntl_dirns[2][2] ); 
						     
	    E->ECO[lev][element].ntl_recip_size[2] = 1.0 / E->ECO[lev][element].ntl_size[2] ;
	   
	    E->ECO[lev][element].ntl_dirns[2][1] *= E->ECO[lev][element].ntl_recip_size[2];
	    E->ECO[lev][element].ntl_dirns[2][2] *= E->ECO[lev][element].ntl_recip_size[2];

	    E->ECO[lev][element].ntl_size[2] *= 0.5;
	    E->ECO[lev][element].ntl_recip_size[2] *= 2.0;
	  }
	  else {
	    xlowmean =  (E->X[lev][1][n1] + 
			 E->X[lev][1][n4] +
			 E->X[lev][1][n5] +
			 E->X[lev][1][n8] )  * 0.25;
	    xhighmean = (E->X[lev][1][n2] + 
			 E->X[lev][1][n3] +
			 E->X[lev][1][n6] +
			 E->X[lev][1][n7] )  * 0.25;
	    zlowmean =  (E->X[lev][2][n1] + 
			 E->X[lev][2][n4] +
			 E->X[lev][2][n5] +
			 E->X[lev][2][n8] )  * 0.25;
	    zhighmean = (E->X[lev][2][n2] + 
			 E->X[lev][2][n3] +
			 E->X[lev][2][n6] +
			 E->X[lev][2][n7] )  * 0.25;
	    ylowmean =  (E->X[lev][3][n1] + 
			 E->X[lev][3][n4] +
			 E->X[lev][3][n5] +
			 E->X[lev][3][n8] )  * 0.25;
	    yhighmean = (E->X[lev][3][n2] + 
			 E->X[lev][3][n3] +
			 E->X[lev][3][n6] +
			 E->X[lev][3][n7] )  * 0.25;

	    E->ECO[lev][element].ntl_dirns[2][1] = xhighmean - xlowmean;
	    E->ECO[lev][element].ntl_dirns[2][2] = zhighmean - zlowmean;
	    E->ECO[lev][element].ntl_dirns[2][3] = yhighmean - ylowmean;

	    E->ECO[lev][element].ntl_size[2] =  sqrt(E->ECO[lev][element].ntl_dirns[2][1] *
						     E->ECO[lev][element].ntl_dirns[2][1] + 
						     E->ECO[lev][element].ntl_dirns[2][2] *
						     E->ECO[lev][element].ntl_dirns[2][2] + 
						     E->ECO[lev][element].ntl_dirns[2][3] *
						     E->ECO[lev][element].ntl_dirns[2][3] ); 
						     
	    E->ECO[lev][element].ntl_recip_size[2] = 1.0 / E->ECO[lev][element].ntl_size[2] ;
	   
	    E->ECO[lev][element].ntl_dirns[2][1] *= E->ECO[lev][element].ntl_recip_size[2];
	    E->ECO[lev][element].ntl_dirns[2][2] *= E->ECO[lev][element].ntl_recip_size[2];
	    E->ECO[lev][element].ntl_dirns[2][3] *= E->ECO[lev][element].ntl_recip_size[2];

	    E->ECO[lev][element].ntl_size[2] *= 0.5;
	    E->ECO[lev][element].ntl_recip_size[2] *= 2.0;
	  }
	  /* 3: y direction */

	  if(3==dims)  {
	    xlowmean =  (E->X[lev][1][n1] + 
			 E->X[lev][1][n2] +
			 E->X[lev][1][n3] +
			 E->X[lev][1][n4] )  * 0.25;
	    xhighmean = (E->X[lev][1][n5] + 
			 E->X[lev][1][n6] +
			 E->X[lev][1][n7] +
			 E->X[lev][1][n8] )  * 0.25;
	    zlowmean =  (E->X[lev][2][n1] + 
			 E->X[lev][2][n2] +
			 E->X[lev][2][n3] +
			 E->X[lev][2][n4] )  * 0.25;
	    zhighmean = (E->X[lev][2][n5] + 
			 E->X[lev][2][n6] +
			 E->X[lev][2][n7] +
			 E->X[lev][2][n8] )  * 0.25;
	    ylowmean =  (E->X[lev][3][n1] + 
			 E->X[lev][3][n2] +
			 E->X[lev][3][n3] +
			 E->X[lev][3][n4] )  * 0.25;
	    yhighmean = (E->X[lev][3][n5] + 
			 E->X[lev][3][n6] +
			 E->X[lev][3][n7] +
			 E->X[lev][3][n8] )  * 0.25;

	    E->ECO[lev][element].ntl_dirns[3][1] = xhighmean - xlowmean;
	    E->ECO[lev][element].ntl_dirns[3][2] = zhighmean - zlowmean;
	    E->ECO[lev][element].ntl_dirns[3][3] = yhighmean - ylowmean;

	    E->ECO[lev][element].ntl_size[3] =  sqrt(E->ECO[lev][element].ntl_dirns[3][1] *
						     E->ECO[lev][element].ntl_dirns[3][1] + 
						     E->ECO[lev][element].ntl_dirns[3][2] *
						     E->ECO[lev][element].ntl_dirns[3][2] + 
						     E->ECO[lev][element].ntl_dirns[3][3] *
						     E->ECO[lev][element].ntl_dirns[3][3] ); 
						     
	    E->ECO[lev][element].ntl_recip_size[3] = 1.0 / E->ECO[lev][element].ntl_size[3] ;
	   
	    E->ECO[lev][element].ntl_dirns[3][1] *= E->ECO[lev][element].ntl_recip_size[3];
	    E->ECO[lev][element].ntl_dirns[3][2] *= E->ECO[lev][element].ntl_recip_size[3];
	    E->ECO[lev][element].ntl_dirns[3][3] *= E->ECO[lev][element].ntl_recip_size[3];

	    E->ECO[lev][element].ntl_size[3] *= 0.5;
	    E->ECO[lev][element].ntl_recip_size[3] *= 2.0;

          /*RAA: 24/10/01, horrid, putrid little temporary fix for front right column of elements 
                by using the info from element 1, so the mesh must be uniform!*/
             if(3==dims && E->mesh.periodic_y && E->mesh.periodic_x && k==ely && i==elx)  {    
	        E->ECO[lev][element].ntl_size[1] = E->ECO[lev][1].ntl_size[1];
	        E->ECO[lev][element].ntl_size[2] = E->ECO[lev][1].ntl_size[2];
	        E->ECO[lev][element].ntl_size[3] = E->ECO[lev][1].ntl_size[3];
	        E->ECO[lev][element].ntl_recip_size[1] = E->ECO[lev][1].ntl_recip_size[1];
	        E->ECO[lev][element].ntl_recip_size[2] = E->ECO[lev][1].ntl_recip_size[2];
	        E->ECO[lev][element].ntl_recip_size[3] = E->ECO[lev][1].ntl_recip_size[3];
	        E->ECO[lev][element].ntl_dirns[1][1] = E->ECO[lev][1].ntl_dirns[1][1];
	        E->ECO[lev][element].ntl_dirns[1][2] = E->ECO[lev][1].ntl_dirns[1][2];
	        E->ECO[lev][element].ntl_dirns[1][3] = E->ECO[lev][1].ntl_dirns[1][3];
	        E->ECO[lev][element].ntl_dirns[2][1] = E->ECO[lev][1].ntl_dirns[2][1];
	        E->ECO[lev][element].ntl_dirns[2][2] = E->ECO[lev][1].ntl_dirns[2][2];
	        E->ECO[lev][element].ntl_dirns[2][3] = E->ECO[lev][1].ntl_dirns[2][3];
	        E->ECO[lev][element].ntl_dirns[3][1] = E->ECO[lev][1].ntl_dirns[3][1];
	        E->ECO[lev][element].ntl_dirns[3][2] = E->ECO[lev][1].ntl_dirns[3][2];
	        E->ECO[lev][element].ntl_dirns[3][3] = E->ECO[lev][1].ntl_dirns[3][3];
             }

	  } /*RAA: end of 'if' (3D)*/
	} /*RAA: end of 'k' loop*/

  /*RAA: data checking */ 
/*   if(E->control.verbose)  
     for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) 
       for (i=1; i<=E->mesh.NEL[lev];i++)  
         fprintf(stderr,"checking element %d, lev: %d,  area: %g,    ntl_size1: %g, ntl_size2: %g, ntl_size3: %g\n",i,lev,E->ECO[lev][i].area,E->ECO[lev][i].ntl_size[1],E->ECO[lev][i].ntl_size[2],E->ECO[lev][i].ntl_size[3]);
*/
/*
   if(E->control.verbose)  
       for (i=1; i<=E->mesh.nel;i++)  
         fprintf(stderr,"checking element %d,  area: %g,  size1: %g, size2: %g, size3: %g\n",i,E->eco[i].area,E->eco[i].size[1],E->eco[i].size[2],E->eco[i].size[3]);
*/

  } /*RAA: end of lev loop*/
 

  return;
}


/* This routine sorts the nodes into "colour" groupings which 
   allow (gs) relaxation of all nodes in one group simultaneously.
   This can happen because all nodes of one colour have no element
   in common with the others. In general, this means there are as
   many colours as nodes per element, I believe, at least for 
   the regular meshes built here */
   


void node_colouring(
		    struct All_variables *E
		    )
{
  int i,j,k,lv;
  int node,node2;
  int nno;
  int nox,noz,noy;
  int nox2,noz2,noy2;

  const int dims = E->mesh.nsd;

  int * node_list; 
  int node_count;

 
  /* For the rectangular mesh of linear elements, this is relatively easy */

  for(lv=E->mesh.levmax;lv>=E->mesh.levmin;lv--) {
   
    if(2==dims) { 
      M_node_cols = 9;
      M_node_col_num[lv] = (int *) Malloc0(M_node_cols * sizeof(int)); /* colours for this mesh */
      M_node_col_list[lv] = (int **) Malloc0(M_node_cols * sizeof(int *));
   
      nox = E->mesh.NOX[lv];
      noz = E->mesh.NOZ[lv];
      nno = E->mesh.NNO[lv];
      
      for(k=0;k<M_node_cols;k++) {
	M_node_col_num[lv][k] = 0;
	M_node_col_list[lv][k] = (int *) Malloc0((1+((nox+2)/3) * ((noz+2)/3) ) * sizeof(int)); 
      }
      

     
      for(i=1;i<=nox;i+=3) 
	for(j=1;j<=noz;j+=3) {
	  node = j + (i-1) * noz;

	  M_node_col_list[lv][0][++M_node_col_num[lv][0]] = node;
	  if(j+1 <= noz)
	    M_node_col_list[lv][1][++M_node_col_num[lv][1]] = node+1;
	  if(j+2 <= noz)
	    M_node_col_list[lv][2][++M_node_col_num[lv][2]] = node+2;
	  if(i+1 <= nox)
	    M_node_col_list[lv][3][++M_node_col_num[lv][3]] = node+noz;
	  if(i+1 <= nox && j+1 <= noz)
	    M_node_col_list[lv][4][++M_node_col_num[lv][4]] = node+noz+1;
	  if(i+1 <= nox && j+2 <= noz)
	    M_node_col_list[lv][5][++M_node_col_num[lv][5]] = node+noz+2;	  
	  if(i+2 <= nox)
	    M_node_col_list[lv][6][++M_node_col_num[lv][6]] = node+2*noz;
	  if(i+2 <= nox && j+1 <= noz)
	    M_node_col_list[lv][7][++M_node_col_num[lv][7]] = node+2*noz+1;
	  if(i+2 <= nox && j+2 <= noz)
	    M_node_col_list[lv][8][++M_node_col_num[lv][8]] = node+2*noz+2;
	}

    }

    /* And the same again for 3 dimensions */
    /*RAA - note that this has not been updated - no Mallocs here, etc - but even in
     * 2D it seems that there is a memory leak due to the lack of 'free' statements. */

    fprintf(stderr,"Level %d: total nodes %d, nodes by colour %d,%d,%d,%d,%d,%d,%d,%d,%d\n",
	    lv,E->mesh.NNO[lv],M_node_col_num[lv][0],M_node_col_num[lv][1],
	    M_node_col_num[lv][2],M_node_col_num[lv][3],M_node_col_num[lv][4],M_node_col_num[lv][5],
	    M_node_col_num[lv][6],M_node_col_num[lv][7],M_node_col_num[lv][8]);

    if(0 && lv<=E->mesh.levmin+1) {
      for(i=1;i<=M_node_col_num[lv][0];i++) 
	fprintf(stderr,"Colour 0, %d is node %d\n",i,M_node_col_list[lv][0][i]);
      for(i=1;i<=M_node_col_num[lv][1];i++) 
	fprintf(stderr,"Colour 1, %d is node %d\n",i,M_node_col_list[lv][1][i]);
      for(i=1;i<=M_node_col_num[lv][2];i++) 
	fprintf(stderr,"Colour 2, %d is node %d\n",i,M_node_col_list[lv][2][i]);
      for(i=1;i<=M_node_col_num[lv][3];i++) 
	fprintf(stderr,"Colour 3, %d is node %d\n",i,M_node_col_list[lv][3][i]);
      for(i=1;i<=M_node_col_num[lv][4];i++) 
	fprintf(stderr,"Colour 4, %d is node %d\n",i,M_node_col_list[lv][4][i]);
      for(i=1;i<=M_node_col_num[lv][5];i++) 
	fprintf(stderr,"Colour 5, %d is node %d\n",i,M_node_col_list[lv][5][i]);
      for(i=1;i<=M_node_col_num[lv][6];i++) 
	fprintf(stderr,"Colour 6, %d is node %d\n",i,M_node_col_list[lv][6][i]);
      for(i=1;i<=M_node_col_num[lv][7];i++) 
	fprintf(stderr,"Colour 7, %d is node %d\n",i,M_node_col_list[lv][7][i]);
      for(i=1;i<=M_node_col_num[lv][8];i++) 
	fprintf(stderr,"Colour 8, %d is node %d\n",i,M_node_col_list[lv][8][i]);
    }
  }  
  return;
}



/*RAA: 05/04/02, the version of mesh_update below is my version from
 * the 3D code, as opposed to the 2D Linux version from Dec 2001. For
 * this reason, there is no use of ACC variables as per that version, 
 * but instead I accomplish this by making sure time_interval is 'on'
 * in the input template, and the linear functions for time dependent
 * loading are handled via 'const' for the intercept and 'slope' for
 * the slope. Some more details are found in global_defs.h, etc. */

void mesh_update(
		 struct All_variables *E,
		 standard_precision timestep, 
		 standard_precision time
		 )

{  
  standard_precision new_X0,new_X1,length_scale_ratio;
  int i,lv;
  int k,found; /*RAA: 19/02/02*/
  int secx0,secx1,secz0,secz1,secy0,secy1;/*RAA: 19/02/02*/
  int seqx0,seqx1,seqz0,seqz1,seqy0,seqy1;/*RAA: 19/02/02*/

  struct RectBc RECT;
  void arbitrary_bc_rectangle();

  const int dims=E->mesh.nsd;

  /*RAA: 19/02/02, initialize these load sequence and trapezoid section markers */
  secx0 = secx1 = secz0 = secz1 = secy0 = secy1 = 0;   /*RAA: 19/02/02*/
  seqx0 = seqx1 = seqz0 = seqz1 = seqy0 = seqy1 = 0;  /*RAA: 19/02/02*/

  /* NB: This is only really appropriate for a cartesian mesh !!! */

 /* RAA: 12/02/02, added lots of lines for the problem of multiple-direction loading 
     and the stuff w/ time interval has no trig. functions, just const. or linear.
     However, I don't know why the linear term is:   timestep*b*(time-timestep)?
     if this is constant accel, shouldn't it be:        timestep*(b/2)*(2*time - timestep)?
     So there is a factor of 2 inserted by my calcs, ie.: timestep*b*(time-(timestep/2)) ?
     Note that this is also done in the linux version, which handles const acc 
     with new terms like BCvelocityaccX0, etc. 
    23/09/02, changed timestep to (timestep/2.0) in 6 places below, see comments above*/

  /* 1. X direction */

 if(!E->mesh.BCvelocity_time_interval) { /*RAA: 12/02/02, added this distinction*/
    new_X0 = E->mesh.layer0[1] + 
       timestep * (  E->mesh.BCvelocityX0 + E->mesh.BCvelocity1X0 * 
		  sin(2*M_PI*E->mesh.BCvelocityomegaX0 * (time - timestep) +
		  E->mesh.BCvelocityphaseX0*M_PI*2.0));

    new_X1 = E->mesh.layer1[1] + 
       timestep * (  E->mesh.BCvelocityX1 + E->mesh.BCvelocity1X1 * 
		  sin(2*M_PI*E->mesh.BCvelocityomegaX1 * (time - timestep) +
	          E->mesh.BCvelocityphaseX1*M_PI*2.0));
 }
 else if(E->mesh.BCvelocity_time_interval) {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityX0_time[0][i][k] && 
	    time <= E->mesh.BCvelocityX0_time[1][i+1][k])  {
            new_X0 = E->mesh.layer0[1] + timestep*(E->mesh.BCvelocityX0_const[i][k] +
	             E->mesh.BCvelocityX0_slope[i][k] * (time - (timestep/2.0))); 
	    seqx0 = k;  /*mark which load sequence, to skip loops below*/
	    secx0 = i;  /*mark which section of the trapezoid, to skip loops below*/
	    found = 1;
	    break;      /*or continue if X0 and X1 are present? */
	 }
         else  
            new_X0 = E->mesh.layer0[1];  /*implies velocity = 0*/
      }
    } while (found == 0 && k <= MAX_LOAD_SEQS);

    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityX1_time[0][i][k] &&
            time <= E->mesh.BCvelocityX1_time[1][i+1][k])  {
            new_X1 = E->mesh.layer1[1] + timestep*(E->mesh.BCvelocityX1_const[i][k] +
	             E->mesh.BCvelocityX1_slope[i][k] * (time - (timestep/2.0))); 
	    seqx1 = k; /*mark which load sequence, to skip loops below*/
	    secx1 = i; /*mark which section of the trapezoid, to skip loops below*/
	    found = 1;
	    break;
	 }
         else  
            new_X1 = E->mesh.layer1[1];  /*implies velocity = 0*/
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);
 }
  

  if(new_X1 < new_X0) {
    fprintf(stderr,"Mesh update has turned the mesh inside out (X) !\n");
    exit(-1);
  }
  
  length_scale_ratio = (new_X1-new_X0) / (E->mesh.layer1[1]-E->mesh.layer0[1]);
  
  for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
    for(i=1;i<=E->mesh.NNO[lv];i++) {
      E->SX[lv][1][i] = E->X[lv][1][i] = (E->X[lv][1][i] - E->mesh.layer0[1]) * length_scale_ratio + new_X0;
    }
  }
 
  E->mesh.layer0[1] = new_X0;
  E->mesh.layer1[1] = new_X1;
  
 /* 2. Z_direction */

 if(!E->mesh.BCvelocity_time_interval) { /*RAA: 12/02/02, added this distinction*/
    new_X0 = E->mesh.layer0[2] + 
             timestep * (  E->mesh.BCvelocityZ0 + E->mesh.BCvelocity1Z0 * 
	     sin(2*M_PI*E->mesh.BCvelocityomegaZ0 * (time - timestep) +
	     E->mesh.BCvelocityphaseZ0*M_PI*2.0));

    new_X1 = E->mesh.layer1[2] + 
             timestep * (  E->mesh.BCvelocityZ1 + E->mesh.BCvelocity1Z1 * 
             sin(2*M_PI*E->mesh.BCvelocityomegaZ1 * (time - timestep) +
	     E->mesh.BCvelocityphaseZ1*M_PI*2.0));
 }
 else if(E->mesh.BCvelocity_time_interval) {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityZ0_time[0][i][k] &&
            time <= E->mesh.BCvelocityZ0_time[1][i+1][k])  {
            new_X0 = E->mesh.layer0[2] + timestep*(E->mesh.BCvelocityZ0_const[i][k] +
	             E->mesh.BCvelocityZ0_slope[i][k] * (time - (timestep/2.0))); 
	    seqz0 = k;
	    secz0 = i;
	    found = 1;
	    break; 
	 }
         else  
            new_X0 = E->mesh.layer0[2];  /*implies velocity = 0*/
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);

    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityZ1_time[0][i][k] &&
            time <= E->mesh.BCvelocityZ1_time[1][i+1][k])  {
            new_X1 = E->mesh.layer1[2] + timestep*(E->mesh.BCvelocityZ1_const[i][k] +
	             E->mesh.BCvelocityZ1_slope[i][k] * (time - (timestep/2.0))); 
	    seqz1 = k;
	    secz1 = i;
	    found = 1;
	    break;
	 }
         else 
            new_X1 = E->mesh.layer1[2];  /*implies velocity = 0*/
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);
 } 
 
  if(new_X1 < new_X0) {
    fprintf(stderr,"Mesh update has turned the mesh inside out (Z) !\n");
    exit(-1);
  }
  
  length_scale_ratio = (new_X1-new_X0) / (E->mesh.layer1[2]-E->mesh.layer0[2]);
  
  for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
    for(i=1;i<=E->mesh.NNO[lv];i++) {
      E->SX[lv][2][i] = E->X[lv][2][i] = (E->X[lv][2][i] - E->mesh.layer0[2]) * length_scale_ratio + new_X0;
    }
  }
 
  E->mesh.layer0[2] = new_X0;
  E->mesh.layer1[2] = new_X1;

  
/* 3. Y_direction */
  
  if(3==E->mesh.nsd) {
    if(!E->mesh.BCvelocity_time_interval) { /*RAA: 12/02/02, added this distinction*/
       new_X0 = E->mesh.layer0[3] + 
          timestep * (  E->mesh.BCvelocityY0 + E->mesh.BCvelocity1Y0 * 
		  sin(2*M_PI*E->mesh.BCvelocityomegaY0 * (time - timestep) +
		  E->mesh.BCvelocityphaseY0*M_PI*2.0));

       new_X1 = E->mesh.layer1[3] + 
          timestep * (  E->mesh.BCvelocityY1 + E->mesh.BCvelocity1Y1 * 
		  sin(2*M_PI*E->mesh.BCvelocityomegaY1 * (time - timestep) +
		  E->mesh.BCvelocityphaseY1*M_PI*2.0));
    }
    else if(E->mesh.BCvelocity_time_interval) {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
       k = 0;
       do {
          found = 0;
          k++;
          for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
             if(time >= E->mesh.BCvelocityY0_time[0][i][k] &&
                time <= E->mesh.BCvelocityY0_time[1][i+1][k])  {
                new_X0 = E->mesh.layer0[3] + timestep*(E->mesh.BCvelocityY0_const[i][k] +
                         E->mesh.BCvelocityY0_slope[i][k] * (time - (timestep/2.0))); 
	        seqy0 = k;
	        secy0 = i;
	        found = 1;
  	        break; 
	     }
             else  
                new_X0 = E->mesh.layer0[3];  /*implies velocity = 0*/
          }
       } while (found == 0 && k <= MAX_LOAD_SEQS);
  
       k = 0;
       do {
          found = 0;
          k++;
          for(i=0; i<3; i++) {  /*find which section of the trapezoid contains the current time*/
             if(time >= E->mesh.BCvelocityY1_time[0][i][k] &&
                time <= E->mesh.BCvelocityY1_time[1][i+1][k])  {
                new_X1 = E->mesh.layer1[3] + timestep*(E->mesh.BCvelocityY1_const[i][k] +
	                 E->mesh.BCvelocityY1_slope[i][k] * (time - (timestep/2.0))); 
	        seqy1 = k;
	        secy1 = i;
	        found = 1;
	        break;
	     }
             else 
               new_X1 = E->mesh.layer1[3];  /*implies velocity = 0*/
          }
       } while (found == 0 && k <= MAX_LOAD_SEQS);
    }

   
    if(new_X1 < new_X0) {
      fprintf(stderr,"Mesh update has turned the mesh inside out (Y) !\n");
      exit(-1);
    }

    length_scale_ratio = (new_X1-new_X0) / (E->mesh.layer1[3]-E->mesh.layer0[3]);
    
    for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
      for(i=1;i<=E->mesh.NNO[lv];i++) {
	E->SX[lv][3][i] = E->X[lv][3][i] = (E->X[lv][3][i] - E->mesh.layer0[3]) * length_scale_ratio + new_X0;
      }
    }
 
    E->mesh.layer0[3] = new_X0;
    E->mesh.layer1[3] = new_X1;
  }

  /* Since the boundary condition may depend upon time, 
     also update the boundary condition values at the
     new boundary location */

  /*RAA: 13/02/02, here's the old part for updating the bc values. Keep it as #if 0
     so that it can be resurrected at a later time if needed - new stuff follows
      -- note that the 'omega' was missing in 4 places below, but i've corrected it*/
#if 0
   if(E->mesh.BCvelocityX0 != 0.0 || E->mesh.BCvelocity1X0 != 0.0 ) {
       RECT.numb=1; /* X-normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[1];
       RECT.mag[0]=E->mesh.BCvelocityX0 + E->mesh.BCvelocity1X0 * 
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaX0*time + 2.0*M_PI*E->mesh.BCvelocityphaseX0) ;
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
   }

  if(E->mesh.BCvelocityX1 != 0.0 || E->mesh.BCvelocity1X1 != 0.0) {
       RECT.numb=1; /* X-normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[1];
       RECT.mag[0]=E->mesh.BCvelocityX1 + E->mesh.BCvelocity1X1 * 
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaX1*time + 2.0*M_PI*E->mesh.BCvelocityphaseX1);

       /* fprintf(stderr,"Time = %g, Updating boundary condition at %g to %g\n",time,E->mesh.layer1[1],
	       E->mesh.BCvelocityX1 + E->mesh.BCvelocity1X1 * 
	       sin(2.0*M_PI*E->mesh.BCvelocityomegaX1*time + 2.0*M_PI*E->mesh.BCvelocityphaseX1)
	       );
       */

       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);

     }


   if(E->mesh.BCvelocityZ0 != 0.0 || E->mesh.BCvelocity1Z0 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[2];
       RECT.mag[0]=E->mesh.BCvelocityZ0 + E->mesh.BCvelocity1Z0 * 
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaZ0*time + 2.0*M_PI*E->mesh.BCvelocityphaseZ0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }

  if(E->mesh.BCvelocityZ1 != 0.0 || E->mesh.BCvelocity1Z1 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[2];
       RECT.mag[0]=E->mesh.BCvelocityZ1 + E->mesh.BCvelocity1Z1 *
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaZ1*time + 2.0*M_PI*E->mesh.BCvelocityphaseZ1);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }


   if(3==dims && (E->mesh.BCvelocityY0 != 0.0 || E->mesh.BCvelocity1Y0 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[3];
       RECT.mag[0]=E->mesh.BCvelocityY0 + E->mesh.BCvelocity1Y0 * 
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaY0*time + 2.0*M_PI*E->mesh.BCvelocityphaseY0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
    }

  if(3==dims && (E->mesh.BCvelocityY1 != 0.0 || E->mesh.BCvelocity1Y1 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[3];
       RECT.mag[0]=E->mesh.BCvelocityY1 + E->mesh.BCvelocity1Y1 * 
	 sin(2.0*M_PI*E->mesh.BCvelocityomegaY1*time + 2.0*M_PI*E->mesh.BCvelocityphaseY1);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
    }
#endif

/*-----------------------------------------------------*/
/*RAA: 13/02/02, here's the new part that checks what the time is before assuming
   that you want to apply a moving boundary to that face.
   The secx0, seqx0, etc are determined above to avoid needless looping here. */
/*But, for the time_interval=on stuff, there was a problem with the velocities
   not making it back to zero, so in order to fix it I re-do some of the 
   velocity assigning steps. It is not optimized and very likely ends up 
   re-calculating things that need not be, but at least it works, and can be
   straightened out at some point. */

 if(!E->mesh.BCvelocity_time_interval) {  
    if(E->mesh.BCvelocityX0 != 0.0 || E->mesh.BCvelocity1X0 != 0.0 ) {
       RECT.numb=1; /* X-normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[1];
       RECT.mag[0]=E->mesh.BCvelocityX0 + E->mesh.BCvelocity1X0 * 
	             sin(2.0*M_PI*E->mesh.BCvelocityomegaX0*time + 
                     2.0*M_PI*E->mesh.BCvelocityphaseX0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
    }
    if(E->mesh.BCvelocityX1 != 0.0 || E->mesh.BCvelocity1X1 != 0.0) {
       RECT.numb=1; /* X-Normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[1];
       RECT.mag[0]=E->mesh.BCvelocityX1 + E->mesh.BCvelocity1X1 * 
	             sin(2.0*M_PI*E->mesh.BCvelocityomegaX1*time +
	             2.0*M_PI*E->mesh.BCvelocityphaseX1);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
    }
    if(E->mesh.BCvelocityZ0 != 0.0 || E->mesh.BCvelocity1Z0 != 0.0) {
         RECT.numb=1; /* Z-normal  */ 
         RECT.norm[0]='Z';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer0[2];
         RECT.mag[0]=E->mesh.BCvelocityZ0 + E->mesh.BCvelocity1Z0 * 
	               sin(2.0*M_PI*E->mesh.BCvelocityomegaZ0*time +
	               2.0*M_PI*E->mesh.BCvelocityphaseZ0);
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }
    if(E->mesh.BCvelocityZ1 != 0.0 || E->mesh.BCvelocity1Z1 != 0.0) {
         RECT.numb=1; /* Z-normal  */ 
         RECT.norm[0]='Z';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer1[2];
         RECT.mag[0]=E->mesh.BCvelocityZ1 + E->mesh.BCvelocity1Z1 *
	               sin(2.0*M_PI*E->mesh.BCvelocityomegaZ1*time +
	               2.0*M_PI*E->mesh.BCvelocityphaseZ1);
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }
    if(3==dims) {
      if(E->mesh.BCvelocityY0 != 0.0 || E->mesh.BCvelocity1Y0 != 0.0) {
         RECT.numb=1; /* Y-normal  */ 
         RECT.norm[0]='Y';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer0[3];
         RECT.mag[0]=E->mesh.BCvelocityY0 + E->mesh.BCvelocity1Y0 * 
	               sin(2.0*M_PI*E->mesh.BCvelocityomegaY0*time +
		       2.0*M_PI*E->mesh.BCvelocityphaseY0);
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
      }
      if(E->mesh.BCvelocityY1 != 0.0 || E->mesh.BCvelocity1Y1 != 0.0) {
         RECT.numb=1; /* Y-normal  */ 
         RECT.norm[0]='Y';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer1[3];
         RECT.mag[0]=E->mesh.BCvelocityY1 + E->mesh.BCvelocity1Y1 * 
	               sin(2.0*M_PI*E->mesh.BCvelocityomegaY1*time +
		       2.0*M_PI*E->mesh.BCvelocityphaseY1);
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
      }
   } 
 } 
 else if(E->mesh.BCvelocity_time_interval) { /*RAA: eg, there are sequentially bi-directional moving boundaries, etc*/  
   /* ------------X0 stuff------------------*/
   if(0.0 != E->mesh.BCvelocityX0_const[secx0][seqx0] || 
             E->mesh.BCvelocityX0_slope[secx0][seqx0] != 0.0) {
           RECT.numb=1; /* X-normal  */ 
           RECT.norm[0]='X';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[1];

	   /*RAA: check for proper time again, just in case*/
           if(time >= E->mesh.BCvelocityX0_time[0][secx0][seqx0] &&
              time <= E->mesh.BCvelocityX0_time[1][secx0+1][seqx0])  
              RECT.mag[0]=E->mesh.BCvelocityX0_const[secx0][seqx0] +
		          E->mesh.BCvelocityX0_slope[secx0][seqx0] * time;
           else  /*Re-set mag to 0.0, but is never accessed, so see below*/
              RECT.mag[0]= 0.0;
    
           arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
   }
   if(0.0 == E->mesh.BCvelocityX0_const[secx0][seqx0] && 
             E->mesh.BCvelocityX0_slope[secx0][seqx0] == 0.0) {
           /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
           RECT.numb=1; /* X-Normal  */ 
           RECT.norm[0]='X';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[1];

           RECT.mag[0]= 0.0;
           arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
   }

   /* ------------X1 stuff------------------*/
   if(0.0 != E->mesh.BCvelocityX1_const[secx1][seqx1] || 
             E->mesh.BCvelocityX1_slope[secx1][seqx1] != 0.0) {
           RECT.numb=1; /* X-Normal  */ 
           RECT.norm[0]='X';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer1[1];

	   /*RAA: check for proper time again, just in case*/
           if(time >= E->mesh.BCvelocityX1_time[0][secx1][seqx1] &&
              time <= E->mesh.BCvelocityX1_time[1][secx1+1][seqx1]) { 
              RECT.mag[0]=E->mesh.BCvelocityX1_const[secx1][seqx1] +
		          E->mesh.BCvelocityX1_slope[secx1][seqx1] * time;
	   }
           else  /*Re-set mag to 0.0*/ /*But is never accessed, so see below*/
              RECT.mag[0]= 0.0;

           arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
   }
   if(0.0 == E->mesh.BCvelocityX1_const[secx1][seqx1] && 
             E->mesh.BCvelocityX1_slope[secx1][seqx1] == 0.0) {
           /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
           RECT.numb=1; /* X-Normal  */ 
           RECT.norm[0]='X';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer1[1];

           RECT.mag[0]= 0.0;
           arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
   }

   /* ------------Z0 stuff------------------*/
   if(0.0 != E->mesh.BCvelocityZ0_const[secz0][seqz0] || 
             E->mesh.BCvelocityZ0_slope[secz0][seqz0] != 0.0) {
           RECT.numb=1; /* Z-normal  */ 
           RECT.norm[0]='Z';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[2];

	   /*RAA: check for proper time again, just in case*/
           if(time >= E->mesh.BCvelocityZ0_time[0][secz0][seqz0] &&
              time <= E->mesh.BCvelocityZ0_time[1][secz0+1][seqz0])  
              RECT.mag[0]=E->mesh.BCvelocityZ0_const[secz0][seqz0] +
		          E->mesh.BCvelocityZ0_slope[secz0][seqz0] * time;
           else  /*Re-set mag to 0.0, but is never accessed, so see below*/
              RECT.mag[0]= 0.0;
    
           arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
   }
   if(0.0 == E->mesh.BCvelocityZ0_const[secz0][seqz0] && 
             E->mesh.BCvelocityZ0_slope[secz0][seqz0] == 0.0) {
    /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
           RECT.numb=1; /* Z-Normal  */ 
           RECT.norm[0]='Z';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[2];

           RECT.mag[0]= 0.0;
           arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
   }

   /* ------------Z1 stuff------------------*/
   if(0.0 != E->mesh.BCvelocityZ1_const[secz1][seqz1] || 
             E->mesh.BCvelocityZ1_slope[secz1][seqz1] != 0.0) {
           RECT.numb=1; /* Z-normal  */ 
           RECT.norm[0]='Z';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer1[2];

	   /*RAA: check for proper time again, just in case*/
           if(time >= E->mesh.BCvelocityZ1_time[0][secz1][seqz1] &&
              time <= E->mesh.BCvelocityZ1_time[1][secz1+1][seqz1])  
              RECT.mag[0]=E->mesh.BCvelocityZ1_const[secz1][seqz1] +
		          E->mesh.BCvelocityZ1_slope[secz1][seqz1] * time;
           else  /*Re-set mag to 0.0, but is never accessed, so see below*/
              RECT.mag[0]= 0.0;
    
           arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }
    if(0.0 == E->mesh.BCvelocityZ1_const[secz1][seqz1] && 
              E->mesh.BCvelocityZ1_slope[secz1][seqz1] == 0.0) {
           /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
           RECT.numb=1; /* Z-Normal  */ 
           RECT.norm[0]='Z';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer1[2];

           RECT.mag[0]= 0.0;
           arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
    }

    if(3==dims) {
       /* ------------Y0 stuff------------------*/
       if (0.0 != E->mesh.BCvelocityY0_const[secy0][seqy0] || 
                  E->mesh.BCvelocityY0_slope[secy0][seqy0] != 0.0) {
           RECT.numb=1; /* Y-normal  */ 
           RECT.norm[0]='Y';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[3];

	   /*RAA: check for proper time again, just in case*/
           if(time >= E->mesh.BCvelocityY0_time[0][secy0][seqy0] &&
              time <= E->mesh.BCvelocityY0_time[1][secy0+1][seqy0])  
              RECT.mag[0]=E->mesh.BCvelocityY0_const[secy0][seqy0] +
		          E->mesh.BCvelocityY0_slope[secy0][seqy0] * time;
           else  /*Re-set mag to 0.0, but is never accessed, so see below*/
              RECT.mag[0]= 0.0;
        
           arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
     }
     if (0.0 == E->mesh.BCvelocityY0_const[secy0][seqy0] && 
                E->mesh.BCvelocityY0_slope[secy0][seqy0] == 0.0) {
           /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
           RECT.numb=1; /* Y-Normal  */ 
           RECT.norm[0]='Y';
           RECT.bb1[0]= -1.0e32;
           RECT.bb2[0]=  1.0e32;
           RECT.aa1[0]= -1.0e32;
           RECT.aa2[0]=  1.0e32;

           RECT.intercept[0]=E->mesh.layer0[3];

           RECT.mag[0]= 0.0;
           arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
     }

     /* ------------Y1 stuff------------------*/
     if(0.0 != E->mesh.BCvelocityY1_const[secy1][seqy1] || 
               E->mesh.BCvelocityY1_slope[secy1][seqy1] != 0.0) {
             RECT.numb=1; /* Y-normal  */ 
             RECT.norm[0]='Y';
             RECT.bb1[0]= -1.0e32;
             RECT.bb2[0]=  1.0e32;
             RECT.aa1[0]= -1.0e32;
             RECT.aa2[0]=  1.0e32;
      
             RECT.intercept[0]=E->mesh.layer1[3];
  
	     /*RAA: check for proper time again, just in case*/
             if(time >= E->mesh.BCvelocityY1_time[0][secy1][seqy1] &&
                time <= E->mesh.BCvelocityY1_time[1][secy1+1][seqy1])  
                RECT.mag[0]=E->mesh.BCvelocityY1_const[secy1][seqy1] +
		            E->mesh.BCvelocityY1_slope[secy1][seqy1] * time;
             else  /*Re-set mag to 0.0, but is never accessed, so see below*/
                RECT.mag[0]= 0.0;
          
             arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
     }
     if(0.0 == E->mesh.BCvelocityY1_const[secy1][seqy1] && 
               E->mesh.BCvelocityY1_slope[secy1][seqy1] == 0.0) {
             /*RAA - new as of 27/09/02 - lets velocities get back to 0.0 */
             RECT.numb=1; /* Y-Normal  */ 
             RECT.norm[0]='Y';
             RECT.bb1[0]= -1.0e32;
             RECT.bb2[0]=  1.0e32;
             RECT.aa1[0]= -1.0e32;
             RECT.aa2[0]=  1.0e32;
  
             RECT.intercept[0]=E->mesh.layer1[3];
  
             RECT.mag[0]= 0.0;
             arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
     }
   }    /*end of if 3D */
 }      /*end of if time_interval is ON */

  /* And now propogate the new measurements through the
     element sizes etc */

  element_locations_from_nodes(E);
  mass_matrix(E);  /*RAA: 21/08/02, fix from LM as part of extension with temp bug*/

  return;
}

/*--------------------------------------------------------------------*/
/*RAA: 20/09/02 - new function, called from pg_solver and pertinent to time
   dependent loading problems with linearly stepped velocity functions.
   I.e., this routine handles 3-branch functions of the form: v(t) = a + bt */

void find_section_of_vel_fn (
                             struct All_variables *E,
			     int direction,
			     standard_precision time,
			     standard_precision *bcvel0,
			     standard_precision *bcvel1
			     )
{  

 int i,k,found;

/* 1. X_direction */
 if (direction==1) {  
    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityX0_time[0][i][k] && 
	    time <= E->mesh.BCvelocityX0_time[1][i+1][k])  {
            *bcvel0 = E->mesh.BCvelocityX0_const[i][k] + E->mesh.BCvelocityX0_slope[i][k] * time; 
          /*  fprintf(E->fp1,"*** bcvel0: %g %g %g\n",*bcvel0,time,E->mesh.BCvelocityX0_const[0][1]); */
	    found = 1;
	    break;      /*breaks from where? do loop or entire function????*/
	 }
         else 
            *bcvel0 = 0.0;
      }
    } while (found == 0 && k <= MAX_LOAD_SEQS);

    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityX1_time[0][i][k] &&
            time <= E->mesh.BCvelocityX1_time[1][i+1][k])  {
             *bcvel1 = E->mesh.BCvelocityX1_const[i][k] + E->mesh.BCvelocityX1_slope[i][k] * time; 
         /*   fprintf(E->fp1,"*** bcvel1: %g %g %g\n",*bcvel1,time,E->mesh.BCvelocityX1_const[0][1]); */
	    found = 1;
	    break;
	 }
         else 
            *bcvel1 = 0.0;
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);
 }

/* 2. Z_direction */
 
 else if (direction==2) {  
    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityZ0_time[0][i][k] &&
            time <= E->mesh.BCvelocityZ0_time[1][i+1][k])  {
            *bcvel0 = E->mesh.BCvelocityZ0_const[i][k] + E->mesh.BCvelocityZ0_slope[i][k] * time;
	    found = 1;
	    break; 
	 }
         else 
            *bcvel0 = 0.0;
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);

    k = 0;
    do {
       found = 0;
       k++;
       for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
         if(time >= E->mesh.BCvelocityZ1_time[0][i][k] &&
            time <= E->mesh.BCvelocityZ1_time[1][i+1][k])  {
            *bcvel1 = E->mesh.BCvelocityZ1_const[i][k] + E->mesh.BCvelocityZ1_slope[i][k] * time;
	    found = 1;
	    break;
	 }
         else 
            *bcvel1 = 0.0;
       }
    } while (found == 0 && k <= MAX_LOAD_SEQS);
 } 
 
/* 3. Y_direction */
  
 else if(3==E->mesh.nsd && direction==3) {  
       k = 0;
       do {
          found = 0;
          k++;
          for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
             if(time >= E->mesh.BCvelocityY0_time[0][i][k] &&
                time <= E->mesh.BCvelocityY0_time[1][i+1][k])  {
                *bcvel0 = E->mesh.BCvelocityY0_const[i][k] + E->mesh.BCvelocityY0_slope[i][k] * time;
	        found = 1;
  	        break; 
	     }
             else 
                *bcvel0 = 0.0;
          }
       } while (found == 0 && k <= MAX_LOAD_SEQS);
  
       k = 0;
       do {
          found = 0;
          k++;
          for(i=0; i<3; i++) {  /*find the section of the trapezoid w/ the current time*/
             if(time >= E->mesh.BCvelocityY1_time[0][i][k] &&
                time <= E->mesh.BCvelocityY1_time[1][i+1][k])  {
                *bcvel1= E->mesh.BCvelocityY1_const[i][k] + E->mesh.BCvelocityY1_slope[i][k] * time;
             /*   fprintf(E->fp1,"*** bcvelY1: %g %g %g\n",*bcvel1,time,E->mesh.BCvelocityY1_const[0][1]); */
	        found = 1;
	        break;
	     }
             else 
                *bcvel1 = 0.0;
          }
       } while (found == 0 && k <= MAX_LOAD_SEQS);
    }
}
