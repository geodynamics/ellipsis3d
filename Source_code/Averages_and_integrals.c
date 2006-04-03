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


void return_horiz_ave(
     struct All_variables *E,
     standard_precision *X,
     standard_precision *H
)
{
  int i,j,k,n,ln;

  standard_precision Have; 
  standard_precision *Z;
  standard_precision return_layer_value(); 

  Z = (standard_precision *) Malloc0((E->mesh.nsf+10)*sizeof(standard_precision));
   
  for(i=1;i<=E->mesh.noz;i++)  /* obtain average at each level */
    { for(j=1;j<=E->mesh.nox;j++)
	for(k=1;k<=E->mesh.noy;k++)
	  { n = i + (j-1)*E->mesh.noz + (k-1)*E->mesh.nox*E->mesh.noz;
	    ln = j + (k-1) * E->mesh.nox;
	    Z[ln] = X[n];
	  }
  
      Have = return_layer_value(E,Z,i,1);
      
      H[i] = Have; 

      /* next row */ 
    }
  free((void *) Z);
  return;
}

void return_horiz_rms(
     struct All_variables *E,
     standard_precision *X,
     standard_precision *H
)
{ int i,j,k,n,ln;
  standard_precision Have; 
  standard_precision *Z;
  standard_precision return_layer_value();

  Z = (standard_precision *) Malloc0((E->mesh.nsf+10)*sizeof(standard_precision));
   
  for(i=1;i<=E->mesh.noz;i++)  /* obtain average at each level */    {
    for(j=1;j<=E->mesh.nox;j++)
      for(k=1;k<=E->mesh.noy;k++) {
	n = i + (j-1)*E->mesh.noz + (k-1)*E->mesh.nox*E->mesh.noz;
	ln = j + (k-1) * E->mesh.nox;
	Z[ln] = X[n]*X[n];
      }
  
    Have = return_layer_value(E,Z,i,1);
      
    H[i] = sqrt(Have); 
    

      /* next row */ 
    }
  
  free((void *) Z);
  
  return;
}

standard_precision return_layer_value(      /* orthogonal meshes !!! */
     struct All_variables *E,
     standard_precision *Z,
     int layer,
     int average
)
{
    int j,k,d,nint,noz,nox,el,elz,elx,ely;
    int node[5],lnode[5];
    standard_precision Have,Hnorm;
    struct Shape_function1 M;
    struct Shape_function1_dA dGamma;
    void get_global_1d_shape_fn();
 

    if (!E->control.ORTHOZ) { 
	return(0.0);
    }

    noz = E->mesh.noz;
    nox = E->mesh.nox;
    elz = E->mesh.elz;
    elx = E->mesh.elx;
    ely = E->mesh.ely;

    if(E->control.SPHERE)
      return(0.0);


    Have = Hnorm = 0.0;

    for(j=1;j<=elx;j++)
	for(k=1;k<=ely;k++) {
	    node[1] = layer +(j-1)*noz + (k-1)*noz*nox;
	    node[2] = layer +(j) * noz + (k-1)*noz*nox;
	    node[4] = layer +(j-1)*noz + k*noz*nox; /* used if 3d */
	    node[3] = layer +(j) * noz + k*noz*nox;
	    lnode[1] = j + (k-1)*nox;
	    lnode[2] = j+1 + (k-1)*nox;
	    lnode[4] = j + k*nox; /* used if 3d */
	    lnode[3] = j+1 + k*nox;

	    if (layer == noz)
		el = layer-1 + (j-1)*elz + (k-1)*elx*elz;
	    else 
		el = layer + (j-1) * elz + (k-1)*elx*elz;
 
	    get_global_1d_shape_fn(E,el,&M,&dGamma,E->mesh.levmax);   
	
	    for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
		for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++){
		    Have += Z[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] * dGamma.vpt[GMVGAMMA(1,nint)];
		    Hnorm += E->M.vpt[GMVINDEX(d,nint)] * dGamma.vpt[GMVGAMMA(1,nint)]; } 

	/* Done one traverse */ 
	}

    if(average && Hnorm != 0.0) 
	Have /= Hnorm;


    return(Have); 
}
	   
standard_precision s_return_layer_value(      /* orthogonal meshes !!! */
     struct All_variables *E,
     standard_precision *Z,
     int layer,
     int average
)
{
    int j,k,d,nint,noz,nox,el,elz,elx,ely;
    int node[5],lnode[5];
    standard_precision Have,Hnorm;
    struct Shape_function1 M;
    struct Shape_function1_dA dGamma;
    void get_global_1d_shape_fn();


    if(E->control.SPHERE || E->control.CYLINDER)
      return(0.0);

    if (!E->control.ORTHOZ) { 
      fprintf(stderr,"Not able to calculate horizontal averages for non-horizontal layers\n");
      return(0.0);
    }

    noz = E->mesh.noz;
    nox = E->mesh.nox;
    elz = E->mesh.elz;
    elx = E->mesh.elx;
    ely = E->mesh.ely;
 
    Have = Hnorm = 0.0;

    for(j=1;j<=elx;j++)
	for(k=1;k<=ely;k++) {
	    node[1] = layer +(j-1)*noz + (k-1)*noz*nox;
	    node[2] = layer +(j) * noz + (k-1)*noz*nox;
	    node[4] = layer +(j-1)*noz + k*noz*nox; /* used if 3d */
	    node[3] = layer +(j) * noz + k*noz*nox;
	    lnode[1] = j + (k-1)*nox;
	    lnode[2] = j+1 + (k-1)*nox;
	    lnode[4] = j + k*nox; /* used if 3d */
	    lnode[3] = j+1 + k*nox;

	    if (layer == noz)
		el = layer-1 + (j-1)*elz + (k-1)*elx*elz;
	    else 
		el = layer + (j-1) * elz + (k-1)*elx*elz;
 
	    get_global_1d_shape_fn(E,el,&M,&dGamma,E->mesh.levmax);   
	
	    for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
		for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++){
		    Have += Z[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] * dGamma.vpt[GMVGAMMA(1,nint)];
		    Hnorm += E->M.vpt[GMVINDEX(d,nint)] * dGamma.vpt[GMVGAMMA(1,nint)]; } 

	/* Done one traverse */ 
	}
 
    if(average && Hnorm != 0.0) 
	Have /= Hnorm;
  
    return(Have); 
}
	   
standard_precision return_bulk_value(
     struct All_variables *E,
     standard_precision *Z,
     int average
)
{  
    void get_global_shape_fn();

    int i,j,k,n,el;
    double volume,integral;
   
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    if(E->control.SPHERE || E->control.CYLINDER)
      return(0.0);

    
    volume=0.0;
    integral=0.0;

    /*
    for(el=1;el<=E->mesh.nel;el++) {
	get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,E->mesh.levmax);
	
	for(j=1;j<=vpts;j++)
	    for(i=1;i<=ends;i++) {
		n = E->ien[el].node[i];
		volume += E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
		integral += Z[n] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
	    }
    }
    */

    for(i=1;i<=E->mesh.nno;i++) {
      volume += E->mass[i];
      integral += E->mass[i] * Z[i];
    }

    if(average && volume != 0.0)
	integral /= volume;

    return((standard_precision)integral);
}

standard_precision return_bulk_value_l(
				       struct All_variables *E,
				       standard_precision *Z,
				       int average,
				       int level
				       )
{  
  void get_global_shape_fn();
  
  int i,j,k,n,el;
  higher_precision volume,integral;
  
  struct Shape_function GN;
  struct Shape_function_dx GNx;
  struct Shape_function_dA dOmega;
  
  const int vpts = vpoints[E->mesh.nsd];
  const int ends = enodes[E->mesh.nsd];
  
  volume=0.0;
  integral=0.0;
  for(el=1;el<=E->mesh.NEL[level];el++) {
    get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,level);
    
    for(j=1;j<=vpts;j++)
      for(i=1;i<=ends;i++) {
	n = E->IEN[level][el].node[i];
	volume += E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
	integral += Z[n] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
      }
  }
  
  if(average && (volume != 0.0)) {
    integral /= volume;
  }
  
  return((standard_precision) integral);
}


void tracer_remove_horiz_ave(
			     struct All_variables *E,
			     standard_precision *field,
			     standard_precision *Have
			     )
{
  void gs_tracers_to_nodes();
  void return_horiz_ave();
  
  standard_precision interpolated_value_1d();

  standard_precision *Node,*NHave,*NZ;

  standard_precision average_value;

  int m,i;

  Node = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  NHave = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  NZ = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));

  /* Obtain nodal value of tracer field, and the average value */
 
  gs_tracers_to_nodes(E,Node,NULL,NULL,NULL,field,E->mesh.levmax,0);
  return_horiz_ave(E,Node,NHave);
 
  for(i=1;i<=E->mesh.noz;i++) {
    NZ[i] = E->x[2][i];
    /* fprintf(stderr,"vertical level %d, horiz ave %g at %g\n",i,NHave[i],NZ[i]); */
  }

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    
    average_value = interpolated_value_1d(E,NHave,NZ,E->mesh.noz,E->tracer.tz[m]); 
    
    /*
    fprintf(stderr,"REmove average from tracer %d (%d) at %g ... %g -= %g -> %g \n",m,
	    E->tracer.property_group[m],E->tracer.tz[m],
	    field[m],average_value,field[m]-average_value); 
    */
    field[m] -= average_value;

  }

  fprintf(stderr,"Have ... 3\n");

  free((void *) Node);
  free((void *) NHave);
  free((void *) NZ);

  return;
}

standard_precision interpolated_value_1d(
     struct All_variables *E,
     standard_precision *YY,
     standard_precision *XX,
     int NX,
     standard_precision X0
)
{

  int i,x0,x1;
  standard_precision interp;
  
  x1 = 1;

  while(X0 > XX[x1] && x1 < NX) 
    x1++;

  x0 = x1-1;

  if((XX[x1]-XX[x0]) != 0.0) 
    interp = YY[x0] + (YY[x1] - YY[x0]) * (X0 - XX[x0])/(XX[x1]-XX[x0]);
  else
    interp = 0.0;

  return(interp);

}


void remove_surf_horiz_press_ave(
  struct All_variables *E,
  standard_precision *PN,
  int level
)
{
  int i,j,k,z;
  int el,node;
  int nodei,nodej; /*RAA: 12/12/01, added this line*/
  standard_precision AveP,area,dx,dy;
  standard_precision minP;
  void sp_to_nodes();
  
  const int elx = E->mesh.ELX[level];
  const int elz = E->mesh.ELZ[level];
  const int ely = E->mesh.ELY[level];
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];

  const int dims = E->mesh.nsd; /*RAA: 29/11/01*/
  /* Run along the top and integrate the pressure/area */
  /* 2D case */ /*RAA: and now 3D, too, yay! */

  minP=1.0e32;
  area=0.0;
  AveP=0.0;
  if(dims==2) { /*RAA: added this line and corresp } */
    for(i=1;i<=nox;i++) {
      node = 1 + (i-1) * noz;
    
      if(PN[node] < minP) { 
        minP = PN[node];
	}
      
      if(i==1) {
	      dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node]);
      }
      else if (i==nox) { 
	      dx = 0.5 * (E->X[level][1][node] - E->X[level][1][node - noz]);
      }
      else  {
	      dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node - noz]);
      }
 
      area += dx;
      AveP += PN[node] * dx;
    }
  }
  else if(dims==3) { /*RAA: 29/11/01, 12/12/01, added this part for 3D, should be ok*/
     for(i=1;i<=nox;i++) 
        for(j=1;j<=noy;j++) {  
           nodei = 1 + (i-1) * noz;
           nodej = nodei + (j-1) * noz * nox;
    
           if(PN[nodej] < minP)
	     minP = PN[nodej];

           if(i==1)        
               dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei]);
           else if (i==nox)
               dx = 0.5 * (E->X[level][1][nodei] - E->X[level][1][nodei - noz]);
           else if (i!=1 && i!=nox) 
               dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei - noz]);

           if(j==1)        
               dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej]);
           else if (j==noy)
               dy = 0.5 * (E->X[level][3][nodej] - E->X[level][3][nodej - noz*nox]);
           else if (j!=1 && j!=noy)
               dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej - noz*nox]);
 
           area += dx*dy;
           AveP += PN[nodej] * dx * dy;
        }
  }
  
  if(area != 0.0) AveP /= area;

  /* Remove this pressure average from the entire field */
  
  for(i=1;i<=E->mesh.NNO[level];i++)  
    PN[i] -= AveP;

  /* Now do this again for each level to find the average pressure difference from the surface */

  if(dims==2) { /*RAA: added this line and corresp } */
    if(level==E->mesh.levmax) {
      for(z=1;z<=noz;z++) {
        area=0.0;
        AveP=0.0;
        for(i=1;i<=nox;i++) {
	  node = z + (i-1) * noz;
	  if(i==1)
	    dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node]);
	  else if (i==nox)
	    dx = 0.5 * (E->X[level][1][node] - E->X[level][1][node - noz]);
	  else
	    dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node - noz]);
   
	  area += dx;
	  AveP += PN[node] * dx;
        }
    
      if(area != 0.0)   E->Have.Pres[z] = AveP/area;
      }
    }
  }
  else if(dims==3) { /*RAA: 12/12/01, added this part for 3D, seems ok */
     if(level==E->mesh.levmax) {
       for(z=1;z<=noz;z++) {
          area=0.0;
          AveP=0.0;
          for(i=1;i<=nox;i++) {
             for(j=1;j<=noy;j++) {  
	        nodei = z + (i-1) * noz;
                nodej = nodei + (j-1) * noz * nox;

               if(i==1)
                  dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei]);
               else if (i==nox)
	          dx = 0.5 * (E->X[level][1][nodei] - E->X[level][1][nodei - noz]);
               else if (i!=1 && i!=nox) 
	          dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei - noz]);
  
	       if(j==1)
	          dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej]);
	       else if (j==noy)
	          dy = 0.5 * (E->X[level][3][nodej] - E->X[level][3][nodej - noz*nox]);
	       else if (j!=1 && j!=noy)
	          dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej - noz*nox]);
  
               area += dx*dy;
               AveP += PN[nodej] * dx * dy;
            }
          }

	  if(area != 0.0)   E->Have.Pres[z] = AveP/area;
       }
    }
  }

  return;
}


void remove_nodal_horiz_press_ave(
  struct All_variables *E,
standard_precision *PN,
int level
)

{
  int i,j,k,z;
  int el,node;
  int nodei,nodej; /*RAA: added this for 3D*/
  standard_precision AveP,area,dx,dy;
  const int dims = E->mesh.nsd; /*RAA: 13/12/01*/
  
  void sp_to_nodes();
  
  const int elx = E->mesh.ELX[level];
  const int elz = E->mesh.ELZ[level];
  const int ely = E->mesh.ELY[level];
  const int nox = E->mesh.NOX[level];
  const int noz = E->mesh.NOZ[level];
  const int noy = E->mesh.NOY[level];

  /* Run along each row and integrate the pressure/area */

  if(dims==2) { /*RAA: added this line and corresp } */
    for(z=1;z<=noz;z++) {
      area=0.0;
      AveP=0.0;
      for(i=1;i<=nox;i++) {
	node = z + (i-1) * noz;
	if(i==1)
	  dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node]);
	else if (i==nox)
	  dx = 0.5 * (E->X[level][1][node] - E->X[level][1][node - noz]);
	else
	  dx = 0.5 * (E->X[level][1][node + noz] - E->X[level][1][node - noz]);
 
	area += dx;
	AveP += PN[node] * dx;
      }
      if(area != 0.0)
	E->Have.Pres[z] = AveP/area;
      
      for(i=1;i<=nox;i++) {
	node = z + (i-1) * noz;
	PN[node] -= E->Have.Pres[z];
      }
    } /*end of 'z' loop*/
  }   /*end of if 2D */
  else if(dims==3) { /*RAA: 13/12/01, added this section */
    for(z=1;z<=noz;z++) {
      area=0.0;
      AveP=0.0;
      for(i=1;i<=nox;i++) {
        for(j=1;j<=noy;j++) {
          nodei = z + (i-1) * noz;
          nodej = nodei + (j-1) * noz * nox;

          if(i==1)
            dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei]);
          else if (i==nox)
            dx = 0.5 * (E->X[level][1][nodei] - E->X[level][1][nodei - noz]);
          else if (i!=1 && i!=nox)
            dx = 0.5 * (E->X[level][1][nodei + noz] - E->X[level][1][nodei - noz]);

          if(j==1)
            dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej]);
          else if (j==noy)
            dy = 0.5 * (E->X[level][3][nodej] - E->X[level][3][nodej - noz*nox]);
          else if (j!=1 && j!=noy)
            dy = 0.5 * (E->X[level][3][nodej + noz*nox] - E->X[level][3][nodej - noz*nox]);

          area += dx*dy;
          AveP += PN[nodej] * dx * dy;
        }
     }
     if(area != 0.0)
       E->Have.Pres[z] = AveP/area;

     for(i=1;i<=nox;i++) 
       for(j=1;j<=noy;j++) {
         nodei = z + (i-1) * noz;
         nodej = nodei + (j-1) * noz * nox;
         PN[nodej] -= E->Have.Pres[z];
       }

    } /*end of 'z' loop*/
  }   /*end of if 3D */

  return;
}
