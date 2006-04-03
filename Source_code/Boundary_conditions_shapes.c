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

void arbitrary_bc_rectangle(
			    struct All_variables *E,
			    struct RectBc *RECT,
			    standard_precision **field,
			    unsigned int **bcbitf,
			    const unsigned int bcmask_on,
			    const unsigned int bcmask_off,
			    const unsigned int surf
			    )
{
  standard_precision weight,radius2,weight2;
  standard_precision aa1,bb1,area,mag;
  int i,j,k,locn;
  int inx,inz;
  int number,node,ndx,ndz;
  int level;
  const int dims = E->mesh.nsd; /*RAA: added this*/
  
  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
    for(number=0;number<RECT->numb;number++) {

      switch(RECT->norm[number]) {
      case 'X':
	for(i=1;i<=E->mesh.NOX[level];i++) {
	  if((i==1  && 
	      fabs(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
	      fabs(E->SX[level][1][1+i*E->mesh.NOZ[level]]-RECT->intercept[number])
	      ) ||
	     (i== E->mesh.NOX[level] &&
	      fabs(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
	      fabs(E->SX[level][1][1+(i-2)*E->mesh.NOZ[level]]-RECT->intercept[number])
	      ) ||
	     ((i!=1 && i != E->mesh.NOX[level]) &&
	      (fabs(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
	       fabs(E->SX[level][1][1+(i-2)*E->mesh.NOZ[level]]-RECT->intercept[number]) ) &&
	      (fabs(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
	       fabs(E->SX[level][1][1+i*E->mesh.NOZ[level]]-RECT->intercept[number])
	       )
	      )
	     ) {
	    locn = i;
	    if(E->control.verbose)
	      fprintf(stderr,"(X) Found bc: number:%d, icpt: %f at %f, node %d, lev: %d\n",number,RECT->intercept[number],E->X[level][1][i],i,level);
	    break;
	  }
	}
	for(k=1;k<=E->mesh.NOY[level];k++)
	  for(j=1;j<=E->mesh.NOZ[level];j++) {
	    node = j+(locn-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
            if(E->NODE[level][node] & ( OFFSIDE  )) /*RAA: 24/01/03, why weren't these 2 lines here before?*/
            continue;
	    aa1=E->SX[level][2][node];
	    bb1=(E->mesh.nsd !=3 ) ? 0.0 : E->SX[level][3][node];
	    inx=(aa1 >= RECT->aa1[number] && aa1 <= RECT->aa2[number]);
	    inz=(E->mesh.nsd != 3) ? 1 : (bb1 >= RECT->bb1[number] && bb1 <= RECT->bb2[number]);
 	    
	    if(inx && inz) {
	      if(field!=NULL) {
		if(surf) { /*RAA: surf is 0 for vel. bcs, 1 for stress bcs*/
		  if(j==1) {
		    area = fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(k==1) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else if(k == E->mesh.NOY[level]) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
	            }
		    if(3==dims && E->mesh.periodic_y) /*RAA: 15/01/01 added this bit, should be verified*/
	              area *= 2.0;
		  }
		  else if(j==E->mesh.NOZ[level]) {
		    area = fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(k==1) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else if(k == E->mesh.NOY[level]) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
	            }
		    if(3==dims && E->mesh.periodic_y) /*RAA: 15/01/01 added this bit, verify it!*/
		      area *= 2.0;
		  }
		  else {
		    area = fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. +
		      fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
		  /* to check in 3D because the node numbering might change */
		  /* if(3==dims)
		     area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]]) / 2. +
		     fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]]) / 2. ;*/
		      if(3==dims) {  /*RAA: 2-1-03, I think the following statements are correct */
		        if(k==1) 
		          area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
		        else if(k == E->mesh.NOY[level]) 
		          area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
		        else
		          area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                  fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
		      }
		  }
		  
		  mag = RECT->mag[number] * area ;
		  field[level][node] = RECT->mag[number] * area ;
		}
		else
		  field[level][node] = RECT->mag[number];
	      }   /*RAA: matches NULL brace*/
	      bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
	      bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
	    }
	  }
	break;
	
      case 'Z':
	locn=0;
	for(j=1;j<=E->mesh.NOZ[level];j++) {
	  if((j==1  && 
	      fabs(E->SX[level][2][j]-RECT->intercept[number]) < 
	      fabs(E->SX[level][2][j+1]-RECT->intercept[number])
	      ) ||
	     (j==E->mesh.NOZ[level] &&
	      fabs(E->SX[level][2][j]-RECT->intercept[number]) < 
	      fabs(E->SX[level][2][j-1]-RECT->intercept[number])
	      ) ||
	     ((j!=1 && j!=E->mesh.NOZ[level]) &&
	      (fabs(E->SX[level][2][j]-RECT->intercept[number]) < 
	       fabs(E->SX[level][2][j+1]-RECT->intercept[number]) ) &&
	      (fabs(E->SX[level][2][j]-RECT->intercept[number]) < 
	       fabs(E->SX[level][2][j-1]-RECT->intercept[number])
	       )
	      )
	     ) {
	    locn = j;
	    if(E->control.verbose)
	      fprintf(stderr,"(Z) Found bc: number:%d, icpt: %f at %f, node %d, lev: %d\n",number,RECT->intercept[number],E->X[level][2][j],j,level);
	    break;
	  }
	}

	/* NEEDS FIXING FOR Y */
	
	for(k=1;k<=E->mesh.NOY[level];k++)
	  for(i=1;i<=E->mesh.NOX[level];i++) {
	    node = locn+(i-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
	    if(E->NODE[level][node] & ( OFFSIDE  )) /*RAA: 24/01/03, why is this not in other 2 cases?*/
	      continue;
	    aa1=E->SX[level][1][node];
	    bb1=(E->mesh.nsd != 3) ? 0.0 : E->SX[level][3][node];
	    inx=(aa1 >= RECT->aa1[number] && aa1 <= RECT->aa2[number]);
	    inz=(E->mesh.nsd != 3) ? 1 : (bb1 >= RECT->bb1[number] && bb1 <= RECT->bb2[number]);
	    
	    if(inx && inz) {
	      if(field!=NULL) {
		if(surf) { /*RAA: surf is 0 for vel. bcs, 1 for stress bcs*/
		  if(i==1) {
		    area =  fabs( E->SX[level][1][node] - E->SX[level][1][node + E->mesh.NOZ[level]]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(k==1) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else if(k == E->mesh.NOY[level]) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		    }
		    if(E->mesh.periodic_x)
		      area *= 2. ;
		  }
		  else if(i==E->mesh.NOX[level]) {
		    area = fabs( E->SX[level][1][node] - E->SX[level][1][node - E->mesh.NOZ[level]]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(k==1) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else if(k == E->mesh.NOY[level]) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
	            }
		    if(E->mesh.periodic_x)
		      area *= 2. ;
		  }
		  else {
		    area = fabs( E->SX[level][1][node] - E->SX[level][1][node - E->mesh.NOZ[level]]) / 2. +
		      fabs( E->SX[level][1][node] - E->SX[level][1][node + E->mesh.NOZ[level]]) / 2. ;
		  /* to check in 3D because the node numbering might change */
		  /* if(3==dims)
		     area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]]) / 2. +
		     fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]]) / 2. ;*/
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(k==1) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else if(k == E->mesh.NOY[level]) 
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2.;
		      else
		        area *= fabs( E->SX[level][3][node] - E->SX[level][3][node - E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. +
		                fabs( E->SX[level][3][node] - E->SX[level][3][node + E->mesh.NOZ[level]*E->mesh.NOX[level]]) / 2. ;
	            }
		  }
		  
		  mag = RECT->mag[number] * area ;			  
		  field[level][node] = RECT->mag[number] * area ;
		}
		else /*RAA: velocity bcs*/
		  field[level][node] = RECT->mag[number];
	      }   /*RAA: matches NULL brace*/
	      bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
	      bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
	    }
	  }
	break;

      case 'Y':
	if(E->mesh.dof != 3)
	  break;
	
	for(k=1;k<=E->mesh.NOY[level];k++) {
      /* if(fabs(E->SX[level][3][1+(k-1)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) < 0.001/E->mesh.NOY[level])  */ /*RAA: 28/5/01, removed current line and substituted lines below for generality*/

          if((k==1  && 
              fabs(E->SX[level][3][1+(k-1)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
              fabs(E->SX[level][3][1+k*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number])
              ) ||
             (k==E->mesh.NOY[level] &&
              fabs(E->SX[level][3][1+(k-1)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
              fabs(E->SX[level][3][1+(k-2)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number])
              ) ||
            ((k!=1 && k != E->mesh.NOY[level]) &&
             (fabs(E->SX[level][3][1+(k-1)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
              fabs(E->SX[level][3][1+(k-2)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) ) && 
             (fabs(E->SX[level][3][1+(k-1)*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) < 
              fabs(E->SX[level][3][1+k*E->mesh.NOX[level]*E->mesh.NOZ[level]]-RECT->intercept[number]) 
              )
             )
            ) {
	    locn = k;
	  if(E->control.verbose)
            fprintf(stderr,"(Y) Found bc: number:%d, icpt: %f at %f, node %d, lev: %d\n",number,RECT->intercept[number],E->X[level][3][k],k,level);
	    break;
	  }
	}

	for(j=1;j<=E->mesh.NOZ[level];j++)
	  for(i=1;i<=E->mesh.NOX[level];i++) {
	    node = j+(i-1)*E->mesh.NOZ[level]+(locn-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
            if(E->NODE[level][node] & ( OFFSIDE  )) /*RAA: 24/01/03, why weren't these 2 lines here before?*/
              continue;
	    aa1=E->SX[level][1][node];
	    bb1=E->SX[level][2][node];
	    inx=(aa1 >= RECT->aa1[number] && aa1 <= RECT->aa2[number]);
	    inz=(bb1 >= RECT->bb1[number] && bb1 <= RECT->bb2[number]);
	    
	    /*RAA: 2/01/03 added surf stuff here but still must check periodics */
	    if(inx && inz) {
	      if(field!=NULL) {
		if(surf) { /*RAA: surf is 0 for vel. bcs, 1 for stress bcs*/
		  if(i==1) {
		    area =  fabs( E->SX[level][1][node] - E->SX[level][1][node + E->mesh.NOZ[level]]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(j==1) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
		      else if(j == E->mesh.NOZ[level]) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. ;
		      else
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. +
		                fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
	            }
		    if(E->mesh.periodic_x)
		      area *= 2. ;
		  }
		  else if(i==E->mesh.NOX[level]) {
		    area = fabs( E->SX[level][1][node] - E->SX[level][1][node - E->mesh.NOZ[level]]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(j==1) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
		      else if(j == E->mesh.NOZ[level]) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. ;
		      else
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. +
		                fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
	            }
		    if(E->mesh.periodic_x)
		      area *= 2. ;
		  }
		  else {
		    area = fabs( E->SX[level][1][node] - E->SX[level][1][node - E->mesh.NOZ[level]]) / 2. +
		           fabs( E->SX[level][1][node] - E->SX[level][1][node + E->mesh.NOZ[level]]) / 2. ;
		    if(3==dims) {  /*RAA: 2-1-03, I think the following needs to be added for 3D*/
		      if(j==1) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
		      else if(j == E->mesh.NOZ[level]) 
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. ;
		      else
		        area *= fabs( E->SX[level][2][node] - E->SX[level][2][node - 1]) / 2. +
		                fabs( E->SX[level][2][node] - E->SX[level][2][node + 1]) / 2. ;
	            }
		  }
		  
		  mag = RECT->mag[number] * area ;			  
		  field[level][node] = RECT->mag[number] * area ;
		}
		else /*RAA: velocity bcs*/
		  field[level][node] = RECT->mag[number];
	      }   /*RAA: matches NULL brace*/
	      bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
	      bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
	    }
	  }
	break;
      }  /*RAA: end of 'switch'*/
    }
      
  return;
}

void arbitrary_bc_circle(
    struct All_variables *E,
    struct CircBc *CIRC,
    standard_precision **field,
    unsigned int **bcbitf,
    unsigned int bcmask_on,
    const unsigned int bcmask_off,
    const unsigned int surf
)
{
    standard_precision weight,radius2,weight2;
    standard_precision aa1,bb1;
    int i,j,k,locn;
    int inx,inz;
    int number,node,ndx,ndz;
    int level;

    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
    for(number=0;number<CIRC->numb;number++) {
	switch(CIRC->norm[number]) {
	case 'X':
	    for(i=1;i<=E->mesh.NOX[level];i++) {
		if(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]] >= CIRC->intercept[number]) {
		    locn = i;
		    break;
		}
	    }
		for(k=1;k<=E->mesh.NOY[level];k++)
		    for(j=1;j<=E->mesh.NOZ[level];j++) {
			node = j+(locn-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
			aa1=E->SX[level][2][node];
			bb1=(E->mesh.nsd != 3) ? CIRC->bb[number] : E->SX[level][3][node];

			radius2=(aa1-CIRC->aa[number])*(aa1-CIRC->aa[number]) + (bb1-CIRC->bb[number])*(bb1-CIRC->bb[number]);
			
			if(radius2 <= CIRC->rad[number]*CIRC->rad[number]) {
			    if(field!=NULL)
				field[level][node] = CIRC->mag[number];
			    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
			}
		    }
	    break;

	case 'Z':
	    for(j=1;j<=E->mesh.NOZ[level];j++) {
		if(E->SX[level][2][j] >= CIRC->intercept[number]) {
		    locn = j;
		    break;
		}
	    }

	    for(k=1;k<=E->mesh.NOY[level];k++)
		for(i=1;i<=E->mesh.NOX[level];i++) {
		    node = locn+(i-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
		    aa1=E->SX[level][1][node];
		    bb1=(E->mesh.nsd != 3) ? CIRC->bb[number]  : E->SX[level][3][node];
		   
		    radius2=(aa1-CIRC->aa[number])*(aa1-CIRC->aa[number]) + (bb1-CIRC->bb[number])*(bb1-CIRC->bb[number]);
		
		    if(radius2 <= CIRC->rad[number]*CIRC->rad[number]) {
			if(field!=NULL)
			    field[level][node] = CIRC->mag[number];
			bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
		    }
		}
	    break;

	case 'Y':
	    if(E->mesh.dof != 3)
		break;
	    
	     for(k=1;k<=E->mesh.NOY[level];k++) {
		if(E->SX[level][3][j] >= CIRC->intercept[number]) {
		    locn = k;
		    break;
		}
	     }
		for(j=1;j<=E->mesh.NOZ[level];j++)
		    for(i=1;i<=E->mesh.NOX[level];i++) {
			node = j+(i-1)*E->mesh.NOZ[level]+(locn-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
			aa1=E->SX[level][1][node];
			bb1=E->SX[level][2][node];
			
			radius2=(aa1-CIRC->aa[number])*(aa1-CIRC->aa[number]) + (bb1-CIRC->bb[number])*(bb1-CIRC->bb[number]);
		
			if(radius2 <= CIRC->rad[number]*CIRC->rad[number]) {
			    if(field!=NULL)
				field[level][node] = CIRC->mag[number];
			    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
			}
		    }
	    break;	    
	}
    }
    return;
}

void arbitrary_bc_harmonic(
  struct All_variables *E,
  struct HarmBc *HARM,
  standard_precision **field,
  unsigned int **bcbitf,
  unsigned int bcmask_on,
  const unsigned int bcmask_off,
    const unsigned int surf
)
{
    standard_precision weight,radius2,weight2;
    standard_precision aa1,bb1;
    int i,j,k,locn,l;
    int inx,inz;
    int number,node,ndx,ndz;
    int level;

    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
    for(number=0;number<HARM->numb;number++) {
	switch(HARM->norm[number]) {
	case 'X':
            /* Find the column of nodes closest to the given intercept. NOTE: Assumes
               rectangular grid. */
	    for(i=1;i<=E->mesh.NOX[level];i++) {
		if(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]] >= HARM->intercept[number]) {
		    locn = i;
		    break;
		}
	    }
		for(k=1;k<=E->mesh.NOY[level];k++)
		    for(j=1;j<=E->mesh.NOZ[level];j++) {
			node = j+(locn-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
			aa1=E->SX[level][2][node];
			bb1=(E->mesh.nsd != 3) ? 0.0 : E->SX[level][3][node];
                        /* Check if we are within coordinate extents of BC. */
			inx=(aa1 >= HARM->aa1[number] && aa1 <= HARM->aa2[number]);
			inz=(E->mesh.nsd != 3) ? 1 : (bb1 >= HARM->bb1[number] && bb1 <= HARM->bb2[number]);
			
			if(inx && inz) {
			    if(field!=NULL) {
				field[level][node]=HARM->off[number];
				for(l=0;l<HARM->harms;l++) 
				    field[level][node] += HARM->amp[l][number] *
					cos((HARM->kaa[l][number]*aa1+HARM->phaa[l][number])*M_PI) *
					cos((HARM->kbb[l][number]*bb1+HARM->phbb[l][number])*M_PI) ;
			    }
			    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
			}
		    }
	    break;

	case 'Z':
            /* Find the row of nodes closest to the given intercept. */
	    for(j=1;j<=E->mesh.NOZ[level];j++) {
		if(E->SX[level][2][j] >= HARM->intercept[number]) {
		    locn = j;
		    break;
		}
	    }

	    for(k=1;k<=E->mesh.NOY[level];k++)
		for(i=1;i<=E->mesh.NOX[level];i++) {
		    node = locn+(i-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
		    aa1=E->SX[level][1][node];
		    bb1=(E->mesh.nsd != 3) ? 0.0 : E->SX[level][3][node];
		    inx=(aa1 >= HARM->aa1[number] && aa1 <= HARM->aa2[number]);
		    inz=(E->mesh.nsd != 3) ? 1 : (bb1 >= HARM->bb1[number] && bb1 <= HARM->bb2[number]);
		   
		    if(inx && inz) {
			if(field != NULL) {
			    field[level][node]=HARM->off[number];
			    for(l=0;l<HARM->harms;l++) 
				field[level][node] += HARM->amp[l][number] *
				    cos((HARM->kaa[l][number]*aa1+HARM->phaa[l][number])*M_PI) *
				    cos((HARM->kbb[l][number]*bb1+HARM->phbb[l][number])*M_PI) ;
			}
			bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
		    }
		}
	    break;

	case 'Y':
	    if(E->mesh.dof!=3) break;
	    
            /* Find the row of nodes closest to the given intercept. */
	     for(k=1;k<=E->mesh.NOY[level];k++) {
		if(E->SX[level][3][j] >= HARM->intercept[number]) {
		    locn = k;
		    break;
		}
	     }
		for(j=1;j<=E->mesh.NOZ[level];j++)
		    for(i=1;i<=E->mesh.NOX[level];i++) {
			node = j+(i-1)*E->mesh.NOZ[level]+(locn-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
			aa1=E->SX[level][1][node];
			bb1=E->SX[level][2][node];
			inx=(aa1 >= HARM->aa1[number] && aa1 <= HARM->aa2[number]);
			inz=(bb1 >= HARM->bb1[number] && bb1 <= HARM->bb2[number]);
			
			if(inx && inz) {
			    if(field != NULL) {
				field[level][node]=HARM->off[number];
				for(l=0;l<HARM->harms;l++) 
				    field[level][node] += HARM->amp[l][number] *
					cos((HARM->kaa[l][number]*aa1+HARM->phaa[l][number])*M_PI) *
					cos((HARM->kbb[l][number]*bb1+HARM->phbb[l][number])*M_PI) ;
			    }  
			    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
			}
		    }
	    break;	    
	}
    }
    return;
}
  
void arbitrary_bc_polynomial(
    struct All_variables *E,
    struct PolyBc *POLY,
    standard_precision **field,
    unsigned int **bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    const unsigned int surf
)
{
    standard_precision weight,radius2,weight2;
    standard_precision aa1,bb1;
    int i,j,k,locn,l;
    int inx,inz;
    int number,node,ndx,ndz;
    int level;
   
    /*RAA: 5/7/01, at some point the 'surf' should be added, as with RECT */
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
     for(number=0;number<POLY->numb;number++) {
	switch(POLY->norm[number]) {
	case 'X':
	    for(i=1;i<=E->mesh.NOX[level];i++) {
		if(E->SX[level][1][1+(i-1)*E->mesh.NOZ[level]] >= POLY->intercept[number]) {
		    locn = i;
		    if(E->control.verbose)
		      fprintf(stderr,"(X) Found bc in poly %d, %f at %f, node %d\n",number,POLY->intercept[number],E->X[level][1][i],i);
		    break;
		}
	    }
		for(k=1;k<=E->mesh.NOY[level];k++)
		    for(j=1;j<=E->mesh.NOZ[level];j++) {
			node = j+(locn-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
			aa1=E->SX[level][2][node];
			bb1=(E->mesh.nsd!=3) ? 0.0 : E->SX[level][3][node];
			inx=(aa1 >= POLY->aa1[number] && aa1 <= POLY->aa2[number]);
			inz=(E->mesh.nsd!=3) ? 1 : (bb1 >= POLY->bb1[number] && bb1 <= POLY->bb2[number]);
			
			if(inx && inz) {
			    if(field!=NULL){
				field[level][node]=0.0;
				for(l=0;l<=POLY->order;l++)  /*RAA, <=, not just < */
				    field[level][node] += 
					POLY->aaa[l][number]*pow(aa1,(double)l) + 
					POLY->abb[l][number]*pow(bb1,(double)l);
			    }
			    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
			}
		    }
	    break;

	case 'Z':
	    for(j=1;j<=E->mesh.NOZ[level];j++) {
		if(E->SX[level][2][j] >= POLY->intercept[number]) {
		    locn = j;
		    if(E->control.verbose)
		      fprintf(stderr,"(Z) Found bc in poly %d, %f at %f, node %d\n",number,POLY->intercept[number],E->X[level][2][j],j);
		    break;
		}
	    }

	    for(k=1;k<=E->mesh.NOY[level];k++)
		for(i=1;i<=E->mesh.NOX[level];i++) {
		    node = locn+(i-1)*E->mesh.NOZ[level]+(k-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
		    aa1=E->SX[level][1][node];
		    bb1=(E->mesh.nsd != 3) ? 0.0 : E->SX[level][3][node];
		    inx=(aa1 >= POLY->aa1[number] && aa1 <= POLY->aa2[number]);
		    inz=(E->mesh.nsd!=3) ? 1 : (bb1 >= POLY->bb1[number] && bb1 <= POLY->bb2[number]);
		   
		    if(inx && inz) {
			if(field!=NULL){
			    field[level][node]=0.0;
			    for(l=0;l<=POLY->order;l++) /*RAA, <=, not just < */
				field[level][node] += 
				    POLY->aaa[l][number]*pow(aa1,(double)l) + 
				    POLY->abb[l][number]*pow(bb1,(double)l); 	
			}
			bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
			bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
		    }
		}
	    break;

	case 'Y':
	    if(E->mesh.dof!=3)
		break;
	    
	     for(k=1;k<=E->mesh.NOY[level];k++) {
		if(E->SX[level][3][j] >= POLY->intercept[number]) {
		    locn = k;
		    if(E->control.verbose)
                      fprintf(stderr,"(Y) Found bc in poly %d, %f at %f, node %d\n",number,POLY->intercept[number],E->X[level][3][k],k);
		    break;
		}
	     }
	    
	    for(j=1;j<=E->mesh.NOZ[level];j++)
		for(i=1;i<=E->mesh.NOX[level];i++) {
		    node = j+(i-1)*E->mesh.NOZ[level]+(locn-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
		    aa1=E->SX[level][1][node];
		    bb1=E->SX[level][2][node];
		    inx=(aa1 >= POLY->aa1[number] && aa1 <= POLY->aa2[number]);
		    inz=(bb1 >= POLY->bb1[number] && bb1 <= POLY->bb2[number]);
		    
		    if(inx && inz) {
			if(field!=NULL){
			    field[level][node]=0.0;
			    for(l=0;l<=POLY->order;l++)   /*RAA, <=, not just < */
				field[level][node] += 
				    POLY->aaa[l][number]*pow(aa1,(double)l) + 
				    POLY->abb[l][number]*pow(bb1,(double)l); 	
			}  
		    bcbitf[level][node] = bcbitf[level][node]  | bcmask_on;
		    bcbitf[level][node] = bcbitf[level][node]  & (~bcmask_off);
		    }
		}   
	    break;	    
	}
     }
    return;
}
