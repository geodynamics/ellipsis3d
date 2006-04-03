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



#if (defined __sunos__)
#include <string.h>
#else
#if (!defined __GNUC__)
#include <strings.h> 
#endif
#endif

/*
#if (! defined __GNUC__)
#include <rpc/xdr.h> 
#endif
*/

#include "element_definitions.h"
#include "global_defs.h"

#include <math.h>


void get_tracer_elts_and_weights(
				 struct All_variables *E,
				 int M,
				 int N
				 )
{
  int level,i,j,k,m,n,q,el,tr,node,imax,ip,ii;
  int subel,near_t,trac,new;
  standard_precision subel_1;
  standard_precision xx,zz,yy,rad,radmin;
  standard_precision xxt,zzt,yyt;
  standard_precision *weight_fn[9];
  standard_precision quads;
  standard_precision eta1,eta2,eta3;
  standard_precision *ETA1[MAX_LEVELS], *ETA2[MAX_LEVELS], *ETA3[MAX_LEVELS];
  standard_precision dOmega;
  standard_precision el_length;
  standard_precision t1,t2,t3,dd;
  standard_precision tt1,tt2,tt3;
  standard_precision W1,W2,W3;
  int get_tracers_element();
  standard_precision get_tracer_jacobian();

  void tracer_within_boundaries();

  standard_precision subelt_dist[25][25];
  int subelt_tracer[25][25];

  void get_element_coords();
  void tracer_copy();
  int all_tracers_elts_and_sfns();

  standard_precision lN[ELNMAX+1];
  standard_precision time,CPU_time();
  
  standard_precision *mem1,*mem2;

  static int been_here = 0;
  static int trans = 0;

  const int SUBD = 9;
  int NUMT0;

  time = CPU_time();

  /* Need to eliminate any lost tracers now ... */
  
  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	  
    all_tracers_elts_and_sfns(E,E->tracer.tracer_elt[level],E->tracer.sfn_values[level],
			      E->tracer.tx,E->tracer.tz,E->tracer.ty,
			      E->tracer.eta1[level],E->tracer.eta2[level],E->tracer.eta3[level],
			      1,E->tracer.NUM_TRACERS,level);
  }
 
  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	  
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {	
      el = E->tracer.tracer_elt[level][m];    
      while (el == -1 && m <= E->tracer.NUM_TRACERS) {
	if(E->control.verbose)
	  fprintf(stderr,"Eliminating tracer %d which could not be found at %g,%g\n",m,E->tracer.tx[m],E->tracer.tz[m]);
	tracer_copy(E,E->tracer.NUM_TRACERS,m);
	E->tracer.NUM_TRACERS--;
	el = get_tracers_element(E,m,&eta1,&eta2,&eta3,level);
      }
    }
  }
 
  NUMT0 = E->tracer.NUM_TRACERS;  /* As we may be changing the number of tracers,
				     we'd best record the number we started with.
				     If we reduce the number of particles, then the loop should
				     stop at the new number. If we increase the number of particles
				     then it should stop at the old number so as not to resplit or anything
				     dangerous */


  /* SPLIT AND MERGE */

  if(E->advection.timesteps != 0) {
    if(E->control.verbose)
      fprintf(stderr,"Tracer split and merge\n");
    
    for(m=1;m<=min(E->tracer.NUM_TRACERS,NUMT0);m++) {
      /* Vow of chastity for certain groups */
      if(E->tracer.particle_reproduction[E->tracer.property_group[m]]) {
	tracer_black_widow(E,m);
      }
    }

    if(E->control.verbose)
      fprintf(stderr,"Tracer split and merge ... done\n");

    fprintf(E->fp,"Merging/Splitting of %d tracers: %g s \n",E->tracer.NUM_TRACERS,CPU_time()-time);
  }

  fprintf(E->fp,"Tracer merge/split: %g s \n",CPU_time()-time);

    /* Get all tracers' element and local coords, shape functions etc */

    /* NOTE - this should only really apply for tracers which have been
       messed with by the black widow */

   if(E->control.verbose)
      fprintf(stderr,"Tracer shape functions\n");

    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	  
      ii=all_tracers_elts_and_sfns(E,E->tracer.tracer_elt[level],E->tracer.sfn_values[level],
				E->tracer.tx,E->tracer.tz,E->tracer.ty,
				E->tracer.eta1[level],E->tracer.eta2[level],E->tracer.eta3[level],
				1,E->tracer.NUM_TRACERS,level);
    }

    if(ii != 0) {
      fprintf(stderr,"%d tracers were not found during shape function routine\n",ii);
    
      /* So we'd better do this again, then */

    /*RAA: 11/4/01, print statement below corrected for 3D*/
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {	  
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) {	
	el = E->tracer.tracer_elt[level][m];    
	while (el == -1 && m <= E->tracer.NUM_TRACERS) {
	  if(E->control.verbose && E->mesh.nsd==2)
	    fprintf(stderr,"** Eliminating tracer %d which could not be found at %g,%g\n",m,E->tracer.tx[m],E->tracer.tz[m]);
	  else if(E->control.verbose && E->mesh.nsd==3)
	    fprintf(stderr,"** Eliminating tracer %d which could not be found at %g,%g,%g\n",m,E->tracer.tx[m],E->tracer.tz[m],E->tracer.ty[m]);
	  tracer_copy(E,E->tracer.NUM_TRACERS,m);
	  E->tracer.NUM_TRACERS--;
	  el = get_tracers_element(E,m,&eta1,&eta2,&eta3,level);
	}
      }
    }
    }


  if(E->control.verbose)
      fprintf(stderr,"Tracer shape functions ... done\n");

  /* List all tracers for each element */

  fprintf(E->fp,"Tracer shape fns: %g s \n",CPU_time()-time);

  /* 1. Count 'em */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)
    for(el=1;el<=E->mesh.NEL[level];el++)
      E->tracer.tr_in_element_number[level][el] = 0;

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
      el = E->tracer.tracer_elt[level][m];
      E->tracer.tr_in_element_number[level][el]++;
    }
  }

  /* 2. Make room for the array */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
    E->tracer.tr_in_element_offset[level][1] = 0;
    for(el=2;el<=E->mesh.NEL[level];el++) {
      E->tracer.tr_in_element_offset[level][el] = 
	E->tracer.tr_in_element_offset[level][el-1] + E->tracer.tr_in_element_number[level][el-1];
      E->tracer.tr_in_element_number[level][el-1]=0; /* you'll see why in a moment */
    }
    E->tracer.tr_in_element_number[level][E->mesh.NEL[level]] = 0;
  }
 
  /* 3. Fill up the array as we go through the tracers one more time */

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
      el = E->tracer.tracer_elt[level][m];
      i = E->tracer.tr_in_element_offset[level][el] + E->tracer.tr_in_element_number[level][el];
      E->tracer.tr_in_element[level][i] = m;
      E->tracer.tr_in_element_number[level][el]++;
    }
  }

  fprintf(E->fp,"Sorting by element of tracers: %g s \n",CPU_time()-time);  
  fprintf(E->fp,"Computed local coords for tracers: %g s \n",CPU_time()-time);

  /* Find the representative volume for the tracers in each element */

  tracer_representative_volumes(E,
				E->tracer.eta1[E->mesh.levmax],
				E->tracer.eta2[E->mesh.levmax],
				E->tracer.eta3[E->mesh.levmax]); 
  fprintf(E->fp,"Reweighing of tracers: %g s \n",CPU_time()-time);

  /* For full accuracy, the volumes themselves do not produce a
     satisfactory integration scheme - it is preferable to allow
     the weightings to differ from the volumes to gain full accuracy */
   
  tracer_constrain_weightings(E,E->tracer.eta1,E->tracer.eta2,E->tracer.eta3);
  fprintf(E->fp,"Adjusted weights tracers: %g s \n",CPU_time()-time);

 /* Jacobians at all levels */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
      eta3 = E->tracer.eta3[level][m];

      E->tracer.tracer_jacobian[level][m] = 
	get_tracer_jacobian(E,m,eta1,eta2,eta3,level);
    }
  }

  /* Node masses from tracers */

  for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
    for(i=1;i<=E->mesh.NNO[level];i++)
      E->Trmass[level][i]=0.0;
    
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      el = E->tracer.tracer_elt[level][m];
        
      /* weight function for tracer from each node point, and
	 corresponding mass-matrix from tracers for each node */

      for(q=1;q<=enodes[E->mesh.nsd];q++) {
	node = E->IEN[level][el].node[q];
	E->Trmass[level][node] += E->tracer.sfn_values[level][m].node[q]
	  * E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];
      }
    }
  }

  /* Opportunity to transmute tracers */


#if 0
  if(E->monitor.elapsed_time > 0.15 && trans==0) {
    /*for(m=1;m<=E->tracer.NUM_TRACERS;m++) {		
      if(E->tracer.tx[m] > 1.0 && E->tracer.tx[m] < 2.5 && E->tracer.tz[m] > 0.5 && E->tracer.tz[m] < 2.5)
	E->tracer.property_group[m] = 2;
    }
    trans++;
    */

     E->mesh.BCvelocityX1=0.00001; 
     E->mesh.BCvelocityX0=0.00000000001; 
  
  }

#endif

  
  fprintf(E->fp,"Reweighing etc of tracers took %g s \n",CPU_time()-time); 

  return;
}

int get_tracers_element(
			struct All_variables *E,
			int m,
			standard_precision *e1,
			standard_precision *e2,
			standard_precision *e3,
			int level
			)
{
  void get_element_coords();
  void tracer_copy();

  standard_precision eta1,eta2,eta3;
  int found,element,count;
  
  const int dims = E->mesh.nsd;

  element = E->tracer.tracer_elt[level][m];

  if(element < 1 || element > E->mesh.NEL[level]) {
    /* fprintf(stderr,"%d: Tracer %d (%g,%g) thinks it is in element %d (%d)\n",level,m,
	    E->tracer.tx[m],E->tracer.tz[m],
	    element, E->mesh.NEL[level]); */
    element = E->mesh.NEL[level]/2;
  }

  count=0;
  do {
    count++;
    found=1;
    get_element_coords(E,element,m,E->tracer.tx,E->tracer.tz,
		       E->tracer.ty,&eta1,&eta2,&eta3,level);
		       
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
      else if(3==dims && (eta3+1.0) < -1.0e-12) { /*RAA: 16-01-03, was missing '-' sign here*/
	element -= E->mesh.ELX[level] * E->mesh.ELZ[level];
	found = 0;
      }

      if(element < 1) 
	element = E->mesh.NEL[level];

      else if(element > E->mesh.NEL[level])
	element = 1;


    } while (found == 0 && count <= E->mesh.NEL[level]);

  if (!found) {

    /* Try being completely methodical */

    /*RAA: 29/7/01, fix for 3D below, with eta3 */
    for(element=1;element <= E->mesh.NEL[level]; element++) {
       get_element_coords(E,element,m,E->tracer.tx,E->tracer.tz,
		       E->tracer.ty,&eta1,&eta2,&eta3,level);
       if(2==dims && eta1 >= -1.0 && eta1 <= 1.0 && eta2 >= -1.0 && eta2 <= 1.0) {
	 found = 1;
	 break;
       }
       if(3==dims && eta1 >= -1.0 && eta1 <= 1.0 && eta2 >= -1.0 && eta2 <= 1.0 && eta3 >= -1.0 && eta3 <= 1.0) {
	 found = 1;
         break;
       }
    }
  }

  if (!found) {
    /*RAA: 5/6/01, put 3D coords here, add verbose, and uncomment*/
    if(E->control.verbose && E->mesh.nsd==2)  
      fprintf(stderr,"Cannot find tracer %d: lev %d (%g,%g) element: %d etas: %g %g\n",m,level,E->tracer.tx[m],E->tracer.tz[m],E->tracer.tracer_elt[E->mesh.levmax][m],eta1,eta2);
    else if(E->control.verbose && E->mesh.nsd==3)  
      fprintf(stderr,"Cannot find tracer %d: lev %d (%g,%g,%g) element: %d etas: %g %g %g\n",m,level,E->tracer.tx[m],E->tracer.tz[m],E->tracer.ty[m],E->tracer.tracer_elt[E->mesh.levmax][m],eta1,eta2,eta3);
    return(-1);  /* not found ... don't change anything else */
  }

  /* send back the local coords if required */

  if(e1 != NULL)
    *e1 = eta1;
  if(e2 != NULL)
    *e2 = eta2;
  if(3==dims && e3 != NULL)
    *e3 = eta3;

  /* store local coords */

  E->tracer.eta1[level][m] = eta1;
  E->tracer.eta2[level][m] = eta2;
  if(3==dims)
     E->tracer.eta3[level][m] = eta3; /*RAA, fixed on 13/3/01 */

  /* store element */

  E->tracer.tracer_elt[level][m] = element;
  return(element);
}

/* Function to search the element space to find tracers 
   a home. Note, the current implementation assumes a
   regular mesh. In an irregular mesh, the element neighbours
   would need to be defined and searched */


int all_tracers_elts_and_sfns(
			      struct All_variables *E,
			      int *el_list,
			      struct TRACER_ELT_WEIGHT *lN,
			      standard_precision *tx,
			      standard_precision *tz,
			      standard_precision *ty,
			      standard_precision *eta1,
			      standard_precision *eta2,
			      standard_precision *eta3,
			      int N1,
			      int N2,
			      int level 
			      )
{ 
  const int elts = E->mesh.NEL[level];
  const int elx = E->mesh.ELX[level];
  const int elz = E->mesh.ELZ[level];
  const int ely = E->mesh.ELY[level];
  const int dims = E->mesh.nsd;

  standard_precision e1,e2,e3;

  int general_tracers_element();
  int get_tracers_element();

  int i,j,k,m,el;
  int *error;
  int check;
  int count;
 
  /* On a workstation, the vectorizing 
     makes for extemely inefficient operation.
     So let's bypass it if we can ! */

  if(!E->control.vector_optimization) {
    k=0;

    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {	
      el = get_tracers_element(E,m,&e1,&e2,&e3,level); /* also determines tracer.eta1,2,3  */
      E->tracer.tracer_elt[level][m] = el;

	/* fprintf(stderr,"Weights: Tracer %d, el %d, eta1,2=%g,%g\n",
			m,el,e1,e2); */

      if(el==-1) { /*RAA: 5/6/01, print coords*/
	if (3==dims)
          fprintf(stderr,"tracer %d (%g,%g,%g) (level %d) is lost ... element: %d\n",m,E->tracer.tx[m],E->tracer.tz[m],E->tracer.ty[m],level,E->tracer.tracer_elt[E->mesh.levmax][m]);
	else if (2==dims)
          fprintf(stderr,"tracer %d (%g,%g) (level %d) is lost ... element: %d\n",m,E->tracer.tx[m],E->tracer.tz[m],level,E->tracer.tracer_elt[E->mesh.levmax][m]);
	k++;
	continue;
      }
	
      /* v_shape_fn(E,el,E->tracer.sfn_values[level][m].node,e1,e2,e3,level); */
      v_shape_fn(E,el,lN[m].node,e1,e2,e3,level);

      if(2==dims) {
	for(i=1;i<=9;i++)
	  E->tracer.sfn_values[level][m].node[i] =  lN[m].node[i];
      }
      else {
	for(i=1;i<=27;i++)
	  E->tracer.sfn_values[level][m].node[i] =  lN[m].node[i];
      }

      /*RAA: 10,4,01, check out the shape function values*/
      /*   if(3==dims && m<=9) 
             if(E->control.verbose) 
               for(i=1;i<=11;i++)
                 fprintf(stderr,"**tracer.sfn_values: i, m, value  %d %d %f\n",i,m,lN[m].node[i]);*/

      eta1[m] = E->tracer.eta1[level][m] = e1;
      eta2[m] = E->tracer.eta2[level][m] = e2;
      if(3==dims)
	eta3[m] = E->tracer.eta3[level][m] = e3; /*RAA, fixed on 13/3/01 */
    }
    return(k);
  }

  printf("...bout to Malloc error in Tracer_interg  \n");

  error = (int *) Malloc0((N2 +1) * sizeof(int));

  check=0;
  for(m=N1;m<=N2;m++) {
    if(el_list[m] < 1 || el_list[m] > elts)
      check++;
  }

  if(check != 0) {
    fprintf(stderr,"Some tracers are lost\n");
  }

  /* Loop until all tracers have been found an element,
     or can safely be assumed to have got lost */

  count=0;
  do{ 
    /* Get the local coordinates for the elements which are
       the initial guesses */

    tr_local_coords(E,el_list,lN,tx,tz,ty,eta1,eta2,eta3,N1,N2,level);
 
    /* Now scan and determine a correction 
       value for each tracers element which can update it */
 
    check=0;
    for(m=N1;m<=N2;m++) {
      error[m] = 0;
    }
    
    for(m=N1;m<=N2;m++) {  /* Correction for X direction */
      if((eta1[m] - 1.0) > 1.0e-12) {
	error[m] += elz;
	check++;
      }
    }

    for(m=N1;m<=N2;m++) {
      if((eta1[m] + 1.0) < -1.0e-12) {
	error[m] -= elz;
	check++;
      }
    }

    for(m=N1;m<=N2;m++) {  /* Correction for Z direction */
      if((eta2[m] - 1.0) > 1.0e-12) {
	error[m] += 1;
	check++;
      }
    }
    
    for(m=N1;m<=N2;m++) {
      if((eta2[m] + 1.0) < -1.0e-12) {
	error[m] -= 1;
	check++;
      }
    }

    if(3==dims) {
      for(m=N1;m<=N2;m++) {  /* Correction for Y direction */
	if((eta3[m] - 1.0) > 1.0e-12) {
	  error[m] += elx*elz;
	  check++;
	}
      }
    
      for(m=N1;m<=N2;m++) {
	if((eta3[m] + 1.0) < -1.0e-12) {
	  error[m] -= elx*elz;
	  check++;
	}
      }
    }

    /* Now apply the offset values to the element
       list for all the tracers */

#pragma loop novrec el_list,error
    for(m=N1;m<=N2;m++) {
      el_list[m] += error[m];
    }
  
    for(m=N1;m<=N2;m++) {  /* e.g. for periodic things */
      if(el_list[m] <1)
	el_list[m] = elts;
      else if (el_list[m] > elts)
	el_list[m] = 1;
    }
    count++;

    /* fprintf(stderr,"Elt search, count = %d, non-found = %d\n",count,check); */

  } while (check != 0  && count <= 2 * (elx + ely + elz) );

  if(check != 0) {  /* need to flag an error condition for some elements */
#pragma loop novrec el_list,error
    for(m=N1;m<=N2;m++) {
      if(error[m] != 0)
	el_list[m] = -1;
    }
    fprintf(stderr,"WARNING: some tracers gone missing\n");
  }

#if 0
   for(m=N1;m<=N2;m++) {
     if(el_list[m] != (k=general_tracers_element(E,el_list[m],tx[m],tz[m],0,&e1,&e2,&e3,level)) ||
	fabs(e1 -eta1[m]) > 1.0e-7 ||
	fabs(e2 -eta2[m]) > 1.0e-7
	)
     fprintf(stderr,"V:%d in %d (%g,%g)\nS:%d in %d (%g,%g) [lev %d: %d/%d  %g,%g]\n",
	     m,el_list[m],eta1[m],eta2[m],m,k,e1,e2,level,error[m],count,tx[m],tz[m]);
 
   } 
#endif

  free((void *) error);
  return(check);
  }

int general_tracers_element(
			    struct All_variables *E,
			    int m,
			    standard_precision x,
			    standard_precision z,
			    standard_precision y,
			    standard_precision *e1,
			    standard_precision *e2,
			    standard_precision *e3,
			    int level
			    )
{
  void get_element_coords();

  standard_precision eta1,eta2,eta3;

  int found,element,count;
  
  const int dims = E->mesh.nsd;

  if(m > 0 && m < E->mesh.NEL[level])
    element=m;  /* initial guess */
  else 
    element = 1;

  count=0;
  do {
    count++;
    found=1;

    /* general_element_coords(E,element,x,z,y,&eta1,&eta2,&eta3,level); */
    get_element_coords(E,element,0,&x,&z,&y,&eta1,&eta2,&eta3,level); /* Should be OK */
    
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
    else if(3==dims && (eta3+1.0) < -1.0e-12) { /*RAA: 16-01-03, was missing '-' sign here*/
      element -= E->mesh.ELX[level] * E->mesh.ELZ[level];
      found = 0;
    }

    if(element < 1) 
      element = 1;
    else if(element > E->mesh.NEL[level])
      element = E->mesh.NEL[level];

    } while (found == 0 && count <= E->mesh.NEL[level]);

   if (0 && count >= E->mesh.NEL[level]) {
     if(2==dims) 
       fprintf(stderr,"%d: Cannot find tracer (%g,%g)\n",level,x,z); 
     else if(3==dims) 
       fprintf(stderr,"%d: Cannot find tracer (%g,%g,%g)\n",level,x,z,y); 
   }

  /* send back the local coords if required */

  if(!found) {
    if(2==dims)  /*RAA: add some statements*/
       fprintf(stderr,"element not found for coords (%g,%g) || etas (%g,%g)\n",x,z,eta1,eta2); 
    else if(3==dims) 
       fprintf(stderr,"element not found for coords (%g,%g,%g) || etas (%g,%g,%g)\n",x,z,y,eta1,eta2,eta3); 
    return(-1);
  }

  if(e1 != NULL)
    *e1 = eta1;
  if(e2 != NULL)
    *e2 = eta2;
  if(3==dims && e3 != NULL)
    *e3 = eta3;

  return(element);
}
