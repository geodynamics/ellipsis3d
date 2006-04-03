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


#if (defined __sunos__) || defined(__GNUC__)
#include <string.h>
#else
#include <strings.h> 
#endif

#include <stdlib.h>  /* for drand48() */

/*
#if (! defined __GNUC__)
#include <rpc/xdr.h> 
#endif
*/

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

void tracer_initial_locations(
			      struct All_variables *E
			      )
{
  void get_tracer_elts_and_weights();
  void nodes_to_tracers();
  standard_precision magcross2d();
  
  int get_eq_phase();

  int *tr_density_element;

  int i,j,k,node,el,e;
  int kkj; /*RAA: added this loop variable */
  int p,q,r,element,counted_el; /*RAA: 30/10/01 these added for per_y & _x*/
  int m,level,tr;
  int fdims,fdatatypes,fmaterials;
  int file_error;

  int num_rect;
  int rect_tracers[40];
  int rect_tracer_colour[40];

  standard_precision rect_XX[40];
  standard_precision rect_tx1[40];
  standard_precision rect_tx2[40];
  standard_precision rect_tx3[40];
  standard_precision rect_tx4[40];
  standard_precision rect_tz1[40];
  standard_precision rect_tz2[40];
  standard_precision rect_tz3[40];
  standard_precision rect_tz4[40];
  standard_precision rect_ty1[40];
  standard_precision rect_ty2[40];
  standard_precision rect_ty3[40];
  standard_precision rect_ty4[40];

  int loc_or_mat; /*RAA: 11/07/02, added 3 variables for initial plastic strains*/
  int p_strain_mat_group[40];
  standard_precision p_strain_mag[40]; 


  standard_precision a1,a2;  /* triangle/tetrahedron  edge vectors */
  standard_precision b1,b2;
  standard_precision c1,c2;
  standard_precision d1,d2;
  standard_precision e1;

  standard_precision p1,p2;  /* test point to vertex vectors */
  standard_precision q1,q2;
  standard_precision r1,r2;

  standard_precision eta1,eta2,eta3;
  standard_precision lN[ELNMAX+1];

  standard_precision ref_area,area_sum,Dz;

  int testcases;
  standard_precision geom[100];
  int imat[100];
  int idx[100];
 
  float *fdata;

  FILE *fp, *fopen();
  char discard[5001];
  char fblockname[20][5];

  standard_precision recip_density;
  standard_precision maxx[4],minx[4];

  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];

  tr_density_element = (int *) Malloc0((1+E->mesh.nel) * sizeof(int));

  if(! E->tracer.TRACERS) 
    return;

 /* Read in locations from a file - Either

     this is a standard format ellipsis  binary file
     (e.g. from a previous run) and
     the requested data types must be read in */

  if((fp=fopen(E->tracer.tracer_input_file,"r")) != NULL) {

   file_error = 0;

    /* read in one line and determine if the file format is correct */
    fgets(discard,4999,fp);  
    if(strstr(discard,"FORMAT=binary") == 0) {
      fprintf(stderr,"Line 1: Particle input data file does not have the correct header information\n");
      file_error=1;
    }
    
    /* Read in next line for the number of data points available */

    if(!file_error) {
      fgets(discard,4999,fp);  
      i=sscanf(discard,"# PARTICLES=%d DIMENSIONS=%d",&(E->tracer.NUM_TRACERS),&fdims);
    
      if(i<1) {
	E->tracer.NUM_TRACERS = 0;
	fprintf(stderr,"Line 2: Particle input data file does not have the correct header information\n");
	file_error=1;
      }
      else {
	tracer_allocate_memory(E);
      }
    }

    if(fdims != dims) {
      fprintf(stderr,"Input particle file has the wrong number of dimensions (%d)\n",fdims);
      file_error=1;
    }

    /* Get the number of data elements */

    if(!file_error) {
      fgets(discard,4999,fp);  
      i=sscanf(discard,"# DATA_ELEMENTS=%d",&(fdatatypes));
    
      if(i<1) {
	E->tracer.NUM_TRACERS = 0;
	fprintf(stderr,"Line 3: Particle input data file does not have the correct header information\n");
	file_error=1;
      }
      fdata = (float *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(float));
    }
    
    /* Get the number of material elements */

    if(!file_error) {
      fgets(discard,4999,fp);  
      i=sscanf(discard,"# MATERIAL_ELEMENTS=%d",&(fmaterials));
    
      if(i<1) {
	E->tracer.NUM_TRACERS = 0;
	fprintf(stderr,"Line 4: Particle input data file does not have the correct header information\n");
	file_error=1;
      } 
      else {
	if(fmaterials > E->tracer.NUM_MATERIALS) {
	  fprintf(stderr,"Warning: The particle file contains more material types than specified\n");
	  fprintf(stderr,"         in the input paramters ... these will have to be truncated !\n");
	}
      }
    }

    /* Read in the data descriptions for the remaining lines of the header */

    if(!file_error) {

      for(i=0;i<=1+dims+fdatatypes;i++) {
	fgets(discard,4999,fp);  
	sscanf(discard,"# %4s=%d",fblockname[i],&k);
	fprintf(stderr,"Available data type: %d is %s (%d)\n",i,fblockname[i],k);
      }

      /* Skip one more line ... the ^L character */
      
      fgets(discard,4999,fp);  
    
      /* Check the test value (geo-PI) to see if the binary ordering is all wrong */

      fread(fdata,sizeof(float),1,fp);
      
      if(fabs(fdata[0]- 3.1) > 1.0e-6) {
	fprintf(stderr,"Binary file might need converting from different architecture (%f)\n",fdata[0]);
	file_error=1;
      }
    }

    /* Now start reading the data one block at a time */


    if(!file_error) {

      /* 1. the particle group numbers */ 
	
      fread(E->tracer.property_group,sizeof(int),E->tracer.NUM_TRACERS+1,fp);

      for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.property_group[i] > E->tracer.NUM_MATERIALS)
	  E->tracer.property_group[i] = E->tracer.NUM_MATERIALS;
      }

      /* 2. the X coordinates */ 
     
      fread(fdata,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
      for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
	E->tracer.tx[i] = (standard_precision) fdata[i];
	}     	

      /* 3. the Z coordinates */ 

      fread(fdata,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
      for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
	E->tracer.tz[i] = (standard_precision) fdata[i];
      }     	
   

      /* 4. the Y coord if 3 D */

      if(3==E->mesh.nsd) {
	fread(fdata,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
	  E->tracer.ty[i] = (standard_precision) fdata[i];
	}
      }

      /* 5. the tracer weights */ 

      /*RAA: 27/07/02 - seg fault within loop below for 3D restarting,
            but, this is almost certainly not where the real problem resides */
      fread(fdata,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
      for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
	E->tracer.tracer_weight_fn[E->mesh.levmax][i] = (standard_precision) fdata[i];
/*        fprintf(stderr,"RAA  num: %d  NUM_TRACERS: %d weight: %g\n",i,E->tracer.NUM_TRACERS,E->tracer.tracer_weight_fn[E->mesh.levmax][i]);*/
      }     	

      /* for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	 fprintf(stderr,"Particle %d (%d): (%g,%g) w=%g\n",i, 
	 E->tracer.property_group[i],E->tracer.tx[i],E->tracer.tz[i],E->tracer.weight[i]);		
	 } */

      /* 6. The remaining data - read and stored if required */ 

      for(k=1;k<=1+dims+fdatatypes;k++) {
	fread(fdata,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
      
	if(strstr(E->control.particle_data.which_data_types_in,fblockname[k])) {
	  for(j=1;j<=E->control.particle_data.numb;j++) {
	    if(strstr(E->control.particle_data.name[j],fblockname[k])) {

	      fprintf(stderr,"Reading %s info from block %d\n",fblockname[k],k);

	      for(i=0;i<=E->tracer.NUM_TRACERS;i++) {
		E->control.particle_data.data[j][i] = (standard_precision) fdata[i];
	      }
	      break;
	    }
	  }
	}
      }
    }
    
    free((void *) fdata);

    fclose(fp);
   
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      for(j=E->mesh.levmin;j<=E->mesh.levmax;j++) {
	E->tracer.tracer_elt[j][E->tracer.NUM_TRACERS] = 1;
      }
      E->tracer.tracer_std_weighting[E->mesh.levmax][i] = E->tracer.tracer_weight_fn[E->mesh.levmax][i]; 
      
      E->tracer.dX11[i] = 0.5 * (E->x[1][E->mesh.nno] - E->x[1][1])  / E->mesh.ELX[E->mesh.levmax]; /* fix this later !!! */
      E->tracer.dX22[i] = 0.5 * (E->x[2][E->mesh.nno] - E->x[2][1])  / E->mesh.ELZ[E->mesh.levmax];  /* fix this later !!! */
    }
  }

  /* Otherwise, if no previous file is available, generate a
     background distribution of tracers on which the initial 
     conditions can be imposed */
  /* RAA: 19/3/01, N.B., for 3D problems, only this part of the tracer
     distribution is corrected (still not finished!), not the previous part which
     is based on initial tracer distributions from a separate file */

  else {

    /* Initialize the tracer density for each element */

    for(el=1;el<=E->mesh.nel;el++) {
      tr_density_element[el] = 1;
    }

    /* Get the rectangular areas for density of tracers. */

    input_int("Tracer_rect",&num_rect,"0");
    input_int_vector("Tracer_rect_density",num_rect,rect_tracers); 
    input_std_precision_vector("Tracer_rect_x1",num_rect,rect_tx1);
    input_std_precision_vector("Tracer_rect_x2",num_rect,rect_tx2);
    input_std_precision_vector("Tracer_rect_z1",num_rect,rect_tz1);
    input_std_precision_vector("Tracer_rect_z2",num_rect,rect_tz2);
    input_std_precision_vector("Tracer_rect_y1",num_rect,rect_ty1);
    input_std_precision_vector("Tracer_rect_y2",num_rect,rect_ty2);
    
    /* And scan them */

    for(el=1;el<=E->mesh.nel;el++) {
      for(j=1;j<=dims;j++) {
	maxx[j] = -1.0e32;
	minx[j] =  1.0e32;
      }
      
      for(j=1;j<=dims;j++)
	for(i=1;i<=ends;i++) {
	  if(maxx[j] < E->x[j][E->ien[el].node[i]]) 
	    maxx[j] = E->x[j][E->ien[el].node[i]];
	  if(minx[j] > E->x[j][E->ien[el].node[i]]) 
	    minx[j] = E->x[j][E->ien[el].node[i]];
	}
      
      for(i=0;i<num_rect;i++) {

	/* does this element overlap this 
	   rectangular region at all ? */

	if(dims==2) { /*RAA */
	  if(minx[1] < rect_tx2[i] && maxx[1] >  rect_tx1[i] &&
	     minx[2] < rect_tz2[i] && maxx[2] >  rect_tz1[i]) {
	 
	    tr_density_element[el] = /* max(tr_density_element[el],min(rect_tracers[i],12)); */
	        min(rect_tracers[i],12);
	  }
	} 
	else {   /*RAA: 3D  21/3/01 */
          if(minx[1] < rect_tx2[i] && maxx[1] >  rect_tx1[i] &&
	     minx[2] < rect_tz2[i] && maxx[2] >  rect_tz1[i] &&
	     minx[3] < rect_ty2[i] && maxx[3] >  rect_ty1[i]) {

	     tr_density_element[el] = /* max(tr_density_element[el],min(rect_tracers[i],12)); */
	         min(rect_tracers[i],12);
	  }
        }
      }  /*RAA: end of 'i' loop */
    }    /*RAA: end of 'el' loop */
    
    E->tracer.NUM_TRACERS=0;
    
/*RAA: 30/5/01, check which nodes are periodic */
/*    if(E->control.verbose) 
       for(e=1;e<=ends;e++) {
         node = E->ien[el].node[e];
         if((E->node[node] & PER_OFFSIDE)) 
          fprintf(stderr,"(Tracer_read), following node is PER_OFFSIDE:local node, node, E->node[node]: %d %d %d\n",e,node,E->node[node]);
       }
    }
*/
    
    /*RAA:  a big loop over the elements, with 2D and 3D separated via 'if'  */
    for(el=1;el<=E->mesh.nel;el++) {
      if(tr_density_element[el] == 0) {
	continue;
      }


      recip_density = 1.0 / tr_density_element[el];
      if(2==E->mesh.nsd) {
        for(i=0;i<tr_density_element[el];i++)
	  for(j=0;j<tr_density_element[el];j++) {
	    E->tracer.NUM_TRACERS++;
	    tracer_allocate_memory(E);

/* Depending on where you want your integration points */
#if 0
	    if(tr_density_element[el] == 2) {
	      eta1 = 0.577350269189626 * (2*i-1) ;
	      eta2 = 0.577350269189626 * (2*j-1) ;
	    }
	    if(tr_density_element[el] == 3) {
	      eta1 = 0.774596669241483 * (i-1) ;
	      eta2 = 0.774596669241483 * (j-1) ;
	    }
	    if(tr_density_element[el] == 4) {
	      eta1 = (0.861136311594053/3.-0.339981043584856) * i*i*i + (4.5*0.339981043584856 - 1.5*0.861136311594053)*i*i
	        +(13*0.861136311594053/6-4.5*0.339981043584856)*i - 0.861136311594053;
	      eta2 = (0.861136311594053/3.-0.339981043584856) * j*j*j + (4.5*0.339981043584856 - 1.5*0.861136311594053)*j*j
	        +(13*0.861136311594053/6-4.5*0.339981043584856)*j - 0.861136311594053;
	    }
#else
	    eta1 = ((2.0 * i + 1.0 + E->control.particle_perturb * (0.5 - drand48())) * recip_density - 1.0);
	    eta2 = ((2.0 * j + 1.0 + E->control.particle_perturb * (0.5 - drand48())) * recip_density - 1.0);
	  
#endif	
	    eta3 = 0.0;  
  
	    E->tracer.tx[E->tracer.NUM_TRACERS] = E->tracer.tz[E->tracer.NUM_TRACERS] = 0.0;
	    v_shape_fn(E,el,lN,eta1,eta2,eta3,E->mesh.levmax);
  
	    for(e=1;e<=ends;e++) {
	      node = E->ien[el].node[e];
	      if((E->node[node] & PER_OFFSIDE)) {
	        if(e==3 || e==4)
		   node += E->mesh.noz * (E->mesh.nox-1);
	      }
	     
	      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];

              /*RAA: check coords*/     
              /*if(E->control.verbose) {
                if (el <= 24) { 
                  fprintf(stderr,"Here is eta1, eta2 : %d %d %f %f \n",  el,E->tracer.NUM_TRACERS,eta1,eta2);
                  fprintf(stderr,"Here is tx, tz : %d %d %f %f \n",  el,E->tracer.NUM_TRACERS,E->tracer.tx[E->tracer.NUM_TRACERS],E->tracer.tz[E->tracer.NUM_TRACERS]);  
                }  
                }  */
	      
	    } /*end of 'e' loop */

	    E->tracer.tracer_elt[E->mesh.levmax][E->tracer.NUM_TRACERS] = el;
  
	    for(level=E->mesh.levmax;level>E->mesh.levmin;level--) {  /* Roughly ... */
	      E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] = 
	        1 + E->tracer.tracer_elt[level][E->tracer.NUM_TRACERS] / ((E->mesh.nsd == 2) ? 4 : 8) ;  
	      if(E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] > E->mesh.NEL[level-1])
	        E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] = E->mesh.NEL[level-1];
	    }
  
	    E->tracer.property_group[E->tracer.NUM_TRACERS] = 0;	 	 

	    /*RAA: removed some if 3D stuff here, etc, will be done separately below*/
#if 0
	  if(tr_density_element[el] == 2) {
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 1.0 ;
	  }
	  if(tr_density_element[el] == 3) {
	    if(eta1 == 0 && eta2 == 0) {
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 64.0 / 81.0 ;
	    }
	    else if(eta1 * eta2 == 0) {
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 40.0 / 81.0 ;
	    }
	    else
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 25.0 / 81.0 ;
	  }
	  if(tr_density_element[el] == 4) {
	    if(fabs(eta1*eta2)>0.73  && fabs(eta1*eta2)<0.76 ) {
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 0.347854845137454 * 0.347854845137454 ;
	    }
	    else if(fabs(eta1*eta2)>0.29 && fabs(eta1*eta2)<0.3 ) {
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 0.347854845137454 * 0.652145154862546 ;
	    }
	    else 
	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
		E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 0.652145154862546*0.652145154862546 ;
	  }
#endif
	  E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
	      E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
	      4.0 * recip_density * recip_density;

	  /*RAA: 09/04/02, used to be 1.0 below, not 0.5, and is 1.0 in old U-Syd 2D version*/
	  E->tracer.dX11[E->tracer.NUM_TRACERS] = 0.5 * E->eco[el].size[1] * recip_density;
	  E->tracer.dX22[E->tracer.NUM_TRACERS] = 0.5 * E->eco[el].size[2] * recip_density;

	  E->tracer.tr_dim[E->tracer.NUM_TRACERS] = 0.5 * 
	    (E->tracer.dX11[E->tracer.NUM_TRACERS]+E->tracer.dX22[E->tracer.NUM_TRACERS]);
	}  /*end of j loop over z-direction*/
      }    /*end of 2D if loop */

      else { /*RAA: 3D */
        for(kkj=0;kkj<tr_density_element[el];kkj++)
      	  for(i=0;i<tr_density_element[el];i++) {
  	    for(j=0;j<tr_density_element[el];j++) {
  	      E->tracer.NUM_TRACERS++;
  	      tracer_allocate_memory(E);

	      eta1 = ((2 * i + 1.0 + E->control.particle_perturb * (0.5 - drand48())) * recip_density - 1.0);
	      eta2 = ((2 * j + 1.0 + E->control.particle_perturb * (0.5 - drand48())) * recip_density - 1.0);
              eta3 = ((2 * kkj + 1.0 + E->control.particle_perturb * (0.5 - drand48())) * recip_density - 1.0);

              E->tracer.tx[E->tracer.NUM_TRACERS] = E->tracer.tz[E->tracer.NUM_TRACERS] = 0.0;
              E->tracer.ty[E->tracer.NUM_TRACERS] = 0.0;
  
	      v_shape_fn(E,el,lN,eta1,eta2,eta3,E->mesh.levmax);

               /*RAA: 31/10/01, added this upper part for both per_x and per_y*/
               /*RAA: This stuff is needed due to ambiguities with edges and 
                      isn't optimum in its writing but should work ok.*/
              if(E->mesh.periodic_y && E->mesh.periodic_x) { 
        	for(p=1;p<=E->mesh.elz;p++)
        	   for(q=1;q<=E->mesh.elx;q++)
		      for(r=1;r<=E->mesh.ely;r++) {
		        element = (r-1)*E->mesh.elz*E->mesh.elx + (q-1)*E->mesh.elz + p;

                        if (element==el && r==1 && q==1) { 
	                   for(e=1;e<=ends;e++) {
	                      node = E->ien[el].node[e];
	                      if((E->node[node] & PER_OFFSIDE)) 
	                        if(e==7 || e==8) 
		                   node += E->mesh.noz * (E->mesh.nox-1); /*per_y ok, too*/
	                      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	                      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	                      E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];
                           } 
                           break;
                        } 
                        else if (element==el && r==1 && q!=1 && q<E->mesh.elx) { 
	                   for(e=1;e<=ends;e++) {
	                      node = E->ien[el].node[e];
	                      if((E->node[node] & PER_OFFSIDE)) 
	                        if(e==5 || e==6 || e==7 || e==8) 
		                   node += E->mesh.noz * E->mesh.nox * (E->mesh.noy-1);
	                      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	                      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	                      E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];
                           } 
                           break;
                        } 
                        else if (element==el && r==1 && q==E->mesh.elx) { 
	                   for(e=1;e<=ends;e++) {
	                      node = E->ien[el].node[e];
	                      if(e==3 || e==4) /* e=3,4 were not counted as PER_OFFSIDE*/  
		                 node += E->mesh.noz * (E->mesh.nox-1);
	                      else if((E->node[node] & PER_OFFSIDE)) {
	                         if(e==7 || e==8) 
		                   node += E->mesh.noz * (E->mesh.nox-1);
                              }
	                      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	                      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	                      E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];
                           } 
                           break; 
                        } 
                        else if (element==el && r==E->mesh.ely) { 
	                   for(e=1;e<=ends;e++) {
	                      node = E->ien[el].node[e]; /*no PER_OFFSIDE statement here*/
	                      if(e==5 || e==6 || e==7 || e==8) 
		                 node += E->mesh.noz * E->mesh.nox * (E->mesh.noy-1);
	                      if(q==E->mesh.elx && (e==3 || e==4))  
		                 node += E->mesh.noz * (E->mesh.nox-1);
	                      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	                      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	                      E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];
                           } 
                           break;
                        } 
                        else if (element==el && r>1 && r<E->mesh.ely) { 
	                   for(e=1;e<=ends;e++) {
	                      node = E->ien[el].node[e];
	                      if((E->node[node] & PER_OFFSIDE)) 
	                        if(e==3 || e==4 || e==7 || e==8) 
		                   node += E->mesh.noz * (E->mesh.nox-1);
	                      E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	                      E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	                      E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];
                           } 
                           break;
                        } 
                      }  /*end of 'r' loop */

              }  /*end of if per_x & per_y */
              else {
	         for(e=1;e<=ends;e++) {
	            node = E->ien[el].node[e];
	            if((E->node[node] & PER_OFFSIDE)) {
                      if(E->mesh.periodic_x && !E->mesh.periodic_y) {  /*RAA: 5/10/01, added this distinction*/
	                if(e==3 || e==4 || e==7 || e==8)   /*RAA: 30/5/01, added e==7 and e==8*/
		           node += E->mesh.noz * (E->mesh.nox-1);
                    }
                    else if(E->mesh.periodic_y && !E->mesh.periodic_x) {  /*RAA: 5/10/01, added lines for per_y*/
	               if(e==5 || e==6 || e==7 || e==8)  
		          node += E->mesh.noz * E->mesh.nox * (E->mesh.noy-1);
                    } 
	         }
	         E->tracer.tx[E->tracer.NUM_TRACERS] += lN[e] * E->x[1][node];
	         E->tracer.tz[E->tracer.NUM_TRACERS] += lN[e] * E->x[2][node];
	         E->tracer.ty[E->tracer.NUM_TRACERS] += lN[e] * E->x[3][node];

	        }  /* end of 'e' loop */
              }  /* end of 'else' part */

/*
                 if(E->control.verbose) {
                   if (el<26) {
                     fprintf(stderr,"Here is el: %d  num: %d; eta1 eta2 eta3: %f %f %f \n",el,E->tracer.NUM_TRACERS,eta1,eta2,eta3);
                   }
                 } 
*/
/*                if(E->control.verbose) {
                    fprintf(stderr,"Here is el: %d  num: %d; tx, tz ty: %f %f %f \n",el,E->tracer.NUM_TRACERS,E->tracer.tx[E->tracer.NUM_TRACERS],E->tracer.tz[E->tracer.NUM_TRACERS],E->tracer.ty[E->tracer.NUM_TRACERS]);  
                } 
*/

	      E->tracer.tracer_elt[E->mesh.levmax][E->tracer.NUM_TRACERS] = el;

            /*RAA: what to do here (below) for 3D?*/
	      for(level=E->mesh.levmax;level>E->mesh.levmin;level--) {  /* Roughly ... */
	        E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] = 
	          1 + E->tracer.tracer_elt[level][E->tracer.NUM_TRACERS] / 8.0; /* RAA: skip the nested if loop here */
	        if(E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] > E->mesh.NEL[level-1])
	          E->tracer.tracer_elt[level-1][E->tracer.NUM_TRACERS] = E->mesh.NEL[level-1];
	      }

	      E->tracer.property_group[E->tracer.NUM_TRACERS] = 0;	 	 

	      E->tracer.tracer_std_weighting[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
	          E->tracer.tracer_weight_fn[E->mesh.levmax][E->tracer.NUM_TRACERS] = 
	          8.0 * recip_density * recip_density * recip_density;

	      /*RAA: not done yet- what about tr_dim in 3D?*/
	      /*RAA: 09/04/02, used to be 1.0 below, not 0.5, and is 1.0 in old U-Syd 3D version*/
	      E->tracer.dX11[E->tracer.NUM_TRACERS] = 0.5 * E->eco[el].size[1] * recip_density;
	      E->tracer.dX22[E->tracer.NUM_TRACERS] = 0.5 * E->eco[el].size[2] * recip_density;
	      E->tracer.dX33[E->tracer.NUM_TRACERS] = 0.5 * E->eco[el].size[3] * recip_density;
	      /*  E->tracer.tr_dim[E->tracer.NUM_TRACERS] = 0.5 * 
	          (E->tracer.dX11[E->tracer.NUM_TRACERS]+E->tracer.dX22[E->tracer.NUM_TRACERS]); */
	      E->tracer.tr_dim[E->tracer.NUM_TRACERS] = 0.333333 * 
	          (E->tracer.dX11[E->tracer.NUM_TRACERS]+E->tracer.dX22[E->tracer.NUM_TRACERS]+E->tracer.dX33[E->tracer.NUM_TRACERS]);
            

            } /*end of j loop over z-direction*/
          }   /*end of i loop over x-direction*/

      }   /*end of 3D 'else' statement*/

    }      /*end of el loop */
      
    fprintf(E->fp,"Allocated %d tracers in %d elements\n",
	    E->tracer.NUM_TRACERS,E->mesh.nel);
    
    /* For the existing distribution of tracer, add color information
       which dictates the material properties */

    input_int("Material_rect",&num_rect,"0");
    input_int_vector("Material_rect_property",num_rect,rect_tracer_colour); 
    input_std_precision_vector("Material_rect_x1",num_rect,rect_tx1);
    input_std_precision_vector("Material_rect_x2",num_rect,rect_tx2);
    input_std_precision_vector("Material_rect_z1",num_rect,rect_tz1);
    input_std_precision_vector("Material_rect_z2",num_rect,rect_tz2);
    input_std_precision_vector("Material_rect_y1",num_rect,rect_ty1);
    input_std_precision_vector("Material_rect_y2",num_rect,rect_ty2);
    
    fprintf(E->fp1,"Making tracer rectangles, dims = %g \n",dims);
    fprintf(stderr,"Making tracer rectangles, dims = %g \n",dims);

    for(i=0;i<num_rect;i++) {   
       if(rect_tx1[i] < E->sx[1][1]) 
	rect_tx1[i] = E->sx[1][1];
      if(rect_tx2[i] < E->sx[1][1]) 
	rect_tx2[i] = E->sx[1][1];
      if(rect_tx1[i] > E->sx[1][E->mesh.nno]) 
	rect_tx1[i] = E->sx[1][E->mesh.nno];
      if(rect_tx2[i] > E->sx[1][E->mesh.nno]) 
	rect_tx2[i] = E->sx[1][E->mesh.nno];
      if(rect_tz1[i] < E->sx[2][1]) 
	rect_tz1[i] = E->sx[2][1];
      if(rect_tz2[i] < E->sx[2][1]) 
	rect_tz2[i] = E->sx[2][1];
      if(rect_tz1[i] > E->sx[2][E->mesh.nno]) 
	rect_tz1[i] = E->sx[2][E->mesh.nno];
      if(rect_tz2[i] > E->sx[2][E->mesh.nno]) 
	rect_tz2[i] = E->sx[2][E->mesh.nno];

      
      if(3==dims) {
	if(rect_ty1[i] < E->sx[3][1]) 
	  rect_ty1[i] = E->sx[3][1];
	if(rect_ty2[i] < E->sx[3][1]) 
	  rect_ty2[i] = E->sx[3][1];
	if(rect_ty1[i] > E->sx[3][E->mesh.nno]) 
	  rect_ty1[i] = E->sx[3][E->mesh.nno];
	if(rect_ty2[i] > E->sx[3][E->mesh.nno]) 
	  rect_ty2[i] = E->sx[3][E->mesh.nno];
      }
    
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) { /* RAA --> added this distinction  */
	if(2==dims) {  
	  if((E->tracer.tx[m] > rect_tx1[i]) && (E->tracer.tx[m] <= rect_tx2[i]) &&
	     (E->tracer.tz[m] > rect_tz1[i]) && (E->tracer.tz[m] <= rect_tz2[i]) ) {
	      E->tracer.property_group[m] = rect_tracer_colour[i];
	  }
	}
	else if(3==dims) {	/* RAA --> added this bit  */
	  if((E->tracer.tx[m] > rect_tx1[i]) && (E->tracer.tx[m] <= rect_tx2[i]) &&
	     (E->tracer.tz[m] > rect_tz1[i]) && (E->tracer.tz[m] <= rect_tz2[i]) &&
	     (E->tracer.ty[m] > rect_ty1[i]) && (E->tracer.ty[m] <= rect_ty2[i]) ) {
	      E->tracer.property_group[m] = rect_tracer_colour[i];
	  }
        }
      }   /*end of for m */
    }     /*end of for i */
      
     input_int("Material_circ",&num_rect,"0");
     input_int_vector("Material_circ_property",num_rect,rect_tracer_colour); 
     input_std_precision_vector("Material_circ_x1",num_rect,rect_tx1);
     input_std_precision_vector("Material_circ_rad",num_rect,rect_tx2);
     input_std_precision_vector("Material_circ_z1",num_rect,rect_tz1);  
     input_std_precision_vector("Material_circ_y1",num_rect,rect_ty1);
   
     fprintf(E->fp1,"dims node_dims %g %g \n",dims,E->mesh.nsd);
     fprintf(stderr,"dims_node_dims %g %g \n",dims,E->mesh.nsd);
     fprintf(stderr,"Tracer_read, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);


     if (2==dims) {
     for(i=0;i<num_rect;i++) {   
       for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	 if(((E->tracer.tx[m]-rect_tx1[i])*(E->tracer.tx[m]-rect_tx1[i]) +
	    (E->tracer.tz[m]-rect_tz1[i])*(E->tracer.tz[m]-rect_tz1[i]))
	    < rect_tx2[i] * rect_tx2[i] ) {
	   E->tracer.property_group[m] = rect_tracer_colour[i];
	   fprintf(E->fp1,"Oh hello, mesh.nsd=2 \n");
	 }
       }
     }
     }
     else if (3==dims) {
	for(i=0;i<num_rect;i++) {   
       	  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	   if(((E->tracer.tx[m]-rect_tx1[i])*(E->tracer.tx[m]-rect_tx1[i]) +
	      (E->tracer.tz[m]-rect_tz1[i])*(E->tracer.tz[m]-rect_tz1[i]) +
		(E->tracer.ty[m]-rect_ty1[i])*(E->tracer.ty[m]-rect_ty1[i])) 	
	       < (rect_tx2[i] * rect_tx2[i] )) {
	     E->tracer.property_group[m] = rect_tracer_colour[i];
	     fprintf(E->fp1,"Goody!, mesh.nsd=3 \n");

	   }
          }
        }
     }
     
     input_int("Material_trgl",&num_rect,"0");
     input_int_vector("Material_trgl_property",num_rect,rect_tracer_colour); 
     input_std_precision_vector("Material_trgl_x1",num_rect,rect_tx1);
     input_std_precision_vector("Material_trgl_x2",num_rect,rect_tx2);
     input_std_precision_vector("Material_trgl_x3",num_rect,rect_tx3);
     input_std_precision_vector("Material_trgl_x4",num_rect,rect_tx4);
     input_std_precision_vector("Material_trgl_z1",num_rect,rect_tz1);
     input_std_precision_vector("Material_trgl_z2",num_rect,rect_tz2);
     input_std_precision_vector("Material_trgl_z3",num_rect,rect_tz3);
     input_std_precision_vector("Material_trgl_z4",num_rect,rect_tz4);
     input_std_precision_vector("Material_trgl_y1",num_rect,rect_ty1);
     input_std_precision_vector("Material_trgl_y2",num_rect,rect_ty2);
     input_std_precision_vector("Material_trgl_y3",num_rect,rect_ty3);
     input_std_precision_vector("Material_trgl_y4",num_rect,rect_ty4);
   
   
     for(i=0;i<num_rect;i++) {   
       /* get the edge vectors and area for each triangle before
	  scanning the tracers for a hit */

       /* 2D */
       a1 =  rect_tx2[i] - rect_tx1[i];
       a2 =  rect_tz2[i] - rect_tz1[i];    
       b1 =  rect_tx3[i] - rect_tx1[i];
       b2 =  rect_tz3[i] - rect_tz1[i];   
       c1 =  rect_tx3[i] - rect_tx2[i];
       c2 =  rect_tz3[i] - rect_tz2[i];
   
       /* 3D coming soon */

       ref_area = magcross2d(a1,a2,b1,b2);  /* 2 x area of triangle */
    
       for(m=1;m<=E->tracer.NUM_TRACERS;m++) {

	 /* 2D */
	 p1 = E->tracer.tx[m] - rect_tx1[i];
	 p2 = E->tracer.tz[m] - rect_tz1[i];
	 q1 = E->tracer.tx[m] - rect_tx2[i];
	 q2 = E->tracer.tz[m] - rect_tz2[i];
	 r1 = E->tracer.tx[m] - rect_tx3[i];
	 r2 = E->tracer.tz[m] - rect_tz3[i];

	 area_sum = 0.0;

	 area_sum += magcross2d(a1,a2,p1,p2); 
	 area_sum += magcross2d(c1,c2,q1,q2); 
	 area_sum += magcross2d(b1,b2,r1,r2); 

	 if(area_sum <= ref_area)
	   E->tracer.property_group[m] = rect_tracer_colour[i];
       }
     }
  }

#if defined (CHEM_TRANS_MODULE_INSTALLED) 
  if(E->control.CHEM_TRANS)
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) 
      if( E->tracer.visc[E->tracer.property_group[i]].mobile_phase_ratio != 1.0) 
	E->tracer.volfraction[i] = 0.5;
      else
	E->tracer.volfraction[i] = 0.0;
      
#endif


    /* Here you can put ad hoc modifications for specific problems !!! */
    input_int("activate_secret_test_problem",&testcases,"0");
    input_std_precision_vector("test_geometry_parameters",99,geom);  /* What if fewer entries ?? */
    input_int_vector("test_geometry_materials",99,imat);  /* What if fewer entries ?? */
    input_int_vector("test_geometry_index",99,idx);

    if(testcases == 1) {
      /* Rayleigh-Taylor Benchmark problem */
      E->control.secret_case_fl[1] = geom[0] ;
      E->control.secret_case_fl[2] = geom[1] ;
      E->control.secret_case_fl[3] = geom[2] ;
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.tz[i] < geom[0] + geom[1] * cos(M_PI * E->tracer.tx[i] / geom[2])) {
	  E->tracer.property_group[i] = 0;
	}
	else
	  E->tracer.property_group[i] = 1;
	/* E->tracer.weight[i] = 50-50*tanh((E->tracer.tz[i] - (0.8 + 0.02 * cos(M_PI * E->tracer.tx[i] / 0.9142)))/0.001); */
      }
    }
    else if(testcases == 7) {
      /* RAA: 12/12/02, similar to above, a sine interface in  a 3 layer problem*/
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.tz[i] < geom[0] + geom[1] * sin(M_PI * E->tracer.tx[i] / geom[2]))
	  E->tracer.property_group[i] = 0;
	else if(E->tracer.tz[i] > geom[3]) /*raa - boundary between mats 1 & 2 */
	  E->tracer.property_group[i] = 2;
	else 
	  E->tracer.property_group[i] = 1;
      }
    }
    else if (testcases == 2) {
      /* Gravity current benchmark */

      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.tx[i] < (geom[3] * E->tracer.tz[i]*E->tracer.tz[i]*E->tracer.tz[i] +
			      geom[2] * E->tracer.tz[i]*E->tracer.tz[i]+
			      geom[1] * E->tracer.tz[i] +
			      geom[0])) 
	  E->tracer.property_group[i] = 1;
	else
	  E->tracer.property_group[i] = 0;
      }
    }
    else if (testcases == 3) {
      /* Suspension flow - lots of different circles or ellipses but centred
	 on a grid  */
      /* Dimensions are 1 x 1 */


      for(i=1;i<=14;i++) {
	for(j=1;j<=7;j++) {
	  a1 = 0.0714286 + 0.142857 * (i-1);
	  a2 = 0.0714286 + 0.142857 * (j-1);

	  e1 = drand48();

	  if(1 || e1 < 0.45) {
	    /* d1 =  M_PI * drand48(); 	
	       d2 = 1.0 + 3.0*drand48();
	       d2 = 1.0 / (1.0e-16+d2*d2);  
	       c1=sin(d1);
	       c2=cos(d1);
	    */

	    d1 = 0.0;
	    d2 = 1.0;
	    c1=0.0;
	    c2=1.0;

      
	    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	      b1 = (E->tracer.tx[m] - a1) * c1 - (E->tracer.tz[m] - a2) * c2;
	      b2 = (E->tracer.tx[m] - a1) * c2 + (E->tracer.tz[m] - a2) * c1;
	  
	      if( (b1 * b1 * d2  + b2 * b2)  < 0.02*0.02)
		E->tracer.property_group[m] = 1;
	    }
	  }
	  else if (e1 < 0.9) {
	    /* d1 = M_PI * drand48();
	       d2 = 1.0 + 3.0*drand48();
	       d2 = 1.0 / (1.0e-16+d2*d2);
	
	       c1=sin(d1);
	       c2=cos(d1); */

	    d1 = 0.0;
	    d2 = 1.0;
	    c1=0.0;
	    c2=1.0;
	
	    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	      b1 = (E->tracer.tx[m] - a1) * c1 - (E->tracer.tz[m] - a2) * c2;
	      b2 = (E->tracer.tx[m] - a1) * c2 + (E->tracer.tz[m] - a2) * c1;
	  
	      if( (b1 * b1 * d2  + b2 * b2)  < 0.02*0.02)
		E->tracer.property_group[m] = 2;
	    }
       
	  }
	}
      }
    }
    else if (testcases == 4) {

      /* Layers + harmonic perturbation */
    
      idx[0]=min(idx[0],19);
   
      /* Set background */
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	if(E->tracer.tz[m] < geom[2])
	  E->tracer.property_group[m] = 0; 
	else
	  E->tracer.property_group[m] = 1; 	
      }
    
      /* Now loop of the number of layers */

      for(i=0;i<idx[0];i++) {
	j=3+i; /* offset into array of geom parameters */
	for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
       
	  /* Layers specified by initial offset & incremental thickness  */
       
	  d1 = geom[2] + sin(E->tracer.tx[m]*geom[0]*M_PI) * geom[1];
	  d2 = d1 + geom[j];
	
	  if(E->tracer.tz[m] >= d1 && E->tracer.tz[m] <= d2)
	    E->tracer.property_group[m] = imat[i];
	}
	geom[2] += geom[j];  /* After the loop, let the thickness be replaced by the depth */
      }
    }
    else if (testcases == 5) {
      /*initial perturbation for buckling */

    if(E->control.ELASTICITY) {
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.property_group[i] == 1) {
	  E->tracer.S11[i] = -geom[3] ;
	  E->tracer.S22[i] = -E->tracer.S11[i] ;
	}
	if(E->tracer.property_group[i] == 2) {
	  E->tracer.S11[i] = -geom[4] ;
	  E->tracer.S22[i] = -E->tracer.S11[i] ;
	}
      }
      E->control.secret_case_fl[5] = 2.*geom[3] ;
    }
    E->control.secret_case_fl[1] = E->control.secret_case_fl[4] = geom[0] ;
    E->control.secret_case_fl[2] = geom[1] ;
    E->control.secret_case_fl[3] = geom[2] ;
    fprintf(stderr,"For buckling problem (particles) w0 = %g k = %g and h = %g stress in beam %g and outer layer = %g \n",
	    geom[0],geom[1],geom[2],geom[3],geom[4]) ;
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      if(E->tracer.property_group[m]>=1) {
	if((E->tracer.tz[m] > 
	    (0.5*E->sx[2][E->mesh.nno]-E->control.secret_case_fl[3]/2.+E->control.secret_case_fl[1]*
	     cos(E->control.secret_case_fl[2]*M_PI/3.0*E->tracer.tx[m]))) &&
	   (E->tracer.tz[m] < 
	    (0.5*E->sx[2][E->mesh.nno]+E->control.secret_case_fl[3]/2.+E->control.secret_case_fl[1]*
	     cos(E->control.secret_case_fl[2]*M_PI/3.0*E->tracer.tx[m])))) {
	  E->tracer.property_group[m] = 1;
/*	  fprintf(stderr,"tracer %d x=%g z=%g belong to the beam \n",m,E->tracer.tx[m],E->tracer.tz[m]) ;*/
	}
	else
	  E->tracer.property_group[m] = 2 ;
      }
    }
   }
    else if (testcases == 6) {
      /* Layering for Paul problems */
      Dz = 0.5*(geom[1]-geom[0])/geom[2] ;
      for(i=1;i<=geom[2];i++) {
	for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
	  if((E->tracer.tz[m] >= geom[0]+Dz*(2*i-2)) && (E->tracer.tz[m] <= geom[0]+Dz*(2*i-1)))
	    E->tracer.property_group[m] = geom[3] ;
	  else if((E->tracer.tz[m] > geom[0]+Dz*(2*i-1)) && (E->tracer.tz[m] <= geom[0]+Dz*2*i)) {
	    E->tracer.property_group[m] = geom[4] ;
	  }
	}
      }
    }

#if 0
    /* Initial stress state */

    for(i=1;i<=E->tracer.NUM_TRACERS;i++)
      if(E->tracer.tz[i] >= 0.2 && E->tracer.tx[i] >= 0.1 && E->tracer.tx[i] <= 0.4) {
	E->tracer.S22[i] = 0.5*(0.2-E->tracer.tz[i])*100  ;
	E->tracer.S11[i] = -0.5*(0.2-E->tracer.tz[i])*100  ;
      }
#endif

    fprintf(stderr,"Allocated/initialized tracer values \n");

/*RAA: 4/6/01, check the tracer locations */
/*  for(j=1;j<=E->tracer.NUM_TRACERS;j++) 
      fprintf(stderr,"Here is #, elnum, tx, tz ty: %d %d %f %f %f \n",j,E->tracer.tracer_elt[E->mesh.levmax][j],E->tracer.tx[j],E->tracer.tz[j],E->tracer.ty[j]);  
*/
    
    /* Figure out which tracers are in which element */
    get_tracer_elts_and_weights(E,1,E->tracer.NUM_TRACERS);
    /* Initialize tracer strain values */

    /*
    fprintf(stderr,"Orthotropy = %d:  %g,%g\n",
		E->control.ORTHOTROPY,E->control.director_perturbkx,E->control.director_perturb);
    */


    if(E->control.ORTHOTROPY) {
    	for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
    		a1 = sin(E->control.director_perturbkx * M_PI * E->tracer.tx[i]) * 
				E->control.director_perturb * M_PI /* (0.5 - drand48()) */;   
	  E->tracer.n1[i] = cos(a1);
	  E->tracer.n2[i] = sin(a1);
	  if(3==dims)
	    E->tracer.n3[i] = 0.0;	
	    
    	}
    }
      /* for(el=1;el<=E->mesh.nel;el++) {
	for(i=0;i<E->tracer.tr_in_element_number[E->mesh.levmax][el];i++) {
	  tr = E->tracer.tr_in_element[E->mesh.levmax][i+E->tracer.tr_in_element_offset[E->mesh.levmax][el]];
	  a1 = cos(E->control.director_perturbkx * M_PI * E->tracer.tx[tr]) * 
			E->control.director_perturb * M_PI ;   
	  E->tracer.n1[tr] = cos(a1);
	  E->tracer.n2[tr] = sin(a1);
	  if(3==dims)
	    E->tracer.n3[tr] = 0.0;
	}
      } */
   
    

    if(read_previous_field(E,E->Pstrain,"plastic_strain","Pstn"))
      nodes_to_tracers(E,E->Pstrain,E->tracer.edotp_integrated,E->mesh.levmax);
    else
      for(i=1;i<=E->tracer.NUM_TRACERS;i++)
	E->tracer.edotp_integrated[i]=0.0;


/*RAA: 11/07/02, added this part for a block of weak material at start*/
    input_int("p_strain_block_rect",&num_rect,"0");
    input_int("p_strain_by_location_or_material",&loc_or_mat,"100");
    input_int("p_strain_material_group",p_strain_mat_group,"0");
    input_std_precision_vector("p_strain_block_rect_x1",num_rect,rect_tx1);
    input_std_precision_vector("p_strain_block_rect_x2",num_rect,rect_tx2);
    input_std_precision_vector("p_strain_block_rect_z1",num_rect,rect_tz1);
    input_std_precision_vector("p_strain_block_rect_z2",num_rect,rect_tz2);
    input_std_precision_vector("p_strain_block_rect_y1",num_rect,rect_ty1);
    input_std_precision_vector("p_strain_block_rect_y2",num_rect,rect_ty2);
    input_std_precision_vector("p_strain_block_rect_mag",num_rect,p_strain_mag); 

    for(i=0;i<num_rect;i++) {   
      if(loc_or_mat==0) {
        for(m=1;m<=E->tracer.NUM_TRACERS;m++) { 
	  if(2==dims) {  
	    if((E->tracer.tx[m] >= rect_tx1[i]) && (E->tracer.tx[m] <= rect_tx2[i]) &&
	       (E->tracer.tz[m] >= rect_tz1[i]) && (E->tracer.tz[m] <= rect_tz2[i]) ) 
               E->tracer.edotp_integrated[m]=p_strain_mag[i];
	  }
	  else if(3==dims) {
	    if((E->tracer.tx[m] >= rect_tx1[i]) && (E->tracer.tx[m] <= rect_tx2[i]) &&
	       (E->tracer.tz[m] >= rect_tz1[i]) && (E->tracer.tz[m] <= rect_tz2[i]) &&
	       (E->tracer.ty[m] >= rect_ty1[i]) && (E->tracer.ty[m] <= rect_ty2[i]) ) 
               E->tracer.edotp_integrated[m]=p_strain_mag[i];
          }
        }
      } /* end of if loc_or_mat 0  */
      else if(loc_or_mat==1) {
        for(m=1;m<=E->tracer.NUM_TRACERS;m++) { 
	  if(p_strain_mat_group[i] == E->tracer.property_group[m]) 
             E->tracer.edotp_integrated[m]=p_strain_mag[i];
	}
      } /* end of if loc_or_mat 0  */
    } /* end of for 'i' */

/*RAA - end of initial plastic strain block stuff */

    
    input_int("Strain_trgl",&num_rect,"0");
    input_std_precision_vector("Strain_trgl_value",num_rect,rect_XX); 
    input_std_precision_vector("Strain_trgl_x1",num_rect,rect_tx1);
    input_std_precision_vector("Strain_trgl_x2",num_rect,rect_tx2);
    input_std_precision_vector("Strain_trgl_x3",num_rect,rect_tx3);
    input_std_precision_vector("Strain_trgl_x4",num_rect,rect_tx4);
    input_std_precision_vector("Strain_trgl_z1",num_rect,rect_tz1);
    input_std_precision_vector("Strain_trgl_z2",num_rect,rect_tz2);
    input_std_precision_vector("Strain_trgl_z3",num_rect,rect_tz3);
    input_std_precision_vector("Strain_trgl_z4",num_rect,rect_tz4);
    input_std_precision_vector("Strain_trgl_y1",num_rect,rect_ty1);
    input_std_precision_vector("Strain_trgl_y2",num_rect,rect_ty2);
    input_std_precision_vector("Strain_trgl_y3",num_rect,rect_ty3);
    input_std_precision_vector("Strain_trgl_y4",num_rect,rect_ty4);
   
   
    for(i=0;i<num_rect;i++) {   
      /* get the edge vectors and area for each triangle before
	 scanning the tracers for a hit */

      /* 2D */
      a1 =  rect_tx2[i] - rect_tx1[i];
      a2 =  rect_tz2[i] - rect_tz1[i];    
      b1 =  rect_tx3[i] - rect_tx1[i];
      b2 =  rect_tz3[i] - rect_tz1[i];   
      c1 =  rect_tx3[i] - rect_tx2[i];
      c2 =  rect_tz3[i] - rect_tz2[i];
   
      /* 3D coming soon */

      ref_area = magcross2d(a1,a2,b1,b2);  /* 2 x area of triangle */
 
      for(m=1;m<=E->tracer.NUM_TRACERS;m++) {

	/* 2D */
	p1 = E->tracer.tx[m] - rect_tx1[i];
	p2 = E->tracer.tz[m] - rect_tz1[i];
	q1 = E->tracer.tx[m] - rect_tx2[i];
	q2 = E->tracer.tz[m] - rect_tz2[i];
	r1 = E->tracer.tx[m] - rect_tx3[i];
	r2 = E->tracer.tz[m] - rect_tz3[i];

	area_sum = 0.0;

	area_sum += magcross2d(a1,a2,p1,p2); 
	area_sum += magcross2d(c1,c2,q1,q2); 
	area_sum += magcross2d(b1,b2,r1,r2); 

	if(area_sum <= ref_area) {
	  E->tracer.yielded[m] = 1.0;
	  E->tracer.edotp_integrated[m] = rect_XX[i];
	}
      }
    }  

    fprintf(stderr,"Computed weights etc \n");

    /* Get temperatures at tracers (BUT what if it was read in 
       already using the tracer file, eh, what're you going to do
       then ?  At the moment this is OK since nodal temperatures
       have the priority but this may be different when advection
       schemes are modified */

    nodes_to_tracers(E,E->T,E->tracer.T,E->mesh.levmax);

    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      E->tracer.Current_phase[m] = E->tracer.Equilibrium_phase[m] = get_eq_phase(E,m); 
    }
    free((void *) tr_density_element);
    return;
  }
