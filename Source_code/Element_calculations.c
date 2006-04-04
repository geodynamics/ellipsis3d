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



#include "config.h"

#include <math.h>

#if HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"


struct SORT {
  standard_precision value;
  int node; 
};

void build_diagonal_of_K(
			 struct All_variables *E,
			 int level
			 )
{
  int i,a,e,k,loc0,found1;
  int compare_function();

  struct SORT *nodesort;
      
  const int neq=E->mesh.NEQ[level];
  const int nno=E->mesh.NNO[level];
  const int nel=E->mesh.NEL[level];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int max_eqn = dofs * max_node_interaction[dims];
  const int max_node = max_node_interaction[dims];

  higher_precision b1,b2,b3,b4,b5,b6;

  standard_precision time,CPU_time();

  nodesort = (struct SORT *) Malloc0((E->mesh.nno+100) * sizeof(struct SORT));
 
/*  E->monitor.max_BI[level] = 0.0;
  E->monitor.min_BI[level] = 1.0e32;*/

  if (E->control.B_is_good[level] == 1) return;
 
  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  )) {
      E->BI1[level][i] = 0.0;
      E->BI2[level][i] = 0.0;
      if(dofs==3)
	E->BI3[level][i] = 0.0;
      else  if(dofs==6) {
	E->BI3[level][i] = 0.0;
	E->BI4[level][i] = 0.0;
	E->BI5[level][i] = 0.0;
	E->BI6[level][i] = 0.0;
      }
      continue;
    }

    found1=0; 
    loc0=i*max_node;
    for(k=0;k<max_node;k++) {
      if(E->Node_map_3[level][loc0+k] == i) {
	b1= E->Node_k11[level][loc0+k];
	b2= E->Node_k22[level][loc0+k];
	if(dofs==3)
	  b3= E->Node_k33[level][loc0+k];
	else if(dofs==6) {
	  b3= E->Node_k33[level][loc0+k];
	  b4= E->Node_k44[level][loc0+k];
	  b5= E->Node_k55[level][loc0+k];
	  b6= E->Node_k66[level][loc0+k];
	}
      }
    }

    if(b1 != 0.0)   /* Boundary conditions */
      E->BI1[level][i] = 1.0/b1;
    else
      E->BI1[level][i] = 0.0;
    if(b2 != 0.0)
      E->BI2[level][i] = 1.0/b2;
    else
      E->BI2[level][i] = 0.0;
    if(dofs==3) {
      if(b3 != 0.0)
	E->BI3[level][i] = 1.0/b3;
      else
	E->BI3[level][i] = 0.0;
    }
    else if(dofs==6) {
      if(b3 != 0.0)   /* Boundary conditions */
	E->BI3[level][i] = 1.0/b3;
      else
	E->BI3[level][i] = 0.0;
      if(b4 != 0.0)
	E->BI4[level][i] = 1.0/b4;
      else
	E->BI4[level][i] = 0.0;
      if(b5 != 0.0)   /* Boundary conditions */
	E->BI5[level][i] = 1.0/b5;
      else
	E->BI5[level][i] = 0.0;
      if(b6 != 0.0)
	E->BI6[level][i] = 1.0/b6;
      else
	E->BI6[level][i] = 0.0;
    }

    E->BI[level][E->ID[level][i].doff[1]] = E->BI1[level][i];
    E->BI[level][E->ID[level][i].doff[2]] = E->BI2[level][i];
    if(dofs==3)
      E->BI[level][E->ID[level][i].doff[3]] = E->BI3[level][i];
    else if(dofs==6) {
      E->BI[level][E->ID[level][i].doff[3]] = E->BI3[level][i];
      E->BI[level][E->ID[level][i].doff[4]] = E->BI4[level][i];
      E->BI[level][E->ID[level][i].doff[5]] = E->BI5[level][i];
      E->BI[level][E->ID[level][i].doff[6]] = E->BI6[level][i];
    }
  }
  E->control.B_is_good[level] = 1;

  /* Create a sort for node ordering based on the values of BI */    

  for(i=1;i<=E->mesh.NNO[level];i++) {
    nodesort[i].value = E->BI1[level][i];
    nodesort[i].node = i;
  }

  qsort(nodesort+1,E->mesh.NNO[level],sizeof(struct SORT),compare_function);
 
  for(i=1;i<=E->mesh.NNO[level];i++) {
    E->index[level][i] = nodesort[i].node;
  }

  free((void *) nodesort);
  return;
}

int compare_function
( struct SORT *s1,
  struct SORT *s2
  )
{
  if(s1->value < s2->value)
    return(-1);
  else if (s1->value > s2->value)
    return(1);
  else
    return(0);
}


void build_diagonal_of_Ahat_level(
				  struct All_variables *E,
				  int level
				  )
{
  higher_precision assemble_dAhatq_entry();
  
  int i,j,k,el,m;
  standard_precision time,time0,CPU_time(),visc,comp,v1;
  
  const int npno = E->mesh.npno;
  const int neq = E->mesh.neq;
  const int dims = E->mesh.nsd;
  
  time=CPU_time();
  
  if(E->control.verbose)
    fprintf(stderr,"Building pressure preconditioner %d\n",level);
  
  for(i=1;i<=E->mesh.NPNO[level];i++)
    E->BPI[level][i]=1.0;

  for(el=1;el<=E->mesh.NEL[level];el++) {
    visc = 0.0;
    comp = 0.0;
    
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
       
      visc += 1.0/(E->tracer.Visc[m]);

       if(3==dims) {
	if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0 && 
	   E->tracer.BulkVisc[m] > 0.6666666 * E->tracer.Visc[m]) 
	  comp += 1.0/(E->tracer.BulkVisc[m] - 0.6666666 * E->tracer.Visc[m]);  
      } 
      else {
	if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0   && 
	   E->tracer.BulkVisc[m] > E->tracer.Visc[m] ) 
	  comp += 1.0/(E->tracer.BulkVisc[m]  - E->tracer.Visc[m]);   
      }
    }

    if(E->tracer.tr_in_element_number[level][el] != 0) {
      visc = E->tracer.tr_in_element_number[level][el]/(visc);
      
      v1 = (higher_precision) assemble_dAhatq_entry(E,el,level);
      
      
      if(comp != 0.0) {
	comp = E->tracer.tr_in_element_number[level][el]/comp;
	/* E->BPI[level][el] = visc * comp / (comp+visc) ; /* */

	    E->BPI[level][el] = 1.0 / ( v1 + 1.0/comp); /**/

      }
      else {
	/* E->BPI[level][el] = visc; /**/
	  E->BPI[level][el] = 1.0/v1;  /**/
      }

    }
    else
      E->BPI[level][el] = 1.0;	

  }
  
  return;
}


void build_diagonal_of_Ahat_level1(
				  struct All_variables *E,
				  int level
				  )
{
  higher_precision assemble_dAhatq_entry();
  
  int i,j,k,el,m;
  standard_precision time,time0,CPU_time(),comp,v1;
  
  const int npno = E->mesh.npno;
  const int neq = E->mesh.neq;
  const int dims = E->mesh.nsd;
  
  time=CPU_time();
  
  if(E->control.verbose)
    fprintf(stderr,"Building pressure preconditioner %d\n",level);
  
  for(i=1;i<=E->mesh.NPNO[level];i++)
    E->BPI[level][i]=1.0;

  for(el=1;el<=E->mesh.NEL[level];el++) {
    comp = 0.0;
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > (2.0/dims)) 
	comp += 1.0/(E->tracer.BulkVisc[m] - (2.0/dims) * E->tracer.Visc[m]);  
    }

    if(E->tracer.tr_in_element_number[level][el] != 0) {
      v1 = (higher_precision) assemble_dAhatq_entry(E,el,level);
      if(comp != 0.0) {
	comp = E->tracer.tr_in_element_number[level][el]/comp;
	E->BPI[level][el] = 1.0 / ( v1 + 1.0/comp);
      }
      else {
	E->BPI[level][el] = 1.0/v1;
      }
    }
    else
      E->BPI[level][el] = 1.0;
  }
  return;
}

#if 0
void build_diagonal_of_Ahat_level(
				  struct All_variables *E,
				  int level
				  )
{
  higher_precision assemble_dAhatq_entry();
  
  int i,j,k,el,m;
  standard_precision time,time0,CPU_time(),visc,comp,v1;
  
  const int npno = E->mesh.npno;
  const int neq = E->mesh.neq;
  const int dims = E->mesh.nsd;
  
  time=CPU_time();
  
  if(E->control.verbose)
    fprintf(stderr,"Building pressure preconditioner %d\n",level);
  
  for(i=1;i<=E->mesh.NPNO[level];i++)
    E->BPI[level][i]=1.0;

  for(el=1;el<=E->mesh.NEL[level];el++) {
    visc = 0.0;
    comp = 0.0;
    
    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      visc += 1.0/(E->tracer.Visc[m]);

      if(3==dims) {
	if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.666666666) 
	  comp += 1.0/(E->tracer.BulkVisc[m] - 0.6666666 * E->tracer.Visc[m]);  
      }
      else {
	if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 1.0) 
	  comp += 1.0/(E->tracer.BulkVisc[m]  - E->tracer.Visc[m]);   
      }
    }

    if(E->tracer.tr_in_element_number[level][el] != 0) {
      visc = E->tracer.tr_in_element_number[level][el]/(visc);
  
      if(comp != 0.0) {
	comp = E->tracer.tr_in_element_number[level][el]/comp;
	E->BPI[level][el] = visc * comp / (comp+visc) ;
      }
      else {
	E->BPI[level][el] = visc;
      }
    }
    else
      E->BPI[level][el] = 1.0;
  }
  return;
}
#endif

higher_precision assemble_dAhatq_entry(
				       struct All_variables *E,
				       int e,
				       int level
				       )
{ 
  int i,j,p,a,b,node,ee,element,lnode;
  higher_precision elt_g[24][1];
  void get_elt_g();

  higher_precision gradQ[81],divU;

  const int nel=E->mesh.nel;
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int npno=E->mesh.npno;
    
  for(i=0;i<81;i++)
    gradQ[i] = 0.0;
  
  divU=0.0;
     
  for(a=1;a<=ends;a++) {
    node=E->IEN[level][e].node[a];
    p = (a-1)*dims;

    if(!(E->NODE[level][node] & BC1 )) {
      j=E->ID[level][node].doff[1];
      gradQ[p] += E->BI[level][j]*E->elt_del[level][e].g[p][0];
    }

    if(!(E->NODE[level][node] & BC2)) {
      j=E->ID[level][node].doff[2];
      gradQ[p+1] += E->BI[level][j]*E->elt_del[level][e].g[p+1][0];
    }
	    
    if(3==dims) { 
      if(!(E->NODE[level][node] & BC3)) {
	j=E->ID[level][node].doff[3];
	gradQ[p+2] += E->BI[level][j]*E->elt_del[level][e].g[p+2][0];
      }
    }
  }

  /* calculate div U from the same thing .... */

  /* only need to run over nodes with non-zero grad P, i.e. the ones in
     the element accessed above, BUT it is only necessary to update the
     value in the original element, because the diagonal is all we use at
     the end ... */

  for(b=1;b<=ends;b++) {
    p = (b-1)*dims;	   
    divU +=E->elt_del[level][e].g[p][0] * gradQ[p];	    
    divU +=E->elt_del[level][e].g[p+1][0] * gradQ[p+1];
    if(3==dims) {
      divU +=E->elt_del[level][e].g[p+2][0] * gradQ[p+2];
    }	
  }
  return(divU);
}

/*==============================================================
  Function to supply the element g matrix for a given element e.
  ==============================================================  */

void get_elt_g(
	       struct All_variables *E,
	       int el,
	       higher_precision elt_del[24][1],
	       int lev
	       )
{
  void get_global_shape_fn();
  void get_global_shape_fn_curvilin();
    
  int p,a,nint,es,d,i,j,k;
  higher_precision ra,ct,si,x[4],rtf[4][9];
  higher_precision dGNdash[3];
  higher_precision recip_radius,temp;
  int lmsize;
  
  higher_precision elt_R[24*24];
  higher_precision elt_del1[24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  
  
  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;
    
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int nn=loc_mat_size[dims][E->control.model];
  
  /* Special case, 4/8 node bilinear v/constant p  element -> 1 pressure point */

  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,2,lev);
  temp=p_point[1].weight[dims-1] * dOmega.ppt[1];
  
  /* for(a=1;a<=ends;a++)
     for(i=1;i<=dims;i++)
     { p = dims*(a-1) + i -1;
     elt_del[p][0] = -GNx.ppt[GNPXINDEX(i-1,a,1)] * temp;
     }
  */

  /* unroll loops for speed */
  for(a=1;a<=ends;a++)  {
    p=dims*(a-1);
    elt_del[p][0] = -GNx.ppt[GNPXINDEX(0,a,1)] * temp;
    elt_del[p+1][0] = -GNx.ppt[GNPXINDEX(1,a,1)] * temp;
    if(3==dims)
      elt_del[p+2][0] = -GNx.ppt[GNPXINDEX(2,a,1)] * temp;
    if(E->control.AXI)
      elt_del[p][0] += -GN.ppt[GNPINDEX(a,1)] * temp / E->ECO[lev][el].centre[1];
  }

  if (E->control.HAVE_SKEWBCS) {
       
    /* Build elt_Rot matrix */
   
    /* Identity */
    for(p=0;p<nn*nn;p++) {
      elt_R[p] = 0.0;
    }
    for(p=0;p<nn;p++)  
      elt_R[p*nn+p] = 1.0;
       
    /* Add in all nodal rotation matrix terms */
    for(a=1;a<=ends;a++) {
      node=E->IEN[lev][el].node[a];
      if(E->NODE[lev][node] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[lev][node];

	px = (a-1) * dims;
	pz = (a-1) * dims + 1;
	py = (a-1) * dims + 2;

	elt_R[px*nn+px] = node_R[0*dims+0];  /* note all nodes independent effect on elt_K */
	elt_R[px*nn+pz] = node_R[0*dims+1];
	elt_R[pz*nn+px] = node_R[1*dims+0];
	elt_R[pz*nn+pz] = node_R[1*dims+1];
	  
	if(dims==3) {
	  elt_R[px*nn+py] = node_R[0*dims+2];
	  elt_R[pz*nn+py] = node_R[1*dims+2];
	  elt_R[py*nn+px] = node_R[2*dims+0];
	  elt_R[py*nn+pz] = node_R[2*dims+1];
	  elt_R[py*nn+py] = node_R[2*dims+2];
	}
      }
    }

    /* pre multiply by RotT */
    for(i=0;i<nn;i++)   {
      sum=0.0;
      for(k=0;k<nn;k++)
	sum += elt_R[k*nn+i] * elt_del[k][0];
	  
      elt_del1[i] = sum;
    }
      
    for(i=0;i<nn;i++)
      elt_del[i][0] = elt_del1[i] ;
  }
  return; 
}


/* This is a wrapper for the various specialized element stiffness
   routines - it just interprets the prevailing settings and calls
   the appropriate function from there */


void get_elt_k(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{
  void get_elt_k_viscous_2D();
  void get_elt_k_viscous_3D();
  void get_elt_k_cosserat_2D();
  void get_elt_k_cosserat_3D();
  void get_elt_k_general();

  const int dims = E->mesh.nsd;

  /* Normally we expect to use a hand-coded
     stiffness matrix formulation. However, for
     debugging, testing and for new problems, the
     perfectly general version is available.

     By putting the tests here we can ensure that
     the code itself is significantly simpler to follow
     and more efficiently executed.
  */

  if(E->control.B_matrix!=0) {
    if(2==dims && 1==E->control.model && !E->control.AXI) { /* 2D viscous, hand-coded */
      get_elt_k_viscous_2D(E,el,elt_k,level);
    }
    else if (3==dims && 1==E->control.model) { /* 3D viscous, hand-coded */
      get_elt_k_viscous_3D(E,el,elt_k,level);
    }
    else if (2==dims && 2==E->control.model) { /* 2D Cosserat, hand-coded */
      get_elt_k_cosserat_2D(E,el,elt_k,level);
    }
    else if (3==dims && 2==E->control.model) { /* 3D Cosserat, hand-coded */
      get_elt_k_cosserat_3D(E,el,elt_k,level);
    }
    else if (2==dims && 1==E->control.model && E->control.AXI) { /* 2D, axisymmetric, hand-coded */
      /* get_elt_k_viscous_axi(E,el,elt_k,level); */
    }
    else { /* use the perfectly general case after all !! */ 
      get_elt_k_general(E,el,elt_k,level);
    }
  }
  else {
    get_elt_k_general(E,el,elt_k,level);
  }

  return;
}

void add_penalty3(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{
  int i,m,pn,qn,a,b;
  standard_precision dOmega,penalty ;
  standard_precision pen_dO;   /* DAS: 21/01/03 */
  standard_precision lNxO[4][ELNMAX+1];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];

    dOmega  = 
      E->tracer.tracer_weight_fn[level][m] 
      * E->tracer.tracer_jacobian[level][m] 
      * E->tracer.Visc[m];

    /* DAS: 21/01/03, factorised penalty and dOmega out of this loop */
    penalty=E->tracer.visc[E->tracer.property_group[m]].Pen_bulk;
    pen_dO = dOmega * penalty;
    
    for(a=1;a<=ends;a++) 
      for(b=a;b<=ends;b++) {
	
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	/* Augmented Lagrangian (penalty) to enhance incompressibility */
	if(penalty!=0.0) {
	  standard_precision    t1a = lNxO[1][a] * pen_dO,
	                        t2a = lNxO[2][a] * pen_dO;

          elt_k[pn] += t1a * lNxO[1][b];
          elt_k[pn+1] += t1a * lNxO[2][b];
          elt_k[pn+nn]  += t2a * lNxO[1][b];
          elt_k[pn+nn+1]  += t2a * lNxO[2][b];
	  if(3==dims) { /*RAA: 10/01/03, put in 3D stuff here*/
	    standard_precision    t3a = lNxO[3][a] * pen_dO;

            elt_k[pn+2] += t1a * lNxO[3][b];
            elt_k[pn+nn+2]  += t2a * lNxO[3][b];
            elt_k[pn+2*nn]  += t3a * lNxO[1][b];
            elt_k[pn+2*nn+1]  += t3a * lNxO[2][b];
            elt_k[pn+2*nn+2]  += t3a * lNxO[3][b];
	  }
	}

	/* And the matrix is symmetric so, for the time being, this is stored */
	elt_k[qn]       = elt_k[pn];
	elt_k[qn+nn]    = elt_k[pn+1];
	elt_k[qn+1]     = elt_k[pn+nn];
	elt_k[qn+nn+1]  = elt_k[pn+nn+1]; 
        if(3==dims)  { /*RAA: 10/01/03, added 3D stuff here*/
	  elt_k[qn+2*nn]  = elt_k[pn+2];
	  elt_k[qn+2*nn+1]= elt_k[pn+nn+2];
	  elt_k[qn+2]     = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]  = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2]= elt_k[pn+2*nn+2];
	}
      }
  }
  return;
}

#if 1
void get_all_elt_k
(  struct All_variables *E,
   higher_precision *elt_k,
   int level  )  {  

  void get_all_elt_k_viscous_2D();
  /* void get_all_elt_k_viscous_3D(); */

  const int dims = E->mesh.nsd;


  /* Normally we expect to use a hand-coded
     stiffness matrix formulation. However, for
     debugging, testing and for new problems, the
     perfectly general version is available.

     By putting the tests here we can ensure that
     the code itself is significantly simpler to follow
     and more efficiently executed.
  */

  if(E->control.B_matrix!=0) {
    if(2==dims && 1==E->control.model && !E->control.AXI) { /* 2D viscous, hand-coded */
      get_all_elt_k_viscous_2D(E,elt_k,level);
    }
    else if (3==dims && 1==E->control.model) { /* 3D viscous, hand-coded */
      /*  get_all_elt_k_viscous_3D(E,elt_k,level); */
    }
    else if (2==dims && 2==E->control.model) { /* 2D Cosserat, hand-coded */
      /* get_elt_k_cosserat_2D(E,el,elt_k,level); */
    }
    else if (3==dims && 2==E->control.model) { /* 3D Cosserat, hand-coded */
      /* get_elt_k_cosserat_3D(E,el,elt_k,level); */
    }
    else if (2==dims && 1==E->control.model && !E->control.AXI) { /* 2D, axisymmetric, hand-coded */
      /* get_elt_k_viscous_axi(E,el,elt_k,level); */
    }
    else { /* use the perfectly general case after all !! */ 
      /* get_elt_k_general(E,el,elt_k,level); */
    }
  }
  else {
    /* get_elt_k_general(E,el,elt_k,level); */
  }

  return;
}
#endif

void get_elt_k_viscous_2D(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2;
  standard_precision dOmega,weight,dOmegaq;
  standard_precision a1,a0,delta;
  standard_precision H,phi;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;
 
  standard_precision penalty;
    
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

  for(p=0;p<nn*nn;p++)
    elt_k[p] = elt_k1[p] = 0.0; /* !!!!!! */

/*	fprintf(stderr,"GEK - 1\n"); */

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  if(E->tracer.tr_in_element_number[level][el] == 0)
    fprintf(E->fp,"Warning: no tracers in %d/%d\n",level,el);

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
 
    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,NULL,level);

    dOmegaq = dOmega  = 
      E->tracer.tracer_weight_fn[level][m] 
      * E->tracer.tracer_jacobian[level][m] 
      * E->tracer.Visc[m];

    if(E->control.CHEM_TRANS  && 
       E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio != 1.0 ) {
      
      H = E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio;
      phi = E->tracer.volfraction[m];

      delta = 1.0 - H / ((phi * H + 1.0 - phi) * (phi + H - phi * H))  ; 

      /*fprintf(stderr,"%d: phi %g , H %g , delta %g\n",m,phi,H,delta); */
      
      a0 = 4.0 * delta *
	E->tracer.n1[m] *  E->tracer.n1[m] *  E->tracer.n2[m] *  E->tracer.n2[m];
      a1 = 2.0 * delta * E->tracer.n1[m] *  E->tracer.n2[m] *
	( E->tracer.n2[m] *  E->tracer.n2[m] - E->tracer.n1[m] *  E->tracer.n1[m] );
      
    }
    else if((E->control.ORTHOTROPY && 
       E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio < 1.0)) { 
      delta = (1.0 - E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio);
      
      a0 = 4.0 * delta *
	E->tracer.n1[m] *  E->tracer.n1[m] *  E->tracer.n2[m] *  E->tracer.n2[m];
      a1 = 2.0 * delta * E->tracer.n1[m] *  E->tracer.n2[m] *
	( E->tracer.n2[m] *  E->tracer.n2[m] - E->tracer.n1[m] *  E->tracer.n1[m] ); 
    }
    
    for(a=1;a<=ends;a++) 
      for(b=a;b<=ends;b++) {
	
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);
	
	 elt_k[pn]      += (2.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b]) * dOmega;
	 elt_k[pn+1]    +=  lNx[2][a] * lNx[1][b] * dOmega;
	 elt_k[pn+nn]   +=  lNx[1][a] * lNx[2][b] * dOmega;
	 elt_k[pn+nn+1] += (2.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b]) * dOmega; 

	 /* Orthotropic - additional terms */

	 if((E->control.ORTHOTROPY && 
	    E->tracer.visc[E->tracer.property_group[m]].Ortho_viscosity_ratio < 1.0) ||
	    (E->control.CHEM_TRANS && 
	     E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio != 1.0)
	    ) {
	   elt_k[pn] += ( -a0 * lNx[1][a] * lNx[1][b] 
			  -a1 * lNx[1][a] * lNx[2][b] 
			  -a1 * lNx[2][a] * lNx[1][b] 
			  +a0 * lNx[2][a] * lNx[2][b] 
			  - delta * lNx[2][a] * lNx[2][b] ) * dOmega;

	   elt_k[pn+1] += ( a0 * lNx[2][a] * lNx[1][b] 
			    +a1 * lNx[2][a] * lNx[2][b] 
			    -a1 * lNx[1][a] * lNx[1][b] 
			    +a0 * lNx[1][a] * lNx[2][b] 
			    - delta * lNx[2][a] * lNx[1][b] ) * dOmega;   /* corrected May 11, 2001 */

	   elt_k[pn+nn] += ( a0 * lNx[1][a] * lNx[2][b] 
			     -a1 * lNx[1][a] * lNx[1][b] 
			     +a1 * lNx[2][a] * lNx[2][b] 
			     +a0 * lNx[2][a] * lNx[1][b] 
			     - delta * lNx[1][a] * lNx[2][b] ) * dOmega; 
	   
	   elt_k[pn+nn+1] += (-a0 * lNx[2][a] * lNx[2][b] 
			     +a1 * lNx[2][a] * lNx[1][b] 
			     +a1 * lNx[1][a] * lNx[2][b] 
			     +a0 * lNx[1][a] * lNx[1][b] 
			     - delta * lNx[1][a] * lNx[1][b] ) * dOmega;
	}

	/* And the matrix is symmetric so, for the time being, this is stored */
	elt_k[qn]       = elt_k[pn];
	elt_k[qn+nn]    = elt_k[pn+1];
	elt_k[qn+1]     = elt_k[pn+nn];
	elt_k[qn+nn+1]  = elt_k[pn+nn+1]; 

   	elt_k1[qn]       = elt_k1[pn];
	elt_k1[qn+nn]    = elt_k1[pn+1];
	elt_k1[qn+1]     = elt_k1[pn+nn];
	elt_k1[qn+nn+1]  = elt_k1[pn+nn+1];
      }
  }

  if(E->control.HAVE_SKEWBCS) {    
    /* Build elt_Rot matrix */
    
    /* Identity */
    for(p=0;p<nn*nn;p++)
      elt_R[p] = 0.0;
    
    for(p=0;p<nn;p++)  
      elt_R[p*nn+p] = 1.0;
    
    /* Add in all nodal rotation matrix terms */
    for(a=1;a<=ends;a++) {/*5*/
      node=E->IEN[level][el].node[a];
      if(E->NODE[level][node] & SKEWBC) {
	node_R = E->curvilinear.NODE_R[level][node];
	
	px = (a-1) * dims;
	pz = (a-1) * dims + 1;
	py = (a-1) * dims + 2;
	
	elt_R[px*nn+px] = node_R[0*dims+0]; 
	/* note all nodes independent effect on elt_K */
	elt_R[px*nn+pz] = node_R[0*dims+1];
	elt_R[pz*nn+px] = node_R[1*dims+0];
	elt_R[pz*nn+pz] = node_R[1*dims+1];
	
      }
    }
    
    /* post multiply by Rot */
    for(i=0;i<nn;i++) 
      for(j=0;j<nn;j++) {
	sum=0.0;
	for(k=0;k<nn;k++)
	  sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	
	elt_k1[i*nn+j] = sum;
      }
    
    /* pre multiply by RotT */
    for(i=0;i<nn;i++)  
      for(j=0;j<nn;j++) {
	sum=0.0;
	for(k=0;k<nn;k++)
	  sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	
	elt_k[i*nn+j] = sum;
      }
  }

  return;
}

void get_all_elt_k_viscous_2D(
			      struct All_variables *E,    
			      higher_precision *Elt_K,
			      int level
			      )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;  
  int el;

  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight,dOmegaq;
  standard_precision t1;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision *elt_k;
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  
  
  standard_precision penalty;
    
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];
  const int elts=E->mesh.NEL[level];
  const int eksize = dofs*ends*dofs*ends;
  const int numtrc = E->tracer.NUM_TRACERS;




  for(p=0;p<elts*eksize;p++)
    Elt_K[p] = 0.0;


  for(el=1;el<=elts;el++) {

    elt_k = Elt_K + (el-1) * eksize;
    get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

    for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
      m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
      eta1 = E->tracer.eta1[level][m];
      eta2 = E->tracer.eta2[level][m];
      if (3==dims)
	eta3 = E->tracer.eta3[level][m];
 
      get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);
  
      dOmegaq = dOmega  = 
	E->tracer.tracer_weight_fn[level][m] 
	* E->tracer.tracer_jacobian[level][m] 
	* E->tracer.Visc[m];
  
    
      /* Algebraic expression for this ?? */
      if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio >= 0.0)
	penalty = max(0.0,min(E->control.AUG_lagrangian,E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio));
      else 
	penalty = E->control.AUG_lagrangian;


      t1 = dOmega; 

#pragma loop novrec elt_k,lNx,lNxO
#pragma loop temp a,b,pn
      /*  for(a=1;a<=ends;a++) 
	  for(b=1;b<=ends;b++) { */
      for(k=0;k<ends*ends;k++) {
	a = k/ends+1;
	b = k%ends+1;
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	/*	
	elt_k[pn]      += t1 * (2.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b] +
				penalty * (lNxO[1][a] * lNxO[1][b]));
	elt_k[pn+1]    += t1 * (lNx[2][a] * lNx[1][b] + 
				penalty * (lNxO[1][a] * lNxO[2][b]));
	elt_k[pn+nn]   += t1 * (lNx[1][a] * lNx[2][b] +
				penalty * (lNxO[2][a] * lNxO[1][b]));
	elt_k[pn+nn+1] += t1 * (2.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b] +
				penalty * (lNxO[2][a] * lNxO[2][b]));
	*/

	/* deviatoric */

	elt_k[pn]      += t1 * (1.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b]);
	elt_k[pn+1]    += t1 * (lNx[2][a] * lNx[1][b] - 1.0 * lNx[1][a] * lNx[2][b]);
	elt_k[pn+nn]   += t1 * (lNx[1][a] * lNx[2][b] - 1.0 * lNx[2][a] * lNx[1][b]);
	elt_k[pn+nn+1] += t1 * (1.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b]);

      }
    }

#if 0
    if(E->control.HAVE_SKEWBCS) {

      /* Build elt_Rot matrix */
   
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;
      
      /* Add in all nodal rotation matrix terms */
      for(a=1;a<=ends;a++) {/*5*/
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];
	
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
    }
#endif  

  }


  return;
}


void get_elt_k_viscous_3D(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  /*standard_precision lNxO[4][ELNMAX+1];*/ /*RAA: 22-01-03, no longer needed*/
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight,dOmegaq;
  /*higher_precision del_elt_k[24*24];*/ /*RAA: 22-01-03, not needed*/
  higher_precision elt_R[24*24];
  /*higher_precision cos_theta,sin_theta;*/ /*RAA: 22-01-03, not needed*/
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  

  /*standard_precision B, Gc ;*/ /*RAA: 22-01-03, not needed*/
  /*standard_precision Gdig, Gnotdig ;*/ /*RAA: 22-01-03, not needed*/
  /*standard_precision penalty;*/ /*RAA: 22-01-03, no longer needed*/
   
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

 
  for(p=0;p<nn*nn;p++)
    elt_k[p] = 0.0;

  /*penalty = E->control.AUG_lagrangian;  */ /*RAA: 22-01-03, no longer needed*/
  /*get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);*/ /*RAA: 22-01-03, no longer needed*/

  if(E->tracer.tr_in_element_number[level][el] == 0) /*RAA: 22-01-03, added from 2D*/
    fprintf(E->fp,"Warning: no tracers in %d/%d\n",level,el);

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    eta3 = E->tracer.eta3[level][m];

    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);

     dOmegaq = dOmega  = 
      E->tracer.tracer_weight_fn[level][m] 
      * E->tracer.tracer_jacobian[level][m] 
      * E->tracer.Visc[m];
      
    for(a=1;a<=ends;a++)
      for(b=a;b<=ends;b++) {
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	/* VISCOUS 3D case */
	
	elt_k[pn]        += (2.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b] + lNx[3][a] * lNx[3][b]) * dOmega;
	elt_k[pn+1]      +=  lNx[2][a] * lNx[1][b] * dOmega;
	elt_k[pn+2]      +=  lNx[3][a] * lNx[1][b] * dOmega;
	elt_k[pn+nn]     +=  lNx[1][a] * lNx[2][b] * dOmega;
	elt_k[pn+nn+1]   += (2.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b] + lNx[3][a] * lNx[3][b]) * dOmega;
	elt_k[pn+nn+2]   +=  lNx[3][a] * lNx[2][b] * dOmega;
        elt_k[pn+2*nn]   +=  lNx[1][a] * lNx[3][b] * dOmega;
        elt_k[pn+2*nn+1] +=  lNx[2][a] * lNx[3][b] * dOmega;
       	elt_k[pn+2*nn+2] += (2.0 * lNx[3][a] * lNx[3][b] + lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b]) * dOmega;
	
	/*RAA: 10/01/03 - removed the penalty stuff that was here --> into add_penalty3( ) */

	  /* And the matrix is symmetric so, for the time being, this is stored */
	  elt_k[qn]       = elt_k[pn];
	  elt_k[qn+nn]    = elt_k[pn+1];
	  elt_k[qn+2*nn]  = elt_k[pn+2];
	  elt_k[qn+1]     = elt_k[pn+nn];
	  elt_k[qn+nn+1]  = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1]= elt_k[pn+nn+2];
	  elt_k[qn+2]     = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]  = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2]= elt_k[pn+2*nn+2];
      }
      
    if(E->control.HAVE_SKEWBCS) {

      /* Build elt_Rot matrix */
   
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;

      /* Add in all nodal rotation matrix terms */
	for(a=1;a<=ends;a++) {/*5*/
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];
	  elt_R[px*nn+py] = node_R[0*dims+2];
	  elt_R[pz*nn+py] = node_R[1*dims+2];
	  elt_R[py*nn+px] = node_R[2*dims+0];
	  elt_R[py*nn+pz] = node_R[2*dims+1];
	  elt_R[py*nn+py] = node_R[2*dims+2];
	  
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
  }
  }
    
  return;
}

void get_elt_k_cosserat_2D(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1],lNO[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2;
  standard_precision dOmega,weight,dOmegaq;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;

  standard_precision one_over_x2,one_over_x;
  standard_precision A11,A12,A22,G11,G12,G22,B31,B32; 
  standard_precision penalty;
    
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

  for(p=0;p<nn*nn;p++)
    elt_k[p] = 0.0;

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  if(E->tracer.tr_in_element_number[level][el] == 0)
    fprintf(stderr,"Warning (error, really): no tracers in %d/%d\n",level,el);

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
 
    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,NULL,level);
  
      dOmegaq = dOmega  = 
	E->tracer.tracer_weight_fn[level][m] 
	* E->tracer.tracer_jacobian[level][m] 
	* E->tracer.Visc[m];
 
    penalty = 0.0 ;
    v_shape_fn(E,el,lN,eta1,eta2,NULL,level);
    v_shape_fn(E,el,lNO,0,0,NULL,level);
    one_over_x = 1.0 / E->tracer.tx[m];
    one_over_x2 = one_over_x * one_over_x;

    A11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A11 - dOmegaq * penalty ;
    A22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A22 - dOmegaq * penalty ;
    A12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A12 - dOmegaq * penalty ;
    G11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G11 ;
    G22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G22 ;
    G12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G12 ;
    B31 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B31 ;
    B32 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B32 ;

    for(a=1;a<=ends;a++) 
      for(b=a;b<=ends;b++) { 
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	/* original */
/*	elt_k[pn]        += lNx[1][a] * lNx[1][b] * A11 + G11 * lNx[2][a] * lNx[2][b] ;
	elt_k[pn+1]      += G12 * lNx[2][a] * lNx[1][b] ;
	elt_k[pn+2]      += (G11 - G12) * lNx[2][a] * lN[b] ;
	elt_k[pn+nn]     += G12 * lNx[1][a] * lNx[2][b] ;
	elt_k[pn+nn+1]   += lNx[2][a] * lNx[2][b] * A22 + G22 * lNx[1][a] * lNx[1][b] ;
	elt_k[pn+nn+2]   += (G12 - G22) * lNx[1][a] * lN[b] ;
	elt_k[pn+2*nn]   += (G11 - G12) * lNx[2][b] * lN[a] ;
	elt_k[pn+2*nn+1] += (G12 - G22) * lNx[1][b] * lN[a] ;
	elt_k[pn+2*nn+2] += (G11 - 2*G12 + G22) * lN[b] * lN[a] + B31 * lNx[1][a] * lNx[1][b] + B32 * lNx[2][a] * lNx[2][b] ;*/
	/* original  + A12 */
/*	elt_k[pn]        += lNx[1][a] * lNx[1][b] * A11 + G11 * lNx[2][a] * lNx[2][b] ;
	elt_k[pn+1]      += G12 * lNx[2][a] * lNx[1][b] + lNx[2][b] * lNx[1][a] * A12 ;
	elt_k[pn+2]      += (G11 - G12) * lNx[2][a] * lN[b] ;
	elt_k[pn+nn]     += G12 * lNx[1][a] * lNx[2][b] +  A12 * lNx[2][a] * lNx[1][b] ;
	elt_k[pn+nn+1]   += lNx[2][a] * lNx[2][b] * A22 + G22 * lNx[1][a] * lNx[1][b] ;
	elt_k[pn+nn+2]   += (G12 - G22) * lNx[1][a] * lN[b] ;
	elt_k[pn+2*nn]   += (G11 - G12) * lNx[2][b] * lN[a] ;
	elt_k[pn+2*nn+1] += (G12 - G22) * lNx[1][b] * lN[a] ;
	elt_k[pn+2*nn+2] += (G11 - 2*G12 + G22) * lN[b] * lN[a] + B31 * lNx[1][a] * lNx[1][b] + B32 * lNx[2][a] * lNx[2][b] ;*/
	/* original + A12 + constante rotation */
	elt_k[pn]        += lNx[1][a] * lNx[1][b] * A11 + G11 * lNx[2][a] * lNx[2][b] ;
	elt_k[pn+1]      += G12 * lNx[2][a] * lNx[1][b] + lNx[2][b] * lNx[1][a] * A12 ;
	elt_k[pn+2]      += (G11 - G12) * lNx[2][a] * lNO[b] ;
	elt_k[pn+nn]     += G12 * lNx[1][a] * lNx[2][b] +  A12 * lNx[2][a] * lNx[1][b] ;
	elt_k[pn+nn+1]   += lNx[2][a] * lNx[2][b] * A22 + G22 * lNx[1][a] * lNx[1][b] ;
	elt_k[pn+nn+2]   += (G12 - G22) * lNx[1][a] * lNO[b] ;
	elt_k[pn+2*nn]   += (G11 - G12) * lNx[2][b] * lNO[a] ;
	elt_k[pn+2*nn+1] += (G12 - G22) * lNx[1][b] * lNO[a] ;
	elt_k[pn+2*nn+2] += (G11 - 2*G12 + G22) * lNO[b] * lNO[a] + 
	  B31 * lNxO[1][a] * lNxO[1][b] + B32 * lNxO[2][a] * lNxO[2][b] ;

	/* Augmented Lagrangian (penalty) to enhance incompressibility */
	if(0 && penalty!=0.0) {
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;
	}
	elt_k[qn]        = elt_k[pn];
	elt_k[qn+nn]     = elt_k[pn+1];
	elt_k[qn+2*nn]   = elt_k[pn+2];
	elt_k[qn+1]      = elt_k[pn+nn];
	elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	elt_k[qn+2]      = elt_k[pn+2*nn];
	elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
      }
  }

    if(E->control.HAVE_SKEWBCS) {

      /* Build elt_Rot matrix */
   
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;
      
      /* Add in all nodal rotation matrix terms */
      for(a=1;a<=ends;a++) {
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];	
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
    }
  return;
}


void get_elt_k_cosserat_3D(
			   struct All_variables *E,
			   int el,
			   higher_precision *elt_k,
			   int level
			   )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight,dOmegaq;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  

  standard_precision B, Gc ;
  standard_precision penalty;
  standard_precision one_over_x2,one_over_x;
  standard_precision A11,A12,A13,A22,A23,A33,G11,G12,G13,G22,G23,G33,B12,B21,B13,B31,B23,B32;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];
 
  for(p=0;p<nn*nn;p++)
    elt_k[p] = 0.0;

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    eta3 = E->tracer.eta3[level][m];

    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);

      dOmegaq = dOmega  = 
	E->tracer.tracer_weight_fn[level][m] 
	* E->tracer.tracer_jacobian[level][m] 
	* E->tracer.Visc[m];
 

    v_shape_fn(E,el,lN,eta1,eta2,eta3,level);
    one_over_x = 1.0 / E->tracer.tx[m];
    one_over_x2 = one_over_x * one_over_x;

    A11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A11 - dOmegaq * penalty ;
    A22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A22 - dOmegaq * penalty ;
    A33 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A33 - dOmegaq * penalty ;
    A12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A12 - dOmegaq * penalty ;
    A13 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A13 - dOmegaq * penalty ;
    A23 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].A23 - dOmegaq * penalty ;
    G11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G11 ;
    G22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G22 ;
    G33 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G33 ;
    G12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G12 ;
    G13 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G13 ;
    G23 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].G23 ;
    B12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B12 ;
    B21 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B21 ;
    B13 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B13 ;
    B31 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B31 ;
    B23 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B23 ;
    B32 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
      E->tracer.coss[E->tracer.property_group[m]].B32 ;

    for(a=1;a<=ends;a++)
      for(b=a;b<=ends;b++) {
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	elt_k[pn]        += 1.0 ;
	elt_k[pn+1]      += 1.0 ;
	elt_k[pn+2]      += 1.0 ;
	elt_k[pn+3]      += 0.0 ;
	elt_k[pn+4]      += 1.0 ;
	elt_k[pn+5]      += 1.0 ;
	elt_k[pn+nn]     += 1.0 ;
	elt_k[pn+nn+1]   += 1.0 ;
	elt_k[pn+nn+2]   += 1.0 ;
	elt_k[pn+nn+3]   += 1.0 ;
	elt_k[pn+nn+4]   += 0.0 ;
	elt_k[pn+nn+5]   += 1.0 ;
	elt_k[pn+2*nn]   += 1.0 ;
	elt_k[pn+2*nn+1] += 1.0 ;
	elt_k[pn+2*nn+2] += 1.0 ;
	elt_k[pn+2*nn+3] += 1.0 ;
	elt_k[pn+2*nn+4] += 1.0 ;
	elt_k[pn+2*nn+5] += 0.0 ;
	elt_k[pn+3*nn]   += 0.0 ;
	elt_k[pn+3*nn+1] += 1.0 ;
	elt_k[pn+3*nn+2] += 1.0 ;
	elt_k[pn+3*nn+3] += 1.0 ;
	elt_k[pn+3*nn+4] += 0.0 ;
	elt_k[pn+3*nn+5] += 0.0 ;
	elt_k[pn+4*nn]   += 1.0 ;
	elt_k[pn+4*nn+1] += 0.0 ;
	elt_k[pn+4*nn+2] += 1.0 ;
	elt_k[pn+4*nn+3] += 0.0 ;
	elt_k[pn+4*nn+4] += 1.0 ;
	elt_k[pn+4*nn+5] += 0.0 ;
	elt_k[pn+5*nn]   += 1.0 ;
	elt_k[pn+5*nn+1] += 1.0 ;
	elt_k[pn+5*nn+2] += 0.0 ;
	elt_k[pn+5*nn+3] += 0.0 ;
	elt_k[pn+5*nn+4] += 0.0 ;
	elt_k[pn+5*nn+5] += 1.0 ;
	
	/* Augmented Lagrangian (penalty) to enhance incompressibility */
	if(0 && penalty!=0.0) {
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2] += penalty * (lNxO[1][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn+2]  += penalty * (lNxO[2][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+2*nn]  += penalty * (lNxO[3][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+2*nn+1]  += penalty * (lNxO[3][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2*nn+2]  += penalty * (lNxO[3][a] * lNxO[3][b]) * dOmegaq;
	}
	elt_k[qn]        = elt_k[pn];
	elt_k[qn+nn]     = elt_k[pn+1];
	elt_k[qn+2*nn]   = elt_k[pn+2];
	elt_k[qn+3*nn]   = elt_k[pn+3];
	elt_k[qn+4*nn]   = elt_k[pn+4];
	elt_k[qn+5*nn]   = elt_k[pn+5];
	elt_k[qn+1]      = elt_k[pn+nn];
	elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	elt_k[qn+3*nn+1] = elt_k[pn+nn+3];
	elt_k[qn+4*nn+1] = elt_k[pn+nn+4];
	elt_k[qn+5*nn+1] = elt_k[pn+nn+5];
	elt_k[qn+2]      = elt_k[pn+2*nn];
	elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	elt_k[qn+3*nn+2] = elt_k[pn+2*nn+3];
	elt_k[qn+4*nn+2] = elt_k[pn+2*nn+4];
	elt_k[qn+5*nn+2] = elt_k[pn+2*nn+5];
	elt_k[qn+3]      = elt_k[pn+3*nn];
	elt_k[qn+nn+3]   = elt_k[pn+3*nn+1];
	elt_k[qn+2*nn+3] = elt_k[pn+3*nn+2];
	elt_k[qn+3*nn+3] = elt_k[pn+3*nn+3];
	elt_k[qn+4*nn+3] = elt_k[pn+3*nn+4];
	elt_k[qn+5*nn+3] = elt_k[pn+3*nn+5];
	elt_k[qn+4]      = elt_k[pn+4*nn];
	elt_k[qn+nn+4]   = elt_k[pn+4*nn+1];
	elt_k[qn+2*nn+4] = elt_k[pn+4*nn+2];
	elt_k[qn+3*nn+4] = elt_k[pn+4*nn+3];
	elt_k[qn+4*nn+4] = elt_k[pn+4*nn+4];
	elt_k[qn+5*nn+4] = elt_k[pn+4*nn+5];
	elt_k[qn+5]      = elt_k[pn+5*nn];
	elt_k[qn+nn+5]   = elt_k[pn+5*nn+1];
	elt_k[qn+2*nn+5] = elt_k[pn+5*nn+2];
	elt_k[qn+3*nn+5] = elt_k[pn+5*nn+3];
	elt_k[qn+4*nn+5] = elt_k[pn+5*nn+4];
	elt_k[qn+5*nn+5] = elt_k[pn+5*nn+5];
      }
      
    if(E->control.HAVE_SKEWBCS) {
      /* Build elt_Rot matrix */
      
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;

      /* Add in all nodal rotation matrix terms */
	for(a=1;a<=ends;a++) {
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];
	  elt_R[px*nn+py] = node_R[0*dims+2];
	  elt_R[pz*nn+py] = node_R[1*dims+2];
	  elt_R[py*nn+px] = node_R[2*dims+0];
	  elt_R[py*nn+pz] = node_R[2*dims+1];
	  elt_R[py*nn+py] = node_R[2*dims+2];
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
    }
  }
  return;
}
#if 0
void get_elt_k1(
	       struct All_variables *E,
	       int el,
	       higher_precision *elt_k,
	       int level
	       )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight,dOmegaq;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  

  standard_precision B, Gc ;
  standard_precision Gdig, Gnotdig ;

  standard_precision penalty;
  standard_precision Visc;
  standard_precision one_over_x2;
  standard_precision one_over_x;

  standard_precision A11,A22,G11,G12,G22,B31,B32;
  
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

  standard_precision DV[16][16];
  standard_precision DP[16][16];
  standard_precision BaT_DV_Bb[7][7];
  standard_precision BaT_DP_Bb[7][7];
  standard_precision BaT_DV[7][16];
  standard_precision Ba[16][7];
  standard_precision Bb[16][7];

  /* FILE *fp ;
  char *name  ;
  name="elt_k_test" ;
  fp = fopen(name,"w") ;
  */

  for(p=0;p<nn*nn;p++)
    elt_k[p] = 0.0;
  penalty = E->control.AUG_lagrangian;  

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  Visc=dOmega=0.0;
  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    if(E->tracer.property_group[m] < 0)
      continue;

    dOmega += E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];
    Visc += E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
  }
  if(dOmega != 0.0)
    Visc /= dOmega;

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    /* if(E->tracer.property_group[m] < 0)
      continue;
    */
  
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    if(3==dims) 
      eta3 = E->tracer.eta3[level][m];

    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);
  
    if(E->control.AXI || E->control.model==2) {
      v_shape_fn(E,el,lN,eta1,eta2,eta3,level);
      one_over_x = 1.0 / E->tracer.tx[m];  
      one_over_x2 = one_over_x * one_over_x;    
    }

    dOmegaq = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
    dOmega  = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
  
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio >= 0.0)
      penalty = max(0.0,min(E->control.AUG_lagrangian,E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio));
    else 
      penalty = E->control.AUG_lagrangian;

    /* dOmega is modified in axisymmetric case */

    if(E->control.AXI) {
      dOmega *= 2 * M_PI * E->tracer.tx[m];
    }

    for(a=1;a<=ends;a++)
      for(b=a;b<=ends;b++) {
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	if(E->control.model==2) {
	  A11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].A11 - dOmegaq * penalty ;
	  A22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].A22 - dOmegaq * penalty ;
	  G11 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].G11 ;
	  G22 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].G22 ;
	  G12 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].G12 ;
	  B31 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].B31 ;
	  B32 = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] *
	    E->tracer.coss[E->tracer.property_group[m]].B32 ;

	  Gc = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] 
	    * E->tracer.coss[E->tracer.property_group[m]].Shear_modulus ;
	  B = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] 
	    * E->tracer.coss[E->tracer.property_group[m]].Bend_stiff ;
	  Gdig = dOmega - Gc ;
	  Gnotdig = dOmega + Gc ;
	}
/*	fprintf(stderr,"A11 = %g , A22 = %g , G11 = %g , G12 = %g , G22 = %g ,B1 = %g and B2 = %g \n",
		A11,A22,G11,G12,G22,B1,B2) ;*/

	/* Supposed to be the most general formulation */

	if(E->control.B_matrix==0) {

	  construct_Bi(E,lNx,a,Ba,one_over_x,lN) ;
	  construct_D_Vel(E,dOmega,DV,B,Gdig,Gnotdig) ;
	  construct_Bi(E,lNx,b,Bb,one_over_x,lN) ;

	  mul_At_B(Ba,DV,BaT_DV) ;
	  mul_A_B(BaT_DV,Bb,BaT_DV_Bb) ;

	  construct_Bi(E,lNxO,a,Ba,one_over_x,lN) ;
	  construct_D_Pres(E,dOmegaq,dOmega,penalty,DP) ;
	  construct_Bi(E,lNxO,b,Bb,one_over_x,lN) ;

	  mul_At_B(Ba,DP,BaT_DV) ;
	  mul_A_B(BaT_DV,Bb,BaT_DP_Bb) ;

	  switch(dofs) {
	  case 2:
	  elt_k[pn]      += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	  elt_k[pn+1]    += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	  elt_k[pn+nn]   += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	  elt_k[pn+nn+1] += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;

	  /* And the matrix is symmetric so, for the time being, this is stored */
	  elt_k[qn]       = elt_k[pn];
	  elt_k[qn+nn]    = elt_k[pn+1];
	  elt_k[qn+1]     = elt_k[pn+nn];
	  elt_k[qn+nn+1]  = elt_k[pn+nn+1];
	    break ;
	  case 3:
	    elt_k[pn]        += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	    elt_k[pn+1]      += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	    elt_k[pn+2]      += BaT_DV_Bb[1][3] + BaT_DP_Bb[1][3] ;
	    elt_k[pn+nn]     += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	    elt_k[pn+nn+1]   += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;
	    elt_k[pn+nn+2]   += BaT_DV_Bb[2][3] + BaT_DP_Bb[2][3] ;
	    elt_k[pn+2*nn]   += BaT_DV_Bb[3][1] + BaT_DP_Bb[3][1] ;
	    elt_k[pn+2*nn+1] += BaT_DV_Bb[3][2] + BaT_DP_Bb[3][2] ;
	    elt_k[pn+2*nn+2] += BaT_DV_Bb[3][3] + BaT_DP_Bb[3][3] ;

	    elt_k[qn]        = elt_k[pn];
	    elt_k[qn+nn]     = elt_k[pn+1];
	    elt_k[qn+2*nn]   = elt_k[pn+2];
	    elt_k[qn+1]      = elt_k[pn+nn];
	    elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	    elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	    elt_k[qn+2]      = elt_k[pn+2*nn];
	    elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	    elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	    break ;
	  case 6:
	    elt_k[pn]        += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	    elt_k[pn+1]      += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	    elt_k[pn+2]      += BaT_DV_Bb[1][3] + BaT_DP_Bb[1][3] ;
	    elt_k[pn+3]      += BaT_DV_Bb[1][4] + BaT_DP_Bb[1][4] ;
	    elt_k[pn+4]      += BaT_DV_Bb[1][5] + BaT_DP_Bb[1][5] ;
	    elt_k[pn+5]      += BaT_DV_Bb[1][6] + BaT_DP_Bb[1][6] ;
	    elt_k[pn+nn]     += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	    elt_k[pn+nn+1]   += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;
	    elt_k[pn+nn+2]   += BaT_DV_Bb[2][3] + BaT_DP_Bb[2][3] ;
	    elt_k[pn+nn+3]   += BaT_DV_Bb[2][4] + BaT_DP_Bb[2][4] ;
	    elt_k[pn+nn+4]   += BaT_DV_Bb[2][5] + BaT_DP_Bb[2][5] ;
	    elt_k[pn+nn+5]   += BaT_DV_Bb[2][6] + BaT_DP_Bb[2][6] ;
	    elt_k[pn+2*nn]   += BaT_DV_Bb[3][1] + BaT_DP_Bb[3][1] ;
	    elt_k[pn+2*nn+1] += BaT_DV_Bb[3][2] + BaT_DP_Bb[3][2] ;
	    elt_k[pn+2*nn+2] += BaT_DV_Bb[3][3] + BaT_DP_Bb[3][3] ;
	    elt_k[pn+2*nn+3] += BaT_DV_Bb[3][4] + BaT_DP_Bb[3][4] ;
	    elt_k[pn+2*nn+4] += BaT_DV_Bb[3][5] + BaT_DP_Bb[3][5] ;
	    elt_k[pn+2*nn+5] += BaT_DV_Bb[3][6] + BaT_DP_Bb[3][6] ;
	    elt_k[pn+3*nn]   += BaT_DV_Bb[4][1] + BaT_DP_Bb[4][1] ;
	    elt_k[pn+3*nn+1] += BaT_DV_Bb[4][2] + BaT_DP_Bb[4][2] ;
	    elt_k[pn+3*nn+2] += BaT_DV_Bb[4][3] + BaT_DP_Bb[4][3] ;
	    elt_k[pn+3*nn+3] += BaT_DV_Bb[4][4] + BaT_DP_Bb[4][4] ;
	    elt_k[pn+3*nn+4] += BaT_DV_Bb[4][5] + BaT_DP_Bb[4][5] ;
	    elt_k[pn+3*nn+5] += BaT_DV_Bb[4][6] + BaT_DP_Bb[4][6] ;
	    elt_k[pn+4*nn]   += BaT_DV_Bb[5][1] + BaT_DP_Bb[5][1] ;
	    elt_k[pn+4*nn+1] += BaT_DV_Bb[5][2] + BaT_DP_Bb[5][2] ;
	    elt_k[pn+4*nn+2] += BaT_DV_Bb[5][3] + BaT_DP_Bb[5][3] ;
	    elt_k[pn+4*nn+3] += BaT_DV_Bb[5][4] + BaT_DP_Bb[5][4] ;
	    elt_k[pn+4*nn+4] += BaT_DV_Bb[5][5] + BaT_DP_Bb[5][5] ;
	    elt_k[pn+4*nn+5] += BaT_DV_Bb[5][6] + BaT_DP_Bb[5][6] ;
	    elt_k[pn+5*nn]   += BaT_DV_Bb[6][1] + BaT_DP_Bb[6][1] ;
	    elt_k[pn+5*nn+1] += BaT_DV_Bb[6][2] + BaT_DP_Bb[6][2] ;
	    elt_k[pn+5*nn+2] += BaT_DV_Bb[6][3] + BaT_DP_Bb[6][3] ;
	    elt_k[pn+5*nn+3] += BaT_DV_Bb[6][4] + BaT_DP_Bb[6][4] ;
	    elt_k[pn+5*nn+4] += BaT_DV_Bb[6][5] + BaT_DP_Bb[6][5] ;
	    elt_k[pn+5*nn+5] += BaT_DV_Bb[6][6] + BaT_DP_Bb[6][6] ;

	    elt_k[qn]        = elt_k[pn];
	    elt_k[qn+nn]     = elt_k[pn+1];
	    elt_k[qn+2*nn]   = elt_k[pn+2];
	    elt_k[qn+3*nn]   = elt_k[pn+3];
	    elt_k[qn+4*nn]   = elt_k[pn+4];
	    elt_k[qn+5*nn]   = elt_k[pn+5];
	    elt_k[qn+1]      = elt_k[pn+nn];
	    elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	    elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	    elt_k[qn+3*nn+1] = elt_k[pn+nn+3];
	    elt_k[qn+4*nn+1] = elt_k[pn+nn+4];
	    elt_k[qn+5*nn+1] = elt_k[pn+nn+5];
	    elt_k[qn+2]      = elt_k[pn+2*nn];
	    elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	    elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	    elt_k[qn+3*nn+2] = elt_k[pn+2*nn+3];
	    elt_k[qn+3*nn+2] = elt_k[pn+2*nn+4];
	    elt_k[qn+5*nn+2] = elt_k[pn+2*nn+5];
	    elt_k[qn+3]      = elt_k[pn+3*nn];
	    elt_k[qn+nn+3]   = elt_k[pn+3*nn+1];
	    elt_k[qn+2*nn+3] = elt_k[pn+3*nn+2];
	    elt_k[qn+3*nn+3] = elt_k[pn+3*nn+3];
	    elt_k[qn+4*nn+3] = elt_k[pn+3*nn+4];
	    elt_k[qn+5*nn+3] = elt_k[pn+3*nn+5];
	    elt_k[qn+4]      = elt_k[pn+4*nn];
	    elt_k[qn+nn+4]   = elt_k[pn+4*nn+1];
	    elt_k[qn+2*nn+4] = elt_k[pn+4*nn+2];
	    elt_k[qn+3*nn+4] = elt_k[pn+4*nn+3];
	    elt_k[qn+4*nn+4] = elt_k[pn+4*nn+4];
	    elt_k[qn+5*nn+4] = elt_k[pn+4*nn+5];
	    elt_k[qn+5]      = elt_k[pn+5*nn];
	    elt_k[qn+nn+5]   = elt_k[pn+5*nn+1];
	    elt_k[qn+2*nn+5] = elt_k[pn+5*nn+2];
	    elt_k[qn+3*nn+5] = elt_k[pn+5*nn+3];
	    elt_k[qn+3*nn+5] = elt_k[pn+5*nn+4];
	    elt_k[qn+5*nn+5] = elt_k[pn+5*nn+5];
	    break ;
	  }
	}
/* End of general formulation */

	/* VISCOUS 2D case .... will be extended to 3D as necessary */

	if(E->control.B_matrix==1) {
	  elt_k[pn]      += (2.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b]) * dOmega;
	  elt_k[pn+1]    +=  lNx[2][a] * lNx[1][b] * dOmega;
	  elt_k[pn+nn]   +=  lNx[1][a] * lNx[2][b] * dOmega;
	  elt_k[pn+nn+1] += (2.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b]) * dOmega;

	  /* Axisymmetric case ... 1 additional term in stiffness matrix  */
	  if(E->control.AXI)
	    elt_k[pn] += 2.0 * (lN[a] * lN[b]) * one_over_x2*dOmega;
	
	  /* Augmented Lagrangian (penalty) to enhance incompressibility */
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;

	  /* Axisymmetric case ... 3 additional terms in constraint matrix */
	  if(E->control.AXI) {
	    elt_k[pn] += penalty * (lNxO[1][a] * lN[b] * one_over_x +
				    lNxO[1][b] * lN[a] * one_over_x +
				    lN[a] * lN[b] * one_over_x2) * dOmegaq ;
	    elt_k[pn+1] += penalty * (lN[a]  * lNxO[2][b] * one_over_x) * dOmegaq ;
	    elt_k[pn+nn]  += penalty * (lN[b]  * lNxO[2][a] * one_over_x) * dOmegaq  ;	
	  }

	  /* And the matrix is symmetric so, for the time being, this is stored */
	  elt_k[qn]       = elt_k[pn];
	  elt_k[qn+nn]    = elt_k[pn+1];
	  elt_k[qn+1]     = elt_k[pn+nn];
	  elt_k[qn+nn+1]  = elt_k[pn+nn+1];
	}

	/* 2D Cosserat case */
	if(E->control.B_matrix==2) {
	  elt_k[pn]        += lNx[1][a] * lNx[1][b] * A11 + G11 * lNx[2][a] * lNx[2][b] ;
	  elt_k[pn+1]      += G12 * lNx[2][a] * lNx[1][b] ;
	  elt_k[pn+2]      += (G11 - G12) * lNx[2][a] * lN[b] ;
	  elt_k[pn+nn]     += G12 * lNx[1][a] * lNx[2][b] ;
	  elt_k[pn+nn+1]   += lNx[2][a] * lNx[2][b] * A22 + G22 * lNx[1][a] * lNx[1][b] ;
	  elt_k[pn+nn+2]   += (G12 - G22) * lNx[1][a] * lN[b] ;
	  elt_k[pn+2*nn]   += (G11 - G12) * lNx[2][b] * lN[a] ;
	  elt_k[pn+2*nn+1] += (G12 - G22) * lNx[1][b] * lN[a] ;
	  elt_k[pn+2*nn+2] += (G11 - 2*G12 + G22) * lN[b] * lN[a] + B1 * lNx[1][a] * lNx[1][b] + B2 * lNx[2][a] * lNx[2][b] ;

	  /* Augmented Lagrangian (penalty) to enhance incompressibility */
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;

	  elt_k[qn]        = elt_k[pn];
	  elt_k[qn+nn]     = elt_k[pn+1];
	  elt_k[qn+2*nn]   = elt_k[pn+2];
	  elt_k[qn+1]      = elt_k[pn+nn];
	  elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	  elt_k[qn+2]      = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	}

	/* VISCOUS 3D case */
	if(E->control.B_matrix==3) {
	  elt_k[pn]      += (2.0 * lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b] + lNx[3][a] * lNx[3][b]) * dOmega;
	  elt_k[pn+1]    +=  lNx[2][a] * lNx[1][b] * dOmega;
	  elt_k[pn+2]    +=  lNx[3][a] * lNx[1][b] * dOmega;
	  elt_k[pn+nn]   +=  lNx[1][a] * lNx[2][b] * dOmega;
	  elt_k[pn+nn+1] += (2.0 * lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b] + lNx[3][a] * lNx[3][b]) * dOmega;
	  elt_k[pn+nn+2] += lNx[3][a] * lNx[2][b] * dOmega;
	  elt_k[pn+2*nn]   +=  lNx[1][a] * lNx[3][b] * dOmega;
	  elt_k[pn+2*nn+1] += lNx[2][a] * lNx[3][b] * dOmega;
	  elt_k[pn+2*nn+2] += (2.0 * lNx[3][a] * lNx[3][b] + lNx[2][a] * lNx[2][b] + lNx[1][a] * lNx[1][b]) * dOmega;
	
	  /* Augmented Lagrangian (penalty) to enhance incompressibility */ 
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2] += penalty * (lNxO[1][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn+2]  += penalty * (lNxO[2][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+2*nn]  += penalty * (lNxO[3][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+2*nn+1]  += penalty * (lNxO[3][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2*nn+2]  += penalty * (lNxO[3][a] * lNxO[3][b]) * dOmegaq;

	  /* And the matrix is symmetric so, for the time being, this is stored */
	  elt_k[qn]       = elt_k[pn];
	  elt_k[qn+nn]    = elt_k[pn+1];
	  elt_k[qn+2*nn]  = elt_k[pn+2];
	  elt_k[qn+1]     = elt_k[pn+nn];
	  elt_k[qn+nn+1]  = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1]= elt_k[pn+nn+2];
	  elt_k[qn+2]     = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]  = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2]= elt_k[pn+2*nn+2];
	}
	/* 3D Cosserat case */
	if(E->control.B_matrix==4) {
	  elt_k[pn]        += 2.0 * lNx[1][a] * lNx[1][b] * dOmega + Gdig * (lNx[2][a] * lNx[2][b]+lNx[3][a] * lNx[3][b]) ;
	  elt_k[pn+1]      += Gnotdig * lNx[2][a] * lNx[1][b] ;
	  elt_k[pn+2]      += Gnotdig * lNx[3][a] * lNx[1][b] ;
	  elt_k[pn+3]      += 0.0 ;
	  elt_k[pn+4]      += (Gdig - Gnotdig) * lNx[3][a] * lN[b] ;
	  elt_k[pn+5]      += (Gdig - Gnotdig) * lNx[2][a] * lN[b] ;
	  elt_k[pn+nn]     += Gnotdig * lNx[1][a] * lNx[2][b] ;
	  elt_k[pn+nn+1]   += 2.0 * lNx[2][a] * lNx[2][b] * dOmega + Gdig *( lNx[1][a] * lNx[1][b]+lNx[3][a] * lNx[3][b]) ;
	  elt_k[pn+nn+2]   += Gnotdig * lNx[3][a] * lNx[2][b] ;
	  elt_k[pn+nn+3]   += (Gdig - Gnotdig) * lNx[3][a] * lN[b] ;
	  elt_k[pn+nn+4]   += 0.0 ;
	  elt_k[pn+nn+5]   += (Gnotdig - Gdig) * lNx[1][a] * lN[b] ;
	  elt_k[pn+2*nn]   += Gnotdig * lNx[3][b] * lNx[1][a] ;
	  elt_k[pn+2*nn+1] += Gnotdig * lNx[3][b] * lNx[2][a] ;
	  elt_k[pn+2*nn+2] += 2.0 * lNx[3][a] * lNx[3][b] * dOmega + Gdig *( lNx[1][a] * lNx[1][b]+lNx[2][a] * lNx[2][b]) ;
	  elt_k[pn+2*nn+3] += (Gnotdig - Gdig) * lNx[2][a] * lN[b] ;
	  elt_k[pn+2*nn+4] += (Gnotdig - Gdig) * lNx[1][a] * lN[b] ;
	  elt_k[pn+2*nn+5] += 0.0 ;
	  elt_k[pn+3*nn]   += 0.0 ;
	  elt_k[pn+3*nn+1] += (Gdig - Gnotdig) * lNx[3][b] * lN[a] ;
	  elt_k[pn+3*nn+2] += (Gnotdig - Gdig) * lNx[2][b] * lN[a] ;
	  elt_k[pn+3*nn+3] += 2.0 * (Gdig - Gnotdig) * lN[b] * lN[a] + B * (lNx[3][a] * lNx[3][b] + lNx[2][a] * lNx[2][b]) ;
	  elt_k[pn+3*nn+4] += 0.0 ;
	  elt_k[pn+3*nn+5] += 0.0 ;
	  elt_k[pn+4*nn]   += (Gdig - Gnotdig) * lNx[3][b] * lN[a] ;
	  elt_k[pn+4*nn+1] += 0.0 ;
	  elt_k[pn+4*nn+2] += (Gnotdig - Gdig) * lNx[1][b] * lN[a] ;
	  elt_k[pn+4*nn+3] += 0.0 ;
	  elt_k[pn+4*nn+4] += 2.0 * (Gdig - Gnotdig) * lN[b] * lN[a] + B * (lNx[1][a] * lNx[1][b] + lNx[3][a] * lNx[3][b]) ;
	  elt_k[pn+4*nn+5] += 0.0 ;
	  elt_k[pn+5*nn]   += (Gdig - Gnotdig) * lNx[2][b] * lN[a] ;
	  elt_k[pn+5*nn+1] += (Gnotdig - Gdig) * lNx[1][b] * lN[a] ;
	  elt_k[pn+5*nn+2] += 0.0 ;
	  elt_k[pn+5*nn+3] += 0.0 ;
	  elt_k[pn+5*nn+4] += 0.0 ;
	  elt_k[pn+5*nn+5] += 2.0 * (Gdig - Gnotdig) * lN[b] * lN[a] + B * (lNx[1][a] * lNx[1][b] + lNx[2][a] * lNx[2][b]) ;

	  /* Augmented Lagrangian (penalty) to enhance incompressibility */
	  elt_k[pn] += penalty * (lNxO[1][a] * lNxO[1][b]) * dOmegaq;
	  elt_k[pn+1] += penalty * (lNxO[1][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2] += penalty * (lNxO[1][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+nn]  += penalty * (lNxO[2][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+nn+1]  += penalty * (lNxO[2][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+nn+2]  += penalty * (lNxO[2][a] * lNxO[3][b]) * dOmegaq;
	  elt_k[pn+2*nn]  += penalty * (lNxO[3][a] * lNxO[1][b]) * dOmegaq;	
	  elt_k[pn+2*nn+1]  += penalty * (lNxO[3][a] * lNxO[2][b]) * dOmegaq;
	  elt_k[pn+2*nn+2]  += penalty * (lNxO[3][a] * lNxO[3][b]) * dOmegaq;

	  elt_k[qn]        = elt_k[pn];
	  elt_k[qn+nn]     = elt_k[pn+1];
	  elt_k[qn+2*nn]   = elt_k[pn+2];
	  elt_k[qn+3*nn]   = elt_k[pn+3];
	  elt_k[qn+4*nn]   = elt_k[pn+4];
	  elt_k[qn+5*nn]   = elt_k[pn+5];
	  elt_k[qn+1]      = elt_k[pn+nn];
	  elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	  elt_k[qn+3*nn+1] = elt_k[pn+nn+3];
	  elt_k[qn+4*nn+1] = elt_k[pn+nn+4];
	  elt_k[qn+5*nn+1] = elt_k[pn+nn+5];
	  elt_k[qn+2]      = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	  elt_k[qn+3*nn+2] = elt_k[pn+2*nn+3];
	  elt_k[qn+4*nn+2] = elt_k[pn+2*nn+4];
	  elt_k[qn+5*nn+2] = elt_k[pn+2*nn+5];
	  elt_k[qn+3]      = elt_k[pn+3*nn];
	  elt_k[qn+nn+3]   = elt_k[pn+3*nn+1];
	  elt_k[qn+2*nn+3] = elt_k[pn+3*nn+2];
	  elt_k[qn+3*nn+3] = elt_k[pn+3*nn+3];
	  elt_k[qn+4*nn+3] = elt_k[pn+3*nn+4];
	  elt_k[qn+5*nn+3] = elt_k[pn+3*nn+5];
	  elt_k[qn+4]      = elt_k[pn+4*nn];
	  elt_k[qn+nn+4]   = elt_k[pn+4*nn+1];
	  elt_k[qn+2*nn+4] = elt_k[pn+4*nn+2];
	  elt_k[qn+3*nn+4] = elt_k[pn+4*nn+3];
	  elt_k[qn+4*nn+4] = elt_k[pn+4*nn+4];
	  elt_k[qn+5*nn+4] = elt_k[pn+4*nn+5];
	  elt_k[qn+5]      = elt_k[pn+5*nn];
	  elt_k[qn+nn+5]   = elt_k[pn+5*nn+1];
	  elt_k[qn+2*nn+5] = elt_k[pn+5*nn+2];
	  elt_k[qn+3*nn+5] = elt_k[pn+5*nn+3];
	  elt_k[qn+4*nn+5] = elt_k[pn+5*nn+4];
	  elt_k[qn+5*nn+5] = elt_k[pn+5*nn+5];
	}
      }
  }

/*  if(el==1)
  for(i=0;i<nn*nn;i++) {
    fprintf(fp,"elt_k[%d] = %g \n",i,elt_k[i]);
    } */

  if(E->control.HAVE_SKEWBCS) {

      /* Build elt_Rot matrix */
   
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;

      /* Add in all nodal rotation matrix terms */
      for(a=1;a<=ends;a++) {/*5*/
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];
	  
	  if(dims==3) {
	    elt_R[px*nn+py] = node_R[0*dims+2];
	    elt_R[pz*nn+py] = node_R[1*dims+2];
	    elt_R[py*nn+px] = node_R[2*dims+0];
	    elt_R[py*nn+pz] = node_R[2*dims+1];
	    elt_R[py*nn+py] = node_R[2*dims+2];
	  }
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
  }
 
  /* fclose(fp) ; */


  return;
}
#endif

void get_elt_k_general(
		       struct All_variables *E,
		       int el,
		       higher_precision *elt_k,
		       int level
		       )
{  
  int a,b,p,q,m,i,j,k,tracers,tr;
  int pn,qn,ad,bd;
  standard_precision lN[ELNMAX+1];
  standard_precision lNx[4][ELNMAX+1];
  standard_precision lNxO[4][ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega,weight,dOmegaq;
  higher_precision del_elt_k[24*24];
  higher_precision elt_R[24*24];
  higher_precision cos_theta,sin_theta;
  higher_precision elt_k1[24*24];
  higher_precision *node_R;
  higher_precision sum;
  int px,pz,py,node;  

  standard_precision B, Gc ;
  standard_precision Gdig, Gnotdig ;

  standard_precision penalty;
  standard_precision Visc;
  standard_precision one_over_x2;
  standard_precision one_over_x;
  
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nn=loc_mat_size[dims][E->control.model];
  const int ends=enodes[dims];

  standard_precision DV[16][16];
  standard_precision DP[16][16];
  standard_precision BaT_DV_Bb[7][7];
  standard_precision BaT_DP_Bb[7][7];
  standard_precision BaT_DV[7][16];
  standard_precision Ba[16][7];
  standard_precision Bb[16][7];

  for(p=0;p<nn*nn;p++)
    elt_k[p] = 0.0;
  penalty = E->control.AUG_lagrangian;  

  get_global_v_x_shape_fn(E,el,lNxO,&dOmega,0.0,0.0,0.0,level);

  Visc=dOmega=0.0;
  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
    if(E->tracer.property_group[m] < 0)
      continue;

    dOmega += E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m];
    Visc += E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
  }
  if(dOmega != 0.0)
    Visc /= dOmega;

  for(i=0;i<E->tracer.tr_in_element_number[level][el];i++) {
    m = E->tracer.tr_in_element[level][i+E->tracer.tr_in_element_offset[level][el]];
  
    eta1 = E->tracer.eta1[level][m];
    eta2 = E->tracer.eta2[level][m];
    if(3==dims) 
      eta3 = E->tracer.eta3[level][m];

    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,level);
  
    if(E->control.AXI || E->control.model==2) {
      v_shape_fn(E,el,lN,eta1,eta2,eta3,level);
      one_over_x = 1.0 / E->tracer.tx[m];  
      one_over_x2 = one_over_x * one_over_x;    
    }

    dOmegaq = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
    dOmega  = E->tracer.tracer_weight_fn[level][m] * E->tracer.tracer_jacobian[level][m] * E->tracer.Visc[m];
  
    if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio >= 0.0)
      penalty = max(0.0,min(E->control.AUG_lagrangian,E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio));
    else 
      penalty = E->control.AUG_lagrangian;

    /* dOmega is modified in axisymmetric case */

    if(E->control.AXI) {
      dOmega *= 2 * M_PI * E->tracer.tx[m];
    }

    for(a=1;a<=ends;a++)
      for(b=a;b<=ends;b++) {
	pn=(dofs*(a-1))*nn+dofs*(b-1);
	qn=(dofs*(b-1))*nn+dofs*(a-1);

	construct_Bi(E,lNx,a,Ba,one_over_x,lN) ;
	construct_D_Vel(E,DV,dOmega,penalty,m,level) ;
	construct_Bi(E,lNx,b,Bb,one_over_x,lN) ;
	
	mul_At_B(Ba,DV,BaT_DV) ;
	mul_A_B(BaT_DV,Bb,BaT_DV_Bb) ;
	
	construct_Bi(E,lNxO,a,Ba,one_over_x,lN) ;
	construct_D_Pres(E,DP,dOmegaq,penalty) ;
	construct_Bi(E,lNxO,b,Bb,one_over_x,lN) ;
	
	mul_At_B(Ba,DP,BaT_DV) ;
	mul_A_B(BaT_DV,Bb,BaT_DP_Bb) ;
	
	switch(dofs) {
	case 2:
	  elt_k[pn]      += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	  elt_k[pn+1]    += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	  elt_k[pn+nn]   += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	  elt_k[pn+nn+1] += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;
	  
	  /* And the matrix is symmetric so, for the time being, this is stored */
	  elt_k[qn]       = elt_k[pn];
	  elt_k[qn+nn]    = elt_k[pn+1];
	  elt_k[qn+1]     = elt_k[pn+nn];
	  elt_k[qn+nn+1]  = elt_k[pn+nn+1];
	  break ;
	case 3:
	  elt_k[pn]        += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	  elt_k[pn+1]      += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	  elt_k[pn+2]      += BaT_DV_Bb[1][3] + BaT_DP_Bb[1][3] ;
	  elt_k[pn+nn]     += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	  elt_k[pn+nn+1]   += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;
	  elt_k[pn+nn+2]   += BaT_DV_Bb[2][3] + BaT_DP_Bb[2][3] ;
	  elt_k[pn+2*nn]   += BaT_DV_Bb[3][1] + BaT_DP_Bb[3][1] ;
	  elt_k[pn+2*nn+1] += BaT_DV_Bb[3][2] + BaT_DP_Bb[3][2] ;
	  elt_k[pn+2*nn+2] += BaT_DV_Bb[3][3] + BaT_DP_Bb[3][3] ;
	  
	  elt_k[qn]        = elt_k[pn];
	  elt_k[qn+nn]     = elt_k[pn+1];
	  elt_k[qn+2*nn]   = elt_k[pn+2];
	  elt_k[qn+1]      = elt_k[pn+nn];
	  elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	  elt_k[qn+2]      = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	  break ;
	case 6:
	  elt_k[pn]        += BaT_DV_Bb[1][1] + BaT_DP_Bb[1][1] ;
	  elt_k[pn+1]      += BaT_DV_Bb[1][2] + BaT_DP_Bb[1][2] ;
	  elt_k[pn+2]      += BaT_DV_Bb[1][3] + BaT_DP_Bb[1][3] ;
	  elt_k[pn+3]      += BaT_DV_Bb[1][4] + BaT_DP_Bb[1][4] ;
	  elt_k[pn+4]      += BaT_DV_Bb[1][5] + BaT_DP_Bb[1][5] ;
	  elt_k[pn+5]      += BaT_DV_Bb[1][6] + BaT_DP_Bb[1][6] ;
	  elt_k[pn+nn]     += BaT_DV_Bb[2][1] + BaT_DP_Bb[2][1] ;
	  elt_k[pn+nn+1]   += BaT_DV_Bb[2][2] + BaT_DP_Bb[2][2] ;
	  elt_k[pn+nn+2]   += BaT_DV_Bb[2][3] + BaT_DP_Bb[2][3] ;
	  elt_k[pn+nn+3]   += BaT_DV_Bb[2][4] + BaT_DP_Bb[2][4] ;
	  elt_k[pn+nn+4]   += BaT_DV_Bb[2][5] + BaT_DP_Bb[2][5] ;
	  elt_k[pn+nn+5]   += BaT_DV_Bb[2][6] + BaT_DP_Bb[2][6] ;
	  elt_k[pn+2*nn]   += BaT_DV_Bb[3][1] + BaT_DP_Bb[3][1] ;
	  elt_k[pn+2*nn+1] += BaT_DV_Bb[3][2] + BaT_DP_Bb[3][2] ;
	  elt_k[pn+2*nn+2] += BaT_DV_Bb[3][3] + BaT_DP_Bb[3][3] ;
	  elt_k[pn+2*nn+3] += BaT_DV_Bb[3][4] + BaT_DP_Bb[3][4] ;
	  elt_k[pn+2*nn+4] += BaT_DV_Bb[3][5] + BaT_DP_Bb[3][5] ;
	  elt_k[pn+2*nn+5] += BaT_DV_Bb[3][6] + BaT_DP_Bb[3][6] ;
	  elt_k[pn+3*nn]   += BaT_DV_Bb[4][1] + BaT_DP_Bb[4][1] ;
	  elt_k[pn+3*nn+1] += BaT_DV_Bb[4][2] + BaT_DP_Bb[4][2] ;
	  elt_k[pn+3*nn+2] += BaT_DV_Bb[4][3] + BaT_DP_Bb[4][3] ;
	  elt_k[pn+3*nn+3] += BaT_DV_Bb[4][4] + BaT_DP_Bb[4][4] ;
	  elt_k[pn+3*nn+4] += BaT_DV_Bb[4][5] + BaT_DP_Bb[4][5] ;
	  elt_k[pn+3*nn+5] += BaT_DV_Bb[4][6] + BaT_DP_Bb[4][6] ;
	  elt_k[pn+4*nn]   += BaT_DV_Bb[5][1] + BaT_DP_Bb[5][1] ;
	  elt_k[pn+4*nn+1] += BaT_DV_Bb[5][2] + BaT_DP_Bb[5][2] ;
	  elt_k[pn+4*nn+2] += BaT_DV_Bb[5][3] + BaT_DP_Bb[5][3] ;
	  elt_k[pn+4*nn+3] += BaT_DV_Bb[5][4] + BaT_DP_Bb[5][4] ;
	  elt_k[pn+4*nn+4] += BaT_DV_Bb[5][5] + BaT_DP_Bb[5][5] ;
	  elt_k[pn+4*nn+5] += BaT_DV_Bb[5][6] + BaT_DP_Bb[5][6] ;
	  elt_k[pn+5*nn]   += BaT_DV_Bb[6][1] + BaT_DP_Bb[6][1] ;
	  elt_k[pn+5*nn+1] += BaT_DV_Bb[6][2] + BaT_DP_Bb[6][2] ;
	  elt_k[pn+5*nn+2] += BaT_DV_Bb[6][3] + BaT_DP_Bb[6][3] ;
	  elt_k[pn+5*nn+3] += BaT_DV_Bb[6][4] + BaT_DP_Bb[6][4] ;
	  elt_k[pn+5*nn+4] += BaT_DV_Bb[6][5] + BaT_DP_Bb[6][5] ;
	  elt_k[pn+5*nn+5] += BaT_DV_Bb[6][6] + BaT_DP_Bb[6][6] ;
	  
	  elt_k[qn]        = elt_k[pn];
	  elt_k[qn+nn]     = elt_k[pn+1];
	  elt_k[qn+2*nn]   = elt_k[pn+2];
	  elt_k[qn+3*nn]   = elt_k[pn+3];
	  elt_k[qn+4*nn]   = elt_k[pn+4];
	  elt_k[qn+5*nn]   = elt_k[pn+5];
	  elt_k[qn+1]      = elt_k[pn+nn];
	  elt_k[qn+nn+1]   = elt_k[pn+nn+1];
	  elt_k[qn+2*nn+1] = elt_k[pn+nn+2];
	  elt_k[qn+3*nn+1] = elt_k[pn+nn+3];
	  elt_k[qn+4*nn+1] = elt_k[pn+nn+4];
	  elt_k[qn+5*nn+1] = elt_k[pn+nn+5];
	  elt_k[qn+2]      = elt_k[pn+2*nn];
	  elt_k[qn+nn+2]   = elt_k[pn+2*nn+1];
	  elt_k[qn+2*nn+2] = elt_k[pn+2*nn+2];
	  elt_k[qn+3*nn+2] = elt_k[pn+2*nn+3];
	  elt_k[qn+3*nn+2] = elt_k[pn+2*nn+4];
	  elt_k[qn+5*nn+2] = elt_k[pn+2*nn+5];
	  elt_k[qn+3]      = elt_k[pn+3*nn];
	  elt_k[qn+nn+3]   = elt_k[pn+3*nn+1];
	  elt_k[qn+2*nn+3] = elt_k[pn+3*nn+2];
	  elt_k[qn+3*nn+3] = elt_k[pn+3*nn+3];
	  elt_k[qn+4*nn+3] = elt_k[pn+3*nn+4];
	  elt_k[qn+5*nn+3] = elt_k[pn+3*nn+5];
	  elt_k[qn+4]      = elt_k[pn+4*nn];
	  elt_k[qn+nn+4]   = elt_k[pn+4*nn+1];
	  elt_k[qn+2*nn+4] = elt_k[pn+4*nn+2];
	  elt_k[qn+3*nn+4] = elt_k[pn+4*nn+3];
	  elt_k[qn+4*nn+4] = elt_k[pn+4*nn+4];
	  elt_k[qn+5*nn+4] = elt_k[pn+4*nn+5];
	  elt_k[qn+5]      = elt_k[pn+5*nn];
	  elt_k[qn+nn+5]   = elt_k[pn+5*nn+1];
	  elt_k[qn+2*nn+5] = elt_k[pn+5*nn+2];
	  elt_k[qn+3*nn+5] = elt_k[pn+5*nn+3];
	  elt_k[qn+3*nn+5] = elt_k[pn+5*nn+4];
	  elt_k[qn+5*nn+5] = elt_k[pn+5*nn+5];
	  break ;
	}
      }
  }
  

  if(E->control.HAVE_SKEWBCS) {
      /* Build elt_Rot matrix */
   
      /* Identity */
      for(p=0;p<nn*nn;p++)
	elt_R[p] = 0.0;

      for(p=0;p<nn;p++)  
	elt_R[p*nn+p] = 1.0;

      /* Add in all nodal rotation matrix terms */
      for(a=1;a<=ends;a++) {/*5*/
	node=E->IEN[level][el].node[a];
	if(E->NODE[level][node] & SKEWBC) {
	  node_R = E->curvilinear.NODE_R[level][node];

	  px = (a-1) * dims;
	  pz = (a-1) * dims + 1;
	  py = (a-1) * dims + 2;

	  elt_R[px*nn+px] = node_R[0*dims+0]; 
         /* note all nodes independent effect on elt_K */
	  elt_R[px*nn+pz] = node_R[0*dims+1];
	  elt_R[pz*nn+px] = node_R[1*dims+0];
	  elt_R[pz*nn+pz] = node_R[1*dims+1];
	  
	  if(dims==3) {
	    elt_R[px*nn+py] = node_R[0*dims+2];
	    elt_R[pz*nn+py] = node_R[1*dims+2];
	    elt_R[py*nn+px] = node_R[2*dims+0];
	    elt_R[py*nn+pz] = node_R[2*dims+1];
	    elt_R[py*nn+py] = node_R[2*dims+2];
	  }
	}
      }

      /* post multiply by Rot */
      for(i=0;i<nn;i++) 
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_k[i*nn+k] * elt_R[k*nn+j];
	    
	  elt_k1[i*nn+j] = sum;
	}

      /* pre multiply by RotT */
      for(i=0;i<nn;i++)  
	for(j=0;j<nn;j++) {
	  sum=0.0;
	  for(k=0;k<nn;k++)
	    sum += elt_R[k*nn+i] * elt_k1[k*nn+j];
	  
	  elt_k[i*nn+j] = sum;
	}
  }
  /* fclose(fp) ; */
  return;
}
