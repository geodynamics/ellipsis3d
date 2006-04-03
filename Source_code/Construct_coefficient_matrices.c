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


#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"


higher_precision *G_eltK[MAX_LEVELS];


void construct_node_maps_3_level(
				 struct All_variables *E,
				 int level
				 )
{
    standard_precision initial_time,CPU_time();

    int el,n,nn,i,j,l;
    int node1,loc1,count,found,element;
    int neq,nno;

    static int been_here = 0;
   
    const int dims=E->mesh.nsd; 
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int max_node = max_node_interaction[dims];

    nno=E->mesh.NNO[level];
        
    for(i=0;i<=(nno+1)*max_node;i++) {
      E->Node_map_3[level][i] = 0;
    }
    
    for(nn=1;nn<=nno;nn++) {
      if(E->NODE[level][nn] & ( OFFSIDE  ))
	continue;

      count=0;
      loc1=nn*max_node;
      
      for(el=1;el<=E->NEI[level].nels[nn];el++)  {
	element = E->NEI[level].element[(nn-1)*ends+el-1]; 
	for(n=1;n<=ends;n++) {
	  node1=E->IEN[level][element].node[n]; /* global node number */
	  if(E->NODE[level][node1] & ( OFFSIDE  ))
	    continue;

	  found=0;
	  for(j=0;j<=count;j++)
	    if(E->Node_map_3[level][loc1+j] == node1) { /* found, index next equation */
	      found++;
	      break;
	    }
	  
	  if(! found) {
	    E->Node_map_3[level][loc1+count] = node1;
	    count++;
	  }
	}
      }
    }
    return;
}


void add_stress_bcs_to_global_F(
			       struct All_variables *E,
			       int level
			       )
{
  int i;
  int w1,w2,w3,w4,w5,w6;

  const int nno=E->mesh.NNO[level];
  const int dofs=E->mesh.dof;

  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  )) {
      continue;
    }
    w1 = (E->NODE[level][i] & SBC1);
    w2 = (E->NODE[level][i] & SBC2);
    w3 = (E->NODE[level][i] & SBC3);
    w4 = (E->NODE[level][i] & SBC4);
    w5 = (E->NODE[level][i] & SBC5);
    w6 = (E->NODE[level][i] & SBC6);
    switch(dofs) {
    case 2:
      if(w1) {
	E->FFb[level][E->ID[level][i].doff[1]] += E->Vb[1][level][i];
      }
      if(w2) {
	E->FFb[level][E->ID[level][i].doff[2]] += E->Vb[2][level][i];
      }
      break ;
    case 3:
      if(w1) {
	E->FFb[level][E->ID[level][i].doff[1]] += E->Vb[1][level][i];
      }
      if(w2) {
	E->FFb[level][E->ID[level][i].doff[2]] += E->Vb[2][level][i];
      }
      if(w3) {
	E->FFb[level][E->ID[level][i].doff[3]] += E->Vb[3][level][i];
      }
      break ;
    case 6:
      if(w1)
	E->FFb[level][E->ID[level][i].doff[1]] += E->Vb[1][level][i];
      if(w2)
	E->FFb[level][E->ID[level][i].doff[2]] += E->Vb[2][level][i];
      if(w3)
	E->FFb[level][E->ID[level][i].doff[3]] += E->Vb[3][level][i];
      if(w4)
	E->FFb[level][E->ID[level][i].doff[4]] += E->Vb[4][level][i];
      if(w5)
	E->FFb[level][E->ID[level][i].doff[5]] += E->Vb[5][level][i];
      if(w6)
	E->FFb[level][E->ID[level][i].doff[6]] += E->Vb[6][level][i];
      break ;
    }
  }
  return;
}


void construct_node_ks_3_level(
			       struct All_variables *E,
			       int level
			       )
{
  int i,j,k,l;
  int node,node1,loc0,found,element,pp,qq;
  int neq,nno,nel;

  /* must change in 3D */
  higher_precision *elt_K;  
  higher_precision *all_elt_K; 

  int w1,w2,w3,ww1,ww2,ww3,w4,w5,w6,ww4,ww5,ww6;
    
  static int been_here[MAX_LEVELS];

  standard_precision initial_time,CPU_time(),elt_k_time,time;
  void get_elt_k();

  standard_precision x1,x2,x3,x4;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int max_node = max_node_interaction[dims];
  const int lms = enodes[dims] * dofs; 
  const int eksize = dofs*ends*dofs*ends;

  initial_time=CPU_time();
  elt_k_time=0.0;

  nel=E->mesh.NEL[level];
  nno=E->mesh.NNO[level];
  neq=E->mesh.NEQ[level];

  /* If there is substantial memory around, store all
     element k's first and then use them - this puts all
     the building in one spot and opens the possibility
     of parallel/vector optimization */

  if(E->control.vector_optimization) {
    if(been_here[level]==0)
      G_eltK[level] = (higher_precision *) Malloc0((nel+1) * eksize * sizeof(higher_precision)); 

    time=CPU_time(); 
    fprintf(stderr,"Getting all element K values in advance for level %d\n",level);
    get_all_elt_k(E,G_eltK[level],level);

    elt_k_time += CPU_time() - time;
  }
  else /* frugal use of memory */ {
    elt_K = (higher_precision *) Malloc0(eksize * sizeof(higher_precision)); 
  }

  if(E->control.verbose)
    fprintf(stderr,"Constructing Node K arrays for level %d\n",level);

  for(i=0;i<neq;i++)
    E->FFb[level][i] = 0.0;
       
  for(i=0;i<=(nno+1)*max_node;i++) {
    E->Node_k11[level][i] = 
      E->Node_k21[level][i] = 
      E->Node_k12[level][i] = 
      E->Node_k22[level][i] = 0.0;
  }
  if(dofs==3) {
    for(i=0;i<=nno*max_node;i++) {
      E->Node_k31[level][i] = 
	E->Node_k32[level][i] = 
	E->Node_k13[level][i] = 
	E->Node_k23[level][i] = 
	E->Node_k33[level][i] = 0.0; 
    }
  }
  else if(dofs==6) {
    for(i=0;i<=nno*max_node;i++) {
      E->Node_k31[level][i] =
	E->Node_k41[level][i] = E->Node_k51[level][i] = E->Node_k61[level][i] = 
	E->Node_k32[level][i] = 
	E->Node_k42[level][i] = E->Node_k52[level][i] = E->Node_k62[level][i] = 
	E->Node_k13[level][i] = E->Node_k23[level][i] = E->Node_k33[level][i] =
	E->Node_k43[level][i] = E->Node_k53[level][i] = E->Node_k63[level][i] =
        E->Node_k14[level][i] = E->Node_k24[level][i] = E->Node_k34[level][i] = 
        E->Node_k44[level][i] = E->Node_k54[level][i] = E->Node_k64[level][i] = 
	E->Node_k15[level][i] = E->Node_k25[level][i] = E->Node_k35[level][i] = 
	E->Node_k45[level][i] = E->Node_k55[level][i] = E->Node_k65[level][i] = 
	E->Node_k16[level][i] = E->Node_k26[level][i] = E->Node_k36[level][i] =
	E->Node_k46[level][i] = E->Node_k56[level][i] = E->Node_k66[level][i] = 0.0; 
      }
  }
    
 
  for(element=1;element<=nel;element++) {
    /* Get or lookup element k matrix */
    time=CPU_time();    
    if(0 && E->control.vector_optimization) {
      elt_K = G_eltK[level] + (element-1) * eksize;
    }
    else /* frugal use of memory */ {
      time=CPU_time();    
      get_elt_k(E,element,elt_K,level);
      elt_k_time += CPU_time() - time;
    }

    add_penalty3(E,element,elt_K,level) ;
    assemble_BC_velocity(E,element,elt_K,level) ;
    assemble_Node_k(E,element,elt_K,level) ;
  }
  E->control.B_is_good[level] = 0;

      
  if(E->control.vector_optimization) {
    /* free((void *) G_eltK); */
  }
  else /* frugal use of memory */ {
    free((void *) elt_K); 
  }

  if(E->control.verbose)
    fprintf(stderr,"Constructed Node K arrays for level %d - elt k build time %g\n",level,elt_k_time);
 
#if 0
  if(level==E->mesh.levmax) {
    for(node=1;node<=E->mesh.NNO[level];node++) {

    if (node != 48 && node != 49 && node != 50)
      continue;

      
      loc0 = node * MAX_NODE_INT_2D;
      /* fprintf(stderr,"node %d\n",node); */
      x4 = x2 = 0.0;
      for(i=0;i<MAX_NODE_INT_2D;i++) {      
	fprintf(stderr,"Interaction with node %03d = %+6e %+6e %+6e\n",E->Node_map_3[level][loc0+i],
		E->Node_k11[level][loc0+i],
		E->Node_k21[level][loc0+i],
		E->Node_k22[level][loc0+i]); 
	
	if(E->Node_map_3[level][loc0+i] != node) {
	  x2 += /*fabs*/(E->Node_k11[level][loc0+i]) + /*fabs*/(E->Node_k21[level][loc0+i]);
	  x4 += /*fabs*/(E->Node_k22[level][loc0+i]) + /*fabs*/(E->Node_k21[level][loc0+i]);
	}
	else {
	  x1 = (E->Node_k11[level][loc0+i]);
	  x3 = (E->Node_k22[level][loc0+i]);
	}
      }
      fprintf(stderr,"1: node %d: Diag = %g, Off diag = %g ... sum = %g\n",node,x1,x2,x1+x2);
      fprintf(stderr,"2: node %d: Diag = %g, Off diag = %g ... sum = %g\n",node,x3,x4,x3+x4);

    }
  }
#endif

  been_here[level] = 1;

  return;
}

void construct_elt_gs_level(
			    struct All_variables *E,
			    int level
			    )
{
  int el;
  void get_elt_g();
  
  for(el=1;el<=E->mesh.NEL[level];el++)
    get_elt_g(E,el,E->elt_del[level][el].g,level);  
  
  return;
}

/*  ==========================================
    construct the lumped mass matrix. The full
    matrix is the FE integration of the density 
    field. The lumped version is the diagonal
    matrix obtained by letting the shape function
    Na be delta(a,b)
    ========================================== */

void mass_matrix(
     struct All_variables *E
)
{ 
  int node,el,a;
  int lv,nint; 

  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];

  for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {  /*1*/
    for(node=1;node<=E->mesh.NNO[lv];node++)
      E->Mass[lv][node] = 0.0;

    /* integral of Na over dOmega */
    for(el=1;el<=E->mesh.NEL[lv];el++) { /*2*/
       get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,lv); 
         
       for(node=1;node<=ends;node++) {
	 for(nint=1;nint<=vpts;nint++) {
      	   E->Mass[lv][E->IEN[lv][el].node[node]] += dOmega.vpt[nint] * E->N.vpt[GNVINDEX(node,nint)] *
	     g_point[nint].weight[dims-1];

	   /* fprintf(stderr,"Mass[%d] += %g * %g * %g\n",E->IEN[lv][el].node[node],dOmega.vpt[nint],E->N.vpt[GNVINDEX(a,nint)],
	      g_point[nint].weight[dims-1]); */
	 }
       }
    }
  }
  /*
  for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
    for(node=1;node<=E->mesh.NNO[lv];node++)
      fprintf(stderr,"%d: Mass[%d] = %g\n",lv,node,E->Mass[lv][node]);
      }*/
  return; 
}

assemble_BC_velocity(
		     struct All_variables *E,
		     int element,
		     higher_precision *elt_K,
		     int level
		     )
{
  int i,j;
  int node,node1,loc0,pp,qq;

  int w1,w2,w3,ww1,ww2,ww3,w4,w5,w6,ww4,ww5,ww6;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int max_node = max_node_interaction[dims];
  const int lms = enodes[dims] * dofs; 

  for(i=1;i<=ends;i++) {  /* i, is the node we are storing to */
    node=E->IEN[level][element].node[i];
    if(E->NODE[level][node] & ( OFFSIDE  ))
      continue;
      
    pp=(i-1)*dofs;
    loc0=node*max_node;

    w1 = !(E->NODE[level][node] & BC1);
    w2 = !(E->NODE[level][node] & BC2);
    w3 = !(E->NODE[level][node] & BC3);
    w4 = !(E->NODE[level][node] & BC4);
    w5 = !(E->NODE[level][node] & BC5);
    w6 = !(E->NODE[level][node] & BC6);
 
    for(j=1;j<=ends;j++) { /* j is the node we are receiving from */
      qq=(j-1)*dofs;
      node1=E->IEN[level][element].node[j];
      if(E->NODE[level][node1] & ( OFFSIDE  ))
	continue;
	
      ww1 = !(E->NODE[level][node1] & BC1);
      ww2 = !(E->NODE[level][node1] & BC2);
      ww3 = !(E->NODE[level][node1] & BC3);
      ww4 = !(E->NODE[level][node1] & BC4);
      ww5 = !(E->NODE[level][node1] & BC5);
      ww6 = !(E->NODE[level][node1] & BC6);


      switch(dofs) {
      case 2:
	if(w1) {  /* No X boundary condition at node */
	  if(!ww1) {
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq] * E->Vb[1][level][node1];
	  }
	  if(!ww2) {
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+1] * E->Vb[2][level][node1];
	  }
	}
	if(w2) {  /* No Z boundary condition at node */
	  if(!ww1) {
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq] * E->Vb[1][level][node1];
	  }
	  if(!ww2) {
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+1] * E->Vb[2][level][node1];
	  }
	}
	break ;
      case 3:
	if(w1) {  /* No X boundary condition at node */
	  if(!ww1) {
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq] * E->Vb[1][level][node1];
	  }
	  if(!ww2) {
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+1] * E->Vb[2][level][node1];
	  }
	  if(!ww3) {
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+2] * E->Vb[3][level][node1];
	  }
	}
	if(w2) {  /* No Z boundary condition at node */
	  if(!ww1) {
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq] * E->Vb[1][level][node1];
	  }
	  if(!ww2) {
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+1] * E->Vb[2][level][node1];
	  }
	  if(!ww3) {
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+2] * E->Vb[3][level][node1];
	  }
	}
	if(w3) {  /* No Y boundary condition at node */
	  if(!ww1) {
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq] * E->Vb[1][level][node1];
	  }
	  if(!ww2) {
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+1] * E->Vb[2][level][node1];
	  }
	  if(!ww3) {
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+2] * E->Vb[3][level][node1];
	  }
	}
	break ;
      case 6:
	if(w1) {  /* No VELX boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[1]] -= elt_K[pp*lms+qq+5] * E->Vb[6][level][node1];
	}
	if(w2) {  /* No VELZ boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[2]] -= elt_K[(pp+1)*lms+qq+5] * E->Vb[6][level][node1];
	}
	if(w3) {  /* No VELY boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[3]] -= elt_K[(pp+2)*lms+qq+5] * E->Vb[6][level][node1];
	}
	if(w4) {  /* No ROTX boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[4]] -= elt_K[(pp+3)*lms+qq+5] * E->Vb[6][level][node1];
	}
	if(w5) {  /* No ROTZ boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[5]] -= elt_K[(pp+4)*lms+qq+5] * E->Vb[6][level][node1];
	}
	if(w6) {  /* No ROTY boundary condition at node */
	  if(!ww1)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq] * E->Vb[1][level][node1];
	  if(!ww2)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq+1] * E->Vb[2][level][node1];
	  if(!ww3)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq+2] * E->Vb[3][level][node1];
	  if(!ww4)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq+3] * E->Vb[4][level][node1];
	  if(!ww5)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq+4] * E->Vb[5][level][node1];
	  if(!ww6)
	    E->FFb[level][E->ID[level][node].doff[6]] -= elt_K[(pp+5)*lms+qq+5] * E->Vb[6][level][node1];
	}
	break ;
      }
    }
  }
  return;
}

assemble_Node_k(
		struct All_variables *E,
		int element,
		higher_precision *elt_K,
		int level
		)
{
  int i,j,k,l;
  int node,node1,loc0,found,pp,qq;
  int neq,nno,nel;

  int w1,w2,w3,ww1,ww2,ww3,w4,w5,w6,ww4,ww5,ww6;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int max_node = max_node_interaction[dims];
  const int lms = enodes[dims] * dofs; 

  for(i=1;i<=ends;i++) {  /* i, is the node we are storing to */
    node=E->IEN[level][element].node[i];
    if(E->NODE[level][node] & ( OFFSIDE  )) {
      continue;
    }
      
    pp=(i-1)*dofs;
    loc0=node*max_node;
      
    w1 = !(E->NODE[level][node] & BC1);
    w2 = !(E->NODE[level][node] & BC2);
    w3 = !(E->NODE[level][node] & BC3);
    w4 = !(E->NODE[level][node] & BC4);
    w5 = !(E->NODE[level][node] & BC5);
    w6 = !(E->NODE[level][node] & BC6);
 
    for(j=1;j<=ends;j++) { /* j is the node we are receiving from */
      qq=(j-1)*dofs;
      node1=E->IEN[level][element].node[j];
      if(E->NODE[level][node1] & ( OFFSIDE  ))
	continue;
	
      ww1 = !(E->NODE[level][node1] & BC1);
      ww2 = !(E->NODE[level][node1] & BC2);
      ww3 = !(E->NODE[level][node1] & BC3);
      ww4 = !(E->NODE[level][node1] & BC4);
      ww5 = !(E->NODE[level][node1] & BC5);
      ww6 = !(E->NODE[level][node1] & BC6);

      found=0;
      for(k=0;k<max_node;k++) {  /* This is quite time consuming */
	if(E->Node_map_3[level][loc0+k] == node1) {
	  found = loc0 + k;
	  break;
	}
      }

      switch(dofs) {
      case 2:
	if(w1) { /* Velocity in direction 1 */
	  if(ww1)
	    E->Node_k11[level][found] += elt_K[pp*lms+qq];
	  if(ww2)
	    E->Node_k12[level][found] += elt_K[pp*lms+qq+1];
	}
	if(w2) { /* Velocity in direction 2 */
	  if(ww1)
	    E->Node_k21[level][found] += elt_K[(pp+1)*lms+qq];
	  if(ww2)
	    E->Node_k22[level][found] += elt_K[(pp+1)*lms+qq+1];
	}
	break ;
      case 3:
	if(w1) { /* Velocity in direction 1 */
	  if(ww1) 
	    E->Node_k11[level][found] += elt_K[pp*lms+qq];
	  if(ww2) 
	    E->Node_k12[level][found] += elt_K[pp*lms+qq+1];
	  if(ww3) 
	    E->Node_k13[level][found] += elt_K[pp*lms+qq+2];
	}
	if(w2) { /* Velocity in direction 2 */
	  if(ww1)
	    E->Node_k21[level][found] += elt_K[(pp+1)*lms+qq];
	  if(ww2) 
	    E->Node_k22[level][found] += elt_K[(pp+1)*lms+qq+1];
	  if(ww3)
	    E->Node_k23[level][found] += elt_K[(pp+1)*lms+qq+2];
	}
	if(w3) { /* Velocity or rotation in direction 3 */
	  if(ww1)
	    E->Node_k31[level][found] += elt_K[(pp+2)*lms+qq];
	  if(ww2)
	    E->Node_k32[level][found] += elt_K[(pp+2)*lms+qq+1];
	  if(ww3)
	    E->Node_k33[level][found] += elt_K[(pp+2)*lms+qq+2];
	}
	break ;
      case 6:
	if(w1) { /* Velocity in direction 1 */
	  if(ww1)
	    E->Node_k11[level][found] += elt_K[pp*lms+qq];       
	  if(ww2) 
	    E->Node_k12[level][found] += elt_K[pp*lms+qq+1];
	  if(ww3) 
	    E->Node_k13[level][found] += elt_K[pp*lms+qq+2];
	  if(ww4) 
	    E->Node_k14[level][found] += elt_K[pp*lms+qq+3];
	  if(ww5) 
	    E->Node_k15[level][found] += elt_K[pp*lms+qq+4];
	  if(ww6) 
	    E->Node_k16[level][found] += elt_K[pp*lms+qq+5];
	}
	if(w2) {    /* Velocity in direction 2 */
	  if(ww1)
	    E->Node_k21[level][found] += elt_K[(pp+1)*lms+qq];       
	  if(ww2) 
	    E->Node_k22[level][found] += elt_K[(pp+1)*lms+qq+1];
	  if(ww3) 
	    E->Node_k23[level][found] += elt_K[(pp+1)*lms+qq+2];
	  if(ww4) 
	    E->Node_k24[level][found] += elt_K[(pp+1)*lms+qq+3];
	  if(ww5) 
	    E->Node_k25[level][found] += elt_K[(pp+1)*lms+qq+4];
	  if(ww6) 
	    E->Node_k26[level][found] += elt_K[(pp+1)*lms+qq+5];
	}
	if(w3) { /* Velocity in direction 3 */
	  if(ww1)
	    E->Node_k31[level][found] += elt_K[(pp+2)*lms+qq];       
	  if(ww2) 
	    E->Node_k32[level][found] += elt_K[(pp+2)*lms+qq+1];
	  if(ww3) 
	    E->Node_k33[level][found] += elt_K[(pp+2)*lms+qq+2];
	  if(ww4) 
	    E->Node_k34[level][found] += elt_K[(pp+2)*lms+qq+3];
	  if(ww5) 
	    E->Node_k35[level][found] += elt_K[(pp+2)*lms+qq+4];
	  if(ww6) 
	    E->Node_k36[level][found] += elt_K[(pp+2)*lms+qq+5];
	}
	if(w4) { /* Rotation in direction 1 */
	  if(ww1)
	    E->Node_k41[level][found] += elt_K[(pp+3)*lms+qq];       
	  if(ww2) 
	    E->Node_k42[level][found] += elt_K[(pp+3)*lms+qq+1];
	  if(ww3) 
	    E->Node_k43[level][found] += elt_K[(pp+3)*lms+qq+2];
	  if(ww4) 
	    E->Node_k44[level][found] += elt_K[(pp+3)*lms+qq+3];
	  if(ww5) 
	    E->Node_k45[level][found] += elt_K[(pp+3)*lms+qq+4];
	  if(ww6) 
	    E->Node_k46[level][found] += elt_K[(pp+3)*lms+qq+5];
	}
	if(w5) { /* Rotation in direction 2 */
	  if(ww1)
	    E->Node_k51[level][found] += elt_K[(pp+4)*lms+qq];       
	  if(ww2) 
	    E->Node_k52[level][found] += elt_K[(pp+4)*lms+qq+1];
	  if(ww3) 
	    E->Node_k53[level][found] += elt_K[(pp+4)*lms+qq+2];
	  if(ww4) 
	    E->Node_k54[level][found] += elt_K[(pp+4)*lms+qq+3];
	  if(ww5) 
	    E->Node_k55[level][found] += elt_K[(pp+4)*lms+qq+4];
	  if(ww6) 
	    E->Node_k56[level][found] += elt_K[(pp+4)*lms+qq+5];
	}
	if(w6) { /* Rotation in direction 3 */
	  if(ww1)
	    E->Node_k61[level][found] += elt_K[(pp+5)*lms+qq];       
	  if(ww2) 
	    E->Node_k62[level][found] += elt_K[(pp+5)*lms+qq+1];
	  if(ww3) 
	    E->Node_k63[level][found] += elt_K[(pp+5)*lms+qq+2];
	  if(ww4) 
	    E->Node_k64[level][found] += elt_K[(pp+5)*lms+qq+3];
	  if(ww5) 
	    E->Node_k65[level][found] += elt_K[(pp+5)*lms+qq+4];
	  if(ww6) 
	    E->Node_k66[level][found] += elt_K[(pp+5)*lms+qq+5];
	}
	break ;
      }
    }
  }



  return;
}
