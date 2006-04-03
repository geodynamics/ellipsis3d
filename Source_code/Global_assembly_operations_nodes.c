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



/* Functions to assemble the element k matrices and the element f vector.
   Note that for the regular grid case the calculation of k becomes repetitive 
   to the point of redundancy. */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <sys/time.h>
#include <sys/resource.h>

/* ================================================================
   Function to assemble the global  F vector.
   ================================================================ */

void assemble_forces(
		     struct All_variables *E,
		     int k,
		     int level
		     )
{
  int a,n,i;
  void add_tracers_to_global_f();

  const int neq=E->mesh.NEQ[level];
  const int nno=E->mesh.NNO[level];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  
  for(a=0;a<neq;a++) {
    E->FF[level][a] = E->FFb[level][a];
  }
 
  /* This should only be computed if tracer properties have
     changed or tracers have moved.

     This may happen if the tracer locations are corrected
     during the velocity iteration */

  add_tracers_to_global_f(E,E->FF[level],k,level);

  return;
}

/* ======================================================
   Assemble Au using stored element K information
   ====================================================== */

void e_assemble_del2_u_6(
     struct All_variables *E,
     standard_precision *V1,
     standard_precision *V2,
     standard_precision *V3,
     standard_precision *V4,
     standard_precision *V5,
     standard_precision *V6,
     standard_precision *Au1,
     standard_precision *Au2,
     standard_precision *Au3,
     standard_precision *Au4,
     standard_precision *Au5,
     standard_precision *Au6,
     int level,
     int strip_bcs
)
{   
  extern higher_precision *G_eltK[MAX_LEVELS];
  

  int el,i,j,k,n;

  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int elts=E->mesh.NEL[level];
  const int neq=E->mesh.NEQ[level];
  const int nno=E->mesh.NNO[level];
  const int eksize = dofs*ends*dofs*ends;
  const int ekusize = dofs*ends;

  standard_precision *all_elt_Au;
  standard_precision *U;
  standard_precision *Elt_K;

  unsigned int *NODE1;

  struct ID *ID1;
  struct IEN *IEN1;


  ID1 = E->ID[level];
  IEN1 = E->IEN[level];
  Elt_K = G_eltK[level];
  NODE1 = E->NODE[level];


  /* We can create the individual element contributions
     independently of one another */

  all_elt_Au = (standard_precision *) Malloc0((ekusize*elts+1) * sizeof(standard_precision));
  U = (standard_precision *) Malloc0((ekusize*elts+1) * sizeof(standard_precision));


  if(2==dofs)
    for(i=1;i<=ends;i++) {
#pragma loop novrec IEN1,U,V1,V2
      for(el=1;el<=elts;el++) {
	U[(el-1)*ekusize+(i-1)*dofs+0] = V1[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+1] = V2[IEN1[el].node[i]];
      }
    }
  else if (3==dofs)
    for(i=1;i<=ends;i++) {
#pragma loop novrec IEN1,U,V1,V2,V3
      for(el=1;el<=elts;el++) {
    	U[(el-1)*ekusize+(i-1)*dofs+0] = V1[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+1] = V2[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+2] = V3[IEN1[el].node[i]];
      }
    }
  else /* 6 */
    for(i=1;i<=ends;i++) {
#pragma loop novrec IEN1,U,V1,V2,V3,V4,V5,V6
      for(el=1;el<=elts;el++) {
	U[(el-1)*ekusize+(i-1)*dofs+0] = V1[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+1] = V2[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+2] = V3[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+3] = V4[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+4] = V5[IEN1[el].node[i]];
	U[(el-1)*ekusize+(i-1)*dofs+5] = V6[IEN1[el].node[i]];
      }
    }

  for(i=0;i<=ekusize*elts;i++)
    all_elt_Au[i] = 0.0;
  


  /* Time is absorbed in this loop ! */
  for(i=0;i<ekusize;i++) {
    for(j=0;j<ekusize;j++) {
#pragma loop novrec Elt_K,all_elt_Au,U 
      for(el=1;el<=elts;el++) {
	all_elt_Au[(el-1)*ekusize+i] += Elt_K[(el-1)*eksize+i*ekusize+j] * U[(el-1)*ekusize+j];
      }
    }
  }


  /* and assemble this back down to the usual places */



  if(2==dofs) {
#pragma loop novrec Au1,Au2
    for(i=1;i<=nno;i++) {
      Au1[i] = 0.0;
      Au2[i] = 0.0;
    }
    for(i=1;i<=ends;i++) {
#pragma loop novrec Au1,Au2,all_elt_Au,NODE1 
      for(el=1;el<=elts;el++) {
	n=IEN1[el].node[i];
	Au1[n] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+0] * ((NODE1[n] & BC1) == 0);
	Au2[n] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+1] * ((NODE1[n] & BC2) == 0);
      }
    }
  } 
  else if (3==dofs) {
    for(i=1;i<=nno;i++) {
      Au1[i] = 0.0;
      Au2[i] = 0.0;
      Au3[i] = 0.0;
    }
     
     for(el=1;el<=elts;el++) {
      for(i=1;i<=ends;i++) {
	Au1[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+0];
	Au2[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+1];
	Au3[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+2];
   }
    }
  }
  else /* 6 */ {
    for(i=1;i<=nno;i++) {
      Au1[i] = 0.0;
      Au2[i] = 0.0;
      Au3[i] = 0.0;
      Au4[i] = 0.0;
      Au5[i] = 0.0;
      Au6[i] = 0.0;
    }
     
    for(el=1;el<=elts;el++) {
      for(i=1;i<=ends;i++) {
	Au1[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+0];
	Au2[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+1];
	Au3[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+2];
	Au4[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+3];
	Au5[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+4];
	Au6[IEN1[el].node[i]] += all_elt_Au[(el-1)*ekusize+(i-1)*dofs+5];
      }
    }
  }

  /* Strip boundary conditions ? */

  
  if(strip_bcs)
    strip_bcs_from_residual_6(E,Au1,Au2,Au3,Au4,Au5,Au6,level);

   
  free((void *) all_elt_Au);
  free((void *) U);


}

/* ======================================================
   Assemble Au using stored, nodal coefficients.
   ====================================================== */

void n_assemble_del2_u_6(
			 struct All_variables *E,
			 standard_precision *V1,
			 standard_precision *V2,
			 standard_precision *V3,
			 standard_precision *V4,
			 standard_precision *V5,
			 standard_precision *V6,
			 standard_precision *Au1,
			 standard_precision *Au2,
			 standard_precision *Au3,
			 standard_precision *Au4,
			 standard_precision *Au5,
			 standard_precision *Au6,
			 int level,
			 int strip_bcs
			 )
{
  int  e,i;
  int *C;
  
  higher_precision *B11,*B12,*B13,*B14,*B15,*B16;
  higher_precision *B21,*B22,*B23,*B24,*B25,*B26;
  higher_precision *B31,*B32,*B33,*B34,*B35,*B36;
  higher_precision *B41,*B42,*B43,*B44,*B45,*B46;
  higher_precision *B51,*B52,*B53,*B54,*B55,*B56;
  higher_precision *B61,*B62,*B63,*B64,*B65,*B66;
  
  const int neq=E->mesh.NEQ[level]; 
  const int nno=E->mesh.NNO[level];
  const int dims=E->mesh.nsd; 
  const int dofs=E->mesh.dof;
  
  if(2==dofs) {
    for(e=0;e<=nno;e++) 
      Au1[e]=Au2[e]=0.0;
    
    for(e=1;e<=nno;e++) {
      if(!(E->NODE[level][e] & ( OFFSIDE  ))) { 	  
	
	C=E->Node_map_3[level]+e*MAX_NODE_INT_2D;
	B11=E->Node_k11[level]+e*MAX_NODE_INT_2D;
	B12=E->Node_k12[level]+e*MAX_NODE_INT_2D;
	B21=E->Node_k21[level]+e*MAX_NODE_INT_2D;
	B22=E->Node_k22[level]+e*MAX_NODE_INT_2D;
	
	for(i=0;i<MAX_NODE_INT_2D;i++) {
	  Au1[C[i]] += B11[i] * V1[e] + B21[i] * V2[e];
	  Au2[C[i]] += B12[i] * V1[e] + B22[i] * V2[e]; 
	}
      }
    }
  }
  else if(3==dofs && 2==dims) /* Cosserat 2D */ {
    for(e=0;e<=nno;e++) 
      Au1[e]=Au2[e]=Au3[e]=0.0;
    
    for(e=1;e<=nno;e++) {
      if(!(E->NODE[level][e] & ( OFFSIDE  ))) {
	
	C=E->Node_map_3[level]+e*MAX_NODE_INT_2D;
	B11=E->Node_k11[level]+e*MAX_NODE_INT_2D;
	B12=E->Node_k12[level]+e*MAX_NODE_INT_2D;
	B13=E->Node_k13[level]+e*MAX_NODE_INT_2D;
	B21=E->Node_k21[level]+e*MAX_NODE_INT_2D;
	B22=E->Node_k22[level]+e*MAX_NODE_INT_2D;
	B23=E->Node_k23[level]+e*MAX_NODE_INT_2D;
	B31=E->Node_k31[level]+e*MAX_NODE_INT_2D;
	B32=E->Node_k32[level]+e*MAX_NODE_INT_2D;
	B33=E->Node_k33[level]+e*MAX_NODE_INT_2D;
	
	for(i=0;i<MAX_NODE_INT_2D;i++) {
	  Au1[C[i]] += B11[i] * V1[e] + B21[i] * V2[e] + B31[i] * V3[e];
	  Au2[C[i]] += B12[i] * V1[e] + B22[i] * V2[e] + B32[i] * V3[e];
	  Au3[C[i]] += B13[i] * V1[e] + B23[i] * V2[e] + B33[i] * V3[e];
	}
      }
    }
  }
  else if(3==dofs && 3==dims) /* Classical 3D */ {
    for(e=0;e<=nno;e++) 
      Au1[e]=Au2[e]=Au3[e]=0.0;
    
    for(e=1;e<=nno;e++)     {
      if(!(E->NODE[level][e] & ( OFFSIDE  ))) {
	
	C=E->Node_map_3[level]+e*MAX_NODE_INT_3D;
	B11=E->Node_k11[level]+e*MAX_NODE_INT_3D;
	B12=E->Node_k12[level]+e*MAX_NODE_INT_3D;
	B13=E->Node_k13[level]+e*MAX_NODE_INT_3D;
	B21=E->Node_k21[level]+e*MAX_NODE_INT_3D;
	B22=E->Node_k22[level]+e*MAX_NODE_INT_3D;
	B23=E->Node_k23[level]+e*MAX_NODE_INT_3D;
	B31=E->Node_k31[level]+e*MAX_NODE_INT_3D;
	B32=E->Node_k32[level]+e*MAX_NODE_INT_3D;
	B33=E->Node_k33[level]+e*MAX_NODE_INT_3D;
	
	for(i=0;i<MAX_NODE_INT_3D;i++) {
	  Au1[C[i]] += B11[i] * V1[e] + B21[i] * V2[e] + B31[i] * V3[e];
	  Au2[C[i]] += B12[i] * V1[e] + B22[i] * V2[e] + B32[i] * V3[e]; 
	  Au3[C[i]] += B13[i] * V1[e] + B23[i] * V2[e] + B33[i] * V3[e]; 
	}
      }
    }
  }  
  else /* if(6==dofs) */ { /* Cosserat 3D */
    for(e=0;e<=nno;e++) 
      Au1[e]=Au2[e]=Au3[e]=Au4[e]=Au5[e]=Au6[e]=0.0;
    
    for(e=1;e<=nno;e++)     {
      if(!(E->NODE[level][e] & ( OFFSIDE  ))) {
	
	C=E->Node_map_3[level]+e*MAX_NODE_INT_3D;
	B11=E->Node_k11[level]+e*MAX_NODE_INT_3D;
	B12=E->Node_k12[level]+e*MAX_NODE_INT_3D;
	B13=E->Node_k13[level]+e*MAX_NODE_INT_3D;
	B14=E->Node_k14[level]+e*MAX_NODE_INT_3D;
	B15=E->Node_k15[level]+e*MAX_NODE_INT_3D;
	B16=E->Node_k16[level]+e*MAX_NODE_INT_3D;
	
	B21=E->Node_k21[level]+e*MAX_NODE_INT_3D;
	B22=E->Node_k22[level]+e*MAX_NODE_INT_3D;
	B23=E->Node_k23[level]+e*MAX_NODE_INT_3D;
	B24=E->Node_k24[level]+e*MAX_NODE_INT_3D;
	B25=E->Node_k25[level]+e*MAX_NODE_INT_3D;
	B26=E->Node_k26[level]+e*MAX_NODE_INT_3D;
	
	B31=E->Node_k31[level]+e*MAX_NODE_INT_3D;
	B32=E->Node_k32[level]+e*MAX_NODE_INT_3D;
	B33=E->Node_k33[level]+e*MAX_NODE_INT_3D;
	B34=E->Node_k34[level]+e*MAX_NODE_INT_3D;
	B35=E->Node_k35[level]+e*MAX_NODE_INT_3D;
	B36=E->Node_k36[level]+e*MAX_NODE_INT_3D;
	
	B41=E->Node_k41[level]+e*MAX_NODE_INT_3D;
	B42=E->Node_k42[level]+e*MAX_NODE_INT_3D;
	B43=E->Node_k43[level]+e*MAX_NODE_INT_3D;
	B44=E->Node_k44[level]+e*MAX_NODE_INT_3D;
	B45=E->Node_k45[level]+e*MAX_NODE_INT_3D;
	B46=E->Node_k46[level]+e*MAX_NODE_INT_3D;
	
	B51=E->Node_k51[level]+e*MAX_NODE_INT_3D;
	B52=E->Node_k52[level]+e*MAX_NODE_INT_3D;
	B53=E->Node_k53[level]+e*MAX_NODE_INT_3D;
	B54=E->Node_k54[level]+e*MAX_NODE_INT_3D;
	B55=E->Node_k55[level]+e*MAX_NODE_INT_3D;
	B56=E->Node_k56[level]+e*MAX_NODE_INT_3D;
	
	B61=E->Node_k61[level]+e*MAX_NODE_INT_3D;
	B62=E->Node_k62[level]+e*MAX_NODE_INT_3D;
	B63=E->Node_k63[level]+e*MAX_NODE_INT_3D;
	B64=E->Node_k64[level]+e*MAX_NODE_INT_3D;
	B65=E->Node_k65[level]+e*MAX_NODE_INT_3D;
	B66=E->Node_k66[level]+e*MAX_NODE_INT_3D;
	
	for(i=0;i<MAX_NODE_INT_3D;i++) {
	  Au1[C[i]] += B11[i] * V1[e] + B21[i] * V2[e] + B31[i] * V3[e] + B41[i] * V4[e] + B51[i] * V5[e] + B61[i] * V6[e];
	  Au2[C[i]] += B12[i] * V1[e] + B22[i] * V2[e] + B32[i] * V3[e] + B42[i] * V4[e] + B52[i] * V5[e] + B62[i] * V6[e]; 
	  Au3[C[i]] += B13[i] * V1[e] + B23[i] * V2[e] + B33[i] * V3[e] + B43[i] * V4[e] + B53[i] * V5[e] + B63[i] * V6[e]; 
	  Au4[C[i]] += B14[i] * V1[e] + B24[i] * V2[e] + B34[i] * V3[e] + B44[i] * V4[e] + B54[i] * V5[e] + B64[i] * V6[e];
	  Au5[C[i]] += B15[i] * V1[e] + B25[i] * V2[e] + B35[i] * V3[e] + B45[i] * V4[e] + B55[i] * V5[e] + B65[i] * V6[e]; 
	  Au6[C[i]] += B16[i] * V1[e] + B26[i] * V2[e] + B36[i] * V3[e] + B46[i] * V4[e] + B56[i] * V5[e] + B66[i] * V6[e]; 
	}
      }
    }
  }
  return; 
}
#if 0
/* ==========================================
   Assemble a div_u vector element by element
   ==========================================  */


void assemble_div_us(
     struct All_variables *E,
     standard_precision *U,
     standard_precision *divU,
     int level
)
{ 
  int e,i,j,k,p,a,b,node;
  higher_precision divmax;
  higher_precision elt_g[24][1];
    
  const int nel=E->mesh.NEL[level];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int npno=E->mesh.NPNO[level];
    
  divmax = 0.0;
  for(i=1;i<=npno;i++)
    divU[i] = 0.0;
   
  for(e=1;e<=nel;e++) {
    for(a=1;a<=ends;a++)   {
      node=E->IEN[level][e].node[a];
      for(i=1;i<=dims;i++)  {
	p = (a-1)*dims +i-1;
	j = E->ID[level][node].doff[i];
	/* for(b=0;b<ploc_mat_size[E->mesh.nsd];b++) */
	divU[e] += E->elt_del[level][e].g[p][0] * U[j];
      }
    }
  }
  return;
}
#endif
/* ==========================================
   Assemble a div_u vector element by element
   ==========================================  */

void assemble_div_us3(
		      struct All_variables *E,
		      standard_precision *U1,
		      standard_precision *U2,
		      standard_precision *U3,
		      standard_precision *divU,
		      int level
		      )
{ 
  int e,i,j,k,p,a,b,node;
  higher_precision divmax;
  higher_precision elt_g[24][1];
  
  const int nel=E->mesh.NEL[level];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd; 
  const int dofs=E->mesh.dof;
  const int npno=E->mesh.NPNO[level];
  
  for(i=1;i<=npno;i++) 
    divU[i] = 0.0;
  
  for(e=1;e<=nel;e++) {
    /* if(E->tracer.tr_in_element_number[level][e] == 0)
       continue; */   
    for(a=1;a<=ends;a++) {
      node=E->IEN[level][e].node[a];

      divU[e] += E->elt_del[level][e].g[(a-1)*dims][0] * U1[node];
      divU[e] += E->elt_del[level][e].g[(a-1)*dims+1][0] * U2[node];
      if(3==dims)
	divU[e] += E->elt_del[level][e].g[(a-1)*dims+2][0] * U3[node];
    } 
  }
  return;
}

/* ==========================================
   Assemble a grad_P vector element by element
   ==========================================  */

void assemble_grad_qs3(
		       struct All_variables *E,
		       standard_precision *Q,
		       standard_precision *gradQ1,
		       standard_precision *gradQ2,
		       standard_precision *gradQ3,
		       int level
		       )
{
  int e,i,j,p,a,node;
 
  const int nel=E->mesh.NEL[level];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;

  if(3==dims/* || dofs > 2*/)
    for(i=1;i<=E->mesh.NNO[level];i++)
      gradQ1[i] = gradQ2[i] = gradQ3[i] = 0.0;
  else
    for(i=1;i<=E->mesh.NNO[level];i++)
      gradQ1[i] = gradQ2[i] = 0.0;
 
  for(e=1;e<=nel;e++) {
    /* if(E->tracer.tr_in_element_number[level][e] < 1 || 0.0==P[e])
       continue;  */
    for(a=1;a<=ends;a++) {
      node=E->IEN[level][e].node[a];	
      gradQ1[node] += E->elt_del[level][e].g[(a-1)*dims][0] * Q[e];
      gradQ2[node] += E->elt_del[level][e].g[(a-1)*dims + 1][0] * Q[e];
      if(3==dims)
	gradQ3[node] += E->elt_del[level][e].g[(a-1)*dims + 2][0] * Q[e];
    }
  }
  strip_bcs_from_residual_6(E,gradQ1,gradQ2,gradQ3,NULL,NULL,NULL,level);
  return;
}
