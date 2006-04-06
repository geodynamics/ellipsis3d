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


/* Gauss seidel routine when R,V are separated in dimensions, and no Au is required */

void gauss_seidel6_nbn(
		       struct All_variables *E,
		       standard_precision *U1,
			standard_precision *U2,
		       standard_precision *U3,
		       standard_precision *U4,
		       standard_precision *U5,
		       standard_precision *U6,
		       standard_precision *F1,
		       standard_precision *F2,
		       standard_precision *F3,
		       standard_precision *F4,
		       standard_precision *F5,
		       standard_precision *F6,
		       standard_precision *au1,
		       standard_precision *au2,
		       standard_precision *au3,
		       standard_precision *au4,
		       standard_precision *au5,
		       standard_precision *au6,
		       standard_precision acc,
		       standard_precision imp,
		       int cycles,
		       int level,
		       int relax
		       )
{     
  int count,i,j,k,l,m,ns,steps,n1,n2,n3,num_act_nodes;
  int *C;
  int count2;
    
  static int been_here = 0;

  standard_precision u1,u2,u3,u4,u5,u6;
  standard_precision du1,du2,du3,du4,du5,du6;
  standard_precision f1,f2,f3,f4,f5,f6;

  standard_precision residual,residual1,delta_res;
  standard_precision highest_B;

  higher_precision *B11,*B12,*B13,*B14,*B15,*B16;
  higher_precision *B21,*B22,*B23,*B24,*B25,*B26;
  higher_precision *B31,*B32,*B33,*B34,*B35,*B36;
  higher_precision *B41,*B42,*B43,*B44,*B45,*B46;
  higher_precision *B51,*B52,*B53,*B54,*B55,*B56;
  higher_precision *B61,*B62,*B63,*B64,*B65,*B66;
    
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int n=ends * dofs;
  const int neq=E->mesh.NEQ[level];
  const int num_nodes=E->mesh.NNO[level];
  const int max_eqn=max_eqn_interaction[dims];
  const int nox=E->mesh.NOX[level];
  const int noz=E->mesh.NOZ[level];
  const int noy=E->mesh.NOY[level];
  const int nno=E->mesh.NNO[level];

  const int inner_cycles=10;

  if(been_here++ == 0 ) {
    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
      E->monitor.total_gs_cycles[i] = 0;
    }
  }

  E->monitor.total_gs_cycles[level]++;

  if(dofs==2) {
    for(i=1;i<=E->mesh.NNO[level];i++) {
      au1[i]=au2[i]=0.0;
      F1[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC1)) == 0);
      F2[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC2)) == 0);
    }
  }
  else if(dofs==3) {
    for(i=1;i<=E->mesh.NNO[level];i++) {
      au1[i]=au2[i]=au3[i]=0.0;
      F1[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC1)) == 0);
      F2[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC2)) == 0);
      F3[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC3)) == 0);
   }
  }
  else /* if(dofs==6) */ {
    for(i=1;i<=E->mesh.NNO[level];i++) {
      au1[i]=au2[i]=au3[i]=au4[i]=au5[i]=au6[i]=0.0;
      F1[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC1)) == 0);
      F2[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC2)) == 0);
      F3[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC3)) == 0);
      F4[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC4)) == 0);
      F5[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC5)) == 0);
      F6[i] *= ((E->NODE[level][i] & ( OFFSIDE  | BC6)) == 0);
   }
  }

  residual1=0.0;
  num_act_nodes=0;

  for(i=1;i<=num_nodes;i++) {
    if(E->NODE[level][i] & ( OFFSIDE  )) 
      continue;
    residual1 += F1[i] * F1[i] + F2[i] * F2[i];
    if(dofs==3) {
      residual1 +=  F3[i] * F3[i];
    }
    else if(dofs==6) {
      residual1 += F3[i] * F3[i] + F4[i] * F4[i] + F5[i] * F5[i] + F6[i] * F6[i];
    }
    num_act_nodes++;
  }

  delta_res = residual = residual1;

  highest_B = E->BI1[level][E->index[level][1]];
 
    U1[0]=U2[0]=0.0; 
    au1[0]=au2[0]=0.0;
    if(3==dofs) {
      U3[0]=0.0;
      au3[0]=0.0;
    }
    else if(6==dofs) {
      U3[0]=U4[0]=U5[0]=U6[0]=0.0; 
      au3[0]=au4[0]=au5[0]=au6[0]=0.0;
    }

  count2 = 1;
  while(((count2++ == 1)  
	 ||  ((residual1 * imp * imp < residual) && (residual > acc * acc * /*num_act_nodes*/num_nodes) && 
	      ((delta_res > acc * acc * imp * imp)) && (residual < residual1 * 4.0) )) 
	&& (count2 <= cycles)) {

    count = 0;

    for(i=1;i<=E->mesh.NNO[level];i++) {
      au1[i]=au2[i]=0.0;
      if(3==dofs) {
	au3[i]=0.0;
      }
      else if(6==dofs) {
	au3[i]=au4[i]=au5[i]=au6[i]=0.0;
      }
    }

    delta_res = residual;

    while (count++ < inner_cycles + count2 / 10) {

      for(n3=1;n3<=noy;n3++)
	for(n2=1;n2<=noz;n2++)
	  for(n1=1;n1<=nox;n1++)   {

      /* for(n1=1;n1<=nno;n1++) {
	 if(relax==1)
	    i=E->index[level][n1]; 
	 else if (relax==2)
	   i=E->index[level][nno-n1+1]; 
	 else 
	 i=n1; */

	  /*RAA: 21/08/202, fix 3 for extension with T bug, as per LM e-mail */
	  /*i=n2 + (n1-1)*noz + (n3-1)*noz*nox; */
	  if(count%2)
	    i=n2 + (n1-1)*noz + (n3-1)*noz*nox; 
	  else
	    i=nno - (n2 + (n1-1)*noz + (n3-1)*noz*nox); 
	  /*RAA: 21/08/02, end fix; so now alternately traverses back and forth
	      during the relaxation, instead of always in one direction */

	  if(E->NODE[level][i] & ( OFFSIDE )) 
	    continue;
	    
	    if(2==dofs){
	      f1=F1[i];
	      f2=F2[i];
	      C=E->Node_map_3[level]+i*MAX_NODE_INT_2D;
	      B11=E->Node_k11[level]+i*MAX_NODE_INT_2D;
	      B12=E->Node_k12[level]+i*MAX_NODE_INT_2D;
	      B21=E->Node_k21[level]+i*MAX_NODE_INT_2D;
	      B22=E->Node_k22[level]+i*MAX_NODE_INT_2D;
	    
	      u1 = u2 = 0.0; 
	      for(j=0;j<MAX_NODE_INT_2D;j++) {
		u1 += B11[j] * U1[C[j]] +  B12[j] * U2[C[j]] ;
		u2 += B21[j] * U1[C[j]] +  B22[j] * U2[C[j]] ;
	      }
	      du1 = (f1 - u1) * E->BI1[level][i];
	      du2 = (f2 - u2) * E->BI2[level][i];
	      U1[i] += du1;
	      U2[i] += du2;

	      if(count==inner_cycles) {
		for(j=0;j<MAX_NODE_INT_2D;j++) {
		  au1[C[j]] += B11[j] * U1[i] + B21[j] * U2[i];
		  au2[C[j]] += B12[j] * U1[i] + B22[j] * U2[i]; 
		}
	      }
	    }
	    else if(3==dofs && 2==dims){
	      f1=F1[i];
	      f2=F2[i];
	      f3=F3[i];
	      C=E->Node_map_3[level]+i*MAX_NODE_INT_2D;
	      B11=E->Node_k11[level]+i*MAX_NODE_INT_2D;
	      B12=E->Node_k12[level]+i*MAX_NODE_INT_2D;
	      B13=E->Node_k13[level]+i*MAX_NODE_INT_2D;
	      B21=E->Node_k21[level]+i*MAX_NODE_INT_2D;
	      B22=E->Node_k22[level]+i*MAX_NODE_INT_2D;
	      B23=E->Node_k23[level]+i*MAX_NODE_INT_2D;
	      B31=E->Node_k31[level]+i*MAX_NODE_INT_2D;
	      B32=E->Node_k32[level]+i*MAX_NODE_INT_2D;
	      B33=E->Node_k33[level]+i*MAX_NODE_INT_2D;
	    
	      u1 = u2 =  u3 = 0.0;
	      for(j=0;j<MAX_NODE_INT_2D;j++) {
		u1 += B11[j] * U1[C[j]] +  B12[j] * U2[C[j]] +  B13[j] * U3[C[j]] ;
		u2 += B21[j] * U1[C[j]] +  B22[j] * U2[C[j]] +  B23[j] * U3[C[j]] ;
		u3 += B31[j] * U1[C[j]] +  B32[j] * U2[C[j]] +  B33[j] * U3[C[j]] ; 
	      }
	      
	      U1[i] += (f1 - u1) * E->BI1[level][i];
	      U2[i] += (f2 - u2) * E->BI2[level][i];
	      U3[i] += (f3 - u3) * E->BI3[level][i];

	      if(count==inner_cycles) {
		for(j=0;j<MAX_NODE_INT_2D;j++) {
		  au1[C[j]] += B11[j] * U1[i] + B21[j] * U2[i] + B31[j] * U3[i];
		  au2[C[j]] += B12[j] * U1[i] + B22[j] * U2[i] + B32[j] * U3[i]; 
		  au3[C[j]] += B13[j] * U1[i] + B23[j] * U2[i] + B33[j] * U3[i]; 
		}
	      }
	    }
	    else if(3==dofs && 3==dims) {
	      f1=F1[i];
	      f2=F2[i];
	      f3=F3[i];
	      C=E->Node_map_3[level]+i*MAX_NODE_INT_3D; 
	      B11=E->Node_k11[level]+i*MAX_NODE_INT_3D;
	      B12=E->Node_k12[level]+i*MAX_NODE_INT_3D;
	      B13=E->Node_k13[level]+i*MAX_NODE_INT_3D;
	      B21=E->Node_k21[level]+i*MAX_NODE_INT_3D;
	      B22=E->Node_k22[level]+i*MAX_NODE_INT_3D;
	      B23=E->Node_k23[level]+i*MAX_NODE_INT_3D;
	      B31=E->Node_k31[level]+i*MAX_NODE_INT_3D;
	      B32=E->Node_k32[level]+i*MAX_NODE_INT_3D;
	      B33=E->Node_k33[level]+i*MAX_NODE_INT_3D;
	      
	      u1 = u2 = u3 = 0.0;
	      for(j=0;j<MAX_NODE_INT_3D;j++) {
		u1 += B11[j] * U1[C[j]] +  B12[j] * U2[C[j]]  +  B13[j] * U3[C[j]] ;
		u2 += B21[j] * U1[C[j]] +  B22[j] * U2[C[j]]  +  B23[j] * U3[C[j]] ;
		u3 += B31[j] * U1[C[j]] +  B32[j] * U2[C[j]]  +  B33[j] * U3[C[j]] ;
	      }

	      U1[i] += (f1 - u1) * E->BI1[level][i];
	      U2[i] += (f2 - u2) * E->BI2[level][i];
	      U3[i] += (f3 - u3) * E->BI3[level][i];
	      
	      if(count==inner_cycles)
		for(j=0;j<MAX_NODE_INT_3D;j++) {
		  au1[C[j]] += B11[j] * U1[i] + B21[j] * U2[i] + B31[j] * U3[i];
		  au2[C[j]] += B12[j] * U1[i] + B22[j] * U2[i] + B32[j] * U3[i]; 
		  au3[C[j]] += B13[j] * U1[i] + B23[j] * U2[i] + B33[j] * U3[i]; 
		}
	    }
	    else /* if (dofs==6) */ {
	      f1=F1[i];
	      f2=F2[i];
	      f3=F3[i];
	      f4=F4[i];
	      f5=F5[i];
	      f6=F6[i];
	      C=E->Node_map_3[level]+i*MAX_NODE_INT_3D;
	      B11=E->Node_k11[level]+i*MAX_NODE_INT_3D;
	      B12=E->Node_k12[level]+i*MAX_NODE_INT_3D;
	      B13=E->Node_k13[level]+i*MAX_NODE_INT_3D;
	      B14=E->Node_k14[level]+i*MAX_NODE_INT_3D;
	      B15=E->Node_k15[level]+i*MAX_NODE_INT_3D;
	      B16=E->Node_k16[level]+i*MAX_NODE_INT_3D;
	      B21=E->Node_k21[level]+i*MAX_NODE_INT_3D;
	      B22=E->Node_k22[level]+i*MAX_NODE_INT_3D;
	      B23=E->Node_k23[level]+i*MAX_NODE_INT_3D;
	      B24=E->Node_k24[level]+i*MAX_NODE_INT_3D;
	      B25=E->Node_k25[level]+i*MAX_NODE_INT_3D;
	      B26=E->Node_k26[level]+i*MAX_NODE_INT_3D;
	      B31=E->Node_k31[level]+i*MAX_NODE_INT_3D;
	      B32=E->Node_k32[level]+i*MAX_NODE_INT_3D;
	      B33=E->Node_k33[level]+i*MAX_NODE_INT_3D;
	      B34=E->Node_k34[level]+i*MAX_NODE_INT_3D;
	      B35=E->Node_k35[level]+i*MAX_NODE_INT_3D;
	      B36=E->Node_k36[level]+i*MAX_NODE_INT_3D;
	      B41=E->Node_k41[level]+i*MAX_NODE_INT_3D;
	      B42=E->Node_k42[level]+i*MAX_NODE_INT_3D;
	      B43=E->Node_k43[level]+i*MAX_NODE_INT_3D;
	      B44=E->Node_k44[level]+i*MAX_NODE_INT_3D;
	      B45=E->Node_k45[level]+i*MAX_NODE_INT_3D;
	      B46=E->Node_k46[level]+i*MAX_NODE_INT_3D;
	      B51=E->Node_k51[level]+i*MAX_NODE_INT_3D;
	      B52=E->Node_k52[level]+i*MAX_NODE_INT_3D;
	      B53=E->Node_k53[level]+i*MAX_NODE_INT_3D;
	      B54=E->Node_k54[level]+i*MAX_NODE_INT_3D;
	      B55=E->Node_k55[level]+i*MAX_NODE_INT_3D;
	      B56=E->Node_k56[level]+i*MAX_NODE_INT_3D;
	      B61=E->Node_k61[level]+i*MAX_NODE_INT_3D;
	      B62=E->Node_k62[level]+i*MAX_NODE_INT_3D;
	      B63=E->Node_k63[level]+i*MAX_NODE_INT_3D;
	      B64=E->Node_k64[level]+i*MAX_NODE_INT_3D;
	      B65=E->Node_k65[level]+i*MAX_NODE_INT_3D;
	      B66=E->Node_k66[level]+i*MAX_NODE_INT_3D;
	      
	      u1 = u2 = u3 = u4 = u5 = u6 = 0.0;
	      for(j=0;j<MAX_NODE_INT_3D;j++) {
		u1 += B11[j]*U1[C[j]]+B12[j]*U2[C[j]]+B13[j]*U3[C[j]]+B14[j]*U4[C[j]]+B15[j]*U5[C[j]]+B16[j]*U6[C[j]] ;
		u2 += B21[j]*U1[C[j]]+B22[j]*U2[C[j]]+B23[j]*U3[C[j]]+B24[j]*U4[C[j]]+B25[j]*U5[C[j]]+B26[j]*U6[C[j]];
		u3 += B31[j]*U1[C[j]]+B32[j]*U2[C[j]]+B33[j]*U3[C[j]]+B34[j]*U4[C[j]]+B35[j]*U5[C[j]]+B36[j]*U6[C[j]];
		u4 += B41[j]*U1[C[j]]+B42[j]*U2[C[j]]+B43[j]*U3[C[j]]+B44[j]*U4[C[j]]+B45[j]*U5[C[j]]+B46[j]*U6[C[j]] ;
		u5 += B51[j]*U1[C[j]]+B52[j]*U2[C[j]]+B53[j]*U3[C[j]]+B54[j]*U4[C[j]]+B55[j]*U5[C[j]]+B56[j]*U6[C[j]];
		u6 += B61[j]*U1[C[j]]+B62[j]*U2[C[j]]+B63[j]*U3[C[j]]+B64[j]*U4[C[j]]+B65[j]*U5[C[j]]+B66[j]*U6[C[j]];
	      }	      
	      U1[i] += (f1 - u1) * E->BI1[level][i];
	      U2[i] += (f2 - u2) * E->BI2[level][i];
	      U3[i] += (f3 - u3) * E->BI3[level][i];
	      U4[i] += (f4 - u4) * E->BI4[level][i];
	      U5[i] += (f5 - u5) * E->BI5[level][i];
	      U6[i] += (f6 - u6) * E->BI6[level][i];
	      
	      if(count==inner_cycles)
		for(j=0;j<MAX_NODE_INT_3D;j++) {
		  au1[C[j]] += B11[j]*U1[i] + B21[j]*U2[i] + B31[j]*U3[i] + B41[j]*U4[i] + B51[j]*U5[i] + B61[j]*U6[i];
		  au2[C[j]] += B12[j]*U1[i] + B22[j]*U2[i] + B32[j]*U3[i] + B42[j]*U4[i] + B52[j]*U5[i] + B62[j]*U6[i]; 
		  au3[C[j]] += B13[j]*U1[i] + B23[j]*U2[i] + B33[j]*U3[i] + B43[j]*U4[i] + B53[j]*U5[i] + B63[j]*U6[i]; 
		  au4[C[j]] += B14[j]*U1[i] + B24[j]*U2[i] + B34[j]*U3[i] + B44[j]*U4[i] + B54[j]*U5[i] + B64[j]*U6[i];
		  au5[C[j]] += B15[j]*U1[i] + B25[j]*U2[i] + B35[j]*U3[i] + B45[j]*U4[i] + B55[j]*U5[i] + B65[j]*U6[i]; 
		  au6[C[j]] += B16[j]*U1[i] + B26[j]*U2[i] + B36[j]*U3[i] + B46[j]*U4[i] + B56[j]*U5[i] + B66[j]*U6[i]; 
		}
	    }
	  }
      residual=0.0;
      for(i=1;i<=num_nodes;i++) {
	if(E->NODE[level][i] & OFFSIDE) 
	  continue;  
	residual += (F1[i] - au1[i]) * (F1[i] - au1[i]) + (F2[i] - au2[i]) * (F2[i] - au2[i]);

	if(3==dofs) {
	  residual += (F3[i] - au3[i]) * (F3[i] - au3[i]);
	}
	else if (6==dofs) {
	  residual += (F3[i]-au3[i])*(F3[i]-au3[i]) + (F4[i]-au4[i])*(F4[i]-au4[i])
	    + (F5[i]-au5[i])*(F5[i]-au5[i]) + (F6[i]-au6[i])*(F6[i] - au6[i]);
	}
      }
    }
    delta_res = fabs(delta_res-residual);
  }

  if(level != E->mesh.levmin && count2 > 100)
    fprintf(stderr,"%d: GS ... %d steps res %g v %g\n",level,count2,residual,residual1);

  return;
}


