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


#include <signal.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"


/* ===================================
   Functions which set up details 
   common to all problems follow ...
   ===================================  */

void allocate_common_vars(
	   struct All_variables *E
)
{ 
  int i,j,l;
  int elx,ely,nxyz;
  
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int max_node = max_node_interaction[dims];
  
  E->mesh.fnodal_malloc_size = (E->mesh.nno+2)*sizeof(standard_precision);

  E->F = (standard_precision *) Malloc0((dofs*E->mesh.nno+1)*sizeof(standard_precision));
  E->Fb = (standard_precision *) Malloc0((dofs*E->mesh.nno+1)*sizeof(standard_precision));
  E->Q = (standard_precision *) Malloc0((E->mesh.npno+1)*sizeof(standard_precision));
  E->mass = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->trmass = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  
  E->Have.T = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  E->Have.Vi = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  E->Have.vrms = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  E->Have.strs = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  E->Have.buoy = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  E->Have.Pres = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
  
  E->T = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->nQ = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->edot = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->strs = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->strd = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->strd1 = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  E->depl = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision)); 
  /*RAA: 24/09/02, depl is melting stuff by C. O'Neill */

  if(2==E->mesh.nsd)
    E->Psi=(standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  
  for(i=1;i<=E->mesh.nno;i++)
    E->T[i] = 0.0;
  
  E->Pstrain = (standard_precision *) Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  
  for(i=1;i<=dofs;i++)  {
    E->Have.V[i] = (standard_precision *) Malloc0((E->mesh.noz+1)*sizeof(standard_precision));
    E->V[i] = (standard_precision *) Malloc0( (E->mesh.nno+2)*sizeof(standard_precision) );
    E->V1[i] = (standard_precision *) Malloc0( (E->mesh.nno+2)*sizeof(standard_precision) );
    E->VB[i] = (standard_precision *)  Malloc0( (E->mesh.nno+2)*sizeof(standard_precision) );
  }
  for(i=1;i<=dims;i++)  {
    E->x[i] = (standard_precision *)  Malloc0((E->mesh.nno+2)*sizeof(standard_precision));
    E->sx[i] = (standard_precision *)  Malloc0((E->mesh.nno+2)*sizeof(standard_precision));
  }

  if(E->control.SPHERE || E->control.CYLINDER) {
    for(i=1;i<=dofs;i++)
      E->sv[i] = (standard_precision *)  Malloc0((E->mesh.nno+2)*sizeof(standard_precision));
    
    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
      E->curvilinear.cosph[i] = (standard_precision *) Malloc0((E->mesh.nox+2)*sizeof(standard_precision));
      E->curvilinear.sinph[i] = (standard_precision *) Malloc0((E->mesh.nox+2)*sizeof(standard_precision));
      if(E->control.SPHERE) {
	E->curvilinear.sint[i] = (standard_precision *) Malloc0((E->mesh.noy+2)*sizeof(standard_precision));
	E->curvilinear.cost[i] = (standard_precision *) Malloc0((E->mesh.noy+2)*sizeof(standard_precision));
      }
    }
  }
  else {
    for(i=1;i<=dofs;i++)
      E->sv[i] = E->V[i];
  }

  E->TB = (standard_precision *)  Malloc0( (E->mesh.nno+2)*sizeof(standard_precision) );
  E->Bpi = (higher_precision *)Malloc0((E->mesh.npno+1)*sizeof(higher_precision));
  E->eco = (struct COORD *) Malloc0((E->mesh.nno+2)*sizeof(struct COORD));
  E->ien = (struct IEN *) Malloc0((E->mesh.nno+2)*sizeof(struct IEN));
  E->ienp = (struct IEN *) Malloc0((E->mesh.nno+2)*sizeof(struct IEN));
  E->id =  (struct ID *) Malloc0((E->mesh.nno+2)*sizeof(struct ID));
  E->elem = (struct ELEMENT *) Malloc0((E->mesh.nno+2)*sizeof(struct ELEMENT));
  
  E->node =  (unsigned int *) Malloc0((E->mesh.nno+2)*sizeof(unsigned int));
  E->tw = (standard_precision  *) Malloc0((E->mesh.nno+2)*sizeof(standard_precision));
  E->Bi = (higher_precision *) Malloc0((E->mesh.nno+2)*sizeof(higher_precision)); 
  E->idd = (struct ID *) Malloc0((E->mesh.NNO[E->mesh.levmin]+1)*sizeof(struct ID));
  
  
  /* Aliases for functions only interested in the highest mg level */
  
  E->ECO[E->mesh.levmax] = E->eco; 
  E->IEN[E->mesh.levmax]=E->ien;
  E->IENP[E->mesh.levmax]=E->ienp;
  E->ID[E->mesh.levmax]=E->id;
  E->NODE[E->mesh.levmax]=E->node;
  E->TW[E->mesh.levmax]=E->tw;
  E->BI[E->mesh.levmax]=E->Bi;
  E->ELEM[E->mesh.levmax]=E->elem;
  E->NQ[E->mesh.levmax]=E->nQ;

  for(i=1;i<=dims;i++) {
    E->X[E->mesh.levmax][i]=E->x[i];
    E->SX[E->mesh.levmax][i]=E->sx[i];
  }
  for(i=1;i<=dofs;i++)  {
    E->VV[E->mesh.levmax][i] = E->V[i];
    E->Vb[i][E->mesh.levmax] = E->VB[i];
  }
  
  
  E->BPI[E->mesh.levmax] = E->Bpi;
  E->QQ[E->mesh.levmax] = E->Q;
  E->FF[E->mesh.levmax] = E->F;
  E->FFb[E->mesh.levmax] = E->Fb;
  E->Mass[E->mesh.levmax] = E->mass;
  E->Trmass[E->mesh.levmax] = E->trmass;
  
  E->Tb[E->mesh.levmax] = E->TB;
  
  /* These things are defined in multigrid and have high level aliases already allocated*/

  for(i=E->mesh.levmin;i<E->mesh.levmax;i++)   {
    E->ECO[i] = (struct COORD *) Malloc0((E->mesh.NNO[i]+2)*sizeof(struct COORD));
    E->IEN[i] = (struct IEN *) Malloc0((E->mesh.NNO[i]+2)*sizeof(struct IEN));
    E->IENP[i] = (struct IEN *) Malloc0((E->mesh.NNO[i]+2)*sizeof(struct IEN));
    E->ID[i]  = (struct ID *) Malloc0((E->mesh.NNO[i]+2)*sizeof(struct ID));
    E->ELEM[i]  = (struct ELEMENT *) Malloc0((E->mesh.NNO[i]+2)*sizeof(struct ELEMENT));
    
    E->NODE[i] = (unsigned int *) Malloc0((E->mesh.NNO[i]+2)*sizeof(unsigned int));
    E->TW[i]  = (standard_precision  *) Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    E->BPI[i] = (higher_precision *) Malloc0((E->mesh.NPNO[i]+2)*sizeof(higher_precision));
   
    E->NQ[i] = (standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    E->QQ[i] = (standard_precision *)Malloc0((E->mesh.NPNO[i]+2)*sizeof(standard_precision));
    E->FF[i] = (standard_precision *)Malloc0((E->mesh.NNO[i]*dofs+2)*sizeof(standard_precision));
    E->FFb[i] = (standard_precision *)Malloc0((E->mesh.NNO[i]*dofs+2)*sizeof(standard_precision));
    E->Tb[i]=(standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    E->Mass[i] = (standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    E->Trmass[i] = (standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    
    for(j=1;j<=dims;j++) {
      E->X[i][j]=(standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
      E->SX[i][j]=(standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    }
    for(j=1;j<=dofs;j++)  {
      E->VV[i][j]=(standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
      E->Vb[j][i]=(standard_precision *)Malloc0((E->mesh.NNO[i]+2)*sizeof(standard_precision));
    }
  }
  
  /* Curvilinear case, boundary condition rotation matrices */
  
  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)
    E->curvilinear.NODE_R[i] = (higher_precision **)Malloc0((E->mesh.NNO[i]+1) * sizeof(higher_precision *));
  
  E->curvilinear.node_r = E->curvilinear.NODE_R[E->mesh.levmax];
  
  /* These things are unique to multigrid */
  
  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)  {
    E->EL[i]  = (struct SUBEL *) Malloc0((E->mesh.NEL[i]+2)*sizeof(struct SUBEL));
    E->NEI[i].nels  = (int *) Malloc0((E->mesh.NNO[i]+2)*sizeof(int));
    E->NEI[i].lnode  = (int *) Malloc0((E->mesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));
    E->NEI[i].element  = (int *) Malloc0((E->mesh.NNO[i]+2)*enodes[E->mesh.nsd]*sizeof(int));  

    E->index[i] =  (int *) Malloc0((E->mesh.NNO[i]+2)*sizeof(int));

    /* E->update_elt_k[i] = (char *) Malloc0((E->mesh.NEL[i]+2)*sizeof(struct char)); */
 
    
    E->elt_del[i]=(struct EG *)Malloc0((E->mesh.NEL[i]+1)*sizeof(struct EG));
    
    E->Node_k11[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
    E->Node_k21[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
    E->Node_k12[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
    E->Node_k22[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
    
    E->BI[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]*dofs+2)*sizeof(higher_precision)); 
    E->BI1[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 
    E->BI2[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 

    if(dofs==3) {
      E->Node_k31[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k32[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k13[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k23[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k33[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      
      E->BI3[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+10)*sizeof(higher_precision)); 
    }
    else if(dofs==6) {
      E->Node_k13[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k14[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k15[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k16[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k23[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k24[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k25[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k26[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k31[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k32[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k33[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k34[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k35[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k36[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k41[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k42[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k43[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k44[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k45[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k46[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k51[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k52[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k53[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k54[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k55[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k56[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k61[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k62[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k63[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k64[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k65[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      E->Node_k66[i] = (higher_precision *) Malloc0((E->mesh.NNO[i] +2) * max_node * sizeof(higher_precision));
      
      E->BI3[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 
      E->BI4[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 
      E->BI5[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 
      E->BI6[i] = (higher_precision *) Malloc0((E->mesh.NNO[i]+1)*sizeof(higher_precision)); 
    }
    E->Node_map_3[i] = (int *) Malloc0((E->mesh.NNO[i]+5) * max_node * sizeof(int));      
    E->control.B_is_good[i] = 0;
  }

  for(l=E->mesh.levmin;l<=E->mesh.levmax;l++) {
    for(i=1;i<=E->mesh.NNO[l];i++) {
      E->NODE[l][i] = (INTX | INTY | INTZ);  /* and any others ... */
      
      E->curvilinear.NODE_R[l][i] = (higher_precision *) NULL;
      
      E->NQ[l][i] = 0.0;
      E->TW[l][i]= 0.0;
      E->Tb[l][i] = 0.0;
      for(j=1;j<=dofs;j++) {
	E->Vb[j][l][i] = 0.0;
	E->VV[l][j][i] = 0.0;
      }
      for(j=1;j<=dims;j++) {
	E->X[l][j][i] = 0.0;
      }
    }

    for(i=1;i<=E->mesh.NPNO[l];i++)
      E->QQ[l][i] = 0.0;
  }
  return;
}

void allocate_observables_storage(
  struct All_variables *E
)
{
  

  E->slice.tpg = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.tpgb = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.tpgk = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.tpgbk = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grv  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geo  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grvb  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grvt  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geob  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geot  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grvk  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geok  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grvbk  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.grvtk  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geobk  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.geotk  = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.shflux = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.bhflux = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.vxsurf[1] = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.vxsurf[2] = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
  E->slice.vxsurf[3] = (standard_precision *)Malloc0((E->mesh.nox * E->mesh.noy+2)*sizeof(standard_precision));
}
