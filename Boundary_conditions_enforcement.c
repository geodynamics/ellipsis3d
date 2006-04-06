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


void s_strip_bcs_from_residual(
			       struct All_variables *E,
			       standard_precision *Res,
			       int level
			       )
{
  int i;
  
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nno=E->mesh.NNO[level];
  
  for(i=1;i<=nno;i++) { 
    if(E->NODE[level][i] & OFFSIDE  )
      continue;
    
    if ( (E->NODE[level][i] & BC1) != 0 )
      Res[ E->ID[level][i].doff[1] ] = 0.0;
    if ( (E->NODE[level][i] & BC2) != 0 )
      Res[ E->ID[level][i].doff[2] ] = 0.0;
    if(dofs==3) {
      if ( (E->NODE[level][i] & BC3) != 0)
	Res[ E->ID[level][i].doff[3] ] = 0.0;
    }
    else if(dofs==6) {
      if ( (E->NODE[level][i] & BC3) != 0)
	Res[ E->ID[level][i].doff[3] ] = 0.0;
      if ( (E->NODE[level][i] & BC4) != 0 )
	Res[ E->ID[level][i].doff[4] ] =  0.0;
      if ( (E->NODE[level][i] & BC5) != 0 )
	Res[ E->ID[level][i].doff[5] ] = 0.0;
      if ( (E->NODE[level][i] & BC6) != 0 )
	Res[ E->ID[level][i].doff[6] ] = 0.0;
    }
  }
  return;
}

void strip_bcs_from_residual_6(
			       struct All_variables *E,
			       standard_precision *U1,
			       standard_precision *U2,
			       standard_precision *U3,
			       standard_precision *U4,
			       standard_precision *U5,
			       standard_precision *U6,
			       int level
			       )
{
  int node;
     
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nno=E->mesh.NNO[level];

  if(3==dofs) {   /* 3D Classical or 2D Cosserat */
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = 0.0; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = 0.0; 
	if (E->NODE[level][node] & BC3)
	  U3[node] = 0.0; 
      }
    }
  }
  
  else if(2==dofs) {   /* 2D Classical */
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = 0.0; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = 0.0; 
      }
    }
  }
  
  else /* if(6==dofs) */ {   /* 3D Cosserat */ 
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = 0.0; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = 0.0; 
	if (E->NODE[level][node] & BC3)
	  U3[node] = 0.0; 
	if (E->NODE[level][node] & BC4)
	  U4[node] = 0.0; 
	if (E->NODE[level][node] & BC5)
	  U5[node] = 0.0; 
	if (E->NODE[level][node] & BC6)
	  U6[node] = 0.0; 
      }
    }
  }
  return;
}


void velocity_conform_bcs_6(
			    struct All_variables *E,
			    standard_precision *U1,
			    standard_precision *U2,
			    standard_precision *U3,
			    standard_precision *U4,
			    standard_precision *U5,
			    standard_precision *U6,
			    int level
			    )
{
  int i;
  
  const int dims=E->mesh.nsd;
  const int dofs=E->mesh.dof;
  const int nno=E->mesh.NNO[level];
  
  for(i=1;i<=nno;i++) {
    if(E->NODE[level][i] & OFFSIDE)
      continue;
    
    if ( (E->NODE[level][i] & BC1) != 0 ) {
      U1[i] = E->Vb[1][level][i];
    }
    if ( (E->NODE[level][i] & BC2) != 0 )
      U2[i] = E->Vb[2][level][i];
    if(dofs==3) {
      if ((E->NODE[level][i] & BC3) != 0)
	U3[i] = E->Vb[3][level][i];
    }
    else if(dofs==6) {
      if ((E->NODE[level][i] & BC3) != 0)
	U3[i] = E->Vb[3][level][i];
      if ((E->NODE[level][i] & BC4) != 0)
	U4[i] = E->Vb[4][level][i];
      if ((E->NODE[level][i] & BC5) != 0)
	U5[i] = E->Vb[5][level][i];
      if ((E->NODE[level][i] & BC6) != 0)
	U6[i] = E->Vb[6][level][i];
    }
  }
  return;
}

void velocities_conform_bcs_6(
			      struct All_variables *E,
			      standard_precision *U1,
			      standard_precision *U2,
			      standard_precision *U3,
			      standard_precision *U4,
			      standard_precision *U5,
			      standard_precision *U6,
			      int level
			      )
{
  int node,d;
  
  const int nno = E->mesh.NNO[level];
  const int dofs = E->mesh.dof;
  
  if(3==dofs) {   /* 3D Classical or 2D Cosserat */
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = E->Vb[1][level][node]; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = E->Vb[2][level][node]; 
	if (E->NODE[level][node] & BC3)
	  U3[node] = E->Vb[3][level][node]; 
      }
    }
  }
  else if(2==dofs) {   /* 2D Classical */
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = E->Vb[1][level][node]; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = E->Vb[2][level][node]; 
      }
    }
  } 
  else /* if(6==dofs) */ {   /* 3D Cosserat */
    for(node=1;node<=nno;node++) {
      if(!(E->NODE[level][node] & OFFSIDE)) {
	if (E->NODE[level][node] & BC1)
	  U1[node] = E->Vb[1][level][node]; 
	if (E->NODE[level][node] & BC2)
	  U2[node] = E->Vb[2][level][node]; 
	if (E->NODE[level][node] & BC3)
	  U3[node] = E->Vb[3][level][node]; 
	if (E->NODE[level][node] & BC4)
	  U4[node] = E->Vb[4][level][node]; 
	if (E->NODE[level][node] & BC5)
	  U5[node] = E->Vb[5][level][node]; 
	if (E->NODE[level][node] & BC6)
	  U6[node] = E->Vb[6][level][node]; 
      }
    }
  }
  return;
}

void temperatures_conform_bcs(
     struct All_variables *E,
     standard_precision *T
)
{
    int node;
    unsigned int type;

   for(node=1;node<=E->mesh.nno;node++)  {
	if(E->node[node] & ( OFFSIDE  ))
	    continue;
      
	if(E->node[node] & TBD)
	  T[node] = E->TB[node];
	
    }

/*RAA: 4/4/02, the stuff below may not be impt, comment it all out*/
/*RAA: 7/6/01, temp bcs are not showing up on nodes which are OFFSIDE or
  PER_OFFSIDE - original (non-corrected) code is commented out above*/
    /*for(node=1;node<=E->mesh.nno;node++)  {  */
 /*    if((E->node[node] & ( OFFSIDE )) & (E->node[node] & ~(TBD)))
         continue;
 */      
       /*if(E->node[node] & (OFFSIDE)) {  */
          /*if (E->node[node] & TBD)  {  */
             /*T[node] = E->TB[node];   */
 /*        fprintf(stderr,"in the top offside loop, big shot!\n"); */
/*RAA, comment not seen, obviously, OFFSIDE shuts off TBD*/
          /*}  */
          /*else if (E->node[node] & ~TBD) {   */
  /*    fprintf(stderr,"in the bottom offside loop, big shot!\n"); */
             /*continue;  */
          /*}  */
       /*}  */
       /*if(E->node[node] & (PER_OFFSIDE)) {  */
          /*if (E->node[node] & TBD)  {  */
             /*T[node] = E->TB[node];   */
  /*     fprintf(stderr,"in the top per_offside loop, short stuff!\n"); */
  /*RAA, comment not seen, obviously, OFFSIDE shuts off TBD*/
          /*}  */
          /*else if (E->node[node] & ~TBD) {   */
   /*     fprintf(stderr,"in the bottom per_offside loop, short stuff!\n"); */
             /*continue;  */
          /*}  */
       /*}  */


/*    else if((E->node[node] & ( OFFSIDE )) & (E->node[node] & TBD))
       T[node] = E->TB[node];
 */     
       /*else if(E->node[node] & TBD)  */
          /*T[node] = E->TB[node];  */
    /*}  */

  return;
}

void check_bc_consistency(
     struct All_variables *E
)
{
  int i,lev;

  const int dofs = E->mesh.dof ;
  const int dims = E->mesh.nsd ;

  for(i=1;i<=E->mesh.nno;i++) {
    if ((E->node[i] & BC1) && (E->node[i] & SBC1))
      printf("Inconsistent x velocity bc at %d\n",i);
    if ((E->node[i] & BC2) && (E->node[i] & SBC2))
      printf("Inconsistent z velocity bc at %d\n",i);
    if(dims==2 && dofs==3) {
      if ((E->node[i] & BC3) && (E->node[i] & SBC3))
	printf("Inconsistent y rotation bc at %d\n",i);
    }
    else if(dims==3 && dofs==3) {
      if ((E->node[i] & BC3) && (E->node[i] & SBC3))
	printf("Inconsistent y velocity bc at %d\n",i);
    }
    else if(6==dofs) {
      if ((E->node[i] & BC3) && (E->node[i] & SBC3))
	printf("Inconsistent y velocity bc at %d\n",i);
      if ((E->node[i] & BC4) && (E->node[i] & SBC4))
	printf("Inconsistent x rotation bc at %d\n",i);
      if ((E->node[i] & BC5) && (E->node[i] & SBC5))
	printf("Inconsistent z rotation bc at %d\n",i);
      if ((E->node[i] & BC6) && (E->node[i] & SBC6))
	printf("Inconsistent y rotation bc at %d\n",i);
    }
    
    if ((E->node[i] & TBD) && (E->node[i] & FBX))
      printf("Inconsistent x temperature bc at %d\n",i);
    if ((E->node[i] & TBD) && (E->node[i] & FBZ))
      printf("Inconsistent z temperature bc at %d\n",i);
    if ((E->node[i] & TBD) && (E->node[i] & FBY))
      printf("Inconsistent y temperature bc at %d\n",i);
  }

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)
    for(i=1;i<=E->mesh.NNO[lev];i++) {
      if ((E->NODE[lev][i] & BC1) && (E->NODE[lev][i]  & SBC1))
	printf("Inconsistent x velocity bc at %d,%d\n",lev,i);
      if ((E->NODE[lev][i]  & BC2) && (E->NODE[lev][i]  & SBC2))
	printf("Inconsistent z velocity bc at %d,%d\n",lev,i);
      
      if(dims==2 && dofs==3) {
	if ((E->NODE[lev][i] & BC3) && (E->NODE[lev][i] & SBC3))
	  printf("Inconsistent y rotation bc at %d,%d\n",lev,i);
      }
      else if(dims==3 && dofs==3) {
	if ((E->NODE[lev][i] & BC3) && (E->NODE[lev][i] & SBC3))
	  printf("Inconsistent y velocity bc at %d,%d\n",lev,i);
      }
      else if(6==dofs) {
	if ((E->NODE[lev][i] & BC3) && (E->NODE[lev][i] & SBC3))
	  printf("Inconsistent y velocity bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i] & BC4) && (E->NODE[lev][i] & SBC4))
	  printf("Inconsistent x rotation bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i] & BC5) && (E->NODE[lev][i] & SBC5))
	  printf("Inconsistent z rotation bc at %d,%d\n",lev,i);
	if ((E->NODE[lev][i] & BC6) && (E->NODE[lev][i] & SBC6))
	  printf("Inconsistent y rotation bc at %d,%d\n",lev,i);
      }
    }
  return;
}
