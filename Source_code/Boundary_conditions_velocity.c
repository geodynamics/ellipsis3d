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
#include <math.h>/* ========================================== */

void velocity_boundary_conditions(
				  struct All_variables *E
				  )
{
    void velocity_refl_vert_bc();
    void velocity_imp_vert_bc();
    void horizontal_bc();
    void velocity_apply_periodic_bcs();
    int lv; 
    int node,d,n1,n2,i,j,n;
 
    void arbitrary_bc_rectangle();
    void arbitrary_bc_circle();
    void arbitrary_bc_harmonic();
    void arbitrary_bc_polynomial();
    void arbitrary_bc_rectangle_file();
    void arbitrary_bc_circle_file();
    void arbitrary_bc_harmonic_file();
    void arbitrary_bc_polynomial_file();
    void free_surface_boundary_conditions();
    void store_node_R();

    struct RectBc RECT;

    const int dims = E->mesh.nsd ;
    const int dofs = E->mesh.dof ;

    /* Default: Free slip top (0) and bottom (1), reflecting sidewalls */
    
    RECT.numb=2;  /* Z-normal - MOBILE */
    RECT.norm[0]=RECT.norm[1]='Z';
    RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
    RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
    RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
    RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
    RECT.intercept[0]=E->mesh.layer0[2];
    RECT.mag[0]=0.0;
    RECT.intercept[1]=E->mesh.layer1[2];
    RECT.mag[1]=0.0;

    if(2==dims) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,SBC1,BC1,1);
      if(3==dofs)
	arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,SBC3,BC3,1);
    }
    else {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,SBC1,BC1,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,SBC3,BC3,1);
      if(6==dofs) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,SBC4,BC4,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,SBC5,BC5,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,SBC6,BC6,1);
      }
    }

    RECT.numb=2; /* X-normal -MOBILE */ 
    RECT.norm[0]=RECT.norm[1]='X';
    RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
    RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
    RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
    RECT.aa2[0]=RECT.aa2[1]=  1.0e32;

    RECT.intercept[0]=E->mesh.layer0[1];
    RECT.mag[0]=0.0;
    RECT.intercept[1]=E->mesh.layer1[1];
    RECT.mag[1]=0.0;

    if(2==dims) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,SBC2,BC2,1);
      if(3==dofs)
	arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,SBC3,BC3,1);
    }
    else {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,SBC2,BC2,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,SBC3,BC3,1);
      if(dofs==6) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,SBC4,BC4,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,SBC5,BC5,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,SBC6,BC6,1);
      }
    }

    if(3==dims) {
      RECT.numb=2; /* Y-normal -MOBILE*/
      RECT.norm[0]=RECT.norm[1]='Y';
      RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
      RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
      RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
      RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
      
      RECT.intercept[0]=E->mesh.layer0[3];
      RECT.mag[0]=0.0;
      RECT.intercept[1]=E->mesh.layer1[3];
      RECT.mag[1]=0.0;
      
      arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,SBC1,BC1,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,SBC2,BC2,1);
      if(6==dofs) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,SBC4,BC4,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,SBC5,BC5,1);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,SBC6,BC6,1);
      }
    }
  
    RECT.numb=2;  /* Z-normal - FIXED */
    RECT.norm[0]=RECT.norm[1]='Z';
    RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
    RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
    RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
    RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
    RECT.intercept[0]=E->mesh.layer0[2];
    RECT.mag[0]=0.0;
    RECT.intercept[1]=E->mesh.layer1[2];
    RECT.mag[1]=0.0;

    arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
/*    if(dofs==3 && 2==dims)
      arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
    else */ if(6==dofs) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,BC4,SBC4,0);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,BC5,SBC5,0);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,BC6,SBC6,0);
    }

    RECT.numb=2; /* X-normal -FIXED */ 
    RECT.norm[0]=RECT.norm[1]='X';
    RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
    RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
    RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
    RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
    
    RECT.intercept[0]=E->mesh.layer0[1];
    RECT.mag[0]=0.0;
    RECT.intercept[1]=E->mesh.layer1[1];
    RECT.mag[1]=0.0;

    arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
/*    if(dofs==3 && 2==dims)
      arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
    else */if(6==dofs) {
      arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,BC4,SBC4,0);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,BC5,SBC5,0);
      arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,BC6,SBC6,0);
    }
    
    if(3==dims) {
      RECT.numb=2; /* Y-normal -FIXED */
      RECT.norm[0]=RECT.norm[1]='Y';
      RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
      RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
      RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
      RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
      
      RECT.intercept[0]=E->mesh.layer0[3];
      RECT.mag[0]=0.0;
      RECT.intercept[1]=E->mesh.layer1[3];
      RECT.mag[1]=0.0;
   
      arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
      if(6==dofs) {
	arbitrary_bc_rectangle(E,&RECT,E->Vb[4],E->NODE,BC4,SBC4,0);
	arbitrary_bc_rectangle(E,&RECT,E->Vb[5],E->NODE,BC5,SBC5,0);
	arbitrary_bc_rectangle(E,&RECT,E->Vb[6],E->NODE,BC6,SBC6,0);
      }
    }

 /* Now read in from any file which is found (top level then inject) */
  read_bc_from_file(E,E->VB[1],E->node,"Velocity_x","Velx",BC1,SBC1);
  read_bc_from_file(E,E->VB[1],E->node,"Stress_x","Strx",SBC1,BC1);
  read_bc_from_file(E,E->VB[2],E->node,"Velocity_z","Velz",BC2,SBC2);
  read_bc_from_file(E,E->VB[2],E->node,"Stress_z","Strz",SBC2,BC2);
  if(dims==3 && 3==dofs) {
    read_bc_from_file(E,E->VB[3],E->node,"Velocity_y","Vely",BC3,SBC3);
    read_bc_from_file(E,E->VB[3],E->node,"Stress_y","Stry",SBC3,BC3);
  }
  else if(2==dims && 3==dofs) {
    read_bc_from_file(E,E->VB[3],E->node,"Rotation_y","Roty",BC3,SBC3);
  }
  else if(6==dofs) {
    read_bc_from_file(E,E->VB[4],E->node,"Rotation_x","Rotx",BC4,SBC4);
    read_bc_from_file(E,E->VB[5],E->node,"Rotation_z","Rotz",BC5,SBC5);
    read_bc_from_file(E,E->VB[6],E->node,"Rotation_y","Roty",BC6,SBC6);
  }

  fprintf(stderr,"RAA: done the part with reading any file which is found\n");

  /* AD HOC adjustment for funnel problem */
#if 0
  for(i=1;i<=E->mesh.noz;i++) {
      E->NODE[E->mesh.levmax][i] = (E->NODE[E->mesh.levmax][i] | BC2);
      E->Vb[2][E->mesh.levmax][i] = 0.0;
    }
#endif

   /* AD HOC adjustment for basal traction problem */

#if 0
  for(i=1;i<=E->mesh.nox;i++) {
      E->NODE[E->mesh.levmax][E->mesh.noz + (i-1) * E->mesh.noz] = (E->NODE[E->mesh.levmax][E->mesh.noz + (i-1) * E->mesh.noz] | BC1);
      E->Vb[1][E->mesh.levmax][E->mesh.noz + (i-1) * E->mesh.noz] = E->x[1][E->mesh.noz + (i-1) * E->mesh.noz];
    }
#endif


#if 0 /* AD HOC adjustment for trapdoors */
   for(i=1;i<=E->mesh.nno;i++) {
     if(E->x[1][i] <= 1.0 && E->x[2][i] >= 0.9) {
       E->Vb[1][E->mesh.levmax][i] =  E->Vb[2][E->mesh.levmax][i] = 0.0;
       E->node[i] |= (BC1 | BC2);
     }
     if(E->x[1][i] >= 2.0 && E->x[2][i] >= 0.9) {
       E->Vb[1][E->mesh.levmax][i] =  E->Vb[2][E->mesh.levmax][i] = 0.0;
       E->node[i] |= (BC1 | BC2);
     }
   
   }
#endif



  for(lv=E->mesh.levmax;lv>E->mesh.levmin;lv--) {
    inject_node_values(E,lv,E->Vb[1][lv],E->Vb[1][lv-1]);
    inject_node_values(E,lv,E->Vb[2][lv],E->Vb[2][lv-1]);
    if(3==dofs) {
      inject_node_values(E,lv,E->Vb[3][lv],E->Vb[3][lv-1]);
      inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
    }
    else if(dofs==6) {
      inject_node_values(E,lv,E->Vb[3][lv],E->Vb[3][lv-1]);
      inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
      inject_node_values(E,lv,E->Vb[4][lv],E->Vb[4][lv-1]);
      inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
      inject_node_values(E,lv,E->Vb[5][lv],E->Vb[5][lv-1]);
      inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
      inject_node_values(E,lv,E->Vb[6][lv],E->Vb[6][lv-1]);
      inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
    }
  }
  if(E->mesh.periodic_x || E->mesh.periodic_y)
    velocity_apply_periodic_bcs(E);

    /* Those were the defaults, now read in general cases from input file */
    fprintf(stderr,"RAA: start the part with reading the input-template bcs\n");

    arbitrary_bc_rectangle_file(E,&(E->mesh.Stxrectbcs),"Stress_x",E->Vb[1],E->NODE,SBC1,BC1,1);
    arbitrary_bc_circle_file(E,&(E->mesh.Stxcircbcs),"Stress_x",E->Vb[1],E->NODE,SBC1,BC1,1);
    arbitrary_bc_harmonic_file(E,&(E->mesh.Stxharmbcs),"Stress_x",E->Vb[1],E->NODE,SBC1,BC1,1);
    arbitrary_bc_polynomial_file(E,&(E->mesh.Stxpolybcs),"Stress_x",E->Vb[1],E->NODE,SBC1,BC1,1); 
 
    arbitrary_bc_rectangle_file(E,&(E->mesh.Stzrectbcs),"Stress_z",E->Vb[2],E->NODE,SBC2,BC2,1);
    arbitrary_bc_circle_file(E,&(E->mesh.Stzcircbcs),"Stress_z",E->Vb[2],E->NODE,SBC2,BC2,1);
    arbitrary_bc_harmonic_file(E,&(E->mesh.Stzharmbcs),"Stress_z",E->Vb[2],E->NODE,SBC2,BC2,1);
    arbitrary_bc_polynomial_file(E,&(E->mesh.Stzpolybcs),"Stress_z",E->Vb[2],E->NODE,SBC2,BC2,1); 
    
    arbitrary_bc_rectangle_file(E,&(E->mesh.Vxrectbcs),"Velocity_x",E->Vb[1],E->NODE,BC1,SBC1,0);
    arbitrary_bc_circle_file(E,&(E->mesh.Vxcircbcs),"Velocity_x",E->Vb[1],E->NODE,BC1,SBC1,0);
    arbitrary_bc_harmonic_file(E,&(E->mesh.Vxharmbcs),"Velocity_x",E->Vb[1],E->NODE,BC1,SBC1,0);
    arbitrary_bc_polynomial_file(E,&(E->mesh.Vxpolybcs),"Velocity_x",E->Vb[1],E->NODE,BC1,SBC1,0); 

    arbitrary_bc_rectangle_file(E,&(E->mesh.Vzrectbcs),"Velocity_z",E->Vb[2],E->NODE,BC2,SBC2,0);
    arbitrary_bc_circle_file(E,&(E->mesh.Vzcircbcs),"Velocity_z",E->Vb[2],E->NODE,BC2,SBC2,0);
    arbitrary_bc_harmonic_file(E,&(E->mesh.Vzharmbcs),"Velocity_z",E->Vb[2],E->NODE,BC2,SBC2,0);
    arbitrary_bc_polynomial_file(E,&(E->mesh.Vzpolybcs),"Velocity_z",E->Vb[2],E->NODE,BC2,SBC2,0); 

    if(2==dims && 3==dofs){
      arbitrary_bc_rectangle_file(E,&(E->mesh.Styrectbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_circle_file(E,&(E->mesh.Stycircbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Styharmbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Stypolybcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Rotation_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Rotation_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Rotation_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Rotation_y",E->Vb[3],E->NODE,BC3,SBC3,0); 
    }
    /* We have to implement new things in 3D, because stress_y has 2 meanings */
    else if(6==dofs) {
      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Rotation_x",E->Vb[4],E->NODE,BC4,SBC4,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Rotation_x",E->Vb[4],E->NODE,BC4,SBC4,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Rotation_x",E->Vb[4],E->NODE,BC4,SBC4,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Rotation_x",E->Vb[4],E->NODE,BC4,SBC4,0); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Rotation_z",E->Vb[5],E->NODE,BC5,SBC5,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Rotation_z",E->Vb[5],E->NODE,BC5,SBC5,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Rotation_z",E->Vb[5],E->NODE,BC5,SBC5,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Rotation_z",E->Vb[5],E->NODE,BC5,SBC5,0); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Rotation_y",E->Vb[6],E->NODE,BC6,SBC6,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Rotation_y",E->Vb[6],E->NODE,BC6,SBC6,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Rotation_y",E->Vb[6],E->NODE,BC6,SBC6,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Rotation_y",E->Vb[6],E->NODE,BC6,SBC6,0); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Styrectbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_circle_file(E,&(E->mesh.Stycircbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Styharmbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Stypolybcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1); 
    }
	                    
/*RAA: 30-12-02, added the stress bit for 3D stress_y */
    else if(3==dims && 3==dofs){ 
      arbitrary_bc_rectangle_file(E,&(E->mesh.Vyrectbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_circle_file(E,&(E->mesh.Vycircbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Vyharmbcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Vypolybcs),"Velocity_y",E->Vb[3],E->NODE,BC3,SBC3,0); 

      arbitrary_bc_rectangle_file(E,&(E->mesh.Styrectbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_circle_file(E,&(E->mesh.Stycircbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_harmonic_file(E,&(E->mesh.Styharmbcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1);
      arbitrary_bc_polynomial_file(E,&(E->mesh.Stypolybcs),"Stress_y",E->Vb[3],E->NODE,SBC3,BC3,1); 
    }

    fprintf(stderr,"\nRAA: finished the part concerning the general cases from the input template\n");
    fprintf(E->fp,"RAA: finished the section concerning the general cases from the input template\n");

    
    /*RAA - check*/
/*     for(i=1;i<=E->mesh.nno;i++) {
       if(3==dims)
         fprintf(E->fp1," nno Vbs1,2,3  %d   %g  %g  %g\n",i,E->Vb[1][0][i],E->Vb[2][0][i],E->Vb[3][0][i]);
       else
        fprintf(E->fp1," nno Vbs1,2  %d   %g  %g\n",i,E->Vb[1][0][i],E->Vb[2][0][i]);
     } */
      
    
    /* Skew boundary conditions */

#if 0
    E->control.HAVE_SKEWBCS = 1;
    for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
      for(i=1;i<E->mesh.NOZ[lv];i++) {
	if(E->X[lv][2][i] > 0.0 * 2.76 && E->x[2][i] < 4.0 + 0.0 * 3.375) {
	  /* E->NODE[lv][i] = (E->NODE[lv][i] | SKEWBC); */
	    E->NODE[lv][i] = (E->NODE[lv][i] | BC2); /**/
	}
      
	/*
	if(fabs(E->X[lv][2][i] - 2.75) < 1.0e-10)
	E->NODE[lv][i] = (E->NODE[lv][i] | BC2); */

	/* if(E->X[lv][2][i] == 3.375)
	   E->NODE[lv][i] = (E->NODE[lv][i] | BC2); */
      }

      /*  Not if rough */ /*
      E->NODE[lv][E->mesh.NOZ[lv]] = E->NODE[lv][E->mesh.NOZ[lv]] & ~BC2 ;
      E->NODE[lv][E->mesh.NOZ[lv]] = E->NODE[lv][E->mesh.NOZ[lv]] | SBC2;
      E->Vb[1][lv][E->mesh.NOZ[lv]] = 0.0;
      E->Vb[2][lv][E->mesh.NOZ[lv]] = 0.0;
      */
    }
    
    /* Allocate memory for any skew boundary conditions set during this 
       procedure */

    for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++)
      for(i=1;i<=E->mesh.NNO[lv];i++) {
	if(E->NODE[lv][i] & SKEWBC)
	  E->curvilinear.NODE_R[lv][i] = (higher_precision *) Malloc0(E->mesh.dof*E->mesh.dof*sizeof(higher_precision));
      }

    /* Compute rotation matrices for these boundary condition nodes */

    for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++)
      for(i=1;i<E->mesh.NOZ[lv];i++) {
	if(E->X[lv][2][i] > 2.75 && E->X[lv][2][i] < 3.375) {
	  store_node_R(E,i,M_PI/4,0.0,lv);
	}
	
	if(0 && E->X[lv][2][i] == 3.375) {
	  store_node_R(E,i,M_PI/8,0.0,lv);
	}
      }
#endif

    /* Here we check to see if any mesh-scaling-with-time
       type boundary conditions are applied. These are only
       valid if the mesh is rectangular ! 

       If these exist, ensure the velocity boundary conditions are compatible
       with the mesh movement implied here.
    */

/*RAA: here is the older version, pre-sequentially bidirectional moving boundaries
     -- lets keep this relatively unperturbed in case it must be resurrected*/
#if 0
    if(E->mesh.BCvelocityX0 != 0.0 || E->mesh.BCvelocity1X0 != 0.0 ) {
       RECT.numb=1; /* X-normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[1];
       RECT.mag[0]=E->mesh.BCvelocityX0 + E->mesh.BCvelocity1X0 * sin(E->mesh.BCvelocityphaseX0*M_PI*2.0) ;
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
       fprintf(stderr,"RAA: x0 moving boundary: %g\n",E->mesh.BCvelocityX0);
    }

    if(E->mesh.BCvelocityX1 != 0.0 || E->mesh.BCvelocity1X1 != 0.0) {
       RECT.numb=1; /* X-Normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[1];
       RECT.mag[0]=E->mesh.BCvelocityX1 + E->mesh.BCvelocity1X1 * sin(E->mesh.BCvelocityphaseX1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
       fprintf(stderr,"RAA: x1 moving boundary: %g\n",E->mesh.BCvelocityX1);
    }


    if(E->mesh.BCvelocityZ0 != 0.0 || E->mesh.BCvelocity1Z0 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[2];
       RECT.mag[0]=E->mesh.BCvelocityZ0 + E->mesh.BCvelocity1Z0 * sin(E->mesh.BCvelocityphaseZ0*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
       fprintf(stderr,"RAA: z0 moving boundary\n");
    }

    if(E->mesh.BCvelocityZ1 != 0.0 || E->mesh.BCvelocity1Z1 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[2];
       RECT.mag[0]=E->mesh.BCvelocityZ1 + E->mesh.BCvelocity1Z1 * sin(E->mesh.BCvelocityphaseZ1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
       fprintf(stderr,"RAA: z1 moving boundary\n");
    }


    if(3==dims && (E->mesh.BCvelocityY0 != 0.0 || E->mesh.BCvelocity1Y0 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[3];
       RECT.mag[0]=E->mesh.BCvelocityY0 + E->mesh.BCvelocity1Y0 * sin(E->mesh.BCvelocityphaseY0*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
       fprintf(stderr,"RAA: y0 moving boundary: %g\n",E->mesh.BCvelocityY0);
    }

    if(3==dims && (E->mesh.BCvelocityY1 != 0.0 || E->mesh.BCvelocity1Y1 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[3];
       RECT.mag[0]=E->mesh.BCvelocityY1 + E->mesh.BCvelocity1Y1 * sin(E->mesh.BCvelocityphaseY1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
       fprintf(stderr,"RAA: y1 moving boundary: %g\n",E->mesh.BCvelocityY1);
    }
#endif

    
/*-----------------------------------------------------*/
/*RAA: 19/02/02, here's the new part that checks what the time is before assuming
   that you want to apply a moving boundary from the start */
/*RAA: for the time interval option, keep it simple, ie. comment out sin( ) stuff
      and the constant &/or linear option extends to the new code in mesh_update()*/    

 if(!E->mesh.BCvelocity_time_interval) { /*RAA:	13/02/02, added this distinction*/
    if(E->mesh.BCvelocityX0 != 0.0 || E->mesh.BCvelocity1X0 != 0.0 ) {
       RECT.numb=1; /* X-normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[1];
       RECT.mag[0]=E->mesh.BCvelocityX0 +
	           E->mesh.BCvelocity1X0 * sin(E->mesh.BCvelocityphaseX0*M_PI*2.0) ;
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
       fprintf(stderr,"\tx0 moving boundary w/o interval: %g\n",E->mesh.BCvelocityX0);
    }

    if(E->mesh.BCvelocityX1 != 0.0 || E->mesh.BCvelocity1X1 != 0.0) {
       RECT.numb=1; /* X-Normal  */ 
       RECT.norm[0]='X';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[1];
       RECT.mag[0]=E->mesh.BCvelocityX1 +
	           E->mesh.BCvelocity1X1 * sin(E->mesh.BCvelocityphaseX1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
       fprintf(stderr,"\tx1 moving boundary w/o interval: %g\n",E->mesh.BCvelocityX1);
    }
 }
 else {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
	 /*assume the time when this is called is 0.0, as within the if statement below*/
    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((0.0 == E->mesh.BCvelocityX0_time[0][i][1]) && (E->mesh.BCvelocityX0_const[i][1] != 0.0)) {
         RECT.numb=1; /* X-normal  */ 
         RECT.norm[0]='X';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;

         RECT.intercept[0]=E->mesh.layer0[1];
         RECT.mag[0]=E->mesh.BCvelocityX0_const[i][1];
    
         arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
         fprintf(stderr,"\tx0 moving boundary w/ interval: %g\n",E->mesh.BCvelocityX0_const[i][1]);
	 break;
      }
    }

    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((0.0 == E->mesh.BCvelocityX1_time[0][i][1]) && (E->mesh.BCvelocityX1_const[i][1] != 0.0)) {
         RECT.numb=1; /* X-Normal  */ 
         RECT.norm[0]='X';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer1[1];
         RECT.mag[0]=E->mesh.BCvelocityX1_const[i][1];
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[1],E->NODE,BC1,SBC1,0);
         fprintf(stderr,"\tx1 moving boundary w/ interval: %g\n",E->mesh.BCvelocityX1_const[i][1]);
         /*RAA: check some stuff*/
         /*for(j=1;j<=E->mesh.nno;j++) {
          fprintf(stderr,"\tlevel 0: x-vel bc at node %d: %g\n",j,E->Vb[1][0][j]);
          }
         */	
	 break;
      }
    }
 }

 if(!E->mesh.BCvelocity_time_interval) { /*RAA: 13/02/02, added this distinction*/
    if(E->mesh.BCvelocityZ0 != 0.0 || E->mesh.BCvelocity1Z0 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[2];
       RECT.mag[0]=E->mesh.BCvelocityZ0 +
	           E->mesh.BCvelocity1Z0 * sin(E->mesh.BCvelocityphaseZ0*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
       fprintf(stderr,"/tz0 moving boundary w/o interval\n");
    }

    if(E->mesh.BCvelocityZ1 != 0.0 || E->mesh.BCvelocity1Z1 != 0.0) {
       RECT.numb=1; /* Z-normal  */ 
       RECT.norm[0]='Z';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[2];
       RECT.mag[0]=E->mesh.BCvelocityZ1 +
	           E->mesh.BCvelocity1Z1 * sin(E->mesh.BCvelocityphaseZ1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
       fprintf(stderr,"\tz1 moving boundary w/o interval \n");
    }
 }
 else {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
	 /*assume the time when this is called is 0.0, as within the if statement below*/
    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((0.0 == E->mesh.BCvelocityZ0_time[0][i][1]) && (E->mesh.BCvelocityZ0_const[i][1] != 0.0)) {
         RECT.numb=1; /* Z-normal  */ 
         RECT.norm[0]='Z';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer0[2];
         RECT.mag[0]=E->mesh.BCvelocityZ0_const[i][1];
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
         fprintf(stderr,"\tz0 moving boundary w/interval\n");
	 break;
      }
    }

    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((0.0 == E->mesh.BCvelocityZ1_time[0][i][1]) && (E->mesh.BCvelocityZ1_const[i][1] != 0.0)) {
         RECT.numb=1; /* Z-normal  */ 
         RECT.norm[0]='Z';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer1[2];
         RECT.mag[0]=E->mesh.BCvelocityZ1_const[i][1];
      
         arbitrary_bc_rectangle(E,&RECT,E->Vb[2],E->NODE,BC2,SBC2,0);
         fprintf(stderr,"\tz1 moving boundary w/ interval\n");
	 break;
      }
    }
 }

 if(!E->mesh.BCvelocity_time_interval) { /*RAA: 13/02/02, added this distinction*/
    if(3==dims && (E->mesh.BCvelocityY0 != 0.0 || E->mesh.BCvelocity1Y0 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer0[3];
       RECT.mag[0]=E->mesh.BCvelocityY0 +
	           E->mesh.BCvelocity1Y0 * sin(E->mesh.BCvelocityphaseY0*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
       fprintf(stderr,"\ty0 moving boundary w/o interval: %g\n",E->mesh.BCvelocityY0);
    }

    if(3==dims && (E->mesh.BCvelocityY1 != 0.0 || E->mesh.BCvelocity1Y1 != 0.0)) {
       RECT.numb=1; /* Y-normal  */ 
       RECT.norm[0]='Y';
       RECT.bb1[0]= -1.0e32;
       RECT.bb2[0]=  1.0e32;
       RECT.aa1[0]= -1.0e32;
       RECT.aa2[0]=  1.0e32;

       RECT.intercept[0]=E->mesh.layer1[3];
       RECT.mag[0]=E->mesh.BCvelocityY1 +
	           E->mesh.BCvelocity1Y1 * sin(E->mesh.BCvelocityphaseY1*M_PI*2.0);
    
       arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
       fprintf(stderr,"\ty1 moving boundary w/o interval: %g\n",E->mesh.BCvelocityY1);
    }
 }
 else {  /*RAA: when the setup has sequentially bi-directional moving boundaries, etc*/  
	 /*assume the time when this is called is 0.0, as within the if statement below*/
    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((3==dims && 0.0 == E->mesh.BCvelocityY0_time[0][i][1]) && (E->mesh.BCvelocityY0_const[i][1] != 0.0)) {
         RECT.numb=1; /* Y-normal  */ 
         RECT.norm[0]='Y';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer0[3];
         RECT.mag[0]=E->mesh.BCvelocityY0_const[i][1];

         arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
         fprintf(stderr,"\ty0 moving boundary w/ interval: %g\n",E->mesh.BCvelocityY0_const[i][1]);
	 break;
      }
    }

    for (i=0; i<3; i++) {  /*loading can start at any of the 3 sections of the trapezoid*/
      if((3==dims && 0.0 == E->mesh.BCvelocityY1_time[0][i][1]) && (E->mesh.BCvelocityY1_const[i][1] != 0.0)) {
         RECT.numb=1; /* Y-normal  */ 
         RECT.norm[0]='Y';
         RECT.bb1[0]= -1.0e32;
         RECT.bb2[0]=  1.0e32;
         RECT.aa1[0]= -1.0e32;
         RECT.aa2[0]=  1.0e32;
  
         RECT.intercept[0]=E->mesh.layer1[3];
         RECT.mag[0]=E->mesh.BCvelocityY1_const[i][1];
    
         arbitrary_bc_rectangle(E,&RECT,E->Vb[3],E->NODE,BC3,SBC3,0);
         fprintf(stderr,"\ty1 moving boundary w/ interval: %g\n",E->mesh.BCvelocityY1_const[i][1]);
	 break;
      }
    }
 }

/*RAA: pause program with infinite loop*/
/*   fprintf(stderr,"Pause program in velocity_boundary_conditions()/n");
   for(i=0; i<=1; j++) {
   }
*/

    return;
}

void velocity_apply_periodic_bcs(
				 struct All_variables *E
				 )
{
  int n1,n2,level;
  int i,j;
  
  fprintf(E->fp,"Periodic boundary conditions\n");
  
/*RAA: 19/10/01, part of line added below for && NOT periodic_y bcs */
  if(E->mesh.periodic_x && !E->mesh.periodic_y) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) { 
      for(i=1;i<=E->mesh.NOZ[level];i++)
	for(j=1;j<=E->mesh.NOY[level];j++) {
	  n1=i+(j-1)*E->mesh.NOX[level]*E->mesh.NOZ[level];
	  n2=n1+(E->mesh.NOX[level]-1)*E->mesh.NOZ[level];
	  
	  E->NODE[level][n1] = E->NODE[level][n1] & ~(BC1);
	  E->NODE[level][n2] = E->NODE[level][n2] & ~(BC1);
	  
	  if(i!=1 && i!=E->mesh.NOZ[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC2);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC2);
	  }
	  
	  if(j!=1 && j!=E->mesh.NOY[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC3);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC3);
	  }
	}
      E->mesh.NEQ[level] -= E->mesh.dof * E->mesh.NOY[level]*E->mesh.NOZ[level];
    }
  }

/*RAA: 23/5/01, code added below for periodic_y bcs (ie, back and front faces)*/
  else if(E->mesh.periodic_y && !E->mesh.periodic_x) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) { 
      for(i=1;i<=E->mesh.NOX[level];i++)
	for(j=1;j<=E->mesh.NOZ[level];j++) {
	  n1=j+(i-1)*E->mesh.NOZ[level];
	  n2=n1+(E->mesh.NOY[level]-1)*(E->mesh.NOZ[level])*(E->mesh.NOX[level]);
	  
	  E->NODE[level][n1] = E->NODE[level][n1] & ~(BC3);
	  E->NODE[level][n2] = E->NODE[level][n2] & ~(BC3);
	  
	  if(i!=1 && i!=E->mesh.NOX[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC1);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC1);
	  }
	  
	  if(j!=1 && j!=E->mesh.NOZ[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC2);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC2);
	  }
	}
      E->mesh.NEQ[level] -= E->mesh.dof * E->mesh.NOZ[level]*E->mesh.NOX[level];
    }
  }
  
/*RAA: 19/10/01, added this part, may not be "smooth" but should be ok*/
  else if(E->mesh.periodic_y && E->mesh.periodic_x) {
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) { 

      /*do the periodic_y part first*/
      for(i=1;i<=E->mesh.NOX[level];i++)
	for(j=1;j<=E->mesh.NOZ[level];j++) {
	  n1=j+(i-1)*E->mesh.NOZ[level];
	  n2=n1+(E->mesh.NOY[level]-1)*(E->mesh.NOZ[level])*(E->mesh.NOX[level]);
	  
	  E->NODE[level][n1] = E->NODE[level][n1] & ~(BC3);
	  E->NODE[level][n2] = E->NODE[level][n2] & ~(BC3);
	  
	  if(i!=1 && i!=E->mesh.NOX[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC1);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC1);
	  }
	  
	  if(j!=1 && j!=E->mesh.NOZ[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC2);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC2);
	  }
	}

      /*now do the periodic_x part*/
      for(i=1;i<=E->mesh.NOZ[level];i++)
	for(j=1;j<=E->mesh.NOY[level];j++) {
	  n1=i+(j-1)*E->mesh.NOX[level]*E->mesh.NOZ[level];
	  n2=n1+(E->mesh.NOX[level]-1)*E->mesh.NOZ[level];
	  
	  E->NODE[level][n1] = E->NODE[level][n1] & ~(BC1);
	  E->NODE[level][n2] = E->NODE[level][n2] & ~(BC1);
	  
	  if(i!=1 && i!=E->mesh.NOZ[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC2);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC2);
	  }
	  
	  if(j!=1 && j!=E->mesh.NOY[level]) {
	    E->NODE[level][n1] = E->NODE[level][n1] & ~(BC3);
	    E->NODE[level][n2] = E->NODE[level][n2] & ~(BC3);
	  }
	}
   /*now update the # of eqtns making sure the front right edge isn't counted twice*/
      E->mesh.NEQ[level] -= (E->mesh.dof * E->mesh.NOZ[level])*(E->mesh.NOY[level] + E->mesh.NOX[level] - 1);
    }
  }


  return;
}
