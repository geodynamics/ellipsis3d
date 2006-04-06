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
/* #include <string.h> */

void temperature_boundary_conditions(
				     struct All_variables *E
				     )
{
  void temperatures_conform_bcs();
  void temperature_apply_periodic_bcs(); 
  void arbitrary_bc_rectangle();
  void arbitrary_bc_rectangle_file();
  void arbitrary_bc_circle_file();
  void arbitrary_bc_harmonic_file();
  void arbitrary_bc_polynomial_file();
  void field_arbitrary_rectangle_file();
  void field_arbitrary_circle_file();
  void field_arbitrary_harmonic_file();

  int i,lv;

  struct RectBc RECT;

  /* Default: Isothermal top (toptbcval) and bottom (bottbcval), reflecting sidewalls. */

  RECT.numb=2; /* X-normal */ 
  RECT.norm[0]=RECT.norm[1]='X';
  RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
  RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
  RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
  RECT.aa2[0]=RECT.aa2[1]=  1.0e32;

  RECT.intercept[0]=E->mesh.layer0[1];
  RECT.mag[0]=0.0;
  RECT.intercept[1]=E->mesh.layer1[1];
  RECT.mag[1]=0.0;
  arbitrary_bc_rectangle(E,&RECT,E->Tb,E->NODE,FBX, TBD | FBZ | FBY, 0); /* Sides */
  
  if(3==E->mesh.nsd) {
      RECT.numb=2; /* Y-normal */ 
      RECT.norm[0]=RECT.norm[1]='Y';
      RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
      RECT.bb2[0]=RECT.bb2[1]=  1.0e32;
      RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
      RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
      
      RECT.intercept[0]=E->mesh.layer0[3];
      RECT.mag[0]=0.0;
      RECT.intercept[1]=E->mesh.layer1[3];
      RECT.mag[1]=0.0;
      arbitrary_bc_rectangle(E,&RECT,E->Tb,E->NODE,FBY, TBD | FBX | FBZ, 0); /* Sides */
  }
    
  RECT.numb=2;  /* Z-normal */
  RECT.norm[0]=RECT.norm[1]='Z';
  RECT.aa1[0]=RECT.aa1[1]= -1.0e32;
  RECT.aa2[0]=RECT.aa2[1]=  1.0e32;
  RECT.bb1[0]=RECT.bb1[1]= -1.0e32;
  RECT.bb2[0]=RECT.bb2[1]=  1.0e32;

  RECT.intercept[0]=E->mesh.layer0[2];
  RECT.mag[0]=E->control.TBCtopval;
  RECT.intercept[1]=E->mesh.layer1[2];
  RECT.mag[1]=E->control.TBCbotval;
  arbitrary_bc_rectangle(E,&RECT,E->Tb,E->NODE,TBD,FBZ | FBX | FBY, 0); 

  /* Now read in from any file which is found (top level). */

  read_bc_from_file(E,E->TB,E->node,"Temperature","Temp",TBD,FBZ | FBX | FBY);
  read_bc_from_file(E,E->TB,E->node,"Heat_flux_z","Hflz",FBZ,TBD | FBX | FBY);

  for(lv=E->mesh.levmax;lv>E->mesh.levmin;lv--) {
    inject_node_values(E,lv,E->Tb[lv],E->Tb[lv-1]);
    inject_node_int_values(E,lv,E->NODE[lv],E->NODE[lv-1]);
  }
  
  /* Read in arbitrary structures from parameter files */
  
  arbitrary_bc_rectangle_file(E,&(E->temperature.Trectbcs),"Heat_flux_z",E->Tb,E->NODE,FBZ,TBD,0);
  arbitrary_bc_circle_file(E,&(E->temperature.Tcircbcs),"Heat_flux_z",E->Tb,E->NODE,FBZ,TBD);
  arbitrary_bc_harmonic_file(E,&(E->temperature.Tharmbcs),"Heat_flux_z",E->Tb,E->NODE,FBZ,TBD);
  arbitrary_bc_polynomial_file(E,&(E->temperature.Tpolybcs),"Heat_flux_z",E->Tb,E->NODE,FBZ,TBD); 

  field_arbitrary_rectangle_file(E,1,&(E->temperature.Trects),"Temp_fixed",E->TB,E->node,TBD,FBZ|FBX|FBY,E->mesh.levmax);
  field_arbitrary_circle_file(E,1,&(E->temperature.Tcircs),"Temp_fixed",E->TB,E->node,TBD,FBZ|FBX|FBY,E->mesh.levmax);
  field_arbitrary_harmonic_file(E,1,&(E->temperature.Tharms),"Temp_fixed",E->TB,E->node,TBD,FBZ|FBX|FBY,E->mesh.levmax);

  for(lv=E->mesh.levmin;lv<E->mesh.levmax;lv++) {
    field_arbitrary_rectangle(E,&(E->temperature.Trects),E->Tb[lv],E->NODE[lv],TBD,FBZ|FBX|FBY,lv);
    field_arbitrary_circle(E,&(E->temperature.Tcircs),E->Tb[lv],E->NODE[lv],TBD,FBZ|FBX|FBY,lv);
    field_arbitrary_harmonic(E,&(E->temperature.Tharms),E->Tb[lv],E->NODE[lv],TBD,FBZ|FBX|FBY,lv);
  }

  arbitrary_bc_rectangle_file(E,&(E->temperature.Trectbcs),"Temp",E->Tb,E->NODE,TBD,FBZ|FBX|FBY,0);
  arbitrary_bc_circle_file(E,&(E->temperature.Tcircbcs),"Temp",E->Tb,E->NODE,TBD,FBZ|FBX|FBY);
  arbitrary_bc_harmonic_file(E,&(E->temperature.Tharmbcs),"Temp",E->Tb,E->NODE,TBD,FBZ|FBX|FBY);
  arbitrary_bc_polynomial_file(E,&(E->temperature.Tpolybcs),"Temp",E->Tb,E->NODE,TBD,FBZ|FBX|FBY); 
 
  if(E->mesh.periodic_x || E->mesh.periodic_y)
    temperature_apply_periodic_bcs(E); 

  field_arbitrary_rectangle_file(E,1,&(E->temperature.Tintz_off),"Tintz_off",(standard_precision *)NULL,E->node,0,INTZ,E->mesh.levmax);
  field_arbitrary_rectangle_file(E,1,&(E->temperature.Tintz_on) ,"Tintz_on" ,(standard_precision *)NULL,E->node,INTZ,0,E->mesh.levmax);

  for(lv=E->mesh.levmin;lv<E->mesh.levmax;lv++) {
    field_arbitrary_rectangle(E,&(E->temperature.Tintz_off),(standard_precision *)NULL,E->NODE[lv],0,INTZ,lv);
    field_arbitrary_rectangle(E,&(E->temperature.Tintz_on) ,(standard_precision *)NULL,E->NODE[lv],INTZ,0,lv);
  }

  temperatures_conform_bcs(E,E->T);
  return;
}

void temperature_apply_periodic_bcs(
				    struct All_variables *E
				    )
{
    int n1,n2,e1,level;
    int i,j;

    void construct_lm();
   
    fprintf(E->fp,"Periodic temperature boundary conditions\n");
   
    if(E->mesh.periodic_x) {
	level = E->mesh.levmax;
	for(i=1;i<=E->mesh.NOZ[level];i++)
	    for(j=1;j<=E->mesh.NOY[level];j++) {
		n1=i+(j-1)*E->mesh.NOX[level]*E->mesh.NOZ[level];
		n2=n1+(E->mesh.NOX[level]-1)*E->mesh.NOZ[level];

		/*RAA: why were these two lines commented out?*/
		/*E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
		  E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD); 
	      
		if((i!=1) && (i!=E->mesh.NOZ[level])) {
		    E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
		    E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD);
		}

		if((3==E->mesh.nsd) && (j!=1 && j!=E->mesh.NOY[level])) {
		    E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
		    E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD);
		}
		 */ /*RAA, I comment out the rest based on some analysis
		         that showed some problems with temp at periodic 
			 nodes - but this should be verified at some point*/
	    }
    }
    /*RAA: 23/5/01 - code added below for periodic bcs in the y-direction*/   
    if(E->mesh.periodic_y) {
      level = E->mesh.levmax;
      for(i=1;i<=E->mesh.NOX[level];i++)
         for(j=1;j<=E->mesh.NOZ[level];j++) {
             n1=j+(i-1)*E->mesh.NOZ[level];
             n2=n1+(E->mesh.NOY[level]-1)*E->mesh.NOZ[level]*E->mesh.NOX[level];
      
     /*      E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
             E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD);  

             if((i!=1) && (i!=E->mesh.NOX[level])) {
                 E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
                 E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD);
             }

             if((3==E->mesh.nsd) && (j!=1 && j!=E->mesh.NOZ[level])) {
                 E->NODE[level][n1] = E->NODE[level][n1] & ~(TBD);
                 E->NODE[level][n2] = E->NODE[level][n2] & ~(TBD);
             }
     */
	}
    }

    return;
}

