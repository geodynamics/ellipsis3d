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

standard_precision Zj0[1000],Zj1[1000];

void set_axi_defaults(
     struct All_variables *E
)
{ 
#include  "Zeros_of_J0_and_J1.h"
  E->control.model = 1 ;
  E->control.AXI = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 2;
}

void set_axicoss_defaults(
     struct All_variables *E
)
{ 
#include  "Zeros_of_J0_and_J1.h"
  E->control.model = 2 ;
  E->control.AXI = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 3;
}

void hankel_transform(
     struct All_variables *E,
     higher_precision data0[],
     higher_precision transform[]
)
{
    int i,j,points ;
    higher_precision k,integrand,data1,data2,norm,horizontal_1scale;
    standard_precision data[2000],Data[2000];
    standard_precision return_layer_value();

    points = E->mesh.nnx[1];
    horizontal_1scale = 1.0 / (E->x[1][E->mesh.nno]);
  
    for(j=1;j<=points;j++)
	data[j] = data0[j]; 

    norm=(higher_precision)return_layer_value(E,data,1,1);
 
    for(j=1;j<=points;j++)
	data[j] -= norm;	/*	Improves the accuracy of the integration to 
						use the axial value as the reference level.	*/
    for(i=0;i<points;i++)  {
	k=Zj1[i]*horizontal_1scale;	/*  Wavenumbers for the problem   */
	norm =0.5* j0(Zj1[i]*horizontal_1scale)*j0(Zj1[i]*horizontal_1scale);

	for(j=1;j<=points;j++) 
	    Data[j] = j0(k*E->x[1][1+(j-1)*E->mesh.noz])*data[j];

	integrand=(higher_precision)return_layer_value(E,Data,1,0);
	for(j=1;j<=points;j++) 
	    Data[j] = j0(k*E->x[1][1+(j-1)*E->mesh.noz])*j0(k*E->x[1][1+(j-1)*E->mesh.noz]);
	norm=(higher_precision)return_layer_value(E,Data,1,0);

	transform[i+1] = integrand/norm;
	}	
    transform[1] = 0.0;	
    return;	
}

void inverse_hankel_transform(
     struct All_variables *E,
     higher_precision data[],
     higher_precision transform[]
)
{
    int i,j;
    higher_precision x,integrand,data1,data2,norm1,norm2;
    higher_precision horizontal_1scale;
    int points;

    points = E->mesh.nnx[1];
    horizontal_1scale = 1.0/E->x[1][E->mesh.nno];

    for(i=0;i<points;i++)    {
	x = E->x[1][1+(i-1)*E->mesh.noz];	
	integrand = 0.0;
	for(j=1;j<points;j++) {
	    norm1 =  0.5*j0(Zj1[j-1]*horizontal_1scale)*j0(Zj1[j-1]*horizontal_1scale);
	    norm2 =  0.5*j0(Zj1[j]*horizontal_1scale)*j0(Zj1[j]*horizontal_1scale);
	    data1=norm1*j0(x*Zj1[j-1]*horizontal_1scale)*data[j];
	    data2=norm2*j0(x*Zj1[j]*horizontal_1scale)*data[j+1];
	    integrand += (data1 + data2) * 0.25 * (Zj1[j]*Zj1[j]-Zj1[j-1]*Zj1[j-1])*
		(horizontal_1scale*horizontal_1scale);		}
	transform[i+1] = integrand;	}
    return;	
}
