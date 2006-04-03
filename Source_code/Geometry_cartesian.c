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



#include "global_defs.h"


void set_2dc_defaults(
		      struct All_variables *E
		      )
{
  E->control.model = 1 ;
  E->control.CART2D = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 2;
}
void set_2dccoss_defaults(
			  struct All_variables *E
			  )
{
  E->control.model = 2 ;
  E->control.CART2D = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 3;
}

void set_2pt5dc_defaults(
			 struct All_variables *E
			 )
{ 
  E->control.model = 1 ;
  E->control.CART2pt5D = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 3;
}
void set_2pt5dccoss_defaults(
			     struct All_variables *E
			     )
{ 
  E->control.model = 2 ;
  E->control.CART2pt5D = 1; 
  E->mesh.nsd = 2;
  E->mesh.dof = 3;
}

void set_3dc_defaults(
		      struct All_variables *E
		      )
{ 
  E->control.model = 1 ;
  E->control.CART3D = 1;
  E->mesh.nsd = 3;
  E->mesh.dof = 3;
}
void set_3dccoss_defaults(
			  struct All_variables *E
			  )
{ 
  E->control.model = 2 ;
  E->control.CART3D = 1;
  E->mesh.nsd = 3;
  E->mesh.dof = 6;
}
