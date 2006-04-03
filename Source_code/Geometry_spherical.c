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



#include "element_definitions.h"
#include "global_defs.h"

void set_cylinder_defaults(
     struct All_variables *E
)
{ 
  E->control.model = 1 ;
  E->control.CYLINDER = 1;
  E->mesh.nsd = 2;
  E->mesh.dof = 2;  
  E->control.ORTHO = 0;
  E->control.ORTHOZ = 0;

}
void set_cylindercoss_defaults(
     struct All_variables *E
)
{ 
  E->control.model = 2 ;
  E->control.CYLINDER = 1;
  E->mesh.nsd = 2;
  E->mesh.dof = 3;
  E->control.ORTHO = 0;
  E->control.ORTHOZ = 0;

}

void set_sphere_defaults(
     struct All_variables *E
)
{ 
  E->control.model = 1 ;
  E->control.SPHERE = 1;
  E->mesh.nsd = 3;
  E->mesh.dof = 3;
  E->control.ORTHO = 0;
  E->control.ORTHOZ = 0;
}
void set_spherecoss_defaults(
     struct All_variables *E
)
{ 
  E->control.model = 2 ;
  E->control.SPHERE = 1;
  E->mesh.nsd = 3;
  E->mesh.dof = 6;
  E->control.ORTHO = 0;
  E->control.ORTHOZ = 0;
}
