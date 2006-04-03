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
#include <math.h>


void set_mg_defaults(
     struct All_variables *E
)
{ 
    void solve_constrained_flow_iterative();
    void mg_allocate_vars();

    E->solver_allocate_vars = mg_allocate_vars;
    E->solve_stokes_problem = solve_constrained_flow_iterative;
    
    return;
}

void mg_allocate_vars(
     struct All_variables *E
)
{  
    /* nothing specific at the moment */
    return;
}
