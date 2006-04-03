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



struct CONVECTION { /* information controlling convection problems */   
  char old_T_file[100];

  int number_of_perturbations;
  standard_precision perturb_mag[33];
  standard_precision perturb_k[33];
  standard_precision perturb_ky[33];

  standard_precision elasticity1;

  struct SOURCES {
     int number;
     standard_precision t_offset;
     standard_precision Q[10];
     standard_precision lambda[10]; 
  }  heat_sources;

} convection;
