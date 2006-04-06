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



struct ADVECTION {
   int ADVECTION;

  standard_precision gamma;
  standard_precision timestep;
  standard_precision elastic_timestep;


  standard_precision timestep_diff;
  standard_precision timestep_adv;
  int diff_ratio;


   standard_precision previous_timestep;
   standard_precision previous_timestep_2;
   standard_precision fine_tune_dt;
   standard_precision fixed_timestep;
   standard_precision max_elapsed_time;

   int min_timesteps;  
   int max_timesteps;
   int max_total_timesteps;
   int timesteps;
   int total_timesteps;
   int temp_iterations;
   int sub_iterations;
   int last_sub_iterations; 
    
  /* PIC MG advection/diffusion information */

  int pic_nodes;
  int *pic_elt_number;

  struct TRACER_ELT_WEIGHT  *pic_sfn_values;

  standard_precision *pic_ls_diag;

  standard_precision *pic_x;
  standard_precision *pic_z;
  standard_precision *pic_y;

  standard_precision *pic_x0;
  standard_precision *pic_z0;
  standard_precision *pic_y0;

  /* Temperature ... any other diffusive fields
     will also have to be computed at nodes and
     the advection term added here (including momentum ?) */

  standard_precision *pic_T;

  int *elements_tracers;
  int *et_start;
  int *et_length;
  
} advection;

