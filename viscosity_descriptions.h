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


            


/* in this file define the contents of the VISC_OPT data structure
   which is used to store information used to create predefined 
   viscosity fields, those determined from prior input, those
   related to temperature/pressure/stress/anything else. */


struct VISC_OPT {
  void (* update_viscosity)();
  
  int SMOOTH;
  int smoothing[MAX_LEVELS];
  
  char STRUCTURE[20];		/* which option to determine viscosity field, one of .... */
  int FROM_SYSTEM;
  int FROM_FILE;
  int FROM_SPECS;

  /* GLOBALLY applicable options  */
  
  int MAX;
  standard_precision max_value;
  int MIN;
  standard_precision min_value;

  /* SYSTEM STATE VISCOSITY PARAMETERS */

  int SDEPV;
  int TDEPV; 
  int CHEMDEPV;
  int H2ODEPV;
  int YIELD;
  int GRAINSIZE;
    			
  int rheologies;

  int RHEOL_T_type[40];                /* 1,2 */
  int DEPL_T_type[40]; /*RAA: 10/07/03, O'Neill melting stuff */
  int MELT_model[40];  /*RAA: 10/07/03, O'Neill melting stuff */


  standard_precision sdepv_max_effect[40]; /* Stress Dependence */ 
  standard_precision sdepv_expt[40];
 
  standard_precision N0[40];
  standard_precision E[40],T0[40];
  standard_precision T[40],Z[40];
  standard_precision Tmin_value[40];  /* Cutoff temperatures */
  standard_precision Tmax_value[40];

  standard_precision Zrange_min[40];  /* Outside these ranges visc -> inf s.t. it does not contribute to composite */
  standard_precision Zrange_max[40];
  standard_precision Trange_min[40];  
  standard_precision Trange_max[40];

  int Phase[40];

  standard_precision yield_stress_B0;
  standard_precision yield_stress_Bz;
  standard_precision yield_stress_Bp;
  standard_precision yield_stress_Bc;
  standard_precision yield_stress_minimum;
  standard_precision yield_stress_maximum;
  standard_precision yield_hysteresis;
  standard_precision sp0,sp1,sp2,sp3; /*RAA: 10/07/03, O'Neill melting stuff */
  standard_precision lp0,lp1,lp2,lp3; /*RAA: 10/07/03, O'Neill melting stuff */

  standard_precision yield_strain_weaken;
  standard_precision yield_strain_onset;
  standard_precision yield_strain_saturate;
  
  /* MODULE BASED VISCOSITY VARIATIONS */
  
  
  standard_precision CHEM_N0[40];
  standard_precision CHEM_E[40];
  standard_precision CHEM_T0[40];
  standard_precision CHEM_T[40];
  standard_precision CHEM_Z[40];
  standard_precision CHEM_T_min_value[40];
  standard_precision CHEM_T_max_value[40];
  standard_precision CHEM_delta_eta[40];
  standard_precision CHEM_yield_stress_B0; 
  standard_precision CHEM_yield_stress_Bz; 
  standard_precision CHEM_yield_stress_Bp; 
 
  standard_precision H2Oeta0[40];
  standard_precision H2OetaN[40];
  standard_precision H_2O_yield_stress_weaken; 
 
  /* Viscosities change with range of phase function, GAMMA */

  standard_precision pt_G_range_min[40];  
  standard_precision pt_G_range_max[40];
  
  char old_file[100];
  /* Specification info */
  
  /* Prespecified viscosity parameters */
  char VISC_OPT[20];
  
  int layers;			/* number of layers with properties .... */
  standard_precision layer_depth[40];
  standard_precision layer_visc[40];
  
  int COSX;
  standard_precision cosx_epsilon;
  standard_precision cosx_k;
  int cosx_exp;
  
  int EXPX;
  standard_precision expx_epsilon;

  struct Rect Strnrects;



  
} viscosity;
