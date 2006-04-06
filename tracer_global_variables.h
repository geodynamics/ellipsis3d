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


struct TRACERS {

  /* Dipstick particles - utterly passive */

  int SAMPLE_PTS;
  standard_precision *sampled_data[MAX_SAMPLE_PTS+3];
  standard_precision sample_x[MAX_SAMPLE_PTS];   
  standard_precision sample_y[MAX_SAMPLE_PTS];
  standard_precision sample_z[MAX_SAMPLE_PTS];
  standard_precision sample_x1[MAX_SAMPLE_PTS];   
  standard_precision sample_y1[MAX_SAMPLE_PTS];
  standard_precision sample_z1[MAX_SAMPLE_PTS];
  standard_precision sample_value[MAX_SAMPLE_PTS];

  int sample_type[MAX_SAMPLE_PTS];
  int sample_direction[MAX_SAMPLE_PTS];
  int sample_in_element[MAX_SAMPLE_PTS];
  int sample_in_plot[MAX_SAMPLE_PTS];
  int sample_lagrangian[MAX_SAMPLE_PTS];
  int sample_normalize[MAX_SAMPLE_PTS];

  standard_precision sample_plotmin[MAX_SAMPLE_PTS];
  standard_precision sample_plotmax[MAX_SAMPLE_PTS];
  standard_precision sample_Red[MAX_SAMPLE_PTS];  /* color per material and per PPM file type */
  standard_precision sample_Green[MAX_SAMPLE_PTS];
  standard_precision sample_Blue[MAX_SAMPLE_PTS];

  /* Active particles - carry history, integration weights and whatnot */

  standard_precision *tx,*tx1;
  standard_precision *tz,*tz1;
  standard_precision *ty,*ty1;

  standard_precision *eta1[MAX_LEVELS],*eta2[MAX_LEVELS],*eta3[MAX_LEVELS];

  standard_precision *tracer_weight_fn[MAX_LEVELS];
  standard_precision *tracer_jacobian[MAX_LEVELS];
  standard_precision *tracer_std_weighting[MAX_LEVELS];

  int *tr_in_element[MAX_LEVELS];
  int *tr_in_element_offset[MAX_LEVELS];
  int *tr_in_element_number[MAX_LEVELS];

  struct TRACER_ELT_WEIGHT {
    standard_precision node[28];  /* RAA: 29/3/01, 28 for 3D, was 10 */
  }  *sfn_values[MAX_LEVELS];

  standard_precision *ls_diag[MAX_LEVELS];

  int *property_group;
  int *tracer_elt[MAX_LEVELS];

  int NUM_MATERIALS;
  int NUM_TRACERS;
  int MEM_TRACERS;
  int TRACERS;

  char tracer_input_file[100];

  /* Physical properties ascribed to tracers */

  standard_precision *T ;
  standard_precision *Tdot ;
  standard_precision *Q;
  standard_precision *Visc ;
  standard_precision *Visc0;
  standard_precision *BulkVisc ;
  standard_precision *Density_change ;

  /* Derived quantities at the tracers */

  standard_precision *weight;
  standard_precision *edot;
  standard_precision *edot_integrated;
  standard_precision *edotp_integrated;
  standard_precision *strs;
  standard_precision *strd;
  standard_precision *strd1;
  standard_precision *strsum;
  standard_precision *strdiff;
  /*RAA: 24/09/02, C. O'Neill melting stuff */
  standard_precision *depl; 
  standard_precision *dFdot;
  int melting[MAX_MATERIALS]; /*RAA: 18/10/02, my flag 4 melting on/off.  C. O'Neill melting stuff */

  standard_precision *Hvisc;

  standard_precision *dX11,*dX12,*dX13;
  standard_precision *dX21,*dX22,*dX23;
  standard_precision *dX31,*dX32,*dX33;

  standard_precision *S11,*S12,*S13; /* Stress terms */
  standard_precision *S21,*S22,*S23;
  standard_precision *S31,*S32,*S33;
  standard_precision *M31,*M32;
  standard_precision *M12,*M13;
  standard_precision *M21,*M23;
  standard_precision *Pt;

  standard_precision *DQ;
  standard_precision *DQ1;

  standard_precision *A1;   /* inertial terms */
  standard_precision *A2;
  standard_precision *A3;

  standard_precision *U1;   /* inertial terms */
  standard_precision *U2;
  standard_precision *U3;
 
  standard_precision *UU1;   /* inertial terms */
  standard_precision *UU2;
  standard_precision *UU3;

  standard_precision *tr_dim;

  standard_precision  *yielded;

  standard_precision *n1,*n2,*n3;        /* Components of director vector for orthotropy */ 



  struct VISCOSITY {
    int rheologies;

    int RHEOL_T_type[MAX_RHEOL_CPTS];
    int DEPL_T_type[MAX_RHEOL_CPTS];  /*RAA: 24/09/02, C. O'Neill melting stuff */
    int MELT_model[MAX_RHEOL_CPTS];  /*RAA: 10/07/03, C. O'Neill melting stuff */
    standard_precision sdepv_expt[MAX_RHEOL_CPTS];
    standard_precision gr_size_expt[MAX_RHEOL_CPTS];
    
    standard_precision N0[MAX_RHEOL_CPTS];
    standard_precision E[MAX_RHEOL_CPTS];
    standard_precision T0[MAX_RHEOL_CPTS];          /* though this ref temperature should apply everywhere */
    standard_precision T[MAX_RHEOL_CPTS],Z[MAX_RHEOL_CPTS];
    standard_precision Tmin_value[MAX_RHEOL_CPTS];  /* Cutoff temperatures */
    standard_precision Tmax_value[MAX_RHEOL_CPTS];

    standard_precision Bulk_visc_ratio;            /* cf Augmented Lagrangian */
    standard_precision Pen_bulk;            /* cf Augmented Lagrangian */
    standard_precision Elas_shear_mod;            /* Elastic shear modulus (constant) */
    standard_precision Elas_Bulk_mod;         /* Fixed ratio for all viscosity laws (for now) */
    standard_precision Num_Elas_Bulk_mod;         /* For numerical viewpoint */

 
    /* Orthotropic materials */

    standard_precision Ortho_viscosity_ratio;     /* Viscosity component applying to "layers" */
 

    int phase[MAX_RHEOL_CPTS];
    
    standard_precision Trange_min[MAX_RHEOL_CPTS];  
    standard_precision Trange_max[MAX_RHEOL_CPTS];

    standard_precision yield_stress_B0;
    standard_precision yield_stress_Bz;
    standard_precision yield_stress_Bp;
    standard_precision yield_stress_Bc;

    standard_precision yield_stress_Edot0;
    standard_precision yield_stress_Edotn;
    standard_precision yield_stress_Edota;

    standard_precision yield_stress_E0;
    standard_precision yield_stress_En;
    standard_precision yield_stress_Ea;
    standard_precision yield_stress_E0dt;

    standard_precision yield_stress_ET;
    standard_precision yield_stress_minimum;
    standard_precision yield_stress_maximum;



#if defined (CHEM_TRANS_MODULE_INSTALLED)
  
    standard_precision mobile_phase_ratio;

#endif



    

  } visc[MAX_MATERIALS];

  struct COSSERAT {
    standard_precision A11,A12,A13,A22,A23,A33 ;
    standard_precision G11,G12,G13,G22,G23,G33 ;
    standard_precision B12,B21,B13,B31,B23,B32 ;
  } coss[MAX_MATERIALS];

  int grain_size_model;  /* which formulation to use for grain growth model ... */



  struct GRAIN_GROWTH {
       
    /* empirical formulation based on piezometric value */

    standard_precision grsz_B;
    standard_precision grsz_epsT;
    standard_precision grsz_stsexp;


    /* theoretical formulation based on nucleation / growth
       balances */

    standard_precision reduction_factor_equil;
    standard_precision reduction_factor_metas;
 
    int GRAIN_T_dep[MAX_MATERIAL_PHASES];
    standard_precision ggrw_a[MAX_MATERIAL_PHASES];
    standard_precision ggrw_m[MAX_MATERIAL_PHASES];
    standard_precision ggrw_Q[MAX_MATERIAL_PHASES];
    standard_precision ggrw_T0[MAX_MATERIAL_PHASES]; /* though this ref temp should apply everywhere */
    standard_precision ggrw_T[MAX_MATERIAL_PHASES]; 
    standard_precision gnuc_a[MAX_MATERIAL_PHASES];
    standard_precision gnuc_m[MAX_MATERIAL_PHASES];
    standard_precision gnuc_meps[MAX_MATERIAL_PHASES];
   
  } grain[MAX_MATERIALS];

  /* Other material properties */

  standard_precision Therm_exp[MAX_MATERIALS];
  standard_precision Therm_diff[MAX_MATERIALS];
  standard_precision Cp[MAX_MATERIALS];
  standard_precision Qt[MAX_MATERIALS];
  standard_precision Hv0[MAX_MATERIALS];

  standard_precision Depl_exp[MAX_MATERIALS];   /*RAA: 24/09/02,  O'Neill melting stuff */
  standard_precision Latent_heat[MAX_MATERIALS];   /*RAA: 10/07/03, O'Neill melting stuff */
  standard_precision sp0[MAX_MATERIALS];   /*RAA: 10/07/03, O'Neill melting stuff */
  standard_precision sp1[MAX_MATERIALS]; 
  standard_precision sp2[MAX_MATERIALS];
  standard_precision sp3[MAX_MATERIALS];
  standard_precision lp0[MAX_MATERIALS];
  standard_precision lp1[MAX_MATERIALS];
  standard_precision lp2[MAX_MATERIALS];
  standard_precision lp3[MAX_MATERIALS];

  standard_precision Red[MAX_MATERIALS][40];  /* color per material and per PPM file type */
  standard_precision Green[MAX_MATERIALS][40];
  standard_precision Blue[MAX_MATERIALS][40];
  standard_precision Opacity[MAX_MATERIALS][40];

  standard_precision Red_hot[MAX_MATERIALS][40];
  standard_precision Green_hot[MAX_MATERIALS][40];
  standard_precision Blue_hot[MAX_MATERIALS][40];
  standard_precision Opacity_hot[MAX_MATERIALS][40];

  standard_precision Red_strained[MAX_MATERIALS][40];
  standard_precision Green_strained[MAX_MATERIALS][40];
  standard_precision Blue_strained[MAX_MATERIALS][40];
  standard_precision Opacity_strained[MAX_MATERIALS][40];


  /* Phase transitions */

  int Phases[MAX_MATERIALS];
  standard_precision Density[MAX_MATERIALS*MAX_MATERIAL_PHASES];
  standard_precision TBlock[MAX_MATERIALS*MAX_MATERIAL_PHASES];

  /* Properties of the phase boundaries */

  standard_precision Clapeyron[MAX_MATERIALS*MAX_MATERIAL_PHASES];
  standard_precision PZ0[MAX_MATERIALS*MAX_MATERIAL_PHASES];
  standard_precision PT0[MAX_MATERIALS*MAX_MATERIAL_PHASES];

  /* Tracking particle phases */

  int *Current_phase;
  int *Equilibrium_phase;

  standard_precision *time_since_phase_change;
  standard_precision *phase_function;

  /* Material grain size at each tracer */

  standard_precision *grain_size;

  /* How keen are particles to eat each other
     and, on the other hand, do we let them reproduce ? */

  standard_precision particle_appetite;
  int particle_reproduction[MAX_MATERIALS];

#if defined (CHEM_TRANS_MODULE_INSTALLED)

  standard_precision *volfraction;
  standard_precision *sigma_n;
  
  standard_precision permeability_factor[MAX_MATERIALS];

#endif

} tracer;
