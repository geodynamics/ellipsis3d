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


#include "config.h"

#include <math.h> 

#if HAVE_STRING_H
#include <string.h>
#endif

#if HAVE_STRINGS_H
#include <strings.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"


void tracer_initialization(
			   struct All_variables *E
			   )
{
  void set_tracer_defaults();
  void tracer_allocate_memory();
  void tracer_initial_locations();

  set_tracer_defaults(E);
  tracer_allocate_memory(E);
  tracer_initial_locations(E);
}


void set_tracer_defaults(
			 struct All_variables *E
			 )
{
  char default_str[256],default_str1[256];
  int m;
  standard_precision xmin,xmax;

  E->tracer.NUM_TRACERS = 0;  

  input_boolean("Tracers",&(E->tracer.TRACERS),"on");
  input_int("Different_materials",&(E->tracer.NUM_MATERIALS),"essential");
  input_std_precision("Tracer_appetite",&(E->tracer.particle_appetite),"0.5,0.0,2.0");

  E->tracer.NUM_MATERIALS--;  /* As the background material is indexed at 0 */

  if (E->tracer.NUM_MATERIALS > MAX_MATERIALS) {
    E->tracer.NUM_MATERIALS = MAX_MATERIALS;
    fprintf(stderr,"Exceeded the hardwired limit of %d different materials\n",MAX_MATERIALS);
  }

  for(m=0;m<=E->tracer.NUM_MATERIALS;m++) {
    
    /* Does this material spawn new particles (assexually) */

    sprintf(default_str,"Material_%d_reproduction",m);
    input_boolean(default_str,&(E->tracer.particle_reproduction[m]),"on");

    /* Phase changes within a single material */

    sprintf(default_str,"Material_%d_phases",m);
    input_int(default_str,&(E->tracer.Phases[m]),"1,1,11");

    sprintf(default_str,"Material_%d_density",m);
    input_std_precision_vector(default_str,E->tracer.Phases[m],&(E->tracer.Density[m*MAX_MATERIAL_PHASES]),1.0);
    sprintf(default_str,"Material_%d_T_block",m);
    input_std_precision_vector(default_str,E->tracer.Phases[m],&(E->tracer.TBlock[m*MAX_MATERIAL_PHASES]),-1.0e32);
    sprintf(default_str,"Material_%d_Clapeyron",m);
    input_std_precision_vector(default_str,E->tracer.Phases[m]-1,&(E->tracer.Clapeyron[m*MAX_MATERIAL_PHASES]),0.0);
    sprintf(default_str,"Material_%d_Z0",m);
    input_std_precision_vector(default_str,E->tracer.Phases[m]-1,&(E->tracer.PZ0[m*MAX_MATERIAL_PHASES]),1.0e32);
    sprintf(default_str,"Material_%d_T0",m);
    input_std_precision_vector(default_str,E->tracer.Phases[m]-1,&(E->tracer.PT0[m*MAX_MATERIAL_PHASES]),0.0);

    sprintf(default_str,"Material_%d_therm_exp",m);
    input_std_precision(default_str,&(E->tracer.Therm_exp[m]),"0.0");
    sprintf(default_str,"Material_%d_therm_diff",m);
    input_std_precision(default_str,&(E->tracer.Therm_diff[m]),"0.0");
    sprintf(default_str,"Material_%d_Cp",m);
    input_std_precision(default_str,&(E->tracer.Cp[m]),"1.0");
    sprintf(default_str,"Material_%d_Qt",m);  
    input_std_precision(default_str,&(E->tracer.Qt[m]),"0.0");

    /* RAA: 24/09/02 -  C. O'Neill - Depletion stuff: expansion, solidus & liquidus */
    sprintf(default_str,"Material_%d_depl_exp",m);  
    input_std_precision(default_str,&(E->tracer.Depl_exp[m]),"0.0");  
    sprintf(default_str,"Material_%d_latent_heat",m);
    input_std_precision(default_str,&(E->tracer.Latent_heat[m]),"0.0");

    sprintf(default_str,"Material_%d_Red",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Red[m],0.0);
    sprintf(default_str,"Material_%d_Green",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Green[m],0.0);
    sprintf(default_str,"Material_%d_Blue",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Blue[m],0.0);
    sprintf(default_str,"Material_%d_Opacity",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Opacity[m],0.0);

    sprintf(default_str, "Material_%d_Red_hot",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Red_hot[m],-1.0);
    sprintf(default_str,"Material_%d_Green_hot",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Green_hot[m],-1.0);
    sprintf(default_str,"Material_%d_Blue_hot",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Blue_hot[m],-1.0);
    sprintf(default_str,"Material_%d_Opacity_hot",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Opacity_hot[m],-1.0);

    sprintf(default_str, "Material_%d_Red_strained",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Red_strained[m],-1.0);
    sprintf(default_str,"Material_%d_Green_strained",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Green_strained[m],-1.0);
    sprintf(default_str,"Material_%d_Blue_strained",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Blue_strained[m],-1.0);
    sprintf(default_str,"Material_%d_Opacity_strained",m);
    input_std_precision_vector(default_str,E->control.PPM_files,E->tracer.Opacity_strained[m],-1.0);
 

#if defined (CHEM_TRANS_MODULE_INSTALLED)
    sprintf(default_str,"Material_%d_permeability",m);
    input_std_precision(default_str,&(E->tracer.permeability_factor[m]),"0.0");
#endif

  }

  /* Here we get the sampling point tracer information */
  
  input_int("Sampling_tracers",&(E->tracer.SAMPLE_PTS),"0");
  input_std_precision_vector("Sampling_x",E->tracer.SAMPLE_PTS,E->tracer.sample_x,0.0);
  input_std_precision_vector("Sampling_z",E->tracer.SAMPLE_PTS,E->tracer.sample_z,0.0);
  input_std_precision_vector("Sampling_y",E->tracer.SAMPLE_PTS,E->tracer.sample_y,0.0);
  input_std_precision_vector("Sampling_plot_min",E->tracer.SAMPLE_PTS,E->tracer.sample_plotmin,1.0e32);
  input_std_precision_vector("Sampling_plot_max",E->tracer.SAMPLE_PTS,E->tracer.sample_plotmax,-1.0e32);
  input_std_precision_vector("Sampling_R",E->tracer.SAMPLE_PTS,E->tracer.sample_Red,0.0);
  input_std_precision_vector("Sampling_G",E->tracer.SAMPLE_PTS,E->tracer.sample_Green,0.0);
  input_std_precision_vector("Sampling_B",E->tracer.SAMPLE_PTS,E->tracer.sample_Blue,0.0);
  input_int_vector("Sampling_field",E->tracer.SAMPLE_PTS,E->tracer.sample_type,0);
  input_int_vector("Sampling_dirn",E->tracer.SAMPLE_PTS,E->tracer.sample_direction,0);
  input_int_vector("Sampling_plot_num",E->tracer.SAMPLE_PTS,E->tracer.sample_in_plot,0);
  input_int_vector("Sampling_lagrangian",E->tracer.SAMPLE_PTS,E->tracer.sample_lagrangian,0);
  input_int_vector("Sampling_normalize",E->tracer.SAMPLE_PTS,E->tracer.sample_normalize,0);

  for(m=0;m<E->tracer.SAMPLE_PTS;m++)
    E->tracer.sample_in_element[m] = 1;

  input_string("particle_input",E->tracer.tracer_input_file,"initialize");
}


void tracer_allocate_memory (
			     struct All_variables *E
			     )
{
  int lv,i;

  standard_precision xmin,xmax,t1;

  int former_highest_tracer;
  int current_highest_tracer;
  int old_allocation;

  static int current_allocation = 0;

  const int granularity = 10000;
  const int dofs=E->mesh.dof ;
  const int dims=E->mesh.nsd ;

  /* First check if we do actually need to 
     allocate more tracers. If so, allocate a
     block of them in one go.
  */

  /*  fprintf(stderr,"There are currently %d tracers allocated and %d required\n", 
	  current_allocation,E->tracer.NUM_TRACERS); */

  if(E->tracer.NUM_TRACERS > current_allocation) {  /* Yes we do */
    old_allocation = current_allocation;
    current_allocation += granularity;
    former_highest_tracer = old_allocation;
    current_highest_tracer = current_allocation;
  }
  else {
    return;
  }

  fprintf(E->fp,"Allocate %d extra tracers to total %d\n",
	  current_highest_tracer - former_highest_tracer, 
	  current_highest_tracer);
 
  if(E->control.verbose)
    fprintf(stderr,"Allocate %d extra tracers to total %d\n",
	  current_highest_tracer - former_highest_tracer, 
	  current_highest_tracer);

  /* Positions */ 

  Realloc0(&(E->tracer.tx),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.tz),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.tx1),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.tz1),old_allocation,current_allocation+10,sizeof(standard_precision)); 

  if(3==dims) {
    Realloc0(&(E->tracer.ty),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.ty1),old_allocation,current_allocation+10,sizeof(standard_precision));
   }

  /* Properties */

  Realloc0(&(E->tracer.T),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Tdot),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Q),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Hvisc),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Visc),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Visc0),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.BulkVisc),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.weight),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Density_change),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.edot),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.strs),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.strd),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.strd1),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.edot_integrated),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.edotp_integrated),old_allocation,current_allocation+10,sizeof(standard_precision));
  /*RAA: 24/09/02, C. O'Neill - melting stuff */
  Realloc0(&(E->tracer.depl),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.dFdot),old_allocation,current_allocation+10,sizeof(standard_precision));

  Realloc0(&(E->tracer.time_since_phase_change),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.phase_function),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.grain_size),old_allocation,current_allocation+10,sizeof(standard_precision));

  Realloc0(&(E->tracer.DQ),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.DQ1),old_allocation,current_allocation+10,sizeof(standard_precision));

  Realloc0(&(E->tracer.A1),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.A2),old_allocation,current_allocation+10,sizeof(standard_precision));
  if(3==dims) {
    Realloc0(&(E->tracer.A3),old_allocation,current_allocation+10,sizeof(standard_precision));
  }


  Realloc0(&(E->tracer.U1),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.U2),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.UU1),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.UU2),old_allocation,current_allocation+10,sizeof(standard_precision));
  if(3==dims) {
    Realloc0(&(E->tracer.U3),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.UU3),old_allocation,current_allocation+10,sizeof(standard_precision));
  }

  Realloc0(&(E->tracer.tr_dim),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.dX11),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.dX12),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.dX21),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.dX22),old_allocation,current_allocation+10,sizeof(standard_precision));
  if(3==dims) {
    Realloc0(&(E->tracer.dX13),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.dX23),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.dX33),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.dX31),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.dX32),old_allocation,current_allocation+10,sizeof(standard_precision));
  }


#if defined (CHEM_TRANS_MODULE_INSTALLED)
  Realloc0(&(E->tracer.volfraction),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.sigma_n),old_allocation,current_allocation+10,sizeof(standard_precision));
#endif



  if(E->control.ORTHOTROPY) {
     Realloc0(&(E->tracer.n1),old_allocation,current_allocation+10,sizeof(standard_precision));
     Realloc0(&(E->tracer.n2),old_allocation,current_allocation+10,sizeof(standard_precision));
     if(3==dims)
        Realloc0(&(E->tracer.n3),old_allocation,current_allocation+10,sizeof(standard_precision));
  }

  if(E->control.ELASTICITY) {

    Realloc0(&(E->tracer.Pt),old_allocation,current_allocation+10,sizeof(standard_precision));
    
    Realloc0(&(E->tracer.S11),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.S22),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.S12),old_allocation,current_allocation+10,sizeof(standard_precision));
    if(3==dims && 3==dofs) {
      Realloc0(&(E->tracer.S33),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S13),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S23),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S21),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S31),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S32),old_allocation,current_allocation+10,sizeof(standard_precision));
    }
    else if(2==dims && 3==dofs) {
      Realloc0(&(E->tracer.S21),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M31),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M32),old_allocation,current_allocation+10,sizeof(standard_precision));
    }
    else if(6==dofs) {
      Realloc0(&(E->tracer.S33),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S21),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S13),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S31),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S23),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.S32),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M12),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M21),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M13),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M31),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M23),old_allocation,current_allocation+10,sizeof(standard_precision));
      Realloc0(&(E->tracer.M32),old_allocation,current_allocation+10,sizeof(standard_precision));
    }
  
  }

  Realloc0(&(E->tracer.property_group),old_allocation,current_allocation+10,sizeof(int));
  Realloc0(&(E->tracer.yielded),old_allocation,current_allocation+10,sizeof(standard_precision));
  Realloc0(&(E->tracer.Current_phase),old_allocation,current_allocation+10,sizeof(int));
  Realloc0(&(E->tracer.Equilibrium_phase),old_allocation,current_allocation+10,sizeof(int));
  
  for(lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++) {
    Realloc0(&(E->tracer.tracer_elt[lv]),old_allocation,current_allocation+10,sizeof(int));
    Realloc0(&(E->tracer.tr_in_element[lv]),old_allocation,current_allocation+10,sizeof(int));
    Realloc0(&(E->tracer.tr_in_element_number[lv]),old_allocation,current_allocation+10,sizeof(int));
    Realloc0(&(E->tracer.tr_in_element_offset[lv]),old_allocation,current_allocation+10,sizeof(int));
    
    Realloc0(&(E->tracer.tracer_weight_fn[lv]),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.tracer_jacobian[lv]),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.tracer_std_weighting[lv]),old_allocation,current_allocation+10,sizeof(standard_precision));
    Realloc0(&(E->tracer.ls_diag[lv]),old_allocation,current_allocation+10,sizeof(standard_precision));
    
    Realloc0(&(E->tracer.eta1[lv]),old_allocation,current_allocation+10,sizeof(standard_precision)); 
    Realloc0(&(E->tracer.eta2[lv]),old_allocation,current_allocation+10,sizeof(standard_precision)); 
    Realloc0(&(E->tracer.eta3[lv]),old_allocation,current_allocation+10,sizeof(standard_precision)); 
    
    Realloc0(&(E->tracer.sfn_values[lv]),old_allocation,current_allocation+10,sizeof(struct TRACER_ELT_WEIGHT));

    /* initialize for element search */
 
    for(i=former_highest_tracer+1;i<=current_highest_tracer;i++) {
      E->tracer.tracer_elt[lv][i] = 1;
    }    
  }

  if(E->control.verbose)
    fprintf(stderr,"Allocated tracer memory ... initializing from %d to %d\n",
	    former_highest_tracer+1,current_highest_tracer );


  
  for(i=former_highest_tracer+1;i<=current_highest_tracer;i++) {
 
    if(E->control.ORTHOTROPY) {
      E->tracer.n1[i] = 0.0;
      E->tracer.n2[i] = 1.0;
      if(3==dims)
	E->tracer.n3[i] = 0.0;
    }

#if defined (CHEM_TRANS_MODULE_INSTALLED)
    if(E->control.CHEM_TRANS)
      E->tracer.sigma_n[i] = 0.0;
#endif
    


    if(E->control.ELASTICITY) {
      E->tracer.S11[i] = 0.0 ;
      E->tracer.S22[i] = 0.0 ;
      E->tracer.S12[i] = 0.0 ;
      if(3==dims && 3==dofs) {
	E->tracer.S33[i] = 0.0 ;
	E->tracer.S13[i] = 0.0 ;
	E->tracer.S23[i] = 0.0 ;
      }
      else if(2==dims && 3==dofs) {
	E->tracer.S21[i] = 0.0 ;
	E->tracer.M31[i] = 0.0 ;
	E->tracer.M32[i] = 0.0 ;
      }
      else if(6==dofs) {
	E->tracer.S33[i] = 0.0 ;
	E->tracer.S21[i] = 0.0 ;
	E->tracer.S13[i] = 0.0 ;
	E->tracer.S31[i] = 0.0 ;
	E->tracer.S23[i] = 0.0 ;
	E->tracer.S32[i] = 0.0 ;
	E->tracer.M31[i] = 0.0 ;
	E->tracer.M21[i] = 0.0 ;
	E->tracer.M12[i] = 0.0 ;
	E->tracer.M13[i] = 0.0 ;
	E->tracer.M23[i] = 0.0 ;
	E->tracer.M32[i] = 0.0 ;
      }     
      E->tracer.Pt[i] = 0.0;
    }

    E->tracer.edot[i] = 0.0;
    E->tracer.strs[i] = 0.0;
    E->tracer.strd[i] = 0.0;
    E->tracer.strd1[i] = 0.0;
    E->tracer.edot_integrated[i] = 0.0;
    E->tracer.edotp_integrated[i] = 0.0;
    E->tracer.T[i] = 0.0;
    E->tracer.Q[i] = 0.0;
    E->tracer.depl[i] = 0.0;   /*RAA: 24/09/02, C. O'Neill - melting stuff */
    E->tracer.time_since_phase_change[i]=1.0e32;
    E->tracer.Density_change[i]=1.0;
    E->tracer.grain_size[i]=1.0;
    E->tracer.Hvisc[i] = 0.0;
    E->tracer.DQ[i]=0.0;
    E->tracer.DQ1[i]=0.0;
    E->tracer.A1[i]=0.0;
    E->tracer.A2[i]=0.0;

    E->tracer.U1[i]=0.0;
    E->tracer.U2[i]=0.0;
    E->tracer.UU1[i]=0.0;
    E->tracer.UU2[i]=0.0;
   
    if(3==dims) {
      E->tracer.U3[i]=0.0;
      E->tracer.UU3[i]=0.0;
    }



    E->tracer.yielded[i]=0.0;

    if(3==dims) 
      E->tracer.A3[i]=0.0;

    E->tracer.Visc[i] = 1.0e32;

    E->tracer.dX11[i] = 0.005;
    E->tracer.dX12[i] = 0.0;
    E->tracer.dX21[i] = 0.0;
    E->tracer.dX22[i] = 0.005;
    
    if(3==dims) {
      E->tracer.dX13[i] = 0.0;
      E->tracer.dX23[i] = 0.0;
      E->tracer.dX33[i] = 0.005;
      E->tracer.dX31[i] = 0.0;
      E->tracer.dX32[i] = 0.0;
     }
  }

  /* And now it is possible to add the pointers to this
     allocated space in the list of available input/output 
     data - since we may now be moving these pointers, this
     must be repeated for each new allocation - and the 
     previous settings erased.
  */

  E->control.particle_data.numb=0;

  add_variable_to_tracer_output(E,E->tracer.T,"Temp");
  add_variable_to_tracer_output(E,E->tracer.Q,"Pres"); /*RAA: 22/05/01, was Pt, but is for elasticity only -led to Seg Faults*/
  add_variable_to_tracer_output(E,E->tracer.Visc,"Visc");
  add_variable_to_tracer_output(E,E->tracer.edot,"Edot");
  add_variable_to_tracer_output(E,E->tracer.edotp_integrated,"StrP");
  add_variable_to_tracer_output(E,E->tracer.edot_integrated,"StrT");
  add_variable_to_tracer_output(E,E->tracer.grain_size,"Grsz");
  add_variable_to_tracer_output(E,E->tracer.strs,"Strs"); /*RAA: 11/01/02, add stress as a variable*/
  add_variable_to_tracer_output(E,E->tracer.strd,"Strd"); /*RAA: 22/05/02, add stressd as a variable*/
  add_variable_to_tracer_output(E,E->tracer.depl,"Depl");  /*RAA: 24/09/02, C. O'Neill - melting stuff */

  if(E->control.ELASTICITY) {

    add_variable_to_tracer_output(E,E->tracer.S11,"S11");
    add_variable_to_tracer_output(E,E->tracer.S22,"S22");
    add_variable_to_tracer_output(E,E->tracer.S33,"S33");
    add_variable_to_tracer_output(E,E->tracer.S12,"S12");
    add_variable_to_tracer_output(E,E->tracer.S21,"S21");
    add_variable_to_tracer_output(E,E->tracer.S13,"S13");
    add_variable_to_tracer_output(E,E->tracer.S31,"S31");
    add_variable_to_tracer_output(E,E->tracer.S23,"S23");
    add_variable_to_tracer_output(E,E->tracer.S32,"S32");
    add_variable_to_tracer_output(E,E->tracer.M12,"M12");
    add_variable_to_tracer_output(E,E->tracer.M21,"M21");
    add_variable_to_tracer_output(E,E->tracer.M13,"M13");
    add_variable_to_tracer_output(E,E->tracer.M31,"M31");
    add_variable_to_tracer_output(E,E->tracer.M23,"M23");
    add_variable_to_tracer_output(E,E->tracer.M32,"M32");
    add_variable_to_tracer_output(E,E->tracer.Pt,"Pt");  /*RAA: may lead to seg fault for 3D elasticity??  check*/

  }

   /* General storage for output from sampling points (currently sampling at 101 locations) */

  for(i=0;i<=E->tracer.SAMPLE_PTS+3;i++) {
    E->tracer.sampled_data[i] = (standard_precision *) Malloc0(103 * sizeof(standard_precision));
  }

  /* load up the coordinate while we're here */

  xmin =  min(E->x[1][1],E->x[1][E->mesh.noz]);
  xmax = max(E->x[1][E->mesh.nno-E->mesh.noz+1],E->x[1][E->mesh.nno]);

  for(i=0;i<=100;i++) { 
    E->tracer.sampled_data[0][i] = 
      xmin + (xmax - xmin) * i / 100.00;
  } 

  xmin = min(E->x[2][1],E->x[2][E->mesh.nno-E->mesh.noz+1]);
  xmax = max(E->x[2][E->mesh.noz],E->x[2][E->mesh.nno]);

  for(i=0;i<=100;i++) {
    E->tracer.sampled_data[1][i] = 
      xmin + (xmax - xmin) * i / 100.00;
  }

/*RAA: 4/01/02, add some lines for 3D, assuming the domain is rectangular*/
  if(3==dims) {
     xmin = min(E->x[3][1],E->x[3][E->mesh.noz]);
     xmax = max(E->x[3][E->mesh.nno-E->mesh.noz+1],E->x[3][E->mesh.nno]);

     for(i=0;i<=100;i++) 
       E->tracer.sampled_data[2][i] = xmin + (xmax - xmin) * i / 100.00;
     }
} 
