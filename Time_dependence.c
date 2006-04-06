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



/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


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


void set_convection_defaults(
			     struct All_variables *E
			     )
{
    void PG_timestep();
    void PIC_timestep();
    void convection_output();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void node_locations();
    void twiddle_thumbs();
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void convection_initial_fields();
    void phase_transition_initialization();
  
    E->next_buoyancy_field = PG_timestep;
    E->special_process_new_velocity = convection_output; 
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;
    E->advection_source_term = NULL;
  
    E->velocity_boundary_conditions=velocity_boundary_conditions;
    E->buoyancy_boundary_conditions=temperature_boundary_conditions;
return;
}

void read_convection_settings(
			      struct All_variables *E
			      )
{ 
    void advection_diffusion_parameters();
    int i;

/* parameters */

    input_int("B_matrix",&(E->control.B_matrix),"1");  /* Defaults to the old, hand-coded viscous flow model */
    input_std_precision("rayleigh",&(E->control.Atemp),"0.0");

    input_int("heating_elements",&(E->convection.heat_sources.number),"0,0,10");
    input_std_precision("heating_time_offset",&(E->convection.heat_sources.t_offset),"0.0");
    input_std_precision_vector("heating_Qs",
			       E->convection.heat_sources.number,E->convection.heat_sources.Q);
    input_std_precision_vector("heating_lambdas",
			       E->convection.heat_sources.number,E->convection.heat_sources.lambda);

    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");

    for(i=0;i<=E->convection.number_of_perturbations;i++)
      E->convection.perturb_mag[i] = E->convection.perturb_k[i] = E->convection.perturb_ky[i] = 0.0;

    input_std_precision_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_std_precision_vector("perturbk",E->convection.number_of_perturbations,E->convection.perturb_k);
    input_std_precision_vector("perturbky",E->convection.number_of_perturbations,E->convection.perturb_ky);

    input_string("prevT",E->convection.old_T_file,"initialize");

    advection_diffusion_parameters(E);

 
    return;
}
/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(
     struct All_variables *E
)
{ 
  /* Nothing here for the time being */
return;
}

void convection_allocate_memory(
     struct All_variables *E
)
{ 
  void advection_diffusion_allocate_memory();
  void convection_setup_output();
   
  advection_diffusion_allocate_memory(E); 
 
  /* Once memory is allocated, the output can be setup */

  convection_setup_output(E);
 
return;
}

/* ============================================ */
void convection_initial_fields(
     struct All_variables *E
)
{ 
    void convection_initial_temperature();
    void tracer_initialization();

    if(E->control.verbose)
      fprintf(stderr,"convection, initial temperature\n");
    convection_initial_temperature(E);
    if(E->control.verbose)
      fprintf(stderr,"convection, initial temperature ... done\n");

#if defined(TRACER_MODULE_INSTALLED)  

    if(E->control.verbose)
      fprintf(stderr,"tracer, initialization\n");
    tracer_initialization(E);
    if(E->control.verbose)
      fprintf(stderr,"tracer, initialization ... done\n");

#endif     

  return;
}
/* =========================================== */

void convection_boundary_conditions(
     struct All_variables *E
)
{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);

    return;
}
/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(
     struct All_variables *E
)
{
    int i,j,k,p,node;
    int el;
    higher_precision temp,base,radius,radius2;
    
    standard_precision eta1,eta2,eta3;
    standard_precision dOmega;

    double drand48();
    FILE *fp;
    void vcopy();
    void temperatures_conform_bcs();
    
    int in1,in2,in3,instance;
    standard_precision x1,y1,z1,weight;

    int read_previous_field();
    void field_arbitrary_rectangle_file();
    void field_arbitrary_circle_file();
    void field_arbitrary_harmonic_file();

    const int dims = E->mesh.nsd;

    /* This is what I want to do:
       Try reading a `new' format temperature,
       if that fails, check that there was no old format file
       and do initial stuff.
       
       This option will persist until the old format files are no longer a problem
       (i.e. when I write an awk script to convert them).

       NB this indulgence is only for temperature fields as they
       take so much work to calculate.
       */
 
    if(read_previous_field(E,E->T,"temperature","Temp")==0)      { /*1*/
      if(strstr(E->convection.old_T_file,"initialize") != NULL) { /*2*/
	base=1.0;

	  fprintf(E->fp,"Linear temperature gradient (%g -> %g)\n",
		  E->control.TBCtopval,E->control.TBCbotval);
	  fflush(E->fp);

	  for(i=1;i<=E->mesh.noz;i++) { /*3*/
	 
	    temp = E->control.TBCtopval+(E->control.TBCbotval-E->control.TBCtopval)
	      * (i-1.0)/(E->mesh.noz-1.0);
	  
	    for(j=1;j<=E->mesh.nox;j++)
	      for(k=1;k<=E->mesh.noy;k++) {
		node = i+(j-1)*E->mesh.noz+(k-1)*E->mesh.noz*E->mesh.nox;
		E->T[node] = temp;
	      }
	  } /*3*/
      } /*2*/
    } /*1*/

 
    /* Then add any other shapes etc on top */

    field_arbitrary_rectangle_file(E,1,&(E->temperature.Trects),"Temp",E->T,NULL,0,0,E->mesh.levmax);
    field_arbitrary_circle_file(E,1,&(E->temperature.Tcircs),"Temp",E->T,NULL,0,0,E->mesh.levmax);
    field_arbitrary_harmonic_file(E,1,&(E->temperature.Tharms),"Temp",E->T,NULL,0,0,E->mesh.levmax);
     
    /* When everything is specified, add some perturbations */
    
    fprintf(E->fp,"Adding in %d perturbations to temperature field\n",E->convection.number_of_perturbations);
    fflush(E->fp);

    for(j=1;j<=E->mesh.nox;j++)
      for(i=1;i<=E->mesh.noz;i++) 	
  	for(k=1;k<=E->mesh.noy;k++) {
	  node = i+(j-1)*E->mesh.noz+(k-1)*E->mesh.noz*E->mesh.nox;
	  x1=E->x[1][node];
	  z1=E->x[2][node];
 	  y1= ((3 == E->mesh.nsd) ? E->x[3][node] : 0.0); 

	  if(3==E->mesh.nsd)
	    if(!E->mesh.periodic_x)
	      for(p=0;p<E->convection.number_of_perturbations;p++) {
		E->T[node] += E->convection.perturb_mag[p] *
		  sin(M_PI*z1)  *
		  cos(E->convection.perturb_k[p] * M_PI * x1) *
		  cos(E->convection.perturb_ky[p] * M_PI * y1) ;
	      }
	    else
	      for(p=0;p<E->convection.number_of_perturbations;p++) {
		E->T[node] += E->convection.perturb_mag[p] *
		  sin(M_PI*z1)*
		  sin(E->convection.perturb_k[p] * M_PI*x1) *
		  cos(E->convection.perturb_ky[p] * M_PI*y1);
		
	    }
	  else
	    if(!E->mesh.periodic_x) {
	      for(p=0;p<E->convection.number_of_perturbations;p++) {
		E->T[node] += E->convection.perturb_mag[p] *
		  sin(M_PI*z1)  *
		  cos(E->convection.perturb_k[p] * M_PI * x1) ;
	      }
	    }
	    else
	      for(p=0;p<E->convection.number_of_perturbations;p++) {
		E->T[node] += E->convection.perturb_mag[p] *
		  sin(M_PI*z1)*
		  sin(E->convection.perturb_k[p] * M_PI*x1) ;
	      }
	}



    temperatures_conform_bcs(E,E->T); 
    return_horiz_ave(E,E->T,E->Have.T); 
    /* thermal_buoyancy(E); */

  /*RAA: 7/6/01, check for proper temps at periodic nodes*/
  /*  if(E->control.verbose) {
        for(p=1;p<=E->mesh.NNO[E->mesh.levmax];p++)
           fprintf(stderr,"(at convection_initial_temperature) Temps at nodes: node  %d %g\n",p,E->T[p]);
      }
  */ 

  return;
}

/* Load output types buffet for convection problem */


void convection_setup_output(
	  struct All_variables *E
)
{
  const int dims = E->mesh.nsd ;
  const int dofs = E->mesh.dof ;
    /* Node output data */


  if(2==dims && 2==dofs) {
    add_variable_to_node_output(E,E->sv[1],"Velx") ;
    add_variable_to_node_output(E,E->sv[2],"Velz") ;
  }
  else if(2==dims && 3==dofs) {
    add_variable_to_node_output(E,E->sv[1],"Velx") ;
    add_variable_to_node_output(E,E->sv[2],"Velz") ;
    add_variable_to_node_output(E,E->sv[3],"Roty") ;
  }
  else if(3==dims && 3==dofs) {
    add_variable_to_node_output(E,E->sv[1],"Velx") ;
    add_variable_to_node_output(E,E->sv[2],"Velz") ;
    add_variable_to_node_output(E,E->sv[3],"Vely") ;
  }
  else {
    add_variable_to_node_output(E,E->sv[1],"Velx") ;
    add_variable_to_node_output(E,E->sv[2],"Velz") ;
    add_variable_to_node_output(E,E->sv[3],"Vely") ;
    add_variable_to_node_output(E,E->sv[4],"Rotx") ;
    add_variable_to_node_output(E,E->sv[5],"Rotz") ;
    add_variable_to_node_output(E,E->sv[6],"Roty") ;
  }

  add_variable_to_node_output(E,E->T,"Temp");
  add_variable_to_node_output(E,E->nQ,"Pres");
  add_variable_to_node_output(E,E->depl,"Depl"); /*RAA; O'Neill's melting stuff*/
  
  if(!E->control.CYLINDER && 3!=E->mesh.dof) 
    add_variable_to_node_output(E,E->Psi,"Strf");
 
 add_variable_to_node_output(E,E->Pstrain,"Pstn");
 
  /* Horizontal Averages */

  add_variable_to_haverage_output(E,E->Have.T,"Temp");
  add_variable_to_haverage_output(E,E->Have.Vi,"Visc");
  add_variable_to_haverage_output(E,E->Have.vrms,"Velo");
  add_variable_to_haverage_output(E,E->Have.V[1],"Urms");
  add_variable_to_haverage_output(E,E->Have.V[2],"Vrms");
  if(3==E->mesh.dof)
    add_variable_to_haverage_output(E,E->Have.V[3],"Wrms");
  
  /* Observables */

  add_variable_to_slice_output(E,E->slice.tpg,"Tpgx");
  add_variable_to_slice_output(E,E->slice.tpgb,"Tpbx");
  add_variable_to_slice_output(E,E->slice.tpgk,"Tpgk");
  add_variable_to_slice_output(E,E->slice.tpgbk,"Tpbk");
  add_variable_to_slice_output(E,E->slice.grv,"Grvx");
  add_variable_to_slice_output(E,E->slice.grvb,"Grbx");
  add_variable_to_slice_output(E,E->slice.grvt,"Grtx");
  add_variable_to_slice_output(E,E->slice.grvk,"Grvk");
  add_variable_to_slice_output(E,E->slice.grvbk,"Grbk");
  add_variable_to_slice_output(E,E->slice.grvtk,"Grtk");
  add_variable_to_slice_output(E,E->slice.geo,"Geox");
  add_variable_to_slice_output(E,E->slice.geob,"Gebx");
  add_variable_to_slice_output(E,E->slice.geot,"Getx"); /*RAA: 6/11/02, bugfix*/
  add_variable_to_slice_output(E,E->slice.geok,"Geok");
  add_variable_to_slice_output(E,E->slice.geotk,"Getk"); /*RAA: 6/11/02, bugfix*/
  add_variable_to_slice_output(E,E->slice.geobk,"Gebk");
  add_variable_to_slice_output(E,E->slice.shflux,"Shfl");
  add_variable_to_slice_output(E,E->slice.bhflux,"Bhfl");
  add_variable_to_slice_output(E,E->slice.vxsurf[1],"Vxsf");
  add_variable_to_slice_output(E,E->slice.vxsurf[2],"Vzsf");
  add_variable_to_slice_output(E,E->slice.vxsurf[3],"Vysf");
  
    /* Timelog stuff */

  add_variable_to_timelog_output(E,&(E->monitor.Nusselt),"Nuss");
  add_variable_to_timelog_output(E,&(E->monitor.F_surface),"Shfl");
  add_variable_to_timelog_output(E,&(E->monitor.F_base),"Bhfl");
  add_variable_to_timelog_output(E,&(E->monitor.Sigma_max),"Stmx");
  add_variable_to_timelog_output(E,&(E->monitor.Vrms_surface),"Svav");
  add_variable_to_timelog_output(E,&(E->monitor.Vrms_base),"Bvav");
  add_variable_to_timelog_output(E,&(E->monitor.Vrms),"Vrms");
  add_variable_to_timelog_output(E,&(E->monitor.Vxrms),"Vxrm");
  add_variable_to_timelog_output(E,&(E->monitor.Vyrms),"Vyrm");
  add_variable_to_timelog_output(E,&(E->monitor.Vzrms),"Vzrm");
  add_variable_to_timelog_output(E,&(E->monitor.Vis2),"Vis2");

  return;
}

void convection_output(
		       struct All_variables *E,
		       int ii
		       )
{ 
  void generic_data_storage();
  standard_precision return_bulk_value();
    
  int i,j,k,p,a1,nint,n,el;
  int node,hnode;        
  int this_node,that_node;
  
  struct Shape_function GN;
  struct Shape_function_dx GNx;
  struct Shape_function_dA dOmega;
  
  const int vpts = vpoints[E->mesh.nsd];
  const int ends = enodes[E->mesh.nsd];
  
  static int been_here=0;
  
  /* if(been_here++==0)
     convection_setup_output(E); */
 
  /*RAA: 13/5/01*/
  if(E->control.verbose)
    fprintf(stderr,"WRITING DATA TO node_data FILE\n");

  generic_data_storage(E,ii);

  /*RAA: 10/4/02*/
  if(E->control.verbose)
    fprintf(stderr,"generic_data_storage finished.\n");

  return;
}
