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



/* Set up the finite element problem to suit: returns with all memory */
/* allocated, temperature, viscosity, node locations and how to use */
/* them all established. 8.29.92 or 29.8.92 depending on your nationality*/

#include "config.h"

#include <math.h>

#if HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <signal.h>

#include "element_definitions.h"
#include "global_defs.h"

int Emergency_stop;

void read_instructions(
		       struct All_variables *E,
		       int argc,
		       char **argv
		       )
{
  int get_process_identifier();
  
  void allocate_common_vars();  
  void common_initial_fields();
  void read_initial_settings();
  void global_default_values();
  void global_derived_values();
  void construct_masks();
  void construct_shape_functions();
  void construct_id();
  void construct_lm();
  void construct_sub_element();
  void element_locations_from_nodes();
  void mass_matrix();
  void construct_node_ks();
  void construct_node_maps();
  void construct_node_maps_3();
  void construct_ien();
  void construct_elem();
  RETSIGTYPE interuption();
  void check_bc_consistency();
  void node_locations();
  void node_colouring();
  void allocate_observables_storage();
  void setup_parser();
  
  char output_file[255],command[255];   

  /* =====================================================
     Global interuption handling routine defined once here
     =====================================================  */
  
  Emergency_stop = 0;
  signal(SIGINT,interuption);
  signal(SIGTERM,interuption);
  
  E->control.PID=get_process_identifier(); 
  
  /* ==================================================
     Initialize from the command line 
     from startup files. (See Parsing.c).
     ==================================================  */

  setup_parser(E,argc,argv);
  global_default_values(E);  

  read_initial_settings(E);  
  
  global_derived_values(E);
   
  (E->problem_derived_values)(E);
  
  allocate_common_vars(E);                 
  allocate_observables_storage(E);

  (E->problem_allocate_vars)(E);

  (E->solver_allocate_vars)(E);
  
  (E->problem_node_positions)(E,0);

  
  check_bc_consistency(E);
  
  construct_ien(E);
  construct_masks(E);		/* order is important here */
  construct_sub_element(E);
  construct_elem(E);
  
  if(E->mesh.nsd!=3)  /*RAA: 4/3/01, coloring is not correct yet for 3D*/
    node_colouring(E);   
  /*  construct_coord_hierarchy(E); */
  
  /* Boundary conditions are now input (potentially
     in a curvilinear system ... compute rotation matrixes for
     appropriate nodes */

  (E->velocity_boundary_conditions)(E);
  (E->buoyancy_boundary_conditions)(E);
  compute_node_R(E);
  
  construct_id(E);
  construct_lm(E);
  
  element_locations_from_nodes(E); 
  construct_shape_functions(E);
  mass_matrix(E);

  (E->problem_initial_fields)(E);    /* temperature/chemistry/melting etc */
  common_initial_fields(E);          /* velocity/pressure/... */
  
  read_viscosity_variables(E);
    
  shutdown_parser(E);
  
  /* Having done that, try and plot initial conditions
     for particles to make sure everything is OK */
  
#if 1 /*RAA: 10/04/02, this was #if 0 for Linux, but put it back in*/
  fprintf(stderr,"Producing initial startup picture\n");
  fprintf(stderr,"IN INSTRUCTIONS, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);
  sprintf(output_file,"%s.start.ppm",E->control.data_file);
   generate_2Ddata_pixmap(E,output_file); 
  fprintf(stderr,"Producing initial startup picture ... done\n");
  fprintf(stderr,"IN VISC HERE, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);
#endif  


  return;
}

void global_default_values(
			   struct All_variables *E
			   )
{
  FILE *fp;
  char command[1000];


  /* ZERO: things which are set by the machinery we are using */

#if defined (__uxp__)
  E->control.vector_optimization = 1;
#else
  E->control.vector_optimization = 0;
#endif

  /* FIRST: values which are not changed routinely by the user */
  E->control.test=0 ;
  E->control.accuracy = 1.0e-4;
  E->control.verbose=0; /* debugging/profiles */

  /* SECOND: values for which an obvious default setting is useful */

  E->control.ORTHO = 1; /* assume short cutting for orthogonal meshes by default */
  E->control.ORTHOZ = 1; /* assume short cutting for orthogonal meshes by default */
 
  E->control.CART2D = 0;
  E->control.CART3D = 0;
  E->control.CART2pt5D = 0;
  E->control.AXI = 0;
  E->control.SPHERE=0;
  E->control.CYLINDER=0;

  E->control.B_matrix=1 ;
  
  E->control.CONJ_GRAD = 0;
  E->control.MULTIGRID = 0;
  E->control.NODE_BY_NODE=0;
  E->control.ELEMENT_BY_ELEMENT=0;
  E->control.COMPRESS = 1;

/* Default: all optional modules set to `off' */
  E->control.MELTING_MODULE = 0;
  E->control.CHEMISTRY_MODULE = 0;
  E->control.CHEM_TRANS=0;

  E->control.ORTHOTROPY=0;

  E->mesh.levmax=0;
  E->mesh.levmin=0;
  E->mesh.nox = 1;  
  E->mesh.noz = 1;
  E->mesh.noy = 1;

  E->mesh.equilibrium_upper_surface=1;
  E->mesh.equilibrium_lower_surface=1;
    
  sprintf(E->viscosity.old_file,"initialize");
 
  E->mesh.periodic_x=0; /* reflection is default*/
  E->mesh.periodic_y=0;
 
  E->data.grav_acc = 9.81;
  
  /* THIRD: you forgot and then went home, let's see if we can help out */
  
  sprintf(E->control.data_file,"ellipsis.tmp.%d",getpid());
  
  E->mesh.layer0[1] =  E->mesh.layer0[2] =  E->mesh.layer0[3] = 1.0;
  E->mesh.layer1[1] =  E->mesh.layer1[2] =  E->mesh.layer1[3] = 1.0;
  E->monitor.elapsed_time=0.0;
  E->monitor.elapsed_time_vsoln=0.0;
  
  E->monitor.time_of_last_store=0.0;
  E->monitor.time_of_last_checkpt = 0.0;
  return;
}

void global_derived_values(
			   struct All_variables *E
			   )
{
  int d,lx,lz,ly,i,nox,noz,noy;
  char logfile[100];
  char debugfile[100]; /*RAA: 12/3/02, added this-removed info for file fred1*/
  char command[1000];
  char *name2, *name3 ;

  FILE *fp,*fp1;
 /* name2 = "fred2" ;
  name3 = "fred3" ; */
/**/
  /* Now the data file names are known,
     copy the input file to keep a record of
     what went on */
  
  sprintf(command,"cp  %s %s.input < /dev/null \n",E->control.input_file_name,E->control.data_file);
  system(command);
  
  /* As early as possible, set up the log file to 
     record information about the progress of the 
     program as it runs */
  
  sprintf(logfile,"%s.log",E->control.data_file);
  sprintf(debugfile,"%s.debug",E->control.data_file); /*RAA: 12/3/02, added this*/
  
  if((fp=fopen(logfile,"w")) == NULL)
    E->fp = stdout;
  else
    E->fp = fp;

  /*RAA: 12/3/02, added these lines for new file- removed fp1 as fred1*/
  if((fp1=fopen(debugfile,"w")) == NULL)
    E->fp1 = stdout;
  else
    E->fp1 = fp1;

  /*if((fp2=fopen(name2,"w")) == NULL)
    E->fp2 = stdout;
  else
    E->fp2 = fp2; */
/**/

  E->mesh.levmax=E->mesh.levels-1;
  
  E->mesh.q_levmin = E->mesh.levmax - (E->mesh.q_levels - 1);
  if(E->mesh.q_levmin < 1)
    E->mesh.q_levmin = 1;
  if(E->mesh.q_levmin > E->mesh.levmax)
    E->mesh.q_levmin = E->mesh.levmax;
  
  E->mesh.nox = E->mesh.mgunitx * (int) pow(2.0,((double)E->mesh.levmax)) + 1;
  E->mesh.noz = E->mesh.mgunitz *(int) pow(2.0,((double)E->mesh.levmax)) + 1; 
  if(E->mesh.nsd == 3 ) {
    E->mesh.noy = E->mesh.mgunity * (int) pow(2.0,((double)E->mesh.levmax)) + 1; 
  }
  else {
    E->mesh.noy = 1;
  }

  E->mesh.nnx[1] = E->mesh.nox;	
  E->mesh.nnx[2] = E->mesh.noz;	
  E->mesh.nnx[3] = E->mesh.noy;	
  E->mesh.nex[1] = E->mesh.elx = E->mesh.nox-1;	
  E->mesh.nex[2] = E->mesh.elz = E->mesh.noz-1;
  E->mesh.nex[3] = E->mesh.ely = max(E->mesh.noy-1,1);
  E->mesh.nno = E->mesh.nel = 1;
  E->mesh.nmx = max(E->mesh.nox,max(E->mesh.noz,E->mesh.noy));
  
  for(d=1;d<=E->mesh.nsd;d++)  {
    E->mesh.nno *= E->mesh.nnx[d];
    E->mesh.nel *= E->mesh.nex[d];
  }

  E->mesh.npno = E->mesh.nel;
  E->mesh.nsf = E->mesh.nox*E->mesh.noy;
 
  for(i=E->mesh.levmax;i>=E->mesh.levmin;i--)  /* set up dimensions for different grids  */  {  
    nox = E->mesh.mgunitx * (int) pow(2.0,(double)i) + 1;
    noz = E->mesh.mgunitz * (int) pow(2.0,(double)i) + 1;
    if(E->mesh.nsd==3)
      noy = E->mesh.mgunity * (int) pow(2.0,(double)i) + 1;
    else 
      noy = 1;

    E->mesh.ELX[i] = nox-1;
    E->mesh.ELZ[i] = noz-1;
    E->mesh.ELY[i] = max(noy-1,1);
    E->mesh.NNO[i] = nox * noz * noy;
    /*    E->mesh.NLNO[i] = nox * noz * noy;*/
    E->mesh.NEL[i] = (nox-1) * (noz-1) * max((noy-1),1);
    E->mesh.NPNO[i] = E->mesh.NEL[i] ;
    E->mesh.NEQ[i] = E->mesh.dof * E->mesh.NNO[i] ;  
 
    E->mesh.NOX[i] = nox;
    E->mesh.NOZ[i] = noz;
    E->mesh.NOY[i] = noy;
    
    E->mesh.NMX[i] = max(nox,max(noz,noy));
  }
  
  if(E->control.print_convergence)
    fprintf(stderr,"Problem has %d x %d x %d nodes\n",E->mesh.nox,E->mesh.noz,E->mesh.noy);
  
  return; 
}

void read_initial_settings(
			   struct All_variables *E
			   )
{
  void set_convection_defaults();
  void set_axi_defaults();
  void set_2dc_defaults();
  void set_3dc_defaults();
  void set_axicoss_defaults();
  void set_2dccoss_defaults();
  void set_3dccoss_defaults();
  void set_direct_defaults();
  void set_mg_defaults();
  
  int i;

 /* Some problems are neatly integrated, others are mutually
    exclusive. Mutually exclusive options are compiled out differently
    and are chosen at this point. This could, however, be done more
    effectively at run time (HINT !!) */
    
#if ( ! defined (GRAIN_GROWTH_MODULE_INSTALLED)) && ( ! defined (CHEMISTRY_MODULE_INSTALLED) )
  E->control.CONVECTION = 1; 
  set_convection_defaults(E);
#elif defined (GRAIN_GROWTH_MODULE_INSTALLED)
  E->control.CONVECTION = 1; 
  set_convection_grain_growth_defaults(E);
#else 
  E->control.CONVECTION = 1; 
  set_convection_chemistry_defaults(E);
#endif  
  
  input_string("Geometry",E->control.GEOMETRY,NULL); 
  if ( strcmp(E->control.GEOMETRY,"cart2d") == 0)  {
    set_2dc_defaults(E);
  }
  else  if ( strcmp(E->control.GEOMETRY,"cart2dcoss") == 0)  {
    set_2dccoss_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"axi") == 0)    {
    set_axi_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"axicoss") == 0)    {
    set_axicoss_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"cart2pt5d") == 0) {
    set_2pt5dc_defaults(E);
    }
  else if ( strcmp(E->control.GEOMETRY,"cart2pt5dcoss") == 0) {
    set_2pt5dccoss_defaults(E);
    }
  else if ( strcmp(E->control.GEOMETRY,"cart3d") == 0) {
    set_3dc_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"cart3dcoss") == 0) {
    set_3dccoss_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"cylinder") == 0)    {
    set_cylinder_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"cylindercoss") == 0)    {
    set_cylindercoss_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"sphere") == 0) {
    set_sphere_defaults(E);
  }
  else if ( strcmp(E->control.GEOMETRY,"spherecoss") == 0) {
    set_spherecoss_defaults(E);
  }
  else {
    fprintf(E->fp,"Unable to determine geometry, assuming cartesian 2d ... \n");
    E->control.CART2D = 1; 
    set_2dc_defaults(E); 
  }

  input_string("Solver",E->control.SOLVER_TYPE,NULL);
  if (1 ||  strcmp(E->control.SOLVER_TYPE,"multigrid") == 0)  {
    E->control.MULTIGRID = 1;
    E->control.NODE_BY_NODE=1;
    set_mg_defaults(E);
  }

  else  {
    fprintf(stderr,"Unable to determine how to solve, specify Solver=VALID_OPTION \n");
    exit(0); 
  }

  /* admin */

    input_string("Spacing_x",E->control.NODE_SPACING[0],"region");
    input_string("Spacing_z",E->control.NODE_SPACING[1],"region");
    input_string("Spacing_y",E->control.NODE_SPACING[2],"region");

    if ( strcmp(E->control.NODE_SPACING[0],"region") == 0)
      E->control.GRID_TYPE[0] = 1;
    else if ( strcmp(E->control.NODE_SPACING[0],"file") == 0)
      E->control.GRID_TYPE[0] = 2;
    else {
      E->control.GRID_TYPE[0] = 1;
    }
    if ( strcmp(E->control.NODE_SPACING[1],"region") == 0)
      E->control.GRID_TYPE[1] = 1;
    else if ( strcmp(E->control.NODE_SPACING[1],"file") == 0)
      E->control.GRID_TYPE[1] = 2;
    else {
      E->control.GRID_TYPE[1] = 1;
    }
    if ( strcmp(E->control.NODE_SPACING[2],"region") == 0)
      E->control.GRID_TYPE[2] = 1;
    else if ( strcmp(E->control.NODE_SPACING[2],"file") == 0)
      E->control.GRID_TYPE[2] = 2;
    else {
      E->control.GRID_TYPE[2] = 1;
    }
    
    /* Information on which files to print, which variables of the flow to calculate and print.
       Default is no information recorded (apart from special things for given applications.
    */

    input_string("datatypes",E->control.node_data.which_data_types,"");
    input_string("averages",E->control.haverage_data.which_data_types,"");
    input_string("timelog",E->control.time_data.which_data_types,"");
    input_string("observables",E->control.slice_data.which_data_types,"");
    input_string("particle_data",E->control.particle_data.which_data_types,"");
    input_string("previous_particle_data",E->control.particle_data.which_data_types_in,"");

    input_string("datafile",E->control.data_file,"initialize");
    input_string("process_command",E->control.output_written_external_command,"");
    input_string("ppm_process_command",E->control.ppm_written_external_command,"");
    input_boolean("AVS",&(E->control.AVS),"off");
    
    input_int("mgunitx",&(E->mesh.mgunitx),"2");
    input_int("mgunitz",&(E->mesh.mgunitz),"2");
    input_int("mgunity",&(E->mesh.mgunity),"2");
    input_int("levels",&(E->mesh.levels),"1");
    input_int("pmg_levels",&(E->mesh.q_levels),"1");
    input_int("activate_secret_output",&(E->control.outputcase),"0");
  
    input_std_precision("aug_lagrangian",&(E->control.AUG_lagrangian),"0.0,0.0,1e50"); 
    input_std_precision("particle_perturb",&(E->control.particle_perturb),"0.0,0.0,1.0"); 
    input_std_precision("director_perturb",&(E->control.director_perturb),"0.0,0.0,1.0"); 
    input_std_precision("director_perturbkx",&(E->control.director_perturbkx),"0.0"); 
 
				/* general mesh structure */
    
    input_boolean("verbose",&(E->control.verbose),"off");
    input_boolean("see_convergence",&(E->control.print_convergence),"off");
    input_boolean("COMPRESS",&(E->control.COMPRESS),"on");
    input_string("gzip",E->control.gzip,"/usr/local/bin/gzip"); /*DAS: 17-01-03 */
    input_boolean("Deformation",&(E->control.deformation),"off");
    input_boolean("ELASTICITY",&(E->control.ELASTICITY),"off");
    input_boolean("ORTHOTROPY",&(E->control.ORTHOTROPY),"off");

  
#if defined (CHEM_TRANS_MODULE_INSTALLED)

    input_boolean("CHEM_TRANS",&(E->control.CHEM_TRANS),"off");

    if(E->control.CHEM_TRANS) {
      E->control.ORTHOTROPY=1;
    }

#endif


    input_boolean("periodicx",&(E->mesh.periodic_x),"off");
    input_boolean("periodic_rm_vx",&(E->mesh.periodic_rm_vx),"off");
    input_boolean("periodicy",&(E->mesh.periodic_y),"off");
    input_boolean("periodic_rm_vy",&(E->mesh.periodic_rm_vy),"off"); /*RAA, 8/10/01 */
   
    input_std_precision("toptbcval",&(E->control.TBCtopval),"0.0");
    input_std_precision("bottbcval",&(E->control.TBCbotval),"1.0");

    input_std_precision("BCmoveX0v",&(E->mesh.BCvelocityX0),"0.0");
    input_std_precision("BCmoveX1v",&(E->mesh.BCvelocityX1),"0.0");
    input_std_precision("BCmoveZ0v",&(E->mesh.BCvelocityZ0),"0.0");
    input_std_precision("BCmoveZ1v",&(E->mesh.BCvelocityZ1),"0.0");
    input_std_precision("BCmoveY0v",&(E->mesh.BCvelocityY0),"0.0");
    input_std_precision("BCmoveY1v",&(E->mesh.BCvelocityY1),"0.0");
    input_std_precision("BCmoveX0v1",&(E->mesh.BCvelocity1X0),"0.0");
    input_std_precision("BCmoveX1v1",&(E->mesh.BCvelocity1X1),"0.0");
    input_std_precision("BCmoveZ0v1",&(E->mesh.BCvelocity1Z0),"0.0");
    input_std_precision("BCmoveZ1v1",&(E->mesh.BCvelocity1Z1),"0.0");
    input_std_precision("BCmoveY0v1",&(E->mesh.BCvelocity1Y0),"0.0");
    input_std_precision("BCmoveY1v1",&(E->mesh.BCvelocity1Y1),"0.0");
    input_std_precision("BCmoveX0f",&(E->mesh.BCvelocityomegaX0),"0.0");
    input_std_precision("BCmoveX1f",&(E->mesh.BCvelocityomegaX1),"0.0");
    input_std_precision("BCmoveZ0f",&(E->mesh.BCvelocityomegaZ0),"0.0");
    input_std_precision("BCmoveZ1f",&(E->mesh.BCvelocityomegaZ1),"0.0");
    input_std_precision("BCmoveY0f",&(E->mesh.BCvelocityomegaY0),"0.0");
    input_std_precision("BCmoveY1f",&(E->mesh.BCvelocityomegaY1),"0.0");
    input_std_precision("BCmoveX0p",&(E->mesh.BCvelocityphaseX0),"0.0");
    input_std_precision("BCmoveX1p",&(E->mesh.BCvelocityphaseX1),"0.0");
    input_std_precision("BCmoveZ0p",&(E->mesh.BCvelocityphaseZ0),"0.0");
    input_std_precision("BCmoveZ1p",&(E->mesh.BCvelocityphaseZ1),"0.0");
    input_std_precision("BCmoveY0p",&(E->mesh.BCvelocityphaseY0),"0.0");
    input_std_precision("BCmoveY1p",&(E->mesh.BCvelocityphaseY1),"0.0");
    input_std_precision("BCmoveX0n",&(E->mesh.BCvelocitypowerX0),"1.0");
    input_std_precision("BCmoveX1n",&(E->mesh.BCvelocitypowerX1),"1.0");
    input_std_precision("BCmoveZ0n",&(E->mesh.BCvelocitypowerZ0),"1.0");
    input_std_precision("BCmoveZ1n",&(E->mesh.BCvelocitypowerZ1),"1.0");
    input_std_precision("BCmoveY0n",&(E->mesh.BCvelocitypowerY0),"1.0");
    input_std_precision("BCmoveY1n",&(E->mesh.BCvelocitypowerY1),"1.0");
    input_std_precision("BCmoveX0a",&(E->mesh.BCvelocityaccX0),"0.0");
    input_std_precision("BCmoveX1a",&(E->mesh.BCvelocityaccX1),"0.0");
    input_std_precision("BCmoveZ0a",&(E->mesh.BCvelocityaccZ0),"0.0");
    input_std_precision("BCmoveZ1a",&(E->mesh.BCvelocityaccZ1),"0.0");
    input_std_precision("BCmoveY0a",&(E->mesh.BCvelocityaccY0),"0.0");
    input_std_precision("BCmoveY1a",&(E->mesh.BCvelocityaccY1),"0.0");
    input_boolean("BCmove_with_time_interval",&(E->mesh.BCvelocity_time_interval),"off"); /*RAA, 12/02/02 */
    input_int("BCmove_time_or_increment",&(E->mesh.BCvelocity_time_or_inc),"0"); /*RAA, 12/02/02 */
     /*RAA: treat increments as std precision in following lines, for simplicity, but not finished yet*/
     /*RAA: below I set up the arrays for 2 load sequences only, in X and Y, but, max-load-seq is 5, I think*/
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ X0- load sequence #1*/
    input_std_precision("BCmoveX0v_startt0_1",&(E->mesh.BCvelocityX0_time[0][0][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht1_1",&(E->mesh.BCvelocityX0_time[1][1][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_startt1_1",&(E->mesh.BCvelocityX0_time[0][1][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht2_1",&(E->mesh.BCvelocityX0_time[1][2][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_startt2_1",&(E->mesh.BCvelocityX0_time[0][2][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht3_1",&(E->mesh.BCvelocityX0_time[1][3][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt0_1",&(E->mesh.BCvelocityX0_const[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet0_1",&(E->mesh.BCvelocityX0_slope[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt1_1",&(E->mesh.BCvelocityX0_const[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet1_1",&(E->mesh.BCvelocityX0_slope[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt2_1",&(E->mesh.BCvelocityX0_const[2][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet2_1",&(E->mesh.BCvelocityX0_slope[2][1]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ X1 - load sequence #1*/
    input_std_precision("BCmoveX1v_startt0_1",&(E->mesh.BCvelocityX1_time[0][0][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht1_1",&(E->mesh.BCvelocityX1_time[1][1][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_startt1_1",&(E->mesh.BCvelocityX1_time[0][1][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht2_1",&(E->mesh.BCvelocityX1_time[1][2][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_startt2_1",&(E->mesh.BCvelocityX1_time[0][2][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht3_1",&(E->mesh.BCvelocityX1_time[1][3][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt0_1",&(E->mesh.BCvelocityX1_const[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet0_1",&(E->mesh.BCvelocityX1_slope[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt1_1",&(E->mesh.BCvelocityX1_const[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet1_1",&(E->mesh.BCvelocityX1_slope[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt2_1",&(E->mesh.BCvelocityX1_const[2][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet2_1",&(E->mesh.BCvelocityX1_slope[2][1]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ X0- load sequence #2*/
    input_std_precision("BCmoveX0v_startt0_2",&(E->mesh.BCvelocityX0_time[0][0][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht1_2",&(E->mesh.BCvelocityX0_time[1][1][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_startt1_2",&(E->mesh.BCvelocityX0_time[0][1][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht2_2",&(E->mesh.BCvelocityX0_time[1][2][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_startt2_2",&(E->mesh.BCvelocityX0_time[0][2][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_finisht3_2",&(E->mesh.BCvelocityX0_time[1][3][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt0_2",&(E->mesh.BCvelocityX0_const[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet0_2",&(E->mesh.BCvelocityX0_slope[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt1_2",&(E->mesh.BCvelocityX0_const[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet1_2",&(E->mesh.BCvelocityX0_slope[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_constt2_2",&(E->mesh.BCvelocityX0_const[2][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX0v_slopet2_2",&(E->mesh.BCvelocityX0_slope[2][2]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ X1 - load sequence #2*/
    input_std_precision("BCmoveX1v_startt0_2",&(E->mesh.BCvelocityX1_time[0][0][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht1_2",&(E->mesh.BCvelocityX1_time[1][1][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_startt1_2",&(E->mesh.BCvelocityX1_time[0][1][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht2_2",&(E->mesh.BCvelocityX1_time[1][2][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_startt2_2",&(E->mesh.BCvelocityX1_time[0][2][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_finisht3_2",&(E->mesh.BCvelocityX1_time[1][3][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt0_2",&(E->mesh.BCvelocityX1_const[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet0_2",&(E->mesh.BCvelocityX1_slope[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt1_2",&(E->mesh.BCvelocityX1_const[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet1_2",&(E->mesh.BCvelocityX1_slope[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_constt2_2",&(E->mesh.BCvelocityX1_const[2][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveX1v_slopet2_2",&(E->mesh.BCvelocityX1_slope[2][2]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ Y0 - load sequence #1*/
    input_std_precision("BCmoveY0v_startt0_1",&(E->mesh.BCvelocityY0_time[0][0][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht1_1",&(E->mesh.BCvelocityY0_time[1][1][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_startt1_1",&(E->mesh.BCvelocityY0_time[0][1][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht2_1",&(E->mesh.BCvelocityY0_time[1][2][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_startt2_1",&(E->mesh.BCvelocityY0_time[0][2][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht3_1",&(E->mesh.BCvelocityY0_time[1][3][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt0_1",&(E->mesh.BCvelocityY0_const[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet0_1",&(E->mesh.BCvelocityY0_slope[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt1_1",&(E->mesh.BCvelocityY0_const[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet1_1",&(E->mesh.BCvelocityY0_slope[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt2_1",&(E->mesh.BCvelocityY0_const[2][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet2_1",&(E->mesh.BCvelocityY0_slope[2][1]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ Y1 - load sequence #1*/
    input_std_precision("BCmoveY1v_startt0_1",&(E->mesh.BCvelocityY1_time[0][0][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht1_1",&(E->mesh.BCvelocityY1_time[1][1][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_startt1_1",&(E->mesh.BCvelocityY1_time[0][1][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht2_1",&(E->mesh.BCvelocityY1_time[1][2][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_startt2_1",&(E->mesh.BCvelocityY1_time[0][2][1]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht3_1",&(E->mesh.BCvelocityY1_time[1][3][1]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt0_1",&(E->mesh.BCvelocityY1_const[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet0_1",&(E->mesh.BCvelocityY1_slope[0][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt1_1",&(E->mesh.BCvelocityY1_const[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet1_1",&(E->mesh.BCvelocityY1_slope[1][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt2_1",&(E->mesh.BCvelocityY1_const[2][1]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet2_1",&(E->mesh.BCvelocityY1_slope[2][1]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ Y0 - load sequence #2*/
    input_std_precision("BCmoveY0v_startt0_2",&(E->mesh.BCvelocityY0_time[0][0][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht1_2",&(E->mesh.BCvelocityY0_time[1][1][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_startt1_2",&(E->mesh.BCvelocityY0_time[0][1][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht2_2",&(E->mesh.BCvelocityY0_time[1][2][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_startt2_2",&(E->mesh.BCvelocityY0_time[0][2][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_finisht3_2",&(E->mesh.BCvelocityY0_time[1][3][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt0_2",&(E->mesh.BCvelocityY0_const[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet0_2",&(E->mesh.BCvelocityY0_slope[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt1_2",&(E->mesh.BCvelocityY0_const[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet1_2",&(E->mesh.BCvelocityY0_slope[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_constt2_2",&(E->mesh.BCvelocityY0_const[2][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY0v_slopet2_2",&(E->mesh.BCvelocityY0_slope[2][2]),"0.0");     /*RAA, 19/02/02 */
     /*RAA: start velocity vs time info for ramp (trapezoid) function w/ Y1 - load sequence #2*/
    input_std_precision("BCmoveY1v_startt0_2",&(E->mesh.BCvelocityY1_time[0][0][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht1_2",&(E->mesh.BCvelocityY1_time[1][1][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_startt1_2",&(E->mesh.BCvelocityY1_time[0][1][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht2_2",&(E->mesh.BCvelocityY1_time[1][2][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_startt2_2",&(E->mesh.BCvelocityY1_time[0][2][2]),"0.0");    /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_finisht3_2",&(E->mesh.BCvelocityY1_time[1][3][2]),"1.0e32");/*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt0_2",&(E->mesh.BCvelocityY1_const[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet0_2",&(E->mesh.BCvelocityY1_slope[0][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt1_2",&(E->mesh.BCvelocityY1_const[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet1_2",&(E->mesh.BCvelocityY1_slope[1][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_constt2_2",&(E->mesh.BCvelocityY1_const[2][2]),"0.0");     /*RAA, 19/02/02 */
    input_std_precision("BCmoveY1v_slopet2_2",&(E->mesh.BCvelocityY1_slope[2][2]),"0.0");     /*RAA, 19/02/02 */


    input_std_precision("offsetx",&(E->mesh.layer0[1]),"0.0");
    input_std_precision("offsetz",&(E->mesh.layer0[2]),"0.0");
    input_std_precision("offsety",&(E->mesh.layer0[3]),"0.0");
    input_std_precision("dimenx",&(E->mesh.layer1[1]),"nodefault");
    input_std_precision("dimenz",&(E->mesh.layer1[2]),"nodefault");
    input_std_precision("dimeny",&(E->mesh.layer1[3]),"nodefault");
 
    input_std_precision("meshx1",&(E->mesh.layer0[1]),"0.0");
    input_std_precision("meshz1",&(E->mesh.layer0[2]),"0.0");
    input_std_precision("meshy1",&(E->mesh.layer0[3]),"0.0");
    input_std_precision("meshx2",&(E->mesh.layer1[1]),"nodefault");
    input_std_precision("meshz2",&(E->mesh.layer1[2]),"nodefault");
    input_std_precision("meshy2",&(E->mesh.layer1[3]),"nodefault");

    input_std_precision("phi0",&(E->mesh.layer0[1]),"nodefault");
    input_std_precision("phi1",&(E->mesh.layer1[1]),"nodefault");
    input_std_precision("r0",&(E->mesh.layer1[2]),"nodefault");
    input_std_precision("r1",&(E->mesh.layer0[2]),"nodefault");
    input_std_precision("theta0",&(E->mesh.layer0[3]),"nodefault");
    input_std_precision("theta1",&(E->mesh.layer1[3]),"nodefault");
    
    input_int("storage_timesteps",&(E->control.storage_timesteps),"50");
    input_int("checkpt_timesteps",&(E->control.checkpt_timesteps),"10");
    input_std_precision("storage_timing",&(E->control.storage_timing),"0.0");
    input_std_precision("checkpt_timing",&(E->control.checkpt_timing),"0.0");

    input_int("PPM_height",&(E->control.PPM_height),"256");
    input_std_precision("PPM_aspect",&(E->control.PPM_aspect),"1.0");
    input_int("PPM_files",&(E->control.PPM_files),"1");
    input_std_precision("PPM_show_strain",&(E->control.PPM_strain),"0.0");
    input_int_vector("PPM_surf_plot",E->control.PPM_files,E->control.PPM_surf_plot,1);
    input_int_vector("PPM_coloring",E->control.PPM_files,E->control.PPM_coloring,1);
    input_int_vector("PPM_coloring_autorange",E->control.PPM_files,E->control.PPM_color_auto,1);
    input_std_precision_vector("PPM_coloring_min",E->control.PPM_files,E->control.PPM_color_min,0.0);
    input_std_precision_vector("PPM_coloring_max",E->control.PPM_files,E->control.PPM_color_max,1.0);
    input_int_vector("PPM_surf_autorange",E->control.PPM_files,E->control.PPM_surf_auto,1);
    input_std_precision_vector("PPM_surf_min",E->control.PPM_files,E->control.PPM_surf_min,0.0);
    input_std_precision_vector("PPM_surf_max",E->control.PPM_files,E->control.PPM_surf_max,1.0);
    input_std_precision("PPM_3D_slice_thickness",&(E->control.PPM_3D_delt),"0.1"); /*RAA: 23/01/03, added this*/

    input_hgr_precision("accuracy",&(E->control.accuracy),"1.0e-4,0.0,1.0");
    input_std_precision("delta_accuracy_factor",&(E->control.delta_accuracy_factor),"1.0,0.01,10.0");

    input_int("mg_cycle",&(E->control.mg_cycle),"1,0,nomax");
    input_int("vel_relaxations",&(E->control.vsolver_relaxations),"2,0,nomax");
    
    input_int("viterations",&(E->control.max_vel_iterations),"251,0,nomax");
    input_int("piterations",&(E->control.p_iterations),"100,0,nomax");

  /* data section */ 

    input_std_precision("gravacc",&(E->data.grav_acc),"9.81");
    input_std_precision("gravtheta",&(E->data.grav_theta),"0.0");

    (E->problem_settings)(E);
    return; 
}
