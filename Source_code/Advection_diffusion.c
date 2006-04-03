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


/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.
*/

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "function_prototypes.h"
#include <stdlib.h>

extern int Emergency_stop;

/*struct el { higher_precision gpt[9]; };*/  /*RAA: structure is not used*/

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */

void advection_diffusion_parameters (
 struct All_variables *E
 )
{ 
 
  int i;
  
  /* Set intial values, defaults & read parameters*/
  
  E->advection.temp_iterations = 2; /* petrov-galerkin iterations: minimum value. */
  E->advection.total_timesteps = 1; 
  E->advection.sub_iterations = 1;
  E->advection.last_sub_iterations = 1;
  E->advection.gamma = 0.5;
  
  input_boolean("ADV",&(E->advection.ADVECTION),"on");
  input_int("minstep",&(E->advection.min_timesteps),"1");
  input_int("maxstep",&(E->advection.max_timesteps),"1000");
  input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000");
  input_std_precision("finetunedt",&(E->advection.fine_tune_dt),"0.9");
  input_std_precision("fixed_timestep",&(E->advection.fixed_timestep),"0.0");
  input_std_precision("elastic_timestep",&(E->advection.elastic_timestep),"0.0");
  input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,2,nomax");
  input_std_precision("maxadvtime",&(E->advection.max_elapsed_time),"1.0e18"); 

  /* Set an initial timestep for use in elasticity
     calculations */

 /*E->advection.timestep=1.0e-16 + E->advection.fixed_timestep;*/
 if(E->advection.fixed_timestep == 0.0)
    E->advection.timestep = 1.e-16 ;
  else
    E->advection.timestep = E->advection.fixed_timestep ;


 if(E->control.ELASTICITY && E->advection.elastic_timestep <= 0.0) {
   fprintf(stderr,"Warning, no elastic timestep specified, turning off elasticity option\n");
   fprintf(E->fp,"Warning, no elastic timestep specified, turning off elasticity option\n");
   
   E->control.ELASTICITY = 0;
 }

 return;
}

void advection_diffusion_allocate_memory
(
 struct All_variables *E
 )
{
  int i;
  
  E->Tdot = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  for(i=1;i<=E->mesh.nno;i++) 
    E->Tdot[i]=0.0;
  
  return;
}

void PG_timestep(
		 struct All_variables *E
		 )
{
  void predictor();
  void corrector();
  void std_timestep();
  void temperatures_conform_bcs();
  
  int get_eq_phase();
  int get_v_estimate();
  int i,j,psc_pass,count,steps,m;
  int keep_going,previous_phase,metastable;
  
  struct SOURCES Qnone;
  
  standard_precision *DTdot,*Fold,*dF,vsselfdot();
  standard_precision *mem1,*mem2;
  standard_precision deltaF,Fmag,Fmag2;
  standard_precision *DSTNdot;
  standard_precision CPU_time(),time;
  
  standard_precision ph_b_distance;
  
  static int loops_since_new_eta = 0;
  static int been_here = 0;
  const int dims = E->mesh.nsd;
  Qnone.number=0;
  
  DTdot = (standard_precision *)Malloc0(E->mesh.fnodal_malloc_size);
  DSTNdot = (standard_precision *)Malloc0(E->mesh.fnodal_malloc_size);
  
  if(E->control.verbose)
    fprintf(stderr,"Advection\n");
  
  if (been_here++ ==0)
    E->advection.timesteps=0;
  
  /* place for stress updating in explicit scheme */
  /* fprintf(stderr,"HERE IS TRACER.PT (13)   %g \n",E->tracer.Pt[436]) ;	 */
  if(E->control.ELASTICITY || E->control.SHEAR_HEATING) {
    if(E->control.verbose)
      fprintf(stderr,"Update stress history \n");
    stress_update(E);
  }
 /* fprintf(stderr,"HERE IS TRACER.PT (14)   %g \n",E->tracer.Pt[436]) ;	*/
/* The positions to be calculated actually reflect the beginning of the next timestep. */
  E->advection.timesteps++; 

  fprintf(stderr,"CONE: Going into std_timestep %g  \n",E->advection.timestep);

  E->advection.previous_timestep_2 = E->advection.previous_timestep + E->advection.timestep;
  E->advection.previous_timestep = E->advection.timestep; 
  std_timestep(E);
  fprintf(stderr,"CONE: Out of std_timestep %g  \n",E->advection.timestep);

  time=CPU_time();
  /* keep_going = get_v_estimate(E); */
   
  /* Diffusion (thermal) */

#if 1 
  for(count=1;count<=E->advection.diff_ratio;count++) {
    predictor(E,E->T,E->Tdot,E->node,OFFSIDE  | TBD);
    
    for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
      pg_solver(E,E->T,E->Tdot,DTdot,E->V,E->convection.heat_sources,1.0,1,E->TB,
		E->node,OFFSIDE  | TBD, FBZ);
      corrector(E,E->T,E->Tdot,DTdot,E->node,OFFSIDE  | TBD); 
    }

    if(E->control.verbose)
      fprintf(stderr,"Sub timestep %d - %g / %g\n",count,E->advection.timestep_diff,E->advection.timestep_adv); 
  }
#endif    

  fprintf(stderr,"TIME UPDATED HERE: monitor_time: %g adv_time: %g \n",E->monitor.elapsed_time,E->advection.timestep);
  E->advection.total_timesteps++; /* i.e. the computed location will apply to the next timestep */
  E->monitor.elapsed_time += E->advection.timestep; 
  count++; 
  fprintf(stderr,"TIME UPDATE DONE: monitor_time: %g adv_time: %g \n",E->monitor.elapsed_time,E->advection.timestep);

  if(E->mesh.periodic_x || E->mesh.periodic_y)  {
    flogical_mesh_to_real(E,E->T,E->mesh.levmax);
    flogical_mesh_to_real(E,E->Tdot,E->mesh.levmax);
  }

  /* TRACER ADVECTION ... */
  /* Update volumetric strain at tracers etc.*/
  /*fprintf(stderr,"HERE IS TRACER.PT (0)   %g \n",E->tracer.Pt[436]) ;	*/
  if(E->control.verbose)
    fprintf(stderr,"Particle based pressure\n");
  tracer_time_dep_terms(E); 

  for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
    E->tracer.UU1[i] =  E->tracer.U1[i];
    E->tracer.UU2[i] =  E->tracer.U2[i];
    if(3==dims)  /*RAA: 04/04/02, added this part*/
      E->tracer.UU3[i] =  E->tracer.U3[i];
  }

  /*RAA:5/4/01, nodes_to_tracers call for 3D added below*/
  nodes_to_tracers(E,E->V[1],E->tracer.U1,E->mesh.levmax);
  nodes_to_tracers(E,E->V[2],E->tracer.U2,E->mesh.levmax);
  if(3==dims)  
    nodes_to_tracers(E,E->V[3],E->tracer.U3,E->mesh.levmax);
 
  
#if defined (CHEM_TRANS_MODULE_INSTALLED)
  if(E->control.CHEM_TRANS) {
    if(2==dims) { /*RAA: added this distinction*/	  
      update_sigma_n(E);
      update_porosity(E);
    }
    else if(3==dims) { /*RAA: added this bit*/	  
      fprintf(stderr,"3D chem transport not properly coded. Adios, muchacho.\n");
      exit(1);
    }
    
  }
#endif



  /* Advection of tracers */

 if(E->control.verbose)
    fprintf(stderr,"Advection of tracers \n");

  tracer_advection(E,E->advection.timestep);
  passive_tracer_advection(E,E->advection.timestep);


  if(E->control.verbose)
    fprintf(stderr,"Particle phase stuff\n");
  
  /* Compute which PHASE each tracer is in */
  
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(E->tracer.Phases[E->tracer.property_group[m]] <= 1)
      continue;
    
    /* 1. Before updating, is the tracer in metastable state ?  */
    
    if(E->tracer.Current_phase[m] != E->tracer.Equilibrium_phase[m])
      metastable = 1;
    else
      metastable = 0;
    
    /* 2. Now find new equilibrium phase, and store previous phase */
    
    E->tracer.Equilibrium_phase[m] = get_eq_phase(E,m);

    if(E->monitor.solution_cycles == 1)
      previous_phase = E->tracer.Equilibrium_phase[m];
    else
      previous_phase = E->tracer.Current_phase[m];
    
    /* 3. If particle can change phase, allow it to do so */
    
    if(((E->monitor.solution_cycles == 1) || /* No knowledge of previous phase if first step */
	(E->tracer.T[m] >
	 E->tracer.TBlock[E->tracer.property_group[m]*MAX_MATERIAL_PHASES+E->tracer.Current_phase[m]]))) {
      E->tracer.Current_phase[m] = E->tracer.Equilibrium_phase[m];
    }    
    
    /* 4. If particle changes, reset this clock */
    
    if(previous_phase != E->tracer.Current_phase[m])
      E->tracer.time_since_phase_change[m] = 0.0;
    else
      E->tracer.time_since_phase_change[m] += E->advection.timestep;
    
    /* 5. Update grainsize for particles crossing boundary ... depends
       on whether they cross equilibrium or metastable boundary */
    
    if(E->viscosity.GRAINSIZE)
      phase_boundary_grain_size(E,m,metastable);
  }
  
  if(E->viscosity.GRAINSIZE)
    grow_grains(E);

  if(E->control.ORTHOTROPY) 
    rotate_director(E);
 

  /* Estimate new velocity solution */
  /* get_v_estimate(E); */ 
  
  temperatures_conform_bcs(E,E->T); 
  nodes_to_tracers(E,E->T,E->tracer.T,E->mesh.levmax);
  
    /* Adjust mesh if required and resize elements */  
  mesh_update(E,E->advection.timestep,E->monitor.elapsed_time);

  /* Update tracer element information */

  get_tracer_elts_and_weights(E,1,E->tracer.NUM_TRACERS); 

    
  E->advection.last_sub_iterations = count;
  
  if(((E->advection.total_timesteps < E->advection.max_total_timesteps) &&
      (E->advection.timesteps < E->advection.max_timesteps) && 
      (E->monitor.elapsed_time < E->advection.max_elapsed_time) ) ||
     (E->advection.total_timesteps < E->advection.min_timesteps) )
    E->control.keep_going = 1;
  else {
    E->control.keep_going = 0;
    if(E->control.verbose) {
      if(E->advection.total_timesteps >= E->advection.max_total_timesteps)
	fprintf(stderr, "DAS: ** maxtotstep reached, stopping.\n");
      else if(E->advection.timesteps >= E->advection.max_timesteps)
	fprintf(stderr, "DAS: ** maxstep reached, stopping.\n");
      else if(E->monitor.elapsed_time >= E->advection.max_elapsed_time)
	fprintf(stderr, "DAS: ** maxadvtime reached, stopping.\n");
      else
	fprintf(stderr, "DAS: ** ARGH! Something unusual just happened!\n");
    }
  }
  
  if(E->control.verbose)
    fprintf(stderr,"Advection ... done\n");
  
  free((void *) DTdot);   
  free((void *) DSTNdot);       /* free memory for vel solver */
  return;  
}


/* Update nodal values of temperature using
   the earlier positions of the nodes */

void nodal_semi_lagrange_advection(
     struct All_variables *E,
     standard_precision *T,
     standard_precision *T1,
     standard_precision timestep 
)  /*RAA: function argument syntax corrected on 4/4/01 */
{
  int node,el;
  int i,j,k;

  standard_precision nodal_interpolated_value();

  standard_precision vx1,vz1,vy1;
  standard_precision vx2,vz2,vy2;
  standard_precision xx1,zz1,yy1;
  standard_precision xx2,zz2,yy2;
  standard_precision xx3,zz3,yy3;
  standard_precision kx1,kz1,ky1;
  standard_precision kx2,kz2,ky2;

  /* For each nodal point, find the advected position */

  for(node=1;node<=E->mesh.nno;node++) {
    vx1 = E->V[1][node];
    vz1 = E->V[2][node];
    xx1 = E->x[1][node];
    zz1 = E->x[2][node];

    kx1 = timestep * vx1;
    kz1 = timestep * vz1;

    /* Midpt ... now get its velocity */
    xx2 = xx1 + kx1 * 0.5;
    zz2 = zz1 + kz1 * 0.5;

    vx2 = nodal_interpolated_value(E,E->V[1],xx2,zz2,NULL);
    vz2 = nodal_interpolated_value(E,E->V[2],xx2,zz2,NULL);

    kx2 = timestep * vx2;
    kz2 = timestep * vz2;

    /* Endpt ... now get its temperature */
 
    xx3 = xx1 + kx2;
    zz3 = zz1 + kz2;

    T1[node] = nodal_interpolated_value(E,T,xx3,zz3,NULL);
    /* fprintf(stderr,"Node %d (%g,%g) used to be at %g,%g where T = %g\n",
       node,xx1,zz1,xx3,zz3,T1[node]); */
  }
}

/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(
  struct All_variables *E,
  standard_precision *field,
  standard_precision *fielddot,
  int *INFO,
  int MASK
)
{ 
    int node;
    standard_precision multiplier;

   multiplier = (1.0-E->advection.gamma) * E->advection.timestep_diff;

    for(node=1;node<=E->mesh.nno;node++)  {
      if(INFO != (int *) NULL) { 
	if( !(INFO[node] & MASK)) 
	  field[node] += multiplier * fielddot[node] ;
      }
      else
	 field[node] += multiplier * fielddot[node] ;

      fielddot[node] = 0.0;
    }
   return; 
}


void corrector(
  struct All_variables *E,
  standard_precision *field,
  standard_precision *fielddot,
  standard_precision *Dfielddot,
  int * INFO,
  int MASK
)

{
  int node;
   standard_precision multiplier;

   multiplier = E->advection.gamma * E->advection.timestep_diff;

   for(node=1;node<=E->mesh.nno;node++) {
     if(INFO != (int *) NULL) { 
       if( !(INFO[node] & MASK))  
	 field[node] += multiplier * Dfielddot[node];
     }
     else
       field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
   }
  return;  
 }

/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style 
   vector optimizations etc and one optimized for 
   workstations.
   =================================================== */

void pg_solver(
  struct All_variables *E,
  standard_precision *T,
  standard_precision *Tdot,
  standard_precision *DTdot,
  standard_precision **V,
  struct SOURCES Q0,
  standard_precision diff,
  int bc,
  standard_precision *TBC,
  unsigned int *FLAGS,
  unsigned int OFFSIDE_MASK,
  unsigned int FLUX_MASK
)
{
    void get_global_shape_fn();
    void pg_element_residual_t();

    int el,a,i,a1;
    higher_precision Eres[9];  /* correction to the (scalar) Tdot field */
    /*RAA: 12/07/02, 3 new variables here from C. Wijns advection fix*/
    standard_precision dist,nodeV;
    standard_precision *realV[4];

    standard_precision vel0,vel1; /*RAA: 23/9/02, 2 new variables here*/
    standard_precision find_section_of_vel_fn(); /*  RAA: 20/09/02, new function*/

    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
 
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];

    /*RAA - as per C.W. - we don't need to worry about rotational components*/
    for(i=1;i<=dims;i++)
      realV[i] = (standard_precision *)Malloc0((E->mesh.nno+1)*sizeof(standard_precision));
  
    for(i=1;i<=E->mesh.nno;i++)
	DTdot[i] = 0.0;


  /*RAA (C.W.) Calculate a true nodal advection velocity subtracting the (constant) mesh velocity,
    w/ time_interval=off & w/ stationary walls, nodeVs will be 0, so this should be fine*/
    if(!E->mesh.BCvelocity_time_interval) {  
      for(i=1;i<=E->mesh.nno;i++) {
        /* X-velocity*/
        dist = (E->x[1][i] - E->mesh.layer0[1]) / (E->mesh.layer1[1] - E->mesh.layer0[1]);
        nodeV = E->mesh.BCvelocityX0 + (E->mesh.BCvelocityX1 - E->mesh.BCvelocityX0) * dist;
        realV[1][i] = V[1][i] - nodeV;
        /* Z-velocity*/
        dist = (E->x[2][i] - E->mesh.layer0[2]) / (E->mesh.layer1[2] - E->mesh.layer0[2]);
        nodeV = E->mesh.BCvelocityZ0 + (E->mesh.BCvelocityZ1 - E->mesh.BCvelocityZ0) * dist;
        realV[2][i] = V[2][i] - nodeV;
        /* Y-velocity*/
        if (dims==3) {
          dist = (E->x[3][i] - E->mesh.layer0[3]) / (E->mesh.layer1[3] - E->mesh.layer0[3]);
          nodeV = E->mesh.BCvelocityY0 + (E->mesh.BCvelocityY1 - E->mesh.BCvelocityY0) * dist;
          realV[3][i] = V[3][i] - nodeV;
        }
      }    
    }
#if 1
/*--------------------------------------*/
    else if(E->mesh.BCvelocity_time_interval) { /*RAA: vel bcs depend linearly on time*/ 
      for(i=1;i<=E->mesh.nno;i++) {
  	vel0 = 0.0;
  	vel1 = 0.0;

        /* X-velocity*/
        dist = (E->x[1][i] - E->mesh.layer0[1]) / (E->mesh.layer1[1] - E->mesh.layer0[1]);
        find_section_of_vel_fn(E,1,E->monitor.elapsed_time,&vel0,&vel1);
	/*RAA:  has the info below already been updated??? apparently not*/
/*      nodeV = E->mesh.BCvelocityX0 + (E->mesh.BCvelocityX1 - E->mesh.BCvelocityX0) * dist; */
        nodeV = vel0 + (vel1 - vel0) * dist;
        realV[1][i] = V[1][i] - nodeV;

        /* Z-velocity*/
        dist = (E->x[2][i] - E->mesh.layer0[2]) / (E->mesh.layer1[2] - E->mesh.layer0[2]);
  	find_section_of_vel_fn(E,2,E->monitor.elapsed_time,&vel0,&vel1);
        nodeV = vel0 + (vel1 - vel0) * dist;
        realV[2][i] = V[2][i] - nodeV;

        /* Y-velocity*/
        if (dims==3) {
          dist = (E->x[3][i] - E->mesh.layer0[3]) / (E->mesh.layer1[3] - E->mesh.layer0[3]);
          find_section_of_vel_fn(E,3,E->monitor.elapsed_time,&vel0,&vel1);
          nodeV = vel0 + (vel1 - vel0) * dist;
       /*   fprintf(E->fp1,"nodeV   vel0,  vel1, dist  %g  %g  %g %g\n",nodeV,vel0,vel1,dist); */
          realV[3][i] = V[3][i] - nodeV;
        }
      }    
    }
/*--------------------------------------*/
#endif

    for(el=1;el<=E->mesh.nel;el++) {
   
      get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,E->mesh.levmax);
      /*RAA - here's the proper call to pg_element_residual w/ realV, as per C.W.'s bug fix */
      pg_element_residual(E,el,GNx,dOmega,realV,T,Tdot,Q0,Eres,1.0,TBC,FLAGS,OFFSIDE_MASK,FLUX_MASK); 

      for(a=1;a<=ends;a++) {
	a1 = E->ien[el].node[a];
	DTdot[a1] += Eres[a]; 
      }
      
    } /* next element */
  
    for(i=1;i<=E->mesh.nno;i++) {
	if(E->node[i] & ( OFFSIDE  )) 
	  continue;
	DTdot[i] /= E->mass[i];         /* lumped mass matrix */
    }

    for(i=1;i<=dims;i++) {
        free((void *) realV[i]);   
    }
    printf("Finished pg_solver, free'd up realV \n");
    return;    
}

/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term. Calculated
   using Petrov-Galerkin Shape functions.
   =========================================  */

void pg_element_residual(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision **vel,
  standard_precision *TT1,
  standard_precision *TT1dot,
  struct SOURCES QQ1,
  higher_precision Eres1[9],
  standard_precision diff1,
  standard_precision *BCC,
  unsigned int *FLAGS1,
  int OFFSIDE_MASK1,
  int FLUX_MASK1
)
{
  int i,j,a,k,m,node,nodes[5],d,aid,back_front,onedfns;
  higher_precision v1[9],v2[9],v3[9];
  higher_precision uc1,uc2,uc3;
  higher_precision size1,size2,size3;

  higher_precision T1xsi1,T1xsi2,T1xsi3;
  higher_precision T1ah1,T1ah2,T1ah3;
  higher_precision Q1,Qi[9];
  higher_precision dT1[9];
  higher_precision T1x1[9],T1x2[9],T1x3[9];
  higher_precision T1,DT1,Ti;

  higher_precision Vsum;
  higher_precision T1xsisum,vT1sum;
  higher_precision unorm;
  higher_precision sfn,pgsfn;
    
  struct Shape_function1 GM;
  struct Shape_function1_dA dGamma;
  
  standard_precision grain_growth_source_term();
  standard_precision diff1_1,QtE,vol;
  void get_global_1d_shape_fn();
    
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const higher_precision recipends=1.0/(higher_precision)ends;
  
  /* Since we're only using this for Temperature
     these days, not chemistry or anything non-diffusive */

   diff1=0.0;
   QtE=0.0;
   vol = 0.0;
   for(j=0;j<E->tracer.tr_in_element_number[E->mesh.levmax][el];j++) {
     m = E->tracer.tr_in_element[E->mesh.levmax][j+E->tracer.tr_in_element_offset[E->mesh.levmax][el]];
  
     diff1 += E->tracer.Therm_diff[E->tracer.property_group[m]] *  
       E->tracer.tracer_weight_fn[E->mesh.levmax][m];
     QtE += E->tracer.Qt[E->tracer.property_group[m]] *  
       E->tracer.tracer_weight_fn[E->mesh.levmax][m] /     /* If Q is supplied as W/kg */
       (E->tracer.Cp[E->tracer.property_group[m]]);

     if(E->control.SHEAR_HEATING) {
       /* fprintf(stderr,"Adding shear heating of %g for tracer %d\n",E->tracer.Hvisc[m],m); */
       
       QtE += E->tracer.Hvisc[m] * E->tracer.tracer_weight_fn[E->mesh.levmax][m];

     }
     
    /*RAA: 09/07/03 - Latent heat of melting stuff - next 2 lines from Craig's 2D routine*/
     QtE -= ((E->tracer.Latent_heat[E->tracer.property_group[m]] *E->tracer.dFdot[m]* E->tracer.T[m] * 
       E->tracer.tracer_weight_fn[E->mesh.levmax][m]) / (E->tracer.Cp[E->tracer.property_group[m]]));
     

     vol += E->tracer.tracer_weight_fn[E->mesh.levmax][m];
   }

   if(vol != 0.0) {
     diff1 /= vol;
     QtE /= vol;
   }

   diff1_1 = 1.0 / (diff1 + 1.0e-32);

   /* And this bit is unchanged */

  for(i=1;i<=vpts;i++)	{ 
    dT1[i]=0.0;
    v1[i] = T1x1[i]=  0.0;
    v2[i] = T1x2[i]=  0.0;
    v3[i] = T1x3[i]=  0.0;
  }
  
  uc1 = uc2 = uc3 = 0.0;

  size1 = (higher_precision) E->eco[el].ntl_size[1]; 
  size2 = (higher_precision) E->eco[el].ntl_size[2]; 
  size3 = (higher_precision) E->eco[el].ntl_size[3]; 
  
  Q1=0.0;
  for(i=0;i<QQ1.number;i++)
    Q1 += QQ1.Q[i] * exp(-QQ1.lambda[i] * (E->monitor.elapsed_time+QQ1.t_offset));
  for(i=1;i<=vpts;i++)
    Qi[i] = Q1;

  /* Track to see if nothing to advect */

  for(j=1;j<=ends;j++) {
    Eres1[j] = 0.0;

    node = E->ien[el].node[j];
    T1 = TT1[node];  
    if(fabs(T1) < 1.0e-16)
      T1 = 0.0;

    if(FLAGS1[node] & (OFFSIDE_MASK1))
      DT1=0.0;
    else
      DT1 = TT1dot[node];
  
    if(3==dims) {  /* Central velocity is expressed in natural coordinates and
		      used to construct the PG shape functions. Elsewhere the
		      cartesian velocity field is used */
      uc1 +=  (E->eco[el].ntl_dirns[1][1] * vel[1][node] +
	       E->eco[el].ntl_dirns[1][2] * vel[2][node] +
	       E->eco[el].ntl_dirns[1][3] * vel[3][node] ) ;
      uc2 +=  (E->eco[el].ntl_dirns[2][1] * vel[1][node] +
	       E->eco[el].ntl_dirns[2][2] * vel[2][node] +
	       E->eco[el].ntl_dirns[2][3] * vel[3][node] ) ;
      uc3 +=  (E->eco[el].ntl_dirns[3][1] * vel[1][node] +
	       E->eco[el].ntl_dirns[3][2] * vel[2][node] +
	       E->eco[el].ntl_dirns[3][3] * vel[3][node] ) ;
	
      Ti=0.0;
      for(i=1;i<=vpts;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
	dT1[i] += DT1 * sfn; 
	v1[i] +=  vel[1][node] * sfn;
	v2[i] +=  vel[2][node] * sfn;
	v3[i] +=  vel[3][node] * sfn;
     	T1x1[i] += GNx.vpt[GNVXINDEX(0,j,i)] * T1;  
	T1x2[i] += GNx.vpt[GNVXINDEX(1,j,i)] * T1; 
	T1x3[i] += GNx.vpt[GNVXINDEX(2,j,i)] * T1; 
	Ti += sfn * T1; 

      if(E->advection_source_term != NULL) 
	Qi[i]=(E->advection_source_term)(E,el,i,TT1,v1[i],v2[i],v3[i],Ti);
      }
    }
  
    else /* 2==dims */ {
      uc1 +=  (E->eco[el].ntl_dirns[1][1] * vel[1][node] +
	       E->eco[el].ntl_dirns[1][2] * vel[2][node] ) ;
      uc2 +=  (E->eco[el].ntl_dirns[2][1] * vel[1][node] +
	       E->eco[el].ntl_dirns[2][2] * vel[2][node] ) ;  
      
      Ti=0.0;
      for(i=1;i<=vpts;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
	dT1[i] += DT1 * sfn;
	v1[i] += vel[1][node] * sfn;
	v2[i] += vel[2][node] * sfn;
	Ti += sfn * T1;
	T1x1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T1; 
	T1x2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T1;   
      
	if(E->advection_source_term != NULL)
	  Qi[i]=(E->advection_source_term)(E,el,i,TT1,v1[i],v2[i],0.0,Ti);
      }
    }
  }

  /* Construct PG shape functions */ 
  
  if(3==dims) {
   if (diff1 != 0.0) {
     T1ah1 = fabs(uc1*recipends) * size1 * diff1_1; 
     T1ah2 = fabs(uc2*recipends) * size2 * diff1_1; 
     T1ah3 = fabs(uc3*recipends) * size3 * diff1_1; 
     T1xsi1 = (T1ah1 > 1.94) ?  T1ah1 - 1.0 : 0.3333333*T1ah1*T1ah1*(1.0-T1ah1*T1ah1*0.06666667); 
     T1xsi2 = (T1ah2 > 1.94) ?  T1ah2 - 1.0 : 0.3333333*T1ah2*T1ah2*(1.0-T1ah2*T1ah2*0.06666667);  
     T1xsi3 = (T1ah3 > 1.94) ?  T1ah3 - 1.0 : 0.3333333*T1ah3*T1ah3*(1.0-T1ah3*T1ah3*0.06666667);  
     T1xsisum = diff1 * (T1xsi1 + T1xsi2 + T1xsi3) * 0.518;
   }
   else {
        T1ah1 = fabs(uc1*recipends) * size1; 
	T1ah2 = fabs(uc2*recipends) * size2; 
	T1ah3 = fabs(uc3*recipends) * size3; 
	T1xsisum = (T1ah1 + T1ah2 + T1ah3) * 0.518;
   }
      for(i=1;i<=VPOINTS3D;i++) {
	unorm = 1.0 / (v1[i] * v1[i] + v2[i] * v2[i] + v3[i] * v3[i] + 1.0e-32);
	vT1sum = dT1[i] /*- Q1*/ + v1[i]*T1x1[i] + v2[i]*T1x2[i] + v3[i]*T1x3[i];

	for(j=1;j<=ENODES3D;j++) {
	  Vsum = (v1[i] * GNx.vpt[GNVXINDEX(0,j,i)] +
		  v2[i] * GNx.vpt[GNVXINDEX(1,j,i)] +
		  v3[i] * GNx.vpt[GNVXINDEX(2,j,i)] ) ;

	  pgsfn = E->N.vpt[GNVINDEX(j,i)] + T1xsisum * Vsum * unorm;

	  Eres1[j] -= dOmega.vpt[i] * (pgsfn  *  vT1sum - E->N.vpt[GNVINDEX(j,i)] * QtE +
					diff1 * (GNx.vpt[GNVXINDEX(0,j,i)]*T1x1[i] +
						 GNx.vpt[GNVXINDEX(1,j,i)]*T1x2[i] +
						 GNx.vpt[GNVXINDEX(2,j,i)]*T1x3[i] ) ); 
	}
      }
    }
  
  else  /* 2==dims */ {
    if (diff1 != 0.0) {
      T1ah1 = fabs(uc1*recipends) * size1 * diff1_1; 
      T1ah2 = fabs(uc2*recipends) * size2 * diff1_1; 
      T1xsi1 = (T1ah1 > 1.94) ?  T1ah1 - 1.0 : 0.3333333*T1ah1*T1ah1*(1.0-T1ah1*T1ah1*0.06666667); 
      T1xsi2 = (T1ah2 > 1.94) ?  T1ah2 - 1.0 : 0.3333333*T1ah2*T1ah2*(1.0-T1ah2*T1ah2*0.06666667);  
      T1xsisum = diff1 * (T1xsi1 + T1xsi2) * 0.518;
    }
    else /* diff == 0.0 */ { 
      T1ah1 = fabs(uc1*recipends) * size1; 
      T1ah2 = fabs(uc2*recipends) * size2; 
      T1xsisum = (T1ah1 + T1ah2) * 0.518;
    }
      for(i=1;i<=VPOINTS2D;i++) {
	unorm = 1.0/(v1[i] * v1[i] + v2[i] * v2[i] + 1.0e-32);
	vT1sum = dT1[i] /* - Q1 */ + v1[i] * T1x1[i] + v2[i] * T1x2[i];

	for(j=1;j<=ENODES2D;j++) {
	  Vsum = (v1[i] * GNx.vpt[GNVXINDEX(0,j,i)] +
		  v2[i] * GNx.vpt[GNVXINDEX(1,j,i)] ) ;
	  pgsfn = E->N.vpt[GNVINDEX(j,i)] + T1xsisum * Vsum * unorm;

	 
 	  /* Add to residual */
	  Eres1[j] -= dOmega.vpt[i] * (1.0 *   pgsfn *  vT1sum - E->N.vpt[GNVINDEX(j,i)] * QtE /*Qi[i]*/ +
					diff1 *  (GNx.vpt[GNVXINDEX(0,j,i)] * T1x1[i] + /* extra term for axi */
						  GNx.vpt[GNVXINDEX(1,j,i)] * T1x2[i] )); 

	}
      }
    }
 
  /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

  if(FLAGS1!=NULL) {
    onedfns=0;
    for(a=1;a<=ends;a++)
      if (FLAGS1[E->ien[el].node[a]] & FLUX_MASK1) {
	if (!onedfns++) 
	  get_global_1d_shape_fn(E,el,&GM,&dGamma,E->mesh.levmax);
	  
	nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
	nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
	nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
	nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];
	  
	for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
	  if (a==nodes[j])
	    aid = j;
	if(aid==0)  
	  printf("%d: mixed up in pg-flux int: looking for %d\n",el,a);
	  
	if (loc[a].plus[1] != 0)
	  back_front = 0;
	else back_front = dims;
	  
	for(j=1;j<=onedvpoints[dims];j++)
	  for(k=1;k<=onedvpoints[dims];k++)
	    Eres1[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] *
	      E->M.vpt[GMVINDEX(aid,j)] * g_1d[j].weight[dims-1] *
	      BCC[E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k,j)];
      }
  } 

  /* fprintf(stderr,"El %d: %g, elt res time %g\n",el,diff1,CPU_time()-time); */

  return; 
}

/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term. Calculated
   using Petrov-Galerking Shape functions.
   =========================================  */

void diff_element_residual(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision *TT1,
  standard_precision *TT1dot,
  struct SOURCES QQ1,
  higher_precision Eres1[9],
  standard_precision diff1,
  standard_precision *BCC,
  unsigned int *FLAGS1,
  int OFFSIDE_MASK1,
  int FLUX_MASK1
)

{
  int i,j,a,k,node,nodes[5],d,aid,back_front,onedfns;

  higher_precision Q1,Qi[9];
  higher_precision dT1[9];
  higher_precision T1x1[9],T1x2[9],T1x3[9];
  higher_precision T1,DT1,Ti;

  higher_precision Vsum;
  higher_precision T1xsisum,vT1sum;
  higher_precision unorm;
  higher_precision sfn,pgsfn;
    
  struct Shape_function1 GM;
  struct Shape_function1_dA dGamma;

  
  standard_precision grain_growth_source_term();
  void get_global_1d_shape_fn();
    
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const higher_precision recipends=1.0/(higher_precision)ends;
  const higher_precision diff1_1 = 1.0/(diff1+1.0e-32);
 

  for(i=1;i<=vpts;i++)	{ 
    dT1[i]=0.0;
    T1x1[i]=  0.0;
    T1x2[i]=  0.0;
    T1x3[i]=  0.0;
  }
  
  
  Q1=0.0;
  for(i=0;i<QQ1.number;i++)
    Q1 += QQ1.Q[i] * exp(-QQ1.lambda[i] * (E->monitor.elapsed_time+QQ1.t_offset));
  for(i=1;i<=vpts;i++)
    Qi[i] = Q1;
  


  for(j=1;j<=ends;j++)       {
    Eres1[j] = 0.0;

    node = E->ien[el].node[j];
    T1 = TT1[node];  
      
    if(FLAGS1[node] & (OFFSIDE_MASK1))
      DT1=0.0;
    else
      DT1 = TT1dot[node];
  
    if(3==dims) { 
	
      Ti=0.0;
      for(i=1;i<=vpts;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
	dT1[i] += DT1 * sfn; 
     	T1x1[i] += GNx.vpt[GNVXINDEX(0,j,i)] * T1;  
	T1x2[i] += GNx.vpt[GNVXINDEX(1,j,i)] * T1; 
	T1x3[i] += GNx.vpt[GNVXINDEX(2,j,i)] * T1; 
	Ti += sfn * T1; 

      if(E->advection_source_term != NULL) 
	Qi[i]=(E->advection_source_term)(E,el,i,TT1,0.0,0.0,0.0,Ti);
      }

    }
  
    else /* 2==dims */ {
          
      Ti=0.0;
      for(i=1;i<=vpts;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
	dT1[i] += DT1 * sfn;
	Ti += sfn * T1;
	T1x1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T1; 
	T1x2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T1;   
      
	if(E->advection_source_term != NULL) {
	  Qi[i]=(E->advection_source_term)(E,el,i,TT1,0.0,0.0,0.0,Ti);
	}
      }
    }
  }

  /* Construct PG shape functions */ 

  
  if(3==dims) {
  
      for(i=1;i<=VPOINTS3D;i++) {
	for(j=1;j<=ENODES3D;j++) {
	  Eres1[j] -= dOmega.vpt[i] * ( E->N.vpt[GNVINDEX(j,i)] *  dT1[i] - Qi[i] +
					diff1 * (GNx.vpt[GNVXINDEX(0,j,i)]*T1x1[i] +
						 GNx.vpt[GNVXINDEX(1,j,i)]*T1x2[i] +
						 GNx.vpt[GNVXINDEX(2,j,i)]*T1x3[i] ) ); 
	}
	
      }
    }

   
  
  else  /* 2==dims */ {
    for(i=1;i<=VPOINTS2D;i++) {
      for(j=1;j<=ENODES2D;j++) {
	/* Add to residual */
	Eres1[j] -= dOmega.vpt[i] * ( E->N.vpt[GNVINDEX(j,i)] * dT1[i] -  Qi[i] +
					diff1 *  (GNx.vpt[GNVXINDEX(0,j,i)] * T1x1[i] + 
						  GNx.vpt[GNVXINDEX(1,j,i)] * T1x2[i] 
						  /* + extra term for axi */)); 

	}
      }
    }

 
 
  /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

  if(0 && FLAGS1!=NULL) {
    onedfns=0;
    for(a=1;a<=ends;a++)
      if (FLAGS1[E->ien[el].node[a]] & FLUX_MASK1) {
	if (!onedfns++) 
	  get_global_1d_shape_fn(E,el,&GM,&dGamma,E->mesh.levmax);
	  
	nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
	nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
	nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
	nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];
	  
	for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
	  if (a==nodes[j])
	    aid = j;
	if(aid==0)  
	  printf("%d: mixed up in pg-flux int: looking for %d\n",el,a);
	  
	if (loc[a].plus[1] != 0)
	  back_front = 0;
	else back_front = dims;
	  
	for(j=1;j<=onedvpoints[dims];j++)
	  for(k=1;k<=onedvpoints[dims];k++)
	    Eres1[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] *
	      E->M.vpt[GMVINDEX(aid,j)] * g_1d[j].weight[dims-1] *
	      BCC[E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k,j)];
      }
  } 
    
  /* fprintf(stderr,"El %d: %g, elt res time %g\n",el,diff1,CPU_time()-time); */

  return; 
}

/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term. Calculated
   using Petrov-Galerking Shape functions.
   =========================================  */

void element_Tdot(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision *TT1,
  standard_precision *dTT1,
  struct SOURCES QQ1,
  higher_precision Eres1[9],
  standard_precision diff1,
  standard_precision *BCC,
  unsigned int *FLAGS1,
  int OFFSIDE_MASK1,
  int FLUX_MASK1
)

{
  int i,j,a,k,node,nodes[5],d,aid,back_front,onedfns;
 
  higher_precision size1,size2,size3;
  
  higher_precision Q1,Qi[9];
  higher_precision dT1[9];
  higher_precision T1x1[9],T1x2[9],T1x3[9];
  higher_precision T1,DT1,Ti;
  higher_precision sfn;
    
  struct Shape_function1 GM;
  struct Shape_function1_dA dGamma;

  void get_global_1d_shape_fn();
    
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const higher_precision recipends=1.0/(higher_precision)ends;
  const higher_precision diff1_1 = 1.0/(diff1+1.0e-32);
 

  for(i=1;i<=vpts;i++)	{ 
    dT1[i]=0.0;
    T1x1[i]=  0.0;
    T1x2[i]=  0.0;
    T1x3[i]=  0.0;
  }
  
  size1 = (higher_precision) E->eco[el].ntl_size[1]; 
  size2 = (higher_precision) E->eco[el].ntl_size[2]; 
  size3 = (higher_precision) E->eco[el].ntl_size[3]; 
  
  Q1=0.0;
  for(i=0;i<QQ1.number;i++)
    Q1 += QQ1.Q[i] * exp(-QQ1.lambda[i] * (E->monitor.elapsed_time+QQ1.t_offset));
  for(i=1;i<=vpts;i++)
    Qi[i] = Q1;
 

  if(3==dims) { 
    for(j=1;j<=ends;j++) {
      Eres1[j] = 0.0;
      node = E->ien[el].node[j];
      T1 = TT1[node];  
      DT1 = dTT1[node];  

   
      Ti=0.0;
      for(i=1;i<=vpts;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
     	T1x1[i] += GNx.vpt[GNVXINDEX(0,j,i)] * T1;  
	T1x2[i] += GNx.vpt[GNVXINDEX(1,j,i)] * T1; 
	T1x3[i] += GNx.vpt[GNVXINDEX(2,j,i)] * T1; 
	Ti += sfn * T1; 
	dT1[i] += sfn * DT1;

      if(E->advection_source_term != NULL) 
	Qi[i]=(E->advection_source_term)(E,el,i,TT1,0.0,0.0,0.0,Ti);
      }
    }


    for(i=1;i<=VPOINTS3D;i++) {
      for(j=1;j<=ENODES3D;j++) {
	Eres1[j] -= dOmega.vpt[i] * (dT1[i] /* - Qi[i]  */ +
				       diff1 * (GNx.vpt[GNVXINDEX(0,j,i)]*T1x1[i] +
						GNx.vpt[GNVXINDEX(1,j,i)]*T1x2[i] +
						GNx.vpt[GNVXINDEX(2,j,i)]*T1x3[i] ) ); 
      }
    }
  }
  
  else /* 2==dims */ {
    for(j=1;j<=ends;j++) {
      Eres1[j] = 0.0;
      node = E->ien[el].node[j];
      T1 = TT1[node];  
      DT1 = dTT1[node];  

      Ti=0.0;
      for(i=1;i<=VPOINTS2D;i++)  {
	sfn = E->N.vpt[GNVINDEX(j,i)];
	Ti += sfn * T1;
	dT1[i] += sfn * DT1;
	T1x1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T1; 
	T1x2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T1;   

	if(E->advection_source_term != NULL) {
	  Qi[i]=(E->advection_source_term)(E,el,i,TT1,0.0,0.0,0.0,Ti);
	}
      }
    }
    
    for(i=1;i<=VPOINTS2D;i++) {
      for(j=1;j<=ENODES2D;j++) {
	Eres1[j] -= dOmega.vpt[i] * (dT1[i] /* - Qi[i] */ +
				     diff1 * (GNx.vpt[GNVXINDEX(0,j,i)]*T1x1[i] +
					      GNx.vpt[GNVXINDEX(1,j,i)]*T1x2[i] ) ); 
      }
    }
  }

  /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

  if(0 && FLAGS1!=NULL) {
    onedfns=0;
    for(a=1;a<=ends;a++)
      if (FLAGS1[E->ien[el].node[a]] & FLUX_MASK1) {
	if (!onedfns++) 
	  get_global_1d_shape_fn(E,el,&GM,&dGamma,E->mesh.levmax);
	  
	nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
	nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
	nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
	nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];
	  
	for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
	  if (a==nodes[j])
	    aid = j;
	if(aid==0)  
	  printf("%d: mixed up in pg-flux int: looking for %d\n",el,a);
	  
	if (loc[a].plus[1] != 0)
	  back_front = 0;
	else back_front = dims;
	  
	for(j=1;j<=onedvpoints[dims];j++)
	  for(k=1;k<=onedvpoints[dims];k++)
	    Eres1[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] *
	      E->M.vpt[GMVINDEX(aid,j)] * g_1d[j].weight[dims-1] *
	      BCC[E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k,j)];
      }
  } 
  /* fprintf(stderr,"El %d: %g, elt res time %g\n",el,diff1,CPU_time()-time); */
  return; 
}

/* =====================================================
   Obtain largest possible timestep (no melt considered)

   For cylindrical/spherical it is necessary to consider
   the natural coordinate system for the element in order
   to determine a suitable timestep.

   =====================================================  */

void std_timestep(
		  struct All_variables *E
		  )
{ 
  static int been_here = 0;
  static standard_precision diff_timestep,root3,root2;
  int i,j,d,m,n,el;
  int diffratio;
  
  standard_precision adv_timestep,advp_timestep,volumetric_timestep;
  standard_precision ts,uc1,uc,size,step;
  standard_precision diffc;
  standard_precision vx,vz,vy,vc; /*RAA: 23/05/02, added vy here*/
  
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const standard_precision recipends=1.0/ends;
  
  diff_timestep = 1.0e32; 
  for(i=1;i<=E->mesh.nel;i++)  {
    if(E->tracer.tr_in_element_number[E->mesh.levmax][i]==0)
      continue;
    
    diffc=1.0e-32;
    for(j=0;j<E->tracer.tr_in_element_number[E->mesh.levmax][i];j++) {
      m = E->tracer.tr_in_element[E->mesh.levmax][j+E->tracer.tr_in_element_offset[E->mesh.levmax][i]];
      if(E->tracer.property_group[m] < 0)
	continue;
      
      diffc=max(diffc, E->tracer.Therm_diff[E->tracer.property_group[m]]);
    }
    
    for(d=1;d<=dims;d++)    {
      ts = 2.0 * E->eco[i].size[d] * E->eco[i].size[d] / diffc;
      if (diff_timestep > ts) diff_timestep = ts;
    }
  }
  
  adv_timestep = 1.0e32;
  for(i=1;i<=E->mesh.nel;i++) {
    uc=0.0;	  
    
    if(3==dims) {
      for(n=1;n<=ENODES3D;n++) {
	uc1=E->V[1][E->ien[i].node[n]]*E->V[1][E->ien[i].node[n]]+
	  E->V[2][E->ien[i].node[n]]*E->V[2][E->ien[i].node[n]]+
	  E->V[3][E->ien[i].node[n]]*E->V[3][E->ien[i].node[n]];
	uc += uc1;
      }
      
      if(uc!= 0.0) {
	uc *= recipends;
	step = E->eco[i].size[1]*E->eco[i].size[2]*E->eco[i].size[3];
	step /= sqrt(uc * (E->eco[i].size[1]*E->eco[i].size[1]*E->eco[i].size[2]*E->eco[i].size[2] +
			   E->eco[i].size[2]*E->eco[i].size[2]*E->eco[i].size[3]*E->eco[i].size[3] +
			   E->eco[i].size[1]*E->eco[i].size[1]*E->eco[i].size[3]*E->eco[i].size[3]));
	
	adv_timestep = min(adv_timestep,step);
      }
    }
    else {
      for(n=1;n<=ENODES2D;n++) {
	uc1=E->V[1][E->ien[i].node[n]]*E->V[1][E->ien[i].node[n]]+
	  E->V[2][E->ien[i].node[n]]*E->V[2][E->ien[i].node[n]];
	uc += uc1;
      }
      
      if(uc != 0.0) {
	uc *= recipends;
	step = E->eco[i].size[1]*E->eco[i].size[2];
	step /= sqrt(uc * (E->eco[i].size[1]*E->eco[i].size[1]+E->eco[i].size[2]*E->eco[i].size[2]));
	adv_timestep = min(adv_timestep,step);
      }
    }
  }

  advp_timestep = 1.0e32;
  
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    el = E->tracer.tracer_elt[E->mesh.levmax][m];
    vx=vz=0.0;
    if(3==dims)  /*RAA: 23/05/02, added these 2 lines*/
      vy=0.0;
    
    for(i=1;i<=enodes[E->mesh.nsd];i++) {
      vx += E->V[1][E->ien[el].node[i]] * E->tracer.sfn_values[E->mesh.levmax][m].node[i];
      vz += E->V[2][E->ien[el].node[i]] * E->tracer.sfn_values[E->mesh.levmax][m].node[i];
      if(3==dims)  /*RAA: 23/05/02, added these 2 lines*/
        vy += E->V[3][E->ien[el].node[i]] * E->tracer.sfn_values[E->mesh.levmax][m].node[i];
    }


    if(2==dims) {  /*RAA: 23/05/02, added this distinction*/
      step = E->eco[el].size[1]*E->eco[el].size[2];
      vc = (vx * vx + vz * vz);
      if(vc != 0.0)
        step /= sqrt(vc * (E->eco[el].size[1]*E->eco[el].size[1]+E->eco[el].size[2]*E->eco[el].size[2]));
    }
    else if(3==dims) {  /*RAA: 23/05/02, added this part*/
      step = E->eco[el].size[1]*E->eco[el].size[2]*E->eco[el].size[3];
      vc = (vx * vx + vz * vz + vy * vy);
      if(vc != 0.0)
        step /= sqrt(vc * (E->eco[el].size[1]*E->eco[el].size[1]+E->eco[el].size[2]*E->eco[el].size[2]+E->eco[el].size[3]*E->eco[el].size[3]));
    }
    
    /* if(step < advp_timestep)
       fprintf(stderr,"Particle %d (%g,%g) dominates the timestep %g\n",m,vx,vz,step); */

    advp_timestep = min(advp_timestep,step);
    
  }

  if(E->advection.fixed_timestep != 0.0) {
    advp_timestep = min(advp_timestep,E->advection.fixed_timestep);
    diff_timestep = min(diff_timestep,E->advection.fixed_timestep);
  }


  /* And the final thing is to check if the
     timestep allows stable elastic behaviour ... an upper limit
     to the timestep is given by elastic_timestep / 3.0  */

/*  if(E->control.ELASTICITY) {
    if(3.0 * advp_timestep >= E->advection.elastic_timestep) {   
      fprintf(E->fp,"Timestep adjusted to elastic: %g <- %g\n",
	      advp_timestep, E->advection.elastic_timestep * 0.333333333333);
      
    advp_timestep = E->advection.elastic_timestep * 0.333333333333;
    }
  }

/*


  /* Adjust timestep according to the fine-tuning
     parameter. C.O'N: removed if 1 stuff*/
  
  advp_timestep *= E->advection.fine_tune_dt;
  diff_timestep *= E->advection.fine_tune_dt;

#if 0
  if(E->advection.fixed_timestep != 0.0) {
    advp_timestep = 
    diff_timestep = E->advection.fixed_timestep;
  }
#endif


  
  /* How many diffusion steps could be crammed into one advection step ? */
  
  diffratio = (int)((advp_timestep/diff_timestep) + 1.0e-10);
  
  if(diffratio < 1)
    diffratio = 1;


  /* If this ratio is bigger than 1, then we adjust the 
     diffusive timestep downwards so that an exact number
     fit into the advective timestep. Then we can do `ratio'
     diffusive steps followed by one advection step 

     
     On the other hand, if it's bigger than 100, then this starts to 
     look a little dangerous. So under that condition, we instead
     drop the advective step to match ... the result is the same
     from the point of view of the diffusion algorithm 
  */
    

  if(diffratio > 1 && diffratio <= 100)
    diff_timestep = advp_timestep / diffratio;
  else if (diffratio > 100) {
    diffratio=100;
    advp_timestep = 100.0 * diff_timestep; 
  }
 
  /* Copy to global variables */
  
  E->advection.timestep_diff = diff_timestep;
  E->advection.timestep_adv = advp_timestep;
  E->advection.diff_ratio = diffratio;
  
  if(E->control.ELASTICITY && E->advection.fixed_timestep == 0.0 && E->advection.timesteps == 1)
    E->advection.timestep = 1.e-16 ;
  else
    E->advection.timestep = advp_timestep;

  /* Stable timestep is the smaller of the values, but
     if we make `ratio' substeps, it should end up 
     being the same value as the advection step */
  
  E->advection.timestep_diff = min(E->advection.timestep,E->advection.timestep_diff);
  
  
  fprintf(E->fp,"Timestep %g --- %g (diff) v %g (adv) ratio %d\n",
	  E->advection.timestep,diff_timestep,advp_timestep,diffratio); 
  printf("Working out timestep \n"); 

  return; 
}
