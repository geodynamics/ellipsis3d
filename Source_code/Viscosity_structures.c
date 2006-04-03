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



#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "function_prototypes.h"
#include <stdlib.h>


void read_viscosity_variables(
     struct All_variables *E
)
{
  void twiddle_thumbs();
  void propogate_info_from_viscosity();
  void get_viscosity_all_elements();
  void get_viscosity_one_level();
 
  int i,j,el,level,rheo,material;
  char default_str[40];
  
  const int vpts = vpoints[E->mesh.nsd];

  /* default values .... */

  for(i=0;i<40;i++) {
    E->viscosity.N0[i]=1.0;
    E->viscosity.T[i] = 0.0;
    E->viscosity.Z[i] = 0.0;
    E->viscosity.E[i] = 0.0;
    E->viscosity.T0[i] = 0.0;
    E->viscosity.layer_depth[i] = 1.0;
  }

  /* read in information */
  
  /* Global parameters */
  input_boolean("VMAX",&(E->viscosity.MAX),"off");
  input_boolean("VMIN",&(E->viscosity.MIN),"off");
  input_boolean("TDEPV",&(E->viscosity.TDEPV),"on");
  input_boolean("SDEPV",&(E->viscosity.SDEPV),"off");
  input_boolean("GRDEPV",&(E->viscosity.GRAINSIZE),"off");
  input_boolean("YIELD",&(E->viscosity.YIELD),"off");
 
  input_std_precision("visc_max",&(E->viscosity.max_value),"nodefault");
  input_std_precision("visc_min",&(E->viscosity.min_value),"nodefault");
      
  E->viscosity.update_viscosity = get_viscosity_all_elements;  /* Function definition */
   
  /* VISCOSITY REFERED TO MATERIAL POINTS ASSOCIATED WITH TRACER
     PARTICLES

     Note: phase change effects are not included as different phases are
     expected to be found as different materials in the rheological definitions.
  */

  for(material=0;material<=E->tracer.NUM_MATERIALS;material++) {

    sprintf(default_str,"Material_%d_rheol_cpts",material);
    input_int(default_str,&(E->tracer.visc[material].rheologies),"1");
    rheo = E->tracer.visc[material].rheologies;

    sprintf(default_str,"Material_%d_Trange_min",material);
    input_std_precision_vector(default_str,rheo,E->tracer.visc[material].Trange_min,-1.0e32);
    sprintf(default_str,"Material_%d_Trange_max",material);
    input_std_precision_vector(default_str,rheo,E->tracer.visc[material].Trange_max, 1.0e32);
      
    sprintf(default_str,"Material_%d_rheol_phase",material);
    input_int_vector(default_str,rheo,E->tracer.visc[material].phase,0);

    sprintf(default_str,"Material_%d_rheol_T_type",material);
    input_int_vector(default_str,rheo,E->tracer.visc[material].RHEOL_T_type,2);

    /*RAA: 24/09/02, C. O'Neill - melting stuff */
    sprintf(default_str,"Material_%d_depl_T_type",material);
    input_int_vector(default_str,rheo,E->tracer.visc[material].DEPL_T_type,1);
    sprintf(default_str,"Material_%d_melt_model",material);
    input_int_vector(default_str,rheo,E->tracer.visc[material].MELT_model,2);

    sprintf(default_str,"Material_%d_viscT1",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].T),1.0);
    sprintf(default_str,"Material_%d_viscZ",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].Z),0.0);
    sprintf(default_str,"Material_%d_viscE",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].E),0.0);
    sprintf(default_str,"Material_%d_viscT0",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].T0),0.0);
    sprintf(default_str,"Material_%d_viscN0",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].N0),1.0);
      
    sprintf(default_str,"Material_%d_viscTmin",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].Tmin_value),0.0);
    sprintf(default_str,"Material_%d_viscTmax",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].Tmax_value),1.0e32);

    sprintf(default_str,"Material_%d_sdepv_expt",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].sdepv_expt),1.0);
    sprintf(default_str,"Material_%d_grsize_expt",material);
    input_std_precision_vector(default_str,rheo,(E->tracer.visc[material].gr_size_expt),0.0);

    sprintf(default_str,"Material_%d_Bulk_visc",material);
    input_std_precision(default_str,&(E->tracer.visc[material].Bulk_visc_ratio),"-1.0");
    sprintf(default_str,"Material_%d_Elas_shear_modulus",material);
    input_std_precision(default_str,&(E->tracer.visc[material].Elas_shear_mod),"0.0");
    sprintf(default_str,"Material_%d_Num_elas_bulk_modulus",material);
    input_std_precision(default_str,&(E->tracer.visc[material].Num_Elas_Bulk_mod),"0.0");
    sprintf(default_str,"Material_%d_Elas_bulk_modulus",material);
    input_std_precision(default_str,&(E->tracer.visc[material].Elas_Bulk_mod),"0.0");
 
    sprintf(default_str,"Material_%d_ortho_visc_ratio",material);
    input_std_precision(default_str,&(E->tracer.visc[material].Ortho_viscosity_ratio),"1.0");

    sprintf(default_str,"Material_%d_mobile_visc_ratio",material);
    input_std_precision(default_str,&(E->tracer.visc[material].mobile_phase_ratio),"1.0");

    sprintf(default_str,"Material_%d_shear_heating_coeff",material);
    input_std_precision(default_str,&(E->tracer.Hv0[material]),"0.0");



    sprintf(default_str,"Material_%d_coss_A11",material);
    input_std_precision(default_str,&(E->tracer.coss[material].A11),"0.0");
    sprintf(default_str,"Material_%d_coss_A22",material);
    input_std_precision(default_str,&(E->tracer.coss[material].A22),"0.0");
    sprintf(default_str,"Material_%d_coss_A12",material);
    input_std_precision(default_str,&(E->tracer.coss[material].A12),"0.0");
    sprintf(default_str,"Material_%d_coss_G11",material);
    input_std_precision(default_str,&(E->tracer.coss[material].G11),"0.0");
    sprintf(default_str,"Material_%d_coss_G12",material);
    input_std_precision(default_str,&(E->tracer.coss[material].G12),"0.0");
    sprintf(default_str,"Material_%d_coss_G22",material);
    input_std_precision(default_str,&(E->tracer.coss[material].G22),"0.0");
    sprintf(default_str,"Material_%d_coss_B31",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B31),"0.0");
    sprintf(default_str,"Material_%d_coss_B32",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B32),"0.0");
    sprintf(default_str,"Material_%d_coss_B12",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B12),"0.0");
    sprintf(default_str,"Material_%d_coss_B21",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B21),"0.0");
    sprintf(default_str,"Material_%d_coss_B13",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B13),"0.0");
    sprintf(default_str,"Material_%d_coss_B23",material);
    input_std_precision(default_str,&(E->tracer.coss[material].B23),"0.0");

    sprintf(default_str,"Material_%d_yield_stress_B0",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_B0),"1.0e32");
    sprintf(default_str,"Material_%d_yield_stress_Bz",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Bz),"0.0");
    sprintf(default_str,"Material_%d_yield_stress_Bp",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Bp),"0.0");
    sprintf(default_str,"Material_%d_yield_stress_Bc",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Bc),"1.0e32");
      
    sprintf(default_str,"Material_%d_yield_stress_ET",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_ET),"1.0");
    sprintf(default_str,"Material_%d_yield_stress_minimum",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_minimum),"1.0e-32");
    sprintf(default_str,"Material_%d_yield_stress_maximum",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_maximum),"1.0e32");
   
    sprintf(default_str,"Material_%d_yield_stress_E0",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_E0),"1.0e32");
    sprintf(default_str,"Material_%d_yield_stress_Ea",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Ea),"1.0,0.0,1.0");
    sprintf(default_str,"Material_%d_yield_stress_En",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_En),"0.0");
    sprintf(default_str,"Material_%d_yield_stress_E0dt",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_E0dt),"0.0");
   
    sprintf(default_str,"Material_%d_yield_stress_Edot0",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Edot0),"0.0");
    sprintf(default_str,"Material_%d_yield_stress_Edota",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Edota),"1.0,0.0,1.0");
    sprintf(default_str,"Material_%d_yield_stress_Edotn",material);
    input_std_precision(default_str,&(E->tracer.visc[material].yield_stress_Edotn),"0.0");

    /*RAA: 24/09/02, C. O'Neill melting stuff  (std -> hgr precision on 10/07/03) */
    sprintf(default_str,"Material_%d_melting",material); 
    input_int(default_str,&(E->tracer.melting[material]),"0"); 
    sprintf(default_str,"Material_%d_sp0",material);
    input_std_precision(default_str,&(E->tracer.sp0[material]),"5.0");
    sprintf(default_str,"Material_%d_sp1",material);
    input_hgr_precision(default_str,&(E->tracer.sp1[material]),"0.0");
    sprintf(default_str,"Material_%d_sp2",material);
    input_hgr_precision(default_str,&(E->tracer.sp2[material]),"0.0");
    sprintf(default_str,"Material_%d_sp3",material);
    input_hgr_precision(default_str,&(E->tracer.sp3[material]),"0.0");

    /*RAA: 24/09/02, C. O'Neill melting stuff  (std -> hgr precision on 10/07/03) */
    sprintf(default_str,"Material_%d_lp0",material);
    input_std_precision(default_str,&(E->tracer.lp0[material]),"10.0");
    sprintf(default_str,"Material_%d_lp1",material);
    input_hgr_precision(default_str,&(E->tracer.lp1[material]),"0.0");
    sprintf(default_str,"Material_%d_lp2",material);
    input_hgr_precision(default_str,&(E->tracer.lp2[material]),"0.0");
    sprintf(default_str,"Material_%d_lp3",material);
    input_hgr_precision(default_str,&(E->tracer.lp3[material]),"0.0");

    /*if(E->tracer.visc[material].Bulk_visc_ratio<0.0)  
      E->tracer.visc[material].Pen_bulk=100.0 ;
    else
      E->tracer.visc[material].Pen_bulk=min(100,E->tracer.visc[material].Bulk_visc_ratio-2.0) ;    */
     E->tracer.visc[material].Pen_bulk=0.0;
  }
  /* If appropriate, initialize grain size information
     which then folds back into the rheological law. */
  if(E->viscosity.GRAINSIZE != 0) grain_growth_initialize(E);

 
#if 0
  if(E->control.ELASTICITY) {
    for(i=0;i<E->tracer.NUM_MATERIALS;i++) {
      if(E->tracer.visc[i].Elas_shear_mod > 0.0) {
	if(E->advection.timestep > VE_TIMETOLERANCE / E->tracer.visc[i].Elas_shear_mod)
	  E->advection.timestep = VE_TIMETOLERANCE / E->tracer.visc[i].Elas_shear_mod;
	 
	fprintf(E->fp,"Elasticity of material %d sets initial timestep to %g\n",i,E->advection.timestep);
      }
    }
  }
#endif
  return;
}

void propogate_info_from_viscosity(
				   struct All_variables *E,
				   int maxlevel
				   )
{
  void build_diagonal_of_K();
  void construct_node_maps_3_level();
  void construct_node_ks_3_level();
  void construct_elt_gs_level(); 
  void add_stress_bcs_to_global_F();
  
  int i,e,lv,jj;
  
  for(i=maxlevel;i>=E->mesh.levmin;i--) {
    /* should be out of the main loop as it is purely geometric*/
    construct_elt_gs_level(E,i);
    construct_node_maps_3_level(E,i);


    construct_node_ks_3_level(E,i);
    add_stress_bcs_to_global_F(E,i) ;
    build_diagonal_of_K(E,i);
    build_diagonal_of_Ahat_level(E,i);
  }
  return;
}

void get_viscosity_all_elements
(
 struct All_variables *E,
 int level
 )
{
  int i,el,j,m,k,l,node;
  int yielded,material;
  int trmin, trmax; /*RAA: 4/02/03, added these variables*/
 
  const int dims = E->mesh.nsd;
  const int dofs = E->mesh.dof;
  const int vpts = vpoints[E->mesh.nsd];

  standard_precision temp;
  standard_precision etal,eta0,etaT,etaS,etaG;
  standard_precision scale,old_visc;
  standard_precision minval,maxval,maxneg;
  standard_precision yield_point;
  /*RAA: 24/09/02, C. O'Neill melting stuff */
  standard_precision Tprime,fdep,xdep,depcase,meltcase,tau1_sign;
  standard_precision *solidus,*liquidus;

  standard_precision elastime1;
  
  standard_precision eta1,eta2,eta3;
  standard_precision lNx[4][ELNMAX+1];
  standard_precision dudx[4][4];
  standard_precision u1,u2,u3;

  standard_precision *pressure,*Dkk,pres,dil;
  standard_precision tau[7][7],D[7][7];

  standard_precision *Node;
  standard_precision *TrD;

  standard_precision vsselfdot6();

  static int been_here = 0;

  Node = (standard_precision *) Malloc0((E->mesh.nno+1) * sizeof(standard_precision));
  TrD = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision));
  pressure = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision));
  Dkk = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision));
  /*RAA: 24/09/02, C. O'Neill melting stuff */
  solidus = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision)); 
  liquidus = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision)); 

  if(E->control.verbose) 
    fprintf(stderr,"Compute all viscosities for level %d\n",level);

  printf("Going into tracer loop!!! \n");

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    /* Supply a dummy viscosity value (0.5) which allows us to get back
       the deformation and stress history in the form we need. If not elastic
       then edot and TrD should come out the same this way. This trick doesn't
       work for the pressure, so we have to fix that later */

    /*E->tracer.Visc[m] = 0.5;

    tracer_deviatoric_stress_and_pressure(E,tau,D,&pres,&dil,NULL,m,level); 
     C.O'N: taking out this dummy visc and tracer_dev call, using instead get_D which I 
     hijakced from FD's version*/

    get_D(E,tau,D,&dil,m,level) ;

    tau1_sign = (tau[1][1] < 0) ? -1 : 1;
  
    if(2==dims && 2==dofs) {
      E->tracer.edot[m] = sqrt(0.5 * (D[1][1] * D[1][1] + D[2][2] * D[2][2] + 
				      2.0 * D[1][2] * D[1][2])); 
      TrD[m] = sqrt(0.5 * (tau[1][1] * tau[1][1] + tau[2][2] * tau[2][2] + 
			   2.0 * tau[1][2] * tau[1][2])); 
      /*RAA: check data*/
      /*
      material=E->tracer.property_group[m]; 
      if (material == 2) 
        fprintf(E->fp1,"tr: %d; edot = %g, D11 D22 D12 = %g  %g  %g\n",m,E->tracer.edot[m],D[1][1],D[2][2],D[1][2]);
	*/
      /*RAA: end of data checking for strain rate*/
    }
    if(3==dims && 3==dofs) {
      E->tracer.edot[m] = sqrt(0.5 * (D[1][1] * D[1][1] + D[2][2] * D[2][2] + D[3][3] * D[3][3] +
				      2.0 * D[1][2] * D[1][2] +
				      2.0 * D[1][3] * D[1][3] +
				      2.0 * D[2][3] * D[2][3] )); 
      TrD[m] = sqrt(0.5 * (tau[1][1] * tau[1][1] + tau[2][2] * tau[2][2] + tau[3][3] * tau[3][3] +
			   2.0 * tau[1][2] * tau[1][2] +
			   2.0 * tau[1][3] * tau[1][3] +
			   2.0 * tau[2][3] * tau[2][3] )); 
    }    
    Dkk[m] = dil;
  }

  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    
    old_visc = E->tracer.Visc[m];
    E->tracer.Visc[m] = 1.0e-32;

    /* Figure out the material */

    material=E->tracer.property_group[m]; 
    if(material > MAX_MATERIALS || material > E->tracer.NUM_MATERIALS) {
      fprintf(stderr,"Tracer %d, color %d is illegal\n",m,E->tracer.property_group[m]);
      material = 0;
    }

    for(l=0;l<E->tracer.visc[material].rheologies;l++) {
      etal = E->tracer.visc[material].N0[l];
     
      if(E->viscosity.TDEPV){ 
	temp = E->tracer.T[m];
      
	if(temp < E->tracer.visc[material].Tmin_value[l]) 
	  temp = E->tracer.visc[material].Tmin_value[l];
	if(temp > E->tracer.visc[material].Tmax_value[l])
	  temp = E->tracer.visc[material].Tmax_value[l];
	
	switch ( E->tracer.visc[material].RHEOL_T_type[l] )   {
	case 1:	/* Frank-Kamenetskii approximation to Arrhenius rheology */
	  eta0 = exp(-E->tracer.visc[material].T[l] * temp);
	  break;    
	case 2:       /* Arrhenius Rheology */
	  eta0 = exp((E->tracer.visc[material].E[l] + E->tracer.visc[material].Z[l] * E->tracer.tz[m])/
		     (E->tracer.visc[material].T[l] * (temp + E->tracer.visc[material].T0[l])));
	  break;
	default:
	  fprintf(stderr,"There is no expression for Rheology(T) type %d\n",l);
	  exit(-1);
	}
	etal *= eta0;
      }

      /*RAA: 24/09/02, C. O'Neill - melting stuff */
      if(E->tracer.melting[material]) {
        switch (E->tracer.visc[material].DEPL_T_type[l] )   { 
          case 1:  
            depcase=1;  
            break; 
          case 2: 
            depcase=2;  
            break;  
          default: 
            depcase=1;  
          /* nfi */  
	      }
        switch (E->tracer.visc[material].MELT_model[l] )   {
          case 1:
	          meltcase=1;
	          break;
          case 2:
	          meltcase=2;
	          break;
          default:
	          meltcase=2;
	          /* nfi */
        }
      } 

      //printf("past first melt stuff \n");

      if(E->viscosity.SDEPV) {
	etaT =  etal;
	if(E->tracer.edot[m] != 0.0) {
	  scale=pow(E->tracer.edot[m],(1.0-E->tracer.visc[material].sdepv_expt[l])/
		    E->tracer.visc[material].sdepv_expt[l] );
	
	  etaS=pow(etaT,1.0/ E->tracer.visc[material].sdepv_expt[l]) * scale;
	  etal = etaS;  
	}
      }
 
      if(E->viscosity.GRAINSIZE) {
	if(E->tracer.grain_size[m] != 0.0) {
	  eta0 = etal;
	  etal = pow(E->tracer.grain_size[m], E->tracer.visc[material].gr_size_expt[l]) * eta0;
	}
      }

      /* For this we assume that dT/dz is ~1 - may need to tone this down !!*/

      if(E->tracer.T[m] < E->tracer.visc[material].Trange_min[l]) {
	scale = pow(31.0,(E->tracer.visc[material].Trange_min[l]-E->tracer.T[m])); 
	etal = max(1.0e32,etal * scale);
      }

      if(E->tracer.T[m]> E->tracer.visc[material].Trange_max[l]) {
	scale = pow(31.0,(E->tracer.T[m]-E->tracer.visc[material].Trange_max[l])); 
	etal = max(1.0e32,etal * scale);
      }
     
      if(E->tracer.Phases[E->tracer.property_group[m]] == 1) {
    	if(etal!=0.0)
	  E->tracer.Visc[m] += 1.0 / etal;
      }
      else
	if(etal!=0.0 && (E->tracer.Current_phase[m] == E->tracer.visc[material].phase[l]))
	  E->tracer.Visc[m] += 1.0 / etal;
      
      if (m == 436) {
	    fprintf(stderr,"IN VISC HERE, PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[m]].Pen_bulk);
	    fprintf(stderr,"IN VISC HERE, WHATS WITH VISC (1) %g \n",E->tracer.Visc[m]);
    }

    } /* Next component of the composite rheology at this tracer */

    E->tracer.Visc[m] = 1.0 / E->tracer.Visc[m] ;


    if (m == 436) {
	    fprintf(stderr,"IN VISC HERE, WHATS WITH VISC (2) %g \n",E->tracer.Visc[m]);
    }

#if 0
    if(E->tracer.tz[m]<=0.6 && E->tracer.tz[m]>=0.4 && E->tracer.tx[m]<=0.1)
      E->tracer.Visc[m] += 1. * (E->tracer.tz[m] - 0.4) * E->tracer.Visc[m] ;
#endif

    /* Compressibility, this parameter is not used for incompressible case, except for elasticity !!!!!!!!!  */
    /* This is actually not the bulk modulus but lambda = bulk -2/dims* visc */
    E->tracer.BulkVisc[m] = 
      (E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio-E->tracer.visc[E->tracer.property_group[m]].Pen_bulk)
      * E->tracer.Visc[m];

#if  defined (CHEM_TRANS_MODULE_INSTALLED) 
    if(E->control.CHEM_TRANS && E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio != 1.0)
      E->tracer.Visc[m] *= 
	( E->tracer.volfraction[m] * E->tracer.visc[E->tracer.property_group[m]].mobile_phase_ratio + 
	  1.0 - E->tracer.volfraction[m]);
#endif

    /* ELASTICITY: If the material is visco-elastic, then we modify the viscosity to
       allow elastic deformation to be accumulated */ 

    E->tracer.Visc0[m] = E->tracer.Visc[m];


    if(E->control.ELASTICITY) {
      /* I - volumetric term - uses shear viscosity to scale bulk,
	 since the interrelationship is not clear at the moment */ 
      if(E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio > 0.0 && 
	 E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod > 0.0) {
	elastime1 =  E->advection.elastic_timestep * E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod;
	E->tracer.BulkVisc[m] = E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio * E->tracer.Visc[m] *
	  elastime1 / (E->tracer.visc[E->tracer.property_group[m]].Bulk_visc_ratio * E->tracer.Visc[m] + elastime1);
      }
      /* Need to be fixed in the incompressible case */

      /* II - Shear term (overwrites previous value) */
      if(E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod > 0.0) {
	elastime1 =  E->advection.elastic_timestep * E->tracer.visc[E->tracer.property_group[m]].Elas_shear_mod;
	E->tracer.Visc[m] = E->tracer.Visc[m] * elastime1 / (E->tracer.Visc[m] + elastime1);
      }
    }
    //printf("past elas stuff \n");
	
    if (m == 436) {
	    fprintf(stderr,"IN VISC HERE, WHATS WITH VISC (3) %g \n",E->tracer.Visc[m]);
    }


    /* Now we can compute the pressure properly */
    
    /*    pressure[m] = E->tracer.DQ1[m]  - (2.0/E->mesh.nsd)*E->tracer.Visc[m] * Dkk[m] ;*/
    if(2==E->mesh.nsd) {
      pressure[m] = E->tracer.DQ1[m]  - (1.0+E->tracer.visc[E->tracer.property_group[m]].Pen_bulk)*E->tracer.Visc[m] * Dkk[m] ;
    }
    else {
      pressure[m] = E->tracer.DQ1[m] - (0.66666666666667+E->tracer.visc[E->tracer.property_group[m]].Pen_bulk) * E->tracer.Visc[m] * Dkk[m];
    }

    if(E->control.ELASTICITY && E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod > 0.0) {    
      pressure[m] += E->tracer.BulkVisc[m] * E->tracer.Pt[m] / 
	(E->tracer.visc[E->tracer.property_group[m]].Elas_Bulk_mod * E->advection.elastic_timestep);
    }
    //printf("past pressure etc\n");

    /* AND shouldn't there be a part from the Kelvin elasticity model ?? */


    /* Depletion stuff - C. O'Neills work on melting*/
#if 0   /*RAA: my start to the melting aspect - but will use Craig's for uniformity  w/ 2D*/
    if(E->tracer.melting[material]) {
        temp=E->tracer.T[m];
        if (depcase==1) {
        solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*pressure[m] + E->tracer.sp2[material]*pressure[m]*pressure[m] + E->tracer.sp3[material]*pressure[m]*pressure[m]*pressure[m];
        liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*E->tracer.tz[m] + E->tracer.lp2[material]*pressure[m]*pressure[m] + E->tracer.lp3[material]*pressure[m]*pressure[m]*pressure[m];
        }
        else {
          solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*E->tracer.tz[m]+ E->tracer.sp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.sp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
          liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*E->tracer.tz[m]+ E->tracer.lp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.lp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
        } 
    
        /*RAA: 17/10/02, define McK & B's  Tprime and X */
        Tprime = (temp - ((solidus[m] + liquidus[m])/2.0))/(liquidus[m] - solidus[m]);
    
        if (Tprime >= 0.5) 
           Tprime = 0.5;
        else if (Tprime <= -0.5) 
           Tprime = -0.5;
    
        xdep = Tprime + (Tprime*Tprime - 0.25)*(0.4256 + 2.988*Tprime) + 0.5;
    
        if (E->monitor.elapsed_time <= 0.0) { /*Pressure initialization to 0 should be noted*/
           xdep = 0.0;
           E->tracer.depl[m] = 0.0;
        }
        /*RAA, don't really need these 4 lines given the above constraints on Tprime*/
        if (xdep > 1.0) 
           xdep = 1.0;
        else if (xdep < 0.0) 
           xdep = 0.0;
    if(E->tracer.melting[material]) {
        temp=E->tracer.T[m];
        if (depcase==1) {
        solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*pressure[m] + E->tracer.sp2[material]*pressure[m]*pressure[m] + E->tracer.sp3[material]*pressure[m]*pressure[m]*pressure[m];
        liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*E->tracer.tz[m] + E->tracer.lp2[material]*pressure[m]*pressure[m] + E->tracer.lp3[material]*pressure[m]*pressure[m]*pressure[m];
        }
        else {
          solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*E->tracer.tz[m]+ E->tracer.sp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.sp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
          liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*E->tracer.tz[m]+ E->tracer.lp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.lp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
        } 

        /*RAA: keep it simple, melt fraction probably should be able to decrease if necessary, 
           but will at some point need to address latent heat issues. */
        E->tracer.depl[m]=xdep;
      /*  if (xdep > 0.0) 
          fprintf(E->fp1,"tracer depth pres temp solidus liquidus  Tprime xdep: %d  %g  %g  %g  %g   %g   %g   %g\n",m,E->tracer.tz[m],pressure[m],E->tracer.T[m],solidus[m],liquidus[m],Tprime,xdep);  */
    
     
    }  /*RAA - end of melting stuff */
#endif
	/* RAA: 11/07/03, start of C. O'Neill's depletion stuff */
       	/*   but I added a melting flag for each material for easy bypassing*/
  if(E->tracer.melting[material]) {

   /* E->tracer.dFdot[m] = 0.0; RAA - initialize, just in case */
   /* E->tracer.depl[m]=0.0;    RAA - initialize, just in case -CRAP!!! What the hell is this in here for????
     	This has just reset all the tracers to zero depletion again!!!*/
    if(E->tracer.melting[material]) {

      temp=E->tracer.T[m];
    
      if (pressure[m]<=0) {
        fdep=0;
      }
      else {
        if (depcase==2){
          solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*E->tracer.tz[m]+ E->tracer.sp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.sp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
          liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*E->tracer.tz[m]+ E->tracer.lp2[material]*E->tracer.tz[m]*E->tracer.tz[m] + E->tracer.lp3[material]*E->tracer.tz[m]*E->tracer.tz[m]*E->tracer.tz[m];
        } 
        else {
          solidus[m]=E->tracer.sp0[material] + E->tracer.sp1[material]*pressure[m] + E->tracer.sp2[material]*pressure[m]*pressure[m] + E->tracer.sp3[material]*pressure[m]*pressure[m]*pressure[m];
          liquidus[m]=E->tracer.lp0[material] + E->tracer.lp1[material]*pressure[m] + E->tracer.lp2[material]*pressure[m]*pressure[m] + E->tracer.lp3[material]*pressure[m]*pressure[m]*pressure[m];
        }
  
        if (temp > solidus[m]) {
          fdep=(temp - solidus[m])/(liquidus[m]-solidus[m]);
          if (fdep > 1) 
	          fdep=1;
        }
        else 
          fdep=0;
      }
      /* McKenzie and Bickle melt curve */
      xdep=fdep;
      fdep=2.0684*xdep - 4.0564*xdep*xdep + 2.988*xdep*xdep*xdep;
      xdep=fdep;
      /*if (E->tracer.tx[m]>0.19 && E->tracer.tx[m]<0.21) */
       /* fprintf(E->fp1,"Depl xdep %g %g \n",E->tracer.depl[m],xdep); */
  
       /* Want to get dF/dt */
       xdep -= E->tracer.depl[m];
       if (xdep > 0) {
         E->tracer.dFdot[m] = xdep/E->advection.timestep;
       }
       else {
         E->tracer.dFdot[m] = 0.0;
       }
  
       if (E->tracer.tx[m]>0.19 && E->tracer.tx[m]<0.21) {
        /* fprintf(E->fp1,"dFdot xdep %g %g \n",E->tracer.dFdot[m],xdep);*/
  /*	printf("Melt case: %g \n",meltcase); */
       }
       if (meltcase==1) {
         E->tracer.depl[m]=fdep;
	/* printf("Why am i here? %g \n",meltcase); */
	 }
       else {
         if (fdep > E->tracer.depl[m]) {
	         E->tracer.depl[m]=fdep;
	/*	 printf("Much better %g %g %g \n",meltcase,fdep,E->tracer.depl[m]); */
	 }
       }
    }  /*RAA - matches 'if melting material' condition */
  } /* RAA: 11/07/03, End of C. O'Neill's depletion stuff */
      
    /* VISCO-PLASTIC PART to simulate brittle failure 
       in the simplest possible manner ... the effective
       viscosity is modified to limit the second invariant of
       the stress tensor (including stored stresses for viscoelastic materials)  */
    
    if(been_here++ && E->viscosity.YIELD && TrD[m] != 0.0) {
      yield_point = E->tracer.visc[material].yield_stress_B0 + 
	E->tracer.visc[material].yield_stress_Bz * E->tracer.tz[m];
      
      /* Pressure dependence */ 
      
      yield_point +=  E->tracer.visc[material].yield_stress_Bp  * pressure[m]; 

      /* Tension cutoff (shouldn't this be a function of strain history ??)*/

      if(pressure[m] < 0.0 && fabs(pressure[m]) > E->tracer.visc[material].yield_stress_Bc)
	yield_point *= 0.001;  /* or something ??!! */
	

      /* Strain dependence - if required */
 
      if(E->tracer.visc[material].yield_stress_E0 != 0.0 && E->tracer.edotp_integrated[m] > 0.0 
	 /* && (E->tracer.DQ1[m] +  E->tracer.DQ[m]) < 0.0 */ ) 
	yield_point -= (1.0-E->tracer.visc[material].yield_stress_Ea) * yield_point *
	  min(1.0,pow(E->tracer.edotp_integrated[m]/E->tracer.visc[material].yield_stress_E0,
		      E->tracer.visc[material].yield_stress_En)); 
      
      /* Strain rate dependence - if required (same form as strain weakening */
 
      if(E->tracer.visc[material].yield_stress_Edot0 != 0.0 && E->tracer.edot[m] != 0.0) 
	yield_point -= (1.0-E->tracer.visc[material].yield_stress_Edota) * yield_point *
	  min(1.0,pow(E->tracer.edot[m]/E->tracer.visc[material].yield_stress_Edot0,
		      E->tracer.visc[material].yield_stress_Edotn)); 
     
      /* Apply truncation at low strength to keep numerically tractable, and 
	 at high strength to mimic semi-brittle effects if required */
      
      if (yield_point < E->tracer.visc[material].yield_stress_minimum)
	yield_point = E->tracer.visc[material].yield_stress_minimum;
   
      else if (yield_point > E->tracer.visc[material].yield_stress_maximum)
	yield_point = E->tracer.visc[material].yield_stress_maximum;
     
      /* Viscosity which limits the stress invariant to
	 be at the yield point */

      etal = yield_point / (TrD[m]);


      /* If this viscosity is smaller than the purely viscous
	 one defined by the composites so far, then we say that this
	 tracer has yielded */
      
      if(etal < E->tracer.Visc[m]) {
	E->tracer.yielded[m] =  1.0 - etal/E->tracer.Visc[m]/* E->tracer.Visc[m]/etal - 1.0 */ ;

	if(etal < E->tracer.Visc[m] * 0.001) /* ???  can we catch this with a global truncation ??? */
	  etal = E->tracer.Visc[m] * 0.001;

	/* Smooth application of transition */
	E->tracer.Visc[m] = etal;

	/* And, if plastic deformation, the relevant viscosity 
	   for computing the shear heating rate is the current one */
	
	E->tracer.Visc0[m] = E->tracer.Visc[m];



      }
      else {
	E->tracer.yielded[m] = 0.0;
      }
    }

    E->tracer.strd[m] = E->tracer.Visc[m] * TrD[m];
    E->tracer.strd1[m] = E->tracer.strd[m] * tau1_sign;
    E->tracer.strs[m] = E->tracer.strd[m] + pressure[m]; /*RAA: 6/02, fixed bug, used to be dims*pressure */

    /* Max/Min values (global) */

    if(E->viscosity.MIN && E->tracer.Visc[m] < E->viscosity.min_value * pow(3.0,(double)(E->mesh.levmax-level)) ) {
      E->tracer.Visc[m] = E->viscosity.min_value /* pow(3.0,(double)(E->mesh.levmax-level))*/; 
    }
    if(E->viscosity.MAX && E->tracer.Visc[m] > E->viscosity.max_value * pow(10.0,-(double)(E->mesh.levmax-level))) {
      E->tracer.Visc[m] = E->viscosity.max_value /* pow(10.0,-(double)(E->mesh.levmax-level))*/;
    }
  }
 
  if (m == 436) {
	    fprintf(stderr,"IN VISC HERE, WHATS WITH VISC (4) %g \n",E->tracer.Visc[m]);
    }

  if(E->control.verbose) 
    fprintf(stderr,"Computed all viscosities for level %d\n",level);

  if(1 /*&&  level == E->mesh.levmax*/) {
    minval =  1.0e32; 
    maxval = -1.0e32;
    for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
      if(E->tracer.Visc[m] > maxval)  {
  	    maxval = E->tracer.Visc[m]; 
	      trmax = m;
      }
      if(E->tracer.Visc[m] < minval) { 
	      minval = E->tracer.Visc[m]; 
	      trmin = m;
      }
    }
    if(E->control.print_convergence)
      fprintf(stderr,"%d tracers, %d: Visc min = %g, Visc max = %g, eg. --> Tracers %d, %d, resp.\n",
	      E->tracer.NUM_TRACERS,level,minval,maxval,trmin,trmax);
    
    fprintf(E->fp,"Visc min = %g, Visc max = %g, eg.--> Tracers %d, %d, resp.\n",minval,maxval,trmin,trmax);
  }

  /*RAA: 11/4/01, check the tracer temp and pressure*/
    /*for(m=1;m<=E->tracer.NUM_TRACERS;m++) {*/
       /*fprintf(E->fp1,"*** tracer temp and press: %d %g %g\n",m,E->tracer.T[m],pressure[m]); */
       /*fprintf(E->fp1,"*** tracer visc and eq. dev. stress(?): %d %g %g\n",m,E->tracer.Visc[m],TrD[m]); */
       /*fprintf(E->fp1,"RAA: Tracer: %d,: strs11,22,33: %g   %g   %g\n",m,E->tracer.S11[m],E->tracer.S22[m],E->tracer.S33[m]);*/
       /*}*/
	  
  /*    for(i=1; i<=20; i++) { 
           fprintf(E->fp1,"RAA: Tracer: %d,: strs: %g, Visc: %g, TrD: %g, Edot: %g,  lev: %d \n",i,E->tracer.strs[i],E->tracer.Visc[i],TrD[i],E->tracer.edot[m],level);
        }
  */
       
  propogate_info_from_viscosity(E,level);
    
  if(E->control.verbose) 
    fprintf(stderr,"Applied all viscosities for level %d\n",level);

  free((void *) Node);
  free((void *) TrD);
  free((void *) pressure);
  free((void *) Dkk);
  free((void *) liquidus);  /*RAA: C. O'Neill - melting stuff */
  free((void *) solidus);   /*RAA: C. O'Neill - melting stuff*/
  printf("Out of visco \n");
  return;
}
