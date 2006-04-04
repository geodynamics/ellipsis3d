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



/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite of files.
*/

#include "config.h"

#include <math.h>

#if HAVE_MALLOC_H
#include <malloc.h>
#endif

#if HAVE_STDLIB_H
#include <stdlib.h> /* for "system" command */
#endif

#if HAVE_STRING_H
#include <string.h>
#endif

#if HAVE_STRINGS_H
#include <strings.h>
#endif

#if HAVE_RPC_RPC_H
#include <rpc/rpc.h>
#endif

#if HAVE_RPC_XDR_H
#include <rpc/xdr.h> 
#endif

#include "element_definitions.h"
#include "global_defs.h"

/* Output files are created, if requested, using the global filename plus a numerical tag.

   Data is stored according to type:
   1. Nodal information (for a given timestep)
   2. Horizontal averages (for a given timestep)
   3. Surface observables (and other slices) (for a given timestep).
   4. Large-scale averages recorded as time sequence.

   The format of the output is a header including individual column labels
   which are required because the information stored can vary. The
   comment character in these files is a `#' terminated by the end of the
   line. 

   NOTE FOR IMPROVEMENT:  THE NODE OUTPUT STRING CAN BE USED FOR THE COLUMN LABELLING 
   INSTEAD OF DOING THE CHECKING TWICE.
*/

void generic_data_storage(
			  struct All_variables *E,
			  int file_number
			  )
{
  int i,j,k,ii,number,node,data_length,sample,l,m,el,element;
  int record_everything,checkpt_and_log,only_timelog;
  char output_file[255],command[255];
  FILE *fp;

  int i_tmp,block,tr;
  float f_tmp;
  double d_tmp;
  float  *data_array;

  static int been_here=0;
  static int time_files=0;
   
  float test;
  
  standard_precision storage_timing,checkpt_timing;
  standard_precision time,time1,CPU_time();
  standard_precision errorva,errorvr,errorra,errorrr,R,Gc,G,sum1,sum2,sum,sumcte,sumerrorav,sumerrorar,sumerrorrv,sumerrorrr;
  standard_precision CC,DD,Co,Si,MM,ave,delta,AA,BB,ta,mu,sigma,omeg,eta1,eta2,eta3,fact,Q,QQ,VV[5];
  standard_precision lN[ELNMAX+1],visc,visc_bulk,theo_GR,axial_stress,stress;
  standard_precision vel[1000],velfft[1000],mag ;
  standard_precision Omega0,Omega1,Omega2,xsie,xsip,kk,a,b,c,k0 ;
  standard_precision maxP,minP,meanP,S11,S12,S22,deflection,upper,lower ;

  const int dofs = E->mesh.dof;
  const int dims = E->mesh.nsd;
  const int nno = E->mesh.nno;
  const int outputcase = E->control.outputcase;

  /* Where did FD's Malloc0's go? - different test cases in FD's version */

  
  /* The data files are written permanently once every "storage_timesteps" steps.
     This is computed in the mean sense .... once the time since the last
     write exceeds "storage_timesteps" * the mean timestep then the data are
     saved for posterity. Otherwise they are just checkpointed 
     to the .00000. files */


 /* fprintf(stderr,"HERE IS TRACER.PT (0)   %g \n",E->tracer.Pt[436]) ;	*/
  printf("In output \n");
  switch(outputcase) {
  case 1: /* Sigma, rot */
    sumerrorav = sumerrorar = sumerrorrv = sumerrorrr = sum1 = sum2 = 0.0 ;
    Gc = 0.5*(E->tracer.coss[0].G11 - E->tracer.coss[0].G12) ;
    G = 0.5 *(E->tracer.coss[0].G11 + E->tracer.coss[0].G12) ;
    MM = E->tracer.coss[0].B31 ;
    R = sqrt(0.5*MM/Gc) ;
    delta = sqrt(R*R*0.5*(1+Gc/G)) ;
    Co = cosh(0.5/delta) ;
    Si = sinh(0.5/delta) ;
    ta = Si / Co ;
    omeg = 0.0 ;
    sigma = 1.0 ;
    AA = sigma ;
    BB = AA/(2*G*Co) ;
    for(i=1;i<=nno;i++) {
      if(fabs(0.5-E->sx[2][i]) >= 0.0000001) {
	errorva = E->V[1][i] - (AA*(0.5-E->sx[2][i])/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-E->sx[2][i])/delta)) ;
	errorvr = errorva / (E->V[1][i] - errorva) ;
      }
      else {
	errorvr = 0.0 ;
	errorva = 0.0 ;
      }
      if(fabs(E->sx[2][i]) >= 0.0000001 && (fabs(1-E->sx[2][i]) >= 0.0000001) ) {
	errorra =  (-E->V[3][i] + (E->V[3][1] + omeg)) - (-AA/(2.*G) + BB*cosh((0.5-E->sx[2][i])/delta));
	errorrr = errorra / (E->V[3][i] - errorra) ;
      }
      else {
	errorrr = 0.0 ;
	errorra = 0.0 ;
      }
      sum1 += fabs(E->V[1][i]) ;
      sum2 += fabs(E->V[3][i] - (E->V[3][1] - omeg)) ;
      sumerrorav += fabs(errorva) ;
      sumerrorar += fabs(errorra) ;
      sumerrorrv += fabs(errorvr) ;
      sumerrorrr += fabs(errorrr) ;
    }
    fprintf(E->fp1,"%g %g %g \n",E->monitor.elapsed_time,E->V[3][25],-(-AA/(2.*G) + BB)) ;
    fprintf(E->fp2,"%g %g %g \n",E->monitor.elapsed_time,E->V[1][13],(AA*(0.5-0.25)/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-0.25)/delta)));
    fprintf(stderr,"Error on V = %g and on R = %g \n",(sumerrorav/sum1)*100,(sumerrorar/sum2)*100) ;
    fprintf(stderr,"Num = %g Ana = %g \n",E->V[1][13],-(-AA/(2.*G) + BB)) ;
    break ;
  case 2: /* vit rot */
    sumerrorav = sumerrorar = sumerrorrv = sumerrorrr = sum1 = sum2 = 0.0 ;
    Gc = 0.5*(E->tracer.coss[0].G11 - E->tracer.coss[0].G12) ;
    G = 0.5 *(E->tracer.coss[0].G11 + E->tracer.coss[0].G12) ;
    MM = E->tracer.coss[0].B31 ;
    R = sqrt(0.5*MM/Gc) ;
    delta = sqrt(R*R*0.5*(1+Gc/G)) ;
    Co = cosh(0.5/delta) ;
    Si = sinh(0.5/delta) ;
    ta = Si / Co ;
    omeg = 0.0 ;
    if(E->control.ELASTICITY) {
      AA = (E->V[1][1]*E->tracer.coss[0].G11+2.*Gc*delta*E->V[3][1]*ta) / (0.5*(1+G/Gc)-(Gc/G)*delta*ta) ;
      BB = (E->V[3][1] + AA*0.5/G) / Co ;
    }
    else {
      BB = (E->V[1][1] - E->V[3][1]) / (Co -2*Gc*delta*Si/(G+Gc)) ;
      AA = 2*G*(BB*Co + E->V[3][1]) ;
    }
    for(i=1;i<=nno;i++) {
      if(fabs(0.5-E->sx[2][i]) >= 0.0000001) {
	errorva = E->V[1][i] - (AA*(0.5-E->sx[2][i])/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-E->sx[2][i])/delta)) ;
	errorvr = errorva / (E->V[1][i] - errorva) ;
      }
      else {
	errorvr = 0.0 ;
	errorva = 0.0 ;
      }
      if(fabs(E->sx[2][i]) >= 0.0000001 && (fabs(1-E->sx[2][i]) >= 0.0000001) ) {
	errorra =  (-E->V[3][i] + (E->V[3][1] + omeg)) - (-AA/(2.*G) + BB*cosh((0.5-E->sx[2][i])/delta));
	errorrr = errorra / (E->V[3][i] - errorra) ;
      }
      else {
	errorrr = 0.0 ;
	errorra = 0.0 ;
      }
      sum1 += fabs(E->V[1][i]) ;
      sum2 += fabs(E->V[3][i] - (E->V[3][1] - omeg)) ;
      sumerrorav += fabs(errorva) ;
      sumerrorar += fabs(errorra) ;
      sumerrorrv += fabs(errorvr) ;
      sumerrorrr += fabs(errorrr) ;
    }
    fprintf(E->fp1,"%g %g %g \n",E->monitor.elapsed_time,E->V[3][25],-(-AA/(2.*G) + BB)) ;
    fprintf(E->fp2,"%g %g %g \n",E->monitor.elapsed_time,E->V[1][13],(AA*(0.5-0.25)/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-0.25)/delta)));
    fprintf(stderr,"Error on V = %g and on R = %g \n",(sumerrorav/sum1)*100,(sumerrorar/sum2)*100) ;
    fprintf(stderr,"Num = %g Ana = %g \n",E->V[1][13],-(-AA/(2.*G) + BB)) ;
    break ;
  case 3: /* couple stress, rot */
    sumerrorav = sumerrorar = sumerrorrv = sumerrorrr = sum1 = sum2 = 0.0 ;
    Gc = 0.5*(E->tracer.coss[0].G11 - E->tracer.coss[0].G12) ;
    G = 0.5 *(E->tracer.coss[0].G11 + E->tracer.coss[0].G12) ;
    MM = E->tracer.coss[0].B31 ;
    R = sqrt(0.5*MM/Gc) ;
    delta = sqrt(R*R*0.5*(1+Gc/G)) ;
    Co = cosh(0.5/delta) ;
    Si = sinh(0.5/delta) ;
    ta = Si / Co ;
    omeg = 0.0 ;
    mu = 1.e2 ;
    BB = mu*delta/(MM*Si) ;
    AA = 2*G*(E->V[1][1]-BB*delta*Si/(G+Gc));
    omeg = -AA/(2*G) + BB*Co ;
    for(i=1;i<=nno;i++) {
      if(fabs(0.5-E->sx[2][i]) >= 0.0000001) {
	errorva = E->V[1][i] - (AA*(0.5-E->sx[2][i])/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-E->sx[2][i])/delta)) ;
	errorvr = errorva / (E->V[1][i] - errorva) ;
      }
      else {
	errorvr = 0.0 ;
	errorva = 0.0 ;
      }
      if(fabs(E->sx[2][i]) >= 0.0000001 && (fabs(1-E->sx[2][i]) >= 0.0000001) ) {
	errorra =  (-E->V[3][i] + (E->V[3][1] + omeg)) - (-AA/(2.*G) + BB*cosh((0.5-E->sx[2][i])/delta));
	errorrr = errorra / (E->V[3][i] - errorra) ;
      }
      else {
	errorrr = 0.0 ;
	errorra = 0.0 ;
      }
      sum1 += fabs(E->V[1][i]) ;
      sum2 += fabs(E->V[3][i] - (E->V[3][1] - omeg)) ;
      sumerrorav += fabs(errorva) ;
      sumerrorar += fabs(errorra) ;
      sumerrorrv += fabs(errorvr) ;
      sumerrorrr += fabs(errorrr) ;
    }
    fprintf(E->fp1,"%g %g %g \n",E->monitor.elapsed_time,E->V[3][25],-(-AA/(2.*G) + BB)) ;
    fprintf(E->fp2,"%g %g %g \n",E->monitor.elapsed_time,E->V[1][13],(AA*(0.5-0.25)/G - 2*Gc*delta*BB/(G+Gc)*sinh((0.5-0.25)/delta)));
    fprintf(stderr,"Error on V = %g and on R = %g \n",(sumerrorav/sum1)*100,(sumerrorar/sum2)*100) ;
    fprintf(stderr,"Num = %g Ana = %g \n",E->V[1][13],-(-AA/(2.*G) + BB)) ;
    break ;
  case 6: /* Egg yolks flow, yummy yummy !!! */
    sum1 = 0.0 ;
    for(i=1;i<=nno;i++) {
      if(fabs(E->sx[1][i] - 1.5) <= 0.000001 && E->sx[2][i] <= 0.7 && E->sx[2][i] >= 0.3) {
	if(fabs(E->sx[2][i] - 3.125) <= 0.000001 || fabs(E->sx[2][i] - 6.875) <= 0.000001)
	  sum1 += E->V[1][i] * (9.*0.625/8.) ;
	else
	  sum1 += E->V[1][i] * 0.625 ;
      }
    }
    sum = 0 ;
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if(E->tracer.property_group[i] == 1 && E->tracer.tx[i] <= 1.5) {
	sum += 1 ;
      }
    }
    break;
  case 5: /* Visco-elastic spreadable square 1D */ /*C.O'N: removed the #if 1 bit, amking them different cases instead */
   
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
/*      fprintf(E->fp3,"%d %g %g ",i,E->tracer.tx[i],E->tracer.tz[i]) ;*/
      if(fabs(E->tracer.tx[i]-0.5) <= 0.1 && fabs(E->tracer.tz[i]-(E->x[2][nno] + E->x[2][1])/2.) <= 0.1) {
	tr = i ;
	break ;
      }
    }
    fprintf(stderr,"Tracer no:  %d  \n",tr) ;
    BB = E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio ;
    fprintf(stderr,"BB:  %g  \n",BB) ;
    if(E->control.ELASTICITY) {
      AA = E->tracer.visc[0].Elas_Bulk_mod ;
      fprintf(stderr,"AA:  %g  \n",AA);
      time = E->monitor.elapsed_time+E->advection.timestep ; 
      fprintf(stderr,"time:  %g  monitor: %g advection %g \n",time,E->monitor.elapsed_time,E->advection.timestep); 
      if(E->V[2][1] == 0) {
	Q =  E->control.secret_case_fl[3]*exp(-AA*(time-E->control.secret_case_fl[2])/BB);
	 fprintf(stderr,"in if, Q:  %g  \n",Q);
	}
      else {
	E->control.secret_case_fl[1] += AA*exp(AA*time/BB)*E->advection.timestep/(1./E->V[2][1]-time) ;
	Q = E->control.secret_case_fl[1] * exp(-AA*time/BB) ;
	E->control.secret_case_fl[3] = Q ;
	E->control.secret_case_fl[2] = time ; 
	fprintf(stderr,"in else, Q:  %g  \n",Q);
      }
      fprintf(E->fp1,"%g  %g  %g \n",time,Q,E->tracer.Pt[tr]) ;
      /*fprintf(stderr,"Viscosity %g  \n",E->tracer.Visc[436]);*/
      fprintf(stderr,"Time %g Ana  %g Num  %g \n",time,Q,E->tracer.Pt[tr]) ;	
      fprintf(stderr,"In output BCmoveZ0 %g Control cases, 1: %g 2:  %g 3:  %g V21 %g Vb %g \n",E->mesh.BCvelocityZ0,E->control.secret_case_fl[1],E->control.secret_case_fl[2],E->control.secret_case_fl[3],E->V[2][1],E->Vb[2][1][1]) ;
      /*fprintf(stderr,"What is PEN_BULK: %g \n",E->tracer.visc[E->tracer.property_group[436]].Pen_bulk);      */
    }
    else {
      time = E->monitor.elapsed_time ;
      Q = (E->tracer.visc[0].Bulk_visc_ratio-1)*E->tracer.Visc[tr]*E->V[2][1]/(E->x[2][nno] - E->x[2][1]) ;
      fprintf(E->fp1,"%g %g ",time,Q) ;
    }
    break ;


    /*Old RAA version */
    /*sum = 0.0 ;
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if(fabs(E->tracer.tx[i]-0.5) <= 0.1 && fabs(E->tracer.tz[i]-(E->x[2][E->mesh.nno] + E->x[2][1])/2.) <= 0.1){
	tr = i ;
      }
    }    
    if(E->control.ELASTICITY) {
      time = E->monitor.elapsed_time+E->advection.timestep ;
      if(E->V[2][1] == 0) {
	Q =  E->control.secret_case_fl[3]*
	  exp(-E->tracer.visc[0].Elas_Bulk_mod*(time-E->control.secret_case_fl[2])/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio));
      }
      else {
	E->control.secret_case_fl[1] += E->tracer.visc[0].Elas_Bulk_mod*E->V[2][1]*E->advection.timestep*
	  exp(E->tracer.visc[0].Elas_Bulk_mod*time/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio))/
	  (1.-E->V[2][1]*time);
	Q = E->control.secret_case_fl[1] * exp(-E->tracer.visc[0].Elas_Bulk_mod*time/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio)) ;
	E->control.secret_case_fl[3] = Q ;
	E->control.secret_case_fl[2] = time ;
      }
      fprintf(E->fp1,"%g  %g  %g \n",time,Q,E->tracer.Pt[tr]) ;
      fprintf(stderr,"Time %g Ana  %g Num  %g \n",time,Q,E->tracer.Pt[tr]) ;	
      fprintf(stderr,"In output BCmoveZ0 %g Control cases, 1: %g 2:  %g 3:  %g V21 %g Vb %g \n",E->mesh.BCvelocityZ0,E->control.secret_case_fl[1],E->control.secret_case_fl[2],E->control.secret_case_fl[3],E->V[2][1],E->Vb[2][1][1]) ;	
    }
    else {
      time = E->monitor.elapsed_time ; */
/*      fprintf(stderr,"time = %g \n",time) ;*/
/*	Q = (E->tracer.visc[0].Bulk_visc_ratio-1)*E->tracer.Visc[tr]*E->V[2][1]/(E->x[2][E->mesh.nno] - E->x[2][1]) ;
    fprintf(stderr,"At x = %g pres n = %g a = %g \n",
	    E->tracer.sample_x[0],E->tracer.Pt[tr],Q) ;
    fprintf(E->fp1,"%g %g %g \n",time,Q,E->tracer.Q[tr]) ;
    }
    break ; */
  case 7: /* Visco-elastic spreadable square 2D */
	/* #else -just removing this -C.O'N */
    fprintf(stderr,"What is tr, firstly? %d \n",tr);
    tr = 0;
    time = E->monitor.elapsed_time + E->advection.timestep ;
    time1 = E->monitor.elapsed_time ;
    fprintf(stderr,"Made it into case 7, time %g Adv_time %g Elapsed_time %g \n",time,E->advection.timestep,E->monitor.elapsed_time);
   /* fprintf(stderr,"HERE IS THE PROBLEM: Elas_Bulk_Modul & Bulk_visc_ratio %g %g \n",E->tracer.visc[0].Elas_Bulk_mod,E->tracer.visc[0].Bulk_visc_ratio);
    fprintf(stderr,"What about Elas shear mod or N0 %g %g \n",E->tracer.visc[0].Elas_shear_mod,E->tracer.visc[0].N0[0] );*/
    for(i=0;i<2;i++) {
      fprintf(stderr,"Trying general_tracers_element for coords %g %g %g \n",E->tracer.sample_x[i],E->tracer.sample_z[i],E->tracer.sample_y[i]);
      /*C.O'N: Why was there no E->tracer.sample_y[i] in this next line? */
      element = general_tracers_element(E,1,E->tracer.sample_x[i],E->tracer.sample_z[i],E->tracer.sample_y[i],&eta1,&eta2,&eta3,E->mesh.levmax) ;
      fprintf(stderr,"Trying v_shape_fn for element %d and coords %g %g %g \n",element,E->tracer.sample_x[i],E->tracer.sample_z[i],E->tracer.sample_y[i]);
      v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax) ;
      VV[i]=0.0 ;
      for(j=1;j<=enodes[E->mesh.nsd];j++) {
	VV[i] += E->V[1][E->ien[element].node[j]] * lN[j];
      }
    }
    fprintf(stderr,"Done general_tracers_element and v_shape\n");
	/* C.O'N: made the thresholds larger (0.01->0.05) so we actually satisfy these conditions */
    fprintf(stderr,"number of tracers %d\n",E->tracer.NUM_TRACERS);
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if( (E->tracer.property_group[i] == 0) && (fabs(E->tracer.tx[i]-0.5) <= 0.05) ){ 
	 fprintf(stderr,"DID I GET HERE??? %d %d\n",tr,i);
	 if (fabs(E->tracer.tz[i]-(E->x[2][E->mesh.nno] + E->x[2][1])/2.) <= 0.05) {
	  fprintf(stderr,"dO i EVEN EVER MAKE IT IN HERE?? %d %d \n",tr,i);    
	  tr = i ;
	  break ;
         }
      }
    }
    fprintf(stderr,"What is tr, by the way? %d \n",tr);
    if(E->control.ELASTICITY) {
	fprintf(stderr,"In elasticity control\n");
    fprintf(stderr,"Testing Elas_Bulk_Modul & Bulk_visc_ratio %g %g \n",E->tracer.visc[0].Elas_Bulk_mod,E->tracer.visc[0].Bulk_visc_ratio);
    fprintf(stderr,"What about Elas shear mod or N0 %g %g \n",E->tracer.visc[0].Elas_shear_mod,E->tracer.visc[0].N0[0] );
    fprintf(stderr,"0: N0  %g Bulk visc ratio  %g \n",E->tracer.visc[0].N0[0],E->tracer.visc[0].Bulk_visc_ratio) ;
 fprintf(stderr,"1: N0  %g Bulk visc ratio  %g \n",E->tracer.visc[1].N0[0],E->tracer.visc[1].Bulk_visc_ratio) ;
 /* Since material 0 exists, and material 1 doesn't, I'm changing all the tracer.visc[1] to [0] -C.O'N.*/
      sum = 0.0 ;
      if(E->tracer.visc[0].Elas_Bulk_mod > 0.0)
	AA = 1./E->tracer.visc[0].Elas_Bulk_mod + 1./E->tracer.visc[0].Elas_shear_mod ;
      else
	AA = 1./E->tracer.visc[0].Elas_shear_mod ;

      if(E->tracer.visc[0].Bulk_visc_ratio > 0.0)
	BB = 1./(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio) + 1./E->tracer.visc[0].N0[0] ;
      else
	BB = 1./E->tracer.visc[0].N0[0] ;

      if(E->V[2][1] == 0) {
	Q =  E->control.secret_case_fl[3]*exp(-BB*(time-E->control.secret_case_fl[2])/AA);
 	fprintf(stderr,"Testing Q(a) (whatever the hell this is) %g\n",Q);

      }
      else {
	E->control.secret_case_fl[1] += 2./AA*exp(BB*time/AA)*E->advection.timestep/(1./E->V[2][1]-time) ;
	Q = E->control.secret_case_fl[1] * exp(-BB*time/AA) ;
	E->control.secret_case_fl[3] = Q ;
	E->control.secret_case_fl[2] = time ;
	fprintf(stderr,"Testing Q(b) (whatever the hell this is) %g\n",Q);
	fprintf(stderr,"Testing secret_case & BB, time & AA %g %g %g %g \n",E->control.secret_case_fl[1],BB,time,AA);
      }
	fprintf(stderr,"Testing Q(c) (whatever the hell this is) %g\n",Q);

	/* C.O'N: changing all the visc[1] to [0] in these statements aswell */
      if(E->tracer.visc[0].Bulk_visc_ratio > 0.0 && E->tracer.visc[0].Elas_Bulk_mod > 0.0) {
	if(E->V[2][1] == 0.0) {
	sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.)*
	  Q*(E->tracer.visc[0].Elas_Bulk_mod/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio) - 
	      E->tracer.visc[0].Elas_shear_mod/E->tracer.visc[0].N0[0])
	  /(E->tracer.visc[0].Elas_shear_mod+E->tracer.visc[0].Elas_Bulk_mod) ;
	}
	else {
	sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.)*
	  (Q*(E->tracer.visc[0].Elas_Bulk_mod/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio) - 
	      E->tracer.visc[0].Elas_shear_mod/E->tracer.visc[0].N0[0])
	   - (E->tracer.visc[0].Elas_Bulk_mod - E->tracer.visc[0].Elas_shear_mod)/(1./E->V[2][1]-time))
	  /(E->tracer.visc[0].Elas_shear_mod+E->tracer.visc[0].Elas_Bulk_mod) ;
	}
      }
      else if(E->tracer.visc[0].Bulk_visc_ratio < 0.0 && E->tracer.visc[0].Elas_Bulk_mod < 0.0) {
	if(E->V[2][1] == 0.0) {
	  sum = 0.0 ;
	}
	else {
	sum = (E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.)/(1./E->V[2][1]-time) ;
	}
      }
      else if(E->tracer.visc[0].Bulk_visc_ratio > 0.0 && E->tracer.visc[0].Elas_Bulk_mod < 0.0) {
	if(E->V[2][1] == 0.0) {
	  sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.) * Q/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio) ;
	}
	else {
	  sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.) 
	    * (Q/(E->tracer.visc[0].N0[0]*E->tracer.visc[0].Bulk_visc_ratio) - 1./(1./E->V[2][1]-time)) ;
	}
      }
      else {
	if(E->V[2][1] == 0.0) {
	sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.) * Q*E->tracer.visc[0].Elas_shear_mod/E->tracer.visc[0].N0[0]
	  /(E->tracer.visc[0].Elas_shear_mod+E->tracer.visc[0].Elas_Bulk_mod) ;
	}
	else {
	sum = -(E->tracer.sample_x[0]-E->sx[1][E->mesh.nno]/2.)*
	  (Q*( E->tracer.visc[0].Elas_shear_mod/E->tracer.visc[0].N0[0])
	   - (E->tracer.visc[0].Elas_Bulk_mod - E->tracer.visc[0].Elas_shear_mod)/(1./E->V[2][1]-time))
	  /(E->tracer.visc[0].Elas_shear_mod+E->tracer.visc[0].Elas_Bulk_mod) ;
	}
      }
	fprintf(stderr,"Out of elastic 1, here is Q %g \n",Q);

    }
    else { /* Viscous */
	    fprintf(stderr," - in viscous\n");

      if(E->tracer.visc[E->tracer.property_group[tr]].Bulk_visc_ratio > 1.0) { /* Compressible viscous */
	sum = (E->tracer.sample_x[0]-0.5)*(E->tracer.visc[1].Bulk_visc_ratio-1.)
	  /(E->tracer.visc[0].Bulk_visc_ratio+1.)/(1./E->V[2][1]-time) ;
	Q = E->tracer.visc[0].N0[0]*(1./(1./E->V[2][1]-time) + sum/(E->tracer.sample_x[0]-0.5)) ;
      }
      else { /* Incompressible viscous */
	sum = (E->tracer.sample_x[0]-0.5)/(1./E->V[2][1]-time) ;
	Q = 2.*E->tracer.visc[0].N0[0]/(1./E->V[2][1]-time) ;
      }
    }
    fprintf(stderr,"Out of elastic 2\n");
/* Changin E->fp1 and fp2 to stderr for now -C.O'N. */
    fprintf(stderr,"Why crash now? \n");
    fprintf(stderr,"Test time %g\n",time1);
    fprintf(stderr,"Testing Q (whatever the hell this is) %g\n",Q);
    fprintf(stderr,"Testing tr %d\n",tr);
    fprintf(stderr,"Testing Pt (which I played with) %g\n",E->tracer.Pt[tr]); 
    fprintf(stderr,"Testing sum %g\n",sum);
    fprintf(stderr,"Testing VV[0] %g VV[1] %g VV[2] %g\n",VV[0],VV[1],VV[2]);
 
    if(E->control.ELASTICITY) {
    	fprintf(E->fp1,"%g %g %g %g %g %g %g \n",time1,Q,E->tracer.Pt[tr],sum,VV[0],VV[1],VV[2]) ;
    	fprintf(stderr,"How'd I go, get here at all? %g %g %g \n",E->tracer.S11[tr],E->tracer.S12[tr],E->tracer.S22[tr]);
    	fprintf(E->fp1,"%g %g %g %g %g \n",time1,E->tracer.S11[tr],E->tracer.S12[tr],E->tracer.S22[tr],E->tracer.Pt[tr]);
    }
    else
      fprintf(E->fp1,"%g %g %g %g %g %g \n",time1,Q,sum,VV[0],VV[1],VV[2]) ;
    fprintf(stderr,"Finished case 7\n");

    break ;
/* #endif -C.O'N: taking out this too */
  case 8: /* buckling of a competent layer embedded in a less competent medium */
    time = E->monitor.elapsed_time + E->advection.timestep ;
    upper = -1e32 ;
    lower = 1.e32 ;
    if(E->control.ELASTICITY) {
      S11 = S12 = S22 = QQ = 0.0 ;
      tr = 0 ;
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.property_group[i] == 1) {
	  upper = max(upper,E->tracer.tz[i]) ;
	  lower = min(lower,E->tracer.tz[i]) ;
	  if (fabs(E->tracer.tx[i])<=0.01) {
	    S11 += E->tracer.S11[i] ;
	    S12 += E->tracer.S12[i] ;
	    S22 += E->tracer.S22[i] ;
	    QQ += E->tracer.Pt[i] ;
	    tr ++ ;
	  }
	}
      }
      S11 /= tr ;
      S12 /= tr ;
      S22 /= tr ;
      QQ /= tr ;
      fprintf(E->fp2,"%g %g %g %g  %g ",time,S11,S12,S22,QQ) ;
      fprintf(stderr,"Beam S11 = %g  S22 = %g ",S11,S22) ;
/*      axial_stress = -S11 ;*/
      S11 = S12 = S22 = QQ = 0.0 ;
      tr = 0 ;
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if((E->tracer.property_group[i] == 2) && (fabs(E->tracer.tx[i])<=0.01) ) {
	  S11 += E->tracer.S11[i] ;
	  S12 += E->tracer.S12[i] ;
	  S22 += E->tracer.S22[i] ;
	  QQ += E->tracer.Pt[i] ;
	  tr ++ ;
	}
      }
      S11 /= tr ;
      S12 /= tr ;
      S22 /= tr ;
      QQ /= tr ;
      fprintf(E->fp2,"%g %g %g %g \n",S11,S12,S22,QQ) ;
      fprintf(stderr,"Outer layer S11 = %g  S22 = %g \n",S11,S22) ;
    }
#if 1
    ii = 1 ;
    for(i=1;i<=E->mesh.NOX[E->mesh.levmax];i++) {
      for(j=(i-1)*E->mesh.NOZ[E->mesh.levmax]+1;j<=i*E->mesh.NOZ[E->mesh.levmax];j++) {
	tr = 0 ;
	for(k=1;k<=E->NEI[E->mesh.levmax].nels[j];k++) {
	  el = E->NEI[E->mesh.levmax].element[(j-1)*ENODES2D+k-1];	  
	  for(l=0;l<E->tracer.tr_in_element_number[E->mesh.levmax][el];l++) {
	    m = E->tracer.tr_in_element[E->mesh.levmax][l+E->tracer.tr_in_element_offset[E->mesh.levmax][el]];
	    if(E->tracer.property_group[m] == 1) {
	      tr++;
	      break;
	    }
	  }
	}
	if(tr==E->NEI[E->mesh.levmax].nels[j]) {
	  vel[ii] = E->V[2][j] ;
	  ii++ ;
	  break ;
	}
      }
    }
#else
  
    j=0 ;
    for(i=1;i<=nno;i++) {
      if(fabs(E->sx[2][i] - E->sx[2][E->mesh.nno]*0.5) <= 0.00001) {
	vel[j] = E->V[2][i] ;
	j++ ;
      }
    }
#endif
    horiz_1d_fft(E,vel,velfft,1) ;
    k = E->control.secret_case_fl[2]+1 ;
    mag = fabs(velfft[k]) ;

    axial_stress = -4*E->V[1][E->mesh.nno-1]*E->tracer.visc[1].N0[0]/E->sx[1][E->mesh.nno];
    fprintf(stderr,"axial stress = %g \n",axial_stress) ;
    kk =  E->control.secret_case_fl[2]/*M_PI/3.0*/ ;
    k0 = pow(6.*E->tracer.visc[2].N0[0]/E->tracer.visc[1].N0[0],1./3.)/E->control.secret_case_fl[3] ;
    Omega0 = axial_stress*k0*E->control.secret_case_fl[3]/(6.*E->tracer.visc[2].N0[0]) ;
/*    Omega0 = -4*pow(6.,-2./3.)*E->V[1][E->mesh.nno-1]*
      pow((E->tracer.visc[1].N0[0]/E->tracer.visc[2].N0[0]),2./3.)/E->sx[1][E->mesh.nno] ;*/
/*    Omega0 = Load*k0/(6.*E->tracer.visc[2].N0[0]) ;*/
    if(E->control.ELASTICITY) {
      kk /= k0 ;
      xsie = Omega0 * E->tracer.visc[2].N0[0] / E->tracer.visc[2].Elas_shear_mod;
      xsip = Omega0 * E->tracer.visc[1].N0[0] / E->tracer.visc[1].Elas_shear_mod;
      a = pow(kk,3.)*xsie-3.*xsie*xsip*kk+2.*xsip ;
      b = pow(kk,3.)-3.*(xsie+xsip)*kk+2. ;
      c = -3.*kk ;
      delta = b*b-4.*a*c;
      if(delta >= 0.0 ) {
	Omega1 = (-b-sqrt(delta)) / (2.*a) ;
	Omega2 = (-b+sqrt(delta)) / (2.*a) ;
      }
      else {
	Omega1 = Omega2 = 0.0 ;
      }
      theo_GR = max(Omega1,Omega2) ;
      if(xsip*xsie*xsie<1.0) {
	fprintf(stderr,"Viscous behaviour \n") ;
      }
      else {
	fprintf(stderr,"Instability expected \n") ;
      }
    fprintf(stderr,"k=%g GR = %g (N)  %g (T)  k0=%g w0=%g xsip = %g and xsie = %g \n",
	    kk,(mag/E->control.secret_case_fl[4])/Omega0,theo_GR,k0,Omega0,xsip,xsie) ;
    }
    else {
      for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
	if(E->tracer.property_group[i] == 1) {
	  upper = max(upper,E->tracer.tz[i]) ;
	  lower = min(lower,E->tracer.tz[i]) ;
	}
      }
      Omega0 = -4*pow(6.,-2./3.)*E->V[1][E->mesh.nno-1]*
	pow((E->tracer.visc[1].N0[0]/E->tracer.visc[2].N0[0]),2./3.)/E->sx[1][E->mesh.nno] ;
      theo_GR = (-(E->V[1][E->mesh.nno-1]*E->control.secret_case_fl[3]*4.*E->tracer.visc[1].N0[0]/E->sx[1][E->mesh.nno])/
		 (4./kk*E->tracer.visc[2].N0[0] + kk*kk*E->tracer.visc[1].N0[0]*
		  E->control.secret_case_fl[3]*E->control.secret_case_fl[3]*E->control.secret_case_fl[3]/3.))/Omega0 ;
      kk /= k0 ;
      fprintf(stderr,"k=%g GR = %g (N)  %g (T)  k0=%g w0=%g \n",
	      kk,(mag/E->control.secret_case_fl[4])/Omega0,theo_GR,k0,Omega0) ;
    }
    deflection = 0.5*(upper - lower - E->control.secret_case_fl[3]);
    fprintf(stderr,"mag = %g var = %g w0 = %g \n",mag,E->control.secret_case_fl[4],Omega0) ;
    fprintf(E->fp1,"%g %g  %g  %g %g %g %g  \n",
	    kk,E->control.secret_case_fl[1],time,theo_GR,(mag/E->control.secret_case_fl[4])/Omega0,
	    (mag/E->control.secret_case_fl[1])/Omega0,(mag/deflection)/Omega0) ;
    E->control.secret_case_fl[4] += mag * E->advection.timestep ;
    break ;
  case 10: /* buckling of a competent layer embedded in a less competent medium */
    S11 = S12 = S22 = 0.0 ;
    tr = 0 ;
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if((E->tracer.property_group[i] == 1) && (fabs(E->tracer.tx[i])<=0.01) ) {
	S11 += E->tracer.S11[i] ;
	S12 += E->tracer.S12[i] ;
	S22 += E->tracer.S22[i] ;
	tr ++ ;
      }
    }
    S11 /= tr ;
    S12 /= tr ;
    S22 /= tr ;
    time = E->monitor.elapsed_time + E->advection.timestep ;

    fprintf(E->fp1,"%g %g %g %g ",time,S11,S12,S22) ;
    S11 = S12 = S22 = 0.0 ;
    tr = 0 ;
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if((E->tracer.property_group[i] == 2) && (fabs(E->tracer.tx[i])<=0.01) ) {
	S11 += E->tracer.S11[i] ;
	S12 += E->tracer.S12[i] ;
	S22 += E->tracer.S22[i] ;
	tr ++ ;
      }
    }
    S11 /= tr ;
    S12 /= tr ;
    S22 /= tr ;
    fprintf(E->fp1,"%g %g %g \n",S11,S12,S22) ;
    break ;
  case 9: /* Rayleigh-Taylor instability */
    for(i=1;i<=nno;i++) {
      if((E->x[1][i] >= 0.9141) && 
	 (fabs(E->x[2][i] - E->control.secret_case_fl[2] * cos(M_PI*E->x[1][i]/0.9142) - E->control.secret_case_fl[1]) <= 0.00001))
	fprintf(stderr,"Growth rate = %g \n",E->V[2][i]/E->control.secret_case_fl[2]) ;
    }
    fprintf(stderr,"vrms = %g \n",E->monitor.Vrms) ;
    fprintf(E->fp2,"%g %g \n",E->monitor.elapsed_time+E->advection.timestep,E->monitor.Vrms) ;
    break ;
  }
  fflush(E->fp1) ;
  fflush(E->fp2) ;
  

  /* If timings have already been supplied, then use these.
     Otherwise compute them */
  /*fprintf(stderr,"HERE IS TRACER.PT (1)   %g \n",E->tracer.Pt[436]) ;	*/
  fprintf(stderr,"The beginning of the rest of my life\n");
  if(E->control.storage_timing != 0.0)
    storage_timing = E->control.storage_timing;
  else if(E->advection.total_timesteps != 0) 
    storage_timing = (E->monitor.elapsed_time/E->advection.total_timesteps) * E->control.storage_timesteps;  
  
  if(E->control.checkpt_timing != 0.0)
    checkpt_timing = E->control.checkpt_timing;
  else if(E->advection.total_timesteps != 0) 
    checkpt_timing = (E->monitor.elapsed_time/E->advection.total_timesteps) * E->control.checkpt_timesteps;  
  
  /* What to do ? */
  fprintf(stderr,"What to do indeed\n");

  record_everything = checkpt_and_log = 0;

    if(E->monitor.solution_cycles == 1) {
      record_everything = 1;
      checkpt_and_log=1;
    }
    else if ((E->monitor.elapsed_time - E->monitor.time_of_last_store) >  storage_timing) {
      record_everything = 1;
      checkpt_and_log = 1;
      E->monitor.time_of_last_store = E->monitor.time_of_last_checkpt = E->monitor.elapsed_time;
    }
    else if ((E->monitor.elapsed_time - E->monitor.time_of_last_checkpt) >  checkpt_timing)  {
      E->monitor.time_of_last_checkpt = E->monitor.elapsed_time;
      checkpt_and_log = 1;
    }

    time1=CPU_time();

    /* Periodic boundary conditions do not calculate anything on the
       wrapped around edge nodes. These must be copied before storing
       the data */
	fprintf(stderr,"Here come periodic boundaries\n");
       /* fprintf(stderr,"HERE IS TRACER.PT (2)   %g \n",E->tracer.Pt[436]) ;	*/
 
    if(E->mesh.periodic_x || E->mesh.periodic_y)  {
	flogical_mesh_to_real(E,E->NQ[E->mesh.levmax],E->mesh.levmax); 
	flogical_mesh_to_real(E,E->T,E->mesh.levmax);
	if(2==dofs) /*RAA: 28/5/01, Psi causes seg fault in 3D*/
	  flogical_mesh_to_real(E,E->Psi,E->mesh.levmax); 
	flogical_mesh_to_real(E,E->V[1],E->mesh.levmax);
	flogical_mesh_to_real(E,E->V[2],E->mesh.levmax);
	if(3==dofs) {
	  flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
	}
	else if(6==dofs) {
	  flogical_mesh_to_real(E,E->V[3],E->mesh.levmax);
	  flogical_mesh_to_real(E,E->V[4],E->mesh.levmax);
	  flogical_mesh_to_real(E,E->V[5],E->mesh.levmax);
	  flogical_mesh_to_real(E,E->V[6],E->mesh.levmax);	
	}
	
   }

    /*fprintf(stderr,"HERE IS TRACER.PT (3)   %g \n",E->tracer.Pt[436]) ;	*/
    if(checkpt_and_log) {
      /* As it is time-consuming to compute the nodal plastic strain
	 value, do this only if it needs to be printed out ... */

      if((strstr(E->control.node_data.which_data_types,"Pstn"))) {
	if(E->advection.timesteps > 0)
	gs_tracers_to_nodes(E,E->Pstrain,NULL,NULL,NULL,E->tracer.edotp_integrated,E->mesh.levmax,0); 
      }
	fprintf(stderr,"Nodal data approacheth\n");

      /* 1. Nodal data (position, velocity etc) */
      if(strlen(E->control.node_data.which_data_types)!=0) {
	E->control.node_data.numb_found=0;
	sprintf(output_file,"%s.%05d.node_data",E->control.data_file, record_everything ? file_number : 0); 	
	if ((fp=fopen(output_file,"w")) != NULL) {
	  fprintf(fp,"# PROCESS_ID=%06d LENGTH_SCALE=%.2e (KM) TIMESCALE=%.2e (MYR)\n",
		  E->control.PID,E->monitor.length_scale,E->monitor.time_scale);
	  fprintf(fp,"# NODESX=%d NODESZ=%d NODESY=%d\n",E->mesh.nox,E->mesh.noz,E->mesh.noy);
	  fprintf(fp,"# ELAPSED_TIME=%.4e VELOCITY_SOLUTIONS=%d TIMESTEPS=%d\n",
		  E->monitor.elapsed_time,E->monitor.solution_cycles,E->advection.total_timesteps);
	  fprintf(fp,"#  Node   |     X      |     Z      |     Y      |");

	  for(number=1;number<=E->control.node_data.numb;number++) {
	    if(strstr(E->control.node_data.which_data_types,E->control.node_data.name[number]) != 0) {
	      E->control.node_data.numb_found++;
	      fprintf(fp,"      %s      |",E->control.node_data.name[number]);
	    }
	  }
	  fprintf(fp,"\n");
	    
	  for(i=1;i<=nno;i++) {
	    fprintf(fp,"  %06d   % .5e % .5e % .5e ",i,E->sx[1][i],E->sx[2][i],(dims==3)?E->sx[3][i]:0.0);
	    for(number=1;number<=E->control.node_data.numb;number++) {
	      if(strstr(E->control.node_data.which_data_types,E->control.node_data.name[number]) != 0)
		fprintf(fp," % .8e ",E->control.node_data.data[number][i]);
	    }
	    fprintf(fp,"\n") ;
	  }
	  fclose(fp);
	  if(E->control.COMPRESS){
	    /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,output_file);*/
	    sprintf(command,"%s -f  %s",E->control.gzip,output_file); /*DAS: 17-01-03*/
	    system(command);
	  }
	}
      }
      /* 2. Horizontal averages */

      if(strlen(E->control.haverage_data.which_data_types)!=0) {
	E->control.haverage_data.numb_found=0;
	sprintf(output_file,"%s.%05d.horiz_ave",E->control.data_file, record_everything ? file_number : 0); 	
	if ((fp=fopen(output_file,"w")) != NULL) {
	  /* print header */
	  fprintf(fp,"# PROCESS_ID=%06d\n",E->control.PID);
	  fprintf(fp,"# NODESZ=%d\n",E->mesh.noz);
	  fprintf(fp,"# ELAPSED_TIME=%.4e VELOCITY_SOLUTIONS=%d TIMESTEPS=%d\n",
		  E->monitor.elapsed_time,E->monitor.solution_cycles,E->advection.total_timesteps);
	  fprintf(fp,"# Layer |    Z     |");
	  
	  for(number=1;number<=E->control.haverage_data.numb;number++) {
	    if(strstr(E->control.haverage_data.which_data_types,E->control.haverage_data.name[number]) != 0) {
	      E->control.haverage_data.numb_found++;
	      fprintf(fp,"    %s    |",E->control.haverage_data.name[number]);
	    }
	  }
	  fprintf(fp,"\n"); 

	  /* Followed by the data ... */

	  for(i=1;i<=E->mesh.noz;i++){
	    fprintf(fp,"  %03d   % .3e",i,E->sx[2][i]);
	    for(number=1;number<=E->control.haverage_data.numb;number++) {
	      if(strstr(E->control.haverage_data.which_data_types,E->control.haverage_data.name[number]) != 0)
		fprintf(fp," % .5e ",E->control.haverage_data.data[number][i]);
	    }
	    fprintf(fp,"\n");
	  }
	    
	  fclose(fp);
	  if(E->control.COMPRESS){
	    /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,output_file);*/
	    sprintf(command,"%s -f  %s",E->control.gzip,output_file); /*DAS: 17-01-03*/
	    system(command);
	  } 
	}
      }

#if 1
    
      /* 3. Surface Observables */

      if(strlen(E->control.slice_data.which_data_types)!=0) {
	E->control.slice_data.numb_found=0;
	sprintf(output_file,"%s.%05d.observables",E->control.data_file, record_everything ? file_number : 0); 
	if ((fp=fopen(output_file,"w")) != NULL) {
	  /* print header */
	  fprintf(fp,"# PROCESS_ID=%06d TPG_SCALE=%.2e (M)\n",E->control.PID,E->monitor.tpgscale);
	  fprintf(fp,"# NODESX=%d NODESY=%d\n",E->mesh.nox,E->mesh.noy);
	  fprintf(fp,"# ELAPSED_TIME=%.4e VELOCITY_SOLUTIONS=%d TIMESTEPS=%d\n",
		  E->monitor.elapsed_time,E->monitor.solution_cycles,E->advection.total_timesteps);
	  fprintf(fp,"# Indx   |   X      |   Y      |");
	    
	  for(number=1;number<=E->control.slice_data.numb;number++) {
	   
	    if(strstr(E->control.slice_data.which_data_types,E->control.slice_data.name[number]) != 0) {
	      E->control.slice_data.numb_found++; 
	      fprintf(fp,"    %s    |",E->control.slice_data.name[number]);
	    }
	  }
	  fprintf(fp,"\n");

	  i=0;
	  for(k=1;k<=E->mesh.noy;k++)
	    for(j=1;j<=E->mesh.nox;j++) {
	      i++;
	      /* node location */
	      fprintf(fp,"  %06d  % .3e % .3e",i,
		      E->sx[1][1+(j-1)*E->mesh.noz],(E->mesh.nsd==3) ? E->sx[3][1+(k-1)*E->mesh.noz*E->mesh.nox]: 0.0) ;
		  
	      for(number=1;number<=E->control.slice_data.numb;number++) {
		if(strstr(E->control.slice_data.which_data_types,E->control.slice_data.name[number]) != 0) 
		  fprintf(fp," % .5e",E->control.slice_data.data[number][i]);
	      }
	      fprintf(fp,"\n");
	    }

	  fclose(fp);
	  if(E->control.COMPRESS){
	    /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,output_file);*/
	    sprintf(command,"%s -f  %s",E->control.gzip,output_file); /*DAS: 17-01-03*/
	    system(command);
	  } 
	}
      }
#endif
    }
    /* Data recorded at each timestep. Note that the organization is different
       in this part, the header is put out in the first timestep (only), then
       the columns of data  */

    if(strlen(E->control.time_data.which_data_types)!=0) {
      E->control.time_data.numb_found=0;
      sprintf(output_file,"%s.timelogs",E->control.data_file); 	
      if ((fp=fopen(output_file,(been_here==0) ? "w" :"a")) != NULL) {
	/* do header ... */
	if(been_here==0) {
	  fprintf(fp,"# PROCESS_ID=%06d LENGTH_SCALE=%.2e (KM) TIMESCALE=%.2e (MYR)\n",
		  E->control.PID,E->monitor.length_scale,E->monitor.time_scale);
	  fprintf(fp,"# Indx:tstep  |  Time    |");
	  
	  for(number=1;number<=E->control.time_data.numb;number++) {
	    if(strstr(E->control.time_data.which_data_types,E->control.time_data.name[number]) != 0) {
	      E->control.time_data.numb_found++;
	      fprintf(fp,"      %s     |",E->control.time_data.name[number]);
	    }
	  }
	  fprintf(fp,"\n");
	}
	fprintf(fp,"  %04d:%05d   %.4e",E->monitor.solution_cycles,E->advection.total_timesteps,E->monitor.elapsed_time);
	for(number=1;number<=E->control.time_data.numb;number++) {
	  if(strstr(E->control.time_data.which_data_types,E->control.time_data.name[number]) != 0)
	    fprintf(fp," % .7e ",*(E->control.time_data.data[number]));
	    }
	fprintf(fp,"\n");
	fclose(fp);
      } 
    }

   /* Data recorded at each timestep. Note that the organization is different
       in this part, the header is put out in the first timestep (only), then
       the columns of data  */

    if(checkpt_and_log)
      if(strlen(E->control.time_data.which_data_types)!=0) {
	E->control.time_data.numb_found=0;
	sprintf(output_file,"%s.checkptlogs",E->control.data_file); 	
	if ((fp=fopen(output_file,(been_here==0) ? "w" :"a")) != NULL) {
	  /* do header ... */
	  if(been_here==0) {
	    fprintf(fp,"# PROCESS_ID=%06d LENGTH_SCALE=%.2e (KM) TIMESCALE=%.2e (MYR)\n",
		    E->control.PID,E->monitor.length_scale,E->monitor.time_scale);
	    fprintf(fp,"# Indx:tstep  |  Time    |");
	    
	    for(number=1;number<=E->control.time_data.numb;number++) {
	      if(strstr(E->control.time_data.which_data_types,E->control.time_data.name[number]) != 0) {
		E->control.time_data.numb_found++;
		fprintf(fp,"      %s     |",E->control.time_data.name[number]);
	      }
	    }
	    fprintf(fp,"\n");
	  }
	  fprintf(fp,"  %04d:%05d   %.4e",E->monitor.solution_cycles,E->advection.total_timesteps,E->monitor.elapsed_time);
	  for(number=1;number<=E->control.time_data.numb;number++) {
	    if(strstr(E->control.time_data.which_data_types,E->control.time_data.name[number]) != 0)
	      fprintf(fp," % .7e ",*(E->control.time_data.data[number]));
	  }
	  fprintf(fp,"\n");
	  fclose(fp);
	} 
    }

    /* AVS .fld file to interpret the information in the node_data file */

    if(checkpt_and_log && E->mesh.nsd == 3 && E->control.AVS)  { 

      sprintf(output_file,"%s.%05d.AVS.fld",E->control.data_file,file_number);
      if ((fp=fopen(output_file,"w")) != NULL) {
	fprintf(fp,"# AVS format file pointing to data for run %s.%d\n",E->control.data_file,record_everything?file_number:0);
	fprintf(fp,"#\n");
	fprintf(fp,"ndim=3\n");
	fprintf(fp,"dim1=%d \ndim2=%d \ndim3=%d\n",E->mesh.noz,E->mesh.nox,E->mesh.noy);
	fprintf(fp,"nspace=3\n");
	fprintf(fp,"veclen=%d\n",E->control.node_data.numb_found); 
	fprintf(fp,"data=float\n");
	fprintf(fp,"field=rectilinear\n");
	fprintf(fp,"coord 2 file=%s.%05d.node_data filetype=ascii skip=4 offset=1 stride=%d\n",
		E->control.data_file,file_number,E->mesh.noz*(E->control.node_data.numb_found+4));
	fprintf(fp,"coord 3 file=%s.%05d.node_data filetype=ascii skip=4 offset=3 stride=%d\n",
		E->control.data_file,file_number,E->mesh.nox*E->mesh.noz*(E->control.node_data.numb_found+4));
	fprintf(fp,"coord 1 file=%s.%05d.node_data filetype=ascii skip=4 offset=2 stride=%d\n",
		E->control.data_file,file_number,E->control.node_data.numb_found+4);
	
	/* all the data that is stored can be read in (assuming AVS is ok about this) */
	for(j=1;j<=E->control.node_data.numb_found;j++)
	  fprintf(fp,"variable %d file=%s.%05d.node_data filetype=ascii skip=4 offset=%d stride=%d\n",
		  j,E->control.data_file,file_number,j+3,E->control.node_data.numb_found+4);
	fclose(fp);
      }
      
      sprintf(output_file,"%s.%05d.AVS_OBS.fld",E->control.data_file,file_number);
      if ((fp=fopen(output_file,"w")) != NULL) {
	fprintf(fp,"# AVS format file pointing to data for run %s.%d\n",E->control.data_file,file_number);
	fprintf(fp,"#\n");
	fprintf(fp,"ndim=2\n");
	fprintf(fp,"dim1=%d \ndim2=%d\n",E->mesh.nox,E->mesh.noy);
	fprintf(fp,"nspace=3\n"); 
	fprintf(fp,"veclen=%d\n",E->control.slice_data.numb_found); 
	fprintf(fp,"data=float\n");
	fprintf(fp,"field=rectilinear\n");
	fprintf(fp,"coord 1 file=%s.%05d.observables filetype=ascii skip=4 offset=1 stride=%d\n",
		E->control.data_file,file_number,(E->control.slice_data.numb_found+3));
	fprintf(fp,"coord 2 file=%s.%05d.observables filetype=ascii skip=4 offset=2 stride=%d\n",
		E->control.data_file,file_number,E->mesh.nox*(E->control.slice_data.numb_found+3));
	fprintf(fp,"coord 3 file=%s.%05d.observables filetype=ascii skip=4 offset=1 stride=0\n",
		E->control.data_file,file_number);
	
	/* all the data that is stored can be read in (assuming AVS is ok about this) */
	for(j=1;j<=E->control.slice_data.numb_found;j++)
	  fprintf(fp,"variable %d file=%s.%05d.observables filetype=ascii skip=4 offset=%d stride=%d\n",
		      j,E->control.data_file,file_number,j+2,E->control.slice_data.numb_found+3);
	
	fclose(fp);
      }	
    }
    
    /* AVS output of T,C,H2O data and coords in an interchange format  */
    

      if(checkpt_and_log && E->mesh.nsd == 3 && E->control.AVS)  { 
#if (HAVE_RPC_RPC_H && HAVE_RPC_XDR_H && HAVE_LIBXDR)
	data_length=1;

	sprintf(output_file,"%s.%05d.AVSI.fld",E->control.data_file,file_number);
	if ((fp=fopen(output_file,"w")) != NULL) {
            XDR xdrs;
            
	    fprintf(fp,"# AVS format file pointing to data for run %s.%d\n",E->control.data_file,file_number);
	    fprintf(fp,"#\n");
	    fprintf(fp,"ndim=3\n");
	    fprintf(fp,"dim1=%d \ndim2=%d \ndim3=%d\n",E->mesh.nox,E->mesh.noy,E->mesh.noz);
	    fprintf(fp,"nspace=3\n");
	    fprintf(fp,"veclen=%d\n",data_length); 
	    fprintf(fp,"data=xdr_float\n");
	    fprintf(fp,"field=rectilinear\n");
	    fprintf(fp,"\f\f");

	    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	  
	     /* Temperature information */

	    for(k=1;k<=E->mesh.noz;k++)
	      for(j=1;j<=E->mesh.noy;j++)
		for(i=1;i<=E->mesh.nox;i++) {
		  node = k + (i-1) * E->mesh.noz + (j-1) * E->mesh.noz * E->mesh.nox;
		    test = E->T[node];
		  xdr_float(&xdrs,&test);
	    }

	    /* X - coordinates */
	    for(i=1;i<=E->mesh.nox;i++) {
	      test = E->sx[1][1+(i-1)*E->mesh.noz];
	      xdr_float(&xdrs,&test);
	    }

	    /* Y - coordinates */
	    if(3==dims)
	      for(i=1;i<=E->mesh.noy;i++) {
		test = E->sx[3][1+(i-1)*E->mesh.noz*E->mesh.nox];
		xdr_float(&xdrs,&test);
	      }
	   
	    /* Z - coordinates */
	    for(i=1;i<=E->mesh.noz;i++) {
	      test = E->sx[2][E->mesh.noz] - E->x[2][i];
	      xdr_float(&xdrs,&test);
	    }

	    /* Finish up */

	    xdr_destroy(&xdrs);
	    fclose(fp);

	    if(E->control.COMPRESS){
	      /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,output_file);*/
	      sprintf(command,"%s -f  %s",E->control.gzip,output_file); /*DAS: 17-01-03*/
	      system(command);
	    }
	}
#endif /* HAVE_RPC_RPC_H && HAVE_RPC_XDR_H && HAVE_LIBXDR */
      }


      /* Particle based data output - binary output */

      if(checkpt_and_log && strlen(E->control.particle_data.which_data_types)!=0) {
	E->control.particle_data.numb_found=0;
	for(number=1;number<=E->control.particle_data.numb;number++) 
	  if(strstr(E->control.particle_data.which_data_types,E->control.particle_data.name[number]) != 0) {
	    E->control.particle_data.numb_found++;
	  }


	sprintf(output_file,"%s.%05d.particles",E->control.data_file, record_everything ? file_number : 0); 

	if ((fp=fopen(output_file,"w")) != NULL) {

	  /* Will we need to copy the data ? */

	  if(sizeof(standard_precision) != sizeof(float)) {
	    data_array = (float *) Malloc0((E->tracer.NUM_TRACERS + 1) * sizeof(float));
	  }
	  else
	    data_array = (float *) Malloc0((10) * sizeof(float));
	 
	  /* HEADER INFORMATION (followed by ^L to 
	     separate ascii information from binary data) */

	  fprintf(fp,"# FORMAT=binary PROCESS_ID=%06d ELAPSED_TIME=%.4e VELOCITY_SOLUTIONS=%d TIMESTEPS=%d\n",
		  E->control.PID,E->monitor.elapsed_time,E->monitor.solution_cycles,
		  E->advection.total_timesteps);

	  fprintf(fp,"# PARTICLES=%d DIMENSIONS=%d GEOMETRY=%s DATA_SIZE=%d \n",
		  E->tracer.NUM_TRACERS,
		  E->mesh.nsd,
		  E->control.GEOMETRY,
		  sizeof(float));

	  fprintf(fp,"# DATA_ELEMENTS=%d \n",E->control.particle_data.numb_found);
	  fprintf(fp,"# MATERIAL_ELEMENTS=%d \n",E->tracer.NUM_MATERIALS);


	  if(2==dims)
	    fprintf(fp,"# Mtrl=0 \n# Xloc=1 \n# Zloc=2 \n# Wght=3");
	  else /*RAA: 25/6/01, this part makes it necessary to create a 3D binary read perl file (extra line for Y)*/ 
	    fprintf(fp,"# Mtrl=0 \n# Xloc=1 \n# Zloc=2 \n# Yloc=3 \n# Wght=4");
	    
	  block = dims+2;
	  for(number=1;number<=E->control.particle_data.numb;number++) {
	    if(strstr(E->control.particle_data.which_data_types,E->control.particle_data.name[number]) != 0) {
	      fprintf(fp,"\n# %s=%d ",E->control.particle_data.name[number],block++);
	    }
	  }
	  fprintf(fp,"\n\f\n");


	  /* Now the body of the file - in binary format */
	  
	  /* 1. A single number "3.1" (AKA the 
	     Geological value of PI) which is used to 
	     check that the format of the file can be
	     read back in */

	  data_array[0] = 3.1;
	  fwrite(data_array,sizeof(float),1,fp);
	 

	  /* Material types */
	  fwrite(E->tracer.property_group,sizeof(int),E->tracer.NUM_TRACERS+1,fp);


	  /* Note, tracer 0 is not used but something needs to go in
	     to the file to pack this number when the tracers are read
	     back in ! */

	  data_array[0] = 0.0;

	  /* X coordinate */
	  if(sizeof(standard_precision) != sizeof(float)) {
	    for(k=1;k<=E->tracer.NUM_TRACERS;k++) {
	      data_array[k] = E->tracer.tx[k];
	    }
	    fwrite(data_array,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }
	  else {
	    fwrite(E->tracer.tx,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }
	  
	  /* Z coordinate */
	  if(sizeof(standard_precision) != sizeof(float)) {
	    for(k=1;k<=E->tracer.NUM_TRACERS;k++) {
	      data_array[k] = E->tracer.tz[k];
	    }
	    fwrite(data_array,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }
	  else {
	    fwrite(E->tracer.tz,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }

	  /* Y coordinate */
	  if(3==E->mesh.nsd)
	    if(sizeof(standard_precision) != sizeof(float)) {
	      for(k=1;k<=E->tracer.NUM_TRACERS;k++) {
		data_array[k] = E->tracer.ty[k];
	      }
	      fwrite(data_array,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	    }
	    else {
	      fwrite(E->tracer.ty,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }
	    
	  /* Weight */
	  
	  if(sizeof(standard_precision) != sizeof(float)) {
	    for(k=1;k<=E->tracer.NUM_TRACERS;k++) {
	      data_array[k] = E->tracer.tracer_weight_fn[E->mesh.levmax][k];
	    }
	    fwrite(data_array,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }
	  else {
	    fwrite(E->tracer.tracer_weight_fn[E->mesh.levmax],sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	  }



	  /* Next check each of the available data types to
	     see if it is required in the output file */

	  for(number=1;number<=E->control.particle_data.numb;number++) {   
	    if(strstr(E->control.particle_data.which_data_types,E->control.particle_data.name[number]) != 0) {

	      if(sizeof(standard_precision) != sizeof(float)) {
		for(k=1;k<=E->tracer.NUM_TRACERS;k++) {
		  data_array[k] = E->control.particle_data.data[number][k];
		}
		fwrite(data_array,sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	      }
	      else {
		fwrite(E->control.particle_data.data[number],sizeof(float),E->tracer.NUM_TRACERS+1,fp);
	      }
	    }
	  } 
	  

	  fclose(fp);

	  if(E->control.COMPRESS){
	    /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,output_file);*/
	    sprintf(command,"%s -f  %s",E->control.gzip,output_file); /*DAS: 17-01-03*/
	    system(command);
	  } 
	}
	
	free((void *) data_array); 
      }

      /* End of binary output */




      /* Particle-sampled profiles ... */

      if(checkpt_and_log && E->tracer.SAMPLE_PTS != 0) {
	sprintf(output_file,"%s.%05d.profiles",E->control.data_file, record_everything ? file_number : 0); 

	if ((fp=fopen(output_file,"w")) != NULL) {
	 
	  /* HEADER INFORMATION (4 lines again for convenience, followed by ^L to 
	     separate ascii information from binary data) */

	  fprintf(fp,"# FORMAT=ascii PROCESS_ID=%06d ELAPSED_TIME=%.4e VELOCITY_SOLUTIONS=%d TIMESTEPS=%d\n",
		  E->control.PID,E->monitor.elapsed_time,E->monitor.solution_cycles,
		  E->advection.total_timesteps);
	  fflush(fp);
	  fprintf(fp,"# PARTICLES=%d  \n",E->tracer.NUM_TRACERS);

	  fprintf(fp,"# DIMENSIONS=%d GEOMETRY=%s SP_SIZE=%s PROFILES=%d\n",
		  E->mesh.nsd,E->control.GEOMETRY,
		  (sizeof(standard_precision) == sizeof(float)) ? "float" : "double",
		  E->tracer.SAMPLE_PTS);

	  fflush(fp);

	  if(2==dims)
	    fprintf(fp,"# Xloc | Zloc |");
	  else 
            fprintf(fp,"#    Xloc     |     Zloc     |     Yloc     |");

	  for(number=0;number<E->tracer.SAMPLE_PTS;number++) {
	    fprintf(fp,"     %1d/%02d     |",number,E->tracer.sample_type[number]);
	  }
	  fprintf(fp,"\n");
	  fflush(fp);

	  for(i=0;i<=100;i++) {

            if(2==dims)  
              fprintf(fp," %.7e  %.7e ", E->tracer.sampled_data[0][i],E->tracer.sampled_data[1][i]);
            else if(3==dims)  /*RAA: 2/1/02, added 'if' above and here for 2D vs 3D coords */
              fprintf(fp," %.7e  %.7e  %.7e ", E->tracer.sampled_data[0][i],E->tracer.sampled_data[1][i],E->tracer.sampled_data[2][i]);
	    for(number=0;number<E->tracer.SAMPLE_PTS;number++) {
	      fprintf(fp," %.7e ",E->tracer.sampled_data[3+number][i]);
	    }
	    fprintf(fp,"\n");
	  }
	  fclose(fp);
	}
      }

      /* Sample point log files */


       for(number=0;number<E->tracer.SAMPLE_PTS;number++){
	sprintf(output_file,"%s.sample%d",E->control.data_file,number); 	
	if ((fp=fopen(output_file,(been_here==0) ? "w" :"a")) != NULL) {
	  /* do header ... */
	  if(been_here==0) {
	    fprintf(fp,"# PROCESS_ID=%06d LENGTH_SCALE=%.2e (KM) TIMESCALE=%.2e (MYR)\n",
		    E->control.PID,E->monitor.length_scale,E->monitor.time_scale);
	    fprintf(fp,"# Indx:tstep  |    Time     |       X       |       Z       |       Y       |     Value  \n");
	  }
	  fprintf(fp,"  %04d:%05d    %.4e     %.4e      %.4e      %.4e      %.4e\n",
		  E->monitor.solution_cycles,
		  E->advection.total_timesteps,
		  E->monitor.elapsed_time,
		  E->tracer.sample_x[number],
		  E->tracer.sample_z[number],
		  (3==dims) ? E->tracer.sample_y[number]: 0.0,
		  E->tracer.sample_value[number]
		  );

	  fclose(fp);
	} 
      }

      /* Grab the file root and call an optional external command
	 to process the data. The root file name is passed as an argument.
	 Ultimately it would be good to offer filename substitution at some
	 given point in the string.
	 
	 The command should not hold up CITCOM for long hence it 
	 is executed in the background of a subshell... of course, if CITCOM
	 gets too far ahead then there will be a problem.  */

    if(checkpt_and_log &&
       (0 < strlen(E->control.output_written_external_command)) ) {
	sprintf(command,"(%s %s %05d 1> /dev/null 2> /dev/null  &); ",E->control.output_written_external_command,
		E->control.data_file,file_number);
	fprintf(E->fp,"executing the following command to process output file ...\n\t %% %s\n",command);
	fprintf(stderr,"Executing %s\n",command);
	system(command);
   }

    /* Generate pixmap for 2D data */ /*RAA: 4/01/02, added nsd 3 below to see a 3D PPM file*/
    if((E->control.verbose || checkpt_and_log) && (E->mesh.nsd == 3 || E->mesh.nsd == 2) && !E->control.CYLINDER && !E->control.SPHERE) {
      sprintf(output_file,"%s.%05d.ppm",E->control.data_file,file_number);
      generate_2Ddata_pixmap(E,output_file);

      if(0 < strlen(E->control.ppm_written_external_command)) {
	sprintf(command,"(%s %05d %05d 1> /dev/null 2> /dev/null  &); ",E->control.ppm_written_external_command,
		E->control.data_file,file_number);
	sprintf(command,"(%s %05d %05d &); ",E->control.ppm_written_external_command,
		E->control.data_file,file_number);
	fprintf(E->fp,"executing the following command to process ppm file ...\n\t %% %s\n",command);
	fprintf(stderr,"Executing %s\n",command);
	system(command);
      }
    }
    been_here++;
    return; 
}

void add_variable_to_tracer_output(
				   struct All_variables *E,
				   standard_precision *VAR,
				   char *name
				   )
{
  static int visits = 0;

  if(0 == visits++)
    E->control.particle_data.numb=0;

  E->control.particle_data.data[++E->control.particle_data.numb] = VAR;
  strcpy(E->control.particle_data.name[E->control.particle_data.numb],name);

  return;
}

void add_variable_to_node_output(
				 struct All_variables *E,
				 standard_precision *VAR,
				 char *name
				 )
{
  static int visits = 0;

  if(0 == visits++)
    E->control.node_data.numb=0;

  E->control.node_data.data[++E->control.node_data.numb] = VAR;
  strcpy(E->control.node_data.name[E->control.node_data.numb],name);
  
  /*  fprintf(E->fp,"Adding variable %s to list for node date file -> %d\n",name,E->control.node_data.numb); */

  return;
}


void add_variable_to_haverage_output(
				     struct All_variables *E,
				     standard_precision *VAR,
				     char *name
				     )
{
  static int visits = 0;

  if(0 == visits++)
    E->control.haverage_data.numb=0;

  E->control.haverage_data.data[++E->control.haverage_data.numb] = VAR;
  strcpy(E->control.haverage_data.name[E->control.haverage_data.numb],name);

  return;
}

void add_variable_to_slice_output(
				  struct All_variables *E,
				  standard_precision *VAR,
				  char *name
				  )
{
  static int visits = 0;

  if(0 == visits++)
    E->control.slice_data.numb=0;

  E->control.slice_data.data[++E->control.slice_data.numb] = VAR;
  strcpy(E->control.slice_data.name[E->control.slice_data.numb],name);

  return;
}

void add_variable_to_timelog_output(
				    struct All_variables *E,
				    standard_precision *VAR,
				    char *name
				    )
{
  static int visits = 0;

  if(0 == visits++)
    E->control.time_data.numb=0;

  E->control.time_data.data[++E->control.time_data.numb] = VAR;
  strcpy(E->control.time_data.name[E->control.time_data.numb],name);

  return;
}
