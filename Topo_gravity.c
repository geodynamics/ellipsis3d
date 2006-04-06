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




#include <stdio.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

#define c_re(a) a.real
#define c_im(a) a.imag
  typedef struct compl {  double real;
			  double imag;    } COMPLEX;

extern standard_precision Zj0[1000],Zj1[1000];

void surface_observables(
			 struct All_variables *E,
			 int ii
			 )
{ 
  static int been_here=0;
  int i,j,m,el,imax,jmax,node,hnode;
  char filename[100];
  void gravity1();
   
  FILE *fp;
   
  standard_precision time,CPU_time();
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega;
  standard_precision lNx[4][ELNMAX+1];
  standard_precision number;
  standard_precision storage_timing,checkpt_timing;
 
  standard_precision *SzzT,*SzzN,*RhoN;
  standard_precision *temptpgk,*temptpgbk,*tempgrvk,*tempgrvbk,*tempgeok,*tempgeobk;
  
  void horiz_fft();
  void horiz_rft();
  void gravity();
  standard_precision return_layer_value();
    
  const int nox = E->mesh.nox;
  const int noy = E->mesh.noy;
  const int noz = E->mesh.noz;
 
  report(E,"Surface Observables");

  return;
  
  /* If timings have already been supplied, then use these.
     Otherwise compute them */
  
  if(E->control.storage_timing != 0.0) storage_timing = E->control.storage_timing;
  else if(E->advection.total_timesteps != 0) 
    storage_timing = (E->monitor.elapsed_time/E->advection.total_timesteps) * E->control.storage_timesteps;  
  
  if(E->control.checkpt_timing != 0.0) checkpt_timing = E->control.checkpt_timing;
  else if(E->advection.total_timesteps != 0) 
    checkpt_timing = (E->monitor.elapsed_time/E->advection.total_timesteps) * E->control.checkpt_timesteps;  
  
  if (((E->monitor.elapsed_time - E->monitor.time_of_last_store) <=  storage_timing) &&
      (E->monitor.solution_cycles != 1))
    return;
  
  /* If we're definitely going ahead, allocate memory */
  
  temptpgk=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  temptpgbk=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  tempgrvk=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  tempgeok=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  tempgrvbk=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  tempgeobk=(standard_precision *)Malloc0((1+E->mesh.equiv_even_nsf)*sizeof(standard_precision));
  
  SzzT=(standard_precision *) Malloc0((1+E->tracer.NUM_TRACERS) * sizeof(standard_precision));
  SzzN=(standard_precision *) Malloc0((1+E->mesh.nno) * sizeof(standard_precision));
  RhoN=(standard_precision *) Malloc0((1+E->mesh.nno) * sizeof(standard_precision));
  
  /* Get surface topography from Szz at particles */
  
  report(E,"Sigma zz (tracers)");
  
    /* 1 get Szz */
  
  for(m=1;m<=E->tracer.NUM_TRACERS;m++) {
    if(E->tracer.property_group[m] < 0)
	continue;
    
    el = E->tracer.tracer_elt[E->mesh.levmax][m];
    eta1 = E->tracer.eta1[E->mesh.levmax][m];
    eta2 = E->tracer.eta2[E->mesh.levmax][m];
    if(3 == E->mesh.nsd) /*RAA: 6/11/02, added this*/
      eta3 = E->tracer.eta3[E->mesh.levmax][m];

    get_global_v_x_shape_fn(E,el,lNx,&dOmega,eta1,eta2,eta3,E->mesh.levmax); 
    
    SzzT[m] = 0.0;
    for(j=1;j<=enodes[E->mesh.nsd];j++)
      SzzT[m] += 2.0 * E->tracer.Visc[m] * E->V[2][E->ien[el].node[j]] * lNx[2][j]; 
  }
  
  /* 2. Send to nodes */
  
  report(E,"Sigma zz to nodes");
  
  gs_tracers_to_nodes(E,SzzN,NULL,NULL,NULL,SzzT,E->mesh.levmax,0);
  
  /* 3. Extract surface topo */

  if(E->data.grav_acc == 0.0)
    report(E,"Topography cannot be calculated because gravitational acceleration is zero");
  else {
    report(E,"Topography from Szz");
  
    for(i=1;i<=E->mesh.nox;i++) {
      node = 1 + (i-1) * E->mesh.noz;
      E->slice.tpg[i] = -SzzN[node] + E->nQ[node] / E->data.grav_acc;
      node += E->mesh.noz - 1;
      E->slice.tpgb[i] = -SzzN[node] + E->nQ[node] / E->data.grav_acc;
    }
  }
  
   /* Get gravity signal due to internal mass distribution and topography */
  
  /* 1. Get nodal densities */
  
   report(E,"Nodal density");
   
   gs_tracers_to_nodes(E,RhoN,NULL,NULL,NULL,E->tracer.weight,E->mesh.levmax,0);
   
   /* 2. Get FFT of topography */
   
   report(E,"FFT's");
   
   horiz_fft(E,E->slice.tpg,E->slice.tpgk);
   horiz_fft(E,E->slice.tpgb,E->slice.tpgbk);
   
   /* 3. Get gravity signal */
   
   report(E,"gravity");

   gravity1(E,RhoN,E->slice.tpg,E->slice.tpgb);
   
   /* Scale topography (currently it is delta M) */
   
   report(E,"gravity 2");
   
   for(i=1;i<=E->mesh.nox;i++) {
     node = 1 + (i-1) * E->mesh.noz;
     if(RhoN[node] != 0.0)
       E->slice.tpg[i] /= (RhoN[node]);
     node += E->mesh.noz - 1;
     if(RhoN[node] != 0.0)
       E->slice.tpgb[i] /= (RhoN[node]);
   }

#if 0
   if(E->control.CART2D || E->control.CART3D) {
     horiz_fft(E,E->slice.tpg,E->slice.tpgk);
     horiz_fft(E,E->slice.tpgb,E->slice.tpgbk);
     E->slice.tpgk[1]=E->slice.tpgbk[1]=0.0;
     /* horiz_rft(E,E->slice.tpgk,E->slice.tpg);      
	horiz_rft(E,E->slice.tpgbk,E->slice.tpgb); */
     
     /* All these fft's get pretty expensive especially in 3D so
	only get gravity for cases which are recorded permanently */
     
     if ((ii < E->control.record_all_until) ||
	 ((ii % E->control.storage_timesteps) == 0)) {
       time=CPU_time();
       gravity(E); 
       report(E,"Obtained  gravity");
       
       /* get magnitude of long wavelength gravity/topography */
       
       for(i=1;i<=nox;i++)
	 for(j=1;j<=noy;j++)  {
	   hnode = i + (j-1)*E->mesh.nox;
	   E->monitor.tpgkmag = max(E->monitor.tpgkmag,fabs(E->slice.tpgk[hnode]));
	   E->monitor.grvkmag = max(E->monitor.grvkmag,fabs(E->slice.grvk[hnode]));
	 }
     }
   }
#endif

    been_here++;

    free((void *) SzzT);
    free((void *) SzzN);
    free((void *) RhoN);
    free(temptpgk);
    free(temptpgbk);
    free(tempgrvk);
    free(tempgeok);
    free(tempgrvbk);
    free(tempgeobk);
    return;
}

/* ======================================================
   Gravity routine accepting specified buoyancy field and 
   topography.
   ====================================================== */
void gravity1(
     struct All_variables *E,
     standard_precision *buoyancy,
     standard_precision *topo,
     standard_precision *topob
)

{  /* Obtain gravity */

  int z,i,j,lnode,node;
  standard_precision kx,ky,k;
  standard_precision *buoy1,*buoy2;
  standard_precision *buoy1ft,*buoy2ft;
  
  standard_precision HaveG,HaveGO;
  standard_precision return_layer_value();

  void horiz_fft();
  void horiz_rft();
  
  buoy1 = (standard_precision *)Malloc0((5+E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
  buoy2 = (standard_precision *)Malloc0((5+E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
  buoy1ft = (standard_precision *)Malloc0((5+E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
  buoy2ft = (standard_precision *)Malloc0((5+E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
   
  for(i=0;i<=E->mesh.nsf+1;i++)
    E->slice.grvk[i] = E->slice.grvbk[i] = 0.0;
  
  for(z=1;z<E->mesh.noz;z++) {
      for(j=1;j<=E->mesh.noy;j++)
	  for(i=1;i<=E->mesh.nox;i++)	  { 
	      lnode = i + (j-1) * E->mesh.nox;
	      node = z + (i-1)*E->mesh.noz + (j-1)*E->mesh.nox*E->mesh.noz;
	      
	      buoy1[lnode] = buoyancy[node];
	      buoy2[lnode] = buoyancy[node+1];
	   
	      /* layer density anomalies loaded */
	  }

      horiz_fft(E,buoy1,buoy1ft);
      horiz_fft(E,buoy2,buoy2ft);
      
      for(j=1;j<=E->mesh.noy;j++)
	for(i=1;i<=E->mesh.nox;i++)  {
	    lnode = i + (j-1) * E->mesh.nox;
	    kx = (i-1)/E->x[1][E->mesh.nno];
	    if(3==E->mesh.nsd)  
		ky = (j-1)/E->x[3][E->mesh.nno];
	    else
		ky = 0.0;
	    k = sqrt((double)(kx*kx+ky*ky));

	    if(!E->control.AXI)
		E->slice.grvk[lnode] += 6.67e-11 * 
		    (buoy1ft[lnode] * exp(-E->x[2][z] * M_PI * k) +
		     buoy2ft[lnode] * exp(-E->x[2][z+1] * M_PI * k))*0.5*(E->x[2][z+1]-E->x[2][z]); 
	    else
		E->slice.grvk[lnode] += 6.67e-11 * 
		    (buoy1ft[lnode] * exp(-E->x[2][z] * Zj1[(int)kx]) +
		     buoy2ft[lnode] * exp(-E->x[2][z+1] * Zj1[(int)kx] ))*0.5*(E->x[2][z+1]-E->x[2][z]);
	  }
      /* Gravity signals for this depth calculated */
    }
  /* Add topography to  gravity signal  */

  for(j=1;j<=E->mesh.noy;j++)
    for(i=1;i<=E->mesh.nox;i++) {
	lnode = i + (j-1) * E->mesh.nox;
	kx = (i-1)/E->x[1][E->mesh.nno];
	if(E->mesh.nsd==3) 
	    ky = (j-1)/E->x[3][E->mesh.nno];
	else
	    ky = 0.0;

	k = 1.0e-32+sqrt((double)kx*kx+ky*ky);
	E->slice.grvtk[lnode] = E->slice.grvk[lnode] + 6.67e-11 * E->slice.tpgk[lnode];
	E->slice.grvbk[lnode] = E->slice.grvtk[lnode] + 6.67e-11 * 
	  E->slice.tpgbk[lnode] * exp(-E->x[2][E->mesh.nno] * (E->control.AXI ? (Zj1[(int)kx]) : M_PI *k));
      
	if(E->control.AXI) { 
	  if(((int)kx) != 0) {
	    E->slice.geobk[lnode] = E->slice.grvbk[lnode] /Zj1[(int)kx];
	    E->slice.geotk[lnode] = E->slice.grvtk[lnode] /Zj1[(int)kx];
	    E->slice.geok[lnode]  = E->slice.grvk[lnode] /Zj1[(int)kx];
	  }
	}
	else
	  if(k > 1.0e-10) {
	    E->slice.geotk[lnode] = E->slice.grvtk[lnode] / (M_PI *k);
	    E->slice.geobk[lnode] = E->slice.grvbk[lnode] / (M_PI *k);
	    E->slice.geok[lnode]  = E->slice.grvk[lnode] / (M_PI *k);
	  }
    }
 
  E->slice.grvk[1] = E->slice.grvbk[1]  =
    E->slice.grvtk[1] = E->slice.geotk[1] =
    E->slice.geok[1] = E->slice.geobk[1] = 0.0;

  horiz_rft(E,E->slice.grvk,E->slice.grv);
  horiz_rft(E,E->slice.geok,E->slice.geo);
  horiz_rft(E,E->slice.grvbk,E->slice.grvb);
  horiz_rft(E,E->slice.geobk,E->slice.geob);
  horiz_rft(E,E->slice.grvtk,E->slice.grvt);
  horiz_rft(E,E->slice.geotk,E->slice.geot);

  HaveG  = return_layer_value(E,E->slice.grv,1,1);
  HaveGO = return_layer_value(E,E->slice.geo,1,1);
  
  /* for(i=1;i<=E->mesh.nsf;i++)
    { E->slice.grv[i]  -= HaveG;
      E->slice.geo[i]  -= HaveGO; }
      
      HaveG  = return_layer_value(E,E->slice.grvb,1,1);
      HaveGO = return_layer_value(E,E->slice.geob,1,1);
  
      for(i=1;i<=E->mesh.nsf;i++)
      { E->slice.grvb[i]  -= HaveG;
      E->slice.geob[i]  -= HaveGO; } */

  free((void *) buoy1);
  free((void *) buoy2);
  free((void *) buoy1ft);
  free((void *) buoy2ft);
  return;
}

 /* ==================================================
    Calculate horizontal fft including the symmetry of
    the reflection boundary conditions at the sides
    ================================================== */


void horiz_fft(
     struct All_variables *E,
     standard_precision *f,
     standard_precision *F
)
{
  void horiz_1d_fft();

  double *line,*lineft,*ff,*FF;
  int i,j,node2d;

  line = (double *)Malloc0((1+max(E->mesh.nox,E->mesh.noy))*sizeof(double));
  lineft = (double *)Malloc0((1+max(E->mesh.nox,E->mesh.noy))*sizeof(double));
  ff  = (double *)Malloc0((1+E->mesh.nsf)*sizeof(double));
  FF = (double *)Malloc0((1+E->mesh.nsf)*sizeof(double));
   

  for(i=1;i<=E->mesh.nsf;i++) {
    ff[i]  = (double) f[i];
  }

  if (E->mesh.nsd==2) 
    horiz_1d_fft(E,ff,FF,1); /* straightforward 1d fft */    


  else {		/*1*/
    for(i=1;i<=E->mesh.noy;i++){		/*2*/
	 for(j=1;j<=E->mesh.nox;j++){
	   node2d = j + (i-1) * E->mesh.nox;
	      line[j] = ff[node2d]; 
		}
	  horiz_1d_fft(E,line,lineft,1);
	  for(j=1;j<=E->mesh.nox;j++) {
		 node2d = j + (i-1) * E->mesh.nox;
	      FF[node2d] = lineft[j]; 
		}
	} 		/*2*/
				/* note change of meanings i,j */
     
       for(i=1;i<=E->mesh.nox;i++){			/*2*/
	   for(j=1;j<=E->mesh.noy;j++) {
		 node2d = i + (j-1) * E->mesh.nox;
	      line[j] = FF[node2d];
		}
	  horiz_1d_fft(E,line,lineft,3);
	  for(j=1;j<=E->mesh.noy;j++) {
	    node2d = i + (j-1) * E->mesh.nox;
	      FF[node2d] = lineft[j]; 
	    }
	 } 						/*2*/
    }						/*1*/
  

  for(i=1;i<=E->mesh.nsf;i++) {
    F[i]  = (standard_precision) FF[i];    
  }

  free((void *) line );
  free((void *) lineft );
  free((void *) ff );
  free((void *) FF );
  return;
}

void horiz_rft(
     struct All_variables *E,
     standard_precision *f,
     standard_precision *F
)  
{ void horiz_1d_rft();
  double *line,*lineft,*ff,*FF;
  int i,j,node2d;

  line = (double *)Malloc0((1+max(E->mesh.nox,E->mesh.noy))*sizeof(double));
  lineft = (double *)Malloc0((1+max(E->mesh.nox,E->mesh.noy))*sizeof(double));
  ff  = (double *)Malloc0((1+E->mesh.nsf)*sizeof(double));
  FF = (double *)Malloc0((1+E->mesh.nsf)*sizeof(double));
     
  for(i=1;i<=E->mesh.nsf;i++)
    ff[i]  = (double) f[i];

  if (E->mesh.nsd==2) /* straightforward 1d fft */
    { horiz_1d_rft(E,ff,FF,1); }

  else
    { for(i=1;i<=E->mesh.noy;i++)
	{ for(j=1;j<=E->mesh.nox;j++)
	    { node2d = j + (i-1) * E->mesh.nox;
	      line[j] = f[node2d]; }
	  horiz_1d_rft(E,line,lineft,1);
	  for(j=1;j<=E->mesh.nox;j++)
	    { node2d = j + (i-1) * E->mesh.nox;
	      FF[node2d] = lineft[j];
	    }
	}
     
      for(i=1;i<=E->mesh.nox;i++) /* note change of meanings i,j */
	{  for(j=1;j<=E->mesh.noy;j++)
	    { node2d = i + (j-1) * E->mesh.nox;
	      line[j] = FF[node2d]; }
	  horiz_1d_rft(E,line,lineft,3);
	  for(j=1;j<=E->mesh.noy;j++)
	    { node2d = i + (j-1) * E->mesh.nox;
	      FF[node2d] = lineft[j]; 
	    }
	 }
      
    }
 
  for(i=1;i<=E->mesh.nsf;i++)
    if(fabs(FF[i]) < 1.0e15)
      F[i]  = (standard_precision) FF[i];
    
  free((void *) line );
  free((void *) lineft );
  free((void *) ff );
  free((void *) FF );
  return;
}


void horiz_1d_fft(
     struct All_variables *E,
     double *f,
     double *F, /* capitals are transforms */
     int dirn
)

{ 
    static double *freg;
    static COMPLEX *packed,*PACKED;
    static int points,been_here = 0;
    int i;
    void horiz_respacing();
    void fft();
    void pack_fft();
    void unpack_fft();
    void hankel_transform();
    
    if (!been_here++) {
	i = E->mesh.nnx[1]+E->mesh.nnx[3];
	freg = (double *)Malloc0((1+i)*sizeof(double));
	packed = (COMPLEX *)Malloc0((10+2*i)*sizeof(COMPLEX));
	PACKED = (COMPLEX *)Malloc0((10+2*i)*sizeof(COMPLEX));
  }
	
    points = 2*E->mesh.nnx[dirn]-2; /* ??? */
  
    if(E->control.AXI)	hankel_transform(E,f,F);

   else {
     horiz_respacing(E,f,freg,dirn);
     pack_fft(E,freg,packed,dirn);
     fft(packed,points,PACKED);
     unpack_fft(E,PACKED,F,dirn);
   }
   return;
}

void horiz_1d_rft(
    struct All_variables *E,
    double *F,
    double *f, /* capitals are transforms */
    int dirn
)
{
  static COMPLEX *packed,*PACKED;
  static int points,been_here = 0;
  int i;
  void horiz_respacing();
  void fft();
  void pack_rft();
  void unpack_rft();
  void inverse_hankel_transform();

  if (!been_here++)  {
      i = E->mesh.nnx[1]+E->mesh.nnx[3];
      packed = (COMPLEX *)Malloc0((3+2*i)*sizeof(COMPLEX));
      PACKED = (COMPLEX *)Malloc0((3+2*i)*sizeof(COMPLEX));
  
  }
 
  points = 2*E->mesh.nnx[dirn]-2;  
 
  if(E->control.AXI)
    inverse_hankel_transform(E,F,f);
  else { 
      pack_rft(E,F,PACKED,dirn);
      rft(PACKED,points,packed);
      unpack_rft(E,packed,f,dirn); 
  }

  return;
}

void pack_fft(
     struct All_variables *E,
     double *f,
     COMPLEX *p,
     int dirn
)
{  
    int i;
    
    for(i=0;i<=E->mesh.nnx[dirn];i++)
	c_im(p[i])  = c_im(p[i+E->mesh.nnx[dirn]]) = 0.0;

    if(1==dirn && !E->mesh.periodic_x)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[2*E->mesh.nnx[dirn]-2-i]) = f[i+1];
    
    if(1==dirn && E->mesh.periodic_x)
	    for(i=0;i<E->mesh.nnx[dirn];i++)
		c_re(p[i]) = c_re(p[E->mesh.nnx[dirn]+i-1]) = f[i+1];
   
    if(3==dirn && !E->mesh.periodic_y)
	    for(i=0;i<E->mesh.nnx[dirn];i++)
		c_re(p[i]) = c_re(p[2*E->mesh.nnx[dirn]-2-i]) = f[i+1];

    if(3==dirn && E->mesh.periodic_y)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[E->mesh.nnx[dirn]+i-1]) = f[i+1];
    
   return; }

void pack_rft(
     struct All_variables *E,
     double *F,
     COMPLEX *p,
     int dirn
)
{ 
    int i;

    for(i=0;i<=E->mesh.nnx[dirn];i++)
	c_im(p[i])  = c_im(p[i+E->mesh.nnx[dirn]]) = 0.0;
    
     
    for(i=0;i<E->mesh.nnx[dirn];i++)
	c_re(p[i]) = c_re(p[2*E->mesh.nnx[dirn]-2-i]) = F[i+1];
    

    if(1==dirn && !E->mesh.periodic_x)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[2*E->mesh.nnx[dirn]-2-i]) = F[i+1];
    
    if(1==dirn && E->mesh.periodic_x)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[E->mesh.nnx[dirn]+i-1]) = F[i+1];
   
    if(3==dirn && !E->mesh.periodic_y)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[2*E->mesh.nnx[dirn]-2-i]) = F[i+1];
    
    if(3==dirn && E->mesh.periodic_y)
	for(i=0;i<E->mesh.nnx[dirn];i++)
	    c_re(p[i]) = c_re(p[E->mesh.nnx[dirn]+i-1]) = F[i+1];
   

   return; }

void unpack_fft(
     struct All_variables *E,
     COMPLEX *p,
     double *F,
     int dirn
)
{ 
    int i,nodes;
    nodes=E->mesh.nnx[dirn];
    
    if(1==dirn && !E->mesh.periodic_x)
      for(i=0;i<nodes;i++) {	
	F[i+1] =  2.0*c_re(p[i]);
      }
    
    if(1==dirn && E->mesh.periodic_x) {
      for(i=0;i<nodes-1;i++) {
	F[i+1] =  2.0*c_re(p[2*i]);
      }
      F[nodes] = F[1];
    }

    if(3==dirn && !E->mesh.periodic_y)
	for(i=0;i<nodes;i++)
	    F[i+1] =  2.0*c_re(p[i]); 	   
    
    if(3==dirn && E->mesh.periodic_y) {
	for(i=0;i<nodes-1;i++)
	    F[i+1] =  2.0*c_re(p[2*i]);
	F[nodes] = F[1];
    } 

    return;
 }


void unpack_rft(
     struct All_variables *E,
     COMPLEX *p,
     double *f,
     int dirn
)
{   
    int i,nodes;

    nodes=E->mesh.nnx[dirn];


    if(1==dirn && !E->mesh.periodic_x)
	for(i=0;i<nodes;i++)
	    f[i+1] = 0.5*c_re(p[i]); 
    
    if(1==dirn && E->mesh.periodic_x) {
	for(i=0;i<nodes-1;i++)
	    f[i+1] =  0.5*c_re(p[2*i]);
	f[nodes] = f[1];
    }
    

    if(3==dirn && !E->mesh.periodic_y)
	for(i=0;i<nodes;i++)
	    f[i+1] = 0.5*c_re(p[i]); 
    
    if(3==dirn && E->mesh.periodic_y) {
	for(i=0;i<nodes-1;i++)
	    f[i+1] =  0.5*c_re(p[2*i]);  
	f[nodes] = f[1];
    }
    return; }

void horiz_respacing(
     struct All_variables *E,
     double *f,
     double *fr,
     int dirn
)
{ standard_precision spacing,currentx;
  int i,bigger,smaller,node_inc;

  /* x dirn increments by noz between sweeps, ydirn increments by nox*noz 
     this gives locations along the axis. (z dirn would be inc of 1) */

  node_inc = E->mesh.noz * ((dirn==3) ? E->mesh.nox : 1);

  spacing = (E->x[dirn][E->mesh.nno] - E->x[dirn][1])/((standard_precision) E->mesh.nnx[dirn] - 1.0);
  fr[1] = f[1];
  fr[E->mesh.nnx[dirn]] = f[E->mesh.nnx[dirn]];

  for(i=2;i<E->mesh.nnx[dirn];i++)
    { currentx = ((standard_precision)i - 1.0) * spacing + E->x[dirn][1];
      smaller = 1;

      while (currentx >= E->x[dirn][1+(smaller-1)*node_inc]) 
	smaller++;
      
      bigger = smaller - 1;

      fr[i] = (f[bigger] +
	       (f[smaller]-f[bigger]) *
	       (currentx - E->x[dirn][1+(bigger-1)*node_inc])/
	       (E->x[dirn][1+(smaller-1)*node_inc] - E->x[dirn][1+(bigger-1)*node_inc]));

      }	/* next point */ 
  return;
}
