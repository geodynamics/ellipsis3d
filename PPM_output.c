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


  /*  */ 

/* Routines to produce ppm format pixmaps of 
   2D citcom data during runs for diagnostic
   purposes and to produce movies.             */


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


#define MAX_CMAP_SECTIONS 99

struct colormap { 
  int sections;
  int background_R;   
  int background_G;   
  int background_B;   
  float low[MAX_CMAP_SECTIONS];
  float opacity_low[MAX_CMAP_SECTIONS];
  int low_R[MAX_CMAP_SECTIONS];
  int low_G[MAX_CMAP_SECTIONS];
  int low_B[MAX_CMAP_SECTIONS];
  float high[MAX_CMAP_SECTIONS];
  float opacity_high[MAX_CMAP_SECTIONS];
  int high_R[MAX_CMAP_SECTIONS];
  int high_G[MAX_CMAP_SECTIONS];
  int high_B[MAX_CMAP_SECTIONS];
};

struct ppm {
  int x;
  int y;
  int *R;
  int *G;
  int *B;
};

/* test routine to create a temperature/viscosity/etc pixmap (no interpolation) */

void generate_2Ddata_pixmap(
			    struct All_variables *E,
			    char *name
			    )
{
  int PPM_plot_number;
  int i,j;
  int newX,newZ,offX,offZ;
  int border,vel_plot_height,vel_plot_width;
  int colorbar_width;
  
  /*RAA: 7/1/02, added these 3 lines for multiple tracer plots (separate y-slices) per ppm file*/
  int newX2,newZ2,offX2,offZ2,newX3,newZ3,offX3,offZ3;
  int border2,vel_plot_height2,vel_plot_width2,yslices;
  standard_precision aspect2;

  char name_plus[256];

  struct colormap cmap;
  struct ppm pixmap;
  
  standard_precision aspect;
  standard_precision *interp_data,*coarse_ndata,*coarse_edata;

  int which_colormap();
  void output_ppm();
  void simple_interp_2D();
  void horiz_respacing_n();
  void add_field_to_pixmap();
  void add_surf_graph_to_pixmap();
  void add_tracers_to_pixmap();
  void add_border_to_pixmap();
  void add_cross_sections_to_pixmap();

#if defined(__uxp__)
  return;
#endif



  vel_plot_width = vel_plot_height = 0;

  /*RAA: 4/01/02, to put multiple plots per ppm file - added the 4 lines below
       but for 3D since x:z  aspect ratios != 1 give Seg Faults, do 1 slice for them*/
  if(E->mesh.nsd == 2 || (E->mesh.nsd==3 && E->x[1][E->mesh.nno]!=E->x[2][E->mesh.nno]))
     yslices = 1;
  else if (E->mesh.nsd==3 && E->x[1][E->mesh.nno]==E->x[2][E->mesh.nno])
     yslices = 3;

  for(PPM_plot_number=0;PPM_plot_number<E->control.PPM_files;PPM_plot_number++) {

    /* Size of the new pixmap. Fix the file size (i.e. area) 
       to be manageable - 150000 pixels give or take*/

    if(E->control.PPM_surf_plot[PPM_plot_number])
      vel_plot_height=max(35,(int) (E->control.PPM_height * 0.15));
    else
      vel_plot_height=0;
    
    for(i=0;i<E->tracer.SAMPLE_PTS;i++) { /*RAA: new stuff here*/
      if(! (E->tracer.sample_in_plot[i] == PPM_plot_number))
	continue;
	 
      if (yslices==1) { /*smaller plot for 1 slice ppms*/
         if (E->tracer.sample_direction[i] == 1) {
	   vel_plot_height=max(35,(int) (E->control.PPM_height * 0.15));
         }
         else if (E->tracer.sample_direction[i] == 2)  {
	   vel_plot_width=max(75,(int) (E->control.PPM_height * 0.3));
         }
         /*RAA: add this part to  help get y sections in a different location*/
         else if (E->mesh.nsd==3 && E->tracer.sample_direction[i] == 3)  {
	   vel_plot_height=max(25,(int) (E->control.PPM_height * 0.15));
         }
      }
      else {     /*more than 1 yslice (for 1:1 aspect ratio w/ 3D problem)*/
         if (E->tracer.sample_direction[i] == 1) {
	   vel_plot_height=max(165,(int) (E->control.PPM_height * 0.15));
         }
         else if (E->tracer.sample_direction[i] == 2)  {
	   vel_plot_width=max(241,(int) (E->control.PPM_height * 0.3));
         }
         /*RAA: add this part to  help get y sections in a different location*/
         else if (E->tracer.sample_direction[i] == 3)  {
	   vel_plot_height=max(110,(int) (E->control.PPM_height * 0.15));
         }
      }
    }

    newZ=min(450,E->control.PPM_height);

    border = 15;
    aspect = (E->x[1][E->mesh.nno]-E->x[1][1]) / (E->x[2][E->mesh.nno]-E->x[2][1]);

    newX = (int) (aspect * newZ);
     
    /* set pixmap dimensions */ 
    pixmap.x = max(newX,E->control.PPM_aspect*newZ) + 3 * border + vel_plot_width;
    pixmap.y = newZ + 3 * border + vel_plot_height;
    
    /* allocate mem to new pixmap */  
    pixmap.R = (int *) Malloc0((pixmap.x * pixmap.y + 1) * sizeof(int));
    pixmap.G = (int *) Malloc0((pixmap.x * pixmap.y + 1) * sizeof(int));
    pixmap.B = (int *) Malloc0((pixmap.x * pixmap.y + 1) * sizeof(int));
    interp_data = (standard_precision *) Malloc0((newX*newZ+1) * sizeof(standard_precision));

    /* set to bg color */
    for(i=0;i<=pixmap.x * pixmap.y;i++) {
      pixmap.R[i] = 255;
      pixmap.G[i] = 255;
      pixmap.B[i] = 255;
    }

    /* Main plot with tracers */

    if(E->control.verbose) 
      fprintf(stderr,"Build tracer pixmap\n");

    offZ = 2 * border + vel_plot_height;
    offX = border;
    add_tracers_to_pixmap(E,E->tracer,&pixmap,newX,newZ,offX,offZ,0,PPM_plot_number,1,yslices);

    /*RAA: this approach obviously not too elegant, but should do for triage*/
    if (E->mesh.nsd==3 && yslices==2) {
      fprintf(stderr,"RAA: Adding another box for 2nd y-slice\n");
      offZ2 = 0.00 * border + 0.05*vel_plot_height;
      offX2 = 10*border;
      add_tracers_to_pixmap(E,E->tracer,&pixmap,newX,newZ,offX2,offZ2,0,PPM_plot_number,2,yslices);
      add_border_to_pixmap(E,&pixmap,newX,newZ,offX2,offZ2,50,50,50); 
    }  
    else if (E->mesh.nsd==3 && yslices==3) {
      fprintf(stderr,"RAA: Adding another box for 2nd y-slice\n");
      offZ2 = -0.00 * border + 0.05*vel_plot_height;
      offX2 = 10*border;
      add_tracers_to_pixmap(E,E->tracer,&pixmap,newX,newZ,offX2,offZ2,0,PPM_plot_number,2,yslices);
      add_border_to_pixmap(E,&pixmap,newX,newZ,offX2,offZ2,50,50,50); 

      fprintf(stderr,"RAA: Adding another box for 3rd y-slice\n");
      offZ3 = offZ2;
      offX3 = 18*border + 7;
      add_tracers_to_pixmap(E,E->tracer,&pixmap,newX,newZ,offX3,offZ3,0,PPM_plot_number,3,yslices);
      add_border_to_pixmap(E,&pixmap,newX,newZ,offX3,offZ3,50,50,50); 
    }

    if(E->control.verbose)
      fprintf(stderr,"Build tracer pixmap ... done\n");

    /* Add outline to plot region of pixmap (RAA: slices 2 and 3 done above)*/
    add_border_to_pixmap(E,&pixmap,newX,newZ,offX,offZ,50,50,50); 
 
    /* Add sample cross-sections to pixmap */

       if(E->control.verbose)
	 fprintf(stderr,"Cross sections to pixmap\n");

    /*RAA: added some stuff here to try and get y-dir plots in different location - NO IMPACT*/
    /*RAA: 22/05/02, but this new stuff may cause a Seg fault, because of i, so comment it out for now.*/
 /*  fprintf(stderr,"Gus, debug - location 0 i: %d \n",i);   */
/*	 fprintf(stderr,"Gus, debug - location 0 i: %d  dir: %d\n",i,E->tracer.sample_direction[i]);  */
/*    if (E->tracer.sample_direction[i] != 3)  {
      offZ = 2 * border + vel_plot_height;
      offX = 2 * border + vel_plot_width;
    }
    else if (E->tracer.sample_direction[i] == 3)  {
      offZ = 3 * border + vel_plot_height;
      offX = 3 * border + vel_plot_width;
    }
*/
    offZ = 2 * border + vel_plot_height;
    offX = 2 * border + vel_plot_width;
    add_cross_sections_to_pixmap(E,&pixmap,border,newX,newZ,
			    offX,offZ,vel_plot_height,vel_plot_width,
                            PPM_plot_number);
   
       if(E->control.verbose)
	 fprintf(stderr,"Cross sections to pixmap ... done \n");

	 
    /* Add  a surface plot pixmap */

    if(0 && (E->control.PPM_surf_plot[PPM_plot_number] <= 5) && (E->control.PPM_surf_plot[PPM_plot_number]> 0)) {
      switch(E->control.PPM_surf_plot[PPM_plot_number]) {
      case 1:
	horiz_respacing_n(E,E->slice.vxsurf[1],interp_data,1,newX); 
	for(i=1;i<=newX;i++)
	  interp_data[i] *= 31.5576e6; /* convert to m/yr */
	break;
      case 2:
	horiz_respacing_n(E,E->slice.tpg,interp_data,1,newX);
	break;
      case 3:
	horiz_respacing_n(E,E->slice.grv,interp_data,1,newX);
	break;
      case 4:
	horiz_respacing_n(E,E->slice.grvt,interp_data,1,newX);
	break;
      case 5:
	horiz_respacing_n(E,E->slice.shflux,interp_data,1,newX);
      }
      add_surf_graph_to_pixmap(E,interp_data,&pixmap,cmap,newX,border,border,
			       vel_plot_height,1,50,50,50,PPM_plot_number);
    }
    sprintf(name_plus,"%s%d",name,PPM_plot_number);

    if(E->control.verbose)
      fprintf(stderr,"Writing ppm file: %s\n",name_plus);

    output_ppm(E,name_plus,pixmap);
 
    /* free up memory (another pixmap may be different size)  */

    free((void *) pixmap.R);
    free((void *) pixmap.G);
    free((void *) pixmap.B);
    free((void *) interp_data);
  }
  return;
}

/* Output ppm file */

void output_ppm(
	   struct All_variables *E,
	   char *name,
	   struct ppm pixmap
	   )
{
  FILE *fp; 
  int i;
  char byte;
  char command[5000];
  
  if((fp=fopen(name,"w")) != NULL) {
    fprintf(fp,"P6\n");  /* Magic # */
    fprintf(fp,"%d\n",pixmap.x);  /* Width */
    fprintf(fp,"%d\n",pixmap.y);  /* Height */
    fprintf(fp,"255\n");  /* Maxval */
    
    
    for(i=1;i<=pixmap.x*pixmap.y;i++) {
      byte=pixmap.R[i];
      putc(byte,fp);
      byte=pixmap.G[i];
      putc(byte,fp);
      byte=pixmap.B[i];
      putc(byte,fp);
    }
    fclose(fp);
  }

#if !defined ( SLAVE_TO_GENETICS ) 

  if(E->control.COMPRESS){
    /*sprintf(command,"%s -f  %s",COMPRESS_BINARY,name);*/
    sprintf(command,"%s -f  %s",E->control.gzip,name); /*DAS: 17-01-03*/
    system(command);
  } 

#else

  fprintf(stderr,"Created image file %s\n",name);

#endif
  return;
}

/* Layers the representation for m*n field onto 
   the existing pixmap. Field is a citcom field
   which has been interpolated to the appropriate
   grid density (but in std coordinate system).
   If the dimensions do not agree
   then the function aborts */

void add_field_to_pixmap(
  struct All_variables *E,
  standard_precision *field,
  struct ppm *pixmap,
  struct colormap cmap,
  int m,
  int n,
  int om,
  int on,
  int edging,
  int eR,
  int eG,
  int eB
)
{
  int i,j;
  int R,G,B;
  int Rs,Gs,Bs;
  int Rh,Gh,Bh;
  int px,py,pp;
  float Op,Oph,Ops;
  
  void supply_color();

  const int mx = pixmap->x;
  const int my = pixmap->y;

  
  for(i=1;i<=m;i++)
    for(j=1;j<=n;j++) {

      /* Where is this pixel in the larger pixmap ? 
	 If it's off the edge, carry on to next pixel */

      px = om + i;
      py = on + j;

      if(px > mx || py > my || px < 0 || py < 0)
	continue;
     
      pp = px + (py-1) * mx;
    
      /* obtain colormap entry */
      supply_color(E,cmap,field[j+(i-1)*n],&R,&G,&B,&Op);

      /* overlay onto existing ppm */
      if(R >= 0)
	pixmap->R[pp] = (int)(Op * R + (1-Op)* pixmap->R[pp]); 
      if(G >= 0)
	pixmap->G[pp] = (int)(Op * G + (1-Op)* pixmap->G[pp]); 
      if(B >= 0)
	pixmap->B[pp] = (int)(Op * B + (1-Op)* pixmap->B[pp]); 

      if(pixmap->R[pp] > 255)
	pixmap->R[pp] = 255;
      if(pixmap->G[pp] > 255)
	pixmap->G[pp] = 255;
      if(pixmap->B[pp] > 255)
	pixmap->B[pp] = 255;

    }

  /* Add border to the pixmap if requested */

  if(edging) {
    for(i=0;i<=m+1;i++) {
      px = om + i;
      py = on;
      if(!(px > mx || py > my || px < 0 || py < 0)) {
	pp = px + (py-1) * mx;
	pixmap->R[pp] = eR;
  	pixmap->G[pp] = eG;
  	pixmap->B[pp] = eB;
      }
      py = on + n + 1;
      if(!(px > mx || py > my || px < 0 || py < 0)) {
	pp = px + (py-1) * mx;
	pixmap->R[pp] = eR;
  	pixmap->G[pp] = eG;
  	pixmap->B[pp] = eB;
      }
    }
    
    for(j=1;j<=n;j++) {
      py = on + j;  
      px = om;
      if(!(px > mx || py > my || px < 0 || py < 0)) {
	pp = px + (py-1) * mx;
	pixmap->R[pp] = eR;
  	pixmap->G[pp] = eG;
  	pixmap->B[pp] = eB;
      }	
      px = om + m + 1;
      if(!(px > mx || py > my || px < 0 || py < 0)) {
	pp = px + (py-1) * mx;
	pixmap->R[pp] = eR;
  	pixmap->G[pp] = eG;
  	pixmap->B[pp] = eB;
      }	
    }
  }
  return;
}

 /* Add border to the pixmap as requested */

void add_border_to_pixmap(
  struct All_variables *E,
  struct ppm *pixmap,
  int m,
  int n,
  int om,
  int on,
  int eR,
  int eG,
  int eB
)
{
  int i,j,k;
  int R,G,B;
  int px,py,pp;
  int resx,resz;
  int n1,n2;
  standard_precision maxX,minX,maxZ,minZ;
  standard_precision X,Z;
  float Op;
  
  void supply_color();

  const int mx = pixmap->x;
  const int my = pixmap->y;

  
  resx = 5 * m / E->mesh.elx;
  resz = 5 * n / E->mesh.elz;


  minX = minZ = 1.0e32;
  maxX = maxZ =-1.0e32;
  for(i=1;i<=E->mesh.nno;i++) {
    if(E->x[1][i] > maxX)
      maxX = E->x[1][i];
    if(E->x[2][i] > maxZ)
      maxZ = E->x[2][i];
    if(E->x[1][i] < minX)
      minX = E->x[1][i];
    if(E->x[2][i] < minZ)
      minZ = E->x[2][i];
  }
  
    for(i=1;i<E->mesh.nox;i++) {
      /* Across the top surface */    

      n1 = 1 + (i-1) * E->mesh.noz;
      n2 = n1 + E->mesh.noz;

      for(k=0;k<=resx;k++) {      
       X = ((E->x[1][n1] + k * (E->x[1][n2]-E->x[1][n1]) / resx) - minX) / (maxX-minX);
       Z = ((E->x[2][n1] + k * (E->x[2][n2]-E->x[2][n1]) / resx) - minZ) / (maxZ-minZ);
       px = om + m * X;
       py = on + n * Z;
       if(!(px > mx || py > my || px < 0 || py < 0)) {
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = eR;
	 pixmap->G[pp] = eG;
	 pixmap->B[pp] = eB;
       }
      }

   /* Across the bottom surface */    

      n1 = E->mesh.noz + (i-1) * E->mesh.noz;
      n2 = n1 + E->mesh.noz;

      for(k=0;k<=resx;k++) {      
       X = ((E->x[1][n1] + k * (E->x[1][n2]-E->x[1][n1]) / resx) - minX) / (maxX-minX);
       Z = ((E->x[2][n1] + k * (E->x[2][n2]-E->x[2][n1]) / resx) - minZ) / (maxZ-minZ);
       px = om + m * X;
       py = 1 + on + n * Z;
       if(!(px > mx || py > my || px < 0 || py < 0)) {
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = eR;
	 pixmap->G[pp] = eG;
	 pixmap->B[pp] = eB;
       }
      }
    }
    
  for(i=1;i<E->mesh.noz;i++) {
      /* Down the left surface */    

      n1 = i;
      n2 = n1 + 1;

      for(k=0;k<=resz;k++) {      
       X = ((E->x[1][n1] + k * (E->x[1][n2]-E->x[1][n1]) / resz)- minX) / (maxX-minX);
       Z = ((E->x[2][n1] + k * (E->x[2][n2]-E->x[2][n1]) / resz)- minZ) / (maxZ-minZ);
       px = om + m * X;
       py = on + n * Z;
       if(!(px > mx || py > my || px < 0 || py < 0)) {
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = eR;
	 pixmap->G[pp] = eG;
	 pixmap->B[pp] = eB;
       }
      }

   /* Down the right surface */    

      n1 = i + (E->mesh.nox-1) * E->mesh.noz;
      n2 = n1 + 1;
 
      for(k=0;k<=resz;k++) {      
       X = ((E->x[1][n1] + k * (E->x[1][n2]-E->x[1][n1]) / resz) - minX) / (maxX-minX);
       Z = ((E->x[2][n1] + k * (E->x[2][n2]-E->x[2][n1]) / resz) - minZ) / (maxZ-minZ);
       px = 1 + om + m * X;
       py = on + n * Z;
       if(!(px > mx || py > my || px < 0 || py < 0)) {
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = eR;
	 pixmap->G[pp] = eG;
	 pixmap->B[pp] = eB;
       }
      }
    }


  return;
}

void add_cross_sections_to_pixmap (
 struct All_variables *E,
 struct ppm *pixmap,
 int border,
 int m,
 int n,
 int om,
 int on,
 int height,
 int width,
 int PPM_plot_number
 )
{  
  int i,j,sample;
  int R,G,B;
  int px,py,pp;
  int dots;
  int element;

  standard_precision xx,z1,z2,zz,x1,x2;
  standard_precision yy; /*RAA: 4/01/02, added this line*/
  standard_precision phi;
  standard_precision *plot_variable;
  standard_precision *Node;
  standard_precision xmin,xmax,zmin,zmax;
  standard_precision ymin,ymax; /*RAA: added this line*/
  standard_precision scale;
  
  const int mx = pixmap->x;
  const int my = pixmap->y;

  standard_precision lN[ELNMAX+1];
  standard_precision eta1,eta2,eta3;
  standard_precision dOmega;


  Node = (standard_precision *)Malloc0((E->mesh.nno+2)*sizeof(standard_precision));


  xmin = min(E->x[1][1],E->x[1][E->mesh.noz]);
  xmax = max(E->x[1][E->mesh.nno-E->mesh.noz+1],E->x[1][E->mesh.nno]);
  zmin = min(E->x[2][1],E->x[2][E->mesh.nno-E->mesh.noz+1]);
  zmax = max(E->x[2][E->mesh.noz],E->x[2][E->mesh.nno]);
  /*RAA: 4 new lines below*/
  if(E->mesh.nsd==3) {
    ymin = min(E->x[3][1],E->x[3][E->mesh.noz]);
    ymax = max(E->x[3][E->mesh.nno-E->mesh.noz+1],E->x[3][E->mesh.nno]);
  }

  /* Run through each of the sampling locations to see if we need to plot anything */


  for(sample=0;sample<E->tracer.SAMPLE_PTS;sample++) {

  if(! (E->tracer.sample_in_plot[sample] == PPM_plot_number))
      continue;
 
    /* Which (interpolated) variable are we going to plot ? */

    switch(E->tracer.sample_type[sample]) {
    case 1:
      plot_variable=E->T;
      break;
    case 2:
      plot_variable=E->V[1];
      break;
    case 3:
      plot_variable=E->V[2];
      break;
    case 4:
      plot_variable=E->nQ;
      break;
    case 5:
      plot_variable=E->edot;
      break;
    case 6:
      plot_variable=E->strd;
      break;
    case 7:
      plot_variable=E->strs;
      break;
    case 8:
      plot_variable=E->V[3];
      break;
    case 9:
      plot_variable=E->V[4];
      break;
    case 10:
      plot_variable=E->V[5];
      break;
    case 11:
      plot_variable=E->V[6];
      break;
    case 12:
      gs_tracers_to_nodes(E,Node,NULL,NULL,NULL,E->tracer.grain_size,E->mesh.levmax,0); 
      plot_variable=Node;
      break;
    case 13:
      gs_tracers_to_nodes(E,Node,NULL,NULL,NULL,E->tracer.sigma_n,E->mesh.levmax,0); 
      plot_variable=Node;
      break;
    case 14:  /*RAA: 24/09/02, C. O'Neill - melting stuff */
      plot_variable=E->depl;  
      break; 
    case 15:
      plot_variable=E->strd1;
      break;
    default:
      fprintf(stderr,
              "This sample type is not supported: %d\n",
              E->tracer.sample_type[sample]);
      exit(1);
    }


   /* Should be able to use the data which was previously computed in the
       processing of the profile output files */

      scale = E->tracer.sample_plotmax[sample] - E->tracer.sample_plotmin[sample];

      if(E->control.print_convergence)
	fprintf(stderr,"Plotting profile %d between %g and %g ... direction %d\n",sample,
		E->tracer.sample_plotmax[sample],E->tracer.sample_plotmin[sample],E->tracer.sample_direction[sample]);

	fprintf(E->fp,"Plotting profile %d between %g and %g\n",sample,
		E->tracer.sample_plotmax[sample],E->tracer.sample_plotmin[sample]);

    /* Which direction */

    if(E->tracer.sample_direction[sample] == 1) {
      /* Obtain a series of sample points */

      zz = E->tracer.sample_z[sample]; 
      if(E->mesh.nsd == 3)
        yy = E->tracer.sample_y[sample]; 
      else
        yy = 0.0;
       /*RAA: 4/01/02, added the 4 lines above*/

      element=1;
      for(i=1;i<=m;i++) {  /* scan across the surface */
	xx = xmin + (xmax - xmin) * (i-1.0) / (m-1);

       /*RAA: 4/01/02, changed 0.0 to yy below*/
	element = general_tracers_element(E,element,xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);
	/* field value at this point */
      
	if(element == -1)  {
          fprintf(stderr,"Problem w/ sample# %d, counter: i= %d, element not found\n",sample,i);
	  continue;
	}

        v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);
	phi=0.0;
	for(j=1;j<=enodes[E->mesh.nsd];j++) {
	  phi += plot_variable[E->ien[element].node[j]] * lN[j];
	}

	/* This tells us where to plot ... */

	if(phi > E->tracer.sample_plotmax[sample] || phi < E->tracer.sample_plotmin[sample])
	  continue;

	px = border + i;
	py = border + height; 
	if(scale != 0.0) 
	  py  -= (int)(height * (phi -  E->tracer.sample_plotmin[sample]) / scale) ;

	pp = px + (py-1) * mx;

	pixmap->R[pp] = (int) 255.0 * E->tracer.sample_Red[sample];
	pixmap->G[pp] = (int) 255.0 * E->tracer.sample_Green[sample];
	pixmap->B[pp] = (int) 255.0 * E->tracer.sample_Blue[sample];
      }

      for(i=1;i<=height;i++) {
	px = border;
	py = border + i;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;

	px = border + m;
	py = border + i;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;
      }
    }

    else if(E->tracer.sample_direction[sample] == 2)   { /* Z */
  /* Obtain a series of sample points */

      xx = E->tracer.sample_x[sample]; 
      if(E->mesh.nsd == 3)
        yy = E->tracer.sample_y[sample]; 
      else
        yy = 0.0;
      /*RAA: 4/01/02, added the 4 lines above*/

      element=1;
      for(i=1;i<=n;i++) {  /* scan across the surface */
	zz = zmin + (zmax - zmin) * (i-1.0) / (n-1);
      
       /*RAA: 4/01/02, changed 0.0 to yy below*/
	element = general_tracers_element(E,element,xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);

	/* field value at this point */
      
	if(element == -1) {
          fprintf(stderr,"Problem w/ sample# %d, counter: i= %d, element not found\n",sample,i);
	  continue;
	}

        v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);
	phi=0.0;
	for(j=1;j<=enodes[E->mesh.nsd];j++) {
	  phi += plot_variable[E->ien[element].node[j]] * lN[j];
	}
	

	/* This tells us where to plot ... */

	if(phi > E->tracer.sample_plotmax[sample] || phi < E->tracer.sample_plotmin[sample])
	  continue;


	py = on + i ;
	px = 2 * border + m;

	if(scale != 0.0)
	  px += (int) (width * (phi -  E->tracer.sample_plotmin[sample]) / scale) ;

	pp = px + (py-1) * mx;

	pixmap->R[pp] = (int) 255.0 * E->tracer.sample_Red[sample];
	pixmap->G[pp] = (int) 255.0 * E->tracer.sample_Green[sample];
	pixmap->B[pp] = (int) 255.0 * E->tracer.sample_Blue[sample];
      }

      for(i=1;i<=width;i++) {
	px = 2 * border + m + i;
	py = 2 * border + height;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;
	
	px = 2 * border + m + i;
	py = 2 * border + height + n;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;
      }
    }

    /*RAA: 7/01/02, added this bit for Y-direction sampling*/
    else if(E->mesh.nsd == 3 && E->tracer.sample_direction[sample] == 3)   { /* Y */
  /* Obtain a series of sample points */

      xx = E->tracer.sample_x[sample]; 
      zz = E->tracer.sample_z[sample]; 

      element=1;
      for(i=1;i<=m;i++) {  /* scan across the surface */
	yy = ymin + (ymax - ymin) * (i-1.0) / (n-1);
      
	element = general_tracers_element(E,element,xx,zz,yy,&eta1,&eta2,&eta3,E->mesh.levmax);

	/* field value at this point */
      
	if(element == -1) {
          fprintf(stderr,"Problem w/ sample# %d, counter: i= %d, element not found\n",sample,i);
	  continue;
        }

        v_shape_fn(E,element,lN,eta1,eta2,eta3,E->mesh.levmax);
	phi=0.0;
	for(j=1;j<=enodes[E->mesh.nsd];j++) {
	  phi += plot_variable[E->ien[element].node[j]] * lN[j];
	}
	

	/* This tells us where to plot ... */ /*RAA: 7/01/02, set up to plot with x-along the top*/

	if(phi > E->tracer.sample_plotmax[sample] || phi < E->tracer.sample_plotmin[sample])
	  continue;

	px = border + i;
	py = border + 1.*height; 
        if(scale != 0.0) 
	  py  -= (int)(height * (phi -  E->tracer.sample_plotmin[sample]) / scale) ;

	pp = px + (py-1) * mx;

	pixmap->R[pp] = (int) 255.0 * E->tracer.sample_Red[sample];
	pixmap->G[pp] = (int) 255.0 * E->tracer.sample_Green[sample];
	pixmap->B[pp] = (int) 255.0 * E->tracer.sample_Blue[sample];
      }

      for(i=1;i<=height;i++) {
	px = border;
	py = border + i;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;

	px = border + m;
	py = border + i;
	pp = px + (py-1) * mx;
	pixmap->R[pp] =	pixmap->G[pp] =	pixmap->B[pp] = 0;
      }
    }


    /* And now plot the sample point itself */

    for(i=-2;i<=2;i++)
      for(j=-2;j<=2;j++) {
	px = border + m * (E->tracer.sample_x[sample] - xmin) / (xmax - xmin);
	py = 2 * border + height + n * (E->tracer.sample_z[sample] - zmin) / (zmax - zmin);
	pp = (px+i) + ((py+j)-1) * mx;

	
	pixmap->R[pp] = (int) 255.0 * E->tracer.sample_Red[sample];
	pixmap->G[pp] = (int) 255.0 * E->tracer.sample_Green[sample];
	pixmap->B[pp] = (int) 255.0 * E->tracer.sample_Blue[sample];
      }  
  }

  free((void *) Node);

  return;
}

void add_surf_graph_to_pixmap(
			      struct All_variables *E,
			      standard_precision *field,
			      struct ppm *pixmap,
			      struct colormap cmap,
			      int m,
			      int om,
			      int on,
			      int height,
			      int axes,
			      int aR,
			      int aG,
			      int aB,
			      int PPM_plot_number
			      )
{
  int i,j;
  int R,G,B;
  int px,py,pp;
  int dots;
  float Op;
  float maxf,minf,scale,scale10,scale11;
  float f1;

 
  void supply_color();

  const int mx = pixmap->x;
  const int my = pixmap->y;
  const int hheight = height / 2;

  /* First, recover the statistics for the data field */

 if(1 || E->control.PPM_surf_auto[PPM_plot_number]) {
  maxf=-1e32;
  minf= 1e32;
  for(i=1;i<=m;i++) {
    if(field[i] > maxf)
      maxf = field[i] ;
    if(field[i] < minf)
      minf =  field[i];
  }
   scale = max(fabs(maxf),fabs(minf));
 }
 /* else */
 if (!E->control.PPM_surf_auto[PPM_plot_number]) {
   scale = max(fabs(E->control.PPM_surf_max[PPM_plot_number]),
	       fabs(E->control.PPM_surf_min[PPM_plot_number]));
 }
  
 if(scale == 0.0) 
    return;
  
 scale10 = pow(10.0,floor(log10(scale*1.001)));
 scale11 = pow(10.0, ceil(log10(scale*0.999)));

 if (E->control.PPM_surf_auto[PPM_plot_number]) {
   scale=scale11;
 }


  /* We'll plot to nearest power of ten */

    for(i=1;i<=m;i++) {
     	f1=field[i];

	if(! E->control.PPM_surf_auto[PPM_plot_number]) {
	 
	  if(f1 < -scale)
	    f1 = -scale;
	  else
	    if(f1 > scale)
	      f1 = scale;
	}
           
      px = om + i;
      py = on + hheight - (int)(hheight * f1 / scale) ;
      
      if(py >= on + hheight) 
	for(j=on + hheight;j<=py;j++) {
	  pp = px + (j-1) * mx;
	  pixmap->B[pp] = 150;
	  pixmap->G[pp] = pixmap->R[pp] = 0;
	}
      else
	for(j=py;j<=on + hheight;j++) {
	  pp = px + (j-1) * mx; 
	  pixmap->R[pp] = 150;
	  pixmap->G[pp] = pixmap->B[pp] = 0;
	}
    }


  /* Fill in pixels to generate the axes */
  if(axes) {
 
  /* X axis */ 
  for(i=1;i<=m;i++) { 
      px = om + i;
      py = on + height / 2;
      pp = px + (py-1) * mx;
      pixmap->R[pp] = aR;
      pixmap->G[pp] = aG;
      pixmap->B[pp] = aB;
      /* py = on + height;
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = aR;
	 pixmap->G[pp] = aG;
	 pixmap->B[pp] = aB;
	 py = on;
	 pp = px + (py-1) * mx;
	 pixmap->R[pp] = aR;
	 pixmap->G[pp] = aG;
	 pixmap->B[pp] = aB;  */
  } 

  /* Y1/Y2 axis */
  for(i=1;i<=height;i++) { 
    px = om;
    py = on + i;
    pp = px + (py-1) * mx;
    pixmap->R[pp] = aR;
    pixmap->G[pp] = aG;
    pixmap->B[pp] = aB;
  
    px = om + m;
    pp = px + (py-1) * mx;
    pixmap->R[pp] = aR;
    pixmap->G[pp] = aG;
    pixmap->B[pp] = aB;
  }

  /* Add dots to indicate Y scaling */
  
  scale10 = pow(10.0,floor(log10(max(fabs(maxf),fabs(minf)))));
  dots = (int) log10(scale10);


  /* dots = nint(log10(scale10)); */
  
  for(i=0;i<5;i++) {
    px = om - i;
    py = on + hheight - (int)(hheight * scale10 / scale)  ;
    pp = px + (py-1) * mx;
    pixmap->R[pp] = aR;
    pixmap->G[pp] = aG;
    pixmap->B[pp] = aB;
    py = on + hheight + (int)(hheight * scale10 / scale) ;
    pp = px + (py-1) * mx;
    pixmap->R[pp] = aR;
    pixmap->G[pp] = aG;
    pixmap->B[pp] = aB;
  }
    

  if (dots > 0) 
    for(i=1;i<=dots;i++) {
      for(j=0;j<=5;j++) {
	px = om + i * 3;
	py = on - j;
	pp = px + (py-1) * mx;
	pixmap->R[pp] = 255;
   	pixmap->G[pp] = pixmap->B[pp] = 0;
      }
    }
  else
    for(i=dots;i<=-1;i++) {
     for(j=0;j<=5;j++) {
       px = om - i * 3;
       py = on + height + j;
       pp = px + (py-1) * mx;
       pixmap->B[pp] = 255;
       pixmap->G[pp] = pixmap->R[pp] = 0;
     }
    }
  }

  return;
}

void add_tracers_to_pixmap(
  struct All_variables *E,
  struct TRACERS tracers,
  struct ppm *pixmap,
  int m,
  int n,
  int om,
  int on,
  int colorbar_width,
  int PPM_plot_number,
  int boxnum,
  int slicenums
)
{
  /*RAA: 7/1/02 added boxnum & slicenums in function call to get more than 1 box per ppm file*/
  int i,j,k,l;
  int R,G,B;
  int Rc,Gc,Bc;
  int Rh,Gh,Bh;
  int Rs,Gs,Bs;
  int px,py,pp;
  int tx,ty;
  int bg,count;
  int are_tracers_shaded_for_T;
  float Op,Opc,Oph,Ops;
  standard_precision scale,CC;

  standard_precision *Colorfield,*WorkTr;
  int *Colormask,*WorkMask;

  void supply_color();

  standard_precision maxT,minT,T,delT_1;
  standard_precision maxEdot,minEdot,aveEdot,refEdot;
  standard_precision a1;
  standard_precision middle_y;  /*RAA: added this variable*/
  const standard_precision delt=E->control.PPM_3D_delt;  /*RAA: added this*/

  static standard_precision running_max_Edot = 0.0;
  static standard_precision running_ave_Edot = 0.0;


  struct colormap cmap;

  const int mx = pixmap->x;
  const int my = pixmap->y;

  WorkTr = (standard_precision *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(standard_precision));
  WorkMask = (int *) Malloc0((E->tracer.NUM_TRACERS+1) * sizeof(int));


  /*RAA: define these new variables*/
  if(2==E->mesh.nsd) {
    /*delt = 0.0;*/
    middle_y = 0.0;
  } 
  else {
    middle_y = (E->x[3][E->mesh.nno] - E->x[3][1])/2.0;
    /*delt = 0.5;*/
  }

  if(boxnum==1) 
    fprintf(stderr,"middle_y and delt for ppms: %g %g\n",middle_y,delt);


  /* Default Mask is NULL */

  Colormask = NULL;

  switch(E->control.PPM_coloring[PPM_plot_number]) {
  case 1:
    Colorfield=E->tracer.T;
    break;
  case 2:
    Colorfield=E->tracer.weight;
    break;
  case 3:
    Colorfield=E->tracer.strd;
    break;
  case 4:  
    Colorfield=E->tracer.DQ1;
    break;
  case 5:
    /* Colorfield=E->tracer.grain_size; */
    for(i=1;i<=tracers.NUM_TRACERS;i++) {
      if(E->tracer.property_group[i] < 0)
	continue;
      if(E->tracer.grain_size[i] >= 0.0)
	WorkTr[i] = log10(E->tracer.grain_size[i]);
      else
	WorkTr[i] = 0.0;
    }
    Colorfield=WorkTr; 

    break;
  case 6:
    Colorfield=E->tracer.DQ;
    break;
  case 7:
    for(i=1;i<=tracers.NUM_TRACERS;i++) {
      if(E->tracer.edot[i] != 0.0)
	WorkTr[i] = (E->tracer.edot[i]);
      else
	WorkTr[i] = 0.0;
    }
    Colorfield=WorkTr; 
    break;
  case 8:
    Colorfield=E->tracer.edotp_integrated;
    break;
  case 9:
    /*RAA: 12/12/02, a few additions for 3D elasticity*/
    if(E->control.ELASTICITY && E->mesh.nsd==2) {
      for(i=1;i<=tracers.NUM_TRACERS;i++) {
	WorkTr[i] = sqrt( 0.5 * (E->tracer.S11[i]*E->tracer.S11[i] + 
				 E->tracer.S12[i]*E->tracer.S12[i] * 2.0 + 
				 E->tracer.S22[i]*E->tracer.S22[i] ) );
      }
      
    }
    else if(E->control.ELASTICITY && E->mesh.nsd==3) {
      for(i=1;i<=tracers.NUM_TRACERS;i++) {
	WorkTr[i] = sqrt( 0.5 * (E->tracer.S11[i]*E->tracer.S11[i] + 
				 E->tracer.S12[i]*E->tracer.S12[i] * 2.0 + 
				 E->tracer.S13[i]*E->tracer.S13[i] * 2.0 + 
				 E->tracer.S23[i]*E->tracer.S23[i] * 2.0 + 
	                         E->tracer.S33[i]*E->tracer.S33[i] + 
				 E->tracer.S22[i]*E->tracer.S22[i] ) );
      }
      
    }
    else
      for(i=1;i<=tracers.NUM_TRACERS;i++) {
	WorkTr[i] = 0.0;
      }

    Colorfield=WorkTr;
    break;

  case 10:
    for(i=1;i<=tracers.NUM_TRACERS;i++) 
      WorkTr[i] =  1.570796;

      if(E->control.ORTHOTROPY) {
	for(i=1;i<=tracers.NUM_TRACERS;i++) {
	  if(E->tracer.visc[E->tracer.property_group[i]].Ortho_viscosity_ratio < 1.0 ||
	     (E->control.CHEM_TRANS && 
	      E->tracer.visc[E->tracer.property_group[i]].mobile_phase_ratio != 1.0)
	     ) {
	    WorkTr[i] = acos(E->tracer.n1[i]);
	    WorkMask[i] = 1;
	  }
	  else
	    WorkMask[i] = 0;
	}
      }

    Colorfield=WorkTr;
    Colormask = WorkMask;
    break;

  case 11:
    for(i=1;i<=tracers.NUM_TRACERS;i++) 
      WorkTr[i] =  0.0;

    if(E->control.CHEM_TRANS) {
	for(i=1;i<=tracers.NUM_TRACERS;i++) {
	  if(E->tracer.visc[E->tracer.property_group[i]].mobile_phase_ratio != 1.0) {
	    WorkTr[i] = E->tracer.sigma_n[i];
	    WorkMask[i] = 1;
	  }
	  else
	    WorkMask[i] = 0;
	}
      }
      
    Colorfield=WorkTr;
    Colormask=WorkMask;
    break;

  case 12:

    for(i=1;i<=tracers.NUM_TRACERS;i++) 
      WorkTr[i] =  0.0;

      if(E->control.CHEM_TRANS) {
	for(i=1;i<=tracers.NUM_TRACERS;i++) {
	  if(E->tracer.visc[E->tracer.property_group[i]].mobile_phase_ratio != 1.0) {
	    WorkTr[i] = E->tracer.volfraction[i];
	    WorkMask[i] = 1;
	  }
	  else
	    WorkMask[i] = 0;
	}
      }

    Colorfield=WorkTr;
    Colormask=WorkMask;
    break;

  case 13:  /*RAA: 24/09/02, C. O'Neill - melting stuff */
    Colorfield=E->tracer.depl; 
    break; 
  default:
    fprintf(stderr,
            "Coloring not implemented for this PPM_coloring choice: %d\n",
            E->control.PPM_coloring[PPM_plot_number]);
    exit(1);
	    
  }
  
  /* Choose colors using specified ranges or automatically
     color using max/min values */


  if(1 || E->control.PPM_color_auto[PPM_plot_number]) {
    maxT=-1.0e32;
    minT=1.0e32;
    k=-1;
    l=-1;
   
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {
      if(Colormask != NULL &&  Colormask[i] == 0)
	continue;
      
      if(maxT < Colorfield[i]) {
	maxT = Colorfield[i];
	k=i;
      }
      if(minT > Colorfield[i]) {
	minT = Colorfield[i];
	l=i;
      }
    }
  }

  fprintf(E->fp,"PPM: Plotting field %d between %g and %g (%d @ %g,%g, %d @ %g,%g) \n",
	  PPM_plot_number,maxT,minT,k,E->tracer.tx[k],E->tracer.tz[k],l,E->tracer.tx[l],E->tracer.tz[l]);

  if(! E->control.PPM_color_auto[PPM_plot_number]) {
    maxT = E->control.PPM_color_max[PPM_plot_number];
    minT = E->control.PPM_color_min[PPM_plot_number];
  }

  if(maxT-minT != 0.0)
    delT_1 = 1.0 / (maxT-minT);
  else
    delT_1 = 1.0;

  /* If yielding, want to colour according to prevailing plastic
     strain rates (which requires knowing typical values) */


  if(E->viscosity.YIELD) {
    maxEdot=-1.0e32;
    minEdot=1.0e32;
  
    count=0;
    aveEdot=0.0;
   
    for(i=1;i<=E->tracer.NUM_TRACERS;i++) {    
      if(E->tracer.yielded[i] == 0.0)
	continue;
      
      count++;
      aveEdot +=E->tracer.edot[i];

      if(maxEdot < E->tracer.edot[i]) {
	maxEdot = E->tracer.edot[i];
      }
      if(minEdot > E->tracer.edot[i]) {
	minEdot = E->tracer.edot[i];
      }
    }
 
    if(count!=0)
      aveEdot /= (standard_precision) count;

    running_max_Edot = 0.333 * maxEdot + 0.667 * running_max_Edot;
    running_ave_Edot = 0.1 * aveEdot + 0.9 * running_ave_Edot;

    if(running_max_Edot <= 0.0)
      running_max_Edot = 1.0; /* for example during printout of initial condition */
    if(running_ave_Edot <= 0.0)
      running_ave_Edot = 1.0; /* for example during printout of initial condition */
 
  }

  refEdot = 0.667 * running_max_Edot + 0.333 * running_ave_Edot;

  if(E->control.print_convergence)
    fprintf(stderr,"Plastic strain rate: Max = %g (%g), Min = %g, Ave/%d = %g (%g) -> %g\n",
	    maxEdot,running_max_Edot,minEdot,count,aveEdot,running_ave_Edot,refEdot);
  



  for(i=1;i<=tracers.NUM_TRACERS;i++)  {

   
    /* Compute location in the pixmap ... SHOULD USE MAXX,MINX ETC*/

    /*RAA: 7/01/01, limit the ppm to slices in y for 3D, added 'if', etc below*/
    if(2==E->mesh.nsd && boxnum==1) {
      tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
      ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
    }
    else if (3==E->mesh.nsd && boxnum==1 && slicenums==1) {
      if(E->tracer.ty[i] > (middle_y - delt) && E->tracer.ty[i] < (middle_y + delt)) { 
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;  
    }
    else if (3==E->mesh.nsd && boxnum==1 && slicenums==2) {
     if(E->tracer.ty[i] > 0.0 && E->tracer.ty[i] < delt) { 
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;
    }
    else if (3==E->mesh.nsd && boxnum==2 && slicenums==2) {
     if(E->tracer.ty[i] > (E->x[3][E->mesh.nno] - delt) && 
        E->tracer.ty[i] < E->x[3][E->mesh.nno]) { 
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;
    }
    else if (3==E->mesh.nsd && boxnum==1 && slicenums==3) {
      if(E->tracer.ty[i] > (middle_y/2. - delt) && E->tracer.ty[i] < (middle_y/2. + delt)) {
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;
    }
    else if (3==E->mesh.nsd && boxnum==2 && slicenums==3) {
      if(E->tracer.ty[i] > (middle_y - delt) && E->tracer.ty[i] < (middle_y + delt)) { 
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;
    }
    else if (3==E->mesh.nsd && boxnum==3 && slicenums==3) {
      if(E->tracer.ty[i] > ((middle_y + middle_y/2.) - delt) && 
         E->tracer.ty[i] < ((middle_y + middle_y/2.) + delt)) { 
        tx = 1 + m *  (tracers.tx[i] - E->x[1][1]) / (E->x[1][E->mesh.nno] - E->x[1][1]);
        ty = 1 + n *  (tracers.tz[i] - E->x[2][1]) / (E->x[2][E->mesh.nno] - E->x[2][1]);
      }
      else
        continue;
    }

    /* Where is this pixel in the larger pixmap ? 
       If it's off the edge, carry on to next pixel */

    px = om + tx;
    py = on + ty;

      if(px > mx || py > my || px < 0 || py < 0)
	continue;
     
      pp = px + (py-1) * mx;
    
      /* obtain colormap entry */
      
      if(E->tracer.property_group[i] < 0) {
	R  = 0.70;
	G  = 0.85;
	B  = 1.00;
	Op = 0.10;
      }
      else {

      Opc = E->tracer.Opacity[E->tracer.property_group[i]][PPM_plot_number];
      Rc = (int) 255.0 * E->tracer.Red[E->tracer.property_group[i]][PPM_plot_number];
      Gc = (int) 255.0 * E->tracer.Green[E->tracer.property_group[i]][PPM_plot_number];
      Bc = (int) 255.0 * E->tracer.Blue[E->tracer.property_group[i]][PPM_plot_number];

      Oph = E->tracer.Opacity_hot[E->tracer.property_group[i]][PPM_plot_number];
      Rh = (int) 255.0 * E->tracer.Red_hot[E->tracer.property_group[i]][PPM_plot_number];
      Gh = (int) 255.0 * E->tracer.Green_hot[E->tracer.property_group[i]][PPM_plot_number];
      Bh = (int) 255.0 * E->tracer.Blue_hot[E->tracer.property_group[i]][PPM_plot_number];

      Ops = E->tracer.Opacity_strained[E->tracer.property_group[i]][PPM_plot_number];
      Rs = (int) 255.0 * E->tracer.Red_strained[E->tracer.property_group[i]][PPM_plot_number];
      Gs = (int) 255.0 * E->tracer.Green_strained[E->tracer.property_group[i]][PPM_plot_number];
      Bs = (int) 255.0 * E->tracer.Blue_strained[E->tracer.property_group[i]][PPM_plot_number];
   
      /* Change shade according to T */
      
      CC = min(maxT,Colorfield[i]);
      CC = max(minT,CC);
      
      if(Oph >= 0.0) {
	R = (int) ((CC-minT) * delT_1 * Rh  + (maxT-CC) * delT_1 * Rc); 
	G = (int) ((CC-minT) * delT_1 * Gh  + (maxT-CC) * delT_1 * Gc); 
	B = (int) ((CC-minT) * delT_1 * Bh  + (maxT-CC) * delT_1 * Bc);
	Op =      ((CC-minT) * delT_1 * Oph + (maxT-CC) * delT_1 * Opc);
      }
      else { /* No value was set */
	R = Rc;
	G = Gc;
	B = Bc;
	Op = Opc;
      }
      
      /* Now change if yielding occurs */

      if(E->viscosity.YIELD && Ops >= 0.0) {           

	a1 = 0.0;

	if(E->tracer.yielded[i] != 0.0  &&        /* Currently failing, color by plastic strain */
	   (E->tracer.visc[E->tracer.property_group[i]].yield_stress_E0 != 0.0)) {

	  /* fprintf(stderr,"1-Tracer %d plots in a yielding state\n",i); */


	  if(E->tracer.edotp_integrated[i] > 0.0) {
	    /* Color according to actual strength change */

	    if(E->control.PPM_strain != 0.0) {
	      a1 =  E->control.PPM_strain * min(1.0,pow(E->tracer.edotp_integrated[i]/
				     E->tracer.visc[E->tracer.property_group[i]].yield_stress_E0,
				     E->tracer.visc[E->tracer.property_group[i]].yield_stress_En));
	      	      
	      if(E->tracer.visc[E->tracer.property_group[i]].yield_stress_Edot0 != 0.0) {
		a1 += (1.0-E->control.PPM_strain) * min(1.0,E->tracer.edot[i]/
							E->tracer.visc[E->tracer.property_group[i]].yield_stress_Edot0);
	      }
	      else {
		a1 += (1.0-E->control.PPM_strain) /* E->tracer.yielded[i]*/ * min(1.0,E->tracer.edot[i]/refEdot);
	      }

	    }
	    /* Or according to current strain rate */
	    else {
	      if(E->tracer.visc[E->tracer.property_group[i]].yield_stress_Edot0 != 0.0) {
		a1 = min(1.0,E->tracer.edot[i]/E->tracer.visc[E->tracer.property_group[i]].yield_stress_Edot0);
	      }
	      else {
		a1 = E->tracer.yielded[i] /* min(1.0,E->tracer.edot[i]/running_max_Edot)*/;
	      }
	    }
	  }
	    
	  else
	    if(E->monitor.solution_cycles == 1)
	      a1 = 0.3;
	}

	else if (E->tracer.yielded[i] != 0.0) {  /* Currently failing but no tracking of plastic strain */
	  a1 = 0.1;
	  /* fprintf(stderr,"3-Tracer %d plots in a yielding state\n",i); */

	}

	if((E->tracer.yielded[i]==0.0) &&     /* Has accumulated plastic strain in earlier timesteps */
	   (E->tracer.edotp_integrated[i] > 0.0) && 
	   (E->tracer.visc[E->tracer.property_group[i]].yield_stress_E0 != 0.0)) {

	  if(E->control.PPM_strain != 0.0) {
	    a1 = 0.2 * E->control.PPM_strain * min(1.0,E->tracer.edotp_integrated[i]/
		     fabs(E->tracer.visc[E->tracer.property_group[i]].yield_stress_E0));
	  }

	}

	R = (int) (a1 * Rs  + (1-a1) * R); 
	G = (int) (a1 * Gs  + (1-a1) * G); 
	B = (int) (a1 * Bs  + (1-a1) * B); 
	Op =      (a1 * Ops + (1-a1) * Op); 

      }
 
      /* Visualize any phase boundaries */
 
	if(E->tracer.time_since_phase_change[i] <=  E->advection.timestep) {
	  R = 255;
	  G = 255;
	  B = 0;
	}

	if(E->tracer.Phases[E->tracer.property_group[i]]!=1) {
	
	  scale = 0.5 * E->tracer.Current_phase[i];

	  R -= (int) (scale * R); 
	  G -= (int) (scale * G); 
	  B -= (int) (scale * B); 
	}

	if( E->tracer.Current_phase[i] != E->tracer.Equilibrium_phase[i]) {
	   R += (int) ((255-R) * 0.5); 
	}
      }

	/* overlay onto existing ppm */
	
	  pixmap->R[pp] = (int)(Op * R + (1-Op)* pixmap->R[pp]); 
	  pixmap->G[pp] = (int)(Op * G + (1-Op)* pixmap->G[pp]); 
	  pixmap->B[pp] = (int)(Op * B + (1-Op)* pixmap->B[pp]); 
	
	/* Slight blurring */

	if(px > om) {
	    pixmap->R[pp-mx] = (int)(Op * 0.66 * R + (1-Op * 0.66)* pixmap->R[pp-mx]); 
	    pixmap->G[pp-mx] = (int)(Op * 0.66 * G + (1-Op * 0.66)* pixmap->G[pp-mx]); 
	    pixmap->B[pp-mx] = (int)(Op * 0.66 * B + (1-Op * 0.66)* pixmap->B[pp-mx]); 
      }

      if(px < om + m) {
	pixmap->R[pp+mx] = (int)(Op * 0.66 * R + (1-Op * 0.66)* pixmap->R[pp+mx]); 
	pixmap->G[pp+mx] = (int)(Op * 0.66 * G + (1-Op * 0.66)* pixmap->G[pp+mx]); 
	pixmap->B[pp+mx] = (int)(Op * 0.66 * B + (1-Op * 0.66)* pixmap->B[pp+mx]); 
      }

      if(py > on) {	
	  pixmap->R[pp-1] = (int)(Op * 0.66 * R + (1-Op * 0.66)* pixmap->R[pp-1]); 
	  pixmap->G[pp-1] = (int)(Op * 0.66 * G + (1-Op * 0.66)* pixmap->G[pp-1]); 
	  pixmap->B[pp-1] = (int)(Op * 0.66 * B + (1-Op * 0.66)* pixmap->B[pp-1]); 
      }

      if(py < on + n) {
	pixmap->R[pp+1] = (int)(Op * 0.66 * R + (1-Op * 0.66)* pixmap->R[pp+1]); 
	pixmap->G[pp+1] = (int)(Op * 0.66 * G + (1-Op * 0.66)* pixmap->G[pp+1]); 
	pixmap->B[pp+1] = (int)(Op * 0.66 * B + (1-Op * 0.66)* pixmap->B[pp+1]); 
      }

      pixmap->R[pp-mx-1] = (int)(Op * 0.4 * R + (1-Op * 0.4)* pixmap->R[pp-mx-1]); 
      pixmap->G[pp-mx-1] = (int)(Op * 0.4 * G + (1-Op * 0.4)* pixmap->G[pp-mx-1]); 
      pixmap->B[pp-mx-1] = (int)(Op * 0.4 * B + (1-Op * 0.4)* pixmap->B[pp-mx-1]); 

      pixmap->R[pp+mx-1] = (int)(Op * 0.4 * R + (1-Op * 0.4)* pixmap->R[pp+mx-1]); 
      pixmap->G[pp+mx-1] = (int)(Op * 0.4 * G + (1-Op * 0.4)* pixmap->G[pp+mx-1]); 
      pixmap->B[pp+mx-1] = (int)(Op * 0.4 * B + (1-Op * 0.4)* pixmap->B[pp+mx-1]); 

      pixmap->R[pp-mx+1] = (int)(Op * 0.4 * R + (1-Op * 0.4)* pixmap->R[pp-mx+1]); 
      pixmap->G[pp-mx+1] = (int)(Op * 0.4 * G + (1-Op * 0.4)* pixmap->G[pp-mx+1]); 
      pixmap->B[pp-mx+1] = (int)(Op * 0.4 * B + (1-Op * 0.4)* pixmap->B[pp-mx+1]); 

      pixmap->R[pp+mx+1] = (int)(Op * 0.4 * R + (1-Op * 0.4)* pixmap->R[pp+mx+1]); 
      pixmap->G[pp+mx+1] = (int)(Op * 0.4 * G + (1-Op * 0.4)* pixmap->G[pp+mx+1]); 
      pixmap->B[pp+mx+1] = (int)(Op * 0.4 * B + (1-Op * 0.4)* pixmap->B[pp+mx+1]); 
  }

  
  for(i=1;i<=m;i++)
    for(j=1;j<=n;j++) {
      px = om + i;
      py = on + j;
      pp = px + (py-1) * mx;
       
      if(pixmap->R[pp] > 255)
	pixmap->R[pp] = 255;
      if(pixmap->G[pp] > 255)
	pixmap->G[pp] = 255;
      if(pixmap->B[pp] > 255)
	pixmap->B[pp] = 255;

      if(pixmap->R[pp] < 0)
	pixmap->R[pp] = 0;
      if(pixmap->G[pp] < 0)
	pixmap->G[pp] = 0;
      if(pixmap->B[pp] < 0)
	pixmap->B[pp] = 0;
    
    }

  free((void *) WorkTr);
  free((void *) WorkMask);
 
  return;
}

/* maps the numerical value onto RGB for the
   colormap which is supplied. */

void supply_color(
  struct All_variables *E,
  struct colormap cmap,
  standard_precision value,
  int *R,
  int *G,
  int *B,
  float *Op
)
{

  int i;
  int sect;
  float frac1,frac2;

  sect=-1;
  for(i=0;i<cmap.sections;i++) {
    if(value >= cmap.low[i] && value < cmap.high[i]) {
      sect = i;
      break;
    }
  }
  
  if(-1 == sect) {  /* not found -> set to bg color */
    *R=cmap.background_R;
    *G=cmap.background_G;
    *B=cmap.background_B;

    /* fprintf(stderr,"Color not found ... value = %g\n",value); */
  }
  else {
    frac1 = (cmap.high[sect] - value)/(cmap.high[sect] - cmap.low[sect]);
    frac2 = (value -  cmap.low[sect])/(cmap.high[sect] - cmap.low[sect]);

    *R = (int)( cmap.low_R[sect] * frac1 + cmap.high_R[sect] * frac2); 
    *G = (int)( cmap.low_G[sect] * frac1 + cmap.high_G[sect] * frac2); 
    *B = (int)( cmap.low_B[sect] * frac1 + cmap.high_B[sect] * frac2); 
    *Op= (float)( cmap.opacity_low[sect] * frac1 + cmap.opacity_high[sect] * frac2);
  }

 return;
}


/* Sorts out the appropriate default colormap for a 
   named data field. If not found, returns a 
   default colormap and -1 */

int which_colormap(
   struct All_variables *E,
   char *fieldname,
   struct colormap *cmap,
   standard_precision *field,
   int m,
   standard_precision *MX,
   standard_precision *MN,
   standard_precision *AV
)
{
  const int default_maps = 4;
  struct colormap defmap[4];
  
  int i,found;
  float maxval,minval,aveval;
  float maxneg,minpos;

 /* max/min/average */

  maxval = -1e32;
  maxneg = -1e32;
  minval =  1e32;
  minpos =  1e32;
  aveval =  0.0;


  for(i=1;i<=m;i++) {
    if(maxval < field[i])
      maxval = field[i];
    if(minval > field[i])
      minval = field[i];
    if(field[i] > 0 && minpos > field[i])
      minpos = field[i];
   if(field[i] < 0 && maxneg < field[i])
      maxneg = field[i];
   aveval  += field[i];
  }

  aveval /= m;

  *MN=minval;
  *MX=maxval;
  *AV=aveval;

  /* We'll do this by a tedious process of typing
     the same thing again and again as it is somewhat easier to
     read than a lookup table (!) */
  
    if(strcmp("Temp",fieldname) == 0) {  /* Temperature */
      found=1;

      cmap->sections=5;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      /* section 0 */
      cmap->low[0]=-0.05;
      cmap->opacity_low[0]=1.0;
      cmap->low_R[0]=0;
      cmap->low_G[0]=50;
      cmap->low_B[0]=200;
      cmap->high[0]=0.5 * aveval;
      cmap->opacity_high[0]=1.0;
      cmap->high_R[0]=0;
      cmap->high_G[0]=90;
      cmap->high_B[0]=255;
      /* section 1 */
      cmap->low[1]=0.5 * aveval;
      cmap->opacity_low[1]=1.0;
      cmap->low_R[1]=0;
      cmap->low_G[1]=90;
      cmap->low_B[1]=255;
      cmap->high[1]=0.9 * aveval;
      cmap->opacity_high[1]=0.9;
      cmap->high_R[1]=90;
      cmap->high_G[1]=120;
      cmap->high_B[1]=140;

      /* section 2 */
      cmap->low[2]=0.9 * aveval;
      cmap->opacity_low[2]=0.9;
      cmap->low_R[2]=90;
      cmap->low_G[2]=120;
      cmap->low_B[2]=140;
      cmap->high[2]=aveval + 0.1 * (maxval - aveval);
      cmap->opacity_high[2]=0.9;
      cmap->high_R[2]=50;
      cmap->high_G[2]=180;
      cmap->high_B[2]=120;

      /* section 3 */
      cmap->low[3]=aveval + 0.1 * (maxval - aveval);
      cmap->opacity_low[3]=0.9;
      cmap->low_R[3]=50;
      cmap->low_G[3]=180;
      cmap->low_B[3]=120;
      cmap->high[3]=aveval + 0.5 * (maxval - aveval);
      cmap->opacity_high[3]=0.9;
      cmap->high_R[3]=250;
      cmap->high_G[3]=140;
      cmap->high_B[3]=50;

      /* section 4 */
      cmap->low[4]=aveval + 0.5 * (maxval - aveval);
      cmap->opacity_low[4]=0.9;
      cmap->low_R[4]=250;
      cmap->low_G[4]=140;
      cmap->low_B[4]=50;
      cmap->high[4]=maxval;
      cmap->opacity_high[4]=1.0;
      cmap->high_R[4]=150;
      cmap->high_G[4]=10;
      cmap->high_B[4]=0;

      /* Need to rescale this to temperature scale
	 and overrun will then be anything above the
	 reference T */


       /* section 5 - overrun*/
      cmap->low[5]=1.01;
      cmap->opacity_low[5]=1.0;
      cmap->low_R[5]=150;
      cmap->low_G[5]=20;
      cmap->low_B[5]=0;
      cmap->high[5]=2.0;
      cmap->opacity_high[5]=1.0;
      cmap->high_R[5]=20;
      cmap->high_G[5]=0;
      cmap->high_B[5]=0;
  }
    else if(strcmp("Visc",fieldname) == 0) {  /* Viscosity (Different If Yielding) */
      found=1;
  
      cmap->sections=1;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      
      /* Section 0 */
      cmap->low[0]=minval;  /* Min(-Maxval*0.01,-Minpos*10);      */         
      cmap->opacity_low[0]=0.0;
      cmap->low_R[0]=230;
      cmap->low_G[0]=20;
      cmap->low_B[0]=190;
      cmap->high[0]=0.0; /* -Minpos; */
      cmap->opacity_high[0]=1.0;
      cmap->high_R[0]=230;
      cmap->high_G[0]=30;
      cmap->high_B[0]=195;
      /* Section 1 */
      cmap->low[1]=-minpos;               
      cmap->opacity_low[1]=0.5;
      cmap->low_R[1]=230;
      cmap->low_G[1]=30;
      cmap->low_B[1]=195;
      cmap->high[1]=maxneg;
      cmap->opacity_high[1]=0.7;
      cmap->high_R[1]=255;
      cmap->high_G[1]=0;
      cmap->high_B[1]=205;
      /* Section 2 */
      cmap->low[2]=maxval * 0.01;               
      cmap->opacity_low[2]=0.0;
      cmap->low_R[2]=120;
      cmap->low_G[2]=120;
      cmap->low_B[2]=120;
      cmap->high[2]=maxval;
      cmap->opacity_high[2]=0.5;
      cmap->high_R[2]=0;
      cmap->high_G[2]=0;
      cmap->high_B[2]=0;

    }
      
    else if(strcmp("Edtp",fieldname) == 0) {  /* Viscosity (different if yielding) */
      found=1;
  
      cmap->sections=1;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      
      /* section 0 */
      cmap->low[0]=minval + 0.01 * (maxval-minval);         
      cmap->opacity_low[0]=0.0;
      cmap->low_R[0]=230;
      cmap->low_G[0]=20;
      cmap->low_B[0]=200;
      cmap->high[0]=maxval; 
      cmap->opacity_high[0]=0.9;
      cmap->high_R[0]=255;
      cmap->high_G[0]=50;
      cmap->high_B[0]=200;
     
       }
    else if(strcmp("Chem",fieldname) == 0) {  /* Chemistry  */
      found=1;

      cmap->sections=1;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      /* section 0 */
      cmap->low[0]=0.5 * maxval + 0.5 * aveval;
      cmap->opacity_low[0]=0.5;
      cmap->low_R[0]=240;
      cmap->low_G[0]=200;
      cmap->low_B[0]=150;
      cmap->high[0]=maxval * 1.01;
      cmap->opacity_high[0]=0.8;
      cmap->high_R[0]=255;
      cmap->high_G[0]=200;
      cmap->high_B[0]=0;
    }
    else if(strcmp("H_2O",fieldname) == 0) {  /* Chemistry 2  */
      found=1;

      cmap->sections=1;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      /* section 0 */
      cmap->low[0]=0.5 * maxval + 0.5 * aveval;
      cmap->opacity_low[0]=0.5;
      cmap->low_R[0]=240;
      cmap->low_G[0]=150;
      cmap->low_B[0]=10;
      cmap->high[0]=maxval * 1.01;
      cmap->opacity_high[0]=0.8;
      cmap->high_R[0]=240;
      cmap->high_G[0]=100;
      cmap->high_B[0]=0;
    }
    else if(strcmp("Strn",fieldname)==0) {  /* `Strain' (if yielding active)  */
      found=1;

      if(maxval < E->viscosity.yield_strain_onset) {
	cmap->sections=1;
	cmap->background_R=-1;
	cmap->background_G=-1;
	cmap->background_B=-1;
	/* section 0 */
	cmap->low[0]=aveval;
	cmap->opacity_low[0]=0.0;
	cmap->low_R[0]=20;
	cmap->low_G[0]=250;
	cmap->low_B[0]=100;
	cmap->high[0]=E->viscosity.yield_strain_onset;
	cmap->opacity_high[0]=0.01;
	cmap->high_R[0]=100;
	cmap->high_G[0]=250;
	cmap->high_B[0]=200;
      }
      else {
	cmap->sections=1;
	cmap->background_R=-1;
	cmap->background_G=-1;
	cmap->background_B=-1;
	/* section 0 */
	cmap->low[0]=min(0.5 *(E->viscosity.yield_strain_onset+maxval), 
			 0.5 *(E->viscosity.yield_strain_onset+E->viscosity.yield_strain_saturate));
	cmap->opacity_low[0]=0.2;
	cmap->low_R[0]=20;
	cmap->low_G[0]=200;
	cmap->low_B[0]=100;
	cmap->high[0]=E->viscosity.yield_strain_saturate*1.01;
	cmap->opacity_high[0]=0.9;
	cmap->high_R[0]=10;
	cmap->high_G[0]=200;
	cmap->high_B[0]=150;
      }
    }

    else {  /* default grayscale map */
      cmap->sections=1;
      cmap->background_R=-1;
      cmap->background_G=-1;
      cmap->background_B=-1;
      /* section 0 */
      cmap->low[0]=minval;
      cmap->opacity_low[0]=1.0;
      cmap->low_R[0]=20;
      cmap->low_G[0]=250;
      cmap->low_B[0]=140;
      cmap->high[0]=maxval;
      cmap->opacity_high[0]=1.0;
      cmap->high_R[0]=255;
      cmap->high_G[0]=255;
      cmap->high_B[0]=255;
    }
  return(found);
}


/* Resampling 2D data to evenly spaced grid using
   simplest conceivable linear interpolation  */

void  simple_interp_2D(
  struct All_variables *E,
  standard_precision *orig_data,
  standard_precision *interp_data,
  int dirn1,
  int new1,
  int dirn2,
  int new2
)
{
  standard_precision *intermediate;
  standard_precision *oldline,*newline;
  int i,j;
  int old1,old2;

  void horiz_respacing_n();

  old1 = E->mesh.nnx[dirn1];
  old2 = E->mesh.nnx[dirn2];

  /* allocate mem for intermediate data */

  intermediate = (standard_precision *) Malloc0((new1*old2+1) * sizeof(standard_precision));
  oldline = (standard_precision *) Malloc0((old1+old2+1) * sizeof(standard_precision));
  newline = (standard_precision *) Malloc0((new1+new2+1) * sizeof(standard_precision));

  /* interpolate in direction 1 first */

  for(j=1;j<=old2;j++) {
    for(i=1;i<=old1;i++) 
      oldline[i] = orig_data[j+(i-1)*old2];
    
    horiz_respacing_n(E,oldline,newline,dirn1,new1);

    for(i=1;i<=new1;i++) 
      intermediate[j+(i-1)*old2] = newline[i];
  }
  
  /* interpolate in direction 2 next */

  for(i=1;i<=new1;i++) {
    for(j=1;j<=old2;j++)
      oldline[j] = intermediate[j+(i-1)*old2];

    horiz_respacing_n(E,oldline,newline,dirn2,new2);
 
    for(j=1;j<=new2;j++) 
      interp_data[j+(i-1)*new2] = newline[j];
  }
  
  free((void *) oldline);
  free((void *) newline);
  free((void *) intermediate);

  return;
}

/* As horiz_respacing but for specified point density */

void horiz_respacing_n(
  struct All_variables *E,
  standard_precision *f,
  standard_precision *fr,
  int dirn,
  int num
)
{
  standard_precision spacing,currentx;
  int i,bigger,smaller,node_inc;

  /* x dirn increments by noz between sweeps, ydirn increments by nox*noz 
     this gives locations along the axis. (z dirn would be inc of 1) */

  switch (dirn) { 
  case 1:  /* X */
    node_inc = E->mesh.noz;
    break;
  case 2: 
    node_inc = 1;
    break;
  case 3:
    node_inc = E->mesh.noz * E->mesh.nox;
    break;

  }

  spacing = (E->x[dirn][E->mesh.nno] - E->x[dirn][1])/((standard_precision) num - 1.0);

  fr[1] = f[1];
  fr[num] = f[E->mesh.nnx[dirn]];

  for(i=2;i<num;i++) {
    currentx = ((standard_precision)i - 1.0) * spacing + E->x[dirn][1];
    smaller = 1;

    while (currentx >= E->x[dirn][1+(smaller-1)*node_inc]) 
	smaller++;
      
    bigger = smaller - 1;

    fr[i] = (f[bigger] +
	     (f[smaller]-f[bigger]) *
	     (currentx - E->x[dirn][1+(bigger-1)*node_inc])/
	     (E->x[dirn][1+(smaller-1)*node_inc] - E->x[dirn][1+(bigger-1)*node_inc]));

    /* fprintf(stderr,"%f lies between %f and %f \n",fr[i],f[bigger],f[smaller]); */

    /* next point */
  }

  return; 
}
