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

#if HAVE_STDLIB_H
#include <stdlib.h>
#endif

#if HAVE_STRING_H
#include <string.h>
#endif

#if HAVE_STRINGS_H
#include <strings.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"


int read_bc_from_file(
		      struct All_variables *E,
		      standard_precision * field,
		      unsigned int *flags,
		      char *name,
		      char *abbr,
		      unsigned int ON,
		      unsigned int OFF
		      )
{
    int input_string();

    standard_precision cross2d();

    char discard[5001];
    char * token;
 
    char *filename;
    char *input_token;
    FILE *fp;
    int fnodesx,fnodesz,fnodesy;
    int i,j,k,node2d,node3d,column,found;
    int interpolate=0;
    double value;

    standard_precision *T;
    standard_precision yloc;
    char instring[500];

    filename=(char *)Malloc0(500*sizeof(char));
    input_token=(char *)Malloc0(1000*sizeof(char));

    /* Define field name, read parameter file to determine file name and column number */

    sprintf(input_token,"%s_bc_file",name);
    if(!input_string(input_token,filename,"initialize")) {
	fprintf(E->fp,"No %s bc information found in input file\n",name);fflush(E->fp);
	return(0);   /* if not found, take no further action, return zero */
    }

    fprintf(E->fp,"%s bc information is in file %s\n",name,filename);fflush(E->fp);
 
    /* Try opening the file, fatal if this fails too */

    if((fp=fopen(filename,"r")) == NULL) {
	fprintf(E->fp,"Unable to open the required file `%s'\n",filename);fflush(E->fp);
	if(E->control.verbose)
	   	fprintf(stderr,"Unable to open the required file `%s'\n",filename);
	return(0);
    }

     /* Read header, get nodes xzy */

    fgets(discard,4999,fp);
    fgets(discard,4999,fp); 
    i=sscanf(discard,"# NODESX=%d NODESZ=%d NODESY=%d",&fnodesx,&fnodesz,&fnodesy);
    if(i<3) {
	fprintf(E->fp,"File %s is not in the correct format\n",filename);fflush(E->fp);
	return(0);
    }

    fgets(discard,4999,fp); /* largely irrelevant line */
    fgets(discard,4999,fp);
    
    /* This last line is the column headers, we need to search for the occurence of abbr to
       find out the column to be read in. */
    if(!strtok(discard,"|")) { 
	fprintf(E->fp,"Unable to decipher the columns in the input file");fflush(E->fp);
	return(0);
    }

    found=0;
    column=1;

    while(found==0 && (token=(char *)strtok(NULL,"|"))) {
	if(strstr(token,abbr)!=0) found=1;
	column++; 
    }

    if(found) {
	fprintf(E->fp,"\t%s (%s) found in column %d\n",name,abbr,column);fflush(E->fp);
    }    
    else {
	fprintf(E->fp,"\t%s (%s) not found in file: %s\n",name,abbr,filename);fflush(E->fp);
	return(0);
    }

    /* A fatal condition - file size is all wrong */
    if(fnodesx != E->mesh.nox || fnodesz != E->mesh.noz || fnodesy != E->mesh.noy ) {
	fprintf(E->fp,"Input data for file `%s' is not %d x %d x %d\n",
		filename,E->mesh.nox,E->mesh.noz,E->mesh.noy);
	fflush(E->fp);
	exit(1);
    }
  
   /* Format for reading the input file  */

    sprintf(input_token," %%d ");
    for(i=2;i<column;i++)
	strcat(input_token," %*s");
    strcat(input_token," %s");

    for(i=1;i<=fnodesx*fnodesz*fnodesy;i++) {
      if(fgets(discard,4999,fp)==NULL) /* EOF encountered */
	return(2);
      sscanf(discard,input_token,&j,instring); 

      if(j<1 || j > E->mesh.nno)
	fprintf(stderr,"Illegal node value in %s bc file line %d: %d\n",name,i+4,j);
	
      if(strcmp(instring,"U")!=0)  { /* This field is not unconstrained (U) at this node */
	sscanf(instring,"%lf",&value);
	/*fprintf(stderr,"Found %s bc at node %d value %g - %s,%s\n",name,j,value,instring,discard); */
	field[j] = (standard_precision) value;
	flags[j] = flags[j] | ON;
	flags[j] = flags[j] & (~OFF);
      }
    }
    fclose(fp); 
    free((void *)filename);
    free((void *)input_token);
    return(1);
}

void arbitrary_bc_rectangle_file(
				 struct All_variables *E,
				 struct RectBc *RECT,
				 char *name,
				 standard_precision **field,
				 unsigned int **bcbitf,
				 const unsigned int bcmask_on,
				 const unsigned int bcmask_off,
				 const unsigned int surf
				 )
{
  char read_string[500];
  void arbitrary_bc_rectangle();


  sprintf(read_string,"%s_bc_rect",name);
  input_int(read_string,&(RECT->numb),"0");
  sprintf(read_string,"%s_bc_rect_aa1",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->aa1);
  sprintf(read_string,"%s_bc_rect_aa2",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->aa2);
  sprintf(read_string,"%s_bc_rect_bb1",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->bb1);
  sprintf(read_string,"%s_bc_rect_bb2",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->bb2);
  sprintf(read_string,"%s_bc_rect_hw",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->halfw);
  sprintf(read_string,"%s_bc_rect_mag",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->mag);
  sprintf(read_string,"%s_bc_rect_icpt",name);
  input_std_precision_vector(read_string,RECT->numb,RECT->intercept);
  sprintf(read_string,"%s_bc_rect_norm",name);
  input_char_vector(read_string,RECT->numb,RECT->norm);

  arbitrary_bc_rectangle(E,RECT,field,bcbitf,bcmask_on,bcmask_off,surf);  
     
  return;
}
  
void arbitrary_bc_circle_file(
    struct All_variables *E,
    struct CircBc *CIRC,
    char *name,
    standard_precision **field,
    unsigned int **bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    const unsigned int surf
)
{
  char read_string[500];
  void arbitrary_bc_circle();

  sprintf(read_string,"%s_bc_circ",name);
  input_int(read_string,&(CIRC->numb),"0");
  sprintf(read_string,"%s_bc_circ_aa",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->aa);
  sprintf(read_string,"%s_bc_circ_bb",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->bb);
  sprintf(read_string,"%s_bc_circ_rad",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->rad);
  sprintf(read_string,"%s_bc_circ_mag",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->mag);
  sprintf(read_string,"%s_circhw",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->halfw);
  sprintf(read_string,"%s_bc_circ_icpt",name);
  input_std_precision_vector(read_string,CIRC->numb,CIRC->intercept);
  sprintf(read_string,"%s_bc_circ_norm",name);
  input_char_vector(read_string,CIRC->numb,CIRC->norm);
 
  arbitrary_bc_circle(E,CIRC,field,bcbitf,bcmask_on,bcmask_off,surf);
 
  return;  
}

void arbitrary_bc_harmonic_file(
  struct All_variables *E,
  struct HarmBc *HARM,
  char *name,
  standard_precision **field,
  unsigned int **bcbitf,
  unsigned int bcmask_on,
  unsigned int bcmask_off,
    const unsigned int surf
)
{
  char read_string[500];
  void arbitrary_bc_harmonic();
  int i;

  sprintf(read_string,"%s_bc_harm",name);
  input_int(read_string,&(HARM->numb),"0");
  sprintf(read_string,"%s_bc_harms",name);
  input_int(read_string,&(HARM->harms),"0,0,19");
  sprintf(read_string,"%s_bc_harm_norm",name);
  input_char_vector(read_string,HARM->numb,HARM->norm);
  sprintf(read_string,"%s_bc_harm_ofst",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->off);
  sprintf(read_string,"%s_bc_harm_icpt",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->intercept);
  sprintf(read_string,"%s_bc_harm_aa1",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->aa1);
  sprintf(read_string,"%s_bc_harm_aa2",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->aa2);
  sprintf(read_string,"%s_bc_harm_bb1",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->bb1);
  sprintf(read_string,"%s_bc_harm_bb2",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->bb2);

  for(i=0;i<HARM->harms;i++) {
      sprintf(read_string,"%s_bc_harm_kaa%02d",name,i);
      input_std_precision_vector(read_string,HARM->numb,HARM->kaa[i]);
      sprintf(read_string,"%s_bc_harm_kbb%02d",name,i);
      input_std_precision_vector(read_string,HARM->numb,HARM->kbb[i]);
      sprintf(read_string,"%s_bc_harm_amp%02d",name,i);
      input_std_precision_vector(read_string,HARM->numb,HARM->amp[i]);
      sprintf(read_string,"%s_bc_harm_phaa%02d",name,i);
      input_std_precision_vector(read_string,HARM->numb,HARM->phaa[i]);
      sprintf(read_string,"%s_bc_harm_phbb%02d",name,i);
      input_std_precision_vector(read_string,HARM->numb,HARM->phbb[i]);
  }
   
  arbitrary_bc_harmonic(E,HARM,field,bcbitf,bcmask_on,bcmask_off,surf);
 
  return;  
}

void arbitrary_bc_polynomial_file(
    struct All_variables *E,
    struct PolyBc *POLY,
    char *name,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    const unsigned int surf
)
{
  char read_string[500], pfx;
  void arbitrary_bc_polynomial();
  int i;

  /*RAA: 5/8/01, eventually, there should be more coefficients for order>1,
       for cross-terms (a6*x*y, etc) */

  sprintf(read_string,"%s_bc_poly",name);
  input_int(read_string,&(POLY->numb),"0");
  sprintf(read_string,"%s_bc_poly_ord",name);
  input_int(read_string,&(POLY->order),"0,0,19");
  sprintf(read_string,"%s_bc_poly_norm",name);
  input_char_vector(read_string,POLY->numb,POLY->norm);
  sprintf(read_string,"%s_bc_poly_icpt",name);
  input_std_precision_vector(read_string,POLY->numb,POLY->intercept);
  sprintf(read_string,"%s_bc_poly_aa1",name);
  input_std_precision_vector(read_string,POLY->numb,POLY->aa1);
  sprintf(read_string,"%s_bc_poly_aa2",name);
  input_std_precision_vector(read_string,POLY->numb,POLY->aa2);
  sprintf(read_string,"%s_bc_poly_bb1",name);
  input_std_precision_vector(read_string,POLY->numb,POLY->bb1);
  sprintf(read_string,"%s_bc_poly_bb2",name);
  input_std_precision_vector(read_string,POLY->numb,POLY->bb2);

  /*RAA: 5/7/01, aaa%02d (abb%02d) used to be oaa.. (obb..),
    should be <= below, not just < !!*/
  /* DAS: 6/1/04, make this work for both _oaa and _aaa forms,
     as long as it's consistent for each polynomial */
  sprintf(read_string,"%s_bc_poly_aaa00",name);
  if(input_std_precision_vector(read_string,POLY->numb,POLY->aaa[i]))
      pfx = 'o';
  else
      pfx = 'a';
  for(i=0;i<=POLY->order;i++) {
      sprintf(read_string,"%s_bc_poly_%caa%02d",name,pfx,i);
      input_std_precision_vector(read_string,POLY->numb,POLY->aaa[i]);
      sprintf(read_string,"%s_bc_poly_%cbb%02d",name,pfx,i);
      input_std_precision_vector(read_string,POLY->numb,POLY->abb[i]);
      }

  arbitrary_bc_polynomial(E,POLY,field,bcbitf,bcmask_on,bcmask_off,surf);
 
  return;  
}
