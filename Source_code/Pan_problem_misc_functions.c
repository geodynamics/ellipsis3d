/*

Copyright (C) 1995 The GeoFramework Consortium

This file is part of Ellipsis3D.

Ellipsis3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, version 2,
as published by the Free Software Foundation.

Ellipsis3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Author:
  Louis Moresi <louis.moresi@sci.monash.edu>

*/


#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

#include <stdlib.h>

#if (!defined __GNUC__)
#include <unistd.h>
#endif

#if (defined __sunos__) || defined(__uxp__) || defined(__GNUC__)
#include <string.h>
#else
#include <strings.h>
#endif

#if defined(__sgi) || defined(__osf__)
#include <sys/types.h>
#endif

extern int Emergency_stop;

/* #include "/home/limbo1/louis/Software/include/dmalloc.h" */

void interuption()
{  
  if (Emergency_stop++) exit(0);
  fprintf(stderr,"Cleaning up before exit\n");
  return; 
}

void twiddle_thumbs(
  struct All_variables *yawn,
  int scratch_groin
)
{ /* Do nothing, just sit back and relax.
     Take it easy for a while, maybe size
     doesn't matter after all. There, there
     that's better. Now ... */
   return;
}

int get_process_identifier()   
{
   int pid;

   pid = (int) getpid();
   return(pid);
}

void report(
  struct All_variables *E,
  char * string
  )
{
   if(E->control.verbose){
      fprintf(stderr,"%s\n",string);
      fflush(stderr);
   }
   return;
}

void record(
  struct All_variables *E,
  char * string
  )
{
   if(E->control.verbose) {
      fprintf(E->fp,"%s\n",string);
      fflush(E->fp);
   }
   return;
}

double SIN_D(
     double x
)
{
#if defined(__osf__)
  return sind(x);
#else
  return sin((x/180.0) * M_PI);
#endif

}

double COT_D(
     double x
)
{
#if defined(__osf__)
  return cotd(x);
#else
  return tan(((90.0-x)/180.0) * M_PI);
#endif

}

/* =================================================
   Useful atan mod for spherical geometry
   ================================================= */

higher_precision myatan(
 higher_precision y,
 higher_precision x
)
 {
 higher_precision fi;

 fi = atan2(y,x);

 if (fi<0.0)
    fi += 2*M_PI;

 return(fi);
 }

/* =================================================
   Sobel sequence (pseudo-random, space-filling) For
   explanation, see Numerical Recipes
   ================================================= */ 

#define MAXBIT 30
#define MAXDIM 6


void sobol(
  int *n,
  standard_precision x[4]
)
{

  int j,k,l;
  unsigned i,im,ipp;
  static standard_precision fac;
  static int been_here=0;
  static unsigned in,ix[MAXDIM+1],*iu[MAXBIT+1];
  static unsigned mdeg[MAXDIM+1] =      {0,1,2,3,3,4,4};
  static unsigned ip[MAXDIM+1]   =      {0,0,1,1,2,1,4};
  static unsigned iv[MAXDIM*MAXBIT+1] = {0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

  if(*n < 0) {
    if(been_here++ ==0) {
      for(j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
	iu[j] = &(iv[k]);

      for(k=1;k<=MAXDIM;k++) {
	for(j=1;j<=mdeg[k];j++)
	  iu[j][k] <<= (MAXBIT-j);

	for(j=mdeg[k]+1;j<=MAXBIT;j++) {
	  ipp=ip[k];
	  i=iu[j-mdeg[k]][k];
	  i^= (i>>mdeg[k]);
	
	  for(l=mdeg[k]-1;l>=1;l--) {
	    if(ipp & 1) 
	      i ^= iu[j-1][k];
	    ipp >>= 1;
	  }
	  iu[j][k]=i;
	}
      }
    }
    
    fac=1.0/(1L << MAXBIT);
    in = 0;
  }

  else {
    im=in;
    for(j=1;j<=MAXBIT;j++) {
      if(!(im & 1)) 
	break;
      im >>= 1;
    }
    if(j > MAXBIT) {
      fprintf(stderr,"Error in Sobel sequence routine !\n");
      exit(-1);
    }

    im=(j-1)*MAXDIM;

    for(k=1;k<=min(*n,MAXDIM);k++) {
      ix[k] ^= iv[im+k];
      x[k] = ix[k]*fac;
    }
    in++;
  }

  return;
}

/* non-runaway malloc */

void * Malloc1
( int bytes,
  char *file,
  int line ) {

    void *ptr;

    /* 
       Check to make sure the memory could be
       allocated and halt if the allocation fails 
    */

    ptr = malloc((size_t)bytes);
    if (ptr == (void *)NULL) {
	fprintf(stderr,"Memory: cannot allocate another %d bytes \n(line %d of file %s)\n",bytes,line,file);
	exit(0);
    }

    return(ptr);
}

/* And a version of realloc which complains
   if it is used to free memory and checks to 
   see if the increased allocation worked */


void Realloc1(
	      void **allocate_me,
	      size_t size0, 
	      size_t size1,
	      size_t entity_size,
	      char *file,
	      int line  
	      )
{  
  void * ptr;
  
  /* 1) size1 == 0, or size0 > size1. 
     This is now made to flag an error - realloc etc don't mind
     this sort of thing and will simply free the memory but this makes inadvertant 
     passing of a zero size difficult to find
  */

  if (size1==0 || size0 > size1) {
    fprintf(stderr,"Error - incorrect use of Realloc %d -> %d bytes at \n(line %d of file %s)\n",
	    size0,size1,line,file);
    exit(1);
  }

  if(size1 == size0) {
    fprintf(stderr,"Warning - using Realloc with same size input/output - no action\n");
    return;
  }


  /* 2) Check if initially a null pointer
     in which case memory will need to be allocated
     rather than reallocated. This would appear to be
     doing the job of realloc but for some reason, this
     doesn't seem to work correctly otherwise !
  */

  if(0 == size0 || *allocate_me == (void *)NULL) {
    /* fprintf(stderr,"This pointer needs initializing ... "); */
    *allocate_me = malloc((size_t) size1 * entity_size);
    /* fprintf(stderr,"%d\n",allocate_me); */   
  } 

  /* 3) Pointer has been allocated and needs
     to be reallocated and the old data copied
     so we use realloc this time */

  else { 
    
    /* fprintf(stderr,"Pointer exists - copy to new location %d -> ? \n",*allocate_me); */

    ptr = *allocate_me;
    *allocate_me = malloc((size_t) size1 * entity_size);

    /* fprintf(stderr,"Pointer exists - new location %d  \n",*allocate_me); */

    memcpy(*allocate_me, ptr, size0 * entity_size);

    /* fprintf(stderr,"Pointer exists - completed copy   \n"); */

    free(ptr);
  }

  if (*allocate_me == (void *)NULL) {
    fprintf(stderr,"Memory: cannot reallocate a block of %d bytes \n(line %d of file %s)\n",size1,line,file);
    exit(0);
  }
  
  /*
  fprintf(stderr,"Successfully realloc'd %d bytes for request at line %d of file %s -> %d \n",
	  size1*entity_size,line,file,*allocate_me);
	  */
}


/* Read in a file containing previous values of a field. The input in the parameter
   file for this should look like: `previous_name_file=string' and `previous_name_column=int' 
   where `name' is substituted by the argument of the function. 

   The file should have the standard CITCOM output format:
     # HEADER LINES etc
     index X Z Y ... field_value1 ...
     index X Z Y ... field_value2 ...
   where index is the node number, X Z Y are the coordinates and
   the field value is in the column specified by the abbr term in the function argument

   If the number of nodes OR the XZY coordinates for the node number (to within a small tolerance)
   are not in agreement with the existing mesh, the data is interpolated. 
   */

int read_previous_field(
    struct All_variables *E,
    standard_precision * field,
    char *name,
    char *abbr
)
{
    void fcopy_interpolating();
    int input_string();

    standard_precision cross2d();

    char discard[5001];
    char *token;
    char *filename;
    char *input_token;
    FILE *fp;
    int fnodesx,fnodesz,fnodesy;
    int i,j,k,node2d,node3d,column,found;
    int interpolate=0;

    standard_precision *X,*Z,*Y,*T;
    standard_precision yloc;
    double in1,in2,in3,in4;

    filename=(char *)Malloc0(500*sizeof(char));
    input_token=(char *)Malloc0(1000*sizeof(char));

    /* Define field name, read parameter file to determine file name and column number */

    sprintf(input_token,"previous_%s_file",name);
    if(!input_string(input_token,filename,"initialize")) {
	fprintf(E->fp,"No previous %s information found in input file\n",name);fflush(E->fp);
	return(0);   /* if not found, take no further action, return zero */
    }

    fprintf(E->fp,"Previous %s information is in file %s\n",name,filename);fflush(E->fp);
 
    /* Try opening the file, fatal if this fails too */

    if((fp=fopen(filename,"r")) == NULL) {
	fprintf(E->fp,"Unable to open the required file `%s'",filename);fflush(E->fp);
	if(E->control.verbose)
	   	fprintf(stderr,"Unable to open the required file `%s'",filename);
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
    
    /* this last line is the column headers, we need to search for the occurence of abbr to
       find out the column to be read in */

    if(!strtok(discard,"|")) { 
	fprintf(E->fp,"Unable to deciphre the columns in the input file");fflush(E->fp);
	return(0);
    }

    found=0;
    column=1;

    while(found==0 && (token=(char *)strtok(NULL,"|"))) {
	if(strstr(token,abbr)!=0)
	    found=1;
	column++;
    }

    if(found) {
	fprintf(E->fp,"\t%s (%s) found in column %d\n",name,abbr,column);fflush(E->fp);
    }    
    else {
	fprintf(E->fp,"\t%s (%s) not found in file: %s\n",name,abbr,filename);fflush(E->fp);
	return(0);
    }

    /* Another fatal condition (not suitable for interpolation: */
    if((3!= E->mesh.nsd) && (fnodesy !=1)) {
	fprintf(E->fp,"Input data for file `%s'  is of inappropriate dimension (not %dD)\n",filename,E->mesh.nsd);fflush(E->fp);
	exit(1);
    }

    X=(standard_precision *)Malloc0((2+fnodesx*fnodesz*fnodesy)*sizeof(standard_precision));
    Z=(standard_precision *)Malloc0((2+fnodesx*fnodesz*fnodesy)*sizeof(standard_precision));
    Y=(standard_precision *)Malloc0((2+fnodesx*fnodesz*fnodesy)*sizeof(standard_precision));
    T=(standard_precision *)Malloc0((2+fnodesx*fnodesz*fnodesy)*sizeof(standard_precision));
    
   /* Format for reading the input file (including coordinates) */

    sprintf(input_token," %%d %%le %%le %%le");
    for(i=5;i<column;i++)
	strcat(input_token," %*f");
    strcat(input_token," %lf");

    for(i=1;i<=fnodesx*fnodesz*fnodesy;i++) {
	fgets(discard,4999,fp);
	sscanf(discard,input_token,&j,&in1,&in2,&in3,&in4);
	X[i] = (standard_precision) in1;
	Z[i] = (standard_precision) in2;
	Y[i] = (standard_precision) in3;
	T[i] = (standard_precision) in4;
    }
    fclose(fp);
    
    /* first we extrude any 2d files into the third dimension (the interpolation can be done afterwards) */
    
    if(3==E->mesh.nsd && 1==fnodesy) {
      fprintf(E->fp,"Extrude the 2d file into the third dimension\n"); 
      fflush(E->fp);
      
      /* Make more room for fields ... */
      T = (standard_precision *)realloc(T,(2+E->mesh.noz*E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
      X = (standard_precision *)realloc(X,(2+E->mesh.noz*E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
      Y = (standard_precision *)realloc(Y,(2+E->mesh.noz*E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
      Z = (standard_precision *)realloc(Z,(2+E->mesh.noz*E->mesh.nox*E->mesh.noy)*sizeof(standard_precision));
      fnodesy=E->mesh.noy;
      for(k=1;k<=E->mesh.noy;k++) {
	yloc = (k-1) * 1.0 / ((standard_precision) E->mesh.noy - 1);
	for(i=1;i<=E->mesh.noz;i++)
	  for(j=1;j<=E->mesh.nox;j++) {
	    node2d = i + (j-1) * E->mesh.noz;
	    node3d = node2d + (k-1) * E->mesh.noz * E->mesh.nox;
	    Y[node3d] = yloc;
	    X[node3d] = X[node2d]; 
	    Z[node3d] = Z[node2d]; 
	    T[node3d] = T[node2d]; 
	  }
      }
    }

    /* check consistency & need for interpolation */

    if(fnodesx != E->mesh.nox || fnodesz != E->mesh.noz || fnodesy != E->mesh.noy)
	interpolate=1;
	
    for(i=1;i<=fnodesx*fnodesz*fnodesy;i++)
	if( fabs(X[i]-E->x[1][i]) > 0.01*fabs(X[i]) ||
	    fabs(Z[i]-E->x[2][i]) > 0.01*fabs(Z[i]) ||
	    ((3==E->mesh.nsd) && fabs(Y[i]-E->x[3][i]) > 0.01*fabs(Y[i]))) {
	    interpolate=i;
	    break;
	}

    if(interpolate!=0) {
	fprintf(E->fp,"\t%s requires interpolation from previous value\n",name);fflush(E->fp);
	fprintf(E->fp,"\tOld nodes = %d/%d/%d and new nodes = %d/%d/%d\n",fnodesx,fnodesz,fnodesy,E->mesh.nox,E->mesh.noz,E->mesh.noy);fflush(E->fp);
	fcopy_interpolating(E,X,Z,Y,fnodesx,fnodesz,fnodesy,T,field);
    }
    else {
	fprintf(E->fp,"\t%s requires no interpolation from previous value\n",name);fflush(E->fp);
	vcopy(field,T,1,E->mesh.nno);
    }

    free((void *)X);
    free((void *)Z);
    free((void *)Y);
    free((void *)T);
    free((void *)filename);
    free((void *)input_token);
    
    return(1);
}

/*
Copy one field to another on a different (but similarly structured) mesh. The
field TT (output) is on the standard mesh for this problem, T (input) is on the
different mesh whose coordinates are described by X,Z,Y 
*/

void fcopy_interpolating(
    struct All_variables *E,
    standard_precision *X,
    standard_precision *Z,
    standard_precision *Y,
    int nx,
    int nz,
    int ny,
    standard_precision *T,
    standard_precision *TT
)
{
    standard_precision cross2d();
    void p_to_nodes();
    void p_to_centres();
    standard_precision CPU_time(),time;

    int i,j,found,not_found;
    int elX,elY,elZ,ex,ez;
    int old_ex,old_ez,old_ey;
    int node1,node2,node3,node4;
    int node5,node6,node7,node8;
    standard_precision inside1,inside2,inside3,inside4;
    standard_precision inside5,inside6,inside7,inside8,inside9,inside10,inside11,inside12;
    standard_precision distance1,distance2,distance3,distance4;
    standard_precision distance5,distance6,distance7,distance8;
    standard_precision d1,d2,d3,d4,d5,d6,d7,d8;

    const int dims=E->mesh.nsd;
         
    /* Run over all the data points (take care to hit only one element), determine
       inside/outside of each element. The criterion for inside-ness is that the
       cross product of the vectors joining the point to the ends of each edge should
       have the same sign as the dot product of the centre to the ends of each edge.
       There are, undoubtedly, better ways to do this !!!

       Because a CITCOM node-ordering is assumed, we can also guess that the best place
       to start looking for a node is where you found the last one !
    */

    old_ex=old_ez=old_ey=1;

    elX = nx-1;
    elZ = nz-1;

     if(E->control.print_convergence) {   
	 time=CPU_time();
	 fprintf(stderr,"Interpolating ...");
     }

    not_found=0;
  
    if(2==dims)
	for(i=1;i<=E->mesh.nno;i++) {
	    found=0;
	    for(ex=old_ex;ex<=elX && found==0;ex++)
		for(ez=1;ez<=elZ && found==0 ;ez++) {
		    node1=ez+(ex-1)*nz;
		    node2=node1+offset[2].vector[1]+offset[2].vector[0]*nz;
		    node3=node1+offset[3].vector[1]+offset[3].vector[0]*nz;
		    node4=node1+offset[4].vector[1]+offset[4].vector[0]*nz;
		    
		    if ((inside1 = cross2d(X[node1]-E->x[1][i],Z[node1]-E->x[2][i],X[node2]-E->x[1][i],Z[node2]-E->x[2][i],3)) <= 0.0 &&
			(inside4 = cross2d(X[node4]-E->x[1][i],Z[node4]-E->x[2][i],X[node1]-E->x[1][i],Z[node1]-E->x[2][i],3)) <= 0.0 &&
			(inside2 = cross2d(X[node2]-E->x[1][i],Z[node2]-E->x[2][i],X[node3]-E->x[1][i],Z[node3]-E->x[2][i],3)) <= 0.0 &&
			(inside3 = cross2d(X[node3]-E->x[1][i],Z[node3]-E->x[2][i],X[node4]-E->x[1][i],Z[node4]-E->x[2][i],3)) <= 0.0) {
			found = node1;
			old_ex=ex;
		    }
		}
	    
	    /* finish the loop if not found */
	    for(ex=1;ex<=old_ex && found==0;ex++)
		for(ez=1;ez<=elZ && found==0 ;ez++) {
		    node1=ez+(ex-1)*nz;
		    node2=node1+offset[2].vector[1]+offset[2].vector[0]*nz;
		    node3=node1+offset[3].vector[1]+offset[3].vector[0]*nz;
		    node4=node1+offset[4].vector[1]+offset[4].vector[0]*nz;
		    
		    if ((inside1 = cross2d(X[node1]-E->x[1][i],Z[node1]-E->x[2][i],X[node2]-E->x[1][i],Z[node2]-E->x[2][i],3)) <= 0.0 &&
			(inside4 = cross2d(X[node4]-E->x[1][i],Z[node4]-E->x[2][i],X[node1]-E->x[1][i],Z[node1]-E->x[2][i],3)) <= 0.0 &&
			(inside2 = cross2d(X[node2]-E->x[1][i],Z[node2]-E->x[2][i],X[node3]-E->x[1][i],Z[node3]-E->x[2][i],3)) <= 0.0 &&
			(inside3 = cross2d(X[node3]-E->x[1][i],Z[node3]-E->x[2][i],X[node4]-E->x[1][i],Z[node4]-E->x[2][i],3)) <= 0.0) {
			found = node1;
			old_ex=ex;
		    }
		}

	    /* and having found the right node location, interpolate the appropriate value to it */

	    if(!found)
		not_found++;
	    else {
		distance1 = ((X[node1]-E->x[1][i])*(X[node1]-E->x[1][i])+(Z[node1]-E->x[2][i])*(Z[node1]-E->x[2][i]));
		distance2 = ((X[node2]-E->x[1][i])*(X[node2]-E->x[1][i])+(Z[node2]-E->x[2][i])*(Z[node2]-E->x[2][i]));
		distance3 = ((X[node3]-E->x[1][i])*(X[node3]-E->x[1][i])+(Z[node3]-E->x[2][i])*(Z[node3]-E->x[2][i]));
		distance4 = ((X[node4]-E->x[1][i])*(X[node4]-E->x[1][i])+(Z[node4]-E->x[2][i])*(Z[node4]-E->x[2][i]));
	    
		d1=distance2*distance3*distance4;
		d2=distance1*distance3*distance4;
		d3=distance2*distance1*distance4;
		d4=distance2*distance3*distance1;
		
		TT[i] = (d1*T[node1]+d2*T[node2]+d3*T[node3]+d4*T[node4])/(d1+d2+d3+d4);
	    }
	}

    else {
	elY = (3==dims)? ny-1 : 1;
	fprintf(stderr,"3D interpolator not yet implemented !!! (ignoring need to interpolate)\n");
   	vcopy(TT,T,1,E->mesh.nno); }
   
    if(E->control.print_convergence)
	fprintf(stderr,". done (%f secs)\n",CPU_time()-time);

    if(not_found)
	fprintf(E->fp,"Warning: unable to interpolate old  data to %d nodes in the new mesh\n",not_found);

    return;
}

/*
Return the out of plane component of the cross product of the two vectors
assuming that one is looking AGAINST the direction of the axis of D,
anti-clockwise angles are positive (are you sure ?), and the axes are
ordered 2,3 or 1,3 or 1,2 .
*/
standard_precision cross2d(
    standard_precision x11,
    standard_precision x12,
    standard_precision x21,
    standard_precision x22,
    int D
)
{
  

   if(1==D)
       return( x11*x22-x12*x21);
   if(2==D) 
       return(-x11*x22+x12*x21);
   if(3==D)
       return( x11*x22-x12*x21);

   /* error condition */

   return(0.0);

}

standard_precision magcross2d(
    standard_precision x11,
    standard_precision x12,
    standard_precision x21,
    standard_precision x22
)
{
  return(fabs(x11*x22-x12*x21));
}

void field_arbitrary_rectangle_file(
    struct All_variables *E,
    int parse_and_apply,
    struct Rect *RECT,
    char *name,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
    char read_string[500];
    standard_precision weight,radius2,weight2;
    standard_precision x1,y1,z1;
    int in1,in2,in3;
    int number,node;
    int combine_option;

    void field_arbitrary_rectangle();

    sprintf(read_string,"%s_rect",name);
    input_int(read_string,&(RECT->numb),"0");
    sprintf(read_string,"%s_rect_x1",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->x1);
    sprintf(read_string,"%s_rect_x2",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->x2);
    sprintf(read_string,"%s_rect_z1",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->z1);
    sprintf(read_string,"%s_rect_z2",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->z2);
    sprintf(read_string,"%s_rect_y1",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->y1);
    sprintf(read_string,"%s_rect_y2",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->y2);
    sprintf(read_string,"%s_rect_hw",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->halfw);
    sprintf(read_string,"%s_rect_mag",name);
    input_std_precision_vector(read_string,RECT->numb,RECT->mag);
    sprintf(read_string,"%s_rect_ovl",name);
    input_char_vector(read_string,RECT->numb,RECT->overlay);
 
    if(parse_and_apply)
	field_arbitrary_rectangle(E,RECT,field,bcbitf,bcmask_on,bcmask_off,level);
 
    return;  
}

void field_arbitrary_rectangle(
    struct All_variables *E,
    struct Rect *RECT,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
    standard_precision weight,radius2,weight2;
    standard_precision x1,y1,z1;
    int in1,in2,in3;
    int number,node;
    int combine_option;
    int change;

 
    for(node=1;node<=E->mesh.NNO[level];node++) {
     change=0;
      x1=E->SX[level][1][node];
      z1=E->SX[level][2][node];
      y1=(E->mesh.nsd!=3) ? 0.0 : E->SX[level][3][node];
	
	for(number=0;number<RECT->numb;number++) {
	  switch (RECT->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }

	    in1=(x1 >= RECT->x1[number] && x1 <= RECT->x2[number]);
	    in2=(z1 >= RECT->z1[number] && z1 <= RECT->z2[number]);
	    in3=(3!=E->mesh.nsd || (y1 >= RECT->y1[number] && y1 <= RECT->y2[number]));
	    
	    if(in1 && in2 && in3) {
		weight = RECT->mag[number];
		radius2=0.0;
		change++;
	    }
	    else {
		radius2 = 
		    (in1 ? 0.0 : min(1.0e-10+(x1-RECT->x1[number])*(x1-RECT->x1[number]),
				     1.0e-10+(x1-RECT->x2[number])*(x1-RECT->x2[number]))) +
		    (in2 ? 0.0 : min(1.0e-10+(z1-RECT->z1[number])*(z1-RECT->z1[number]),
				     1.0e-10+(z1-RECT->z2[number])*(z1-RECT->z2[number]))) +
		    (in3 ? 0.0 : min(1.0e-10+(y1-RECT->y1[number])*(y1-RECT->y1[number]),
				     1.0e-10+(y1-RECT->y2[number])*(y1-RECT->y2[number])));
		
		weight += RECT->mag[number]*exp(-radius2/(1.0e-16 + 0.5 * RECT->halfw[number]*RECT->halfw[number]));
		if(radius2 < RECT->halfw[number]* RECT->halfw[number])
		  change++;
	    }
		
	    if(field!=NULL)
	      switch (RECT->overlay[number]) {
	      case 'R':
		weight2=1.0/(1.0+radius2/(1.0e-16+RECT->halfw[number]*RECT->halfw[number]));
		if(weight2 > 1.0e-10 && weight > 1.0e-10) {
		  field[node]=(1.0-weight2)*field[node]+weight2*weight;
		  
		  /* fprintf(stderr,"Node %d (%g,%g,%g -> %d,%d,%d) is fixed at T=%g ... %g,%g,%g,%g\n",
		     node,x1,z1,y1,in1,in2,in3,field[node],weight,weight2,radius2,RECT->mag[number]); */
		}
		break;
	      case 'M':
		field[node]*=weight;
		break;
	      case 'A':
		field[node]+=weight;
		break;
	      default:
		fprintf(E->fp,"RECTANGLE: %d can't work out how to combine new/old fields\n",number);
	      break;
	    }
	    
	    if(change && bcbitf!=NULL) {
	      bcbitf[node] = (bcbitf[node] | bcmask_on);
	      bcbitf[node] = (bcbitf[node] & (~bcmask_off));
	    }
	}
    }
    return;  
}

void field_arbitrary_circle_file(
    struct All_variables *E,
    int parse_and_apply,
    struct Circ *CIRC,
    char *name,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
    
    char read_string[500];
    standard_precision weight,radius2,weight2;
    standard_precision x1,y1,z1;
    int in1;
    int number,node;
    int combine_option;
    void field_arbitrary_circle();
    
    sprintf(read_string,"%s_circ",name);
    input_int(read_string,&(CIRC->numb),"0");
    sprintf(read_string,"%s_circ_x",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->x);
    sprintf(read_string,"%s_circ_z",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->z);
    sprintf(read_string,"%s_circ_y",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->y);
    sprintf(read_string,"%s_circ_rad",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->rad);
    sprintf(read_string,"%s_circ_mag",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->mag);
    sprintf(read_string,"%s_circ_hw",name);
    input_std_precision_vector(read_string,CIRC->numb,CIRC->halfw);
    sprintf(read_string,"%s_circ_ovl",name);
    input_char_vector(read_string,CIRC->numb,CIRC->overlay);
 
    if(parse_and_apply)
	field_arbitrary_circle(E,CIRC,field,bcbitf,bcmask_on,bcmask_off,level);
   
    return;  
}

void field_arbitrary_circle(
    struct All_variables *E,
    struct Circ *CIRC,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
    char read_string[500];
    standard_precision weight,radius2,weight2;
    standard_precision x1,y1,z1;
    int in1;
    int number,node;
    int combine_option;
    int change;
   
    for(node=1;node<=E->mesh.NNO[level];node++) {
      change=0;
	x1=E->SX[level][1][node];
	z1=E->SX[level][2][node];
	y1=(E->mesh.nsd!=3) ? 0.0 : E->SX[level][3][node];
	
	for(number=0;number<CIRC->numb;number++) {
	    switch (CIRC->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }
	    
	    radius2 =
		(x1-CIRC->x[number])*(x1-CIRC->x[number]) +
		(z1-CIRC->z[number])*(z1-CIRC->z[number]) +
		((E->mesh.nsd != 3) ? 0.0 : (y1-CIRC->y[number])*(y1-CIRC->y[number]));
 
	    if(radius2 <= CIRC->rad[number]*CIRC->rad[number]) {
		weight = CIRC->mag[number];
		radius2=0.0;
	    }
	    else {
		radius2 -= CIRC->rad[number] * CIRC->rad[number];
		weight += CIRC->mag[number]*exp(-2.0*radius2/(1.0e-16+CIRC->halfw[number]*CIRC->halfw[number]));
	    }

		switch (CIRC->overlay[number]) {
		case 'R':
		  if(radius2 > CIRC->halfw[number]*CIRC->halfw[number])
		    break;
		  weight2=1.0/(1.0+radius2/(1.0e-16+CIRC->halfw[number]*CIRC->halfw[number]));
		  if(field != NULL)
		    field[node]=(1.0-weight2)*field[node]+weight2*weight;
		  change++;
		  break;
		case 'M':
		  if(field != NULL) 
		    field[node]*=weight;
		  change++;
		  break;
		case 'A':
		  if(field != NULL)
		    field[node]+=weight;
		    change++;
		    break;
		default:
		    fprintf(E->fp,"CIRCLE: %d can't work out how to combine new/old fields\n",number);
		    break;
		}

	    if(change && bcbitf!=NULL) {
		bcbitf[node] = (bcbitf[node] | bcmask_on);
		bcbitf[node] = (bcbitf[node] & (~bcmask_off));
	    }
	}
    }
    return;  
}

void field_arbitrary_harmonic_file(
    struct All_variables *E,
    int parse_and_apply,
    struct Harm *HARM,
    char *name,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
  char read_string[500];
  void field_arbitrary_harmonic();
  int i;

  sprintf(read_string,"%s_harm",name);
  input_int(read_string,&(HARM->numb),"0");
  sprintf(read_string,"%s_harms",name);
  input_int(read_string,&(HARM->harms),"0,0,19");
  sprintf(read_string,"%s_harm_off",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->off);
  sprintf(read_string,"%s_harm_x1",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->x1);
  sprintf(read_string,"%s_harm_x2",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->x2);
  sprintf(read_string,"%s_harm_z1",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->z1);
  sprintf(read_string,"%s_harm_z2",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->z2);
  sprintf(read_string,"%s_harm_y1",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->y1);
  sprintf(read_string,"%s_harm_y2",name);
  input_std_precision_vector(read_string,HARM->numb,HARM->y2);
  sprintf(read_string,"%s_harm_ovl",name);
  input_char_vector(read_string,HARM->numb,HARM->overlay);
 
  for(i=0;i<HARM->harms;i++) {
      sprintf(read_string,"%s_harm_kx%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->kx[i]);
      sprintf(read_string,"%s_harm_kz%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->kz[i]);
      sprintf(read_string,"%s_harm_ky%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->ky[i]);
      sprintf(read_string,"%s_harm_ka%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->ka[i]);
      sprintf(read_string,"%s_harm_phx%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->phx[i]);
      sprintf(read_string,"%s_harm_phz%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->phz[i]);
      sprintf(read_string,"%s_harm_phy%02d",name,i+1);
      input_std_precision_vector(read_string,HARM->numb,HARM->phy[i]);
  }

  if(parse_and_apply)
	field_arbitrary_harmonic(E,HARM,field,bcbitf,bcmask_on,bcmask_off,level);

  return;  
}

void field_arbitrary_harmonic(
    struct All_variables *E,
    struct Harm *HARM,
    standard_precision *field,
    unsigned int *bcbitf,
    unsigned int bcmask_on,
    unsigned int bcmask_off,
    int level
)
{
    
    standard_precision weight,radius2,weight2;
    standard_precision x1,y1,z1;
    int in1,in2,in3;
    int number,node,l;
    int combine_option;
    int change;
 
    for(node=1;node<=E->mesh.NNO[level];node++) {
      change=0;
      x1=E->SX[level][1][node];
      z1=E->SX[level][2][node];
      y1=(E->mesh.nsd!=3) ? 0.0 : E->SX[level][3][node];
	
	for(number=0;number<HARM->numb;number++) {
	   
	    switch (HARM->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }

	    in1=(x1 >= HARM->x1[number] && x1 <= HARM->x2[number]);
	    in2=(z1 >= HARM->z1[number] && z1 <= HARM->z2[number]);
	    in3=(3!=E->mesh.nsd || y1 >= HARM->y1[number] && y1 <= HARM->y2[number]);
	    
	    if(in1 && in2 && in3) {
		weight = HARM->off[number];
		for(l=0;l<HARM->harms;l++) {
		  weight += HARM->ka[l][number] *
		    cos((HARM->kx[l][number]*x1+HARM->phx[l][number])*M_PI) *
		    cos((HARM->kz[l][number]*z1+HARM->phz[l][number])*M_PI) *
		    ((3!=E->mesh.nsd) ? 1.0 : cos((HARM->ky[l][number]*y1+HARM->phy[l][number])*M_PI)) ;
			change++;
		}

	   	if(field != NULL)
		    switch (HARM->overlay[number]) {
		    case 'R':
			field[node]=weight;
			break;
		    case 'M':
			field[node]*=weight;
			break;
		    case 'A':
			field[node]+=weight;
			break;
		    default:
			fprintf(E->fp,"POLYNOMIAL: %d can't work out how to combine new/old fields\n",number);
			break;
		    }

		if(change && bcbitf!=NULL) {
		    bcbitf[node] = (bcbitf[node] | bcmask_on);
		    bcbitf[node] = (bcbitf[node] & (~bcmask_off));
		}
	    }
	}
    }
    return;  
}


#if defined( No_drand48 )

double drand48(void) 
{
	return(0.1);   /* Just to get things to compile */
	
}



#endif
