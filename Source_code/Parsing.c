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




/* Routines which read filenames from the command line and
   then parse the contents as parameters for citcom */


#include "config.h"

#include <math.h>
#include <stdio.h>

#if HAVE_STDLIB_H
#include <stdlib.h>
#endif

#if HAVE_STRING_H
#include <string.h>
#endif

#if HAVE_STRINGS_H
#include <strings.h>
#endif

#include "global_defs.h"





#define MAXLINE		2048	/* max length of line in input file */
#define MAXNAME		128	/* max length of name */
#define MAXVALUE	2048	/* max length of value */
#define MAXFILENAME	128	/* max length of par file name */
#define MAXVECTOR	40	/* max # of elements for unspecified vectors */

/* abbreviations: */
#define AL 		struct arglist
#define PROGNAME	ext_par.progname
#define FLAGS		ext_par.argflags
#define ARGLIST		ext_par.arglist
#define ARGHEAD		ext_par.arghead
#define ARGBUF		ext_par.argbuf
#define NLIST		ext_par.nlist
#define NBUF		ext_par.nbuf
#define LISTMAX		ext_par.listmax
#define BUFMAX		ext_par.bufmax
#define LISTFILE	ext_par.listout

#define LISTINC		128	/* increment size for arglist */
#define BUFINC		2048	/* increment size for argbuf */

struct ext_par		/* global variables for getpar */
{
  char *progname;
  int argflags;
  struct arglist *arglist;
  struct arglist *arghead;
  char *argbuf;
  int nlist;
  int nbuf;
  int listmax;
  int bufmax;
  FILE *listout;
}	ext_par;

struct arglist		/* structure of list set up by setpar */
{
    int argname_offset;
    int argval_offset;
    int hash;
};

int VERBOSE = 0;
int DESCRIBE = 0;
int BEGINNER = 0;

void setup_parser(
     struct All_variables *E,
     int ac,
     char **av
)
{
    void unique_copy_file();
    
    FILE * fp;
    char *pl,*pn,*pv;
    char t1, t2, line[MAXLINE], name[MAXNAME], value[MAXVALUE];
    int i,j,k;
    char command[600];
  
    /* should get file length & cpp &c before any further parsing */

    /* for now, read one filename from the command line, we'll parse that ! */
    
    if (ac < 2)   {
	fprintf(stderr,"Usage: citcom PARAMETERFILE\n");
	exit(10);
    }
  
    if ((fp = fopen(av[1],"r")) == NULL)  {
      fprintf(stderr,"File: %s is unreadable\n",av[1]);
      exit(11);
    }

    strcpy(E->control.input_file_name,av[1]);
 
  /* now the parameter file is open, read into memory */

  while( fgets(line,MAXLINE,fp) != NULL )  {
    pl= line;
    /* loop over entries on each line */
  loop:	
    while(*pl==' ' || *pl=='\t')
      pl++;
    if(*pl=='\0'|| *pl=='\n') 
      continue; /* end of line */
    if(*pl=='#') 
      continue; /* end of interpretable part of line */
    
    /* get name */
    pn= name;
    while(*pl != '=' && *pl != '\0' && *pl != ' '
	  && *pl != '\n'		/* FIX by Glenn Nelson */
	  && *pl != '\t') 
      *pn++ = *pl++;
    *pn = '\0';
    if(*pl == '=') pl++;
      
    /* get value */
    *value= '\0';
    pv= value;
    if(*pl=='"' || *pl=='\'')
      t1= t2= *pl++; 
    else
    { t1= ' ';
    t2= '\t';
    }
    while(*pl!=t1 && *pl!=t2 &&
	  *pl!='\0' && *pl!='\n') *pv++= *pl++;
    *pv= '\0';
    if(*pl=='"' || *pl=='\'')
      pl++;
    add_to_parameter_list(name,value);
    
    goto loop;
  }

  fclose(fp);

  ARGHEAD= ARGLIST;

  /* Now we can use our routines to check & set their own flags ! */

  input_boolean("VERBOSE",&i,"off");
  input_boolean("DESCRIBE",&j,"off");
  input_boolean("BEGINNER",&k,"off");
  VERBOSE=i;
  DESCRIBE=j;
  BEGINNER=k;
  
}

void shutdown_parser(
     struct All_variables *E
)
{
	if(ARGLIST != NULL)
	  free(ARGLIST);
	if(ARGBUF  != NULL) 
	  free(ARGBUF);
	ARGBUF=  NULL;
	ARGLIST= NULL;
}


add_to_parameter_list(	/* add an entry to arglist, expanding memory */
     register char *name,
     register char *value	/* if necessary */
)
{
  struct arglist *alptr;
  int len;
  register char *ptr;
  int oldlistmax;
  int oldbufmax;


  /* check arglist memory */
  if(NLIST >= LISTMAX)  {
    if(ARGLIST == NULL) {
      LISTMAX += LISTINC;
      ARGLIST= (AL *)Malloc0(LISTMAX * sizeof(AL));
    }
    else {	
      oldlistmax = LISTMAX;  
      LISTMAX += LISTINC;
      Realloc0(&(ARGLIST),oldlistmax,LISTMAX,sizeof(AL));
    }
  }
  /* check argbuf memory */
  len= strlen(name) + strlen(value) + 2; /* +2 for terminating nulls */
  if(NBUF+len >= BUFMAX) {
    if(ARGBUF == NULL) {
      BUFMAX += BUFINC;
      ARGBUF= (char *)Malloc0(BUFMAX);
    }
      else {
	oldbufmax = BUFMAX;
	BUFMAX += BUFINC;
	Realloc0(&(ARGBUF),oldbufmax,BUFMAX,sizeof(char));
    }
  }
  if(ARGBUF == NULL || ARGLIST == NULL)
   fprintf(stderr,"cannot allocate memory\n");

  /* add name */
  alptr= ARGLIST + NLIST;
  alptr->hash= compute_parameter_hash_table(name);
  alptr->argname_offset = NBUF;
  ptr= ARGBUF + NBUF;
  do 
    *ptr++ = *name; 
  while(*name++);
  
  /* add value */
  NBUF += len;
  alptr->argval_offset= ptr - ARGBUF;
  do
    *ptr++ = *value;
  while(*value++);
  NLIST++;

  return(1);
}

compute_parameter_hash_table(
     register char *s
)
{ register int h;
  
  h= s[0];
  if(s[1])
    h |= (s[1])<<8;
  else
    return(h);
  if(s[2])
    h |= (s[2])<<16;
  else 
    return(h);
  if(s[3])
    h |= (s[3])<<24;
  return(h);
}

int input_int(
     char *name,
     int *value,
     char *interpret
)
{
    int interpret_control_string();
    struct arglist *alptr; 
    int h, found;
    char  *str;
    
  int exists,essential;
  double Default,minvalue,maxvalue;

  if(DESCRIBE)
    fprintf(stderr,"input_int: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
  exists = interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
 
  if(Default != 1.01e32)
    *value = (int)Default;
 
  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      sscanf(str,"%d",value);
      found=1;
      break;
    } 

  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
  if((minvalue!=1.01e32) && (*value < (int) minvalue))
     { *value = (int) minvalue;
     }
  if((maxvalue!=1.01e32) && (*value > (int) maxvalue))
    {  *value = (int) maxvalue;
    }

  if(VERBOSE) {
    if (found)
      fprintf(stderr,"%25s: (int) = %d \n",name,*value); 
    else
      if (Default != 1.01e32)
	fprintf(stderr,"%25s: (int) = not found (%d) \n",name,(int)(Default)); 
      else {
	fprintf(stderr,"%25s: (int) = not found (no default) \n",name); 
	if(BEGINNER) {
	  fprintf(stderr,"\t\t Previously set value gives ...");
	  fprintf(stderr,"%d\n",*value);
	}
      } 
  }
  return(found);
}

int input_string(  /* in the case of a string default=NULL forces input */
     char *name,
     char *value,
     char *Default
)
{ 
    char *sptr;
  struct arglist *alptr; 
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
  int essential;

 
  if(DESCRIBE)
    fprintf(stderr,"input_string: searching for '%s' with default '%s'\n",
	    name,(Default == NULL) ? "no default" : Default);
 
  h=compute_parameter_hash_table(name);
  essential=found=0;

    
    if (Default != NULL)   /* Cannot use "Essential" as this is a valid input */
	strcpy(value,Default);  
    else
	essential=1;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      strcpy(value,str);
      found=1;
      break;
    } 
  
  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
 
  if(VERBOSE)
    fprintf(stderr,"%25s: (string) = %s (%s)\n",name,
	    (found ? value : "not found"),
	    (Default != NULL ?  Default : "no default")); 

  return(found);
}

int input_boolean(  /* supports name=on/off too */
     char *name,
     int *value,
     char *interpret
)
{ char *sptr;
  struct arglist *alptr; 
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  int essential;
  double Default,minvalue,maxvalue;

  if(DESCRIBE)
    fprintf(stderr,"input_boolean: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
 
  interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
  
  if(Default != 1.01e32)
    *value = (int)(Default);
 
  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
 
 
  if(!found)
    { if(VERBOSE)
	if (Default != 1.01e32)
	  fprintf(stderr,"%25s: (boolean int) = not found (%d) \n",name,(int)(Default)); 
	else
	 { fprintf(stderr,"%25s: (boolean int) = not found (no default) \n",name); 
	   if(BEGINNER)
	     { fprintf(stderr,"\t\t Previously set value gives ...");
	       fprintf(stderr,"%d\n",*value);
	     }
	 }
	 
      return(0);
    }
 
  if((strstr(str,"on")!=NULL) || (strstr(str,"ON")!=NULL))
    *value=1;
  else if ((strstr(str,"off") != NULL) || (strstr(str,"OFF")!=NULL))
    *value=0;
  else /* assume some numerical value */
    *value=atoi(str);

  if(VERBOSE)
    fprintf(stderr,"%25s: (boolean int) = %d \n",name,*value); 
  
  return(found);
}

int input_std_precision(
			char *name,
			standard_precision *value,
			char *interpret
			)
{ 
  char *sptr;
  struct arglist *alptr;
  
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
  int exists,essential;
  double Default,minvalue,maxvalue,readvalue;

 
  if(DESCRIBE)
    fprintf(stderr,"input_std_precision: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);

  exists=interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
 
  if(Default != 1.01e32)
    *value = (standard_precision) Default;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--) {
    if(alptr->hash != h)
      continue;
    if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
      continue;
    str= ARGBUF + alptr->argval_offset;
    
    sscanf(str,"%lf",&readvalue);
    *value = (standard_precision) readvalue;
    found=1;
    break;
  }
 
  if(essential && !found) {
    fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
    exit(12);
  }

  if((minvalue!=1.01e32) && (*value < (standard_precision) minvalue))
    *value = (standard_precision) minvalue;
  if((maxvalue!=1.01e32) && (*value > (standard_precision) maxvalue))
    *value = (standard_precision) maxvalue;

  if(VERBOSE) {
    if (found)
      fprintf(stderr,"%25s: (standard_precision) = %f \n",name,*value); 
    else
      if (Default != 1.01e32)
	fprintf(stderr,"%25s: (standard_precision) = not found (%f) \n",name,Default); 
      else {
	fprintf(stderr,"%25s: (standard_precision) = not found (no default) \n",name); 
	if(BEGINNER)   {
	  fprintf(stderr,"\t\t Previously set value gives ...");
	  fprintf(stderr,"%g\n",*value);
	}
      }
  }
  return(found);
}
  
int input_hgr_precision(
     char *name,
     higher_precision *value,
     char *interpret
)
{ char *sptr;
  struct arglist *alptr;
  double read_value;
  
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  int exists,essential;
  double Default,minvalue,maxvalue;


  if(DESCRIBE)
   fprintf(stderr,"input_double: searching for '%s' with default/range '%s'\n",
	   name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
 
  exists=interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
 
  if(Default != 1.01e32)
    *value = Default;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      sscanf(str,"%lf",&read_value);
      *value = (higher_precision) read_value;
      found=1;
      break;
    } 
 
  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
  if((minvalue!=1.01e32) && (*value <  minvalue))
    *value =  minvalue;
  if((maxvalue!=1.01e32) && (*value >  maxvalue))
    *value =  maxvalue;

  if(VERBOSE)
   { if (found)
       fprintf(stderr,"%25s: (higher_precision) = %g \n",name,*value); 
     else
       if (Default != 1.01e32)
	  fprintf(stderr,"%25s: (higher_precision) = not found (%g) \n",name,Default); 
       else
	  { fprintf(stderr,"%25s: (higher_precision) = not found (no default)\n",name); 
	    if(BEGINNER)
	       { fprintf(stderr,"\t\t Previously set value gives ...");
		 fprintf(stderr,"%g\n",*value);
	       }
	  }
   }
  

  return(found);
}


int input_int_vector(
		     char *name,
		     int number,
		     int *value, /* comma-separated list of ints */
		     int default_value
)
{ 
  char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if(DESCRIBE)
    fprintf(stderr,"input_int_vector: searching for %s (%d times)\n",name,number);

  for(i=0;i<number;i++)
    value[i] = default_value;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)
    { sprintf(control_string,"");
      for(i=0;i<h;i++)
	strcat(control_string,"%*f,");
      strcat(control_string,"%d");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (vector) = %s\n",name,str); 

  return(found);
}

int input_char_vector(
     char *name,
     int number,
     char *value /* comma-separated list of ints */
)
{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if(DESCRIBE)
    fprintf(stderr,"input_char_vector: searching for %s (%d times)\n",name,number);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)
    { sprintf(control_string,"");
      for(i=0;i<h;i++)
	strcat(control_string,"%*c,");
      strcat(control_string,"%c");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (vector) = %s\n",name,str); 

  return(found);
}

int input_std_precision_vector(
     char *name,
     int number,
     standard_precision *value, /* comma-separated list of standard_precisions */
     standard_precision default_value
)
{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
  double readvalue;
  
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if(0==number)
      return(0);

  if(DESCRIBE)
    fprintf(stderr,"input_std_precision_vector: searching for %s (%d times)\n",name,number);

 for(i=0;i<number;i++)
    value[i] = default_value;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)    {
    sprintf(control_string,"");
    for(i=0;i<h;i++)
      strcat(control_string,"%*f,");
    strcat(control_string,"%lf");
    if(sscanf(str,control_string,&readvalue) == 1)
      value[h] = (standard_precision) readvalue;
  }

  if(VERBOSE)
   fprintf(stderr,"%25s: (standard_precision vector) = %s\n",name,str); 

  return(found);
}

int input_hgr_precision_vector(
     char *name,
     int number,
     higher_precision *value, /* comma-separated list of double s */
     higher_precision default_value
)
{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
  double readvalue;
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
 
  if(DESCRIBE)
    fprintf(stderr,"input_hgr_precision_vector: searching for %s (%d times)\n",name,number);
 
  for(i=0;i<number;i++)
    value[i] = default_value;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 

  if(!found) return(0);

 /* now interpret vector */
  
  for(h=0;h<number;h++)  {
    sprintf(control_string,"");
    for(i=0;i<h;i++)
      strcat(control_string,"%*f,");
    strcat(control_string,"%lf");
    sscanf(str,control_string,&(readvalue));
    value[h] = readvalue;
  }
  

  if(VERBOSE)
   fprintf(stderr,"%25s: (higher_precision vector) = %s\n",name,str); 

  return(found);
}

/* =================================================== */

int interpret_control_string(
			     char *interpret,
			     int *essential,
			     double *Default,
			     double *minvalue,
			     double *maxvalue
			     )
{ 
  char *substring;
  char interpcopy[10000];

  *Default=*maxvalue=*minvalue=1.01e32;
  *essential=0;

  strcpy(interpcopy,interpret);
     
  if ((substring=(char *)strtok(interpcopy,",")) == NULL)
    return(0);  /* nothing to interpret */
    
  if (strstr(substring,"essential")!=NULL)
    *essential=1; /* no default possible, must read a value */
  else
    if (strstr(substring,"nodefault")==NULL) {
      if((strstr(substring,"on")!=NULL) || (strstr(substring,"ON")!=NULL))
	*Default = 1.0;
      else 
	if ((strstr(substring,"off") != NULL) || (strstr(substring,"OFF")!=NULL))
	  *Default = 0.0; 
	else
	  sscanf(substring,"%lf",Default);  /* read number as a default value */
    }


  if (!(substring=(char *)strtok(NULL,","))) /* minvalue */
    { /* no minimum, no maximum */
      return(1);
    }

  if (strstr(substring,"nomin")==NULL)  { 
      sscanf(substring,"%lf",minvalue);
    }
  
  if (!(substring=(char *)strtok(NULL,","))) /* maxvalue */
    { /* no maximum */
      if (DESCRIBE)
	fprintf(stderr,"minimum but no maximum\n");
      return(1);
    }

  if (strstr(substring,"nomax")==NULL)  {
    sscanf(substring,"%lf",maxvalue);
  }
  return(1);
}
