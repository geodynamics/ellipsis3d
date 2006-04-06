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



#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

/*=============================================================
  Functions to allocate/remove space for variable sized vector.
  =============================================================  */

standard_precision lpdot(
  struct All_variables *E,
  standard_precision *A,
  standard_precision *B,
  int level
)
    
{   standard_precision prod;
    int e;
    const int n=E->mesh.NPNO[level];

    prod = 0.0;
    for(e=1;e<=n;e++)
      prod += A[e] * B[e] ;
    return(prod); 
}

standard_precision lpselfdot(
  struct All_variables *E,
  standard_precision *A,
  int level
)
{   standard_precision prod;
    int e;
    const int n=E->mesh.NPNO[level];

    prod = 0.0;
    for(e=1;e<=n;e++)
      prod += A[e] * A[e] ;
    return(prod);  }

/* Beware of the alias if A=B on vector machines, use vselfdot instead  */

standard_precision vsdot(
     struct All_variables *E,
     standard_precision *A,
     standard_precision *B,
     int level
)
{   
    higher_precision prod,mprod[1];
    int i,incx=1;
    
    char trans='N';
    higher_precision alpha=1.0;
    
    const int n = E->mesh.NEQ[level];

    prod = 0.0;
    for(i=0;i<n;i++)
	prod += A[i] * B[i] ;

    return(prod);  }


standard_precision vsselfdot(
     struct All_variables *E,
     standard_precision *A,
     int level
)
{ 
    standard_precision prod;
    int i,n;

  n = E->mesh.NEQ[level];

  prod = 0.0;
  for(i=0;i<n;i++) prod += A[i] * A[i] ;
  
  return(prod);
}


standard_precision vsselfdot3(
     struct All_variables *E,
     standard_precision *A1,
     standard_precision *A2,
     standard_precision *A3,
     int level
)
{ 
    standard_precision prod;
    int i;

  const int n = E->mesh.NNO[level];
  const int dofs = E->mesh.dof;

  prod = 0.0;

  switch(dofs) {
  case 2:
     for(i=1;i<=n;i++) {
       if(E->NODE[level][i] & ( OFFSIDE  )) continue;
       prod += A1[i] * A1[i] + A2[i] * A2[i] ;
     }
    break ;
  case 3:
    for(i=1;i<=n;i++) {
      if(E->NODE[level][i] & ( OFFSIDE  )) continue;
      prod += A1[i] * A1[i] + A2[i] * A2[i] +  A3[i] * A3[i] ;
    }
    break;
  }
  return(prod);
}

standard_precision vsselfdot6(
     struct All_variables *E,
     standard_precision *A1,
     standard_precision *A2,
     standard_precision *A3,
     standard_precision *A4,
     standard_precision *A5,
     standard_precision *A6,
     int level
)
{ 
    standard_precision prod;
    int i;

  const int n = E->mesh.NNO[level];
  const int dofs = E->mesh.dof;

  prod = 0.0;

  switch(dofs) {
  case 2:
     for(i=1;i<=n;i++) {
       if(E->NODE[level][i] & ( OFFSIDE  )) continue;
       prod += A1[i] * A1[i] + A2[i] * A2[i] ;
     }
    break ;
  case 3:
    for(i=1;i<=n;i++) {
      if(E->NODE[level][i] & ( OFFSIDE  )) continue;
      prod += A1[i] * A1[i] + A2[i] * A2[i] +  A3[i] * A3[i] ;
    }
    break ;
  case 6:
    for(i=1;i<=n;i++) {
      if(E->NODE[level][i] & ( OFFSIDE  )) continue;
      prod += A1[i] * A1[i] + A2[i] * A2[i] +  A3[i] * A3[i] + A4[i] * A4[i] + A5[i] * A5[i] +  A6[i] * A6[i] ;
    }
    break ;
  }
  return(prod);
}

standard_precision fselfdot(
      standard_precision *A,
     int n1,
     int n2
)

{  standard_precision prod;
   int i;
 
   prod = 0.0;
   for(i=n1;i<=n2;i++)
     prod += A[i] * A[i] ;

   return(prod);
}


higher_precision hdot(
     higher_precision *A,
     higher_precision *B,
     int n1,
     int n2
)

{  higher_precision prod;
   int i;
   
   prod = 0.0;
   for(i=n1;i<=n2;i++)
     prod += A[i] * B[i] ;
  
   return(prod);  }

higher_precision hselfdot(
  higher_precision *A,
  int n1,
  int n2
)
{  higher_precision prod;
   int i;
 
   prod = 0.0;
   for(i=n1;i<=n2;i++)
     prod += A[i] * A[i] ;

   return(prod);
}

higher_precision ddot(
     higher_precision *A,
     higher_precision *B,
     int n1,
     int n2
)
{ 
    higher_precision prod;
    int i;
    
    prod = 0.0;
    for(i=n1;i<=n2;i++)
	prod += A[i] * B[i] ;
  
    return(prod);  
}

higher_precision dselfdot(
     higher_precision *A,
     int n1,
     int n2
)
{
    higher_precision prod;
    int i;
 
    prod = 0.0;
    for(i=n1;i<=n2;i++)
	prod += A[i] * A[i] ;
  
   return(prod);
}

void vcopy(
     standard_precision *A,
     standard_precision *B,
     int a,
     int b
)
{   int i;

    for(i=a;i<=b;i++)
      A[i] = B[i]; 

    return; }
