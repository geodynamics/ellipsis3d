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


int epsilon[4][4] = {   /* Levi-Cita epsilon */
  {0, 0, 0, 0},
  {0, 1,-1, 1},
  {0,-1, 1,-1},
  {0, 1,-1, 1} };

static standard_precision cost_per_level[MAX_LEVELS]; /* this will accumulate data over the run */
static int total_cycles[MAX_LEVELS];






higher_precision cofactor(   /* n <= 3 in this case */
     higher_precision A[4][4],
     int i,
     int j,
     int n
)

{ 
    int k,l,p,q;
    higher_precision determinant();
    static int been_here = 0;
    static higher_precision B[4][4][4]; /* because of recursive behaviour of det/cofac, need to use
			       new copy of B at each 'n' level of this routine */
 
  if (n>3) printf("Error, no cofactors for matrix more than 3x3\n");

  p=q=1;

  for(k=1;k<=n;k++)
    { if(k==i) continue;
       for(l=1;l<=n;l++)
	  { if (l==j) continue;
	    B[n][p][q]=A[k][l];
	    q++ ;
	  }
     q=1;p++;  
    }
       
  return(epsilon[i][j]*determinant(B[n],n-1));
}

/* Fast (conditional) determinant for 3x3 or 2x2 ... otherwise calls general routine */

higher_precision determinant(
     higher_precision A[4][4],
     int n
)
{ higher_precision gen_determinant();

  switch (n)
    { case 1: 
	return(A[1][1]);
	
      case 2:
	return(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	
      case 3:
	return(A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])-
	       A[1][2]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+
	       A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]));
	
      default:
	return(gen_determinant(A,n));
      }
}

/* recursive function to determine matrix determinant but n<=3 only, for this type of A */

higher_precision gen_determinant(
     higher_precision A[4][4],
     int n
)
{
    higher_precision det;
    higher_precision cofactor();

    int i;
    
    if(n==1) return(A[1][1]); /* need a way to break the recursion */
  
    det=0.0; 
    for(i=1;i<=n;i++)
	det += A[1][i]*cofactor(A,1,i,n);
    
    return(det);
}


