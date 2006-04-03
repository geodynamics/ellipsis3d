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

/* Procedure which multiplies the m*n matrix A by the n*p matrix B 
The result is the m*p matrix C*/

void mul_A_B(
  standard_precision A[7][16],
  standard_precision B[16][7],
  standard_precision C[7][7]
)
{
  int i,j,k ;
  for(i=1;i<7;i++){
    for(j=1;j<7;j++){
      C[i][j] = 0.0 ;
      for(k=1;k<16;k++){
	C[i][j] += A[i][k] * B[k][j] ;
      }
    }
  }
  return ;
}

/* Procedure which multiplies the transpose of the m*n matrix A by the m*p matrix B 
The result is the n*p matrix C*/

void mul_At_B(
   standard_precision A[16][7],
   standard_precision B[16][16],
   standard_precision C[7][16]
)
{
  int i,j,k ;

  for(i=1;i<7;i++){
    for(j=1;j<16;j++){
      C[i][j] = 0.0 ;
    }
  }
  for(i=1;i<7;i++){
    for(j=1;j<16;j++){
      for(k=1;k<16;k++){
	C[i][j] += A[k][i] * B[k][j] ;
      }
    }
  }
  return ;
}
