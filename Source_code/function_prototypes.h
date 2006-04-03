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


/* Function prototypes for citcom */

void get_element_viscosity(
  struct All_variables *E,
  int el,
  int level,
  int non_linear_repetition,
  standard_precision *U1,
  standard_precision *U2,
  standard_precision *U3,
  standard_precision *Pres
 );

void get_system_viscosity(
  struct All_variables *E
  );

void pg_element_residual(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision **vel,
  standard_precision *TT1,
  standard_precision *TT1dot,
  struct SOURCES QQ1,
  higher_precision Eres1[9],
  standard_precision diff1,
  standard_precision *BCC,
  unsigned int *FLAGS1,
  int OFFSIDE_MASK1,
  int FLUX_MASK1
  );

void pg_element_residual2(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision **vel,
  standard_precision *TT1,
  standard_precision *TT2,
  standard_precision *TT1dot,
  standard_precision *TT2dot,
  struct SOURCES QQ1,
  struct SOURCES QQ2,
  higher_precision Eres1[9],
  higher_precision Eres2[9],
  standard_precision diff1,
  standard_precision diff2,
  standard_precision *BCC,
/*  standard_precision *BC2,*/
  unsigned int *FLAGS1,
  unsigned int *FLAGS2,
  int OFFSIDE_MASK1,
  int OFFSIDE_MASK2,
  int FLUX_MASK1,
  int FLUX_MASK2
  );

void pg_element_residual3(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision **vel,
  standard_precision *TT1,
  standard_precision *TT2,
  standard_precision *TT3,
  standard_precision *TT1dot,
  standard_precision *TT2dot,
   standard_precision *TT3dot,
  struct SOURCES QQ1,
  struct SOURCES QQ2,
  struct SOURCES QQ3,
  higher_precision Eres1[9],
  higher_precision Eres2[9],
  higher_precision Eres3[9],
  standard_precision diff1,
  standard_precision diff2,
  standard_precision diff3,
  standard_precision *BCC,
/*  standard_precision *BC2,
  standard_precision *BC3,*/
  unsigned int *FLAGS1,
  unsigned int *FLAGS2,
  unsigned int *FLAGS3,
  int OFFSIDE_MASK1,
  int OFFSIDE_MASK2,
  int OFFSIDE_MASK3,
  int FLUX_MASK1,
  int FLUX_MASK2,
  int FLUX_MASK3
  );

void pg_solver(
  struct All_variables *E,
  standard_precision *T,
  standard_precision *Tdot,
  standard_precision *DTdot,
  standard_precision **V,
  struct SOURCES Q0,
  standard_precision diff,
  int bc,
  standard_precision *TBC,
  unsigned int *FLAGS,
  unsigned int OFFSIDE_MASK,
  unsigned int FLUX_MASK );

void pg_solver2(
  struct All_variables *E,
  standard_precision *T1,
  standard_precision *T2,
  standard_precision *T1dot,
  standard_precision *T2dot,
  standard_precision *DT1dot,
  standard_precision *DT2dot,
  standard_precision **V,
  struct SOURCES Q1,
  struct SOURCES Q2,
  standard_precision diff1,
  standard_precision diff2,
  int bc1,
  int bc2,
/*  standard_precision *TBC1,
  standard_precision *TBC2,*/
  unsigned int *FLAGS1,
  unsigned int *FLAGS2,
  unsigned int OFFSIDE_MASK1,
  unsigned int OFFSIDE_MASK2,
  unsigned int FLUX_MASK1,
  unsigned int FLUX_MASK2
);

void pg_solver3(
  struct All_variables *E,
  standard_precision *T1,
  standard_precision *T2,
  standard_precision *T3,
  standard_precision *T1dot,
  standard_precision *T2dot,
  standard_precision *T3dot,
  standard_precision *DT1dot,
  standard_precision *DT2dot,
  standard_precision *DT3dot,
  standard_precision **V,
  struct SOURCES Q1,
  struct SOURCES Q2,
  struct SOURCES Q3,
  standard_precision diff1,
  standard_precision diff2,
  standard_precision diff3,
  int bc1,
  int bc2,
  int bc3,
/*  standard_precision *TBC1,
  standard_precision *TBC2,
  standard_precision *TBC3,*/
  unsigned int *FLAGS1,
  unsigned int *FLAGS2,
  unsigned int *FLAGS3,
  unsigned int OFFSIDE_MASK1,
  unsigned int OFFSIDE_MASK2,
  unsigned int OFFSIDE_MASK3,
  unsigned int FLUX_MASK1,
  unsigned int FLUX_MASK2,
  unsigned int FLUX_MASK3
);


void pic_solver(
  struct All_variables *E,
  standard_precision *T,
  standard_precision *Tdot,
  standard_precision *DTdot,
  standard_precision **V,
  struct SOURCES Q0,
  standard_precision diff,
  int bc,
  standard_precision *TBC,
  unsigned int *FLAGS,
  unsigned int OFFSIDE_MASK,
  unsigned int FLUX_MASK
);
void diff_element_residual(
  struct All_variables *E,
  int el,
  struct Shape_function_dx GNx,
  struct Shape_function_dA dOmega,
  standard_precision *TT1,
  standard_precision *TT1dot,
  struct SOURCES QQ1,
  higher_precision Eres1[9],
  standard_precision diff1,
  standard_precision *BCC,
  unsigned int *FLAGS1,
  int OFFSIDE_MASK1,
  int FLUX_MASK1
);
