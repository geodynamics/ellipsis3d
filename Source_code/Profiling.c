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



/* Profiling functions .... return elapsed CPU time etc. 
   These functions seem the most likely to get broken by
   different architectures/operating systems &c */

#if defined(_UNICOS) || defined(__hpux) || (defined(__SVR4) && !defined(__sunos__))
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#define TIMES_STYLE_TIME
#elif defined(__osf__) || defined(__aix__) || defined(__sunos__) || defined(__sgi)  || defined(__uxp__)
#include <sys/time.h> 
#include <sys/resource.h>
#define RUSAGE_STYLE_TIME
#endif


#include "global_defs.h"

/* ===============================================
   Function to return currently elapsed CPU usage
   =============================================== */

standard_precision CPU_time()
     
{ 
#if defined(RUSAGE_STYLE_TIME)
  struct rusage rusage;
  double time;

  getrusage(RUSAGE_SELF,&rusage);
  time = rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec ;
#elif defined(TIMES_STYLE_TIME)
  struct tms time_now;
  time_t utime;
  long sometime;
  
  standard_precision time;
  static standard_precision initial_time;
  static int visit = 0;

  if (visit==0)
  { sometime=times(&time_now);
    initial_time = (standard_precision) time_now.tms_utime / (standard_precision) CLK_TCK;
    visit++;
  }

  sometime=times(&time_now);
  time = (standard_precision) time_now.tms_utime / (standard_precision) CLK_TCK - initial_time;

#else  /* stupid, break nothing "timer" */
  static standard_precision time;
  
  time += 0.0001;

#endif

   return((standard_precision)time);

 }
