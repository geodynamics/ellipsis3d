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

#include "config.h"

#if HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#if HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#if HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#if HAVE_TIME_H
#include <time.h>
#endif

#include "global_defs.h"


#define GETRUSAGE 0
#define TIMES 1
#define STUPID 2

#if (!defined(TIME_FUNC) && HAVE_GETRUSAGE && HAVE_SYS_RESOURCE_H)
#define TIME_FUNC GETRUSAGE
#endif

#if (!defined(TIME_FUNC) && HAVE_TIMES && HAVE_SYS_TIMES_H)
#define TIME_FUNC TIMES
#endif

#ifndef TIME_FUNC
#define TIME_FUNC STUPID
#endif


/* ===============================================
   Function to return currently elapsed CPU usage
   =============================================== */

standard_precision CPU_time()
     
{
    standard_precision time;
    
    switch (TIME_FUNC) {
        
    case GETRUSAGE:
#if (HAVE_GETRUSAGE && HAVE_SYS_RESOURCE_H)
    {
        struct rusage rusage;
        
        getrusage(RUSAGE_SELF,&rusage);
        time = (standard_precision)(rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec);
    }
#endif
    break;
    
    case TIMES:
#if (HAVE_TIMES && HAVE_SYS_TIMES_H)
    {
        struct tms time_now;
        time_t utime;
        long sometime;
        
        static standard_precision initial_time;
        static int visit = 0;
        
        if (visit==0) {
            sometime=times(&time_now);
            initial_time = (standard_precision) time_now.tms_utime / (standard_precision) CLK_TCK;
            visit++;
        }
        
        sometime=times(&time_now);
        time = (standard_precision) time_now.tms_utime / (standard_precision) CLK_TCK - initial_time;
    }
#endif
    break;
    
    case STUPID: /* stupid, break nothing "timer" */
    default:
    {
        static standard_precision counter;
        counter += 0.0001;
        time = counter;
    }
    break;
    
    } /* switch (TIME_FUNC) */
    
    return time;
}
