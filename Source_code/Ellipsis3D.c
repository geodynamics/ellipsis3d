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

/*#include <sys/time.h>
#include <sys/resource.h>*/

extern int Emergency_stop;

main(
     int argc,
     char **argv
)
{	/* Functions called by main*/
  void read_instructions();
  void solve_derived_velocities();
  void vcopy();

  standard_precision dot();
  standard_precision CPU_time();

  int i,k;
  standard_precision time,initial_time,start_time;
  standard_precision total_work;
 
  struct All_variables E;
 
   /*==================================================
    =                   INITIALIZE                   =
    ================================================== */
  
  E.monitor.solution_cycles=1;
  start_time = time = CPU_time();

  read_instructions(&E,argc,argv);

  fprintf(E.fp,"Input parameters taken from file '%s'\n",argv[1]);

  E.control.keep_going=1;
  E.monitor.cpu_time_on_vp_it = 0.0;
  E.monitor.cpu_time_on_forces = 0.0;
  E.monitor.cpu_time_on_mg_maps = 0.0;

  fprintf(E.fp,"Initialization complete after %g seconds\n\n",CPU_time()-time);
  fflush(E.fp);
  initial_time = CPU_time()-time;

  /*==================================================
    =                    SOLVE                       =
    ================================================== 
  */
  while ( E.control.keep_going   &&  (Emergency_stop == 0) )   {
      time = CPU_time();
      (E.solve_stokes_problem)(&E); /* solve_constrained_flow_iterative */
      fflush(E.fp);

      E.monitor.cpu_time_on_vp_it += CPU_time()-time;
      process_new_velocity(&E,E.monitor.solution_cycles); 
      (E.special_process_new_velocity)(&E,E.monitor.solution_cycles); /* convection output */

       fflush(E.fp);

      (E.next_buoyancy_field)(&E);/*pg timestep*/
      fflush(E.fp);
  
      process_new_buoyancy_field(&E,E.monitor.solution_cycles);  
      (E.special_process_new_buoyancy)(&E,E.monitor.solution_cycles);  /* hook into user supplied stuff */
      fflush(E.fp);

      /* (E.problem_update_node_positions)(&E,E.monitor.solution_cycles);  /* caution ! */
      (E.problem_update_bcs)(&E,E.monitor.solution_cycles);
      fflush(E.fp);

      E.monitor.cpu_time_elapsed = CPU_time();
      fprintf(E.fp,"Cpu time elapsed after %d cycles  = %g seconds\n\n\n",
	      E.monitor.solution_cycles, E.monitor.cpu_time_elapsed); fflush(E.fp);
      E.monitor.solution_cycles++; 

      if(E.monitor.solution_cycles>E.control.print_convergence)
	  E.control.print_convergence=0;

    
  }
  
  if (Emergency_stop) {
    (E.special_process_new_buoyancy)(&E,0); 
  }


  fprintf(E.fp,"Initialization overhead = %f\n",initial_time);
  fprintf(E.fp,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((standard_precision)(E.monitor.solution_cycles-1)));
  fprintf(E.fp,"Average cpu time taken for problem solving per v step =%f\n",
	 (E.monitor.cpu_time_elapsed-E.monitor.cpu_time_on_vp_it)/((standard_precision)(E.monitor.solution_cycles-1)));
  fclose(E.fp);
  fclose(E.fp1);
  fclose(E.fp2);


  total_work=0.0;
  for(i=E.mesh.levmin+1;i<=E.mesh.levmax;i++)
    total_work += (double)(E.control.vsolver_relaxations * (E.mesh.levmax-i+1) ) *
      E.monitor.total_gs_cycles[i]*E.mesh.NEQ[i]*max_node_interaction[E.mesh.nsd]*E.mesh.dof;

  for(i=E.mesh.levmin+1;i<=E.mesh.levmax;i++)
    fprintf(stderr,"Work at level %d was %g [%d] or %g %%\n",i,(double)(E.control.vsolver_relaxations * (E.mesh.levmax-i+1)) *
	    E.monitor.total_gs_cycles[i]*E.mesh.NEQ[i]*max_node_interaction[E.mesh.nsd]*E.mesh.dof,
	    E.monitor.total_gs_cycles[i],
	    100.0 * (E.control.vsolver_relaxations * (E.mesh.levmax-i+1)) *E.monitor.total_gs_cycles[i]*E.mesh.NEQ[i]*
	    max_node_interaction[E.mesh.nsd]*E.mesh.dof/total_work);

  fprintf(stderr,"Total work done by solver = %e\n",total_work);
 
  return(1);  
}  /*  End of element_driver section  */
