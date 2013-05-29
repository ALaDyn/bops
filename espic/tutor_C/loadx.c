
/**  =======================================
 **
 **   Particle loading for uniform density
 **
     ======================================= **/

#include "es.h"

/*  local variables  */

float dpx;        /* particle spacing */ 
int i;            /* counter */

/*  perturbation test */

float a_perturb = 0.8, x0, x1;

void loadx(void)
{
      if (bc_particle >= 2)
      {

         /**  for reflective bcs to work properly, must place particle
              boundaries half a mesh spacing inside field boundaries **/

         wall_left = dx/2.;
         wall_right = grid_length-dx/2.;
	 plasma_end = wall_right;          /* actually want min(end,wr) */

      }
      else 
      {

              /*  periodic boundaries */
         plasma_start = 0.;
	 plasma_end = grid_length;
         wall_left = 0.;
         wall_right = grid_length;
      }

      xload = plasma_end - plasma_start;   /*  length for particle loading */
      dpx = xload/ne;                   /*  particle spacing */

 /*  pseudo-particle charge normalised to give ncrit=1 (rhoc=-1) */

      qe = -rho0*dpx;

 /*  pseudo-particle mass (need for kinetic energy diagnostic) */

      e_mass = qe / q_over_me;

 /*  set up initial positions according to chosen density profile */

      switch ( profile ) {
      
        default:        /* (0 or 1) */               

	  /* particles */	
	  for (i = 1; i <= ne; i++) {  x[i] = plasma_start + dpx*(i-0.5); }
	  break;

        case 2:

	  /* oscillation test: perturb positions */ 
	  printf("\n Plasma oscillation test:\n - perturbing initial positions\
 with amplitude %f\n", a_perturb );
	  
	  for (i = 1; i <= ne; i++) 
	    { 
	      x[i] = plasma_start + dpx*(i-0.5);
	      x[i] += a0*cos(x[i]);
	    }
	  break;
      }
}
