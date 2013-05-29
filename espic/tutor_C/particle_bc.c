/**  ===================================
 **
 **   Particle boundary conditions
 **
 ** ===================================  
 **/

#include "es.h"

void particle_bc(void)
{
   
  int iseed1 = 28631;        /* random number seed */
  int iseed2 = 1631;         /* random number seed */
  int i, i_inject;

  float xl = grid_length;     /*  shorthand */  
  float wl = wall_left;
  float wr = wall_right;
  float vm, theta, rs;

  /**  loop over all particles to see if any have
       left simulation region: if so, we put them back again
       according to the switch 'bc_particle' **/

  for ( i=1; i<=ne; i++)
    {


      switch (bc_particle) {

      case 1:              /*  periodic */
	   
	    if ( x[i] < 0.0 ) { x[i] += xl; }
	    if ( x[i] >= xl ) { x[i] -= xl; }
	    break;

      case 2:          /* reflective at x = wall_left and x = wall_right */

	    if ( x[i] <= wl )
	      {
	       x[i] = 2*wl - x[i];
	       vx[i] = - vx[i];
	      }
 

	    if ( x[i] >= wr )
	      {
	       x[i] = 2*wr - x[i];
	       vx[i] = - vx[i];
	      }
	    break;

      case 3:

	/**   reflect at LH boundary
	      and re-inject with thermal velocity 
	      at RH (high density) boundary **/


	    if ( x[i] <= wl )
	      {
		x[i] = 2*wl - x[i]; 
		vx[i] = -vx[i];
	      }
       

	    if ( x[i] >= wr )
	      {
		x[i] = 2*wr - x[i];

		/*  find new random velocity in cold Maxwellian */

		rs = (float) rand() / random_max;   /* random (0-1) */

		i_inject = ne*rs + 1;
		vm = vte*sqrt(-2.*log((i_inject-0.5)/ne));

		rs = (float) rand() / random_max;   /* random (0-1) */
		theta = pi*rs;   /* 1st quadrant */
		vx[i] = -vm*sin(theta);
	      }	
	    break;

      default:           /*  default to periodic */
	   
	    if ( x[i] < 0.0 ) { x[i] += xl; }
	    if ( x[i] >= xl ) { x[i] -= xl; }
	    break;
      }

    }

}
