/**  ===================================
 **
 **         push.c
 **
 **   Electrostatic particle pusher
 **
 ** ===================================  
 **/

#include "es.h"

void push(void)
{

  /* locals */

   float xa;            /* fractional particle position */
   float  b1, b2;       /* linear weights */
   float  exi;          /* interpolated e-field */
   int i, j1, j2;          /* grid indices */

   for( i=1; i<=ne; i++)
     {

       /*   interpolate field Ex from grid to particle */

      xa = x[i]/dx;
      j1 = xa;
      j2 = j1 + 1;
      b2 = xa - j1;
      b1 = 1.0 - b2;
      exi = b1*Ex[j1] + b2*Ex[j2];
      
      /*  update velocities */

      vx[i] = vx[i] + q_over_me*dt*exi; 
     }



   /*  update positions (2nd half of leap-frog) */

   for( i=1; i<=ne; i++)
     {
      x[i] = x[i] + dt*vx[i];
     }

}
