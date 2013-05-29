/** ======================================
 **
 **   Electron density
 **
 **  ======================================
**/

#include "es.h"

void density(void)
{
   float re, xa, f1, f2;
   int i, j, j1, j2, iwl, iwr;

   re = qe/dx;     /*  charge weighting factor */
     
   for (j=0; j<=nx+1; j++) 
     { rhoe[j]=0.0; }  /* zero density array */

   /*  gather density onto grid from particles */

   for (i=1; i<=ne; i++)
     {
      xa = x[i]/dx;
      j1 = xa;
      j2 = j1 + 1;
      f2 = xa - j1;
      f1 = 1.0 - f2;
      rhoe[j1] = rhoe[j1] + re*f1;
      rhoe[j2] = rhoe[j2] + re*f2;
     }


   if (bc_field == 1)
     {
       /*   periodic boundaries */ 
      rhoe[0] += rhoe[nx];
      rhoe[nx] = rhoe[0];
     }

   else     /* bc_field =2 */
     {
     /*   reflective - 1st and last (ghost) cells folded back onto physical grid */
      iwl = wall_left/dx;
      rhoe[iwl+1] += rhoe[iwl];
      rhoe[iwl] = 0.0;	

      iwr = wall_right/dx;
      rhoe[iwr] += rhoe[iwr+1];
      rhoe[iwr+1] = rhoe[iwr];
     }     

}
