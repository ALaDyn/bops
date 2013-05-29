
/**  ============================================
 **
 **   Direct field integration - slab geometry
 **
 **  ============================================ 
 **/

#include "es.h"


void field(void)
{
   int j;
   float edc;
   float rhot[501];     /*  net density */


  /*  Add neutral background to get net density */
 
   for (j=1; j<=nx+1; j++)
     { rhot[j] = rhoe[j] + rho0; }
      
  /*   integrate div.E=rho directly (trapezium approx) */

  /*   end point - ex=0 mirror at right wall */


  /*  integrate from right to left */

   Ex[nx+1] = 0.0;
   edc = 0.0;

   for (j=nx; j>=1; j--)
      {
	Ex[j] = Ex[j+1] - 0.5*( rhot[j] + rhot[j+1] )*dx;
        edc = edc + Ex[j];
      }


   if (bc_field == 1)
      {
	/*  periodic fields:  subtract off DC component */ 
	/*  -- need this for consistency with charge conservation */
      
	for (j=1; j<=nx; j++)
	  { Ex[j] -= edc/nx; }
          
	/*  end points periodic */
	Ex[0] = Ex[nx];
      }

   /*   potential - same again */
   /*   - integration const. is arbitrary - phi not used for anything at present */

      phi[nx+1] = 0.0;

      for (j=nx; j>=1; j--)
	{
         phi[j] = phi[j+1] + 0.5*( Ex[j] + Ex[j+1] )*dx;
	}
}
 
