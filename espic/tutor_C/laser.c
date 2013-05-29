/**  ===================================
 **
 **         laser.c
 **
 **   Laser routine
 **
 ** ===================================  
 **/
#include "es.h"			/* common variables */

float xj, earg, xlend;
int j, jend;

void laser()
{


   for (j=1; j<=nx; j++)
     {
       if (ilas == 1)             /* homogeneous laser field */	
	 {
	   xlend = plasma_start + 1.0;  /* penetrates 1 skin depth */
	   jend = xlend/dx;
	   if (j <= jend)
	     {
	       Elaser[j] = a0*sin(w0*itime*dt); 
	     }
	   else
	     { Elaser[j] = 0.; }
 	 }

       else if (ilas == 2)        /* ponderomotive force */
	 {
	   xj = dx*j;
	   earg = (xj - xcentre)*(xj - xcentre)/tdel/tdel;
	   Elaser[j] = -(xj - xcentre)*a0*a0/2./tdel/tdel*exp(-earg);
	 }

       else                       /* no laser */
	 { Elaser[j]=0.;}
     }

/* advance centroid of laser */

  xcentre = xcentre + dt;
}
