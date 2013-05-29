/**  ============================================
 **
 **   Time-histories
 **
 ** ============================================
 **/

#include <stdio.h>
#include <math.h>
#include "es.h"

FILE *history_file;     /* file for writing out time histories */            

void histories(void)
{

  int i;
  float ukin;

  /*   kinetic energy */
 
  ukin = 0.0;

  for ( i=1; i<=ne; i++ )
    {
      ukin += 0.5*e_mass*vx[i]*vx[i];
    }

  /*  write energies out to file */

   fprintf( history_file, "%f  %f\n", itime*dt, ukin );

}
