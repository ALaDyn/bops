/**  ============================================
 **
 **   Time-histories
 **
 ** ============================================
 **/


#include "es.h"

FILE *history_file;     /* file for writing out time histories */            

void histories(void)
{

  int i;
  float ukin, upot, utot;

  /*   kinetic energy */
 
  ukin = 0.0;
  for ( i=1; i<=ne; i++ )  { ukin += 0.5*e_mass*vx[i]*vx[i]; }

  /*   potential energy */
 
  upot = 0.0;
  for ( i=1; i<=nx; i++ )   { upot += 0.5*dx*Ex[i]*Ex[i]; }

  /*  total energy */

   utot = upot + ukin;

  /*  write energies out to file */

   fprintf( history_file, "%f  %f  %f  %f\n", itime*dt, ukin, upot, utot );

}
