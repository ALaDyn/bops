
/**  ============================================
 **
 **    Diagnostics - graphics and printed output
 **
 **  ============================================
 **/

#include "es.h"

void diagnostics(void)
{


   /*  write run information to terminal */

   if (itime % iout == 0) 
     {
       printf("timestep: %d\n", itime);
     }

   /*  do graphics snapshots */

   if (itime % igraph == 0) 
     {
       plots();
      }

   /*  write out time-dep. quantities to file */
 
   if (itime % ihist == 0) 
     {
      histories(); 
     }
}
