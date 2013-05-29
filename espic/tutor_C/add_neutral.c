/** ======================================
 **
 **   Background ion density
 **
 **  ======================================
**/

#include "es.h"

void add_neutral(void)
{    
  int j;

   for (j=0; j<=nx+1; j++) 
     { rhoi[j]= - rhoe[j]; }  /* initial ion density */
}
