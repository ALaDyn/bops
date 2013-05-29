/**  ======================================
 **
 **  Set up particle velocity distribution
 **
 ** ====================================== 
 **/


#include <math.h>
#include <time.h>

#include "es.h"

void loadv(void)
{

/* local working variables */

  float vm, rs, theta, temp1;
  int i, iseed, idum1, n1, j, time_in_seconds;


/* random number seeds */
  iseed = 76523; idum1 = 137;


/*       Rand. particles with clock time */

/*   time_in_seconds = time(time_in_seconds); */
/*   srandom(time_in_seconds); */

  for (i = 1; i<=ne; i++)
   {
      /*  inverted 2v-distribution - amplitude */

    vm = vte*sqrt( (-2.*log((i-0.5)/ne)) );

    rs = (float) rand() / random_max;        /* random angle */
    theta = 2*pi*rs;
  
    vx[i] = vm*sin(theta); /* x-component of v */
   }



/*   exclude odd one out */

      n1 = ne;
      if (ne%2 !=0 ) {
	n1--;
	  }

/**  scramble particle indicies to remove correlations
       between x and vx  **/

      for (i = 1; i<=n1; i++)
	{
	 j = n1* (float) rand()/random_max + 1;

	 temp1 = vx[i];        /* switch i,j */
	 vx[i] = vx[j];
	 vx[j] = temp1;
	}				  
}
