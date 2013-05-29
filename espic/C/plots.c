
/**  =========================================
 **
 **   Graphical snapshots
 **
 **  =========================================
 **/

#include "es.h"


 
FILE  *plot_file;     /* file for dumping plot data */            

int isnap = 0;                /* counts number of calls to routine 'plots' */

void plots(void)
{

   float xgrid[NX_MAX];             /* grid work array */

   char cfile[40];               /* plot filename string */

   int i,j;

   for (j=0; j<=nx; j++)
      {
	xgrid[j] = j*dx;   /* set up x-values for grid plots */
      }



/*   electron density */

   sprintf(cfile,"density%d.data", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     {
       fprintf( plot_file, "%f  %f\n" ,xgrid[j],-rhoe[j] );
     }

   fclose(plot_file);



/*   electrostatic field */

   sprintf(cfile,"field%d.data", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     {
       fprintf( plot_file, "%f  %f\n" ,xgrid[j],Ex[j] );
     }

   fclose(plot_file);


   isnap++;               /* increment snapshot counter */
}     

