
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
   float fv[NX_MAX];           /* velocity distribution */

   char cfile[40];               /* plot filename string */

   int nvx = 100;        /* # velocity points */
   float vmax = 5*vte;   /* max velocity */
   float vax, dv;

   int i,j, iv;

   for (j=0; j<=nx; j++)
      {
	xgrid[j] = j*dx;   /* set up x-values for grid plots */
      }


 /*  electron px-x phase space */

 /*  compute filename including plot id character '0-9' */

   sprintf(cfile,"phase%d.data", isnap );   
   plot_file = fopen(cfile, "w");

   for (i=1; i<=ne; i++)
     {
       fprintf( plot_file, "%f  %f\n" ,x[i],vx[i] );
     }

   fclose(plot_file);



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

/*   laser field */

   sprintf(cfile,"laser%d.data", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     {
       fprintf( plot_file, "%f  %f\n" ,xgrid[j],Elaser[j] );
     }

   fclose(plot_file);


/*   potential */

   sprintf(cfile,"potential%d.data", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     { fprintf( plot_file, "%f  %f\n" ,xgrid[j], phi[j] );}

   fclose(plot_file);


/*  distribution function */

   sprintf(cfile,"vel_dist%d.data", isnap ); 
   plot_file = fopen(cfile, "w");

/* first compute f(vx) on grid */

   for (j=1; j<=nvx; j++) {fv[j]=0;}      /* zero distn. fn. array */

   dv = 2*vmax/nvx;    /* bin separation */

   for (i=1; i<=ne; i++)
     {
       vax= ( vx[i] + vmax )/dv;       /* normalised velocity */
       iv = vax+1;         /* index */	 
       if (iv <= nvx & iv > 0) { fv[iv]++; }  /* increment dist. fn if within range */
      }


   for (j=1; j<=nvx; j++)
     {       
       fprintf( plot_file, "%f  %f\n" ,j*dv-vmax, fv[j] );
     }

   fclose(plot_file);


   isnap++;               /* increment snapshot counter */
}     

