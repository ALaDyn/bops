/************************************************
 **
 **              init.c
 **
 **  initialise main particle and grid variables
 **
 *************************************************/



#include "es.h"


void init(void)
{

  /*  local variables  */

    int ncell;
    float xdodx;

/* read input data */

#include "inputs.h"



    /**   ------------------
     **   derived parameters 
     **   ------------------ **/

    dx = grid_length/nx;     /*  mesh size */

    omega_p = sqrt(rho0);    /* plasma frequency */
    x_debye = vte/omega_p;   /* Debye length */
    xdodx = x_debye/dx;      /* ratio */
    ncell = ne/nx;           /* # particles per cell */

    xcentre = -tdel;        /* initialise laser centroid */

/*      printf("time: %f, x: %f, y: %f\n", time, x, y); */  

    printf("\n Input parameters: \n\n");
    printf("# particles = %d\n", ne);
    printf("# mesh points = %d\n", nx);
    printf("# particles/cell = %d\n", ncell);
    printf("grid length = %f\n\n", grid_length);

    printf("thermal velocity = %f\n",vte);
    printf("mesh size = %f\n",dx);
    printf("Debye length = %f\n\n",x_debye);

    printf("timestep = %f\n",dt);
    printf("# timesteps = %d\n", ntrun);
    printf("run time = %f\n",dt*ntrun);


}
