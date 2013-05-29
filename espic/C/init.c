/************************************************
 **
 **              init.c
 **
 **  initialise main particle and grid variables
 **
 *************************************************/


#include <stdio.h>
#include <math.h>

#include "es.h"


void init(void)
{

  /*  local variables  */

    int ncell;
    float xdodx;

   /* set up simulation parameters */

    ntrun = 100;     /* # timesteps  */

    nx = 200;        /*  # grid points */
    ne = 5000;       /*  # electrons   */
    ni = 0;          /* # ions (fixed) */

     
    grid_length = 2*pi;   /* size of spatial grid */
    plasma_start = 0.;    /* LH plasma edge */

    dt = 0.1;             /* normalised timestep */
     
    q_over_me=-1.0;       /* electron charge:mass ratio */

    rho0 = 1.0;           /* background ion density */
    vte = 0.1;            /* thermal velocity */

    ilas = 0;             /*  laser switch:   0 = off  */
                          /*	              1 = uniform sine wave    */
    w0=1.0;               /*  laser frequency   */
    a0 = 0.1;             /*  laser amplitude */

    bc_field = 1;      /*  field boundary conditions:  1 = periodic, */
                       /*                              2 = reflective */

    bc_particle = 1;   /*  particle BCs:  1 = periodic, */
                       /*                 2 = reflective, */
                       /*                 3 = thermal */

    ihist = 1;         /* frequency of time-history output */
    igraph = 100;      /* freq. of graphical snapshots */
    
    iout = 10;         /* freq. of printed diags. */

    itime = 0;         /* initialise time counter */



    /**   ------------------
     **   derived parameters 
     **   ------------------ **/

    dx = grid_length/nx;     /*  mesh size */

    omega_p = sqrt(rho0);    /* plasma frequency */
    x_debye = vte/omega_p;   /* Debye length */
    xdodx = x_debye/dx;      /* ratio */
    ncell = ne/nx;           /* # particles per cell */

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
