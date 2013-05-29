/**  
 **       es.h
 **
 ** Header file for pices.c 
 ** containing global arrays and variables
 **
 **/

#include<stdio.h>
#include<math.h>
   
#define pi 3.14159265
#define random_max 2147483647   /* max. allowed integer (32-bit) */
#define NP_MAX 20001             /* max. # particles */
#define NX_MAX 501              /* max. # grid points */

/* particle arrays */ 

 float   x[NP_MAX];         /* position */
 float  vx[NP_MAX];         /* velocity */ 

/* grid arrays */

 float  rhoe[NX_MAX];        /* electron density */
 float  rhoi[NX_MAX];        /* ion density */
 float  Ex[NX_MAX];           /* electric field */ 
 float  phi[NX_MAX];          /* potential */
 float Elaser[NX_MAX];           /* laser field or pond force  */


/* plasma constants, parameters */

 float qe;                /* macro-electron charge */
 float q_over_me;         /* charge/mass ratio */
 float e_mass;            /* macro-electron mass */

 float xload;             /* plasma length */
 float plasma_start;      /* plasma start */
 float plasma_end;        /* plasma end */
 float wall_left;         /* LH particle boundary */
 float wall_right;        /* RH particle boundary */
 float grid_length;       /* grid length */

 float rho0;              /* background (ion) density */
 float vte;               /* thermal velocity */
 float omega_p;           /* plasma frequency */
 float x_debye;           /* Debye-length */ 

 int ne, ni;              /* # simulation electrons, ions */

/* grid variables */

 int  ntrun;             /* # timesteps */
 int  itime;             /* current timestep */
 float  dt;              /* timestep */

 int  nx;                /* # grid points */
 float  dx;              /* mesh size */

 int  bc_particle;       /* particle boundary condition switch */
 int  bc_field;          /* field boundary condition */

/* laser parameters */

 float  a0;              /* laser amplitude (= vosc/c) */
 float tpulse;           /* pulse length */
 float tdel;           /* pulse delay */
 float xcentre;          /* pulse centroid */
 float  w0;              /* laser frequency */
 int  ilas;              /* laser switch */

/* diagnostic parameters */
 
 int  iout;             /* frequeny of printed output */
 int  igraph;           /* frequeny of graphical output */
 int  ihist;           /* frequeny of time-history output */

 int  profile;          /* density profile switch */
