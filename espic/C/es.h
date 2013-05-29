/**  
 **       es.h
 **
 ** Header file for pices.c 
 ** containing global arrays and variables
 **
 **/

#include<stdio.h>
#include<math.h>
#include<time.h>
   
#define pi 3.14159265
#define random_max 2147483647   /* max. allowed integer (32-bit) */
#define NP_MAX 50001             /* max. # particles */
#define NX_MAX 5001              /* max. # grid points */

/* particle arrays */ 

 float   x[NP_MAX];         /* position */
 float  vx[NP_MAX];         /* velocity */ 

/* grid arrays */

 float  rhoe[NX_MAX];        /* electron density */
 float  Ex[NX_MAX];           /* electric field */ 
 float  phi[NX_MAX];          /* potential */
 float  a[NX_MAX];           /* laser amplitude */


/* plasma constants, parameters */

 float qe;                /* macro-electron charge */
 float vte;               /* thermal velocity */
 float q_over_me;         /* charge/mass ratio */
 float e_mass;            /* macro-electron mass */
 float xload;             /* plasma length */
 float plasma_start;      /* LH edge of plasma */
 float plasma_end;        /* RH edge of plasma */
 float wall_left;         /* LH particle boundary */
 float wall_right;        /* RH particle boundary */
 float rho0;              /* background (ion) density */
 float grid_length;       /* grid length */
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
 float  w0;              /* laser frequency */
 float  trise;           /* laser rise-time */
 int  ilas;              /* laser switch */

/* diagnostic parameters */
 
 int  iout;             /* frequeny of printed output */
 int  igraph;           /* frequeny of graphical output */
 int  ihist;           /* frequeny of time-history output */

