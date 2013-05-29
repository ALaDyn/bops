/*************************************************** 

    Electrostatic PIC code for Heraeus School
    laser-plasma simulation practical 

    1d1v,  slab geometry, C version
  
    Paul Gibbon, September 2000


*************************************************/

#include "es.h"			/* common variables */

        /*  external routines  */

extern void init(void);               /* initialisation */ 
extern void loadx(void);              /* particle loading - positions */
extern void loadv(void);              /* particle loading - velocities */
extern void density(void);            /* density gather */
extern void field(void);              /* Poisson solver */
extern void push(void);               /* particle pusher */
extern void boundaries(void);        /* particle boundary conditions */
extern void diagnostics(void);        /* diagnostic routine */

FILE *history_file;     /* file for writing out time histories */            

int main() 
{

  int i;

  history_file = fopen("hist.dat", "w");

  init();			/* do initialisation */
  loadx();			/* load particles onto grid */

  loadv();			/* define velocity distribution */

  density();			/* compute initial density from particles */

  field();			/* compute initial electric field */

  diagnostics();		/* output initial conditions */

 
  for (itime=1; itime<=ntrun; itime++) 
    {

      push();			/* push particles */
      boundaries();		/* particle boundary conditions */

      density();		/* compute density */

      field();		        /* compute electric field (Poisson) */

      diagnostics();		/* output snapshots and time-histories */

    }


  close(history_file);          /*   close time-hist files */


}
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

    nx = 100;        /*  # grid points */
    ne = 1000;       /*  # electrons   */
    ni = 0;          /* # ions (fixed) */

     
    grid_length = 10.0;   /* size of spatial grid */
    plasma_start = 0.;    /* LH plasma edge */

    dt = 0.2;             /* normalised timestep */
     
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
    igraph = 50;      /* freq. of graphical snapshots */
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

/**  =======================================
 **
 **   Particle loading for uniform density
 **
     ======================================= **/

#include "es.h"

/*  local variables  */

float dpx;        /* particle spacing */ 
int i;            /* counter */

void loadx(void)
{
      if (bc_particle >= 2)
      {

         /**  for reflective bcs to work properly, must place particle
              boundaries half a mesh spacing inside field boundaries **/

         wall_left = dx/2.;
         wall_right = grid_length-dx/2.;
	 plasma_end = wall_right;        /* actually want min(end,wr) */
      }
      else 
      {

              /*  periodic boundaries */

         wall_left = 0.;
         wall_right = grid_length;
         plasma_start = 0.;
	 plasma_end = grid_length;
      }

      xload = plasma_end - plasma_start;  /*  length for particle loading */
      dpx = xload/ne;                   /*  particle spacing */

 /*  pseudo-particle charge normalised to give ncrit=1 (rhoc=-1) */

      qe = -rho0*dpx;

 /*  pseudo-particle mass (need for kinetic energy diagnostic) */

      e_mass = qe / q_over_me;

 /*  set up initial positions */

      for (i = 1; i <= ne; i++)
      { 
          x[i] = plasma_start + dpx*(i-0.5);
      }
}
/**  ======================================
 **
 **  Set up particle velocity distribution
 **
 ** ====================================== 
 **/


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

    rs = (float) random() / random_max;        /* random angle */
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
	 j = n1* (float) random()/random_max + 1;

	 temp1 = vx[i];        /* switch i,j */
	 vx[i] = vx[j];
	 vx[j] = temp1;
	}				  
}
/** ======================================
 **
 **   Electron density
 **
 **  ======================================
**/

#include "es.h"

void density(void)
{
   float re, xa, f1, f2;
   int i, j, j1, j2, iwl, iwr;

   re = qe/dx;     /*  charge weighting factor */
     
   for (j=0; j<=nx+1; j++) 
     { rhoe[j]=0.0; }  /* zero density array */

   /*  gather density onto grid from particles */

   for (i=1; i<=ne; i++)
     {
      xa = x[i]/dx;
      j1 = xa;
      j2 = j1 + 1;
      f2 = xa - j1;
      f1 = 1.0 - f2;
      rhoe[j1] = rhoe[j1] + re*f1;
      rhoe[j2] = rhoe[j2] + re*f2;
     }


   if (bc_field == 1)
     {
       /*   periodic boundaries */ 
      rhoe[0] += rhoe[nx];
      rhoe[nx] = rhoe[0];
     }

   else     /* bc_field =2 */
     {
     /*   reflective - 1st and last (ghost) cells folded back onto physical grid */
      iwl = wall_left/dx;
      rhoe[iwl+1] += rhoe[iwl];
      rhoe[iwl] = 0.0;	

      iwr = wall_right/dx;
      rhoe[iwr] += rhoe[iwr+1];
      rhoe[iwr+1] = rhoe[iwr];
     }     

}

/**  ============================================
 **
 **   Direct field integration - slab geometry
 **
 **  ============================================ 
 **/

#include "es.h"


void field(void)
{
   int j;
   float edc;
   float rhot[501];     /*  net density */


  /*  Add neutral background to get net density */
 
   for (j=1; j<=nx+1; j++)
     { rhot[j] = rhoe[j] + rho0; }
      
  /*   integrate div.E=rho directly (trapezium approx) */

  /*   end point - ex=0 mirror at right wall */


  /*  integrate from right to left */

   Ex[nx+1] = 0.0;
   edc = 0.0;

   for (j=nx; j>=1; j--)
      {
	Ex[j] = Ex[j+1] - 0.5*( rhot[j] + rhot[j+1] )*dx;
        edc = edc + Ex[j];
      }


   if (bc_field == 1)
      {
	/*  periodic fields:  subtract off DC component */ 
	/*  -- need this for consistency with charge conservation */
      
	for (j=1; j<=nx; j++)
	  { Ex[j] -= edc/nx; }
          
	/*  end points periodic */
	Ex[0] = Ex[nx];
      }

   /*   potential - same again */
   /*   - integration const. is arbitrary - phi not used for anything at present */

      phi[nx+1] = 0.0;

      for (j=nx; j>=1; j--)
	{
         phi[j] = phi[j+1] + 0.5*( Ex[j] + Ex[j+1] )*dx;
	}
}
 

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
/**  ===================================
 **
 **         push.c
 **
 **   Electrostatic particle pusher
 **
 ** ===================================  
 **/

#include "es.h"

void push(void)
{

  /* locals */

   float xa;            /* fractional particle position */
   float  b1, b2;       /* linear weights */
   float  exi;          /* interpolated e-field */
   int i, j1, j2;          /* grid indices */

   for( i=1; i<=ne; i++)
     {

       /*   interpolate field Ex from grid to particle */

      xa = x[i]/dx;
      j1 = xa;
      j2 = j1 + 1;
      b2 = xa - j1;
      b1 = 1.0 - b2;
      exi = b1*Ex[j1] + b2*Ex[j2];
      
      /*  update velocities */

      vx[i] = vx[i] + q_over_me*dt*exi; 
     }



   /*  update positions (2nd half of leap-frog) */

   for( i=1; i<=ne; i++)
     {
      x[i] = x[i] + dt*vx[i];
     }

}
/**  ===================================
 **
 **   Particle boundary conditions
 **
 ** ===================================  
 **/

#include "es.h"

void boundaries(void)
{
   
  int iseed1 = 28631;        /* random number seed */
  int iseed2 = 1631;         /* random number seed */
  int i, i_inject;

  float xl = grid_length;     /*  shorthand */  
  float wl = wall_left;
  float wr = wall_right;
  float vm, theta, rs;

  /**  loop over all particles to see if any have
       left simulation region: if so, we put them back again
       according to the switch 'bc_particle' **/

  for ( i=1; i<=ne; i++)
    {


      switch (bc_particle) {

      case 1:              /*  periodic */
	   
	    if ( x[i] < 0.0 ) { x[i] += xl; }
	    if ( x[i] >= xl ) { x[i] -= xl; }
	    break;

      case 2:          /* reflective at x = wall_left and x = wall_right */

	    if ( x[i] <= wl )
	      {
	       x[i] = 2*wl - x[i];
	       vx[i] = - vx[i];
	      }
 

	    if ( x[i] >= wr )
	      {
	       x[i] = 2*wr - x[i];
	       vx[i] = - vx[i];
	      }
	    break;

      case 3:

	/**   reflect at LH boundary
	      and re-inject with thermal velocity 
	      at RH (high density) boundary **/


	    if ( x[i] <= wl )
	      {
		x[i] = 2*wl - x[i]; 
		vx[i] = -vx[i];
	      }
       

	    if ( x[i] >= wr )
	      {
		x[i] = 2*wr - x[i];

		/*  find new random velocity in cold Maxwellian */

		rs = (float) random() / random_max;   /* random (0-1) */

		i_inject = ne*rs + 1;
		vm = vte*sqrt(-2.*log((i_inject-0.5)/ne));

		rs = (float) random() / random_max;   /* random (0-1) */
		theta = pi*rs;   /* 1st quadrant */
		vx[i] = -vm*sin(theta);
	      }	
	    break;

      default:           /*  default to periodic */
	   
	    if ( x[i] < 0.0 ) { x[i] += xl; }
	    if ( x[i] >= xl ) { x[i] -= xl; }
	    break;
      }

    }

}
/**  ============================================
 **
 **   Time-histories
 **
 ** ============================================
 **/

#include <stdio.h>
#include <math.h>
#include "es.h"

FILE *history_file;     /* file for writing out time histories */            

void histories(void)
{

  int i;
  float ukin;

  /*   kinetic energy */
 
  ukin = 0.0;

  for ( i=1; i<=ne; i++ )
    {
      ukin += 0.5*e_mass*vx[i]*vx[i];
    }

  /*  write energies out to file */

   fprintf( history_file, "%f  %f\n", itime*dt, ukin );

}

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

   sprintf(cfile,"density%d.dat", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     {
       fprintf( plot_file, "%f  %f\n" ,xgrid[j],-rhoe[j] );
     }

   fclose(plot_file);



/*   electrostatic field */

   sprintf(cfile,"field%d.dat", isnap ); 
   plot_file = fopen(cfile, "w");

   for (j=0; j<=nx; j++)
     {
       fprintf( plot_file, "%f  %f\n" ,xgrid[j],Ex[j] );
     }

   fclose(plot_file);


   isnap++;               /* increment snapshot counter */
}     

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
#define NP_MAX 5001             /* max. # particles */
#define NX_MAX 501              /* max. # grid points */

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

