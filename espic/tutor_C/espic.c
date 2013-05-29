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
extern void add_neutral(void);           /* initial ion density */
extern void field(void);              /* Poisson solver */
extern void push(void);               /* particle pusher */
extern void particle_bc(void);        /* particle boundary conditions */
extern void diagnostics(void);        /* diagnostic routine */

FILE *history_file;     /* file for writing out time histories */            

int main() 
{

  int i;

  history_file = fopen("hist.data", "w");

  init();			/* do initialisation */
  loadx();			/* load particles onto grid */
  particle_bc();		/* particle boundary conditions */
                     
  loadv();			/* define velocity distribution */

  density();			/* compute initial density from particles */
  add_neutral();                   /* background ion density */
  laser();                      /* initialise laser field */
  field();			/* compute initial electric field */

  diagnostics();		/* output initial conditions */

 
  for (itime=1; itime<=ntrun; itime++) 
    {

      push();			/* Push particles */
      particle_bc();		/* particle boundary conditions */

      density();		/* compute density */

      laser();                  /* update laser field */
      field();		        /* compute electric field (Poisson) */

      diagnostics();		/* output snapshots and time-histories */

    }


  close(history_file);          /*   close time-hist files */


}
