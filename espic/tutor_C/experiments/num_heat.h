   /* set up simulation parameters */

    ntrun = 500;     /* # timesteps  */

    nx = 10;        /*  # grid points */
    ne = 5000;       /*  # electrons   */
    ni = 0;          /* # ions (fixed) */

     
    grid_length = 10.0;   /* size of spatial grid */

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

    ihist = 5;         /* frequency of time-history output */
    igraph = 500;      /* freq. of graphical snapshots */
    iout = 50;         /* freq. of printed diags. */

    itime = 0;         /* initialise time counter */

