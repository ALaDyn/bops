   /* set up simulation parameters */

    ntrun = 200;     /* # timesteps  */

    nx = 100;        /*  # grid points */
    ne = 5000;       /*  # electrons   */
    ni = 0;          /* # ions (fixed) */

     
    grid_length = 20.;   /* size of spatial grid */
    plasma_start = 10.;   /* plasma edge */
    plasma_end = grid_length;


    dt = 0.1;             /* normalised timestep */
     
    q_over_me=-1.0;       /* electron charge:mass ratio */

    rho0 = 1.0;           /* background ion density */
    vte = 0.04;            /* thermal velocity */

    ilas = 1;             /*  laser switch:   0 = off  */
                          /*	              1 = uniform sine wave    */
                          /*	              2 = pond. force v_g=c    */
    w0=0.3;               /*  laser frequency   */
    a0 = 0.2;             /*  laser amplitude */
    tpulse = 1.;
    tdel = 2.;

    bc_field = 2;      /*  field boundary conditions:  1 = periodic, */
                       /*                              2 = reflective */

    bc_particle = 3;   /*  particle BCs:  1 = periodic, */
                       /*                 2 = reflective, */
                       /*                 3 = thermal */

    profile = 1;   /* density profile switch */

    ihist = 5;         /* frequency of time-history output */
    igraph = 50;      /* freq. of graphical snapshots */
    iout = 50;         /* freq. of printed diags. */

    itime = 0;         /* initialise time counter */



