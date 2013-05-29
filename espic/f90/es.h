!  ===================================
!
!    es.h  Header file for ES PIC code
!
!  ===================================  

parameter(pi=3.141592654,c=3.e8, &
	  nxm=500,               &
	  npm=5000,              &           !  constants
	  nvxm=200,              &
	  ntm=300                & 
	  )

!  Variable declarations which deviate from defaults (integer i-n, otherwise real)

integer :: bc_particle, bc_field


!  ** particle arrays **

common/part/               &         
          x(0:npm), vx(0:npm)

!  ** grid arrays **

common/grid/  & 
          rhoe(0:nxm), Ex(0:nxm), phi(0:nxm), a(0:nxm)


!  ** particle constants/ parameters **

common/phys/           &  
         qe,   vte,   q_over_me,  e_mass, ne,   &
         xload, wall_left, wall_right, rho0, grid_length


!  ** grid variables **

common/grid/            &  
         nt,    nx,  nxo2,   dx,   xn0,   rdx,   dkx,   &
         dt, itime,    yi,   nvx,   nvy,   dvx,   dvy,   vxm,  &
         bc_particle,  bc_field, ntrun


!  ** laser parameters **

common/las/             &  
         a0,    w0, trise,  ilas,  tdel,    tp, xlam

!  ** diagnostic parameters **

common/diag/             &  
       iout,  ihist,   igraph
