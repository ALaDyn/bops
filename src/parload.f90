!
!    PARLOAD
!
!     Configure particle distn

!     $Revision: 26 $
!
!     ========================================================

subroutine parload
  use bopsvars
  implicit none

  integer ip, l, i

  integer ipow2(0:20)
  real*8 uth0, u2

  data ipow2/1,2,4,8,16,32,64,128,256,512,1024,2048,4096, &
       8192,16384,32768,65536,131072,262144,524288,1048576/


  !  electron density profile
  call denprof(1,ne,qe)

!  ion density profile
!#### Anupam & Bin 2009/2010:! change ni with ni_tot including protons
! CHECK: notice qi is improper here, even for three layers target it is ok ! 

  if (ni_tot.gt.0) call denprof(ion1,ni_tot,qi)

! assign particle charges, masses, and labels

! electrons
  do i=1,ne
    q(i) = qe
    m(i) = me
    species(i) = 1
  end do

!#### Anupam & Bin 2009/2010: set up for multispecies_target
! add protons for proton layers
! protons
  if (np .gt. 0) then
    do i=ne+1,ne+np
      q(i) = qp
      m(i) = mp
      species(i) =3
    end do  
  end if
!#### 

! ions
 if ((ni_tot-np) .gt. 0) then
  do i=ne+np+1,ne+ni_tot        ! here ni_tot-np=ni
    q(i) = qi
    m(i) = mi
    species(i) = 2
  end do
 end if
!#### Anupam & Bin 2009/2010

!  relativistic thermal distribution in vx, vy, vz
!  - load 1st octant in each direction

!#### Anupam & Bin 2009/2010: add relativistic thermal distribution for proton layers 

! Proton layer on front
  if (n_pp1.gt.0) then
     !  electrons
     call urelmax(ux(1:n_pp1),n_pp1,vte)
     call urelmax(uy(1:n_pp1),n_pp1,vte)
     call urelmax(uz(1:n_pp1),n_pp1,vte)
     call fill3v(1,n_pp1)
     ! protons  
     call urelmax(ux(ne+1:ne+n_pp1),n_pp1,vtp)
     call urelmax(uy(ne+1:ne+n_pp1),n_pp1,vtp)
     call urelmax(uz(ne+1:ne+n_pp1),n_pp1,vtp)
     !  fill in remaining 7 octants
     call fill3v(ne+1,n_pp1)
  end if

! Proton layer on back
  if (n_pp2.gt.0) then
    ! electrons
     call urelmax(ux(1+n_pp1:np),n_pp2,vte)
     call urelmax(uy(1+n_pp1:np),n_pp2,vte)
     call urelmax(uz(1+n_pp1:np),n_pp2,vte)
     call fill3v(1+n_pp1,n_pp2)
    ! protons
     call urelmax(ux(ne+1+n_pp1:ne+np),n_pp2,vtp)
     call urelmax(uy(ne+1+n_pp1:ne+np),n_pp2,vtp)
     call urelmax(uz(ne+1+n_pp1:ne+np),n_pp2,vtp)
     !  fill in remaining 7 octants
     call fill3v(ne+1+n_pp1,n_pp2)
  end if
!#### Anupam & Bin 2009/2010

! Electrons in heavy ion substrate layer
  if ((ne-np) .gt. 0) then
   call urelmax(ux(np+1:ne),ne-np,vte)
   call urelmax(uy(np+1:ne),ne-np,vte)
   call urelmax(uz(np+1:ne),ne-np,vte)
  !  fill in remaining 7 octants
   call fill3v(np+1,ne-np)
  end if 
!  heavy ions
 if ((ni_tot-np) .gt. 0) then 
  call urelmax(ux(ne+np+1:ne+ni_tot),ni_tot-np,vti)
  call urelmax(uy(ne+np+1:ne+ni_tot),ni_tot-np,vti)
  call urelmax(uz(ne+np+1:ne+ni_tot),ni_tot-np,vti)
  !  fill in remaining 7 octants
  call fill3v(ne+np+1,ni_tot-np)
 end if
!  scramble positions and velocities

  if (n_pp1.gt.0) call scramb(1,n_pp1)
  if (n_pp2.gt.0) call scramb(n_pp1+1,n_pp2)
  if ((ne-np).gt.0) call scramb(np+1,ne-np)  
  if (n_pp1.gt.0) call scramb(ne+1,n_pp1)
  if (n_pp2.gt.0) call scramb(ne+n_pp1+1,n_pp2) 
  if ((ni_tot-np).gt.0) call scramb(ne+np+1,ni_tot-np)

!  ni is now # heavy ions, ni_tot is total # ions

!  calculate relativistic gamma factor
  write (*,*) npart, ne, ni, np
  do ip=1,npart
     gamma(ip) = sqrt(1.0+ ux(ip)**2 + uy(ip)**2 + uz(ip)**2)
  end do

  !  initial kinetic energy in lab frame
  uth0=0.
  do l=1,ne
     u2=ux(l)**2+uy(l)**2+uz(l)**2
     uth0=uth0+me*u2/(gamma(l)+1.0)
  end do

!#### Anupam & Bin 2009/2010: Changed for calculating protons separately
  if (np .gt. 0) then
    do l=ne+1,ne+np
       u2=ux(l)**2+uy(l)**2+uz(l)**2
       uth0=uth0+mp*u2/(gamma(l)+1.0)           ! changed for protons
    end do
  end if 
 if ((ni_tot-np).gt.0) then
  do l=ne+np+1,ne+ni_tot
     u2=ux(l)**2+uy(l)**2+uz(l)**2
     uth0=uth0+mi*u2/(gamma(l)+1.0)
  end do
 end if

  !  relativistic boost along y-axis for oblique incidence
  if (the0.gt.0) call boost


  ! get initial densities for first plot
  call eden

  if (ni.eq.0 .or. target_config==0) then
     call addneut
  else
     call iden(ion1)  !  ion1 is index of 1st ion (heavy)
  endif

!  Initial current densities
  call jeden
  call jiden(ion1)

end subroutine parload





















