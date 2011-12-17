
!  =================================================================
!
!  Particle boundary conditions (ibound=ipbc in input deck)
!
!     ibound = 1     periodic
!              2     reflective
!              3     absorb/reemit at both sides
!              4     absorb ions at LHB, electrons only if charged
!             5     absorb electrons, reflect ions at RHB
!  $ Revision: $
!
!  TODO:  write out energy spectrum when escaped particle arrays full
!	  need to adapt book-keeping to include protons
! 
!  =================================================================

subroutine pbcs(i1,n,vt,mass,ibound)
  use bopsvars
  implicit none

  real*8 uxi,uyi,uzi,uxo,uyo,uzo,gami,gamout,gamin,g
  real*8 uxtot, uytot, gtot, vt, mass
  integer :: idum,l, n, i1, ibound, iswap, ip


  data idum,uxtot,uytot,gtot/-3,0.,0.,0./
  save idum,uxtot,uytot,gtot

  do l=1,n
     ip=i1+l-1

     !  ************************************
     !   periodic
     !  ************************************

     if (ibound.eq.1) then
        if (xn(ip).lt.0.) xn(ip)=xn(ip)+xl
        if (xn(ip).ge.xl) xn(ip)=xn(ip)-xl

        !  ************************************
        !   reflective at x=wl and x=wr (left and right walls)
        !    - gamma, uy invarient
        !  ************************************

     else if (ibound.eq.2) then
        if (xn(ip).le.wl) then
           xn(ip)=2*wl-xn(ip)
           ux(ip)=-ux(ip)
        endif
        if (xn(ip).ge.wr) then
           xn(ip)=2*wr-xn(ip)
           ux(ip)=-ux(ip)
        endif

        !  ************************************
        !  absorption and re-emission
        !  ************************************

     else if (ibound.eq.3) then
        !  first do LH boundary:
        if (xn(ip).le.wl) then
           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           !  find new random position in cold maxwellian
           if (i1.eq.1) then
              call reinj(vt,n,1,uxi,uyi,uzi,g)
           else
              call reinji(vt,n,1,uxi,uyi,uzi,g)
           endif

           gamin=g

           !  add net escape energy at boundary - lab frame
           elhb=elhb + (gamout-gamin)*mass/gam0**2

           !  boost back to simulation frame
           call lab2sim(vy0,gam0,uyi,g,uy(ip),gamma(ip))
           ux(ip)=uxi
           uz(ip)=uzi
           !  re-inject at left wall with vthermal
           xn(ip)=wl+uxi/gamma(ip)*dt

           !  now do RH boundary
        else if (xn(ip).ge.wr) then

           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           !  find new random position in cold maxwellian (lab frame)
           !  chosen to conserve flux in ux
           if (i1.eq.1) then
              call reinj(vt,n,-1,uxi,uyi,uzi,g)
           else
              call reinji(vt,n,-1,uxi,uyi,uzi,g)
           endif

           gamin=g

           !  add escapee to list:
           !  store net energy loss from system (lab frame)

           if (i1.eq.1) then
              !  electron
              nesc=nesc+1
              if (nesc.ge.npart) nesc=1
              iesc(nesc)=ip
              uesc(nesc)=gamout-1.d0
              uinj(nesc)=gamin-1.d0
              !  write out exit data: label, time in fs, energy in keV
!              write(90,'(i8,f12.3,f15.4)') ip,dt*itime/tcfs,uesc(nesc)*511
           else
              !  ion
              nesci=nesci+1
              if (nesci.ge.npart) nesci=1
              iesci(nesci)=ip
              uesci(nesci)=gamout-gamin
           endif

           !  add net escape energy  - lab frame
           erhb=erhb+(gamout-gamin)*mass/gam0**2

           !  boost back to simulation frame
           call lab2sim(vy0,gam0,uyi,g,uy(ip),gamma(ip))
           ux(ip)=uxi
           uz(ip)=uzi
           !  reinject at right wall
           xn(ip)=wr+uxi/gamma(ip)*dt
        endif

        !  ************************************
        !  Absorb selectively at LHB, reemit at RHB
        !  ************************************

     else if (ibound.eq.4) then


        !  first do LH boundary for ions:
        if (xn(ip).le.wl.and.ip.gt.ne) then

           !  boost back to lab frame
           call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyo,gamout)

           !  write out ion momenta to file
!           write(82,'(3(1pe14.4))') ux(ip),uyo,uz(ip)


           !  add net escape energy at boundary - lab frame
           elhb=elhb + gamout*mass/gam0**2

           !  absorb ion: swap with last one in list
           iswap=ion1+ni-1
           xn(ip)=xn(iswap)
           xo(ip)=xo(iswap)
           ux(ip)=ux(iswap)
           uy(ip)=uy(iswap)
           uz(ip)=uz(iswap)
           gamma(ip)=gamma(iswap)
           ni=ni-1
           npart=npart-1

        else if (xn(ip).le.wl.and.ip.le.ne) then

           !  Absorb electron if plasma charged ...
           !   (ie: ion has escaped previously)

           if (ni.gt.0.and.ne.gt.Z*(ni-np)+np) then

              !  boost back to lab frame
              call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyo,gamout)

              !  write out electron momenta to file
!              write(80,'(3(1pe14.4))') ux(ip),uyo,uz(ip)

              !  add net escape energy at boundary - lab frame
              elhb=elhb + gamout*mass/gam0**2

              iswap=ne
              xn(ip)=xn(iswap)
              xo(ip)=xo(iswap)
              ux(ip)=ux(iswap)
              uy(ip)=uy(iswap)
              uz(ip)=uz(iswap)
              gamma(ip)=gamma(iswap)
              ne=ne-1
              npart=npart-1

           else
              ! .. otherwise reflect

              !  boost back to lab frame
              call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyo,gamout)

              !  write out electron momenta to file
!              write(80,'(3(1pe14.4))') ux(ip),uyo,uz(ip)
              xn(ip)=2*wl-xn(ip)
              ux(ip)=-ux(ip)
           endif

           !  RH boundary:

           !  absorb and re-emit electrons
        else if (xn(ip).ge.wr.and.ip.le.ne) then

           !  boost back to lab frame
           call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyo,gamout)

           !  write out energetic electron momenta to file
!           write(81,'(3(1pe14.4))') ux(ip),uyo,uz(ip)


           !  find new random position in cold maxwellian (lab frame)
           !  chosen to conserve flux in ux

           call reinj(vt,n,-1,uxi,uyi,uzi,g)
           gamin=g

           !  add escapee to list:
           !  store net energy loss from system (lab frame)

           !  electron
           nesc=nesc+1
           if (nesc.ge.10*npart) nesc=1
           iesc(nesc)=ip
           !  store outgoing and reinjected kinetic energies
           uesc(nesc)=gamout-1.d0
           uinj(nesc)=gamin-1.d0
           !  write out exit data: label, time in w0^-1, energy in keV
!           write(90,'(i8,f12.3,f15.4)') ip,dt*itime/tcfs,uesc(nesc)*511

           !  add net escape energy  -  lab frame
           erhb=erhb+(gamout-gamin)*mass/gam0**2

!  Debug reinjection
!	write(*,'(i6,a,4(1pe12.3))') itime,'reinject rhb: ',erhb,uxi,uyi,uzi 

           !  boost back to simulation frame
           call lab2sim(vy0,gam0,uyi,g,uy(ip),gamma(ip))
           ux(ip)=uxi
           uz(ip)=uzi
           !  reinject at right wall
           !          xn(ip)=wr+uxi/gamma(ip)*dt
           xn(ip)=2*wr-xn(ip)

           !  reflect ions
        else if (xn(ip).ge.wr.and.ip.gt.ne) then
           xn(ip)=2*wr-xn(ip)
           ux(ip)=-ux(ip)
        endif

        !  ************************************

        !  Absorb selectively at LHB,
        !  absorb electrons, reflect ions at RHB

        !  ************************************

     else if (ibound.eq.5) then


        !  first do LH boundary for ions:
        if (xn(ip).le.wl.and.ip.gt.ne) then

           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           !  add net escape energy at boundary - lab frame
           elhb=elhb + gamout*mass/gam0**2

           !  absorb ion: swap with last one in list
           iswap=ion1+ni-1
           xn(ip)=xn(iswap)
           xo(ip)=xo(iswap)
           ux(ip)=ux(iswap)
           uy(ip)=uy(iswap)
           uz(ip)=uz(iswap)
           gamma(ip)=gamma(iswap)
           ni=ni-1
           npart=npart-1

        else if (xn(ip).le.wl.and.ip.le.ne) then

           !  Absorb electron if plasma charged ...
           !   (ie: ion has escaped previously)

           if (ni.gt.0.and.ne.gt.Z*(ni-np)+np) then

              !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
              gamout=gam0*( gamma(ip) - vy0*uy(ip) )

              !  add net escape energy at boundary - lab frame
              elhb=elhb + gamout*mass/gam0**2

              iswap=ne
              xn(ip)=xn(iswap)
              xo(ip)=xo(iswap)
              ux(ip)=ux(iswap)
              uy(ip)=uy(iswap)
              uz(ip)=uz(iswap)
              gamma(ip)=gamma(iswap)
              ne=ne-1
              npart=npart-1

           else
              ! .. otherwise reflect
              xn(ip)=2*wl-xn(ip)
              ux(ip)=-ux(ip)
           endif

           !  RH boundary:


        else if (xn(ip).ge.wr.and.ip.le.ne) then

           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           if (ux(ip)/gamma(ip).lt.6*vt) then

              !  find new random position in cold maxwellian (lab frame)
              !  chosen to conserve flux in ux

              call reinj(vt,n,-1,uxi,uyi,uzi,g)
              gamin=g

              !  add escapee to list:
              !  store net energy loss from system (lab frame)

              !  electron
              nesc=nesc+1
              if (nesc.ge.npart) nesc=1
              iesc(nesc)=ip
              uesc(nesc)=gamout-1.d0
              uinj(nesc)=gamin-1.d0
              !  write out exit data: label, time in w0^-1, energy in keV
!              write(90,'(i8,f12.3,f15.4)') ip,dt*itime/tcfs,uesc(nesc)*511

              !  boost back to simulation frame
              call lab2sim(vy0,gam0,uyi,g,uy(ip),gamma(ip))
              ux(ip)=uxi
              uz(ip)=uzi
              xn(ip)=wr-uxi/gamma(ip)*dt

           else

              !  absorb fast electrons - allow charge build-up!

              nesc=nesc+1
              if (nesc.ge.npart) nesc=1
              iesc(nesc)=ip
              uesc(nesc)=gamout-1.d0
              gamin=1.

              iswap=ne
              xn(ip)=xn(iswap)
              xo(ip)=xo(iswap)
              ux(ip)=ux(iswap)
              uy(ip)=uy(iswap)
              uz(ip)=uz(iswap)
              gamma(ip)=gamma(iswap)
              ne=ne-1
              npart=npart-1
           endif

           !  add net escape energy  -  lab frame
           erhb=erhb+(gamout-gamin)*mass/gam0**2

           !  reflect ions
        else if (xn(ip).ge.wr.and.ip.gt.ne) then
           xn(ip)=2*wr-xn(ip)
           ux(ip)=-ux(ip)
        endif

        !  ************************************
        !  Foil boundaries: absorb selectively at both LHB, RHB
        !  ************************************

     else if (ibound.eq.6) then


        !  first do LH boundary for ions:
        if (xn(ip).le.wl.and.ip.gt.ne) then

           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           !  add net escape energy at boundary - lab frame
           elhb=elhb + gamout*mass/gam0**2

           !  absorb ion: swap with last one in list
           iswap=ion1+ni-1
           xn(ip)=xn(iswap)
           xo(ip)=xo(iswap)
           ux(ip)=ux(iswap)
           uy(ip)=uy(iswap)
           uz(ip)=uz(iswap)
           gamma(ip)=gamma(iswap)
           ni=ni-1
           npart=npart-1

        else if (xn(ip).le.wl.and.ip.le.ne) then

           !  Absorb electron if plasma charged ...
           !   (ie: ion has escaped previously)

           if (ni.gt.0.and.ne.gt.Z*(ni-np)+np) then

              !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
              gamout=gam0*( gamma(ip) - vy0*uy(ip) )

              !  add net escape energy at boundary - lab frame
              elhb=elhb + gamout*mass/gam0**2

              iswap=ne
              xn(ip)=xn(iswap)
              xo(ip)=xo(iswap)
              ux(ip)=ux(iswap)
              uy(ip)=uy(iswap)
              uz(ip)=uz(iswap)
              gamma(ip)=gamma(iswap)
              ne=ne-1
              npart=npart-1

           else
              ! .. otherwise reflect
              xn(ip)=2*wl-xn(ip)
              ux(ip)=-ux(ip)
           endif

           !  absorb ions at RH boundary:

        else if (xn(ip).ge.wr.and.ip.gt.ne) then

           !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
           gamout=gam0*( gamma(ip) - vy0*uy(ip) )

           !  add net escape energy at boundary - lab frame
           erhb=erhb + gamout*mass/gam0**2

           !  absorb ion: swap with last one in list
           iswap=ion1+ni-1
           xn(ip)=xn(iswap)
           xo(ip)=xo(iswap)
           ux(ip)=ux(iswap)
           uy(ip)=uy(iswap)
           uz(ip)=uz(iswap)
           gamma(ip)=gamma(iswap)
           ni=ni-1
           npart=npart-1

        else if (xn(ip).ge.wr.and.ip.le.ne) then
           !  Absorb electron if plasma charged ...
           !   (ie: ion has escaped previously)

           if (ni.gt.0.and.ne.gt.Z*(ni-np)+np) then

              !  escaped particle energy (lab frame) - normalised to mc**2 (=gamma-1)
              gamout=gam0*( gamma(ip) - vy0*uy(ip) )

              !  add net escape energy at boundary - lab frame
              erhb=erhb + gamout*mass/gam0**2

              iswap=ne
              xn(ip)=xn(iswap)
              xo(ip)=xo(iswap)
              ux(ip)=ux(iswap)
              uy(ip)=uy(iswap)
              uz(ip)=uz(iswap)
              gamma(ip)=gamma(iswap)
              ne=ne-1
              npart=npart-1

           else
              ! .. otherwise reflect
              xn(ip)=2*wr-xn(ip)
              ux(ip)=-ux(ip)
           endif

        endif

     endif

  end do
end subroutine pbcs
