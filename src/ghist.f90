
!     =====================
!
!     Time history graphs
!
!     =====================

      subroutine ghist
      use bopsvars
      implicit none

      integer i, icycdel, iu
      real*8 dtav, tmax
      real*8 tmin, tstep, abstep, abmin, abmax
      real*8 etat(0:nt/itav+1),etah(0:nt/itav+1),etai(0:nt/itav+1)
      real*8 :: work1(nt/itav+1), work2(nt/itav+1), xt(0:max(nt/itav,ntc)+1)
      character chx*15

!  limits (if needed)
      tmin=0.
      tmax=nint(nt*dt)/tconv
      tstep=tmax/5.
      abmin=0.
      abmax=1.
      abstep=abmax/5.

      if (ntc.eq.0) return
      ntc=ntc-1
      if (iunits.eq.2) then
         chx = '     t/fs      '
      else if (iunits.eq.1) then
         chx = '     w_pt      '
      else
         chx = '    w_0t       '
      endif

      if (.not.lcycave) then

         do i=1,ntc
            xt(i)=dt*(i-1)*itc+(itim0-1)*dt
            xt(i)=xt(i)/tconv
         end do


         call grxy(xt,ues,ntc,15100,1,1 &
              ,chx,'      Ues      ','uesp            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,uth,ntc,15200,1,1 &
              ,chx,'      Uth      ','uthm            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,usys,ntc,15600,1,1 &
              ,chx,'      Usys     ','usys            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,erh,ntc,15500,1,1 &
              ,chx,'     Urhb      ','urhb            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,elh,ntc,12100,1,1 &
              ,chx,'     Ulhb      ','ulhb            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,uem,ntc,15000,1,1 &
              ,chx,'      Uem      ','ueme            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,utot,ntc,15900,1,1 &
              ,chx,'      Utot     ','utot            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     max electron density
         call grxy(xt,edenmax,ntc,17100,1,1 &
              ,chx,'       ne      ','nehi            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)
!     max ion density
         call grxy(xt,idenmax,ntc,17000,1,1 &
              ,chx,'       ni      ','nihi            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


      else

!
!     time-averaged quantities
         ntav=ntav-1
         call i0prnt('  itav',itav)
         call i0prnt('  ntav',ntav)
         call r0form('   tav',tav,'f12.2')
         dtav=dt*itav

         do i=1,ntav
            work1(i)=uinc(i)
            work2(i)=xi0(i)
            xt(i)=(i-1)*dtav/tconv
         end do

         call grxy(xt,work1,ntav,15700,1,1 &
              ,chx,'     <P_L>     ','uinc            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)



         call grxy(xt,uback,ntav,15800,1,1 &
              ,chx,'     <P_R>     ','ubac            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,utrans,ntav,15750,1,1 &
              ,chx,'     <P_T>     ','utra            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

         call grxy(xt,eta,ntav,16000,1,-1 &
              ,chx,'  Abs (1-R-T)  ','absr            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     Electrostatic energy
         call grxy(xt,uesa,ntav,15100,1,1 &
              ,chx,'     <Ues>     ','uesp            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Electromagnetic energy - TE fields
         call grxy(xt,uema,ntav,15000,1,1 &
              ,chx,'     <U(TE)>   ','ueme            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Electromagnetic energy - TM fields
         call grxy(xt,utma,ntav,15050,1,1 &
              ,chx,'     <U(TM)>   ','ueme            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Total thermal energy in plasma
         call grxy(xt,utha,ntav,15200,1,1 &
              ,chx,'     <Utherm>  ','uthm            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Electron thermal energy
         call grxy(xt,uthea,ntav,15300,1,1 &
              ,chx,'      <Uthe>   ','uthe            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Ion thermal energy
         call grxy(xt,uthia,ntav,15400,1,1 &
              ,chx,'      <Uthi>   ','uthi            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     Proton thermal energy
         call grxy(xt,uthpa,ntav,15450,1,1 &
              ,chx,'      <Uthp>   ','uthp            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Total escapee energy - RH boundary
         call grxy(xt,urha,ntav,15500,1,1 &
              ,chx,'      <Urhb>   ','urhb            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Total escapee energy - LH boundary
         call grxy(xt,ulha,ntav,15600,1,1 &
              ,chx,'      <Ulhb>   ','ulhb            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     Total plasma energy
         call grxy(xt,usysa,ntav,15700,1,1 &
              ,chx,'      <Usys>   ','usys            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     Total energy
         call grxy(xt,utota,ntav,15900,1,1 &
              ,chx,'     <Utot>    ','utot            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     absorption from rate of increase of <Utot>, <Uhote>, <Uion>
         icycdel=(xm1/tav+0.5)
         do i=1,ntav
            iu=i-icycdel
            if (uinc(iu).ne.0) then
               etat(i)=gam0*(utota(i+1)-utota(i))/dtav/uinc(iu)
               etah(i)=gam0*(urha(i+1)-urha(i))/dtav/uinc(iu)
               etai(i)=gam0*(uthia(i+1)-uthia(i))/dtav/uinc(iu)
            else
               etat(i)=0.
               etah(i)=0.
               etai(i)=0.
            endif
         end do


!     Absorption from rate of increase of total energy
         call grxy(xt,etat,ntav,16100,1,-1 &
              ,chx,'    dUtot/dt   ','abut            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     Absorption from rate of increase of total energy
         call grxy(xt,etah,ntav,16200,1,-1 &
              ,chx,'    dUhot/dt   ','abuh            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     Absorption from rate of increase of total energy
         call grxy(xt,etai,ntav,16300,1,-1 &
              ,chx,'   dUion/dt    ','abui            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)

!     max ion density
         call grxy(xt,rhoim,ntav,17000,1,1 &
              ,chx,'       <ni>    ','nihi            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)


!     position of critical ion density
         call grxy(xt,xcrni,ntav,17200,1,1 &
              ,chx,'       xci     ','xcni            ' &
              ,tmin,tmax,tstep,abmin,abmax,abstep)



      endif

      end







