
! ===========================================
!
!  Energy distribution of plasma ions
!    - Tihot, Tphot plots
!
! ===========================================

subroutine uidist
  use bopsvars
  implicit none

  integer i, ip, l, i0
  real*8 fmin, fpmin
  real*8 u,du, up, dup
  real*8 :: work1(0:nvx), work2(0:nvx)
  real*8 :: work1_p(0:nvx), work2_p(0:nvx)

  if (ni_tot.eq.0) return

  !  energy interval in mc^2: plot runs from -uionmax to +uionmax
! keep same energy scales for ions, protons (but different momenta)
!  uionmax=max(uimax,upmax)
  du=uimax/0.511/miome/nvx*2  ! normalized energy intervals
!  dup=(upmax-upmin)/0.511/mpome/nvx*2
  dup=upmax/0.511/mpome/nvx*2

  if (mod(itime,igr).eq.0) then
     if (itime.gt.0) then

        !  plot distn fn and re-zero; Te in KeV

        do i=1,nvx
           xv(i)=(i*du-du/2)*.511*miome-uimax
           work2(i) = fion(i)
        end do
        call grxy(xv,work2,nvx,6400+idc,1,2 &
             ,'      U(MeV)   ','     fi(U)     ','fuip'//ctime(1:12) )
! forward proton spectrum 
	if (np.ne.0) then
          do i=1,nvx
             xv(i)=(i*dup-dup/2)*.511*mpome-upmax
             work2(i) = fproton(i)
          end do
          call grxy(xv,work2,nvx,6600+idc,1,2 &
             ,'      U(MeV)   ','     fp(U)     ','fupp'//ctime(1:12) )
	endif
     endif

     do i=1,nvx
        fion(i)=0.
        fproton(i)=0.
     end do

  else if (mod(itime,igr).ge.igr-nuav) then
     do i=1,nvx
        work1(i)=0.
        work1_p(i)=0.
     end do

     !  relativistic thermal energy distn (lab frame)

     do l=1,ni_tot
        ip=l+ion1-1
        u=gam0*(gamma(ip)-vy0*uy(ip)) - 1.0
! heavy ions
	if (species(ip).eq.2) then
          if (ux(ip).gt.0) then
           i0=u/du+nvx/2
          else
           i0=-u/du+nvx/2
          endif
          if (i0.le.nvx .and. i0.ge.0) work1(i0)=work1(i0)+1

! protons
        else if (species(ip).eq.3) then
	! range -upmax -> upmax
           if (ux(ip).gt.0) then
            i0=u/dup+nvx/2
           else
            i0=-u/dup+nvx/2
           endif
!	  else   ! range 0-umax
!            i0=u/dup+1
!          endif
          if (i0.le.nvx .and. i0.ge.0) work1_p(i0)=work1_p(i0)+1

        endif
     end do


     fmin=0.1
     fpmin=0.1

     do i=1,nvx
        if (work1(i).ne.0) then
           fion(i)=fion(i)+work1(i)/nuav
        else
           fion(i)=fmin
        endif

        if (work1_p(i).ne.0) then
           fproton(i)=fproton(i)+work1_p(i)/nuav
        else
           fproton(i)=fpmin
        endif
     end do

  endif
end subroutine uidist
