
!     ==========

subroutine splot
  use bopsvars
  implicit none

  real*8 uxmi, uzma, uzmi, uymi, xlimit
  integer i, iplot


  do iplot=1,nsp
     if (itime.eq.isp(iplot)) then
        !  find max/min of momenta in lab/boost frame
        uxma=ux(1)
        uxmi=ux(1)
        uzma=uz(1)
        uzmi=uz(1)

        if (ioboost.eq.1) then
           uyma=uy(1)
        else
           uyma=gam0*( uy(1) - vy0*gamma(1) )
        endif

        uymi=uyma

        do i=1,ne
           uux(i) = ux(i)
           uuz(i) = uz(i)

           if (ioboost.eq.1) then
              !  boost frame
              uuy(i) = uy(i)
              xnlab(i) = xn(i)
           else
              !  lab frame
              uuy(i) = gam0*( uy(i) - vy0*gamma(i) )
              xnlab(i) = gam0*xn(i)
           endif

           uxma=max(uxma,uux(i))
           uxmi=min(uxmi,uux(i))
           uyma=max(uyma,uuy(i))
           uymi=min(uymi,uuy(i))
           uzma=max(uzma,uuz(i))
           uzmi=min(uzmi,uuz(i))
        end do
        xlimit = (xsol+xl)/2*gam0
        !  x, ux, uy plot in rayshade format
        call gr3d(xnlab,uux,uuz,1,ne,1,iplot &
             ,'movie/elec'//ctime(1:6) &
             ,0,xllab,-1.,1.,-1.,1.,xlimit,1.5*vte)

     endif
  end do
end subroutine splot
