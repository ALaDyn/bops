!  =================================
!
!      open/close particle o/p files
!
!  =================================

           subroutine fileop
!  reflected laser amplitude
      open(70,file='ey_back.t')
      open(71,file='ez_back.t')

!  momenta files
      open(80,file='e_mom_lhb.dat')
      open(81,file='e_mom_rhb.dat')
      open(82,file='ion_mom_lhb.dat')

!  Output files
        open(15,file='bops.out')
        open(50,file='bops.oddata')
        open(20,file='bops.header')
        open(90,file='exit_energies.dat')
           end

          subroutine fileclo
      parameter(nfile=9)
      integer ifile(nfile)
      data ifile/15,20,50,70,71,80,81,82,90/

      do i=1,nfile
       close(ifile(i))
      end do
           end
