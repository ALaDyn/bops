      program od2gle

c  (c)  Paul Gibbon, May 1995

c  Converts ODPLOT graphics data to GLE format
c  
c  Now operates entirely in data directory
c  Fetches id files from base directory cbase

      real x(500000),y(500000)

      character*15 chx,chy,ctitle*30,cbase*80,cdir*30,cyn*1,cidfull*80
     :,climx*1,chead*1,climy*1,cval*19,cdate*24,ccode*80
      character*40 cp,csym(0:9)
     :,cfile*80,csnap*2,cmem*80,cpage*2,cid*80
      logical lexist
      logical :: debug=.false.
      integer :: id, mid
      data csym/'dot','cross','wcircle','wsquare','wtriangle'
     :         ,'fcircle','fsquare','ftriangle','fdiamond','dot'/

      write(6,'(a)') 'od2gle graphics postprocessor (c) P. Gibbon, 1994'
      write(6,*) 'Base directory (absolute path):'
      read(5,'(a)') cbase
      write(6,'(a)') cbase
      write(6,*) 'Run name:'
      read(5,'(a)') ccode
      write(6,'(a)') ccode
      write(6,*) 'File with ID list'
      read(5,'(a)') cid
      write(6,'(a)') cid
      write(6,*) '# plots per page:'
      read(5,*) nbox
      write(6,*) 'header page (y/n)'
      read(5,'(a)') chead
      itypo=1

      ibox=1
      ipage=1
      lcode=lench(ccode)
      lcbase=lench(cbase)
      lcid=lench(cid)
      cidfull = cbase(1:lcbase)//'/'//cid(1:lcid)//'.id' 
      lcidfull = lench(cidfull)
!      write (6,'(a)') cidfull(1:lcidfull)
      
c  fortran streams for read file (header and data)
      iread=50
      idata=51
      ihead=52

c  id list file 
      irlist=40

c  fortran stream for cplot command file
      iwrcom=60
      iwrhead=61

      if (chead.eq.'y') then
c  open the header file
        open(ihead,file='header')
c  open the command file
        open(iwrhead,file='header.gle')
        write (iwrhead,'(a/a/a/a)') 
     :  '! header','size 18 26 ','begin translate .5 1',' box 17 20'
c  run date/time
c        read(ihead,'(a/a/a)') cyn,cdate,cyn
c        write (iwrhead,'(a/a)') 'set font rmi','set hei 0.4'
c        write (iwrhead,'(a)') 'amove 2. 17.'
c        write (iwrhead,'(a)') 'text '//cdate
        write (iwrhead,'(a/a)') 'set font rm','set hei 0.55'

        rx = 2.
        ry = 18.
        
   5    read(ihead,'(a)',end=6) cval
        if (cval(1:6).eq.'      ') goto 6
        write (iwrhead,'(a,2f13.2/a)') 'amove ',rx,ry,'set just left'
        write (iwrhead,'(a)') 'text '//cval(1:6)
        write (iwrhead,'(a/a)') 'rmove 5 0','set just right'
        write (iwrhead,'(a)') 'text '//cval(8:19)


        ry=ry-0.8
        if (ry.le.2) then
          rx=rx+8.
          ry=18.
        endif
        goto 5
   6    continue
      write (iwrhead,'(a)') 'end translate'
      close(iwrhead)
      close(ihead)
      endif

       


c  open the read file
      open(iread,file='bops.oddata')


c  open id list file
      open(irlist,file=cidfull(1:lcidfull))
      

      if (nbox.eq.9) then
        irow=3
      else if (nbox.eq.25) then
        irow=5
      else if (nbox.eq.30 .or. nbox.eq.16) then
	irow=4
      else if (nbox.eq.10.or.nbox.eq.4.or.nbox.eq.8) then
        irow=2
      else 
        irow=1
      endif
 
 
c  open the command file
      open(iwrcom,file=ccode(1:lcode)//'-1.gle')

  10  continue

c  page string
      call chr(1.*ipage,0,cp,lcpage)
      cpage=cp(1:lcpage)

c  read next id
      read (irlist,*,end=999) idr
      rewind (iread)  ! rewind oddata file

c  skip to new line if end of sequence:  ID=0 in .id file
      if (idr.eq.0) then
	ibinc=mod(ibox-1,irow)
        if (ibinc.ne.0) then
          ibox=ibox+irow-ibinc
        endif
        goto 100

c  new page:  ID=-999
      else if (idr.eq.-999) then
        ibox=nbox+1
        goto 100
      endif

c  read plot data in odplot-format
!  20  read (iread,'(i6)',end=10) id
  20  read (iread,'(2i6,4i4,2a15,a)',end=10) id,n,itype,idash,isym,ispage
     :,chx,chy,ctitle
!  20  read (iread,*) id!

      if (debug) write (6,'(i6,a6)') id,ctitle(1:4)

!      read (iread,'(2a15,a16)') chx,chy,ctitle
!  20  read (iread,101) id,n,itype,idash,isym,ispage
!     :,chx,chy,ctitle
      if (id.lt.0) then
         if (debug) write (6,*) 'reading ranges'
         read (iread,'(6(1pe12.3))') xmin,xmax,dx,ymin,ymax,dy
      endif
 101  format (2i6,4i4,2a15,a)

c  check if id found
      if (abs(id).ne.abs(idr)) then
         if (debug) write(6,*) 'ID ',idr,'not found in bops.oddata' 
         goto 20
      else
         if (debug) write(6,*) 'ID ',idr,'found - making graph' 
      endif
 
c  ID found, so read in data
c  Open data file gibbon/run#/cmember
c   - ensure directory exists before running postprocessor

c  convert snapshot id to character
c   - convention is mod(id,10)=1,2,3 etc for snapshots, 0 for time histories
 
!      call chr(1.0*abs(mod(id,10)),0,csnap,l1)
!  2-digit timestamp
      mid=abs(id)
      csnap(2:2) = achar(mod(mid,10) + 48)  
      csnap(1:1) = achar(mod(mid/10,10) + 48)
      cmem=ctitle(1:4)
      cfile = cmem(1:lench(cmem))//csnap(1:2)//'.xy'

      inquire(file=cfile,exist=lexist)

      if (lexist) then
        open(idata,file=cfile)
c  old data extension
c        cfile=cmem(1:lench(cmem))//csnap(1:1)//'.xydata'
c        open(idata,file=cfile)
         if (debug) write (*,*) 'Reading data for ', cfile
         rewind(idata)
         read(idata,*,end=50) (x(i),y(i),i=1,n)
         goto 55
 50      write (*,*) 'End of data ',i-1,' out of ',n,' found'

        
 55      close(idata)
      else
         write (*,*) 'Data file',cfile,' not found'
         goto 10
      endif


c  check for constant data
      ymi=y(1)
      yma=y(1)

      do i=1,n
        ymi=amin1(ymi,y(i))
        yma=amax1(yma,y(i))
      end do

      if (id.gt.0 .and. yma.eq.ymi) then
        id=-id
        ymax=2*yma
        if (ymax.eq.0) ymax=1.
        if (ymax.gt.0) then
         ymin=yma/2.
        else if (ymax.lt.0) then
         ymin=yma*2
        endif
        xmax=x(n)
        xmin=x(1)
        if (xmax.eq.xmin) xmax=xmax+1.
      endif

c  check for zeros if lin-log
      if (cyn.eq.'y' .and. itypo.ne.itype) then
        do i=1,n
          if (itypo.eq.2) then
            y(i)=amax1(y(i),1.e-20)
          endif
        end do
        itype=itypo
      endif

      if (climy.eq.'y') then
        ymin=ymino
        ymax=ymaxo
      endif

      if (climx.eq.'y') then
        xmin=xmino
        xmax=xmaxo
      endif

c  write to cplot dataset


c  initialise graphics filter - skip if plot on same page
      if (ibox.eq.1) then  
        write (6,'(a)') 'Page '//cpage
	if (nbox.eq.25.or.nbox.eq.9) then
          write (iwrcom,1005) 
     :  '! new page','size 26 18 ','set font rm'
	else
          write (iwrcom,1005) 
     :  '! new page','size 18 26 ','set font rm'
	endif
1005   format(a/a/a)
      endif  

c  echo file name to diagnostic output
      write (6,'(a,i8)') cmem(1:4)//csnap(1:2),id



 
c  calculate window positions

c  single plot

      if (nbox.eq.1) then
        xo=3
        yo=10
        xs=10
        ys=7
        hlab = 0.6
        hax = 0.5
        wline = 0.02

      else if (nbox.eq.2) then
        j=2-ibox
        xo=5
        yo=5+8*j
        xs=8
        ys=5
        hlab = 0.6
        hax = 0.5
        wline = 0.02

      else if (nbox.eq.3) then
        j=3-ibox
        xo=2
        yo=1+8*j
        xs=14
        ys=8
        hlab = 0.5
        hax = 0.5
        wline = 0.025

      else if (nbox.eq.10) then
        j=4-(ibox-1)/2
        i=mod(ibox+1,2)
        xo=1.+8.5*i
        yo=.5+5.2*j
        xs=7.5
        ys=4.5
        hlab = 0.4
        wline = 0.01
        hax = 0.3

      else if (nbox.eq.8) then
        j=3-(ibox-1)/2
        i=mod(ibox+1,2)
        xo=1+8.5*i
        yo=1+6*j
        xs=9
        ys=5.5
        hlab = 0.4
        wline = 0.02
        hax = 0.4


c  9 per page

      else if (nbox.eq.9) then
        i=mod(ibox+2,3)
        j=2-(ibox-1)/3
        xo=1+7*i
        yo=1+5.5*j
        xs=7
        ys=5
        hlab = 0.3
        hax = 0.25
        wline = 0.01

c  25 per page

      else if (nbox.eq.25) then
        i=mod(ibox+4,5)
        j=4-(ibox-1)/5
        xo=1+4.5*i
        yo=1+3.5*j
        xs=4.2
        ys=3.
        hlab = 0.25
        hax = 0.16
        wline = 0.01

c  30 per page

      else 
        i=mod(ibox+3,4)
        j=6-(ibox-1)/4
        xo=2+4.*i
        yo=1+3.4*j
        xs=3.5
        ys=2.5
        hlab = 0.25
        hax = 0.16
        wline = 0.01
      endif

c  shift origin

      write (iwrcom,'(a,2(1pe13.4))') 'begin translate',xo,yo
      write (iwrcom,'(a)') 'begin graph'
      write (iwrcom,'(a,2(1pe13.4))') 'size',xs,ys

c  data fetch and font commands

      write (iwrcom,'(a)') 'data '//cfile(1:lench(cfile))


c  box limits:
c   axis scaling is manually chosen if id<0

      write (iwrcom,'(a)') 'nobox'

      if (id.lt.0) then
	write (iwrcom,1009) 'xaxis min ',xmin,' max ',xmax
     :                     ,'yaxis min ',ymin,' max ',ymax
 1009   format(a,1pe13.4,a,1pe13.4/a,1pe13.4,a,1pe13.4) 

      else if (climx.eq.'y') then
	write (iwrcom,1020) 'xaxis min ',xmin,' max ',xmax 

      else if (climy.eq.'y') then
	write (iwrcom,1020) 'yaxis min ',ymin,' max ',ymax
 1020   format(a,1pe13.4,a,1pe13.4) 
      endif


c  Axis labels

      write (iwrcom,'(a,1pe13.4)') 'xaxis hei ',hax
      write (iwrcom,'(a,1pe13.4)') 'yaxis hei ',hax
      write (iwrcom,'(a,1pe13.4)') 'xlabels hei ',hlab
      write (iwrcom,'(a,1pe13.4)') 'ylabels hei ',hlab



c  Axis type (other than lin-lin)

      if (itype.eq.2) then
c  lin-log
        write (iwrcom,'(a)') 'yaxis log'
      else if (itype.eq.3) then
c  log-lin
        write (iwrcom,'(a)') 'xaxis log'
      else if (itype.eq.4) then
c  log-log
        write (iwrcom,'(a/a)') 'xaxis log','yaxis log'
      endif

c  Tick marks
      write (iwrcom,'(a)') 'xaxis nticks 5'

c  length
      
      write (iwrcom,'(a/a)') 'xticks length .07','yticks length .07' 

c  Text labels and units

      write (iwrcom,'(a7,a17,a13,a20/2a/2a)') 'title ','"'//cfile(1:10)
     ://', '
     :,ctitle(5:16)//'"'
     :,' hei 0.25 font psh'
     :                           ,'xtitle ','"'//chx//'"'
     :                           ,'ytitle ','"'//chy//'"'


c  Line style
      if (idash.gt.0) then
        write (iwrcom,1001) 'd1 line'
     :                    ,' lwidth ',wline
     :                    ,' lstyle ',idash
      endif

      isym = mod(isym,10)
      if (isym.eq.9) then
	 write (iwrcom,'(a)') 'd1 marker dot msize 0.05'
      else if (isym.ne.0) then
        write (iwrcom,'(a,f6.1)') 'd1 marker '//csym(isym)
     ://'msize ',0.2
      endif
1001  format(a,a,1pe13.4,a,i6)


      
      write (iwrcom,'(a)') 'end graph'
      write (iwrcom,'(a///)') 'end translate'

c  increment box #
      ibox=ibox+1

 100  continue

c  new page

      if (ibox.gt.nbox) then
c  close command file
        close(iwrcom)
        ibox=1
        ipage=ipage+1
c  page string
        call chr(1.*ipage,0,cp,lcpage)
        cpage=cp(1:lcpage)
c  open new command file
        open(iwrcom,file=ccode(1:lcode)//'-'//cpage(1:lcpage)//'.gle')
      endif

      goto 10
 999  write (6,*) 'End of data'

      close(iwrcom)
      stop
 900  print *,'Error in data at i=',i
      print *,x(i),y(i)
      end

c ==============

c  converts real to character

      subroutine chr(z,ndp,ch,l)
      character ch*40,chnum(0:9)*1
      data chnum/'0','1','2','3','4','5','6','7','8','9'/
      ch=' '
      ia=log10(amax1(1.,abs(z)))
      e=float(ia)
      ie=1
c  turn number into integer
      do 50 i=1,ndp
  50    ie=ie*10
      iz=int(abs(z)*ie+0.5)
      is=sign(1.,z)
      l=0
c  negative number
      if (is.lt.0) then
	l=l+1
	ch(l:l)='-'
      endif
c  make ie as big as integer number
      do 60 i=1,ia
  60    ie=ie*10
c  do number
      do 100 idp=ia,0,-1
	l=l+1
	idigit=iz/ie
	ch(l:l)=chnum(idigit)
	iz=iz-idigit*ie
  100   ie=ie/10
c  decimal point
      if (ndp.gt.0) then
	l=l+1
	ch(l:l)='.'
      endif
c  do decimal places
      do 200 idp=-1,-ndp,-1
	l=l+1
	idigit=iz/ie
	ch(l:l)=chnum(idigit)
	iz=iz-idigit*ie
 200    ie=ie/10
      end

      include 'lench.f'
