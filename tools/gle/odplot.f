      program odsun

c  (c)  Paul Gibbon, June 1990

c  adapted for ibm 3090 CMS, March 1992

c  uses GDDM GKS library

c  IBM 3179 terminal is   gopwk(1,istream,1)
c  gdf metfile is         gopwk(1,istream,5)
c  gksm metafile is       gopwk(1,istream,3)


      real x(10000),y(10000)
      logical loride,lrnorm,lfound,lbox,lhc,lrw,lmanu,lsegmt
      character chs*80,copt*2,carg*20,chform*5,chfout*6,chfout2*9
     :         ,clabx*15,claby*15,ctxt*20
      chform='(a)'
      chfout='(1x,a)'
      chfout2='(1x,a,i3)'
      ncall=0
      ifile=20
      lsegmt=.false.
      lrnorm=.false.
      loride=.false.
      lrw=.true.
      lmanu=.false.
      rewind (ifile)
      id=999999
      icmd=5
c
      nplot=1
c  lin-lin axes
      itype=1
c  line type
      idash=1
c  plotting symbol
      isym=0
c  new axes for each graph
      ignum=1
      iseg=1
c  # boxes
      nbox=1
      ibox=0
c  movie
      n1=1
      n2=1
      icol=2
      call grin
      call setup
      call setws(1)
      call grocgm(nplot)

c  create header segment
c     call gcrsg(100)
c     call headpg
c     call gclsg

      lbox=.false.
      lhc=.false.
c     call gmsg(1,'ODPLOT  Paul Gibbon 1990')
      call help

   10   continue
	carg=' '
c Read in command
	if (icmd.eq.5) write (6,chfout) 'Command:'
c       call grqst(1,1,jstat,lch,chs)
	jstat=1
	read (icmd,chform,end=999,err=999) chs
	call csiz(chs,lch)
	copt=chs(1:2)
	carg=chs(3:lch)
c        write (*,'(3(1x,a/),i8)') chs,copt,carg,lch

c  store command
	if (copt.ne.'ex') write (40,'(a)') chs
c  echo to screen if read from file
	if (icmd.ne.5) write (6,chfout) chs

c  List data ids
	if (copt.eq.'ld'.or.jstat.ne.1) then
	  call idop(ifile)

c  Batch mode - preview everything
	else if (copt.eq.'pv') then
	  istart=ich(carg)
	  call bat(ifile,nbox,1,istart,nplot)

c  Batch mode - plot everything with hardcopy |
	else if (copt.eq.'ba') then
	  istart=ich(carg)
	  call bat(ifile,nbox,2,istart,nplot)

c  Read commands from file
	else if (copt.eq.'cm') then
	  icmd=ich(carg)
	  if (icmd.lt.10) goto 99
	  if (icmd.ne.5) write (6,102) 'Reading command list',icmd
	  rewind(icmd)
 102      format(1x,a25,i4)

c  Make hardcopy
	else if (copt.eq.'hc') then
	  call grsgwk(nplot+1)
	  call grccgm(nplot)
	  nplot=nplot+1
	  iseg=iseg+1
	  call grocgm(nplot)
	  write (6,chfout2) 'Hardcopy ',nplot-1
	  call gclrwk(1,0)
	  ibox=0

c  Show graph
	else if (copt.eq.'sh') then
c  close segment and display
	  call grsgwk(1)
c  pause
	  call grqch(1,2,istat,ichar)

c  Clear graphs and reset
	else if (copt.eq.'cl') then
c  clear segments
	  call gclrwk(1,0)
	  do i=1,iseg
	    call gdsg(i)
	  end do
	  itype=1
	  idash=1
	  isym=0
	  icol=2
	  iseg=1
	  lsegmt=.false.
	  write (6,chfout) 'cleared'
	  rewind (30)

c Exit
	else if (copt.eq.'ex') then
c  close gks
	  write (40,'(a)') 'cm5'
	  call grccgm(nplot)
	  call grst
	  call grprint(nplot-1)
	  stop

c  Select axis type: lin-lin, lin-log, log-lin, log-log
	else if (copt.eq.'ax') then
	  itype=ich(carg)

c  Symbol
	else if (copt.eq.'sy') then
	  isym=ich(carg)
	  if (isym.eq.9) idash=0

c  Line-type
	else if (copt.eq.'da') then
	  idash=ich(carg)

c  Colour
	else if (copt.eq.'co') then
	  icol=ich(carg)

c  Line Options overide toggle
	else if (copt.eq.'or') then
	  loride=.not.loride
	  if (loride) write (6,chfout) 'Manual line options'
	  if (.not.loride) write (6,chfout) 'Auto line options'

c  Number of boxes per page
	else if (copt.eq.'nb') then
	  nbox=ich(carg)
	  ibox=0

c  Box number for next graph
	else if (copt.eq.'ib') then
	  ibox=ich(carg)-1

c  File stream to read data from (default=20)
	else if (copt.eq.'fn'.or.copt.eq.'fi') then
	  ifile=ich(carg)
	  ncall=0
	  if (ifile.gt.99) then
	    write (6,chfout) 'File # no good - choose one < 100'
	    ifile=20
	  endif

c  Rewind switch
	else if (copt.eq.'rw') then
	  lrw=.not.lrw
	  if (lrw) write (6,chfout) 'Rewind on'
	  if (.not.lrw) write (6,chfout) 'Rewind off'

c  Text label
	else if (copt.eq.'tx') then
	  write (6,chfout) 'text,x,y,theta,allign'
	  read (5,'(a)') ctxt
	  read (5,*) xpos,ypos,ang,iall
	  call lab(xpos,ypos,ang,0.4,ctxt,20,iall)

c  Manual axis scaling
	else if (copt.eq.'ms') then
	  lmanu=.true.
	  jstat=1
c         call gmsg(1,'xmin:')
c         call grqvl(1,1,jstat,xmin)
c         call gmsg(1,'xmax:')
c         call grqvl(1,1,jstat,xmax)
c         call gmsg(1,'dx:')
c         call grqvl(1,1,jstat,dx)
c         call gmsg(1,'ymin:')
c         call grqvl(1,1,jstat,ymin)
c         call gmsg(1,'ymax:')
c         call grqvl(1,1,jstat,ymax)
c         call gmsg(1,'dy:')
c         call grqvl(1,1,jstat,dy)
	  write (6,chfout) 'Manual scaling:'
	  write (6,chfout) 'xmin,xmax,dx,ymin,ymax,dy:'
	  read (icmd,*,end=999,err=999) xmin,xmax,dx,ymin,ymax,dy
	  write (40,*) xmin,xmax,dx,ymin,ymax,dy
c         if (jstat.ne.1) lmanu=.false.

c  Automatic axis scaling
	else if (copt.eq.'as') then
	  write (6,chfout) 'Auto scaling'
	  lmanu=.false.

	else if (copt.eq.'su') then
	  write (6,chfout) 'Superposition on'
	  ignum=3

	else if (copt.eq.'so') then
	  write (6,chfout) 'Superposition off'
	  ignum=1

	else if (copt.eq.'rn') then
	  write (6,chfout) 'Renormalisation: x,y'
	  lrnorm=.true.
	  read (icmd,*) xnorm,ynorm

c  Help page
	else if (copt.eq.'h') then
	  call help

c  Draw graph on new (vn) or old (vo) axes
	else if (copt.eq.'vn'.or.copt.eq.'vo') then
c  create graph segment and set flag true
	  n1=1
	  idq=ich(carg)
	  if (lrw) rewind ifile
	  lfound=.false.
c  data heading
  15        read (ifile,'(2i6,4i4,2a15)',end=20)
     :           idr,n,ityp,ida,isy,ispage,clabx,claby
	    id=abs(idr)

c  Auto line options
	    if (.not.loride) then
	      itype=ityp
	      idash=ida
	      isym=isy
	    endif

	    if (id.eq.idq) then
	      lfound=.true.
	    endif
c  preset scaling
	    if (idr.lt.0) read (ifile,*,end=20)
     : xmi,xma,dax,ymi,yma,day

c  data
	    read (ifile,*)  (x(i),y(i),i=1,n)
	    if (.not.lfound)  goto 15
  20      if (.not.lfound) then
	    write (6,chfout) 'Dataset not found'
	    goto 99
	  endif

c  manual overide
	  if (lmanu) then
	    xmi=xmin
	    xma=xmax
	    dax=dx
	    ymi=ymin
	    yma=ymax
	    day=dy
	  endif

	  if (copt.eq.'vo'.and.lbox) then
	    if (ignum.ne.3) ignum=2
	    ibox=ibo
	  else if (ignum.eq.3) then
	    ibox=ibo
	  else
	    ignum=1
	    if (ibox.eq.0)  then
c             call gclrwk(1,0)
c             call grccgm(nplot)
c             call grocgm(nplot)
	    endif
	    lbox=.true.
	  endif

c renormalisation
	  if (lrnorm) then
	    do i=1,n
	      x(i)=x(i)*xnorm
	      y(i)=y(i)*ynorm
	      lrnorm=.false.
	    end do
	  endif

c check for zeros if log axis selected
	  if (itype.eq.2) then
	    yb=0.
	    ybn0=1.e24
	    yt=0.
	    do i=1,n
	      yb=amin1(yb,y(i))
	      if (y(i).ne.0) ybn0=amin1(ybn0,y(i))
	      yt=amax1(yt,y(i))
	    end do
	    if (yb.eq.0) then
	      yb=ybn0
	      do i=1,n
		y(i)=amax1(yb,y(i))
	      end do
	    endif
	  endif

	  call gcrsg(iseg)
	  if (idr.lt.0.or.lmanu) then
	    call odpm(x,y,n,id,itype,idash,isym,ignum,icol,nbox,ibox+1
     :,xmi,xma,dax,ymi,yma,day,clabx,claby)
	  else
	    call odp(x,y,n,id,itype,idash,isym,ignum,icol,nbox,ibox+1
     :,clabx,claby)
c            call movi(x,y,n,id,itype,idash,isym,ignum,icol,nbox,ibox+1
c     :,xmin,xmax,dx,ymin,ymax,dy)
	  endif
	  call gclsg
	  iseg=iseg+1
	  ibo=ibox
	  ibox=mod(ibox+1,nbox)

c  Movie mode
	else if (copt.eq.'mv') then
	  nmove=ich(carg)
	  n2=n1+nmove-1
	  call movie(ibo+1,n1,n2,icol)
	  n1=n2+1

c  Time-slice mode
	else if (copt.eq.'3d') then
c         call gclrwk(1,0)
	  call grccgm(nplot)
	  call grocgm(nplot)
	  rewind (ifile+1)
	  read(ifile+1,*) n,ngraph,yrange
	  rewind ifile
	  print *,'points ',n,' no. graphs ',ngraph,' yrange ',yrange
c         read (ifile,*,end=250) ((f(i,j),i=1,n),j=1,ngraph)
	  do 250 ignum=1,ngraph
	    j1=ngraph-ignum+1
	    do 300 i=1,n
c             y(i)=f(i,j1)+yrange/ngraph*j1
 300          x(i)=i
c           call movi(x,y,n,id,itype,-1,isym,ignum,icol,1,1
c    :,0.,float(n),n/5.,0.,1.5*yrange,yrange/5)
 250      continue
	else
	  write (6,chfout) 'Unknown command'
	endif
  99  goto 10
 999  rewind(5)
      goto 10
 101  format(a)
      end
c
c     =========
c
      subroutine idop(ifile)
      real x(10000),y(10000),xmin(200),xmax(200),ymin(200),ymax(200)
      integer idlist(200),ncall
      character*15 cx(200),cy(200),clabx*15,claby*15
      write (6,103) 'ID list:','id','xmin','xmax','ymin','ymax'
 103  format(1x,a/a5,4a10)
c     if (ncall.gt.0) goto 99
      ncall=ncall+1
      rewind ifile
      nid=1
  10  read (ifile,'(2i6,4i4,2a15)',end=99)
     : id,n,i1,i2,i3,i4,clabx,claby
	if (id.lt.0) read(ifile,*) xd1,xd2,xd3,yd1,yd2,yd3
	read (ifile,*) (x(i),y(i),i=1,n)
	cx(nid)=clabx
	cy(nid)=claby
	idlist(nid)=id
	xmax(nid)=x(1)
	xmin(nid)=x(1)
	ymax(nid)=y(1)
	ymin(nid)=y(1)
	do i=1,n
	  xmax(nid)=amax1(xmax(nid),x(i))
	  xmin(nid)=amin1(xmin(nid),x(i))
	  ymax(nid)=amax1(ymax(nid),y(i))
	  ymin(nid)=amin1(ymin(nid),y(i))
	end do
	nid=nid+1
      goto 10
  99  write (6,101) (i,idlist(i),xmin(i),xmax(i),ymin(i),ymax(i)
     :,cy(i),cx(i),i=1,nid-1)
 101  format ((1x,i3,i5,4(1pe10.2),2a15))
 102  format (1x,3i6,6(1pe10.2))
      end
c
c     =========
c
      subroutine bat(ifile,nbox,ihard,istart,nplot)
      real x(10000),y(10000)
      character*15 clabx*15,claby*15,chfout*12
      write (6,103) 'Batch processing ..'
 103  format(1x,a)
      rewind ifile
      icol=2
      nid=0
      ibox=1
      iseg=1
      chfout='(1x,a,1x,i4)'
  10  nid=nid+1

c  create segment
      call gcrsg(iseg)
       read (ifile,'(2i6,4i4,2a15)',end=20)
     :    id,n,itype,idash,isym,ispage,clabx,claby
       if (id.lt.0) read (ifile,*) xmin,xmax,dx,ymin,ymax,dy
       read (ifile,*,end=20) (x(i),y(i),i=1,n)

c  skip graphs up to start #
       if (nid.ge.istart) then

c  check limits for log plots
	xdmi=x(1)
	ydmi=y(1)
	do 100 i=1,n
	  xdmi=amin1(xdmi,x(i))
 100      ydmi=amin1(ydmi,y(i))
	if (itype.eq.3 .and. xdmi.le.0) itype=1
	if (itype.eq.2 .and. ydmi.le.0) itype=1
	if (itype.eq.4 .and. (ydmi.le.0. .or. xdmi.le.0.)) itype=1

c  same axes switch
c       if (ispage.eq.2 .and. ibox.gt.1) ibox=ibox-1

c  create graph segment
	write (6,'(3i6,a15)') nid,id,ibox,claby
	if (id.lt.0) then
	  call odpm(x,y,n,id,itype,idash,isym,ispage,icol,nbox,ibox
     :,xmin,xmax,dx,ymin,ymax,dy,clabx,claby)
	else
	  call odp(x,y,n,id,itype,idash,isym,ignum,icol,nbox,ibox
     :,clabx,claby)
	endif

c  close segment
	call gclsg
	ibox=ibox+1
	if (ibox.gt.nbox) then
c  close old segment and create new one
	  if (ihard.eq.2) then
	    call grsgwk(nplot+1)
	    call grccgm(nplot)
	    write (6,chfout) 'Hardcopy ',nplot
	    nplot=nplot+1
	    call grocgm(nplot)
	  else
	    call grsgwk(1)
c  pause
	    call grqch(1,2,istat,ichar)
	    if (ichar.eq.5) then
	      call grsgwk(nplot+1)
	      call grccgm(nplot)
	      write (6,chfout) 'Hardcopy ',nplot
	      nplot=nplot+1
	      call grocgm(nplot)
	    endif
	    call gclrwk(1,0)
	    if (ichar.eq.3) goto 20

	  endif
	  ibox=1
	endif
       endif
       iseg=iseg+1
      goto 10
  20  continue
      if (ibox.ne.1) then
	if (ihard.eq.2) then
	  call grsgwk(nplot+1)
	  call grccgm(nplot)
	  write (6,chfout) 'Hardcopy ',nplot
	  nplot=nplot+1
	  call grocgm(nplot)
	else
	  call grsgwk(1)
	  call grqch(1,2,istat,ichar)
	  call gclrwk(1,0)
	endif
      endif
      end
c
      integer function ich(cw)
      character cw*20
      rewind 3
      write (3,'(a)') cw
      backspace 3
      read (3,*,end=99) i
      ich=i
      return
 99   ich=1
      end
c
      subroutine help
      write (6,101)
 101  format(
     :1x,'h        -  help'/
     :1x,'vn(id)   -  plot (id) on new axes'/
     :1x,'vo(id)   -  plot (id) on old axes'/
     :1x,'ld       -  list data set'/
     :1x,'fn(m)    -  change data set to fort.(m)'/
     :1x,'hc       -  make hardcopy'/
     :1x,'ax(m)    -  axis type  1=lin-lin,'
     :,'  2=lin-log, 3=log-lin, 4=log-log'/
     :1x,'da(m)    -  line type (1-5)'/
     :1x,'sy(m)    -  symbol type (1-9)'/
     :1x,'co(m)    -  line colour'/
     :1x,'nb(m)    -  number of boxes per page'/
     :1x,'ib(m)    -  plot on box (m)'/
     :1x,'ms       -  manual scaling'/
     :1x,'as       -  automatic scaling'/
     :1x,'fn(m)    -  use data on file fort.m'/
     :1x,'cm(m)    -  use plot commands from file fort.m'/
     :1x,'ba       -  batch mode: hardcopy everything in datafile'/
     :1x,'pv       -  preview everything in datafile'/
     :1x,'rw       -  rewind on/off'/
     :1x,'rn       -  renormalise x,y axes'/
     :1x,'or       -  toggle to override graph options'/
     :1x,'ex       -  exit'
     :)
      end
c
      subroutine grprint(nplot)
	rewind 3
	write(3,*) nplot
      end

c
      subroutine headpg
      character*8 cword
      open(30,file='/od header *')
      j=1
      do while (j.lt.20)
	read(30,'(a8,f12.4,i6)',end=999) cword,x,ndp
	call hnum(cword,x,ndp)
	j=j+1
      end do
 999  return
      end
c
c  ====================
c
      subroutine hnum(cword,z,ndp)
      include 'odgks common'
      character cword*8
      save rx,ry
      data rx,ry/0.2,0.95/
      xw=rx*xl
      yw=ry*yl
      height=xl/40.
      call lab(xw,yw,0.,height,cword,8,3)
      call numbr(xw+cht,yw,height,z,1,ndp)
      ry=ry-0.05
      if (ry.le.0.86) then
	rx=rx+0.3
	ry=0.95
      endif
      end
c
c  ====================
c     include 'odgks for *'
