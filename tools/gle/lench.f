
      function lench(c1)
      character*80 c1,c0*1
      l=80
  10  l=l-1
      c0=c1(l:l)
      if ((c0.eq.' ').and.(l.gt.1)) goto 10
      if (c0.eq.' ') l=0
      lench=l
      end

