      parameter (npm=10000,nbm=10,nms=10,nidm=20,nmovm=200)
      implicit logical(l)
      character*12 ci*80,c1*1,c2*2,c3*3,cw*80
      real nsiz
      common/odarr/ xori(nbm),   yori(nbm),   axl(2,nbm),
     :   xd(npm),   yd(npm),     xp(npm),      yp(npm)
     :,  xg(2),     yg(2)
     :,space(nms), smark(nms),ra(2,nbm), vor(2,nbm)
     :,  xmo(npm,nbm), ymo(npm,nbm)

      common/odvar/ xl,yl,cht,ibc,np,xmg,ymg,xsg,ysg,lmanex,lmanc
     :           ,lrnorm,nbox,ityp,isym,idash,icolor,lag,lbox,lhc
     :           ,lautid,lbadid,nlab,ndash,itp,lora,lprl,nsiz,ndpo
     :           ,ldpo,gid,igraph,lautax,ibacol,xaxgap,yaxgap
     :           ,xa1,xa2,dax,ya1,ya2,day,tics,tick,iaxcol
      common/llog/ lbadlog
      common/cv/ ci,c1,c2,c3,cw

