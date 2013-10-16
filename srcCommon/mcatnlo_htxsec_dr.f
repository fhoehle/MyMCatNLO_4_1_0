      subroutine frealht(s,x,yi,cth2,tk,uk,q1q,q2q,xinv,
     #                   jproc,xmatout)
c Returns the real matrix contribution, times the FKS damping factor: 
c    xii**2*(1-yi**2)=4*tk*uk/s**2    
c Only initial-state kinematics is used here; final-state collinear
c singularities will not occur for Ht. The normalization is such that
c    dsigma_real = g_S^2 (|a^2|+|b^2|) |V_tj|^2 freal dphi_3
c with g_S^2, |a^2|+|b^2| and CKM factors inserted in the main code. 
c
c The matrix elements are passed through the following arrays
c    xmat(out/in)(i):
c    jproc=1 --> gg processes
c    jproc=2 --> q(bar)q(bar) processes
c    jproc=3 --> gq(bar),q(bar)g processes
c The meaning of i(==idr) depends on jproc -- see note
c
      implicit none
      real * 8 s,x,yi,cth2,tk,uk,q1q,q2q,xinv(5),xmatout(7)
      integer jproc
      include 'stpcblks.h'
      real * 8 xii,xmatin(7)
      integer i,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
c Hard-coded choice for Ht mode; for s- and t-channel see freal()
      parameter (jch=3)
c
c      write(6,*) 'jproc= ',jproc
      xii=1-x
      if(jproc.eq.1)then
        call frealht_gg(s,xii,yi,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jch)
          xmatout(i)=xmatin(i)
        enddo
      elseif(jproc.eq.2)then
        call frealht_qq(s,xii,yi,cth2,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jch)
          xmatout(i)=xmatin(i)
        enddo
      elseif(jproc.eq.3)then
        call frealht_qg(s,xii,yi,cth2,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jch)
          xmatout(i)=xmatin(i)
        enddo
      else
        write(*,*)'Unknown process in frealht',jproc
        stop
      endif
c
      return
      end


      subroutine frealht_gg(xs,xxii,xyi,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for gg --> tHbbar. See frealht for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xtk,xuk,xq1q,xq2q,xxinv(5),xmatout(7)
      real * 8 s,xii,yi,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vtr,xnc,s_red,t_red,x_ap,ap_kern,xfact,xinv(5),xmatin(7)
      integer i,jch,ione,ithree,icode,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      parameter (tiny=1.d-6)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      parameter (jch=3)
      parameter (ione=1)
      parameter (ithree=3)
c
      s=xs
      xii=xxii
      yi=xyi
      tk=xtk
      uk=xuk
      q1q=xq1q
      q2q=xq2q
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      xnorm=vtr/(4*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(ione,jch)
        xmatout(i)=0.d0
      enddo
      if(xii.lt.tiny)then
c Soft limit
        continue
      elseif(abs(1-yi).lt.tiny)then
c Collinear + limit
        icode=2
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q2q
        call fbornht(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ione,jch)
          if(idrlimcp(jch,ione,i).ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jch,ione,i))
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        icode=2
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornht(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ione,jch)
          if(idrlimcm(jch,ione,i).ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jch,ione,i))
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call htmat_gg(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(ione,jch)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
c
      return
      end


      subroutine frealht_qq(xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for qq --> tHbbar. See frealht for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,xxinv(5),xmatout(7)
      real * 8 s,xii,yi,cth2,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vtr,xnc,s_red,t_red,x_ap,ap_kern,qin_kern,xfact,xinv(5),
     # xmatin(7),xaziin(7)
      integer i,jch,itwo,ithree,ipone,imone,icode,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      parameter (tiny=1.d-6)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      parameter (jch=3)
      parameter (itwo=2)
      parameter (ithree=3)
      parameter (ipone=1)
      parameter (imone=-1)
c
      s=xs
      xii=xxii
      yi=xyi
      cth2=xcth2
      tk=xtk
      uk=xuk
      q1q=xq1q
      q2q=xq2q
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      xnorm=vtr/(4*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(itwo,jch)
        xmatout(i)=0.d0
      enddo
      if(xii.lt.tiny)then
c Soft limit
        continue
      elseif(abs(1-yi).lt.tiny)then
c Collinear + limit
        icode=3
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q2q
        call fbornht(s_red,t_red,ithree,xmatin)
        call faziht(s_red,t_red,cth2,itwo,ipone,xaziin)
        do i=1,idrmax(itwo,jch)
          if(idrlimcp(jch,itwo,i).ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jch,itwo,i)) +
     #                 4*(1+yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        icode=3
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornht(s_red,t_red,ithree,xmatin)
        call faziht(s_red,t_red,cth2,itwo,imone,xaziin)
        do i=1,idrmax(itwo,jch)
          if(idrlimcm(jch,itwo,i).ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jch,itwo,i)) +
     #                 4*(1-yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call htmat_qq(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(itwo,jch)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
c
      return
      end


      subroutine frealht_qg(xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for qg --> tHbbar. See frealht for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,xxinv(5),
     # xmatout(7)
      real * 8 s,xii,yi,cth2,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vca,vtr,xnc,s_red,t_red,x_ap,ap_kern,qin_kern,xfact,xinv(5),
     # xmatin(7),softfc(3),xaziin(7)
      integer i,jch,ithree,ipone,imone,icode,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer iapcp(3,7),iapcm(3,7)
      common/ciap/iapcp,iapcm
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      parameter (tiny=1.d-6)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vca=3.d0)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      parameter (jch=3)
      parameter (ithree=3)
      parameter (ipone=1)
      parameter (imone=-1)
c
      s=xs
      xii=xxii
      yi=xyi
      cth2=xcth2
      tk=xtk
      uk=xuk
      q1q=xq1q
      q2q=xq2q
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      xnorm=vtr/(4*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(ithree,jch)
        xmatout(i)=0.d0
        softfc(i)=0d0
      enddo
      if(xii.lt.tiny)then
        s_red = s
        t_red = q1q
c Soft limit
c The following assumes that 
c   p_i.k = sqrt{s}*(1-x)*xinv(i)    i=1,2
c   k_i.k = sqrt{s}*(1-x)*xinv(i+3)  i=1,2
c and that xinv(1)=sqrt(s)*(1-yi)/4, xinv(2)=sqrt(s)*(1+yi)/4; 
c softfc here is the combination of the eikonals, times the damping 
c factor (1-yi**2)*xii**2
         softfc(1) = 16/s_red * (
     #    (-vcf)*xm12*(1-yi**2)/(16*xinv(4)**2)  + 
     #    (vcf-vca/2.)*(-q1q+xm12)*(1.+yi)/(4*xinv(4)*sqrt(s)) + 
     #    (vca/2.)*(-q2c+xm12)*(1.-yi)/(4*xinv(4)*sqrt(s)) + 
     #    (vca/2.)*1.)
         softfc(3) = 16/s_red * (
     #    (-vcf)*xm12*(1-yi**2)/(16*xinv(4)**2)  + 
     #    (vcf-vca/2.)*(-q2c+xm12)*(1.-yi)/(4*xinv(4)*sqrt(s)) + 
     #    (vca/2.)*(-q1q+xm12)*(1.+yi)/(4*xinv(4)*sqrt(s)) + 
     #    (vca/2.)*1.)
        call fbornht(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ithree,jch)
             xmatout(i)=softfc(i)*xmatin(i) 
        enddo
      elseif(abs(1-yi).lt.tiny)then
c Collinear + limit
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q2q
        call fbornht(s_red,t_red,ithree,xmatin)
        call faziht(s_red,t_red,cth2,ithree,ipone,xaziin)
        do i=1,idrmax(ithree,jch)
          icode=iapcp(jch,i)
          if(idrlimcp(jch,ithree,i).ne.0.and.icode.ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jch,ithree,i)) +
     #                 4*(1+yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornht(s_red,t_red,ithree,xmatin)
        call faziht(s_red,t_red,cth2,ithree,imone,xaziin)
        do i=1,idrmax(ithree,jch)
          icode=iapcm(jch,i)
          if(idrlimcm(jch,ithree,i).ne.0.and.icode.ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jch,ithree,i)) +
     #                 4*(1-yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call htmat_qg(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(ithree,jch)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
      return
      end


      subroutine f2svht(s,t,jproc,xmatout)
c Soft-virtual contribution. The normalization is such that
      implicit none
      include 'stpcblks.h'
      real*8 s,t,xmatout(7)
      real*8 pi,vca,vtr,u,beta0,scalefac,svfac,svnonfac,
     # xfac(3),xnonfac(3),xmatin(7),vcf
      real*8 xfac2(3),xnonfac2(3)
      real*8 svfac2,svnonfac_stef
      parameter (pi=3.14159265358979312D0)
      parameter (vca=3.d0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtr=0.5d0)
      integer idrmax(1:3,3)
      common/cidrmax/idrmax
      integer jproc,jch,ithree,i
      parameter (jch=3)
      parameter (ithree=3)
c
      if(jproc.lt.1.or.jproc.gt.3)then
        write(*,*)'Error in f2svht: jproc=',jproc
        stop
      endif
      do i=1,7
        xmatout(i)=0.d0
        if(i.le.3)then
          xfac(i)=0.d0
          xnonfac(i)=0.d0
        endif
      enddo
      if(jproc.eq.3)then
        u=xm12+xm22-s-t
        beta0=(11*vca-4*vtr*nl)/(12*pi)
        scalefac=(2*pi*beta0+2.*vcf*3./2.)*log(xmur2/xmuf2h1) ! Second term from Yukawa coupling
c Soft-virtual function have been derived for idr=3. See finiteWt.m
        xfac(1)=svfac(s,u,t,xm12,xm22,xmuf2h1)
        xnonfac(1)=svnonfac(s,u,t,xm12,xm22,xmuf2h1)
        xfac(3)=svfac(s,t,u,xm12,xm22,xmuf2h1)
        xnonfac(3)=svnonfac(s,t,u,xm12,xm22,xmuf2h1)

c        xfac2(1)=svfac2(s,u,t,xm12,xm22,xmuf2h1,5)
c        xnonfac2(1)=svnonfac_stef(s,u,t,xm12,xm22,xmuf2h1,5)
c        xfac2(3)=svfac2(s,t,u,xm12,xm22,xmuf2h1,5)
c        xnonfac2(3)=svnonfac_stef(s,t,u,xm12,xm22,xmuf2h1,5)
        
        call fbornht(s,t,ithree,xmatin)
        do i=1,idrmax(jproc,jch)
          xmatout(i)=(xfac(i)+scalefac)*xmatin(i)+xnonfac(i)/(2*s)
c          write(*,*) 'Chris: ',xfac(i)*xmatin(i)+xnonfac(i)/2/s
c          write(*,*) 'Stefano: ',xfac2(i)*xmatin(i)+xnonfac2(i)/2/s
        enddo
      endif
      return
      end


      subroutine f2prht(xs,xt,xx,xxc,xyic,xxlmude,jproc,xmatout)
c Collinear reminder contribution of FKS. The normalization is such that
c    dsigma_pr = g_S^2 alpha_S/(2*pi) g_W^2 |V_tj|^2 f2prht/xii dxii dphi_2
c with g_S^4, g_W^2 and CKM factors inserted in the main code. See frealht
c for the conventions used for the array xmatout. Entries:
c   xs, xt = Mandelstam invariants
c   xx = 1-xii, with xii that of the event
c   xxc = 1-xii, with xii that of the event or of the counterevent
c   xyic = 1  ==> selects the first 3 lines in 5.7  (collinear +)
c   xyic = -1 ==> selects the second 3 lines in 5.7  (collinear -)
c   if xxc # 1, takes the regular part of K 
c   if xxc = 1  takes the delta-function part of K 
      implicit none
      real * 8 xs,xt,xx,xxc,xyic,xxlmude,xmatout(7)
      include 'stpcblks.h'
      real * 8 one,s,t,xii,xiic,yic,xlmude,x_ap,s_red,xdfct1,ap_kern,
     # xdfct2,apprime_kern,xdfct3p,xdfct3l,xdfct5,xkplus,xklog,xkreg,
     # xkdelta,tmp,xmatin(7)
      character * 2 scheme
      real * 8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
      integer jch,ithree,i,jproc,icode,ilim
      integer idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      integer iapcp(3,7),iapcm(3,7)
      common/ciap/iapcp,iapcm
      parameter (one=1.d0)
      parameter (jch=3)
      parameter (ithree=3)
c
      s=xs
      t=xt
      xii=1-xx
      xiic=1-xxc
      yic=xyic
      xlmude=xxlmude
c
      x_ap=1-xiic
      s_red=s*x_ap
c The appropriate value for t is set in the main code,
c here it is the same as entered in the call to f2prht
      call fbornht(s_red,t,ithree,xmatin)
      if(yic.eq.1.d0)then
        scheme=schhad1
      elseif(yic.eq.-1.d0)then
        scheme=schhad2
      else
        write(6,*)'Error in f2prht',yic
        stop
      endif

c     Have checked the following Altarelli-Parisi identifications:

c
      do i=1,idrmax(jproc,jch)
        if(jproc.eq.1)then
          icode=2
        elseif(jproc.eq.2)then
          icode=3
        elseif(jproc.eq.3)then
          if(yic.eq.1.d0)then
            icode=iapcp(jch,i)
          elseif(yic.eq.-1.d0)then
            icode=iapcm(jch,i)
          endif
        else
          write(*,*)'Unknown process in f2prht',jproc
          stop
        endif

c     Have checked the following IDR identifications:
        if(yic.eq.1.d0)then
          ilim=idrlimcp(jch,jproc,i)
        elseif(yic.eq.-1.d0)then
          ilim=idrlimcm(jch,jproc,i)
        endif
c


        if(icode.ne.0.and.ilim.ne.0)then
          xdfct1=ap_kern(x_ap,abs(icode))
          xdfct2=apprime_kern(x_ap,abs(icode))
          xdfct3p=0.d0
          xdfct3l=0.d0
          xdfct5=0.d0
c
          if(scheme.eq.'DI')then
            xdfct3p=xkplus(x_ap,abs(icode))
            xdfct3l=xklog(x_ap,abs(icode))
            if(xiic.ne.0.d0)then
              xdfct5=xkreg(x_ap,abs(icode))
            else
              xdfct5=xkdelta(abs(icode))
     #              +xkplus(one,abs(icode))*log(xicut)
     #              +xklog(one,abs(icode))*log(xicut)**2/2.d0
c This part contributes to sig2pr(soft), which is integrated in xi
c over the range (0,xicut). This implies the presence of a jacobian
c equal to xicut in the soft term, which has to be removed by hand
c in this case
              xdfct5=xdfct5/xicut
            endif
          elseif(scheme.ne.'MS')then
            write(6,*)'Error in f2prht, y=',yic
            write(6,*)'Factorization scheme ',scheme,' not known'
          endif
c 1/xi is the main code

c     Have checked the following against FKS eq. (5.7):
          tmp=xdfct1*(xlmude+2*log(xii))-xdfct2
     #       -xdfct3p-xdfct3l*log(xii) 
     #       -xii*xdfct5
        else
          tmp=0.d0
        endif
c
        if(ilim.ne.0)then
          if(icode.eq.0)then
            write(6,*)'Inconsistency in f2prht'
            write(6,*)ilim,icode,i,jproc
            stop
          endif
          xmatout(i)=tmp*xmatin(ilim)
        else
          if(icode.ne.0.and.jproc.eq.3)then
            write(6,*)'Inconsistency in f2prht'
            write(6,*)ilim,icode,i,jproc
            stop
          endif
          xmatout(i)=0.d0
        endif
      enddo
      return
      end


      subroutine fbornht(s,t,jproc,xmatout)
c Returns the partonic Born contribution, with the following normalization
c    dsigma_born = g_S^2 (|a_tj|^2+|b_tj|^2) |V_tj|^2 fborn dphi_2
c with g_S^2, a_tj, b_tj and CKM factors inserted in the main code. See frealht
c for the conventions used for the array xmatout
c
c With p1 denoting the momentum of the parton coming from the left
c (= +z-component) and k1 the momentum of the top quark, and k2 that
c of the W, we define the invariants as follows
c   s=(p1+p2)^2=(k1+k2)^2
c   t=(p1-k1)^2=(p2-k2)^2
c   u=(p1-k2)^2=(p2-k1)^2
      implicit none
      include 'stpcblks.h'
      real*8 s,t,u,vtr,xnc,xnorm,elborn,xmatout(7)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      integer jproc,i,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (jch=3)
c
      if(jproc.lt.1.or.jproc.gt.3)then
        write(*,*)'Error in fbornht: jproc=',jproc
        stop
      endif
c
      do i=1,idrmax(jproc,jch)
        xmatout(i) = 0d0
      enddo
c
      if(jproc.eq.3)then
        xnorm=vtr/(4*xnc)  ! Factor 4 from spin averaging. 1 / (2 N_C) from colour factor and colour average.
        u=xm12+xm22-s-t
        xmatout(1)=xnorm*elborn(s,u)/(2*s)
        xmatout(3)=xnorm*elborn(s,t)/(2*s)
      endif
c
      return
      end



      subroutine htmat_gg(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,
     # t23,tpq,tp3,tq3,H_gg
      integer i,ione,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (ione=1)
      parameter (jch=3)
c
      do i=1,idrmax(ione,jch)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c Write a proper comment here for the mapping of invariants
c
      t12 = s
      t1p = q1q-xm12
      t1q = -q1c+xm22
      t13 = tk
      t2p = q2c-xm12
      t2q = -q2q+xm22
      t23 = uk
      tpq = -s2+xm12+xm22
      tp3 = w1-xm12
      tq3 = -w2+xm22
c
      xmatout(1)=H_gg(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
c the two gluon are exchanged in the call below
c      xmatout(1)=H_gg(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
c
      return
      end


      subroutine htmat_qq(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,
     # t23,tpq,tp3,tq3,Nc,HQqqq_1,HQqqq_2
      real*8 HQqqq_1_nores,HQqqq_2_nores
      integer i,itwo,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (itwo=2)
      parameter (jch=3)
      parameter(Nc = 3d0)

      do i=1,idrmax(itwo,jch)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c Write a proper comment here for the mapping of invariants
c
      t12 = s
      t1p = q1q-xm12
      t1q = -q1c+xm22
      t13 = tk
      t2p = q2c-xm12
      t2q = -q2q+xm22
      t23 = uk
      tpq = -s2+xm12+xm22
      tp3 = w1-xm12
      tq3 = -w2+xm22

c


      if(xm22.gt.xm12) then
      xmatout(1)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .      HQqqq_1(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
      xmatout(2)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .      HQqqq_1(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
      xmatout(3)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .      HQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
      xmatout(4)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .      HQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
      xmatout(5)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .      HQqqq_1(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)+
     .      HQqqq_1(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
     .  +    0.5*(Nc/4.-1./(4.*Nc))*(
     .      HQqqq_2(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)) )
      xmatout(6)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .      HQqqq_1(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)+
     .      HQqqq_1(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
     .  +    0.5*(Nc/4.-1./(4.*Nc))*(
     .      HQqqq_2(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)) )
      xmatout(7)= (Nc**2-1.)/(2.*Nc)*(-1)*(
     .      HQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)+
     .      HQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
     .  +    0.5*(Nc/4.-1./(4.*Nc))*(
     .      HQqqq_2(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)) )
      else
      xmatout(1)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  HQqqq_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
      xmatout(2)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  HQqqq_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
      xmatout(3)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  HQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
      xmatout(4)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  HQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
      xmatout(5)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .  HQqqq_1(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)+
     .  HQqqq_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
     . +0.5*(Nc/4.-1./(4.*Nc))*(
     .  HQqqq_2_nores(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)) )
      xmatout(6)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .  HQqqq_1(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)+
     .  HQqqq_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
     . +0.5*(Nc/4.-1./(4.*Nc))*(
     .  HQqqq_2_nores(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)) )
      xmatout(7)= (Nc**2-1.)/(2.*Nc)*(-1)*(
     .  HQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)+
     .  HQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
     .  +  0.5*(Nc/4.-1./(4.*Nc))*(
     .  HQqqq_2(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)) )         
      endif

c
      return
      end


      subroutine htmat_qg(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,
     # t23,tpq,tp3,tq3,H_qg
      integer i,ithree,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (ithree=3)
      parameter (jch=3)
c
      do i=1,idrmax(ithree,jch)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c Write a proper comment here for the mapping of invariants
c
      t12 = s
      t1p = q1q-xm12
      t1q = -q1c+xm22
      t13 = tk
      t2p = q2c-xm12
      t2q = -q2q+xm22
      t23 = uk
      tpq = -s2+xm12+xm22
      tp3 = w1-xm12
      tq3 = -w2+xm22
c
      xmatout(1)=H_qg(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
      xmatout(3)=H_qg(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)


      return
      end


      subroutine faziht(s,t,cth2,jproc,icoll,xmatout)
c Returns the azimuthal-dependent part of collinear limits; the definitions
c and normalization are as in fbornht
      implicit none
      include 'stpcblks.h'
      real*8 s,t,cth2,xmatout(7)
      real*8 vtr,xnc,xnorm,u,azidep,azi1,azi3
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      integer jproc,icoll,i,jch,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (jch=3)
c
      if(jproc.lt.2.or.jproc.gt.3)then
        write(*,*)'Error in faziht: jproc=',jproc
        stop
      endif
c
      do i=1,idrmax(jproc,jch)
        xmatout(i) = 0d0
      enddo
c
      xnorm=vtr/(4*xnc)
      u=xm12+xm22-s-t
      azi1=xnorm*azidep(s,t,u,xm12,xm22,cth2)/(2*s)
      azi3=xnorm*azidep(s,u,t,xm12,xm22,cth2)/(2*s)
      if(icoll.eq.1)then
        if(jproc.eq.2)then
          xmatout(4)=azi3
          xmatout(6)=azi3
          xmatout(7)=azi3
        elseif(jproc.eq.3)then
          xmatout(3)=azi3
        endif
      elseif(icoll.eq.-1)then
        if(jproc.eq.2)then
          xmatout(3)=azi1
          xmatout(5)=azi1
          xmatout(7)=azi1
        elseif(jproc.eq.3)then
          xmatout(1)=azi1
        endif
      else
        write(*,*)'Error in faziht: icoll=',icoll
        stop
      endif
c
      return
      end

c
c
c End of wrapper routines
c
c
c
c
c Begin of actual matrix element routines
c
c


      function elborn(s,t)
c Original Born function; checked against two independent computations
      implicit none
      real*8 elborn,s,t
      include 'stpcblks.h'
      real*8 m2,t1,q2,MbMb,mH2,MbMb2
c
      m2 = xm12
      t1 = t-m2
      mH2 = xm22
c
      MbMb =
     &  - 8 - 8*s**(-1)*t1**(-1)*mH2**2 + 16*s**(-1)*t1**(-1)*m2*mH2 -
     & 8*s**(-1)*t1**(-1)*m2**2 + 8*s**(-1)*mH2 - 8*s**(-1)*m2 - 4*
     & s**(-1)*t1 + 8*t1**(-2)*m2*mH2 - 8*t1**(-2)*m2**2 + 8*t1**(-1)*
     & mH2 - 8*t1**(-1)*m2 - 4*s*t1**(-1)

      elborn = MbMb
      return
      end


      function azidep(s,t,u,m2,q2,cth2)
c Returns the azimuthal-dependent part for the y=-1 collinear singularity
c of any real process that factorizes
c  b(p1)+g(p2) --> t(k1)+H(k2) 
c The definitions of the invariants and the normalizations are as in bgborn()
      implicit none
      real * 8 azidep,s,t,u,m2,q2,cth2
c

      azidep =
     & 4*q2*s**(-1) - 4*m2*s**(-1) + 4/(u - m2)*q2 - 4/(u - m2)*q2**2*
     & s**(-1) - 4/(u - m2)*m2 + 8/(u - m2)*m2*q2*s**(-1) - 4/(u - m2)*
     & m2**2*s**(-1) - 8/(u - m2)/(u - m2)*q2*s**(-1)*u**2*cth2**2 - 8
     & /(u - m2)/(u - m2)*q2*u*cth2**2 + 8/(u - m2)/(u - m2)*q2**2*
     & s**(-1)*u*cth2**2 + 8/(u - m2)/(u - m2)*m2*s**(-1)*u**2*cth2**2
     &  + 8/(u - m2)/(u - m2)*m2*u*cth2**2 + 4/(u - m2)/(u - m2)*m2*q2
     &  - 8/(u - m2)/(u - m2)*m2*q2**2*s**(-1)*cth2**2 - 8/(u - m2)/(u
     &  - m2)*m2**2*s**(-1)*u*cth2**2 - 4/(u - m2)/(u - m2)*m2**2 + 8/(
     & u - m2)/(u - m2)*m2**2*q2*s**(-1)*cth2**2


c      azidep = -8*(2*cth2**2-1)*(mw2-m2)*(2*mw2+m2)*(m2*mw2-t*u)/(mw2*s*
c     1   (u-m2)**2)
      return 
      end 



*---------------------------------------------------------------
      real*8 function H_gg(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2, Nc
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)
      parameter(Nc = 3d0)

      real*8 HQqgg_1, HQqgg_2,HQqgg_1_nores,HQqgg_2_nores
      real*8 m2, q2, s12, s13, s23, sQ1, sQ2, sQ3
      
      H_gg=0.


      if(q2.gt.m2) then
      H_gg = 0.5*
     .  (HQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  HQqgg_1(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2) +
     .  (HQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  HQqgg_1(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2))/(Nc**2-1d0)+
     .  HQqgg_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)*(-1d0/(Nc**2-1d0)))
      else
      H_gg = 0.5*
     .  (HQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  HQqgg_1_nores(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2) +
     .  (HQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  HQqgg_1_nores(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2))/(Nc**2-1d0)+
     .  HQqgg_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)*(-1d0/(Nc**2-1d0)))
      endif
      

c      write(*,*) 'ggm ggnores= ',HQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3),
c     &HQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)

      return
      end



*---------------------------------------------------------------
      real*8 function H_qg(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2, Nc
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)
      parameter(Nc = 3d0)

      real*8 HQqgg_1, HQqgg_2,HQqgg_1_nores,HQqgg_2_nores
      real*8 m2, q2, s12, s13, s23, sQ1, sQ2, sQ3

      H_qg=0.

C Overall minus sign due to fermion crossing!
      H_qg = -0.5*Nc*(
     .  HQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  HQqgg_1(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2) ) 
     . + HQqgg_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)*( 1d0/(2*Nc) )


      return
      end


*---------------------------------------------------------------
      real*8 function HQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqgg1
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2
      real*8 lam,Dlam
      integer gaugeflag

      integer i


      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3


      HQqgg1 =
     &  + D2**(-2)*s12**(-1)*s23 * ( 8*q2 - 8*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2) * ( 24*q2 - 24*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2)*s12*s23**(-1) * ( 32*q2 - 32*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2)*s12**2*s23**(-2) * ( 16*q2 - 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*sQ3**(-1) * (  - 8*q2**2 +
     &    16*m**2*q2 - 8*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1) * (  - 8*q2 + 8*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*sQ2*sQ3**(-1) * ( 8*q2 - 8*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*s23*sQ3**(-1) * ( 8*q2 - 24*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ3**(-1) * (  - 16*q2**2 +
     &    32*m**2*q2 - 16*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1) * (  - 16*q2 + 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ3 * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*sQ3**(-1) * ( 16*q2 - 48*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*sQ2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-2) * (  - 16*q2 + 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-2)*sQ3 * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-2)*sQ2*sQ3**(-1) * (  - 16*
     &    q2 + 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-2)*sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1)*sQ3**(-1) * ( 8*q2 - 24*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1) * (  - 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1)*sQ2*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**2*s23**(-2) * (  - 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**2*s23**(-2)*sQ2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-2)*s13*sQ2*sQ3**(-1) * ( 8*q2
     &     - 8*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-2)*s13*sQ2 * (  - 4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-2)*s13*sQ2**2*sQ3**(-1) * (
     &     - 4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-2)*s13**2*sQ2*sQ3**(-1) * (
     &     - 4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ3**(-1) * ( 24*q2**2 - 48
     &    *m**2*q2 + 24*m**4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1) * (  - 20*q2 + 20*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ3 * ( 12 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2*sQ3**(-1) * (  - 12*q2
     &     + 12*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2 * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ3**(-2) * (  - 32*
     &    m**2*q2 + 32*m**4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ3**(-1) * (  - 20*q2
     &     + 20*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13 * ( 24 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-2) * ( 16*
     &    m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-1) * ( 4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13**2*sQ3**(-1) * ( 12 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ3**(-2) * (  - 16*m**2*q2 + 16*m**4
     &     )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ3**(-1) * (  - 28*q2 + 44*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1) * ( 24 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*s13*sQ3**(-2) * ( 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s13*sQ3**(-1) * ( 24 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23*sQ3**(-1) * ( 12 )
      HQqgg1 = HQqgg1 + s23**(-2) * (  - 4*q2 + 4*m**2 )
      HQqgg1 = HQqgg1 + s23**(-2)*sQ2*sQ3**(-1) * ( 16*q2 - 16*m**2 )
      HQqgg1 = HQqgg1 + s23**(-2)*sQ2 * (  - 8 )
      HQqgg1 = HQqgg1 + s23**(-2)*sQ2**2*sQ3**(-1) * (  - 8 )
      HQqgg1 = HQqgg1 + s23**(-2)*s13*sQ2*sQ3**(-1) * (  - 16 )
      HQqgg1 = HQqgg1 + s23**(-1)*sQ3**(-1) * (  - 24*q2 + 48*m**2 )
      HQqgg1 = HQqgg1 + s23**(-1) * ( 20 )
      HQqgg1 = HQqgg1 + s23**(-1)*sQ2*sQ3**(-1) * (  - 12 )
      HQqgg1 = HQqgg1 + s23**(-1)*s13*sQ3**(-2) * ( 32*m**2 )
      HQqgg1 = HQqgg1 + s23**(-1)*s13*sQ3**(-1) * ( 12 )
      HQqgg1 = HQqgg1 + sQ3**(-2) * ( 16*m**2 )
      HQqgg1 = HQqgg1 + sQ3**(-1) * ( 12 )
      HQqgg1 = HQqgg1 + s12*s23**(-2) * ( 12 )
      HQqgg1 = HQqgg1 + s12*s23**(-2)*sQ2*sQ3**(-1) * (  - 8 )
      HQqgg1 = HQqgg1 + s12*s23**(-1)*sQ3**(-1) * ( 8 )


      HQqgg_1_nores = HQqgg1


      return
      end


*---------------------------------------------------------------
      real*8 function HQqgg_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
c      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqgg2
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2


      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      HQqgg2 =
     &  + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*sQ2*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     & sQ2**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2**2 * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1*sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     & sQ2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2**2 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ2 * (  -
     &    24 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ1 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ2 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ1 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ3**2 * (
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ2**2 * (
     &     - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ1*sQ3 * (
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ1*sQ2 * (
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ3**2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ2*sQ3 * ( 40
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ1*sQ3 * ( 80
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2 * ( 80
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ2 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ1 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*sQ2*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*sQ2**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*sQ1*sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*sQ1*sQ2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s13*sQ3 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s13*sQ1 * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23*sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23*sQ2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s23*sQ1 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ3**2 * (
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ2*sQ3
     &  * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ2**2 * (
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ1*sQ3
     &  * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ1*sQ2
     &  * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ2*sQ3 * ( 40
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ2**2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ1*sQ3 * ( 80
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2 * ( 80
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ3 * ( 16
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ2 * ( 16
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ1 * ( 64
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ3**2 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*sQ2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*sQ1 * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ1 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ3**2 * (
     &     - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ2**2 * (
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ1*sQ3 * (
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ1*sQ2 * (
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ3**2 * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1*sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ1 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ3 * (  -
     &    24 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ1 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     & sQ3**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     & sQ2*sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     & sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     & sQ2*sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     & sQ2**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     & sQ2**3*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     & sQ1*sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     & sQ1*sQ2**2*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     & sQ2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     & sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     & sQ3**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     & sQ2**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     & sQ2**3*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     & sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     & sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     & sQ2**2*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2 * ( 72 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2**2*sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ2**2*sQ3**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1*sQ2*sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     & sQ1*sQ2*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     & sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     & sQ2*sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     & sQ2*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**3*sQ3**(-1)
     &  * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**(-1)*
     & sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**2*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1)*
     & sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ2**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ3**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13 * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1)
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)*
     & sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     & sQ3**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     & sQ3**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ3**(-1)
     &  * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     & sQ3**(-2) * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-2)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2*sQ2**(-1)*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2*sQ2*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**3*sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**3*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23 * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2*sQ3**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)*
     & sQ3 * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ3**(-1)
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1 * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ2**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**3*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ3**(-1)
     &  * (  - 56*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ2*
     & sQ3**(-2) * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ2*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ2**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ3**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     & sQ3**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ3**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1 * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**2*sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**2*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**3*sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ3**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ2**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*
     & sQ3**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ2**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1**2*sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ3**(-2)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**4*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ3**2 * (  -
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ2*sQ3 * (  -
     &    40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ2**3*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ3 * (  -
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ2 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ2**2*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ2*sQ3 * (  - 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ2**2 * (  - 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ1*sQ3 * (  - 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ1*sQ2 * (  - 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**(-1)*sQ3
     &  * ( 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**(-1)*
     & sQ3**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13 * ( 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ3 * (  - 72 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2 * (  - 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**2*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ2**(-1)*sQ3**3
     &  * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**2*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*
     & sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ2**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ3**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-2)*sQ3**2 * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3 * ( 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1) * ( 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ3 * (  - 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2*sQ3**(-1) * ( 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2 * (  - 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-1) * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-2)*sQ3 * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1) * (  - 48*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ3**(-1) * (  - 48*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1) * ( 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1)*sQ3 * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ3**(-1) * ( 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13 * (  - 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-1) * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-2)*
     & sQ3**2 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*sQ3
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*
     & sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-2)*
     & sQ3 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ3**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1 * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**3*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-2)*sQ3 * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1) * ( 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1)*sQ3 * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ3**(-1) * ( 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23 * (  - 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-1) * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1**2*sQ2**(-1)*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-1)*sQ3**(-1)
     &  * ( 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-1)*sQ3**(-1)
     &  * ( 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ2**(-2)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2**(-1)*
     & sQ3**3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2*sQ3 * (  -
     &    40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2**2 * (  -
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ2**(-1)*
     & sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ2 * (  -
     &    32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**(-1)*
     & sQ3**2 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**2*
     & sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     & sQ2**(-1)*sQ3**3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     & sQ3**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     & sQ2**2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     & sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     & sQ2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**(-2)*
     & sQ3**2 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1) * ( 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * (  - 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2*sQ3**(-1)
     &  * ( 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2 * (  - 72 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**2*
     & sQ3**(-1) * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-1)*
     & sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-1)*
     & sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2*
     & sQ3**(-1) * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-2)*
     & sQ3**2 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*
     & sQ3**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-2)*
     & sQ3 * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-2)*sQ3 * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1) * ( 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1)*sQ3 * (  - 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ3**(-1) * ( 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12 * (  - 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2*sQ3**(-1) * (  - 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-1)
     &  * (  - 56*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-1)*
     & sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1)*sQ3**(-1)
     &  * ( 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*
     & sQ2**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*
     & sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     & sQ2**(-1)*sQ3**3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     & sQ3**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     & sQ2*sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     & sQ1*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2**(-2)*sQ3**2 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2**(-1)*sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2**(-1)*sQ3**2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ3 * ( 72 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1*sQ2**(-2)*sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1*sQ2**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     & sQ1 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-1)*
     & sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)
     &  * (  - 80*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)*
     & sQ3 * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     & sQ2**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     & sQ2**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     & sQ2**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     & sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     & sQ2**(-2)*sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     & sQ2**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     & sQ2**(-1)*sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)
     &  * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-2)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-1)
     &  * ( 16 )


      HQqgg_2_nores =  HQqgg2
      
      return
      end





*---------------------------------------------------------------
      real*8 function HQqqq_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqqq1
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C HQqqq1 corresponds to the g graph squared (two fermion loops):
C HQqqq2 corresponds to the interference graph
C Resonant diagrams are removed


      HQqqq1 =
     +  + D2**(-2)*s12*s23**(-1) * (  - 16*m**2 + 16*q2 )
      HQqqq1 = HQqqq1 + D2**(-2)*s12**2*s23**(-2) * (  - 16*m**2 + 16*
     +    q2 )
      HQqqq1 = HQqqq1 + D2**(-2) * (  - 8*m**2 + 8*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-2)*sQ2 * (  - 8 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-2)*sQ3 * ( 8 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-2) * ( 16*m**2 - 16*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-1) * (  - 16 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12**2*s23**(-2) * (  - 16 )
      HQqqq1 = HQqqq1 + D2**(-1)*s23**(-1)*sQ2 * (  - 8 )
      HQqqq1 = HQqqq1 + D2**(-1)*s23**(-1) * ( 8*m**2 - 8*q2 )
      HQqqq1 = HQqqq1 + D2**(-1) * (  - 8 )
      HQqqq1 = HQqqq1 + s12*s23**(-2) * ( 16 )
      HQqqq1 = HQqqq1 + s23**(-2)*sQ2 * ( 8 )
      HQqqq1 = HQqqq1 + s23**(-1) * ( 8 )


      HQqqq_1_nores = HQqqq1

      return
      end


*---------------------------------------------------------------
      real*8 function HQqqq_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqqq2
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C HQqqq1 corresponds to the g graph squared (two fermion loops):
C HQqqq2 corresponds to the interference graph
C Resonant diagrams are removed


      HQqqq2 =
     +  + D2**(-2)*s12**2*s23**(-1)*s13**(-1) * ( 16*m**2 - 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * ( 16*m**2 - 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s13**(-1) * (  - 16*m2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)
     +  * ( 16*m**2 - 16*m2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3 * ( 16*m**2 - 
     +    16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1) * (  - 16*m**2*q2
     +     + 8*m**4 + 8*q2**2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1) * (  - 24*m**2
     +     + 24*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s13**(-1) * (  - 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*s23*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s23**(-1)*s13**(-1) * (  - 16*m**2
     +     + 16*m2 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3 * (  - 8*m**2
     +     + 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s13**(-1) * ( 8*m**2 - 8*q2 )
      HQqqq2 = HQqqq2 + s23**(-1)*s13**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + s23**(-1)*s13**(-1) * ( 8*m**2 - 8*q2 )
      HQqqq2 = HQqqq2 + s23**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + s13**(-1) * ( 16 )



      HQqqq_2_nores = HQqqq2

      return
      end







*---------------------------------------------------------------
      real*8 function HQqqq_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqqq1
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

c      write (*,*) 'q2=',q2
c      write (*,*) 'm2=',m2
c      write (*,*) 's12= ',s12
c      write (*,*) 's13= ',s13
c      write (*,*) 's23= ',s23
c      write (*,*) 'sQ1= ',sQ1
c      write (*,*) 'sQ2= ',sQ2
c      write (*,*) 'sQ3= ',sQ3

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C HQqqq1 corresponds to the f+g graphs squared (two fermion loops):
C They contribute at 1st and 3rd color order
C HQqqq2 corresponds to the interference graphs. They contribute 
C at 2nd and 4th color order.
 
      HQqqq1 =
     +  + D2**(-2)*s12*s23**(-1) * (  - 16*m**2 + 16*q2 )
      HQqqq1 = HQqqq1 + D2**(-2)*s12**2*s23**(-2) * (  - 16*m**2 + 16*
     +    q2 )
      HQqqq1 = HQqqq1 + D2**(-2) * (  - 8*m**2 + 8*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ3 * (  - 32*
     +    m**2 + 32*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 32 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-1) * (  - 48*m**2
     +     + 32*m2 + 16*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s12**2*s23**(-2) * (  - 32*
     +    m**2 + 32*m2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * (  - 16 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3 * (  - 16*m**2
     +     + 16*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 16 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1)*s23**(-1) * ( 32*m**2*q2 - 16
     +    *m**4 - 16*q2**2 )
      HQqqq1 = HQqqq1 + D2**(-1)*D4**(-1) * (  - 16*m**2 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-2) * ( 32*m**2 - 32*q2 )
      HQqqq1 = HQqqq1 + D2**(-1)*s12*s23**(-1) * (  - 16 )
      HQqqq1 = HQqqq1 + D2**(-1)*s23**(-1)*sQ2 * (  - 8 )
      HQqqq1 = HQqqq1 + D2**(-1)*s23**(-1)*sQ3 * ( 8 )
      HQqqq1 = HQqqq1 + D2**(-1) * (  - 8 )
      HQqqq1 = HQqqq1 + D4**(-2)*s12*s23**(-2)*sQ3 * (  - 8*m**2 + 8*m2
     +     )
      HQqqq1 = HQqqq1 + D4**(-2)*s23**(-2)*s13*sQ3 * ( 8*m**2 - 8*m2 )
      HQqqq1 = HQqqq1 + D4**(-2)*s23**(-2)*sQ3**2 * (  - 16*m**2 + 16*
     +    q2 )
      HQqqq1 = HQqqq1 + D4**(-2)*s23**(-1)*s13 * ( 8*m**2 - 8*m2 )
      HQqqq1 = HQqqq1 + D4**(-2)*s23**(-1)*sQ3 * (  - 16*m**2 + 16*q2 )
      HQqqq1 = HQqqq1 + D4**(-2)*s23**(-1) * ( 16*m**2*q2 - 16*m**4 )
      HQqqq1 = HQqqq1 + D4**(-2) * (  - 8*m**2 + 8*q2 )
      HQqqq1 = HQqqq1 + D4**(-1)*s12*s23**(-2) * ( 32*m**2 - 32*m2 )
      HQqqq1 = HQqqq1 + D4**(-1)*s12*s23**(-1) * ( 8 )
      HQqqq1 = HQqqq1 + D4**(-1)*s23**(-2)*s13 * (  - 8*m**2 + 8*m2 )
      HQqqq1 = HQqqq1 + D4**(-1)*s23**(-2)*sQ3 * ( 32*m**2 - 32*q2 )
      HQqqq1 = HQqqq1 + D4**(-1)*s23**(-1)*s13 * (  - 8 )
      HQqqq1 = HQqqq1 + D4**(-1)*s23**(-1)*sQ3 * (  - 16 )
      HQqqq1 = HQqqq1 + D4**(-1) * (  - 8 )
      HQqqq1 = HQqqq1 + s23**(-2) * (  - 16*m**2 + 16*q2 )




      HQqqq_1 = HQqqq1

      return
      end


*---------------------------------------------------------------
      real*8 function HQqqq_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqqq2
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C HQqqq1 corresponds to the f+g graphs squared (two fermion loops):
C They contribute at 1st and 3rd color order
C HQqqq2 corresponds to the interference graphs. They contribute 
C at 2nd and 4th color order.


      HQqqq2 =
     +  + D2**(-2)*s12**2*s23**(-1)*s13**(-1) * ( 16*m**2 - 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * ( 16*m**2 - 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12*s13**(-1) * (  - 16*m2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)
     +  * ( 16*m**2 - 16*m2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3 * ( 16*m**2 - 
     +    16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*s13**(-1) * (  - 16*m**2*q2
     +     + 8*m**4 + 8*q2**2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2*
     + sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2
     +  * ( 8*m**2 - 8*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * (  - 8*m**2 + 8*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*sQ2 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s23*s13**(-1) * ( 32 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1)*sQ2 * ( 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1) * ( 8*m**2 + 
     +    16*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12 * ( 24 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2 * ( 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)
     +  * ( 32*m**2 - 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1) * ( 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s13**(-1) * ( 40 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s12**3*s23**(-1)*s13**(-1)
     +  * ( 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s23*s13**(-1)*sQ2 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s23 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s23**2*s13**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ3 * (  - 16*m**2
     +     + 16*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1) * ( 16*m**2*q2 - 8
     +    *m**4 - 8*q2**2 )
      HQqqq2 = HQqqq2 + D2**(-1)*D4p**(-1)*sQ2 * ( 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1) * (  - 24*m**2
     +     + 24*q2 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s23**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1)*s12*s13**(-1) * (  - 16 )
      HQqqq2 = HQqqq2 + D2**(-1)*s23*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D2**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * ( 8*m**2 + 8*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1) * (  - 16*m**2
     +     + 16*m2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s23*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s13**(-1)*sQ3 * (  - 16
     +     )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12*s13**(-1) * (  - 8*m**2
     +     + 8*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s12 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ3**2
     +  * (  - 16*m**2 + 16*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*sQ3 * (  - 8*m**2
     +     - 8*m2 + 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s23*s13**(-1) * (  - 8*m**2
     +     + 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1)*sQ3 * (  - 24*m**2
     +     + 24*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1) * (  - 8*m**2*m2
     +     - 16*m**2*q2 + 8*m**4 + 8*m2*q2 + 8*q2**2 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*D4p**(-1) * ( 8*m**2 - 16*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s23**(-1)*s13**(-1) * (  - 16*m**2
     +     + 16*m2 )
      HQqqq2 = HQqqq2 + D4**(-1)*s12*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3 * (  - 8*m**2
     +     + 8*q2 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3**2 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s23**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4**(-1)*s13**(-1) * ( 8*m**2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2 * (  - 16
     +     )
      HQqqq2 = HQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 8
     +     )
      HQqqq2 = HQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1) * (  - 32*
     +    m**2 + 16*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s12*s13**(-1) * (  - 24 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s12**2*s23**(-1)*s13**(-1) * (  - 16
     +     )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ2*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ2 * (  - 8*m**2
     +     + 8*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ3 * ( 16*m**2
     +     - 16*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*s13 * ( 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*sQ2 * ( 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23**(-1) * ( 8*m**2 + 8*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s23*s13**(-1) * (  - 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s13**(-1)*sQ3 * (  - 8 )
      HQqqq2 = HQqqq2 + D4p**(-1)*s13**(-1) * ( 8*m**2 - 8*m2 - 8*q2 )
      HQqqq2 = HQqqq2 + D4p**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + s23**(-1)*s13**(-1)*sQ3 * ( 8 )
      HQqqq2 = HQqqq2 + s23**(-1)*s13**(-1) * ( 8*m**2 - 8*q2 )
      HQqqq2 = HQqqq2 + s23**(-1) * ( 8 )
      HQqqq2 = HQqqq2 + s13**(-1) * ( 16 )


      HQqqq_2 = HQqqq2

      return
      end





*---------------------------------------------------------------
      real*8 function HQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqgg1
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      HQqgg1 =
     &  + D2**(-2)*s12**(-1)*s23 * ( 8*q2 - 8*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2) * ( 24*q2 - 24*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2)*s12*s23**(-1) * ( 32*q2 - 32*m**2 )
      HQqgg1 = HQqgg1 + D2**(-2)*s12**2*s23**(-2) * ( 16*q2 - 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1)*sQ3**(-1) * ( 16*
     &    q2**3 - 48*m**2*q2**2 + 48*m**4*q2 - 16*m**6 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1) * (  - 32*q2**2 +
     &    64*m**2*q2 - 32*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1)*sQ3 * ( 24*q2 - 24*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1)*sQ3**2 * (  - 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ3**(-1) * ( 8
     &    *m**2*q2 - 8*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**(-1)*s23 * (  - 8*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s23**(-1) * (  - 64*q2**2 +
     &    128*m**2*q2 - 64*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3 * ( 32*q2 - 32*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 16 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*sQ3**(-1) * (  - 32*q2**2 +
     &    32*m**2*q2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1) * (  - 16*q2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*sQ3 * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s23*sQ3**(-1) * (  - 8*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ3 * ( 32*q2
     &     - 32*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3**(-1) * (
     &     - 32*m**2*q2 + 32*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12*s23**(-1) * ( 32*q2 - 32*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 32 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12*sQ3**(-1) * ( 24*q2 - 24*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12 * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * (  - 16 )
      HQqgg1 = HQqgg1 + D2**(-1)*D4**(-1)*s12**2*sQ3**(-1) * (  - 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*sQ3**(-1) * (  - 16*q2**2 +
     &    32*m**2*q2 - 16*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*sQ2*sQ3**(-1) * ( 8*q2 - 8*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**(-1)*s23*sQ3**(-1) * ( 8*q2 - 24*
     &    m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ3**(-1) * ( 16*q2**2 - 32*
     &    m**2*q2 + 16*m**4 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ3 * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ2*sQ3**(-1) * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ2 * (  - 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*sQ3**(-1) * ( 16*q2 - 24*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*sQ2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-2) * (  - 32*q2 + 32*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1)*sQ3**(-1) * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1) * (  - 16 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12*s23**(-1)*sQ2*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + D2**(-1)*s12**2*s23**(-1)*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D4**(-2)*s23**(-2)*sQ3**2 * ( 16*q2 - 16*m**2 )
      HQqgg1 = HQqgg1 + D4**(-2)*s23**(-1) * ( 32*m**2*q2 - 32*m**4 )
      HQqgg1 = HQqgg1 + D4**(-2)*s23**(-1)*sQ3 * ( 32*q2 - 32*m**2 )
      HQqgg1 = HQqgg1 + D4**(-2)*sQ3**(-2) * ( 32*m**4*q2 - 32*m**6 )
      HQqgg1 = HQqgg1 + D4**(-2)*sQ3**(-1) * ( 32*m**2*q2 - 32*m**4 )
      HQqgg1 = HQqgg1 + D4**(-2) * ( 24*q2 - 24*m**2 )
      HQqgg1 = HQqgg1 + D4**(-2)*s23*sQ3**(-1) * ( 8*q2 - 8*m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1) * ( 16*q2**2 - 32*
     &    m**2*q2 + 16*m**4 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*sQ3 * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*sQ3**2 * ( 8 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ3**(-1) * (
     &     - 32*m**2*q2 + 32*m**4 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*s13 * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ3 * ( 16 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23**(-1)*s13**2 * ( 8 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*sQ3**(-2) * (  - 32*m**2*
     &    q2**2 + 64*m**4*q2 - 32*m**6 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*sQ3**(-1) * (  - 16*q2**2 +
     &    32*m**2*q2 - 16*m**4 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1) * ( 16*q2 - 8*m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s13*sQ3**(-1) * ( 8*q2 - 8*
     &    m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s13 * ( 8 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12**(-1)*s23*sQ3**(-1) * ( 8*q2 - 8*
     &    m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s23**(-2)*sQ3 * (  - 32*q2 + 32*m**2 )
      HQqgg1 = HQqgg1 + D4**(-1)*s23**(-1)*sQ3 * (  - 16 )
      HQqgg1 = HQqgg1 + D4**(-1)*s23**(-1)*s13 * (  - 8 )
      HQqgg1 = HQqgg1 + D4**(-1)*sQ3**(-2) * ( 32*m**2*q2 - 32*m**4 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12*s23**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + D4**(-1)*s12*sQ3**(-2) * (  - 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ3**(-1) * ( 16*q2**2 - 32
     &    *m**2*q2 + 16*m**4 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1) * (  - 16*q2 + 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ3 * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2*sQ3**(-1) * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2 * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ3**(-1) * (  - 16*q2
     &     + 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13 * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23**(-1)*s13**2*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ3**(-2) * ( 32*m**2*q2 - 32*m**4 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ3**(-1) * (  - 16*q2 + 8*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ2*sQ3**(-2) * (  - 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*sQ2*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*s13*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23*sQ3**(-2) * (  - 16*m**2 )
      HQqgg1 = HQqgg1 + s12**(-1)*s23*sQ3**(-1) * ( 8 )
      HQqgg1 = HQqgg1 + s23**(-2) * ( 16*q2 - 16*m**2 )
      HQqgg1 = HQqgg1 + s23**(-1)*sQ3**(-1) * (  - 16*q2 + 16*m**2 )
      HQqgg1 = HQqgg1 + s23**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s23**(-1)*sQ2*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s23**(-1)*s13*sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + sQ3**(-2) * (  - 32*m**2 )
      HQqgg1 = HQqgg1 + sQ3**(-1) * ( 16 )
      HQqgg1 = HQqgg1 + s12*s23**(-1)*sQ3**(-1) * ( 8 )

      
      HQqgg_1 = HQqgg1


      return
      end



*---------------------------------------------------------------
      real*8 function HQqgg_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 HQqgg2
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2



      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      
      HQqgg2 =
     &  + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ3**2 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2**2 * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ2 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ1 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ3**2 * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2**2 * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ2 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ1 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ2**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13*sQ1*sQ2*sQ3**(-1) * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ2**(-1)*sQ3 * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s13**2*sQ2*sQ3**(-1) * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ2**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*sQ1*sQ2*sQ3**(-1) * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ2**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s23*s13*sQ2*sQ3**(-1) * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ2**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*sQ1*sQ2*sQ3**(-1) * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ2**(-2) * ( 64*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ2**(-1)*sQ3**(-1)
     &  * ( 128*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ2**(-1) * ( 64*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ2**(-1)*sQ3 * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ3**(-2) * ( 64*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ3**(-1) * ( 64*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s13*sQ2*sQ3**(-1) * (  -
     &    16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ2**(-1)*sQ3 * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12*s23*sQ2*sQ3**(-1) * (  -
     &    8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ2**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ2**(-1)*sQ3**(-1)
     &  * ( 64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ2**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ2**(-1)*sQ3 * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ3**(-2) * ( 32*m**4
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ3**(-1) * ( 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-2)*s12**2*sQ2*sQ3**(-1) * (  - 8
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**3*sQ3**(-1)
     &  * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**(-1)*
     & sQ3**2 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ3 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2 * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1)*
     & sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2 * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ2**(-1)
     &  * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ3**(-1)
     &  * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)*
     & sQ3 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1 * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ3**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-2)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-1)
     &  * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     & sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2*sQ2**(-1)*sQ3 * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**2*sQ2*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**3*sQ2**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     & sQ1**3*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23 * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2*sQ3**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     & sQ3**(-2) * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)*
     & sQ3 * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ3**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1 * ( 96 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ2**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     & sQ3**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**3*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ3**(-1)
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ2*
     & sQ3**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ2**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ3**(-2) * (  - 96*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     & sQ3**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     & sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3 * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ3**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1 * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1*sQ2*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**2*sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**2*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     & sQ1**3*sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ3**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ2*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ2**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ3**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*
     & sQ3**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ2**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     & sQ1**2*sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ3**(-2)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**4*s13**(-1)*
     & sQ1*sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ2**(-1)*sQ3**3
     &  * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3**2 * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ3 * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2 * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**2*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*
     & sQ3 * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2 * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ2**(-1)
     &  * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ3**(-1)
     &  * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-2)*sQ3**2 * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3 * (  - 48*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2*sQ3**(-1) * (  - 48*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-2)*sQ3 * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1) * (  - 128*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ3**(-1) * (  - 128*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1 * ( 128 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-1) * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ3**(-1) * ( 80 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-2)*sQ3 * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1)*sQ3**(-1) * (
     &    64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1) * (  - 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1)*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ3**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ3**(-1) * (  - 64*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * (  - 96*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-2) * (  - 96*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ2**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ2**(-1)*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-2)*
     & sQ3**2 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*sQ3
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*
     & sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1) * (  - 8*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3 * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1 * ( 96 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ2**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ2**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     & sQ3**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**3*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-2)*sQ3 * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1)*sQ3**(-1) * (
     &    64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1) * (  - 40*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1)*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ3**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ3**(-1) * (  - 40*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * (  - 128*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1) * ( 96 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-1) * ( 96 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1**2*sQ2**(-1)*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-1)*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ3**(-2) * (  - 96*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-1)
     &  * (  - 24*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     & sQ3**(-1) * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ3**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ2**(-2)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-2)*
     & sQ3**2 * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*
     & sQ3**2 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)
     &  * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)*
     & sQ3 * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1 * ( 64 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ2**(-2) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ2**(-1) * ( 48 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     & sQ3**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-2)*sQ3 * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1)*sQ3**(-1) * (
     &    64*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1) * (  - 64*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1)*sQ3 * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ3**(-2) * ( 32*m**4 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ3**(-1) * (  - 32*m**2
     &     )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2*sQ3**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-2) * (  - 96*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1)*sQ3**(-1)
     &  * (  - 96*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ3**(-2) * (  - 32*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ3**(-1) * ( 40 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-1)*sQ3**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ3**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-1)
     &  * (  - 40*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1) * (  - 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 96*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     & sQ3**(-1) * ( 56 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1**2*
     & sQ2**(-1)*sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-2) * (  - 96*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1)*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ3**(-2) * (  - 48*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ1*sQ2**(-1)*
     & sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*
     & sQ2**(-2) * (  - 48*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 24 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-2)*
     & sQ3 * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)*
     & sQ3 * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ2**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ2**(-1) * ( 32 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     & sQ3**(-1) * ( 16 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-2) * (  - 64*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1)*sQ3**(-1)
     &  * (  - 16*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ3**(-2) * (  - 16*
     &    m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     & sQ2**(-2) * (  - 64*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     & sQ2**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*sQ1*
     & sQ2**(-1)*sQ3**(-1) * ( 8 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-2)
     &  * (  - 32*m**2 )
      HQqgg2 = HQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-1)
     &  * ( 8 )




      HQqgg_2 =  HQqgg2

      return
      end




c
c Begin of soft-virtual functions. See f2svht and finiteHt.m
c
c
      function svfac(s,t,u,m2,mh2,mu2)
c     Note that the normalisation of this function is such 
c     that it needs a 1./4./2./NC (supplied to the Born
c     contribution in fbornht).

      implicit none
      real*8 svfac,s,t,u,m2,mh2,mu2
      real*8 cf,ca,pi
      parameter (cf=4.d0/3.d0)
      parameter (ca=3.d0)
      parameter (pi=3.14159265358979312D0)
      real*8 lnsmbeo2m2,lnspbeo2m2,beta,lnu,lns,lnt,dilg4,dilg3,
     # dilg2,dilg1,lnmu,lnxi,ro,t1,u1,ddilog
      real*8 xicut,deltai,deltao
      real*8 zeta2,Eop0
      common/parsub/xicut,deltai,deltao
      real*8 tmp,svfac2,svfac3

      zeta2=pi**2./6.
      ro=2*(m2+mh2)/s-(m2-mh2)**2/s**2
      beta=sqrt(1-ro) ! Note that E/p0=(s+m^2-mH^2)/beta/s
      Eop0=(s+m2-mh2)/beta/s
      t1=t-m2
      u1=u-m2
      lnxi=log(xicut)
      lns=log(s/m2)  
      lnt=log(-t1/m2)
      lnu=log(-u1/m2)
      lnmu=log(mu2/m2)
      lnsmbeo2m2=log((m2-mh2+s-beta*s)/(2*m2))
      lnspbeo2m2=log((m2-mh2+s+beta*s)/(2*m2))
      dilg1=ddilog((m2+mh2-s-beta*s-2*u)/(m2-mh2+s-beta*s))
      dilg2=ddilog((m2+mh2-s-beta*s-2*u)/(2*m2-2*u))
      dilg3=ddilog((m2+mh2-s+beta*s-2*u)/
     #             (-m2+mh2+(-1+beta)*s))
      dilg4=ddilog((m2+mh2-s+beta*s-2*u)/(2*mh2-2*s-2*u))


      svfac =
     &  + CA * (  - 2*lnxi*lnmu + 2*lnxi**2 - dilg4 - dilg3 + dilg2 + 
     &    dilg1 + lnu*lnsmbeo2m2 - 2*lnu*lnxi - 1./2.*lnu**2 - lnt*
     &    lnsmbeo2m2 + 2*lnt*lnxi + 1./2.*lnt**2 + 2*lns*lnxi - lns*lnu
     &     + lns*lnt + 1./2.*lns**2 - 3./2.*zeta2 )
      svfac = svfac + CF * (  - 3./2.*lnmu + lnspbeo2m2*Eop0 - 1./
     &    2.*lnspbeo2m2**2 - lnsmbeo2m2*Eop0 + lnsmbeo2m2*
     &    lnspbeo2m2 + 1./2.*lnsmbeo2m2**2 - 2*lnxi - 2*lnxi*lnmu + 2*
     &    lnxi**2 + 2*dilg3 - 2*dilg2 - 2*lnu*lnsmbeo2m2 + 4*lnu*lnxi
     &     + lnu**2 - lns + 2*lns*lnu - 1./2.*lns**2 - 3./2.*zeta2 )


c      tmp=svfac2(s,t,u,m2,mh2,mu2,3)

c      write(*,*) 'svfac, svfac2= ',svfac,tmp

      return 
      end 

      function svfac2(s,t,u,m2,mh2,mu2,nf)
c     This is Stefano's result
      implicit none
      real*8 svfac,s,t,u,m2,mh2,mu2,svfac2
      integer nf
      real*8 cf,ca,pi,zeta2,nc
      parameter (cf=4.d0/3.d0)
      parameter (ca=3.d0)
c      parameter (pi=3.14159265358979312D0)
c      parameter (zeta2=pi**2/6.d0)
      parameter (nc=3.d0)
c Set the couplings to one; physical values are assigned elsewhere
      real*8 gs,gw
      parameter (gs=1.d0)
      parameter (gw=1.d0)
      real*8 u1,t1,lnxi,lnmu,lns,lnt,lnu,ro,beta,lnsmbeo2m2,
     # lnspbeo2m2,ddilog,dilg1,dilg2,dilg3,dilg4
      real*8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
c
      u1=u-m2
      t1=t-m2
      lnxi=log(xicut)
      lnmu=log(mu2/m2)
      lns=log(s/m2)
      lnt=log(-t1/m2)
      lnu=log(-u1/m2)
      ro=2*(m2+mh2)/s-(m2-mh2)**2/s**2
      beta=sqrt(1-ro)
      lnsmbeo2m2=log((m2-mh2+s-beta*s)/(2*m2))
      lnspbeo2m2=log((m2-mh2+s+beta*s)/(2*m2))
      dilg1=ddilog((m2+mh2-s-beta*s-2*u)/(m2-mh2+s-beta*s))
      dilg2=ddilog((m2+mh2-s-beta*s-2*u)/(2*m2-2*u))
      dilg3=ddilog((m2+mh2-s+beta*s-2*u)/
     #             (-m2+mh2+(-1+beta)*s))
      dilg4=ddilog((m2+mh2-s+beta*s-2*u)/(2*mh2-2*s-2*u))
c

      pi=0.

      SVFAC = -CF*(LNSMBEO2M2-LNSPBEO2M2)*(S-MH2+M2)/(BETA*S)+(-CA*(PI**
     1   2+2*LNU**2-4*LNSMBEO2M2*LNU+4*LNS*LNU-2*LNT**2+4*LNSMBEO2M2*LNT
     2   -4*LNS*LNT-2*LNS**2+4*DILG4+4*DILG3-4*DILG2-4*DILG1)-CF*(PI**2-
     3   4*LNU**2+8*LNSMBEO2M2*LNU-8*LNS*LNU+2*LNSPBEO2M2**2-4*LNSMBEO2M
     4   2*LNSPBEO2M2-2*LNSMBEO2M2**2+2*LNS**2+4*LNS+6*LNMU-8*DILG3+8*DI
     5   LG2))/4.0D0+2*(CF+CA)*LNXI**2-2*((CA-2*CF)*LNU-CA*LNT-CA*LNS+(C
     6   F+CA)*LNMU+CF)*LNXI

      svfac2=svfac

      return 
      end 



      function svnonfac(s,t,u,m2,mh2,mu2)
c     This is the component of the virtual corrections which is 
c     not proportional to the Born amplitude. Normalisation is such
c     that an overall gs^2*(|a|^2+|b|^2) is needed.

      implicit none
      real*8 svnonfac,s,t,u,m2,mh2,mu2
      real*8 svnonfac_stef
      real*8 cf,ca,pi,zeta2,nc
      parameter (cf=4.d0/3.d0)
      parameter (ca=3.d0)
      parameter (pi=3.14159265358979312D0)
      parameter (zeta2=pi**2/6.d0)
      parameter (nc=3.d0)
c Set the couplings to one; physical values are assigned elsewhere
      real*8 gs,gw,m
      parameter (gs=1.d0)
      parameter (gw=1.d0)
      real*8 cs4,limdel,li2u1,bs6f,li2t1,lit,
     &svnonfac2,svnonfac3,lndel2,Del,lndel,liq,liu
      real*8 q2;
      real*8 elborn,tmp,SVNON
      real*8 LiCS41,LiCS42,LiCS43,LiCS44,LiCS45,lnCS,cs2fb
c
      real*8 lam,t1,u1,delq,lns,lnt,lnu,lnq,xl2tot1,xl2uou1,xl2tom,
     # xl2uom,xl2mwom,xl2omdqom,xl2omdqot1,xl2omdqou1,as1f,bs1f,bs2f,
     # bs3f,bs4f,bs5f,cs1f0,cs1f,cs2f0,cs2f,cs3f0,cs3f,cs4f,cs5f,cs6f,
     # cs7f,cs8f,ds1f0,ds1f,ds2f0,ds2f,ds3f0,ds3f,ddilog,lnmu,tmp5,
     # tmp4,tmp3,tmp2,tmp1,tmp0

      real*8 rflag

c      write(*,*) 's,t,u,m2,mh2= ',s,t,u,m2,mh2
c      stop

c      lam=s**2.+m2**2.+mh2**2.-2.*s*m2-2.*s*mh2-2.*mh2*m2
c      t1=t-m2
c      u1=u-m2
c      delq=(mh2-m2)  ! WAS -(mh2-m2)
c      q2=MH2;

c      lns=log(s/m2)  
c      lnt=log(-t1/m2)
c      lnu=log(-u1/m2)
c      lnq=log(mh2/m2)
c      Del=delq;
c      liu=ddilog(u/m2)
c      lit=ddilog(t/m2)
      
c      liq=ddilog(mh2/m2)
c      limdel=ddilog(-delq/m2) 
c      li2u1=ddilog(1.-u1/delq) 
c      li2t1=ddilog(1.-t1/delq)

      rflag=1. ! =1 for on-shell Yukawa, 0 for MSbar scheme

      liu=ddilog(u/m2)
      q2=mh2
      lam=s**2.+m2**2.+mh2**2.-2.*s*m2-2.*s*mh2-2.*mh2*m2
      t1=t-m2
      u1=u-m2
      delq=mh2-m2
      lns=log(s/m2)
      lnt=log(-t1/m2)
      lnu=log(-u1/m2)
      lnmu=log(mu2/m2)
      if(mh2.gt.m2) THEN
         lnq=log(delq/m2)
      else
         lnq=log(-delq/m2)
      endif
      xl2tot1=ddilog(t/t1)
      xl2uou1=ddilog(u/u1)
      xl2tom=ddilog(t/m2)
      xl2uom=ddilog(u/m2)
      xl2mwom=ddilog(mh2/m2)
      xl2omdqom=ddilog(1-delq/m2)
      xl2omdqot1=ddilog(1-delq/t1)
      xl2omdqou1=ddilog(1-delq/u1)
c As before, but without the 1/(16*pi**2). Rename it cs4f
      CS4f=-1./sqrt(lam)*(2*ddilog(1.-sqrt(lam)/s)+ddilog((del
     &q**2.+delq*(s+sqrt(lam))-2*s*mh2)/(delq**2.+delq*(s-sqrt(lam))-2*s
     &*mh2))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sqrt(lam)))/(d
     &elq**2.+delq*(s-sqrt(lam))-2*s*mh2))-2*ddilog((s-sqrt(lam)-delq)/(
     &s+sqrt(lam)-delq))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sq
     &rt(lam)))/(delq**2.+delq*(s+sqrt(lam))-2*s*mh2))-1./2.*(log((s-sqr
     &t(lam)-delq)/(s+sqrt(lam)-delq)))**2.+zeta2)

c Other finite parts
c Equation numbers refer to Ellis & Zanderighi
      as1f = m2   ! Y
      BS1F  = 2-LNS  ! Y
      BS2F  = 2 ! Y
      BS3F  = 2-LNU*U1/U  ! Y
      BS4F  = 2-DELQ*LNQ/MH2  ! Y
      BS5F  = 2-LNT*T1/T   !  Y
      CS1F0 = (LNS**2-PI**2)/S/2.0D0  ! Eq. (4.5) 
      CS1F = CS1F0-PI**2/S/1.2D1 ! From r_G
      CS2F0 = (-XL2UOU1+PI**2/1.2D1+LNU**2/2.0D0)/U1 ! Eq. (4.11)
      CS2fb = 1./u1*(lnu**2+Liu+1./4.*zeta2)
      CS2F = CS2F0-PI**2/U1/2.4D1 ! From r_G
c      write(*,*) 'CS2f, CS2fb= ',CS2f,CS2fb
c      write(*,*) '(CS2f-CS2fb)*u1/zeta2= ',(CS2f-CS2fb)*u1/zeta2
      CS3F0 = (-XL2TOT1+PI**2/1.2D1+LNT**2/2.0D0)/T1 ! Y
      CS3F = CS3F0-PI**2/T1/2.4D1  ! Y
      CS5F  = (-XL2UOM+XL2MWOM-LNU**2+LNQ**2)/(DELQ-U1) 
      CS5f=CS5f-PI**2/(DELQ-U1) ! Added by Chris
      CS6F  = (-XL2TOM+XL2MWOM-LNT**2+LNQ**2)/(DELQ-T1)
      CS6f=CS6f-PI**2/(DELQ-T1) ! Added by Chris

      cs7f = (zeta2-xl2tom)/t1
      cs8f = (xl2mwom-xl2uom)/(u1-delq)
      DS1F0 = (-2*XL2OMDQOU1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNU-LNQ**2)/(S*U
     1   1)
      DS1F = DS1F0-PI**2/(S*U1)/8.0D0
      DS1F=DS1F+PI**2/S/U1 ! Added by Chris
      DS2F0 = (-2*XL2OMDQOT1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNT-LNQ**2)/(S*T
     1   1)
      DS2F = DS2F0-PI**2/(S*T1)/8.0D0
      DS2F=DS2F+PI**2/S/T1
      DS3F0 = (-2*XL2OMDQOU1-2*XL2OMDQOT1-PI**2/1.2D1+2*LNT*LNU-LNQ**2)/
     1   (T1*U1)
      DS3F = DS3F0-PI**2/(T1*U1)/2.4D1
      DS3F=DS3F+PI**2/T1/U1
      if (mh2.gt.m2) then ! High Higgs mass
c      lndel=log(delq/m2)
c      lndel2=lndel**2

      LiCS41=ddilog(1.-sqrt(lam)/s)
      LiCS42=ddilog((delq**2+delq*(s+sqrt(lam))-2.*s*q2)
     &/(delq**2+delq*(s-sqrt(lam))-2.*s*q2))
      LiCS43=ddilog((delq**2+delq*(s+sqrt(lam))-2.*q2*(s+sqrt(lam)))
     &/(delq**2+delq*(s-sqrt(lam))-2.*q2*s))
      LiCS44=ddilog((s-sqrt(lam)-delq)/(s+sqrt(lam)-delq))
      LiCS45=ddilog((delq**2+delq*(s+sqrt(lam))-2.*q2*(s+sqrt(lam)))
     &/(delq**2+delq*(s+sqrt(lam))-2.*q2*s))
      lnCS=log((s-sqrt(lam)-delq)/(s+sqrt(lam)-delq))

      else ! Low Higgs mass      

c As before, but without the 1/(16*pi**2). Rename it cs4f
      cs4f=-1./sqrt(lam)*(2*ddilog(1.-sqrt(lam)/s)+ddilog((del
     &q**2.+delq*(s+sqrt(lam))-2*s*mh2)/(delq**2.+delq*(s-sqrt(lam))-2*s
     &*mh2))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sqrt(lam)))/(d
     &elq**2.+delq*(s-sqrt(lam))-2*s*mh2))-2*ddilog((s-sqrt(lam)-delq)/(
     &s+sqrt(lam)-delq))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sq
     &rt(lam)))/(delq**2.+delq*(s+sqrt(lam))-2*s*mh2))-1./2.*(log((s-sqr
     &t(lam)-delq)/(s+sqrt(lam)-delq)))**2.+zeta2)
c Other finite parts
      as1f = m2
      BS1F  = 2-LNS
      BS2F  = 2
      BS3F  = 2-LNU*U1/U
      BS4F  = 2-DELQ*LNQ/Mh2
      BS5F  = 2-LNT*T1/T
      CS1F0 = (LNS**2-PI**2)/S/2.0D0
      CS1F = CS1F0-PI**2/S/1.2D1
      CS2F0 = (-XL2UOU1+PI**2/1.2D1+LNU**2/2.0D0)/U1
      CS2F = CS2F0-PI**2/U1/2.4D1
      CS3F0 = (-XL2TOT1+PI**2/1.2D1+LNT**2/2.0D0)/T1
      CS3F = CS3F0-PI**2/T1/2.4D1
      CS5F  = (-XL2UOM+XL2MWOM-LNU**2+LNQ**2)/(DELQ-U1)
      CS6F  = (-XL2TOM+XL2MWOM-LNT**2+LNQ**2)/(DELQ-T1)
      cs7f = (zeta2-xl2tom)/t1
      cs8f = (xl2mwom-xl2uom)/(u1-delq)
      DS1F0 = (-2*XL2OMDQOU1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNU-LNQ**2)/(S*U
     1   1)
      DS1F = DS1F0-PI**2/(S*U1)/8.0D0
      DS2F0 = (-2*XL2OMDQOT1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNT-LNQ**2)/(S*T
     1   1)
      DS2F = DS2F0-PI**2/(S*T1)/8.0D0
      DS3F0 = (-2*XL2OMDQOU1-2*XL2OMDQOT1-PI**2/1.2D1+2*LNT*LNU-LNQ**2)/
     1   (T1*U1)
      DS3F = DS3F0-PI**2/(T1*U1)/2.4D1

      endif

      m=sqrt(m2)


      svnonfac =
     &  + CA*NC**(-1) * ( 1./4.*s**(-1)*lam**(-1)*mh2**3*BS4f - 1./4.*
     &    s**(-1)*lam**(-1)*mh2**3*BS2f - 1./2.*s**(-1)*mh2 - 1./4.*
     &    s**(-1)*mh2*BS4f + 1./4.*s**(-1)*mh2*BS2f - 1./2.*s**(-1)*
     &    mh2**2*CS8f + 1./2.*s**(-1)*mh2**2*CS7f - 1./2.*s**(-1)*
     &    mh2**2*CS5f + 1./2.*s**(-1)*mh2**2*CS3f - 1./2.*s**(-1)*
     &    mh2**3*DS3f - 1./4.*s**(-1)*t1*lam**(-1)*mh2**2*BS4f + 1./4.*
     &    s**(-1)*t1*lam**(-1)*mh2**2*BS2f + 1./4.*s**(-1)*t1*BS4f - 1./
     &    4.*s**(-1)*t1*BS2f + 1./2.*s**(-1)*t1*mh2*CS8f - 1./2.*
     &    s**(-1)*t1*mh2*CS7f + 1./2.*s**(-1)*t1*mh2*CS5f - 1./2.*
     &    s**(-1)*t1*mh2*CS3f + s**(-1)*t1*mh2**2*DS3f - 1./4.*s**(-1)*
     &    t1**2*CS8f + 1./4.*s**(-1)*t1**2*CS7f - 1./4.*s**(-1)*t1**2*
     &    CS5f + 1./4.*s**(-1)*t1**2*CS3f - 3./4.*s**(-1)*t1**2*mh2*
     &    DS3f + 1./4.*s**(-1)*t1**3*DS3f - 1./8.*t1**(-1)*lam**(-1)*
     &    mh2**3*BS4f + 1./8.*t1**(-1)*lam**(-1)*mh2**3*BS2f - 1./2.*
     &    t1**(-1)*mh2 - 1./2.*t1**(-1)*mh2*BS5f + 13./8.*t1**(-1)*mh2*
     &    BS4f )
      svnonfac = svnonfac + CA*NC**(-1) * (  - t1**(-1)*mh2*BS3f - 1./8.
     &    *t1**(-1)*mh2*BS2f - 1./2.*t1**(-1)*mh2**2*CS8f - 1./2.*
     &    t1**(-1)*mh2**2*CS5f + t1**(-1)*mh2**2*CS1f - 1./2.*t1**(-1)*
     &    mh2**3*DS1f - 5./8.*lam**(-1)*mh2**2*BS4f + 5./8.*lam**(-1)*
     &    mh2**2*BS2f + 1./4.*lam**(-1)*mh2**3*CS4f - 7./8.*BS4f + BS3f
     &     - 1./8.*BS2f + mh2*CS8f - 1./2.*mh2*CS7f + mh2*CS5f - 1./4.*
     &    mh2*CS4f - 1./2.*mh2*CS3f - mh2*CS1f + mh2**2*DS3f + 1./2.*
     &    mh2**2*DS2f + mh2**2*DS1f + 1./2.*t1*lam**(-1)*mh2*BS4f - 1./
     &    2.*t1*lam**(-1)*mh2*BS2f - 1./4.*t1*lam**(-1)*mh2**2*CS4f - 3.
     &    /4.*t1*CS8f + 1./2.*t1*CS7f - 3./4.*t1*CS5f + 1./4.*t1*CS4f
     &     + 1./2.*t1*CS3f + 1./2.*t1*CS1f - 3./2.*t1*mh2*DS3f - 1./2.*
     &    t1*mh2*DS2f - 3./4.*t1*mh2*DS1f + 3./4.*t1**2*DS3f + 1./4.*
     &    t1**2*DS2f + 1./4.*t1**2*DS1f + 3./8.*s*t1**(-1)*lam**(-1)*
     &    mh2**2*BS4f - 3./8.*s*t1**(-1)*lam**(-1)*mh2**2*BS2f - 1./8.*
     &    s*t1**(-1)*lam**(-1)*mh2**3*CS4f + 1./2.*s*t1**(-1)*BS5f - 13.
     &    /8.*s*t1**(-1)*BS4f )
      svnonfac = svnonfac + CA*NC**(-1) * ( s*t1**(-1)*BS3f + 1./8.*s*
     &    t1**(-1)*BS2f + 1./2.*s*t1**(-1)*mh2*CS8f + 1./2.*s*t1**(-1)*
     &    mh2*CS5f + 1./8.*s*t1**(-1)*mh2*CS4f - s*t1**(-1)*mh2*CS1f + 
     &    s*t1**(-1)*mh2**2*DS1f + 1./2.*s*lam**(-1)*mh2*BS4f - 1./2.*s
     &    *lam**(-1)*mh2*BS2f - 5./8.*s*lam**(-1)*mh2**2*CS4f - 3./4.*s
     &    *CS8f + 1./4.*s*CS7f - 3./4.*s*CS5f + 1./8.*s*CS4f + 1./4.*s*
     &    CS3f + s*CS1f - 3./4.*s*mh2*DS3f - 1./2.*s*mh2*DS2f - 3./2.*s
     &    *mh2*DS1f - 1./4.*s*t1*lam**(-1)*BS4f + 1./4.*s*t1*lam**(-1)*
     &    BS2f + 1./2.*s*t1*lam**(-1)*mh2*CS4f + 3./4.*s*t1*DS3f + 1./2.
     &    *s*t1*DS2f + 3./4.*s*t1*DS1f - 3./8.*s**2*t1**(-1)*lam**(-1)*
     &    mh2*BS4f + 3./8.*s**2*t1**(-1)*lam**(-1)*mh2*BS2f + 3./8.*
     &    s**2*t1**(-1)*lam**(-1)*mh2**2*CS4f - 1./4.*s**2*t1**(-1)*
     &    CS8f - 1./4.*s**2*t1**(-1)*CS5f - 1./8.*s**2*t1**(-1)*CS4f + 
     &    1./2.*s**2*t1**(-1)*CS1f - 3./4.*s**2*t1**(-1)*mh2*DS1f - 1./
     &    8.*s**2*lam**(-1)*BS4f + 1./8.*s**2*lam**(-1)*BS2f + 1./2.*
     &    s**2*lam**(-1)*mh2*CS4f )
      svnonfac = svnonfac + CA*NC**(-1) * ( 1./4.*s**2*DS3f + 1./4.*
     &    s**2*DS2f + 3./4.*s**2*DS1f - 1./4.*s**2*t1*lam**(-1)*CS4f + 
     &    1./8.*s**3*t1**(-1)*lam**(-1)*BS4f - 1./8.*s**3*t1**(-1)*
     &    lam**(-1)*BS2f - 3./8.*s**3*t1**(-1)*lam**(-1)*mh2*CS4f + 1./
     &    4.*s**3*t1**(-1)*DS1f - 1./8.*s**3*lam**(-1)*CS4f + 1./8.*
     &    s**4*t1**(-1)*lam**(-1)*CS4f - 1./4.*m**2*s**(-1)*lam**(-1)*
     &    mh2**2*BS4f + 1./4.*m**2*s**(-1)*lam**(-1)*mh2**2*BS2f + 1./2.
     &    *m**2*s**(-1) - 1./4.*m**2*s**(-1)*BS4f + 1./4.*m**2*s**(-1)*
     &    BS2f + m**2*s**(-1)*mh2*CS8f - m**2*s**(-1)*mh2*CS7f + 1./4.*
     &    m**2*s**(-1)*mh2*CS6f + m**2*s**(-1)*mh2*CS5f - m**2*s**(-1)*
     &    mh2*CS3f - 1./4.*m**2*s**(-1)*mh2*CS2f + 3./2.*m**2*s**(-1)*
     &    mh2**2*DS3f + 1./2.*m**2*s**(-1)*t1*lam**(-1)*mh2*BS4f - 1./2.
     &    *m**2*s**(-1)*t1*lam**(-1)*mh2*BS2f - 1./4.*m**2*s**(-1)*t1*
     &    CS8f + 1./4.*m**2*s**(-1)*t1*CS7f - 1./4.*m**2*s**(-1)*t1*
     &    CS6f - 1./2.*m**2*s**(-1)*t1*CS5f + 1./2.*m**2*s**(-1)*t1*
     &    CS3f )
      svnonfac = svnonfac + CA*NC**(-1) * ( 1./4.*m**2*s**(-1)*t1*CS2f
     &     - 7./4.*m**2*s**(-1)*t1*mh2*DS3f + 1./2.*m**2*s**(-1)*t1**2*
     &    DS3f - 1./2.*m**2*t1**(-2)*mh2*BS5f + 1./2.*m**2*t1**(-2)*mh2
     &    *BS2f + 3./8.*m**2*t1**(-1)*lam**(-1)*mh2**2*BS4f - 3./8.*
     &    m**2*t1**(-1)*lam**(-1)*mh2**2*BS2f + m**2*t1**(-1) + m**2*
     &    t1**(-1)*BS5f - 13./8.*m**2*t1**(-1)*BS4f + m**2*t1**(-1)*
     &    BS3f - 3./8.*m**2*t1**(-1)*BS2f + m**2*t1**(-1)*mh2*CS8f - 1./
     &    2.*m**2*t1**(-1)*mh2*CS7f + m**2*t1**(-1)*mh2*CS5f + 1./4.*
     &    m**2*t1**(-1)*mh2*CS4f - 1./4.*m**2*t1**(-1)*mh2*CS2f - 2*
     &    m**2*t1**(-1)*mh2*CS1f + m**2*t1**(-1)*mh2**2*DS3f + 3./2.*
     &    m**2*t1**(-1)*mh2**2*DS1f - 3./4.*m**2*lam**(-1)*mh2*BS4f + 3.
     &    /4.*m**2*lam**(-1)*mh2*BS2f - 1./4.*m**2*lam**(-1)*mh2**2*
     &    CS4f - 1./4.*m**2*CS8f - 3./4.*m**2*CS5f - 3./4.*m**2*CS4f + 
     &    1./2.*m**2*CS3f + 1./2.*m**2*CS2f + m**2*CS1f - 3*m**2*mh2*
     &    DS3f - m**2*mh2*DS2f - 2*m**2*mh2*DS1f + 1./2.*m**2*t1*
     &    lam**(-1)*BS4f )
      svnonfac = svnonfac + CA*NC**(-1) * (  - 1./2.*m**2*t1*lam**(-1)*
     &    BS2f + 1./2.*m**2*t1*lam**(-1)*mh2*CS4f + 5./4.*m**2*t1*DS3f
     &     + 1./2.*m**2*t1*DS2f + 3./4.*m**2*t1*DS1f + 1./2.*m**2*s*
     &    t1**(-2) + 1./2.*m**2*s*t1**(-2)*BS5f - 1./2.*m**2*s*t1**(-2)
     &    *BS2f - 1./4.*m**2*s*t1**(-1)*lam**(-1)*mh2*BS4f + 1./4.*m**2
     &    *s*t1**(-1)*lam**(-1)*mh2*BS2f + 3./8.*m**2*s*t1**(-1)*
     &    lam**(-1)*mh2**2*CS4f - 1./4.*m**2*s*t1**(-1)*CS5f - 3./8.*
     &    m**2*s*t1**(-1)*CS4f + 1./4.*m**2*s*t1**(-1)*CS2f + 3./4.*
     &    m**2*s*t1**(-1)*CS1f - m**2*s*t1**(-1)*mh2*DS3f - m**2*s*
     &    t1**(-1)*mh2*DS2f - 7./4.*m**2*s*t1**(-1)*mh2*DS1f + 1./2.*
     &    m**2*s*lam**(-1)*BS4f - 1./2.*m**2*s*lam**(-1)*BS2f - 3./4.*
     &    m**2*s*lam**(-1)*mh2*CS4f + 3./4.*m**2*s*DS3f + 1./2.*m**2*s*
     &    DS2f + 5./4.*m**2*s*DS1f + 1./2.*m**2*s*t1*lam**(-1)*CS4f - 1.
     &    /8.*m**2*s**2*t1**(-1)*lam**(-1)*BS4f + 1./8.*m**2*s**2*
     &    t1**(-1)*lam**(-1)*BS2f - 1./4.*m**2*s**2*t1**(-1)*lam**(-1)*
     &    mh2*CS4f )
      svnonfac = svnonfac + CA*NC**(-1) * ( 1./2.*m**2*s**2*t1**(-1)*
     &    DS1f + 1./2.*m**2*s**2*lam**(-1)*CS4f - 1./8.*m**2*s**3*
     &    t1**(-1)*lam**(-1)*CS4f - 1./4.*m**4*s**(-1)*lam**(-1)*mh2*
     &    BS4f + 1./4.*m**4*s**(-1)*lam**(-1)*mh2*BS2f - 1./2.*m**4*
     &    s**(-1)*CS8f + 1./2.*m**4*s**(-1)*CS7f - 1./4.*m**4*s**(-1)*
     &    CS6f - 1./2.*m**4*s**(-1)*CS5f + 1./2.*m**4*s**(-1)*CS3f + 1./
     &    4.*m**4*s**(-1)*CS2f - 3./2.*m**4*s**(-1)*mh2*DS3f - 1./4.*
     &    m**4*s**(-1)*t1*lam**(-1)*BS4f + 1./4.*m**4*s**(-1)*t1*
     &    lam**(-1)*BS2f + 3./4.*m**4*s**(-1)*t1*DS3f + 1./2.*m**4*
     &    t1**(-2)*BS5f - 1./2.*m**4*t1**(-2)*BS2f - 3./8.*m**4*
     &    t1**(-1)*lam**(-1)*mh2*BS4f + 3./8.*m**4*t1**(-1)*lam**(-1)*
     &    mh2*BS2f - 1./2.*m**4*t1**(-1)*CS8f + 1./2.*m**4*t1**(-1)*
     &    CS7f - 1./2.*m**4*t1**(-1)*CS5f - 1./4.*m**4*t1**(-1)*CS4f + 
     &    1./4.*m**4*t1**(-1)*CS2f + m**4*t1**(-1)*CS1f - 2*m**4*
     &    t1**(-1)*mh2*DS3f - 3./2.*m**4*t1**(-1)*mh2*DS1f - 5./8.*m**4
     &    *lam**(-1)*BS4f )
      svnonfac = svnonfac + CA*NC**(-1) * ( 5./8.*m**4*lam**(-1)*BS2f
     &     - 1./4.*m**4*lam**(-1)*mh2*CS4f + 2*m**4*DS3f + 1./2.*m**4*
     &    DS2f + m**4*DS1f - 1./4.*m**4*t1*lam**(-1)*CS4f + 1./2.*m**4*
     &    s*t1**(-3)*BS5f - 1./2.*m**4*s*t1**(-3)*BS2f + 1./2.*m**4*s*
     &    t1**(-2)*CS7f - 1./8.*m**4*s*t1**(-1)*lam**(-1)*BS4f + 1./8.*
     &    m**4*s*t1**(-1)*lam**(-1)*BS2f - 3./8.*m**4*s*t1**(-1)*
     &    lam**(-1)*mh2*CS4f + m**4*s*t1**(-1)*DS3f + m**4*s*t1**(-1)*
     &    DS2f + 3./4.*m**4*s*t1**(-1)*DS1f - 5./8.*m**4*s*lam**(-1)*
     &    CS4f - 1./8.*m**4*s**2*t1**(-1)*lam**(-1)*CS4f + 1./4.*m**6*
     &    s**(-1)*lam**(-1)*BS4f - 1./4.*m**6*s**(-1)*lam**(-1)*BS2f + 
     &    1./2.*m**6*s**(-1)*DS3f + 1./8.*m**6*t1**(-1)*lam**(-1)*BS4f
     &     - 1./8.*m**6*t1**(-1)*lam**(-1)*BS2f + m**6*t1**(-1)*DS3f + 
     &    1./2.*m**6*t1**(-1)*DS1f + 1./4.*m**6*lam**(-1)*CS4f + 1./8.*
     &    m**6*s*t1**(-1)*lam**(-1)*CS4f + 1./2./(u1 + s)*t1**(-1)*
     &    mh2**2*BS5f - 1./2./(u1 + s)*t1**(-1)*mh2**2*BS4f - 1./2./(u1
     &     + s)*mh2*BS5f )
      svnonfac = svnonfac + CA*NC**(-1) * ( 1./2./(u1 + s)*mh2*BS4f - 1.
     &    /2./(u1 + s)*s*t1**(-1)*mh2*BS5f + 1./2./(u1 + s)*s*t1**(-1)*
     &    mh2*BS4f + 1./2./(u1 + s)*s*BS5f - 1./2./(u1 + s)*s*BS4f - 
     &    1/(u1 + s)*m**2*t1**(-1)*mh2*BS5f + 1/(u1 + s)*m**2*t1**(-1)*
     &    mh2*BS4f + 1./2./(u1 + s)*m**2*BS5f - 1./2./(u1 + s)*m**2*
     &    BS4f + 1./2./(u1 + s)*m**2*s*t1**(-1)*BS5f - 1./2./(u1 + s)*
     &    m**2*s*t1**(-1)*BS4f + 1./2./(u1 + s)*m**4*t1**(-1)*BS5f - 1./
     &    2./(u1 + s)*m**4*t1**(-1)*BS4f - 1/(t1 + s)*mh2*BS4f + 1/(t1
     &     + s)*mh2*BS3f + 1/(t1 + s)*t1*BS4f - 1/(t1 + s)*t1*BS3f - 
     &    1/(t1 + s)*s*t1**(-1)*mh2*BS4f + 1/(t1 + s)*s*t1**(-1)*mh2*
     &    BS3f + 2/(t1 + s)*s*BS4f - 2/(t1 + s)*s*BS3f + 1/(t1 + s)*
     &    s**2*t1**(-1)*BS4f - 1/(t1 + s)*s**2*t1**(-1)*BS3f + 1/(t1 + 
     &    s)*m**2*BS4f - 1/(t1 + s)*m**2*BS3f + 1/(t1 + s)*m**2*s*
     &    t1**(-1)*BS4f - 1/(t1 + s)*m**2*s*t1**(-1)*BS3f )
      svnonfac = svnonfac + CF*NC**(-1) * ( 4 + s**(-1)*t**(-1)*mh2*
     &    AS1f + s**(-1)*t**(-1)*t1*mh2*BS5f - 1./2.*s**(-1)*t**(-1)*
     &    t1**2*BS5f - s**(-1)*t1**(-2)*mh2**2*AS1f + 4*s**(-1)*
     &    t1**(-1)*mh2**2 + 4*s**(-1)*t1**(-1)*mh2**2*rflag - 3*s**(-1)
     &    *t1**(-1)*mh2**2*lnmu - 2*s**(-1)*t1**(-1)*mh2**2*BS2f + 
     &    s**(-1)*t1**(-1)*mh2**3*CS6f + s**(-1)*t1**(-1)*mh2**3*CS4f
     &     - 7./2.*s**(-1)*mh2 - 4*s**(-1)*mh2*rflag + 3*s**(-1)*mh2*
     &    lnmu - s**(-1)*mh2*BS5f + 2*s**(-1)*mh2*BS2f + s**(-1)*mh2**2
     &    *CS8f - s**(-1)*mh2**2*CS7f - 2*s**(-1)*mh2**2*CS6f + s**(-1)
     &    *mh2**2*CS5f - s**(-1)*mh2**2*CS4f + s**(-1)*mh2**3*DS3f + 2*
     &    s**(-1)*t1 + 2*s**(-1)*t1*rflag - 3./2.*s**(-1)*t1*lnmu + 1./
     &    2.*s**(-1)*t1*BS5f - s**(-1)*t1*BS2f - s**(-1)*t1*mh2*CS8f + 
     &    s**(-1)*t1*mh2*CS7f + 3./2.*s**(-1)*t1*mh2*CS6f - s**(-1)*t1*
     &    mh2*CS5f + 1./2.*s**(-1)*t1*mh2*CS4f - 2*s**(-1)*t1*mh2**2*
     &    DS3f + 1./2.*s**(-1)*t1**2*CS8f - 1./2.*s**(-1)*t1**2*CS7f - 
     &    1./2.*s**(-1)*t1**2*CS6f )
      svnonfac = svnonfac + CF*NC**(-1) * ( 1./2.*s**(-1)*t1**2*CS5f + 
     &    3./2.*s**(-1)*t1**2*mh2*DS3f - 1./2.*s**(-1)*t1**3*DS3f - 1./
     &    2.*t**(-1)*t1**(-1)*mh2*AS1f + 1./2.*t**(-1)*AS1f - 1./2.*
     &    t**(-1)*mh2*BS5f - 1./2.*t**(-1)*t1*BS5f + 3./2.*t1**(-2)*mh2
     &    *AS1f - t1**(-1)*AS1f - 7./2.*t1**(-1)*mh2 - 4*t1**(-1)*mh2*
     &    rflag + 3*t1**(-1)*mh2*lnmu - 1./2.*t1**(-1)*mh2*BS5f - 
     &    t1**(-1)*mh2*BS4f + 2*t1**(-1)*mh2*BS3f + 2*t1**(-1)*mh2*BS2f
     &     + t1**(-1)*mh2**2*CS8f - t1**(-1)*mh2**2*CS6f + t1**(-1)*
     &    mh2**2*CS5f - 2*t1**(-1)*mh2**2*CS4f - t1**(-1)*mh2**2*CS1f
     &     + t1**(-1)*mh2**3*DS1f + 4*rflag - 3*lnmu + 1./2.*BS5f + 2*
     &    BS4f - 2*BS3f - 2*BS2f - 2*mh2*CS8f + mh2*CS7f + 2*mh2*CS6f
     &     - 2*mh2*CS5f + 2*mh2*CS4f + mh2*CS1f - 2*mh2**2*DS3f - 2*
     &    mh2**2*DS1f + 3./2.*t1*CS8f - t1*CS7f - t1*CS6f + 3./2.*t1*
     &    CS5f - 1./2.*t1*CS4f - 1./2.*t1*CS1f + 3*t1*mh2*DS3f + 3./2.*
     &    t1*mh2*DS1f - 3./2.*t1**2*DS3f - 1./2.*t1**2*DS1f + 2*s*
     &    t1**(-1) )
      svnonfac = svnonfac + CF*NC**(-1) * ( 2*s*t1**(-1)*rflag - 3./2.*
     &    s*t1**(-1)*lnmu + 2*s*t1**(-1)*BS4f - 2*s*t1**(-1)*BS3f - s*
     &    t1**(-1)*BS2f - s*t1**(-1)*mh2*CS8f + 1./2.*s*t1**(-1)*mh2*
     &    CS6f - s*t1**(-1)*mh2*CS5f + 3./2.*s*t1**(-1)*mh2*CS4f + s*
     &    t1**(-1)*mh2*CS1f - 2*s*t1**(-1)*mh2**2*DS1f + 3./2.*s*CS8f
     &     - 1./2.*s*CS7f - 1./2.*s*CS6f + 3./2.*s*CS5f - s*CS4f - s*
     &    CS1f + 3./2.*s*mh2*DS3f + 3*s*mh2*DS1f - 3./2.*s*t1*DS3f - 3./
     &    2.*s*t1*DS1f + 1./2.*s**2*t1**(-1)*CS8f + 1./2.*s**2*t1**(-1)
     &    *CS5f - 1./2.*s**2*t1**(-1)*CS4f - 1./2.*s**2*t1**(-1)*CS1f
     &     + 3./2.*s**2*t1**(-1)*mh2*DS1f - 1./2.*s**2*DS3f - 3./2.*
     &    s**2*DS1f - 1./2.*s**3*t1**(-1)*DS1f + m**2*s**(-1)*t**(-1)*
     &    t1**(-1)*mh2*AS1f - 1./2.*m**2*s**(-1)*t**(-1)*AS1f + m**2*
     &    s**(-1)*t**(-1)*mh2*BS5f - 3./2.*m**2*s**(-1)*t**(-1)*t1*BS5f
     &     + 2*m**2*s**(-1)*t1**(-2)*mh2*AS1f + 5*m**2*s**(-1)*t1**(-2)
     &    *mh2**2 - 2*m**2*s**(-1)*t1**(-2)*mh2**2*BS2f + 1./2.*m**2*
     &    s**(-1)*t1**(-1)*lam**(-1)*mh2**3*BS4f )
      svnonfac = svnonfac + CF*NC**(-1) * (  - 1./2.*m**2*s**(-1)*
     &    t1**(-1)*lam**(-1)*mh2**3*BS2f - 1./2.*m**2*s**(-1)*t1**(-1)*
     &    AS1f - 13*m**2*s**(-1)*t1**(-1)*mh2 - 8*m**2*s**(-1)*t1**(-1)
     &    *mh2*rflag + 6*m**2*s**(-1)*t1**(-1)*mh2*lnmu + 3./2.*m**2*
     &    s**(-1)*t1**(-1)*mh2*BS4f + 9./2.*m**2*s**(-1)*t1**(-1)*mh2*
     &    BS2f - 3*m**2*s**(-1)*t1**(-1)*mh2**2*CS6f - 3*m**2*s**(-1)*
     &    t1**(-1)*mh2**2*CS4f - 2*m**2*s**(-1)*lam**(-1)*mh2**2*BS4f
     &     + 2*m**2*s**(-1)*lam**(-1)*mh2**2*BS2f + 7./2.*m**2*s**(-1)
     &     + 4*m**2*s**(-1)*rflag - 3*m**2*s**(-1)*lnmu + m**2*s**(-1)*
     &    BS5f - 2*m**2*s**(-1)*BS2f - 2*m**2*s**(-1)*mh2*CS8f + 2*m**2
     &    *s**(-1)*mh2*CS7f + 7./2.*m**2*s**(-1)*mh2*CS6f - 2*m**2*
     &    s**(-1)*mh2*CS5f + 2*m**2*s**(-1)*mh2*CS4f + 1./2.*m**2*
     &    s**(-1)*mh2*CS2f - 3*m**2*s**(-1)*mh2**2*DS3f + m**2*s**(-1)*
     &    t1*lam**(-1)*mh2*BS4f - m**2*s**(-1)*t1*lam**(-1)*mh2*BS2f + 
     &    1./2.*m**2*s**(-1)*t1*CS8f - 1./2.*m**2*s**(-1)*t1*CS7f - 
     &    m**2*s**(-1)*t1*CS6f )
      svnonfac = svnonfac + CF*NC**(-1) * ( m**2*s**(-1)*t1*CS5f - 1./2.
     &    *m**2*s**(-1)*t1*CS4f - 1./2.*m**2*s**(-1)*t1*CS2f + 7./2.*
     &    m**2*s**(-1)*t1*mh2*DS3f - m**2*s**(-1)*t1**2*DS3f + 1./2.*
     &    m**2*t**(-1)*t1**(-1)*AS1f - 1./2.*m**2*t**(-1)*BS5f + 2*m**2
     &    *t1**(-3)*mh2*AS1f - 5./2.*m**2*t1**(-2)*AS1f - 23./2.*m**2*
     &    t1**(-2)*mh2 - 4*m**2*t1**(-2)*mh2*rflag + 3*m**2*t1**(-2)*
     &    mh2*lnmu + m**2*t1**(-2)*mh2*BS5f + 4*m**2*t1**(-2)*mh2*BS2f
     &     - 2*m**2*t1**(-2)*mh2**2*CS6f - m**2*t1**(-1)*lam**(-1)*
     &    mh2**2*BS4f + 3./2.*m**2*t1**(-1)*lam**(-1)*mh2**2*BS2f - 1./
     &    2.*m**2*t1**(-1)*lam**(-1)*mh2**2*BS1f + 1./2.*m**2*t1**(-1)*
     &    lam**(-1)*mh2**3*CS4f + 8*m**2*t1**(-1) + 4*m**2*t1**(-1)*
     &    rflag - 3*m**2*t1**(-1)*lnmu - 1./2.*m**2*t1**(-1)*BS5f + 
     &    m**2*t1**(-1)*BS4f - 2*m**2*t1**(-1)*BS3f - 7./2.*m**2*
     &    t1**(-1)*BS2f + 1./2.*m**2*t1**(-1)*BS1f - 2*m**2*t1**(-1)*
     &    mh2*CS8f + m**2*t1**(-1)*mh2*CS7f + 4*m**2*t1**(-1)*mh2*CS6f
     &     - 2*m**2*t1**(-1)*mh2*CS5f )
      svnonfac = svnonfac + CF*NC**(-1) * ( m**2*t1**(-1)*mh2*CS4f + 1./
     &    2.*m**2*t1**(-1)*mh2*CS2f + 2*m**2*t1**(-1)*mh2*CS1f - 2*m**2
     &    *t1**(-1)*mh2**2*DS3f - 3*m**2*t1**(-1)*mh2**2*DS1f + 2*m**2*
     &    lam**(-1)*mh2*BS4f - 4*m**2*lam**(-1)*mh2*BS2f + 2*m**2*
     &    lam**(-1)*mh2*BS1f - 2*m**2*lam**(-1)*mh2**2*CS4f + 1./2.*
     &    m**2*CS8f - 2*m**2*CS6f + 3./2.*m**2*CS5f - m**2*CS2f - m**2*
     &    CS1f + 6*m**2*mh2*DS3f + 4*m**2*mh2*DS1f + m**2*t1*lam**(-1)*
     &    BS4f + m**2*t1*lam**(-1)*BS2f - 2*m**2*t1*lam**(-1)*BS1f + 
     &    m**2*t1*lam**(-1)*mh2*CS4f - 5./2.*m**2*t1*DS3f - 3./2.*m**2*
     &    t1*DS1f + 2*m**2*s*t1**(-2) - 1./2.*m**2*s*t1**(-2)*BS5f - 
     &    m**2*s*t1**(-2)*BS2f + 1./2.*m**2*s*t1**(-1)*lam**(-1)*mh2*
     &    BS4f - 3./2.*m**2*s*t1**(-1)*lam**(-1)*mh2*BS2f + m**2*s*
     &    t1**(-1)*lam**(-1)*mh2*BS1f - 3./2.*m**2*s*t1**(-1)*lam**(-1)
     &    *mh2**2*CS4f - 1./2.*m**2*s*t1**(-1)*CS6f + 1./2.*m**2*s*
     &    t1**(-1)*CS5f - 1./2.*m**2*s*t1**(-1)*CS4f - 1./2.*m**2*s*
     &    t1**(-1)*CS2f )
      svnonfac = svnonfac + CF*NC**(-1) * (  - 1./2.*m**2*s*t1**(-1)*
     &    CS1f + 2*m**2*s*t1**(-1)*mh2*DS3f + 7./2.*m**2*s*t1**(-1)*mh2
     &    *DS1f + 2*m**2*s*lam**(-1)*BS2f - 2*m**2*s*lam**(-1)*BS1f + 4
     &    *m**2*s*lam**(-1)*mh2*CS4f - 3./2.*m**2*s*DS3f - 5./2.*m**2*s
     &    *DS1f - m**2*s*t1*lam**(-1)*CS4f + 1./2.*m**2*s**2*t1**(-1)*
     &    lam**(-1)*BS2f - 1./2.*m**2*s**2*t1**(-1)*lam**(-1)*BS1f + 3./
     &    2.*m**2*s**2*t1**(-1)*lam**(-1)*mh2*CS4f - m**2*s**2*t1**(-1)
     &    *DS1f - 2*m**2*s**2*lam**(-1)*CS4f - 1./2.*m**2*s**3*t1**(-1)
     &    *lam**(-1)*CS4f - 1./2.*m**4*s**(-1)*t**(-1)*t1**(-1)*AS1f - 
     &    m**4*s**(-1)*t**(-1)*BS5f - m**4*s**(-1)*t1**(-2)*AS1f - 10*
     &    m**4*s**(-1)*t1**(-2)*mh2 + 4*m**4*s**(-1)*t1**(-2)*mh2*BS2f
     &     - 3./2.*m**4*s**(-1)*t1**(-1)*lam**(-1)*mh2**2*BS4f + 3./2.*
     &    m**4*s**(-1)*t1**(-1)*lam**(-1)*mh2**2*BS2f + 9*m**4*s**(-1)*
     &    t1**(-1) + 4*m**4*s**(-1)*t1**(-1)*rflag - 3*m**4*s**(-1)*
     &    t1**(-1)*lnmu - 3./2.*m**4*s**(-1)*t1**(-1)*BS4f - 5./2.*m**4
     &    *s**(-1)*t1**(-1)*BS2f )
      svnonfac = svnonfac + CF*NC**(-1) * ( 3*m**4*s**(-1)*t1**(-1)*mh2
     &    *CS6f + 3*m**4*s**(-1)*t1**(-1)*mh2*CS4f + 4*m**4*s**(-1)*
     &    lam**(-1)*mh2*BS4f - 4*m**4*s**(-1)*lam**(-1)*mh2*BS2f + m**4
     &    *s**(-1)*CS8f - m**4*s**(-1)*CS7f - 3./2.*m**4*s**(-1)*CS6f
     &     + m**4*s**(-1)*CS5f - m**4*s**(-1)*CS4f - 1./2.*m**4*s**(-1)
     &    *CS2f + 3*m**4*s**(-1)*mh2*DS3f - m**4*s**(-1)*t1*lam**(-1)*
     &    BS4f + m**4*s**(-1)*t1*lam**(-1)*BS2f - 3./2.*m**4*s**(-1)*t1
     &    *DS3f - 2*m**4*t1**(-3)*AS1f - 10*m**4*t1**(-3)*mh2 + 4*m**4*
     &    t1**(-3)*mh2*BS2f + 33./2.*m**4*t1**(-2) + 4*m**4*t1**(-2)*
     &    rflag - 3*m**4*t1**(-2)*lnmu - m**4*t1**(-2)*BS5f - 2*m**4*
     &    t1**(-2)*BS4f - 4*m**4*t1**(-2)*BS2f + 4*m**4*t1**(-2)*mh2*
     &    CS6f + m**4*t1**(-1)*lam**(-1)*mh2*BS4f - 2*m**4*t1**(-1)*
     &    lam**(-1)*mh2*BS2f + m**4*t1**(-1)*lam**(-1)*mh2*BS1f - 3./2.
     &    *m**4*t1**(-1)*lam**(-1)*mh2**2*CS4f + m**4*t1**(-1)*CS8f - 
     &    m**4*t1**(-1)*CS7f - 3*m**4*t1**(-1)*CS6f + m**4*t1**(-1)*
     &    CS5f )
      svnonfac = svnonfac + CF*NC**(-1) * ( m**4*t1**(-1)*CS4f - 1./2.*
     &    m**4*t1**(-1)*CS2f - m**4*t1**(-1)*CS1f + 4*m**4*t1**(-1)*mh2
     &    *DS3f + 3*m**4*t1**(-1)*mh2*DS1f + 2*m**4*lam**(-1)*BS4f - 2*
     &    m**4*lam**(-1)*BS1f + 4*m**4*lam**(-1)*mh2*CS4f - 4*m**4*DS3f
     &     - 2*m**4*DS1f - m**4*t1*lam**(-1)*CS4f - m**4*s*t1**(-3)*
     &    BS5f + m**4*s*t1**(-3)*BS2f - m**4*s*t1**(-2)*CS7f + 1./2.*
     &    m**4*s*t1**(-1)*lam**(-1)*BS4f + 1./2.*m**4*s*t1**(-1)*
     &    lam**(-1)*BS2f - m**4*s*t1**(-1)*lam**(-1)*BS1f + 2*m**4*s*
     &    t1**(-1)*lam**(-1)*mh2*CS4f - 2*m**4*s*t1**(-1)*DS3f - 3./2.*
     &    m**4*s*t1**(-1)*DS1f - 1./2.*m**4*s**2*t1**(-1)*lam**(-1)*
     &    CS4f + 5*m**6*s**(-1)*t1**(-2) - 2*m**6*s**(-1)*t1**(-2)*BS2f
     &     + 3./2.*m**6*s**(-1)*t1**(-1)*lam**(-1)*mh2*BS4f - 3./2.*
     &    m**6*s**(-1)*t1**(-1)*lam**(-1)*mh2*BS2f - m**6*s**(-1)*
     &    t1**(-1)*CS6f - m**6*s**(-1)*t1**(-1)*CS4f - 2*m**6*s**(-1)*
     &    lam**(-1)*BS4f + 2*m**6*s**(-1)*lam**(-1)*BS2f - m**6*s**(-1)
     &    *DS3f )
      svnonfac = svnonfac + CF*NC**(-1) * ( 10*m**6*t1**(-3) - 4*m**6*
     &    t1**(-3)*BS2f - 2*m**6*t1**(-2)*CS6f + 1./2.*m**6*t1**(-1)*
     &    lam**(-1)*BS2f - 1./2.*m**6*t1**(-1)*lam**(-1)*BS1f + 3./2.*
     &    m**6*t1**(-1)*lam**(-1)*mh2*CS4f - 2*m**6*t1**(-1)*DS3f - 
     &    m**6*t1**(-1)*DS1f - 2*m**6*lam**(-1)*CS4f - 1./2.*m**6*s*
     &    t1**(-1)*lam**(-1)*CS4f - 1./2.*m**8*s**(-1)*t1**(-1)*
     &    lam**(-1)*BS4f + 1./2.*m**8*s**(-1)*t1**(-1)*lam**(-1)*BS2f
     &     - 1./2.*m**8*t1**(-1)*lam**(-1)*CS4f + 1/(u1 + s)*t1**(-1)*
     &    mh2**2*BS5f - 1/(u1 + s)*t1**(-1)*mh2**2*BS4f - 1/(u1 + s)*
     &    mh2*BS5f + 1/(u1 + s)*mh2*BS4f + 2/(t1 + s)*mh2*BS4f - 2/(t1
     &     + s)*mh2*BS3f - 2/(t1 + s)*t1*BS4f + 2/(t1 + s)*t1*BS3f + 2
     &    /(t1 + s)*s*t1**(-1)*mh2*BS4f - 2/(t1 + s)*s*t1**(-1)*mh2*
     &    BS3f - 4/(t1 + s)*s*BS4f + 4/(t1 + s)*s*BS3f - 2/(t1 + s)*
     &    s**2*t1**(-1)*BS4f + 2/(t1 + s)*s**2*t1**(-1)*BS3f - 2/(t1 + 
     &    s)*m**2*BS4f + 2/(t1 + s)*m**2*BS3f - 2/(t1 + s)*m**2*s*
     &    t1**(-1)*BS4f )
      svnonfac = svnonfac + CF*NC**(-1) * ( 2/(t1 + s)*m**2*s*t1**(-1)*
     &    BS3f )

      return 
      end 



      function svnonfac_stef(s,t,u,m2,mh2,mu2,nf)
c     This is Stefano's result
      implicit none
      real*8 svnonfac,s,t,u,m2,mh2,mu2
      real*8 svnonfac_stef
      integer nf
      real*8 cf,ca,pi,zeta2,nc
      parameter (cf=4.d0/3.d0)
      parameter (ca=3.d0)
      parameter (pi=3.14159265358979312D0)
      parameter (zeta2=pi**2/6.d0)
      parameter (nc=3.d0)
c Set the couplings to one; physical values are assigned elsewhere
      real*8 gs,gw
      parameter (gs=1.d0)
      parameter (gw=1.d0)
      real*8 lam,u1,t1,lnmu,delq,lns,lnt,lnu,lnq,xl2tot1,xl2uou1,
     # xl2tom,xl2uom,xl2mwom,xl2omdqom,xl2omdqot1,xl2omdqou1,cs4f,
     # ddilog,as1f,BS1F,BS2F,BS3F,BS4F,BS5F,CS1F0,CS1F,CS2F0,CS2F,
     # CS3F0,CS3F,CS5F,CS6F,cs7f,cs8f,DS1F0,DS1F,DS2F0,DS2F,DS3F0,
     # DS3F,m
c
      m=sqrt(m2)
      lam=s**2.+m2**2.+mh2**2.-2.*s*m2-2.*s*mh2-2.*mh2*m2
      u1=u-m2
      t1=t-m2
      lnmu=log(mu2/m2)
      call getsymbols(s,t,u,m2,mh2,u1,t1,delq,lns,lnt,lnu,lnq,
     #                xl2tot1,xl2uou1,xl2tom,xl2uom,xl2mwom,
     #                xl2omdqom,xl2omdqot1,xl2omdqou1)
c Finite parts of scalar integrals
      cs4f=-1./sqrt(lam)*(2*ddilog(1.-sqrt(lam)/s)+ddilog((del
     &q**2.+delq*(s+sqrt(lam))-2*s*mh2)/(delq**2.+delq*(s-sqrt(lam))-2*s
     &*mh2))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sqrt(lam)))/(d
     &elq**2.+delq*(s-sqrt(lam))-2*s*mh2))-2*ddilog((s-sqrt(lam)-delq)/(
     &s+sqrt(lam)-delq))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mh2*(s+sq
     &rt(lam)))/(delq**2.+delq*(s+sqrt(lam))-2*s*mh2))-1./2.*(log((s-sqr
     &t(lam)-delq)/(s+sqrt(lam)-delq)))**2.+zeta2)
c 
      as1f = m2
      BS1F  = 2-LNS
      BS2F  = 2
      BS3F  = 2-LNU*U1/U
      BS4F  = 2-DELQ*LNQ/MH2
      BS5F  = 2-LNT*T1/T
      CS1F0 = (LNS**2-PI**2)/S/2.0D0
      CS1F = CS1F0-PI**2/S/1.2D1
      CS2F0 = (-XL2UOU1+PI**2/1.2D1+LNU**2/2.0D0)/U1
      CS2F = CS2F0-PI**2/U1/2.4D1
      CS3F0 = (-XL2TOT1+PI**2/1.2D1+LNT**2/2.0D0)/T1
      CS3F = CS3F0-PI**2/T1/2.4D1
      if(m2.lt.mh2)then
        CS5F = (-XL2UOM+XL2MWOM-PI**2-LNU**2+LNQ**2)/(DELQ-U1)
        CS6F = (-XL2TOM+XL2MWOM-PI**2-LNT**2+LNQ**2)/(DELQ-T1)
      else
        CS5F  = (-XL2UOM+XL2MWOM-LNU**2+LNQ**2)/(DELQ-U1)
        CS6F  = (-XL2TOM+XL2MWOM-LNT**2+LNQ**2)/(DELQ-T1)
      endif
      cs7f = (zeta2-xl2tom)/t1
      cs8f = (xl2mwom-xl2uom)/(u1-delq)
      if(m2.lt.mh2)then
        DS1F0 = (-2*XL2OMDQOU1+7.0D0*PI**2/1.2D1+2*LNS*LNU-LNQ**2)/
     #          (S*U1)
        DS1F = DS1F0-PI**2/(S*U1)/8.0D0
        DS2F0 = (-2*XL2OMDQOT1+7.0D0*PI**2/1.2D1+2*LNS*LNT-LNQ**2)/
     #          (S*T1)
        DS2F = DS2F0-PI**2/(S*T1)/8.0D0
        DS3F0 = (-2*XL2OMDQOU1-2*XL2OMDQOT1+1.1D1*PI**2/1.2D1+
     #           2*LNT*LNU-LNQ**2)/(T1*U1)
        DS3F = DS3F0-PI**2/(T1*U1)/2.4D1
      else
        DS1F0 = (-2*XL2OMDQOU1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNU-LNQ**2)/
     #          (S*U1)
        DS1F = DS1F0-PI**2/(S*U1)/8.0D0
        DS2F0 = (-2*XL2OMDQOT1+(-5.0D0)*PI**2/1.2D1+2*LNS*LNT-LNQ**2)/
     #          (S*T1)
        DS2F = DS2F0-PI**2/(S*T1)/8.0D0
        DS3F0 = (-2*XL2OMDQOU1-2*XL2OMDQOT1-PI**2/1.2D1+2*LNT*LNU-
     #           LNQ**2)/(T1*U1)
        DS3F = DS3F0-PI**2/(T1*U1)/2.4D1
      endif
c

      SVNONFAC = -BS5F*CA*CF*MH2*S/((T-M2)*(U+S-M2))/8.0D0+BS4F*CA*CF*MH
     1   2*S/((T-M2)*(U+S-M2))/8.0D0+BS5F*CA*CF*M**2*S/((T-M2)*(U+S-M2))
     2   /8.0D0
      SVNONFAC = -BS4F*CA*CF*M**2*S/((T-M2)*(U+S-M2))/8.0D0+BS5F*CF**2*M
     1   H2**2/((T-M2)*(U+S-M2))/4.0D0-BS4F*CF**2*MH2**2/((T-M2)*(U+S-M2
     2   ))/4.0D0+BS5F*CA*CF*MH2**2/((T-M2)*(U+S-M2))/8.0D0-BS4F*CA*CF*M
     3   H2**2/((T-M2)*(U+S-M2))/8.0D0-BS5F*CA*CF*M2*MH2/((T-M2)*(U+S-M2
     4   ))/8.0D0+BS4F*CA*CF*M2*MH2/((T-M2)*(U+S-M2))/8.0D0-BS5F*CA*CF*M
     5   **2*MH2/((T-M2)*(U+S-M2))/8.0D0+BS4F*CA*CF*M**2*MH2/((T-M2)*(U+
     6   S-M2))/8.0D0+BS5F*CA*CF*M**2*M2/((T-M2)*(U+S-M2))/8.0D0-BS4F*CA
     7   *CF*M**2*M2/((T-M2)*(U+S-M2))/8.0D0+BS5F*CA*CF*S/(U+S-M2)/8.0D0
     8   -BS4F*CA*CF*S/(U+S-M2)/8.0D0-BS5F*CF**2*MH2/(U+S-M2)/4.0D0+BS4F
     9   *CF**2*MH2/(U+S-M2)/4.0D0-BS5F*CA*CF*MH2/(U+S-M2)/8.0D0+BS4F*CA
     :   *CF*MH2/(U+S-M2)/8.0D0+BS5F*CA*CF*M2/(U+S-M2)/1.6D1-BS4F*CA*CF*
     ;   M2/(U+S-M2)/1.6D1+BS5F*CA*CF*M**2/(U+S-M2)/1.6D1-BS4F*CA*CF*M**
     <   2/(U+S-M2)/1.6D1-BS4F*CF**2*(T-M2)/(T+S-M2)/2.0D0+BS3F*CF**2*(T
     =   -M2)/(T+S-M2)/2.0D0+BS4F*CA*CF*(T-M2)/(T+S-M2)/4.0D0-BS3F*CA*CF
     >   *(T-M2)/(T+S-M2)/4.0D0-BS4F*CF**2*S**2/((T-M2)*(T+S-M2))/2.0D0+
     ?   BS3F*CF**2*S**2/((T-M2)*(T+S-M2))/2.0D0+BS4F*CA*CF*S**2/((T-M2)
     @   *(T+S-M2))/4.0D0-BS3F*CA*CF*S**2/((T-M2)*(T+S-M2))/4.0D0+BS4F*C
     1   F**2*MH2*S/((T-M2)*(T+S-M2))/2.0D0-BS3F*CF**2*MH2*S/((T-M2)*(T+
     2   S-M2))/2.0D0-BS4F*CA*CF*MH2*S/((T-M2)*(T+S-M2))/4.0D0+BS3F*CA*C
     3   F*MH2*S/((T-M2)*(T+S-M2))/4.0D0-BS4F*CF**2*M2*S/((T-M2)*(T+S-M2
     4   ))/2.0D0+BS3F*CF**2*M2*S/((T-M2)*(T+S-M2))/2.0D0+SVNONFAC
      SVNONFAC = BS4F*CA*CF*M2*S/((T-M2)*(T+S-M2))/4.0D0-BS3F*CA*CF*M2*S
     1   /((T-M2)*(T+S-M2))/4.0D0-BS4F*CF**2*S/(T+S-M2)+BS3F*CF**2*S/(T+
     2   S-M2)+BS4F*CA*CF*S/(T+S-M2)/2.0D0-BS3F*CA*CF*S/(T+S-M2)/2.0D0+B
     3   S4F*CF**2*MH2/(T+S-M2)/2.0D0-BS3F*CF**2*MH2/(T+S-M2)/2.0D0-BS4F
     4   *CA*CF*MH2/(T+S-M2)/4.0D0+BS3F*CA*CF*MH2/(T+S-M2)/4.0D0-CF**2*D
     5   S3F*(T-M2)**3/S/8.0D0+CA*CF*DS3F*(T-M2)**3/S/1.6D1-BS5F*CF**2*(
     6   T-M2)**2/(S*T)/8.0D0+3.0D0*CF**2*DS3F*MH2*(T-M2)**2/(8.0D0*S)+(
     7   -3.0D0)*CA*CF*DS3F*MH2*(T-M2)**2/(1.6D1*S)-CF**2*DS3F*M2*(T-M2)
     8   **2/S/4.0D0+CA*CF*DS3F*M2*(T-M2)**2/S/8.0D0+CF**2*CS8F*(T-M2)**
     9   2/S/8.0D0-CA*CF*CS8F*(T-M2)**2/S/1.6D1+SVNONFAC
      SVNONFAC = -CF**2*CS7F*(T-M2)**2/S/8.0D0+CA*CF*CS7F*(T-M2)**2/S/1.
     1   6D1-CF**2*CS6F*(T-M2)**2/S/8.0D0+CF**2*CS5F*(T-M2)**2/S/8.0D0-C
     2   A*CF*CS5F*(T-M2)**2/S/1.6D1+CA*CF*CS3F*(T-M2)**2/S/1.6D1+(-3.0D
     3   0)*CF**2*DS3F*(T-M2)**2/8.0D0+3.0D0*CA*CF*DS3F*(T-M2)**2/1.6D1+
     4   CA*CF*DS2F*(T-M2)**2/1.6D1-CF**2*DS1F*(T-M2)**2/8.0D0+CA*CF*DS1
     5   F*(T-M2)**2/1.6D1+BS5F*CF**2*MH2*(T-M2)/(S*T)/4.0D0+(-3.0D0)*BS
     6   5F*CF**2*M2*(T-M2)/(8.0D0*S*T)-BS5F*CF**2*(T-M2)/T/8.0D0+SVNONF
     7   AC
      SVNONFAC = -BS4F*CF**2*M2/(T+S-M2)/2.0D0-CA*CF*CS4F*S**2*(T-M2)/LA
     1   M/1.6D1+CA*CF*CS4F*MH2*S*(T-M2)/LAM/8.0D0-CF**2*CS4F*M2*S*(T-M2
     2   )/LAM/4.0D0+CA*CF*CS4F*M2*S*(T-M2)/LAM/8.0D0-BS4F*CA*CF*S*(T-M2
     3   )/LAM/1.6D1+BS2F*CA*CF*S*(T-M2)/LAM/1.6D1+(-3.0D0)*CF**2*DS3F*S
     4   *(T-M2)/8.0D0+3.0D0*CA*CF*DS3F*S*(T-M2)/1.6D1+CA*CF*DS2F*S*(T-M
     5   2)/8.0D0+(-3.0D0)*CF**2*DS1F*S*(T-M2)/8.0D0+3.0D0*CA*CF*DS1F*S*
     6   (T-M2)/1.6D1-BS4F*CA*CF*MH2**2*(T-M2)/(LAM*S)/1.6D1+BS2F*CA*CF*
     7   MH2**2*(T-M2)/(LAM*S)/1.6D1-CF**2*DS3F*MH2**2*(T-M2)/S/2.0D0+CA
     8   *CF*DS3F*MH2**2*(T-M2)/S/4.0D0+BS4F*CF**2*M2*MH2*(T-M2)/(LAM*S)
     9   /4.0D0-BS2F*CF**2*M2*MH2*(T-M2)/(LAM*S)/4.0D0+BS4F*CA*CF*M2*MH2
     :   *(T-M2)/(LAM*S)/8.0D0-BS2F*CA*CF*M2*MH2*(T-M2)/(LAM*S)/8.0D0+SV
     ;   NONFAC
      SVNONFAC = BS3F*CF**2*M2/(T+S-M2)/2.0D0+7.0D0*CF**2*DS3F*M2*MH2*(T
     1   -M2)/(8.0D0*S)+(-7.0D0)*CA*CF*DS3F*M2*MH2*(T-M2)/(1.6D1*S)-CF**
     2   2*CS8F*MH2*(T-M2)/S/4.0D0+CA*CF*CS8F*MH2*(T-M2)/S/8.0D0+CF**2*C
     3   S7F*MH2*(T-M2)/S/4.0D0-CA*CF*CS7F*MH2*(T-M2)/S/8.0D0+3.0D0*CF**
     4   2*CS6F*MH2*(T-M2)/(8.0D0*S)-CF**2*CS5F*MH2*(T-M2)/S/4.0D0+CA*CF
     5   *CS5F*MH2*(T-M2)/S/8.0D0+CF**2*CS4F*MH2*(T-M2)/S/8.0D0-CA*CF*CS
     6   3F*MH2*(T-M2)/S/8.0D0-BS4F*CF**2*M2**2*(T-M2)/(LAM*S)/4.0D0+BS2
     7   F*CF**2*M2**2*(T-M2)/(LAM*S)/4.0D0-BS4F*CA*CF*M2**2*(T-M2)/(LAM
     8   *S)/1.6D1+BS2F*CA*CF*M2**2*(T-M2)/(LAM*S)/1.6D1+(-3.0D0)*CF**2*
     9   DS3F*M2**2*(T-M2)/(8.0D0*S)+3.0D0*CA*CF*DS3F*M2**2*(T-M2)/(1.6D
     :   1*S)+CF**2*CS8F*M2*(T-M2)/S/8.0D0-CA*CF*CS8F*M2*(T-M2)/S/1.6D1+
     ;   SVNONFAC
      SVNONFAC = BS4F*CA*CF*M2/(T+S-M2)/4.0D0-CF**2*CS7F*M2*(T-M2)/S/8.0
     1   D0+CA*CF*CS7F*M2*(T-M2)/S/1.6D1-CF**2*CS6F*M2*(T-M2)/S/4.0D0-CA
     2   *CF*CS6F*M2*(T-M2)/S/1.6D1+CF**2*CS5F*M2*(T-M2)/S/4.0D0-CA*CF*C
     3   S5F*M2*(T-M2)/S/8.0D0-CF**2*CS4F*M2*(T-M2)/S/8.0D0+CA*CF*CS3F*M
     4   2*(T-M2)/S/8.0D0-CF**2*CS2F*M2*(T-M2)/S/8.0D0+CA*CF*CS2F*M2*(T-
     5   M2)/S/1.6D1+(-3.0D0)*CF**2*LNMU*(T-M2)/(8.0D0*S)+BS5F*CF**2*(T-
     6   M2)/S/8.0D0-BS2F*CF**2*(T-M2)/S/4.0D0+CF**2*(T-M2)/S/2.0D0+BS4F
     7   *CA*CF*(T-M2)/S/1.6D1-BS2F*CA*CF*(T-M2)/S/1.6D1-CA*CF*CS4F*MH2*
     8   *2*(T-M2)/LAM/1.6D1+CF**2*CS4F*M2*MH2*(T-M2)/LAM/4.0D0+CA*CF*CS
     9   4F*M2*MH2*(T-M2)/LAM/8.0D0+SVNONFAC
      SVNONFAC = -BS3F*CA*CF*M2/(T+S-M2)/4.0D0+BS4F*CA*CF*MH2*(T-M2)/LAM
     1   /8.0D0-BS2F*CA*CF*MH2*(T-M2)/LAM/8.0D0+3.0D0*CF**2*DS3F*MH2*(T-
     2   M2)/4.0D0+(-3.0D0)*CA*CF*DS3F*MH2*(T-M2)/8.0D0-CA*CF*DS2F*MH2*(
     3   T-M2)/8.0D0+3.0D0*CF**2*DS1F*MH2*(T-M2)/8.0D0+(-3.0D0)*CA*CF*DS
     4   1F*MH2*(T-M2)/1.6D1-CF**2*CS4F*M2**2*(T-M2)/LAM/4.0D0-CA*CF*CS4
     5   F*M2**2*(T-M2)/LAM/1.6D1+BS4F*CF**2*M2*(T-M2)/LAM/4.0D0+BS2F*CF
     6   **2*M2*(T-M2)/LAM/4.0D0-BS1F*CF**2*M2*(T-M2)/LAM/2.0D0+BS4F*CA*
     7   CF*M2*(T-M2)/LAM/8.0D0-BS2F*CA*CF*M2*(T-M2)/LAM/8.0D0+(-5.0D0)*
     8   CF**2*DS3F*M2*(T-M2)/8.0D0+5.0D0*CA*CF*DS3F*M2*(T-M2)/1.6D1+CA*
     9   CF*DS2F*M2*(T-M2)/8.0D0+(-3.0D0)*CF**2*DS1F*M2*(T-M2)/8.0D0+3.0
     :   D0*CA*CF*DS1F*M2*(T-M2)/1.6D1+SVNONFAC
      SVNONFAC = 3.0D0*CF**2*CS8F*(T-M2)/8.0D0+(-3.0D0)*CA*CF*CS8F*(T-M2
     1   )/1.6D1-CF**2*CS7F*(T-M2)/4.0D0+CA*CF*CS7F*(T-M2)/8.0D0-CF**2*C
     2   S6F*(T-M2)/4.0D0+3.0D0*CF**2*CS5F*(T-M2)/8.0D0+(-3.0D0)*CA*CF*C
     3   S5F*(T-M2)/1.6D1-CF**2*CS4F*(T-M2)/8.0D0+CA*CF*CS4F*(T-M2)/1.6D
     4   1+CA*CF*CS3F*(T-M2)/8.0D0-CF**2*CS1F*(T-M2)/8.0D0+CA*CF*CS1F*(T
     5   -M2)/8.0D0+AS1F*CF**2*M2*MH2/(S*T*(T-M2))/4.0D0-AS1F*CF**2*M2**
     6   2/(S*T*(T-M2))/8.0D0-AS1F*CF**2*MH2/(T*(T-M2))/8.0D0+AS1F*CF**2
     7   *M2/(T*(T-M2))/8.0D0+CA*CF*CS4F*S**4/(LAM*(T-M2))/3.2D1+(-3.0D0
     8   )*CA*CF*CS4F*MH2*S**3/(3.2D1*LAM*(T-M2))-CF**2*CS4F*M2*S**3/(LA
     9   M*(T-M2))/8.0D0-CA*CF*CS4F*M2*S**3/(LAM*(T-M2))/3.2D1+BS4F*CA*C
     :   F*S**3/(LAM*(T-M2))/3.2D1-BS2F*CA*CF*S**3/(LAM*(T-M2))/3.2D1-CF
     ;   **2*DS1F*S**3/(T-M2)/8.0D0+CA*CF*DS1F*S**3/(T-M2)/1.6D1+3.0D0*C
     <   A*CF*CS4F*MH2**2*S**2/(3.2D1*LAM*(T-M2))+3.0D0*CF**2*CS4F*M2*MH
     =   2*S**2/(8.0D0*LAM*(T-M2))-CA*CF*CS4F*M2*MH2*S**2/(LAM*(T-M2))/1
     >   .6D1+SVNONFAC
      SVNONFAC = (-3.0D0)*BS4F*CA*CF*MH2*S**2/(3.2D1*LAM*(T-M2))+3.0D0*B
     1   S2F*CA*CF*MH2*S**2/(3.2D1*LAM*(T-M2))+3.0D0*CF**2*DS1F*MH2*S**2
     2   /(8.0D0*(T-M2))+(-3.0D0)*CA*CF*DS1F*MH2*S**2/(1.6D1*(T-M2))-CF*
     3   *2*CS4F*M2**2*S**2/(LAM*(T-M2))/8.0D0-CA*CF*CS4F*M2**2*S**2/(LA
     4   M*(T-M2))/3.2D1+BS2F*CF**2*M2*S**2/(LAM*(T-M2))/8.0D0-BS1F*CF**
     5   2*M2*S**2/(LAM*(T-M2))/8.0D0-BS4F*CA*CF*M2*S**2/(LAM*(T-M2))/3.
     6   2D1+BS2F*CA*CF*M2*S**2/(LAM*(T-M2))/3.2D1-CF**2*DS1F*M2*S**2/(T
     7   -M2)/4.0D0+CA*CF*DS1F*M2*S**2/(T-M2)/8.0D0+CF**2*CS8F*S**2/(T-M
     8   2)/8.0D0-CA*CF*CS8F*S**2/(T-M2)/1.6D1+CF**2*CS5F*S**2/(T-M2)/8.
     9   0D0-CA*CF*CS5F*S**2/(T-M2)/1.6D1-CF**2*CS4F*S**2/(T-M2)/8.0D0-C
     :   A*CF*CS4F*S**2/(T-M2)/3.2D1-CF**2*CS1F*S**2/(T-M2)/8.0D0+CA*CF*
     ;   CS1F*S**2/(T-M2)/8.0D0-CA*CF*CS4F*MH2**3*S/(LAM*(T-M2))/3.2D1+(
     <   -3.0D0)*CF**2*CS4F*M2*MH2**2*S/(8.0D0*LAM*(T-M2))+3.0D0*CA*CF*C
     =   S4F*M2*MH2**2*S/(3.2D1*LAM*(T-M2))+3.0D0*BS4F*CA*CF*MH2**2*S/(3
     >   .2D1*LAM*(T-M2))+(-3.0D0)*BS2F*CA*CF*MH2**2*S/(3.2D1*LAM*(T-M2)
     ?   )-CF**2*DS1F*MH2**2*S/(T-M2)/2.0D0+CA*CF*DS1F*MH2**2*S/(T-M2)/4
     @   .0D0+CF**2*CS4F*M2**2*MH2*S/(LAM*(T-M2))/2.0D0+(-3.0D0)*CA*CF*C
     1   S4F*M2**2*MH2*S/(3.2D1*LAM*(T-M2))+BS4F*CF**2*M2*MH2*S/(LAM*(T-
     2   M2))/8.0D0+(-3.0D0)*BS2F*CF**2*M2*MH2*S/(8.0D0*LAM*(T-M2))+BS1F
     3   *CF**2*M2*MH2*S/(LAM*(T-M2))/4.0D0-BS4F*CA*CF*M2*MH2*S/(LAM*(T-
     4   M2))/1.6D1+BS2F*CA*CF*M2*MH2*S/(LAM*(T-M2))/1.6D1+CF**2*DS3F*M2
     5   *MH2*S/(T-M2)/2.0D0-CA*CF*DS3F*M2*MH2*S/(T-M2)/4.0D0-CA*CF*DS2F
     6   *M2*MH2*S/(T-M2)/4.0D0+7.0D0*CF**2*DS1F*M2*MH2*S/(8.0D0*(T-M2))
     7   +(-7.0D0)*CA*CF*DS1F*M2*MH2*S/(1.6D1*(T-M2))+SVNONFAC
      SVNONFAC = -CF**2*CS8F*MH2*S/(T-M2)/4.0D0+CA*CF*CS8F*MH2*S/(T-M2)/
     1   8.0D0+CF**2*CS6F*MH2*S/(T-M2)/8.0D0-CF**2*CS5F*MH2*S/(T-M2)/4.0
     2   D0+CA*CF*CS5F*MH2*S/(T-M2)/8.0D0+3.0D0*CF**2*CS4F*MH2*S/(8.0D0*
     3   (T-M2))+CA*CF*CS4F*MH2*S/(T-M2)/3.2D1+CF**2*CS1F*MH2*S/(T-M2)/4
     4   .0D0-CA*CF*CS1F*MH2*S/(T-M2)/4.0D0-CF**2*CS4F*M2**3*S/(LAM*(T-M
     5   2))/8.0D0+CA*CF*CS4F*M2**3*S/(LAM*(T-M2))/3.2D1+BS4F*CF**2*M2**
     6   2*S/(LAM*(T-M2))/8.0D0+BS2F*CF**2*M2**2*S/(LAM*(T-M2))/8.0D0-BS
     7   1F*CF**2*M2**2*S/(LAM*(T-M2))/4.0D0-BS4F*CA*CF*M2**2*S/(LAM*(T-
     8   M2))/3.2D1+BS2F*CA*CF*M2**2*S/(LAM*(T-M2))/3.2D1-CF**2*DS3F*M2*
     9   *2*S/(T-M2)/2.0D0+CA*CF*DS3F*M2**2*S/(T-M2)/4.0D0+CA*CF*DS2F*M2
     :   **2*S/(T-M2)/4.0D0+(-3.0D0)*CF**2*DS1F*M2**2*S/(8.0D0*(T-M2))+3
     ;   .0D0*CA*CF*DS1F*M2**2*S/(1.6D1*(T-M2))-CF**2*CS6F*M2*S/(T-M2)/8
     <   .0D0+CF**2*CS5F*M2*S/(T-M2)/8.0D0-CA*CF*CS5F*M2*S/(T-M2)/1.6D1-
     =   CF**2*CS4F*M2*S/(T-M2)/8.0D0+(-3.0D0)*CA*CF*CS4F*M2*S/(3.2D1*(T
     >   -M2))-CF**2*CS2F*M2*S/(T-M2)/8.0D0+CA*CF*CS2F*M2*S/(T-M2)/1.6D1
     ?   -CF**2*CS1F*M2*S/(T-M2)/8.0D0+3.0D0*CA*CF*CS1F*M2*S/(1.6D1*(T-M
     @   2))+(-3.0D0)*CF**2*LNMU*S/(8.0D0*(T-M2))+BS4F*CF**2*S/(T-M2)/2.
     1   0D0-BS3F*CF**2*S/(T-M2)/2.0D0-BS2F*CF**2*S/(T-M2)/4.0D0+CF**2*S
     2   /(T-M2)/2.0D0+BS5F*CA*CF*S/(T-M2)/8.0D0+(-1.3D1)*BS4F*CA*CF*S/(
     3   3.2D1*(T-M2))+BS3F*CA*CF*S/(T-M2)/4.0D0+BS2F*CA*CF*S/(T-M2)/3.2
     4   D1+SVNONFAC
      SVNONFAC = BS4F*CF**2*M2*MH2**3/(LAM*S*(T-M2))/8.0D0-BS2F*CF**2*M2
     1   *MH2**3/(LAM*S*(T-M2))/8.0D0+CF**2*CS6F*MH2**3/(S*(T-M2))/4.0D0
     2   +CF**2*CS4F*MH2**3/(S*(T-M2))/4.0D0+(-3.0D0)*BS4F*CF**2*M2**2*M
     3   H2**2/(8.0D0*LAM*S*(T-M2))+3.0D0*BS2F*CF**2*M2**2*MH2**2/(8.0D0
     4   *LAM*S*(T-M2))+(-3.0D0)*CF**2*CS6F*M2*MH2**2/(4.0D0*S*(T-M2))+(
     5   -3.0D0)*CF**2*CS4F*M2*MH2**2/(4.0D0*S*(T-M2))+(-3.0D0)*CF**2*LN
     6   MU*MH2**2/(4.0D0*S*(T-M2))-BS2F*CF**2*MH2**2/(S*(T-M2))/2.0D0+C
     7   F**2*MH2**2/(S*(T-M2))+3.0D0*BS4F*CF**2*M2**3*MH2/(8.0D0*LAM*S*
     8   (T-M2))+(-3.0D0)*BS2F*CF**2*M2**3*MH2/(8.0D0*LAM*S*(T-M2))+3.0D
     9   0*CF**2*CS6F*M2**2*MH2/(4.0D0*S*(T-M2))+3.0D0*CF**2*CS4F*M2**2*
     :   MH2/(4.0D0*S*(T-M2))+3.0D0*CF**2*LNMU*M2*MH2/(2.0D0*S*(T-M2))+3
     ;   .0D0*BS4F*CF**2*M2*MH2/(8.0D0*S*(T-M2))+9.0D0*BS2F*CF**2*M2*MH2
     <   /(8.0D0*S*(T-M2))+(-1.3D1)*CF**2*M2*MH2/(4.0D0*S*(T-M2))-BS4F*C
     =   F**2*M2**4/(LAM*S*(T-M2))/8.0D0+BS2F*CF**2*M2**4/(LAM*S*(T-M2))
     >   /8.0D0-CF**2*CS6F*M2**3/(S*(T-M2))/4.0D0-CF**2*CS4F*M2**3/(S*(T
     ?   -M2))/4.0D0+(-3.0D0)*CF**2*LNMU*M2**2/(4.0D0*S*(T-M2))+(-3.0D0)
     @   *BS4F*CF**2*M2**2/(8.0D0*S*(T-M2))+(-5.0D0)*BS2F*CF**2*M2**2/(8
     1   .0D0*S*(T-M2))+9.0D0*CF**2*M2**2/(4.0D0*S*(T-M2))-AS1F*CF**2*M2
     2   /(S*(T-M2))/8.0D0+CF**2*CS4F*M2*MH2**3/(LAM*(T-M2))/8.0D0-BS4F*
     3   CA*CF*MH2**3/(LAM*(T-M2))/3.2D1+BS2F*CA*CF*MH2**3/(LAM*(T-M2))/
     4   3.2D1+CF**2*DS1F*MH2**3/(T-M2)/4.0D0-CA*CF*DS1F*MH2**3/(T-M2)/8
     5   .0D0+(-3.0D0)*CF**2*CS4F*M2**2*MH2**2/(8.0D0*LAM*(T-M2))-BS4F*C
     6   F**2*M2*MH2**2/(LAM*(T-M2))/4.0D0+3.0D0*BS2F*CF**2*M2*MH2**2/(8
     7   .0D0*LAM*(T-M2))-BS1F*CF**2*M2*MH2**2/(LAM*(T-M2))/8.0D0+3.0D0*
     8   BS4F*CA*CF*M2*MH2**2/(3.2D1*LAM*(T-M2))+(-3.0D0)*BS2F*CA*CF*M2*
     9   MH2**2/(3.2D1*LAM*(T-M2))+SVNONFAC
      SVNONFAC = -CF**2*DS3F*M2*MH2**2/(T-M2)/2.0D0+CA*CF*DS3F*M2*MH2**2
     1   /(T-M2)/4.0D0+(-3.0D0)*CF**2*DS1F*M2*MH2**2/(4.0D0*(T-M2))+3.0D
     2   0*CA*CF*DS1F*M2*MH2**2/(8.0D0*(T-M2))+CF**2*CS8F*MH2**2/(T-M2)/
     3   4.0D0-CA*CF*CS8F*MH2**2/(T-M2)/8.0D0-CF**2*CS6F*MH2**2/(T-M2)/4
     4   .0D0+CF**2*CS5F*MH2**2/(T-M2)/4.0D0-CA*CF*CS5F*MH2**2/(T-M2)/8.
     5   0D0-CF**2*CS4F*MH2**2/(T-M2)/2.0D0-CF**2*CS1F*MH2**2/(T-M2)/4.0
     6   D0+CA*CF*CS1F*MH2**2/(T-M2)/4.0D0+3.0D0*CF**2*CS4F*M2**3*MH2/(8
     7   .0D0*LAM*(T-M2))+BS4F*CF**2*M2**2*MH2/(LAM*(T-M2))/4.0D0-BS2F*C
     8   F**2*M2**2*MH2/(LAM*(T-M2))/2.0D0+BS1F*CF**2*M2**2*MH2/(LAM*(T-
     9   M2))/4.0D0+(-3.0D0)*BS4F*CA*CF*M2**2*MH2/(3.2D1*LAM*(T-M2))+3.0
     :   D0*BS2F*CA*CF*M2**2*MH2/(3.2D1*LAM*(T-M2))+CF**2*DS3F*M2**2*MH2
     ;   /(T-M2)-CA*CF*DS3F*M2**2*MH2/(T-M2)/2.0D0+3.0D0*CF**2*DS1F*M2**
     <   2*MH2/(4.0D0*(T-M2))+(-3.0D0)*CA*CF*DS1F*M2**2*MH2/(8.0D0*(T-M2
     =   ))-CF**2*CS8F*M2*MH2/(T-M2)/2.0D0+CA*CF*CS8F*M2*MH2/(T-M2)/4.0D
     >   0+CF**2*CS7F*M2*MH2/(T-M2)/4.0D0-CA*CF*CS7F*M2*MH2/(T-M2)/8.0D0
     ?   +CF**2*CS6F*M2*MH2/(T-M2)-CF**2*CS5F*M2*MH2/(T-M2)/2.0D0+CA*CF*
     @   CS5F*M2*MH2/(T-M2)/4.0D0+CF**2*CS4F*M2*MH2/(T-M2)/4.0D0+CA*CF*C
     1   S4F*M2*MH2/(T-M2)/1.6D1+CF**2*CS2F*M2*MH2/(T-M2)/8.0D0-CA*CF*CS
     2   2F*M2*MH2/(T-M2)/1.6D1+CF**2*CS1F*M2*MH2/(T-M2)/2.0D0-CA*CF*CS1
     3   F*M2*MH2/(T-M2)/2.0D0+3.0D0*CF**2*LNMU*MH2/(4.0D0*(T-M2))-BS5F*
     4   CF**2*MH2/(T-M2)/8.0D0-BS4F*CF**2*MH2/(T-M2)/4.0D0+BS3F*CF**2*M
     5   H2/(T-M2)/2.0D0+SVNONFAC
      SVNONFAC = BS2F*CF**2*MH2/(T-M2)/2.0D0+(-7.0D0)*CF**2*MH2/(8.0D0*(
     1   T-M2))-BS5F*CA*CF*MH2/(T-M2)/8.0D0+1.3D1*BS4F*CA*CF*MH2/(3.2D1*
     2   (T-M2))-BS3F*CA*CF*MH2/(T-M2)/4.0D0-BS2F*CA*CF*MH2/(T-M2)/3.2D1
     3   -CA*CF*MH2/(T-M2)/8.0D0-CF**2*CS4F*M2**4/(LAM*(T-M2))/8.0D0+BS2
     4   F*CF**2*M2**3/(LAM*(T-M2))/8.0D0-BS1F*CF**2*M2**3/(LAM*(T-M2))/
     5   8.0D0+BS4F*CA*CF*M2**3/(LAM*(T-M2))/3.2D1-BS2F*CA*CF*M2**3/(LAM
     6   *(T-M2))/3.2D1-CF**2*DS3F*M2**3/(T-M2)/2.0D0+CA*CF*DS3F*M2**3/(
     7   T-M2)/4.0D0-CF**2*DS1F*M2**3/(T-M2)/4.0D0+CA*CF*DS1F*M2**3/(T-M
     8   2)/8.0D0+CF**2*CS8F*M2**2/(T-M2)/4.0D0-CA*CF*CS8F*M2**2/(T-M2)/
     9   8.0D0-CF**2*CS7F*M2**2/(T-M2)/4.0D0+CA*CF*CS7F*M2**2/(T-M2)/8.0
     :   D0+(-3.0D0)*CF**2*CS6F*M2**2/(4.0D0*(T-M2))+CF**2*CS5F*M2**2/(T
     ;   -M2)/4.0D0-CA*CF*CS5F*M2**2/(T-M2)/8.0D0+CF**2*CS4F*M2**2/(T-M2
     <   )/4.0D0-CA*CF*CS4F*M2**2/(T-M2)/1.6D1-CF**2*CS2F*M2**2/(T-M2)/8
     =   .0D0+CA*CF*CS2F*M2**2/(T-M2)/1.6D1-CF**2*CS1F*M2**2/(T-M2)/4.0D
     >   0+CA*CF*CS1F*M2**2/(T-M2)/4.0D0+(-3.0D0)*CF**2*LNMU*M2/(4.0D0*(
     ?   T-M2))-BS5F*CF**2*M2/(T-M2)/8.0D0+BS4F*CF**2*M2/(T-M2)/4.0D0-BS
     @   3F*CF**2*M2/(T-M2)/2.0D0+(-7.0D0)*BS2F*CF**2*M2/(8.0D0*(T-M2))+
     1   BS1F*CF**2*M2/(T-M2)/8.0D0+2*CF**2*M2/(T-M2)+BS5F*CA*CF*M2/(T-M
     2   2)/4.0D0+(-1.3D1)*BS4F*CA*CF*M2/(3.2D1*(T-M2))+BS3F*CA*CF*M2/(T
     3   -M2)/4.0D0+SVNONFAC
      SVNONFAC = (-3.0D0)*BS2F*CA*CF*M2/(3.2D1*(T-M2))+CA*CF*M2/(T-M2)/4
     1   .0D0-AS1F*CF**2/(T-M2)/4.0D0-CF**2*CS7F*M2**2*S/(T-M2)**2/4.0D0
     2   +CA*CF*CS7F*M2**2*S/(T-M2)**2/8.0D0-BS5F*CF**2*M2*S/(T-M2)**2/8
     3   .0D0-BS2F*CF**2*M2*S/(T-M2)**2/4.0D0+CF**2*M2*S/(T-M2)**2/2.0D0
     4   +BS5F*CA*CF*M2*S/(T-M2)**2/8.0D0-BS2F*CA*CF*M2*S/(T-M2)**2/8.0D
     5   0+CA*CF*M2*S/(T-M2)**2/8.0D0-BS2F*CF**2*M2*MH2**2/(S*(T-M2)**2)
     6   /2.0D0+5.0D0*CF**2*M2*MH2**2/(4.0D0*S*(T-M2)**2)-AS1F*CF**2*MH2
     7   **2/(S*(T-M2)**2)/4.0D0+BS2F*CF**2*M2**2*MH2/(S*(T-M2)**2)+(-5.
     8   0D0)*CF**2*M2**2*MH2/(2.0D0*S*(T-M2)**2)+AS1F*CF**2*M2*MH2/(S*(
     9   T-M2)**2)/2.0D0-BS2F*CF**2*M2**3/(S*(T-M2)**2)/2.0D0+5.0D0*CF**
     :   2*M2**3/(4.0D0*S*(T-M2)**2)-AS1F*CF**2*M2**2/(S*(T-M2)**2)/4.0D
     ;   0-CF**2*CS6F*M2*MH2**2/(T-M2)**2/2.0D0+CF**2*CS6F*M2**2*MH2/(T-
     <   M2)**2+3.0D0*CF**2*LNMU*M2*MH2/(4.0D0*(T-M2)**2)+BS5F*CF**2*M2*
     =   MH2/(T-M2)**2/4.0D0+BS2F*CF**2*M2*MH2/(T-M2)**2+(-2.3D1)*CF**2*
     >   M2*MH2/(8.0D0*(T-M2)**2)-BS5F*CA*CF*M2*MH2/(T-M2)**2/8.0D0+BS2F
     ?   *CA*CF*M2*MH2/(T-M2)**2/8.0D0+3.0D0*AS1F*CF**2*MH2/(8.0D0*(T-M2
     @   )**2)-CF**2*CS6F*M2**3/(T-M2)**2/2.0D0+(-3.0D0)*CF**2*LNMU*M2**
     1   2/(4.0D0*(T-M2)**2)-BS5F*CF**2*M2**2/(T-M2)**2/4.0D0-BS4F*CF**2
     2   *M2**2/(T-M2)**2/2.0D0-BS2F*CF**2*M2**2/(T-M2)**2+3.3D1*CF**2*M
     3   2**2/(8.0D0*(T-M2)**2)+BS5F*CA*CF*M2**2/(T-M2)**2/8.0D0-BS2F*CA
     4   *CF*M2**2/(T-M2)**2/8.0D0+(-5.0D0)*AS1F*CF**2*M2/(8.0D0*(T-M2)*
     5   *2)-BS5F*CF**2*M2**2*S/(T-M2)**3/4.0D0+SVNONFAC
      SVNONFAC = BS2F*CF**2*M2**2*S/(T-M2)**3/4.0D0+BS5F*CA*CF*M2**2*S/(
     1   T-M2)**3/8.0D0-BS2F*CA*CF*M2**2*S/(T-M2)**3/8.0D0+BS2F*CF**2*M2
     2   **2*MH2/(T-M2)**3+(-5.0D0)*CF**2*M2**2*MH2/(2.0D0*(T-M2)**3)+AS
     3   1F*CF**2*M2*MH2/(T-M2)**3/2.0D0-BS2F*CF**2*M2**3/(T-M2)**3+5.0D
     4   0*CF**2*M2**3/(2.0D0*(T-M2)**3)-AS1F*CF**2*M2**2/(T-M2)**3/2.0D
     5   0+BS5F*CF**2*M2*MH2/(S*T)/4.0D0+AS1F*CF**2*MH2/(S*T)/4.0D0-BS5F
     6   *CF**2*M2**2/(S*T)/4.0D0-AS1F*CF**2*M2/(S*T)/8.0D0-BS5F*CF**2*M
     7   H2/T/8.0D0-BS5F*CF**2*M2/T/8.0D0+AS1F*CF**2/T/8.0D0+SVNONFAC-CA
     8   *CF*CS4F*S**3/LAM/3.2D1+CA*CF*CS4F*MH2*S**2/LAM/8.0D0-CF**2*CS4
     9   F*M2*S**2/LAM/2.0D0+CA*CF*CS4F*M2*S**2/LAM/8.0D0-BS4F*CA*CF*S**
     :   2/LAM/3.2D1+BS2F*CA*CF*S**2/LAM/3.2D1-CF**2*DS3F*S**2/8.0D0+CA*
     ;   CF*DS3F*S**2/1.6D1+CA*CF*DS2F*S**2/1.6D1+(-3.0D0)*CF**2*DS1F*S*
     <   *2/8.0D0+3.0D0*CA*CF*DS1F*S**2/1.6D1+(-5.0D0)*CA*CF*CS4F*MH2**2
     =   *S/(3.2D1*LAM)+CF**2*CS4F*M2*MH2*S/LAM+(-3.0D0)*CA*CF*CS4F*M2*M
     >   H2*S/(1.6D1*LAM)+BS4F*CA*CF*MH2*S/LAM/8.0D0-BS2F*CA*CF*MH2*S/LA
     ?   M/8.0D0+3.0D0*CF**2*DS3F*MH2*S/8.0D0+(-3.0D0)*CA*CF*DS3F*MH2*S/
     @   1.6D1-CA*CF*DS2F*MH2*S/8.0D0+3.0D0*CF**2*DS1F*MH2*S/4.0D0+(-3.0
     1   D0)*CA*CF*DS1F*MH2*S/8.0D0+(-5.0D0)*CA*CF*CS4F*M2**2*S/(3.2D1*L
     2   AM)+BS2F*CF**2*M2*S/LAM/2.0D0
      SVNONFAC = SVNONFAC-BS1F*CF**2*M2*S/LAM/2.0D0+BS4F*CA*CF*M2*S/LAM/
     1   8.0D0-BS2F*CA*CF*M2*S/LAM/8.0D0+(-3.0D0)*CF**2*DS3F*M2*S/8.0D0+
     2   3.0D0*CA*CF*DS3F*M2*S/1.6D1+CA*CF*DS2F*M2*S/8.0D0+(-5.0D0)*CF**
     3   2*DS1F*M2*S/8.0D0+5.0D0*CA*CF*DS1F*M2*S/1.6D1+3.0D0*CF**2*CS8F*
     4   S/8.0D0+(-3.0D0)*CA*CF*CS8F*S/1.6D1-CF**2*CS7F*S/8.0D0+CA*CF*CS
     5   7F*S/1.6D1-CF**2*CS6F*S/8.0D0+3.0D0*CF**2*CS5F*S/8.0D0+(-3.0D0)
     6   *CA*CF*CS5F*S/1.6D1-CF**2*CS4F*S/4.0D0+CA*CF*CS4F*S/3.2D1+CA*CF
     7   *CS3F*S/1.6D1-CF**2*CS1F*S/4.0D0+CA*CF*CS1F*S/4.0D0+BS4F*CA*CF*
     8   MH2**3/(LAM*S)/1.6D1-BS2F*CA*CF*MH2**3/(LAM*S)/1.6D1+CF**2*DS3F
     9   *MH2**3/S/4.0D0-CA*CF*DS3F*MH2**3/S/8.0D0-BS4F*CF**2*M2*MH2**2/
     :   (LAM*S)/2.0D0+BS2F*CF**2*M2*MH2**2/(LAM*S)/2.0D0-BS4F*CA*CF*M2*
     ;   MH2**2/(LAM*S)/1.6D1+BS2F*CA*CF*M2*MH2**2/(LAM*S)/1.6D1+(-3.0D0
     <   )*CF**2*DS3F*M2*MH2**2/(4.0D0*S)+3.0D0*CA*CF*DS3F*M2*MH2**2/(8.
     =   0D0*S)+CF**2*CS8F*MH2**2/S/4.0D0-CA*CF*CS8F*MH2**2/S/8.0D0-CF**
     >   2*CS7F*MH2**2/S/4.0D0+CA*CF*CS7F*MH2**2/S/8.0D0-CF**2*CS6F*MH2*
     ?   *2/S/2.0D0+CF**2*CS5F*MH2**2/S/4.0D0-CA*CF*CS5F*MH2**2/S/8.0D0-
     @   CF**2*CS4F*MH2**2/S/4.0D0+CA*CF*CS3F*MH2**2/S/8.0D0
      SVNONFAC = SVNONFAC+BS4F*CF**2*M2**2*MH2/(LAM*S)-BS2F*CF**2*M2**2*
     1   MH2/(LAM*S)-BS4F*CA*CF*M2**2*MH2/(LAM*S)/1.6D1+BS2F*CA*CF*M2**2
     2   *MH2/(LAM*S)/1.6D1+3.0D0*CF**2*DS3F*M2**2*MH2/(4.0D0*S)+(-3.0D0
     3   )*CA*CF*DS3F*M2**2*MH2/(8.0D0*S)-CF**2*CS8F*M2*MH2/S/2.0D0+CA*C
     4   F*CS8F*M2*MH2/S/4.0D0+CF**2*CS7F*M2*MH2/S/2.0D0-CA*CF*CS7F*M2*M
     5   H2/S/4.0D0+7.0D0*CF**2*CS6F*M2*MH2/(8.0D0*S)+CA*CF*CS6F*M2*MH2/
     6   S/1.6D1-CF**2*CS5F*M2*MH2/S/2.0D0+CA*CF*CS5F*M2*MH2/S/4.0D0+CF*
     7   *2*CS4F*M2*MH2/S/2.0D0-CA*CF*CS3F*M2*MH2/S/4.0D0+CF**2*CS2F*M2*
     8   MH2/S/8.0D0-CA*CF*CS2F*M2*MH2/S/1.6D1+3.0D0*CF**2*LNMU*MH2/(4.0
     9   D0*S)-BS5F*CF**2*MH2/S/4.0D0+BS2F*CF**2*MH2/S/2.0D0+(-7.0D0)*CF
     :   **2*MH2/(8.0D0*S)-BS4F*CA*CF*MH2/S/1.6D1+BS2F*CA*CF*MH2/S/1.6D1
     ;   -CA*CF*MH2/S/8.0D0-BS4F*CF**2*M2**3/(LAM*S)/2.0D0+BS2F*CF**2*M2
     <   **3/(LAM*S)/2.0D0+BS4F*CA*CF*M2**3/(LAM*S)/1.6D1-BS2F*CA*CF*M2*
     =   *3/(LAM*S)/1.6D1-CF**2*DS3F*M2**3/S/4.0D0+CA*CF*DS3F*M2**3/S/8.
     >   0D0+CF**2*CS8F*M2**2/S/4.0D0-CA*CF*CS8F*M2**2/S/8.0D0-CF**2*CS7
     ?   F*M2**2/S/4.0D0+CA*CF*CS7F*M2**2/S/8.0D0+(-3.0D0)*CF**2*CS6F*M2
     @   **2/(8.0D0*S)-CA*CF*CS6F*M2**2/S/1.6D1+CF**2*CS5F*M2**2/S/4.0D0
     1   -CA*CF*CS5F*M2**2/S/8.0D0
      SVNONFAC = SVNONFAC-CF**2*CS4F*M2**2/S/4.0D0+CA*CF*CS3F*M2**2/S/8.
     1   0D0-CF**2*CS2F*M2**2/S/8.0D0+CA*CF*CS2F*M2**2/S/1.6D1+(-3.0D0)*
     2   CF**2*LNMU*M2/(4.0D0*S)+BS5F*CF**2*M2/S/4.0D0-BS2F*CF**2*M2/S/2
     3   .0D0+7.0D0*CF**2*M2/(8.0D0*S)-BS5F*CA*CF*M2/S/1.6D1+BS2F*CA*CF*
     4   M2/S/1.6D1+CA*CF*M2/S/8.0D0+BS5F*CA*CF*M**2/S/1.6D1-BS4F*CA*CF*
     5   M**2/S/1.6D1+CA*CF*CS4F*MH2**3/LAM/1.6D1-CF**2*CS4F*M2*MH2**2/L
     6   AM/2.0D0-CA*CF*CS4F*M2*MH2**2/LAM/1.6D1+(-5.0D0)*BS4F*CA*CF*MH2
     7   **2/(3.2D1*LAM)+5.0D0*BS2F*CA*CF*MH2**2/(3.2D1*LAM)-CF**2*DS3F*
     8   MH2**2/2.0D0+CA*CF*DS3F*MH2**2/4.0D0+CA*CF*DS2F*MH2**2/8.0D0-CF
     9   **2*DS1F*MH2**2/2.0D0+CA*CF*DS1F*MH2**2/4.0D0+CF**2*CS4F*M2**2*
     :   MH2/LAM-CA*CF*CS4F*M2**2*MH2/LAM/1.6D1+BS4F*CF**2*M2*MH2/LAM/2.
     ;   0D0-BS2F*CF**2*M2*MH2/LAM+BS1F*CF**2*M2*MH2/LAM/2.0D0+(-3.0D0)*
     <   BS4F*CA*CF*M2*MH2/(1.6D1*LAM)+3.0D0*BS2F*CA*CF*M2*MH2/(1.6D1*LA
     =   M)+3.0D0*CF**2*DS3F*M2*MH2/2.0D0+(-3.0D0)*CA*CF*DS3F*M2*MH2/4.0
     >   D0-CA*CF*DS2F*M2*MH2/4.0D0+CF**2*DS1F*M2*MH2-CA*CF*DS1F*M2*MH2/
     ?   2.0D0-CF**2*CS8F*MH2/2.0D0+CA*CF*CS8F*MH2/4.0D0+CF**2*CS7F*MH2/
     @   4.0D0-CA*CF*CS7F*MH2/8.0D0
      SVNONFAC = SVNONFAC+CF**2*CS6F*MH2/2.0D0-CF**2*CS5F*MH2/2.0D0+CA*C
     1   F*CS5F*MH2/4.0D0+CF**2*CS4F*MH2/2.0D0-CA*CF*CS4F*MH2/1.6D1-CA*C
     2   F*CS3F*MH2/8.0D0+CF**2*CS1F*MH2/4.0D0-CA*CF*CS1F*MH2/4.0D0-CF**
     3   2*CS4F*M2**3/LAM/2.0D0+CA*CF*CS4F*M2**3/LAM/1.6D1+BS4F*CF**2*M2
     4   **2/LAM/2.0D0-BS1F*CF**2*M2**2/LAM/2.0D0+(-5.0D0)*BS4F*CA*CF*M2
     5   **2/(3.2D1*LAM)+5.0D0*BS2F*CA*CF*M2**2/(3.2D1*LAM)-CF**2*DS3F*M
     6   2**2+CA*CF*DS3F*M2**2/2.0D0+CA*CF*DS2F*M2**2/8.0D0-CF**2*DS1F*M
     7   2**2/2.0D0+CA*CF*DS1F*M2**2/4.0D0+CF**2*CS8F*M2/8.0D0-CA*CF*CS8
     8   F*M2/1.6D1-CF**2*CS6F*M2/2.0D0+3.0D0*CF**2*CS5F*M2/8.0D0+(-3.0D
     9   0)*CA*CF*CS5F*M2/1.6D1+(-3.0D0)*CA*CF*CS4F*M2/1.6D1+CA*CF*CS3F*
     :   M2/8.0D0-CF**2*CS2F*M2/4.0D0+CA*CF*CS2F*M2/8.0D0-CF**2*CS1F*M2/
     ;   4.0D0+CA*CF*CS1F*M2/4.0D0+(-3.0D0)*CF**2*LNMU/4.0D0+BS5F*CF**2/
     <   8.0D0+BS4F*CF**2/2.0D0-BS3F*CF**2/2.0D0-BS2F*CF**2/2.0D0+CF**2+
     =   (-7.0D0)*BS4F*CA*CF/3.2D1+BS3F*CA*CF/4.0D0-BS2F*CA*CF/3.2D1
      SVNONFAC = 4*SVNONFAC/(CF*NC)

      svnonfac_stef=svnonfac
      return 
      end 


      subroutine getsymbols(s,t,u,xm12,xm22,u1,t1,delq,lns,lnt,lnu,lnq,
     #                      xl2tot1,xl2uou1,xl2tom,xl2uom,xl2mwom,
     #                      xl2omdqom,xl2omdqot1,xl2omdqou1)
      implicit none
      real*8 s,t,u,xm12,xm22,u1,t1,delq,lns,lnt,lnu,lnq,xl2tot1,
     # xl2uou1,xl2tom,xl2uom,xl2mwom,xl2omdqom,xl2omdqot1,xl2omdqou1,
     # ddilog
c
      t1=t-xm12
      u1=u-xm12
      delq=xm22-xm12
      lns=log(s/xm12)
      lnt=log(-t1/xm12)
      lnu=log(-u1/xm12)
      if(delq.gt.0.d0)then
        lnq=log(delq/xm12)
      else
        lnq=log(-delq/xm12)
      endif
      xl2tot1=ddilog(t/t1)
      xl2uou1=ddilog(u/u1)
      xl2tom=ddilog(t/xm12)
      xl2uom=ddilog(u/xm12)
      xl2mwom=ddilog(xm22/xm12)
      xl2omdqom=ddilog(1-delq/xm12)
      xl2omdqot1=ddilog(1-delq/t1)
      xl2omdqou1=ddilog(1-delq/u1)
c
      return
      end


      subroutine getmudep(xlgmuom,SDP,SSP,SF0,xSDP,xSSP,xSF0)
      implicit none
      real*8 xlgmuom,SDP,SSP,SF0,xSDP,xSSP,xSF0
c
      xSDP = SDP
      xSSP = SSP + SDP*xlgmuom
      xSF0 = SF0 + SSP*xlgmuom + SDP*xlgmuom**2/2.d0
c
      return
      end



c
c
c Utility function for azimuthal correlations
c
c
      function qin_kern(x,index)
c This function returns the quantity (1-x)*Q_{a*b}(x), where
c Q_{a*b} are the kernels given in FKS (B.42)-(B.45), and the splitting 
c partons {ab} are defined with the following conventions
c
c         index          ab
c
c           1            gg
c           2            qg
c           3            gq
c           4            qq
c
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        qin_kern=-4*vca*(1-x)**2/x
      elseif(index.eq.2)then
        qin_kern=0.d0
      elseif(index.eq.3)then
        qin_kern=-4*vcf*(1-x)**2/x
      elseif(index.eq.4)then
        qin_kern=0.d0
      else
        write(6,*)'Error in qin_kern: wrong index value'
        stop
      endif
      return
      end
c
c
c From the jet package, Altarelli-Parisi kernels and change of scheme
c
c
      function ap_kern(x,index)
c This function returns the quantity (1-x)*P_{ab}(x), where
c P_{ab} are the Altarelli-Parisi kernels, and the splitting partons
c {ab} are defined with the following conventions
c
c         index          ab
c
c           1            gg
c           2            qg
c           3            gq
c           4            qq
c
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        ap_kern=2*vca*(x+(1-x)**2/x+x*(1-x)**2)
      elseif(index.eq.2)then
        ap_kern=vtf*(1-x)*(x**2+(1-x)**2)
      elseif(index.eq.3)then
        ap_kern=vcf*(1-x)*(1+(1-x)**2)/x
      elseif(index.eq.4)then
        ap_kern=vcf*(1+x**2)
      else
        write(6,*)'Error in ap_kern: wrong index value'
        stop
      endif
      return
      end


      function apprime_kern(x,index)
c This function returns the quantity (1-x)*P_{ab}^{prime}(x), where
c P_{ab}^{prime} is the ep-dependent part of the Altarelli-Parisi kernels, 
c and the codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        apprime_kern=0.d0
      elseif(index.eq.2)then
        apprime_kern=-2*vtf*x*(1-x)**2
      elseif(index.eq.3)then
        apprime_kern=-vcf*(1-x)*x
      elseif(index.eq.4)then
        apprime_kern=-vcf*(1-x)**2
      else
        write(6,*)'Error in apprime_kern: wrong index value'
        stop
      endif
      return
      end


      function xkdelta(index)
c This function returns the quantity K^{(d)}_{ab}, relevant for
c the MS --> DIS change in the factorization scheme. 
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkdelta=0.d0
      elseif(index.eq.2)then
        xkdelta=0.d0
      elseif(index.eq.3)then
        xkdelta=vcf*(9.d0/2.d0+pi**2/3.d0)
      elseif(index.eq.4)then
        xkdelta=-vcf*(9.d0/2.d0+pi**2/3.d0)
      else
        write(6,*)'Error in xkdelta: wrong index value'
        stop
      endif
      return
      end


      function xkplus(x,index)
c This function returns the quantity K^{(+)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkplus=0.d0
      elseif(index.eq.2)then
        xkplus=0.d0
      elseif(index.eq.3)then
        xkplus=-vcf*(-3.d0/2.d0-(1+x**2)*log(x)+(1-x)*(3+2*x))
      elseif(index.eq.4)then
        xkplus=vcf*(-3.d0/2.d0-(1+x**2)*log(x)+(1-x)*(3+2*x))
      else
        write(6,*)'Error in xkplus: wrong index value'
        stop
      endif
      return
      end


      function xklog(x,index)
c This function returns the quantity K^{(l)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xklog=0.d0
      elseif(index.eq.2)then
        xklog=0.d0
      elseif(index.eq.3)then
        xklog=-vcf*(1+x**2)
      elseif(index.eq.4)then
        xklog=vcf*(1+x**2)
      else
        write(6,*)'Error in xklog: wrong index value'
        stop
      endif
      return
      end


      function xkreg(x,index)
c This function returns the quantity K^{(reg)}_{ab}(x), relevant for
c the MS --> DIS change in the factorization scheme. Notice that
c there's NO multiplicative (1-x) factor like in the previous functions.
c The codes for the splitting partons {ab} are defined above
      implicit real * 8 (a-h,o-z)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
      parameter (xnc=3.d0)
      common/nl/nl
c
      if(index.eq.1)then
        xkreg=-2*nl*vtf*( (x**2+(1-x)**2)*log((1-x)/x)+8*x*(1-x)-1 )
      elseif(index.eq.2)then
        xkreg=vtf*( (x**2+(1-x)**2)*log((1-x)/x)+8*x*(1-x)-1 )
      elseif(index.eq.3)then
        xkreg=0.d0
      elseif(index.eq.4)then
        xkreg=0.d0
      else
        write(6,*)'Error in xkreg: wrong index value'
        stop
      endif
      return
      end
c
c
c MadEvent routines
c
c
      function xmadevht(iborn,jproc,idr,s,tk,uk,xmom)
c Wrapper for MadEvent functions. Inputs are 
c   iborn = 0(born), 1(real)
c   jproc = 1(gg), 2(qq), 3(qg)
c   idr   = 1..7 (depends on jproc)
c   s     = parton cm energy squared
c   tk,uk = FNR invariants
c   xmom  = 4-momenta obtained from invar
c Output is the matrix element squared in GeV^-2, times the flux factor,
c times 4*tk*uk/s**2 in the case of real matrix elements.
c MadEvent routines use the parameters as defined in setmeFRpar();
c for kinematic-dependent couplings, and if simple factorized
c expressions can't be found (which appears to happen only if a Higgs 
c is involved in the reaction), the call to setmeFRpar() must be included
c in this routine. 
c In the present case, all the Born (real) formulae are proportional 
c to g^2 gw^4 (A^2+B^2) [g^4 gw^8 (A^2+B^2)] (neglecting a CKM factor, which
c is trivial for unweighting purposes), and this call can be placed somewhere
c else (and done only once). The main code assumes that g^n gw^4 (A^2+B^2) 
c be stripped off matrix element routines, which is equivalent to setting
c g, gw, and (A^2+B^2) equal to one. Since setmeFRpar() is called, 
c for consistency with other processes, with g=1 and e=1, gw=1 implies 
c that the results of MadEvent must be multiplied by sthw2^2=e^4/gw^4.
c Furthermore, the call to setmeFRpar() does not set (A^2+B^2)=1, and thus
c we divide here by (A^2+B^2).
c NOTE: in the case of A and B depending upon flavour, in a way non 
c factorizable in the CKM matrix, their values must be set prior to each
c call to this routine, after choosing the flavour of the down-type quark
c in the tHD vertex
      implicit none
      integer iborn,jproc,idr
      real*8 xmadevht,s,tk,uk,xmom(10,4),xfact,tmp(1)
      real*8 pme0dc(0:3,7),pme1dc(0:3,8)
      real*8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
c This is a^2+b^2, and is set by setmeFRpar()
      real*8 apbsq
      common/capbsq/apbsq
      integer ipart,icomp
c Components: MC@NLO conventions   -> 1=px, 2=py, 3=pz, 4=E
c             MadEvent conventions -> 0=E, 1=px, 2=py, 3=pz
      integer mapcomp(0:3)
      data mapcomp/4,1,2,3/
c Labelling conventions for subprocess: 
c MC@NLO   -> x(1)y(2)  ->  z(3)t(4)[->l+(6)nu(7)b(8)]H-(5)
c MadEvent -> b(1)g(2)  ->  l+(3)nu(4)b(5) H-(6) 
c MadEvent -> g(1)g(2)  ->  l+(3)nu(4)b(5) H-(6) bx(7)
c MadEvent -> u(1)ux(2) ->  l+(3)nu(4)b(5) H-(6) bx(7)
c MadEvent -> b(1)u(2)  ->  l+(3)nu(4)b(5) H-(6) u(7)
c MadEvent -> b(1)bx(2) ->  l+(3)nu(4)b(5) H-(6) bx(7)
c MadEvent -> b(1)b(2)  ->  l+(3)nu(4)b(5) H-(6) b(7)
c MadEvent -> b(1)g(2)  ->  l+(3)nu(4)b(5) H-(6) g(7)
c In the case of bb initial state, the two final state b quarks are treated
c as different particles, by excluding the diagrams which can be otained
c from other diagrams with the formal exchange 5<->7
      integer mapdir(7),mapref(7)
c mapping for direct events
      data mapdir/1,2,6,7,8,5,3/
c mapping for reflected events
      data mapref/2,1,6,7,8,5,3/
c
      if(iborn.eq.0)then
        if(jproc.ne.3)then
          write(*,*)'Error #1 in xmadevht: iborn, jproc=',iborn,jproc
          stop
        endif
        xfact=1.d0
        do ipart=1,6
          do icomp=0,3
            if(idr.eq.1)then
              pme0dc(icomp,ipart)=xmom(mapdir(ipart),mapcomp(icomp))
            elseif(idr.eq.3)then
              pme0dc(icomp,ipart)=xmom(mapref(ipart),mapcomp(icomp))
            else
              write(*,*)'Error #2 in xmadevht: idr=',idr
              stop
            endif
          enddo
        enddo
      elseif(iborn.eq.1)then
        if( (jproc.eq.1.and.idr.ne.1) .or.
     #      (jproc.eq.3.and.(idr.ne.1.and.idr.ne.3)) )then
          write(*,*)'Error #3 in xmadevht: jproc,idr=',jproc,idr
          stop
        endif
        xfact=4*tk*uk/s**2
        do ipart=1,7
          do icomp=0,3
            if( jproc.eq.1 .or.
     #          (jproc.eq.2.and.(idr.eq.1.or.idr.eq.3.or.
     #                          idr.eq.5.or.idr.eq.7)) .or.
     #          (jproc.eq.3.and.idr.eq.1) )then
              pme1dc(icomp,ipart)=xmom(mapdir(ipart),mapcomp(icomp))
            else
              pme1dc(icomp,ipart)=xmom(mapref(ipart),mapcomp(icomp))
            endif
          enddo
        enddo
      else
        write(*,*)'xmadevht: unknown iborn value',iborn
        stop
      endif
      if(iborn.eq.0)then
        call bgth_dc(pme0dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.1)then
        call ggthb_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.1.or.idr.eq.2))then
        call uuxthbx_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.3.or.idr.eq.4))then
        call buthu_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.5.or.idr.eq.6))then
        call bbxthbx_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2.and.idr.eq.7)then
        call bbthb_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.3)then
        call bgthg_dc(pme1dc,tmp)
      else
        write(*,*)'xmadevht: Error #5'
        stop
      endif
c Insert sin(theta_W) for e^4 -> gw^4 (due to decay), and divide by A^2+B^2
c (due to production)
      xmadevht=sthw2**2/apbsq*xfact*tmp(1)/(2*s)
      return
      end


      subroutine setmeFRpar(xiwmass,xiwwidth,xizmass,xizwidth,
     #                      xitmass,xitwidth,xihmass,xihwidth,
     #                      xibmass,xisin2w,xiee2,xig,xia,xib)
c This routine has been deeply changed wrt setmeFRpar to fit the MG/ME4 
c (with FeynRules models) structure.
c The basic idea is the same, i.e. couplings are defined here.

      implicit none

c These are input parameters from the MC@NLO code.
      real * 8 xiwmass,xiwwidth,xihmass,xihwidth,xizmass,xizwidth,
     # xitmass,xitwidth,xibmass,xisin2w,xiee2,xig,xia,xib
c
      real*8 apbsq
      common/capbsq/apbsq

c These are variable definitions from MG
      include 'MEFRinput.inc'
      include 'MEFRcoupl.inc'


c This is where the matching between MC@NLO and MG variables is done

      mhc= xihmass
      whc= xihwidth
      mh = 0d0
      mz = xizmass
      wz = xizwidth
      mt = xitmass
      wt = xitwidth
      mb = xibmass
      mw = xiwmass
      ww = xiwwidth
      
      a = xia
      b = xib
      apbsq = a**2+b**2
     
      cabi = 0d0
      aewm1 = 3.5449077018110318**2/xiee2
      as = xig**2/3.5449077018110318**2

c Begin of file intparam_definition.inc. Get rid of unwanted quantities,
c and define others to give priority to inputs to this routine

c This file has been generated automatically by FeynRules

c Version: 1.3.22   Date: 14. 11. 2008,    11:25

c Internal parameters definition :
      aEW = 1/aEWM1
c Original definition was sw2 = 1 - MW**2/MZ**2. We want to keep a
c complete freedom in the choice of scheme
      sw2 = xisin2w
      ee = 3.5449077018110318*sqrt(aEW)
      cw = sqrt(1. - 1.*sw2)
      sw = sqrt(sw2)
      gw = ee/sw
      g1 = ee/cw
      G = 3.5449077018110318*sqrt(aS)
      v = (2*MW*sw)/ee
      lam = MH**2/(2.*v**2)
      muH = sqrt(lam*v**2)
      ytau = (sqrt(2.)*ymtau)/v
      yc = (sqrt(2.)*ymc)/v
      yt = (sqrt(2.)*ymt)/v
      yb = (sqrt(2.)*ymb)/v
      CKM11 = Cos(cabi)
      CKM12 = Sin(cabi)
      CKM21 = -Sin(cabi)
      CKM22 = Cos(cabi)
      Sqrt2 = sqrt(2.)
      SqrtPi = sqrt(3.141592653589793)
      Psw2 = sw**2
      Pgw2 = gw**2
      Pcw2 = cw**2
      CONJCKM21 = conjg(CKM21)
      CONJCKM11 = conjg(CKM11)
      CONJCKM22 = conjg(CKM22)
      CONJCKM12 = conjg(CKM12)


c Definition of the EW coupling used in the write out of aqed
      gal(1) = ee
      gal(2) = ee


c Definition of DUM symbols
      DUM0 = 0
      DUM1 = 1

c End of file intparam_definition.inc

c Couplings coming from FR Lagrangian 1		
      mgvx11 = -(gw*sw)
      mgvx3 = g
      mgvx15 = cw*gw
      mgvx14 = pgw2*psw2
      mgvx21 = cw*pgw2*sw
      ggt1 = sqrt(g**2)
      mgvx16 = pgw2
      mgvx24 = pcw2*pgw2

c Couplings coming from FR Lagrangian 2
      mgvx2 = -6*lam*v
      mgvx1 = -6*lam
      mgvx13 = (ee**2*v)/(2.*psw2)
      mgvx23 = ee**2*v + (ee**2*pcw2*v)/(2.*psw2) + (ee**2*psw2*v)
     &/(2.*pcw2)
      mgvx12 = ee**2/(2.*psw2)
      mgvx22 = ee**2 + (ee**2*pcw2)/(2.*psw2) + (ee**2*psw2)/(2.*pcw2)
      mgvx18(1) = -(ytau/sqrt2)
      mgvx18(2) = -(ytau/sqrt2)
      mgvx28(1) = ee
      mgvx28(2) = ee
      mgvx34(1) = -(ee/(sqrt2*sw))
      mgvx34(2) = 0
      mgvx53(1) = (cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      mgvx53(2) = -((ee*sw)/cw)
      mgvx59(1) = -(cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      mgvx59(2) = 0
      mgvx25(1) = ee/3.
      mgvx25(2) = ee/3.
      mgvx31(1) = (-2*ee)/3.
      mgvx31(2) = (-2*ee)/3.
      mgvx17(1) = -(yb/sqrt2)
      mgvx17(2) = -(yb/sqrt2)
      mgvx19(1) = -(yc/sqrt2)
      mgvx19(2) = -(yc/sqrt2)
      mgvx20(1) = -(yt/sqrt2)
      mgvx20(2) = -(yt/sqrt2)
      mgvx37(1) = -((ckm21*ee)/(sqrt2*sw))
      mgvx37(2) = 0
      mgvx38(1) = -((ckm22*ee)/(sqrt2*sw))
      mgvx38(2) = 0
      mgvx40(1) = -((ckm11*ee)/(sqrt2*sw))
      mgvx40(2) = 0
      mgvx41(1) = -((ckm12*ee)/(sqrt2*sw))
      mgvx41(2) = 0
      mgvx46(1) = -((conjckm21*ee)/(sqrt2*sw))
      mgvx46(2) = 0
      mgvx47(1) = -((conjckm11*ee)/(sqrt2*sw))
      mgvx47(2) = 0
      mgvx48(1) = -((conjckm22*ee)/(sqrt2*sw))
      mgvx48(2) = 0
      mgvx49(1) = -((conjckm12*ee)/(sqrt2*sw))
      mgvx49(2) = 0
      mgvx50(1) = (cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      mgvx50(2) = -(ee*sw)/(3.*cw)
      mgvx56(1) = -(cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      mgvx56(2) = (2*ee*sw)/(3.*cw)
      gg(1) = -g
      gg(2) = -g
      hmcoup(2) = -(a+b)
      hmcoup(1) = -(a-b)
      hpcoup(2) = -(a+b)
      hpcoup(1) = -(a-b)

      return
      end


      subroutine switchmom(p1,p,ic,jc,nexternal)
c**************************************************************************
c     Changes stuff for crossings
c**************************************************************************
      implicit none
      integer nexternal
      integer jc(nexternal),ic(nexternal)
      real*8 p1(0:3,nexternal),p(0:3,nexternal)
      integer i,j
c-----
c Begin Code
c-----
      do i=1,nexternal
         do j=0,3
            p(j,ic(i))=p1(j,i)
         enddo
      enddo
      do i=1,nexternal
         jc(i)=1
      enddo
      jc(ic(1))=-1
      jc(ic(2))=-1
      end


C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BGTH_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b h-  
C  
C Crossing   1 is b g -> e+ ve b h-  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  6)
      integer    nincoming
      parameter (nincoming=  2)
      
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATBGTH_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   3),IHEL=1, 6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   5),IHEL=1, 6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,   7),IHEL=1, 6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,   9),IHEL=1, 6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  11),IHEL=1, 6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  13),IHEL=1, 6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  15),IHEL=1, 6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  17),IHEL=1, 6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  19),IHEL=1, 6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  21),IHEL=1, 6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  23),IHEL=1, 6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  25),IHEL=1, 6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  27),IHEL=1, 6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  29),IHEL=1, 6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  31),IHEL=1, 6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 6) / 1, 1, 1, 1, 1,-1/
      DATA (  IC(IHEL,  1),IHEL=1, 6) / 1, 2, 3, 4, 5, 6/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBGTH_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBGTH_DC
      REAL*8 FUNCTION MATBGTH_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b h-  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  6)
      integer    nincoming
      parameter (nincoming=  2)
      
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[ 5, 1, 2]                                                
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),MB ,NHEL(1   ),+1*IC(1   ),W(1,1   ))          
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     7   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,7   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     8   ))                                                          
      CALL FVOXXX(W(1,8   ),W(1,2   ),GG ,MT      ,WT      ,W(1,9   ))     
      CALL IOSXXX(W(1,1   ),W(1,9   ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,MB      ,ZERO    ,W(1,10  ))     
      CALL IOSXXX(W(1,10  ),W(1,8   ),W(1,6   ),hmcoup ,AMP(2   ))         
      JAMP(   1) = +AMP(   1)+AMP(   2)
      MATBGTH_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBGTH_DC =MATBGTH_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE GGTHB_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> e+ ve b h- b~  
C  
C Crossing   1 is g g -> e+ ve b h- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)

      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATGGTHB_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATGGTHB_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATGGTHB_DC
      REAL*8 FUNCTION MATGGTHB_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> e+ ve b h- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)

      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  20, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 5, 7, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 5, 7, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL IXXXXX(P(0,7   ),MB ,NHEL(7   ),-1*IC(7   ),W(1,7   ))          
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     8   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     9   ))                                                          
      CALL FVIXXX(W(1,7   ),W(1,1   ),GG ,MB      ,ZERO    ,W(1,10  ))     
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,10  ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVOXXX(W(1,11  ),W(1,1   ),GG ,MT      ,WT      ,W(1,12  ))     
      CALL IOSXXX(W(1,7   ),W(1,12  ),W(1,6   ),hmcoup ,AMP(2   ))         
      CALL FVIXXX(W(1,7   ),W(1,2   ),GG ,MB      ,ZERO    ,W(1,13  ))     
      CALL FVOXXX(W(1,9   ),W(1,1   ),GG ,MT      ,WT      ,W(1,14  ))     
      CALL IOSXXX(W(1,13  ),W(1,14  ),W(1,6   ),hmcoup ,AMP(3   ))         
      CALL FVOXXX(W(1,14  ),W(1,2   ),GG ,MT      ,WT      ,W(1,15  ))     
      CALL IOSXXX(W(1,7   ),W(1,15  ),W(1,6   ),hmcoup ,AMP(4   ))         
      CALL FVIXXX(W(1,10  ),W(1,2   ),GG ,MB      ,ZERO    ,W(1,16  ))     
      CALL IOSXXX(W(1,16  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(5   ))         
      CALL FVIXXX(W(1,13  ),W(1,1   ),GG ,MB      ,ZERO    ,W(1,17  ))     
      CALL IOSXXX(W(1,17  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(6   ))         
      CALL JVVXXX(W(1,1   ),W(1,2   ),MGVX3 ,ZERO    ,ZERO    ,W(1,        
     &     18  ))                                                          
      CALL FVOXXX(W(1,9   ),W(1,18  ),GG ,MT      ,WT      ,W(1,19  ))     
      CALL IOSXXX(W(1,7   ),W(1,19  ),W(1,6   ),hmcoup ,AMP(7   ))         
      CALL FVIXXX(W(1,7   ),W(1,18  ),GG ,MB      ,ZERO    ,W(1,20  ))     
      CALL IOSXXX(W(1,20  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(8   ))         
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels may not be the same as in the case of production without decay

c Chris
      if(MHC.le.MT) then
         AMP(   2)=0
         AMP(   4)=0
         AMP(   7)=0
      endif

      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      MATGGTHB_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATGGTHB_DC =MATGGTHB_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BGTHG_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b h- g  
C  
C Crossing   1 is b g -> e+ ve b h- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATBGTHG_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBGTHG_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBGTHG_DC
      REAL*8 FUNCTION MATBGTHG_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b h- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  10,NEIGEN=  2) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  23, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 5, 1, 2, 7]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 5, 1, 7, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),MB ,NHEL(1   ),+1*IC(1   ),W(1,1   ))          
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL VXXXXX(P(0,7   ),ZERO ,NHEL(7   ),+1*IC(7   ),W(1,7   ))        
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     8   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     9   ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,7   ),GG ,MB      ,ZERO    ,W(1,10  ))     
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,10  ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL JVVXXX(W(1,7   ),W(1,2   ),MGVX3 ,ZERO    ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVOXXX(W(1,9   ),W(1,12  ),GG ,MT      ,WT      ,W(1,13  ))     
      CALL IOSXXX(W(1,1   ),W(1,13  ),W(1,6   ),hmcoup ,AMP(2   ))         
      CALL FVOXXX(W(1,5   ),W(1,7   ),GG ,MB      ,ZERO    ,W(1,14  ))     
      CALL FVOXXX(W(1,14  ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     15  ))                                                          
      CALL FVOXXX(W(1,15  ),W(1,2   ),GG ,MT      ,WT      ,W(1,16  ))     
      CALL IOSXXX(W(1,1   ),W(1,16  ),W(1,6   ),hmcoup ,AMP(3   ))         
      CALL FVOXXX(W(1,9   ),W(1,7   ),GG ,MT      ,WT      ,W(1,17  ))     
      CALL FVOXXX(W(1,17  ),W(1,2   ),GG ,MT      ,WT      ,W(1,18  ))     
      CALL IOSXXX(W(1,1   ),W(1,18  ),W(1,6   ),hmcoup ,AMP(4   ))         
      CALL FVOXXX(W(1,11  ),W(1,7   ),GG ,MT      ,WT      ,W(1,19  ))     
      CALL IOSXXX(W(1,1   ),W(1,19  ),W(1,6   ),hmcoup ,AMP(5   ))         
      CALL FVIXXX(W(1,10  ),W(1,2   ),GG ,MB      ,ZERO    ,W(1,20  ))     
      CALL IOSXXX(W(1,20  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(6   ))         
      CALL FVIXXX(W(1,1   ),W(1,12  ),GG ,MB      ,ZERO    ,W(1,21  ))     
      CALL IOSXXX(W(1,21  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(7   ))         
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,MB      ,ZERO    ,W(1,22  ))     
      CALL IOSXXX(W(1,22  ),W(1,15  ),W(1,6   ),hmcoup ,AMP(8   ))         
      CALL FVIXXX(W(1,22  ),W(1,7   ),GG ,MB      ,ZERO    ,W(1,23  ))     
      CALL IOSXXX(W(1,23  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(9   ))         
      CALL IOSXXX(W(1,22  ),W(1,17  ),W(1,6   ),hmcoup ,AMP(10  ))         
C SF: eliminate graphs with gluon emission from the b resulting from top
      AMP(   3)=0
      AMP(   8)=0
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   5)+AMP(   6)-AMP(   7)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   4)+AMP(   7)+AMP(   8)
     &             +AMP(   9)+AMP(  10)
      MATBGTHG_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBGTHG_DC =MATBGTHG_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       

C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BBXTHBX_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> e+ ve b h- b~  
C  
C Crossing   1 is b b~ -> e+ ve b h- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATBBXTHBX_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBBXTHBX_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBBXTHBX_DC
      REAL*8 FUNCTION MATBBXTHBX_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> e+ ve b h- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  15, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 5, 7]T[ 2, 1]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 5, 1]T[ 2, 7]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),MB ,NHEL(1   ),+1*IC(1   ),W(1,1   ))          
      CALL OXXXXX(P(0,2   ),MB ,NHEL(2   ),-1*IC(2   ),W(1,2   ))          
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL IXXXXX(P(0,7   ),MB ,NHEL(7   ),-1*IC(7   ),W(1,7   ))          
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     8   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     9   ))                                                          
      CALL JIOXXX(W(1,7   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL FVOXXX(W(1,9   ),W(1,10  ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,1   ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVIXXX(W(1,1   ),W(1,10  ),GG ,MB      ,ZERO    ,W(1,12  ))     
      CALL IOSXXX(W(1,12  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(2   ))         
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVOXXX(W(1,9   ),W(1,13  ),GG ,MT      ,WT      ,W(1,14  ))     
      CALL IOSXXX(W(1,7   ),W(1,14  ),W(1,6   ),hmcoup ,AMP(3   ))         
      CALL FVIXXX(W(1,7   ),W(1,13  ),GG ,MB      ,ZERO    ,W(1,15  ))     
      CALL IOSXXX(W(1,15  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(4   ))         
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels may not be the same as in the case of production without decay
c Chris:
      if(MHC.le.MT) then
         AMP(   3)=0
      endif

      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      MATBBXTHBX_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBBXTHBX_DC =MATBBXTHBX_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BUTHU_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b u -> e+ ve b h- u  
C  
C Crossing   1 is b u -> e+ ve b h- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATBUTHU_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBUTHU_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBUTHU_DC
      REAL*8 FUNCTION MATBUTHU_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b u -> e+ ve b h- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  12, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 7, 1]T[ 5, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),MB ,NHEL(1   ),+1*IC(1   ),W(1,1   ))          
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL OXXXXX(P(0,7   ),ZERO ,NHEL(7   ),+1*IC(7   ),W(1,7   ))        
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     8   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     9   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL FVOXXX(W(1,9   ),W(1,10  ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,1   ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVIXXX(W(1,1   ),W(1,10  ),GG ,MB      ,ZERO    ,W(1,12  ))     
      CALL IOSXXX(W(1,12  ),W(1,9   ),W(1,6   ),hmcoup ,AMP(2   ))         
      JAMP(   1) = -AMP(   1)-AMP(   2)
      MATBUTHU_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBUTHU_DC =MATBUTHU_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       

C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BBTHB_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b -> e+ ve b h- b  
C  
C Crossing   1 is b b -> e+ ve b h- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATBBTHB_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  72/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBBTHB_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
c SF: insert a factor 2, since identical final-state b's have been treated
c SF as different particles
      ANS(IPROC)=2*ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBBTHB_DC
      REAL*8 FUNCTION MATBBTHB_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b -> e+ ve b h- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  22, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 7, 1]T[ 5, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 7, 2]T[ 5, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),MB ,NHEL(1   ),+1*IC(1   ),W(1,1   ))          
      CALL IXXXXX(P(0,2   ),MB ,NHEL(2   ),+1*IC(2   ),W(1,2   ))          
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL OXXXXX(P(0,7   ),MB ,NHEL(7   ),+1*IC(7   ),W(1,7   ))          
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     8   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,7   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     10  ))                                                          
      CALL FVOXXX(W(1,10  ),W(1,9   ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,2   ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVIXXX(W(1,2   ),W(1,9   ),GG ,MB      ,ZERO    ,W(1,12  ))     
      CALL IOSXXX(W(1,12  ),W(1,10  ),W(1,6   ),hmcoup ,AMP(2   ))         
      CALL FVOXXX(W(1,5   ),W(1,8   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     13  ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,14  ))     
      CALL FVOXXX(W(1,13  ),W(1,14  ),GG ,MT      ,WT      ,W(1,15  ))     
      CALL IOSXXX(W(1,1   ),W(1,15  ),W(1,6   ),hmcoup ,AMP(3   ))         
      CALL FVIXXX(W(1,1   ),W(1,14  ),GG ,MB      ,ZERO    ,W(1,16  ))     
      CALL IOSXXX(W(1,16  ),W(1,13  ),W(1,6   ),hmcoup ,AMP(4   ))         
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,17  ))     
      CALL FVOXXX(W(1,10  ),W(1,17  ),GG ,MT      ,WT      ,W(1,18  ))     
      CALL IOSXXX(W(1,1   ),W(1,18  ),W(1,6   ),hmcoup ,AMP(5   ))         
      CALL FVIXXX(W(1,1   ),W(1,17  ),GG ,MB      ,ZERO    ,W(1,19  ))     
      CALL IOSXXX(W(1,19  ),W(1,10  ),W(1,6   ),hmcoup ,AMP(6   ))         
      CALL JIOXXX(W(1,1   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,20  ))     
      CALL FVOXXX(W(1,13  ),W(1,20  ),GG ,MT      ,WT      ,W(1,21  ))     
      CALL IOSXXX(W(1,2   ),W(1,21  ),W(1,6   ),hmcoup ,AMP(7   ))         
      CALL FVIXXX(W(1,2   ),W(1,20  ),GG ,MB      ,ZERO    ,W(1,22  ))     
      CALL IOSXXX(W(1,22  ),W(1,13  ),W(1,6   ),hmcoup ,AMP(8   ))         
C SF: eliminate graphs in which b quark #7 results from the top decay
      AMP(   1)=0
      AMP(   2)=0
      AMP(   5)=0
      AMP(   6)=0
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)
      JAMP(   2) = +AMP(   5)+AMP(   6)+AMP(   7)+AMP(   8)
      MATBBTHB_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBBTHB_DC =MATBBTHB_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END


C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE UUXTHBX_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> e+ ve b h- b~  
C  
C Crossing   1 is u u~ -> e+ ve b h- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATUUXTHBX_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1, 7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1, 7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1, 7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1, 7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1, 7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1, 7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1, 7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1, 7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1, 7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1, 7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 7) / 1, 2, 3, 4, 5, 6, 7/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATUUXTHBX_DC(P ,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATUUXTHBX_DC
      REAL*8 FUNCTION MATUUXTHBX_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> e+ ve b h- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      integer    nincoming
      parameter (nincoming=  2)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  12, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
C SF Changed first entry of W() from 18 to 6 to be consistent with 
C SF mcatnlo_helas2.f, based on pre-2007 MadEvent
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      include "MEFRcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 2, 7]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),MB ,NHEL(5   ),+1*IC(5   ),W(1,5   ))          
      CALL SXXXXX(P(0,6   ),+1*IC(6   ),W(1,6   ))                         
      CALL IXXXXX(P(0,7   ),MB ,NHEL(7   ),-1*IC(7   ),W(1,7   ))          
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,3   ),W(1,4   ),MGVX34 ,MW      ,WW      ,W(1,       
     &     9   ))                                                          
      CALL FVOXXX(W(1,5   ),W(1,9   ),MGVX34 ,MT      ,WT      ,W(1,       
     &     10  ))                                                          
      CALL FVOXXX(W(1,10  ),W(1,8   ),GG ,MT      ,WT      ,W(1,11  ))     
      CALL IOSXXX(W(1,7   ),W(1,11  ),W(1,6   ),hmcoup ,AMP(1   ))         
      CALL FVIXXX(W(1,7   ),W(1,8   ),GG ,MB      ,ZERO    ,W(1,12  ))     
      CALL IOSXXX(W(1,12  ),W(1,10  ),W(1,6   ),hmcoup ,AMP(2   ))         
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels may not be the same as in the case of production without decay

c Chris:

      if(MHC.le.MT) then
         AMP(   1)=0
      endif

      JAMP(   1) = -AMP(   1)-AMP(   2)
      MATUUXTHBX_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATUUXTHBX_DC =MATUUXTHBX_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
