      subroutine frealwt(s,x,yi,cth2,tk,uk,q1q,q2q,xinv,
     #                   jproc,xmatout)
c Returns the real matrix contributions, times the FKS damping factor: 
c    xii**2*(1-yi**2)=4*tk*uk/s**2    
c Only initial-state kinematics is used here; final-state collinear
c singularities will not occur for Wt. The normalization is such that
c    dsigma_real = g_S^4 g_W^2 |V_tj|^2 (1/xii)_+ [(1/(1-yi))+ + (1/(1+yi))+] 
c                  frealwt (dphi_3/(2*xii))
c with g_S^4, g_W^2 and CKM factors inserted in the main code. 
c
c The matrix elements are passed through the following arrays
c    xmat(out/in)(i):
c    jproc=1 --> gg processes
c    jproc=2 --> q(bar)q(bar) processes
c    jproc=3 --> gq(bar),q(bar)g processes
c The meaning of i(==idr) depends on jproc -- see notes
c
      implicit none
      real * 8 s,x,yi,cth2,tk,uk,q1q,q2q,xinv(5),xmatout(7)
      integer jproc
      include 'stpcblks.h'
      real * 8 xii,xmatin(7)
      integer i,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
c Hard-coded choice for Wt mode; for s- and t-channel see freal()
      parameter (jwt=3)
c
      xii=1-x
      if(jproc.eq.1)then
        call frealwt_gg(s,xii,yi,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jwt)
          xmatout(i)=xmatin(i)
        enddo
      elseif(jproc.eq.2)then
        call frealwt_qq(s,xii,yi,cth2,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jwt)
          xmatout(i)=xmatin(i)
        enddo
      elseif(jproc.eq.3)then
        call frealwt_qg(s,xii,yi,cth2,tk,uk,q1q,q2q,xinv,xmatin)
        do i=1,idrmax(jproc,jwt)
          xmatout(i)=xmatin(i)
        enddo
      else
        write(*,*)'Unknown process in frealwt',jproc
        stop
      endif
c
      return
      end


      subroutine frealwt_gg(xs,xxii,xyi,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for gg --> tWbbar. See frealwt for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xtk,xuk,xq1q,xq2q,xxinv(5),xmatout(7)
      real * 8 s,xii,yi,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vtr,xnc,s_red,t_red,x_ap,ap_kern,xfact,xinv(5),xmatin(7)
      integer i,jwt,ione,ithree,icode,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      parameter (tiny=1.d-6)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      parameter (jwt=3)
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
      xnorm=vtr/(32*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(ione,jwt)
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
        call fbornwt(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ione,jwt)
          if(idrlimcp(jwt,ione,i).ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jwt,ione,i))
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        icode=2
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornwt(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ione,jwt)
          if(idrlimcm(jwt,ione,i).ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jwt,ione,i))
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call wtmat_gg(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(ione,jwt)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
c
      return
      end


      subroutine frealwt_qq(xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for qq --> tWbbar. See frealwt for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,xxinv(5),xmatout(7)
      real * 8 s,xii,yi,cth2,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vtr,xnc,s_red,t_red,x_ap,ap_kern,qin_kern,xfact,xinv(5),
     # xmatin(7),xaziin(7)
      integer i,jwt,itwo,ithree,ipone,imone,icode,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      parameter (tiny=1.d-6)
      parameter (pi=3.14159265358979312D0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      parameter (jwt=3)
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
      xnorm=vtr/(32*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(itwo,jwt)
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
        call fbornwt(s_red,t_red,ithree,xmatin)
        call faziwt(s_red,t_red,cth2,itwo,ipone,xaziin)
        do i=1,idrmax(itwo,jwt)
          if(idrlimcp(jwt,itwo,i).ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jwt,itwo,i)) +
     #                 4*(1+yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        icode=3
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornwt(s_red,t_red,ithree,xmatin)
        call faziwt(s_red,t_red,cth2,itwo,imone,xaziin)
        do i=1,idrmax(itwo,jwt)
          if(idrlimcm(jwt,itwo,i).ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jwt,itwo,i)) +
     #                 4*(1-yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call wtmat_qq(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(itwo,jwt)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
c
      return
      end


      subroutine frealwt_qg(xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,
     #                      xxinv,xmatout)
c Real matrix elements for qg --> tWbbar. See frealwt for comments
      implicit none
      include 'stpcblks.h'
      real * 8 xs,xxii,xyi,xcth2,xtk,xuk,xq1q,xq2q,xxinv(5),
     # xmatout(7)
      real * 8 s,xii,yi,cth2,tk,uk,q1q,q2q,q1c,q2c,xnorm,tiny,pi,vcf,
     # vca,vtr,xnc,s_red,t_red,x_ap,ap_kern,qin_kern,xfact,xinv(5),
     # xmatin(7),softfc(3),xaziin(7)
      integer i,jwt,ithree,ipone,imone,icode,idrmax(1:3,3)
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
      parameter (jwt=3)
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
      xnorm=vtr/(32*xnc)
      do i=1,5
        xinv(i)=xxinv(i)
      enddo
      do i=1,idrmax(ithree,jwt)
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
        call fbornwt(s_red,t_red,ithree,xmatin)
        do i=1,idrmax(ithree,jwt)
             xmatout(i)=softfc(i)*xmatin(i) 
        enddo
      elseif(abs(1-yi).lt.tiny)then
c Collinear + limit
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q2q
        call fbornwt(s_red,t_red,ithree,xmatin)
        call faziwt(s_red,t_red,cth2,ithree,ipone,xaziin)
        do i=1,idrmax(ithree,jwt)
          icode=iapcp(jwt,i)
          if(idrlimcp(jwt,ithree,i).ne.0.and.icode.ne.0)
     #      xmatout(i)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcp(jwt,ithree,i)) +
     #                 4*(1+yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      elseif(abs(1+yi).lt.tiny)then
c Collinear - limit
        x_ap=1-xii
        s_red=s*x_ap
        t_red=q1q
        call fbornwt(s_red,t_red,ithree,xmatin)
        call faziwt(s_red,t_red,cth2,ithree,imone,xaziin)
        do i=1,idrmax(ithree,jwt)
          icode=iapcm(jwt,i)
          if(idrlimcm(jwt,ithree,i).ne.0.and.icode.ne.0)
     #      xmatout(i)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*
     #                 xmatin(idrlimcm(jwt,ithree,i)) +
     #                 4*(1-yi)/s*qin_kern(x_ap,abs(icode))*
     #                 xaziin(i)
        enddo
      else
c The kinematical configuration is not in one of the singular regions:
c use the full expression for the matrix element.
c Factor multiplying the matrix element, eq.(4.37)
        xfact=(1-yi**2)*xii**2
        call wtmat_qg(s,tk,uk,q1q,q2q,xmatin)
        do i=1,idrmax(ithree,jwt)
          xmatout(i)=xnorm*xfact*xmatin(i)/(2*s)
        enddo
      endif
      return
      end


      subroutine f2svwt(s,t,jproc,xmatout)
c Soft-virtual contribution. The normalization is such that
      implicit none
      include 'stpcblks.h'
      real*8 s,t,xmatout(7)
      real*8 pi,vca,vtr,u,beta0,scalefac,svfac,svnonfac,
     # xfac(3),xnonfac(3),xmatin(7)
      parameter (pi=3.14159265358979312D0)
      parameter (vca=3.d0)
      parameter (vtr=0.5d0)
      integer idrmax(1:3,3)
      common/cidrmax/idrmax
      integer jproc,jwt,ithree,i
      parameter (jwt=3)
      parameter (ithree=3)
c
      if(jproc.lt.1.or.jproc.gt.3)then
        write(*,*)'Error in f2svwt: jproc=',jproc
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
        scalefac=2*pi*beta0*log(xmur2/xmuf2h1)
c Soft-virtual function have been derived for idr=3. See finiteWt.m
        xfac(1)=svfac(s,u,t,xm12,xm22,xmuf2h1,nl)
        xnonfac(1)=svnonfac(s,u,t,xm12,xm22,xmuf2h1,nl)
        xfac(3)=svfac(s,t,u,xm12,xm22,xmuf2h1,nl)
        xnonfac(3)=svnonfac(s,t,u,xm12,xm22,xmuf2h1,nl)
        call fbornwt(s,t,ithree,xmatin)
        do i=1,idrmax(jproc,jwt)
          xmatout(i)=(xfac(i)+scalefac)*xmatin(i)+xnonfac(i)/(2*s)
        enddo
      endif
      return
      end


      subroutine f2prwt(xs,xt,xx,xxc,xyic,xxlmude,jproc,xmatout)
c Collinear reminder contribution of FKS. The normalization is such that
c    dsigma_pr = g_S^2 alpha_S/(2*pi) g_W^2 |V_tj|^2 f2prwt/xii dxii dphi_2
c with g_S^4, g_W^2 and CKM factors inserted in the main code. See frealwt
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
      integer jwt,ithree,i,jproc,icode,ilim
      integer idrmax(1:3,3)
      common/cidrmax/idrmax
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      integer iapcp(3,7),iapcm(3,7)
      common/ciap/iapcp,iapcm
      parameter (one=1.d0)
      parameter (jwt=3)
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
c here it is the same as entered in the call to f2prwt
      call fbornwt(s_red,t,ithree,xmatin)
      if(yic.eq.1.d0)then
        scheme=schhad1
      elseif(yic.eq.-1.d0)then
        scheme=schhad2
      else
        write(6,*)'Error in f2prwt',yic
        stop
      endif
c
      do i=1,idrmax(jproc,jwt)
        if(jproc.eq.1)then
          icode=2
        elseif(jproc.eq.2)then
          icode=3
        elseif(jproc.eq.3)then
          if(yic.eq.1.d0)then
            icode=iapcp(jwt,i)
          elseif(yic.eq.-1.d0)then
            icode=iapcm(jwt,i)
          endif
        else
          write(*,*)'Unknown process in f2prwt',jproc
          stop
        endif
        if(yic.eq.1.d0)then
          ilim=idrlimcp(jwt,jproc,i)
        elseif(yic.eq.-1.d0)then
          ilim=idrlimcm(jwt,jproc,i)
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
            write(6,*)'Error in f2prwt, y=',yic
            write(6,*)'Factorization scheme ',scheme,' not known'
          endif
c 1/xi is the main code
          tmp=xdfct1*(xlmude+2*log(xii))-xdfct2
     #       -xdfct3p-xdfct3l*log(xii) 
     #       -xii*xdfct5
        else
          tmp=0.d0
        endif
c
        if(ilim.ne.0)then
          if(icode.eq.0)then
            write(6,*)'Inconsistency in f2prwt'
            write(6,*)ilim,icode,i,jproc
            stop
          endif
          xmatout(i)=tmp*xmatin(ilim)
        else
          if(icode.ne.0.and.jproc.eq.3)then
            write(6,*)'Inconsistency in f2prwt'
            write(6,*)ilim,icode,i,jproc
            stop
          endif
          xmatout(i)=0.d0
        endif
      enddo
      return
      end


      subroutine fbornwt(s,t,jproc,xmatout)
c Returns the partonic Born contribution, with the following normalization
c    dsigma_born = g_S^2 g_W^2 |V_tj|^2 fborn dphi_2
c with g_S^2, g_W^2 and CKM factors inserted in the main code. See frealwt
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
      integer jproc,i,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (jwt=3)
c
      if(jproc.lt.1.or.jproc.gt.3)then
        write(*,*)'Error in fbornwt: jproc=',jproc
        stop
      endif
c
      do i=1,idrmax(jproc,jwt)
        xmatout(i) = 0d0
      enddo
c
      if(jproc.eq.3)then
        xnorm=vtr/(32*xnc)
        u=xm12+xm22-s-t
        xmatout(1)=xnorm*elborn(s,u)/(2*s)
        xmatout(3)=xnorm*elborn(s,t)/(2*s)
      endif
c
      return
      end


      subroutine faziwt(s,t,cth2,jproc,icoll,xmatout)
c Returns the azimuthal-dependent part of collinear limits; the definitions
c and normalization are as in fbornwt
      implicit none
      include 'stpcblks.h'
      real*8 s,t,cth2,xmatout(7)
      real*8 vtr,xnc,xnorm,u,azidep,azi1,azi3
      parameter (vtr=0.5d0)
      parameter (xnc=3.d0)
      integer jproc,icoll,i,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (jwt=3)
c
      if(jproc.lt.2.or.jproc.gt.3)then
        write(*,*)'Error in faziwt: jproc=',jproc
        stop
      endif
c
      do i=1,idrmax(jproc,jwt)
        xmatout(i) = 0d0
      enddo
c
      xnorm=vtr/(32*xnc)
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
        write(*,*)'Error in faziwt: icoll=',icoll
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
      subroutine wtmat_gg(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,
     # t23,tpq,tp3,tq3,W_gg_nores
      integer i,ione,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (ione=1)
      parameter (jwt=3)
c
      do i=1,idrmax(ione,jwt)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c See the notes for the mapping of invariants. W_gg calls WQqgg_1 and
c WQqgg_2, which are computed with
c   W(q) --> t(p) bbar(k1) g(k2) g(k3)
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
      xmatout(1)=W_gg_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
c
      return
      end


      subroutine wtmat_qq(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,t23,tpq,
     # tp3,tq3,Nc,WQqqq_1,WQqqq_2,WQqqq_1_nores,WQqqq_2_nores
      integer i,itwo,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (itwo=2)
      parameter (jwt=3)
      parameter(Nc = 3d0)

      do i=1,idrmax(itwo,jwt)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c See the notes for the mapping of invariants. WQqqq_1 and WQqqq_2 are
c computed with 
c   W(q) --> t(p) bbar(k1) ubar(k2) u(k3)
c   W(q) --> t(p) bbar(k1) bbar(k2) b(k3)
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
      xmatout(1)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  WQqqq_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
      xmatout(2)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  WQqqq_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
      xmatout(3)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  WQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
      xmatout(4)=(Nc**2-1.)/(2.*Nc)*(-1)*
     .  WQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
      xmatout(5)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .  WQqqq_1(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)+
     .  WQqqq_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
     . +0.5*(Nc/4.-1./(4.*Nc))*(
     .  WQqqq_2_nores(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)) )
      xmatout(6)=(Nc**2-1.)/(2.*Nc)*(-1)*(
     .  WQqqq_1(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)+
     .  WQqqq_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
     . +0.5*(Nc/4.-1./(4.*Nc))*(
     .  WQqqq_2_nores(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)) )
      xmatout(7)= (Nc**2-1.)/(2.*Nc)*(-1)*(
     .  WQqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)+
     .  WQqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
     .  +  0.5*(Nc/4.-1./(4.*Nc))*(
     .  WQqqq_2(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)) )
c
      return
      end


      subroutine wtmat_qg(s,tk,uk,q1q,q2q,xmatout)
      implicit none
      real*8 s,tk,uk,q1q,q2q,xmatout(7)
      include 'stpcblks.h'
      real*8 s2,q1c,q2c,w1,w2,t12,t1p,t1q,t13,t2p,t2q,
     # t23,tpq,tp3,tq3,W_qg
      integer i,ithree,jwt,idrmax(1:3,3)
      common/cidrmax/idrmax
      parameter (ithree=3)
      parameter (jwt=3)
c
      do i=1,idrmax(ithree,jwt)
        xmatout(i) = 0d0
      enddo
c
      s2 = s+tk+uk
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c
c See the notes for the mapping of invariants. W_qg calls WQqgg_1 and
c WQqgg_2, which are computed with
c   W(q) --> t(p) bbar(k1) g(k2) g(k3)
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
      xmatout(1)=W_qg(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
      xmatout(3)=W_qg(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
c
      return
      end


      function elborn(s,t)
c Original Born function; checked against two independent computations
      implicit none
      real*8 elborn,s,t
      include 'stpcblks.h'
      real*8 m2,t1,q2,MbMb
c
      m2 = xm12
      t1 = t-m2
      q2 = xm22
c
      MbMb =
     &  + s**(-1)*t1**(-1) * (  - 32*q2**2 + 48*m2*q2 - 16*m2**3*
     &    q2**(-1) )
      MbMb = MbMb + s**(-1) * ( 32*q2 - 16*m2 - 16*m2**2*q2**(-1) )
      MbMb = MbMb + s**(-1)*t1 * (  - 16 - 8*m2*q2**(-1) )
      MbMb = MbMb + t1**(-2) * ( 32*m2*q2 - 16*m2**2 - 16*m2**3*
     &    q2**(-1) )
      MbMb = MbMb + t1**(-1) * ( 32*q2 - 16*m2 - 16*m2**2*q2**(-1) )
      MbMb = MbMb + s*t1**(-1) * (  - 16 - 8*m2*q2**(-1) )
      MbMb = MbMb - 16*m2*q2**(-1)           
c
      elborn = MbMb
      return
      end


      function azidep(s,t,u,m2,mw2,cth2)
c Returns the azimuthal-dependent part for the y=-1 collinear singularity
c of any real process that factorizes
c  b(p1)+g(p2) --> t(k1)+W(k2) 
c The definitions of the invariants and the normalizations are as in fbornwt.
c See also the comments in born.mac and aziborn.mac
      implicit none
      real * 8 azidep,s,t,u,m2,mw2,cth2
c
      azidep = -8*(2*cth2**2-1)*(mw2-m2)*(2*mw2+m2)*(m2*mw2-t*u)/(mw2*s*
     1   (u-m2)**2)
      return 
      end 


*---------------------------------------------------------------
      real*8 function W_gg_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2, Nc
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)
      parameter(Nc = 3d0)

      real*8 WQqgg_1_nores, WQqgg_2_nores
      real*8 m2, q2, s12, s13, s23, sQ1, sQ2, sQ3

      W_gg_nores = 0.5*
     .  (WQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  WQqgg_1_nores(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2) +
     .  (WQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  WQqgg_1_nores(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2))/(Nc**2-1d0)+
     .  WQqgg_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)*(-1d0/(Nc**2-1d0)))
      return
      end


*---------------------------------------------------------------
      real*8 function W_qg(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2, Nc
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)
      parameter(Nc = 3d0)

      real*8 WQqgg_1, WQqgg_2
      real*8 m2, q2, s12, s13, s23, sQ1, sQ2, sQ3

C Overall minus sign due to fermion crossing!
      W_qg = -0.5*Nc*(
     .  WQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3) +
     .  WQqgg_1(q2,m2,s13,s12,s23,sQ1,sQ3,sQ2) ) 
     . + WQqgg_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)*( 1d0/(2*Nc) )

      return
      end


*---------------------------------------------------------------
      real*8 function WQqqq_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqqq1
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C WQqqq1 corresponds to the f+g graphs squared (two fermion loops):
C They contribute at 1st and 3rd color order
C WQqqq2 corresponds to the interference graphs. They contribute 
C at 2nd and 4th color order.

      WQqqq1 =0d0
      WQqqq1 = WQqqq1 + D2**(-2)*s12*s23**(-1) * (  - 32*m2 - 32*m2**2*
     +    q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-2)*s12**2*s23**(-2) * (  - 32*m2 - 32*
     +    m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-2) * (  - 16*m2 - 16*m2**2*q2**(-1) + 32*
     +    q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ3 * (  - 64*
     +    m2 - 64*m2**2*q2**(-1) + 128*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 64*m2*
     +    q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s12*s23**(-1) * (  - 32*m2 - 
     +    32*m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * (  - 64 - 
     +    32*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3 * (  - 32*m2 - 
     +    32*m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 64 - 
     +    32*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s23**(-1) * ( 96*m2*q2 - 32*
     +    m2**3*q2**(-1) - 64*q2**2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*s13 * ( 64 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1)*sQ3 * (  - 64 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*D4**(-1) * (  - 32*m2**2*q2**(-1) - 64
     +    *q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*s12*s23**(-2) * ( 64*m2 + 64*m2**2*
     +    q2**(-1) - 128*q2 )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*s12*s23**(-1) * (  - 32*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*s23**(-1)*sQ2 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D2**(-1)*s23**(-1)*sQ3 * ( 32 + 16*m2*q2**(-1)
     +     )
     +
      WQqqq1 = WQqqq1 + D2**(-1) * ( 32 - 16*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D4**(-2)*s23**(-2)*sQ3**2 * (  - 32*m2 - 32*
     +    m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D4**(-2)*s23**(-1)*sQ3 * (  - 32*m2 - 32*m2**2*
     +    q2**(-1) + 64*q2 )
     +
      WQqqq1 = WQqqq1 + D4**(-2)*s23**(-1) * ( 64*m2*q2 - 32*m2**2 - 32
     +    *m2**3*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D4**(-2) * (  - 16*m2 - 16*m2**2*q2**(-1) + 32*
     +    q2 )
     +
      WQqqq1 = WQqqq1 + D4**(-1)*s12*s23**(-1) * ( 32 + 16*m2*q2**(-1)
     +     )
     +
      WQqqq1 = WQqqq1 + D4**(-1)*s23**(-2)*sQ3 * ( 64*m2 + 64*m2**2*
     +    q2**(-1) - 128*q2 )
     +
      WQqqq1 = WQqqq1 + D4**(-1)*s23**(-1)*s13 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D4**(-1)*s23**(-1)*sQ3 * (  - 32*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + D4**(-1) * (  - 32 - 16*m2*q2**(-1) )
     +
      WQqqq1 = WQqqq1 + s23**(-2) * (  - 32*m2 - 32*m2**2*q2**(-1) + 64
     +    *q2 )

      WQqqq_1 = WQqqq1

      return
      end


*---------------------------------------------------------------
      real*8 function WQqqq_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqqq2
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C WQqqq1 corresponds to the f+g graphs squared (two fermion loops):
C They contribute at 1st and 3rd color order
C WQqqq2 corresponds to the interference graphs. They contribute 
C at 2nd and 4th color order.


      WQqqq2 =0d0
      WQqqq2 = WQqqq2 + D2**(-2)*s12**2*s23**(-1)*s13**(-1) * ( 32*m2
     +     + 32*m2**2*q2**(-1) - 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * ( 32*m2 + 32*m2**2*q2**(-1) - 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ3 * (  - 32*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s12*s13**(-1) * (  - 32*m2**2
     +    *q2**(-1) - 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3 * ( 32*m2 + 32*
     +    m2**2*q2**(-1) - 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1)*s13**(-1) * (  - 48*m2*q2 + 
     +    16*m2**3*q2**(-1) + 32*q2**2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2*
     + sQ3 * (  - 160 - 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2
     +  * (  - 128*m2 - 32*m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ2**2*sQ3 * ( 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ2**2 * (  - 32 + 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ3
     +  * (  - 16*m2 - 48*m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 32 - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1) * ( 
     +    64*m2*q2 - 32*m2**2 - 32*m2**3*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*sQ2*sQ3 * ( 32
     +    *q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*sQ2 * (  - 32
     +     + 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1)*sQ3 * (  - 32
     +     - 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s23**(-1) * (  - 64*m2
     +     - 32*m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1)*sQ2*sQ3 * ( 96
     +    *q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1)*sQ2 * (  - 96
     +     + 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1)*sQ3 * (  - 128
     +     - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*s13**(-1) * (  - 160*m2
     +     + 64*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12*sQ3 * (  - 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12 * ( 32 - 16*m2*q2**(-1)
     +     )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2*sQ3 * ( 96*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2 * (  - 32 + 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * (  - 128 - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)
     +  * (  - 96*m2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1)*sQ3 * ( 32*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s23**(-1) * ( 32*m2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s13**(-1)*sQ3 * ( 64*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**2*s13**(-1) * ( 32*m2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ3 * ( 64*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s12**3*s23**(-1)*s13**(-1)
     +  * ( 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2*sQ3
     +  * (  - 16*m2 - 16*m2**2*q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2 * ( 
     +    48*m2*q2 - 16*m2**3*q2**(-1) - 32*q2**2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2**2*
     + sQ3 * (  - 32 - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2**2
     +  * (  - 16*m2 - 16*m2**2*q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*sQ2*sQ3 * (  - 32
     +     - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s23**(-1)*sQ2 * (  - 16*m2
     +     - 16*m2**2*q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ2*sQ3 * (  - 96
     +     - 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ2 * (  - 48*m2
     +     - 48*m2**2*q2**(-1) + 96*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ3 * (  - 48*m2
     +     - 48*m2**2*q2**(-1) + 96*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1)*sQ3**2 * (  - 32
     +     - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*s13**(-1) * ( 96*m2*q2 - 32*
     +    m2**3*q2**(-1) - 64*q2**2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1)*sQ3 * ( 32 + 16*m2*q2**(-1)
     +     )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*D4p**(-1) * ( 16*m2 + 16*m2**2*
     +    q2**(-1) - 32*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1)*sQ2*sQ3 * ( 32
     +    *q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1)*sQ2 * (  - 32
     +     + 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * (  - 64
     +     - 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s23**(-1)*s13**(-1) * (  - 112*m2
     +     - 80*m2**2*q2**(-1) + 96*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s13**(-1)*sQ3 * ( 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12*s13**(-1) * (  - 32 + 16*m2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12**2*s23**(-1)*s13**(-1)*sQ3 * ( 32*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s12**2*s23**(-1)*s13**(-1) * ( 32 + 16
     +    *m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s23**(-1)*s13**(-1)*sQ2*sQ3 * (  - 32
     +     - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s23**(-1)*s13**(-1)*sQ2 * (  - 16*m2
     +     - 16*m2**2*q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s13**(-1)*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D2**(-1)*s13**(-1) * (  - 16*m2 - 16*m2**2*
     +    q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 64 - 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s23**(-1)*sQ3 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s23*s13**(-1) * (  - 32*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s13**(-1)*sQ2 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s13**(-1)*sQ3 * (  - 64*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12*s13**(-1) * ( 16*m2 - 32
     +    *m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12 * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2*sQ3
     +  * (  - 32*m2 - 16*m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ2*
     + sQ3**2 * (  - 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ3 * ( 
     +    32*m2*q2 - 16*m2**2 - 16*m2**3*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13**(-1)*sQ3**2
     +  * (  - 64*m2**2*q2**(-1) + 64*q2 )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*s13*sQ3 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*sQ2*sQ3 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*sQ3 * (  - 16*m2
     +     - 32*m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23**(-1)*sQ3**2 * (  - 32*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23*s13**(-1)*sQ3 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s23*s13**(-1) * ( 64*m2 - 32
     +    *m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1)*sQ2*sQ3 * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1)*sQ2 * ( 48*m2 - 16
     +    *m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1)*sQ3 * ( 16*m2 - 64
     +    *m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1)*sQ3**2 * (  - 32*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*s13**(-1) * (  - 16*m2*q2 + 
     +    48*m2**2 - 32*m2**3*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1)*sQ3 * (  - 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*D4p**(-1) * ( 48*m2 - 16*m2**2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * ( 32 - 
     +    16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s12*s13**(-1) * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3 * (  - 48*m2
     +     - 32*m2**2*q2**(-1) + 32*q2 )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s23**(-1)*s13**(-1)*sQ3**2 * (  - 32
     +     - 48*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s23**(-1)*sQ3 * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s13**(-1)*sQ3 * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4**(-1)*s13**(-1) * ( 48*m2 - 16*m2**2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2*sQ3 * ( 
     +     - 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ2 * (  - 16
     +    *m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 * ( 96 + 
     +    32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s23**(-1)*s13**(-1) * ( 48*m2 )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s13**(-1)*sQ3 * ( 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12*s13**(-1) * (  - 64 + 32*m2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12**2*s23**(-1)*s13**(-1)*sQ3 * ( 
     +     - 32*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s12**2*s23**(-1)*s13**(-1) * (  - 16*
     +    m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ2*sQ3 * ( 64 + 
     +    64*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ2 * ( 32*m2 + 
     +    16*m2**2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1)*sQ3 * ( 64*m2**2*
     +    q2**(-1) - 64*q2 )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*s13**(-1) * (  - 32*m2*q2
     +     + 16*m2**2 + 16*m2**3*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*s13 * ( 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*sQ2 * ( 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1)*sQ3 * ( 32*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23**(-1) * ( 16*m2 + 32*m2**2*
     +    q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s23*s13**(-1) * (  - 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + D4p**(-1)*s13**(-1) * (  - 32*m2 - 32*m2**2*
     +    q2**(-1) + 64*q2 )
     +
      WQqqq2 = WQqqq2 + D4p**(-1) * ( 16*m2*q2**(-1) )
     +
      WQqqq2 = WQqqq2 + s12*s23**(-1)*s13**(-1) * (  - 32 )
     +
      WQqqq2 = WQqqq2 + s23**(-1)*s13**(-1)*sQ3 * ( 32 + 48*m2*q2**(-1)
     +     )
     +
      WQqqq2 = WQqqq2 + s23**(-1)*s13**(-1) * ( 48*m2 + 32*m2**2*
     +    q2**(-1) - 32*q2 )
     +
      WQqqq2 = WQqqq2 + s23**(-1) * ( 16*m2*q2**(-1) )

      WQqqq_2 = WQqqq2

      return
      end


*---------------------------------------------------------------
      real*8 function WQqgg_1(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqgg
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2


      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      WQqgg = 0d0
      WQqgg = WQqgg + D2**(-2)*s12**(-1)*s23 * (  - 16*m**2 - 16*m**4*
     +    q2**(-1) + 32*q2 )
     +
      WQqgg = WQqgg + D2**(-2)*s12*s23**(-1) * (  - 64*m**2 - 64*m**4*
     +    q2**(-1) + 128*q2 )
     +
      WQqgg = WQqgg + D2**(-2)*s12**2*s23**(-2) * (  - 32*m**2 - 32*
     +    m**4*q2**(-1) + 64*q2 )
     +
      WQqgg = WQqgg + D2**(-2) * (  - 48*m**2 - 48*m**4*q2**(-1) + 96*
     +    q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ3**(-1) * (  - 
     +    96*m**2*q2 + 48*m**4 - 16*m**6*q2**(-1) + 64*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*s23 * ( 32*m**2 - 16*
     +    m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ3**(-1) * ( 
     +     - 32*m**2 + 32*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*s23**2 * (  - 32 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*sQ3**(-1) * (  - 160*
     +    m**2*q2**2 + 96*m**4*q2 + 32*m**6 - 32*m**8*q2**(-1) + 64*
     +    q2**3 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*sQ3 * (  - 48*m**2 - 
     +    48*m**4*q2**(-1) + 96*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1)*sQ3**2 * (  - 32 - 16
     +    *m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**(-1) * ( 192*m**2*q2 - 64*
     +    m**6*q2**(-1) - 128*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ3 * (  - 64*
     +    m**2 - 64*m**4*q2**(-1) + 128*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3**(-1) * (  - 
     +    128*m**2*q2 + 64*m**4 + 64*m**6*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * ( 64*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12*s23**(-1) * (  - 64*m**2 - 
     +    64*m**4*q2**(-1) + 128*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12*sQ3**(-1) * (  - 176*m**2
     +     - 48*m**4*q2**(-1) + 96*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12 * (  - 96 + 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * (  - 64 - 32
     +    *m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s12**2*sQ3**(-1) * (  - 32 - 16
     +    *m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23**(-1)*sQ3 * (  - 64*m**2 - 
     +    64*m**4*q2**(-1) + 128*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 64 - 32
     +    *m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23**(-1) * ( 384*m**2*q2 - 128
     +    *m**6*q2**(-1) - 256*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23*sQ3**(-1) * (  - 96*m**2 - 
     +    16*m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23 * (  - 192 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*s23**2*sQ3**(-1) * (  - 32 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*sQ3**(-1) * ( 64*m**2*q2 + 64*
     +    m**4 - 128*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1)*sQ3 * (  - 96 + 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*D4**(-1) * ( 96*m**2 - 192*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*s23*sQ3**(-1) * (  - 16*m**2
     +     - 48*m**4*q2**(-1) - 32*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*sQ2*sQ3**(-1) * (  - 16*m**2
     +     - 16*m**4*q2**(-1) + 32*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*sQ3**(-1) * ( 96*m**2*q2 - 32*
     +    m**6*q2**(-1) - 64*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-2) * ( 64*m**2 + 64*m**4*
     +    q2**(-1) - 128*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-1)*sQ2*sQ3**(-1) * ( 32*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-1)*sQ3**(-1) * ( 32*m**2 + 32
     +    *m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-1) * (  - 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*s12**2*s23**(-1)*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ2*sQ3**(-1) * ( 32*m**2 + 32
     +    *m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ2 * (  - 32 - 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ3**(-1) * (  - 96*m**2*q2 + 
     +    32*m**6*q2**(-1) + 64*q2**2 )
     +
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ3 * ( 32 + 16*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + D2**(-1)*s23*sQ3**(-1) * ( 32 )
     +
      WQqgg = WQqgg + D2**(-1)*sQ2*sQ3**(-1) * (  - 32 + 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D2**(-1)*sQ3**(-1) * (  - 64*m**2 - 48*m**4*
     +    q2**(-1) + 64*q2 )
     +
      WQqgg = WQqgg + D2**(-1) * ( 64 )
     +
      WQqgg = WQqgg + D4**(-2)*s23**(-2)*sQ3**2 * (  - 32*m**2 - 32*
     +    m**4*q2**(-1) + 64*q2 )
     +
      WQqgg = WQqgg + D4**(-2)*s23**(-1)*sQ3 * (  - 64*m**2 - 64*m**4*
     +    q2**(-1) + 128*q2 )
     +
      WQqgg = WQqgg + D4**(-2)*s23**(-1) * ( 128*m**2*q2 - 64*m**4 - 64
     +    *m**6*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-2)*s23*sQ3**(-1) * (  - 16*m**2 - 16*m**4*
     +    q2**(-1) + 32*q2 )
     +
      WQqgg = WQqgg + D4**(-2)*sQ3**(-2) * ( 128*m**4*q2 - 64*m**6 - 64
     +    *m**8*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-2)*sQ3**(-1) * ( 128*m**2*q2 - 64*m**4 - 64
     +    *m**6*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-2) * (  - 48*m**2 - 48*m**4*q2**(-1) + 96*
     +    q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ3**(-1) * ( 
     +     - 128*m**2*q2 + 64*m**4 + 64*m**6*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ3 * ( 32*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*s13 * ( 32*m**2 + 32
     +    *m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*s13**2 * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*sQ3 * ( 32*m**2 + 32
     +    *m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1)*sQ3**2 * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23**(-1) * (  - 96*m**2*q2 + 
     +    32*m**6*q2**(-1) + 64*q2**2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23*sQ3**(-1) * ( 48*m**2 - 16
     +    *m**4*q2**(-1) - 32*q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s23 * ( 32 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s13*sQ3**(-1) * (  - 16*m**2
     +     - 16*m**4*q2**(-1) + 32*q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*s13 * (  - 32 + 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*sQ3**(-2) * (  - 128*m**2*
     +    q2**2 + 192*m**4*q2 - 64*m**8*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1)*sQ3**(-1) * ( 96*m**2*q2 - 32*
     +    m**6*q2**(-1) - 64*q2**2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12**(-1) * (  - 16*m**4*q2**(-1) + 64*
     +    q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s12*s23**(-1) * ( 32 + 16*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + D4**(-1)*s12*sQ3**(-2) * (  - 64*m**2 - 32*m**4*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s23**(-2)*sQ3 * ( 64*m**2 + 64*m**4*
     +    q2**(-1) - 128*q2 )
     +
      WQqgg = WQqgg + D4**(-1)*s23**(-1)*s13 * (  - 32 - 16*m**2*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*s23**(-1)*sQ3 * (  - 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*sQ3**(-2) * ( 128*m**2*q2 - 64*m**4 - 64
     +    *m**6*q2**(-1) )
     +
      WQqgg = WQqgg + D4**(-1)*sQ3**(-1) * ( 128*m**2 )
     +
      WQqgg = WQqgg + D4**(-1) * ( 64 )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-1) * ( 32*m**2
     +    *q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ3**(-1) * ( 32*m**2 + 
     +    32*m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13 * ( 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13**2*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2*sQ3**(-1) * ( 32*m**2 + 
     +    32*m**4*q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2 * ( 64 + 32*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ3**(-1) * (  - 96*m**2*q2
     +     + 32*m**6*q2**(-1) + 64*q2**2 )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ3 * ( 32 + 16*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s12**(-1)*s23**(-1) * ( 32*m**2 + 32*m**4*
     +    q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + s12**(-1)*s23*sQ3**(-2) * (  - 64*m**2 - 32*m**4*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*s23*sQ3**(-1) * ( 64 + 16*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s12**(-1)*s13*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s12**(-1)*sQ2*sQ3**(-2) * (  - 64*m**2 - 32*m**4*
     +    q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*sQ2*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s12**(-1)*sQ3**(-2) * ( 128*m**2*q2 - 64*m**4 - 
     +    64*m**6*q2**(-1) )
     +
      WQqgg = WQqgg + s12**(-1)*sQ3**(-1) * ( 16*m**4*q2**(-1) - 64*q2
     +     )
     +
      WQqgg = WQqgg + s12**(-1) * ( 64 + 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + s12*s23**(-1)*sQ3**(-1) * ( 32 + 16*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s23**(-2) * (  - 32*m**2 - 32*m**4*q2**(-1) + 64*
     +    q2 )
     +
      WQqgg = WQqgg + s23**(-1)*s13*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1)
     +     )
     +
      WQqgg = WQqgg + s23**(-1)*sQ2*sQ3**(-1) * ( 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + s23**(-1)*sQ3**(-1) * ( 32*m**2 + 32*m**4*
     +    q2**(-1) - 64*q2 )
     +
      WQqgg = WQqgg + s23**(-1) * ( 32*m**2*q2**(-1) )
     +
      WQqgg = WQqgg + sQ3**(-2) * (  - 64*m**4*q2**(-1) )
     +
      WQqgg = WQqgg + sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )


      WQqgg_1 = WQqgg


      return
      end


*---------------------------------------------------------------
      real*8 function WQqgg_2(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 W
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2



      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      
      W = 0d0
      W = W + D2**(-2)*s12**(-1)*s13*sQ1 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-2)*s12**(-1)*s13*sQ2 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-2)*s12**(-1)*s13*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-2)*s12*s13**(-1)*sQ1 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-2)*s12*s13**(-1)*sQ2 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-2)*s12*s13**(-1)*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)*sQ3
     +  * (  - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)
     +  * (  - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*sQ3**(-1)
     +  * ( 32 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1)
     +  * (  - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1 * (  - 32 - 
     +    16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*sQ2**(-1)
     +  * ( 32 + 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*sQ3**(-1)
     +  * ( 96 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1)
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2 * (  - 32 - 
     +    16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13 * ( 32*m2 + 16*
     +    m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*sQ2**(-1)
     +  * ( 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*sQ3**(-1)
     +  * ( 32 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*sQ3**(-1)
     +  * ( 32 )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ3**(-1) * ( 
     +     - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**2 * (  - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-1) * ( 
     +    32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**(-1)*sQ3
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2*sQ3**(-1)
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2 * ( 32 - 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**2*sQ3**(-1)
     +  * ( 32 - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ3 * (  - 32*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1 * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1)*sQ3
     +  * (  - 32 - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2*sQ3**(-1)
     +  * ( 32 + 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ2**(-1) * ( 
     +    64 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ3**(-1) * ( 
     +    64 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)*sQ3
     +  * ( 32 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)
     +  * (  - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2*sQ3**(-1)
     +  * (  - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1)
     +  * (  - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1 * (  - 32 - 
     +    16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*sQ2**(-1)
     +  * ( 96 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*sQ3**(-1)
     +  * ( 32 + 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ3 * (  - 32 - 
     +    16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13**(-1) * ( 32*m2 + 16*
     +    m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-1) * ( 48*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*s13*sQ3**(-1) * ( 48*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1) * ( 64 + 160*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ1*sQ3**(-1) * ( 32 + 112*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ2**(-1) * (  - 64*m2 - 
     +    160*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ2 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ3**(-1) * (  - 128*m2**2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12*sQ3 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12 * (  - 32 - 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*sQ2**(-1)
     +  * ( 32 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*sQ3**(-1)
     +  * ( 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)*sQ3
     +  * ( 32 )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1) * ( 
     +     - 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * (  - 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1) * ( 64*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**2*sQ3**(-1) * ( 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-1) * ( 
     +    32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3**2
     +  * ( 32 - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2*sQ3**(-1)
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2 * (  - 32*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ3 * ( 32 - 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1 * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*sQ3
     +  * ( 32 + 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2*sQ3**(-1)
     +  * (  - 32 - 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ2**(-1) * ( 
     +    64 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ3**(-1) * ( 
     +    64 + 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1) * ( 32 + 112*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-1) * ( 64 + 160*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ2**(-1) * (  - 128*m2**2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ2 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ3**(-1) * (  - 64*m2 - 
     +    160*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13*sQ3 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13 * (  - 32 - 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*D4**(-1)*s13**2*sQ2**(-1) * ( 16*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*s13**2*sQ3**(-1) * ( 64*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1)*sQ3 * (  - 64 )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1) * (  - 256*m2 - 
     +    128*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-1) * (  - 64 )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1*sQ3**(-1) * (  - 256*m2 - 
     +    128*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1 * (  - 128 )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-1) * ( 128 + 128*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ1**2*sQ3**(-1) * ( 128 + 128*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3 * (  - 32*m2 + 48
     +    *m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ2*sQ3**(-1) * (  - 32*m2 + 48
     +    *m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ2*sQ3 * (  - 128*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ2 * ( 32 - 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ2**2 * (  - 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ3 * ( 32 - 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1)*sQ3**2 * (  - 64*q2**(-1) )
     +
      W = W + D2**(-1)*D4**(-1) * (  - 64*m2 + 96*m2**2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1) * (  - 32 + 
     +    32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1) * ( 32 + 64*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1) * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*s13 * (  - 32 - 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*s12**(-1)*s13**2*sQ3**(-1) * ( 32*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1*sQ2**(-1) * (  - 32*m2 - 
     +    16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1*sQ2*sQ3**(-1) * ( 64 + 64*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1*sQ3**(-1) * (  - 32*m2 - 
     +    16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1 * ( 32 + 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1**2*sQ2**(-1) * ( 96 + 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ1**2*sQ3**(-1) * ( 96 + 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**(-1)*sQ2**2*sQ3**(-1) * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1) * ( 32 + 64*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1) * (  - 32 + 
     +    32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3 * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12*s13**(-1) * (  - 32 - 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*s12*sQ2**(-1) * ( 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12*sQ3**(-1) * ( 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s12**2*s13**(-1)*sQ2**(-1) * ( 32*m2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 64 + 64*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1*sQ2**(-1) * (  - 32*m2 - 
     +    16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1*sQ2*sQ3**(-1) * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1*sQ3**(-1) * (  - 32*m2 - 
     +    16*m2**2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1 * ( 32 + 48*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1**2*sQ2**(-1) * ( 96 + 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ1**2*sQ3**(-1) * ( 96 + 48*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13**(-1)*sQ2**(-1)*sQ3**2 * ( 32 + 16*
     +    m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13*sQ2**(-1) * ( 16*m2*q2**(-1) )
     +
      W = W + D2**(-1)*s13*sQ3**(-1) * ( 32*m2*q2**(-1) )
     +
      W = W + D2**(-1)*sQ1*sQ2**(-1) * ( 32 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*sQ1*sQ3**(-1) * ( 32 + 80*m2*q2**(-1) )
     +
      W = W + D2**(-1)*sQ2**(-1)*sQ3 * (  - 32 + 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*sQ2**(-1) * (  - 64*m2 - 96*m2**2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*sQ2*sQ3**(-1) * (  - 32 + 16*m2*q2**(-1)
     +     )
     +
      W = W + D2**(-1)*sQ2 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1)*sQ3**(-1) * (  - 64*m2 - 96*m2**2*
     +    q2**(-1) )
     +
      W = W + D2**(-1)*sQ3 * ( 64*q2**(-1) )
     +
      W = W + D2**(-1) * (  - 64 + 32*m2*q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ2**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ2**(-1)*sQ3**(-1) * ( 256*m2**2 + 
     +    128*m2**3*q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ2**(-1)*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ2**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ2*sQ3**(-1) * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ3**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s12*sQ3**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ2**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ2**(-1)*sQ3**(-1) * ( 256*m2**2 + 
     +    128*m2**3*q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ2**(-1)*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ2**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ2*sQ3**(-1) * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ3**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*s13*sQ3**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ2**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ2**(-1)*sQ3**(-1) * ( 256*m2**2 + 
     +    128*m2**3*q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ2**(-1)*sQ3 * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ2**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ2*sQ3**(-1) * (  - 32 - 16*m2*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ3**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-2)*sQ1*sQ3**(-1) * ( 128*m2 + 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m2 - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-2) * (  - 256*
     +    m2 - 128*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1) * ( 64 - 32*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1) * (  - 16*m2
     +    *q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13*sQ3**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13 * (  - 32 )
     +
      W = W + D4**(-1)*s12**(-1)*s13**2*sQ3**(-2) * (  - 128*m2
     +     - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*s13**2*sQ3**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1*sQ2**(-1)*sQ3 * (  - 32 )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1*sQ2**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1*sQ2*sQ3**(-1) * ( 32 - 32*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1*sQ3**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1 * (  - 32*m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m2 - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1**2*sQ3**(-2) * (  - 128*m2
     +     - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**(-1)*sQ1**2*sQ3**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-2) * (  - 256*
     +    m2 - 128*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m2 - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1) * ( 64 - 32*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3 * (  - 16*m2
     +    *q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1)*sQ2**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*s13**(-1) * (  - 32 )
     +
      W = W + D4**(-1)*s12*sQ2**(-2) * (  - 128*m2 - 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*s12*sQ2**(-1)*sQ3**(-1) * (  - 64*m2 - 
     +    32*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*sQ2**(-1) * (  - 32*m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12*sQ3**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*s12*sQ3**(-1) * (  - 16*m2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**2*s13**(-1)*sQ2**(-2) * (  - 128*m2
     +     - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s12**2*s13**(-1)*sQ2**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 32 - 32*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1*sQ2**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1*sQ2*sQ3**(-1) * (  - 32 )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1*sQ3**(-1) * ( 64*m2 + 32*
     +    m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1 * (  - 32*m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-2) * (  - 128*m2
     +     - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m2 - 64*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13**(-1)*sQ1**2*sQ3**(-1) * ( 32 - 16*
     +    m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13*sQ2**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*s13*sQ2**(-1)*sQ3**(-1) * (  - 64*m2 - 
     +    32*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*s13*sQ2**(-1) * (  - 16*m2*q2**(-1) )
     +
      W = W + D4**(-1)*s13*sQ3**(-2) * (  - 128*m2 - 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*s13*sQ3**(-1) * (  - 32*m2*q2**(-1) )
     +
      W = W + D4**(-1)*sQ1*sQ2**(-2) * (  - 128*m2 - 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ1*sQ2**(-1)*sQ3**(-1) * (  - 256*m2 - 
     +    128*m2**2*q2**(-1) )
     +
      W = W + D4**(-1)*sQ1*sQ2**(-1) * (  - 64*m2*q2**(-1) )
     +
      W = W + D4**(-1)*sQ1*sQ3**(-2) * (  - 128*m2 - 64*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ1*sQ3**(-1) * (  - 64*m2*q2**(-1) )
     +
      W = W + D4**(-1)*sQ2**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ2**(-1)*sQ3**(-1) * ( 256*m2**2 + 128*
     +    m2**3*q2**(-1) )
     +
      W = W + D4**(-1)*sQ2**(-1)*sQ3 * (  - 32 - 32*m2*q2**(-1)
     +     )
     +
      W = W + D4**(-1)*sQ2**(-1) * ( 128*m2 + 128*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ2*sQ3**(-1) * (  - 32 - 32*m2*q2**(-1)
     +     )
     +
      W = W + D4**(-1)*sQ2 * (  - 64*q2**(-1) )
     +
      W = W + D4**(-1)*sQ3**(-2) * ( 128*m2**2 + 64*m2**3*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ3**(-1) * ( 128*m2 + 128*m2**2*
     +    q2**(-1) )
     +
      W = W + D4**(-1)*sQ3 * (  - 64*q2**(-1) )
     +
      W = W + D4**(-1) * (  - 32*m2*q2**(-1) )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3**(-1) * ( 
     +     - 16*m2*q2 - 16*m2**2 + 32*m2**3*q2**(-1) )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 32 )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1*sQ2**(-1) * ( 32*m2 + 16*
     +    m2**2*q2**(-1) )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1*sQ2*sQ3**(-1) * ( 32 )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1*sQ3**(-1) * ( 32*m2 + 16*
     +    m2**2*q2**(-1) )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1 * ( 64 )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1)
     +  * ( 32*m2 + 16*m2**2*q2**(-1) )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1**2*sQ2**(-1) * ( 64 )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1**2*sQ3**(-1) * ( 64 )
     +
      W = W + s12**(-1)*s13**(-1)*sQ1**3*sQ2**(-1)*sQ3**(-1)
     +  * ( 64 + 16*m2*q2**(-1) )
     +
      W = W + s12**(-1)*s13*sQ1*sQ2**(-1)*sQ3**(-1) * ( 32 )
     +
      W = W + s12**(-1)*s13*sQ3**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + s12**(-1)*s13*sQ3**(-1) * ( 32 )
     +
      W = W + s12**(-1)*sQ1*sQ2**(-1)*sQ3**(-1) * (  - 96*m2 - 
     +    48*m2**2*q2**(-1) )
     +
      W = W + s12**(-1)*sQ1*sQ2**(-1) * ( 64 + 16*m2*q2**(-1) )
     +
      W = W + s12**(-1)*sQ1*sQ3**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + s12**(-1)*sQ1*sQ3**(-1) * ( 64 )
     +
      W = W + s12**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1) * ( 64 )
     +
      W = W + s12**(-1)*sQ2*sQ3**(-1) * ( 32 )
     +
      W = W + s12**(-1)*sQ3**(-2) * (  - 64*m2*q2 + 32*m2**2 + 
     +    32*m2**3*q2**(-1) )
     +
      W = W + s12**(-1)*sQ3**(-1) * ( 16*m2**2*q2**(-1) + 32*q2
     +     )
     +
      W = W + s12**(-1) * (  - 32 )
     +
      W = W + s12*s13**(-1)*sQ1*sQ2**(-1)*sQ3**(-1) * ( 32 )
     +
      W = W + s12*s13**(-1)*sQ2**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + s12*s13**(-1)*sQ2**(-1) * ( 32 )
     +
      W = W + s13**(-1)*sQ1*sQ2**(-2) * (  - 64*m2 - 32*m2**2*
     +    q2**(-1) )
     +
      W = W + s13**(-1)*sQ1*sQ2**(-1)*sQ3**(-1) * (  - 96*m2 - 
     +    48*m2**2*q2**(-1) )
     +
      W = W + s13**(-1)*sQ1*sQ2**(-1) * ( 64 )
     +
      W = W + s13**(-1)*sQ1*sQ3**(-1) * ( 64 + 16*m2*q2**(-1) )
     +
      W = W + s13**(-1)*sQ1**2*sQ2**(-1)*sQ3**(-1) * ( 64 )
     +
      W = W + s13**(-1)*sQ2**(-2) * (  - 64*m2*q2 + 32*m2**2 + 
     +    32*m2**3*q2**(-1) )
     +
      W = W + s13**(-1)*sQ2**(-1)*sQ3 * ( 32 )
     +
      W = W + s13**(-1)*sQ2**(-1) * ( 16*m2**2*q2**(-1) + 32*q2
     +     )
     +
      W = W + s13**(-1) * (  - 32 )
     +
      W = W + sQ1*sQ2**(-1)*sQ3**(-1) * ( 64 )
     +
      W = W + sQ2**(-2) * ( 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + sQ2**(-1)*sQ3**(-1) * (  - 128*m2 )
     +
      W = W + sQ2**(-1) * (  - 32 )
     +
      W = W + sQ3**(-2) * ( 64*m2 - 32*m2**2*q2**(-1) )
     +
      W = W + sQ3**(-1) * (  - 32 )


      WQqgg_2 =  W

      return
      end


*---------------------------------------------------------------
      real*8 function WQqqq_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqqq
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C WQqqq1 corresponds to the g graph squared (two fermion loops):
C WQqqq2 corresponds to the interference graph
C Resonant diagrams are removed

      WQqqq =
     +  + D2**(-2) * (  - 32*m2*s12*s23**(-1) - 32*m2*s12**2*s23**(-2)
     +     - 16*m2 - 32*m2**2*q2**(-1)*s12*s23**(-1) - 32*m2**2*
     +    q2**(-1)*s12**2*s23**(-2) - 16*m2**2*q2**(-1) + 64*q2*s12*
     +    s23**(-1) + 64*q2*s12**2*s23**(-2) + 32*q2 )
      WQqqq = WQqqq + D2**(-1) * (  - 32 - 16*m2*q2**(-1)*s12*
     +    s23**(-2)*sQ2 + 16*m2*q2**(-1)*s12*s23**(-2)*sQ3 + 32*m2*
     +    q2**(-1)*s12*s23**(-1) + 32*m2*q2**(-1)*s12**2*s23**(-2) - 16
     +    *m2*q2**(-1)*s23**(-1)*sQ2 + 16*m2*q2**(-1) + 32*m2*s12*
     +    s23**(-2) + 16*m2*s23**(-1) + 32*m2**2*q2**(-1)*s12*s23**(-2)
     +     + 16*m2**2*q2**(-1)*s23**(-1) - 64*q2*s12*s23**(-2) - 32*q2*
     +    s23**(-1) - 32*s12*s23**(-2)*sQ2 + 32*s12*s23**(-2)*sQ3 - 64*
     +    s12*s23**(-1) - 64*s12**2*s23**(-2) - 32*s23**(-1)*sQ2 )
      WQqqq = WQqqq - 32*m2*q2**(-1)*s12*s23**(-2) + 16*m2*q2**(-1)*
     +    s23**(-2)*sQ2 - 16*m2*q2**(-1)*s23**(-1) - 16*q2**(-1)*s12*
     +    s23**(-2)*sQ3 - 16*q2**(-1)*s23**(-2)*s13*sQ2 + 64*s12*
     +    s23**(-2) + 32*s23**(-2)*sQ2 + 32*s23**(-1)


      WQqqq_1_nores = WQqqq

      return
      end


*---------------------------------------------------------------
      real*8 function WQqqq_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqqq2
      real*8 D2, D4, D4p
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2

      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3
      D4p = sQ1+sQ3+s13


C WQqqq1 corresponds to the g graph squared (two fermion loops):
C WQqqq2 corresponds to the interference graph
C Resonant diagrams are removed

      WQqqq2 =
     +  + D2**(-2) * ( 32*m2*s12**2*s23**(-1)*s13**(-1) + 32*m2**2*
     +    q2**(-1)*s12**2*s23**(-1)*s13**(-1) - 64*q2*s12**2*s23**(-1)*
     +    s13**(-1) )
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1) * (  - 32*m2*q2**(-1)*s12*
     +    s13**(-1)*sQ3 + 16*m2*q2**(-1)*s12**2*s13**(-1) + 16*m2*
     +    q2**(-1)*s13**(-1)*sQ3**2 - 48*m2*q2*s13**(-1) + 32*m2*s12*
     +    s23**(-1)*s13**(-1)*sQ3 + 32*m2*s13**(-1)*sQ3 + 32*m2**2*
     +    q2**(-1)*s12*s23**(-1)*s13**(-1)*sQ3 - 32*m2**2*q2**(-1)*s12*
     +    s13**(-1) + 32*m2**2*q2**(-1)*s13**(-1)*sQ3 + 16*m2**3*
     +    q2**(-1)*s13**(-1) - 64*q2*s12*s23**(-1)*s13**(-1)*sQ3 - 64*
     +    q2*s12*s13**(-1) - 64*q2*s13**(-1)*sQ3 + 32*q2**2*s13**(-1)
     +     + 32*s12**2*s13**(-1) )
      WQqqq2 = WQqqq2 + D2**(-1)*D4**(-1) * ( 32*s13**(-1)*sQ3**2 )
      WQqqq2 = WQqqq2 + D2**(-1) * (  - 16*m2*q2**(-1)*s12*s23**(-1)*
     +    s13**(-1)*sQ3 - 16*m2*q2**(-1)*s12**2*s23**(-1)*s13**(-1) - 
     +    48*m2*s12*s23**(-1)*s13**(-1) - 48*m2**2*q2**(-1)*s12*
     +    s23**(-1)*s13**(-1) + 96*q2*s12*s23**(-1)*s13**(-1) - 32*s12*
     +    s23**(-1)*s13**(-1)*sQ3 + 32*s12**2*s23**(-1)*s13**(-1) )
      WQqqq2 = WQqqq2 + D4**(-1) * (  - 16*m2*q2**(-1)*s12*s23**(-1)*
     +    s13**(-1)*sQ3 - 16*m2*q2**(-1)*s23**(-1)*s13**(-1)*sQ3**2 - 
     +    16*m2*s23**(-1)*s13**(-1)*sQ3 - 16*m2**2*q2**(-1)*s23**(-1)*
     +    s13**(-1)*sQ3 + 32*q2*s23**(-1)*s13**(-1)*sQ3 + 32*s12*
     +    s23**(-1)*s13**(-1)*sQ3 - 32*s23**(-1)*s13**(-1)*sQ3**2 )
      WQqqq2 = WQqqq2 + 16*m2*q2**(-1)*s12*s23**(-1)*s13**(-1) + 16*m2*
     +    q2**(-1)*s23**(-1)*s13**(-1)*sQ3 + 16*m2*s23**(-1)*s13**(-1)
     +     + 16*m2**2*q2**(-1)*s23**(-1)*s13**(-1) - 32*q2*s23**(-1)*
     +    s13**(-1) - 32*s12*s23**(-1)*s13**(-1) + 32*s23**(-1)*
     +    s13**(-1)*sQ3

      WQqqq_2_nores = WQqqq2

      return
      end


*---------------------------------------------------------------
      real*8 function WQqgg_1_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqgg
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2


      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      WQqgg = 0d0


      WQqgg =
     +  + D2**(-2)*s12**(-1)*s23 * (  - 16*m**2 - 16*m**4*q2**(-1) + 32
     +    *q2 )
      WQqgg = WQqgg + D2**(-2)*s12*s23**(-1) * (  - 64*m**2 - 64*m**4*
     +    q2**(-1) + 128*q2 )
      WQqgg = WQqgg + D2**(-2)*s12**2*s23**(-2) * (  - 32*m**2 - 32*
     +    m**4*q2**(-1) + 64*q2 )
      WQqgg = WQqgg + D2**(-2) * (  - 48*m**2 - 48*m**4*q2**(-1) + 96*
     +    q2 )
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*s23*sQ3**(-1) * (  - 64*m**2
     +     - 32*m**4*q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*s23 * (  - 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*sQ2*sQ3**(-1) * (  - 16*m**2
     +     - 16*m**4*q2**(-1) + 32*q2 )
      WQqgg = WQqgg + D2**(-1)*s12**(-1)*sQ3**(-1) * ( 48*m**2*q2 - 16*
     +    m**6*q2**(-1) - 32*q2**2 )
      WQqgg = WQqgg + D2**(-1)*s12**(-1) * ( 16*m**2 + 16*m**4*q2**(-1)
     +     - 32*q2 )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-2)*sQ2*sQ3**(-1) * ( 32*m**2
     +     + 32*m**4*q2**(-1) - 64*q2 )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-2)*sQ2**2*sQ3**(-1) * ( 32 + 
     +    16*m**2*q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-2)*sQ3 * ( 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-2) * ( 32*m**2 + 32*m**4*
     +    q2**(-1) - 64*q2 )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-1)*sQ3**(-1) * (  - 112*m**2
     +     - 16*m**4*q2**(-1) + 32*q2 )
      WQqgg = WQqgg + D2**(-1)*s12*s23**(-1) * (  - 96 + 48*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12**2*s23**(-2)*sQ2*sQ3**(-1) * ( 32 - 
     +    16*m**2*q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s12**2*s23**(-2) * (  - 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ3**(-1) * ( 96*m**2*q2 - 32*
     +    m**6*q2**(-1) - 64*q2**2 )
      WQqgg = WQqgg + D2**(-1)*s23**(-1)*sQ3 * ( 32 + 16*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + D2**(-1)*s23**(-1) * ( 32*m**2 + 32*m**4*q2**(-1)
     +     - 64*q2 )
      WQqgg = WQqgg + D2**(-1)*sQ3**(-1) * (  - 128*m**2 - 64*m**4*
     +    q2**(-1) )
      WQqgg = WQqgg + D2**(-1) * (  - 96 + 48*m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-2)*s13*sQ2*sQ3**(-1) * (  - 16*
     +    m**2 - 16*m**4*q2**(-1) + 32*q2 )
      WQqgg = WQqgg + s12**(-1)*s23**(-2)*s13*sQ2 * (  - 16 - 8*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-2)*s13*sQ2**2*sQ3**(-1) * (  - 
     +    16 - 8*m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-2)*s13**2*sQ2*sQ3**(-1) * (  - 
     +    16 + 8*m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-2) * ( 64*m**2
     +     + 32*m**4*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ2*sQ3**(-1) * (  - 32
     +     + 16*m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ3**(-2) * (  - 128*m**2
     +    *q2 + 64*m**4 + 64*m**6*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13*sQ3**(-1) * ( 24*m**2 + 
     +    56*m**4*q2**(-1) - 80*q2 )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13 * (  - 32 + 64*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*s13**2*sQ3**(-1) * ( 48 + 40*
     +    m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2*sQ3**(-1) * ( 24*m**2 + 
     +    24*m**4*q2**(-1) - 48*q2 )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2 * ( 32 + 16*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ2**2*sQ3**(-1) * ( 16 + 8*
     +    m**2*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ3**(-1) * (  - 144*m**2*q2
     +     + 48*m**6*q2**(-1) + 96*q2**2 )
      WQqgg = WQqgg + s12**(-1)*s23**(-1)*sQ3 * ( 48 + 24*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12**(-1)*s23**(-1) * ( 40*m**2 + 40*m**4*
     +    q2**(-1) - 80*q2 )
      WQqgg = WQqgg + s12**(-1)*s23*sQ3**(-1) * ( 32 + 32*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12**(-1)*s13*sQ3**(-2) * ( 64*m**2 + 32*m**4*
     +    q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*s13*sQ3**(-1) * ( 48 + 56*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12**(-1)*sQ2*sQ3**(-1) * ( 16 + 24*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12**(-1)*sQ3**(-2) * (  - 64*m**2*q2 + 32*m**4
     +     + 32*m**6*q2**(-1) )
      WQqgg = WQqgg + s12**(-1)*sQ3**(-1) * ( 80*m**2 + 80*m**4*
     +    q2**(-1) - 64*q2 )
      WQqgg = WQqgg + s12**(-1) * ( 48 + 40*m**2*q2**(-1) )
      WQqgg = WQqgg + s12*s23**(-2)*sQ2*sQ3**(-1) * (  - 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s12*s23**(-2)*sQ2 * ( 16*q2**(-1) )
      WQqgg = WQqgg + s12*s23**(-2) * ( 48 - 24*m**2*q2**(-1) )
      WQqgg = WQqgg + s12*s23**(-1)*sQ2*sQ3**(-1) * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s12*s23**(-1)*sQ3**(-1) * ( 32 - 32*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s12*s23**(-1) * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s12*sQ3**(-1) * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*s13*sQ2*sQ3**(-1) * (  - 64 + 32*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*s13*sQ2 * ( 16*q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*s13*sQ2**2*sQ3**(-1) * ( 24*q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*s13*sQ3 * ( 8*q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*sQ2*sQ3**(-1) * (  - 32*m**2 - 32*m**4*
     +    q2**(-1) + 64*q2 )
      WQqgg = WQqgg + s23**(-2)*sQ2 * (  - 32 - 16*m**2*q2**(-1) )
      WQqgg = WQqgg + s23**(-2)*sQ2**2*sQ3**(-1) * (  - 32 - 16*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s23**(-2) * ( 8*m**2 + 8*m**4*q2**(-1) - 16*q2 )
      WQqgg = WQqgg + s23**(-1)*s13*sQ2*sQ3**(-2) * (  - 32*m**2*
     +    q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*s13*sQ2*sQ3**(-1) * ( 16*q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*s13*sQ3**(-2) * ( 128*m**2 - 64*m**4*
     +    q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*s13*sQ3**(-1) * ( 48 - 40*m**2*q2**(-1)
     +     )
      WQqgg = WQqgg + s23**(-1)*s13 * ( 16*q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*sQ2*sQ3**(-1) * ( 24 )
      WQqgg = WQqgg + s23**(-1)*sQ2 * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*sQ2**2*sQ3**(-1) * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s23**(-1)*sQ3**(-1) * ( 208*m**2 + 32*m**4*
     +    q2**(-1) - 96*q2 )
      WQqgg = WQqgg + s23**(-1) * ( 128 - 16*m**2*q2**(-1) )
      WQqgg = WQqgg + s23*sQ3**(-1) * (  - 8*q2**(-1) )
      WQqgg = WQqgg + s13*sQ3**(-2) * (  - 32*m**2*q2**(-1) )
      WQqgg = WQqgg + s13*sQ3**(-1) * ( 8*q2**(-1) )
      WQqgg = WQqgg + sQ2*sQ3**(-1) * (  - 16*q2**(-1) )
      WQqgg = WQqgg + sQ3**(-2) * ( 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg = WQqgg + sQ3**(-1) * ( 72 - 8*m**2*q2**(-1) )
      WQqgg = WQqgg - 8*q2**(-1)


      WQqgg_1_nores = WQqgg


      return
      end


*---------------------------------------------------------------
      real*8 function WQqgg_2_nores(q2,m2,s12,s13,s23,sQ1,sQ2,sQ3)
*---------------------------------------------------------------
      implicit none

      real*8 pi, zeta2
      parameter(pi = 3.14159265359d0)
      parameter(zeta2 = 1.644934066848226d0)

      real*8 WQqgg2
      real*8 D2, D4
      real*8 s12, s13, s23, sQ1, sQ2, sQ3
      real*8 m, m2, q2


      m = dsqrt(m2)
      D2 = s12+s13+s23
      D4 = s23+sQ2+sQ3

      

      WQqgg2 =
     +  + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*sQ2*sQ3 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     + sQ2**2 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**4*
     + sQ2*sQ3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-2)*s13**4*
     + sQ2**2 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1*sQ2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1*sQ3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2*sQ3 * (  - 32 - 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2**2 * (  - 96 - 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ3**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2*sQ3 * ( 112*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2 * (  - 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2**2 * ( 112*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ3 * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23**(-1)*s13**4*
     + sQ2 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ1 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ2 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13*sQ3 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13**2*sQ2 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s23*s13**2 * ( 64*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3 * (  - 
     +    32 - 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ2**2 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13*sQ3**2 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ1 * ( 64
     +     + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ2*sQ3 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ2 * (  - 
     +    96 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ2**2 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**2*sQ3 * ( 64
     +     + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**3*sQ2 * ( 112*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**(-1)*s13**3 * ( 64*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ1*sQ2
     +  * ( 256 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ1*sQ3
     +  * ( 256 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ2*sQ3
     +  * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ2**2 * ( 
     +    64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13*sQ3**2 * ( 
     +    64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13**2*sQ2*sQ3
     +  * ( 320*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13**2*sQ2 * ( 
     +    384*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13**2*sQ2**2
     +  * ( 192*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13**2*sQ3 * ( 
     +    384*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-2)*s13**2*sQ3**2
     +  * ( 128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ1*sQ2
     +  * (  - 64*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ1*sQ3
     +  * (  - 64*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ1 * ( 256
     +     + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ2*sQ3
     +  * ( 544*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ2 * ( 64
     +     + 800*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ2**2 * ( 
     +    272*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ3 * ( 64
     +     + 800*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13*sQ3**2 * ( 
     +    272*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13**2*sQ2 * ( 
     +    192*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13**2*sQ3 * ( 
     +    128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*s13**2 * ( 384*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2 * ( 320
     +     + 160*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ1*sQ3 * ( 320
     +     + 160*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ2*sQ3 * ( 160
     +     + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23**(-1)*sQ2**2 * ( 160
     +     + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ1 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ2 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*s13**(-1)*sQ3 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*sQ1 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*sQ2 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23*sQ3 * ( 128*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s23 * ( 256*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2*sQ3 * (  - 
     +    32 - 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ2**2 * ( 32
     +     + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13**(-1)*sQ3**2 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13*sQ1 * (  - 64*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13*sQ2 * ( 272*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13*sQ3 * ( 272*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*s13 * ( 768*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ1*sQ2 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ1*sQ3 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ1 * ( 320 + 160*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ2*sQ3 * ( 176*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ2 * ( 160 + 336*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ2**2 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ3 * ( 256*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12*sQ3**2 * ( 128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*s13*sQ2*sQ3
     +  * ( 320*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*s13*sQ2 * ( 
     +    384*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*s13*sQ2**2
     +  * ( 128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*s13*sQ3 * ( 
     +    384*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*s13*sQ3**2
     +  * ( 192*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ1*sQ2 * ( 
     +    128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ1*sQ3 * ( 
     +    128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ2**2 * ( 
     +    64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-2)*sQ3**2 * ( 
     +     - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1*sQ2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1*sQ3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2*sQ3 * (  - 32 - 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 96 - 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13*sQ2 * ( 
     +    128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13*sQ3 * ( 
     +    192*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*s13 * ( 384*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ1*sQ2 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ1*sQ3 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ1 * ( 128
     +     + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ2*sQ3 * ( 
     +    384*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ2 * ( 64
     +     + 480*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ2**2 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ3 * (  - 
     +    64 + 416*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23**(-1)*sQ3**2 * ( 
     +    304*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23*s13**(-1)*sQ3 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s23*s13**(-1) * ( 64*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ1 * ( 64
     +     + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ2*sQ3 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ2 * ( 64
     +     + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ3 * (  - 
     +    96 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*s13**(-1)*sQ3**2 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*sQ1 * (  - 32*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*sQ2 * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2*sQ3 * ( 304*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**2 * ( 448*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ2*sQ3 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ3**2 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*sQ2*sQ3 * ( 
     +    160*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*sQ2 * ( 128*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*sQ2**2 * ( 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*sQ3 * ( 128*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-2)*sQ3**2 * ( 
     +    128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ2*sQ3 * ( 112*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ2 * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ3 * (  - 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ3**2 * ( 112*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*sQ2 * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1)*sQ3 * ( 128*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s23**(-1) * ( 128*m**2
     +    *q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s13**(-1)*sQ3 * ( 112*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**3*s13**(-1) * ( 64*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**4*s23**(-2)*s13**(-1)*
     + sQ2*sQ3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**4*s23**(-2)*s13**(-1)*
     + sQ3**2 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s12**4*s23**(-1)*s13**(-1)*
     + sQ3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ1*sQ2 * ( 
     +    128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ1*sQ3 * ( 
     +    128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ2**2 * ( 
     +     - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**2*sQ3**2 * ( 
     +    64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**3*sQ2*sQ3 * ( 
     +    160*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**3*sQ2 * ( 128*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**3*sQ2**2 * ( 
     +    128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**3*sQ3 * ( 128*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-2)*s13**3*sQ3**2 * ( 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2 * ( 320
     +     + 160*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ1*sQ3 * ( 320
     +     + 160*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ2*sQ3 * ( 160
     +     + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13*sQ3**2 * ( 160
     +     + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ1*sQ2 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ1*sQ3 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ1 * ( 128
     +     + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ2*sQ3 * ( 
     +    384*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ2 * (  - 
     +    64 + 416*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ2**2 * ( 
     +    304*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ3 * ( 64
     +     + 480*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**2*sQ3**2 * ( 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**3*sQ2 * ( 128*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**3*sQ3 * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23**(-1)*s13**3 * ( 128*m**2
     +    *q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*s13*sQ1 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*s13*sQ2 * ( 128*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*s13*sQ3 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*s13 * ( 256*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*sQ1 * ( 128 + 64*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*sQ2 * ( 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s23*sQ3 * ( 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ1*sQ2 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ1*sQ3 * (  - 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ1 * ( 320 + 160*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ2*sQ3 * ( 176*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ2 * ( 256*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ2**2 * ( 128*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ3 * ( 160 + 336*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13*sQ3**2 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13**2*sQ1 * (  - 32*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13**2*sQ2 * ( 304*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13**2*sQ3 * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*s13**2 * ( 448*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*sQ1*sQ2 * ( 128 + 64*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*sQ1*sQ3 * ( 128 + 64*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*sQ2*sQ3 * ( 64 + 32*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*sQ2**2 * ( 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-2)*D4**(-1)*sQ3**2 * ( 32 + 16*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     + sQ1*sQ2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     + sQ1*sQ2**2*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     + sQ2*sQ3 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     + sQ2**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**2*
     + sQ2**3*sQ3**(-1) * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     + sQ2*sQ3 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     + sQ2 * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     + sQ2**2*sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-2)*s13**3*
     + sQ2**2 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     + sQ2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     + sQ2**2*sQ3**(-1) * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*sQ1*
     + sQ3 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     + sQ2**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     + sQ2**3*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13*
     + sQ3**2 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1*sQ2*sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1*sQ2*sQ3**(-1) * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1*sQ3**(-1) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ1 * ( 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2*sQ3**(-1) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2*sQ3 * (  - 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2 * ( 192 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2**2*sQ3**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2**2*sQ3**(-1) * ( 128 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ2**2 * (  - 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2*
     + sQ3 * (  - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**2
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2*sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2*sQ3**(-1) * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ2 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3*
     + sQ3**(-1) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**(-1)*s13**3
     +  * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**2 * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     + sQ2 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     + sQ2**2*sQ3**(-1) * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*sQ1*
     + sQ3 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     + sQ1**2*sQ2**(-1)*sQ3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     + sQ1**2*sQ2*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     + sQ1**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     + sQ1**3*sQ2**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**(-1)*
     + sQ1**3*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     + sQ2**(-1)*sQ3**(-1) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     + sQ2**(-1) * ( 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     + sQ3**(-2) * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1*
     + sQ3**(-1) * ( 128 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ1**2*
     + sQ2**(-1)*sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ2*
     + sQ3**(-2) * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ2*
     + sQ3**(-1) * ( 224 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13*sQ3**(-1)
     +  * (  - 224*m**2 - 112*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13 * ( 64 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     + sQ3**(-2) * (  - 256*m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*s13**2*
     + sQ3**(-1) * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)*
     + sQ3 * ( 96 + 112*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2**(-1)
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     + sQ3**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ2*
     + sQ3**(-1) * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1*sQ3**(-1)
     +  * (  - 160*m**2 - 80*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1 * ( 160 + 
     +    144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     + sQ2**(-1) * ( 128 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**2*
     + sQ3**(-1) * ( 64 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ1**3*
     + sQ2**(-1)*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2*sQ3**(-1)
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2 * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     + sQ3**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23*sQ2**2*
     + sQ3**(-1) * ( 96 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23 * (  - 32*m**2
     +     - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1*sQ2**(-1)*sQ3 * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1*sQ2**(-1) * (  - 32*m**2 - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1*sQ2*sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1*sQ3**(-1) * (  - 32*m**2 - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1 * ( 64 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1**2*sQ2**(-1) * ( 64 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1**2*sQ3**(-1) * ( 64 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13**(-1)*
     + sQ1**3*sQ2**(-1)*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*
     + sQ3**(-2) * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*s13*
     + sQ3**(-1) * ( 96 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     + sQ2**(-1)*sQ3**(-1) * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     + sQ2**(-1) * ( 32 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     + sQ3**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1*
     + sQ3**(-1) * ( 64 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ1**2*
     + sQ2**(-1)*sQ3**(-1) * ( 32 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ2*
     + sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ2*
     + sQ3**(-1) * ( 96 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2*sQ3**(-1)
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**2 * ( 32 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     + sQ1*sQ2**(-1)*sQ3**(-1) * (  - 32*m**2 - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     + sQ1*sQ2**(-1) * ( 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     + sQ1*sQ3**(-1) * ( 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*s13**(-1)*
     + sQ1**2*sQ2**(-1)*sQ3**(-1) * ( 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ3**(-2)
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**3*sQ3**(-1)
     +  * ( 32 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s23**4*s13**(-1)*
     + sQ1*sQ2**(-1)*sQ3**(-1) * ( 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)*
     + sQ3 * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2**(-1)
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     + sQ3**(-2) * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ2*
     + sQ3**(-1) * ( 256 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1*sQ3**(-1)
     +  * (  - 256*m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1 * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     + sQ2**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ1**2*
     + sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2*sQ3**(-1)
     +  * (  - 160*m**2 - 80*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2 * ( 160 + 
     +    48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     + sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ2**2*
     + sQ3**(-1) * ( 192 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13*sQ3 * (  - 96
     +     - 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13 * (  - 32*m**2
     +     - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     + sQ2**(-1) * ( 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     + sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ1*
     + sQ3**(-1) * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     + sQ3**(-2) * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2*
     + sQ3**(-1) * ( 224 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ2 * (  - 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2*sQ3**(-1)
     +  * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**2 * ( 128 - 64
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-2)
     +  * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*s13**3*sQ3**(-1)
     +  * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**(-1)*
     + sQ3**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2 * ( 128 + 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ2**2*
     + sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1*sQ3 * ( 128 + 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2**(-1)*
     + sQ3 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2*sQ2*
     + sQ3**(-1) * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**2 * ( 128 + 64
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ2**(-1)
     +  * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ1**3*sQ3**(-1)
     +  * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**2 * ( 32 + 16*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**(-1)*sQ2**3*sQ3**(-1)
     +  * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**(-1)*
     + sQ3**2 * (  - 32 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**(-1)*
     + sQ3**3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2 * (  - 
     +    32 - 176*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**2*
     + sQ3**(-1) * (  - 32 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**2 * ( 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ2**3*
     + sQ3**(-1) * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ3 * (  - 
     +    32 - 176*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*s13*sQ3**2 * ( 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ2**(-1)*
     + sQ3**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ2 * (  - 
     +    128 - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ1*sQ3 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2**(-1)*
     + sQ3**3 * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2*sQ3 * (  - 
     +    160 - 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-2)*sQ2**2 * (  - 
     +    128 - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     + sQ2 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*sQ1*
     + sQ3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ2**(-1)*sQ3**3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ2**2 * (  - 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13**(-1)*
     + sQ3**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-2)*
     + sQ3 * (  - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-2)*
     + sQ3**2 * (  - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-1)*
     + sQ3 * (  - 64 - 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-1)*
     + sQ3**2 * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**(-1)
     +  * (  - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2*
     + sQ3**(-2) * (  - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2*
     + sQ3**(-1) * (  - 64 - 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**2*
     + sQ3**(-2) * (  - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ2**2*
     + sQ3**(-1) * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ3**(-1)
     +  * (  - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13*sQ3 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*s13 * (  - 64
     +     - 544*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-2)*
     + sQ3 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-1)*
     + sQ3 * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1*sQ2**(-1)
     +  * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ1 * (  - 64
     +     - 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**(-2)*
     + sQ3**2 * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**(-1)*
     + sQ3**2 * ( 32 + 144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**(-1)*
     + sQ3**3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2*sQ3**(-1)
     +  * ( 128*m**2 + 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2*sQ3 * (  - 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2 * (  - 96
     +     - 336*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ2**2 * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3 * (  - 160*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1)*sQ3**2 * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**(-1) * ( 192*m**2 + 
     +    96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     + sQ2**(-2) * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**(-1) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     + sQ2**(-1) * ( 128 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1*
     + sQ3**(-1) * ( 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ1**2*
     + sQ2**(-1)*sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-2)*
     + sQ3 * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-1)*
     + sQ3 * ( 224 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1)*sQ2**(-1)
     +  * (  - 224*m**2 - 112*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13**(-1) * ( 64 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13*sQ2**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13*sQ2**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13*sQ3**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*s13*sQ3**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ1*sQ2**(-1)*
     + sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ1*sQ2**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ1*sQ3**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-2)*sQ3 * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-2) * (  - 128*
     +    m**2 - 256*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m**2 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1)*sQ3 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ2**(-1) * ( 32 - 
     +    144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ3**(-2) * (  - 64*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23*sQ3**(-1) * (  - 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*
     + sQ2**(-2) * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*s13**(-1)*
     + sQ2**(-1) * ( 96 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*sQ1*sQ2**(-1)*
     + sQ3**(-1) * (  - 16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*sQ2**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s23**2*sQ2**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-2)*
     + sQ3 * (  - 192*m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)*
     + sQ3 * ( 256 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2**(-1)
     +  * (  - 256*m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ2*
     + sQ3**(-1) * ( 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1*sQ3**(-1)
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1 * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     + sQ2**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ1**2*
     + sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-2)*
     + sQ3**2 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*sQ3
     +  * (  - 160*m**2 - 80*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2**(-1)*
     + sQ3**2 * ( 192 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ2 * (  - 96
     +     - 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1)*sQ3 * ( 160 + 
     +    48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13**(-1) * (  - 32*m**2
     +     - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-2)*sQ3 * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-2) * (  - 128*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-1)*sQ3 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2**(-1) * (  - 32
     +     - 176*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2*sQ3**(-2) * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ2*sQ3**(-1) * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ3**(-2) * (  - 128*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13*sQ3**(-1) * (  - 32
     +     - 176*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*s13 * ( 96*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-2) * (  - 128*
     +    m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1)*sQ3 * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2**(-1) * ( 96 + 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ2*sQ3**(-1) * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1*sQ3**(-1) * ( 32 + 48
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ1 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-2)*sQ3 * (  - 192*
     +    m**2 - 288*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-2)*sQ3**2 * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1)*sQ3 * ( 32 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1)*sQ3**2 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2**(-1) * (  - 192*m**2
     +     - 224*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2*sQ3**(-2) * (  - 64*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12*sQ2*sQ3**(-1) * (  - 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12 * ( 64 - 480*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     + sQ1*sQ2**(-1)*sQ3**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     + sQ1*sQ3 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     + sQ2**(-1)*sQ3**3 * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     + sQ2*sQ3 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*s13**(-1)*
     + sQ3**2 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*sQ2**(-1)*
     + sQ3**2 * ( 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*sQ2**(-1)*
     + sQ3**3 * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*sQ2*sQ3 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-2)*sQ2 * (  - 
     +    128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1*sQ2**(-2)*sQ3 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1*sQ2**(-1)*sQ3 * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1*sQ2**(-1) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ1 * ( 64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2**(-2)*sQ3**2 * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2**(-1)*sQ3 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2**(-1)*sQ3**2 * ( 128 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2*sQ3 * (  - 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ2 * (  - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3 * ( 192 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)*
     + sQ3**2 * (  - 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*s13**(-1)
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-2)*
     + sQ3 * (  - 128*m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-2)*
     + sQ3**2 * (  - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-1)*
     + sQ3 * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-1)*
     + sQ3**2 * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2**(-1)
     +  * (  - 128*m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2*
     + sQ3**(-1) * (  - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ2 * ( 16*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1)*sQ3 * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23**(-1) * ( 32 - 272
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     + sQ2**(-2) * (  - 256*m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*s13**(-1)*
     + sQ2**(-1) * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*sQ2**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s23*sQ2**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     + sQ2**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     + sQ2**(-1) * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ1*
     + sQ3**(-1) * ( 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-2)*
     + sQ3 * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)*
     + sQ3 * ( 224 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ2**(-1)
     +  * (  - 320*m**2 - 160*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1)*sQ3 * (  - 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*s13**(-1) * ( 128 - 64
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-2)*sQ3 * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-2) * (  - 128*
     +    m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1)*sQ3 * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ2**(-1) * ( 32 - 112
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2*sQ3**(-1) * (  - 64*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**2 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ2**(-1)*sQ3**2 * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ2*sQ3 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ3 * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-2)*s13**(-1)*
     + sQ3**2 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ2**(-2)*sQ3 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ2**(-1)*sQ3 * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ2**(-1) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)*
     + sQ3 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s23**(-1)*s13**(-1)
     +  * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-2)
     +  * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s12**3*s13**(-1)*sQ2**(-1)
     +  * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ2 * (  - 
     +    64 - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ2**2*
     + sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ1*sQ3 * (  - 
     +    128 - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ2*sQ3 * (  - 
     +    160 - 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ2**3*
     + sQ3**(-1) * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13*sQ3**2 * (  - 
     +    128 - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13**2*sQ2*sQ3 * ( 
     +     - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13**2*sQ2**2*
     + sQ3**(-1) * ( 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13**2*sQ2**3*
     + sQ3**(-1) * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-2)*s13**2*sQ3 * (  - 
     +    128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2*
     + sQ3**(-2) * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ2*
     + sQ3**(-1) * ( 128 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1*sQ3**(-1)
     +  * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ1 * (  - 64
     +     - 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**(-1)*sQ3
     +  * ( 128*m**2 + 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2*sQ3 * (  - 
     +    80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2 * (  - 160*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**2*
     + sQ3**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**2*
     + sQ3**(-1) * ( 32 + 144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**2 * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ2**3*
     + sQ3**(-1) * ( 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ3 * (  - 96
     +     - 336*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13*sQ3**2 * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13 * ( 192*m**2 + 
     +    96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2**(-1)*
     + sQ3 * (  - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2*
     + sQ3**(-2) * (  - 128*m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2*
     + sQ3**(-1) * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2 * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2**2*
     + sQ3**(-2) * (  - 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ2**2*
     + sQ3**(-1) * ( 80*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ3**(-1)
     +  * (  - 128*m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2*sQ3 * ( 16*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*s13**2 * ( 32 - 272
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ1*sQ2 * (  - 128
     +     - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ1*sQ3 * (  - 128
     +     - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ2*sQ3 * (  - 256
     +     - 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ2**2 * (  - 128
     +     - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**(-1)*sQ3**2 * (  - 128
     +     - 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-2)*
     + sQ3 * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)*
     + sQ3 * ( 128 + 64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2**(-1)
     +  * (  - 160*m**2 - 80*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ2*
     + sQ3**(-1) * ( 96 + 112*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1*sQ3**(-1)
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1 * ( 160 + 
     +    144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     + sQ2**(-1) * ( 64 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**2*
     + sQ3**(-1) * ( 128 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ1**3*
     + sQ2**(-1)*sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-2)*
     + sQ3**2 * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*sQ3
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ2**(-1)*
     + sQ3**2 * ( 96 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1)*sQ3 * ( 64 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**(-1) * (  - 32*m**2
     +     - 16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ1*sQ2**(-1)*
     + sQ3**(-1) * ( 32 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ1*sQ2**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ1*sQ3**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-2) * (  - 64*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m**2 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2**(-1) * (  - 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2*sQ3**(-2) * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ2*sQ3**(-1) * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ3**(-2) * (  - 128*
     +    m**2 - 256*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13*sQ3**(-1) * ( 32 - 
     +    144*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**2*sQ3**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*s13**2*sQ3**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-2) * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1)*sQ3 * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2**(-1) * ( 128 + 
     +    128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ2*sQ3**(-1) * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-2) * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1*sQ3**(-1) * ( 128 + 
     +    128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1 * (  - 64*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ1**2*sQ2**(-1)*
     + sQ3**(-1) * ( 64 + 96*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-2)*sQ3 * (  - 128*
     +    m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1)*sQ3 * ( 64 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2**(-1) * (  - 96*m**2
     +     - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-2) * (  - 128*
     +    m**2 - 128*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ2*sQ3**(-1) * ( 64 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23*sQ3**(-1) * (  - 96*m**2
     +     - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23 * ( 64 - 64*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     + sQ2**(-2) * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**(-1) * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     + sQ2**(-1) * ( 64 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1*
     + sQ3**(-1) * ( 32 + 128*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ1**2*
     + sQ2**(-1)*sQ3**(-1) * ( 32 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-2)*
     + sQ3 * (  - 128*m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-1)*
     + sQ3 * ( 96 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1)*sQ2**(-1)
     +  * (  - 96*m**2 - 48*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13**(-1) * ( 32 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13*sQ1*sQ2**(-1)*
     + sQ3**(-1) * (  - 16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13*sQ3**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*s13*sQ3**(-1) * ( 32*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ1*sQ2**(-1)*
     + sQ3**(-1) * ( 32 + 80*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ1*sQ2**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ1*sQ3**(-1) * (  - 
     +    32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-2) * (  - 64*
     +    m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-1)*sQ3**(-1)
     +  * (  - 128*m**2 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ2**(-1) * ( 32 - 32*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ3**(-2) * (  - 64*
     +    m**2 - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**2*sQ3**(-1) * ( 32 - 32*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ1*
     + sQ2**(-1)*sQ3**(-1) * ( 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ2**(-2)
     +  * (  - 64*m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**3*s13**(-1)*sQ2**(-1)
     +  * ( 32 )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s23**3*sQ1*sQ2**(-1)*
     + sQ3**(-1) * (  - 16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**(-1)*
     + sQ3**2 * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2 * ( 128 + 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ2**2*
     + sQ3**(-1) * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1*sQ3 * ( 128 + 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2**(-1)*
     + sQ3 * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2*sQ2*
     + sQ3**(-1) * ( 96 + 48*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**2 * ( 128 + 64
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ2**(-1)
     +  * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ1**3*sQ3**(-1)
     +  * ( 64 + 32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ2**(-1)*sQ3**3
     +  * ( 32 + 16*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**(-1)*sQ3**2 * ( 32 + 16*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1)*sQ3 * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2**(-1) * ( 32 + 48
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ2*sQ3**(-1) * (  - 
     +    16*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-2) * (  - 128*
     +    m**2 - 64*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1*sQ3**(-1) * ( 96 + 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ1 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-2)*sQ3 * (  - 64*
     +    m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**(-1)*sQ3 * (  - 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-2) * (  - 192*
     +    m**2 - 288*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2*sQ3**(-1) * ( 32 - 16
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**2*sQ3**(-2) * (  - 
     +    32*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ2**2*sQ3**(-1) * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13*sQ3**(-1) * (  - 192*m**2
     +     - 224*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13 * ( 64 - 480*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ2**(-1) * (  - 64*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ2*sQ3**(-2) * (  - 
     +    64*m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ2*sQ3**(-1) * ( 80*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-2) * (  - 128*
     +    m**2 - 192*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2*sQ3**(-1) * ( 32 - 112
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*s13**2 * ( 48*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-2)*sQ3 * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1)*sQ3 * ( 96 + 48
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2**(-1) * (  - 192*m**2
     +     - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-2) * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2*sQ3**(-1) * ( 96 + 48
     +    *m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ2 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ3**(-1) * (  - 192*m**2
     +     - 96*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1*sQ3 * (  - 32*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1 * ( 64 + 32*m**2*q2**(-1)
     +     )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ2**(-1) * ( 64 + 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ1**2*sQ3**(-1) * ( 64 + 96*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-2)*sQ3**2 * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3 * ( 32*m**2 + 
     +    16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2**(-1)*sQ3**2 * ( 32 + 16*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2*sQ3**(-1) * ( 32*m**2 + 
     +    16*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2 * (  - 96 - 48*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-2) * (  - 64*
     +    m**2 - 32*m**4*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ2**2*sQ3**(-1) * ( 32 + 16*
     +    m**2*q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1)*sQ3 * (  - 96 - 48*m**2*
     +    q2**(-1) )
      WQqgg2 = WQqgg2 + D2**(-1)*D4**(-1) * ( 192*m**2 + 96*m**4*
     +    q2**(-1) )


      WQqgg_2_nores =  WQqgg2
      return
      end
c
c
c Begin of soft-virtual functions. See f2svwt and finiteWt.m
c
c
      function svfac(s,t,u,m2,mw2,mu2,nf)
      implicit none
      real*8 svfac,s,t,u,m2,mw2,mu2
      integer nf
      real*8 cf,ca,pi
      parameter (cf=4.d0/3.d0)
      parameter (ca=3.d0)
      parameter (pi=3.14159265358979312D0)
      real*8 lnsmbeo2m2,lnspbeo2m2,beta,lnu,lns,lnt,dilg4,dilg3,
     # dilg2,dilg1,lnmu,lnxi,ro,t1,u1,ddilog
      real*8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
c
      ro=2*(m2+mw2)/s-(m2-mw2)**2/s**2
      beta=sqrt(1-ro)
      t1=t-m2
      u1=u-m2
      lnxi=log(xicut)
      lns=log(s/m2)  
      lnt=log(-t1/m2)
      lnu=log(-u1/m2)
      lnmu=log(mu2/m2)
      lnsmbeo2m2=log((m2-mw2+s-beta*s)/(2*m2))
      lnspbeo2m2=log((m2-mw2+s+beta*s)/(2*m2))
      dilg1=ddilog((m2+mw2-s-beta*s-2*u)/(m2-mw2+s-beta*s))
      dilg2=ddilog((m2+mw2-s-beta*s-2*u)/(2*m2-2*u))
      dilg3=ddilog((m2+mw2-s+beta*s-2*u)/
     #             (-m2+mw2+(-1+beta)*s))
      dilg4=ddilog((m2+mw2-s+beta*s-2*u)/(2*mw2-2*s-2*u))
c
      svfac = -cf*(lnsmbeo2m2-lnspbeo2m2)*(s-mw2+m2)/(beta*s)
      svfac = svfac+(-ca*(pi**2+2*lnu**2-4*lnsmbeo2m2*lnu+4*lns*lnu-2*ln
     1   t**2+4*lnsmbeo2m2*lnt-4*lns*lnt-2*lns**2+4*dilg4+4*dilg3-4*dilg
     2   2-4*dilg1)-cf*(pi**2-4*lnu**2+8*lnsmbeo2m2*lnu-8*lns*lnu+2*lnsp
     3   beo2m2**2-4*lnsmbeo2m2*lnspbeo2m2-2*lnsmbeo2m2**2+2*lns**2+4*ln
     4   s+6*lnmu-8*dilg3+8*dilg2))/4.0d0-lnmu*(11*ca-2*nf)/6.0d0+2*(cf+
     5   ca)*lnxi**2-2*((ca-2*cf)*lnu-ca*lnt-ca*lns+(cf+ca)*lnmu+cf)*lnx
     6   i
      return 
      end 


      function svnonfac(s,t,u,m2,mw2,mu2,nf)
      implicit none
      real*8 svnonfac,s,t,u,m2,mw2,mu2
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
      real*8 lam,t1,u1,delq,lns,lnt,lnu,lnq,xl2tot1,xl2uou1,xl2tom,
     # xl2uom,xl2mwom,xl2omdqom,xl2omdqot1,xl2omdqou1,as1f,bs1f,bs2f,
     # bs3f,bs4f,bs5f,cs1f0,cs1f,cs2f0,cs2f,cs3f0,cs3f,cs4f,cs5f,cs6f,
     # cs7f,cs8f,ds1f0,ds1f,ds2f0,ds2f,ds3f0,ds3f,ddilog,lnmu,tmp5,
     # tmp4,tmp3,tmp2,tmp1,tmp0
c
      lam=s**2.+m2**2.+mw2**2.-2.*s*m2-2.*s*mw2-2.*mw2*m2
      t1=t-m2
      u1=u-m2
      delq=mw2-m2
      lns=log(s/m2)
      lnt=log(-t1/m2)
      lnu=log(-u1/m2)
      lnmu=log(mu2/m2)
      lnq=log(-delq/m2)
      xl2tot1=ddilog(t/t1)
      xl2uou1=ddilog(u/u1)
      xl2tom=ddilog(t/m2)
      xl2uom=ddilog(u/m2)
      xl2mwom=ddilog(mw2/m2)
      xl2omdqom=ddilog(1-delq/m2)
      xl2omdqot1=ddilog(1-delq/t1)
      xl2omdqou1=ddilog(1-delq/u1)
c As before, but without the 1/(16*pi**2). Rename it cs4f
      cs4f=-1./sqrt(lam)*(2*ddilog(1.-sqrt(lam)/s)+ddilog((del
     &q**2.+delq*(s+sqrt(lam))-2*s*mw2)/(delq**2.+delq*(s-sqrt(lam))-2*s
     &*mw2))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mw2*(s+sqrt(lam)))/(d
     &elq**2.+delq*(s-sqrt(lam))-2*s*mw2))-2*ddilog((s-sqrt(lam)-delq)/(
     &s+sqrt(lam)-delq))-ddilog((delq**2.+delq*(s+sqrt(lam))-2*mw2*(s+sq
     &rt(lam)))/(delq**2.+delq*(s+sqrt(lam))-2*s*mw2))-1./2.*(log((s-sqr
     &t(lam)-delq)/(s+sqrt(lam)-delq)))**2.+zeta2)
c Other finite parts
      as1f = m2
      BS1F  = 2-LNS
      BS2F  = 2
      BS3F  = 2-LNU*U1/U
      BS4F  = 2-DELQ*LNQ/MW2
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
c
      tmp5 = 8*bs5f*mw2*s/((t-m2)*(u+s-m2))-6*bs4f*mw2*s/((t-m2)*(u+s-m2
     1   ))-2*mw2*s/((t-m2)*(u+s-m2))-2*as1f*s/((t-m2)*(u+s-m2))-14*bs5f
     2   *mw2**2/((t-m2)*(u+s-m2))+12*bs4f*mw2**2/((t-m2)*(u+s-m2))+2*mw
     3   2**2/((t-m2)*(u+s-m2))+8*bs5f*m2*mw2/((t-m2)*(u+s-m2))-6*bs4f*m
     4   2*mw2/((t-m2)*(u+s-m2))-2*m2*mw2/((t-m2)*(u+s-m2))+2*as1f*mw2/(
     5   (t-m2)*(u+s-m2))-2*as1f*m2/((t-m2)*(u+s-m2))-2*bs5f*s/(u+s-m2)+
     6   10*bs5f*mw2/(u+s-m2)-8*bs4f*mw2/(u+s-m2)-2*mw2/(u+s-m2)+2*bs5f*
     7   m2/(u+s-m2)-2*bs4f*m2/(u+s-m2)-2*as1f/(u+s-m2)+2*bs5f*mw2*s/(u+
     8   s-m2)**2-2*bs4f*mw2*s/(u+s-m2)**2-4*bs3f*m2*mw2/((t-m2)*(u-m2))
     9   +ds3f*m2*(t-m2)**3/(mw2*s)+4*ds3f*(t-m2)**3/s+bs5f*m2*(t-m2)**2
     :   /(mw2*s*t)
      tmp5 = tmp5-as1f*(t-m2)**2/(mw2*s*t)+2*bs5f*(t-m2)**2/(s*t)-12*ds3
     1   f*mw2*(t-m2)**2/s+2*ds3f*m2**2*(t-m2)**2/(mw2*s)-cs8f*m2*(t-m2)
     2   **2/(mw2*s)+cs7f*m2*(t-m2)**2/(mw2*s)+cs6f*m2*(t-m2)**2/(mw2*s)
     3   -cs5f*m2*(t-m2)**2/(mw2*s)+7*ds3f*m2*(t-m2)**2/s-4*cs8f*(t-m2)*
     4   *2/s+4*cs7f*(t-m2)**2/s+4*cs6f*(t-m2)**2/s-2*cs5f*(t-m2)**2/s
      tmp5 = 4*bs2f*m2*mw2/((t-m2)*(u-m2))+4*bs3f*m2**2/((t-m2)*(u-m2))+
     1   tmp5-2*cs2f*(t-m2)**2/s+3*ds3f*m2*(t-m2)**2/mw2+ds1f*m2*(t-m2)*
     2   *2/mw2+8*ds3f*(t-m2)**2+2*ds1f*(t-m2)**2-4*bs5f*mw2*(t-m2)/(s*t
     3   )+3*bs5f*m2**2*(t-m2)/(mw2*s*t)-3*as1f*m2*(t-m2)/(mw2*s*t)+6*bs
     4   5f*m2*(t-m2)/(s*t)+bs5f*m2*(t-m2)/(mw2*t)-as1f*(t-m2)/(mw2*t)+2
     5   *bs5f*(t-m2)/t+3*ds3f*m2*s*(t-m2)/mw2+3*ds1f*m2*s*(t-m2)/mw2+6*
     6   ds3f*s*(t-m2)+6*ds1f*s*(t-m2)
      tmp5 = -4*bs2f*m2**2/((t-m2)*(u-m2))+tmp5+12*ds3f*mw2**2*(t-m2)/s-
     1   18*ds3f*m2*mw2*(t-m2)/s+8*cs8f*mw2*(t-m2)/s-8*cs7f*mw2*(t-m2)/s
     2   -12*cs6f*mw2*(t-m2)/s+4*cs5f*mw2*(t-m2)/s-2*cs4f*mw2*(t-m2)/s+6
     3   *cs2f*mw2*(t-m2)/s+3*ds3f*m2**3*(t-m2)/(mw2*s)-cs8f*m2**2*(t-m2
     4   )/(mw2*s)+cs7f*m2**2*(t-m2)/(mw2*s)+2*cs6f*m2**2*(t-m2)/(mw2*s)
     5   -2*cs5f*m2**2*(t-m2)/(mw2*s)+cs4f*m2**2*(t-m2)/(mw2*s)+cs2f*m2*
     6   *2*(t-m2)/(mw2*s)-bs5f*m2*(t-m2)/(mw2*s)-2*bs2f*m2*(t-m2)/(mw2*
     7   s)-m2*(t-m2)/(mw2*s)/2.0d0+2*as1f*(t-m2)/(mw2*s)
      tmp5 = 8*bs3f*m2/(u-m2)+tmp5+3*ds3f*m2**2*(t-m2)/s-4*cs8f*m2*(t-m2
     1   )/s+4*cs7f*m2*(t-m2)/s+7*cs6f*m2*(t-m2)/s-2*cs5f*m2*(t-m2)/s+cs
     2   4f*m2*(t-m2)/s-4*cs2f*m2*(t-m2)/s-2*bs5f*(t-m2)/s-6*bs4f*(t-m2)
     3   /s+4*bs2f*(t-m2)/s-16*ds3f*mw2*(t-m2)-6*ds1f*mw2*(t-m2)+5*ds3f*
     4   m2**2*(t-m2)/mw2+3*ds1f*m2**2*(t-m2)/mw2-3*cs8f*m2*(t-m2)/mw2+2
     5   *cs7f*m2*(t-m2)/mw2+2*cs6f*m2*(t-m2)/mw2-3*cs5f*m2*(t-m2)/mw2-1
     6   9*cs4f*m2*(t-m2)/mw2
      tmp5 = -8*bs2f*m2/(u-m2)+tmp5+cs1f*m2*(t-m2)/mw2-6*bs4f*(t-m2)/mw2
     1   +6*bs3f*(t-m2)/mw2+10*bs2f*(t-m2)/mw2-8*bs1f*(t-m2)/mw2+8*ds3f*
     2   m2*(t-m2)+3*ds1f*m2*(t-m2)-8*cs8f*(t-m2)+4*cs7f*(t-m2)+4*cs6f*(
     3   t-m2)-6*cs5f*(t-m2)+2*cs4f*(t-m2)-10*cs2f*(t-m2)+2*cs1f*(t-m2)-
     4   4*as1f*m2*mw2/(s*t*(t-m2))+2*as1f*m2**2/(s*t*(t-m2))+2*as1f*mw2
     5   /(t*(t-m2))-2*as1f*m2/(t*(t-m2))+ds1f*m2*s**3/(mw2*(t-m2))+4*ds
     6   1f*s**3/(t-m2)-12*ds1f*mw2*s**2/(t-m2)+2*ds1f*m2**2*s**2/(mw2*(
     7   t-m2))-cs8f*m2*s**2/(mw2*(t-m2))-cs5f*m2*s**2/(mw2*(t-m2))
      tmp5 = tmp5+cs4f*m2*s**2/(mw2*(t-m2))+cs1f*m2*s**2/(mw2*(t-m2))-4*
     1   bs4f*s**2/(mw2*(t-m2))+4*bs3f*s**2/(mw2*(t-m2))+7*ds1f*m2*s**2/
     2   (t-m2)-2*cs8f*s**2/(t-m2)-4*cs5f*s**2/(t-m2)+4*cs4f*s**2/(t-m2)
     3   -2*cs2f*s**2/(t-m2)+4*cs1f*s**2/(t-m2)+12*ds1f*mw2**2*s/(t-m2)-
     4   8*ds3f*m2*mw2*s/(t-m2)-18*ds1f*m2*mw2*s/(t-m2)+4*cs8f*mw2*s/(t-
     5   m2)-2*cs6f*mw2*s/(t-m2)+8*cs5f*mw2*s/(t-m2)-12*cs4f*mw2*s/(t-m2
     6   )+6*cs2f*mw2*s/(t-m2)-8*cs1f*mw2*s/(t-m2)+4*ds3f*m2**3*s/(mw2*(
     7   t-m2))+3*ds1f*m2**3*s/(mw2*(t-m2))+cs6f*m2**2*s/(mw2*(t-m2))-cs
     8   5f*m2**2*s/(mw2*(t-m2))+2*cs4f*m2**2*s/(mw2*(t-m2))+cs2f*m2**2*
     9   s/(mw2*(t-m2))+cs1f*m2**2*s/(mw2*(t-m2))-8*bs4f*m2*s/(mw2*(t-m2
     :   ))+8*bs3f*m2*s/(mw2*(t-m2))+2*bs2f*m2*s/(mw2*(t-m2))+(-1.1d1)*m
     ;   2*s/(2.0d0*mw2*(t-m2))+4*ds3f*m2**2*s/(t-m2)+3*ds1f*m2**2*s/(t-
     <   m2)-6*cs8f*m2*s/(t-m2)+cs6f*m2*s/(t-m2)-4*cs5f*m2*s/(t-m2)-cs4f
     =   *m2*s/(t-m2)-4*cs2f*m2*s/(t-m2)+4*cs1f*m2*s/(t-m2)-8*bs5f*s/(t-
     >   m2)
      tmp5 = tmp5-8*bs4f*s/(t-m2)+8*bs3f*s/(t-m2)+4*bs2f*s/(t-m2)+2*s/(t
     1   -m2)-4*cs6f*mw2**3/(s*(t-m2))-4*cs4f*mw2**3/(s*(t-m2))+10*cs6f*
     2   m2*mw2**2/(s*(t-m2))+10*cs4f*m2*mw2**2/(s*(t-m2))-12*bs4f*mw2**
     3   2/(s*(t-m2))+8*bs2f*mw2**2/(s*(t-m2))-6*mw2**2/(s*(t-m2))-6*cs6
     4   f*m2**2*mw2/(s*(t-m2))-6*cs4f*m2**2*mw2/(s*(t-m2))+12*bs4f*m2*m
     5   w2/(s*(t-m2))-14*bs2f*m2*mw2/(s*(t-m2))+21*m2*mw2/(s*(t-m2))+3*
     6   as1f*mw2/(s*(t-m2))+2*cs6f*m2**4/(mw2*s*(t-m2))+2*cs4f*m2**4/(m
     7   w2*s*(t-m2))+2*bs2f*m2**3/(mw2*s*(t-m2))-21*m2**3/(mw2*s*(t-m2)
     8   )+5*as1f*m2**2/(mw2*s*(t-m2))-2*cs6f*m2**3/(s*(t-m2))-2*cs4f*m2
     9   **3/(s*(t-m2))+4*bs2f*m2**2/(s*(t-m2))+6*m2**2/(s*(t-m2))-6*as1
     :   f*m2/(s*(t-m2))-4*ds1f*mw2**3/(t-m2)+8*ds3f*m2*mw2**2/(t-m2)+10
     ;   *ds1f*m2*mw2**2/(t-m2)-4*cs8f*mw2**2/(t-m2)+4*cs6f*mw2**2/(t-m2
     <   )-4*cs5f*mw2**2/(t-m2)+12*cs4f*mw2**2/(t-m2)-4*cs2f*mw2**2/(t-m
     =   2)+4*cs1f*mw2**2/(t-m2)-12*ds3f*m2**2*mw2/(t-m2)-6*ds1f*m2**2*m
     >   w2/(t-m2)+6*cs8f*m2*mw2/(t-m2)
      tmp5 = tmp5-4*cs7f*m2*mw2/(t-m2)-14*cs6f*m2*mw2/(t-m2)+6*cs5f*m2*m
     1   w2/(t-m2)-10*cs4f*m2*mw2/(t-m2)+6*cs2f*m2*mw2/(t-m2)-6*cs1f*m2*
     2   mw2/(t-m2)+12*bs5f*mw2/(t-m2)+16*bs4f*mw2/(t-m2)-16*bs3f*mw2/(t
     3   -m2)-8*bs2f*mw2/(t-m2)+2*mw2/(t-m2)+4*ds3f*m2**4/(mw2*(t-m2))+2
     4   *ds1f*m2**4/(mw2*(t-m2))-2*cs8f*m2**3/(mw2*(t-m2))+2*cs7f*m2**3
     5   /(mw2*(t-m2))+6*cs6f*m2**3/(mw2*(t-m2))-2*cs5f*m2**3/(mw2*(t-m2
     6   ))-55*cs4f*m2**3/(mw2*(t-m2))+cs2f*m2**3/(mw2*(t-m2))+2*cs1f*m2
     7   **3/(mw2*(t-m2))+bs5f*m2**2/(mw2*(t-m2))-4*bs4f*m2**2/(mw2*(t-m
     8   2))+4*bs3f*m2**2/(mw2*(t-m2))+48*bs2f*m2**2/(mw2*(t-m2))-36*bs1
     9   f*m2**2/(mw2*(t-m2))-24*m2**2/(mw2*(t-m2))-as1f*m2/(mw2*(t-m2))
     :   -2*ds1f*m2**3/(t-m2)-8*cs8f*m2**2/(t-m2)+2*cs7f*m2**2/(t-m2)+4*
     ;   cs6f*m2**2/(t-m2)-cs4f*m2**2/(t-m2)-3*cs2f*m2**2/(t-m2)+3*bs5f*
     <   m2/(t-m2)-14*bs4f*m2/(t-m2)+12*bs3f*m2/(t-m2)-2*bs2f*m2/(t-m2)-
     =   4*as1f/(t-m2)+2*cs7f*m2**3*s/(mw2*(t-m2)**2)
      tmp5 = tmp5+bs5f*m2**2*s/(mw2*(t-m2)**2)+2*bs2f*m2**2*s/(mw2*(t-m2
     1   )**2)-4*m2**2*s/(mw2*(t-m2)**2)+4*cs7f*m2**2*s/(t-m2)**2+2*bs5f
     2   *m2*s/(t-m2)**2+4*bs2f*m2*s/(t-m2)**2-8*m2*s/(t-m2)**2+8*bs2f*m
     3   2*mw2**2/(s*(t-m2)**2)-20*m2*mw2**2/(s*(t-m2)**2)+4*as1f*mw2**2
     4   /(s*(t-m2)**2)-12*bs2f*m2**2*mw2/(s*(t-m2)**2)+30*m2**2*mw2/(s*
     5   (t-m2)**2)-6*as1f*m2*mw2/(s*(t-m2)**2)+4*bs2f*m2**4/(mw2*s*(t-m
     6   2)**2)-10*m2**4/(mw2*s*(t-m2)**2)+2*as1f*m2**3/(mw2*s*(t-m2)**2
     7   )+8*cs6f*m2*mw2**2/(t-m2)**2-12*cs6f*m2**2*mw2/(t-m2)**2-4*bs5f
     8   *m2*mw2/(t-m2)**2+12*bs4f*m2*mw2/(t-m2)**2-16*bs2f*m2*mw2/(t-m2
     9   )**2+36*m2*mw2/(t-m2)**2-6*as1f*mw2/(t-m2)**2+4*cs6f*m2**4/(mw2
     :   *(t-m2)**2)+2*bs5f*m2**3/(mw2*(t-m2)**2)+8*bs2f*m2**3/(mw2*(t-m
     ;   2)**2)-41*m2**3/(mw2*(t-m2)**2)+9*as1f*m2**2/(mw2*(t-m2)**2)+2*
     <   bs5f*m2**2/(t-m2)**2+8*bs2f*m2**2/(t-m2)**2-25*m2**2/(t-m2)**2+
     =   3*as1f*m2/(t-m2)**2+2*bs5f*m2**3*s/(mw2*(t-m2)**3)-2*bs2f*m2**3
     >   *s/(mw2*(t-m2)**3)+4*bs5f*m2**2*s/(t-m2)**3-4*bs2f*m2**2*s/(t-m
     ?   2)**3-16*bs2f*m2**2*mw2/(t-m2)**3+40*m2**2*mw2/(t-m2)**3-8*as1f
     @   *m2*mw2/(t-m2)**3
      tmp5 = tmp5+8*bs2f*m2**4/(mw2*(t-m2)**3)-20*m2**4/(mw2*(t-m2)**3)+
     1   4*as1f*m2**3/(mw2*(t-m2)**3)+8*bs2f*m2**3/(t-m2)**3-20*m2**3/(t
     2   -m2)**3+4*as1f*m2**2/(t-m2)**3-4*bs5f*m2*mw2/(s*t)-4*as1f*mw2/(
     3   s*t)+2*bs5f*m2**3/(mw2*s*t)-2*as1f*m2**2/(mw2*s*t)+4*bs5f*m2**2
     4   /(s*t)+2*as1f*m2/(s*t)+2*bs5f*mw2/t+bs5f*m2**2/(mw2*t)-as1f*m2/
     5   (mw2*t)-bs5f*m2/t+as1f/t+ds3f*m2*s**2/mw2+3*ds1f*m2*s**2/mw2+2*
     6   ds3f*s**2+8*ds1f*s**2-6*ds3f*mw2*s-16*ds1f*mw2*s+3*ds3f*m2**2*s
     7   /mw2+5*ds1f*m2**2*s/mw2-3*cs8f*m2*s/mw2+cs7f*m2*s/mw2+cs6f*m2*s
     8   /mw2-3*cs5f*m2*s/mw2+2*cs4f*m2*s/mw2+2*cs1f*m2*s/mw2-10*bs4f*s/
     9   mw2+10*bs3f*s/mw2+3*ds3f*m2*s+8*ds1f*m2*s-6*cs8f*s+2*cs7f*s+2*c
     :   s6f*s-8*cs5f*s
      tmp5 = tmp5+4*cs4f*s-10*cs2f*s+4*cs1f*s-4*ds3f*mw2**3/s+10*ds3f*m2
     1   *mw2**2/s-4*cs8f*mw2**2/s+4*cs7f*mw2**2/s+12*cs6f*mw2**2/s-4*cs
     2   5f*mw2**2/s+4*cs4f*mw2**2/s-4*cs2f*mw2**2/s-6*ds3f*m2**2*mw2/s+
     3   6*cs8f*m2*mw2/s-6*cs7f*m2*mw2/s-18*cs6f*m2*mw2/s+6*cs5f*m2*mw2/
     4   s-6*cs4f*m2*mw2/s+6*cs2f*m2*mw2/s+4*bs5f*mw2/s+12*bs4f*mw2/s-8*
     5   bs2f*mw2/s+4*mw2/s+2*ds3f*m2**4/(mw2*s)-2*cs8f*m2**3/(mw2*s)+2*
     6   cs7f*m2**3/(mw2*s)+3*cs6f*m2**3/(mw2*s)-2*cs5f*m2**3/(mw2*s)+2*
     7   cs4f*m2**3/(mw2*s)+cs2f*m2**3/(mw2*s)-2*bs5f*m2**2/(mw2*s)-2*bs
     8   2f*m2**2/(mw2*s)-10*m2**2/(mw2*s)+5*as1f*m2/(mw2*s)-2*ds3f*m2**
     9   3/s+3*cs6f*m2**2/s-3*cs2f*m2**2/s-4*bs5f*m2/s-2*bs2f*m2/s+6*m2/
     :   s
      tmp5 = tmp5-3*as1f/s+8*ds3f*mw2**2+8*ds1f*mw2**2-20*ds3f*m2*mw2-12
     1   *ds1f*m2*mw2+4*cs8f*mw2-4*cs7f*mw2-8*cs6f*mw2+12*cs5f*mw2-8*cs4
     2   f*mw2+16*cs2f*mw2-4*cs1f*mw2+8*ds3f*m2**3/mw2+4*ds1f*m2**3/mw2-
     3   cs8f*m2**2/mw2+4*cs6f*m2**2/mw2-3*cs5f*m2**2/mw2-68*cs4f*m2**2/
     4   mw2+2*cs2f*m2**2/mw2+2*cs1f*m2**2/mw2-bs5f*m2/mw2-10*bs4f*m2/mw
     5   2+10*bs3f*m2/mw2+50*bs2f*m2/mw2-40*bs1f*m2/mw2-6*m2/mw2-4*as1f/
     6   mw2+4*ds3f*m2**2-10*cs8f*m2+6*cs7f*m2+4*cs6f*m2-6*cs5f*m2+4*cs4
     7   f*m2-12*cs2f*m2+2*cs1f*m2-6*bs5f+12*bs3f-8*bs1f+6
      tmp5 = -cf**2*nc*tmp5/pi**2
      tmp4 = -2*bs5f*mw2/(u+s-m2)+2*bs4f*mw2/(u+s-m2)+2*bs3f*m2*mw2/((t-
     1   m2)*(u-m2))-2*bs2f*m2*mw2/((t-m2)*(u-m2))-ds3f*m2*(t-m2)**3/(mw
     2   2*s)/2.0d0-2*ds3f*(t-m2)**3/s-bs5f*m2*(t-m2)**2/(mw2*s*t)/2.0d0
      tmp4 = tmp4+as1f*(t-m2)**2/(mw2*s*t)/2.0d0+6*ds3f*mw2*(t-m2)**2/s-
     1   ds3f*m2**2*(t-m2)**2/(mw2*s)+cs8f*m2*(t-m2)**2/(mw2*s)/2.0d0-cs
     2   7f*m2*(t-m2)**2/(mw2*s)/2.0d0+cs5f*m2*(t-m2)**2/(mw2*s)/2.0d0-c
     3   s3f*m2*(t-m2)**2/(mw2*s)/2.0d0+(-7.0d0)*ds3f*m2*(t-m2)**2/(2.0d
     4   0*s)+2*cs8f*(t-m2)**2/s-2*cs7f*(t-m2)**2/s-cs6f*(t-m2)**2/s+cs5
     5   f*(t-m2)**2/s-cs3f*(t-m2)**2/s
      tmp4 = -2*bs3f*m2**2/((t-m2)*(u-m2))+2*bs2f*m2**2/((t-m2)*(u-m2))+
     1   tmp4+cs2f*(t-m2)**2/s+(-3.0d0)*ds3f*m2*(t-m2)**2/(2.0d0*mw2)-ds
     2   2f*m2*(t-m2)**2/mw2/2.0d0-ds1f*m2*(t-m2)**2/mw2/2.0d0-4*ds3f*(t
     3   -m2)**2-ds2f*(t-m2)**2-ds1f*(t-m2)**2+(-3.0d0)*bs5f*m2**2*(t-m2
     4   )/(2.0d0*mw2*s*t)+3.0d0*as1f*m2*(t-m2)/(2.0d0*mw2*s*t)-bs5f*m2*
     5   (t-m2)/(mw2*t)/2.0d0+as1f*(t-m2)/(mw2*t)/2.0d0+(-3.0d0)*ds3f*m2
     6   *s*(t-m2)/(2.0d0*mw2)-ds2f*m2*s*(t-m2)/mw2+(-3.0d0)*ds1f*m2*s*(
     7   t-m2)/(2.0d0*mw2)-3*ds3f*s*(t-m2)
      tmp4 = -4*bs3f*m2/(u-m2)+tmp4-3*ds1f*s*(t-m2)-6*ds3f*mw2**2*(t-m2)
     1   /s+9*ds3f*m2*mw2*(t-m2)/s-4*cs8f*mw2*(t-m2)/s+4*cs7f*mw2*(t-m2)
     2   /s+3*cs6f*mw2*(t-m2)/s-2*cs5f*mw2*(t-m2)/s+2*cs3f*mw2*(t-m2)/s-
     3   3*cs2f*mw2*(t-m2)/s+(-3.0d0)*ds3f*m2**3*(t-m2)/(2.0d0*mw2*s)+cs
     4   8f*m2**2*(t-m2)/(mw2*s)/2.0d0-cs7f*m2**2*(t-m2)/(mw2*s)/2.0d0+c
     5   s6f*m2**2*(t-m2)/(mw2*s)/2.0d0+cs5f*m2**2*(t-m2)/(mw2*s)-cs3f*m
     6   2**2*(t-m2)/(mw2*s)-cs2f*m2**2*(t-m2)/(mw2*s)/2.0d0+bs5f*m2*(t-
     7   m2)/(mw2*s)/2.0d0+1.1d1*m2*(t-m2)/(6.0d0*mw2*s)-as1f*(t-m2)/(mw
     8   2*s)/2.0d0
      tmp4 = 4*bs2f*m2/(u-m2)+tmp4+(-3.0d0)*ds3f*m2**2*(t-m2)/(2.0d0*s)+
     1   2*cs8f*m2*(t-m2)/s-2*cs7f*m2*(t-m2)/s-2*cs6f*m2*(t-m2)/s+cs5f*m
     2   2*(t-m2)/s-cs3f*m2*(t-m2)/s+2*cs2f*m2*(t-m2)/s+1.6d1*(t-m2)/(3.
     3   0d0*s)+8*ds3f*mw2*(t-m2)+2*ds2f*mw2*(t-m2)+3*ds1f*mw2*(t-m2)+(-
     4   5.0d0)*ds3f*m2**2*(t-m2)/(2.0d0*mw2)-ds2f*m2**2*(t-m2)/mw2+(-3.
     5   0d0)*ds1f*m2**2*(t-m2)/(2.0d0*mw2)+3.0d0*cs8f*m2*(t-m2)/(2.0d0*
     6   mw2)-cs7f*m2*(t-m2)/mw2+3.0d0*cs5f*m2*(t-m2)/(2.0d0*mw2)-cs3f*m
     7   2*(t-m2)/mw2-cs1f*m2*(t-m2)/mw2
      tmp4 = tmp4+3*bs4f*(t-m2)/mw2-3*bs3f*(t-m2)/mw2-4*ds3f*m2*(t-m2)-d
     1   s2f*m2*(t-m2)+(-3.0d0)*ds1f*m2*(t-m2)/2.0d0+4*cs8f*(t-m2)-2*cs7
     2   f*(t-m2)-2*cs6f*(t-m2)+3*cs5f*(t-m2)+5*cs2f*(t-m2)-2*cs1f*(t-m2
     3   )-ds1f*m2*s**3/(mw2*(t-m2))/2.0d0-2*ds1f*s**3/(t-m2)+6*ds1f*mw2
     4   *s**2/(t-m2)-ds1f*m2**2*s**2/(mw2*(t-m2))+cs8f*m2*s**2/(mw2*(t-
     5   m2))/2.0d0+cs5f*m2*s**2/(mw2*(t-m2))/2.0d0-cs1f*m2*s**2/(mw2*(t
     6   -m2))+2*bs4f*s**2/(mw2*(t-m2))-2*bs3f*s**2/(mw2*(t-m2))+(-7.0d0
     7   )*ds1f*m2*s**2/(2.0d0*(t-m2))+cs8f*s**2/(t-m2)+2*cs5f*s**2/(t-m
     8   2)-cs4f*s**2/(t-m2)+cs2f*s**2/(t-m2)-3*cs1f*s**2/(t-m2)-6*ds1f*
     9   mw2**2*s/(t-m2)+4*ds3f*m2*mw2*s/(t-m2)
      tmp4 = tmp4+4*ds2f*m2*mw2*s/(t-m2)+9*ds1f*m2*mw2*s/(t-m2)-2*cs8f*m
     1   w2*s/(t-m2)-4*cs5f*mw2*s/(t-m2)+3*cs4f*mw2*s/(t-m2)-3*cs2f*mw2*
     2   s/(t-m2)+6*cs1f*mw2*s/(t-m2)-2*ds3f*m2**3*s/(mw2*(t-m2))-2*ds2f
     3   *m2**3*s/(mw2*(t-m2))+(-3.0d0)*ds1f*m2**3*s/(2.0d0*mw2*(t-m2))+
     4   cs5f*m2**2*s/(mw2*(t-m2))/2.0d0+cs4f*m2**2*s/(mw2*(t-m2))/2.0d0
     5   -cs2f*m2**2*s/(mw2*(t-m2))/2.0d0+(-3.0d0)*cs1f*m2**2*s/(2.0d0*m
     6   w2*(t-m2))+4*bs4f*m2*s/(mw2*(t-m2))-4*bs3f*m2*s/(mw2*(t-m2))+1.
     7   1d1*m2*s/(6.0d0*mw2*(t-m2))-2*ds3f*m2**2*s/(t-m2)-2*ds2f*m2**2*
     8   s/(t-m2)+(-3.0d0)*ds1f*m2**2*s/(2.0d0*(t-m2))+3*cs8f*m2*s/(t-m2
     9   )+2*cs5f*m2*s/(t-m2)+2*cs4f*m2*s/(t-m2)+2*cs2f*m2*s/(t-m2)-3*cs
     :   1f*m2*s/(t-m2)+4*bs4f*s/(t-m2)-4*bs3f*s/(t-m2)+1.6d1*s/(3.0d0*(
     ;   t-m2))+2.2d1*mw2**2/(3.0d0*s*(t-m2))+(-4.4d1)*m2*mw2/(3.0d0*s*(
     <   t-m2))+2.2d1*m2**2/(3.0d0*s*(t-m2))+2*ds1f*mw2**3/(t-m2)-4*ds3f
     =   *m2*mw2**2/(t-m2)-5*ds1f*m2*mw2**2/(t-m2)+2*cs8f*mw2**2/(t-m2)+
     >   2*cs5f*mw2**2/(t-m2)-2*cs4f*mw2**2/(t-m2)+2*cs2f*mw2**2/(t-m2)-
     ?   4*cs1f*mw2**2/(t-m2)
      tmp4 = tmp4+6*ds3f*m2**2*mw2/(t-m2)+3*ds1f*m2**2*mw2/(t-m2)-3*cs8f
     1   *m2*mw2/(t-m2)+2*cs7f*m2*mw2/(t-m2)-3*cs5f*m2*mw2/(t-m2)+3*cs4f
     2   *m2*mw2/(t-m2)-3*cs2f*m2*mw2/(t-m2)+6*cs1f*m2*mw2/(t-m2)-8*bs4f
     3   *mw2/(t-m2)+8*bs3f*mw2/(t-m2)+(-1.6d1)*mw2/(3.0d0*(t-m2))-2*ds3
     4   f*m2**4/(mw2*(t-m2))-ds1f*m2**4/(mw2*(t-m2))+cs8f*m2**3/(mw2*(t
     5   -m2))-cs7f*m2**3/(mw2*(t-m2))+cs5f*m2**3/(mw2*(t-m2))+cs4f*m2**
     6   3/(mw2*(t-m2))/2.0d0-cs2f*m2**3/(mw2*(t-m2))/2.0d0-2*cs1f*m2**3
     7   /(mw2*(t-m2))-bs5f*m2**2/(mw2*(t-m2))+2*bs4f*m2**2/(mw2*(t-m2))
     8   -2*bs3f*m2**2/(mw2*(t-m2))+bs2f*m2**2/(mw2*(t-m2))-2*m2**2/(mw2
     9   *(t-m2))+ds1f*m2**3/(t-m2)+4*cs8f*m2**2/(t-m2)-cs7f*m2**2/(t-m2
     :   )+(-3.0d0)*cs4f*m2**2/(2.0d0*(t-m2))+3.0d0*cs2f*m2**2/(2.0d0*(t
     ;   -m2))+2*bs5f*m2/(t-m2)+4*bs4f*m2/(t-m2)-6*bs3f*m2/(t-m2)+2.5d1*
     <   m2/(3.0d0*(t-m2))-cs7f*m2**3*s/(mw2*(t-m2)**2)-bs5f*m2**2*s/(mw
     =   2*(t-m2)**2)+bs2f*m2**2*s/(mw2*(t-m2)**2)-m2**2*s/(mw2*(t-m2)**
     >   2)-2*cs7f*m2**2*s/(t-m2)**2-2*bs5f*m2*s/(t-m2)**2
      tmp4 = tmp4+2*bs2f*m2*s/(t-m2)**2-2*m2*s/(t-m2)**2+2*bs5f*m2*mw2/(
     1   t-m2)**2-2*bs2f*m2*mw2/(t-m2)**2+(-2.2d1)*m2*mw2/(3.0d0*(t-m2)*
     2   *2)-bs5f*m2**3/(mw2*(t-m2)**2)+bs2f*m2**3/(mw2*(t-m2)**2)-bs5f*
     3   m2**2/(t-m2)**2+bs2f*m2**2/(t-m2)**2+2.2d1*m2**2/(3.0d0*(t-m2)*
     4   *2)-bs5f*m2**3*s/(mw2*(t-m2)**3)+bs2f*m2**3*s/(mw2*(t-m2)**3)-2
     5   *bs5f*m2**2*s/(t-m2)**3+2*bs2f*m2**2*s/(t-m2)**3-bs5f*m2**3/(mw
     6   2*s*t)+as1f*m2**2/(mw2*s*t)-bs5f*m2**2/(mw2*t)/2.0d0+as1f*m2/(m
     7   w2*t)/2.0d0-ds3f*m2*s**2/mw2/2.0d0-ds2f*m2*s**2/mw2/2.0d0+(-3.0
     8   d0)*ds1f*m2*s**2/(2.0d0*mw2)-ds3f*s**2-ds2f*s**2-4*ds1f*s**2+3*
     9   ds3f*mw2*s+2*ds2f*mw2*s+8*ds1f*mw2*s+(-3.0d0)*ds3f*m2**2*s/(2.0
     :   d0*mw2)-ds2f*m2**2*s/mw2+(-5.0d0)*ds1f*m2**2*s/(2.0d0*mw2)+3.0d
     ;   0*cs8f*m2*s/(2.0d0*mw2)-cs7f*m2*s/mw2/2.0d0+3.0d0*cs5f*m2*s/(2.
     <   0d0*mw2)-cs3f*m2*s/mw2/2.0d0-2*cs1f*m2*s/mw2+5*bs4f*s/mw2-5*bs3
     =   f*s/mw2+(-3.0d0)*ds3f*m2*s/2.0d0-ds2f*m2*s
      tmp4 = tmp4-4*ds1f*m2*s+3*cs8f*s-cs7f*s+4*cs5f*s-2*cs4f*s-cs3f*s+5
     1   *cs2f*s-2*cs1f*s+2*ds3f*mw2**3/s-5*ds3f*m2*mw2**2/s+2*cs8f*mw2*
     2   *2/s-2*cs7f*mw2**2/s-2*cs6f*mw2**2/s+2*cs5f*mw2**2/s-2*cs3f*mw2
     3   **2/s+2*cs2f*mw2**2/s+3*ds3f*m2**2*mw2/s-3*cs8f*m2*mw2/s+3*cs7f
     4   *m2*mw2/s+3*cs6f*m2*mw2/s-3*cs5f*m2*mw2/s+3*cs3f*m2*mw2/s-3*cs2
     5   f*m2*mw2/s+(-1.6d1)*mw2/(3.0d0*s)-ds3f*m2**4/(mw2*s)+cs8f*m2**3
     6   /(mw2*s)-cs7f*m2**3/(mw2*s)+cs6f*m2**3/(mw2*s)/2.0d0+cs5f*m2**3
     7   /(mw2*s)-cs3f*m2**3/(mw2*s)-cs2f*m2**3/(mw2*s)/2.0d0+bs5f*m2**2
     8   /(mw2*s)-m2**2/(mw2*s)-as1f*m2/(mw2*s)+ds3f*m2**3/s+(-3.0d0)*cs
     9   6f*m2**2/(2.0d0*s)+3.0d0*cs2f*m2**2/(2.0d0*s)+1.9d1*m2/(3.0d0*s
     :   )-4*ds3f*mw2**2
      tmp4 = tmp4-2*ds2f*mw2**2-4*ds1f*mw2**2+10*ds3f*m2*mw2+3*ds2f*m2*m
     1   w2+6*ds1f*m2*mw2-2*cs8f*mw2+2*cs7f*mw2+2*cs6f*mw2-6*cs5f*mw2+2*
     2   cs4f*mw2+2*cs3f*mw2-8*cs2f*mw2+4*cs1f*mw2-4*ds3f*m2**3/mw2-ds2f
     3   *m2**3/mw2-2*ds1f*m2**3/mw2+cs8f*m2**2/mw2/2.0d0+3.0d0*cs5f*m2*
     4   *2/(2.0d0*mw2)+cs4f*m2**2/mw2-cs3f*m2**2/mw2-cs2f*m2**2/mw2-2*c
     5   s1f*m2**2/mw2+bs5f*m2/mw2/2.0d0+5*bs4f*m2/mw2-5*bs3f*m2/mw2+1.1
     6   d1*m2/(3.0d0*mw2)-as1f/mw2/2.0d0-2*ds3f*m2**2+5*cs8f*m2-3*cs7f*
     7   m2-2*cs6f*m2+3*cs5f*m2-cs3f*m2+6*cs2f*m2-2*cs1f*m2+2*bs5f+4*bs4
     8   f-6*bs3f+3.3333333333333335
      tmp4 = -ca*cf*nc*tmp4/pi**2
      tmp4 = tmp5+tmp4
      tmp3 = 6*bs4f*(t-m2)**2/mw2-6*bs3f*(t-m2)**2/mw2+16*bs4f*s*(t-m2)/
     1   mw2-16*bs3f*s*(t-m2)/mw2
      tmp3 = tmp3+8*cs8f*mw2*(t-m2)+10*bs4f*m2*(t-m2)/mw2-10*bs3f*m2*(t-
     1   m2)/mw2+4*bs4f*(t-m2)-4*bs3f*(t-m2)+4*bs4f*s**3/(mw2*(t-m2))-4*
     2   bs3f*s**3/(mw2*(t-m2))+8*bs4f*m2*s**2/(mw2*(t-m2))-8*bs3f*m2*s*
     3   *2/(mw2*(t-m2))+8*bs4f*s**2/(t-m2)-8*bs3f*s**2/(t-m2)-16*bs4f*m
     4   w2*s/(t-m2)+16*bs3f*mw2*s/(t-m2)+4*bs4f*m2**2*s/(mw2*(t-m2))-4*
     5   bs3f*m2**2*s/(mw2*(t-m2))+8*bs4f*m2*s/(t-m2)-8*bs3f*m2*s/(t-m2)
     6   -8*bs4f*m2*mw2/(t-m2)+8*bs3f*m2*mw2/(t-m2)+14*bs4f*s**2/mw2-14*
     7   bs3f*s**2/mw2+8*cs8f*mw2*s+18*bs4f*m2*s/mw2-18*bs3f*m2*s/mw2+12
     8   *bs4f*s-12*bs3f*s+8*cs8f*m2*mw2-16*bs4f*mw2+16*bs3f*mw2+8*mw2+4
     9   *bs4f*m2**2/mw2-4*bs3f*m2**2/mw2+8*bs4f*m2-8*bs3f*m2
      tmp3 = -cf**2*nc*tmp3/(pi**2*(t+s-m2))
      tmp3 = tmp4+tmp3
      tmp2 = -3*bs4f*(t-m2)**2/mw2+3*bs3f*(t-m2)**2/mw2-8*bs4f*s*(t-m2)/
     1   mw2+8*bs3f*s*(t-m2)/mw2
      tmp2 = tmp2-4*cs8f*mw2*(t-m2)-5*bs4f*m2*(t-m2)/mw2+5*bs3f*m2*(t-m2
     1   )/mw2-2*bs4f*(t-m2)+2*bs3f*(t-m2)-2*bs4f*s**3/(mw2*(t-m2))+2*bs
     2   3f*s**3/(mw2*(t-m2))-4*bs4f*m2*s**2/(mw2*(t-m2))+4*bs3f*m2*s**2
     3   /(mw2*(t-m2))-4*bs4f*s**2/(t-m2)+4*bs3f*s**2/(t-m2)+8*bs4f*mw2*
     4   s/(t-m2)-8*bs3f*mw2*s/(t-m2)-2*bs4f*m2**2*s/(mw2*(t-m2))+2*bs3f
     5   *m2**2*s/(mw2*(t-m2))-4*bs4f*m2*s/(t-m2)+4*bs3f*m2*s/(t-m2)+4*b
     6   s4f*m2*mw2/(t-m2)-4*bs3f*m2*mw2/(t-m2)-7*bs4f*s**2/mw2+7*bs3f*s
     7   **2/mw2-4*cs8f*mw2*s-9*bs4f*m2*s/mw2+9*bs3f*m2*s/mw2-6*bs4f*s+6
     8   *bs3f*s-4*cs8f*m2*mw2+8*bs4f*mw2-8*bs3f*mw2-4*mw2-2*bs4f*m2**2/
     9   mw2+2*bs3f*m2**2/mw2-4*bs4f*m2+4*bs3f*m2
      tmp2 = -ca*cf*nc*tmp2/(pi**2*(t+s-m2))
      tmp2 = tmp3+tmp2
      tmp1 = 20*cs4f*m2*s**2*(t-m2)/mw2-10*bs2f*s**2*(t-m2)/mw2+8*bs1f*s
     1   **2*(t-m2)/mw2-38*cs4f*m2**2*s*(t-m2)/mw2+22*bs2f*m2*s*(t-m2)/m
     2   w2-16*bs1f*m2*s*(t-m2)/mw2-2*as1f*s*(t-m2)/mw2-52*cs4f*m2*s*(t-
     3   m2)+6*bs4f*s*(t-m2)+20*bs2f*s*(t-m2)
      tmp1 = tmp1-22*bs1f*s*(t-m2)+2*s*(t-m2)-6*bs4f*m2*mw2*(t-m2)/s+6*b
     1   s2f*m2*mw2*(t-m2)/s+6*bs4f*m2**2*(t-m2)/s-6*bs2f*m2**2*(t-m2)/s
     2   +16*cs4f*m2*mw2*(t-m2)-8*bs4f*mw2*(t-m2)-10*bs2f*mw2*(t-m2)+16*
     3   bs1f*mw2*(t-m2)-2*mw2*(t-m2)+22*cs4f*m2**3*(t-m2)/mw2-16*bs2f*m
     4   2**2*(t-m2)/mw2+12*bs1f*m2**2*(t-m2)/mw2+2*as1f*m2*(t-m2)/mw2-3
     5   8*cs4f*m2**2*(t-m2)-6*bs4f*m2*(t-m2)+28*bs2f*m2*(t-m2)-16*bs1f*
     6   m2*(t-m2)
      tmp1 = tmp1+2*m2*(t-m2)-2*as1f*(t-m2)+4*cs4f*m2*s**3/(t-m2)-4*cs4f
     1   *m2*mw2*s**2/(t-m2)+56*cs4f*m2**3*s**2/(mw2*(t-m2))-42*bs2f*m2*
     2   *2*s**2/(mw2*(t-m2))+36*bs1f*m2**2*s**2/(mw2*(t-m2))+4*as1f*m2*
     3   s**2/(mw2*(t-m2))-8*cs4f*m2**2*s**2/(t-m2)+2*bs4f*m2*s**2/(t-m2
     4   )-4*bs2f*m2*s**2/(t-m2)+2*bs1f*m2*s**2/(t-m2)+2*m2*s**2/(t-m2)-
     5   4*cs4f*m2**2*mw2*s/(t-m2)-12*bs4f*m2*mw2*s/(t-m2)+4*bs2f*m2*mw2
     6   *s/(t-m2)+8*bs1f*m2*mw2*s/(t-m2)-2*m2*mw2*s/(t-m2)-106*cs4f*m2*
     7   *4*s/(mw2*(t-m2))+82*bs2f*m2**3*s/(mw2*(t-m2))-68*bs1f*m2**3*s/
     8   (mw2*(t-m2))-10*as1f*m2**2*s/(mw2*(t-m2))-106*cs4f*m2**3*s/(t-m
     9   2)-2*bs4f*m2**2*s/(t-m2)+92*bs2f*m2**2*s/(t-m2)-72*bs1f*m2**2*s
     :   /(t-m2)+2*m2**2*s/(t-m2)-14*as1f*m2*s/(t-m2)+54*cs4f*m2**3*mw2/
     ;   (t-m2)-6*bs4f*m2**2*mw2/(t-m2)-38*bs2f*m2**2*mw2/(t-m2)+36*bs1f
     <   *m2**2*mw2/(t-m2)+6*as1f*m2*mw2/(t-m2)+54*cs4f*m2**5/(mw2*(t-m2
     =   ))-44*bs2f*m2**4/(mw2*(t-m2))+36*bs1f*m2**4/(mw2*(t-m2))
      tmp1 = tmp1+6*as1f*m2**3/(mw2*(t-m2))-108*cs4f*m2**4/(t-m2)+6*bs4f
     1   *m2**3/(t-m2)+82*bs2f*m2**3/(t-m2)-72*bs1f*m2**3/(t-m2)-12*as1f
     2   *m2**2/(t-m2)+72*cs4f*m2**2*s**2/mw2-48*bs2f*m2*s**2/mw2+40*bs1
     3   f*m2*s**2/mw2+4*as1f*s**2/mw2-4*bs4f*s**2+4*bs1f*s**2+16*cs4f*m
     4   2*mw2*s+4*bs4f*mw2*s-4*bs1f*mw2*s-136*cs4f*m2**3*s/mw2+96*bs2f*
     5   m2**2*s/mw2-76*bs1f*m2**2*s/mw2-12*as1f*m2*s/mw2-144*cs4f*m2**2
     6   *s+24*bs4f*m2*s+96*bs2f*m2*s-100*bs1f*m2*s+8*m2*s-12*as1f*s+72*
     7   cs4f*m2**2*mw2+20*bs4f*m2*mw2-64*bs2f*m2*mw2+32*bs1f*m2*mw2+8*a
     8   s1f*mw2+72*cs4f*m2**4/mw2-56*bs2f*m2**3/mw2+44*bs1f*m2**3/mw2+8
     9   *as1f*m2**2/mw2-144*cs4f*m2**3+4*bs4f*m2**2+96*bs2f*m2**2-76*bs
     :   1f*m2**2-16*as1f*m2
      tmp1 = -cf**2*nc*tmp1/(lam*pi**2)
      tmp1 = tmp2+tmp1
      tmp1 = tmp1-cf**2*nc*(12*cs4f*m2*s**3*(t-m2)-12*cs4f*m2*mw2*s**2*(
     1   t-m2)+12*cs4f*m2**2*s**2*(t-m2)-12*bs4f*m2*s**2*(t-m2)-12*bs2f*
     2   m2*s**2*(t-m2)+24*bs1f*m2*s**2*(t-m2)-12*bs4f*m2*mw2*s*(t-m2)+1
     3   2*bs2f*m2*mw2*s*(t-m2)+12*bs4f*m2**2*s*(t-m2)-12*bs2f*m2**2*s*(
     4   t-m2)+12*cs4f*m2**2*s**4/(t-m2)-12*cs4f*m2**2*mw2*s**3/(t-m2)+1
     5   2*cs4f*m2**3*s**3/(t-m2)-12*bs4f*m2**2*s**3/(t-m2)-12*bs2f*m2**
     6   2*s**3/(t-m2)+24*bs1f*m2**2*s**3/(t-m2)-12*bs4f*m2**2*mw2*s**2/
     7   (t-m2)+12*bs2f*m2**2*mw2*s**2/(t-m2)+12*bs4f*m2**3*s**2/(t-m2)-
     8   12*bs2f*m2**3*s**2/(t-m2)+48*cs4f*m2**2*s**3-24*bs4f*m2*s**3+24
     9   *bs1f*m2*s**3+24*bs4f*m2*mw2*s**2-24*bs1f*m2*mw2*s**2+24*bs4f*m
     :   2**2*s**2-48*bs2f*m2**2*s**2+24*bs1f*m2**2*s**2)/(lam**2*pi**2)
      tmp1 = tmp1+cf**2*lnmu*nc*(3*m2*(t-m2)/(mw2*s)+2*(t-m2)/s+3*m2*s/(
     1   mw2*(t-m2))+2*s/(t-m2)+12*mw2**2/(s*(t-m2))-14*m2*mw2/(s*(t-m2)
     2   )+10*m2**3/(mw2*s*(t-m2))-8*m2**2/(s*(t-m2))-12*mw2/(t-m2)+10*m
     3   2**2/(mw2*(t-m2))+2*m2/(t-m2)-12*m2*mw2/(t-m2)**2+10*m2**3/(mw2
     4   *(t-m2)**2)+2*m2**2/(t-m2)**2-12*mw2/s+10*m2**2/(mw2*s)+2*m2/s+
     5   6*m2/mw2-8)/pi**2/2.0d0
      tmp1 = tmp1+cf**2*nc*(3*m2*(t-m2)/(mw2*s)+2*(t-m2)/s+3*m2*s/(mw2*(
     1   t-m2))+2*s/(t-m2)+12*mw2**2/(s*(t-m2))-14*m2*mw2/(s*(t-m2))+10*
     2   m2**3/(mw2*s*(t-m2))-8*m2**2/(s*(t-m2))-12*mw2/(t-m2)+10*m2**2/
     3   (mw2*(t-m2))+2*m2/(t-m2)-12*m2*mw2/(t-m2)**2+10*m2**3/(mw2*(t-m
     4   2)**2)+2*m2**2/(t-m2)**2-12*mw2/s+10*m2**2/(mw2*s)+2*m2/s+6*m2/
     5   mw2-8)/pi**2/2.0d0
      tmp1 = tmp1+ca*cf*lnmu*lnu*nc*(2*m2*(t-m2)/(mw2*s)+4*(t-m2)/s+2*m2
     1   *s/(mw2*(t-m2))+4*s/(t-m2)+8*mw2**2/(s*(t-m2))-12*m2*mw2/(s*(t-
     2   m2))+4*m2**3/(mw2*s*(t-m2))-8*mw2/(t-m2)+4*m2**2/(mw2*(t-m2))+4
     3   *m2/(t-m2)-8*m2*mw2/(t-m2)**2+4*m2**3/(mw2*(t-m2)**2)+4*m2**2/(
     4   t-m2)**2-8*mw2/s+4*m2**2/(mw2*s)+4*m2/s+4*m2/mw2)/pi**2/2.0d0
      tmp1 = tmp1+ca*cf*lnu*nc*(2*m2*(t-m2)/(mw2*s)+4*(t-m2)/s+2*m2*s/(m
     1   w2*(t-m2))+4*s/(t-m2)+8*mw2**2/(s*(t-m2))-12*m2*mw2/(s*(t-m2))+
     2   4*m2**3/(mw2*s*(t-m2))-8*mw2/(t-m2)+4*m2**2/(mw2*(t-m2))+4*m2/(
     3   t-m2)-8*m2*mw2/(t-m2)**2+4*m2**3/(mw2*(t-m2)**2)+4*m2**2/(t-m2)
     4   **2-8*mw2/s+4*m2**2/(mw2*s)+4*m2/s+4*m2/mw2)/pi**2/2.0d0-ca*cf*
     5   lnmu*nc*(1.1d1*m2*(t-m2)/(6.0d0*mw2*s)+1.1d1*(t-m2)/(3.0d0*s)+1
     6   .1d1*m2*s/(6.0d0*mw2*(t-m2))+1.1d1*s/(3.0d0*(t-m2))+2.2d1*mw2**
     7   2/(3.0d0*s*(t-m2))-11*m2*mw2/(s*(t-m2))+1.1d1*m2**3/(3.0d0*mw2*
     8   s*(t-m2))+(-2.2d1)*mw2/(3.0d0*(t-m2))+1.1d1*m2**2/(3.0d0*mw2*(t
     9   -m2))+1.1d1*m2/(3.0d0*(t-m2))+(-2.2d1)*m2*mw2/(3.0d0*(t-m2)**2)
     :   +1.1d1*m2**3/(3.0d0*mw2*(t-m2)**2)+1.1d1*m2**2/(3.0d0*(t-m2)**2
     ;   )+(-2.2d1)*mw2/(3.0d0*s)+1.1d1*m2**2/(3.0d0*mw2*s)+1.1d1*m2/(3.
     <   0d0*s)+1.1d1*m2/(3.0d0*mw2))/pi**2
      tmp1 = tmp1+ca*cf*lnmu*nc*(5.0d0*m2*(t-m2)/(3.0d0*mw2*s)+(-2.0d0)*
     1   (t-m2)/(3.0d0*s)+5.0d0*m2*s/(3.0d0*mw2*(t-m2))+(-2.0d0)*s/(3.0d
     2   0*(t-m2))+2.0d1*mw2**2/(3.0d0*s*(t-m2))-6*m2*mw2/(s*(t-m2))+2.2
     3   d1*m2**3/(3.0d0*mw2*s*(t-m2))-8*m2**2/(s*(t-m2))+(-2.0d1)*mw2/(
     4   3.0d0*(t-m2))+2.2d1*m2**2/(3.0d0*mw2*(t-m2))+(-2.0d0)*m2/(3.0d0
     5   *(t-m2))+(-2.0d1)*m2*mw2/(3.0d0*(t-m2)**2)+2.2d1*m2**3/(3.0d0*m
     6   w2*(t-m2)**2)+(-2.0d0)*m2**2/(3.0d0*(t-m2)**2)+(-2.0d1)*mw2/(3.
     7   0d0*s)+2.2d1*m2**2/(3.0d0*mw2*s)+(-2.0d0)*m2/(3.0d0*s)+1.0d1*m2
     8   /(3.0d0*mw2)-8)/pi**2/2.0d0
      tmp1 = tmp1+ca*cf*nc*(5.0d0*m2*(t-m2)/(3.0d0*mw2*s)+(-2.0d0)*(t-m2
     1   )/(3.0d0*s)+5.0d0*m2*s/(3.0d0*mw2*(t-m2))+(-2.0d0)*s/(3.0d0*(t-
     2   m2))+2.0d1*mw2**2/(3.0d0*s*(t-m2))-6*m2*mw2/(s*(t-m2))+2.2d1*m2
     3   **3/(3.0d0*mw2*s*(t-m2))-8*m2**2/(s*(t-m2))+(-2.0d1)*mw2/(3.0d0
     4   *(t-m2))+2.2d1*m2**2/(3.0d0*mw2*(t-m2))+(-2.0d0)*m2/(3.0d0*(t-m
     5   2))+(-2.0d1)*m2*mw2/(3.0d0*(t-m2)**2)+2.2d1*m2**3/(3.0d0*mw2*(t
     6   -m2)**2)+(-2.0d0)*m2**2/(3.0d0*(t-m2)**2)+(-2.0d1)*mw2/(3.0d0*s
     7   )+2.2d1*m2**2/(3.0d0*mw2*s)+(-2.0d0)*m2/(3.0d0*s)+1.0d1*m2/(3.0
     8   d0*mw2)-8)/pi**2/2.0d0-ca*cf*lnu*nc*(m2*(t-m2)/(mw2*s)+4*(t-m2)
     9   /s+m2*s/(mw2*(t-m2))+4*s/(t-m2)+4*mw2**2/(s*(t-m2))-8*m2*mw2/(s
     :   *(t-m2))+4*m2**2/(s*(t-m2))-4*mw2/(t-m2)+4*m2/(t-m2)-4*m2*mw2/(
     ;   t-m2)**2+4*m2**2/(t-m2)**2-4*mw2/s+4*m2/s+2*m2/mw2+4)/pi**2
      tmp1 = tmp1-cf*lnmu*nc*nf*(-m2*(t-m2)/(mw2*s)/3.0d0+(-2.0d0)*(t-m2
     1   )/(3.0d0*s)-m2*s/(mw2*(t-m2))/3.0d0+(-2.0d0)*s/(3.0d0*(t-m2))+(
     2   -4.0d0)*mw2**2/(3.0d0*s*(t-m2))+2*m2*mw2/(s*(t-m2))+(-2.0d0)*m2
     3   **3/(3.0d0*mw2*s*(t-m2))+4.0d0*mw2/(3.0d0*(t-m2))+(-2.0d0)*m2**
     4   2/(3.0d0*mw2*(t-m2))+(-2.0d0)*m2/(3.0d0*(t-m2))+4.0d0*m2*mw2/(3
     5   .0d0*(t-m2)**2)+(-2.0d0)*m2**3/(3.0d0*mw2*(t-m2)**2)+(-2.0d0)*m
     6   2**2/(3.0d0*(t-m2)**2)+4.0d0*mw2/(3.0d0*s)+(-2.0d0)*m2**2/(3.0d
     7   0*mw2*s)+(-2.0d0)*m2/(3.0d0*s)+(-2.0d0)*m2/(3.0d0*mw2))/pi**2-c
     8   f*nc*nf*(-m2*(t-m2)/(mw2*s)/3.0d0+(-4.0d0)*(t-m2)/(3.0d0*s)-m2*
     9   s/(mw2*(t-m2))/3.0d0+(-4.0d0)*s/(3.0d0*(t-m2))+(-4.0d0)*mw2**2/
     :   (3.0d0*s*(t-m2))+8.0d0*m2*mw2/(3.0d0*s*(t-m2))+(-4.0d0)*m2**2/(
     ;   3.0d0*s*(t-m2))+4.0d0*mw2/(3.0d0*(t-m2))+(-4.0d0)*m2/(3.0d0*(t-
     <   m2))+4.0d0*m2*mw2/(3.0d0*(t-m2)**2)+(-4.0d0)*m2**2/(3.0d0*(t-m2
     =   )**2)+4.0d0*mw2/(3.0d0*s)+(-4.0d0)*m2/(3.0d0*s)+(-2.0d0)*m2/(3.
     >   0d0*mw2)-1.3333333333333333)/pi**2
      tmp1 = tmp1+cf*lnmu*nc*nf*((-2.0d0)*m2*(t-m2)/(3.0d0*mw2*s)+(-4.0d
     1   0)*(t-m2)/(3.0d0*s)+(-2.0d0)*m2*s/(3.0d0*mw2*(t-m2))+(-4.0d0)*s
     2   /(3.0d0*(t-m2))+(-8.0d0)*mw2**2/(3.0d0*s*(t-m2))+4*m2*mw2/(s*(t
     3   -m2))+(-4.0d0)*m2**3/(3.0d0*mw2*s*(t-m2))+8.0d0*mw2/(3.0d0*(t-m
     4   2))+(-4.0d0)*m2**2/(3.0d0*mw2*(t-m2))+(-4.0d0)*m2/(3.0d0*(t-m2)
     5   )+8.0d0*m2*mw2/(3.0d0*(t-m2)**2)+(-4.0d0)*m2**3/(3.0d0*mw2*(t-m
     6   2)**2)+(-4.0d0)*m2**2/(3.0d0*(t-m2)**2)+8.0d0*mw2/(3.0d0*s)+(-4
     7   .0d0)*m2**2/(3.0d0*mw2*s)+(-4.0d0)*m2/(3.0d0*s)+(-4.0d0)*m2/(3.
     8   0d0*mw2))/pi**2/2.0d0+cf*nc*nf*((-2.0d0)*m2*(t-m2)/(3.0d0*mw2*s
     9   )+(-4.0d0)*(t-m2)/(3.0d0*s)+(-2.0d0)*m2*s/(3.0d0*mw2*(t-m2))+(-
     :   4.0d0)*s/(3.0d0*(t-m2))+(-8.0d0)*mw2**2/(3.0d0*s*(t-m2))+4*m2*m
     ;   w2/(s*(t-m2))+(-4.0d0)*m2**3/(3.0d0*mw2*s*(t-m2))+8.0d0*mw2/(3.
     <   0d0*(t-m2))+(-4.0d0)*m2**2/(3.0d0*mw2*(t-m2))+(-4.0d0)*m2/(3.0d
     =   0*(t-m2))+8.0d0*m2*mw2/(3.0d0*(t-m2)**2)+(-4.0d0)*m2**3/(3.0d0*
     >   mw2*(t-m2)**2)+(-4.0d0)*m2**2/(3.0d0*(t-m2)**2)+8.0d0*mw2/(3.0d
     ?   0*s)+(-4.0d0)*m2**2/(3.0d0*mw2*s)+(-4.0d0)*m2/(3.0d0*s)+(-4.0d0
     @   )*m2/(3.0d0*mw2))/pi**2/2.0d0
      tmp1 = tmp1-ca*cf*lnt*nc*(-m2*(t-m2)/(mw2*s)-4*(t-m2)/s-m2*s/(mw2*
     1   (t-m2))-4*s/(t-m2)-4*mw2**2/(s*(t-m2))+8*m2*mw2/(s*(t-m2))-4*m2
     2   **2/(s*(t-m2))+4*mw2/(t-m2)-4*m2/(t-m2)+4*m2*mw2/(t-m2)**2-4*m2
     3   **2/(t-m2)**2+4*mw2/s-4*m2/s-2*m2/mw2-4)/pi**2-ca*cf*lns*nc*(-m
     4   2*(t-m2)/(mw2*s)-4*(t-m2)/s-m2*s/(mw2*(t-m2))-4*s/(t-m2)-4*mw2*
     5   *2/(s*(t-m2))+8*m2*mw2/(s*(t-m2))-4*m2**2/(s*(t-m2))+4*mw2/(t-m
     6   2)-4*m2/(t-m2)+4*m2*mw2/(t-m2)**2-4*m2**2/(t-m2)**2+4*mw2/s-4*m
     7   2/s-2*m2/mw2-4)/pi**2
      tmp1 = tmp1+ca*cf*lnmu*lnt*nc*(-2*m2*(t-m2)/(mw2*s)-4*(t-m2)/s-2*m
     1   2*s/(mw2*(t-m2))-4*s/(t-m2)-8*mw2**2/(s*(t-m2))+12*m2*mw2/(s*(t
     2   -m2))-4*m2**3/(mw2*s*(t-m2))+8*mw2/(t-m2)-4*m2**2/(mw2*(t-m2))-
     3   4*m2/(t-m2)+8*m2*mw2/(t-m2)**2-4*m2**3/(mw2*(t-m2)**2)-4*m2**2/
     4   (t-m2)**2+8*mw2/s-4*m2**2/(mw2*s)-4*m2/s-4*m2/mw2)/pi**2/2.0d0+
     5   ca*cf*lnt*nc*(-2*m2*(t-m2)/(mw2*s)-4*(t-m2)/s-2*m2*s/(mw2*(t-m2
     6   ))-4*s/(t-m2)-8*mw2**2/(s*(t-m2))+12*m2*mw2/(s*(t-m2))-4*m2**3/
     7   (mw2*s*(t-m2))+8*mw2/(t-m2)-4*m2**2/(mw2*(t-m2))-4*m2/(t-m2)+8*
     8   m2*mw2/(t-m2)**2-4*m2**3/(mw2*(t-m2)**2)-4*m2**2/(t-m2)**2+8*mw
     9   2/s-4*m2**2/(mw2*s)-4*m2/s-4*m2/mw2)/pi**2/2.0d0
      tmp1 = tmp1+ca*cf*lnmu*lns*nc*(-2*m2*(t-m2)/(mw2*s)-4*(t-m2)/s-2*m
     1   2*s/(mw2*(t-m2))-4*s/(t-m2)-8*mw2**2/(s*(t-m2))+12*m2*mw2/(s*(t
     2   -m2))-4*m2**3/(mw2*s*(t-m2))+8*mw2/(t-m2)-4*m2**2/(mw2*(t-m2))-
     3   4*m2/(t-m2)+8*m2*mw2/(t-m2)**2-4*m2**3/(mw2*(t-m2)**2)-4*m2**2/
     4   (t-m2)**2+8*mw2/s-4*m2**2/(mw2*s)-4*m2/s-4*m2/mw2)/pi**2/2.0d0+
     5   ca*cf*lns*nc*(-2*m2*(t-m2)/(mw2*s)-4*(t-m2)/s-2*m2*s/(mw2*(t-m2
     6   ))-4*s/(t-m2)-8*mw2**2/(s*(t-m2))+12*m2*mw2/(s*(t-m2))-4*m2**3/
     7   (mw2*s*(t-m2))+8*mw2/(t-m2)-4*m2**2/(mw2*(t-m2))-4*m2/(t-m2)+8*
     8   m2*mw2/(t-m2)**2-4*m2**3/(mw2*(t-m2)**2)-4*m2**2/(t-m2)**2+8*mw
     9   2/s-4*m2**2/(mw2*s)-4*m2/s-4*m2/mw2)/pi**2/2.0d0
      tmp1 = tmp1+cf**2*lnmu*lnu*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m
     1   2*s/(mw2*(t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(
     2   t-m2))-8*m2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2)
     3   )-8*m2/(t-m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2*
     4   *2/(t-m2)**2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/2.
     5   0d0+cf**2*lnu*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/(mw2*(
     6   t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2))-8*m
     7   2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*m2/(t-
     8   m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(t-m2)*
     9   *2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/2.0d0
      tmp1 = tmp1-cf**2*lnmu**2*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2
     1   *s/(mw2*(t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t
     2   -m2))-8*m2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))
     3   -8*m2/(t-m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**
     4   2/(t-m2)**2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/8.0
     5   d0-ca*cf*lnmu**2*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/(mw
     6   2*(t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2))-
     7   8*m2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*m2/
     8   (t-m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(t-m
     9   2)**2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/8.0d0
      tmp1 = tmp1-cf**2*lnmu*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/
     1   (mw2*(t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2
     2   ))-8*m2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*
     3   m2/(t-m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(
     4   t-m2)**2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/4.0d0-
     5   ca*cf*lnmu*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/(mw2*(t-m
     6   2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2))-8*m2**
     7   3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*m2/(t-m2)
     8   +16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(t-m2)**2+
     9   16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/4.0d0
      tmp1 = tmp1-cf**2*nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/(mw2*
     1   (t-m2))-8*s/(t-m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2))-8*
     2   m2**3/(mw2*s*(t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*m2/(t
     3   -m2)+16*m2*mw2/(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(t-m2)
     4   **2+16*mw2/s-8*m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/4.0d0-ca*cf
     5   *nc*(-4*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-4*m2*s/(mw2*(t-m2))-8*s/(t
     6   -m2)-16*mw2**2/(s*(t-m2))+24*m2*mw2/(s*(t-m2))-8*m2**3/(mw2*s*(
     7   t-m2))+16*mw2/(t-m2)-8*m2**2/(mw2*(t-m2))-8*m2/(t-m2)+16*m2*mw2
     8   /(t-m2)**2-8*m2**3/(mw2*(t-m2)**2)-8*m2**2/(t-m2)**2+16*mw2/s-8
     9   *m2**2/(mw2*s)-8*m2/s-8*m2/mw2)/pi**2/4.0d0
      tmp1 = tmp1-cf**2*lnu*nc*(-2*m2*(t-m2)/(mw2*s)-8*(t-m2)/s-2*m2*s/(
     1   mw2*(t-m2))-8*s/(t-m2)-8*mw2**2/(s*(t-m2))+16*m2*mw2/(s*(t-m2))
     2   -8*m2**2/(s*(t-m2))+8*mw2/(t-m2)-8*m2/(t-m2)+8*m2*mw2/(t-m2)**2
     3   -8*m2**2/(t-m2)**2+8*mw2/s-8*m2/s-4*m2/mw2-8)/pi**2-ca*cf*nc*(2
     4   *cs4f*m2*s**3/(t-m2)-2*cs4f*m2*mw2*s**2/(t-m2)+2*cs4f*m2**2*s**
     5   2/(t-m2)-2*bs4f*m2*s**2/(t-m2)-2*bs2f*m2*s**2/(t-m2)+4*bs1f*m2*
     6   s**2/(t-m2)-2*bs4f*m2*mw2*s/(t-m2)+2*bs2f*m2*mw2*s/(t-m2)+2*bs4
     7   f*m2**2*s/(t-m2)-2*bs2f*m2**2*s/(t-m2)+4*cs4f*m2*s**2-2*bs4f*s*
     8   *2+2*bs1f*s**2+2*bs4f*mw2*s-2*bs1f*mw2*s+2*bs4f*m2*s-4*bs2f*m2*
     9   s+2*bs1f*m2*s)/(lam*pi**2)-cf**2*(8*bs4f*mw2**2-8*bs3f*mw2**2)*
     :   nc/(pi**2*(t+s-m2)**2)-ca*cf*(4*bs3f*mw2**2-4*bs4f*mw2**2)*nc/(
     ;   pi**2*(t+s-m2)**2)
      tmp1 = gs**2*gw**2*pi**2*tmp1/(nc*(nc**2-1))/4.0d0
      tmp0 = -8*m2*(t-m2)/(mw2*s)-16*(t-m2)/s-8*m2*s/(mw2*(t-m2))-16*s/(
     1   t-m2)-32*mw2**2/(s*(t-m2))+48*m2*mw2/(s*(t-m2))-16*m2**3/(mw2*s
     2   *(t-m2))+32*mw2/(t-m2)-16*m2**2/(mw2*(t-m2))-16*m2/(t-m2)+32*m2
     3   *mw2/(t-m2)**2-16*m2**3/(mw2*(t-m2)**2)-16*m2**2/(t-m2)**2+32*m
     4   w2/s-16*m2**2/(mw2*s)-16*m2/s-16*m2/mw2
      tmp0 = (lnmu**2/5.12d2+lnmu/2.56d2+0.00390625)*tmp0/pi**2
      tmp0 = tmp0+(lnmu/2.56d2+0.00390625)*(8*m2*(t-m2)/(mw2*s)+32*(t-m2
     1   )/s+8*m2*s/(mw2*(t-m2))+32*s/(t-m2)+32*mw2**2/(s*(t-m2))-64*m2*
     2   mw2/(s*(t-m2))+32*m2**2/(s*(t-m2))-32*mw2/(t-m2)+32*m2/(t-m2)-3
     3   2*m2*mw2/(t-m2)**2+32*m2**2/(t-m2)**2-32*mw2/s+32*m2/s+16*m2/mw
     4   2+32)/pi**2+(-16*(t-m2)/s-16*s/(t-m2)-32)/pi**2/2.56d2
      tmp0 = 8*cf*(cf+ca)*gs**2*gw**2*pi**2*tmp0/(nc**2-1)
      tmp0 = tmp1+tmp0
      svnonfac = (8*m2*(t-m2)/(mw2*s)+32*(t-m2)/s+8*m2*s/(mw2*(t-m2))+32
     1   *s/(t-m2)+32*mw2**2/(s*(t-m2))-64*m2*mw2/(s*(t-m2))+32*m2**2/(s
     2   *(t-m2))-32*mw2/(t-m2)+32*m2/(t-m2)-32*m2*mw2/(t-m2)**2+32*m2**
     3   2/(t-m2)**2-32*mw2/s+32*m2/s+16*m2/mw2+32)/pi**2/2.56d2
      svnonfac = (lnmu/2.56d2+0.00390625)*(-8*m2*(t-m2)/(mw2*s)-16*(t-m2
     1   )/s-8*m2*s/(mw2*(t-m2))-16*s/(t-m2)-32*mw2**2/(s*(t-m2))+48*m2*
     2   mw2/(s*(t-m2))-16*m2**3/(mw2*s*(t-m2))+32*mw2/(t-m2)-16*m2**2/(
     3   mw2*(t-m2))-16*m2/(t-m2)+32*m2*mw2/(t-m2)**2-16*m2**3/(mw2*(t-m
     4   2)**2)-16*m2**2/(t-m2)**2+32*mw2/s-16*m2**2/(mw2*s)-16*m2/s-16*
     5   m2/mw2)/pi**2+svnonfac
      svnonfac = 8*cf*gs**2*gw**2*(-nf/3.0d0+ca*(lnu-lnt-lns+1.833333333
     1   3333333)+cf*(5-4*lnu)/2.0d0)*pi**2*svnonfac/(nc**2-1)
      svnonfac = tmp0+svnonfac
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
      function xmadevwt(iborn,jproc,idr,s,tk,uk,xmom)
c Wrapper for MadEvent functions. Inputs are 
c   iborn = 0(born), 1(real)
c   jproc = 1(gg), 2(qq), 3(qg)
c   idr   = 1..7 (depends on jproc)
c   s     = parton cm energy squared
c   tk,uk = FNR invariants
c   xmom  = 4-momenta obtained from invar
c Output is the matrix element squared in GeV^-2, times the flux factor,
c times 4*tk*uk/s**2 in the case of real matrix elements.
c MadEvent routines use the SM parameters as defined in setmepar();
c for kinematic-dependent couplings, and if simple factorized
c expressions can't be found (which appears to happen only if a Higgs 
c is involved in the reaction), the call to setmepar() must be included
c in this routine. 
c In the present case, all the Born (real) formulae are proportional 
c to g^2 gw^8 (g^4 gw^8), and this call can be placed somewhere
c else (and done only once). The main code assumes that g^n gw^8 be
c stripped off matrix element routines, which is equivalent to setting
c g and gw equal to one. Since setmepar() is called, for consistency with
c other processes, with g=1 and e=1, gw=1 implies that the results of
c MadEvent must be multiplied by sthw2^4=e^8/gw^8
      implicit none
      integer iborn,jproc,idr
      real*8 xmadevwt,s,tk,uk,xmom(10,4),xfact,tmp(1)
      real*8 pme0dc(0:3,7),pme1dc(0:3,8)
      real*8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      integer ipart,icomp
c Components: MC@NLO conventions   -> 1=px, 2=py, 3=pz, 4=E
c             MadEvent conventions -> 0=E, 1=px, 2=py, 3=pz
      integer mapcomp(0:3)
      data mapcomp/4,1,2,3/
c Labelling conventions for subprocess: 
c MC@NLO   -> x(1)y(2)  ->  z(3)t(4)[->l+(6)nu(7)b(8)]W-(5)[->l-(9)nub(10)]
c MadEvent -> b(1)g(2)  ->  l+(3)nu(4)b(5) l-(6)nub(7)
c MadEvent -> g(1)g(2)  ->  l+(3)nu(4)b(5) l-(6)nub(7) bx(8)
c MadEvent -> u(1)ux(2) ->  l+(3)nu(4)b(5) l-(6)nub(7) bx(8)
c MadEvent -> b(1)u(2)  ->  l+(3)nu(4)b(5) l-(6)nub(7) u(8)
c MadEvent -> b(1)bx(2) ->  l+(3)nu(4)b(5) l-(6)nub(7) bx(8)
c MadEvent -> b(1)b(2)  ->  l+(3)nu(4)b(5) l-(6)nub(7) b(8)
c MadEvent -> b(1)g(2)  ->  l+(3)nu(4)b(5) l-(6)nub(7) g(8)
c In the case of bb initial state, the two final state b quarks are treated
c as different particles, by excluding the diagrams which can be otained
c from other diagrams with the formal exchange 5<->8
      integer mapdir(8),mapref(8)
c mapping for direct events
      data mapdir/1,2,6,7,8,9,10,3/
c mapping for reflected events
      data mapref/2,1,6,7,8,9,10,3/
c
      if(iborn.eq.0)then
        if(jproc.ne.3)then
          write(*,*)'Error #1 in xmadevwt: iborn, jproc=',iborn,jproc
          stop
        endif
        xfact=1.d0
        do ipart=1,7
          do icomp=0,3
            if(idr.eq.1)then
              pme0dc(icomp,ipart)=xmom(mapdir(ipart),mapcomp(icomp))
            elseif(idr.eq.3)then
              pme0dc(icomp,ipart)=xmom(mapref(ipart),mapcomp(icomp))
            else
              write(*,*)'Error #2 in xmadevwt: idr=',idr
              stop
            endif
          enddo
        enddo
      elseif(iborn.eq.1)then
        if( (jproc.eq.1.and.idr.ne.1) .or.
     #      (jproc.eq.3.and.(idr.ne.1.and.idr.ne.3)) )then
          write(*,*)'Error #3 in xmadevwt: jproc,idr=',jproc,idr
          stop
        endif
        xfact=4*tk*uk/s**2
        do ipart=1,8
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
        write(*,*)'xmadevwt: unknown iborn value',iborn
        stop
      endif
      if(iborn.eq.0)then
        call bgtw_dc(pme0dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.1)then
        call ggtwb_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.1.or.idr.eq.2))then
        call uuxtwbx_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.3.or.idr.eq.4))then
        call butwu_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2 .and.
     #       (idr.eq.5.or.idr.eq.6))then
        call bbxtwbx_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.2.and.idr.eq.7)then
        call bbtwb_dc(pme1dc,tmp)
      elseif(iborn.eq.1.and.jproc.eq.3)then
        call bgtwg_dc(pme1dc,tmp)
      else
        write(*,*)'xmadevwt: Error #5'
        stop
      endif
c Insert sin(theta_W) for e^8 -> gw^8
      xmadevwt=sthw2**4*xfact*tmp(1)/(2*s)
      return
      end


      subroutine setmepar(xiwmass,xiwwidth,xizmass,xizwidth,
     #                    xitmass,xitwidth,xibmass,xisin2w,xiee2,xig)
c Fills HELAS common blocks for masses and couplings. The electron charge
c squared and the masses may eventually be passed through a common block
c on a event-by-event basis. This code is mainly taken from coupsm-ORIGINAL.F 
c of the HELAS package. Here, we limit ourselves to setting the following
c parameters:
c
c       real    gw                : weak coupling constant
c       real    gwwa              : dimensionless WWA  coupling
c       real    gwwz              : dimensionless WWZ  coupling
c       complex gal(2)            : coupling with A of charged leptons
c       complex gau(2)            : coupling with A of up-type quarks
c       complex gad(2)            : coupling with A of down-type quarks
c       complex gwf(2)            : coupling with W-,W+ of fermions
c       complex gzn(2)            : coupling with Z of neutrinos
c       complex gzl(2)            : coupling with Z of charged leptons
c       complex gzu(2)            : coupling with Z of up-type quarks
c       complex gzd(2)            : coupling with Z of down-type quarks
c       complex gg(2)             : QCD gqq coupling (L,R)
c
c through the following parameters, given in input
c
c       real    zmass,wmass       : weak boson masses
c       real    zwidth,wwidth     : weak boson width
c       real    tmass,twidth      : top mass and width
c       real    bmass             : bottom mass
c       real    sin2w             : square of sine of the weak angle
c       real    ee2               : positron charge squared
c       real    g                 : QCD 3-,4-gluon coupling
c
      implicit none
      real * 8 xiwmass,xiwwidth,xizmass,xizwidth,xitmass,xitwidth,
     #         xibmass,xisin2w,xiee2,xig
      include "MEcoupl.inc"
      double precision zero,half,one,two,three,pi,ee2,sw,cw,ez,ey,sc2,v
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
      parameter (three=3.d0)
      parameter (pi=3.14159265358979312D0)
c
      wmass = xiwmass
      wwidth= xiwwidth
      zmass = xizmass
      zwidth= xizwidth
      tmass = xitmass
      twidth= xitwidth
      bmass = xibmass
      sin2w = xisin2w
      ee2   = xiee2
      g     = xig
c
      amass=0.d0
      awidth=1.d-99
c
      ee=sqrt(ee2)
      alpha=ee2/(4*pi)
c
      sw  = sqrt( sin2w )
      cw  = sqrt( One - sin2w )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*zmass*sqrt(sc2)/ee
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )
c
      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )
c
      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )
c
c QCD coupling
c
      gg(1) = dcmplx( -g, Zero )
      gg(2) = gg(1)
c
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
      SUBROUTINE BGTW_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b mu- vm~  
C  
C Crossing   1 is b g -> e+ ve b mu- vm~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  7)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 128, NCROSS=  1)
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
      REAL*8 MATBGTW_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,7) / 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,7) / 1, 2, 3, 4, 5, 6, 7/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBGTW_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATBGTW_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBGTW_DC
      REAL*8 FUNCTION MATBGTW_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b mu- vm~  
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
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  13, NCOLOR=   1) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[5,1,2]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,8   ))    
      CALL FVOXXX(W(1,5   ),W(1,8   ),GWF ,TMASS   ,TWIDTH  ,W(1,9   ))    
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL JIOXXX(W(1,1   ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,11  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,11  ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL JIOXXX(W(1,12  ),W(1,9   ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(2   ))            
      JAMP(   1) = -AMP(   1)-AMP(   2)
      MATBGTW_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBGTW_DC =MATBGTW_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       

C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE GGTWB_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> e+ ve b mu- vm~ b~  
C  
C Crossing   1 is g g -> e+ ve b mu- vm~ b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATGGTWB_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATGGTWB_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATGGTWB_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATGGTWB_DC
      REAL*8 FUNCTION MATGGTWB_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> e+ ve b mu- vm~ b~  
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
      parameter (nexternal=  8)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  29, NCOLOR=   2) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[5,8,2,1]                                                 
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[5,8,1,2]                                                 
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL IXXXXX(P(0,8   ),BMASS ,NHEL(8   ),-1*IC(8   ),W(1,8   ))       
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL FVOXXX(W(1,5   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,8   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,11  ))     
      CALL FVOXXX(W(1,10  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL JIOXXX(W(1,11  ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,12  ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,14  ))     
      CALL JIOXXX(W(1,8   ),W(1,14  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,15  ),GWF ,AMP(2   ))            
      CALL FVIXXX(W(1,8   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,16  ))     
      CALL FVOXXX(W(1,10  ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,17  ))     
      CALL JIOXXX(W(1,16  ),W(1,17  ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,17  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,19  ))     
      CALL JIOXXX(W(1,8   ),W(1,19  ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,20  ),GWF ,AMP(4   ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,21  ))     
      CALL JIOXXX(W(1,21  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,22  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,22  ),GWF ,AMP(5   ))            
      CALL FVIXXX(W(1,16  ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,23  ))     
      CALL JIOXXX(W(1,23  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,24  ),GWF ,AMP(6   ))            
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,25  ))      
      CALL FVIXXX(W(1,8   ),W(1,25  ),GG ,BMASS   ,ZERO    ,W(1,26  ))     
      CALL JIOXXX(W(1,26  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,27  ),GWF ,AMP(7   ))            
      CALL FVOXXX(W(1,10  ),W(1,25  ),GG ,TMASS   ,TWIDTH  ,W(1,28  ))     
      CALL JIOXXX(W(1,8   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,29  ),GWF ,AMP(8   ))            
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels are not the same as in the case of production without decay
      AMP(   2)=0
      AMP(   4)=0
      AMP(   8)=0
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   5)-AMP(   7)-AMP(   8)
      JAMP(   2) = +AMP(   3)+AMP(   4)+AMP(   6)+AMP(   7)+AMP(   8)
      MATGGTWB_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATGGTWB_DC =MATGGTWB_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BGTWG_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b mu- vm~ g  
C  
C Crossing   1 is b g -> e+ ve b mu- vm~ g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATBGTWG_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /   10/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBGTWG_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATBGTWG_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBGTWG_DC
      REAL*8 FUNCTION MATBGTWG_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b g -> e+ ve b mu- vm~ g  
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
      parameter (nexternal=  8)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  34, NCOLOR=   2) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[5,1,2,8]                                                 
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[5,1,8,2]                                                 
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL VXXXXX(P(0,8   ),ZERO ,NHEL(8   ),+1*IC(8   ),W(1,8   ))        
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL FVOXXX(W(1,5   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,10  ))    
      CALL FVIXXX(W(1,1   ),W(1,8   ),GG ,BMASS   ,ZERO    ,W(1,11  ))     
      CALL FVOXXX(W(1,10  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL JIOXXX(W(1,11  ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL JVVXXX(W(1,8   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,14  ))      
      CALL FVOXXX(W(1,10  ),W(1,14  ),GG ,TMASS   ,TWIDTH  ,W(1,15  ))     
      CALL JIOXXX(W(1,1   ),W(1,15  ),GWF ,WMASS   ,WWIDTH  ,W(1,16  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,16  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,8   ),GG ,BMASS   ,ZERO    ,W(1,17  ))     
      CALL FVOXXX(W(1,17  ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,18  ))    
      CALL FVOXXX(W(1,18  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,19  ))     
      CALL JIOXXX(W(1,1   ),W(1,19  ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,20  ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,10  ),W(1,8   ),GG ,TMASS   ,TWIDTH  ,W(1,21  ))     
      CALL FVOXXX(W(1,21  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,22  ))     
      CALL JIOXXX(W(1,1   ),W(1,22  ),GWF ,WMASS   ,WWIDTH  ,W(1,23  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,23  ),GWF ,AMP(4   ))            
      CALL FVOXXX(W(1,12  ),W(1,8   ),GG ,TMASS   ,TWIDTH  ,W(1,24  ))     
      CALL JIOXXX(W(1,1   ),W(1,24  ),GWF ,WMASS   ,WWIDTH  ,W(1,25  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,25  ),GWF ,AMP(5   ))            
      CALL FVIXXX(W(1,11  ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,26  ))     
      CALL JIOXXX(W(1,26  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,27  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,27  ),GWF ,AMP(6   ))            
      CALL FVIXXX(W(1,1   ),W(1,14  ),GG ,BMASS   ,ZERO    ,W(1,28  ))     
      CALL JIOXXX(W(1,28  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,29  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,30  ))     
      CALL JIOXXX(W(1,30  ),W(1,18  ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,31  ),GWF ,AMP(8   ))            
      CALL JIOXXX(W(1,30  ),W(1,21  ),GWF ,WMASS   ,WWIDTH  ,W(1,32  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,32  ),GWF ,AMP(9   ))            
      CALL FVIXXX(W(1,30  ),W(1,8   ),GG ,BMASS   ,ZERO    ,W(1,33  ))     
      CALL JIOXXX(W(1,33  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,34  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,34  ),GWF ,AMP(10  ))            
C SF: eliminate graphs with gluon emission from the b resulting from top
      AMP(   3)=0
      AMP(   8)=0
      JAMP(   1) = -AMP(   1)+AMP(   2)-AMP(   5)-AMP(   6)+AMP(   7)
      JAMP(   2) = -AMP(   2)-AMP(   3)-AMP(   4)-AMP(   7)-AMP(   8)
     &             -AMP(   9)-AMP(  10)
      MATBGTWG_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBGTWG_DC =MATBGTWG_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       

C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BBXTWBX_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> e+ ve b mu- vm~ b~  
C  
C Crossing   1 is b b~ -> e+ ve b mu- vm~ b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATBBXTWBX_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBBXTWBX_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATBBXTWBX_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBBXTWBX_DC
      REAL*8 FUNCTION MATBBXTWBX_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> e+ ve b mu- vm~ b~  
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
      parameter (nexternal=  8)
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[5,8]T[2,1]                                               
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[5,1]T[2,8]                                               
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL IXXXXX(P(0,8   ),BMASS ,NHEL(8   ),-1*IC(8   ),W(1,8   ))       
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL FVOXXX(W(1,5   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,10  ))    
      CALL JIOXXX(W(1,8   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,11  ))     
      CALL FVOXXX(W(1,10  ),W(1,11  ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL JIOXXX(W(1,1   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,1   ),W(1,11  ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL JIOXXX(W(1,14  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,15  ),GWF ,AMP(2   ))            
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,16  ))     
      CALL FVIXXX(W(1,8   ),W(1,16  ),GG ,BMASS   ,ZERO    ,W(1,17  ))     
      CALL JIOXXX(W(1,17  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,18  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,18  ),GWF ,AMP(3   ))            
      CALL FVOXXX(W(1,10  ),W(1,16  ),GG ,TMASS   ,TWIDTH  ,W(1,19  ))     
      CALL JIOXXX(W(1,8   ),W(1,19  ),GWF ,WMASS   ,WWIDTH  ,W(1,20  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,20  ),GWF ,AMP(4   ))            
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels are not the same as in the case of production without decay
      AMP(   4)=0
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      MATBBXTWBX_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBBXTWBX_DC =MATBBXTWBX_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BUTWU_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b u -> e+ ve b mu- vm~ u  
C  
C Crossing   1 is b u -> e+ ve b mu- vm~ u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATBUTWU_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C   
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBUTWU_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATBUTWU_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel 
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBUTWU_DC
      REAL*8 FUNCTION MATBUTWU_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b u -> e+ ve b mu- vm~ u  
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
      parameter (nexternal=  8)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  15, NCOLOR=   1) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[8,1]T[5,2]                                               
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL OXXXXX(P(0,8   ),ZERO ,NHEL(8   ),+1*IC(8   ),W(1,8   ))        
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL FVOXXX(W(1,5   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,10  ))    
      CALL JIOXXX(W(1,2   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,11  ))     
      CALL FVOXXX(W(1,10  ),W(1,11  ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL JIOXXX(W(1,1   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,1   ),W(1,11  ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL JIOXXX(W(1,14  ),W(1,10  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,15  ),GWF ,AMP(2   ))            
      JAMP(   1) = +AMP(   1)+AMP(   2)
      MATBUTWU_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBUTWU_DC =MATBUTWU_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE BBTWB_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b -> e+ ve b mu- vm~ b  
C  
C Crossing   1 is b b -> e+ ve b mu- vm~ b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATBBTWB_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATBBTWB_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATBBTWB_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
c SF: insert a factor 2, since identical final-state b's have been treated
c SF as different particles
      ANS(IPROC)=2*ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATBBTWB_DC
      REAL*8 FUNCTION MATBBTWB_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b -> e+ ve b mu- vm~ b  
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
      parameter (nexternal=  8)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  31, NCOLOR=   2) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[8,1]T[5,2]                                               
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[8,2]T[5,1]                                               
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL OXXXXX(P(0,8   ),BMASS ,NHEL(8   ),+1*IC(8   ),W(1,8   ))       
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,9   ))    
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL FVOXXX(W(1,8   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,11  ))    
      CALL FVOXXX(W(1,11  ),W(1,10  ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL JIOXXX(W(1,2   ),W(1,12  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL FVIXXX(W(1,2   ),W(1,10  ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL JIOXXX(W(1,14  ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,15  ),GWF ,AMP(2   ))            
      CALL FVOXXX(W(1,5   ),W(1,9   ),GWF ,TMASS   ,TWIDTH  ,W(1,16  ))    
      CALL JIOXXX(W(1,2   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,17  ))     
      CALL FVOXXX(W(1,16  ),W(1,17  ),GG ,TMASS   ,TWIDTH  ,W(1,18  ))     
      CALL JIOXXX(W(1,1   ),W(1,18  ),GWF ,WMASS   ,WWIDTH  ,W(1,19  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,19  ),GWF ,AMP(3   ))            
      CALL FVIXXX(W(1,1   ),W(1,17  ),GG ,BMASS   ,ZERO    ,W(1,20  ))     
      CALL JIOXXX(W(1,20  ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,21  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,21  ),GWF ,AMP(4   ))            
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,22  ))     
      CALL FVOXXX(W(1,11  ),W(1,22  ),GG ,TMASS   ,TWIDTH  ,W(1,23  ))     
      CALL JIOXXX(W(1,1   ),W(1,23  ),GWF ,WMASS   ,WWIDTH  ,W(1,24  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,24  ),GWF ,AMP(5   ))            
      CALL FVIXXX(W(1,1   ),W(1,22  ),GG ,BMASS   ,ZERO    ,W(1,25  ))     
      CALL JIOXXX(W(1,25  ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,26  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,26  ),GWF ,AMP(6   ))            
      CALL JIOXXX(W(1,1   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,27  ))     
      CALL FVOXXX(W(1,16  ),W(1,27  ),GG ,TMASS   ,TWIDTH  ,W(1,28  ))     
      CALL JIOXXX(W(1,2   ),W(1,28  ),GWF ,WMASS   ,WWIDTH  ,W(1,29  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,29  ),GWF ,AMP(7   ))            
      CALL FVIXXX(W(1,2   ),W(1,27  ),GG ,BMASS   ,ZERO    ,W(1,30  ))     
      CALL JIOXXX(W(1,30  ),W(1,16  ),GWF ,WMASS   ,WWIDTH  ,W(1,31  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,31  ),GWF ,AMP(8   ))            
C SF: eliminate graphs in which b quark #8 results from the top decay
      AMP(   1)=0
      AMP(   2)=0
      AMP(   5)=0
      AMP(   6)=0
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)
      JAMP(   2) = -AMP(   5)-AMP(   6)-AMP(   7)-AMP(   8)
      MATBBTWB_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATBBTWB_DC =MATBBTWB_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
C SF: The following routine is SMATRIX generated by MadEvent, suitably modified
      SUBROUTINE UUXTWBX_DC(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> e+ ve b mu- vm~ b~  
C  
C Crossing   1 is u u~ -> e+ ve b mu- vm~ b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
C SF: replace the following include with the explicit inclusion of
c nexternal.inc, originally included by genps.inc
c      Include "genps.inc"
      integer    nexternal
      parameter (nexternal=  8)
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB= 256, NCROSS=  1)
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
      REAL*8 MATUUXTWBX_DC
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2

      character*79         hel_buff
C SF: comment out all common blocks
c      common/to_helicity/  hel_buff

      integer          isum_hel
      logical                    multi_channel
C SF: comment out all common blocks
c      common/to_matrix/isum_hel, multi_channel
C SF: comment out all instances of mapconfig, used by multi_channel
c      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
c      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
C SF: comment out all instances of amp2, used by multi_channel
c      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1,8) /-1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1,8) /-1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1,8) /-1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1,8) /-1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1,8) /-1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1,8) /-1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1,8) /-1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1,8) /-1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1,8) /-1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1,8) /-1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1,8) /-1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1,8) /-1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1,8) /-1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1,8) /-1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1,8) /-1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1,8) /-1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  34),IHEL=1,8) /-1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  35),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1,8) /-1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1,8) /-1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  40),IHEL=1,8) /-1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  41),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1,8) /-1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1,8) /-1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  46),IHEL=1,8) /-1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  47),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1,8) /-1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  50),IHEL=1,8) /-1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  51),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  52),IHEL=1,8) /-1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  53),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  54),IHEL=1,8) /-1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1,8) /-1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  57),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  58),IHEL=1,8) /-1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  59),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  60),IHEL=1,8) /-1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  61),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1,8) /-1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  63),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  64),IHEL=1,8) /-1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  65),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  66),IHEL=1,8) /-1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  68),IHEL=1,8) /-1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  69),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  70),IHEL=1,8) /-1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  71),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  72),IHEL=1,8) /-1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1,8) /-1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  75),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  76),IHEL=1,8) /-1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  77),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  78),IHEL=1,8) /-1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  79),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1,8) /-1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  81),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  82),IHEL=1,8) /-1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  83),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  84),IHEL=1,8) /-1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  86),IHEL=1,8) /-1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  87),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  88),IHEL=1,8) /-1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  89),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  90),IHEL=1,8) /-1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1,8) /-1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  93),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  94),IHEL=1,8) /-1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  95),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  96),IHEL=1,8) /-1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  97),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1,8) /-1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  99),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 100),IHEL=1,8) /-1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 101),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 102),IHEL=1,8) /-1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 104),IHEL=1,8) /-1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 105),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 106),IHEL=1,8) /-1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 107),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 108),IHEL=1,8) /-1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 109),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 110),IHEL=1,8) /-1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 111),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 112),IHEL=1,8) /-1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 113),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 114),IHEL=1,8) /-1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 115),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 116),IHEL=1,8) /-1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 117),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 118),IHEL=1,8) /-1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 119),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 120),IHEL=1,8) /-1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 121),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 122),IHEL=1,8) /-1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 123),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 124),IHEL=1,8) /-1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 125),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 126),IHEL=1,8) /-1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 127),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 128),IHEL=1,8) /-1, 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 129),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 130),IHEL=1,8) / 1,-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 131),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 132),IHEL=1,8) / 1,-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 133),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 134),IHEL=1,8) / 1,-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 135),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 136),IHEL=1,8) / 1,-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 137),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 138),IHEL=1,8) / 1,-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 139),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 140),IHEL=1,8) / 1,-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 141),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 142),IHEL=1,8) / 1,-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 143),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 144),IHEL=1,8) / 1,-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 145),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 146),IHEL=1,8) / 1,-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 147),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 148),IHEL=1,8) / 1,-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 149),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 150),IHEL=1,8) / 1,-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 151),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 152),IHEL=1,8) / 1,-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 153),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 154),IHEL=1,8) / 1,-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 155),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 156),IHEL=1,8) / 1,-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 157),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 158),IHEL=1,8) / 1,-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 159),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 160),IHEL=1,8) / 1,-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 161),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 162),IHEL=1,8) / 1,-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 163),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 164),IHEL=1,8) / 1,-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 165),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 166),IHEL=1,8) / 1,-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 167),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 168),IHEL=1,8) / 1,-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 169),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 170),IHEL=1,8) / 1,-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 171),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 172),IHEL=1,8) / 1,-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 173),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 174),IHEL=1,8) / 1,-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 175),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 176),IHEL=1,8) / 1,-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 177),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 178),IHEL=1,8) / 1,-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 179),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 180),IHEL=1,8) / 1,-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 181),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 182),IHEL=1,8) / 1,-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 183),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 184),IHEL=1,8) / 1,-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 185),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 186),IHEL=1,8) / 1,-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 187),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 188),IHEL=1,8) / 1,-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 189),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 190),IHEL=1,8) / 1,-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 191),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 192),IHEL=1,8) / 1,-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 193),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 194),IHEL=1,8) / 1, 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 195),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 196),IHEL=1,8) / 1, 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 197),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 198),IHEL=1,8) / 1, 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 199),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 200),IHEL=1,8) / 1, 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 201),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 202),IHEL=1,8) / 1, 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 203),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 204),IHEL=1,8) / 1, 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 205),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 206),IHEL=1,8) / 1, 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 207),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 208),IHEL=1,8) / 1, 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 209),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 210),IHEL=1,8) / 1, 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 211),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 212),IHEL=1,8) / 1, 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 213),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 214),IHEL=1,8) / 1, 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 215),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 216),IHEL=1,8) / 1, 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 217),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 218),IHEL=1,8) / 1, 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 219),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 220),IHEL=1,8) / 1, 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 221),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 222),IHEL=1,8) / 1, 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 223),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 224),IHEL=1,8) / 1, 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 225),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 226),IHEL=1,8) / 1, 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 227),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 228),IHEL=1,8) / 1, 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 229),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 230),IHEL=1,8) / 1, 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 231),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 232),IHEL=1,8) / 1, 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 233),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 234),IHEL=1,8) / 1, 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 235),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 236),IHEL=1,8) / 1, 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 237),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 238),IHEL=1,8) / 1, 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 239),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 240),IHEL=1,8) / 1, 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL, 241),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL, 242),IHEL=1,8) / 1, 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL, 243),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL, 244),IHEL=1,8) / 1, 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL, 245),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL, 246),IHEL=1,8) / 1, 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL, 247),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL, 248),IHEL=1,8) / 1, 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL, 249),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL, 250),IHEL=1,8) / 1, 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL, 251),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL, 252),IHEL=1,8) / 1, 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL, 253),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 254),IHEL=1,8) / 1, 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 255),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 256),IHEL=1,8) / 1, 1, 1, 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1,8) / 1, 2, 3, 4, 5, 6, 7, 8/
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
       
C SF: comment out all instances of multi_channel
c      IF (multi_channel) THEN
c          DO IHEL=1,NGRAPHS
c              amp2(ihel)=0d0
c              jamp2(ihel)=0d0
c          ENDDO
c          DO IHEL=1,int(jamp2(0))
c              jamp2(ihel)=0d0
c          ENDDO
c      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
              IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATUUXTWBX_DC(P ,NHEL(1,IHEL),JC(1))            
                 ANS(IPROC)=ANS(IPROC)+T
                  IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                      GOODHEL(IHEL,IPROC)=.TRUE.
                      NGOOD = NGOOD +1
                      IGOOD(NGOOD) = IHEL
C                WRITE(*,*) ngood,IHEL,T
                  ENDIF
              ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=MATUUXTWBX_DC(P ,NHEL(1,IHEL),JC(1))            
           ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
C SF: comment out all instances of multi_channel
c      IF (MULTI_CHANNEL) THEN
c          XTOT=0D0
c          DO IHEL=1,MAPCONFIG(0)
c              XTOT=XTOT+AMP2(MAPCONFIG(IHEL))
c          ENDDO
c          ANS(IPROC)=ANS(IPROC)*AMP2(MAPCONFIG(ICONFIG))/XTOT
c      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
C SF: the original name MATRIX has been replaced by MATUUXTWBX_DC
      REAL*8 FUNCTION MATUUXTWBX_DC(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> e+ ve b mu- vm~ b~  
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
      parameter (nexternal=  8)
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  15, NCOLOR=   1) 
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
C SF: comment out all instances of amp2, used by multi_channel
c      Double Precision amp2(maxamps), jamp2(0:maxamps)
c      common/to_amps/  amp2,       jamp2
C SF: The original coupl.inc has been renamed MEcoupl.inc
      include "MEcoupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[5,1]T[2,8]                                               
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL IXXXXX(P(0,3   ),ZERO ,NHEL(3   ),-1*IC(3   ),W(1,3   ))        
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL OXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL IXXXXX(P(0,7   ),ZERO ,NHEL(7   ),-1*IC(7   ),W(1,7   ))        
      CALL IXXXXX(P(0,8   ),BMASS ,NHEL(8   ),-1*IC(8   ),W(1,8   ))       
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL JIOXXX(W(1,3   ),W(1,4   ),GWF ,WMASS   ,WWIDTH  ,W(1,10  ))    
      CALL FVOXXX(W(1,5   ),W(1,10  ),GWF ,TMASS   ,TWIDTH  ,W(1,11  ))    
      CALL FVIXXX(W(1,8   ),W(1,9   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL JIOXXX(W(1,12  ),W(1,11  ),GWF ,WMASS   ,WWIDTH  ,W(1,13  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,13  ),GWF ,AMP(1   ))            
      CALL FVOXXX(W(1,11  ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,14  ))     
      CALL JIOXXX(W(1,8   ),W(1,14  ),GWF ,WMASS   ,WWIDTH  ,W(1,15  ))    
      CALL IOVXXX(W(1,7   ),W(1,6   ),W(1,15  ),GWF ,AMP(2   ))            
C SF: eliminate graphs with interference with ttbar. Note that the
C SF labels are not the same as in the case of production without decay
      AMP(   2)=0
      JAMP(   1) = +AMP(   1)+AMP(   2)
      MATUUXTWBX_DC = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATUUXTWBX_DC =MATUUXTWBX_DC+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C SF: comment out all instances of amp2, used by multi_channel
c      Do I = 1, NGRAPHS
c          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
c      Enddo
c      Do I = 1, NCOLOR
c          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
c      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
