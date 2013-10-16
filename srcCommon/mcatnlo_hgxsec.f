      subroutine f2body(xs,xxii,xyi,xt,xu,res)
c Returns the real matrix elements times xii**2*(1-yi**2)=4*t*u/s**2.
c The normalization is such that
c   sigma_2body=gs**6*gf*(f2body/(xii**2*(1-yi**2)))*dphi_2
      implicit none
      real * 8 xs,xxii,xyi,xt,xu,res(4)
      include 'hgscblks.h'
      real * 8 s,xii,yi,t,u,tmp(4)
      integer jproc0,j
      common/cjproc/jproc0
c
      s=xs
      xii=xxii
      yi=xyi
      t=xt
      u=xu
      if(jproc0.eq.1)then
        call f2body_gg(s,xii,yi,t,u,tmp)
      elseif(jproc0.eq.2)then
        call f2body_qq(s,xii,yi,t,u,tmp)
      elseif(jproc0.eq.3)then
        call f2body_qg(s,xii,yi,t,u,tmp)
      else
        write(*,*)'Unknown process in f2body',jproc0
        stop
      endif
      do j=1,4
        res(j)=tmp(j)
      enddo
      return
      end

        
      subroutine f2body_gg(xs,xxii,xyi,xt,xu,res)
c gg --> Hg, real matrix element
      implicit none
      real * 8 xs,xxii,xyi,xt,xu,res(4)
      include 'hgscblks.h'
      real * 8 s,xii,yi,t,u,tiny,v2,pi,xnorm,x_ap,s_red,
     # ap_kern,vca,xmat,abdv_Rgg,tmp(4)
      integer icode,ione,itwo,j
      parameter (tiny=1.d-6)
c v2=1/sqrt(2)
      parameter (v2=0.70710678118654757d0)
      parameter (pi=3.14159265358979312D0)
      parameter (vca=3.d0)
      parameter (ione=1)
      parameter (itwo=2)
      integer ibornuse
      common/cibornuse/ibornuse
c
      s=xs
      xii=xxii
      yi=xyi
      t=xt
      u=xu
      do j=1,4
        res(j)=0.d0
      enddo
      if(xii.lt.tiny)then
        s_red=s*(1-xii)
        call f1born(s_red,ione,ibornuse,'sf',tmp)
        res(1)=16*vca/s_red*tmp(1)
      elseif(yi.gt.1-tiny)then
        x_ap=1-xii
        s_red=s*x_ap
        icode=1
        call f1born(s_red,ione,ibornuse,'c+',tmp)
        res(1)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*tmp(1)
      elseif(yi.lt.-1+tiny)then
        x_ap=1-xii
        s_red=s*x_ap
        icode=1
        call f1born(s_red,ione,ibornuse,'c-',tmp)
        res(1)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*tmp(1)
      else
        if(ibornuse.eq.2)then
c From eq.(3.18) of NPB359(91)283
          xnorm=32/(256.d0*3.d0*pi*v2*(4*pi)**3)
          xmat=4*(xmh2**4+s**4+t**4+u**4)/s**3
          res(1)=xnorm*xmat/(2*s)
        elseif(ibornuse.eq.3)then
          res(1)=abdv_Rgg(s,t,u,xii,yi,xmh2,xmb2,xmt2)
        endif
      endif
      return
      end


      subroutine f2body_qq(xs,xxii,xyi,xt,xu,res)
c qq --> Hg, real matrix element
      implicit none
      real * 8 xs,xxii,xyi,xt,xu,res(4)
      include 'hgscblks.h'
      real * 8 s,xii,yi,t,u,tiny,v2,pi,xnorm,xmat
      integer j
      parameter (tiny=1.d-6)
c v2=1/sqrt(2)
      parameter (v2=0.70710678118654757d0)
      parameter (pi=3.14159265358979312D0)
      integer ibornuse
      common/cibornuse/ibornuse
c
      s=xs
      xii=xxii
      yi=xyi
      t=xt
      u=xu
      do j=1,4
        res(j)=0.d0
      enddo
      if(xii.lt.tiny)then
        continue
      elseif(yi.gt.1-tiny)then
        continue
      elseif(yi.lt.-1+tiny)then
        continue
      else
        if(ibornuse.eq.2)then
c From eq.(3.1) of NPB359(91)283
          xnorm=16/(36.d0*9.d0*pi*v2*(4*pi)**3)
          xmat=4*t*u*(t**2+u**2)/s**3
          res(1)=xnorm*xmat/(2*s)
          res(3)=res(1)
        elseif(ibornuse.eq.3)then
          call abdv_Rqq(s,t,u,xii,yi,xmh2,xmb2,xmt2,res)
        endif
      endif
      return
      end


      subroutine f2body_qg(xs,xxii,xyi,xt,xu,res)
c qg --> Hq, real matrix element
      implicit none
      real * 8 xs,xxii,xyi,xt,xu,res(4)
      include 'hgscblks.h'
      real * 8 s,xii,yi,t,u,tiny,v2,pi,xnorm,x_ap,s_red,
     # ap_kern,xmatd,xmatr,tmp(4)
      integer icode,ione,itwo,j
      parameter (tiny=1.d-6)
c v2=1/sqrt(2)
      parameter (v2=0.70710678118654757d0)
      parameter (pi=3.14159265358979312D0)
      parameter (ione=1)
      parameter (itwo=2)
      integer ibornuse
      common/cibornuse/ibornuse
c
      s=xs
      xii=xxii
      yi=xyi
      t=xt
      u=xu
      do j=1,4
        res(j)=0.d0
      enddo
      if(xii.lt.tiny)then
        continue
      elseif(yi.gt.1-tiny)then
        x_ap=1-xii
        s_red=s*x_ap
        icode=3
        call f1born(s_red,ione,ibornuse,'c+',tmp)
        do j=1,2
          res(j)=4*(1+yi)/s*ap_kern(x_ap,abs(icode))*tmp(j)
        enddo
      elseif(yi.lt.-1+tiny)then
        x_ap=1-xii
        s_red=s*x_ap
        icode=3
        call f1born(s_red,ione,ibornuse,'c-',tmp)
        do j=3,4
          res(j)=4*(1-yi)/s*ap_kern(x_ap,abs(icode))*tmp(j)
        enddo
      else
        if(ibornuse.eq.2)then
c From eq.(3.3) of NPB359(91)283
          xnorm=-16/(96.d0*9.d0*pi*v2*(4*pi)**3)
          xmatd=4*u*(s**2+u**2)/s**2
          xmatr=4*t*(s**2+t**2)/s**2
          do j=1,2
            res(j)=xnorm*xmatd/(2*s)
          enddo
          do j=3,4
            res(j)=xnorm*xmatr/(2*s)
          enddo
        elseif(ibornuse.eq.3)then
          call abdv_Rqg(s,t,u,xii,yi,xmh2,xmb2,xmt2,res)
        endif
      endif
      return
      end


      subroutine f1born(xs,jproc,iborn,c2,res)
c Born matrix element, times flux factor times normalizations and averages.
c   sigma_born=gs**4*gf*f1born*dphi_1
c The results are given exact in M_top (iborn=1) or in the M_top --> inf 
c limit (iborn=2). We assume that only the top quark is flowing in the loop
      implicit none
      character * 2 c2
      real * 8 xs,res(4)
      integer jproc,iborn,j
      include 'hgscblks.h'
      real * 8 s,v2,pi,tiny,tmp,xnorm,tauq,etapl,etamn,abdv_f1born
      complex * 16 zic,tmpc
      parameter (v2=0.70710678118654757d0)
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-8)
      parameter (zic=(0.d0,1.d0))
c
      s=xs
      if(gah.eq.0.d0.and.abs(s-xmh2).gt.tiny)then
        write(*,*)'Fatal error in f1born',s,xmh2,c2
        stop
      endif
      if(iborn.eq.1.or.iborn.eq.2)then
        xnorm=1/(2.d0*pi)*1/(pi*(4.d0*pi)**2)*xmh2/(256*v2)
      else
        xnorm=1.d0
      endif
      if(jproc.eq.1)then
        if(iborn.eq.1)then
c From eq.(2.2) of NPB359(91)283
          tauq=4*xmt2/xmh2
          if(tauq.gt.1.d0)then
            tmp=tauq*(1+(1-tauq)*(asin(1/sqrt(tauq)))**2)
            tmp=tmp**2
          else
            etapl=1+sqrt(1-tauq)
            etamn=1-sqrt(1-tauq)
            tmpc=tauq*(1-(1-tauq)*(log(etapl/etamn)-zic*pi)**2/4.d0)
            tmp=abs(tmpc)**2
          endif
        elseif(iborn.eq.2)then
          tmp=4.d0/9.d0
        elseif(iborn.eq.3)then
c Exact top and bottom mass dependence
          tmp=abdv_f1born(xmh2,xmb2,xmt2)
        else
          write(*,*)'Unknown option in f1born',iborn
          stop
        endif
        do j=1,4
          res(j)=xnorm*tmp
        enddo
      else
        do j=1,4
          res(j)=0.d0
        enddo
      endif
      return
      end


      subroutine f1sv(xs,jproc,res)
c Returns sig_2pv of FKS. It is derived from the subroutine xmatel_2pv_contr
c of that package. The Ellis-Sexton scale is set equal to the factorization
c scale, thus there's no contribution from Q. The normalization is
c   sigma_f1sv=gs**6*gf/(8*pi**2)*f1sv*dphi_1
      implicit none
      real * 8 xs,res(4)
      include 'hgscblks.h'
      real * 8 tiny,pi,zero,s,eikcon,fincon,fincon2,vca,xicut,delta,
     # xmat,abdv_f1sv_nonpi,tmp(4)
      integer jproc,ione,itwo,j
      common/parsub/xicut,delta
      parameter (tiny=1.d-8)
      parameter (pi=3.14159265358979312D0)
      parameter (zero=0.d0)
      parameter (vca=3.d0)
      parameter (ione=1)
      parameter (itwo=2)
      integer ibornuse
      common/cibornuse/ibornuse
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c
      s=xs
      if(gah.eq.0.d0.and.abs(s-xmh2).gt.tiny)then
        write(*,*)'Fatal error in f1sv',s,xmh2
        stop
      endif
      if(abs(xmuf2h1-xmuf2h2).gt.tiny .or.
     #   abs(xmuf2h1-xmur2).gt.tiny)then
        write(*,*)'No such scale choice'
        stop
      endif
      fincon2=0.d0
      if(jproc.eq.1)then
        eikcon=2*vca*( 0.5d0*log(xicut**2*s/xmuf2h1)**2-
     #                 pi**2/6.d0 )
        if(ibornuse.eq.2)then
          fincon=-3*log(xmh2/xmuf2h1)**2+3*pi**2+11
        elseif(ibornuse.eq.3)then
          if(iabdvnotop.eq.0)then
            fincon=-3*log(xmh2/xmuf2h1)**2+3*pi**2+
     #             2*abdv_f1sv_nonpi(xmh2,xmb2,xmt2,xmur2)
          elseif(iabdvnotop.eq.1)then
            fincon=-3*log(xmh2/xmuf2h1)**2+3*pi**2
            fincon2=2*abdv_f1sv_nonpi(xmh2,xmb2,xmt2,xmur2)-
     #              2*abdv_f1sv_nonpi(xmh2,zero,xmt2,xmur2)
          endif
        endif
        xmat=eikcon+fincon
      else
        xmat=0.d0
      endif
      if(xmat.ne.0.d0)call f1born(s,ione,ibornuse,'sv',tmp)
      do j=1,4
        res(j)=xmat*tmp(j)+fincon2
      enddo
      return
      end


      subroutine f2b_coll(xs,xxii,xxiic,xyic,xxlmude,res)
c Returns sig_2pr of FKS. It is derived from the subroutine xmatel_coll
c of that package, which returns xmtel and xmtel_sc, the latter being
c the contribution of the delta term and of the regular part of the
c change of scheme, as defined in NPB357(91)409. These contributions
c are not associated to plus prescriptions: therefore, they are multiplied
c here by xii, since a factor 1/xii appears in the main code. The possible
c numerical inaccuracies motivated the definition of the subroutine in
c the jet package -- here we just ignore the problem. The normalization is
c   sigma_2pr=gs**6*gf/(8*pi**2)*f2b_coll/xii*dphi_1
      implicit none
      real * 8 xs,xxii,xxiic,xyic,xxlmude,res(4)
      include 'hgscblks.h'
      real * 8 s,xii,xiic,yic,xlmude,x_ap,s_red,one,xicut,delta,
     # xdfct1,xdfct2,xdfct3p,xdfct3l,xdfct5,xrfct1,xrfct2,xrfct3p,
     # xrfct3l,xrfct5,ap_kern,apprime_kern,xkplus,xklog,xkreg,
     # xkdelta,xfct4(4)
      common/parsub/xicut,delta
      parameter (one=1.d0)
      character * 2 scheme
      integer jproc0,icoded,icoder,j
      common/cjproc/jproc0
      integer ione,itwo
      parameter (ione=1)
      parameter (itwo=2)
      integer ibornuse
      common/cibornuse/ibornuse
c
      s=xs
      xii=xxii
      xiic=xxiic
      yic=xyic
      xlmude=xxlmude
c
      x_ap=1-xiic
      s_red=s*x_ap
      if(yic.eq.1.d0)then
        scheme=schhad1
        if(jproc0.eq.1)then
          icoded=1
          icoder=1
        elseif(jproc0.eq.2)then
          icoded=0
          icoder=0
        elseif(jproc0.eq.3)then
          icoded=3
          icoder=0
        else
          write(*,*)'Unknown process in f2b_coll',jproc0
          stop
        endif
      elseif(yic.eq.-1.d0)then
        scheme=schhad2
        if(jproc0.eq.1)then
          icoded=1
          icoder=1
        elseif(jproc0.eq.2)then
          icoded=0
          icoder=0
        elseif(jproc0.eq.3)then
          icoded=0
          icoder=3
        else
          write(*,*)'Unknown process in f2b_coll',jproc0
          stop
        endif
      else
        write(6,*)'Error in f2b_coll',yic
        stop
      endif
      if(icoded.ne.0.or.icoder.ne.0)then
        call f1born(s_red,ione,ibornuse,'pr',xfct4)
      else
        do j=1,4
          res(j)=0.d0
        enddo
      endif


      if(icoded.ne.0)then
        xdfct1=ap_kern(x_ap,abs(icoded))
        xdfct2=apprime_kern(x_ap,abs(icoded))
        xdfct3p=0.d0
        xdfct3l=0.d0
        xdfct5=0.d0
c
        if(scheme.eq.'DI')then
        xdfct3p=xkplus(x_ap,abs(icoded))
          xdfct3l=xklog(x_ap,abs(icoded))
          if(xiic.ne.0.d0)then
            xdfct5=xkreg(x_ap,abs(icoded))
          else
            xdfct5=xkdelta(abs(icoded))
     #            +xkplus(one,abs(icoded))*log(xicut)
     #            +xklog(one,abs(icoded))*log(xicut)**2/2.d0
c This part contributes to sig2pr(soft), which is integrated in xi
c over the range (0,xicut). This implies the presence of a jacobian
c equal to xicut in the soft term, which has to be removed by hand
c in this case
            xdfct5=xdfct5/xicut
          endif
        elseif(scheme.ne.'MS')then
          write(6,*)'Error in f2b_coll, y=',yic
          write(6,*)'Factorization scheme ',scheme,' not known'
        endif
c
        do j=1,2
          res(j)=( xdfct1*(xlmude+2*log(xii))-xdfct2
     #            -xdfct3p-xdfct3l*log(xii) )*xfct4(j)
     #           -xii*xdfct5*xfct4(j)
        enddo
      else
        do j=1,2
          res(j)=0.d0
        enddo
      endif
      if(icoder.ne.0)then
        xrfct1=ap_kern(x_ap,abs(icoder))
        xrfct2=apprime_kern(x_ap,abs(icoder))
        xrfct3p=0.d0
        xrfct3l=0.d0
        xrfct5=0.d0
c
        if(scheme.eq.'DI')then
        xrfct3p=xkplus(x_ap,abs(icoder))
          xrfct3l=xklog(x_ap,abs(icoder))
          if(xiic.ne.0.d0)then
            xrfct5=xkreg(x_ap,abs(icoder))
          else
            xrfct5=xkdelta(abs(icoder))
     #            +xkplus(one,abs(icoder))*log(xicut)
     #            +xklog(one,abs(icoder))*log(xicut)**2/2.d0
c This part contributes to sig2pr(soft), which is integrated in xi
c over the range (0,xicut). This implies the presence of a jacobian
c equal to xicut in the soft term, which has to be removed by hand
c in this case
            xrfct5=xrfct5/xicut
          endif
        elseif(scheme.ne.'MS')then
          write(6,*)'Error in f2b_coll, y=',yic
          write(6,*)'Factorization scheme ',scheme,' not known'
        endif
c
        do j=3,4
          res(j)=( xrfct1*(xlmude+2*log(xii))-xrfct2
     #            -xrfct3p-xrfct3l*log(xii) )*xfct4(j)
     #           -xii*xrfct5*xfct4(j)
        enddo
      else
        do j=3,4
          res(j)=0.d0
        enddo
      endif
      return
      end

c
c From the jet package, Altarelli-Parisi kernels and change of scheme
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
c ABDV routines (for exact mass dependence)
c
c
      function abdv_f1born(xmh2,xmb2,xmt2)
c Eq.(11) of 1111.2854, with the same normalization as f1born,
c and the sum of amplitudes of eq.(22)
      implicit none
      real*8 abdv_f1born,xmh2,xmb2,xmt2
      complex*16 abdv_H1l_bare,tmp,tmpt
      real*8 sqrt2,pi,xnorm,dummy
      parameter (sqrt2=1.4142135623730950488d0)
      parameter (pi=3.1415926535897932385d0)
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c
      xnorm=xmh2/(256*sqrt2*pi**2*(4*pi)**2)
      tmp=dcmplx(0.d0,0.d0)
      tmpt=dcmplx(0.d0,0.d0)
      if(xmb2.ne.0.d0)tmp=tmp+abdv_H1l_bare(xmh2,xmb2)
      if(xmt2.ne.0.d0)then
        tmpt=abdv_H1l_bare(xmh2,xmt2)
        tmp=tmp+tmpt
      endif
      if(iabdvnotop.eq.0)then
        dummy=xnorm*abs(tmp)**2
      elseif(iabdvnotop.eq.1)then
        dummy=xnorm*(abs(tmp)**2-abs(tmpt)**2)
      endif
      abdv_f1born=dummy
      return
      end


      function abdv_H1l_bare(xmh2,xmq2)
c Eq.(22) of 1111.2854, for one heavy quark flavour of mass squared xmq2
      implicit none
      complex*16 abdv_H1l_bare
      real*8 xmh2,xmq2
      complex*16 xq,abdv_sqrt
      real*8 yq,vtf,lambdaq
      parameter (vtf=0.5d0)
      parameter (lambdaq=1d0)
c
      yq=xmq2/xmh2
      xq=(abdv_sqrt(1-4*yq)-1)/(abdv_sqrt(1-4*yq)+1)
      abdv_H1l_bare=-4*vtf*lambdaq*yq*(2-(1-4*yq)*log(xq)**2/2.d0)
      return
      end


      function abdv_f1sv_nonpi(xmh2,xmb2,xmt2,xmu2)
c This is 2*Re(H^{(2l)}/H^{(1l)}), which appears in eq.(13) of 1111.2854,
c when iabdvnotop=0 (all contributions). When iabdvnotop=1 (ie, exclude top^2)
c it returns 2*xnorm*Re(conjug[H^{(2l)}] H^{(1l)}). The two are equivalent,
c being based on the complex-number identity 
c  2*Re(a*conjug[b]) = 2*Re(b/a)*|a|^2
c and because the return value of this function is not multiplied
c by the Born squared in f1sv when iabdvnotop=1 (hence the necessity of
c the normalization xnorm, included in Born when iabdvnotop=0).
c For the normalization of this function and its use in the virtual
c contribution to an FKS cross section, see virt_abdv.mac
      implicit none
      real*8 abdv_f1sv_nonpi,xmh2,xmb2,xmt2,xmu2
      complex*16 tmpn,tmpd,abdv_H1l_bare,abdv_H2l_bare
      real*8 sqrt2,pi,xnorm,dummy
      parameter (sqrt2=1.4142135623730950488d0)
      parameter (pi=3.1415926535897932385d0)
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c Use on-shell scheme for virtuals (see abdv_ghalf2lcr)
      integer ischeme
      parameter (ischeme=1)
c
      tmpn=dcmplx(0.d0,0.d0)
      tmpd=dcmplx(0.d0,0.d0)
      if(xmb2.ne.0.d0)then
        tmpn=tmpn+abdv_H2l_bare(xmh2,xmb2,xmu2,ischeme)
        tmpd=tmpd+abdv_H1l_bare(xmh2,xmb2)
      endif
      if(xmt2.ne.0.d0)then
        tmpn=tmpn+abdv_H2l_bare(xmh2,xmt2,xmu2,ischeme)
        tmpd=tmpd+abdv_H1l_bare(xmh2,xmt2)
      endif
      if(iabdvnotop.eq.0)then
        dummy=2*dreal(tmpn/tmpd)
        xnorm=1.d0
      elseif(iabdvnotop.eq.1)then
        dummy=2*dreal(conjg(tmpn)*tmpd)
        xnorm=xmh2/(256*sqrt2*pi**2*(4*pi)**2)
      endif
      abdv_f1sv_nonpi=xnorm*dummy
      return
      end


      function abdv_H2l_bare(xmh2,xmq2,xmu2,ischeme)
c Eq.(24) of 1111.2854, for one heavy quark flavour of mass squared xmq2.
c Note that the +h.c. in eq.(24) is incorrect, at least if one uses it
c in eq.(13) as done here. This has been confirmed by Degrassi
      implicit none
      complex*16 abdv_H2l_bare
      real*8 xmh2,xmq2,xmu2
      integer ischeme
      complex*16 xq,abdv_sqrt,tmp,abdv_ghalf2lcr,abdv_ghalf2lca
      real*8 yq,vca,vcf,vtf,lambdaq
      parameter (vca=3.d0)
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=0.5d0)
      parameter (lambdaq=1d0)
c
      yq=xmq2/xmh2
      xq=(abdv_sqrt(1-4*yq)-1)/(abdv_sqrt(1-4*yq)+1)
      tmp=vcf*abdv_ghalf2lcr(xq,xmq2,xmu2,ischeme) + 
     #    vca*abdv_ghalf2lca(xq)
      abdv_H2l_bare=vtf*lambdaq*tmp
      return
      end


      function abdv_ghalf2lca(x)
c Eq.(28) of hep-ph/0611266
      implicit none
      complex*16 abdv_ghalf2lca,x
      complex*16 tmp,abdv_H2slanted
      real*8 zeta3
      parameter (zeta3=1.2020569031595942854d0)
      integer izero,ione
      parameter (izero=0)
      parameter (ione=1)
      double complex HPL3
      external HPL3
c
      tmp=3+x*(1+8*x+3*x**2)/(x-1)**3*HPL3(izero,izero,izero,x)-
     #    2*(1+x)**2/(x-1)**2*abdv_H2slanted(x)+zeta3-
     #    HPL3(ione,izero,izero,x)
      abdv_ghalf2lca=4*x/(x-1)**2 * tmp
      return
      end


      function abdv_ghalf2lcr(x,xmq2,xmu2,ischeme)
c Eq.(11) [for MSbar] or eq.(15) [for on-shell] of hep-ph/0611266,
c since g_{1/2}^{(2l,CR)=F_{1/2}^{(2l)} as written after eq.(27).
c Call with
c  ischeme=1 -> on-shell
c  ischeme=2 -> MSbar
      implicit none
      complex*16 abdv_ghalf2lcr,x
      real*8 xmq2,xmu2
      integer ischeme
      complex*16 tmp,abdv_fhalf2la,abdv_fhalf2lb
c
      if(ischeme.eq.1)then
        tmp=abdv_fhalf2la(x)+4/3.d0*abdv_fhalf2lb(x)
      elseif(ischeme.eq.2)then
        tmp=abdv_fhalf2la(x)+log(xmq2/xmu2)*abdv_fhalf2lb(x)
      else
        write(*,*)'Error in abdv_ghalf2lcr: unknown scheme',ischeme
        stop
      endif
      abdv_ghalf2lcr=tmp
      return
      end


      function abdv_fhalf2la(x)
c Eq.(12) of hep-ph/0611266
      implicit none
      complex*16 abdv_fhalf2la,x
      complex*16 tmp,abdv_H1slanted
      real*8 pi,zeta2,zeta3
      parameter (pi=3.1415926535897932385d0)
      parameter (zeta2=pi**2/6.d0)
      parameter (zeta3=1.2020569031595942854d0)
      integer izero,ione,imone
      parameter (izero=0)
      parameter (ione=1)
      parameter (imone=-1)
      double complex HPL1,HPL2,HPL3
      external HPL1,HPL2,HPL3
c
      tmp=36/(x-1)**2-4*(1-14*x+x**2)/(x-1)**4*zeta3-
     #    4*(1+x)/(x-1)**3*HPL1(izero,x)-
     #    8*(1+9*x+x**2)/(x-1)**4*HPL2(izero,izero,x)+
     #    2*(3+25*x-7*x**2+3*x**3)/(x-1)**5*
     #      HPL3(izero,izero,izero,x)+
     #    4*(1+2*x+x**2)/(x-1)**4*
     #      ( zeta2*HPL1(izero,x)+
     #        4*HPL3(izero,imone,izero,x)-
     #        HPL3(izero,ione,izero,x) )+
     #    4*(5-6*x+5*x**2)/(x-1)**4*HPL3(ione,izero,izero,x)-
     #    8*(1+x+x**2+x**3)/(x-1)**5*abdv_H1slanted(x)
      abdv_fhalf2la=x*tmp
      return
      end


      function abdv_fhalf2lb(x)
c Eq.(12) of hep-ph/0611266
      implicit none
      complex*16 abdv_fhalf2lb,x
      complex*16 tmp
      integer izero
      parameter (izero=0)
      double complex HPL1,HPL2
      external HPL1,HPL2
c
      tmp=-12/(x-1)**2-6*(1+x)/(x-1)**3*HPL1(izero,x)+
     #    6*(1+6*x+x**2)/(x-1)**4*HPL2(izero,izero,x)
      abdv_fhalf2lb=x*tmp
      return
      end


      function abdv_H1slanted(x)
c Eq.(14) of hep-ph/0611266
      implicit none
      complex*16 abdv_H1slanted,x
      complex*16 tmp
      real*8 pi,zeta2,zeta3
      parameter (pi=3.1415926535897932385d0)
      parameter (zeta2=pi**2/6.d0)
      parameter (zeta3=1.2020569031595942854d0)
      integer izero,ione,imone
      parameter (izero=0)
      parameter (ione=1)
      parameter (imone=-1)
      double complex HPL1,HPL2,HPL4
      external HPL1,HPL2,HPL4
c
      tmp=9*zeta2**2/10.d0+2*zeta3*HPL1(izero,x)+
     #    zeta2*HPL2(izero,izero,x)+
     #    HPL4(izero,izero,izero,izero,x)/4.d0+
     #    7/2.d0*HPL4(izero,ione,izero,izero,x)-
     #    2*HPL4(izero,imone,izero,izero,x)+
     #    4*HPL4(izero,izero,imone,izero,x)-
     #    HPL4(izero,izero,ione,izero,x)
      abdv_H1slanted=tmp
      return
      end


      function abdv_H2slanted(x)
c Eq.(30) of hep-ph/0611266
      implicit none
      complex*16 abdv_H2slanted,x
      complex*16 tmp
      real*8 pi,zeta2,zeta3
      parameter (pi=3.1415926535897932385d0)
      parameter (zeta2=pi**2/6.d0)
      parameter (zeta3=1.2020569031595942854d0)
      integer izero,ione,imone
      parameter (izero=0)
      parameter (ione=1)
      parameter (imone=-1)
      double complex HPL1,HPL2,HPL3,HPL4
      external HPL1,HPL2,HPL3,HPL4
c
      tmp=4*zeta2**2/5.d0+2*zeta3+3*zeta3/2.d0*HPL1(izero,x)+
     #    3*zeta3*HPL1(ione,x)+zeta2*HPL2(ione,izero,x)+
     #    (1+2*zeta2)/4.d0*HPL2(izero,izero,x)-
     #    2*HPL3(ione,izero,izero,x)+
     #    HPL4(izero,izero,imone,izero,x)+
     #    HPL4(izero,izero,izero,izero,x)/4.d0+
     #    2*HPL4(ione,izero,imone,izero,x)-
     #    HPL4(ione,izero,izero,izero,x)
      abdv_H2slanted=tmp
      return
      end


      function abdv_Rgg(s,t,u,xi,y,xmh2,xmb2,xmt2)
c Eq.(25) of 1111.2854, with normalization consistent with f2body_gg,
c hence including the prefactors of eqs.(14) and (15) in 1111.2854
      implicit none
      real*8 abdv_Rgg,s,t,u,xi,y,xmh2,xmb2,xmt2
      integer i
      complex*16 abdv_A2_bare,abdv_A4_bare,tmp(4),tmpt(4)
      real*8 xs,xt,xu,xres,sqrt2,pi,xnorm,xfks
      parameter (sqrt2=1.4142135623730950488d0)
      parameter (pi=3.1415926535897932385d0)
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c
      xnorm=3*xmh2**4/(sqrt2*pi*(4*pi)**3)
c xkfs=xi**2*(1-y**2)/(s*t*u); these factors appear in eqs.(14) and (15)
c of 1111.2854
      xfks=4*(1-xi)/(xmh2*s**2)
c Mandelstam invariants as defined in 1111.2854, after eq.(17).
c NOTE THE EXCHANGE u<->t WRT WHAS IS DONE IN INVAR
      xs=s
      xt=u
      xu=t
      do i=1,4
        tmp(i)=dcmplx(0.d0,0.d0)
        tmpt(i)=dcmplx(0.d0,0.d0)
      enddo
      if(xmb2.ne.0.d0)then
        tmp(1)=tmp(1)+abdv_A2_bare(xs,xt,xu,xmh2,xmb2)
        tmp(2)=tmp(2)+abdv_A2_bare(xu,xs,xt,xmh2,xmb2)
        tmp(3)=tmp(3)+abdv_A2_bare(xt,xu,xs,xmh2,xmb2)
        tmp(4)=tmp(4)+abdv_A4_bare(xs,xt,xu,xmh2,xmb2)
      endif
      if(xmt2.ne.0.d0)then
        tmpt(1)=tmpt(1)+abdv_A2_bare(xs,xt,xu,xmh2,xmt2)
        tmpt(2)=tmpt(2)+abdv_A2_bare(xu,xs,xt,xmh2,xmt2)
        tmpt(3)=tmpt(3)+abdv_A2_bare(xt,xu,xs,xmh2,xmt2)
        tmpt(4)=tmpt(4)+abdv_A4_bare(xs,xt,xu,xmh2,xmt2)
        tmp(1)=tmp(1)+tmpt(1)
        tmp(2)=tmp(2)+tmpt(2)
        tmp(3)=tmp(3)+tmpt(3)
        tmp(4)=tmp(4)+tmpt(4)
      endif
      xres=0.d0
      do i=1,4
        if(iabdvnotop.eq.0)then
          xres=xres+abs(tmp(i))**2
        elseif(iabdvnotop.eq.1)then
          xres=xres+abs(tmp(i))**2-abs(tmpt(i))**2
        endif
      enddo
      abdv_Rgg=xnorm*xfks*xres/(2*s)
      return
      end


      subroutine abdv_Rqq(s,t,u,xi,y,xmh2,xmb2,xmt2,res)
c Eq.(21) of 1111.2854, with normalization consistent with f2body_qq
      implicit none
      real*8 s,t,u,xi,y,xmh2,xmb2,xmt2,res(4)
      integer i
      complex*16 abdv_Aqq
      real*8 xs,xt,xu,sqrt2,pi,zero,xnorm,xfks
      parameter (sqrt2=1.4142135623730950488d0)
      parameter (pi=3.1415926535897932385d0)
      parameter (zero=0.d0)
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c
      xnorm=8*xmh2**2/(sqrt2*pi*(4*pi)**3*9)
c At variance with eq.(21) of 1111.2854, f2body_qq contains the FKS
c damping factor
      xfks=xi**2*(1-y**2) * (t**2+u**2)/(s*(t+u)**2)
c Mandelstam invariants as defined in 1111.2854, after eq.(17).
c NOTE THE EXCHANGE u<->t WRT WHAS IS DONE IN INVAR
      xs=s
      xt=u
      xu=t
      if(iabdvnotop.eq.0)then
        res(1)=abs(abdv_Aqq(xs,xt,xu,xmh2,xmb2,xmt2))**2
      elseif(iabdvnotop.eq.1)then
        res(1)=abs(abdv_Aqq(xs,xt,xu,xmh2,xmb2,xmt2))**2-
     #         abs(abdv_Aqq(xs,xt,xu,xmh2,zero,xmt2))**2
      endif
      res(2)=0.d0
      res(3)=res(1)
      res(4)=0.d0
      do i=1,4
        res(i)=res(i) * xnorm*xfks/(2*s)
      enddo
      return
      end


      subroutine abdv_Rqg(s,t,u,xi,y,xmh2,xmb2,xmt2,res)
c Eqs.(16) and (17) of 1111.2854, with normalization consistent with
c f2body_qg. Use has been made of eq.(30)
      implicit none
      real*8 s,t,u,xi,y,xmh2,xmb2,xmt2,res(4)
      integer i
      complex*16 abdv_Aqq
      real*8 xs,xt,xu,sqrt2,pi,zero,xnorm,xfksd,xfksr
      parameter (sqrt2=1.4142135623730950488d0)
      parameter (pi=3.1415926535897932385d0)
      parameter (zero=0.d0)
      integer iabdvnotop
      common/ciabdvnotop/iabdvnotop
c
      xnorm=-xmh2**2/(sqrt2*pi*(4*pi)**3*3)
c xkfsd=xi**2*(1-y**2)/t*(s**2+u**2)/(s+u)**2
c xkfsr=xi**2*(1-y**2)/u*(s**2+t**2)/(s+t)**2
c These factors appear in eqs.(16) and (17) of 1111.2854. 
c Note that according to 1111.2854, t<->u wrt our conventions,
c which is the reason for the forms used here
      xfksd=-2*xi*(1+y)/s * (s**2+u**2)/(s+u)**2
      xfksr=-2*xi*(1-y)/s * (s**2+t**2)/(s+t)**2
c Mandelstam invariants as defined in 1111.2854, after eq.(17).
c NOTE THE EXCHANGE u<->t WRT WHAS IS DONE IN INVAR
      xs=s
      xt=u
      xu=t
      if(iabdvnotop.eq.0)then
        res(1)=abs(abdv_Aqq(xu,xs,xt,xmh2,xmb2,xmt2))**2
      elseif(iabdvnotop.eq.1)then
        res(1)=abs(abdv_Aqq(xu,xs,xt,xmh2,xmb2,xmt2))**2-
     #         abs(abdv_Aqq(xu,xs,xt,xmh2,zero,xmt2))**2
      endif
      res(2)=res(1)
      if(iabdvnotop.eq.0)then
        res(3)=abs(abdv_Aqq(xt,xs,xu,xmh2,xmb2,xmt2))**2
      elseif(iabdvnotop.eq.1)then
        res(3)=abs(abdv_Aqq(xt,xs,xu,xmh2,xmb2,xmt2))**2-
     #         abs(abdv_Aqq(xt,xs,xu,xmh2,zero,xmt2))**2
      endif
      res(4)=res(3)
      do i=1,2
        res(i)=res(i) * xnorm*xfksd/(2*s)
      enddo
      do i=3,4
        res(i)=res(i) * xnorm*xfksr/(2*s)
      enddo
      return
      end


      function abdv_Aqq(s,t,u,xmh2,xmb2,xmt2)
c Eq.(29) of 1111.2854
      implicit none
      complex*16 abdv_Aqq
      real*8 s,t,u,xmh2,xmb2,xmt2
      complex*16 tmp,abdv_dhalf
      real*8 sq,tq,uq,yq,vtf,lambdaq
      parameter (vtf=0.5d0)
      parameter (lambdaq=1d0)
c
      tmp=dcmplx(0.d0,0.d0)
      if(xmb2.ne.0.d0)then
        yq=xmb2/xmh2
        sq=s/xmb2
        tq=t/xmb2
        uq=u/xmb2
        tmp=tmp+yq*abdv_dhalf(sq,tq,uq,xmh2,xmb2)
      endif
      if(xmt2.ne.0.d0)then
        yq=xmt2/xmh2
        sq=s/xmt2
        tq=t/xmt2
        uq=u/xmt2
        tmp=tmp+yq*abdv_dhalf(sq,tq,uq,xmh2,xmt2)
      endif
      abdv_Aqq=vtf*lambdaq*tmp
      return
      end


      function abdv_A2_bare(s,t,u,xmh2,xmq2)
c Eq.(26) of 1111.2854, for one heavy quark flavour of mass squared xmq2
      implicit none
      complex*16 abdv_A2_bare
      real*8 s,t,u,xmh2,xmq2
      complex*16 abdv_bhalf,tmp
      real*8 sq,tq,uq,yq,vtf,lambdaq
      parameter (vtf=0.5d0)
      parameter (lambdaq=1d0)
c
      yq=xmq2/xmh2
      sq=s/xmq2
      tq=t/xmq2
      uq=u/xmq2
      tmp=dcmplx(0.d0,0.d0)
      tmp=tmp+abdv_bhalf(sq,tq,uq,xmh2,xmq2)
      tmp=tmp+abdv_bhalf(sq,uq,tq,xmh2,xmq2)
      abdv_A2_bare=vtf*lambdaq*yq**2*tmp
      return
      end


      function abdv_A4_bare(s,t,u,xmh2,xmq2)
c Eq.(27) of 1111.2854, for one heavy quark flavour of mass squared xmq2
      implicit none
      complex*16 abdv_A4_bare
      real*8 s,t,u,xmh2,xmq2
      complex*16 abdv_chalf,tmp
      real*8 sq,tq,uq,yq,vtf,lambdaq
      parameter (vtf=0.5d0)
      parameter (lambdaq=1d0)
c
      yq=xmq2/xmh2
      sq=s/xmq2
      tq=t/xmq2
      uq=u/xmq2
      tmp=dcmplx(0.d0,0.d0)
      tmp=tmp+abdv_chalf(sq,tq,uq,xmh2,xmq2)
      tmp=tmp+abdv_chalf(tq,uq,sq,xmh2,xmq2)
      tmp=tmp+abdv_chalf(uq,sq,tq,xmh2,xmq2)
      abdv_A4_bare=vtf*lambdaq*yq**2*tmp
      return
      end


      function abdv_bhalf(s,t,u,xmh2,xmq2)
c Eq.(22) of 0709.4227. Uses eq.(26) for B_{1/2}
      implicit none
      complex*16 abdv_bhalf
      real*8 s,t,u,xmh2,xmq2
      complex*16 xhalf,xs,xt,abdv_sqrt,abdv_H3,tmp
      real*8 yhalf
      integer izero
      parameter (izero=0)
      double complex HPL1,HPL2
      external HPL1,HPL2
c
      yhalf=xmq2/xmh2
      xhalf=(abdv_sqrt(1-4*yhalf)-1)/(abdv_sqrt(1-4*yhalf)+1)
      xs=(abdv_sqrt(1-4/s)-1)/(abdv_sqrt(1-4/s)+1)
      xt=(abdv_sqrt(1-4/t)-1)/(abdv_sqrt(1-4/t)+1)
c This is B_{1/2}
      tmp=s*(t-s)/(s+t)+2*(t*u**2+2*s*t*u)/(s+u)**2*
     #      ( abdv_sqrt(1-4*yhalf)*HPL1(izero,xhalf)-
     #        abdv_sqrt(1-4/t)*HPL1(izero,xt) ) -
     #    (1+t*u/s)*HPL2(izero,izero,xhalf) + HPL2(izero,izero,xs) -
     #    2*( 2*s**2/(s+u)**2-1-t*u/s )*
     #      (HPL2(izero,izero,xhalf)-HPL2(izero,izero,xt))+
     #    0.5d0*(t*u/s+3)*abdv_H3(s,u,t)-
     #    abdv_H3(t,s,u)
c Now all the other terms in eq.(22)
      tmp=tmp+s/4.d0*(HPL2(izero,izero,xhalf)-HPL2(izero,izero,xs))-
     #    (s/2-s**2/(s+u))*
     #        (HPL2(izero,izero,xhalf)-HPL2(izero,izero,xt))-
     #    s/8.d0*abdv_H3(s,u,t)+
     #    s/4.d0*abdv_H3(t,s,u)
      abdv_bhalf=tmp
      return
      end


      function abdv_chalf(s,t,u,xmh2,xmq2)
c Eq.(24) of 0709.4227. Uses eq.(27) for C_{1/2}
      implicit none
      complex*16 abdv_chalf
      real*8 s,t,u,xmh2,xmq2
      complex*16 xhalf,xs,xt,abdv_sqrt,abdv_H3,tmp
      real*8 yhalf
      integer izero
      parameter (izero=0)
      double complex HPL2
      external HPL2
c
      yhalf=xmq2/xmh2
      xhalf=(abdv_sqrt(1-4*yhalf)-1)/(abdv_sqrt(1-4*yhalf)+1)
      xs=(abdv_sqrt(1-4/s)-1)/(abdv_sqrt(1-4/s)+1)
      xt=(abdv_sqrt(1-4/t)-1)/(abdv_sqrt(1-4/t)+1)
c This is C_{1/2}
      tmp=-2*s-2*(HPL2(izero,izero,xhalf)-HPL2(izero,izero,xs))-
     #    abdv_H3(u,s,t)
c Now all the other terms in eq.(24)
      tmp=tmp+
     #    1/(2*yhalf)*(HPL2(izero,izero,xhalf)-HPL2(izero,izero,xs))+
     #    1/(4*yhalf)*abdv_H3(s,u,t)
      abdv_chalf=tmp
      return
      end


      function abdv_dhalf(s,t,u,xmh2,xmq2)
c Eq.(31) of 0709.4227. Uses eq.(33) for D_{1/2}
      implicit none
      complex*16 abdv_dhalf
      real*8 s,t,u,xmh2,xmq2
      complex*16 xhalf,xs,xt,abdv_sqrt,tmp
      real*8 yhalf
      integer izero
      parameter (izero=0)
      double complex HPL1,HPL2
      external HPL1,HPL2
c
      yhalf=xmq2/xmh2
      xhalf=(abdv_sqrt(1-4*yhalf)-1)/(abdv_sqrt(1-4*yhalf)+1)
      xs=(abdv_sqrt(1-4/s)-1)/(abdv_sqrt(1-4/s)+1)
      xt=(abdv_sqrt(1-4/t)-1)/(abdv_sqrt(1-4/t)+1)
c This is D_{1/2}
      tmp=4+4*s/(t+u)*
     #      ( abdv_sqrt(1-4*yhalf)*HPL1(izero,xhalf)-
     #        abdv_sqrt(1-4/s)*HPL1(izero,xs) ) +
     #    8/(t+u)*(HPL2(izero,izero,xhalf)-HPL2(izero,izero,xs))
c Now all the other terms in eq.(31)
      tmp=tmp-
     #    2*(HPL2(izero,izero,xhalf)-HPL2(izero,izero,xs))
      abdv_dhalf=tmp
      return
      end


      function abdv_H3(a,b,c)
c See below eq.(28) in 0709.4227
      implicit none
      complex*16 abdv_H3,tmp,W3fun
      real*8 a,b,c,sum
c
      sum=a+b+c
      tmp=-W3fun(b,a,c,sum)
      abdv_H3=tmp
      return
      end


      function W3fun(s,t,u,v)
c Eq.(A.17) of NPB297(1988)221
      implicit none
      complex*16 W3fun,I3fun,tmp
      real*8 s,t,u,v
c
      tmp=I3fun(s,t,u,v)-I3fun(s,t,u,s)-I3fun(s,t,u,u)
      W3fun=tmp
      return
      end


      function I3fun(s,t,u,v)
c Eq.(A.21) of NPB297(1988)221, and subsequent formulas. Note that,
c according to footnote 3 of 0709.4227, one has to put mf=1 here.
c Note that ddilogs are special cases of HPLs (see eq.(2.6) of 1106.5739 --
c in particular, Li2(z) = HPL2(0,1;z))
      implicit none
      complex*16 I3fun
      real*8 s,t,u,v
      complex*16 tmp,zic
      parameter (zic=(0.d0,1.d0))
      real*8 beta,gamma,a,r,cosphi,costh,sinphi,sinth,phi,th,ddlog2a,
     # getcos,getsinfromcos,xmf2,pi,zero
      parameter (xmf2=1.d0)
      parameter (pi=3.1415926535897932385d0)
      parameter (zero=0.d0)
      integer izero,ione
      parameter (izero=0)
      parameter (ione=1)
      double complex HPL2
      external HPL2
c
      beta=1/2.d0*(1+sqrt(1+4*t*xmf2/(u*s)))
      if(v.gt.0.d0.and.v.lt.(4*xmf2))then
        gamma=0.d0
        a=sqrt(4*xmf2/v-1)
        r=sqrt( (a**2+1)/(a**2+(2*beta-1)**2) )
        cosphi=r*(a**2+2*beta-1)/(1+a**2)
        cosphi=getcos(cosphi,'cll1  ')
        costh=r*(a**2-2*beta+1)/(1+a**2)
        costh=getcos(costh,'cll2  ')
c phi and theta are between 0 and pi (see eq.(A.22) of NPB297(1988)221)
        sinphi=getsinfromcos(cosphi,'cll1  ')
        sinth=getsinfromcos(costh,'cll2  ')
        phi=atan2(sinphi,cosphi)
        th=atan2(sinth,costh)
      else
        gamma=1/2.d0*(1+sqrt(1-4*xmf2/v))
        a=0.d0
        r=0.d0
        cosphi=0.d0
        costh=0.d0
        phi=0.d0
        th=0.d0
      endif
c
      if(v.lt.0.d0)then
        tmp=-HPL2(izero,ione,dcmplx(gamma/(gamma+beta-1),zero))+
     #      HPL2(izero,ione,dcmplx((gamma-1)/(gamma+beta-1),zero))+
     #      HPL2(izero,ione,dcmplx((beta-gamma)/beta,zero))-
     #      HPL2(izero,ione,dcmplx((beta-gamma)/(beta-1),zero))+
     #      0.5d0*(log(beta)**2-log(beta-1)**2)+
     #      log(gamma)*log((gamma+beta-1)/beta)+
     #      log(gamma-1)*log((beta-1)/(gamma+beta-1))
      elseif(v.gt.0.d0.and.v.lt.(4*xmf2))then
        tmp=2*ddlog2a(r,th)-2*ddlog2a(r,phi)+
     #      (phi-th)*(phi+th-pi)
      elseif(v.gt.(4*xmf2))then
        tmp=-HPL2(izero,ione,dcmplx(gamma/(gamma+beta-1),zero))+
     #      HPL2(izero,ione,dcmplx((gamma-1)/(gamma+beta-1),zero))+
     #      HPL2(izero,ione,dcmplx(gamma/(gamma-beta),zero))-
     #      HPL2(izero,ione,dcmplx((gamma-1)/(gamma-beta),zero))+
     #      log(gamma/(1-gamma))*log((gamma+beta-1)/(beta-gamma))-
     #      zic*pi*log((gamma+beta-1)/(beta-gamma))
      else
        tmp=0.d0
      endif
      I3fun=tmp * 2/(2*beta-1)
      return
      end


      function ddlog2a(x,th)
c This is Li2(x,theta) of NPB297(1988)221, after eq.(A.22).
c Its analytic form has been computed with mathematica:
c   Integrate[Log[1-2*a*z+z^2]/z,{z,0,x},Assumptions->a>-1&&a<1]
c Note that the resulting logarithmic part is zero after analytical
c simplification, and hence it does not appear below
      implicit none
      real*8 ddlog2a,x,th,a,s1ma
      complex*16 tmp,zic
      parameter (zic=(0.d0,1.d0))
      integer izero,ione
      parameter (izero=0)
      parameter (ione=1)
      double complex HPL2,zz
      external HPL2
c
      a=cos(th)
      s1ma=sqrt(1-a**2)
      tmp=-HPL2(izero,ione,-zic*x/(s1ma-zic*a))
     #    -HPL2(izero,ione,zic*x/(s1ma+zic*a))
      ddlog2a=-0.5d0*tmp
      return
      end


      function abdv_log(x)
c Utility function: log of a complex variable
      implicit none
      complex*16 abdv_log,x
      real*8 re,im,rho2,theta
c
      re=dreal(x)
      im=dimag(x)
      rho2=re**2+im**2
      theta=atan2(im,re)
      abdv_log=dcmplx(0.5d0*log(rho2),theta)
      return
      end


      function abdv_sqrt(x)
c Utility function: square root of a real variable of either sign
      implicit none
      complex*16 abdv_sqrt,tmp
      real*8 x,sx
c
      sx=sqrt(abs(x))
      if(x.ge.0.d0)then
        tmp=dcmplx(sx,0.d0)
      else
        tmp=dcmplx(0.d0,sx)
      endif
      abdv_sqrt=tmp
      return
      end


      function getcos(cosx,strid)
c Prevents numerical inaccuracies: returns cosx (meant to be a cosine) 
c or \pm 1
      implicit none
      real*8 getcos,cosx,tmp,tiny
      parameter (tiny=1.d-6)
      character*6 strid
c
      if(abs(cosx).ge.(1+tiny))then
        write(*,*)'Error in getcos',cosx
        write(*,*)' Called: ',strid
        stop
      elseif(abs(cosx).lt.(1+tiny).and.abs(cosx).ge.1.d0)then
        tmp=dsign(1.d0,cosx)
      else
        tmp=cosx
      endif
      getcos=tmp
      return
      end


      function getsinfromcos(cosx,strid)
c Prevents numerical inaccuracies: given cosx (meant to be a cosine) 
c returns the corresponding sine (with positive sign)
      implicit none
      real*8 getsinfromcos,cosx,tmp,tiny
      parameter (tiny=1.d-6)
      character*6 strid
c
      if(abs(cosx).ge.(1+tiny))then
        write(*,*)'Error in getsinfromcos',cosx
        write(*,*)' Called: ',strid
        stop
      elseif(abs(cosx).lt.(1+tiny).and.abs(cosx).ge.1.d0)then
        tmp=0.d0
      else
        tmp=sqrt(1-cosx**2)
      endif
      getsinfromcos=tmp
      return
      end
