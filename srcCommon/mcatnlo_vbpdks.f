c
c
c Begin of dkswrap.f (WZ)
c
c
c
      function dkswzborn(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
c Returns spin-averaged Born, with the same normalization as vbborn().
c Coincides with vbborn() in the case of SM couplings
      implicit none
      real*8 dkswzborn,zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q
      character*2 str
      real*8 convfac,sig,phil,cthl,phir,cthr,v1a,v1b,v1c,v3a,
     #  v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'born'
      call dkswrap(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c dkswrap has units GeV, convert to TeV before returning. Born amplitude
c squared has dimension 1/energy^2
      dkswzborn = convfac*sig
      return 
      end


      function dkswzborn2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13)
c Returns Born with spin correlations. 
c This is the fully undecayed matrix elements
      implicit none 
      real*8 dkswzborn2,q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13
      character*2 str
      real*8 convfac,sig
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'born'
      call dkswrap(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c dkswrap has units GeV, convert to TeV before returning. Born amplitude
c squared has dimension 1/energy^2
      dkswzborn2 = convfac*sig
      return 
      end


      function dkswzreal(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
c Returns spin-averaged NLO real corrections, times FKS damping factor,
c with the same normalization as vbfpp(). Coincides with vbfpp() in the 
c case of SM couplings
      implicit none
      real*8 dkswzreal,zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 convfac,sig,phil,cthl,phir,cthr,v1a,v1b,v1c,v3a,
     #  v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'real'
      call dkswrap(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c Apply FKS damping factor (4*tk*uk), and DKS prefactor ((xnc**2-1)/xnc).
c dkswrap has units GeV, convert to TeV before returning. Tree-level NLO
c amplitudes squared have dimension 1/energy^4
      dkswzreal = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sig
      return 
      end


      function dkswzreal2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13)
c Returns NLO real corrections with spin correlations, times FKS damping 
c factor. This is the fully undecayed matrix elements
      implicit none 
      real*8 dkswzreal2,q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 convfac,sig
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'real'
      call dkswrap(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c Apply FKS damping factor (4*tk*uk), and DKS prefactor ((xnc**2-1)/xnc).
c dkswrap has units GeV, convert to TeV before returning. Tree-level NLO
c amplitudes squared have dimension 1/energy^4
      dkswzreal2 = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sig
      return 
      end


      function dkswzvirt(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
c Returns spin-averaged finite virtual contribution, according to DKS
c conventions
      implicit none
      real*8 dkswzvirt,zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      complex*16 imag
      real*8 zero,qes2,muscloop
      integer abspart
      common/nparton/imag,zero,qes2,muscloop,abspart
      real*8 convfac,sig,phil,cthl,phir,cthr,v1a,v1b,v1c,v3a,
     #  v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'virt'
c Scale for DKS virtual computation; must be in GeV. Ellis-Sexton scale is
c set to zero at the beginning of the run (setdkspar)
      muscloop = sqrt(xmufct2*convfac)
      call dkswrap(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c Apply DKS prefactor ((xnc**2-1)/xnc); dkswrap has units GeV, convert 
c to TeV before returning. Virtual NLO amplitudes times Born amplitudes
c have dimension 1/energy^2
      dkswzvirt = (xnc**2-1)/xnc*convfac*sig
      return 
      end


      function dkswzvirt2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13)
c Returns finite virtual contribution with spin correlations. 
c This is the fully undecayed matrix elements
      implicit none
      real*8 dkswzvirt2,q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      complex*16 imag
      real*8 zero,qes2,muscloop
      integer abspart
      common/nparton/imag,zero,qes2,muscloop,abspart
      real*8 convfac,sig
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'virt'
c Scale for DKS virtual computation; must be in GeV. Ellis-Sexton scale is
c set to zero at the beginning of the run (setdkspar)
      muscloop = sqrt(xmufct2*convfac)
      call dkswrap(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,sig,me)
c Apply DKS prefactor ((xnc**2-1)/xnc); dkswrap has units GeV, convert 
c to TeV before returning. Virtual NLO amplitudes times Born amplitudes
c have dimension 1/energy^2
      dkswzvirt2 = (xnc**2-1)/xnc*convfac*sig
      return 
      end

      
      function dkswzsv(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
c Returns spin-averaged soft-virtual contribution, with the same normalization 
c as vb2(). Coincides with vb2() in the case of SM couplings
      implicit none
      real*8 dkswzsv,zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 betfac,delta
      common/betfac/betfac,delta
      real*8 born,btil,vcf,pi,sig,dksv,dkswzborn,dkswzvirt
      parameter (pi=3.14159265358979312D0)
c
      born=dkswzborn(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
      btil=betfac*sqrt(1-(sqrt(zmw2)+sqrt(zmz2))**2/s0)
      vcf=(xnc**2-1)/(2*xnc)
      sig=vcf*born*( log(s0/xmufct2)*(6+16*log(btil))+
     #               32*log(btil)**2-2*pi**2/3.d0-2.d0+
     #               2*log(s0/xmufct2)**2-
     #               6*log(s0/xmufct2) )
      dksv=dkswzvirt(zmw2,zmz2,s0,x,y,cth1,cth2,str,tk,uk,q1q,q2q) 
      dkswzsv=2*dksv+sig
      return
      end


      subroutine dkswrap(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  avg,sig,me)
c Main routine for calling DKS functions, using MC@NLO kinematics
c WZ process
      implicit none
      real*8 q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,sig
      character*2 str
      character*4 me
      logical avg
      real*8 xmom_cm(9,4)
      common/cxmomcm/xmom_cm
      real*8 xkm(3,3),xkm2(3,3),sw,ze2,gup,gdown,ez,gw
      common/weakcoup/xkm,xkm2,sw,ze2,gup,gdown,ez,gw
      character*2 prc,prdct
      common/process/prc,prdct
      character*3 drstr
      common/cdrstr/drstr
      logical bornonly,qqonly,gqonly
      common/subcont/ bornonly,qqonly,gqonly
      real * 8 wdtwmsb,wdtzmsb
      common/barewidths/wdtwmsb,wdtzmsb
      real * 8 vcZ2,acZ2,pwdtZ2
      common/Z2cplg/vcZ2,acZ2,pwdtZ2
      real * 8 zmw,zmz,zmw2,zmz2
      common/zmass/zmw,zmz,zmw2,zmz2
      real * 8 ga1,ga1mn,ga1pl,bw1fmpl,bw1fmmn,bw1delf,xm1low2,xm1upp2
      common/bw1cmm/ga1,ga1mn,ga1pl,bw1fmpl,bw1fmmn,bw1delf,
     #              xm1low2,xm1upp2
      real * 8 ga2,ga2mn,ga2pl,bw2fmpl,bw2fmmn,bw2delf,xm2low2,xm2upp2
      common/bw2cmm/ga2,ga2mn,ga2pl,bw2fmpl,bw2fmmn,bw2delf,
     #              xm2low2,xm2upp2
      real*8 s,q1c,q2c
      complex*16 a,b,zs
      common/qinne/a(7,7),b(7,7),zs(7,7),s(7,7)
      real*8 xi,ximax,convfac,pi,k1pk2sq,facavg,factor,dkswzborn,
     #  dkswzvirt,dkswzqqreal,dkswzqgreal,dkswzagreal
      real*8 sigmat(2,4),sigmatb(2,4),sigup(2,4),sigdw(2,4)
      real*8 sigmat7(6,4),sigmatb7(6,4),sigup7(6,4),sigdw7(6,4)
      real*8 pv1(4),pv2(4)
      real*8 pfac
      integer izero,i1,i2,ispinmax,ispin,ii
      parameter (convfac=1.d6)
      parameter (pi=3.14159265358979312D0)
      parameter (izero=0)
c Set xi smaller than ximax in order to compute DKS real matrix elements
c (for unknown reasons, if this is not the case DKS skip the computation).
c In the case of 6-body ME's, if xi>ximax some CPU time can be saved,
c since collinear reminders are not computed
      parameter (xi=0.1d0)
      parameter (ximax=1.d0)
c
c When averaging over spins, vector bosons must be on shell
      if(avg)then
        if( abs(q12-zmw2).gt.1.d-6*zmw2 .or.
     #      abs(q22-zmz2).gt.1.d-6*zmz2 )then
          write(*,*)
     #      'Error #1 in dkswrap: bosons must be on shell ',q12,q22
          stop
        endif
      endif
c reset the result array for 6-body ME's
      do i1 = 1,2
        do i2 = 1,4
          sigup(i1,i2) = 0.d0
          sigdw(i1,i2) = 0.d0
        enddo
      enddo
c reset the result array for 7-body ME's
      do i1 = 1,6
        do i2 = 1,4
          sigup7(i1,i2) = 0.d0
          sigdw7(i1,i2) = 0.d0
        enddo
      enddo
c
      if(avg)then
        ispinmax=36
        facavg=36.d0
        pfac=1.d0
      else
        ispinmax=1
        facavg=1.d0
c For each decayed vector boson, insert a factor
c   BR * (8pi) * (2pi) * Gamma/(pi*M^3) * q^4/BW(q,M)
c Use BR * Gamma = partial decay width, and compute such width
c at the LO in the SM
        pfac=wdtwmsb*16*pi*q12**2/zmw**3/
     #       ((q12-zmw2)**2+(zmw*ga1)**2) *
     #       pwdtZ2*16*pi*q22**2/zmz**3/
     #       ((q22-zmz2)**2+(zmz*ga2)**2)
      endif
c     
      do ispin=1,ispinmax
c steer the DKS 
        if (me.eq.'born') then
          bornonly = .true.
        else
          bornonly = .false.
        endif
c pass s to the DKS routines TMART6NEW and TMATR7NEW
        do ii=1,4
          pv1(ii) = xmom_cm(1,ii) 
          pv2(ii) = xmom_cm(2,ii)
        enddo
        sigmat(1,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
        sigmatb(1,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
c pass sww
        do ii=1,4
          pv1(ii) = xmom_cm(4,ii) 
          pv2(ii) = xmom_cm(5,ii)
        enddo
        sigmat(2,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
        sigmatb(2,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
c needs to be populated 4 times because the way sigmat is passed to 
c sigtyp in wzanDKS.f
        do ii=1,4
          sigmat7(1,ii) = sigmat(1,1)
          sigmat7(2,ii) = sigmat(2,1)
          sigmatb7(1,ii) = sigmatb(1,1)
          sigmatb7(2,ii) = sigmatb(2,1)
        enddo
c get DKS kinematics from MC@NLO
        if(avg)then
          call getangles(ispin,phil,cthl,phir,cthr)
          call invar(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #      str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,izero)
        endif
        if(drstr.eq.'dir')then
          call fillqinne(s0,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #                   q12,q22,x)
        elseif(drstr.eq.'ref')then
          q1c=q12+q22-s0-tk-q1q
          q2c=q12+q22-s0-uk-q2q
          call reflex_xmom_cm()
          call fillqinne(s0,uk,tk,q2c,q1c,v1b,v1a,v1c,v3b,v3a,v3c,v13,
     #                   q12,q22,x)
        else
          write(*,*)'Error in dkswrap: unknown drstr',drstr
          stop
        endif
        if (me.eq.'born'.or.me.eq.'virt') then 
           call tmatr6new(xi,ximax,sigmat,sigmatb)
        elseif (me.eq.'real') then
           call tmatr7new(xi,ximax,sigmat7,sigmatb7)
        else
          write(*,*)'Error in dkswrap: unknown me ',me
          stop
        endif
c
        do i1 = 1,6
          do i2 = 1,4
            if (i1.le.2) then
              sigup(i1,i2) = sigup(i1,i2) + sigmat(i1,i2)/facavg
              sigdw(i1,i2) = sigdw(i1,i2) + sigmatb(i1,i2)/facavg
            endif
            sigup7(i1,i2) = sigup7(i1,i2) + sigmat7(i1,i2)/facavg
            sigdw7(i1,i2) = sigdw7(i1,i2) + sigmatb7(i1,i2)/facavg
          enddo
        enddo
c end loop over lepton configurations (computation of the spin-averaged
c cross sections)
      enddo
c The DKS code used here is normalized as follows (see setdkspar())
c   M^{tree}(eq.(3.19)) = BR * BR * |V_ud|^2 * e^4 * DKS_code
c with eq.(3.19) in hep-ph/9803250. When using the narrow-width
c approximation, M^{tree} fulfills eq.(3.9). Since each angular
c integral can be obtained as an averaged sum over spin, times (4*pi)
c (see Zolly's note), it follows that
c   M_undecayed = |V_ud|^2 * e^4 * (4*pi)^2 (spin-averaged DKS_code)
c The MC@NLO code is normalized as follows:
c   M_undecayed = |V_ud|^2 * e^4 / (8*sinthw**2) * MC@NLO_code
c with the factor |V_ud|^2 / (8*sinthw**2) included in strfun() [the
c denominator is 1/gw**2]. Therefore
c   (spin-averaged DKS_code) = 1/(4*pi)^2 * 1/(8*sinthw**2) MC@NLO_code
c We remove this factor here, in order to have the DKS result returned
c by dkswrap with the same normalization as the MC@NLO code
      factor = (1.d0/(4*pi))**2*1.d0/(8*sw**2)
      if (me.eq.'born'.or.me.eq.'virt') then 
        dkswzborn = sigup(2,2)
        dkswzvirt = sigup(2,1) 
      elseif (me.eq.'real') then
        dkswzqqreal = sigup7(2,1)
        dkswzqgreal = sigdw7(6,1)
        dkswzagreal = sigdw7(3,1)
      endif
c For Born, returns always that of the qq process, regardless of the value
c of prc, since it is needed for the computation of the soft and collinear
c limits
      if (me.eq.'born') then
        sig = dkswzborn/factor*pfac
      elseif (me.eq.'virt') then
        if(prc.ne.'qq')then
          write(*,*)'Error in dkswrap: no such process for SV ',prc
          stop
        endif
        sig = dkswzvirt/factor*pfac
      elseif (me.eq.'real') then
        if(prc.eq.'qq')then
          sig = dkswzqqreal/factor*pfac
        elseif(prc.eq.'qg')then
          sig = dkswzqgreal/factor*pfac
        elseif(prc.eq.'ag')then
          sig = dkswzagreal/factor*pfac
        else
          write(*,*)'Error in dkswrap: no such process for real ',prc
          stop
        endif
      endif
c Reinstate original kinematics if reflected
      if(drstr.eq.'ref')
     #     call invar(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #        str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,izero)
c
      return
      end


      subroutine getangles(ispin,phil,cthl,phir,cthr)
c Given 1<=ispin<=36, returns the angles phil,cthl,phir,cthr which 
c have to be given in input to invar() in order to compute the averaged
c cross section.
      implicit none
      integer ispin,iv1conf,iv2conf
      real*8 phil,cthl,phir,cthr,pi
      parameter (pi=3.1415926535897932d0)
c
c decide which angles should be used
c for leptons from v1 and v2
c v?conf convention:
c i l1 l2
c 1 x -x
c 2 -x x
c 3 y -y
c 4 -y y
c 5 z -z
c 6 -z z
c- cos th from z axis
c- phi=0 at +x axis
c
      iv1conf = mod(ispin-1,6)+ 1  !  1..6
      iv2conf = (ispin-1)/6 + 1  !  1..6
c
      if (iv1conf.eq.1) then
         phil = 0.
         cthl = 0.
      else if (iv1conf.eq.2) then
         phil = pi
         cthl = 0.
      else if (iv1conf.eq.3) then
         phil = -pi/2.d0
         cthl = 0.
      else if (iv1conf.eq.4) then
         phil = +pi/2.d0
         cthl = 0.
      else if (iv1conf.eq.5) then
         phil = 0.
         cthl = 1.
      else if (iv1conf.eq.6) then
         phil = 0.
         cthl = -1.
      endif
c
      if (iv2conf.eq.1) then
         phir = 0.
         cthr = 0.
      else if (iv2conf.eq.2) then
         phir = pi
         cthr = 0.
      else if (iv2conf.eq.3) then
         phir = -pi/2.d0
         cthr = 0.
      else if (iv2conf.eq.4) then
         phir = +pi/2.d0
         cthr = 0.
      else if (iv2conf.eq.5) then
         phir = 0.
         cthr = 1.
      else if (iv2conf.eq.6) then
         phir = 0.
         cthr = -1.
      endif
c
      return
      end


      subroutine reflex_xmom_cm()
      implicit none
      integer i,j
      real*8 xmom_cm(9,4)
      common/cxmomcm/xmom_cm
c
      do i=3,9
        do j=1,3
          xmom_cm(i,j)=-xmom_cm(i,j)
        enddo
      enddo
      return
      end
c
c
c End of dkswrap.f (WZ)
c
c
c
c
c Begin of dkswrap.f (WW)
c
c
      subroutine dkswwborn(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,xsecup,xsecdo) 
c Returns spin-averaged Born, with the same normalization as vbborn().
c Coincides with vbborn() in the case of SM couplings
      implicit none
      real*8 zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q,xsecup,xsecdo
      character*2 str
      real*8 convfac,sigup,sigdo,phil,cthl,phir,cthr,v1a,v1b,v1c,
     #  v3a,v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'born'
      call dkswrap2(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  avg,sigup,sigdo,me)
c dkswrap2 has units GeV, convert to TeV before returning. Born amplitude
c squared has dimension 1/energy^2
      xsecup = convfac*sigup
      xsecdo = convfac*sigdo
      return 
      end


      subroutine dkswwborn2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  xsecup,xsecdo)
c Returns Born with spin correlations. 
c This is the fully undecayed matrix elements
      implicit none 
      real*8 q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,xsecup,xsecdo
      character*2 str
      real*8 convfac,sigup,sigdo
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'born'
      call dkswrap2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,
     #  sigup,sigdo,me)
c dkswrap2 has units GeV, convert to TeV before returning. Born amplitude
c squared has dimension 1/energy^2
      xsecup = convfac*sigup
      xsecdo = convfac*sigdo
      return 
      end


      subroutine dkswwreal(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,xsecup,xsecdo) 
c Returns spin-averaged NLO real corrections, times FKS damping factor,
c with the same normalization as vbfpp(). Coincides with vbfpp() in the 
c case of SM couplings
      implicit none
      real*8 zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q,xsecup,xsecdo
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 convfac,sigup,sigdo,phil,cthl,phir,cthr,v1a,v1b,v1c,v3a,
     #  v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'real'
      call dkswrap2(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  avg,sigup,sigdo,me)
c Apply FKS damping factor (4*tk*uk), and DKS prefactor ((xnc**2-1)/xnc).
c dkswrap2 has units GeV, convert to TeV before returning. Tree-level NLO
c amplitudes squared have dimension 1/energy^4
      xsecup = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sigup
      xsecdo = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sigdo
      return 
      end


      subroutine dkswwreal2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  xsecup,xsecdo)
c Returns NLO real corrections with spin correlations, times FKS damping 
c factor. This is the fully undecayed matrix elements
      implicit none 
      real*8 q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,xsecup,xsecdo
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 convfac,sigup,sigdo
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'real'
      call dkswrap2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,
     #  sigup,sigdo,me)
c Apply FKS damping factor (4*tk*uk), and DKS prefactor ((xnc**2-1)/xnc).
c dkswrap2 has units GeV, convert to TeV before returning. Tree-level NLO
c amplitudes squared have dimension 1/energy^4
      xsecup = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sigup
      xsecdo = 4*tk*uk*(xnc**2-1)/xnc*convfac**2*sigdo
      return 
      end


      subroutine dkswwvirt(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,xsecup,xsecdo) 
c Returns spin-averaged finite virtual contribution, according to DKS
c conventions
      implicit none
      real*8 zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q,xsecup,xsecdo
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      complex*16 imag
      real*8 zero,qes2,muscloop
      integer abspart
      common/nparton/imag,zero,qes2,muscloop,abspart
      real*8 convfac,sigup,sigdo,phil,cthl,phir,cthr,v1a,v1b,v1c,v3a,
     #  v3b,v3c,v13
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .true.
      me = 'virt'
c Scale for DKS virtual computation; must be in GeV. Ellis-Sexton scale is
c set to zero at the beginning of the run (setdkspar)
      muscloop = sqrt(xmufct2*convfac)
      call dkswrap2(zmw2,zmz2,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  avg,sigup,sigdo,me)
c Apply DKS prefactor ((xnc**2-1)/xnc); dkswrap2 has units GeV, convert 
c to TeV before returning. Virtual NLO amplitudes times Born amplitudes
c have dimension 1/energy^2
      xsecup = (xnc**2-1)/xnc*convfac*sigup
      xsecdo = (xnc**2-1)/xnc*convfac*sigdo
      return 
      end


      subroutine dkswwvirt2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  xsecup,xsecdo)
c Returns finite virtual contribution with spin correlations. 
c This is the fully undecayed matrix elements
      implicit none
      real*8 q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,xsecup,xsecdo
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      complex*16 imag
      real*8 zero,qes2,muscloop
      integer abspart
      common/nparton/imag,zero,qes2,muscloop,abspart
      real*8 convfac,sigup,sigdo
      character*4 me
      logical avg
      parameter (convfac=1.d6)
c
      avg = .false.
      me = 'virt'
c Scale for DKS virtual computation; must be in GeV. Ellis-Sexton scale is
c set to zero at the beginning of the run (setdkspar)
      muscloop = sqrt(xmufct2*convfac)
      call dkswrap2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,avg,
     #  sigup,sigdo,me)
c Apply DKS prefactor ((xnc**2-1)/xnc); dkswrap2 has units GeV, convert 
c to TeV before returning. Virtual NLO amplitudes times Born amplitudes
c have dimension 1/energy^2
      xsecup = (xnc**2-1)/xnc*convfac*sigup
      xsecdo = (xnc**2-1)/xnc*convfac*sigdo
      return 
      end

      
      subroutine dkswwsv(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,xsecup,xsecdo) 
c Returns spin-averaged soft-virtual contribution, with the same normalization 
c as vb2(). Coincides with vb2() in the case of SM couplings
      implicit none
      real*8 zmw2,zmz2,s0,x,y,cth1,cth2,tk,uk,q1q,q2q,xsecup,xsecdo
      character*2 str
      real * 8 sh,xmufct2,xmuren2,as,xnc
      common/fixvar/sh,xmufct2,xmuren2,as,xnc
      real*8 betfac,delta
      common/betfac/betfac,delta
      real*8 btil,vcf,pi,dksvup,dksvdo,bornup,borndo,sigfact,
     # sigup,sigdo
      parameter (pi=3.14159265358979312D0)
c
      call dkswwborn(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,bornup,borndo) 
      btil=betfac*sqrt(1-(sqrt(zmw2)+sqrt(zmz2))**2/s0)
      vcf=(xnc**2-1)/(2*xnc)
      sigfact=vcf*( log(s0/xmufct2)*(6+16*log(btil))+
     #              32*log(btil)**2-2*pi**2/3.d0-2.d0+
     #              2*log(s0/xmufct2)**2-
     #              6*log(s0/xmufct2) )
      sigup=sigfact*bornup
      sigdo=sigfact*borndo
      call dkswwvirt(zmw2,zmz2,s0,x,y,cth1,cth2,str,
     #  tk,uk,q1q,q2q,dksvup,dksvdo)
      xsecup=2*dksvup+sigup
      xsecdo=2*dksvdo+sigdo
      return
      end


      subroutine dkswrap2(q12,q22,s0,x,y,cth1,cth2,phil,cthl,
     #  phir,cthr,str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #  avg,xsigup,xsigdo,me)
c Main routine for calling DKS functions, using MC@NLO kinematics
c WW process
      implicit none
      real*8 q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #  tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,xsigup,xsigdo
      character*2 str
      character*4 me
      logical avg
      real*8 xmom_cm(9,4)
      common/cxmomcm/xmom_cm
      real*8 xkm(3,3),xkm2(3,3),sw,ze2,gup,gdown,ez,gw
      common/weakcoup/xkm,xkm2,sw,ze2,gup,gdown,ez,gw
      character*2 prc,prdct
      common/process/prc,prdct
      character*3 drstr
      common/cdrstr/drstr
      logical bornonly,qqonly,gqonly
      common/subcont/ bornonly,qqonly,gqonly
      real * 8 wdtwmsb,wdtzmsb
      common/barewidths/wdtwmsb,wdtzmsb
      real * 8 vcZ2,acZ2,pwdtZ2
      common/Z2cplg/vcZ2,acZ2,pwdtZ2
      real * 8 zmw,zmz,zmw2,zmz2
      common/zmass/zmw,zmz,zmw2,zmz2
      real * 8 ga1,ga1mn,ga1pl,bw1fmpl,bw1fmmn,bw1delf,xm1low2,xm1upp2
      common/bw1cmm/ga1,ga1mn,ga1pl,bw1fmpl,bw1fmmn,bw1delf,
     #              xm1low2,xm1upp2
      real * 8 ga2,ga2mn,ga2pl,bw2fmpl,bw2fmmn,bw2delf,xm2low2,xm2upp2
      common/bw2cmm/ga2,ga2mn,ga2pl,bw2fmpl,bw2fmmn,bw2delf,
     #              xm2low2,xm2upp2
      real*8 s,q1c,q2c
      complex*16 a,b,zs
      common/qinne/a(7,7),b(7,7),zs(7,7),s(7,7)
      real*8 xi,ximax,convfac,pi,k1pk2sq,facavg,factor,dkswwbornup,
     #  dkswwborndo,dkswwvirtup,dkswwvirtdo,dkswwqqrealup,dkswwqqrealdo,
     #  dkswwqgrealup,dkswwqgrealdo,dkswwagrealup,dkswwagrealdo
      real*8 sigmat(2,4),sigmatb(2,4),sigup(2,4),sigdw(2,4)
      real*8 sigmat7(6,4),sigmatb7(6,4),sigup7(6,4),sigdw7(6,4)
      real*8 pv1(4),pv2(4)
      real*8 pfac
      integer izero,i1,i2,ispinmax,ispin,ii
      parameter (convfac=1.d6)
      parameter (pi=3.14159265358979312D0)
      parameter (izero=0)
c Set xi smaller than ximax in order to compute DKS real matrix elements
c (for unknown reasons, if this is not the case DKS skip the computation).
c In the case of 6-body ME's, if xi>ximax some CPU time can be saved,
c since collinear reminders are not computed
      parameter (xi=0.1d0)
      parameter (ximax=1.d0)
c
c When averaging over spins, vector bosons must be on shell
      if(avg)then
        if( abs(q12-zmw2).gt.1.d-6*zmw2 .or.
     #      abs(q22-zmz2).gt.1.d-6*zmz2 )then
          write(*,*)
     #      'Error #1 in dkswrap2: bosons must be on shell ',q12,q22
          stop
        endif
      endif
c reset the result array for 6-body ME's
      do i1 = 1,2
        do i2 = 1,4
          sigup(i1,i2) = 0.d0
          sigdw(i1,i2) = 0.d0
        enddo
      enddo
c reset the result array for 7-body ME's
      do i1 = 1,6
        do i2 = 1,4
          sigup7(i1,i2) = 0.d0
          sigdw7(i1,i2) = 0.d0
        enddo
      enddo
c
      if(avg)then
        ispinmax=36
        facavg=36.d0
        pfac=1.d0
      else
        ispinmax=1
        facavg=1.d0
c For each decayed vector boson, insert a factor
c   BR * (8pi) * (2pi) * Gamma/(pi*M^3) * q^4/BW(q,M)
c Use BR * Gamma = partial decay width, and compute such width
c at the LO in the SM
        pfac=wdtwmsb*16*pi*q12**2/zmw**3/
     #       ((q12-zmw2)**2+(zmw*ga1)**2) *
     #       wdtwmsb*16*pi*q22**2/zmz**3/
     #       ((q22-zmz2)**2+(zmz*ga2)**2)
      endif
c     
      do ispin=1,ispinmax
c steer the DKS 
        if (me.eq.'born') then
          bornonly = .true.
        else
          bornonly = .false.
        endif
c pass s to the DKS routines TMART6NEW and TMATR7NEW
        do ii=1,4
          pv1(ii) = xmom_cm(1,ii) 
          pv2(ii) = xmom_cm(2,ii)
        enddo
        sigmat(1,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
        sigmatb(1,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
c pass sww
        do ii=1,4
          pv1(ii) = xmom_cm(4,ii) 
          pv2(ii) = xmom_cm(5,ii)
        enddo
        sigmat(2,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
        sigmatb(2,1) = k1pk2sq(pv1(1),pv2(1),1) * convfac
c needs to be populated 4 times because the way sigmat is passed to 
c sigtyp in wwpackage.f (same as for WZ)
        do ii=1,4
          sigmat7(1,ii) = sigmat(1,1)
          sigmat7(2,ii) = sigmat(2,1)
          sigmatb7(1,ii) = sigmatb(1,1)
          sigmatb7(2,ii) = sigmatb(2,1)
        enddo
c get DKS kinematics from MC@NLO
        if(avg)then
          call getangles(ispin,phil,cthl,phir,cthr)
          call invar(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #      str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,izero)
        endif
        if(drstr.eq.'dir')then
          call fillqinne(s0,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,
     #                   q12,q22,x)
        elseif(drstr.eq.'ref')then
          q1c=q12+q22-s0-tk-q1q
          q2c=q12+q22-s0-uk-q2q
          call reflex_xmom_cm()
          call fillqinne(s0,uk,tk,q2c,q1c,v1b,v1a,v1c,v3b,v3a,v3c,v13,
     #                   q12,q22,x)
        else
          write(*,*)'Error in dkswrap2: unknown drstr',drstr
          stop
        endif
        if (me.eq.'born'.or.me.eq.'virt') then 
           call tmatr6new_ww(xi,ximax,sigmat,sigmatb)
        elseif (me.eq.'real') then
           call tmatr7new_ww(xi,ximax,sigmat7,sigmatb7)
        else
          write(*,*)'Error in dkswrap2: unknown me ',me
          stop
        endif
c
        do i1 = 1,6
          do i2 = 1,4
            if (i1.le.2) then
              sigup(i1,i2) = sigup(i1,i2) + sigmat(i1,i2)/facavg
              sigdw(i1,i2) = sigdw(i1,i2) + sigmatb(i1,i2)/facavg
            endif
            sigup7(i1,i2) = sigup7(i1,i2) + sigmat7(i1,i2)/facavg
            sigdw7(i1,i2) = sigdw7(i1,i2) + sigmatb7(i1,i2)/facavg
          enddo
        enddo
c end loop over lepton configurations (computation of the spin-averaged
c cross sections)
      enddo
c The relative normalizations of the DKS and MC@NLO code appear to be
c the same as in the case of WZ production, except for the fact that
c for WW MC@NLO does not include a factor 1/(8*sinthw**2) in the main 
c code as it does in the case of WZ production. Hence, while for WZ we had
c        factor = (1.d0/(4*pi))**2*1.d0/(8*sw**2)
c here only the pi-dependent term is necessary
      factor = (1.d0/(4*pi))**2*1.d0
      if (me.eq.'born'.or.me.eq.'virt') then 
        dkswwbornup = sigdw(1,2)
        dkswwborndo = sigdw(2,2)
        dkswwvirtup = sigdw(1,1) 
        dkswwvirtdo = sigdw(2,1) 
      elseif (me.eq.'real') then
        dkswwqqrealup = sigdw7(1,1)
        dkswwqqrealdo = sigdw7(2,1)
        dkswwqgrealup = sigdw7(3,1)
        dkswwqgrealdo = sigdw7(4,1)
        dkswwagrealup = sigdw7(5,1)
        dkswwagrealdo = sigdw7(6,1)
      endif
c For Born, returns always that of the qq process, regardless of the value
c of prc, since it is needed for the computation of the soft and collinear
c limits
      if (me.eq.'born') then
        xsigup = dkswwbornup/factor*pfac
        xsigdo = dkswwborndo/factor*pfac
      elseif (me.eq.'virt') then
        if(prc.ne.'qq')then
          write(*,*)'Error in dkswrap2: no such process for SV ',prc
          stop
        endif
        xsigup = dkswwvirtup/factor*pfac
        xsigdo = dkswwvirtdo/factor*pfac
      elseif (me.eq.'real') then
        if(prc.eq.'qq')then
          xsigup = dkswwqqrealup/factor*pfac
          xsigdo = dkswwqqrealdo/factor*pfac
        elseif(prc.eq.'qg')then
          xsigup = dkswwqgrealup/factor*pfac
          xsigdo = dkswwqgrealdo/factor*pfac
        elseif(prc.eq.'ag')then
          xsigup = dkswwagrealup/factor*pfac
          xsigdo = dkswwagrealdo/factor*pfac
        else
          write(*,*)'Error in dkswrap2: no such process for real ',prc
          stop
        endif
      endif
c Reinstate original kinematics if reflected
      if(drstr.eq.'ref')
     #     call invar(q12,q22,s0,x,y,cth1,cth2,phil,cthl,phir,cthr,
     #        str,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,izero)
c
      return
      end
c
c
c End of dkswrap.f (WW)
c
c
c
c
c Helper functions
c
c
      real*8 function mysqrt(a)
      real*8 a
      if (a.gt.0.) then
         mysqrt = sqrt(a)
      else
         write(6,*) 'mcatnlo_vbpdks:mysqrt() neg argument, returning 0'
      endif
      return 
      end
c
c
c End of helper functions
c
c
c
c
c Begin of fillqinne.f
c
c
      subroutine fillqinne(s0,tk,uk,q1q,q2q,v1a,v1b,v1c,
     #                     v3a,v3b,v3c,v13,xm12,xm22,x)
c Fills the DKS common block
c  COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
c and uses A, B, ZS, S to compute the matrix elements. By definition
c  A(i,j) = <i,j>
c  B(i,j) = [i,j]
c  S(i,j) = p_i.p_j
c  ZS(i,j)= ( S(i,j),0 )   [<==> Re(ZS(i,j))=S(i,j), Im(ZS(i,j))=0]
      implicit none
      real*8 s0,tk,uk,q1q,q2q,v1a,v1b,v1c,v3a,v3b,v3c,v13,xm12,xm22,x
c variables to transform the p1, p2 momenta to the convention of DKS that
c all momenta should be outgoing
      complex*16 cc,rc,zt,zt1,zt2,zt3
      integer i,j,i_xmom,j_xmom,iqinxmo(7)
      logical flfail,checks
      complex*16 xna_2,xnb_2
      complex*16 a,b,zs
      real*8 s
      common/qinne/a(7,7),b(7,7),zs(7,7),s(7,7)
      real*8 schk(7,7)
      real*8 cpsi,spsi
      common/cspsi/cpsi,spsi
      real*8 xmom_cm(9,4)
      common/cxmomcm/xmom_cm
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop, abspart
      real*8 pi,convfac,sconvfac
      integer wcharge
      common/wcharge/wcharge
      parameter (pi=3.1415926535897932d0)
      parameter( convfac = 1.d6 )
      parameter ( sconvfac = 1.d3)
      parameter ( cc = (0.d0,1.d0) )
      parameter ( rc = (1.d0,0.d0) )
      parameter ( checks = .false. )
c
c Map MC@NLO momenta into DKS momenta, using
c   p_DKS(i) = p_MC@NLO(IQINXMO(i))
c The conventions for particle identities of DKS are given in terms of
c those relevant to processes initiated by a quark,antiquark pair; DKS
c routines internally do the necessary permutations in order to compute
c all necessary physical processes. The conventions are as follows:
c [actually, in hep-ph/9803250 the quark and antiquark are exchanged;
c thus, we take here the complex conjugates. These assignments are likely
c related to those in dkswrap*, in particular the choice of sigup or sigdw]
c
c  d(1)ubar(2) -> W^- Z g(7) -> e^-(3)nubar_e(4) mu^+(5) mu^-(6) g(7)
c  u(1)dbar(2) -> W^+ Z g(7) -> e^+(3)nu_e(4)    mu^-(5) mu^+(6) g(7)
c  q(1)qbar(2) -> W^+ W^- g(7) -> e^+(3)nu_e(4) mu^-(5)nubar_mu(6) g(7)
c
c while the conventions of MC@NLO are
c
c  d(1)ubar(2) -> W^-(4) Z(5) g(3)   -> e^-(6)nubar_e(7) mu^-(8) mu^+(9) g(3)
c  u(1)dbar(2) -> W^+(4) Z(5) g(3)   -> e^+(6)nu_e(7)    mu^-(8) mu^+(9) g(3)
c  q(1)qbar(2) -> W^+(4) W^-(5) g(3) -> e^+(6)nu_e(7) mu^-(8)nubar_mu(9) g(3)
c
      DATA IQINXMO(1) /1/
      DATA IQINXMO(2) /2/
      DATA IQINXMO(3) /6/
      DATA IQINXMO(4) /7/
      DATA IQINXMO(5) /0/
      DATA IQINXMO(6) /0/
      DATA IQINXMO(7) /3/
c
      if (wcharge.eq.1.or.wcharge.eq.2) then
         IQINXMO(5)=8
         IQINXMO(6)=9
      elseif (wcharge.eq.-1) then
         IQINXMO(5)=9
         IQINXMO(6)=8
      else
         write(*,*)'Error in fillqinne: unknown W charge',wcharge
         stop
      endif
c
c Set cosine and sine of psi (used by xnaxnb2)
      call get_psi1(flfail)
      if(flfail)then
        cpsi = cos(pi)
        spsi = sin(pi)
      endif
c
c Fill s(i,j) array for DKS; i and j follow DKS conventions, and thus
c incoming partons (i,j=1,2) have negative energies. This has been taken
c into account when expressing them in terms of invariants returned by invar()
      s(1,2) = s0*convfac
      s(1,3) = v1a*convfac
      s(1,4) = -(xm12-q1q+v1a)*convfac
      s(1,7) = tk*convfac
      s(2,3) = v1b*convfac
      s(2,4) = -(-xm22+q2q+s0+uk+v1b)*convfac
      s(2,7) = uk*convfac
      s(3,4) = xm12*convfac 
      s(3,7) = (-xm12-v1a-v1b-v1c)*convfac
      s(4,7) = (xm12-q1q+q2q-tk+v1a+v1b+v1c)*convfac
      s(5,6) = xm22*convfac
      if (wcharge.eq.1.or.wcharge.eq.2) then
        s(1,6) = -(-xm12+q1q+s0+tk+v3a)*convfac
        s(1,5) = v3a*convfac
        s(2,6) = -(xm22-q2q+v3b)*convfac
        s(2,5) = v3b*convfac
        s(3,6) = (v1c-v13)*convfac
        s(3,5) = v13*convfac
        s(4,6) = (-xm12-xm22+s0+tk+uk+v13-v1c-v3c)*convfac
        s(4,5) = (v3c-v13)*convfac
        s(6,7) = (xm22+q1q-q2q-uk+v3a+v3b+v3c)*convfac
        s(5,7) = (-xm22-v3a-v3b-v3c)*convfac
      elseif (wcharge.eq.-1) then
        s(1,5) = -(-xm12+q1q+s0+tk+v3a)*convfac
        s(1,6) = v3a*convfac
        s(2,5) = -(xm22-q2q+v3b)*convfac
        s(2,6) = v3b*convfac
        s(3,5) = (v1c-v13)*convfac
        s(3,6) = v13*convfac
        s(4,5) = (-xm12-xm22+s0+tk+uk+v13-v1c-v3c)*convfac
        s(4,6) = (v3c-v13)*convfac
        s(5,7) = (xm22+q1q-q2q-uk+v3a+v3b+v3c)*convfac
        s(6,7) = (-xm22-v3a-v3b-v3c)*convfac
      endif
c Symmetrize and set diagonal elements equal to zero
      do i=1,6
        s(i,i)=0.d0
        do j=i+1,7
          s(j,i) = s(i,j)
        enddo
      enddo
c Check consistency of mapping
      if(checks)then
        do i=1,7
          do j=1,7
             schk(i,j) = (
     &              xmom_cm(IQINXMO(i),4)*xmom_cm(IQINXMO(j),4) -
     &              xmom_cm(IQINXMO(i),1)*xmom_cm(IQINXMO(j),1) -
     &              xmom_cm(IQINXMO(i),2)*xmom_cm(IQINXMO(j),2) -
     &              xmom_cm(IQINXMO(i),3)*xmom_cm(IQINXMO(j),3) ) *
     &              2.d0 * convfac
             if (i.le.2) schk(i,j) = -schk(i,j)
             if (j.le.2) schk(i,j) = -schk(i,j)
             if(abs(schk(i,j)).gt.1.d0)then
               if(abs(s(i,j)/schk(i,j)-1.d0).gt.1.d-6)then
                 write(*,*)'fillqinne: inconsistency in mapping'
                 write(*,*)i,j,s(i,j),schk(i,j)
                 stop
               endif
             else
               if(abs(schk(i,j)-s(i,j)).gt.1.d-6)then
                 write(*,*)'fillqinne: inconsistency in mapping'
                 write(*,*)i,j,s(i,j),schk(i,j)
                 stop
               endif
             endif
          enddo
        enddo
      endif

c Fill zs, a, and b arrays for DKS
      do i=1,7
        do j=1,i
c symmetry          
          zs(i,j) = dcmplx(s(i,j),0)
          if (i.ne.j) zs(j,i) = zs(i,j)
        enddo
      enddo
c calculate inner and outer products to be filled in QINNE/A(i,j) 
c and QINNE/B(i,j). Momenta are stored in xmom_cm
c      cc = dcmplx(0.d0,1.d0)
c      rc = dcmplx(1.d0,0.d0)
      zt1 = dcmplx(1.d0,0.d0)
      zt2 = dcmplx(0.d0,1.d0)
      zt3 = dcmplx(-1.d0,0.d0)
      do i = 1,7
        do j = 1,i
          i_xmom = iqinxmo(i)
          j_xmom = iqinxmo(j)
c convention ala DKS to multiply beam particle entries
c          zt = rc
c          if (i.le.2) zt = zt * cc
c          if (j.le.2) zt = zt * cc

          if (j.le.2.and.i.le.2) then 
            zt = dcmplx(-1.d0,0.d0)
          else if (i.le.2) then
            zt = dcmplx(0.d0,1.d0)
          else if (j.le.2) then
            zt = dcmplx(0.d0,1.d0)
          else
            zt = dcmplx(1.d0,0.d0)
          endif

c
c inner products <p1,p2>=xna_2 and [p1,p2]=xnb_2
          call xnaxnb2(
     &           xmom_cm( i_xmom, 4) ,
     &           xmom_cm( i_xmom, 1) ,
     &           xmom_cm( i_xmom, 2) ,
     &           xmom_cm( i_xmom, 3) ,
     &           xmom_cm( j_xmom, 4) ,
     &           xmom_cm( j_xmom, 1) ,
     &           xmom_cm( j_xmom, 2) ,
     &           xmom_cm( j_xmom, 3) ,
     &           xna_2,   
     &           xnb_2)


c Convert to DKS conventions
          a(i,j) =   
     &         xna_2
     &         *dcmplx(sconvfac)
     &         *zt
          b(i,j) = 
     &         xnb_2
     &         *dcmplx(sconvfac)
     &         *zt

c anti-symmetry 
          if (i.ne.j) then
            a(j,i)=-a(i,j)
            b(j,i)=-b(i,j)
          endif

        enddo
      enddo
c
      return
      end


c      function iqinxmo(i)
c Helper function to calculate the correct index when going 
c from QINNE block to XMOM block convention
c
c DKS:
c the POSITION of i1,i2,i3,i4,i5,i6,i7 in the function call
c corresponds to
c ubar(1)^-u(2)^+l(3)^-nubar(4)^+lpbar(5)^+nup(6)^-g(7)^+
c
c MCNLO:
c xmom_cm(i,j) is the j component of the four vector of the particle # i,
c given in the partonic CM frame. j=4 is the energy. i=1,2 are the incoming
c partons, 3 is the outgoing parton, 4 is V1, 5 is V2, 6 and 7 are the leptons
c originating from the decay of V1, 8 and 9 those originating from the 
c decay of V2. See xmadevww() for the fermion/antifermion assignments
c in the case of WW production. Momentum conservation is 
c (1+2)-(3+4+5)=0 or (1+2)-(3+6+7+8+9)=0 
c Subprocess: MC@NLO   -> a(1)b(2) -> c(3)W+(4)[->e+(6)n(7)]W-(5)[->e-(8)n(9)]
c
c The 4-momenta of the particles involved in the hard reaction are stored
c in xmom_cm(ipart,icomp), with icomp=1,2,3,4 corresponding to px,py,pz,E,
c and ipart=1,..,9 with the following identifications:
c ipart=1 -> a; ipart=2 -> b; ipart=3 -> c; ipart=4 -> V1; ipart=5 -> V2;
c ipart=6 -> l; ipart=7 -> lbar; ipart=8 -> r; ipart=9 -> rbar.
c             MadEvent -> a(1)b(2) -> e+(3)n(4)e-(5)n(6)c(7)
c      implicit none
c      integer i, iqinxmo
c
c      if (i.eq.1) iqinxmo = 1
c      if (i.eq.2) iqinxmo = 2
c      if (i.eq.3) iqinxmo = 6
c      if (i.eq.4) iqinxmo = 7
c      if (i.eq.5) iqinxmo = 9
c      if (i.eq.6) iqinxmo = 8
c      if (i.eq.7) iqinxmo = 3
c
c      return
c      end




      subroutine xnaxnb2(a0,a1,a2,a3,b0,b1,b2,b3,xna_2,xnb_2)
c xna_2
c Computes the inner product <p_a,p_b>. The momenta p_a=(a0,a1,a2,a3)
c and p_b=(b0,b1,b2,b3) given in input are rotated in the plane 2-3
c according to the definition adopted in the routine invar_out6
c (taken from gagadiff.for of the gamma*gamma* package)
c
c xnb_2
c Computes the inner product [p_a,p_b]. The momenta p_a=(a0,a1,a2,a3)
c and p_b=(b0,b1,b2,b3) given in input are rotated in the plane 2-3
c according to the definition adopted in the routine invar_out6
c (taken from gagadiff.for of the gamma*gamma* package)
      implicit none
      real*8 aplus,aminus,bplus,bminus,xmoda,xmodb,ca,sa,cb,sb
      real*8 a0,a1,a2,a3,b0,b1,b2,b3,tinyx,a3p,a2p,b3p,b2p
      real*8 mysqrt 
      complex * 16 zic,tmp1,tmp2,tmpa,tmpb,sambp,sbmap
      complex*16 xna_2,xnb_2
      parameter (zic=(0.d0,1.d0))
      parameter (tinyx=1.d-10)
c
      if(abs(a0).lt.tinyx.or.abs(b0).lt.tinyx)then
        tmpa=dcmplx(0.d0)
        tmpb=dcmplx(0.d0)
      else
        call rotate_2(a3,a2,a3p,a2p)
        call rotate_2(b3,b2,b3p,b2p)
        aplus=a0+a3p
        aminus=a0-a3p
        bplus=b0+b3p
        bminus=b0-b3p
        xmoda=sqrt(a1**2+a2p**2)
        xmodb=sqrt(b1**2+b2p**2)
        ca=a1/xmoda
        sa=a2p/xmoda
        cb=b1/xmodb
        sb=b2p/xmodb
c xna_2
        sambp = dcmplx(mysqrt(aminus*bplus))
        sbmap = dcmplx(mysqrt(bminus*aplus))

        tmp1= sambp*(dcmplx(ca)+zic*dcmplx(sa))
        tmp2= sbmap*(dcmplx(cb)+zic*dcmplx(sb))
        tmpa=tmp1-tmp2
c xnb_2
        tmp1= sbmap*(dcmplx(cb)-zic*dcmplx(sb))
        tmp2= sambp*(dcmplx(ca)-zic*dcmplx(sa))
        tmpb=tmp1-tmp2


      endif
      xna_2=tmpa
      xnb_2=tmpb
      return
      end


      subroutine rotate_2(x0,y0,x1,y1)
      implicit none
      real*8 x0,y0,x1,y1
      common/cspsi/cpsi,spsi
      real*8 cpsi,spsi
c
      x1=x0*cpsi+y0*spsi
      y1=y0*cpsi-x0*spsi
      return
      end




      subroutine get_psi1(flfail)
c This routine finds the optimal angle psi, which defines a rotation
c in the 2-3 plane, in order to define the inner products in a 
c numerically safe way. The routine tries to avoid that k+=k0+k3p
c is smaller than 5.d-2 (the four-momenta have here k0=1), where k3p 
c is obtained by rotating the original four-momentum of an angle psi 
c in the 2-3 plane. The first try is with psi=pi/2; if this fails,
c the code tries psi=pi/4, pi/8, pi/16, .... (or pi-pi/4, pi-pi/8, ...,
c depending on the sign of k3 of the particle that has the smallest k+),
c until either succeeds in having k+>5.d-2 for all the particles, or
c maxnum steps are performed. In the latter case, the stability
c condition is relaxed, and the code just tries to avoid that k+=0
c for some momentum. Up to ten steps are performed, using 
c psi=3*pi/8-n*pi/160 (or psi=5*pi/8+n*pi/160), n=0,...,9.
c If for one of these angles all k+#0, then psi is found. Otherwise,
c everything has failed, and the routine returns flfail=.true.
c When this happens, the only way out is to throw the event away,
c in order to avoid divisions by zero when computing inner products.
c The present routine is derived from the one of the gamma*gamma* package.
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      parameter (pio2=1.5707963267948966d0)
      parameter (psi01=1.1780972450961724d0)
      parameter (psi02=1.9634954084936207d0)
      parameter (epalpha=0.019634954084936207d0)
      parameter (unstable=5.d-2)
      parameter (maxnum=10)
      logical flag,flag1,flfail
      real*8 xmom_cm(9,4)
      common/cxmomcm/xmom_cm
      common/cspsi/cpsi,spsi
c
      flfail=.false.
      flag1=.true.
      icount=0
      alpha=pio2
      dalpha=pio2
      yplusmin=3.d0
      dowhile(flag1.and.icount.lt.maxnum)
c First attempt: all k+ larger than unstable=5.d-2
        flag=.true.
        icount=icount+1
        iunst=-100
        ca=cos(alpha)
        sa=sin(alpha)
        do i=1,7
          yplus=xmom_cm(i,4)+ca*xmom_cm(i,3)+sa*xmom_cm(i,2)
          if(yplus.lt.unstable.and.yplus.lt.yplusmin)then
            flag=.false.
            iunst=i
            yplusmin=yplus
          endif
        enddo
        if(.not.flag)then
c Decrease alpha and try again
          dalpha=dalpha/2.d0
          if(xmom_cm(iunst,3).gt.0.d0)then
            alpha=dalpha
          else
            alpha=pi-dalpha
          endif
        else
c Psi found: skip the following search and exit
          flag1=.false.
        endif
      enddo
      if(flag1)then
c Second attempt: all k+ different from zero
        jcount=0
        flag=.true.
        if(xmom_cm(iunst,3).gt.0.d0)then
          xsgn=-1.d0
          psi0=psi01
        else
          xsgn=1.d0
          psi0=psi02
        endif
        dowhile(flag.and.jcount.lt.10)
          alpha=psi0+xsgn*epalpha*jcount
          jcount=jcount+1
          ca=cos(alpha)
          sa=sin(alpha)
          icheck=1
          do i=1,7
            yplus=xmom_cm(i,4)+ca*xmom_cm(i,3)+sa*xmom_cm(i,2)
            if(yplus.eq.0.d0)icheck=0
          enddo
          if(icheck.eq.1)flag=.false.
        enddo
        if(flag)then
c Failed: at least one k+ equal to zero. Set flfail to true and exit
          flfail=.true.
        else
          cpsi=ca
          spsi=sa
        endif
      else
        cpsi=ca
        spsi=sa
      endif
      return
      end
c
c
c End of fillqinne.f
c
c
c
c
c Begin of packanomDKS.f
c
c
*****************   functions_fm_A.f  *************
*  this package uses common blocks /initv/ /quinne/ /pmom/ 
*  
*
* this  package contains the functions
*       FUNCTION fmwwall(npart,iso,ii)    used common blocks:  /initv,//pmom,//qinnne/
*       FUNCTION k1pk2sq(k1,k2,isign)
*       FUNCTION AmplA(npart,ii,N1,N2,N3,N4,N5,N6,n7) used common blocks: /qinne/
*       FUNCTION AmplB(npart,ii,N1,N2,N3,N4,N5,N6,n7) used common blocks: /qinne/
*       FUNCTION AmplS(npart,ii,N1,N2,N3,N4,N5,N6,n7) used common blocks: /qinne/
*       FUNCTION ann4(N1,N2,N3,N4)                    used common blocks: /qinne/      
*       FUNCTION AloopA(i1,i2,i3,i4,i5,i6)  used common blocks: /qinne/,/nparton/
*       FUNCTION AloopB(i1,i2,i3,i4,i5,i6)  used common blocks: /qinne/,/nparton/
*       FUNCTION AtreeA6(N1,N2,N3,N4,N5,N6) used common blocks: /qinne/ 
*       FUNCTION Atreeb6(N1,N2,N3,N4,N5,N6) used common blocks: /qinne/  
*       FUNCTION RAND2(randy)
*       FUNCTION TuWW(i1,i2,i3,i4,i5,i6)    used common blocks: /qinne/  
*       FUNCTION LOGuS56uWW(i1,i2,i3,i4,i5,i6) used common blocks: /qinne/    
*       FUNCTION LOGuS34uWW(i1,i2,i3,i4,i5,i6) used common blocks: /qinne/   
*       FUNCTION LOGuT34uWW(i1,i2,i3,i4,i5,i6) used common blocks: /qinne/
*       FUNCTION LOGuT56uWW(i1,i2,i3,i4,i5,i6) used common blocks: /qinne/
*       FUNCTION POLYuAuWW(i1,i2,i3,i4,i5,i6)  used common blocks: /qinne/
*       FUNCTION LTILDEuWW(i1,i2,i3,i4,i5,i6)  used common blocks: /qinne/
*       FUNCTION F6uAuWW(i1,i2,i3,i4,i5,i6)    used common blocks: /qinne/
*       FUNCTION SF(k1,k2)                                                  
*       FUNCTION SDOTF(k1)                                                  
*       FUNCTION SDOT3(i1,i2,i3)               used common blocks: /qinne/                           
*       FUNCTION SDOT4(i1,i2,i3,i4)            used common blocks: /qinne/                              
*       FUNCTION GRAM3(s1,s2,s3)                                            
*       FUNCTION COSTH(k1,k2)                  used common blocks: /qinne/                             
*       FUNCTION AM3(i1,i2,i3)                 used common blocks: /qinne/
*       FUNCTION AM4(i1,i2,i3,i4)              used common blocks: /qinne/
*       FUNCTION AM5(i1,i2,i3,i4,i5)           used common blocks: /qinne/
*       FUNCTION BM3(i1,i2,i3)                 used common blocks: /qinne/
*       FUNCTION BM4(i1,i2,i3,i4)              used common blocks: /qinne/
*       FUNCTION AM6(i1,i2,i3,i4,i5,i6)        used common blocks: /qinne/   
*       FUNCTION BM6(i1,i2,i3,i4,i5,i6)        used common blocks: /qinne/
*       FUNCTION LI2(x)                  
*       FUNCTION LN(arg1, arg2)              used common blocks: /nparton/        
*       FUNCTION DILOG(arg1,arg2)                                           
*       FUNCTION L0(arg1,arg2)                                              
*       FUNCTION L1(arg1,arg2)                                              
*       FUNCTION LSM12MHNEW(ss,tt,m1sq,m2sq)  used common blocks: /nparton/                             
*       FUNCTION I33M(s1,s2,s3)                                            
*       FUNCTION NA(k1,k2)                        
*       FUNCTION NNB(k1,k2)
*       FUNCTION Tiny(x)
*       FUNCTION NB(k1,k2)
*       FUNCTION NNA(k1,k2)
*       SUBROUTINE CRASH(functionname)
*       SUBROUTINE fillquin(pcm,a,b,s,zs)
*************************************************************************************************

*
      FUNCTION fmwwall(npart,iso,ii)
            IMPLICIT NONE
            INTEGER iso,ii,npart
            REAL*8 fmwwall
*     amplitude squared for qqbar'->w^+w^-g and pos. hel. gluon 
*     m(i)_tree_n=|ampla
*     (n,ii,1,2,3,4,5,6,7)+cl(i)amplb(n,ii,1,2,3,4,5,6,7)|^2+
*               |cr(i)*amplb(n,ii2,1,3,4,5,6,7)|^2
*     where n=6,7 in case of n=6 the seventh vector is not used
* ..............  common blocks ..........................
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*                                                                 4)pmom
      real*8 plaball,pcmall,pwwg
      common/pmom/ plaball(4,7,4),pcmall(4,7,4),
     & pwwg(5,3,4)             
*                                                    5) quinne               
      complex*16 a,b,zs
      real*8 s
      common/qinne/a(7,7),b(7,7),zs(7,7),s(7,7)
*.............. end of common block ....................
*
      real*8 s12,sww,fmtr1,fmtr2,
     & constww,k1pk2sq
      real*8 alfa
      integer i,j
      parameter( alfa=1.d0/128.d0)
      complex *16 fcl(2),fcr(2)
      complex*16 a1,a2,a3,z1u,z2u,ampla,amplb
      complex*16 b1,b2,b3,zc1u,zc2u
* -------------------------------------------------------------
      s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
      if(ii.eq.1.and.npart.eq.7)then
      sww =k1pk2sq(pwwg(1,1,ii),pwwg(1,2,ii),1)
      else
      sww=s12
      endif
*
      call qinn(npart,pcmall(1,1,ii),a,b,zs,s)
      if(s(3,4).eq.0.d0.and.ii.eq.3)then
      write(1,1)ii,s(3,4),((pcmall(i,j,3),i=1,4),j=1,6)
 1    format('ii,s(3,4):',i4,g12.5,/,6(4g17.7,/))
      endif
* couplings introduced by lance (eq. 17,19)
* ----------------------------------------
      fcl(1)=dcmplx(2*qu*swsq+sww*(1-2*qu*swsq)/(sww-mz**2))
      fcl(2)=dcmplx(-2*qd*swsq+sww*(1+2*qd*swsq)/(sww-mz**2))
      fcr(1)=2*dcmplx(-qu*swsq*mz**2/(sww-mz**2))
      fcr(2)=2*dcmplx(qd*swsq*mz**2/(sww-mz**2))
*
* calculate the contribution to the ampliltude squared
* -----------------
      fmtr1=0
      fmtr2=0
      fmwwall=0.d0
* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(iso.eq.1)then
* for positive helicity gluons
      a1=ampla(npart,ii,1,2,3,4,5,6,7)
      a2=amplb(npart,ii,1,2,3,4,5,6,7)
      a3=amplb(npart,ii,2,1,3,4,5,6,7)
      b1=a1
      b2=a2
      b3=a3
        if(npart.eq.6.and.ii.eq.1)then
* loop amplitudes
        b1=ampla(npart,ii+1,1,2,3,4,5,6,7)
        b2=amplb(npart,ii+1,1,2,3,4,5,6,7)
        b3=amplb(npart,ii+1,2,1,3,4,5,6,7)
        endif
      zc1u=(b1)+fcl(1)*(b2)
      zc2u=fcr(1)*(b3)
      z1u=a1+fcl(1)*a2
      z2u=fcr(1)*a3
       fmtr1=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
       fmtr2=0.d0
* for negative helicity gluons 'flip1' if npart=7
      if(npart.eq.7)then
        a1=ampla(npart,ii,2,1,5,6,3,4,7)
        a2=amplb(npart,ii,2,1,5,6,3,4,7)
        a3=amplb(npart,ii,1,2,5,6,3,4,7)
        b1=a1
        b2=a2
        b3=a3
        zc1u=(b1)+fcl(1)*(b2)
        zc2u=fcr(1)*(b3)
        z1u=a1+fcl(1)*a2
        z2u=fcr(1)*a3
        fmtr2=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
        endif
      fmwwall=fmtr2+fmtr1
      endif
* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(iso.eq.2)then
* positive helicity gluons
      a1=ampla(npart,ii,1,2,6,5,4,3,7)
      a2=amplb(npart,ii,1,2,6,5,4,3,7)
      a3=amplb(npart,ii,2,1,6,5,4,3,7)
      b1=a1
      b2=a2
      b3=a3
        if(npart.eq.6.and.ii.eq.1)then
* loop contributions
        b1=ampla(npart,ii+1,1,2,6,5,4,3,7)
        b2=amplb(npart,ii+1,1,2,6,5,4,3,7)
        b3=amplb(npart,ii+1,2,1,6,5,4,3,7)
        endif
      zc1u=(b1)+fcl(2)*(b2)
      zc2u=fcr(2)*(b3)
      z1u=a1+fcl(2)*a2
      z2u=fcr(2)*a3
      fmtr1=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
      fmtr2=0.d0
* for negative helicity gluons 'flip1' if npart=7
      if(npart.eq.7)then
        a1=ampla(npart,ii,2,1,4,3,6,5,7)
        a2=amplb(npart,ii,2,1,4,3,6,5,7)
        a3=amplb(npart,ii,1,2,4,3,6,5,7)
        b1=a1
        b2=a2
        b3=a3
        zc1u=b1+fcl(2)*(b2)
        zc2u=fcr(2)*(b3)
        z1u=a1+fcl(2)*a2
        z2u=fcr(2)*a3
        fmtr2=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
      endif
      fmwwall=fmtr2+fmtr1
      endif
      if(iso.ne.1.and.iso.ne.2)then
       write(6,*)'invalid iso in fmtree=',iso
       stop
       endif
* final over all constant
* -----------------------
       constww=(6*alfa/swsq)**2*mw**4/(32*s12*nc)
       if(npart.eq.7)constww=2*constww
       if(npart.eq.6.and.ii.eq.1)constww=2*cf*constww
       fmwwall=constww*fmwwall
       return
       end
*
************** (11) (K1+/-K2)**2 ********
*
       FUNCTION k1pk2sq(k1,k2,isign)
       implicit none
       real*8  k1pk2sq,k1(4),k2(4),tes
       integer i,isign
*
        tes=0
        do i=1,3
        tes = tes + (k1(i)+isign*k2(i))**2
        enddo
        k1pk2sq=(k1(4)+isign*k2(4))**2-tes 
        return
        end 

*
******  (15) AmplA as defined in eq. (4) of Lance note *****
      FUNCTION AmplA(npart,ii,N1,N2,N3,N4,N5,N6,n7)
* npart=6,ii=1: loop
* npart=6,ii=2: born
* npart=7,ii=1,2,3,4:gluon bremsstrahlung: full,soft,collp,collm  
      implicit none 
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 AmplA,aloopa,atreea6
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag,aa1,aa2,am3,t134,t256
* //////////////////////////////////
*
      imag = (0.0 ,1.d0 )
      if(ii.eq.1.and.npart.eq.6)then
         AmplA=aloopa(n1,n2,n3,n4,n5,n6)
      endif
      if(ii.eq.1.and.npart.eq.7) then
         t134=zs(n1,n3)+zs(n1,n4)+zs(n3,n4)
         t256=zs(n2,n5)+zs(n2,n6)+zs(n5,n6)
         aa1=a(n1,n3)*b(n3,n4)*b(n2,n5)*(am3(n6,n2,n7)+am3(n6,n5,n7) )
     &        /t256
         aa2=(am3(n6,n1,n4)+am3(n6,n3,n4))*(am3(n1,n2,n5)+am3(n1,n7,n5))
     &        /a(n7,n2)
         AmplA=imag*a(n1,n3)*(aa1+aa2)
     &        /(a(n1,n7)*zs(n3,n4)*zs(n5,n6)*t134)
      endif 
      if(ii.ne.1)then
         ampla=atreea6(n1,n2,n3,n4,n5,n6)
      endif
      return
      end 

*
*
******  (16) AmplB as defined in eq. (4) of Lance note *****
      FUNCTION AmplB(npart,ii,N1,N2,N3,N4,N5,N6,n7)
      implicit none 
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 AmplB,aloopb,atreeb6
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag,aa1,aa2,aa3,t127,am3,ann4
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      if(ii.eq.1.and.npart.eq.6)then
       AmplB=aloopb(n1,n2,n3,n4,n5,n6)
      endif
      if(ii.eq.1.and.npart.eq.7)then
      t127=zs(n1,n2)+zs(n1,n7)+zs(n2,n7)
      aa1=a(n3,n6)*b(n4,n5)*(ann4(n1,n5,n2,n1)+ann4(n1,n5,n7,n1)+
     & ann4(n1,n6,n2,n1)+ann4(n1,n6,n7,n1))
      aa2=a(n1,n3)*(am3(n1,n2,n4)+am3(n1,n7,n4))*
     & (am3(n6,n3,n5)+am3(n6,n4,n5))
      aa3=a(n1,n6)*(am3(n1,n2,n5)+am3(n1,n7,n5))*
     & (am3(n3,n5,n4)+am3(n3,n6,n4))
      AmplB=imag*(-aa1+aa2-aa3)/
     & (t127*a(n1,n7)*a(n7,n2)*zs(n3,n4)*zs(n5,n6))
      endif
      if(ii.ne.1)then
       amplb=atreeb6(n1,n2,n3,n4,n5,n6)
      endif
      return
      end 
*

******  (15') AmplS as defined as
      FUNCTION AmplS(npart,ii,N1,N2,N3,N4,N5,N6,n7)
* npart=6,ii=1: loop
* npart=6,ii=2: born
* npart=7,ii=1,2,3,4:gluon bremsstrahlung: full,soft,collp,collm  
      implicit none 
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 AmplS,AmplA
*...
       AmplS= AmplA(npart,ii,N1,N2,N3,N4,N5,N6,n7)
     &  + AmplA(npart,ii,N1,N2,N6,N5,N4,N3,n7)
       return
       end
******  Ann4 auxiliary functions for products of inner products
*       in AmplB
      FUNCTION ann4(N1,N2,N3,N4)
* <n1|n2 n3|n4>
      implicit none 
      integer n1,n2,n3,n4
      complex*16 ann4
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      ann4=a(n1,n2)*b(n2,n3)*a(n3,n4)
      return
      end
      FUNCTION bnn4(N1,N2,N3,N4)
* <n1|n2 n3|n4>
      implicit none 
      integer n1,n2,n3,n4
      complex*16 bnn4
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      bnn4=b(n1,n2)*a(n2,n3)*b(n3,n4)
      return
      end
*
       FUNCTION AloopA(i1,i2,i3,i4,i5,i6)
       implicit none
       integer i1,i2,i3,i4,i5,i6                                   
       complex*16 Aloopa,atreea6, Vua,F6uauww                                         
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*                                                                 1) nparton
       complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart
*
       real*8 pi,s12,musq 
       parameter(    pi = 3.14159265358979323846d0 ) 
       s12 = S(i1,i2)    
       musq=muscloop**2                                                         
                                                                                
*            !! Add DR->CDR shifts later in formulas below                         
                                                                                  
       Vua = - 0.5*(log(musq/s12)**2 - pi**2) - 1.5d0*log(musq/s12)      
     &        - 4.d0 - imag*abspart*pi*(log(musq/s12) + 1.5d0)            


       Aloopa =atreea6(i1,i2,i3,i4,i5,i6)*Vua +
     & F6uauww(i1,i2,i3,i4,i5,i6)

       end                                                                           
*                                                                                  
       FUNCTION AloopB(i1,i2,i3,i4,i5,i6)
       implicit none
       integer  i1,i2,i3,i4,i5,i6                                   
       complex*16 Aloopb, Vub                                         
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      real*8 s12,musq,pi
      parameter(          pi = 3.14159265358979323846d0 )                                                                           
*                                                                 1) nparton
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart
*
      complex*16  atreeb6
       musq=muscloop**2
       s12 = S(i1,i2)                                                             
                                                                                  
*            !! Add DR->CDR shifts later in formulas below                         
                                                                                  
       Vub = - 0.5*(log(musq/s12)**2 - pi**2) - 1.5*log(musq/s12)      
     &        - 3.5d0 - imag*abspart*pi*(log(musq/s12) + 1.5)            
       Aloopb = atreeb6(i1,i2,i3,i4,i5,i6)*Vub
        end

******   AtreeA6 as defined in eq. (4) of Lance note *****
      FUNCTION AtreeA6(N1,N2,N3,N4,N5,N6)
      implicit none 
      integer n1,n2,n3,n4,n5,n6
      complex*16 AtreeA6
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag,aa1,aa2,t134
* //////////////////////////////////
*
      imag = (0.0 ,1.d0 )
       aa1 = A(N1,N3)*B(N2,N5)*(- A(N6,N3)*B(N3,N4)
     X   - A(N6,N1)*B(N1,N4) )
      t134=zs(n1,n3)+zs(n1,n4)+zs(n3,n4)
       aa2=imag/(t134*zs(n3,n4)*zs(n5,n6)) 
      atreea6=aa2*aa1/imag
      return
      end

******  (16) Atreeb6 as defined in eq. (4) of Lance note *****
      FUNCTION Atreeb6(N1,N2,N3,N4,N5,N6)
      implicit none 
      integer n1,n2,n3,n4,n5,n6
      complex*16 Atreeb6
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag,aa1,aa2,aa3
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      Aa1 =imag/(zs(n1,n2)*zs(n3,n4)*zs(n5,n6))
      aa2=A(n1,n3)*B(n2,n5)*(A(n6,n2)*B(n2,n4)+A(n6,n5)*b(n5,n4))
      aa3=A(n1,n6)*B(n2,n4)*(A(n3,n1)*B(n1,n5)+A(n3,n6)*b(n6,n5))
      atreeb6=aa1*(aa2+aa3)/imag
      return
      end 


*
*
************ (17) random number with seed ************
*
      FUNCTION RAND2(randy)
*    This is the usual "random" 
      implicit none
      DOUBLE PRECISION MINV,RAND2
      integer m,a,Qran,r,hi,lo,randy
      PARAMETER(M=2147483647,A=16807,Qran=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = RANDY/Qran
      LO = MOD(RANDY,Qran)
      RANDY = A*LO - R*HI
      IF(RANDY.LE.0) RANDY = RANDY + M
      RAND2 = RANDY*MINV
      END 
*
* ----------------------------------
*
       FUNCTION TuWW(i1,i2,i3,i4,i5,i6)
       implicit none
       integer i1,i2,i3,i4,i5,i6                                   
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       real*8 sdot3,gram3,del,del3u12, del3u34, del3u56                         
       complex*16  tuww,am3,am4                                                
       del3u12 = S(i2,i1) - S(i5,i6) - S(i3,i4)                                  
       del3u34 = - S(i2,i1) - S(i5,i6) + S(i3,i4)                                
       del3u56 = - S(i2,i1) + S(i5,i6) - S(i3,i4)               
       del = gram3(S(i5,i6), S(i2,i1), S(i3,i4))                           
       tuww =                         
     &   + 1.5*S(i2,i1)*del3u12*(Sdot3(i5,i6,i2)-Sdot3(i5,i6,i1))   
     &    *AM4(i6,i2,i1,i5)*AM4(i3,i2,i1,i4) /AM4(i2,i5,i6,i1)/del**2            
     &   - 0.5*(3*S(i2,i1)+2*Sdot3(i5,i6,i2))*AM4(i6,i2,i1,i5)                  
     &           *AM4(i3,i2,i1,i4)/AM4(i2,i5,i6,i1)/del       
     & + A(i6,i3)*B(i5,i4)*S(i2,i1)*Sdot3(i5,i6,i2)/AM4(i2,i5,i6,i1)/del         
     &   - B(i5,i6)*A(i3,i4)*AM4(i6,i2,i1,i4)**2/AM4(i2,i5,i6,i1)/del   
     &   + Sdot3(i5,i6,i2)/AM4(i2,i5,i6,i1)**2/del*(                             
     &     + B(i5,i1)*A(i2,i3)*                                
     &             ( AM3(i6,i5,i4)*del3u56 - AM3(i6,i3,i4)*del3u34 )  
     &     - A(i6,i2)*B(i1,i4)*                            
     &             ( AM3(i3,i6,i5)*del3u56 - AM3(i3,i4,i5)*del3u34 ) )   
     &   + A(i6,i2)*B(i1,i4)*AM4(i3,i6,i2,i5)/AM4(i2,i5,i6,i1)**2     
     &   - 0.5*A(i6,i2)**2*B(i1,i4)**2*                                          
     &            (Sdot3(i5,i6,i2)*del3u12+2*S(i5,i6)*S(i3,i4))                  
     &      /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)**3                               
       tuww =   tuww    - 2*A(i6,i1)*B(i2,i4)                                    
     &      *( AM3(i6,i5,i4)*del3u56 - AM3(i6,i3,i4)*del3u34 )                   
     &      /A(i5,i6)/B(i3,i4)/del                                               
     &   + 2*AM4(i6,i5,i2,i4)/AM4(i2,i5,i6,i1)/del*(           
     &      + ( AM3(i6,i5,i2)*AM3(i2,i1,i4)*del3u56                       
     &        - AM3(i6,i2,i1)*AM3(i1,i3,i4)*del3u34                
     &        + AM4(i6,i5,i2,i4)*S(i2,i1)*del3u12 ) /A(i5,i6)/B(i3,i4)      
     &      + 2*AM4(i3,i6,i2,i5)*S(i2,i1) )                                      
     &   - 2*AM4(i6,i5,i2,i4)*B(i5,i1)*A(i2,i3)/AM4(i2,i5,i6,i1)**2     
     &   + A(i6,i2)*B(i1,i4)*AM4(i6,i5,i2,i4)*del3u12                      
     &      /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)**2                              
     &   + 0.5/AM4(i2,i5,i6,i1)*(                                                
     &     + 3*AM3(i6,i2,i4)*AM3(i6,i1,i4)/A(i5,i6)/B(i3,i4)                     
     &     + A(i6,i1)*B(i5,i4)*B(i1,i4)/B(i3,i4)                          
     &     - A(i6,i2)*A(i6,i3)*B(i2,i4)/A(i5,i6)                                 
     &     + AM3(i3,i2,i5)*AM3(i3,i1,i5)/B(i5,i6)/A(i3,i4)                
     &     + B(i5,i2)*A(i6,i3)*A(i2,i3)/A(i3,i4)                        
     &     - B(i5,i1)*B(i5,i4)*A(i1,i3)/B(i5,i6)                 
     &     + 4*B(i5,i4)*A(i6,i3) )                      
     &   + 0.5/AM4(i1,i5,i6,i2)*(                                               
     &     + A(i6,i1)**2*B(i2,i4)**2/A(i5,i6)/B(i3,i4)           
     &     - B(i5,i2)**2*A(i1,i3)**2/B(i5,i6)/A(i3,i4) )                       
        end                                                  
*
*                                                         
       FUNCTION LOGuS56uWW(i1,i2,i3,i4,i5,i6)
       implicit none                          
       integer i1,i2,i3,i4,i5,i6                                   
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       real*8 gram3,del,del3u12, del3u34,sdot3
       complex*16 logus56uww,am3,am4,bm3,am6
*          !! coefficient of log(s/s56) = log(s/Mw) 
*          !! coefficient of log(s/s34) can be obtained by "flip" 
       del3u12 = S(i2,i1) - S(i5,i6) - S(i3,i4)      
       del3u34 = - S(i2,i1) - S(i5,i6) + S(i3,i4)    
       del = gram3(S(i5,i6), S(i2,i1), S(i3,i4))     
       logus56uww =         
     &     1.5 * del3u34 * AM4(i6,i2,i1,i5)*AM4(i3,i2,i1,i4)         
     &   * (Sdot3(i5,i6,i2)-Sdot3(i5,i6,i1))/AM4(i2,i5,i6,i1)/del**2   
     &   + 1.5 * B(i5,i4)*AM6(i6,i2,i1,i5,i6,i3)/AM4(i2,i5,i6,i1)/del 
     &   + 0.5 * AM3(i3,i6,i5)*AM6(i6,i3,i4,i2,i1,i3)/A(i3,i4)   
     &                   /AM4(i2,i5,i6,i1)/del         
     &   + 0.5 * Sdot3(i5,i6,i2)/AM4(i2,i5,i6,i1)/del          
     &     * (-2 * A(i6,i3)*B(i5,i4) + B(i5,i6)*A(i6,i3)**2/A(i3,i4)  
     &          + A(i5,i6)*B(i5,i4)**2/B(i3,i4) )                
     &   + ( AM4(i3,i6,i2,i5)/A(i3,i4)                                         
     &     - B(i5,i6)*A(i6,i2)*B(i1,i4)/AM4(i2,i5,i6,i1) )  
     &      * (A(i6,i3)*del3u12 - 2*AM3(i6,i5,i4)*A(i4,i3)) 
     &                   /AM4(i2,i5,i6,i1)/del               
     &+A(i6,i2)*B(i1,i4)*Sdot3(i5,i6,i2)/B(i3,i4)
     &  /AM4(i2,i5,i6,i1)**2/del       
     &      * (B(i5,i4)*del3u12 - 2*BM3(i5,i6,i3)*B(i3,i4))  
     &   + 4/AM4(i2,i5,i6,i1)/del                                              
     &*(AM3(i3,i6,i5)*AM4(i6,i5,i2,i4) + AM3(i6,i5,i4)*AM4(i3,i6,i1,i5))  
     &   + 2 * del3u12/AM4(i2,i5,i6,i1)/del   
     &       * ( A(i6,i3) * AM4(i3,i6,i1,i5)/A(i3,i4)                      
     &         - B(i5,i4) * AM4(i6,i5,i2,i4)/B(i3,i4) )  
*        !!! Now add the terms (with the flip) coming from the log*log terms 
       logus56uww =  logus56uww +                                
     &     conjg( 0.5 * A(i4,i1)*B(i2,i6) * AM4(i4,i3,i1,i6)        
     &               /A(i3,i4)/B(i5,i6)/AM4(i1,i3,i4,i2)**2              
     &          - 0.75 * AM4(i4,i3,i1,i6)**2                                   
     &      /A(i3,i4)/B(i5,i6)/Sdot3(i3,i4,i1)/AM4(i1,i3,i4,i2) )     
*         !! change of argument of log from s56/s12 -> s12/s56   
       logus56uww = - logus56uww
        end


       FUNCTION LOGuS34uWW(i1,i2,i3,i4,i5,i6)    
       implicit none           
       integer i1,i2,i3,i4,i5,i6                                   
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      complex*16 logus34uww,logus56uww, Abak(7,7), Bbak(7,7)  
       integer i,j
       do i=1,7
       do j=1,7                                                                           
       Abak(i,j) = A(i,j) 
       Bbak(i,j) = B(i,j)                                                         
       A(i,j) = Bbak(i,j) 
       B(i,j) = Abak(i,j)                                                         
       enddo
       enddo                                                  
       logus34uww = - logus56uww(i2,i1,i5,i6,i3,i4)                               
       do i=1,7
       do j=1,7                                                                           
       A(i,j) = Abak(i,j)
       B(i,j) = Bbak(i,j)                                                         
       enddo
       enddo
        end                                                                          
                                                                                  
       FUNCTION LOGuT34uWW(i1,i2,i3,i4,i5,i6)
       implicit none
       integer i1,i2,i3,i4,i5,i6                                   
       complex*16 logut34uww                                          
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      complex*16 am4
      real*8 sdot3
                                                                                  
*         !! coefficient of log(-t/s34)  t == s134                                 
*         !! coefficient of log(-t/s56)  can be obtained with "flip"               
                                                                                  
       logut34uww = - 0.5 * B(i4,i1)**2*A(i1,i6)**2                              
     &          /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)                                
     &    /(1.d0 - S(i3,i4)/Sdot3(i5,i6,i2))**2 / Sdot3(i5,i6,i2)               
     &    + 2 * B(i4,i1)*A(i1,i6) * AM4(i6,i5,i2,i4)                               
     &     /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)                                     
     &           /(1. - Sdot3(i5,i6,i2)/S(i3,i4)) / S(i3,i4)                       
     &   - A(i6,i2)*A(i6,i1)*B(i4,i1)**2 * Sdot3(i5,i6,i2)                         
     &         /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)**2                              
     &           /(1. - Sdot3(i5,i6,i2)/S(i3,i4)) / S(i3,i4)                       
     &   - 0.5 * A(i6,i2)*B(i4,i1) * AM4(i6,i5,i2,i4)                              
     &      /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)**2                                 
     &   - 0.75 * AM4(i6,i5,i2,i4)**2                                              
     &        /A(i5,i6)/B(i3,i4)/Sdot3(i5,i6,i2)/AM4(i2,i5,i6,i1)                    
        end                                                                          
*                                                                                  
       FUNCTION LOGuT56uWW(i1,i2,i3,i4,i5,i6)
       implicit none
       integer i1,i2,i3,i4,i5,i6 
       complex*16 logut56uww, logut34uww,Abak(7,7), Bbak(7,7)                    
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)

       integer i,j
       do i=1,7
       do j=1,7                                                                           
       Abak(i,j) = A(i,j)
       Bbak(i,j) = B(i,j)                                                         
       A(i,j) = Bbak(i,j)
       B(i,j) = Abak(i,j)                                                         
       enddo
       enddo
                                                                                  
       logut56uww = - logut34uww(i2,i1,i5,i6,i3,i4)                               
       do i=1,7
       do j=1,7                                                                           
                                                                                  
       A(i,j) = Abak(i,j)
       B(i,j) = Bbak(i,j)                                                         
       enddo
       enddo
       end
*
       FUNCTION POLYuAuWW(i1,i2,i3,i4,i5,i6)
       implicit none
       integer i1,i2,i3,i4,i5,i6                                   
       real*8 sdot3,del,gram3                                              
       complex*16 polyuauww                                           
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      complex*16 am4,bm4
                                                                                  
       del = gram3(S(i5,i6), S(i2,i1), S(i3,i4))                                  
                                                                                  
       polyuauww = 0.5 /AM4(i2,i5,i6,i1)/del *                                    
     &     (Sdot3(i5,i6,i2)*(S(i5,i6)+S(i3,i4)-S(i2,i1)) -                         
     &           2*S(i5,i6)*S(i3,i4)) *                                           
     & (B(i4,i5)**2/B(i5,i6)/B(i3,i4) + A(i3,i6)**2/A(i5,i6)/A(i3,i4))          
     &    + (S(i5,i2)+S(i6,i2)-S(i5,i1)-S(i6,i1))*B(i4,i5)*A(i3,i6)               
     &         /AM4(i2,i5,i6,i1)/del                                                
                                                                                  
*       !! Add terms coming from L1 and its flip                                   
                                                                                  
       polyuauww =  polyuauww  +   0.5 * B(i4,i1)**2*A(i1,i6)**2                 
     & /A(i5,i6)/B(i3,i4)/AM4(i2,i5,i6,i1)/(Sdot3(i5,i6,i2)
     & - S(i3,i4))   -   0.5 * A(i6,i2)**2*B(i2,i4)**2                 
     &  /B(i3,i4)/A(i5,i6)/BM4(i1,i3,i4,i2)/
     & (Sdot3(i3,i4,i1) - S(i5,i6))          
                                                                                  
*       !! Add additional antisymmetric term                                       
                                                                                  
       polyuauww =  polyuauww  - 0.5*( B(i4,i5)**2/B(i5,i6)/B(i3,i4) 
     &           + A(i3,i6)**2/A(i5,i6)/A(i3,i4))/AM4(i2,i5,i6,i1)    
     &     - 0.5* AM4(i6,i5,i2,i4)*AM4(i3,i6,i1,i5)/                   
     &               S(i5,i6)/S(i4,i3)/AM4(i2,i5,i6,i1)                             
                                                                                  
*     !!$  !! Add DR -> CDR shift terms (0.5) and constant terms from VuauWW (4)   
*     !!$                                                                          
*     !!$  polyuauww =  polyuauww - 4.5*a0u6ua(i1,i2,i3,i4,i5,i6)                  
        end                                      
*
       FUNCTION LTILDEuWW(i1,i2,i3,i4,i5,i6)
       implicit none           
       integer i1,i2,i3,i4,i5,i6                                   
       complex*16 ltildeuww,am4                                           
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      real*8 sdot3 
*         !! coefficient of Lstildeu{-1}^{2mh}                                     
       ltildeuww = B(i5,i2)**2*A(i1,i3)**2                                       
     &        /B(i5,i6)/A(i3,i4)/Sdot3(i5,i6,i2)/AM4(i1,i5,i6,i2)      
     &    - AM4(i2,i5,i6,i4)**2 * AM4(i6,i5,i2,i1)**2                      
     &        /A(i5,i6)/B(i3,i4)/Sdot3(i5,i6,i2)/AM4(i2,i5,i6,i1)**3             
        end                                                                      
*     end funCTION LTILDEuWW                                                     

       FUNCTION F6uAuWW(i1,i2,i3,i4,i5,i6) 
       implicit none                                            
       integer i1,i2,i3,i4,i5,i6                                   
       complex*16 f6uauww                                             
       real*8 sdot3,s12, s34, s56, t134,I33m                                  
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      complex*16 ln,Ltildeuww,tuww,logus56uww,logus34uww,
     &logut34uww,logut56uww,Lsm12mhnew,polyuauww
       s12 = S(i1,i2)
       s34 = S(i3,i4) 
       s56 = S(i5,i6)  
       t134 = Sdot3(i1,i3,i4)
* 
       f6uauww =  
     &Ltildeuww(i1,i2,i3,i4,i5,i6)*Lsm12mhnew(s12,t134,s34,s56)
     &   + tuww(i1,i2,i3,i4,i5,i6)*I33m(s12,s34,s56)
     &   + logus56uww(i1,i2,i3,i4,i5,i6)*ln(s12,s56) 
     &   + logus34uww(i1,i2,i3,i4,i5,i6)*ln(s12,s34)      
     &   + logut34uww(i1,i2,i3,i4,i5,i6)*ln(t134,s34)   
     &   + logut56uww(i1,i2,i3,i4,i5,i6)*ln(t134,s56)                       
     &   + polyuauww(i1,i2,i3,i4,i5,i6)         
        end
*
       FUNCTION SF(k1,k2)                                                  
       implicit none                                                       
       real*8 k1(4),k2(4)                         
       real*8 SF,dotdot                                      
       dotdot = k1(4)*k2(4) - k1(1)*k2(1) - k1(2)*k2(2) - k1(3)*k2(3)     
       SF =  2*dotdot                                                     
       end                                                                    
                                                                           
       FUNCTION SDOTF(k1)                                                  
       implicit none                                                       
       real*8 k1(4)                               
       real*8 SdotF,SF                                              
       SdotF = 0.5d0*SF(k1,k1)                                          
       end                                                                    

                                                                           
                                                                           
       FUNCTION SDOT3(i1,i2,i3)                                          
       implicit none                                                       
       integer i1,i2,i3                                     
       real*8  Sdot3                                            
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
         Sdot3 = S(i1,i2) + S(i2,i3) + S(i1,i3)                             
       end  
*
*
       FUNCTION SDOT4(i1,i2,i3,i4)                                          
       implicit none                                                       
       integer i1,i2,i3, i4                                             
       real*8  Sdot4                                            
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
         Sdot4 = S(i1,i2) + S(i1,i3) + S(i1,i4) + S(i2,i3)
     &        + S(i2,i4) + S(i3,i4)                                     
       end  

                                                                           
                                                                           
       FUNCTION GRAM3(s1,s2,s3)                                            
       implicit none                                                       
       real*8 s1,s2,s3                            
       real*8  gram3                                           
       gram3 = s1**2 + s2**2 + s3**2 - 2*(s1*s2 + s2*s3+ s3*s1)            
       end                                                                    


                                                                           
                                                                           
                                                                           
       FUNCTION COSTH(k1,k2)                                               
       implicit none                                       
       real*8 k1(4),k2(4)                         
       real*8  costh, abs1, abs2                               
       abs1 = dsqrt(k1(1)**2 + k1(2)**2 + k1(3)**2)                         
       abs2 = dsqrt(k2(1)**2 + k2(2)**2 + k2(3)**2)                         
       costh = k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3)                     
       costh = costh/abs1/abs2                                             
       end                                                                           

                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION AM3(i1,i2,i3)
       implicit none
       integer  i1,i2,i3                                    
       complex*16  AM3                                              
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       AM3 = A(i1,i2)*B(i2,i3)                                              
       end                                                                    
                                                                           
       FUNCTION AM4(i1,i2,i3,i4)
       implicit none
       integer  i1,i2,i3,i4                                    
       complex*16  AM4                                              
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       AM4 = A(i1,i2)*B(i2,i4) + A(i1,i3)*B(i3,i4)                        
       end                                                                    
*


       FUNCTION AM5(i1,i2,i3,i4,i5)
       implicit none
       integer  i1,i2,i3,i4,i5                                    
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       complex*16  AM5                                              
       AM5 = A(i1,i2)*B(i2,i5) + A(i1,i3)*B(i3,i5) + A(i1,i4)*B(i4,i5)    
       end                                                                 
                                                                           
                                                                           
                                                                           
                                                                             
       FUNCTION BM3(i1,i2,i3)
       implicit none
       integer  i1,i2,i3                                     
       complex*16  BM3                                              
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       BM3 = B(i1,i2)*A(i2,i3)                                              
       end                                                                    

                                                                           
                                                                             
       FUNCTION BM4(i1,i2,i3,i4)                                         
       implicit none
       integer  i1,i2,i3,i4                                        
       complex*16  BM4                       
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       BM4 = B(i1,i2)*A(i2,i4) + B(i1,i3)*A(i3,i4)                        
       end                                                                    
                                                                           
                                                                             
       FUNCTION BM5(i1,i2,i3,i4,i5)                                         
       implicit none
       integer  i1,i2,i3,i4,i5                                    
       complex*16  BM5                                              
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
       BM5 = B(i1,i2)*A(i2,i5) + B(i1,i3)*A(i3,i5) + B(i1,i4)*A(i4,i5)    
       end                                                                    
                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION AM6(i1,i2,i3,i4,i5,i6)                                     
       implicit none
       integer  i1,i2,i3,i4,i5,i6                           
       complex*16  AM6                                             
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
                                                                           
       AM6 = A(i1,i2)*B(i2,i4)*A(i4,i6)                 
     &            + A(i1,i3)*B(i3,i4)*A(i4,i6)                 
     &              + A(i1,i2)*B(i2,i5)*A(i5,i6)               
     &           + A(i1,i3)*B(i3,i5)*A(i5,i6)                    
         end 

                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION BM6(i1,i2,i3,i4,i5,i6)                                     
       implicit none                             
       integer  i1,i2,i3,i4,i5,i6                           
       complex*16  BM6                                             
       COMPLEX*16 A,B,ZS
       real*8 s           
       COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
                                                                           
       BM6 = B(i1,i2)*A(i2,i4)*B(i4,i6)             
     &           + B(i1,i3)*A(i3,i4)*B(i4,i6)
     &             + B(i1,i2)*A(i2,i5)*B(i5,i6)
     &            + B(i1,i3)*A(i3,i5)*B(i5,i6)            
                end

                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
        FUNCTION LI2(x)                                                     
        implicit none                                                           
*      !! Dilogarithm for arguments x < = 1.0                             
        real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO  
        real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD                 
        real*8 Li2                                             
        integer  i                                                       
                                                                           
        DATA ZERO /0.0d0/, ONE /1.0d0/                               
        DATA HALF /0.5d0/, MALF /-0.5d0/                             
        DATA MONE /-1.0d0/, MTWO /-2.0d0/                            
        DATA PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/    
                                                                           
        DATA C( 0) / 0.4299669356081370d0/                              
        DATA C( 1) / 0.4097598753307711d0/                              
        DATA C( 2) /-0.0185884366501460d0/                              
        DATA C( 3) / 0.0014575108406227d0/                              
        DATA C( 4) /-0.0001430418444234d0/                              
        DATA C( 5) / 0.0000158841554188d0/                              
        DATA C( 6) /-0.0000019078495939d0/                              
        DATA C( 7) / 0.0000002419518085d0/                              
        DATA C( 8) /-0.0000000319334127d0/                              
        DATA C( 9) / 0.0000000043454506d0/                              
        DATA C(10) /-0.0000000006057848d0/                              
        DATA C(11) / 0.0000000000861210d0/                              
        DATA C(12) /-0.0000000000124433d0/                              
        DATA C(13) / 0.0000000000018226d0/                              
        DATA C(14) /-0.0000000000002701d0/                              
        DATA C(15) / 0.0000000000000404d0/                              
        DATA C(16) /-0.0000000000000061d0/                              
        DATA C(17) / 0.0000000000000009d0/                              
        DATA C(18) /-0.0000000000000001d0/                              
                                                                           
        if(x .gt. 1.00000000001d0) then                                    
          call crash('LI2')                                                
        elseif(x .gt. 1.0d0) then                                          
          x = 1.d0                                                      
        endif                                                              
                                                                           
        IF(X .EQ. ONE) THEN                                                
         LI2OLD=PI6
         LI2=LI2OLD                                                       
         RETURN                                                            
        ELSE IF(X .EQ. MONE) THEN                                          
         LI2OLD=MALF*PI6
         LI2=LI2OLD                                                  
         RETURN                                                            
        END IF                                                             
        T=-X                                                               
        IF(T .LE. MTWO) THEN                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
        ELSE IF(T .LT. MONE) THEN                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
        ELSE IF(T .LE. MALF) THEN                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
        ELSE IF(T .LT. ZERO) THEN                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
        ELSE IF(T .LE. ONE) THEN                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
        ELSE                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
        END IF                                                             
                                                                           
        H=Y+Y-ONE                                                          
        ALFA=H+H                                                           
        B1=ZERO                                                            
        B2=ZERO                                                            
        DO  I = 18,0,-1                                                    
          B0=C(I)+ALFA*B1-B2                                               
          B2=B1                                                            
          B1=B0                                                            
        ENDDO                                                              
        LI2OLD=-(S*(B0-H*B2)+A)                                           
              ! Artificial conversion                                      
*        Li2 = DReal(LI2OLD)   
         LI2=LI2OLD
        end                                        

                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION LN(arg1, arg2)                                             
       implicit none                                                           
*     !! log(-arg1/-arg2) for positive and negative arguments           
       complex*16   Ln                                           
*                                                                 1) nparton
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart
       real*8  arg1, arg2,pi                           
       parameter(pi = 3.14159265358979323846d0 )            
       if(arg1/arg2 .gt. 0.d0) then                                        
         ln = dcmplx(dlog(arg1/arg2),0.d0)                                               
       elseif(arg1.gt.0.d0 ) then                                              
         ln = dlog(dabs(arg1/arg2)) - imag* abspart*pi                       
       else                                                                
         ln = dlog(dabs(arg1/arg2))+imag*abspart*pi                       
       endif                                                               
       end                                                                    

                                                                           
                                                                           
                                                                           
       FUNCTION DILOG(arg1,arg2)                                           
       implicit none                                                           
*         !! dilogarithm of 1 - arg1/arg2 for pos and neg arguments         
                                                                           
       complex*16   Dilog,ln
       real*8 arg1,arg2,li2
       real*8 pi
       parameter(pi = 3.14159265358979323846d0 )    
*
       if(arg1/arg2 .gt. 0.d0) then                                        
         dilog = li2(1.d0 - arg1/arg2)
       else                                                                
         dilog = pi**2/6 - li2(arg1/arg2)                      
     &  - ln(arg1,arg2)*dlog(1.d0 - arg1/arg2)                      
       endif                                                               
       end                                                                    

                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION L0(arg1,arg2)                                              
       implicit none                                 
       complex*16   L0,ln                                           
       real*8 arg1, arg2                           
       real*8  arg                                              
                                                                           
       arg = arg1/arg2                                                     
       if( abs(1.d0 - arg) .gt. 0.00001d0) then                         
          L0 = ln(arg1,arg2)/(1.d0 - arg)                               
       else                                                                
*   !! in this case we never get an absorptive part                
          L0 = -1.d0 + 0.5d0*(arg - 1) - (arg - 1)**2/3.d0        
       endif                                                               
       end                                                                    

                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION L1(arg1,arg2)                                              
       implicit none                                 
       complex*16   L1,ln                                           
       real*8 arg1, arg2                           
       real*8  arg                                              
                                                                           
       arg = arg1/arg2                                                     
       if( abs(1.d0 - arg) .gt. 0.00001d0) then                         
         L1 = (ln(arg1,arg2) + 1 - arg)/(1.d0 - arg)**2                 
       else                                                                
         L1 = - 0.5d0 + (arg - 1)/3.d0 - 0.25d0*(arg - 1)**2      
       endif                                                               
       end                                                                     

                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION LSM12MHNEW(ss,tt,m1sq,m2sq)                              
       implicit none                                 
       real*8  ss,tt,m1sq,m2sq
       complex*16    lsm12mhnew,dilog
*                                                                 1) nparton
       complex* 16 imag
       real*8 zero,QES2,muscloop
       integer abspart
       common/nparton/imag,zero,QES2,muscloop,abspart
*
       real*8 pi
       parameter(pi = 3.14159265358979323846d0 )
       if( min(ss,m1sq,m2sq).lt. 0.d0 ) call crash('LSM12MHNEW')
       lsm12mhnew = - Dilog(m1sq,tt) - Dilog(m2sq,tt)    
     &     + 0.5d0*log(ss/m1sq)*log(ss/m2sq)
       if(tt .gt. 0.d0) then                                               
         lsm12mhnew = lsm12mhnew - 0.5d0*log(ss/tt)**2                   
       else                                                                
         lsm12mhnew = lsm12mhnew - 0.5d0*(log(-ss/tt)**2 - pi**2)        
         if(abspart.eq.1) lsm12mhnew =
     & lsm12mhnew + imag*pi*(log(ss) - log(-tt))                  
       endif                                                               
       end                                                                    

                                                                           
                                                                           
                                                                           
                                                                           
                                                                           
       FUNCTION I33M(s1,s2,s3)                                            
       implicit none                                 
       real*8 s1,s2,s3                            
       real*8  i33m, gram, xx, yy, lambda, rho,
     &            eta1, eta2, eta3,gram3,li2                               
       real*8 pi
       parameter(pi = 3.14159265358979323846d0 )
       if( min(s1,s2,s3).lt. 0.d0 ) call crash('I33M')                   
       gram = gram3(s1,s2,s3)                                              
                                                                           
       if(gram .gt. 0.d0) then                                             
         eta3 = max(s1,s2,s3)                                              
         eta2 = min(s1,s2,s3)                                              
         eta1 = s1 + s2 + s3 - eta2 - eta3                                 
         xx = eta1/eta3
         yy = eta2/eta3                                   
         lambda = dsqrt( (1 - xx - yy)**2 - 4*xx*yy )                       
         rho = 2.d0/(1 - xx - yy + lambda)                              
         i33m = - ( 2 * ( li2(-rho*xx) + li2(-rho*yy) )   
     &  + log(rho*xx)*log(rho*yy)                
     &  + log(yy/xx)*log((1+rho*yy)/(1+rho*xx))   
     &        + Pi**2/3   )/eta3/lambda                             
       else                                                                
         print*,' I33M not programmed for these arguments '               
         print*,' s1,s2,s3 : ', s1,s2,s3                                   
         stop                                                              
       endif                                                               
       end                                                                    

                                                                           

* 
       FUNCTION NA(k1,k2)                        
        implicit none                                                       
        real *8  k1(4),k2(4),k1m(4),k2m(4)                        
        complex *16  na,nna
         complex* 16 imag
         parameter( imag = (0.0d0 ,1.0d0 ))
         integer i
*                                                                            
        do i=1,4
        k1m(i)=-k1(i)    
        k2m(i)=-k2(i)
        enddo    
        if(k1(4). gt. 0. .and. k2(4) .gt. 0. ) then                     
         na=nna(k1,k2)
        elseif(k1(4).gt. 0. .and. k2(4) .lt.0. ) then                 
          na = imag*nna(k1,k2m)                                          
        elseif(k1(4) .lt. 0. .and. k2(4) .gt. 0. ) then                 
          na = imag*nna(k1m,k2)                                          
        elseif(k1(4).lt. 0. .and. k2(4) .lt. 0. ) then                 
          na = -nna(k1m,k2m)                                             
        else                                                                
          call crash('NA')                     
        endif             
*                                                                            
        END 
*
*
       FUNCTION NNB(k1,k2)
c                                                                              
      implicit none                                                      
      real *8  k1(4),k2(4)                       
      complex *16 nnb,sq1,sq2
      real*8 Tiny                              
         complex* 16 imag
         parameter( imag = (0.0 ,1.0 ))
         external Tiny
c                                                                         
      if(k1(4) .gt. 0.  .and. k2(4) .gt. 0.  ) then                    
        sq1 = (k1(3) + imag*k1(2) + Tiny(k1(4)))/
     &    (Dsqrt(k1(3)**2 + k1(2)**2) + Tiny(k1(4)))                  
        sq2 = (k2(3) + imag*k2(2) + Tiny(k2(4)))/
     &        (Dsqrt(k2(3)**2 + k2(2)**2) + Tiny(k2(4))) 
        nnb = - Dsqrt(k1(4) - k1(1))*Dsqrt(k2(4) + k2(1))/sq1
     &         + Dsqrt(k1(4) + k1(1))*Dsqrt(k2(4) - k2(1))/sq2           
      else                                                               
        call crash('NNB')                                                 
      endif     
*
      END
   
       FUNCTION Tiny(x)
      real*8 Tiny,x
      if(dabs(x).lt.1d-12)x=1d-12
      tiny=(1d-24)*dabs(x)
      end
                                                          
       FUNCTION NB(k1,k2)
      implicit none                                       
      real *8 k1(4),k2(4),k1m(4),k2m(4)                                 
      complex *16 nb,nnb                          
         complex* 16 imag
         parameter( imag = (0.0 ,1.0 ))
      integer i
*
        do i=1,4
        k1m(i)=-k1(i)    
        k2m(i)=-k2(i)
        enddo    
      if(k1(4) .gt. 0.  .and. k2(4) .gt. 0.  ) then       
       nb=nnb(k1,k2)                                      
      elseif(k1(4) .gt. 0.  .and. k2(4) .lt. 0.  ) then   
        nb = imag*nnb(k1,k2m)                             
      elseif(k1(4) .lt. 0.  .and. k2(4) .gt.  0.  ) then  
        nb = imag*nnb(k1m,k2)                             
      elseif(k1(4) .lt. 0.  .and. k2(4) .lt. 0.  ) then   
        nb = -nnb(k1m,k2m)                                
      else                                                
        call crash('NB')                                  
      endif                                               
                                                          
      END
*
       FUNCTION NNA(k1,k2)
c
        implicit none
        real *8  k1(4),k2(4)
        complex *16 nna,sq1,sq2 
         complex* 16 imag
         parameter( imag = (0.0 ,1.0 ))
         real*8 Tiny
         external Tiny
c
       if(k1(4) .gt. 0.  .and. k2(4) .gt. 0.  ) then
c 
       sq1=dcmplx(k1(3)+Tiny(k1(4)),k1(2))
       sq1 = sq1/(Dsqrt(k1(3)**2 + k1(2)**2) + Tiny(k1(4)))

       sq2=dcmplx(k2(3)+Tiny(k2(4)),k2(2))
       sq2 = sq2/
     &             (Dsqrt(k2(3)**2 + k2(2)**2) + Tiny(k2(4)))
       nna =  dsqrt((k1(4) - k1(1))*(k2(4) + k2(1)))*sq1 
       nna=nna- dsqrt((k1(4) + k1(1))*(k2(4) - k2(1)))*sq2
c
      else
       call crash('NNA')
      endif
      End

      SUBROUTINE CRASH(functionname)
      implicit none
      character*50  functionname
      write(6,*) 'Program crashes because of a call to the function ', 
     &          functionname
      stop
      END
 
********************* () fillquin ************************************
* and  functions na,nb,nna,nnb, tiny and crash
* adrian's type verion for calculating the inner products
* it agrees with quinn 
        SUBROUTINE fillquin(pcm,a,b,s,zs)
        implicit none
        real*8 pcm(4,7),s(7,7)
        complex*16 a(7,7),b(7,7),zs(7,7),na,nb
        integer i,j
*   
        do i=1,6
        a(i,i)=(0.d0,0.d0)
        b(i,i)=(0.d0,0.d0)
        zs(i,i)=(0.d0,0.d0)
        s(i,i)=0.d0
        enddo
        do i=1,5
        do j=1+i,6
        a(i,j)=na(pcm(1,i),pcm(1,j))
        b(i,j)=nb(pcm(1,i),pcm(1,j))
        zs(i,j)=-a(i,j)*b(i,j)
        s(i,j)=dreal(zs(i,j))
        enddo
        enddo
*
        do i=1,5
        do j=i+1,6
        b(j,i)=-b(i,j)
        a(j,i)=-a(i,j)
        zs(j,i)=zs(i,j)
        s(j,i)=s(i,j)
        enddo
        enddo

        return
        end
* package phase_pack.f ww,zz,wz production
* read-only
* routines :
* SUBROUTINE THREEBODYWWGFIX(etot,qm,mx,my,zm,y,f, pr3,weight,rancosww,ranphiww)
* SUBROUTINE TWOBOfix(etot,XM,YM,PR,WEIGHT,rancosww,ranphiww)
* SUBROUTINE BOOSTG(N,WMASS,PBOO,PIN,POUT)                         
* SUBROUTINE GDECWfix(PR,PL,XM,WEIG,cran,fran)                   
* SUBROUTINE CHANGWall(NPART,plaball,pcmall)    
* SUBROUTINE QINN(NPART,PPCM,H,CH,ZD,d)
******************************************************* 
      SUBROUTINE THREEBODYWWGFIX(
     & etot,qm,mx,my,zm,y,f, pr3,weight,rancosww,ranphiww)

*
       implicit none
       real*8 etot,qm,y,f,mx,my,zm,pr3(5,3),weight,rancosww,
     & ranphiww
*
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart


      real*8 pi
      parameter ( pi = 3.14159265358979323846d0 )
      integer i,j,nout
      real*8 s,eq,ptotq,sc,pr2(5,2),pr(5,2),
     & pin(4,2),pout(4,2),weight2
*
      if((etot-qm).lt.etot*1d-24)then
      qm=etot
      endif
      SC=DSQRT(1.D0-y**2)
      s=etot**2
      eq=(s+qm**2-zm**2)/(2D0*etot)                                   
      PR2(4,1)=eq 
      ptotq=etot-eq
      PR2(5,1)=ptotq

      PR2(4,2)=ptotq
      PR2(3,1)=-PR2(5,1)*y                                                
      PR2(5,2)=PR2(5,1)                                                  
      PR2(2,1)=-PR2(5,1)*SC*DCOS(F)                                       
      PR2(1,1)=-PR2(5,1)*SC*DSIN(F)                                       
      DO  I=1,3 
      PR2(I,2)=-PR2(I,1)
      enddo 
      do i=1,5
      pr3(i,3)=pr2(i,2)                                                
      enddo
*
      call TWOBOfix(qm,mx,my,PR,WEIGHT2,rancosww,ranphiww) 
      do i=1,2
      do j=1,4
      pin(j,i)=pr(j,i)
      enddo                       
      enddo                       
*           PRINT 2007,((PR(I,J),I=1,5),J=1,2),qm            
* 2007      FORMAT(' Pr',/,2(2X,5E16.8,/),'  qm:',E15.7,/)
* boost the W-bosons from the rest frame of the pair to the moving frame
      nout=2
      CALL BOOSTG(nout,Qm,PR2(1,1),PIN(1,1),Pout(1,1))
      do i=1,2
      do j=1,4
      pr3(j,i)=pout(j,i)
      enddo          
      pr3(5,i)=dsqrt(pout(1,i)**2+pout(2,i)**2+pout(3,i)**2)             
      enddo                       
*      PRINT 2008,((PR3(I,J),I=1,5),J=1,3)
* 2008      FORMAT(' Pr3',/,3(2X,5E15.7,/))
* boost the W-bosons from the rest frame of the pair to the moving frame
      weight=weight2
      return
      end      
*
      SUBROUTINE TWOBOfix(etot,XM,YM,PR,WEIGHT,rancosww,ranphiww)
      IMPLICIT none
      real*8 etot,xm,ym,pr(5,2),weight,rancosww,ranphiww
*
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart

*
      real*8 pi,xm1,xm2
      parameter ( pi = 3.14159265358979323846d0 )
      integer i
      real*8 f,c,sc
*
      xm1=xm**2
      xm2=ym**2                                                      
      F=2D0*PI*ranphiww
      c=2.d0*rancosww-1.d0                                                   
      SC=DSQRT(1.D0-C*C)                                               
      PR(4,1)=(etot**2+XM1-XM2)/(2D0*ETOT)                                   
      PR(5,1)=DSQRT(PR(4,1)*PR(4,1)-XM1)                               
      PR(4,2)=ETOT-PR(4,1)                                             
      PR(3,1)=PR(5,1)*C                                                
      PR(5,2)=PR(5,1)                                                  
      PR(2,1)=PR(5,1)*SC*DCOS(F)                                       
      PR(1,1)=PR(5,1)*SC*DSIN(F)                                       
      DO  I=1,3 
      PR(I,2)=-PR(I,1)
      enddo                                                 
      WEIGHT=PI*PR(5,1)/ETOT                                   
      RETURN                                                           
      END                                                              
*
*
      SUBROUTINE BOOSTG(N,WMASS,PBOO,PIN,POUT)                         
      implicit none                                          
      integer n,i,j
      real*8 wmass,pboo,pin,pout,q0,qm,e,f
      DIMENSION PBOO(5),PIN(4,20),POUT(4,20)                           
*
      Q0=PBOO(4)
      QM=WMASS 
      DO 1 I=1,N                                                       
      E=PIN(4,I)                                                       
      POUT(4,I)=(Q0*E+PBOO(1)*PIN(1,I)+PBOO(2)*PIN(2,I)+               
     X          PBOO(3)*PIN(3,I))/QM                                   
      F=(E+POUT(4,I))/(Q0+QM)                                          
      DO 2 J=1,3                                                       
2     POUT(J,I)=PIN(J,I)+PBOO(J)*F                                     
1     CONTINUE 
      RETURN                                                           
      END                                                              
*                                                                      
C                                                                      
       SUBROUTINE GDECWfix(PR,PL,XM,WEIG,cran,fran)                   
       implicit none   
       real*8 PR(5),PL(4,2),xm,weig,cran,fran
*
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart

*
      real*8 pin(4,2)
      integer nout,i
      real*8 pi,c,f,sc                                
      parameter ( pi = 3.14159265358979323846d0 )
      c=2.d0*cran-1
      f=2.d0*pi*fran
      SC=DSQRT(1D0-C*C)  
      PIN(4,1)=0.5*XM                                                  
      PIN(3,1)=PIN(4,1)*C                                              
      PIN(2,1)=PIN(4,1)*SC*DCOS(F)                                     
      PIN(1,1)=PIN(4,1)*SC*DSIN(F)                                     
      PIN(4,2)=PIN(4,1)                                                
      DO  I=1,3                                                       
      PIN(I,2)=-PIN(I,1)                                               
      enddo
      nout=2
      CALL BOOSTG(nout,XM,PR(1),PIN(1,1),PL(1,1))
      weig=0.5d0*pi                         
      RETURN                                                           
      END                                                              
*
      SUBROUTINE CHANGWall(NPART,plaball,pcmall)    
* 
      implicit none
      integer npart
      real*8 plaball(4,7,4),pcmall(4,7,4)
      REAL*8 P(4,7)
      integer ii,i,j,k,init                            
      DATA INIT/0/                                                     
* 
      do ii=1,4
      DO 10 I=1,NPART                                                  
      P(1,I)=PLABall(1,I,ii)                                                 
      P(2,I)=-PLABall(2,I,ii)                                                 
      P(3,I)=PLABall(3,I,ii)                                                 
      P(4,I)=PLABall(4,I,ii)                                                 
  10  CONTINUE                                                         
      DO 20 J=1,4                                                      
      DO 30 K=3,NPART                                              
  30  PCMALL(J,K,ii)=P(J,K)                                                
      PCMALL(J,1,ii)=-P(J,1)                                           
      PCMALL(J,2,ii)=-P(J,2)                                             
  20  CONTINUE     
      enddo                                                    
      IF(INIT.gt.-2) RETURN                                             
* test print
* ---------
*      INIT=Init+1
      do ii=1,4
      write(6,*)'in change, ii=',ii
      PRINT 66,((PLABALL(I,J,ii),I=1,4),J=1,NPART)                           
  66  FORMAT(7(2X,4G13.7,/))                                           
      PRINT 77,((PCMALL(I,J,ii),I=1,4),J=1,NPART)                            
  77  FORMAT(7(2X,4G13.7,/))
      enddo
      RETURN                                                           
      END                                                              
C                                                                      
*
      SUBROUTINE QINN(NPART,PPCM,H,CH,ZD,d)
      implicit none
*
C PCM SHOULD BE DEFINED SUCH THAT ALL 4-MOMENTA ARE OUTGOING.          
C CONVENTION FOR PCM AND P IS THAT DIRECTION 1 =BEAM, COMPONENT        
C 4 = ENERGY AND COMPONENT 2 AND 3 ARE TRANSVERSE COMPONENTS.          
C THUS INCOMING MOMENTA SHOULD CORRESPON TO OUTGOING MOMENTA           
C OF NEGATIVE ENERGY. 
*
* in new notation:H->A,CH->B, ZD->xomplex S, D-> S 
*
c
C PCM IS FILLED BY PHASE SPACE MONTE CARLO.                            
      COMPLEX*16 PT5,ZT,Z1,ZI,ZP,ZQ,ZPS,ZQS,H,CH,ZD                  
      COMPLEX*16 ZDPM,ZDMP                                             
      real*8 d,ppcm,p,wrn,eps,p1,q1,pp,qp,qm,q2,p2,
     & pt,qt,dpm,pti,pm,dmp,qti
      DIMENSION H(7,7),CH(7,7),zd(7,7)                                  
      DIMENSION PPCM(4,7),D(7,7)                                              
      DIMENSION P(4,7)                                                 
      DIMENSION WRN(7)
      integer npart,i,j,l,ij,ii,jj,ip1,ipp1                         
*
      EPS=1.d-40                                                    
      ZI=DCMPLX(0.D0,1.D0)                                             
      Z1=DCMPLX(1.D0,0.D0)
                                                   
C FOLLOWING DO LOOP IS TO CONVERT TO OUR STANDARD INDEXING             
      DO 1 L=1,NPART                                                   
      DO 1 IJ=1,4                                                      
 1       P(IJ,L)=PPCM(IJ,L)                                               

      DO 2 II=1,NPART                                                      
      WRN(II)=1.D0                                                     
      IF(P(4,II).LT.0.D0) WRN(II)=-1.D0                                
      DO 2 JJ=1,4                                                      
      P(JJ,II)=WRN(II)*P(JJ,II)                                        
    2 CONTINUE                                                         
C THE ABOVE CHECKS FOR MOMENTA WITH NEGATIVE ENERGY,INNER PRODUCTS     
C ARE EXPRESSED DIFFERENTLY FOR DIFFERENT CASES                        
      DO 11 I=1,NPART-1                                                
      IP1=I+1                                                          
      DO 11 J=IP1,NPART                                                
C      PRINT *,I,J                                                     
      Q1=P(4,I)+P(1,I)                                                 
      QP=0.d0                                                           
      IF(Q1.GT.EPS*p(4,i))QP=Q1                                        
      Q2=P(4,I)-P(1,I)                                                 
      QM=0.d0                                                           
      IF(Q2.GT.EPS*p(4,i))QM=Q2                                        
      P1=P(4,J)+P(1,J)                                                 
      PP=0.d0                                                            
      IF(P1.GT.EPS*p(4,i))PP=P1                                        
      P2=P(4,J)-P(1,J)                                                 
      PM=0.d0                                                            
      IF(P2.GT.EPS*p(4,j))PM=P2                                        
      DMP=dsqrt(PM*QP)                                                        
      ZDMP=DCMPLX(DMP,0.D0)                                            
      DPM=dsqrt(PP*QM)                                                        
      ZDPM=DCMPLX(DPM,0.D0)                                            
C NOTE THAT IN OUR INNER PRODUCT NOTATION WE ARE COMPUTING <P,Q>       
      PT=DSQRT(P(2,J)**2+P(3,J)**2)                                    
      QT=DSQRT(P(2,I)**2+P(3,I)**2)                                    
      IF(PT.GT.EPS*p(4,j)) GO TO 99                                           
      ZP=Z1                                                            
      GO TO 98                                                         
   99 PTI=1.D0/PT                                                      
      ZP=DCMPLX(PTI*P(3,J),PTI*P(2,J))                                 
   98 ZPS=DCONJG(ZP)                                                   
      IF(QT.GT.EPS*p(4,i)) GO TO 89                                           
      ZQ=Z1                                                            
      GO TO 88                                                         
   89 QTI=1.D0/QT                                                      
      ZQ=DCMPLX(QTI*P(3,I),QTI*P(2,I))                                 
   88 ZQS=DCONJG(ZQ)                                                   
      ZT=Z1                                                            
      IF(WRN(I).LT.0) ZT=ZT*ZI                                         
      IF(WRN(J).LT.0) ZT=ZT*ZI                                         
      H(J,I)=(ZDMP*ZP-ZDPM*ZQ)*ZT                                      
      CH(I,J)=(ZDMP*ZPS-ZDPM*ZQS)*ZT                                   
      ZD(j,i)=H(J,I)*CH(I,J)                                                
      PT5=DCMPLX(.5D0,0.D0)                                            
      D(J,I)=Dreal(ZD(j,i))                                                    
*      write(1,2003),J,I,H(J,I),CH(I,J),D(J,I)                            
*2003  FORMAT(2X,'I,J,H(J,I),CH(I,J),D(J,I):   ',
*     & 2I2,/,4(2X,4d19.9,/))                
   11 CONTINUE                                                         
      DO  I=1,NPART-1                                                
      IPP1=I+1                                                         
      DO  J=IPP1,NPART                                               
      H(I,J)=-H(J,I)                                                   
      CH(J,I)=-CH(I,J)                                                 
      zd(i,j)=zd(j,i)
      D(I,J)=D(J,I) 
      enddo
      enddo
      RETURN                                                           
      END                                                              
c
c
c
c End of packanomDKS.f
c
c
c
c
c Begin of wzanDKS.f
c
c
          SUBROUTINE TMATR7NEW(xi,ximax,sigmat,sigmatb)
* sigmat(1 or 2,j): uubar->ww, sigmat(2,j):ddbar->ww 
* sigmatb(1 or 2,j):ubaru->ww, sigmatb(2,j):dbard->ww
* sigmat(iso,j):j=1,2,3,4: loop, born, collp,collm
          implicit none
          real*8 xi,ximax,sigmat(6,4),sigmatb(6,4)
          integer npart
*                                                   7)pmom
      REAL*8 plaball,pcmall,pwwg
      COMMON/PMOM/ PLABall(4,7,4),PCMall(4,7,4), pwwg(5,3,4)             
*                                                    10)qinne

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*
      real*8 constwz,ctgtetw
      common/constGG/constwz,ctgtetw
*
      integer i,ii
      real*8 norm(6)
*
       npart=7
*
      do i=1,6
      if(i.lt.3)then
      norm(i)=constwz/nc**2
      else
      norm(i)=constwz/(nc*(nc**2-1))
      endif
      enddo
* -------------------------------------------------------------
* full
      ii=1
      if(xi.lt.ximax)then
       call sigtyp(ii,npart,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
      endif
* soft
      ii=2
      call sigtyp(ii,npart-1,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
* collp and collm
      do ii=3,4
      if(xi.lt.ximax)then
      call sigtyp(ii,npart-1,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
      endif
      enddo
      return
      end    

*
*
          SUBROUTINE TMATR6NEW(xi,ximax,sigmat,sigmatb)
* sigmat( 1,j): 0, sigmat(2,j):dubar->ww 
* sigmatb(1 or 2,j):ubaru->ww, sigmatb(2,j):dbard->ww
* sigmat(iso,j):j=1,2,3,4: loop, born, collp,collm
          implicit none
          real*8 xi,ximax,sigmat(2,4),sigmatb(2,4)
          integer npart
*                                                   7)pmom
      REAL*8 plaball,pcmall,pwwg
      COMMON/PMOM/ PLABall(4,7,4),PCMall(4,7,4), pwwg(5,3,4)             
*                                                    10)qinne

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*
      real*8 ubd4l,sww,s12
*                                                          ??) subcont
      logical bornonly,qqonly,gqonly
      common/subcont/ bornonly,qqonly,gqonly
*
      real*8 constwz,ctgtetw
      common/constGG/constwz,ctgtetw
*
      complex*16 fcl
      integer ii
      real*8 normqq

*
      npart=6

c aohcms      CALL CHANGWALL(NPART+1,plaball,pcmall)

* -------------------------------------------------------------
           ii=1
c aohcms           call qinn(npart,pcmall(1,1,ii+1),a,b,zs,s)                       
c aohcms           s12 =k1pk2sq(plaball(1,1,ii+1),plaball(1,2,ii+1),1)
           s12 = sigmat(1,1) ! special s12 passing from dks wrapper sh = (p_had1 + p_had2)^2 
           sww=s12
           
* ------------------------------------------------------
* couplings introduced by lance (eq. 17,19)
* ----------------------------------------
      fcl=-dcmplx(ctgtetw*sww/(sww-mw**2))
* if uubar intitial state iso=1 if ddbar initital state iso=2
* for ubaru and dbard we have the same but in sigmatb
* if loop ii=1
* if born or soft ii=2
      normqq=constwz/(s12*nc**2)

c     ii = 1 : loop
c     ii = 2 : born
      do ii=1,2
      sigmat(1,ii)=0.0
      sigmatb(1,ii)=0.0
        if((bornonly.eqv..true.).and.(ii.eq.1))goto 1
*d_1-ubar_2
      sigmat(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &    fcl,2,1,3,4,5,6,7)*normqq

*ubar_1-d_2
      sigmatb(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &    fcl,1,2,3,4,5,6,7)*normqq
 1    continue
      enddo
      if(bornonly.eqv..true.)goto 2
      if(xi.lt.ximax)then
* ii=3 
           ii=3
c aohcms           call qinn(npart,pcmall(1,1,ii),a,b,zs,s)                       
c aohcms           s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c aohcms           sww=s12
      fcl=-dcmplx(ctgtetw*sww/(sww-mw**2))
      normqq=constwz/(s12*nc**2)
*
      sigmat(1,ii)=0.0
      sigmatb(1,ii)=0.0
*d_1-ubar_2
      sigmat(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &     fcl,2,1,3,4,5,6,7)*normqq
*ubar_1-d_2
      sigmatb(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &     fcl,1,2,3,4,5,6,7)*normqq

      ii=4
c aohcms      call qinn(npart,pcmall(1,1,ii),a,b,zs,s)                       
c aohcms      s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c aohcms      sww=s12
      fcl=-dcmplx(ctgtetw*sww/(sww-mw**2))
*
      normqq=constwz/(s12*nc**2)
       sigmat(1,ii)=0.0
       sigmatb(1,ii)=0.0
       sigmat(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &     fcl,2,1,3,4,5,6,7)*normqq
       sigmatb(2,ii)=ubd4l(npart,ii,sww,el,er,dl,ul,
     &     fcl,1,2,3,4,5,6,7)*normqq
      endif
*
 2    continue

      return
      end    
*
      subroutine sigtyp(ii,np,norm,pcmall,plaball,pwwg,
     & sig, sigb)
      implicit none
      integer ii,np
      real*8 norm(6),pcmall(4,7,4),plaball(4,7,4),
     & pwwg(5,3,4)
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*                                                    10)qinne
      real*8 constwz,ctgtetw
      common/constGG/constwz,ctgtetw

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*
      real*8 ubd4l,normqq(6),
     & sig(6),sigb(6),s12,sww
      complex *16 fcl
      integer i
* -------------------------
c aohcms hack
      s12 = sig(1)         ! special s12 passing from dks wrappper sh = (p_had1 + p_had2)^2 

      if (ii.eq.1.and.np.eq.7) then
         sww = sig(2)      ! special sww passing from dks wrapper sh = (p_v1 + p_v2)^2 
      else
         sww = s12
      endif
c aohcms hack end

c commented out
c      s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c        if(ii.eq.1.and.np.eq.7)then
c         do i=1,4
c            pw1(i)=plaball(i,3,ii)+plaball(i,4,ii)
c            pw2(i)=plaball(i,5,ii)+plaball(i,6,ii)
c         enddo
c          sww =k1pk2sq(pw1(1),pw2(1),1)
c        else
c           sww=s12
c        endif
c       call qinn(np,pcmall(1,1,ii),a,b,zs,s)               
c end commented out        
* -------------------------------
      do i=1,6
      normqq(i)=norm(i)/s12
      sig(i)=0
      sigb(i)=0
      enddo
*
      fcl=-dcmplx(ctgtetw*sww/(sww-mw**2))
*
* for any ii
*u_1-dbar_2
      sig(1)=0.0
*d_1-ubar_2
      sig(2)=ubd4l(np,ii,sww,el,er,dl,ul,
     &     fcl,2,1,3,4,5,6,7)*normqq(2)
*
*ubar_1-d_2
      sigb(2)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,1,2,3,4,5,6,7)*normqq(1)
*dbar_1-d_2
      sigb(1)=0.0
*
      if(ii.eq.1)then
*g_1-ubar_2
      sig(3)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,2,7,3,4,5,6,1)*normqq(3)
*g_1-dbar_2
      sig(4)=0.
*g_1-u_2
      sig(5)=0
*g_1-d_2
      sig(6)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,7,2,3,4,5,6,1)*normqq(6)
*u_1-g_2
      sigb(5)=0.0
*d_1-g_2
      sigb(6)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,7,1,3,4,5,6,2)*normqq(4)
*ubar_1-g_2
      sigb(3)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,1,7,3,4,5,6,2)*normqq(5)
*dbar_1-g_2
      sigb(4)=0.0
      endif
* in the collinear limit only quark initial states
      if(ii.eq.3)then
*g_1-ubar_2-->(dbar)d_1-ubar_2
      sig(3)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,2,1,3,4,5,6,7)*normqq(1)
*g_1-dbar_2-->(ubar)u_1-dbar_2
      sig(4)=0.0
*g_1-u_2-->(d)dbar_1-u_2
      sig(5)=0.0
*g_1-d_2-->(u)ubar_1-d_2
      sig(6)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,1,2,3,4,5,6,7)*normqq(1)

      endif 
*
      if(ii.eq.4)then
*ubar_1-g_2->ubar_1-d_2(dbar) 
      sigb(3)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,1,2,3,4,5,6,7)*normqq(1)
*dbar_1-g_2->d_1-ubar_2(u)
      sigb(4)=0
*u_1-g_2-> u_1-dbar_2(d) 
      sigb(5)=0.0
*d_1-g_2->d_1-ubar_2(u)
      sigb(6)=ubd4l(np,ii,sww,el,er,dl,ul,
     &    fcl,2,1,3,4,5,6,7)*normqq(1)
      endif
      return
      end    
*
      FUNCTION ubd4l(npart,ii,sww,el,er,dl,ul,fcl,
     &   i1,i2,i3,i4,i5,i6,i7)
            IMPLICIT NONE
            INTEGER npart,ii,i1,i2,i3,i4,i5,i6,i7
            REAL*8 el,er,dl,ul,ubd4l,sww
            COMPLEX*16 fcl
* the POSITION of i1,i2,i3,i4,i5,i6,i7 in the function call
* corresponds to
* ubar(1)^-u(2)^+l(3)^-nubar(4)^+lpbar(5)^+nup(6)^-g(7)^+
* where all particles are outgoing and the upper label refers
* to helicity (see eq. (2.6) in DKS.
* the ACTUAL VALUE of i1,i2,i3,i4,i5,i6,i7 is the 
* the momentum label in the physical process
* if npart=6, the gluon momenta=0 and the 
* label i7 is dummy in a sense 
* from the primitive amplitude we calculate
* (see
*     m(i)_tree_n=el**2*|dl*ampla(n,ii,1,2,3,4,5,6,7)+
* ul*ampla(n,ii,1,2,6,5,4,3,7)+fcl*
* amplb(n,ii,1,2,3,4,5,6,7)|^2+(el<->er,5<->6)
* =el**2*|amplwz(dl,ul,ctgtetw,s12,mw2,
*     &      npart,ii,N1,N2,N3,N4,N5,N6,n7)|^2
* +er**2*|amplwz(dl,ul,ctgtetw,s12,mw2,
*     &      npart,ii,N1,N2,N3,N4,N6,N5,n7)|^2
*
      real*8 gluneghel
      complex*16 a1,a2,amplwz
      complex*16 b1,b2
      integer gluhel
*
      ubd4l=0
      gluneghel=0.d0
      gluhel=1
* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

* for ubarq
      a1=amplwz(dl,ul,fcl,sww,
     & npart,ii,i1,i2,i3,i4,i5,i6,i7,gluhel)
      a2=amplwz(dl,ul,fcl,sww,
     & npart,ii,i1,i2,i3,i4,i6,i5,i7,gluhel)
      b1=a1
      b2=a2

      if(npart.eq.6.and.ii.eq.1)then
*     loop amplitudes
         b1=amplwz(dl,ul,fcl,sww,
     &        npart,ii+1,i1,i2,i3,i4,i5,i6,i7,gluhel)
         b2=amplwz(dl,ul,fcl,sww,
     &        npart,ii+1,i1,i2,i3,i4,i6,i5,i7,gluhel)
      endif
      ubd4l=el**2*dreal(a1*dconjg(b1))
     &     + er**2*dreal(a2*dconjg(b2))
      gluneghel=0
      if(npart.eq.7.and.ii.eq.1)then
         gluhel=-1
*     negative helicity gluons can also contribute 
         a1=amplwz(dl,ul,fcl,sww,
     &        npart,ii,i1,i2,i3,i4,i5,i6,i7,gluhel)
         a2=amplwz(dl,ul,fcl,sww,
     &        npart,ii,i1,i2,i3,i4,i6,i5,i7,gluhel)
         b1=a1
         b2=a2
         gluneghel=el**2*dreal(a1*dconjg(b1))
     &        + er**2*dreal(a2*dconjg(b2))
      endif
      ubd4l=ubd4l+gluneghel
      return 
      end
*     

*


      FUNCTION Amplwz(dl,ul,fcl,sww,
     &      npart,ii,N1,N2,N3,N4,N5,N6,n7,gluhel)
* npart=6,ii=1: loop
* npart=6,ii=2: born
* npart=7,ii=1,2,3,4:gluon bremsstrahlung: full,soft,collp,collm  
      implicit none 
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7,gluhel
      complex*16 Amplwz,AmplA,AmplBwzan,fcl
      complex*16 FAmplA,FAmplBwzan
      real*8 dl,ul,sww
      integer wcharge
      common/wcharge/wcharge
*...
       if(wcharge.eq.-1)then
         if(gluhel.eq.1)then
       Amplwz=dcmplx(dl)* AmplA(npart,ii,N1,N2,N3,N4,N5,N6,n7)
     &  + dcmplx(ul)*AmplA(npart,ii,N1,N2,N6,N5,N4,N3,n7)
     &  +fcl* AmplBwzan(npart,ii,sww,N1,N2,N3,N4,N5,N6,n7)
          endif
         if(gluhel.eq.-1)then
       Amplwz=-dcmplx(dl)* FAmplA(npart,ii,N2,N1,N5,N6,N3,N4,n7)
     &  - dcmplx(ul)*FAmplA(npart,ii,N2,N1,N4,N3,N6,N5,n7)
     &  +fcl* FAmplBwzan(npart,ii,sww,N2,N1,N4,N3,N6,N5,n7)
          endif
       endif
       if(wcharge.eq.1)then
        if(gluhel.eq.1)then
       Amplwz=dcmplx(ul)* AmplA(npart,ii,N1,N2,N4,N3,N6,N5,n7)
     &  + dcmplx(dl)*AmplA(npart,ii,N1,N2,N5,N6,N3,N4,n7)
     &  -fcl* AmplBwzan(npart,ii,sww,N1,N2,N4,N3,N6,N5,n7)
        endif
        if(gluhel.eq.-1)then
       Amplwz=-dcmplx(ul)* FAmplA(npart,ii,N2,N1,N6,N5,N4,N3,n7)
     &  - dcmplx(dl)*FAmplA(npart,ii,N2,N1,N3,N4,N5,N6,n7)
     &  -fcl* FAmplBwzan(npart,ii,sww,N2,N1,N3,N4,N5,N6,n7)
        endif
       endif

       return
       end
*
      FUNCTION fAmplA(npart,ii,N1,N2,N3,N4,N5,N6,n7)
* npart=6,ii=1: loop
* npart=6,ii=2: born
* npart=7,ii=1,2,3,4:gluon bremsstrahlung: full,soft,collp,collm  
      implicit none 
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 fAmplA
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag,aa1,aa2,bm3,t134,t256
* //////////////////////////////////
*
      imag = (0.0 ,1.d0 )
      if(ii.eq.1.and.npart.eq.7) then
      t134=zs(n1,n3)+zs(n1,n4)+zs(n3,n4)
      t256=zs(n2,n5)+zs(n2,n6)+zs(n5,n6)
      aa1=b(n1,n3)*a(n3,n4)*a(n2,n5)*(bm3(n6,n2,n7)+bm3(n6,n5,n7) )
     &  /t256
      aa2=(bm3(n6,n1,n4)+bm3(n6,n3,n4))*(bm3(n1,n2,n5)+bm3(n1,n7,n5))
     &/b(n7,n2)
      fAmplA=imag*b(n1,n3)*(aa1+aa2)
     &        /(b(n1,n7)*zs(n3,n4)*zs(n5,n6)*t134)
      endif 
      return
      end 
*

      FUNCTION AmplBwzan(npart,ii,sww,N1,N2,N3,N4,N5,N6,n7)
      implicit none 
      real*8 sww
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 AmplBwzan,aloopbwzan,atree6bwzan,atree7bwzan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      if(ii.eq.1.and.npart.eq.6)then
       AmplBwzan=aloopbwzan(sww,n1,n2,n3,n4,n5,n6)
      endif
      if(ii.eq.1.and.npart.eq.7)then
      AmplBwzan=atree7bwzan(sww,n1,n2,n3,n4,n5,n6,n7)
      endif
      if(ii.ne.1)then
       amplbwzan=atree6bwzan(sww,n1,n2,n3,n4,n5,n6)
      endif
      return
      end 

      FUNCTION FAmplBwzan(npart,ii,sww,N1,N2,N3,N4,N5,N6,n7)
      implicit none 
      real*8 sww
      integer npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 fAmplBwzan,fatree7bwzan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
      complex*16 imag
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      if(ii.ne.1.and.npart.ne.7)then
      write(6,*)'wrong ii,npart:',ii,npart
      stop
      endif
      if(ii.eq.1.and.npart.eq.7)then
      fAmplBwzan=fatree7bwzan(sww,n1,n2,n3,n4,n5,n6,n7)
      endif
      return
      end 
*
*

      FUNCTION Atree7Bwzan(sww,N1,N2,N3,N4,N5,N6,n7)
      implicit none 
      real*8 sww
      integer n1,n2,n3,n4,n5,n6,n7
      complex*16 Atree7Bwzan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*..................................                   )anomcpz
      real*8 g1zcp,kapzcp,lamzcp,mwcp,formlam
      common/anomcpz/g1zcp,kapzcp,lamzcp,mwcp,formlam
* .............  end of common block .............
      real*8 s127,gcp1,gcp2,gcp3,gcp4,formf
      complex *16 dgcp1,dgcp2,dgcp3,dgcp4,imag,t127     
      complex*16 aa0,aa1,aa2,aa30,aa31,aa32,am4,ann4
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      t127=zs(n1,n2)+zs(n1,n7)+zs(n2,n7)
      s127=s(n1,n2)+s(n1,n7)+s(n2,n7)
      formf= 1d0/(1d0+sww/formlam**2)**2
*
      aa0=2*t127*a(n1,n7)*a(n7,n2)*zs(n3,n4)*zs(n5,n6)
      aa1=a(n3,n6)*b(n4,n5)*(
     &  ann4(n1,n5,n2,n1)+ann4(n1,n5,n7,n1)+
     &  ann4(n1,n6,n2,n1)+ann4(n1,n6,n7,n1))
      aa2=a(n1,n6)*(a(n1,n2)*b(n2,n5)+a(n1,n7)*b(n7,n5))
     &  *(a(n3,n5)*b(n5,n4)+a(n3,n6)*b(n6,n4))
      aa30=a(n6,n3)*b(n3,n5)+a(n6,n4)*b(n4,n5)
      aa31=a(n1,n3)*am4(n1,n2,n7,n4)
      aa32=am4(n1,n5,n6,n4)*
     & ( ann4(n3,n5,n2,n1)+ann4(n3,n5,n7,n1)
     & + ann4(n3,n6,n2,n1)+ann4(n3,n6,n7,n1))
      gcp1=2d0+(g1zcp+kapzcp+lamzcp*s127/mwcp**2)*formf
      gcp2=2d0+(g1zcp+kapzcp+lamzcp)*formf
      gcp3=-2*(g1zcp*formf+1d0)
      gcp4=lamzcp/mwcp**2*formf
      dgcp1=dcmplx(gcp1)
      dgcp2=dcmplx(gcp2)
      dgcp3=dcmplx(gcp3)
      dgcp4=dcmplx(gcp4)
      Atree7Bwzan=-imag/aa0*
     & (dgcp1*aa1+dgcp2*aa2+aa30*(dgcp3*aa31+dgcp4*aa32))
      return
      end 


      FUNCTION FAtree7Bwzan(sww,N1,N2,N3,N4,N5,N6,n7)
      implicit none 
      real*8 sww
      integer n1,n2,n3,n4,n5,n6,n7
      complex*16 fAtree7Bwzan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*..................................                   )anomcpz
      real*8 g1zcp,kapzcp,lamzcp,mwcp,formlam
      common/anomcpz/g1zcp,kapzcp,lamzcp,mwcp,formlam
* .............  end of common block .............
      real*8 s127,gcp1,gcp2,gcp3,gcp4,formf
      complex *16 dgcp1,dgcp2,dgcp3,dgcp4,imag,t127     
      complex*16 aa0,aa1,aa2,aa30,aa31,aa32,bm4,bnn4
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      formf= 1d0/(1d0+sww/formlam**2)**2
      t127=zs(n1,n2)+zs(n1,n7)+zs(n2,n7)
      s127=s(n1,n2)+s(n1,n7)+s(n2,n7)
      aa0=2*t127*b(n1,n7)*b(n7,n2)*zs(n3,n4)*zs(n5,n6)
      aa1=b(n3,n6)*a(n4,n5)*(
     &  bnn4(n1,n5,n2,n1)+bnn4(n1,n5,n7,n1)+
     &  bnn4(n1,n6,n2,n1)+bnn4(n1,n6,n7,n1))
      aa2=b(n1,n6)*(b(n1,n2)*a(n2,n5)+b(n1,n7)*a(n7,n5))
     &  *(b(n3,n5)*a(n5,n4)+b(n3,n6)*a(n6,n4))
      aa30=b(n6,n3)*a(n3,n5)+b(n6,n4)*a(n4,n5)
      aa31=b(n1,n3)*bm4(n1,n2,n7,n4)
      aa32=bm4(n1,n5,n6,n4)*
     & ( bnn4(n3,n5,n2,n1)+bnn4(n3,n5,n7,n1)
     & + bnn4(n3,n6,n2,n1)+bnn4(n3,n6,n7,n1))
      gcp1=2d0+(g1zcp+kapzcp+lamzcp*s127/mwcp**2)*formf
      gcp2=2d0+(g1zcp+kapzcp+lamzcp)*formf
      gcp3=-2d0*(g1zcp*formf+1d0)
      gcp4=lamzcp/mwcp**2*formf
      dgcp1=dcmplx(gcp1)
      dgcp2=dcmplx(gcp2)
      dgcp3=dcmplx(gcp3)
      dgcp4=dcmplx(gcp4)
      FAtree7Bwzan=-imag/aa0*
     & (dgcp1*aa1+dgcp2*aa2+aa30*(dgcp3*aa31+dgcp4*aa32))
      return
      end 


******  (16a2) Atree6Bwzan as defined in eq. (7) anom. note *****
      FUNCTION Atree6Bwzan(sww,N1,N2,N3,N4,N5,N6)
      implicit none 
      real*8 sww
      integer n1,n2,n3,n4,n5,n6
      complex*16 Atree6Bwzan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*..................................                   )anomcpz
      real*8 g1zcp,kapzcp,lamzcp,mwcp,formlam
      common/anomcpz/g1zcp,kapzcp,lamzcp,mwcp,formlam
* .............  end of common block .............
      complex*16 imag,aa0,aa1,aa2,aa30,aa31,aa32
      real*8 gcp1,gcp2,gcp3,gcp4,formf     
      complex*16 dgcp1,dgcp2,dgcp3,dgcp4     
* //////////////////////////////////
*
      formf= 1d0/(1d0+sww/formlam**2)**2
      imag = (0.0 ,1.0 )
      aa0 =-imag/(2*zs(n1,n2)*zs(n3,n4)*zs(n5,n6))
      aa1=A(n3,n6)*B(n4,n5)*(A(n1,n5)*B(n5,n2)+A(n1,n6)*B(n6,n2))
      aa2=A(n1,n6)*B(n2,n5)*(A(n3,n1)*B(n1,n4)+A(n3,n2)*b(n2,n4))
      aa30=A(n6,n3)*B(n3,n5)+A(n6,n4)*b(n4,n5)
      aa31=A(n1,n3)*B(n2,n4)
      aa32=(A(n3,n5)*B(n5,n2)+A(n3,n6)*B(n6,n2))
     &    *(A(n1,n5)*B(n5,n4)+A(n1,n6)*b(n6,n4))
      gcp1=2d0+(g1zcp+kapzcp + lamzcp*s(n1,n2)/mwcp**2)*formf
      gcp2=2d0+(g1zcp+kapzcp + lamzcp)*formf
      gcp3=2*(g1zcp*formf+1d0)
      gcp4=lamzcp/mwcp**2*formf
      dgcp1=dcmplx(gcp1,0d0)
      dgcp2=dcmplx(gcp2,0d0)
      dgcp3=dcmplx(gcp3,0d0)
      dgcp4=dcmplx(gcp4,0d0)
      atree6Bwzan=aa0*(dgcp1*aa1+dgcp2*aa2+
     & aa30*(dgcp3*aa31+dgcp4*aa32))/imag
      return
      end 

*
       FUNCTION AloopBwzan(sww,i1,i2,i3,i4,i5,i6)
       implicit none
       real*8 sww
       integer  i1,i2,i3,i4,i5,i6                                   
       complex*16 Aloopbwzan, Vub                                         
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      real*8 s12,musq,pi
      parameter(          pi = 3.14159265358979323846d0 )                                                                           
*                                                                 1) nparton
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart
*
      complex*16  atree6bwzan
       musq=muscloop**2
       s12 = S(i1,i2)                                                             
                                                                                  
*            !! Add DR->CDR shifts later in formulas below                         
                                                                                  
       Vub = - 0.5*(log(musq/s12)**2 - pi**2) - 1.5*log(musq/s12)      
     &        - 3.5d0 - imag*abspart*pi*(log(musq/s12) + 1.5)            

       Aloopbwzan = atree6bwzan(sww,i1,i2,i3,i4,i5,i6)*Vub
 
      end
c
c
c End of wzanDKS.f
c
c
c
c
c Begin of wwpackage.f
c
c
c
c In order to avoid conflicts with the WZ part (implemented previously)
c we have renamed:
c    TMATR6NEW --> TMATR6NEW_WW
c    TMATR7NEW --> TMATR7NEW_WW
c    sigtyp    --> sigtyp_WW
c 
c
          SUBROUTINE TMATR7NEW_WW(xi,ximax,sigmat,sigmatb)
* sigmat(1 or 2,j): uubar->ww, sigmat(2,j):ddbar->ww 
* sigmatb(1 or 2,j):ubaru->ww, sigmatb(2,j):dbard->ww
* sigmat(iso,j):j=1,2,3,4: loop, born, collp,collm
          implicit none
          real*8 xi,ximax,sigmat(6,4),sigmatb(6,4)
          integer npart
*                                                   7)pmom
      REAL*8 plaball,pcmall,pwwg
      COMMON/PMOM/ PLABall(4,7,4),PCMall(4,7,4), pwwg(5,3,4)             
*                                                    10)qinne

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
cSF constww is computed in setdkspar, hence this extra common block
      real*8 constww
      common/DKSconstww/constww
*
      integer i,ii
cSF Previously, constww and alfa were defined here (alfa also set)
      real*8 norm(6)
*
       npart=7
cSF constww is computed in setdkspar
c$$$       constww=nc*(6*alfa/swsq)**2*mw**4/(32.d0)
      do i=1,6
      if(i.lt.3)then
      norm(i)=constww/nc**2
      else
      norm(i)=constww/(nc*(nc**2-1))
      endif
      enddo
cSF Eliminated for consistency with what is done for WZ
c$$$      CALL CHANGWALL(NPART,plaball,pcmall)
* -------------------------------------------------------------
* full
      ii=1
      if(xi.lt.ximax)then
      call sigtyp_WW(ii,npart,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
      endif
* soft
      ii=2
      call sigtyp_WW(ii,npart-1,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
* collp and collm
      do ii=3,4
      if(xi.lt.ximax)then
      call sigtyp_WW(ii,npart-1,norm,pcmall,plaball,pwwg,
     & sigmat(1,ii), sigmatb(1,ii))
      endif
      enddo
      return
      end    

*
*
          SUBROUTINE TMATR6NEW_WW(xi,ximax,sigmat,sigmatb)
* sigmat(1 or 2,j): uubar->ww, sigmat(2,j):ddbar->ww 
* sigmatb(1 or 2,j):ubaru->ww, sigmatb(2,j):dbard->ww
* sigmat(iso,j):j=1,2,3,4: loop, born, collp,collm
          implicit none
          real*8 xi,ximax,sigmat(2,4),sigmatb(2,4)
          integer npart
*                                                   7)pmom
      REAL*8 plaball,pcmall,pwwg
      COMMON/PMOM/ PLABall(4,7,4),PCMall(4,7,4), pwwg(5,3,4)             
*                                                    10)qinne

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)

*
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*                                                          ??) subcont
      logical bornonly,qqonly,gqonly
      common/subcont/ bornonly,qqonly,gqonly
*
      real*8 ubu4l,dbd4l,sww,s12
      complex*16 fcl(2,2),fcr(2,2)
      integer ii
      real*8 normqq
cSF constww is computed in setdkspar, hence this extra common block
cSF constww and alfa were previously defined here
      real*8 constww
      common/DKSconstww/constww

*
       npart=6
cSF Remove call as done in WZ
c$$$      CALL CHANGWALL(NPART+1,plaball,pcmall)
*
           ii=1
cSF Remove call as done in WZ
c$$$           call qinn(npart,pcmall(1,1,ii+1),a,b,zs,s)                      
c$$$           s12 =k1pk2sq(plaball(1,1,ii+1),plaball(1,2,ii+1),1)
c$$$ special s12 passing from dks wrapper sh = (p_had1 + p_had2)^2 
           s12 = sigmat(1,1)
           sww=s12
* ------------------------------------------------------
* couplings introduced by lance (eq. 17,19)
* ----------------------------------------
      fcl(1,1)=dcmplx(2*qu*swsq)
      fcl(1,2)=dcmplx(-2*qd*swsq)
      fcr(1,1)=dcmplx(2*qu*swsq)
      fcr(1,2)=dcmplx(-2*qd*swsq)
      fcl(2,1)=dcmplx(sww*(1-2*qu*swsq)/(sww-mz**2))
      fcl(2,2)=dcmplx(sww*(1+2*qd*swsq)/(sww-mz**2))
      fcr(2,1)=dcmplx(-2*qu*swsq*sww/(sww-mz**2))
      fcr(2,2)=dcmplx(2*qd*swsq*sww/(sww-mz**2))

     
* if uubar intitial state iso=1 if ddbar initital state iso=2
* in sigmat
* for ubaru and dbard we have the same but in sigmatb
* if loop ii=1
* if born or soft ii=2
      normqq=constww/(s12*nc**2)
      do ii=1,2
        if((bornonly.eqv..true.).and.(ii.eq.1))goto 1
          sigmat(1,ii)=ubu4l(npart,ii,sww,
     &        fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq
          sigmat(2,ii)=dbd4l(npart,ii,sww,
     &       fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq
          sigmatb(1,ii)=ubu4l(npart,ii,sww,
     &       fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq
          sigmatb(2,ii)=dbd4l(npart,ii,sww,
     &    fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq
 1    continue
      enddo
      if(bornonly.eqv..true.)goto 2
      if(xi.lt.ximax)then
* ii=3 
           ii=3
cSF Remove call as done in WZ
c$$$           call qinn(npart,pcmall(1,1,ii),a,b,zs,s)                       
c$$$           s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c$$$           sww=s12
*
      fcr(2,2)=dcmplx(2*qd*swsq*sww/(sww-mz**2))
      fcl(1,1)=dcmplx(2*qu*swsq)
      fcl(1,2)=dcmplx(-2*qd*swsq)
      fcr(1,1)=dcmplx(2*qu*swsq)
      fcr(1,2)=dcmplx(-2*qd*swsq)
      fcl(2,1)=dcmplx(sww*(1-2*qu*swsq)/(sww-mz**2))
      fcl(2,2)=dcmplx(sww*(1+2*qd*swsq)/(sww-mz**2))
      fcr(2,1)=dcmplx(-2*qu*swsq*sww/(sww-mz**2))
      fcr(2,2)=dcmplx(2*qd*swsq*sww/(sww-mz**2))
      normqq=constww/(s12*nc**2)
      sigmat(1,ii)=ubu4l(npart,ii,sww,
     &     fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq
      sigmat(2,ii)=dbd4l(npart,ii,sww,
     &     fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq
      sigmatb(1,ii)=ubu4l(npart,ii,sww,
     &     fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq
      sigmatb(2,ii)=dbd4l(npart,ii,sww,
     &     fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq
      ii=4
cSF Remove call as done in WZ
c$$$      call qinn(npart,pcmall(1,1,ii),a,b,zs,s)                       
c$$$      s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c$$$      sww=s12
*
      fcl(1,1)=dcmplx(2*qu*swsq)
      fcl(1,2)=dcmplx(-2*qd*swsq)
      fcr(1,1)=dcmplx(2*qu*swsq)
      fcr(1,2)=dcmplx(-2*qd*swsq)
      fcl(2,1)=dcmplx(sww*(1-2*qu*swsq)/(sww-mz**2))
      fcl(2,2)=dcmplx(sww*(1+2*qd*swsq)/(sww-mz**2))
      fcr(2,1)=dcmplx(-2*qu*swsq*sww/(sww-mz**2))
      fcr(2,2)=dcmplx(2*qd*swsq*sww/(sww-mz**2))
      normqq=constww/(s12*nc**2)
      sigmat(1,ii)=ubu4l(npart,ii,sww,
     &     fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq
      sigmat(2,ii)=dbd4l(npart,ii,sww,
     &     fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq
      sigmatb(1,ii)=ubu4l(npart,ii,sww,
     &     fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq
      sigmatb(2,ii)=dbd4l(npart,ii,sww,
     &     fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq
      endif
*
 2    continue
      return
      end    
*
      subroutine sigtyp_WW(ii,np,norm,pcmall,plaball,pwwg,
     & sig, sigb)
      implicit none
      integer ii,np
      real*8 norm(6),pcmall(4,7,4),plaball(4,7,4),
     & pwwg(5,3,4)
*                                               2)initv
      real *8 pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
      common/initv/ pi,nc,tf,cf,swsq, s2w,
     & el,er,ul, ur,dl,dr, mz ,mw, gz,gw,mtop,
     & qu,qd
*                                                    10)qinne

      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*
      real*8 ubu4l,dbd4l,normqq(6),
     & sig(6),sigb(6),s12,sww
      complex *16 fcl(2,2),fcr(2,2)
      integer i
* -------------------------
cSF hack as done for WZ
      s12 = sig(1)
      if (ii.eq.1.and.np.eq.7) then
         sww = sig(2)
      else
         sww = s12
      endif
cSF end hack
c$$$The hack replaces what is commented out below
c$$$      s12 =k1pk2sq(plaball(1,1,ii),plaball(1,2,ii),1)
c$$$        if(ii.eq.1.and.np.eq.7)then
c$$$          sww =k1pk2sq(pwwg(1,1,ii),pwwg(1,2,ii),1)
c$$$        else
c$$$           sww=s12
c$$$        endif
c$$$       call qinn(np,pcmall(1,1,ii),a,b,zs,s)                       
c$$$
* -------------------------------
      do i=1,6
      normqq(i)=norm(i)/s12
      sig(i)=0
      sigb(i)=0
      enddo
*
      fcl(1,1)=dcmplx(2*qu*swsq)
      fcl(1,2)=dcmplx(-2*qd*swsq)
      fcr(1,1)=dcmplx(2*qu*swsq)
      fcr(1,2)=dcmplx(-2*qd*swsq)
      fcl(2,1)=dcmplx(sww*(1-2*qu*swsq)/(sww-mz**2))
      fcl(2,2)=dcmplx(sww*(1+2*qd*swsq)/(sww-mz**2))
      fcr(2,1)=dcmplx(-2*qu*swsq*sww/(sww-mz**2))
      fcr(2,2)=dcmplx(2*qd*swsq*sww/(sww-mz**2))

* for any ii
*uubar
      sig(1)=ubu4l(np,ii,sww,
     &     fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq(1)
*ddbar
      sig(2)=dbd4l(np,ii,sww,
     &     fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq(2)
*
*ubar-u
      sigb(1)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq(1)
*dbar-d
      sigb(2)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq(2)
*
      if(ii.eq.1)then
*g-ubar
      sig(3)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),2,7,3,4,5,6,1)*normqq(3)
*g-dbar
      sig(4)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),2,7,3,4,5,6,1)*normqq(4)
*g-u
      sig(5)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),7,2,3,4,5,6,1)*normqq(5)
*g-d
      sig(6)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),7,2,3,4,5,6,1)*normqq(6)
*u-g
      sigb(5)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),7,1,3,4,5,6,2)*normqq(3)
*d-g
      sigb(6)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),7,1,3,4,5,6,2)*normqq(4)
*ubar-g
      sigb(3)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),1,7,3,4,5,6,2)*normqq(5)
*dbar-g
      sigb(4)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),1,7,3,4,5,6,2)*normqq(6)
      endif
* in the collinear limit only quark initial states
      if(ii.eq.3)then
*g-ubar-->(ubar)u-ubar
      sig(3)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq(1)
      sig(4)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq(1)
*g-u-->(u)ubar-u
      sig(5)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq(1)
      sig(6)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq(1)

      endif 
*
      if(ii.eq.4)then
      sigb(5)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),2,1,3,4,5,6,7)*normqq(1)
      sigb(6)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),2,1,3,4,5,6,7)*normqq(1)
      sigb(3)=ubu4l(np,ii,sww,
     &    fcl(1,1),fcr(1,1),1,2,3,4,5,6,7)*normqq(1)
      sigb(4)=dbd4l(np,ii,sww,
     &    fcl(1,2),fcr(1,2),1,2,3,4,5,6,7)*normqq(1)
      endif
      return
      end    
*
      FUNCTION ubu4l(npart,ii,sww,
     &   fcl,fcr,i1,i2,i3,i4,i5,i6,i7)
            IMPLICIT NONE
            INTEGER npart,ii,i1,i2,i3,i4,i5,i6,i7
            REAL*8 sww,ubu4l
            COMPLEX*16 fcl(2),fcr(2)
* the POSITION of i1,i2,i3,i4,i5,i6,i7 in the function call
* corresponds to
* ubar(1)^-u(2)^+l(3)^-nubar(4)^+lpbar(5)^+nup(6)^-g(7)^+
* where all particles are outgoing and the upper label refers
* to helicity (see eq. (2.6) in DKS.
* the ACTUAL VALUE of i1,i2,i3,i4,i5,i6,i7 is the 
* the momentum label in the physical process
* if npart=6, the gluon momenta=0 and the 
* label i7 is dummy in a sense 
* from the primitive amplitude we calculate
* (see
*     m(i)_tree_n=|ampla
*     (n,ii,1,2,3,4,5,6,7)+cl(i)amplbwwan(n,ii,1,2,3,4,5,6,7)|^2+
*               |cr(i)*amplbwwan(n,ii2,1,3,4,5,6,7)|^2
*
      real*8 gluneghel
      complex*16 a1,a21,a22,a31,a32,z1u,z2u,ampla,amplbwwan
      complex*16 b1,b21,b22,b31,b32,zc1u,zc2u
*
      ubu4l=0
      gluneghel=0.d0
* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
* for up quarks
      a1=ampla(npart,ii,i1,i2,i3,i4,i5,i6,i7)
      a21=amplbwwan(sww,1,npart,ii,i1,i2,i3,i4,i5,i6,i7)
      a22=amplbwwan(sww,2,npart,ii,i1,i2,i3,i4,i5,i6,i7)
      a31=amplbwwan(sww,1,npart,ii,i2,i1,i3,i4,i5,i6,i7)
      a32=amplbwwan(sww,2,npart,ii,i2,i1,i3,i4,i5,i6,i7)
      b1=a1
      b21=a21
      b31=a31
      b22=a22
      b32=a32
        if(npart.eq.6.and.ii.eq.1)then
* loop amplitudes
          b1=ampla(npart,ii+1,i1,i2,i3,i4,i5,i6,i7)
          b21=amplbwwan(sww,1,npart,ii+1,i1,i2,i3,i4,i5,i6,i7)
          b31=amplbwwan(sww,1,npart,ii+1,i2,i1,i3,i4,i5,i6,i7)
          b22=amplbwwan(sww,2,npart,ii+1,i1,i2,i3,i4,i5,i6,i7)
          b32=amplbwwan(sww,2,npart,ii+1,i2,i1,i3,i4,i5,i6,i7)
        endif
      if(ii.eq.1.and.npart.eq.6)then
      endif
      zc1u=(b1)+fcl(1)*(b21)+fcl(2)*(b22)
      zc2u=fcr(1)*(b31)+fcr(2)*(b32)
      z1u=a1+fcl(1)*a21+fcl(2)*a22
      z2u=fcr(1)*a31+fcr(2)*a32
      ubu4l=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
      gluneghel=0
      if(npart.eq.7.and.ii.eq.1)then
* negative helicity gluons can also contribute
        a1=ampla(npart,ii,i2,i1,i5,i6,i3,i4,i7)
        a21=amplbwwan(sww,1,npart,ii,i2,i1,i5,i6,i3,i4,i7)
        a31=amplbwwan(sww,1,npart,ii,i1,i2,i5,i6,i3,i4,i7)
        a22=amplbwwan(sww,2,npart,ii,i2,i1,i5,i6,i3,i4,i7)
        a32=amplbwwan(sww,2,npart,ii,i1,i2,i5,i6,i3,i4,i7)
        b1=a1
        b21=a21
        b22=a22
        b31=a31
        b32=a32
        zc1u=(b1)+fcl(1)*(b21)+fcl(2)*(b22)
        zc2u=fcr(1)*(b31)+fcr(2)*(b32)
        z1u=a1+fcl(1)*a21+fcl(2)*a22
        z2u=fcr(1)*a31+fcr(2)*a32
        gluneghel=dreal(z1u*dconjg(zc1u))+dreal(z2u*dconjg(zc2u))
       endif
        ubu4l=ubu4l+gluneghel
      return 
      end
*
      FUNCTION dbd4l(npart,ii,sww,
     &   fcl,fcr,i1,i2,i3,i4,i5,i6,i7)
            IMPLICIT NONE
            INTEGER npart,ii,i1,i2,i3,i4,i5,i6,i7
            REAL*8 sww,dbd4l,ubu4l
            COMPLEX*16 fcl(2),fcr(2)
      dbd4l=ubu4l(
     & npart,ii,sww,fcl,fcr,i1,i2,i6,i5,i4,i3,i7)
       return
       end

*
      FUNCTION AmplBwwan(sww,ia,npart,ii,N1,N2,N3,N4,N5,N6,n7)
      implicit none
       real*8 sww  
      integer ia,npart,ii,n1,n2,n3,n4,n5,n6,n7
      complex*16 AmplBwwan,aloopbwwan,atreeb6wwan,
     & atree7bwwan,imag 
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
* .............  end of common block .............
*
      imag = (0.0 ,1.0 )
      if(ii.eq.1.and.npart.eq.6)then
       AmplBwwan=aloopbwwan(sww,ia,n1,n2,n3,n4,n5,n6)
      endif
      if(ii.eq.1.and.npart.eq.7)then
* calculate atreeb7
      AmplBwwan=atree7bwwan(sww,ia,n1,n2,n3,n4,n5,n6,n7)
      endif
      if(ii.ne.1)then
       amplbwwan=atreeb6wwan(sww,ia,n1,n2,n3,n4,n5,n6)
      endif
      return
      end 
* .............
       FUNCTION AloopBwwan(sww,ia,i1,i2,i3,i4,i5,i6)
       implicit none
       real*8 sww 
       integer  ia,i1,i2,i3,i4,i5,i6                                   
       complex*16 Aloopbwwan, Vub                                         
*                                                                 5)QUINNE
* innner products for given pcmall momenta
      COMPLEX*16 A,B,ZS
      real*8 s           
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
      real*8 s12,musq,pi
      parameter(          pi = 3.14159265358979323846d0 )                                                                           
*                                                                 1) nparton
      complex* 16 imag
      real*8 zero,QES2,muscloop
      integer abspart
      common/nparton/imag,zero,QES2,muscloop,abspart
*
      complex*16  atreeb6wwan
       musq=muscloop**2
       s12 = S(i1,i2)                                                             
                                                                                  
*            !! Add DR->CDR shifts later in formulas below                         
                                                                                  
       Vub = - 0.5*(log(musq/s12)**2 - pi**2) - 1.5*log(musq/s12)      
     &        - 3.5d0 - imag*abspart*pi*(log(musq/s12) + 1.5)            
       Aloopbwwan = atreeb6wwan(sww,ia,i1,i2,i3,i4,i5,i6)*Vub
        end


**
**------------------
      FUNCTION Atree7Bwwan(sww,ia,N1,N2,N3,N4,N5,N6,n7)
** -------------
      implicit none 
      integer ia,n1,n2,n3,n4,n5,n6,n7
      complex*16 Atree7Bwwan
      real*8 sww
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*..................................                   )anom coup
      real*8 g1vcp,kapvcp,lamvcp,mwcp,formlam
      common/anomcp/g1vcp(2),kapvcp(2),lamvcp(2),mwcp,formlam
* .............  end of common block .............
      complex *16 gcp1,gcp2,gcp3     
* .............  end of common block .............
      complex*16 imag,aa0,aa1,aa20,aa21,aa22,t127,am3,ann4
* //////////////////////////////////
*
      imag = (0.0 ,1.0 )
      t127=zs(n1,n2)+zs(n1,n7)+zs(n2,n7)
      aa0=2*t127*a(n1,n7)*a(n7,n2)*zs(n3,n4)*zs(n5,n6)
      aa1=a(n1,n3)*(am3(n1,n2,n4)+am3(n1,n7,n4))*
     &             (am3(n6,n3,n5)+am3(n6,n4,n5))
     & -a(n1,n6)*(am3(n1,n2,n5)+am3(n1,n7,n5))*
     &             (am3(n3,n5,n4)+am3(n3,n6,n4))
      aa20=
     & ( ann4(n1,n3,n2,n1)+ann4(n1,n3,n7,n1)
     &  + ann4(n1,n4,n2,n1)+ann4(n1,n4,n7,n1))
      aa21=a(n3,n6)*b(n4,n5)
      aa22=
     & ( am3(n3,n4,n5)+am3(n3,n6,n5))*
     & ( am3(n6,n3,n4)+am3(n6,n5,n4))
      gcp1=2+(g1vcp(ia)+kapvcp(ia)+lamvcp(ia))/(1+sww/formlam**2)**2
      gcp2=2+2*g1vcp(ia)/(1+sww/formlam**2)**2
      gcp3=lamvcp(ia)/(mwcp**2*(1+sww/formlam**2)**2)
      Atree7bwwan=imag/aa0*
     & (gcp1*aa1+aa20*(gcp2*aa21+gcp3*aa22))
      return
      end 


      FUNCTION Atreeb6wwan(sww,ia,N1,N2,N3,N4,N5,N6)
      implicit none
      real*8 sww 
      integer ia,n1,n2,n3,n4,n5,n6
      complex*16 Atreeb6wwan
*............    common blocks ....................
*                                                    5) QUINNE               
      COMPLEX*16 A,B,ZS
      real*8 s
      COMMON/QINNE/A(7,7),B(7,7),ZS(7,7),S(7,7)
*..................................                  .. )anom coup
      real*8 g1vcp,kapvcp,lamvcp,mwcp,formlam
      common/anomcp/g1vcp(2),kapvcp(2),lamvcp(2),mwcp,formlam
* .............  end of common block .............
      complex*16 imag,aa0,aa1,aa20,aa21,aa22
      complex *16 gcp1,gcp2,gcp3     
*
      imag = (0.0 ,1.0 )
      aa0 =imag/(2*zs(n1,n2)*zs(n3,n4)*zs(n5,n6))
      aa1=A(n1,n3)*B(n2,n4)*(A(n6,n1)*B(n1,n5)+A(n6,n2)*B(n2,n5))
     &   +A(n1,n6)*B(n2,n5)*(A(n3,n5)*B(n5,n4)+A(n3,n6)*b(n6,n4))
      aa20=A(n1,n3)*B(n3,n2)+A(n1,n4)*b(n4,n2)
      aa21=A(n3,n6)*B(n4,n5)
      aa22=(A(n3,n1)*B(n1,n5)+A(n3,n2)*B(n2,n5))
     &    *(A(n6,n1)*B(n1,n4)+A(n6,n2)*b(n2,n4))
      gcp1=2+(g1vcp(ia)+kapvcp(ia) + lamvcp(ia))/(1+sww/formlam**2)**2
      gcp2=2+2*g1vcp(ia)/(1+sww/formlam**2)**2
      gcp3=lamvcp(ia)/(mwcp**2*(1+sww/formlam**2)**2)
      atreeb6wwan=aa0*(gcp1*aa1+aa20*(gcp2*aa21+gcp3*aa22))/imag
      return
      end 
c
c
c End of wwpackage.f
c
c

