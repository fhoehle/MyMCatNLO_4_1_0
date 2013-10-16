      PROGRAM MCATNLO_HTMAIN_DR_PP
      implicit none
      include 'stpcblks.h'
      real * 8 value(20),xmass(-5:21),xmomshifts(4),vickm(1:6,1:6),
     #  xlep1mass(2),xlep2mass(2),alsfht(3),besfht(3),alclht(3),
     #  beclht(3),alazi(3),beazi(3)
      real * 8 pi,xalsfht,xbesfht,xalclht,xbeclht,xalazi,xbeazi,
     #  al_spcfun,be_spcfun,wgtaev,wgtbev,xicut,deltai,deltao,xicutss,
     #  deltas,deltac,etacut,xmw,gaw,xmt,twidth,tga1,tga1mn,tga1pl,
     #  tbw1fmpl,tbw1fmmn,tbw1delf,ym1low2,ym1upp2,ga2,ga2mn,ga2pl,
     #  bw2fmpl,bw2fmmn,bw2delf,xm2low2,xm2upp2,xm012,ga1,bw1delf,
     #  bw1fmmn,xm1low2,xm1upp2,brrtop1,brrv2msb,xmone,tmas,xpdflam4,
     #  xpdflam5,ecm,xren,xfh,xrenmc,xfhmc,gammay1,ym1low,ym1upp,
     #  gammax1,xm1low,xm1upp,gammax2,xm2low,xm2upp,tmp,ac1,ac2,rohlim,
     #  xtotal,ytotal,dtot,xares,yares,xbres,ybres,av3a,d3a,av3nega,
     #  d3nega,av3b,d3b,av3negb,d3negb,ctime,avtot,evfrac,evprcfrac,
     #  dummy,tbw1mdpl,tbw1mdmn,bw2mdpl,bw2mdmn,bw1mdpl,bw1mdmn,
     #  bw1fmpl,scptveto,xbrrtoplep,xbrrtophad,xbrrwlep,xbrrwhad,
     #  hwidth,wgtmax,tanbeta,Ainput,Binput,xsecpp,xerrpp
      integer ih1,ih2,ndns1,ndns2,iinput,iwgtnorm,iprespl,ifxdaem,
     #  ichkmom,ichkpid,isubttype,nsamp,idec,jwidth,iwidth,il1hw,il2hw,
     #  izero,ione,inloscale,imcscale,ia1ora2,loproc,maproc,ichmin,
     #  ichmax,iwrong,iwrong1,neventsuw,nqeventsuw,ifailuw,ncntuws,
     #  nqcntuws,nmaxuw,nqmaxuw,ifuntype,ionshell,ideconsh,ifk88seed,
     #  ifk88ih,ifk88ndns,ipdfih,ipdfgroup,ipdfndns,mode,nlf,lo,
     #  iverbose,ibswrite,ifk88istrl,iprdct0hw,i,j,itmpih,idpdfset,
     #  itmpndns,maxevt,iseed0,iseed,maxtrials,iproc,it1,it2,iseld,
     #  ncl3,ndim,nwild,mx_of_evta,mx_of_evtb,itd1,itd2,ibscall,
     #  ntotal,ndiff,nevts,ntrls,iunita,iunitb,ioutput,itot,ii,
     #  iunit,ittintf,inonbtop,ievffmt,i2hdmtype
      character * 2 scheme
      character * 4 part1,part2
      character * 20 parm(20),gname
      character * 80 fname,fnamea,fnameb
      character * 80 fname1,fnamev
      character * 80 pref,prefn,prefev,prefnev
      character * 70 strin,strout,lhapdf
      logical evgen
      external sig5azw_ht,sig5bzw_ht
      parameter (pi=3.14159265358979312D0)
      parameter (xmone=-1.d0)
      parameter (izero=0)
      parameter (ione=1)
c
c common /strfun0/ is only in strfunht:
c ndns = pdf type
c ih1,ih2 = beam type (0=(p+n)/2, 1=p, -1=pbar, 2=n, -2=nbar)
      common/strfun0/ih1,ih2,ndns1,ndns2
c quark and gluon masses, used by Herwig. PDF labeling convention
      common/parmass/xmass
c CKM matrix elements entered by the user
      common/cvickm/vickm
c alsfht and besfht are the parameters entering gfunsoftht
      common/cgfunsfht/alsfht,besfht
c alclht and beclht are the parameters entering gfuncoll
      common/cgfunclht/alclht,beclht
c alazi and beazi are the parameters entering gfunazi
      common/cgfunazi/alazi,beazi
c al_spcfun, be_spcfun are the parameters entering spcdamp
      common/cspcpar/al_spcfun,be_spcfun
c scptveto is used to compute mu_F and/or mu_R when scale factors are negative
      common/cscptveto/scptveto
c iwgtnorm=0 for weight=+1/-1, iwgtnorm=1 otherwise
      common/ciwgtnorm/iwgtnorm
c number of events generated
      common/cmaxevt/maxevt
c ievffmt=0 for MC@NLO event file format, ievffmt=1 for LHEF format
      common/cievffmt/ievffmt
c wgtaev and wgtbev are the norms of weights for H and S events respectively
      common/cwgtev/wgtaev,wgtbev
c iprespl=0 ==> preserves rapidity
c iprespl=1 ==> preserves longitudinal momentum
      common/ciprespl/iprespl
c ifxdaem=0 ==> uses running alpha_EM(M^2)
c ifxdaem=1 ==> uses alpha_EM=1/137.0359895
      common/cifxdaem/ifxdaem
c ichkmom=0 --> enables checks on kinematics
      common/cichkmom/ichkmom
c ichkpid=0 --> enables checks on parton identities
      common/cichkpid/ichkpid
c----------------------------------------------------------
c Variables that control the integrations
c
      common/cisubttype/isubttype
      common/parsub/xicut,deltai,deltao
      common/xisave/xicutss
      common/pmerge/deltas,deltac
      common/samp/nsamp
c etacut is the maximum allowed for [2*kt(gluon)/sqrt(shat)]^2
      common/cetacut/etacut
c----------------------------------------------------------
c Decay variables
c Decay of the top and W: idec=0    -->   top decays
c                         idec=1    -->   top doesn't decay
      common/cidec/idec
c Mass ranges: jwidth=0    -->   top on shell
c              jwidth=1    -->   top off shell
      common/cjwidth/jwidth
c Mass range of W from top: iwidth=0    -->   W on shell
c                           iwidth=1    -->   W off shell
      common/ciwidth/iwidth
c Type of W decays; il1hw is entered with the following conventions,
c which control the decay of the W's
c  IL=0     ==> all W decays (quark+leptons)
c  IL=1,2,3 ==> W -> e\nu_e, mu\nu_mu, tau\nu_tau
c  IL=4     ==> W -> e\nu_e + mu\nu_mu
c  IL=5     ==> W -> all quarks
c  IL=6     ==> W -> e\nu_e + mu\nu_mu + all quarks (ie all decays except tau)
c  IL=7     ==> the top does not decay
c il1hw is relevant to the W from the top, il2hw is dummy here and set equal
c to a negative number to force the program to crash if used
      common/cilhw/il1hw,il2hw
c W mass and width (W mass squared is in cmass)
      common/cwparam/xmw,gaw
c top mass and width; top mass and its square are also stored in cmass 
      common/ctparam/xmt,twidth
c Top mass range
      common/tbw1/tga1,tga1mn,tga1pl,tbw1fmpl,tbw1fmmn,tbw1delf,
     #            ym1low2,ym1upp2
c Hard W mass range
      common/bw2cmm/ga2,ga2mn,ga2pl,bw2fmpl,bw2fmmn,bw2delf,
     #              xm2low2,xm2upp2
c Decay W mass range
      common/cbw1/xm012,ga1,bw1delf,bw1fmmn,xm1low2,xm1upp2
c top branching ratios, for lepton and hadron decays; apart from testing 
c purposes, these should be about 0.111 and 0.333 respectively,
c ie for W->e nu_e and W->udbar+usbar+ubar
      common/xibrratios/xbrrtoplep,xbrrtophad
c W branching ratios, for lepton and hadron decays; these are the
c analogues of xbrrtoplep,xbrrtophad
      common/xwibrratios/xbrrwlep,xbrrwhad
c reweight factors when top decays, that include branching ratios
c brrv2msb is set to one, and may be used for Higgs decays if done here
      common/brratios/brrtop1,brrv2msb
c mass of particles from W decays: not necessarily leptons!
      common/clepmass/xlep1mass,xlep2mass
c Decay of the top: inonbtop=0    -->   t->Wb only
c                   inonbtop=1    -->   t->W+any down-type quark
      common/cinonbtop/inonbtop
c tHb coupling stuff
c 2HDM model: i2hdmtype=1  -->  type I
c             i2hdmtype=2  -->  type II
      common/ci2hdmtype/i2hdmtype
c parameters used in coupling
      common/ctHbcoupl/tanbeta,Ainput,Binput
c----------------------------------------------------------
c inloscale controls the reference scale in the NLO computation
      common/cinloscale/inloscale
c imcscale controls the reference scale in the MC subtraction terms
      common/cimcscale/imcscale
c----------------------------------------------------------
c The following refer to the computation of MC subtraction terms
c ia1ora2=1 -> full invariants, ia1ora2=2 -> simplified invariants
      common/cia1ora2/ia1ora2
c----------------------------------------------------------
c Subprocesses: 'gg', 'qq', 'qg', corresponding to jproc=jproc0=1,2,3
c In the integration routines, loproc<=jproc<=maproc
      common/cwchproc/loproc,maproc
c Production channels: ichmin<=ich<=ichmax, ich=1,2,3 => s-, t-, Ht-channel
      common/cichrange/ichmin,ichmax
c ittintf=0 -> all processes; ittintf=1 -> no processes with ttb interference
      common/cittintf/ittintf
c Number of failures in flavour determination
      common/ciwrong/iwrong,iwrong1
c Common blocks for statistics relevant to secondary unweighting
      common/c1iunwgt/neventsuw,nqeventsuw,ifailuw
      common/c2iunwgt/ncntuws,nqcntuws,nmaxuw,nqmaxuw
c Average shifts in momenta, due to quark and lepton masses
      common/cshifts/xmomshifts
c----------------------------------------------------------
c ifuntype=1 for sig5a, ifuntype=2 for sig5b
      common/cifuntype/ifuntype
c Flag to put partons on shell, according to Herwig list of masses
      common/cionshell/ionshell
c Flag to put top decay products on shell
      common/cideconsh/ideconsh
c----------------------------------------------------------
c Common blocks for general MC@NLO routines
c common block for internal rnd number generation, independent of bases
      common/cifk88seed/ifk88seed
c common block fk88ipdfs is filled by our interface to MLMPDF
      common/fk88ipdfs/ifk88ih,ifk88ndns
c common block w50511 and w50512 are filled by PDFLIB 
      common/w50511/ipdfih,ipdfgroup,ipdfndns,mode,nlf,lo,tmas
      common/w50512/xpdflam4,xpdflam5
C
C------------------------------------------------------------------------
C                             START                                     -
C------------------------------------------------------------------------
c iinput=1 ==> all inputs are given by the user
      iinput=0
c iverbose=1 ==> writes more on standard output
      iverbose=0
c ichkmom=0 ==> enables checks on kinematics
      ichkmom=1
c ichkpid=0 ==> enables checks on parton identities
      ichkpid=0
c if linked to PDFLIB, these quantities stay negative
      ifk88ih=-100
      ifk88ndns=-100
c forces the code to get Lambda value if not obtained from PDFLIB/MLMPDF
      xpdflam5=-1.d0
C Set system dependent parameters
      call sysdep
c----- vegas prints nothing
c      call nopr(0)
c Bases writes data file
      ibswrite=1
c-----
c Open the file collecting all the input parameter. This file is meant 
c to be converted in a command file in a subsequent run
      open(unit=11,file='Htplog',status=newver)
c
      write(*,*)' '
      write(*,*)
     # 'Enter prefix for name of BASES files'
      read (*,*) pref
      write(11,*) ''''//pref(1:ifk88istrl(pref))//'''',
     # '  ! prefix for BASES files'
      write(*,*)' '
      write(*,*)
     # 'Enter prefix for name of event files'
      read (*,*) prefev
      write(11,*) ''''//prefev(1:ifk88istrl(prefev))//'''',
     # '  ! prefix for event files'
c----------------------------------------------------------
c Parameters of the run
      write(*,*)' '
      write(*,*)
     # 'Enter Ecm(GeV),fren[NLO],ffact[NLO],fren[MC],ffact[MC]'
      write(*,*)' fren=mu_ren/mu0'
      write(*,*)' ffact=mu_fac/mu0'
      write(*,*)' mu_ren=renormalization scale'
      write(*,*)' mu_fac=factorization scale'
      write(*,*)' mu0=reference scale'
      read(*,*) ecm,xren,xfh,xrenmc,xfhmc
      write(11,'(5(1x,d10.4),1x,a)') ecm,xren,xfh,xrenmc,xfhmc
     #     ,'! Ecm, fren, ffact, frenmc, ffactmc'
      sh = ecm**2
c Will allow more flexibility in future versions
      xrenmc = xren
      xfhmc = xfh
      xren2 = xren**2
      xf2h1 = xfh**2
      xf2h2 = xfh**2
      xren2mc = xrenmc**2
      xf2h1mc = xfhmc**2
      xf2h2mc = xfhmc**2
c----------------------------------------------------------
c Select process; enter HERWIG code, and converts to NLO codes
      write(*,*)' '
      write(*,*)'Enter -(1)2000-IT for all channels'
      write(*,*)'      -(1)2010-IT for s-channel'
      write(*,*)'      -(1)2020-IT for t-channel'
      write(*,*)'      -(1)2030-IT for Ht-channel'
      write(*,*)'  with IT=0 for t+tbar production'
      write(*,*)'       IT=1 for tbar production only'
      write(*,*)'       IT=4 for t production only'
      read(*,*) iprdct0hw
      write(11,'(1x,i6,27x,a)') iprdct0hw,
     #  '! 0/1/4 -> t+tb/tb/t'
      iprdct0hw=mod(-iprdct0hw,10000)
      call getnloiproc(iprdct0hw)
      if(ichmin.ne.3)then
        write(*,*)'Use this code for Ht-channel only'
        stop
      endif
c
      write(*,*)' '
      write(*,*)'Enter top mass and width (GeV)'
      read(*,*)xm1,twidth
      write(11,'(2(1x,d10.4),12x,a)') xm1,twidth,'! M_top, Gamma_top'
      xm12 = xm1**2
      xmt = xm1
c
      write(*,*)' '
      write(*,*)'Enter Higgs mass and width (GeV)'
      read(*,*)xm2,hwidth
      write(11,'(2(1x,d10.4),12x,a)') xm2,hwidth,'! M_H, Gamma_H'
      xm22 = xm2**2
c
      write(*,*)' '
      write(*,*)'Enter W mass and width (GeV)'
      read(*,*)xmw,gaw
      write(11,'(2(1x,d10.4),12x,a)') xmw,gaw,'! M_W, Gamma_W'
      xmW2 = xmw**2
c Top decay parameter
      write(*,*)' '
      write(*,*)'Enter IL=0..6 for t->W(->d1_IL d2_IL) b'
      write(*,*)'      IL=7 for undecayed top'
      read(*,*) il1hw
      write(11,'(1x,i2,31x,a)') il1hw,
     #  '! 0..6 -> t dec, 7 -> t undec'
      write(*,*)' '
      il2hw=-1
      if(il1hw.eq.7)then
        idec=1
      elseif( il1hw.ge.0.and.il1hw.le.6 )then
        idec=0
      else
        write(*,*) 'Unknown option:',il1hw
        stop
      endif
      if(idec.eq.0)then
        write(*,*)' '
        write(*,*)'Enter GammaX, M_T(min), M_T(max) for top'
        write(*,*)
     #   '  If GammaX>0, the top mass is chosen in the range'
        write(*,*)'      M0-GammaX*width < M_T < M0+GammaX*width'
        write(*,*)'  and M_T(min), M_T(max) are ignored'
        write(*,*)
     #   '  If GammaX<0, the top mass is chosen in the range'
        write(*,*)'            M_T(min) < M_T < M_T(max)'
        write(*,*)
     #   '  If GammaX=0, the top mass is set equal to the pole mass'
        read(*,*)gammay1,ym1low,ym1upp
        write(11,'(3(1x,d10.4),1x,a)') gammay1,ym1low,ym1upp,
     #   '! GammaX, M_T(min), M_T(max)'
        if(gammay1.lt.0.and.ym1low.ge.ym1upp)then
          write(*,*)'Enter a non-zero range'
          stop
        endif
c W from the top decay
        write(*,*)' '
        write(*,*)'Enter GammaX, M_V1(min), M_V1(max) for W from top'
        write(*,*)
     #   '  If GammaX>0, the boson mass is chosen in the range'
        write(*,*)'      M0-GammaX*width < M_W < M0+GammaX*width'
        write(*,*)'  and M_V1(min), M_V1(max) are ignored'
        write(*,*)
     #   '  If GammaX<0, the boson mass is chosen in the range'
        write(*,*)'            M_V1(min) < M_W < M_V1(max)'
        write(*,*)
     #   '  If GammaX=0, the boson mass is set equal to the pole mass'
        read(*,*)gammax1,xm1low,xm1upp
        write(11,'(3(1x,d10.4),1x,a)') gammax1,xm1low,xm1upp,
     #   '! GammaX, M_V1(min), M_V1(max)'
        if(gammax1.lt.0.and.xm1low.ge.xm1upp)then
          write(*,*)'Enter a non-zero range'
          stop
        endif
c Enter here gammax2,xm2low,xm2upp for Higgs with finite width
        gammax2=0.d0
        xm2low=0.d0
        xm2upp=0.d0
c
        write(*,*)' '
        write(*,*) 'Enter 0 for t->Wb only'
        write(*,*) '      1 for t->W+any down-type quark'
        read(*,*) inonbtop
        write(11,'(1x,i2,31x,a)') inonbtop,'! 0=t->Wb, 1=t->W+any d'
        if(inonbtop.ne.0.and.inonbtop.ne.1) then
          write(*,*)'Option not implemented: inonbtop=',inonbtop
          stop
        endif
c
        write(*,*)' '
        write(*,*)'Enter top -> leptons branching ratio'
        read(*,*)xbrrtoplep
        write(11,'(1x,d10.4,23x,a)') xbrrtoplep,
     #    '! t -> leptons branching ratio'
        write(*,*)' '
        write(*,*)'Enter top -> hadrons branching ratio'
        read(*,*)xbrrtophad
        write(11,'(1x,d10.4,23x,a)') xbrrtophad,
     #    '! t -> hadrons branching ratio'
        write(*,*)' '
        write(*,*)'Enter W -> leptons branching ratio'
        read(*,*)xbrrwlep
        write(11,'(1x,d10.4,23x,a)') xbrrwlep,
     #    '! W -> leptons branching ratio'
        write(*,*)' '
        write(*,*)'Enter W -> hadrons branching ratio'
        read(*,*),xbrrwhad
        write(11,'(1x,d10.4,23x,a)') xbrrwhad,
     #    '! W -> hadrons branching ratio'
c Redefine width or branching ratios if need be
        call reset_twdbr2(xmt,twidth,xmw,gaw,
     #                    xbrrtoplep,xbrrtophad,
     #                    xbrrwlep,xbrrwhad,
     #                    gammax1,gammax2)
      else
        if(twidth.lt.0.d0.or.gaw.lt.0.d0)then
          write(*,*)'Option not implemented if IL1CODE=7'
          write(*,*)' Enter positive top and W widths'
          stop
        endif
        twidth=0.d0
        gaw=0.d0
        inonbtop=-1
      endif
c CKM matrix elements
      do i=1,6
        do j=1,6
          vickm(i,j)=0.d0
        enddo
      enddo
      write(*,*)' '
      write(*,*)'Enter |V_ud|, |V_us|, |V_ub|'
      write(*,*)' all equal to zero to use PDG values'
      read(*,*)vickm(1,2),vickm(1,3),vickm(1,5)
      write(11,'(3(1x,d10.4),1x,a)')vickm(1,2),vickm(1,3),vickm(1,5),
     #      '! |V_ud|,|V_us|,|V_ub|'
      write(*,*)'Enter |V_cd|, |V_cs|, |V_cb|'
      read(*,*)vickm(4,2),vickm(4,3),vickm(4,5)
      write(11,'(3(1x,d10.4),1x,a)')vickm(4,2),vickm(4,3),vickm(4,5),
     #      '! |V_cd|,|V_cs|,|V_cb|'
      write(*,*)'Enter |V_td|, |V_ts|, |V_tb|'
      read(*,*)vickm(6,2),vickm(6,3),vickm(6,5)
      write(11,'(3(1x,d10.4),1x,a)')vickm(6,2),vickm(6,3),vickm(6,5),
     #      '! |V_td|,|V_ts|,|V_tb|'
c Quark and gluon masses (must be consistent with Herwig)
      do i=-5,21
        xmass(i)=0.d0
      enddo
      write(*,*)' '
      write(*,*)'Enter d, u, s, c, b, glu (Herwig) masses'
      read(*,*)xmass(1),xmass(2),xmass(3),xmass(4),xmass(5),xmass(21)
      write(11,'(6(1x,d10.4),1x,a)') xmass(1),xmass(2),xmass(3),
     #  xmass(4),xmass(5),xmass(21),'! quark and gluon masses'
      do i=-5,-1
        xmass(i)=xmass(-i)
      enddo
c 2HDM type and parameters relevant to tHb coupling
      write(*,*)' '
      write(*,*)'Enter 1 for a type I 2HDM'
      write(*,*)'      2 for a type II 2HDM'
      read(*,*) i2hdmtype
      write(11,'(1x,i2,31x,a)') i2hdmtype,'! 1=typeI, 2=typeII'
      if(i2hdmtype.ne.1.and.i2hdmtype.ne.2) then
        write(*,*)'Option not implemented: i2hdmtype=',i2hdmtype
        stop
      endif
c
      write(*,*)' '
      write(*,*)'Enter tan(beta), A and B'
      write(*,*)'   tan(beta) will be ignored in a type I 2HDM'
      write(*,*)'   A and B will be ignored in a type II 2HDM'
      read(*,*)tanbeta,Ainput,Binput
      write(11,'(3(1x,d10.4),1x,a)') tanbeta,Ainput,Binput,
     #   '! tan(beta), A, B'
c PDFs
      write(*,*)' '
      write(*,*)
     #  'Enter beam type for beam1 and beam2 (p, pbar, n, nbar):'
      read(*,*) part1,part2
      write(11,'(1x,a,2x,a,19x,a)') ''''//part1//'''',
     #  ''''//part2//'''','! hadron types'
      strin=part1
      call fk88low_to_upp(strin,strout)
      part1=strout
      strin=part2
      call fk88low_to_upp(strin,strout)
      part2=strout
      if( (part1.ne.'P   ' .and. part1.ne.'PBAR' .and. 
     #     part1.ne.'N   ' .and. part1.ne.'NBAR') .or.
     #    (part2.ne.'P   ' .and. part2.ne.'PBAR' .and. 
     #     part2.ne.'N   ' .and. part2.ne.'NBAR') )then
        write(*,*)'This code only works for hadronic collisions'
        stop
      endif
      call getihpart(part1,itmpih)
      ih1=itmpih
      call getihpart(part2,itmpih)
      ih2=itmpih
c
      write(*,*)' '
      write(*,*)'Enter group name and id number for PDF set'
      read(*,*)gname,idpdfset
      write(11,'(1x,a,1x,i6,21x,a)') 
     # ''''//gname(1:ifk88istrl(gname))//'''',
     # idpdfset,'! PDF group and id number'
      strin=gname
      call fk88low_to_upp(strin,strout)
      if(strout.eq.'LHAPDF'.or.strout.eq.'LHAEXT')then
        lhapdf='FREEZE'
        if(strout.eq.'LHAEXT')lhapdf='EXTRAPOLATE'
        call setlhacblk(lhapdf)
        parm(1)='DEFAULT'
      else
        lhapdf='NOTLKD'
        parm(1)=gname
      endif
      value(1)=idpdfset
      call pdfset(parm,value)
      if(ipdfih.ne.1)then
        write(*,*)'PDFLIB could not understand the input'
        write(*,*)'Hadron type:',ipdfih
        stop
      endif
      if(ifk88ih.eq.-100.and.ifk88ndns.eq.-100)then
        if(lhapdf.eq.'NOTLKD')then
c the code is linked to PDFLIB; get the MLM pdf id number from
c ipdfih, ipdfgroup, and ipdfndns returned by PDFLIB in /w50511/
          call pdftomlm(ipdfih,ipdfgroup,ipdfndns,itmpih,itmpndns)
          ndns1=itmpndns
          ndns2=itmpndns
        elseif(lhapdf.eq.'FREEZE'.or.lhapdf.eq.'EXTRAPOLATE')then
c the code is linked to LHAPDF, which doesn't fill /w50511/
          call pdftomlm(ione,izero,idpdfset,itmpih,itmpndns)
          ndns1=itmpndns
          ndns2=itmpndns
        else
          write(*,*) 'Unknown lhapdf value: ',lhapdf
          stop
        endif
      else
c the code is linked to the interface to MLMPDF
        ndns1=ifk88ndns
        ndns2=ifk88ndns
      endif
c
      scheme='**'
      if(lhapdf.eq.'NOTLKD')then
        xlam=xpdflam5
      else
        call getlamlha(xlam,xpdflam5)
      endif
c
      if(xlam.gt.0) then
         write(*,*)' '
         write(*,*) 'Enter Lambda_QCD_5, < 0 for default'
         read(*,*) tmp
         write(11,'(1x,d10.4,23x,a)') tmp,'! Lambda_5, < 0 for default'
         if(tmp.gt.0) xlam=tmp
      else
         dowhile(xlam.le.0)
            write(*,*)' '
            write(*,*)'Enter Lambda_5_2loop'
            read(*,*) xlam
c            if (xlam.le.0) call prntsf
         enddo
         write(11,'(1x,d10.4,23x,a)') xlam,'! lambda'
      endif
      write(*,*) 'Lambda_5=',xlam,' GeV'
c
 22   if(scheme.ne.'DI'.and.scheme.ne.'MS') then
         write(*,*)' '
         write(*,'(1x,a)') 'Enter scheme: ''DI'' or ''MS'''
         read(*,*) scheme
         if(scheme.ne.'DI'.and.scheme.ne.'MS') then
c            call prntsf
            goto 22
         endif
         write(11,'(1x,a,29x,a)') ''''//scheme//'''','! scheme'
      endif
      write(*,*) 'Scheme=',scheme
      schhad1=scheme
      schhad2=scheme
c-----------------------------------------------------------------
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_soft [bg]'
        write(*,*)' Defaults are: alpha=1, beta=-0.1'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalsfht,xbesfht
        write(11,'(2(1x,d10.4),12x,a)') xalsfht,xbesfht,
     #    '! alpha, beta [soft;bg]'
      else
        xalsfht=1.d0
        xbesfht=-0.1d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_coll [bg]'
        write(*,*)' Defaults are: alpha=0, beta=0'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalclht,xbeclht
        write(11,'(2(1x,d10.4),12x,a)') xalclht,xbeclht,
     #    '! alpha, beta [coll;bg]'
      else
        xalclht=0.d0
        xbeclht=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_azi [bg]'
        write(*,*)' Defaults are: alpha=-1, beta=-0.1'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalazi,xbeazi
        write(11,'(2(1x,d10.4),12x,a)') xalazi,xbeazi,
     #    '! alpha, beta [azi;bg]'
      else
        xalazi=-1.d0
        xbeazi=-0.1d0
      endif
      alsfht(3)=xalsfht
      besfht(3)=xbesfht
      alclht(3)=xalclht
      beclht(3)=xbeclht
      alazi(3)=xalazi
      beazi(3)=xbeazi
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_soft [gg,qq]'
        write(*,*)' Defaults are: alpha=-1, beta=0'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalsfht,xbesfht
        write(11,'(2(1x,d10.4),12x,a)') xalsfht,xbesfht,
     #    '! alpha, beta [soft;gg,qq]'
      else
        xalsfht=-1.d0
        xbesfht=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_coll [gg,qq]'
        write(*,*)' Defaults are: alpha=0, beta=0'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalclht,xbeclht
        write(11,'(2(1x,d10.4),12x,a)') xalclht,xbeclht,
     #    '! alpha, beta [coll;gg,qq]'
      else
        xalclht=0.d0
        xbeclht=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #   'Enter alpha and beta for G_azi [gg,qq]'
        write(*,*)' Defaults are: alpha=-1, beta=-0.1'
        write(*,*)' Allowed ranges: alpha>=1, 0<|beta|<=1'
        read(*,*) xalazi,xbeazi
        write(11,'(2(1x,d10.4),12x,a)') xalazi,xbeazi,
     #    '! alpha, beta [azi;gg,qq]'
      else
        xalazi=-1.d0
        xbeazi=-0.1d0
      endif
      alsfht(1)=xalsfht
      besfht(1)=xbesfht
      alclht(1)=xalclht
      beclht(1)=xbeclht
      alazi(1)=xalazi
      beazi(1)=xbeazi
      alsfht(2)=xalsfht
      besfht(2)=xbesfht
      alclht(2)=xalclht
      beclht(2)=xbeclht
      alazi(2)=xalazi
      beazi(2)=xbeazi
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha and beta for the function SPC_damp'
        write(*,*)' Defaults are: alpha=1, beta=0.5'
        write(*,*)' Allowed ranges: alpha>=1, 0<beta<=1'
        read(*,*) al_spcfun,be_spcfun
        write(11,'(2(1x,d10.4),12x,a)') al_spcfun,be_spcfun,
     #    '! alpha, beta [spin corr]'
      else
        al_spcfun=1.d0
        be_spcfun=0.5d0
      endif
c-----------------------------------------------------------------
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     # 'Enter accuracies for grid setup and for integral evaluation'
        read(*,*)ac1,ac2
        write(11,'(2(2x,d10.4),10x,a)') ac1,ac2,'! ac1,ac2'
      else
        ac1=0.2d0
        ac2=0.05d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
     #    'For the computation of the MEs in the MC subtraction terms'
        write(*,*)'Enter 1 to use full 2->3 invariants'
        write(*,*)'      2 to use simplified invariants'
        write(*,*)' The default is 1'
        read(*,*) ia1ora2
        write(11,'(1x,i2,31x,a)') ia1ora2,
     #    '! 1 for full, 2 for simplified invariants'
      else
        ia1ora2=1
      endif
c 
      write(*,*)' '
      write(*,*) 'Enter mass scale related to pt-veto'
      write(*,*) ' It is used to compute mu_fac if ffact<0'
      write(*,*) ' and mu_ren if fren<0'
      read(*,*) scptveto
      write(11,'(1x,d10.4,23x,a)') scptveto,'! ptveto[mu_F,mu_R]'
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'For the computation of alpha_S in NLO terms'
        write(*,*)'Enter 1 to set mu_0^2=(M_t^2+pt_t^2+M_H^2+pt_H^2)/2'
        write(*,*)'      2 to set mu_0^2=M_t^2'
        write(*,*)'      3 to set mu_0^2=(M_t^2+M_H^2)/2'
        write(*,*)'      4 to set mu_0^2=[(M_t+M_H)/2]^2'
        write(*,*)' The default is 4'
        read(*,*) inloscale
        write(11,'(1(1x,i8),25x,a)') inloscale,
     #    '! 1->mu_0=mt+<pt>, 2->mu_0=mt'
      else
        inloscale=4
      endif
      if(xfh.lt.0.and.xren.gt.0)then
        inloscale=inloscale+10
      elseif(xfh.gt.0.and.xren.lt.0)then
        inloscale=inloscale+20
      elseif(xfh.lt.0.and.xren.lt.0)then
        inloscale=inloscale+30
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'For the computation of alpha_S in MC terms'
        write(*,*)'Enter 1 to set mu_0^2=(M_t^2+pt_t^2+M_H^2+pt_H^2)/2'
        write(*,*)'      2 to set mu_0^2=M_t^2'
        write(*,*)'      3 to set mu_0^2=(M_t^2+M_H^2)/2'
        write(*,*)'      4 to set mu_0^2=[(M_t+M_H)/2]^2'
        write(*,*)' The default is 4'
        read(*,*) imcscale
        write(11,'(1(1x,i8),25x,a)') imcscale,
     #    '! 1->mu_0=mt+<pt>, 2->mu_0=mt'
      else
        imcscale=4
      endif
      if(xfhmc.lt.0.and.xrenmc.gt.0)then
        imcscale=imcscale+10
      elseif(xfhmc.gt.0.and.xrenmc.lt.0)then
        imcscale=imcscale+20
      elseif(xfhmc.lt.0.and.xrenmc.lt.0)then
        imcscale=imcscale+30
      endif
c
c Set constants [need Lambda_QCD in setpar_2hdm()]
      call setpar_2hdm()
      call setpar()
c
      write(*,*)' '
      write(*,*)'Enter the maximum number of events to generate;'
      write(*,*)'enter 0 to skip the event generation step'
      read(*,*)maxevt
      write(11,'(1(1x,i8),25x,a)') maxevt,'! number of events'
      evgen=.true.
      if(maxevt.eq.0)then
        evgen=.false.
        maxevt=100000
      endif
      write(*,*)' '
      write(*,*)'Enter 0 to have +1/-1 event weights'
      write(*,*)'      1 to normalize the weights, in such a way that'
      write(*,*)'        their sum is equal to the total rate'
      read(*,*)iwgtnorm
      write(11,'(1(1x,i8),25x,a)') iwgtnorm,
     #  '! 0 => wgt=+1/-1, 1 otherwise'
c iseed0 is the seed for the integration step
      iseed0=12345
      write(*,*)' '
      write(*,*)'Enter the seed for random numbers; it will be used'
      write(*,*)'to generate events. Enter 0 for default'
      read(*,*)iseed
      write(11,'(1(1x,i8),25x,a)') iseed,'! seed for rnd numbers'
      if(iseed.lt.0)then
        stop
      elseif(iseed.eq.0)then
        iseed=iseed0
      endif
c initialization of internal random number generation (same seed as
c used in bases/spring, this should not be a problem as they are different
c random number generators)
      ifk88seed=iseed
c Here, assume that the unweighting efficiency is larger than 10%
      maxtrials=10*maxevt
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to use standard subtraction'
        write(*,*)'      1 to use zeta subtraction'
        read(*,*)isubttype
        write(11,'(1(1x,i8),25x,a)') isubttype,
     #    '! 0=subt, 1=zeta subt'
      else
        isubttype=0
      endif
      if(iinput.eq.1)then
        if(isubttype.eq.0)then
          write(*,*)' '
          write(*,*)'enter xicut, deltaI, deltaO'
          read(*,*)xicut,deltai,deltao
          write(11,'(3(1x,d10.4),1x,a)') xicut,deltai,deltao,
     #      '! xicut,deltaI,deltaO'
        else
          write(*,*)' '
          write(*,*)'Enter zi ( [ 2*kt(gluon)/sqrt(shat) ]^2 < zi )'
          read(*,*) etacut
          write(11,'(1x,d10.4,23x,a)') etacut,'! zi'
          xicut = 1.d0
c Set deltaI=1 as prescribed by the zeta subtraction (see first MC@NLO paper)
          deltai = 1.d0
          deltao = 1.d0
        endif
        xicutss = xicut
        deltas = 0
        deltac = 0
      else
        if(isubttype.eq.0)then
          xicut = 1.d0
          deltai = 2.d0
          deltao = 2.d0
          xicutss = xicut
          deltas = 0
          deltac = 0
        else
          etacut = 1.d0
          deltai = 1.d0
          deltao = 1.d0
        endif
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to preserve rapidity'
        write(*,*)'      1 to preserve longitudinal momentum'
        read(*,*)iprespl
        write(11,'(1(1x,i8),25x,a)') iprespl,'! 0=y, 1=k_3 preserved'
      else
        iprespl=0
      endif
      write(*,*)' '
      write(*,*)'Enter 0 to use running alpha_EM'
      write(*,*)'      1 to use alpha_EM=1/137.0359895'
      read(*,*)ifxdaem
      write(11,'(1(1x,i8),25x,a)') ifxdaem,
     #  '! 0=running, 1=fixed alpha_EM'
      if(ifxdaem.ne.0.and.ifxdaem.ne.1)then
        write(*,*)'No such option for alpha_em'
        stop
      endif
      if(ifxdaem.eq.0)ze2=0.d0
c---------------------------------------------------------------
c Initialize parameters, such as labelling for parton processes
      call parsetpar()
c Select partonic subprocess
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*) 'Enter 1 for gg, 2 for qq, 3 for qg, 0 for all'
        read(*,*) iproc
        write(11,'(1x,i2,31x,a)') iproc,'! 1=gg, 2=qq, 3=qg, 0=all'
      else
        iproc=0
      endif
      if(iproc.lt.0.or.iproc.gt.3) then
        write(*,*)'Option not implemented: iproc=',iproc
        stop
      endif
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*) 
     #    'Enter 1 to exclude processes with ttbar interference'
        write(*,*) '      0 otherwise'
        read(*,*) ittintf
        write(11,'(1x,i2,31x,a)') ittintf,'! 1=no prc/w ttb, 0=all'
      else
        ittintf=0
      endif
      if(ittintf.ne.0.and.ittintf.ne.1) then
        write(*,*)'Option not implemented: ittintf=',ittintf
        stop
      endif
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to leave the partons massless'
        write(*,*)'      1(11) to put partons(FS) on mass shell'
        write(*,*)
     #   '      2(12) to put partons(FS) on mass shell, shat constant'
        read(*,*) ionshell
        write(11,'(1x,i1,32x,a)') ionshell,
     #   '! 0=massless, 1(11)=massive, 2(12)=massive, shat const'
      else
        ionshell=0
      endif
      if( ionshell.ne.0.and.ionshell.ne.1.and.ionshell.ne.2.and.
     #    ionshell.ne.11.and.ionshell.ne.12 ) then
        write(*,*) 'Error: enter 0, 1, 2, 11 or 12'
        stop
      endif
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to leave the decay products massless'
        write(*,*)'      2 to put them on their mass shell'
        write(*,*)'      12 to put only leptons on their mass shell'
        read(*,*) ideconsh
        write(11,'(1x,i1,32x,a)') ideconsh,
     #    '! 0=massless, 2(12)=massive decay(leptons) products'
      else
        ideconsh=12
      endif
      if(ideconsh.ne.0.and.ideconsh.ne.2.and.ideconsh.ne.12) then
        write(*,*) 'Error: enter 0, 2 or 12'
        stop
      endif
c
      write(*,*)' '
      write(*,*)'Enter 0 for MC@NLO-format event file'
      write(*,*)'      1 for LHEF-format event file'
      read(*,*)ievffmt
      write(11,'(1(1x,i8),25x,a)') ievffmt,
     #  '! 0 => MC@NLO format, 1 => LHEF'
      if(ievffmt.ne.0.and.ievffmt.ne.1)then
        write(*,*)'Unknown event file format',ievffmt
        stop
      endif
      write(*,*)' '
      write(*,*)'Enter number of iterations'
      write(*,*)'for grid setup and for integral evaluation;'
      write(*,*)'set either one to 0 to skip the integration step'
      read(*,*) it1,it2
      write(11,'(2(1x,i4),24x,a)') it1,it2,'! itmx1,itmx2'
      iseld=1
      if(it1.eq.0.or.it2.eq.0)iseld=0
c---------------------------------------------------------------
c Integration parameters
c
      if(iinput.eq.1)then
        if(iseld.eq.1)then
          write(*,*)' '
          write(*,*)
     #     'Enter number of calls for bases (<0 for default)'
          read(*,*)ncl3
          write(11,'(1x,i9,24x,a)')ncl3,'! # of calls for bases'
        endif
      else
        ncl3=-1
      endif
      if(ncl3.lt.0)ncl3=120000
c---- close logfile
      close(11)
c----------------------------------------------------------------
c  *********************  START INTEGRATION *********************
c----------------------------------------------------------------
      ifuntype=0
      loproc = 1
      maproc = 3
      if(iproc.ne.0) then
        loproc=iproc
        maproc=iproc
      endif
c When top decays, compute the relevant parameters
      if(idec.eq.0)then
        if( gammay1.ne.0.d0.and.twidth.eq.0.d0 )then
          write(*,*)'Non-zero mass range require non-zero width'
          stop
        endif
        if(gammay1.eq.0)then
          jwidth=0
          ym1low2=-1.d0
          ym1upp2=-1.d0
          tbw1delf=0.d0
          xm2low2=-1.d0
          xm2upp2=-1.d0
          bw2delf=0.d0
        elseif(gammay1.ne.0)then
          write(*,*)'Off-shell effects not yet implemented'
          stop
          jwidth=1
          tga1=twidth
          ga2=gaw
          if(gammay1.ge.0)then
            ym1low2=(max( 10.d0,xlep1mass(1)+xlep2mass(1)+xmass(5),
     #                    xmt-gammay1*tga1 ))**2
            ym1upp2=(min(ecm-10.d0,xmt+gammay1*tga1))**2
          else
            ym1low2=(max( 10.d0,xlep1mass(1)+xlep2mass(1)+xmass(5),
     #                    ym1low) )**2
            ym1upp2=(min(ecm-10.d0,ym1upp))**2
          endif
          if(ym1low2.gt.ym1upp2)then
            write(*,*)'Error in top mass range'
            write(*,*)ym1low2,ym1upp2
            stop
          endif
c Parameters for the skewed Breit Wigner function
          tga1mn=tga1
          tga1pl=1.15d0*tga1
          tbw1mdpl=ym1upp2-xm12
          tbw1mdmn=xm12-ym1low2
          tbw1fmpl=tga1pl/tga1*atan(tbw1mdpl/(xmt*tga1pl))
          tbw1fmmn=tga1mn/tga1*atan(tbw1mdmn/(xmt*tga1mn))
          tbw1delf=(tbw1fmpl+tbw1fmmn)/pi
c
          if(gammax2.ge.0)then
            xm2low2=(max( 10.d0,xlep1mass(2)+xlep2mass(2),
     #                    xmw-gammax2*ga2 ))**2
            xm2upp2=(min(ecm-10.d0,xmw+gammax2*ga2))**2
          else
            xm2low2=(max( 10.d0,xlep1mass(2)+xlep2mass(2),
     #                    xm2low) )**2
            xm2upp2=(min(ecm-10.d0,xm2upp))**2
          endif
          if(xm2low2.gt.xm2upp2)then
            write(*,*)'Error in hard W mass range'
            write(*,*)xm2low2,xm2upp2
            stop
          endif
c Parameters for the skewed Breit Wigner function
          ga2mn=ga2
          ga2pl=1.15d0*ga2
          bw2mdpl=xm2upp2-xmw2
          bw2mdmn=xmw2-xm2low2
          bw2fmpl=ga2pl/ga2*atan(bw2mdpl/(xmw*ga2pl))
          bw2fmmn=ga2mn/ga2*atan(bw2mdmn/(xmw*ga2mn))
          bw2delf=(bw2fmpl+bw2fmmn)/pi
        else
          write(*,*)'Unknown option: gammay1=',gammay1
          stop
        endif
c
        if( gammax1.ne.0.d0.and.gaw.eq.0.d0 )then
          write(*,*)'Non-zero W mass range requires non-zero width'
          stop
        endif
        xm012=xmw2
        if(gammax1.eq.0)then
          iwidth=0
          xm1low2=-1.d0
          xm1upp2=-1.d0
          bw1delf=0.d0
        else
          iwidth=1
          ga1=gaw
          if(gammax1.ge.0)then
            xm1low2=(max( 1.d-1,xlep1mass(1)+xlep2mass(1),
     #                    xmw-gammax1*ga1 ))**2
            xm1upp2=(xmw+gammax1*ga1)**2
          else
            xm1low2=(max( 1.d-1,xlep1mass(1)+xlep2mass(1),
     #                    xm1low ))**2
            xm1upp2=xm1upp**2
          endif
          if(jwidth.eq.0)then
            xm1upp2=min((xmt-xmass(5))**2-1.d-1,xm1upp2)
          else
            xm1upp2=min((sqrt(ym1upp2)-xmass(5))**2-1.d-1,xm1upp2)
          endif
          if(xm1low2.gt.xm1upp2)then
            write(*,*)'Error in decay W mass range'
            write(*,*)xm1low2,xm1upp2
            stop
          endif
c Parameters for the Breit Wigner
          bw1mdpl=xm1upp2-xmw2
          bw1mdmn=xmw2-xm1low2
          bw1fmpl=atan(bw1mdpl/(xmw*ga1))
          bw1fmmn=atan(bw1mdmn/(xmw*ga1))
          bw1delf=(bw1fmpl+bw1fmmn)/pi
        endif
c Initialize parameters relevant to decay; MadEvent is initialized 
c in setpardec()
        call setpardec()
      endif
c
      prefn = pref
      prefnev = prefev
c tau generated according to a flat distribution in (1/tau)**nsamp
      nsamp = 1
c
      ndim=6
      nwild=6
      rohlim=(xm1+xm2)**2/sh
      xicut=xicutss*(1-rohlim)
c Perform the integration step
      if(iseld.eq.1)then
        xtotal=0.d0
        ytotal=0.d0
        dtot=0.d0
        xares=0.d0
        yares=0.d0
        xbres=0.d0
        ybres=0.d0
        av3a=0.d0
        d3a=0.d0
        av3nega=0.d0
        d3nega=0.d0
        av3b=0.d0
        d3b=0.d0
        av3negb=0.d0
        d3negb=0.d0
        mx_of_evta=0
        mx_of_evtb=0
c
        ifuntype=1
        call fk88strcat(prefn,'_a',fnamea)
        call run_bases(sig5azw_ht,fnamea,ndim,nwild,ncl3,it1,it2,
     #    ac1,ac2,av3a,d3a,av3nega,d3nega,ctime,itd1,itd2,iseed0,
     #    ibswrite,ibscall)
        write(*,*)'   '
        write(*,*)'|integral[a]|:',av3a,' +- ',d3a
        write(*,*)' integral[a] :',av3nega,' +- ',d3nega
        xares=av3a
        yares=av3nega
        xtotal=xtotal+xares
        ytotal=ytotal+yares
        dtot=dtot+d3nega**2
c
        ifuntype=2
        call fk88strcat(prefn,'_b',fnameb)
        call run_bases(sig5bzw_ht,fnameb,ndim,nwild,ncl3,it1,it2,
     #    ac1,ac2,av3b,d3b,av3negb,d3negb,ctime,itd1,itd2,iseed0,
     #    ibswrite,ibscall)
        write(*,*)'   '
        write(*,*)'|integral[b]|:',av3b,' +- ',d3b
        write(*,*)' integral[b] :',av3negb,' +- ',d3negb
        xbres=av3b
        ybres=av3negb
        xtotal=xtotal+xbres
        ytotal=ytotal+ybres
        dtot=dtot+d3negb**2
c
        avtot=ytotal
        dtot=sqrt(dtot)
        call fk88strcat(prefn,'.integrals',fname)
        open(unit=21,file=fname,
     #       form='formatted',status='unknown')
        write(21,240)xares
        write(21,240)xbres
        write(21,240)yares
        write(21,240)ybres
        close(21)
 240    format(1x,d14.8)
      endif
c Sanity check
      if(isubttype.eq.1.and.deltai.ne.1.d0)then
        write(*,*)'Fatal error: xicut, deltaI=',xicut,deltai
        stop
      endif
      if(iseld.eq.0)then
c Read integrals from disk only if the integration step has been skipped
        call fk88strcat(prefn,'.integrals',fname)
        open(unit=21,file=fname,
     #       form='formatted',status='old')
        read(21,240)xares
        read(21,240)xbres
        read(21,240)yares
        read(21,240)ybres
        close(21)
      endif
c
c Generates events when evgen=.true.; if evgen=.false., maxevt=100000 in
c order to estimate the number of negative weights
      if(maxevt.ne.0)then
        ntotal=0
        xtotal=0.d0
        ytotal=0.d0
        xtotal=xtotal+xares+xbres
        ytotal=ytotal+yares+ybres
        avtot=ytotal
        if(iseld.eq.0)dtot=0.d0
c For future upgrades, define the weights of H and S events; this is 
c necessary when the relative number of H and S events is not generated
c according to total rates
        if(iwgtnorm.eq.0)then
          wgtaev=1.d0
          wgtbev=1.d0
        else
          wgtaev=xtotal/dfloat(maxevt)
          wgtbev=xtotal/dfloat(maxevt)
        endif
c Maximum weight cannot be set here owing to possible branching ratios
        wgtmax=-1.d8
        mx_of_evta=int(maxevt*xares/xtotal)
        mx_of_evtb=int(maxevt*xbres/xtotal)
        ntotal=ntotal+mx_of_evta+mx_of_evtb
        ndiff=maxevt-ntotal
        if(ndiff.gt.0)mx_of_evta=mx_of_evta+ndiff
        if(ndiff.lt.0)then
          write(6,*)'Fatal error:',maxevt,ntotal
          stop
        endif
        if(evgen)then
          write(*,*)'  '
          write(*,*)
     #  'The following number of events will be generated'
          write(*,*)'# events[a]:',mx_of_evta
          write(*,*)'# events[b]:',mx_of_evtb
        endif
        write(*,*)'  '
        write(*,*)
     #  'Estimated fractions of events with negative weights'
        evfrac=0.d0
        evprcfrac=(xares-yares)/
     #            (xares+yares)
        evprcfrac=evprcfrac/(1+evprcfrac)
        evfrac=evfrac+evprcfrac*mx_of_evta
        write(*,*)'Events[a]: w<0/all:',evprcfrac
        evprcfrac=(xbres-ybres)/
     #            (xbres+ybres)
        evprcfrac=evprcfrac/(1+evprcfrac)
        write(*,*)'Events[b]: w<0/all:',evprcfrac
        evfrac=evfrac+evprcfrac*mx_of_evtb
        evfrac=evfrac/dfloat(maxevt)
        write(*,*)'Events[all]: w<0/all:',evfrac
c
        if(.not.evgen)goto 111
        fname=prefnev
        call fk88strcat(fname,'_a.events',fname1)
        open(unit=22,file=fname1,
     #       form='formatted',status='unknown')
        write(22,250)mx_of_evta
        close(22)
        call fk88strcat(fname,'_b.events',fname1)
        open(unit=22,file=fname1,
     #       form='formatted',status='unknown')
        write(22,250)mx_of_evtb
        close(22)
c
        fname=prefn
        fnamev=prefnev
c
        iwrong=0
        iwrong1=0
        neventsuw=0
        nqeventsuw=0
        ifailuw=0
        ncntuws=0
        nqcntuws=0
        nmaxuw=0
        nqmaxuw=0
        do i=1,4
          xmomshifts(i)=0.d0
        enddo
        ifuntype=1
        call fk88strcat(fname,'_a',fnamea)
        call fk88strcat(fnamev,'_a.events',fname1)
        open(unit=22,file=fname1,
     #       form='formatted',status='old',access='append')
        call run_spring(sig5azw_ht,fnamea,mx_of_evta,maxtrials,
     #                  nevts,ntrls,ndim,nwild,iseed)
        close(22)
        if(iverbose.eq.1)then
          write(*,*)'   '
          write(*,*)'Events[a]'
          write(*,*)'Trials:',ntrls
          write(*,*)'Events generated:',nevts
          write(*,*)'Unlike sign events(1):',iwrong
          write(*,*)'Unlike sign events(2):',iwrong1
          write(*,*)'Unlike sign(1)/all events:',
     #              iwrong/dfloat(nevts)
          write(*,*)'Unlike sign(2)/all events:',
     #              iwrong1/dfloat(nevts)
          if(idec.eq.0)then
            if(neventsuw.ne.mx_of_evta)then
              write(*,*)'Error in spin correlations [a]'
              stop
            endif
            write(*,*)'   '
            write(*,*)'Secondary unweighting for spin correlations'
            write(*,*)'Failures',ifailuw
            write(*,*)'Average trials',ncntuws/dfloat(neventsuw)
            write(*,*)'Maximum trials',nmaxuw
            write(*,*)'Efficiency',neventsuw/dfloat(ncntuws)
            if(iwidth.eq.1)then
              write(6,*)'Maximum trials [Q]',nqmaxuw
              write(6,*)'Efficiency [Q]',
     #                  nqeventsuw/dfloat(nqcntuws)
            endif
          endif
          write(*,*)'   '
          write(*,*)'Average momentum shifts due to masses'
          do i=1,4
            if(idec.eq.0)then
              write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(8*nevts)
            else
              write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(5*nevts)
            endif
          enddo
        endif
c
        iwrong=0
        iwrong1=0
        neventsuw=0
        nqeventsuw=0
        ifailuw=0
        ncntuws=0
        nqcntuws=0
        nmaxuw=0
        nqmaxuw=0
        do i=1,4
          xmomshifts(i)=0.d0
        enddo
        ifuntype=2
        call fk88strcat(fname,'_b',fnameb)
        call fk88strcat(fnamev,'_b.events',fname1)
        open(unit=22,file=fname1,
     #       form='formatted',status='old',access='append')
        call run_spring(sig5bzw_ht,fnameb,mx_of_evtb,maxtrials,
     #                  nevts,ntrls,ndim,nwild,iseed)
        close(22)
        if(iverbose.eq.1)then
          write(*,*)'   '
          write(*,*)'Events[b]'
          write(*,*)'Trials:',ntrls
          write(*,*)'Events generated:',nevts
          write(*,*)'Unlike sign events(1):',iwrong
          write(*,*)'Unlike sign events(2):',iwrong1
          write(*,*)'Unlike sign(1)/all events:',
     #              iwrong/dfloat(nevts)
          write(*,*)'Unlike sign(2)/all events:',
     #              iwrong1/dfloat(nevts)
          if(idec.eq.0)then
            if(neventsuw.ne.mx_of_evtb)then
              write(*,*)'Error in spin correlations [b]'
              stop
            endif
            write(*,*)'   '
            write(*,*)'Secondary unweighting for spin correlations'
            write(*,*)'Failures',ifailuw
            write(*,*)'Average trials',ncntuws/dfloat(neventsuw)
            write(*,*)'Maximum trials',nmaxuw
            write(*,*)'Efficiency',neventsuw/dfloat(ncntuws)
            if(iwidth.eq.1)then
              write(6,*)'Maximum trials [Q]',nqmaxuw
              write(6,*)'Efficiency [Q]',
     #                  nqeventsuw/dfloat(nqcntuws)
            endif
          endif
          write(*,*)'   '
          write(*,*)'Average momentum shifts due to masses'
          do i=1,4
            if(idec.eq.0)then
              write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(7*nevts)
            else
              write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(4*nevts)
            endif
          enddo
        endif
c write a single event file
        iunita=21
        call fk88strcat(prefnev,'_a.events',fname1)
        open(unit=iunita,file=fname1,form='formatted',status='old')
        read(iunita,250)mx_of_evta
        iunitb=22
        call fk88strcat(prefnev,'_b.events',fname1)
        open(unit=iunitb,file=fname1,form='formatted',status='old')
        read(iunitb,250)mx_of_evtb
c
        call fk88strcat(prefnev,'.events',fname1)
        ioutput=30
        open(unit=ioutput,file=fname1,form='formatted',
     #       status='unknown')
        if(ievffmt.eq.1)then
c Rescale maximum weight for consistency with what done with XWGTUP
          wgtmax=wgtmax*maxevt
c Write LHEF header
          write(ioutput,'(a)')
     # '<LesHouchesEvents version="1.0">'
          write(ioutput,'(a)')
     # '  <!--'
        endif
c Write all the quantities which identify the run
        write(ioutput,801)
     #    ecm,xren,xfh,xrenmc,xfhmc,
     #    ': CM energy, muR/mu0[NLO], muF/mu0[NLO], '//
     #    'muR/mu0[MC], muF/mu0[MC]'
        write(ioutput,802)abs(iprdct0hw),': 0/1/4 -> t+tb/tb/t'
        write(ioutput,803)xm1,twidth,': M_top, Gamma_top'
        write(ioutput,803)xm2,hwidth,': M_H, Gamma_H'
        write(ioutput,815)il1hw,': IL1 (0..7)'
        write(ioutput,804)xmass(1),xmass(2),
     #                    xmass(3),xmass(4),
     #                    xmass(5),xmass(21),
     #                    ': quark and gluon masses'
        write(ioutput,805)part1,part2,': colliding particles'
        write(ioutput,806)gname(1:8),idpdfset,
     #    ': PDF group and id number'
        write(ioutput,807)xlam,scheme,': Lambda_5, scheme'
        write(ioutput,811)'H++',': Herwig++ (v3.3 and higher)'
        write(ioutput,250)maxevt
        if(ievffmt.eq.1)then
c Close LHEF header
          if(idec.eq.0)then
            xsecpp=avtot*brrtop1*brrv2msb
            xerrpp=dtot*brrtop1*brrv2msb
          else
            xsecpp=avtot
            xerrpp=dtot
          endif
          write(ioutput,'(a)')
     # '  -->'
          write(ioutput,'(a)')
     # '  <header>'
          write(ioutput,'(a)')
     # '  </header>'
          call write_lhef_init(ioutput,iprdct0hw,ecm,
     #                         xsecpp,xerrpp,part1,part2)
        endif
        itot=maxevt
        do ii=1,maxevt
          call whichone(iseed,itot,mx_of_evta,mx_of_evtb,iunit)
          call retrieve_events(iunit,ii,dummy)
          call store_events(ioutput,xmone)
        enddo
        if(ievffmt.eq.1)then
c Close LHEF file
          write(ioutput,'(a)')
     # '</LesHouchesEvents>'
        endif
        call crosscheck(itot,mx_of_evta,mx_of_evtb)
        close(iunita)
        close(iunitb)
        close(ioutput)
 111    continue
      endif
      if(idec.eq.0)then
        write(*,*) '   '
        write(*,*)'Branching ratios used in the computation:'
        write(*,*)' BR(t -> e nu [b+d+s])= ',xbrrtoplep
        write(*,*)' BR(t -> [udbar+usbar+ubbar] [b+d+s])= ',
     #            xbrrtophad
        write(*,*) '   '
        write(*,*)'Normalization factor due to decays:',
     #            brrtop1*brrv2msb
      endif 
      write(*,*) '   '
      write(*,*) 'Total for fully inclusive'
      write(*,200)ih1,ih2,ndns1,ndns2,nl,xlam
      write(*,202) 
      write(*,270)xm1
      write(*,201) 'tot'
      write(*,300)ecm,xfh,xren,avtot,dtot
 200  format(' had1=',i2,'  had2=',i2,'  strf1=',i6,'  strf2=',i6,
     #  '  nl=',i2,'  lambda5=',d10.4)
 201  format(' ecm or ebeam  xf   xr   ',a,
     # '        err    ')
 202  format(' M_top')
 270  format(1x,1pd9.3)
 300  format((1x,1pd9.3),4x,2(1x,0pf4.2),2(1x,0pd10.4))
 250  format(1x,i8)
 801  format(5(1x,d10.4),1x,a)
 802  format(1x,i6,1x,a)
 803  format(2(1x,d10.4),1x,a)
 804  format(6(1x,d10.4),1x,a)
 805  format(2(1x,a4),1x,a)
 806  format(1x,a8,1x,i6,1x,a)
 807  format(1x,d10.4,1x,a2,1x,a)
 810  format(2(1x,i2),1x,a)
 811  format(1x,a3,1x,a)
 813  format(3(1x,d10.4),1x,a)
 814  format(1x,d10.4,1x,a)
 815  format(1x,i2,1x,a)
      end


      subroutine getset(str,ndns,ih)
      implicit real * 8 (a-h,o-z)
      character * (*) str
 2    write(*,*) str
      write(*,*)
     # '   (< 0 for a display of the features of the various sets'
      read(*,*) ndns
      if(ndns.lt.0) then
        call prntsf
        go to 2
      endif
      end


      subroutine toend(iunit)
      ios = 0    
      dowhile(ios.eq.0)
         read(unit=iunit,fmt='(1x)',iostat=ios)
      enddo                        
      backspace(iunit)
      end


      subroutine getihpart(part,ih)
c Converts particle naming conventions, for Herwig to MLM
      implicit real * 8 (a-h,o-z)
      character * 4 part
c
      ih=-100
      if(part.eq.'P   ')then
        ih=1
      elseif(part.eq.'PBAR')then
        ih=-1
      elseif(part.eq.'N   ')then
        ih=2
      elseif(part.eq.'NBAR')then
        ih=-2
      elseif(part.eq.'GAMA')then
        ih=4
      elseif(part.eq.'E-  ')then
        ih=5
      else
        write(*,*)'Error in getihpart'
        write(*,*)'No such particle in MLM:',part
        stop
      endif
      return
      end


      subroutine strfunht(x1,x2,sf)
c Return parton densities through the matrix
c  sf(idr,jproc,itype,ittbar), with the following conventions:
c   idr -> identifies the partonic process given jproc and ich
c   jproc=1,2,3 -> gg, q(bar)q(bar), q(bar)g processes respectively
c   itype -> identifies the individual contribution to a given jproc
c   ittbar=1,2 -> t production (IT=4), tbar production (IT=1), IT being
c                 HERWIG v6.50? labeling convention
c Note in the following that CKM matrix elements are weighted by the matrix of
c Higgs couplings!
c ckm2(i,j)=|CKM matrix elements|^2, with  i=1,4,6 --> up,charm,top
c                                          j=2,3,5 --> down,strange,bottom
c and the following combination must be also defined (here in setpar)
c ruckm   = |V_ud|^2+|V_us|^2+|V_ub|^2
c rcckm   = |V_cd|^2+|V_cs|^2+|V_cb|^2
c rtckm   = |V_td|^2+|V_ts|^2+|V_tb|^2
c rducckm = |V_ud|^2+|V_cd|^2
c rsucckm = |V_us|^2+|V_cs|^2
c rbucckm = |V_ub|^2+|V_cb|^2
      implicit none
      real*4 fh1x1(-5:5),fh2x2(-5:5),smuf2h1,smuf2h2
      real * 8 pi,x1,x2,sf(7,3,24,2)
      integer ih1,ih2,ndns1,ndns2,ii,jproc,itype,ittbar,jht
      parameter(pi=3.14159265358979312D0)
      parameter (jht=3)
      include 'stpcblks.h'
      common/strfun0/ih1,ih2,ndns1,ndns2
      real*8 zel(1:6),zel2(1:6)
      common/charges/zel,zel2
      real*8 cab2(1:6,1:6)
      common/ccab2/cab2
      real*8 rucab,rccab,rtcab,rduccab,rsuccab,rbuccab
      common/ccabfct/rucab,rccab,rtcab,rduccab,rsuccab,rbuccab
      integer ichmin,ichmax
      common/cichrange/ichmin,ichmax
      integer ittmin,ittmax
      common/cittrange/ittmin,ittmax
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer idrmax(3,3)
      common/cidrmax/idrmax
      integer ittintf
      common/cittintf/ittintf
      integer ipdfscale
      common/cipdfscale/ipdfscale
c ipdfscale=1 --> use NLO factorization scale
c ipdfscale=2 --> use MC factorization scale
c
      do jproc=1,3
        do ii=1,idrmax(jproc,jht)
          do itype=1,jtypemax(jproc,ii)
            do ittbar=ittmin,ittmax
              sf(ii,jproc,itype,ittbar)=0.d0
            enddo
          enddo
        enddo
      enddo
c
      if(ipdfscale.eq.1)then
        smuf2h1=sngl(xmuf2h1)
        smuf2h2=sngl(xmuf2h2)
      elseif(ipdfscale.eq.2)then
        smuf2h1=sngl(xmumcf2h1)
        smuf2h2=sngl(xmumcf2h2)
      else
        write(*,*)'Fatal error in strfun: unknown ipdfscale',ipdfscale
        stop
      endif
c
      call mlmpdf(ndns1,ih1,smuf2h1,sngl(x1),fh1x1,5)
      call mlmpdf(ndns2,ih2,smuf2h2,sngl(x2),fh2x2,5)
c
c jproc=1
      if(ittmin.eq.1)then
c t production
        sf(1,1,1,1)=dble(fh1x1( 0) * fh2x2( 0))*rtcab
      endif
      if(ittmax.eq.2)then
c tbar production
        sf(1,1,1,2)=dble(fh1x1( 0) * fh2x2( 0))*rtcab
      endif
c jproc=2
      if(ittmin.eq.1)then
c t production
        sf(1,2,1,1)=dble(fh1x1( 1) * fh2x2(-1))*rtcab
        sf(1,2,2,1)=dble(fh1x1( 4) * fh2x2(-4))*rtcab
        sf(1,2,3,1)=dble(fh1x1( 2) * fh2x2(-2))*(rtcab-cab2(6,2))
        sf(1,2,4,1)=dble(fh1x1( 3) * fh2x2(-3))*(rtcab-cab2(6,3))
c
        sf(2,2,1,1)=dble(fh1x1(-1) * fh2x2( 1))*rtcab
        sf(2,2,2,1)=dble(fh1x1(-4) * fh2x2( 4))*rtcab
        sf(2,2,3,1)=dble(fh1x1(-2) * fh2x2( 2))*(rtcab-cab2(6,2))
        sf(2,2,4,1)=dble(fh1x1(-3) * fh2x2( 3))*(rtcab-cab2(6,3))
c
        sf(3,2, 1,1)=dble(fh1x1( 2) * fh2x2( 1))*cab2(6,2)
        sf(3,2, 2,1)=dble(fh1x1( 2) * fh2x2( 4))*cab2(6,2)
        sf(3,2, 3,1)=dble(fh1x1( 2) * fh2x2(-1))*cab2(6,2)
        sf(3,2, 4,1)=dble(fh1x1( 2) * fh2x2(-4))*cab2(6,2)
        sf(3,2, 5,1)=dble(fh1x1( 3) * fh2x2( 1))*cab2(6,3)
        sf(3,2, 6,1)=dble(fh1x1( 3) * fh2x2( 4))*cab2(6,3)
        sf(3,2, 7,1)=dble(fh1x1( 3) * fh2x2(-1))*cab2(6,3)
        sf(3,2, 8,1)=dble(fh1x1( 3) * fh2x2(-4))*cab2(6,3)
        sf(3,2, 9,1)=dble(fh1x1( 5) * fh2x2( 1))*cab2(6,5)
        sf(3,2,10,1)=dble(fh1x1( 5) * fh2x2( 4))*cab2(6,5)
        sf(3,2,11,1)=dble(fh1x1( 5) * fh2x2(-1))*cab2(6,5)
        sf(3,2,12,1)=dble(fh1x1( 5) * fh2x2(-4))*cab2(6,5)
        sf(3,2,13,1)=dble(fh1x1( 2) * fh2x2( 3))*cab2(6,2)
        sf(3,2,14,1)=dble(fh1x1( 2) * fh2x2( 5))*cab2(6,2)
        sf(3,2,15,1)=dble(fh1x1( 2) * fh2x2(-3))*cab2(6,2)
        sf(3,2,16,1)=dble(fh1x1( 2) * fh2x2(-5))*cab2(6,2)
        sf(3,2,17,1)=dble(fh1x1( 3) * fh2x2( 2))*cab2(6,3)
        sf(3,2,18,1)=dble(fh1x1( 3) * fh2x2( 5))*cab2(6,3)
        sf(3,2,19,1)=dble(fh1x1( 3) * fh2x2(-2))*cab2(6,3)
        sf(3,2,20,1)=dble(fh1x1( 3) * fh2x2(-5))*cab2(6,3)
        sf(3,2,21,1)=dble(fh1x1( 5) * fh2x2( 2))*cab2(6,5)
        sf(3,2,22,1)=dble(fh1x1( 5) * fh2x2( 3))*cab2(6,5)
        sf(3,2,23,1)=dble(fh1x1( 5) * fh2x2(-2))*cab2(6,5)
        sf(3,2,24,1)=dble(fh1x1( 5) * fh2x2(-3))*cab2(6,5)
c
        sf(4,2, 1,1)=dble(fh1x1( 1) * fh2x2( 2))*cab2(6,2)
        sf(4,2, 2,1)=dble(fh1x1( 4) * fh2x2( 2))*cab2(6,2)
        sf(4,2, 3,1)=dble(fh1x1(-1) * fh2x2( 2))*cab2(6,2)
        sf(4,2, 4,1)=dble(fh1x1(-4) * fh2x2( 2))*cab2(6,2)
        sf(4,2, 5,1)=dble(fh1x1( 1) * fh2x2( 3))*cab2(6,3)
        sf(4,2, 6,1)=dble(fh1x1( 4) * fh2x2( 3))*cab2(6,3)
        sf(4,2, 7,1)=dble(fh1x1(-1) * fh2x2( 3))*cab2(6,3)
        sf(4,2, 8,1)=dble(fh1x1(-4) * fh2x2( 3))*cab2(6,3)
        sf(4,2, 9,1)=dble(fh1x1( 1) * fh2x2( 5))*cab2(6,5)
        sf(4,2,10,1)=dble(fh1x1( 4) * fh2x2( 5))*cab2(6,5)
        sf(4,2,11,1)=dble(fh1x1(-1) * fh2x2( 5))*cab2(6,5)
        sf(4,2,12,1)=dble(fh1x1(-4) * fh2x2( 5))*cab2(6,5)
        sf(4,2,13,1)=dble(fh1x1( 3) * fh2x2( 2))*cab2(6,2)
        sf(4,2,14,1)=dble(fh1x1( 5) * fh2x2( 2))*cab2(6,2)
        sf(4,2,15,1)=dble(fh1x1(-3) * fh2x2( 2))*cab2(6,2)
        sf(4,2,16,1)=dble(fh1x1(-5) * fh2x2( 2))*cab2(6,2)
        sf(4,2,17,1)=dble(fh1x1( 2) * fh2x2( 3))*cab2(6,3)
        sf(4,2,18,1)=dble(fh1x1( 5) * fh2x2( 3))*cab2(6,3)
        sf(4,2,19,1)=dble(fh1x1(-2) * fh2x2( 3))*cab2(6,3)
        sf(4,2,20,1)=dble(fh1x1(-5) * fh2x2( 3))*cab2(6,3)
        sf(4,2,21,1)=dble(fh1x1( 2) * fh2x2( 5))*cab2(6,5)
        sf(4,2,22,1)=dble(fh1x1( 3) * fh2x2( 5))*cab2(6,5)
        sf(4,2,23,1)=dble(fh1x1(-2) * fh2x2( 5))*cab2(6,5)
        sf(4,2,24,1)=dble(fh1x1(-3) * fh2x2( 5))*cab2(6,5)
c
        sf(5,2,1,1)=dble(fh1x1( 2) * fh2x2(-2))*cab2(6,2)
        sf(5,2,2,1)=dble(fh1x1( 3) * fh2x2(-3))*cab2(6,3)
        sf(5,2,3,1)=dble(fh1x1( 5) * fh2x2(-5))*cab2(6,5)
c
        sf(6,2,1,1)=dble(fh1x1(-2) * fh2x2( 2))*cab2(6,2)
        sf(6,2,2,1)=dble(fh1x1(-3) * fh2x2( 3))*cab2(6,3)
        sf(6,2,3,1)=dble(fh1x1(-5) * fh2x2( 5))*cab2(6,5)
c
        sf(7,2,1,1)=dble(fh1x1( 2) * fh2x2( 2))*cab2(6,2)
        sf(7,2,2,1)=dble(fh1x1( 3) * fh2x2( 3))*cab2(6,3)
        sf(7,2,3,1)=dble(fh1x1( 5) * fh2x2( 5))*cab2(6,5)
      endif
      if(ittmax.eq.2)then
c tbar production
        sf(1,2,1,2)=dble(fh1x1(-1) * fh2x2( 1))*rtcab
        sf(1,2,2,2)=dble(fh1x1(-4) * fh2x2( 4))*rtcab
        sf(1,2,3,2)=dble(fh1x1(-2) * fh2x2( 2))*(rtcab-cab2(6,2))
        sf(1,2,4,2)=dble(fh1x1(-3) * fh2x2( 3))*(rtcab-cab2(6,3))
c
        sf(2,2,1,2)=dble(fh1x1( 1) * fh2x2(-1))*rtcab
        sf(2,2,2,2)=dble(fh1x1( 4) * fh2x2(-4))*rtcab
        sf(2,2,3,2)=dble(fh1x1( 2) * fh2x2(-2))*(rtcab-cab2(6,2))
        sf(2,2,4,2)=dble(fh1x1( 3) * fh2x2(-3))*(rtcab-cab2(6,3))
c
        sf(3,2, 1,2)=dble(fh1x1(-2) * fh2x2(-1))*cab2(6,2)
        sf(3,2, 2,2)=dble(fh1x1(-2) * fh2x2(-4))*cab2(6,2)
        sf(3,2, 3,2)=dble(fh1x1(-2) * fh2x2( 1))*cab2(6,2)
        sf(3,2, 4,2)=dble(fh1x1(-2) * fh2x2( 4))*cab2(6,2)
        sf(3,2, 5,2)=dble(fh1x1(-3) * fh2x2(-1))*cab2(6,3)
        sf(3,2, 6,2)=dble(fh1x1(-3) * fh2x2(-4))*cab2(6,3)
        sf(3,2, 7,2)=dble(fh1x1(-3) * fh2x2( 1))*cab2(6,3)
        sf(3,2, 8,2)=dble(fh1x1(-3) * fh2x2( 4))*cab2(6,3)
        sf(3,2, 9,2)=dble(fh1x1(-5) * fh2x2(-1))*cab2(6,5)
        sf(3,2,10,2)=dble(fh1x1(-5) * fh2x2(-4))*cab2(6,5)
        sf(3,2,11,2)=dble(fh1x1(-5) * fh2x2( 1))*cab2(6,5)
        sf(3,2,12,2)=dble(fh1x1(-5) * fh2x2( 4))*cab2(6,5)
        sf(3,2,13,2)=dble(fh1x1(-2) * fh2x2(-3))*cab2(6,2)
        sf(3,2,14,2)=dble(fh1x1(-2) * fh2x2(-5))*cab2(6,2)
        sf(3,2,15,2)=dble(fh1x1(-2) * fh2x2( 3))*cab2(6,2)
        sf(3,2,16,2)=dble(fh1x1(-2) * fh2x2( 5))*cab2(6,2)
        sf(3,2,17,2)=dble(fh1x1(-3) * fh2x2(-2))*cab2(6,3)
        sf(3,2,18,2)=dble(fh1x1(-3) * fh2x2(-5))*cab2(6,3)
        sf(3,2,19,2)=dble(fh1x1(-3) * fh2x2( 2))*cab2(6,3)
        sf(3,2,20,2)=dble(fh1x1(-3) * fh2x2( 5))*cab2(6,3)
        sf(3,2,21,2)=dble(fh1x1(-5) * fh2x2(-2))*cab2(6,5)
        sf(3,2,22,2)=dble(fh1x1(-5) * fh2x2(-3))*cab2(6,5)
        sf(3,2,23,2)=dble(fh1x1(-5) * fh2x2( 2))*cab2(6,5)
        sf(3,2,24,2)=dble(fh1x1(-5) * fh2x2( 3))*cab2(6,5)
c
        sf(4,2, 1,2)=dble(fh1x1(-1) * fh2x2(-2))*cab2(6,2)
        sf(4,2, 2,2)=dble(fh1x1(-4) * fh2x2(-2))*cab2(6,2)
        sf(4,2, 3,2)=dble(fh1x1( 1) * fh2x2(-2))*cab2(6,2)
        sf(4,2, 4,2)=dble(fh1x1( 4) * fh2x2(-2))*cab2(6,2)
        sf(4,2, 5,2)=dble(fh1x1(-1) * fh2x2(-3))*cab2(6,3)
        sf(4,2, 6,2)=dble(fh1x1(-4) * fh2x2(-3))*cab2(6,3)
        sf(4,2, 7,2)=dble(fh1x1( 1) * fh2x2(-3))*cab2(6,3)
        sf(4,2, 8,2)=dble(fh1x1( 4) * fh2x2(-3))*cab2(6,3)
        sf(4,2, 9,2)=dble(fh1x1(-1) * fh2x2(-5))*cab2(6,5)
        sf(4,2,10,2)=dble(fh1x1(-4) * fh2x2(-5))*cab2(6,5)
        sf(4,2,11,2)=dble(fh1x1( 1) * fh2x2(-5))*cab2(6,5)
        sf(4,2,12,2)=dble(fh1x1( 4) * fh2x2(-5))*cab2(6,5)
        sf(4,2,13,2)=dble(fh1x1(-3) * fh2x2(-2))*cab2(6,2)
        sf(4,2,14,2)=dble(fh1x1(-5) * fh2x2(-2))*cab2(6,2)
        sf(4,2,15,2)=dble(fh1x1( 3) * fh2x2(-2))*cab2(6,2)
        sf(4,2,16,2)=dble(fh1x1( 5) * fh2x2(-2))*cab2(6,2)
        sf(4,2,17,2)=dble(fh1x1(-2) * fh2x2(-3))*cab2(6,3)
        sf(4,2,18,2)=dble(fh1x1(-5) * fh2x2(-3))*cab2(6,3)
        sf(4,2,19,2)=dble(fh1x1( 2) * fh2x2(-3))*cab2(6,3)
        sf(4,2,20,2)=dble(fh1x1( 5) * fh2x2(-3))*cab2(6,3)
        sf(4,2,21,2)=dble(fh1x1(-2) * fh2x2(-5))*cab2(6,5)
        sf(4,2,22,2)=dble(fh1x1(-3) * fh2x2(-5))*cab2(6,5)
        sf(4,2,23,2)=dble(fh1x1( 2) * fh2x2(-5))*cab2(6,5)
        sf(4,2,24,2)=dble(fh1x1( 3) * fh2x2(-5))*cab2(6,5)
c
        sf(5,2,1,2)=dble(fh1x1(-2) * fh2x2( 2))*cab2(6,2)
        sf(5,2,2,2)=dble(fh1x1(-3) * fh2x2( 3))*cab2(6,3)
        sf(5,2,3,2)=dble(fh1x1(-5) * fh2x2( 5))*cab2(6,5)
c
        sf(6,2,1,2)=dble(fh1x1( 2) * fh2x2(-2))*cab2(6,2)
        sf(6,2,2,2)=dble(fh1x1( 3) * fh2x2(-3))*cab2(6,3)
        sf(6,2,3,2)=dble(fh1x1( 5) * fh2x2(-5))*cab2(6,5)
c
        sf(7,2,1,2)=dble(fh1x1(-2) * fh2x2(-2))*cab2(6,2)
        sf(7,2,2,2)=dble(fh1x1(-3) * fh2x2(-3))*cab2(6,3)
        sf(7,2,3,2)=dble(fh1x1(-5) * fh2x2(-5))*cab2(6,5)
      endif
c jproc=3
      if(ittmin.eq.1)then
c t production
        sf(1,3,1,1)=dble(fh1x1( 2) * fh2x2( 0))*cab2(6,2)
        sf(1,3,2,1)=dble(fh1x1( 3) * fh2x2( 0))*cab2(6,3)
        sf(1,3,3,1)=dble(fh1x1( 5) * fh2x2( 0))*cab2(6,5)
c
        sf(3,3,1,1)=dble(fh1x1( 0) * fh2x2( 2))*cab2(6,2)
        sf(3,3,2,1)=dble(fh1x1( 0) * fh2x2( 3))*cab2(6,3)
        sf(3,3,3,1)=dble(fh1x1( 0) * fh2x2( 5))*cab2(6,5)
      endif
      if(ittmax.eq.2)then
c tbar production 
        sf(1,3,1,2)=dble(fh1x1(-2) * fh2x2( 0))*cab2(6,2)
        sf(1,3,2,2)=dble(fh1x1(-3) * fh2x2( 0))*cab2(6,3)
        sf(1,3,3,2)=dble(fh1x1(-5) * fh2x2( 0))*cab2(6,5)
c
        sf(3,3,1,2)=dble(fh1x1( 0) * fh2x2(-2))*cab2(6,2)
        sf(3,3,2,2)=dble(fh1x1( 0) * fh2x2(-3))*cab2(6,3)
        sf(3,3,3,2)=dble(fh1x1( 0) * fh2x2(-5))*cab2(6,5)
      endif
c Set to zero luminosity for processes with ttbar interference
      if(ittintf.eq.1)then
        jproc=1
        ii=1
        do itype=1,jtypemax(jproc,ii)
          do ittbar=ittmin,ittmax
            sf(ii,jproc,itype,ittbar)=0.d0
          enddo
        enddo
c
        jproc=2
        ii=1
        do itype=1,jtypemax(jproc,ii)
          do ittbar=ittmin,ittmax
            sf(ii,jproc,itype,ittbar)=0.d0
          enddo
        enddo
        ii=2
        do itype=1,jtypemax(jproc,ii)
          do ittbar=ittmin,ittmax
            sf(ii,jproc,itype,ittbar)=0.d0
          enddo
        enddo
        ii=5
        do itype=1,jtypemax(jproc,ii)
          do ittbar=ittmin,ittmax
            sf(ii,jproc,itype,ittbar)=0.d0
          enddo
        enddo
        ii=6
        do itype=1,jtypemax(jproc,ii)
          do ittbar=ittmin,ittmax
            sf(ii,jproc,itype,ittbar)=0.d0
          enddo
        enddo
      elseif(ittintf.ne.0)then
        write(*,*)'Error in strfunht: unknown ittintf',ittintf
      endif
c
      return
      end
c
c
c NLO cross section
c
c
      function sig5azw_ht(xx)
c Integrand function for H events
      implicit none
      real * 8 sig5azw_ht,xx
      real * 8 pi,tiny
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-5)
      dimension xx(6)
      include 'stpcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      real * 8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
      integer iprespl
      common/ciprespl/iprespl
      integer nsamp
      common/samp/nsamp
      integer ifxdaem
      common/cifxdaem/ifxdaem
      real * 8 xjac,rohlim,zzz,x,ttt,th,yi,csi,rx,rohlimx,taumax,
     #  ximax0,ximin0,tmp,ymax,ymin,xxa1,xxa2,xxc,xxymax,xxymin,
     #  s,xmi2,xalfaem,rox,cth1,th2,cth2,tot5a_ht
c
c xx(1) --> tau, xx(2)-->ycm, xx(3) --> x, xx(4) --> y, xx(5) --> cth1,
c xx(6) --> cth2
c
      xjac = 1.d0
      rohlim=(sqrt(xm12)+sqrt(xm22))**2/sh
c
c To improve convergence in the soft regions
c
      zzz = tiny+(1-tiny)*xx(3)**2
      xjac = xjac * xx(3) * 2
      x = 1 - zzz*(1-rohlim)
      xjac = xjac * (1-rohlim)
c
c To improve convergence in the initial state collinear regions
c
      zzz = 1-2*xx(4)
      xjac = xjac * 2
      ttt = tiny+(1-tiny)*zzz**2
      xjac = xjac * 2 * abs(zzz)
      if(zzz.gt.0) then
         th = ttt * pi/2
      else
         th = pi-ttt*pi/2
      endif
      xjac = xjac * pi/2
      yi    = cos(th)
      xjac = xjac * sin(th)
c
      csi = sqrt((1-(1-x)*(1+yi)/2.d0)/(1-(1-x)*(1-yi)/2.d0))
      rx = sqrt(x)
      rohlimx = rohlim/x
      taumax = 1/x
      ximax0 = rohlimx**(-nsamp)
      ximin0 = taumax**(-nsamp)
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
      if(iprespl.eq.0)then
        ymax= -log(tau)/2 + log(1/(csi*rx))
        ymin=  log(tau)/2 - log(csi/rx)
      else
        xxa1 = (1+x-yi*(1-x))/2.d0
        xxa2 = (1+x+yi*(1-x))/2.d0
        xxc = (1-x*tau)/sqrt(tau)
        xxymax = (xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        xxymin = (-xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        ymax = max(log(xxymax),-log(tau)/2.d0)
        ymin = min(log(xxymin),log(tau)/2.d0)
      endif
      ycm = ymin + xx(2)*(ymax-ymin)
      xjac= xjac * (ymax-ymin)
      s=tau*sh
c Hard coded choice for scale of running e.m. coupling: mtop
      xmi2=xm12
      if(ifxdaem.eq.0)ze2=4*pi*xalfaem(xmi2)
c
      rox = 2*(xm12+xm22)/(s*x)-(xm12-xm22)**2/(s*x)**2
c zzchvar: a change of variables xx(5) --> cth1
      call zzchvar(xx(5),cth1,xjac,rox)
c
      th2 = xx(6) * 2 * pi
      xjac = xjac * 2* pi
      cth2 = cos(th2)
c
      sig5azw_ht = tot5a_ht(s,x,yi,cth1,cth2,xjac)
      return
      end


      function tot5a_ht(s,x,yi,cth1,cth2,xjac)
      implicit none
      real * 8 tot5a_ht,tot5as_ht,tot5az_ht,s,x,yi,cth1,cth2,
     #  xjac,tmp
      integer isubttype
      common/cisubttype/isubttype
c
      if(isubttype.eq.0)then
        tmp=tot5as_ht(s,x,yi,cth1,cth2,xjac)
      elseif(isubttype.eq.1)then
        tmp=tot5az_ht(s,x,yi,cth1,cth2,xjac)
      else
        write(*,*)'Fatal error in tot5a_ht:',isubttype
        stop
      endif
      tot5a_ht=tmp
      return
      end


      function tot5as_ht(xs,xx,xyi,xcth1,xcth2,xjac)
c Implements standard subtraction
      implicit none
      real * 8 tot5as_ht,xs,xx,xyi,xcth1,xcth2,xjac
      real * 8 pi,pi2,hc2
      parameter (pi=3.14159265358979312D0)
      parameter (pi2 = pi*pi)
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
c TeV to picobarn conversion factor: sigma (pb) = hc2 * sigma (TeV^-2)
c sigma (pb) = 10^6 * sigma (mub)
      parameter (hc2=3.8937966d2)
      integer jht
      parameter (jht=3)
      integer ione
      parameter (ione=1)
      character * 2 str
      parameter (str='p1')
      include 'stpcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      real * 8 bsfsgn
      common/cbssgn/bsfsgn
      real * 8 xevsign
      common/cxevsign/xevsign
      real * 8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      real * 8 ps,px,pyi,pcth1,pcth2
      common/cpsave/ps,px,pyi,pcth1,pcth2
      real * 8 vv(7,3,24,2),vvs(7,3,24,2)
      common/cvv/vv
      common/cvvs/vvs
      real * 8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
      real*8 emwgt_prc(3)
      common/cemwgtprc/emwgt_prc
      integer ittmin,ittmax
      common/cittrange/ittmin,ittmax
      integer idrmap(8,1:3,2)
      common/cidrmap/idrmap
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer idrmax(3,3)
      common/cidrmax/idrmax
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      integer ipdfscale
      common/cipdfscale/ipdfscale
      integer idec
      common/cidec/idec
      logical fx1x2,flagmc,lzone(7,3,6)
      real * 8 xinv(5)
      real * 8 sf(7,3,24,2)
      real * 8 vecre(7),vecmc(7,3,6),xmcz(7,3,6)
      real * 8 vecmccl(7),vecmcsf(7),vecmcsc(7)
      real * 8 gfsf(3),gfcl(3)
      real*8 s,x,yi,cth1,cth2,sx,xii,ro,beta,rox,betax,xphsp_ns,
     #  xphsp_s,x1,x2,x1t,x2t,tk,uk,q1q,q2q,zg2_nlo,zgmu2_nlo,
     #  xnorm,zg2_mc,zgmu2_mc,gfactsf,gfactcl,zhwfct,x1soft,x2soft,
     #  betafact,x1x2j,x1x2jac,ytmp,xfact,xtmp,xsum,dummy,xint
      integer jproc,j,itype,itt,jp,ileg,ie0sq,iret,i2b,itoosoftkin
c
      s = xs
      x = xx
      yi = xyi
      cth1 = xcth1
      cth2 = xcth2
      sx = x*s
      xii = 1-x
      ro = 2*(xm12+xm22)/s-(xm12-xm22)**2/s**2
      beta = sqrt(1-ro)
      rox = 2*(xm12+xm22)/sx-(xm12-xm22)**2/sx**2
      betax = sqrt(1-rox)
c The normalization of born and soft phase spaces already accounts
c for event projection
      xphsp_ns = xjac * betax * s/(2*1024*pi**4)
      xphsp_s = xjac * betax * sx/(2*1024*pi**4)
c
      do jproc=1,3
        do j=1,7
          do itype=1,24
            do itt=1,2
              vv(j,jproc,itype,itt)=0.d0
              vvs(j,jproc,itype,itt)=0.d0
            enddo
          enddo
        enddo
      enddo
c
      x1 = sqrt(tau) * exp(ycm)
      x2 = tau/x1
c
      if(x1.lt.1.and.x2.lt.1)then
        call invar_in(xm12,xm22,s,x,yi,cth1,cth2,str,
     #                tk,uk,q1q,q2q,xinv)
        zg2_nlo = zgmu2_nlo() 
        ipdfscale=1
        call strfunht(x1,x2,sf)
C This corresponds to Eq. (4.37) in FKS
        xnorm = zg2_nlo**2 * xphsp_ns
        xnorm = xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
        do jproc=loproc,maproc
          call frealht(s,x,yi,cth2,tk,uk,q1q,q2q,xinv,jproc,vecre)
          do j=1,idrmax(jproc,jht)
            do itype=1,jtypemax(jproc,j)
              do itt=ittmin,ittmax
                jp=idrmap(j,jproc,itt)
                vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #            sf(j,jproc,itype,itt)*xnorm*vecre(jp)
              enddo
            enddo
          enddo
        enddo
c
c MC subtraction terms: pure MC
c
        zg2_mc = zgmu2_mc()
        ipdfscale=2
        xnorm = zg2_mc**2 * xphsp_ns
        xnorm = xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
        do jproc=loproc,maproc
          call xmcsubthtpp(jproc,x1,x2,xm12,xm22,s,x,yi,cth1,cth2,
     #                     gfactsf,gfactcl,flagmc,lzone,xmcz,vecmc)
          gfsf(jproc)=gfactsf
          gfcl(jproc)=gfactcl
          if(flagmc)then
            do j=1,idrmax(jproc,jht)
              do ileg=1,3
                do ie0sq=1,6
                  if(lzone(j,ileg,ie0sq))then
                    if(ileg.eq.1)then
                      zhwfct=xmcz(j,ileg,ie0sq)
                      x1t=x1soft(x1,x2,x,yi)/zhwfct
                      x2t=x2soft(x1,x2,x,yi)
                      betafact=1.d0
                      fx1x2=x1t.lt.1.and.x2t.lt.1
                    elseif(ileg.eq.2)then
                      zhwfct=xmcz(j,ileg,ie0sq)
                      x1t=x1soft(x1,x2,x,yi)
                      x2t=x2soft(x1,x2,x,yi)/zhwfct
                      betafact=1.d0
                      fx1x2=x1t.lt.1.and.x2t.lt.1
                    else
                      zhwfct=1.d0
                      x1t=x1
                      x2t=x2
                      betafact=beta/betax
                      fx1x2=x1t.lt.1.and.x2t.lt.1
                    endif
                    if(fx1x2)then
                      call strfunht(x1t,x2t,sf)
                      x1x2j=x1x2jac(x1,x2,x,yi,ileg)/zhwfct
                      do itype=1,jtypemax(jproc,j)
                        do itt=ittmin,ittmax
                          jp=idrmap(j,jproc,itt)
                          vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)-
     #                      sf(j,jproc,itype,itt)*xnorm*x1x2j*
     #                      betafact*vecmc(j,ileg,ie0sq)
                        enddo
                      enddo
                    endif
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
c
c MC subtraction term: collinear ME
c
        do jproc=loproc,maproc
          if(gfsf(jproc).lt.1.d0)then
            if(yi.gt.0.d0)then
              ytmp=1.d0
              x1t=x1soft(x1,x2,x,yi)/x
              x2t=x2soft(x1,x2,x,yi)
              xfact=1.d0/( xii*(1-yi) )
            else
              ytmp=-1.d0
              x1t=x1soft(x1,x2,x,yi)
              x2t=x2soft(x1,x2,x,yi)/x
              xfact=1.d0/( xii*(1+yi) )
            endif
            if(x1t.lt.1.and.x2t.lt.1)then
              x1x2j = x1x2jac(x1,x2,x,yi,ione)/x
              call invar_in(xm12,xm22,s,x,ytmp,cth1,cth2,str,
     #                      tk,uk,q1q,q2q,xinv)
              zg2_nlo = zgmu2_nlo() 
              ipdfscale=1
              call strfunht(x1t,x2t,sf)
              xnorm = xfact * x1x2j * zg2_nlo**2 * xphsp_ns
              xnorm = xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
              call frealht(s,x,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecmccl)
              do j=1,idrmax(jproc,jht)
                do itype=1,jtypemax(jproc,j)
                  do itt=ittmin,ittmax
                    jp=idrmap(j,jproc,itt)
                    vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)-
     #                sf(j,jproc,itype,itt)*xnorm*vecmccl(jp)
                  enddo
                enddo
              enddo
            endif
          endif
c
c MC subtraction term: soft ME
c
          if(gfsf(jproc).lt.1.d0)then
            xtmp=1.d0
            x1t=x1soft(x1,x2,x,yi)
            x2t=x2soft(x1,x2,x,yi)
            if(x1t.lt.1.and.x2t.lt.1)then
              x1x2j = x1x2jac(x1,x2,x,yi,ione)
              call invar_in(xm12,xm22,sx,xtmp,yi,cth1,cth2,str,
     #                      tk,uk,q1q,q2q,xinv)
              zg2_nlo = zgmu2_nlo() 
              ipdfscale=1
              call strfunht(x1t,x2t,sf)
              xnorm = x1x2j * zg2_nlo**2 * xphsp_s 
              xnorm = xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
              xnorm = xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
              call frealht(sx,xtmp,yi,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecmcsf)
              do j=1,idrmax(jproc,jht)
                do itype=1,jtypemax(jproc,j)
                  do itt=ittmin,ittmax
                    jp=idrmap(j,jproc,itt)
                    vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)-
     #                sf(j,jproc,itype,itt)*xnorm*vecmcsf(jp)
                  enddo
                enddo
              enddo
c
c MC subtraction term: soft-collinear ME
c
              if(gfsf(jproc).lt.1.d0)then
                if(yi.gt.0.d0)then
                  ytmp=1.d0
                  xfact=1.d0/( xii*(1-yi) )
                else
                  ytmp=-1.d0
                  xfact=1.d0/( xii*(1+yi) )
                endif
                call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                        tk,uk,q1q,q2q,xinv)
                xnorm = xfact * x1x2j * zg2_nlo**2 * xphsp_s
                xnorm = - xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
                call frealht(sx,xtmp,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                       jproc,vecmcsc)
                do j=1,idrmax(jproc,jht)
                  do itype=1,jtypemax(jproc,j)
                    do itt=ittmin,ittmax
                      jp=idrmap(j,jproc,itt)
                      vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)-
     #                  sf(j,jproc,itype,itt)*xnorm*vecmcsc(jp)
                    enddo
                  enddo
                enddo
              endif
            endif
          endif
        enddo
      endif
c
      call checkvv(xsum,dummy,iret)
      if(iret.eq.1)then
        call invar_in(xm12,xm22,s,x,yi,cth1,cth2,str,
     #                tk,uk,q1q,q2q,xinv)
        if(idec.eq.0)then
          ps=s
          px=x
          pyi=yi
          pcth1=cth1
          pcth2=cth2
        endif
c Cross section in pb (momenta are in GeV)
        xint=1.d6*hc2*xsum
        xevsign=1.d0
        if(xint.lt.0.d0)xevsign=-1.d0
        i2b=itoosoftkin()
        if(i2b.eq.1)then
          xtmp=1.d0
          ytmp=1.d0
          call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                  tk,uk,q1q,q2q,xinv)
          if(idec.eq.0)then
            ps=sx
            px=xtmp
            pyi=ytmp
            pcth1=cth1
            pcth2=cth2
          endif
        endif
c Store x1, x2 and Q2 for PDF reweighting using event file
        ux1 = x1 
        ux2 = x2
        dummy = zgmu2_nlo()
        uq2 = xmuf2h1
      else
        xint=0.d0
        xevsign=1.d0
      endif
c
      bsfsgn=xevsign
      tot5as_ht=abs(xint)
c
      return
      end


C Empty now, but for zeta-subtraction, to be done,
C if the number of negative weight events for single top
C happens to be too large
      function tot5az_ht(xs,xx,xyi,xcth1,xcth2,xjac)
      implicit none
      real * 8 tot5az_ht,xs,xx,xyi,xcth1,xcth2,xjac
      return
      end


      function sig5bzw_ht(xx)
c Integrand function for S events
      implicit none
      real * 8 sig5bzw_ht,xx
      real * 8 pi,tiny
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-5)
      dimension xx(6)
      include 'stpcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      real * 8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
      integer iprespl
      common/ciprespl/iprespl
      integer nsamp
      common/samp/nsamp
      integer ifxdaem
      common/cifxdaem/ifxdaem
      real * 8 xjac,rohlim,zzz,x,ttt,th,yi,csi,rx,rohlimx,taumax,
     #  ximax0,ximin0,tmp,ymax,ymin,xxa1,xxa2,xxc,xxymax,xxymin,
     #  s,xmi2,xalfaem,rox,cth1,th2,cth2,tot5b_ht
c
c xx(1) --> tau, xx(2)-->ycm, xx(3) --> x, xx(4) --> y, xx(5) --> cth1,
c xx(6) --> cth2
c
      xjac = 1.d0
      rohlim=(sqrt(xm12)+sqrt(xm22))**2/sh
c
c To improve convergence in the soft regions
c
      zzz = tiny+(1-tiny)*xx(3)**2
      xjac = xjac * xx(3) * 2
      x = 1 - zzz*(1-rohlim)
      xjac = xjac * (1-rohlim)
c
c To improve convergence in the initial state collinear regions
c
      zzz = 1-2*xx(4)
      xjac = xjac * 2
      ttt = tiny+(1-tiny)*zzz**2
      xjac = xjac * 2 * abs(zzz)
      if(zzz.gt.0) then
         th = ttt * pi/2
      else
         th = pi-ttt*pi/2
      endif
      xjac = xjac * pi/2
      yi    = cos(th)
      xjac = xjac * sin(th)
c
      csi = sqrt((1-(1-x)*(1+yi)/2.d0)/(1-(1-x)*(1-yi)/2.d0))
      rx = sqrt(x)
      rohlimx = rohlim/x
      taumax = 1/x**2
      ximax0 = rohlimx**(-nsamp)
      ximin0 = taumax**(-nsamp)
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
      if(iprespl.eq.0)then
        ymax= -log(tau)/2 + log(1/(csi*rx))
        ymin=  log(tau)/2 - log(csi/rx)
      else
        xxa1 = (1+x-yi*(1-x))/2.d0
        xxa2 = (1+x+yi*(1-x))/2.d0
        xxc = (1-x*tau)/sqrt(tau)
        xxymax = (xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        xxymin = (-xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        ymax = max(log(xxymax),-log(tau)/2.d0)
        ymin = min(log(xxymin),log(tau)/2.d0)
      endif
      ycm = ymin + xx(2)*(ymax-ymin)
      xjac= xjac * (ymax-ymin)
      s=tau*sh
c Hard coded choice for scale of running e.m. coupling: mtop
      xmi2=xm12
      if(ifxdaem.eq.0)ze2=4*pi*xalfaem(xmi2)
c
      rox = 2*(xm12+xm22)/(s*x)-(xm12-xm22)**2/(s*x)**2
c zzchvar: a change of variables xx(5) --> cth1
      call zzchvar(xx(5),cth1,xjac,rox)
c
      th2 = xx(6) * 2 * pi
      xjac = xjac * 2* pi
      cth2 = cos(th2)
c
      sig5bzw_ht = tot5b_ht(s,x,yi,cth1,cth2,xjac)
      return
      end


      function tot5b_ht(s,x,yi,cth1,cth2,xjac)
      implicit none
      real * 8 tot5b_ht,tot5bs_ht,tot5bz_ht,s,x,yi,cth1,cth2,
     #  xjac,tmp
      integer isubttype
      common/cisubttype/isubttype
c
      if(isubttype.eq.0)then
        tmp=tot5bs_ht(s,x,yi,cth1,cth2,xjac)
      elseif(isubttype.eq.1)then
        tmp=tot5bz_ht(s,x,yi,cth1,cth2,xjac)
      else
        write(*,*)'Fatal error in tot5b_ht:',isubttype
        stop
      endif
      tot5b_ht=tmp
      return
      end


      function tot5bs_ht(xs,xx,xyi,xcth1,xcth2,xjac)
      implicit none
      real * 8 tot5bs_ht,xs,xx,xyi,xcth1,xcth2,xjac
      real * 8 pi,pi2,hc2
      parameter (pi=3.14159265358979312D0)
      parameter (pi2 = pi*pi)
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
c TeV to picobarn conversion factor: sigma (pb) = hc2 * sigma (TeV^-2)
c sigma (pb) = 10^6 * sigma (mub)
      parameter (hc2=3.8937966d2)
      integer jht
      parameter (jht=3)
      integer ione
      parameter (ione=1)
      character * 2 str
      parameter (str='p1')
      include 'stpcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      real * 8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
      real * 8 bsfsgn
      common/cbssgn/bsfsgn
      real * 8 xevsign
      common/cxevsign/xevsign
      real * 8 ps,px,pyi,pcth1,pcth2
      common/cpsave/ps,px,pyi,pcth1,pcth2
      real * 8 vv(7,3,24,2),vvs(7,3,24,2)
      common/cvv/vv
      common/cvvs/vvs
      real * 8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
      real*8 emwgt_prc(3)
      common/cemwgtprc/emwgt_prc
      integer ittmin,ittmax
      common/cittrange/ittmin,ittmax
      integer idrmap(8,1:3,2)
      common/cidrmap/idrmap
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer idrmax(3,3)
      common/cidrmax/idrmax
      real * 8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      integer ipdfscale
      common/cipdfscale/ipdfscale
      integer idec
      common/cidec/idec
      logical fx1x2,flagmc,lzone(7,3,6)
      real * 8 xinv(5)
      real * 8 sf(7,3,24,2)
      real * 8 xintsvc(7,3),xbornvc(7,3)
      real * 8 xcsvc(7,3),xsvvc(7,3)
      real * 8 vecre(7),veccl(7)
      real * 8 vecbrn(7),vec2sv(7)
      real * 8 vecmc(7,3,6),xmcz(7,3,6)
      real * 8 vecmccl(7),vecmcsf(7),vecmcsc(7)
      real * 8 gfsf(3),gfcl(3)
      real * 8 s,x,yi,cth1,cth2,sx,xii,ro,beta,rox,betax,xphsp_ns,
     #  xphsp_s,xphspb,x1,x2,x1t,x2t,tk,uk,q1q,q2q,zg2_mc,
     #  zgmu2_mc,xnorm,gfactsf,gfactcl,zhwfct,x1soft,x2soft,
     #  betafact,x1x2j,x1x2jac,ytmp,xfact,zg2_nlo,zgmu2_nlo,xtmp,
     #  xlmude,xnormc,xnormb,xnormsv,xsum,dummy,xint
      integer jproc,j,itype,itt,ileg,ie0sq,jp,iret
c
      s = xs
      x = xx
      yi = xyi
      cth1 = xcth1
      cth2 = xcth2
      sx = x*s
      xii = 1-x
      ro = 2*(xm12+xm22)/s-(xm12-xm22)**2/s**2
      beta = sqrt(1-ro)
      rox = 2*(xm12+xm22)/sx-(xm12-xm22)**2/sx**2
      betax = sqrt(1-rox)
c The normalization of born and soft phase spaces already accounts
c for event projection
      xphsp_ns = xjac * betax * s/(2*1024*pi**4)
      xphsp_s = xjac * betax * sx/(2*1024*pi**4)
      xphspb = xjac * betax/(32*pi2)
c
      do jproc=1,3
        do j=1,7
          do itype=1,24
            do itt=1,2
              vv(j,jproc,itype,itt)=0.d0
              vvs(j,jproc,itype,itt)=0.d0
            enddo
          enddo
        enddo
      enddo
c
      x1 = sqrt(tau) * exp(ycm)
      x2 = tau/x1
c
      if(x1.lt.1.and.x2.lt.1)then
        call invar_in(xm12,xm22,s,x,yi,cth1,cth2,str,
     #                tk,uk,q1q,q2q,xinv)
c
c MC subtraction terms: pure MC
c
        zg2_mc = zgmu2_mc()
        ipdfscale=2
        xnorm = zg2_mc**2 * xphsp_ns
        xnorm = xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
        do jproc=loproc,maproc
          call xmcsubthtpp(jproc,x1,x2,xm12,xm22,s,x,yi,cth1,cth2,
     #                     gfactsf,gfactcl,flagmc,lzone,xmcz,vecmc)
          gfsf(jproc)=gfactsf
          gfcl(jproc)=gfactcl
          if(flagmc)then
            do j=1,idrmax(jproc,jht)
              do ileg=1,3
                do ie0sq=1,6
                  if(lzone(j,ileg,ie0sq))then
                    if(ileg.eq.1)then
                      zhwfct=xmcz(j,ileg,ie0sq)
                      x1t=x1soft(x1,x2,x,yi)/zhwfct
                      x2t=x2soft(x1,x2,x,yi)
                      betafact=1.d0
                      fx1x2=x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1
                    elseif(ileg.eq.2)then
                      zhwfct=xmcz(j,ileg,ie0sq)
                      x1t=x1soft(x1,x2,x,yi)
                      x2t=x2soft(x1,x2,x,yi)/zhwfct
                      betafact=1.d0
                      fx1x2=x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1
                    else
                      zhwfct=1.d0
                      x1t=x1soft(x1,x2,x,yi)
                      x2t=x2soft(x1,x2,x,yi)
                      betafact=1.d0
                      fx1x2=x1t.lt.1.and.x2t.lt.1.and.(x**2*tau).lt.1
                    endif
                    if(fx1x2)then
                      call strfunht(x1t,x2t,sf)
                      x1x2j=x1x2jac(x1,x2,x,yi,ione)/zhwfct
                      do itype=1,jtypemax(jproc,j)
                        do itt=ittmin,ittmax
                          jp=idrmap(j,jproc,itt)
                          vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #                      sf(j,jproc,itype,itt)*xnorm*x1x2j*
     #                      betafact*vecmc(j,ileg,ie0sq)
                        enddo
                      enddo
                    endif
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
c
c MC subtraction term: collinear ME
c
        do jproc=loproc,maproc
          if(gfsf(jproc).lt.1.d0)then
            if(yi.gt.0.d0)then
              ytmp=1.d0
              x1t=x1soft(x1,x2,x,yi)/x
              x2t=x2soft(x1,x2,x,yi)
              xfact=1.d0/( xii*(1-yi) )
            else
              ytmp=-1.d0
              x1t=x1soft(x1,x2,x,yi)
              x2t=x2soft(x1,x2,x,yi)/x
              xfact=1.d0/( xii*(1+yi) )
            endif
            if(x1t.lt.1.and.x2t.lt.1)then
              x1x2j = x1x2jac(x1,x2,x,yi,ione)/x
              call invar_in(xm12,xm22,s,x,ytmp,cth1,cth2,str,
     #                      tk,uk,q1q,q2q,xinv)
              zg2_nlo = zgmu2_nlo() 
              ipdfscale=1
              call strfunht(x1t,x2t,sf)
              xnorm = xfact * x1x2j * zg2_nlo**2 * xphsp_ns
              xnorm = xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
              call frealht(s,x,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecmccl)
              do j=1,idrmax(jproc,jht)
                do itype=1,jtypemax(jproc,j)
                  do itt=ittmin,ittmax
                    jp=idrmap(j,jproc,itt)
                    vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #                sf(j,jproc,itype,itt)*xnorm*vecmccl(jp)
                  enddo
                enddo
              enddo
            endif
          endif
c
c MC subtraction term: soft ME
c
          if(gfsf(jproc).lt.1.d0)then
            xtmp=1.d0
            x1t=x1soft(x1,x2,x,yi)
            x2t=x2soft(x1,x2,x,yi)
            if(x1t.lt.1.and.x2t.lt.1)then
              x1x2j = x1x2jac(x1,x2,x,yi,ione)
              call invar_in(xm12,xm22,sx,xtmp,yi,cth1,cth2,str,
     #                      tk,uk,q1q,q2q,xinv)
              zg2_nlo = zgmu2_nlo() 
              ipdfscale=1
              call strfunht(x1t,x2t,sf)
              xnorm = x1x2j * zg2_nlo**2 * xphsp_s 
              xnorm = xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
              xnorm = xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
              call frealht(sx,xtmp,yi,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecmcsf)
              do j=1,idrmax(jproc,jht)
                do itype=1,jtypemax(jproc,j)
                  do itt=ittmin,ittmax
                    jp=idrmap(j,jproc,itt)
                    vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #                sf(j,jproc,itype,itt)*xnorm*vecmcsf(jp)
                  enddo
                enddo
              enddo
c
c MC subtraction term: soft-collinear ME
c
              if(gfsf(jproc).lt.1.d0)then
                if(yi.gt.0.d0)then
                  ytmp=1.d0
                  xfact=1.d0/( xii*(1-yi) )
                else
                  ytmp=-1.d0
                  xfact=1.d0/( xii*(1+yi) )
                endif
                call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                        tk,uk,q1q,q2q,xinv)
                xnorm = xfact * x1x2j * zg2_nlo**2 * xphsp_s
                xnorm = - xnorm * (1-gfsf(jproc)) * emwgt_prc(jproc)
                call frealht(sx,xtmp,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                       jproc,vecmcsc)
                do j=1,idrmax(jproc,jht)
                  do itype=1,jtypemax(jproc,j)
                    do itt=ittmin,ittmax
                      jp=idrmap(j,jproc,itt)
                      vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #                  sf(j,jproc,itype,itt)*xnorm*vecmcsc(jp)
                    enddo
                  enddo
                enddo
              endif
            endif
          endif
        enddo
      endif
c
c Counter-event :
c
      ipdfscale=1
      if(yi.gt.1-deltai) then
        ytmp = 1.d0
c The arguments of the pdf, see (A.43) in FW; x1, x2 are called z1, z2 there 
        x1t = x1soft(x1,x2,x,yi)/x
        x2t = x2soft(x1,x2,x,yi)
        if(x1t.lt.1.and.x2t.lt.1)then
          x1x2j = x1x2jac(x1,x2,x,yi,ione)/x
          call invar_in(xm12,xm22,s,x,ytmp,cth1,cth2,str,
     #                  tk,uk,q1q,q2q,xinv)
          zg2_nlo = zgmu2_nlo()
          call strfunht(x1t,x2t,sf)
          xnorm = x1x2j * zg2_nlo**2 * xphsp_ns
          xnorm = xnorm * 1.d0/xii*( - 1/(1-yi) ) 
c The following term comes from (5.7) in FKS; the 1/xi_c term has 
c the same kinematics as the collinear counter-term in the real event, 
c and is therefore conveniently included here, via f2pr
          xlmude = log(s/xmuf2h1)+log(deltai/2)
          xnormc = x1x2j * zg2_nlo**2 * xphspb /(8*pi2 * deltai)
          xnormc = xnormc/xii 
          do jproc=loproc,maproc
            call frealht(s,x,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                   jproc,vecre)
            call f2prht(s,q2q,x,x,ytmp,xlmude,jproc,veccl)
            do j=1,idrmax(jproc,jht)
              do itype=1,jtypemax(jproc,j)
                do itt=ittmin,ittmax
                  jp=idrmap(j,jproc,itt)
                  vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #              sf(j,jproc,itype,itt)*
     #              ( xnorm*vecre(jp) + xnormc*veccl(jp) )
                enddo
              enddo
            enddo
          enddo
        endif
      endif
c
      if(yi.lt.-1+deltai) then
        ytmp = -1.d0
        x1t = x1soft(x1,x2,x,yi)
        x2t = x2soft(x1,x2,x,yi)/x
        if(x1t.lt.1.and.x2t.lt.1)then
          x1x2j = x1x2jac(x1,x2,x,yi,ione)/x
          call invar_in(xm12,xm22,s,x,ytmp,cth1,cth2,str,
     #                  tk,uk,q1q,q2q,xinv)
          zg2_nlo = zgmu2_nlo()
          call strfunht(x1t,x2t,sf)
          xnorm = x1x2j * zg2_nlo**2 * xphsp_ns
          xnorm = xnorm * 1.d0/xii*( - 1/(1+yi) ) 
          xlmude = log(s/xmuf2h2)+log(deltai/2)
          xnormc = x1x2j * zg2_nlo**2 * xphspb /(8*pi2 * deltai)
          xnormc = xnormc/xii 
          do jproc=loproc,maproc
            call frealht(s,x,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                   jproc,vecre)
            call f2prht(s,q1q,x,x,ytmp,xlmude,jproc,veccl)
            do j=1,idrmax(jproc,jht)
              do itype=1,jtypemax(jproc,j)
                do itt=ittmin,ittmax
                  jp=idrmap(j,jproc,itt)
                  vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #              sf(j,jproc,itype,itt)*
     #              ( xnorm*vecre(jp) + xnormc*veccl(jp) )
                enddo
              enddo
            enddo
          enddo
        endif
      endif
c
c     Soft part of the counter-event:
c
      if(xii.lt.xicut) then
        xtmp = 1.d0
        x1t = x1soft(x1,x2,x,yi)
        x2t = x2soft(x1,x2,x,yi)
        if(x1t.lt.1.and.x2t.lt.1)then
          x1x2j = x1x2jac(x1,x2,x,yi,ione)
          call invar_in(xm12,xm22,sx,xtmp,yi,cth1,cth2,str,
     #                  tk,uk,q1q,q2q,xinv)
          zg2_nlo = zgmu2_nlo()
          call strfunht(x1t,x2t,sf)
          xnorm = x1x2j * zg2_nlo**2 * xphsp_s 
          xnorm = - xnorm * 1.d0/xii*( 1/(1-yi) + 1/(1+yi) )
          xnormb = x1x2j * zg2_nlo * xphspb /(2*xicut)
          xnormsv = x1x2j * zg2_nlo**2 * xphspb / 
     #             (8*pi2 * 2*xicut)
          do jproc=loproc,maproc
            call frealht(sx,xtmp,yi,cth2,tk,uk,q1q,q2q,xinv,
     #                   jproc,vecre)
            call fbornht(sx,q1q,jproc,vecbrn)
            call f2svht(sx,q1q,jproc,vec2sv)
            do j=1,idrmax(jproc,jht)
              xintsvc(j,jproc)=xnorm*vecre(j)
              xbornvc(j,jproc)=xnormb*vecbrn(j)
              xsvvc(j,jproc)=xnormsv*vec2sv(j)
              xcsvc(j,jproc)=0.d0
            enddo
          enddo
c
          if(yi.gt.1-deltai) then
            ytmp = 1.d0
            call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                    tk,uk,q1q,q2q,xinv)
            xnorm = x1x2j * zg2_nlo**2 * xphsp_s
            xnorm = - xnorm * 1.d0/xii*( - 1/(1-yi) ) 
            xlmude = log(sx/xmuf2h1)+log(deltai/2)
            xnormc = x1x2j * zg2_nlo**2 * xphspb /
     #               (8*pi2 * deltai)
            xnormc = -xnormc/xii
            do jproc=loproc,maproc
              call frealht(sx,xtmp,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecre)
              call f2prht(sx,q2q,x,xtmp,ytmp,xlmude,jproc,veccl)
              do j=1,idrmax(jproc,jht)
                xintsvc(j,jproc)=xintsvc(j,jproc)+
     #                           xnorm*vecre(j)
                xcsvc(j,jproc)=xnormc*veccl(j)
              enddo
            enddo
          endif
c
          if(yi.lt.-1+deltai) then
            ytmp = -1.d0
            call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                    tk,uk,q1q,q2q,xinv)
            xnorm = x1x2j * zg2_nlo**2 * xphsp_s
            xnorm = - xnorm * 1.d0/xii*( - 1/(1+yi) ) 
            xlmude = log(sx/xmuf2h2)+log(deltai/2)
            xnormc = x1x2j * zg2_nlo**2 * xphspb /
     #               (8*pi2 * deltai)
            xnormc = -xnormc/xii
            do jproc=loproc,maproc
              call frealht(sx,xtmp,ytmp,cth2,tk,uk,q1q,q2q,xinv,
     #                     jproc,vecre)
              call f2prht(sx,q1q,x,xtmp,ytmp,xlmude,jproc,veccl)
              do j=1,idrmax(jproc,jht)
                xintsvc(j,jproc)=xintsvc(j,jproc)+
     #                           xnorm*vecre(j)
                xcsvc(j,jproc)=xcsvc(j,jproc)+
     #                         xnormc*veccl(j)
              enddo
            enddo
          endif
          do jproc=loproc,maproc
            do j=1,idrmax(jproc,jht)
              do itype=1,jtypemax(jproc,j)
                do itt=ittmin,ittmax
                  jp=idrmap(j,jproc,itt)
                  vv(j,jproc,itype,itt)=vv(j,jproc,itype,itt)+
     #              sf(j,jproc,itype,itt)*
     #              ( xintsvc(jp,jproc)+xbornvc(jp,jproc)+
     #                xsvvc(jp,jproc)+xcsvc(jp,jproc) )
                enddo
              enddo
            enddo
          enddo
        endif
      endif
c
      call checkvv(xsum,dummy,iret)
      if(iret.eq.1)then
        xtmp = 1.d0
        ytmp = 1.d0
        call invar_in(xm12,xm22,sx,xtmp,ytmp,cth1,cth2,str,
     #                tk,uk,q1q,q2q,xinv)
        x1t = x1soft(x1,x2,x,yi)
        x2t = x2soft(x1,x2,x,yi)
        ycm = 0.5d0*log(x1t/x2t)
        tau=x*tau
        if(idec.eq.0)then
          ps=sx
          px=xtmp
          pyi=ytmp
          pcth1=cth1
          pcth2=cth2
        endif
c Cross section in pb (momenta are in GeV)
        xint=1.d6*hc2*xsum
        xevsign=1.d0
        if(xint.lt.0.d0)xevsign=-1.d0
c Store x1, x2 and Q2 for PDF reweighting using event file
        ux1 = x1t
        ux2 = x2t
        dummy = zgmu2_nlo()
        uq2 = xmuf2h1
      else
        xint=0.d0
        xevsign=1.d0
      endif
c
      bsfsgn=xevsign
      tot5bs_ht=abs(xint)
c
      return
      end


C Empty now, but for zeta-subtraction, to be done,
C if the number of negative weight events for single top
C happens to be too large
      function tot5bz_ht(xs,xx,xyi,xcth1,xcth2,xjac)
      implicit none
      real * 8 tot5bz_ht,xs,xx,xyi,xcth1,xcth2,xjac
      tot5bz_ht=0.d0
      return
      end
c
c
c End of cross section routines
c
c
c
c
c Begin of utility functions for zeta subtraction
c
c
      function svn(ro)
      implicit none
      real*8 svn,ro,tmp,be4,ybar,etacut
      common/cetacut/etacut
c
      tmp=0.d0
      if(ro.lt.1.d0-sqrt(etacut))then
        be4=(1-ro)**2
        ybar=sqrt(1-etacut/be4)
        tmp=-(1-ro)*ybar+sqrt(etacut)*asin(ybar)
      endif
      svn=tmp
      return
      end


      function f1fun(ro)
      implicit real * 8 (a-z)
      common/cetacut/etacut
c
      tmp=0.d0
      if(ro.lt.1.d0-sqrt(etacut))then
        be4=(1-ro)**2
        ybar=sqrt(1-etacut/be4)
        tmp=log((1+ybar)/(1-ybar))*( log(etacut/be4)
     #        -log(1-ybar**2)/2.d0-log(2.d0) )
     #     +ddilog((1+ybar)/2.d0)-ddilog((1-ybar)/2.d0) 
        tmp=tmp/4.d0
      endif
      f1fun=tmp
      return
      end


      function bdelta(x)
      implicit none
      real*8 bdelta,x,tmp,etacut
      common/cetacut/etacut
c
      tmp=0.d0
      if(x.lt.1.d0-dsqrt(etacut))tmp=sqrt(1-etacut/(1-x)**2)
      bdelta=tmp
      return
      end
c
c
c End of utility functions for zeta subtraction
c
c
c
c
c Begin of event-generation routines
c
c
      subroutine sprfin()
c This routine is called by run_spring; the entry is dummy, all the 
c parameters must be passed through common blocks
      implicit none
      integer iunit
      parameter (iunit=22)
      real*8 xone
      parameter (xone=1.d0)
      real*8 ycm,tau
      common/x1x2/ycm,tau
      integer i0,jproc0,itype0,ich0,itt0
      common/cidproc/i0,jproc0,itype0,ich0,itt0
      integer idec
      common/cidec/idec
      integer idrmax(3,3)
      common/cidrmax/idrmax
      integer iret
      real*8 ycm0
c
      call xout(iret)
      if(iret.eq.1)then
        if(i0.lt.1.or.i0.gt.idrmax(jproc0,ich0).or.ich0.ne.3)then
          write(*,*)'Fatal error in sprfin'
          stop
        endif
        if(idec.eq.0)call getspincorr()
        ycm0=ycm
        call getx1x2(tau,ycm0)
        call getmom(tau,ycm0)
        call store_events(iunit,xone)
      endif
      return
      end


      subroutine getx1x2(tau,ycm)
      implicit none
      real*8 tau,ycm,x1,x2,stau,ey
      common/cx1x2/x1,x2
c
      stau=sqrt(tau)
      ey=exp(ycm)
      x1=stau*ey
      x2=stau/ey
      return
      end


      subroutine getmom(xtau,xycm)
      implicit none
      real*8 xtau,xycm
      include 'stpcblks.h'
      real*8 pi
      parameter (pi=3.14159265358979312D0)
      integer i,imax,itype
      real*8 tau,ycm,theta,cth,sth,fk88random,sqsh,ycmnew
      real*8 x1,x2
      common/cx1x2/x1,x2
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_lb(10,4)
      common/cxmomlb/xmom_lb
      real*8 xmom_prime(10,4)
      common/cxmomprime/xmom_prime
      integer ionshell
      common/cionshell/ionshell
      integer ideconsh
      common/cideconsh/ideconsh
      integer ichkmom
      common/cichkmom/ichkmom
      integer idec
      common/cidec/idec
      integer ifk88seed
      common/cifk88seed/ifk88seed
c
      imax=5
      if(idec.eq.0)imax=8
      tau=xtau
      ycm=xycm
      call getx1x2(tau,ycm)
c perform a random rotation in the transverse plane
      theta=2*pi*fk88random(ifk88seed)
      cth=cos(theta)
      sth=sin(theta)
      do i=3,imax
        call transrot(cth,sth,xmom_cm(i,1),xmom_cm(i,2))
      enddo
      if(ichkmom.eq.0)call checkmom(xmom_cm,sh,0.d0,3,2)
c determine colour connections
      call getcolconn()
c put partons on Herwig mass shell
      if(ionshell.eq.0.and.ideconsh.eq.0)then
c keep the parton massless
        sqsh=sqrt(sh)
        xmom_lb(1,1)=0.d0
        xmom_lb(1,2)=0.d0
        xmom_lb(1,3)=x1*sqsh/2.d0
        xmom_lb(1,4)=x1*sqsh/2.d0
        xmom_lb(2,1)=0.d0
        xmom_lb(2,2)=0.d0
        xmom_lb(2,3)=-x2*sqsh/2.d0
        xmom_lb(2,4)=x2*sqsh/2.d0
        do i=3,imax
          call boost(-ycm,
     #         xmom_cm(i,1),xmom_cm(i,2),
     #         xmom_cm(i,3),xmom_cm(i,4),
     #         xmom_lb(i,1),xmom_lb(i,2),xmom_lb(i,3),xmom_lb(i,4))
        enddo
        call setxmss()
      else
c put the partons on Herwig mass shell
        call put_on_shell(ycm,ycmnew)
        do i=1,imax
          call boost(-ycmnew,
     #         xmom_prime(i,1),xmom_prime(i,2),
     #         xmom_prime(i,3),xmom_prime(i,4),
     #         xmom_lb(i,1),xmom_lb(i,2),xmom_lb(i,3),xmom_lb(i,4))
        enddo
      endif
      if(ichkmom.eq.0)then
        itype=3-idec
        call checkmom(xmom_lb,sh,-ycmnew,2,itype)
      endif
      call momnewformat()
      return
      end


      subroutine momnewformat()
c Replaces the energy with the mass in the fourth component of xmom_lb,
c to comply with the new format of the event file. Must be called as the
c last step before storing events on the temporary event files.
c If the energy is zero, the fourth component is left unchanged,
c since the LH interface uses it to distinguish between S and H events.
      implicit none
      real*8 xmom_lb(10,4)
      common/cxmomlb/xmom_lb
      real*8 xmss(1:10)
      common/procmass/xmss
      integer i
c
      do i=1,10
        if(xmom_lb(i,4).ne.0.d0)xmom_lb(i,4)=xmss(i)
      enddo
      return
      end


      subroutine setxmss()
c Fills the common block xmss. Used only if put_on_shell is not called;
c thus, set all masses equal to zero, except those of the primary top and W
      implicit none
      include 'stpcblks.h'
      integer i
      real*8 tq12
      common/ctvirt/tq12
      real*8 q12,q22
      common/cvirt/q12,q22
      real*8 xmss(1:10)
      common/procmass/xmss
      integer idec
      common/cidec/idec
c
      do i=1,10
        xmss(i)=0.d0
      enddo
      if(idec.eq.0)then
        xmss(4) = sqrt(tq12)
        xmss(5) = xm2
      elseif(idec.eq.1)then
        xmss(4) = xm1
        xmss(5) = xm2
      endif
      return
      end


      subroutine getcolconn()
c Determines colour connections
      implicit none
      include 'stpcblks.h'
      real*8 fk88random,crnd,xmomprod,s,tk,uk,q1q,q2q,q2c,w1,t12,t13,
     # t23,t1p,t2p,tp3,wg9,wg10,pr9,wg11,wg12,pr11,wg13,wg14,pr13,
     # wg15,wg16,pr15,wg17,wg18,pr17,wg19,wg20,pr19,Hqqqq_1,Hqqgg_1,
     # Hqqqq_1_nores,Hqqgg_1_nores
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer i0,jproc0,itype0,ich0,itt0
      common/cidproc/i0,jproc0,itype0,ich0,itt0
      integer i1hpro
      common/ci1hpro/i1hpro
      integer iccode
      common/ciccode/iccode
      integer idec
      common/cidec/idec
      integer ifk88seed
      common/cifk88seed/ifk88seed
c
      if(xmom_cm(3,4).eq.0.d0)then
c 2-body kinematics; may be the crossing of qq or gg events, so don't
c use jproc0 and i0 to determine the nature of the process
        if( (i1hpro.eq.402.and.itt0.eq.1) .or.
     #      (i1hpro.eq.404.and.itt0.eq.2) )then
          iccode=1
        elseif( (i1hpro.eq.405.and.itt0.eq.1) .or.
     #          (i1hpro.eq.406.and.itt0.eq.2) )then
          iccode=2
        else
          write(*,*)'Error #1 in getcolconn',
     #              i0,jproc0,itype0,ich0,itt0
          stop
        endif
      else
c 3-body kinematics; safe to use jproc0 and i0
        if( (jproc0.eq.2.and.(i0.eq.5.or.i0.eq.6.or.i0.eq.7)).or.
     #      (jproc0.eq.1.and.i0.eq.1).or.jproc0.eq.3 )then
          s=2*xmomprod(xmom_cm,1,2)
          tk=-2*xmomprod(xmom_cm,1,3)
          uk=-2*xmomprod(xmom_cm,2,3)
          q1q=-2*xmomprod(xmom_cm,1,4)+xm12
          q2q=-2*xmomprod(xmom_cm,2,5)+xm22
          q2c = xm12 + xm22 - s - uk - q2q
          w1  = xm12 - q1q + q2q - tk
c
          t12 = s
          t1p = q1q-xm12
          t13 = tk
          t2p = q2c-xm12
          t23 = uk
          tp3 = w1-xm12
c
          crnd=fk88random(ifk88seed)
          if(jproc0.eq.2.and.i0.eq.5)then
            wg9=Hqqqq_1(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)
            wg10=Hqqqq_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
            pr9=wg9/(wg9+wg10)
            if(crnd.le.pr9)then
              iccode=9
            else
              iccode=10
            endif
          elseif(jproc0.eq.2.and.i0.eq.6)then
            wg11=Hqqqq_1(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)
            wg12=Hqqqq_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
            pr11=wg11/(wg11+wg12)
            if(crnd.le.pr11)then
              iccode=11
            else
              iccode=12
            endif
          elseif(jproc0.eq.2.and.i0.eq.7)then
            wg13=Hqqqq_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
            wg14=Hqqqq_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
            pr13=wg13/(wg13+wg14)
            if(crnd.le.pr13)then
              iccode=13
            else
              iccode=14
            endif
          elseif(jproc0.eq.1.and.i0.eq.1)then
            wg15=Hqqgg_1_nores(xm22,xm12,t13,t23,t12,tp3,t1p,t2p)
            wg16=Hqqgg_1_nores(xm22,xm12,t23,t13,t12,tp3,t2p,t1p)
            pr15=wg15/(wg15+wg16)
            if(crnd.le.pr15)then
              iccode=15
            else
              iccode=16
            endif
          elseif(jproc0.eq.3.and.i0.eq.1)then
            wg17=Hqqgg_1(xm22,xm12,t12,t13,t23,t1p,t2p,tp3)
            wg18=Hqqgg_1(xm22,xm12,t13,t12,t23,t1p,tp3,t2p)
            pr17=wg17/(wg17+wg18)
            if(crnd.le.pr17)then
              iccode=17
            else
              iccode=18
            endif
          elseif(jproc0.eq.3.and.i0.eq.3)then
            wg19=Hqqgg_1(xm22,xm12,t12,t23,t13,t2p,t1p,tp3)
            wg20=Hqqgg_1(xm22,xm12,t23,t12,t13,t2p,tp3,t1p)
            pr19=wg19/(wg19+wg20)
            if(crnd.le.pr19)then
              iccode=19
            else
              iccode=20
            endif
          endif
        elseif(jproc0.eq.2.and.i0.eq.1)then
          iccode=3
        elseif(jproc0.eq.2.and.i0.eq.2)then
          iccode=4
        elseif(jproc0.eq.2.and.i0.eq.3.and.
     #         mod((itype0+mod(itype0,2))/2,2).eq.1)then
          iccode=5
        elseif(jproc0.eq.2.and.i0.eq.3.and.
     #         mod((itype0+mod(itype0,2))/2,2).eq.0)then
          iccode=6
        elseif(jproc0.eq.2.and.i0.eq.4.and.
     #         mod((itype0+mod(itype0,2))/2,2).eq.1)then
          iccode=7
        elseif(jproc0.eq.2.and.i0.eq.4.and.
     #         mod((itype0+mod(itype0,2))/2,2).eq.0)then
          iccode=8
        else
          write(*,*)'Error #2 in getcolconn',
     #              i0,jproc0,itype0,ich0,itt0
          stop
        endif
      endif
      return
      end


      function xmomprod(xmom,i,j)
      implicit none
      real*8 xmomprod,dotprod,xmom(10,4)
      integer i,j
c
      xmomprod=dotprod(xmom(i,1),xmom(i,2),xmom(i,3),xmom(i,4),
     #                 xmom(j,1),xmom(j,2),xmom(j,3),xmom(j,4))
      return
      end


      subroutine boost(y,a1,a2,a3,a4,b1,b2,b3,b4)
      implicit none
      real*8 y,a1,a2,a3,a4,b1,b2,b3,b4
c
      b1=a1
      b2=a2
      b3=a3*cosh(y)-a4*sinh(y)
      b4=a4*cosh(y)-a3*sinh(y)
      return
      end


      subroutine transrot(cth,sth,xpt1,xpt2)
      implicit none
      real*8 cth,sth,xpt1,xpt2,pt1,pt2
c
      pt1=xpt1
      pt2=xpt2
      xpt1=pt1*cth+pt2*sth
      xpt2=-pt1*sth+pt2*cth
      return
      end


      subroutine put_on_shell(ycm,ycmnew)
      implicit none
      include 'stpcblks.h'
      integer i2b,i,j,il,in,ib,ii,it
      real*8 xmss(1:10),xtmp(1:4),xk1tmp(1:4),ytmp1(1:4),ytmp2(1:4),
     #  xavg3(1:3),wvec(1:4),wvec2(1:4)
      real*8 ycm,ycmnew,pi,one,delta_thrs,shat,xkp2prime_norm2,
     #  xkp2prime_norm,xkprime_0,xsign,xnorm_3,delta,gamma,xmprime,
     #  xk1prime_norm,fakemass,xk1tmp_norm,xkprime_norm,xavgnorm,
     #  qw2,qw,xnormsq,xbwnorm,xlepnorm,tmplmass
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      parameter (delta_thrs=0.5d-3)
c Actual top mass squared if top is off shell
      real*8 tq12
      common/ctvirt/tq12
c Actual W masses squared: q12 for the W from the decay, q22 for the hard W
      real*8 q12,q22
      common/cvirt/q12,q22
      common/procmass/xmss
      real*8 xmass(-5:21)
      common/parmass/xmass
c top pole mass and width; top mass squared is stored in cmass; xmt must
c be used only in those parts of the code relevant to top decay
      real*8 xmt,twidth
      common/ctparam/xmt,twidth
c W pole mass and width
      real*8 xmw,gaw
      common/cwparam/xmw,gaw
c Lepton masses
      real*8 xlep1mass(2),xlep2mass(2)
      common/clepmass/xlep1mass,xlep2mass
c x1 and x2 are the Bjorken variables; x1 is relevant to the parton
c coming from the left
      real*8 x1,x2
      common/cx1x2/x1,x2
c xmom_cm(i,j) is the j component of the four vector of the particle # i,
c given in the partonic CM frame. j=4 is the energy. i=1,2 are the incoming
c partons, 3 is the outgoing FKS parton, 4 is the top or antitop, 5 is the H.
c When the top decays, 6=l, 7=nu, 8=b are the decay products of the top, 
c 9 and 10 are unused (this may change if Higgs will be decayed here).
c Momentum conservation is (1+2)-(3+4+5)=0 or (1+2)-(3+5+6+7+8)=0. 
c In the following, 4+5 will be referred to as the pair
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
c new momenta (put on shell) are stored here
      real*8 xmom_prime(10,4)
      common/cxmomprime/xmom_prime
c ipX is the parton code relevant to parton # X. PDG conventions are
c used: 1=d, 2=u, 3=s, 4=c, 5=b, 21=g. Only observables particles are
c listed here. Therefore, ip3=final-state light parton, ip4=top/tbar,
c ip5=H when the top doesn't decay, and ip3=final-state light parton,
c ip4=H, and (ip5,ip6,ip7)=top decay products when the top decays
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
c At least one of the variables ionshell and ideconsh must be non zero here.
c  ionshell=0 -> all partons of the hard process are massless
c  ionshell=1,2 -> all partons of the hard process are massive
c  ionshell=11,12 -> initial-state partons are massless, FKS parton is massive
c  ideconsh=0 -> all decay products are massless
c  ideconsh=2 -> all decay products are massive
c  ideconsh=12 -> all decay products are massless except leptons
      integer ionshell
      common/cionshell/ionshell
      integer ideconsh
      common/cideconsh/ideconsh
      integer ichkmom
      common/cichkmom/ichkmom
      integer idec
      common/cidec/idec
      integer iwidth
      common/ciwidth/iwidth
c 
      if(ionshell.eq.1.or.ionshell.eq.2)then
        xmss(1) = xmass(ip1)
        xmss(2) = xmass(ip2)
      else
        xmss(1) = 0.d0
        xmss(2) = 0.d0
      endif
      if(ionshell.ne.0)then
        xmss(3) = xmass(ip3)
      else
        xmss(3) = 0.d0
      endif
      if(idec.eq.0)then
        xmss(4) = sqrt(tq12)
        xmss(5) = xm2
        if(ideconsh.eq.0)then
          do i=6,8
            xmss(i) = 0.d0
          enddo
        elseif(ideconsh.eq.2)then
          xmss(6) = xlep1mass(1)
          xmss(7) = xlep2mass(1)
          xmss(8) = xmass(ip7)
          xmss(9) = -1.d10
          xmss(10) = -1.d10
        elseif(ideconsh.eq.12)then
          if( abs(ip4).ge.11.and.abs(ip4).le.16 )then
            xmss(6) = xlep1mass(1)
            xmss(7) = xlep2mass(1)
          else
            xmss(6) = 0.d0
            xmss(7) = 0.d0
          endif
          xmss(8) = 0.d0
          xmss(9) = -1.d10
          xmss(10) = -1.d10
        else
          write(*,*)'Error in put_on_shell: unknown ideconsh',ideconsh
          stop
        endif
      elseif(idec.eq.1)then
        xmss(4) = xm1
        xmss(5) = xm2
        do i=6,10
          xmss(i) = -1.d10
        enddo
      else
        write(6,*) 'Error in put_on_shell: idec=',idec
        stop
      endif
c i2b=0 --> 3-body kinematics; i2b=1 --> 2-body kinematics
      i2b = 0
      if(xmom_cm(3,4).lt.1.d-14)i2b=1
      if(ionshell.eq.0.or.ionshell.eq.1.or.ionshell.eq.11)then
c don't change the 3-momenta of partons 1,2 and 3, if possible
        do i=1,3
          do j=1,3
            xmom_prime(i,j)=xmom_cm(i,j)
          enddo
        enddo
        shat=(xmom_cm(1,4)+xmom_cm(2,4))**2
      elseif(ionshell.eq.2.or.ionshell.eq.12)then
c don't change the 3-momentum of parton 3, and shat, if possible
        do j=1,3
          xmom_prime(3,j)=xmom_cm(3,j)
        enddo
        do i=1,2
          do j=1,2
            xmom_prime(i,j)=xmom_cm(i,j)
          enddo
        enddo
        shat=(xmom_cm(1,4)+xmom_cm(2,4))**2
        call getxmss(shat,ycm,
     #               xmom_cm(1,3),xmss(1),
     #               xmom_cm(2,3),xmss(2),
     #               xmom_prime(1,3),xmom_prime(2,3))
      else
        write(*,*)'Error in put_on_shell: unknown ionshell'
        stop
      endif
      xkprime_0=0.d0
      do i=1,3
        xsign=1.d0
        if(i.eq.3)xsign=-1.d0
        if(i.eq.3.and.i2b.eq.1)then
          xmom_prime(i,4)=0.d0
        else
          call getenergy(xmom_prime(i,1),xmom_prime(i,2),
     #                   xmom_prime(i,3),xmss(i),xmom_prime(i,4))
        endif
        xkprime_0=xkprime_0+xsign*xmom_prime(i,4)
      enddo
c compute the modulus of the 3-momentum of the pair, which is equal
c to that of parton 3 in the CM frame. The energy doesn't play any role
      call fillvec(xmom_cm(3,1),xmom_cm(3,2),
     #             xmom_cm(3,3),xmom_cm(3,4),xtmp)
      xkprime_norm=xnorm_3(xtmp)
c delta is the would-be invariant mass of the pair, minus the sum
c of the masses of the top and the non-FKS parton
      delta=sqrt(xkprime_0**2-xkprime_norm**2)-xmss(4)-xmss(5)
      if(delta.lt.delta_thrs)then
c parton 3-momenta cannot be kept fixed: the total available energy
c is not sufficient; modify 3-momenta of the incoming partons
        gamma=sqrt( (xmss(4)+xmss(5)+delta_thrs)**2+xkprime_norm**2 )+
     #        xmom_prime(3,4)
        if(gamma.lt.(xmss(1)+xmss(2)))then
          write(6,*)'Error #0 in put_on_shell'
          write(6,*)gamma,xmom_prime(3,4)
          stop
        endif
        xkp2prime_norm2=( gamma**2-2*(xmss(1)**2+xmss(2)**2)+
     #                    (xmss(1)**2-xmss(2)**2)**2/gamma**2 )/4.d0
        xkp2prime_norm=sqrt(xkp2prime_norm2)
        xmom_prime(1,3)=sign(1.d0,xmom_cm(1,3))*xkp2prime_norm
        xmom_prime(1,4)=gamma/2.d0*(1+(xmss(1)**2-xmss(2)**2)/gamma**2)
        xmom_prime(2,3)=sign(1.d0,xmom_cm(2,3))*xkp2prime_norm
        xmom_prime(2,4)=gamma/2.d0*(1-(xmss(1)**2-xmss(2)**2)/gamma**2)
        xkprime_0=xmom_prime(1,4)+xmom_prime(2,4)-xmom_prime(3,4)
        shat=(xmom_prime(1,4)+xmom_prime(2,4))**2 -
     #       (xmom_prime(1,3)+xmom_prime(2,3))**2
      endif
c now the parton 3-momenta have been defined in such a way
c that the momenta of the top and the non-FKS parton can be transformed.
c xtmp is the 4-momentum of the pair. CM frame stays the same,
c so does the boost
      ycmnew=ycm
      xtmp(1)=-xtmp(1)
      xtmp(2)=-xtmp(2)
      xtmp(3)=-xtmp(3)
      xtmp(4)=xkprime_0
      xmprime=sqrt(xkprime_0**2-xkprime_norm**2)
      xk1prime_norm=xmprime**2-2*(xmss(4)**2+xmss(5)**2)+
     #              (xmss(4)**2-xmss(5)**2)**2/xmprime**2
      xk1prime_norm=sqrt(xk1prime_norm)/2.d0
      do j=1,3
        xavg3(j)=0.d0
      enddo
      do i=4,5
        xsign=1.d0
        if(i.eq.5)xsign=-1.d0
        call fillvec(xmom_cm(i,1),xmom_cm(i,2),
     #               xmom_cm(i,3),xmom_cm(i,4),ytmp1)
        call xhwulof(xtmp,xmprime,
     #               ytmp1,xmss(i),
     #               xk1tmp,fakemass)
        if(abs(fakemass-xmss(i)).gt.1.d-4)then
          write(6,*)'Error #1 in put_on_shell'
          write(6,*)i,xmss(i),fakemass
          stop
        endif
        xk1tmp_norm=xnorm_3(xk1tmp)
c xavg is the direction along which the top and non-FKS parton momenta are 
c placed in the new pair rest frame. It is arbitrarily defined by averaging 
c (hence the 1/2 in the definition) the directions of the original top and
c non-FKS parton momenta. It may not have modulus 1, so normalize it
        do j=1,3
          xavg3(j)=xavg3(j)+xsign*xk1tmp(j)/(2*xk1tmp_norm)
        enddo
      enddo
      xavgnorm=sqrt(xavg3(1)**2+xavg3(2)**2+xavg3(3)**2)
      do j=1,3
        xavg3(j)=xavg3(j)/xavgnorm
      enddo
      do i=4,5
        xsign=1.d0
        if(i.eq.5)xsign=-1.d0
        do j=1,3
          xk1tmp(j)=xsign*xk1prime_norm*xavg3(j)
        enddo
        xk1tmp(4)=xmprime/2.d0*
     #            (1+xsign*(xmss(4)**2-xmss(5)**2)/xmprime**2)
        call xhwulob(xtmp,xmprime,
     #               xk1tmp,xmss(i),
     #               ytmp2,fakemass)
        if(abs(fakemass-xmss(i)).gt.1.d-4)then
          write(6,*)'Error #2 in put_on_shell'
          write(6,*)i,xmss(i),fakemass
          stop
        endif
        call getvec(ytmp2,xmom_prime(i,1),xmom_prime(i,2),
     #                    xmom_prime(i,3),xmom_prime(i,4))
      enddo
c 
      if(idec.eq.0)then
c Top decay; adapted from QQ production
        it=4
        il=6
        in=7
        ib=8
        call fillvec(xmom_prime(it,1),xmom_prime(it,2),
     #               xmom_prime(it,3),xmom_prime(it,4),xtmp)
c First deal with the Wb pair; define W momentum, and compute W mass
c (when iwidth=1, W is off shell)
        call vecsum(xmom_cm(il,1),xmom_cm(il,2),
     #              xmom_cm(il,3),xmom_cm(il,4),one,
     #              xmom_cm(in,1),xmom_cm(in,2),
     #              xmom_cm(in,3),xmom_cm(in,4),one,wvec)
        qw2=xnormsq(wvec)
        qw=sqrt(qw2)
        if( ichkmom.eq.0 .and. 
     #      ( (iwidth.eq.0.and.abs(qw/xmw-1.d0).gt.1.d-4) .or.
     #        (iwidth.eq.1.and.it.eq.4.and.
     #         abs(qw/sqrt(q12)-1.d0).gt.1.d-4) ) )then
          write(6,*)'Error #3 in put_on_shell'
          write(6,*)qw,it,il,in
          stop
        endif
        if( ichkmom.eq.0 .and. iwidth.eq.1 .and.
     #      qw.gt.xmss(it) )then
          write(6,*)'Error #4 in put_on_shell'
          write(6,*)qw,it,il,in
          stop
        endif
        xbwnorm=xmss(it)**2-2*(xmss(ib)**2+qw2)+
     #          (xmss(ib)**2-qw2)**2/xmss(it)**2
        xbwnorm=sqrt(xbwnorm)/2.d0
        do j=1,3
          xavg3(j)=0.d0
        enddo
        xsign=1.d0
        call xhwulof(xtmp,xmss(it),wvec,qw,xk1tmp,fakemass)
        xk1tmp_norm=xnorm_3(xk1tmp)
        do j=1,3
          xavg3(j)=xavg3(j)+xsign*xk1tmp(j)/(2*xk1tmp_norm)
        enddo
        xsign=-1.d0
        call fillvec(xmom_cm(ib,1),xmom_cm(ib,2),
     #               xmom_cm(ib,3),xmom_cm(ib,4),ytmp1)
        call xhwulof(xtmp,xmss(it),ytmp1,xmss(ib),xk1tmp,fakemass)
        xk1tmp_norm=xnorm_3(xk1tmp)
        do j=1,3
          xavg3(j)=xavg3(j)+xsign*xk1tmp(j)/(2*xk1tmp_norm)
        enddo
        xavgnorm=sqrt(xavg3(1)**2+xavg3(2)**2+xavg3(3)**2)
        do j=1,3
          xavg3(j)=xavg3(j)/xavgnorm
        enddo
        xsign=1.d0
        do j=1,3
          xk1tmp(j)=xsign*xbwnorm*xavg3(j)
        enddo
        xk1tmp(4)=xmss(it)/2.d0*
     #            (1+xsign*(qw2-xmss(ib)**2)/xmss(it)**2)
        call xhwulob(xtmp,xmss(it),xk1tmp,qw,wvec2,fakemass)
        xsign=-1.d0
        do j=1,3
          xk1tmp(j)=xsign*xbwnorm*xavg3(j)
        enddo
        xk1tmp(4)=xmss(it)/2.d0*
     #            (1+xsign*(qw2-xmss(ib)**2)/xmss(it)**2)
        call xhwulob(xtmp,xmss(it),xk1tmp,xmss(ib),ytmp2,fakemass)
        call getvec(ytmp2,xmom_prime(ib,1),xmom_prime(ib,2),
     #                    xmom_prime(ib,3),xmom_prime(ib,4))
c Next deal with the lepton pair; W has momentum wvec2
        xlepnorm=qw2-2*(xmss(il)**2+xmss(in)**2)+
     #           (xmss(il)**2-xmss(in)**2)**2/qw2
        xlepnorm=sqrt(xlepnorm)/2.d0
        do j=1,3
          xavg3(j)=0.d0
        enddo
        do i=1,2
          if(i.eq.1)then
            xsign=1.d0
            ii=il
          else
            xsign=-1.d0
            ii=in
          endif
          tmplmass=xmss(ii)
          call fillvec(xmom_cm(ii,1),xmom_cm(ii,2),
     #                 xmom_cm(ii,3),xmom_cm(ii,4),ytmp1)
          call xhwulof(wvec2,qw,ytmp1,tmplmass,xk1tmp,fakemass)
          xk1tmp_norm=xnorm_3(xk1tmp)
          do j=1,3
            xavg3(j)=xavg3(j)+xsign*xk1tmp(j)/(2*xk1tmp_norm)
          enddo
        enddo
        xavgnorm=sqrt(xavg3(1)**2+xavg3(2)**2+xavg3(3)**2)
        do j=1,3
          xavg3(j)=xavg3(j)/xavgnorm
        enddo
        do i=1,2
          if(i.eq.1)then
            xsign=1.d0
            ii=il
          else
            xsign=-1.d0
            ii=in
          endif
          tmplmass=xmss(ii)
          do j=1,3
            xk1tmp(j)=xsign*xlepnorm*xavg3(j)
          enddo
          xk1tmp(4)=qw/2.d0*
     #      (1+xsign*(xmss(il)**2-xmss(in)**2)/qw2)
          call xhwulob(wvec2,qw,xk1tmp,tmplmass,ytmp2,fakemass)
          call getvec(ytmp2,xmom_prime(ii,1),xmom_prime(ii,2),
     #                      xmom_prime(ii,3),xmom_prime(ii,4))
        enddo
c Insert here H decay if need be. Here set momenta consistently with
c what done in genhdmom
        do j=1,4
          xmom_prime(9,j)=xmom_prime(5,j)
          xmom_prime(10,j)=0.d0
        enddo
      else
        do i=6,10
          do j=1,4
            xmom_prime(i,j)=0.d0
          enddo
        enddo
      endif
c
      if(ichkmom.eq.0)then
        if(idec.eq.0)then
          call checktdec2(xmom_prime,4,6,7,8)
          call checkmom(xmom_prime,shat,0.d0,4,1)
        else
          call checkmom(xmom_prime,shat,0.d0,4,2)
        endif
        if(xmass(1).eq.0.and.xmass(2).eq.0.and.xmass(3).eq.0.and.
     #     xmass(4).eq.0.and.xmass(5).eq.0.and.xmass(21).eq.0.and.
     #     xlep1mass(1).eq.0.and.xlep2mass(1).eq.0)then
          call checkonsh(1)
        else
          call checkonsh(2)
        endif
      endif
      return
      end


      subroutine getxmss(shat,ycm,p13cm,xm1,p23cm,xm2,p13,p23)
c After putting the momenta on shell, the two incoming partons may
c travel in the same direction. This routine prevents this to happen,
c redefining Herwig masses if necessary
      implicit none
      real*8 shat,ycm,p13cm,xm1,p23cm,xm2,p13,p23
      real*8 tiny,fact,sqs,xm1s,xm2s,xkp2prime_norm2,xkp2prime_norm,
     #  ytmp,e1,e2,p13p,p23p,s1p,s2p,xif,sol
      integer iflag,idone,ileg
      integer ionshell
      common/cionshell/ionshell
      parameter (fact=0.98d0)
      parameter (tiny=1.d-6)
c
      sqs=sqrt(shat)
      xm1s=xm1
      xm2s=xm2
      ytmp=-ycm
      idone=0
 100  continue
      xkp2prime_norm2=( shat-2*(xm1**2+xm2**2)+
     #                  (xm1**2-xm2**2)**2/shat )/4.d0
      xkp2prime_norm=sqrt(xkp2prime_norm2)
      if(sign(1.d0,p13cm).ne.1.d0.or.sign(1.d0,p23cm).ne.-1.d0)then
        write(*,*)'Error # 0 in getxmss'
        stop
      endif
      p13=xkp2prime_norm
      p23=-xkp2prime_norm
      e1=sqrt(shat)/2.d0*(1+(xm1**2-xm2**2)/shat)
      e2=sqrt(shat)/2.d0*(1-(xm1**2-xm2**2)/shat)
      p13p=p13*cosh(ytmp)-e1*sinh(ytmp)
      p23p=p23*cosh(ytmp)-e2*sinh(ytmp)
      s1p=sign(1.d0,p13p)
      s2p=sign(1.d0,p23p)
      iflag=0
      if(s1p.eq.1.d0 .and. s2p.eq.-1.d0)then
        iflag=1
      elseif(s1p.eq.-1.d0 .and. s2p.eq.-1.d0)then
        if(ytmp.lt.0.d0)then
          write(*,*)'Wrong y sign, # 1'
          stop
        endif
        ileg=1
        xif=xm2**2/shat
      elseif(s1p.eq.1.d0 .and. s2p.eq.1.d0)then
        if(ytmp.gt.0.d0)then
          write(*,*)'Wrong y sign, # 2'
          stop
        endif
        ileg=2
        xif=xm1**2/shat
      else
        write(*,*)'Error # 1 in getxmss'
        stop
      endif
      if(iflag.eq.0)then
        if(ionshell.eq.0.or.ionshell.eq.11.or.ionshell.eq.12)then
          write(*,*)'Error # 3 in getxmss',xm1,xm2,p13cm,p23cm,p13,p23
          stop
        endif
        sol=xif+cosh(2*ytmp)-
     #      sqrt(2.d0)*cosh(ytmp)*sqrt(cosh(2*ytmp)-1+2*xif)
        if(sol.le.0.d0.or.idone.eq.1)then
c The procedure failed; pass the massless event to Herwig, and let Herwig
c deal with it
          xm1=0.d0
          xm2=0.d0
          p13=sqs/2.d0
          p23=-sqs/2.d0
          return
        endif
        if(ileg.eq.1)then
          xm1=fact*sqrt(sol*shat)
          if(xm1.gt.xm1s)then
            write(*,*)'Mass # 1 too large in getxmss'
            stop
          endif
        elseif(ileg.eq.2)then
          xm2=fact*sqrt(sol*shat)
          if(xm2.gt.xm2s)then
            write(*,*)'Mass # 2 too large in getxmss'
            stop
          endif
        else
          write(*,*)'Error # 2 in getxmss'
          stop
        endif
        idone=1
        goto 100
      endif
      return
      end


      subroutine fillvec(p1,p2,p3,p4,ytmp)
      implicit none
      real*8 p1,p2,p3,p4,ytmp(1:4)
c
      ytmp(1)=p1
      ytmp(2)=p2
      ytmp(3)=p3
      ytmp(4)=p4
      return
      end


      subroutine getvec(ytmp,p1,p2,p3,p4)
      implicit none
      real*8 ytmp(1:4),p1,p2,p3,p4
c
      p1=ytmp(1)
      p2=ytmp(2)
      p3=ytmp(3)
      p4=ytmp(4)
      return
      end

c-----------------------------------------------------------------------
      subroutine xhwulob(ps,ps5,pi,pi5,pf,pf5)
c     transforms pi (given in rest frame of ps) into pf (in lab)
c     n.b. p(1,2,3,4,5) = (px,py,pz,e,m)
c-----------------------------------------------------------------------
      real*8 pf4,fn,ps(4),ps5,pi(4),pi5,pf(4),pf5
      if (ps(4).eq.ps5) then
        pf(1)= pi(1)
        pf(2)= pi(2)
        pf(3)= pi(3)
        pf(4)= pi(4)
      else
        pf4  = (pi(1)*ps(1)+pi(2)*ps(2)
     &         +pi(3)*ps(3)+pi(4)*ps(4))/ps5
        fn   = (pf4+pi(4)) / (ps(4)+ps5)
        pf(1)= pi(1) + fn*ps(1)
        pf(2)= pi(2) + fn*ps(2)
        pf(3)= pi(3) + fn*ps(3)
        pf(4)= pf4
      end if
      pf5= pi5
      end


c-----------------------------------------------------------------------
      subroutine xhwulof(ps,ps5,pi,pi5,pf,pf5)
c     transforms pi (given in lab) into pf (in rest frame of ps)
c     n.b. p(1,2,3,4,5) = (px,py,pz,e,m)
c-----------------------------------------------------------------------
      real*8 pf4,fn,ps(4),ps5,pi(4),pi5,pf(4),pf5
      if (ps(4).eq.ps5) then
        pf(1)= pi(1)
        pf(2)= pi(2)
        pf(3)= pi(3)
        pf(4)= pi(4)
      else
        pf4  = (pi(4)*ps(4)-pi(3)*ps(3)
     &         -pi(2)*ps(2)-pi(1)*ps(1))/ps5
        fn   = (pf4+pi(4)) / (ps(4)+ps5)
        pf(1)= pi(1) - fn*ps(1)
        pf(2)= pi(2) - fn*ps(2)
        pf(3)= pi(3) - fn*ps(3)
        pf(4)= pf4
      end if
      pf5= pi5
      end


C-----------------------------------------------------------------------
      SUBROUTINE HWULB4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(1)*PS(1)+PI(2)*PS(2)
     &         +PI(3)*PS(3)+PI(4)*PS(4))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) + FN*PS(1)
        PF(2)= PI(2) + FN*PS(2)
        PF(3)= PI(3) + FN*PS(3)
        PF(4)= PF4
      END IF
      END


C-----------------------------------------------------------------------
      SUBROUTINE HWULF4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN LAB) INTO PF (IN REST FRAME OF PS)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(4)*PS(4)-PI(3)*PS(3)
     &         -PI(2)*PS(2)-PI(1)*PS(1))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) - FN*PS(1)
        PF(2)= PI(2) - FN*PS(2)
        PF(3)= PI(3) - FN*PS(3)
        PF(4)= PF4
      END IF
      END


      subroutine getenergy(p1,p2,p3,xm,en)
      implicit none
      real*8 p1,p2,p3,xm,en
c
      en=sqrt(p1**2+p2**2+p3**2+xm**2)
      return
      end


      function dotprod(p1,p2,p3,p4,q1,q2,q3,q4)
      implicit none
      real*8 dotprod,p1,p2,p3,p4,q1,q2,q3,q4
c
      dotprod=p4*q4-p1*q1-p2*q2-p3*q3
      return
      end


      function xnormsq(p)
c Computes p.p, assuming the energy is the fourth component
      implicit none
      real*8 xnormsq,p(1:4),dotprod
c
      xnormsq=dotprod(p(1),p(2),p(3),p(4),p(1),p(2),p(3),p(4))
      return
      end


      function xnorm_3(p)
c Evaluates the norm of the spatial component of a four-momentum
c The result is positive by definition, regardless of the 4-metric
      implicit none
      real*8 xnorm_3,p(1:4),tmp
c
      tmp=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
      xnorm_3=sqrt(tmp)
      return
      end


      subroutine vecsum(p1,p2,p3,p4,pfact,q1,q2,q3,q4,qfact,r)
c Weighted sum of the four-vectors p and q. The result is r
      implicit none
      real*8 p1,p2,p3,p4,pfact,q1,q2,q3,q4,qfact,r(1:4)
c
      r(1)=pfact*p1+qfact*q1
      r(2)=pfact*p2+qfact*q2
      r(3)=pfact*p3+qfact*q3
      r(4)=pfact*p4+qfact*q4
      return
      end


      subroutine checkonsh(itype)
c Checks that put_on_shell is harmless if masses are zero (itype=1),
c or computes (itype=2) the average of the shifts due to the masses
      real*8 tiny
      parameter (tiny=1.d-4)
      integer itype
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_prime(10,4)
      common/cxmomprime/xmom_prime
      real*8 xmomshifts(4)
      common/cshifts/xmomshifts
      integer i,j,imax,iflag
      integer idec
      common/cidec/idec
c
      if(itype.ne.1.and.itype.ne.2)then
        write(*,*)'Unknown option in checkonsh'
        stop
      endif
      iflag=1
      imax=5
      if(idec.eq.0)imax=8
      do i=1,imax
        do j=1,4
          if(abs(xmom_cm(i,j)).lt.1.d0)then
            xtmp=xmom_cm(i,j)-xmom_prime(i,j)
          else
            xtmp=(xmom_cm(i,j)-xmom_prime(i,j))/xmom_cm(i,j)
          endif
          if(abs(xtmp).gt.tiny)iflag=0
          xmomshifts(j)=xmomshifts(j)+abs(xtmp)
        enddo
      enddo
      if(iflag.eq.0.and.itype.eq.1)then
        write(*,*)'Error in checkonsh'
        write(*,*)'  '
        write(*,*)'xmom_cm:'
        do i=1,imax
          write(*,'(4(d14.8,1x))') (xmom_cm(i,j),j=1,4)
        enddo
        write(*,*)'  '
        write(*,*)'xmom_prime:'
        do i=1,imax
          write(*,'(4(d14.8,1x))') (xmom_prime(i,j),j=1,4)
        enddo
        stop
      endif
      return
      end


      subroutine xout(iret)
c This routine is called by sprfin; it determines, on statistical
c basis, which partonic process has been generated.
c It also counts the number of unlike sign events (iwrong), and the number
c of these events (iwrong1) for which the relative difference between
c unlike signs exceeds 5%. If all the entries of vv are equal to zero,
c iret is set equal to 0 (by checkvv), and no operation is performed
      implicit none
      integer iret,iretvv,iretvvs,iproc,iproclo,iprocma,i,itype,iz,io,
     #  itt,iwh,iflag,ihpro,i1,i2,i3,i4,i5,i1hproo,ip1o,ip2o,ip3o
      real*8 wwx(7,3,24,2),xsum,xsumabs,xsumvvs,xsumabsvvs,xstsign,
     #  xg,wh,rmax,fk88random
      parameter (iz=0)
      parameter (io=1)
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      integer ifuntype
      common/cifuntype/ifuntype
      integer ichmin,ichmax
      common/cichrange/ichmin,ichmax
      integer ittmin,ittmax
      common/cittrange/ittmin,ittmax
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer idrmax(3,3)
      common/cidrmax/idrmax
      real*8 vv(7,3,24,2)
      common/cvv/vv
      real*8 vvs(7,3,24,2)
      common/cvvs/vvs
      real*8 emsca_prc(3)
      common/cemscaprc/emsca_prc
      real*8 emsca
      common/cemsca/emsca
      integer iwrong,iwrong1
      common/ciwrong/iwrong,iwrong1
      integer i0,jproc0,itype0,ich0,itt0
      common/cidproc/i0,jproc0,itype0,ich0,itt0
      integer ivbhpro(7,3,24,2)
      common/civbhpro/ivbhpro
      integer idp1(7,3,24,2),idp2(7,3,24,2),idp3(7,3,24,2)
      common/cidpart/idp1,idp2,idp3
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer ichkpid
      common/cichkpid/ichkpid
      integer idec
      common/cidec/idec
      integer iallzero
      parameter (iallzero=1)
      integer jht
      parameter (jht=3)
c
      if(ichmin.ne.3.or.ichmax.ne.3)then
        write(*,*)'This version of xout must be called for Ht channel'
        stop
      endif
      i0=0
      jproc0=0
      itype0=0
      ich0=jht
      itt0=0
      iret=0
      call checkvv(xsum,xsumabs,iretvv)
      call checkvvs(xsumvvs,xsumabsvvs,iretvvs)
      if(iretvv.eq.0.and.iretvvs.eq.1)then
        write(6,*)'Fatal error in xout:',iretvv,iretvvs
        stop
      endif
      if(iretvv.eq.1)then
        iret=iretvv
        if(ifuntype.eq.1)then
          iproclo=loproc
          iprocma=maproc
        elseif(ifuntype.eq.2)then
          iproclo=loproc
          iprocma=maproc
        else
          write(*,*)'Fatal error in xout: ifuntype=',ifuntype
          stop
        endif
        if(iretvvs.eq.1)then
          xsum=xsumvvs
          xsumabs=xsumabsvvs
          do iproc=iproclo,iprocma
            do i=1,idrmax(iproc,jht)
              do itype=1,jtypemax(iproc,i)
                do itt=ittmin,ittmax
                  wwx(i,iproc,itype,itt)=
     #              vvs(i,iproc,itype,itt)
                enddo
              enddo
            enddo
          enddo
        else
          do iproc=iproclo,iprocma
            do i=1,idrmax(iproc,jht)
              do itype=1,jtypemax(iproc,i)
                do itt=ittmin,ittmax
                  wwx(i,iproc,itype,itt)=
     #              vv(i,iproc,itype,itt)
                enddo
              enddo
            enddo
          enddo
        endif
        xstsign=sign(1.d0,xsum)
        xg=fk88random(ifk88seed)
        wh=0.d0
        iwh=0
        iflag=0
        rmax=0.d0
        do iproc=iproclo,iprocma
          do i=1,idrmax(iproc,jht)
            do itype=1,jtypemax(iproc,i)
              do itt=ittmin,ittmax
                if(iwh.eq.0)then
                  wh=wh+abs(wwx(i,iproc,itype,itt))/xsumabs
                  if(wh.gt.xg)then
                    i0=i
                    jproc0=iproc
                    itype0=itype
                    itt0=itt
                    iwh=1
                  endif
                endif
                if(wwx(i,iproc,itype,itt).ne.0.d0)then
                  if(xstsign.ne.
     #               sign(1.d0,wwx(i,iproc,itype,itt)))then
                    if(iflag.eq.0)then
                      iwrong=iwrong+1
                      iflag=1
                    endif
                    rmax=max(rmax,abs(wwx(i,iproc,itype,itt)))
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(iflag.eq.1)then
          if(xsum.ne.0.d0)rmax=rmax/xsum
          if(rmax.gt.0.05d0.or.xsum.eq.0.d0)iwrong1=iwrong1+1
        endif
        if(i0.eq.0.or.jproc0.eq.0.or.itype0.eq.0.or.
     #     ich0.eq.0.or.itt0.eq.0)then
          write(*,*)'Fatal error in xout',i0,jproc0,itype0,ich0,itt0
          stop
        endif
        ihpro=ivbhpro(i0,jproc0,itype0,itt0)
        i1=idp1(i0,jproc0,itype0,itt0)
        i2=idp2(i0,jproc0,itype0,itt0)
        i3=idp3(i0,jproc0,itype0,itt0)
        if(abs(i3).gt.100)call ckmunwgt(i3)
        if(itt0.eq.1)then
          i4=6
          i5=-37
        elseif(itt0.eq.2)then
          i4=-6
          i5=37
        else
          write(*,*)'Fatal error: wrong top identity',itt0
          stop
        endif
        if(ichkpid.eq.0)
     #    call parcheckfin(ihpro,i1,i2,i3,i4,i5,iallzero,iz,
     #                     i0,jproc0,itype0,ich0,itt0)
        call parcrossing(i0,jproc0,ihpro,i1,i2,i3,i4,
     #                   i1hproo,ip1o,ip2o,ip3o)
        i1hpro=i1hproo
        ip1=ip1o
        ip2=ip2o
        ip3=ip3o
        ip4=i4
        ip5=i5
        if(ichkpid.eq.0)
     #    call parcheckfin(i1hpro,ip1,ip2,ip3,ip4,ip5,iallzero,io,
     #                     i0,jproc0,itype0,ich0,itt0)
c The top decays. In such case, ip3 is a final-state light parton, 
c ip4 is the Higgs, and ip5,ip6,ip7 are top decay products identities. The
c routine getpdecids determines the identities of the decay products, 
c and assigns them the relevant masses (to be used by put_on_shell)
        if(idec.eq.0)then
          ip4=ip5
          call getpdecids()
        endif
        emsca=emsca_prc(jproc0)
      endif
      return
      end


      subroutine getpdecids()
c Determine the identities of top decay products. Proceed assuming
c top (ie not antitop) production; charge conjugate at the end if 
c an antitop was actually generated
      implicit none
      integer iemu,iemutau,ipone,imone
      parameter (iemu=1113)
      parameter (iemutau=111315)
      parameter (ipone=1)
      parameter (imone=-1)
      real*8 xg,fk88random
c frac12 is the fraction of decays W->e+mu/W->e+mu+all quarks
c frac123 is the fraction of decays W->e+mu+tau/W->e+mu+tau+all quarks
      real*8 frac12,frac123
      common/cfracs123/frac12,frac123
      real*8 xlep1mass(2),xlep2mass(2)
      common/clepmass/xlep1mass,xlep2mass
      real*8 xmass(-5:21)
      common/parmass/xmass
      real*8 pdglepmss(11:16)
      common/cpdglepmss/pdglepmss
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer inonbtop
      common/cinonbtop/inonbtop
      integer il1hw,il2hw
      common/cilhw/il1hw,il2hw
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
      integer ip1s,ip2s,ip3s,ip4s,ip5s,ip6s,ip7s
      common/ci1parts/ip1s,ip2s,ip3s,ip4s,ip5s,ip6s,ip7s
      integer i0,jproc0,itype0,ich0,itt0
      common/cidproc/i0,jproc0,itype0,ich0,itt0
c
      if(il1hw.eq.7)then
        write(*,*)'Error #1 in getpdecids()',il1hw
        stop
      endif
      if(inonbtop.eq.0)then
c t->Wb
        ip7=ip7s
      elseif(inonbtop.eq.1)then
c t->W+any down-type quark
        call dectopuwgt(ipone,ip7)
      else
        write(*,*)'Unknown option in getpdecids()',inonbtop
        stop
      endif
c W+ decay (for top decay)
      if(il1hw.eq.0)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac123)then
          call declepuwgt(-iemutau,ip5,ip6)
        else
          call decqrkuwgt(ipone,ip5,ip6)
        endif
      elseif(il1hw.ge.1.and.il1hw.le.3)then
        ip5=ip5s
        ip6=ip6s
      elseif(il1hw.eq.4)then
        call declepuwgt(-iemu,ip5,ip6)
      elseif(il1hw.eq.5)then
        call decqrkuwgt(ipone,ip5,ip6)
      elseif(il1hw.eq.6)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac12)then
          call declepuwgt(-iemu,ip5,ip6)
        else
          call decqrkuwgt(ipone,ip5,ip6)
        endif
      else
        write(*,*)'Error #2 in getpdecids()',il1hw
        stop
      endif
      if(abs(ip5).le.5.and.abs(ip6).le.5)then
        xlep1mass(1)=xmass(ip5)
        xlep2mass(1)=xmass(ip6)
      elseif( abs(ip5).ge.11.and.abs(ip5).le.16 .and.
     #        abs(ip6).ge.11.and.abs(ip6).le.16 )then
        xlep1mass(1)=pdglepmss(abs(ip5))
        xlep2mass(1)=pdglepmss(abs(ip6))
      else
        write(*,*)'Error #4 in getpdecids()',ip5,ip6
        stop
      endif
c A tbar was produced: charge-conjugate
      if(itt0.eq.2)then
        ip5=-ip5
        ip6=-ip6
        ip7=-ip7
      elseif(itt0.ne.1)then
        write(*,*)'Error #6 in getpdecids()',itt0
        stop
      endif
c
      return
      end


      subroutine declepuwgt(ip4s,ip4,ip5)
c Determines on statistical basis the identity of the leptons emerging
c from W decays
      implicit none
      integer ip4s,ip4,ip5,ii
      real*8 xg,fk88random
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer ichlw0(1:3)
      data ichlw0/11,13,15/
      integer ineuw0(1:3)
      data ineuw0/12,14,16/
c
      xg=fk88random(ifk88seed)
      if(abs(ip4s).eq.1113)then
        ii=1+int(2*xg)
        ip4=sign(1,ip4s)*ichlw0(ii)
        ip5=-sign(1,ip4s)*ineuw0(ii)
      elseif(abs(ip4s).eq.111315)then
        ii=1+int(3*xg)
        ip4=sign(1,ip4s)*ichlw0(ii)
        ip5=-sign(1,ip4s)*ineuw0(ii)
      else
        write(*,*)'Error in declepuwgt',ip4s
        stop
      endif
      return
      end


      subroutine decqrkuwgt(ipmone,ip4,ip5)
c Determines on statistical basis the identity of the quarks emerging
c from W+ (ipmone=1) or W- (ipmone=-1) decays
      implicit none
      integer ipmone,ip4,ip5,iwh,iresu,iresd,ii,jj
      real*8 xg,fk88random,xden,wh
      real*8 ckm2(1:6,1:6)
      common/cckm2/ckm2
      real*8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      integer imapp(0:5)
      common/cimapp/imapp
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer imapd(3)
      data imapd/2,3,5/
      integer imapu(3)
      data imapu/1,4,6/
c 
c All decays: (u+any downbar)+(c+any downbar)
      xden=ruckm+rcckm
      xg=fk88random(ifk88seed)
      wh=0.d0
      iwh=0
      iresu=0
      iresd=0
      do ii=1,2
        do jj=1,3
          if(iwh.eq.0)then
            wh=wh+ckm2(imapu(ii),imapd(jj))/xden
            if(wh.gt.xg)then
              iresu=imapu(ii)
              iresd=imapd(jj)
              iwh=1
            endif
          endif
        enddo
      enddo
c
      if(iresu.eq.0.or.iresd.eq.0)then
        write(*,*)'Error #1 in decqrkuwgt:',ipmone,ip4,ip5
        stop
      else
        if(ipmone.eq.1)then
          ip4=-imapp(iresd)
          ip5=imapp(iresu)
        elseif(ipmone.eq.-1)then
          ip4=imapp(iresd)
          ip5=-imapp(iresu)
        else
          write(*,*)'Error #2 in decqrkuwgt:',ipmone,ip4,ip5
          stop
        endif
      endif
      return
      end


      subroutine dectopuwgt(ipmone,ip6)
c Determines on statistical basis the identity of the down-type quark
c emerging from top (ipmone=1) or antitop (ipmone=-1) decays
      implicit none
      integer ipmone,ip6,iwh,iresd,jj
      real*8 xden,xg,fk88random,wh
      real*8 ckm2(1:6,1:6)
      common/cckm2/ckm2
      real*8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      integer imapp(0:5)
      common/cimapp/imapp
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer imapd(3)
      data imapd/2,3,5/
c
      xden=rtckm
      xg=fk88random(ifk88seed)
      wh=0.d0
      iwh=0
      iresd=0
      do jj=1,3
        if(iwh.eq.0)then
          wh=wh+ckm2(6,imapd(jj))/xden
          if(wh.gt.xg)then
            iresd=imapd(jj)
            iwh=1
          endif
        endif
      enddo
c
      if(iresd.eq.0)then
        write(*,*)'Error #1 in dectopuwgt:',ipmone,ip6
        stop
      else
        if(ipmone.eq.1)then
          ip6=imapp(iresd)
        elseif(ipmone.eq.-1)then
          ip6=-imapp(iresd)
        else
          write(*,*)'Error #2 in dectopuwgt:',ipmone,ip6
          stop
        endif
      endif
      return
      end


      subroutine parcrossing(i0,jproc0,ihpro,i1,i2,i3,i4,
     #                       i1hproo,ip1o,ip2o,ip3o)
      implicit none
      integer i0,jproc0,ihpro,i1,i2,i3,i4,i1hproo,ip1o,ip2o,ip3o,
     # ihpropl(401:409),ihpromn(401:409)
      real*8 xg,fk88random
      integer ifuntype
      common/cifuntype/ifuntype
      integer ifk88seed
      common/cifk88seed/ifk88seed
      data ihpropl/406,0,405,0,0,0,0,405,406/
      data ihpromn/402,0,404,0,0,0,0,402,404/
c
      i1hproo=0
      ip1o=0
      ip2o=0
      ip3o=0
      if( ifuntype.eq.1 .or. 
     #    (ifuntype.eq.2.and.jproc0.eq.3) )then
        i1hproo=ihpro
        ip1o=i1
        ip2o=i2
        ip3o=i3
      elseif(ifuntype.eq.2.and.(jproc0.eq.1.or.jproc0.eq.2))then
        if( (jproc0.eq.1.and.(ihpro.ne.407.or.i0.ne.1)) .or.
     #      (jproc0.eq.2.and.(i0.eq.1.or.i0.eq.2)) .or.
     #      ihpro.eq.402 .or. ihpro.eq.404 .or.
     #      ihpro.eq.405 .or. ihpro.eq.406 )then
          write(*,*)'Error #1 in parcrossing:',
     #              i0,jproc0,ihpro,i1,i2,i3,i4
          stop
        endif
        xg=fk88random(ifk88seed)
        if(jproc0.eq.1)then
          if(xg.le.0.5d0)then
            ip1o=i1
            ip2o=-i3
            ip3o=21
            if(i4.eq.6)then
              i1hproo=405
            elseif(i4.eq.-6)then
              i1hproo=406
            endif
          else
            ip1o=-i3
            ip2o=i2
            ip3o=21
            if(i4.eq.6)then
              i1hproo=402
            elseif(i4.eq.-6)then
              i1hproo=404
            endif
          endif
        elseif( jproc0.eq.2 .and.
     #          (i0.eq.3.or.i0.eq.5.or.(i0.eq.7.and.xg.le.0.5d0)) )then
          ip1o=i1
          ip2o=21
          ip3o=21
          i1hproo=ihpromn(ihpro)
        elseif( jproc0.eq.2 .and.
     #          (i0.eq.4.or.i0.eq.6.or.(i0.eq.7.and.xg.gt.0.5d0)) )then
          ip1o=21
          ip2o=i2
          ip3o=21
          i1hproo=ihpropl(ihpro)
        else
          write(*,*)'Error #2 in parcrossing:',
     #              i0,jproc0,ihpro,i1,i2,i3,i4
          stop
        endif
      else
        write(*,*)'parcrossing: do not know what to do'
        write(*,*)ifuntype,jproc0
        stop
      endif
      return
      end


      subroutine ckmunwgt(i3)
c Determines parton identity of the final state light parton for the partonic
c processes whose weights involve sums over CKM matrix elements
      implicit none
      integer i3,idisc,ii,ires,iwh
      real * 8 fk88random,xden,xg,wh
      real * 8 ckm2(1:6,1:6)
      common/cckm2/ckm2
      real * 8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer imapp(0:5)
      common/cimapp/imapp
      integer imapd(3)
      data imapd/2,3,5/
c
      if(abs(i3).eq.6235)then
        idisc=0
        xden=rtckm
      elseif(abs(i3).eq.625)then
        idisc=3
        xden=rtckm-ckm2(6,3)
      elseif(abs(i3).eq.635)then
        idisc=2
        xden=rtckm-ckm2(6,2)
      else
        write(*,*)'Error in ckmunwgt:',i3
        stop
      endif
c
      xg=fk88random(ifk88seed)
      wh=0.d0
      iwh=0
      ires=0
      do ii=1,3
        if(imapd(ii).ne.idisc.and.iwh.eq.0)then
          wh=wh+ckm2(6,imapd(ii))/xden
          if(wh.gt.xg)then
            ires=imapd(ii)
            iwh=1
          endif
        endif
      enddo
      i3=sign(1,i3)*imapp(ires)
c
      return
      end


      subroutine checkvv(xsum,xsumabs,iret)
c iret=0 -> all vv entries are equal to zero
c iret=1 -> there is at least one entry which is not zero
c xsum is the sum of all the entries of vv
c xsumabs is the sum of the absolute value of all the entries of vv
      implicit none
      integer jproc,iret,i,itype,itt
      real * 8 vv(7,3,24,2)
      common/cvv/vv
      real * 8 xsum,xsumabs
c
      xsum=0.d0
      xsumabs=0.d0
      iret=0
      do jproc=1,3
        do i=1,7
          do itype=1,24
            do itt=1,2
              if(vv(i,jproc,itype,itt).ne.0.d0)iret=1
              xsum=xsum+vv(i,jproc,itype,itt)
              xsumabs=xsumabs+abs(vv(i,jproc,itype,itt))
            enddo
          enddo
        enddo
      enddo
      return
      end


      subroutine checkvvs(xsum,xsumabs,iret)
c identical to checkvv, except for the fact that works on vvs instead of vv,
c and jproc is not fixed
      implicit none
      integer jproc,iret,i,itype,itt
      real * 8 vvs(7,3,24,2)
      common/cvvs/vvs
      real * 8 xsum,xsumabs
c
      xsum=0.d0
      xsumabs=0.d0
      iret=0
      do jproc=1,3
        do i=1,7
          do itype=1,24
            do itt=1,2
              if(vvs(i,jproc,itype,itt).ne.0.d0)iret=1
              xsum=xsum+vvs(i,jproc,itype,itt)
              xsumabs=xsumabs+abs(vvs(i,jproc,itype,itt))
            enddo
          enddo
        enddo
      enddo
      return
      end


      subroutine getspincorr()
c Determines the lepton momenta, by performing an unweighting using
c the exact real and Born lepton matrix elements
      implicit none
      real*8 pi,tolerance,bdredfact
      parameter (pi=3.14159265358979312D0)
      parameter (tolerance=1.d-2)
c Divide the bound by bdredfact to compensate for any peculiar behaviour of
c tH cross section. May become a process-dependent correction if need be.
      parameter (bdredfact=1.d0)
      character*2 str
      parameter (str='p1')
      include 'stpcblks.h'
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      real*8 ps,px,pyi,pcth1,pcth2
      common/cpsave/ps,px,pyi,pcth1,pcth2
      real*8 tq12
      common/ctvirt/tq12
      real*8 q12,q22
      common/cvirt/q12,q22
      real*8 xm012,ga1,bw1delf,bw1fmmn,xm1low2,xm1upp2
      common/cbw1/xm012,ga1,bw1delf,bw1fmmn,xm1low2,xm1upp2
      real*8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      real*8 xmt,twidth
      common/ctparam/xmt,twidth
      real*8 xmw,gaw
      common/cwparam/xmw,gaw
      integer i0,jproc0,itype0,ich0,itt0
      common/cidproc/i0,jproc0,itype0,ich0,itt0
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
      integer ifuntype
      common/cifuntype/ifuntype
      integer iwidth
      common/ciwidth/iwidth
      integer ichkmom
      common/cichkmom/ichkmom
      integer neventsuw,nqeventsuw,ifailuw
      common/c1iunwgt/neventsuw,nqeventsuw,ifailuw
      integer ncntuws,nqcntuws,nmaxuw,nqmaxuw
      common/c2iunwgt/ncntuws,nqcntuws,nmaxuw,nqmaxuw
      integer idec
      common/cidec/idec
      integer ifk88seed
      common/cifk88seed/ifk88seed
      real*8 xtmp,prob,spcdamp,rrnd,fk88random,yitmp,e1,f1,g1,h1,
     # phitq1,cthtq1,phitq2,cthtq2,o,xbwmass3,
     # rat1,qphsp,q1,tk,uk,q1q,q2q,xdec,xmadevht,unxdec,xsngltht,
     # dmfactb1,dmfact1,phspfact1,xboundb,rat,
     # xinv(5),xtq(4),xbq(4),xtl(4),xtn(4)
      integer iborn,iproj,icross,jjprc,jidr,icntuw,iqcntuw
c
      if(ich0.ne.3)then
        write(*,*)'Error #0 in getspincorr'
        stop
      endif
      if(ichkmom.eq.0)call spccheck(1)
c Set top mass to be used by put_on_shell. Generate it in this
c routine when off-shell effects will be implemented. q22 would be
c the analogous quantity to be used for Higgs off-shell effects
      tq12=xm12
      if(ifuntype.eq.2)then
        if(px.ne.1.d0.or.xmom_cm(3,4).ne.0.d0)then
          write(*,*)'Error #1 in getspincorr'
          stop
        else
          iborn=0
          iproj=0
          if(jproc0.eq.3)then
            icross=0
          else
            icross=1
          endif
          xtmp=px
        endif
      endif
      if(ifuntype.eq.1)then
        prob=spcdamp(px,pyi)
        rrnd=fk88random(ifk88seed)
        if(rrnd.ge.prob)then
c Close to the soft/collinear limits: use Born kinematics for the unweighting
          iborn=0
          iproj=1
          icross=1
          xtmp=1.d0
        else
c Away from the soft/collinear limit: use real kinematics for the unweighting
          iborn=1
          iproj=0
          icross=0
          xtmp=px
        endif
      endif
c yitmp will not be equal to pyi only when real configurations which needed
c to be crossed failed to do so because of a corresponding zero collinear limit
      yitmp=pyi
c When iproj=0, the Born and real kinematics are used to perform unweighting
c for S and H events respectively. When iproj=1, the real kinematics is close 
c to the soft/collinear limits, and the Born is used to unweight. In the case 
c of the gg and qq processes, the Born is chosen according to whether the 
c real configuration is closer to the y->1 (pyi>0 ==> idrlimcp is used) or 
c to the y->-1 (pyi<0 ==> idrlimcm is used) limit. This strategy cannot be
c be adopted in the case of S events due to the gg and qq contributions, since
c for such event we always have py=1. For S events, we rather define jidr
c using the value of i1hpro defined in xout(), where the conversion 
c gg/qq->qg or gg/qq->gq has already been made. Any manipulations on 
c parton identities must also be carried out here. 
c In summary, we have what follows:
c
c          icross=0                               icross=1
c  H events, not projected onto S        H events, projected onto S
c  S events, qq and gg contributions     S events, qg contribution
c          iproj=0                                iproj=1
c  H events, not projected onto S        H events, projected onto S
c  S events
c          iborn=0                                iborn=1
c  Use Born matrix elements              Use real matrix elements 
c
c It may happen that, for H events and when icross=1, the corresponding 
c collinear limit is zero; in such a case, when don't do the crossing,
c and use the real matrix elements, preventing them to be too close to
c the soft/collinear limits. We also set iproj=2 in such a case
      jjprc=0
      jidr=0
      if(icross.eq.0)then
        jjprc=jproc0
        jidr=i0
      elseif(icross.eq.1)then
        jjprc=3
        if(ifuntype.eq.1)then
          if(pyi.ge.0.d0)then
            jidr=idrlimcp(ich0,jproc0,i0)
          else
            jidr=idrlimcm(ich0,jproc0,i0)
          endif
        else
          if(i1hpro.eq.402.or.i1hpro.eq.404)then
            jidr=1
          elseif(i1hpro.eq.405.or.i1hpro.eq.406)then
            jidr=3
          else
            write(*,*)'Error #10 in getspincorr'
            stop
          endif
        endif
        if(jidr.eq.0.and.ifuntype.eq.1)then
          xtmp=px
          yitmp=pyi
          iborn=1
          jjprc=jproc0
          jidr=i0
          if(px.gt.0.995d0)then
            xtmp=0.995d0
            iproj=2
          endif
          if(abs(pyi).gt.0.995d0)then
            yitmp=0.995d0*sign(1.d0,pyi)
            iproj=2
          endif
        endif
      else
        write(*,*)'Error #4 in getspincorr'
        stop
      endif
      if(jjprc.eq.0.or.jidr.eq.0)then
        write(*,*)'Error #11 in getspincorr'
        stop
      endif
c
      neventsuw=neventsuw+1
      icntuw=0
 100  icntuw=icntuw+1
      e1=fk88random(ifk88seed)
      f1=fk88random(ifk88seed)
      g1=fk88random(ifk88seed)
      h1=fk88random(ifk88seed)
      phitq1=2*pi*e1
      cthtq1=-1.d0+2*f1
      phitq2=2*pi*g1
      cthtq2=-1.d0+2*h1
 300  continue
      if(iwidth.eq.1)then
        iqcntuw=0
 200    iqcntuw=iqcntuw+1
        o=fk88random(ifk88seed)
c First distribute q according to the matrix element upper bound,
c which can be done exactly the upper bound being a Breit Wigner
        q12=xbwmass3(o,xm012,ga1,bw1delf,bw1fmmn)
c Then reject some of the values generated according to the phase-space
c q-dependent factor. A 1->1+(1->2) phase-space decomposition has been used.
c Much better here than after computing matrix elements. The following
c form works since qphsp is a function decreasing with q2
        rat1=qphsp(q12,xm12)/qphsp(xm1low2,xm12)
        rrnd=fk88random(ifk88seed)
        if(rat1.lt.rrnd)goto 200
        nqcntuws=nqcntuws+iqcntuw
        if(iqcntuw.gt.nqmaxuw)nqmaxuw=iqcntuw
        nqeventsuw=nqeventsuw+1
        q1=sqrt(q12)
      else
        q12=xm012
        q1=sqrt(q12)
      endif
c No complications here due to off-shell top and hard W; can use the same 
c kinematics for decayed and undecayed matrix elements
      call invar_in(xm12,xm22,ps,xtmp,yitmp,pcth1,pcth2,str,
     #              tk,uk,q1q,q2q,xinv)
      call gentopdmom(xmt,q1,cthtq1,phitq1,cthtq2,phitq2,
     #                xtq,xbq,xtl,xtn,1)
      if(ichkmom.eq.0)call checktdec1(xmt,xtq,xbq,xtl,xtn,1)
      call genhdmom()
      if(ichkmom.eq.0)call checkmom(xmom_cm,ps,0.d0,10,1)
      xdec=xmadevht(iborn,jjprc,jidr,ps,tk,uk,xmom_cm)
      unxdec=xsngltht(iborn,jjprc,jidr,xm12,xm22,ps,xtmp,yitmp,
     #                pcth1,pcth2,tk,uk,q1q,q2q,xinv)
      dmfactb1=256*xm12**2/16.d0
      dmfact1=1/(64.d0*sthw2**2)*
     #        1.d0/((q12-xm012)**2+xm012*ga1**2)
c e^4 -> gw^4; single-top cross sections in mcatnlo_stxsec factorize gw^4
      dmfact1=dmfact1*sthw2**2
      phspfact1=1.d0/(xm12*twidth**2)
      xboundb=dmfactb1*dmfact1*phspfact1
      rat=xdec/((1+tolerance)*unxdec*xboundb)
      rat=rat*bdredfact
      if(rat.gt.1.d0)then
        ifailuw=ifailuw+1
        goto 300
      endif
      rrnd=fk88random(ifk88seed)
      if(rat.lt.rrnd)goto 100
      ncntuws=ncntuws+icntuw
      if(icntuw.gt.nmaxuw)nmaxuw=icntuw
c The event is accepted; regenerate kinematics if Born was used for 
c unweighting H events (to get xmom_cm corresponding to a real emission
c configuration), or if x or y values had been freezed
      if(iproj.eq.0)then
        if(px.ne.xtmp)then
          write(*,*)'Error #6 in getspincorr',px,xtmp
          stop
        endif
      elseif(iproj.eq.1.or.iproj.eq.2)then
        call invar_in(xm12,xm22,ps,px,pyi,pcth1,pcth2,str,
     #                tk,uk,q1q,q2q,xinv)
        call gentopdmom(xmt,q1,cthtq1,phitq1,cthtq2,phitq2,
     #                  xtq,xbq,xtl,xtn,1)
        call genhdmom()
        if(ichkmom.eq.0)call checkmom(xmom_cm,ps,0.d0,20,1)
      else
        write(*,*)'Error #7 in getspincorr'
        stop
      endif
      if(ichkmom.eq.0)call spccheck(2)
      return
      end


      function spcdamp(x,y)
c This function is defined in such a way that
c    spcdamp=0  if  tt=0
c  0<spcdamp<1  if  0<tt<1
c    spcdamp=1  if  tt>1
c and tt is a measure in the (x,y) plane, such that tt=0 in the soft
c and collinear limits (x=1, or y=1, or y=-1), growing monotonically
c away from these limits. In terms of invariants, tt=4*tk*uk/((1-xlim)*s)**2,
c which can easily be generalized for any kind of emissions. 
c Since when spcdamp=1 the real matrix elements are used in the 
c unweighting, xlim has been defined in such a way that, if be_spcfun=1,
c spcdamp is equal to 1 in a region similar to the dead zone. This is
c by no means necessary, and the dependence upon xlim in tt can be
c eliminated altogether. Call this function with al_spcfun>=1, 0<be_spcfun<=1
      implicit none
      real * 8 spcdamp,x,y,xmin,tt,xlim
      parameter (xmin=0.69519410160110384d0)
      real * 8 al_spcfun,be_spcfun
      common/cspcpar/al_spcfun,be_spcfun
c
      xlim=1.d0-be_spcfun+xmin*be_spcfun
      tt=(1-x)**2*(1-y**2)/(1-xlim)**2
      if(tt.lt.0.d0)then
        write(*,*)'Error in spcdamp',tt
        stop
      endif
      if(tt.gt.1.d0)tt=1.d0
      spcdamp=tt**(2*al_spcfun)/
     #       (tt**(2*al_spcfun)+(1-tt)**(2*al_spcfun))
      return
      end


      subroutine spccheck(iflag)
c Stores hard-process four-momenta (in xmom_save) at the beginning of 
c getspincorr(), and checks as the last step of getspincorr() that
c the manipulations carried out there did not change them, by comparing
c xmom_save and xmom_cm
      implicit none
      real*8 tiny,xmom_save(5,4)
      parameter (tiny=1.d-4)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer iflag,i,j,itmp
      save xmom_save
c
      if(iflag.eq.1)then
        itmp=0
        do i=1,5
          do j=1,4
            xmom_save(i,j)=xmom_cm(i,j)
          enddo
        enddo
      elseif(iflag.eq.2)then
        itmp=0
        do i=1,5
          do j=1,4
            if(abs(xmom_save(i,j)-xmom_cm(i,j)).gt.tiny)itmp=1
          enddo
        enddo
      else
        write(*,*)'Wrong call to spccheck'
        stop
      endif
      if(itmp.eq.1)then
        write(*,*)'The check in spccheck failed'
        write(*,*)'Original momenta:'
        do i=1,5
          write(*,900)(xmom_save(i,j),j=1,4)
        enddo
        write(*,*)'  '
        write(*,*)'New momenta:'
        do i=1,5
          write(*,900)(xmom_cm(i,j),j=1,4)
        enddo
        stop
      endif
 900  format(4(1x,d14.8))
      return
      end


      function xsngltht(iborn,jproc,idr,xm12,xm22,s,x,yi,cth1,cth2,
     #                tk,uk,q1q,q2q,xinv)
c Wrapper for the undecayed matrix elements of the original code.
c For Born matrix elements, q1q is t (consistently with invar_in)
      implicit none
      real*8 xsngltht,xm12,xm22,s,x,yi,cth1,cth2,tk,uk,q1q,q2q,xinv(5)
      integer iborn,jproc,idr
      real*8 xmatin(7)
c
      if(iborn.eq.0)then
        call fbornht(s,q1q,jproc,xmatin)
      else
        call frealht(s,x,yi,cth2,tk,uk,q1q,q2q,xinv,jproc,xmatin)
      endif
      xsngltht=xmatin(idr)
      return
      end


      subroutine gentopdmom(xmt,xmw,cth1,phi1,cth2,phi2,
     #                      xtq,xbq,xel,xnu,iqrk)
c Generates the four-momenta of the decay products of the top (iqrk=1).
c These four-momenta are returned in the top rest frame (xbq, xel, xnu; 
c the trivial top momentum is returned as well, xtq). The four-momenta 
c are also boosted to the frame in which the top has momentum xmom_cm(4,*),
c  and the common block xmomcm is filled according to the identifications
c   l+ --> xmom_cm(6,*), nu --> xmom_cm(7,*), b --> xmom_cm(8,*)
c consistently with the labelling conventions used in MC@NLO:
c   x(1)y(2) -> z(3)t(4)[->l+(6)nu(7)b(8)]W-(5)[->l-(9)nub(10)]
c The inputs of the routine are cth1,phi1,cth2,phi2, which are cosines of
c polar angles and azimuthal angles, with
c   (cth1,phi1) --> direction of W in the top rest frame
c   (cth2,phi2) --> direction of l in the W rest frame
c This routine has been derived from that of the ttbar package; the option
c iqrk=2 has been forbidden
      implicit none
      real*8 xmt,xmw,cth1,phi1,cth2,phi2,xtq(4),xbq(4),xel(4),xnu(4)
      integer iqrk
      real*8 tiny,xmt2,xmw2,sth1,sth2,ew,eb,pwx,pwy,pwz,pbx,pby,pbz,
     # eel,enu,pex,pey,pez,pnx,pny,pnz,xmq,tmp(5),tmp1(4),tmp2(4)
      parameter (tiny=1.d-4)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer itop,iel,inu,ib
c
      if(iqrk.eq.1)then
        itop=4
        iel=6
        inu=7
        ib=8
      else
        write(6,*)'gentopdmom called improperly',iqrk
      endif
      xmq=xmom_cm(itop,4)**2-xmom_cm(itop,1)**2-
     #    xmom_cm(itop,2)**2-xmom_cm(itop,3)**2
      xmq=sqrt(xmq)
      if(abs(xmt-xmq).gt.tiny*xmt)then
        write(6,*)'Subroutine gentopdmom'
        write(6,*)'Top is off the mass shell',xmt,xmq
        stop
      endif
c
      xmt2=xmt**2
      xmw2=xmw**2
      sth1=sqrt(1-cth1**2)
      sth2=sqrt(1-cth2**2)
c
      xtq(1)=0.d0
      xtq(2)=0.d0
      xtq(3)=0.d0
      xtq(4)=xmt
c W and b momenta, top rest frame
      ew=(xmt2+xmw2)/(2*xmt)
      eb=(xmt2-xmw2)/(2*xmt)
      pwx=eb*sth1*cos(phi1)
      pwy=eb*sth1*sin(phi1)
      pwz=eb*cth1
      pbx=-pwx
      pby=-pwy
      pbz=-pwz
      xbq(1)=pbx
      xbq(2)=pby
      xbq(3)=pbz
      xbq(4)=eb
c l+ and nu momenta, W rest frame
      eel=xmw/2.d0
      enu=eel
      pex=eel*sth2*cos(phi2)
      pey=eel*sth2*sin(phi2)
      pez=eel*cth2
      pnx=-pex
      pny=-pey
      pnz=-pez
c Boost lepton momenta to top rest frame
      tmp(1)=pwx
      tmp(2)=pwy
      tmp(3)=pwz
      tmp(4)=ew
      tmp(5)=xmw
c Boost l+
      tmp1(1)=pex
      tmp1(2)=pey
      tmp1(3)=pez
      tmp1(4)=eel
      call hwulb4(tmp,tmp1,tmp2)
      xel(1)=tmp2(1)
      xel(2)=tmp2(2)
      xel(3)=tmp2(3)
      xel(4)=tmp2(4)
c Boost nu
      tmp1(1)=pnx
      tmp1(2)=pny
      tmp1(3)=pnz
      tmp1(4)=enu
      call hwulb4(tmp,tmp1,tmp2)
      xnu(1)=tmp2(1)
      xnu(2)=tmp2(2)
      xnu(3)=tmp2(3)
      xnu(4)=tmp2(4)
c Boost all momenta to cm frame
      tmp(1)=xmom_cm(itop,1)
      tmp(2)=xmom_cm(itop,2)
      tmp(3)=xmom_cm(itop,3)
      tmp(4)=xmom_cm(itop,4)
      tmp(5)=xmt
c
      call filltopdec(tmp,xel,iel)
      call filltopdec(tmp,xnu,inu)
      call filltopdec(tmp,xbq,ib)
c
      return
      end


      subroutine filltopdec(tmp,tmp1,ipart)
c Utility routine for gentopdmom; performs the boost and fills xmom_cm 
c for top decay products
      implicit none
      real*8 tmp(5),tmp1(4),tmp2(4)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer ipart
c
      call hwulb4(tmp,tmp1,tmp2)
      xmom_cm(ipart,1)=tmp2(1)
      xmom_cm(ipart,2)=tmp2(2)
      xmom_cm(ipart,3)=tmp2(3)
      xmom_cm(ipart,4)=tmp2(4)
      return
      end


      subroutine genwdmom(xmw,cth1,phi1,xwb,xel,xnu,iwb)
c Generates the four-momenta of the decay products of the W (iwb=2).
c These four-momenta are returned in the W rest frame (xel, xnu; 
c the trivial W momentum is returned as well, xwb). The four-momenta 
c are also boosted to the frame in which the W has momentum xmom_cm(5,*),
c and the common block xmomcm is filled according to the identifications
c   l- --> xmom_cm(9,*), nub --> xmom_cm(10,*)
c consistently with the labelling conventions used in MC@NLO:
c   x(1)y(2) -> z(3)t(4)[->l+(6)nu(7)b(8)]W-(5)[->l-(9)nub(10)]
c The inputs of the routine are cth1,phi1, which are cosines of
c polar angles and azimuthal angles, with
c   (cth1,phi1) --> direction of l in the W rest frame
c This routine must be called with iwb=2, consistently with the convention
c that in Wt production W's #1 and 2 are those from the top decay and
c the hard process respectively
      implicit none
      real*8 xmw,cth1,phi1,xwb(4),xel(4),xnu(4)
      integer iwb
      real*8 tiny,sth1,eel,xmq,tmp(5)
      parameter (tiny=1.d-4)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer iw,iel,inu
c
      if(iwb.eq.2)then
        iw=5
        iel=9
        inu=10
      else
        write(6,*)'genwdmom called improperly',iwb
      endif
      xmq=xmom_cm(iw,4)**2-xmom_cm(iw,1)**2-
     #    xmom_cm(iw,2)**2-xmom_cm(iw,3)**2
      xmq=sqrt(xmq)
      if(abs(xmw-xmq).gt.tiny*xmw)then
        write(6,*)'Subroutine genwmom'
        write(6,*)'W is off the mass shell',xmw,xmq
        stop
      endif
c
      sth1=sqrt(1-cth1**2)
c
      xwb(1)=0.d0
      xwb(2)=0.d0
      xwb(3)=0.d0
      xwb(4)=xmw
c l- and nub momenta, W rest frame
      eel=xmw/2.d0
      xel(1)=eel*sth1*cos(phi1)
      xel(2)=eel*sth1*sin(phi1)
      xel(3)=eel*cth1
      xel(4)=eel
      xnu(1)=-xel(1)
      xnu(2)=-xel(2)
      xnu(3)=-xel(3)
      xnu(4)=eel
c Boost all momenta to cm frame
      tmp(1)=xmom_cm(iw,1)
      tmp(2)=xmom_cm(iw,2)
      tmp(3)=xmom_cm(iw,3)
      tmp(4)=xmom_cm(iw,4)
      tmp(5)=xmw
c
      call fillwdec(tmp,xel,iel)
      call fillwdec(tmp,xnu,inu)
c
      return
      end


      subroutine fillwdec(tmp,tmp1,ipart)
c Utility routine for genwdmom; performs the boost and fills xmom_cm 
c for W decay products. Identical to filltopdec
      implicit none
      real*8 tmp(5),tmp1(4),tmp2(4)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer ipart
c
      call hwulb4(tmp,tmp1,tmp2)
      xmom_cm(ipart,1)=tmp2(1)
      xmom_cm(ipart,2)=tmp2(2)
      xmom_cm(ipart,3)=tmp2(3)
      xmom_cm(ipart,4)=tmp2(4)
      return
      end


      subroutine genhdmom()
c Generates fake momenta of H decay products, by setting one of the momenta
c equal to the Higgs momentum, and the other to zero. This enforces momentum
c conservation and ensures backwards compatibility with Wt case
      implicit none
      integer i
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
c
      do i=1,4
        xmom_cm(9,i)=xmom_cm(5,i)
        xmom_cm(10,i)=0.d0
      enddo
      return
      end


      function qphsp(q12,xmt2)
c Non-trivial factor of the t->bW phase space, in the t rest frame; q12 is
c the W mass squared
      implicit none
      real*8 qphsp,q12,xmt2,tmp
c
      tmp=0.d0
      if(q12.gt.0.d0.and.q12.lt.xmt2)tmp=(xmt2-q12)/(2.d0*xmt2)
      qphsp=tmp
      return
      end


      function itoosoftkin()
c Returns 1 when a three-body kinematics can be safely approximated
c with a two-body kinematics. It is useful when three-body NLO configurations
c are obtained, which cannot be produced through showering
      implicit none
      integer itoosoftkin,itmp
c
      itmp=0
      itoosoftkin=itmp
      return
      end
c
c
c End of event-generation routines
c
c
c
c
c
c Begin of phase-space routines
c
c
      subroutine invar_in(xm12,xm22,s,x,y,cth1,cth2,str,
     #                    tk,uk,q1q,q2q,xinv)
c This routine has been obtained by modifying the analogous routine
c in the VH code. 
c The names of the invariants are taken from Nucl.Phys.B383:3-44,1992 [FNR] 
c (q1q is q_1 of the paper, q2q is q_2, q1c is \hat{q}_1, q2c is \hat{q}_2).
c
c The hard process is
c   a(p1)+b(p2) --> t(k1)+c(k2)+d(k)
c where a, b, c and d are light partons, t is top quark with k1^2=xm12.
c The quarks t and c are attached to the W-vertex, k2^2=xm22=0. The process 
c can be described by the same invariants as in FNR [eqs.(2.6) and (2.7)].
c
c In terms of the
c invariants, the dot products are 
c
c    p1.p2 = s/2
c    p1.k  = -tk/2
c    p2.k  = -uk/2
c    p1.k1 = -(q1q-xm12)/2
c    p2.k2 = -(q2q-xm22)/2
c    k1.k2 = (s2-xm12-xm22)/2
c    p2.k1 = -(q2c-xm12)/2
c    p1.k2 = -(q1c-xm22)/2
c    k.k1  = (w1-xm12)/2
c    k.k2  = (w2-xm22)/2
c
c The four momenta are given in the t-c rest frame as follows
c     p1 = p10*(1,0,spsi2,cpsi2)
c     p2 = p20*(1,0,spsi ,cpsi )
c     k  = k0*(1,0,spsi1,cpsi1).
c     k1 = (k10, bx*sth2*sth1, bx*cth2*sth1, bx*cth1)
c     k2 = (k20,-bx*sth2*sth1,-bx*cth2*sth1,-bx*cth1).
c The argument str should be set to 'p1': then p1 = p10 (1,0,0,1) (psi2 =0), 
c with psi and psi1 determined using momentum conservation; according to the 
c work done for Drell Yan, the other options for str have been disabled.
c
c The four momenta of the partons in the c.m. frame of the incoming
c partons are stored in xmom_cm(ipart,icomp), with the conventions:
c   icomp=1 -> px, icomp=2 -> py, icomp=3 -> pz, icomp=4 -> E;
c   ipart=1 -> p1, ipart=2 -> p2, ipart=3 -> k, ipart=4 -> k1, ipart=5 -> k2.
c
c Notice that  bx = sqrt(s2)/2 * beta_x[FNR paper]
c
c
      implicit none
      real * 8 xm12,xm22,s,x,y,cth1,cth2,tk,uk,q1q,q2q,xinv(5)
      character * 2 str
c
      real * 8 ptv1,ptv2,ptvg,y1,y2,yg
      common/perpen/ptv1(2),ptv2(2),ptvg(2)
      common/ycmvar/y1,y2,yg
c
      real * 8 s2,drs2,p10,p20,k0,k10,k20,bx,sth1,cpsi,
     # spsi,cpsi2,spsi2,cpsi1,spsi1,xktsq,xkt1sq,xkt2sq,
     # xkt,xkt1,xkt2,tmp,sqs,tiny,zero,sth2,q1c,q2c,w1,w2,
     # e1lab,pl1lab,e2lab,pl2lab,beta,xcpsi1,xspsi1
      parameter (tiny=1.d-8)
      parameter (zero=0.d0)
      real*8 xmom_cm(10,4)
      common/cxmomcm/xmom_cm
      integer ichkmom
      common/cichkmom/ichkmom
c     
      sqs=sqrt(s)
      if(sqs.lt.(sqrt(xm12)+sqrt(xm22)))then
        write(*,*)'Error in invar: energy too small',s,xm12,xm22
        stop
      endif
      tk=-s*(1-x)*(1-y)/2.d0
      uk=-s*(1-x)*(1+y)/2.d0
      s2 = tk+uk+s
      drs2 = 2*sqrt(s2)
      p10 = (s+tk)/drs2
      p20 = (s+uk)/drs2
      k0  = -(tk+uk)/drs2
      k10 = (s2+xm12-xm22)/drs2
      k20 = (s2+xm22-xm12)/drs2
      bx=sqrt(s2**2+xm22**2+xm12**2-2*(s2*xm22+s2*xm12+xm22*xm12))/drs2
      sth1 = sqrt(1-cth1**2)
      sth2 = sqrt(1-cth2**2)
      if(str.eq.'p1') then
         cpsi2 = 1
         spsi2 = 0
         cpsi = 1-8*x/((1+y+x*(1-y))*(1-y+x*(1+y)))
         spsi = 4*(1-x)*sqrt(x*(1-y**2))/
     #          ((1+y+x*(1-y))*(1-y+x*(1+y)))
         cpsi1 = (1+y-x*(1-y))/(1+y+x*(1-y))
         spsi1 = sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
      else
         write(6,*) 'Error in invar: str=',str
         stop
      endif
      q1q = xm12 - 2*p10*(k10-bx*(cth2*sth1*spsi2+cth1*cpsi2))
      q2q = xm22 - 2*p20*(k20+bx*(cth2*sth1*spsi +cth1*cpsi ))
      q1c = xm12 + xm22 - s - tk - q1q
      q2c = xm12 + xm22 - s - uk - q2q
      w1  = xm12 - q1q + q2q - tk
      w2  = xm22 - q2q + q1q - uk
c Here define xinv, according to
c   p_i.k = sqrt{s}*(1-x)*xinv(i)    i=1,2
c   k_i.k = sqrt{s}*(1-x)*xinv(i+3)  i=1,2
c and for consistency xinv(3)=0; xinv thus factor out analytically the 
c xi-dependence. These quantities are used to compute the soft limits
c of matrix elements and S_in 
      if(x.ne.1.d0)then
        beta=4*bx/drs2
        xcpsi1 = cpsi1
        xspsi1 = spsi1
      else
        beta=sqrt(1-2*(xm12+xm22)/s+(xm12-xm22)**2/s**2)
        xcpsi1 = y
        xspsi1 = sqrt(1-y**2)
      endif
      xinv(1) = sqs*(1-y)/4.d0
      xinv(2) = sqs*(1+y)/4.d0
      xinv(3) = 0d0
      xinv(4) = sqs/4.d0*(1-(xm22-xm12)/(x*s)-
     #    beta*(cth2*sth1*xspsi1+cth1*xcpsi1))
      xinv(5) = sqs/4.d0*(1+(xm22-xm12)/(x*s)+
     #    beta*(cth2*sth1*xspsi1+cth1*xcpsi1))
c
c Recall: y1,y2, yg are rapidities in the partonic cm frame
      if(abs(q1q-xm12).lt.tiny) then
        y1  = 1.d8
      elseif(abs(q2c-xm12).lt.tiny) then
        y1  = -1.d8
      else
        y1 = .5d0*log( (xm12-q2c)/(xm12-q1q) )
      endif
      if(abs(q1c-xm22).lt.tiny) then
        y2  = 1.d8
      elseif(abs(q2q-xm22).lt.tiny) then
        y2  = -1.d8
      else
        y2 = .5d0*log( (xm22-q2q)/(xm22-q1c) )
      endif
      if(abs(tk).lt.tiny) then
        yg  = 1.d8
      elseif(abs(uk).lt.tiny) then
        yg  = -1.d8
      else
        yg  = .5d0*log( uk/tk )
      endif
c-----------------------------------------------------------------
c xktsq, xkt1sq e xkt2sq are the square of transverse momenta of d, t,
c and c respectively. The axis orientation is such that t is always
c along the x direction. The component of p_T(t) along the y direction
c is always positive or zero
c
      xktsq = uk*tk/s
      if(xktsq.lt.tiny) then
         ptv1(1) = bx*sth1
         ptv1(2) = 0.d0
         ptv2(1) = -ptv1(1)
         ptv2(2) = 0.d0
         ptvg(1) = 0.d0
         ptvg(2) = 0.d0
         xkt1 = ptv1(1)
         xkt2 = xkt1
      else
         xkt1sq = (xm12-q2c)*(xm12-q1q)/s - xm12
         xkt2sq = (xm22-q2q)*(xm22-q1c)/s - xm22
         xkt = sqrt(xktsq)
         xkt1 = sqrt(xkt1sq)
         xkt2 = sqrt(xkt2sq)
         ptv1(1) = xkt1
         ptv1(2) = 0.d0
         ptv2(1) = (xktsq-xkt1sq-xkt2sq)/(2.d0*xkt1)
         tmp = xkt2sq-ptv2(1)**2
         if(tmp.gt.0.d0)then
            ptv2(2) = sqrt(tmp)
         else
            ptv2(2) = 0.d0
         endif
         ptvg(1) = (xkt2sq-xkt1sq-xktsq)/(2.d0*xkt1)
         tmp = xktsq-ptvg(1)**2
         if(tmp.gt.0.d0)then
            ptvg(2) = -sqrt(tmp)
         else
            ptvg(2) = 0.d0
         endif
      endif
      if(ichkmom.eq.0)call checkptcon(ptv1,ptv2,ptvg)
c
c xmom_cm(1,mu) = p1(mu)
      xmom_cm(1,1)=0.d0
      xmom_cm(1,2)=0.d0
      xmom_cm(1,3)=sqs/2.d0
      xmom_cm(1,4)=sqs/2.d0
c xmom_cm(2,mu) = p2(mu)
      xmom_cm(2,1)=0.d0
      xmom_cm(2,2)=0.d0
      xmom_cm(2,3)=-sqs/2.d0
      xmom_cm(2,4)=sqs/2.d0
c xmom_cm(3,mu) = k(mu)
      if(tk.eq.0.d0.and.uk.eq.0.d0)then
        xmom_cm(3,1)=0.d0
        xmom_cm(3,2)=0.d0
        xmom_cm(3,3)=0.d0
        xmom_cm(3,4)=0.d0
      elseif(tk.eq.0)then
        xmom_cm(3,1)=0.d0
        xmom_cm(3,2)=0.d0
        xmom_cm(3,3)=-uk/(2*sqs)
        xmom_cm(3,4)=xmom_cm(3,3)
      elseif(uk.eq.0)then
        xmom_cm(3,1)=0.d0
        xmom_cm(3,2)=0.d0
        xmom_cm(3,3)=tk/(2*sqs)
        xmom_cm(3,4)=-xmom_cm(3,3)
      else
        xmom_cm(3,1)=ptvg(1)
        xmom_cm(3,2)=ptvg(2)
        xmom_cm(3,3)=sqs/2.d0*(1-x)*y
        xmom_cm(3,4)=sqs/2.d0*(1-x)
      endif
c xmom_cm(4,mu) = k1(mu)
      e1lab=(2*xm12-q1q-q2c)/(2*sqs)
      pl1lab=(q1q-q2c)/(2*sqs)
      xmom_cm(4,1)=ptv1(1)
      xmom_cm(4,2)=ptv1(2)
      xmom_cm(4,3)=pl1lab
      xmom_cm(4,4)=e1lab
c xmom_cm(5,mu) = k2(mu)
      e2lab=(2*xm22-q1c-q2q)/(2*sqs)
      pl2lab=(q1c-q2q)/(2*sqs)
      xmom_cm(5,1)=ptv2(1)
      xmom_cm(5,2)=ptv2(2)
      xmom_cm(5,3)=pl2lab
      xmom_cm(5,4)=e2lab
c
      if(ichkmom.eq.0) call checkmom(xmom_cm,s,0.d0,1,2)
      return
      end


      subroutine checkmom(xmom,smax,ybst,iflag,itype)
      implicit none
      real * 8 xmom(10,4)
      real * 8 smax,ybst,xpmax
      real*8 x1,x2
      common/cx1x2/x1,x2
      real * 8 tiny,vtiny,xsum(4),xsuma(4),xsign,xrat(4)
      parameter (tiny=5.d-3)
      parameter (vtiny=1.d-4)
      integer iflag,itype,i,j,jj,jflag,jeflag,jmax
c
      if(itype.eq.1)then
        jmax=10
      elseif(itype.eq.2)then
        jmax=5
      elseif(itype.eq.3)then
        jmax=8
      else
        write(6,*)'Wrong option in checkmom'
        stop
      endif
      jflag=0
      jeflag=0
      xpmax=sqrt(smax)/2.d0*(1+vtiny)
      do i=1,4
        xsum(i)=0.d0
        xsuma(i)=0.d0
        do j=1,jmax
          if((itype.eq.1.and.j.ne.4.and.j.ne.5).or.itype.eq.2.or.
     #       (itype.eq.3.and.j.ne.4.and.j.ne.9.and.j.ne.10))then
            if(i.ne.4.and.xmom(j,i).gt.xpmax)jeflag=1
            xsign=1.d0
            if(j.eq.1.or.j.eq.2)xsign=-1.d0
            xsum(i)=xsum(i)+xmom(j,i)*xsign
            xsuma(i)=xsuma(i)+abs(xmom(j,i))
          endif
        enddo
        if(xsuma(i).lt.1.d0)then
          xrat(i)=abs(xsum(i))
        else
          xrat(i)=abs(xsum(i))/xsuma(i)
        endif
        if(xrat(i).gt.tiny.and.jflag.eq.0)then
          write(*,*)'Momentum is not conserved'
          write(*,*)'iflag,i=',iflag,i
          write(*,*)'smax,y=',smax,ybst
          write(*,*)'x1,x2=',x1,x2
          do j=1,10
            write(*,'(4(d14.8,1x))') (xmom(j,jj),jj=1,4)
          enddo
          jflag=1
        endif
      enddo
      if(jflag.eq.1)then
        write(*,'(4(d14.8,1x))') (xsum(jj),jj=1,4)
        write(*,'(4(d14.8,1x))') (xrat(jj),jj=1,4)
        stop
      endif
      if(jeflag.eq.1)then
        write(*,*)'Momentum component larger than sqrt(s)/2'
        write(*,*)'iflag=',iflag
        write(*,*)'s,pmax,y=',smax,xpmax,ybst
        write(*,*)'x1,x2=',x1,x2
        do j=1,10
          write(*,'(4(d14.8,1x))') (xmom(j,jj),jj=1,4)
        enddo
        stop
      endif
      return
      end


      subroutine checktdec1(xmt,xtq,xbq,xel,xnu,itop)
c Checks momentum conservation in top decay
      implicit none
      real*8 xmt,tiny,diff,xm1,xm2,xtq(4),xbq(4),xel(4),xnu(4)
      parameter (tiny=1.d-8)
      integer itop,i
c
      xm1=0.d0
      xm2=0.d0
      do i=1,4
        if(i.le.3)then
          xm1=xm1+(xbq(i)+xel(i)+xnu(i))**2
          xm2=xm2+xtq(i)**2
        else
          xm1=(xbq(i)+xel(i)+xnu(i))**2-xm1
          xm2=xtq(i)**2-xm2
        endif
        diff=xtq(i)-xbq(i)-xel(i)-xnu(i)
        if(abs(diff).gt.tiny*xmt)then
          write(6,*)'Subroutine checktdec1'
          write(6,*)'Momentum is not conserved in decay',i,itop
          stop
        endif
      enddo
      xm1=sqrt(xm1)
      xm2=sqrt(xm2)
      if(abs(xm1-xmt).gt.tiny*xmt.or.abs(xm2-xmt).gt.tiny*xmt)then
        write(6,*)'Subroutine checktdec1'
        write(6,*)'Top is off the mass shell',xm1,xm2,xmt
        stop
      endif
      return
      end


      subroutine checktdec2(xmom,idec,iprod1,iprod2,iprod3)
c Checks momentum conservation in top decay, after manipulations 
c in put_on_shell()
      implicit none
      real * 8 xmom(10,4)
      real * 8 tiny,vtiny,xm1,xm2,xsum(4),xsuma(4),xrat(4)
      parameter (tiny=5.d-3)
      parameter (vtiny=1.d-8)
      integer idec,iprod1,iprod2,iprod3,jflag,i,jj
c
      jflag=0
      xm1=0.d0
      xm2=0.d0
      do i=1,4
        if(i.le.3)then
          xm1=xm1+(xmom(iprod1,i)+xmom(iprod2,i)+xmom(iprod3,i))**2
          xm2=xm2+xmom(idec,i)**2
        else
          xm1=(xmom(iprod1,i)+xmom(iprod2,i)+xmom(iprod3,i))**2-xm1
          xm2=xmom(idec,i)**2-xm2
        endif
        xsum(i)=xmom(idec,i)-xmom(iprod1,i)-
     #          xmom(iprod2,i)-xmom(iprod3,i)
        xsuma(i)=abs(xmom(idec,i))+abs(xmom(iprod1,i))+
     #           abs(xmom(iprod2,i))+abs(xmom(iprod3,i))
        if(xsuma(i).lt.1.d0)then
          xrat(i)=abs(xsum(i))
        else
          xrat(i)=abs(xsum(i))/xsuma(i)
        endif
        if(xrat(i).gt.tiny.and.jflag.eq.0)then
          write(*,*)'Subroutine checktdec2'
          write(*,*)'Momentum is not conserved'
          write(*,*)idec,iprod1,iprod2,iprod3
          write(*,'(4(d14.8,1x))') (xmom(idec,jj),jj=1,4)
          write(*,'(4(d14.8,1x))') (xmom(iprod1,jj),jj=1,4)
          write(*,'(4(d14.8,1x))') (xmom(iprod2,jj),jj=1,4)
          write(*,'(4(d14.8,1x))') (xmom(iprod3,jj),jj=1,4)
          jflag=1
        endif
      enddo
      if(jflag.eq.1)then
        write(*,'(4(d14.8,1x))') (xsum(jj),jj=1,4)
        write(*,'(4(d14.8,1x))') (xrat(jj),jj=1,4)
        stop
      endif
      xm1=sqrt(xm1)
      xm2=sqrt(xm2)
      if(abs(xm1-xm2).gt.tiny*xm1)then
        write(6,*)'Subroutine checktdec2'
        write(6,*)'Top is off the mass shell',xm1,xm2
        stop
      endif
      return
      end


      subroutine checkwdec1(xmw,xwb,xel,xnu,iw)
c Checks momentum conservation in W decay
      implicit none
      real*8 xmw,tiny,diff,xm1,xm2,xwb(4),xel(4),xnu(4)
      parameter (tiny=1.d-8)
      integer iw,i
c
      xm1=0.d0
      xm2=0.d0
      do i=1,4
        if(i.le.3)then
          xm1=xm1+(xel(i)+xnu(i))**2
          xm2=xm2+xwb(i)**2
        else
          xm1=(xel(i)+xnu(i))**2-xm1
          xm2=xwb(i)**2-xm2
        endif
        diff=xwb(i)-xel(i)-xnu(i)
        if(abs(diff).gt.tiny*xmw)then
          write(6,*)'Subroutine checkwdec1'
          write(6,*)'Momentum is not conserved in decay',i,iw
          stop
        endif
      enddo
      xm1=sqrt(xm1)
      xm2=sqrt(xm2)
      if(abs(xm1-xmw).gt.tiny*xmw.or.abs(xm2-xmw).gt.tiny*xmw)then
        write(6,*)'Subroutine checkwdec1'
        write(6,*)'Top is off the mass shell',xm1,xm2,xmw
        stop
      endif
      return
      end


      subroutine checkwdec2(xmom,idec,iprod1,iprod2)
c Checks momentum conservation in W decay, after manipulations 
c in put_on_shell()
      implicit none
      real * 8 xmom(10,4)
      real * 8 tiny,vtiny,xm1,xm2,xsum(4),xsuma(4),xrat(4)
      parameter (tiny=5.d-3)
      parameter (vtiny=1.d-8)
      integer idec,iprod1,iprod2,jflag,i,jj
c
      jflag=0
      xm1=0.d0
      xm2=0.d0
      do i=1,4
        if(i.le.3)then
          xm1=xm1+(xmom(iprod1,i)+xmom(iprod2,i))**2
          xm2=xm2+xmom(idec,i)**2
        else
          xm1=(xmom(iprod1,i)+xmom(iprod2,i))**2-xm1
          xm2=xmom(idec,i)**2-xm2
        endif
        xsum(i)=xmom(idec,i)-xmom(iprod1,i)-xmom(iprod2,i)
        xsuma(i)=abs(xmom(idec,i))+abs(xmom(iprod1,i))+
     #           abs(xmom(iprod2,i))
        if(xsuma(i).lt.1.d0)then
          xrat(i)=abs(xsum(i))
        else
          xrat(i)=abs(xsum(i))/xsuma(i)
        endif
        if(xrat(i).gt.tiny.and.jflag.eq.0)then
          write(*,*)'Subroutine checkwdec2'
          write(*,*)'Momentum is not conserved'
          write(*,*)idec,iprod1,iprod2
          write(*,'(4(d14.8,1x))') (xmom(idec,jj),jj=1,4)
          write(*,'(4(d14.8,1x))') (xmom(iprod1,jj),jj=1,4)
          write(*,'(4(d14.8,1x))') (xmom(iprod2,jj),jj=1,4)
          jflag=1
        endif
      enddo
      if(jflag.eq.1)then
        write(*,'(4(d14.8,1x))') (xsum(jj),jj=1,4)
        write(*,'(4(d14.8,1x))') (xrat(jj),jj=1,4)
        stop
      endif
      xm1=sqrt(xm1)
      xm2=sqrt(xm2)
      if(abs(xm1-xm2).gt.tiny*xm1)then
        write(6,*)'Subroutine checkwdec2'
        write(6,*)'W is off the mass shell',xm1,xm2
        stop
      endif
      return
      end


      subroutine checkptcon(ptvl1,ptvl2,ptvg)
      implicit none
      real*8 ptvl1(2),ptvl2(2),ptvg(2),tiny,pt1,pt2,ptmax
      parameter (tiny=1.d-5)
      integer jj
c
      ptmax=max(abs(ptvl1(1)),abs(ptvl2(1)),abs(ptvg(1)),
     #          abs(ptvl1(2)),abs(ptvl2(2)),abs(ptvg(2)))
      pt1=ptvl1(1)+ptvl2(1)+ptvg(1)
      pt2=ptvl1(2)+ptvl2(2)+ptvg(2)
      if(pt1.gt.ptmax*tiny.or.pt2.gt.ptmax*tiny)then
        write(*,*)'Transverse momentum is not conserved'
        write(*,'(4(d14.8,1x))') (ptvl1(jj),jj=1,2)
        write(*,'(4(d14.8,1x))') (ptvl2(jj),jj=1,2)
        write(*,'(4(d14.8,1x))') (ptvg(jj),jj=1,2)
        stop
      endif
      return
      end


      function bwfunc(s,xm02,gah)
c Returns the Breit Wigner function, normalized in such a way that
c its integral in the range (-inf,inf) is one
      implicit none
      real*8 bwfunc,s,xm02,gah
      real*8 pi,xm0
      parameter (pi=3.1415926535897932d0)
c
      xm0=sqrt(xm02)
      bwfunc=xm0*gah/(pi*((s-xm02)**2+xm02*gah**2))
      return
      end


      function xbwmass3(t,xm02,ga,bwdelf,bwfmmn)
c Returns the boson mass squared, given 0<t<1, the nominal mass (xm0),
c and the mass range (implicit in bwdelf and bwfmmn). This function
c is the inverse of F(M^2), where
c   F(M^2)=\int_{xmlow2}^{M^2} ds BW(sqrt(s),M0,Ga)
c   BW(M,M0,Ga)=M0 Ga/pi 1/((M^2-M0^2)^2+M0^2 Ga^2
c and therefore eats up the Breit-Wigner when changing integration 
c variable M^2 --> t
      implicit none
      real*8 xbwmass3,t,xm02,ga,bwdelf,bwfmmn
      real*8 pi,xm0
      parameter (pi=3.1415926535897932d0)
c
      xm0=sqrt(xm02)
      xbwmass3=xm02+xm0*ga*tan(pi*bwdelf*t-bwfmmn)
      return
      end


      subroutine zzchvar(parth1,cth1,xjac,ro)
c
c Given 0<parth1<1 returns -1<cth1<1
c and multiplies xjac times the d cth1 / d parth1 jacobian
c
      implicit none
      real * 8 parth1,cth1,xjac,ro,bb,xlgbb,yy,expyy
      bb = 1-ro**2/16
      xlgbb = log((1+bb)/(1-bb))
      yy = ( parth1 * 2 - 1 ) * xlgbb
      xjac = xjac * 2 * xlgbb
      expyy = exp(-yy)
      cth1 = (1-expyy)/(1+expyy)/bb
      xjac = xjac * 2 * expyy/(1+expyy)**2 / bb
      return
      end
c
c
c End of phase-space routines
c
c
      FUNCTION FK88RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,FK88RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      FK88RANDOM = SEED*MINV
      END
c
c
c Initialization
c
c
      subroutine setpar()
      implicit none
      include 'stpcblks.h'
      real * 8 pi,aem,cthw2,sthw2
      real * 8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      parameter (pi=3.14159265358979312D0)
      real * 8 ckm(1:6,1:6),ckm2(1:6,1:6),vickm(1:6,1:6)
      common/cckm2/ckm2
      common/cvickm/vickm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cweinan/sthw2,cthw2
      integer i,j,k,l,idrmax(1:3,3)
      common/cidrmax/idrmax
      integer iapcp(3,7),iapcm(3,7)
      common/ciap/iapcp,iapcm
      integer idrlimcp(3,1:3,8),idrlimcm(3,1:3,8)
      common/cidrlims/idrlimcp,idrlimcm
      integer itypemax(2:3)
      common/citypemax/itypemax
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer ialwsplit(3,3,6)
      common/cialwsplit/ialwsplit
      integer ie0ht(3,7,3,6)
      common/cie0ht/ie0ht
      integer imcmult(3,7,3)
      common/cimcmult/imcmult
c
c Fermi constant, from PDG2002
      gf=1.16639d-5
c alpha_em
      aem=1/137.0359895d0
c electron charge squared
      ze2=4*pi*aem
c sin and cos squared of theta_W; MSbar scheme, from PDG2003
      sthw2=0.23113d0
      cthw2=1-sthw2
c ckm(i,j)=|CKM matrix elements|, with  i=1,4,6 --> up,charm,top
c                                       j=2,3,5 --> down,strange,bottom
      if(vickm(1,2).eq.0.d0.and.vickm(1,3).eq.0.d0.and.
     #   vickm(1,5).eq.0.d0)then
        do i=1,6
          do j=1,6
            ckm(i,j)=0.d0
          enddo
        enddo
c Values from PDG 2003.
c Centers of the ranges given in eq.(11.2), supposedly taking unitarity
c into account; with the following entries, it holds better than 0.1%
        ckm(1,2)=0.9748d0
        ckm(1,3)=0.2225d0
        ckm(1,5)=0.0036d0
        ckm(4,2)=0.2225d0
        ckm(4,3)=0.9740d0
        ckm(4,5)=0.041d0
        ckm(6,2)=0.009d0
        ckm(6,3)=0.0405d0
        ckm(6,5)=0.9992d0
      else
        do i=1,6
          do j=1,6
            ckm(i,j)=vickm(i,j)
          enddo
        enddo
      endif
      do i=1,6
        do j=1,6
          ckm2(i,j)=ckm(i,j)**2
        enddo
      enddo
c Combines Higgs couplings to CKM matrix
      call combine_coupl()
c Combinations used in strfunht; need them also for unweighting
      ruckm=ckm2(1,2)+ckm2(1,3)+ckm2(1,5)
      rcckm=ckm2(4,2)+ckm2(4,3)+ckm2(4,5)
      rtckm=ckm2(6,2)+ckm2(6,3)+ckm2(6,5)
      rducckm=ckm2(1,2)+ckm2(4,2)
      rsucckm=ckm2(1,3)+ckm2(4,3)
      rbucckm=ckm2(1,5)+ckm2(4,5)
c Fills the array idrmap used in the main code; this uses the information
c on t production to get tbar production matrix elements
      call idrfill()
c Fills the array idrmax(jproc,ich) of maximum values for idr; depends 
c on process type and channel
      do i=2,3
        do j=1,2
          idrmax(i,j)=4
        enddo
      enddo
      idrmax(3,2)=8
      idrmax(1,3)=1
      idrmax(2,3)=7
      idrmax(3,3)=3
c idrlimcp(ich,jproc,idr) returns the idr code relevant to the Born matrix 
c element that factorizes when the y->1 collinear limit is taken in the real 
c matrix element identified by (ich,jproc,idr); idrlimcm has a similar
c meaning for the y->-1 limit. If the limit is not singular, idrlimcp=0.
      do k=1,3
        do i=1,3
          do j=1,8
            idrlimcp(k,i,j)=0
            idrlimcm(k,i,j)=0
          enddo
        enddo
      enddo
c s- and t-channels
      do k=1,2
        do j=1,idrmax(2,k)
          idrlimcp(k,2,j)=j
          idrlimcm(k,2,j)=j
        enddo
      enddo
      idrlimcm(1,3,1)=1
      idrlimcm(1,3,2)=3
      idrlimcp(1,3,3)=3
      idrlimcp(1,3,4)=1
      idrlimcm(2,3,1)=1
      idrlimcm(2,3,3)=2
      idrlimcm(2,3,5)=3
      idrlimcm(2,3,7)=4
      idrlimcp(2,3,2)=1
      idrlimcp(2,3,4)=2
      idrlimcp(2,3,6)=3
      idrlimcp(2,3,8)=4
c Ht-channel, gg initial state
      idrlimcp(3,1,1)=1
      idrlimcm(3,1,1)=3
c Ht-channel, qq initial state
      idrlimcm(3,2,3)=1
      idrlimcp(3,2,4)=3
      idrlimcm(3,2,5)=1
      idrlimcp(3,2,6)=3
      idrlimcm(3,2,7)=1
      idrlimcp(3,2,7)=3
c Ht-channel, qg initial state
      idrlimcp(3,3,1)=1
      idrlimcm(3,3,1)=1
      idrlimcp(3,3,3)=3
      idrlimcm(3,3,3)=3
c itypemax(jproc) is the maximum value the index itype can get -- used
c for s- and t-channels, see the routine strfun
      itypemax(2)=9
      itypemax(3)=3
c jtypemax(jproc,idr) is the maximum value the index jtype can get -- used
c for Ht-channel, see the routine strfunht
      jtypemax(1,1)=1
      jtypemax(2,1)=4
      jtypemax(2,2)=4
      jtypemax(2,3)=24
      jtypemax(2,4)=24
      jtypemax(2,5)=3
      jtypemax(2,6)=3
      jtypemax(2,7)=3
      jtypemax(3,1)=3
      jtypemax(3,3)=3
c iapcp(jproc,idr) returns the Altarelli-Parisi code relevant to the
c collinear y=1 splitting for the Ht channel. Set to zero if unused or
c if the splitting is not allowed; iapcm has the same meaning for the 
c y=-1 splitting
      do i=1,3
        do j=1,7
          iapcp(i,j)=0
          iapcm(i,j)=0
        enddo
      enddo
      iapcp(1,1)=2
      iapcm(1,1)=2
      iapcm(2,3)=3
      iapcp(2,4)=3
      iapcm(2,5)=3
      iapcp(2,6)=3
      iapcp(2,7)=3
      iapcm(2,7)=3
      iapcp(3,1)=4
      iapcm(3,1)=1
      iapcp(3,3)=1
      iapcm(3,3)=4
c ialwsplit(jproc,ileg,ie0sq)=1 if there is at least one non-zero
c contribution from the MC subtraction terms which corresponds
c to (jproc,ileg,ie0sq) -- no information on idr is needed, since
c this quantity is relevant to shower variables and scale 
      do i=1,3
        do j=1,3
          do k=1,6
            ialwsplit(i,j,k)=0
          enddo
        enddo
      enddo
      ialwsplit(1,1,1)=1
      ialwsplit(1,2,1)=1
      ialwsplit(2,1,1)=1
      ialwsplit(2,1,2)=1
      ialwsplit(2,2,1)=1
      ialwsplit(2,2,4)=1
      ialwsplit(3,1,1)=1
      ialwsplit(3,1,2)=1
      ialwsplit(3,2,1)=1
      ialwsplit(3,2,4)=1
      ialwsplit(3,3,2)=1
      ialwsplit(3,3,4)=1
c ie0ht(jproc,idr,ileg,ie0sq)=ie0sq if the corresponding combination in
c the MC subtraction terms is an allowed branching. Set to zero otherwise.
      do i=1,3
        do l=1,7
          do j=1,3
            do k=1,6
              ie0ht(i,l,j,k)=0
            enddo
          enddo
        enddo
      enddo
      ie0ht(1,1,1,1)=1
      ie0ht(1,1,2,1)=1
      ie0ht(2,3,2,1)=1
      ie0ht(2,3,2,4)=4
      ie0ht(2,4,1,1)=1
      ie0ht(2,4,1,2)=2
      ie0ht(2,5,2,1)=1
      ie0ht(2,5,2,4)=4
      ie0ht(2,6,1,1)=1
      ie0ht(2,6,1,2)=2
      ie0ht(2,7,1,1)=1
      ie0ht(2,7,1,2)=2
      ie0ht(2,7,2,1)=1
      ie0ht(2,7,2,4)=4
      ie0ht(3,1,1,2)=2
      ie0ht(3,1,2,1)=1
      ie0ht(3,1,2,4)=4
      ie0ht(3,1,3,2)=2
      ie0ht(3,3,1,1)=1
      ie0ht(3,3,1,2)=2
      ie0ht(3,3,2,4)=4
      ie0ht(3,3,3,4)=4
c imcmult(jproc,idr,ileg) returns the pre-factor for the MC subtraction
c terms that takes into account the fact that gluons and quarks have two
c and one colour partners respectively. In the case of no branching, this
c quantity is set to zero to cause the program to crash if improperly used
      do i=1,3
        do l=1,7
          do j=1,3
            imcmult(i,l,j)=0
          enddo
        enddo
      enddo
      imcmult(1,1,1)=1
      imcmult(1,1,2)=1
      imcmult(2,3,2)=2
      imcmult(2,4,1)=2
      imcmult(2,5,2)=2
      imcmult(2,6,1)=2
      imcmult(2,7,1)=2
      imcmult(2,7,2)=2
      imcmult(3,1,1)=1
      imcmult(3,1,2)=2
      imcmult(3,1,3)=1
      imcmult(3,3,1)=2
      imcmult(3,3,2)=1
      imcmult(3,3,3)=1
c
      return
      end


      subroutine setpar_2hdm()
      implicit none
      real * 8 xmw
c Values from PDG 2003
      parameter (xmw=80.425d0)
      integer i,j
c aij() is the scalar Higgs coupling matrix, with CKM factored out
      real*8 aij(1:6,1:6)  ! 
c aij() is the pseudo-scalar Higgs coupling matrix, with CKM factored out
      real*8 bij(1:6,1:6)  ! 
      real*8 apbij2(1:6,1:6)
      common/cxab2/apbij2
      real*8 anorm,bnorm
      common/higgs/anorm,bnorm
      integer i2hdmtype
      common/ci2hdmtype/i2hdmtype
      real*8 tanbeta,Ainput,Binput
      common/ctHbcoupl/tanbeta,Ainput,Binput
      real*8 xmt,twidth
      common/ctparam/xmt,twidth
      real*8 xmass(-5:21)
      common/parmass/xmass
      real*8 xmb
      include 'stpcblks.h'
c
c Number of light flavours
      nl=5
c
      if(i2hdmtype.eq.1)then
        anorm=Ainput
        bnorm=Binput
      elseif(i2hdmtype.eq.2)then
        xmb=xmass(5)
        call typeIIcoupl(xmb,xmt,tanbeta,anorm,bnorm)
      else
        write(*,*)'Error in setpar_2hdm: i2hdmtype=',i2hdmtype
        stop
      endif
c
      do i=1,6
        do j=1,6
          aij(i,j)=0.d0
          bij(i,j)=0.d0
          apbij2(i,j)=0.d0
        enddo
      enddo
c
      aij(1,2)=anorm
      aij(1,3)=anorm
      aij(1,5)=anorm  
      aij(4,2)=anorm
      aij(4,3)=anorm
      aij(4,5)=anorm
      aij(6,2)=anorm
      aij(6,3)=anorm
      aij(6,5)=anorm
c
      bij(1,2)=bnorm
      bij(1,3)=bnorm
      bij(1,5)=bnorm  
      bij(4,2)=bnorm
      bij(4,3)=bnorm
      bij(4,5)=bnorm
      bij(6,2)=bnorm
      bij(6,3)=bnorm
      bij(6,5)=bnorm

c The following only applies for t->bH decays only.
      do i=1,6
        do j=1,6
          apbij2(i,j)=aij(i,j)**2+bij(i,j)**2
        enddo
      enddo
c
      return
      end


      subroutine typeIIcoupl(xmb,xmt,tanbeta,xia,xib)
c Calculates couplings A and B in a two-Higgs doublet model,
c from eq. (1.96) of hep-ph/0503173. The masses entering the
c Yukawa coupling are run from the pole masses used as inputs.
      implicit none
      include 'stpcblks.h'
      real*8 xmb,xmb_run,xmt,xmt_run,tanbeta,xia,xib,xmu2,
     # xmu2a,mbrun,mtrun
      real*8 v
      parameter (v=246.22796) ! Higgs VeV
      integer inloscale
      common/cinloscale/inloscale
      real * 8 scptveto
      common/cscptveto/scptveto
c
c Calculate running masses. The scale is set according to what done
c in zgmu2_nlo. Needs further refinements in the case of dynamical
c scale (such choice is overruled here)
      if(mod(inloscale,10).eq.2)then
        xmu2 = xm12
      elseif(mod(inloscale,10).eq.1.or.mod(inloscale,10).eq.3)then
        xmu2 = (xm12+xm22)/2.d0
      elseif(mod(inloscale,10).eq.4)then
        xmu2 = ( (sqrt(xm12)+sqrt(xm22))/2.d0 )**2
      else
        write(*,*)'Unknown option in typeIIcoupl',inloscale
        stop
      endif
      if(inloscale.le.20)then
        xmu2a = xmu2
      else
        xmu2a = scptveto**2
      endif
      xmur2  = xmu2a*xren2
c
      xmb_run=mbrun(sqrt(xmur2),xmb,xmt)
      xmt_run=mtrun(sqrt(xmur2),xmb,xmt)

c Now calculate Yukawa couplings:
      xia=(xmb_run*tanbeta+xmt_run/tanbeta)/v/sqrt(2.)
      xib=(xmb_run*tanbeta-xmt_run/tanbeta)/v/sqrt(2.)
      
      return
      end


      function mbrun(mu,mbpole,mtpole)
c Calculates the msbar running b mass for mu>mb, as in hep-ph/9704448
      implicit none
      include 'stpcblks.h'
      real*8 x1,x2,x3,alfas,mu,mbrun,pi,mbpole,mtpole,c,c2,mbmt
      parameter (pi=3.14159265358979312D0)
c
      x1=alfas(mu**2.,xlam,nl)/pi
      x2=alfas(mbpole**2.,xlam,nl)/pi

      if(mu.ge.mtpole) then
         x3=alfas(mtpole**2.,xlam,nl)/pi
         mbmt=4.23*c2(x3)/c2(x2) ! 4.23 is mb(mb) in hep-ph/9704448
         mbrun=mbmt*c(x1)/c(x3)
      else if (mu.ge.mbpole) then
         mbrun=4.23*c2(x1)/c2(x2)
      endif

      return
      end


      function mtrun(mu,mbpole,mtpole)
c Calculates the msbar running t mass for mu>mb, as in hep-ph/9704448
      implicit none
      include 'stpcblks.h'
      real*8 x1,x2,alfas,mu,mtrun,pi,mbpole,mtpole,c,c2,mtmt
      parameter (pi=3.14159265358979312D0)
c
      x1=alfas(mu**2.,xlam,nl)/pi
      x2=alfas(mtpole**2.,xlam,nl)/pi

      if(mu.ge.mtpole) then
         mtmt=167.4 ! Value of mt(mt) according in hep-ph/99704448
         mtrun=mtmt*c(x1)/c(x2)
      else if (mu.ge.mbpole) then
         mtrun=167.4*c2(x1)/c2(x2)
      endif

      return
      end


      function c(x)
c Function entering the running b mass for mu>mt
      real * 8 c,x

      c=(7./2.*x)**(4./7.)*(1.+1.398*x+1.793*x**2)

      return
      end

      function c2(x)
c Function entering the running b mass for mb<mu<mt
      real * 8 c2,x

      c2=(23./6.*x)**(12./23.)*(1.+1.175*x+1.501*x**2)

      return
      end


      subroutine combine_coupl()
      implicit none
      integer i,j
      real*8 ckm2(1:6,1:6)
      common/cckm2/ckm2
      real*8 apbij2(1:6,1:6)
      common/cxab2/apbij2
      real*8 rucab,rccab,rtcab,rduccab,rsuccab,rbuccab
      common/ccabfct/rucab,rccab,rtcab,rduccab,rsuccab,rbuccab
      real*8 cab2(1:6,1:6)
      common/ccab2/cab2
c
      do i=1,6
        do j=1,6
          cab2(i,j)=ckm2(i,j)*apbij2(i,j)
        enddo
      enddo
c
      rucab=cab2(1,2)+cab2(1,3)+cab2(1,5)
      rccab=cab2(4,2)+cab2(4,3)+cab2(4,5)
      rtcab=cab2(6,2)+cab2(6,3)+cab2(6,5)
      rduccab=cab2(1,2)+cab2(4,2)
      rsuccab=cab2(1,3)+cab2(4,3)
      rbuccab=cab2(1,5)+cab2(4,5)
c
      return
      end


      subroutine idrfill()
c Fills the array idrmap
      implicit none
      integer idr,jproc,itt
      integer idrmap(8,1:3,2)
      common/cidrmap/idrmap
c
c t production: trivial
      itt=1
      do idr=1,8
        do jproc=1,3
          idrmap(idr,jproc,itt)=idr
        enddo
      enddo
c tbar production: charge conjugation
      itt=2
      do idr=1,8
        do jproc=1,3
          idrmap(idr,jproc,itt)=idr
        enddo
      enddo
c
      return
      end


      subroutine parsetpar()
      implicit none
      integer jproc,i,itype,itt,ichconj
      integer imapp(0:5)
      common/cimapp/imapp
      integer ivbhpro(7,3,24,2)
      common/civbhpro/ivbhpro
      integer idp1(7,3,24,2),idp2(7,3,24,2),idp3(7,3,24,2)
      common/cidpart/idp1,idp2,idp3
      integer jtypemax(3,7)
      common/cjtypemax/jtypemax
      integer ichkpid
      common/cichkpid/ichkpid
c
c imapp(i) returns the PDG id number (1=d, 2=u, 3=s, 4=c, 5=b, 21=g)
c given our id number (1=u, 2=d, 3=s, 4=c, 5=b, 0=g)
      imapp(0)=21
      imapp(1)=2
      imapp(2)=1
      imapp(3)=3
      imapp(4)=4
      imapp(5)=5
c
c ivbhpro returns the process number associated to the entries; this is
c identical to i1hpro (see the routine store_events)
      do i=1,7
        do jproc=1,3
          do itype=1,24
            do itt=1,2
              ivbhpro(i,jproc,itype,itt)=0
            enddo
          enddo
        enddo
      enddo
c
c t production
c
c jproc=1
      ivbhpro(1,1,1,1)=407
c jproc=2
      do itype=1,jtypemax(2,1)
        ivbhpro(1,2,itype,1)=401
      enddo
      do itype=1,jtypemax(2,2)
        ivbhpro(2,2,itype,1)=403
      enddo
      do itype=1,jtypemax(2,3)
        if(mod((itype+mod(itype,2))/2,2).eq.1)then
          ivbhpro(3,2,itype,1)=408
        else
          ivbhpro(3,2,itype,1)=401
        endif
      enddo
      do itype=1,jtypemax(2,4)
        if(mod((itype+mod(itype,2))/2,2).eq.1)then
          ivbhpro(4,2,itype,1)=408
        else
          ivbhpro(4,2,itype,1)=403
        endif
      enddo
      do itype=1,jtypemax(2,5)
        ivbhpro(5,2,itype,1)=401
      enddo
      do itype=1,jtypemax(2,6)
        ivbhpro(6,2,itype,1)=403
      enddo
      do itype=1,jtypemax(2,7)
        ivbhpro(7,2,itype,1)=408
      enddo
c jproc=3
      do itype=1,jtypemax(3,1)
        ivbhpro(1,3,itype,1)=402
      enddo
      do itype=1,jtypemax(3,3)
        ivbhpro(3,3,itype,1)=405
      enddo
c
c tbar production; charge conjugation of those relevant to t production,
c using the map 
c 401->403; 402->404; 403->401; 404->402; 
c 405->406; 406->405; 408->409; 409->408;
c 407->407
c
c jproc=1
      ivbhpro(1,1,1,2)=407
c jproc=2
      do itype=1,jtypemax(2,1)
        ivbhpro(1,2,itype,2)=403
      enddo
      do itype=1,jtypemax(2,2)
        ivbhpro(2,2,itype,2)=401
      enddo
      do itype=1,jtypemax(2,3)
        if(mod((itype+mod(itype,2))/2,2).eq.1)then
          ivbhpro(3,2,itype,2)=409
        else
          ivbhpro(3,2,itype,2)=403
        endif
      enddo
      do itype=1,jtypemax(2,4)
        if(mod((itype+mod(itype,2))/2,2).eq.1)then
          ivbhpro(4,2,itype,2)=409
        else
          ivbhpro(4,2,itype,2)=401
        endif
      enddo
      do itype=1,jtypemax(2,5)
        ivbhpro(5,2,itype,2)=403
      enddo
      do itype=1,jtypemax(2,6)
        ivbhpro(6,2,itype,2)=401
      enddo
      do itype=1,jtypemax(2,7)
        ivbhpro(7,2,itype,2)=409
      enddo
c jproc=3
      do itype=1,jtypemax(3,1)
        ivbhpro(1,3,itype,2)=404
      enddo
      do itype=1,jtypemax(3,3)
        ivbhpro(3,3,itype,2)=406
      enddo
c
c idpX returns the flavour of parton number X (1=coming from the left,
c 2=coming from the right, 3=FKS parton) in the process associated 
c with the entries. The labelling scheme of PDG has been used.
c For some partonic subprocesses in single-top production, the identity
c of the final-state light parton is determined statistically
c on an event-by-event basis (in ckmunwgt called by xout). In such cases,
c special non-PDG labels are used here:
c  625  -> a down-type parton; corresponding weight: rtckm-ckm2(6,3)
c  635  -> a down-type parton; corresponding weight: rtckm-ckm2(6,2)
c  6235 -> a down-type parton; corresponding weight: rtckm
      do i=1,7
        do jproc=1,3
          do itype=1,24
            do itt=1,2
              idp1(i,jproc,itype,itt)=0
              idp2(i,jproc,itype,itt)=0
              idp3(i,jproc,itype,itt)=0
            enddo
          enddo
        enddo
      enddo
c
c t production
c
c jproc=1
      idp1(1,1,1,1)=imapp(0)
      idp2(1,1,1,1)=imapp(0)
      idp3(1,1,1,1)=-6235
c jproc=2
      idp1(1,2,1,1)=imapp( 1)
      idp1(1,2,2,1)=imapp( 4)
      idp1(1,2,3,1)=imapp( 2)
      idp1(1,2,4,1)=imapp( 3)
c
      idp2(1,2,1,1)=-imapp(1)
      idp2(1,2,2,1)=-imapp(4)
      idp2(1,2,3,1)=-imapp(2)
      idp2(1,2,4,1)=-imapp(3)
c
      idp3(1,2,1,1)=-6235
      idp3(1,2,2,1)=-6235
      idp3(1,2,3,1)=-635
      idp3(1,2,4,1)=-625
c
      do itype=1,4
        idp1(2,2,itype,1)=idp2(1,2,itype,1)
        idp2(2,2,itype,1)=idp1(1,2,itype,1)
        idp3(2,2,itype,1)=idp3(1,2,itype,1)
      enddo
c
      idp1(3,2, 1,1)=imapp( 2)
      idp1(3,2, 2,1)=imapp( 2)
      idp1(3,2, 3,1)=imapp( 2)
      idp1(3,2, 4,1)=imapp( 2)
      idp1(3,2, 5,1)=imapp( 3)
      idp1(3,2, 6,1)=imapp( 3)
      idp1(3,2, 7,1)=imapp( 3)
      idp1(3,2, 8,1)=imapp( 3)
      idp1(3,2, 9,1)=imapp( 5)
      idp1(3,2,10,1)=imapp( 5)
      idp1(3,2,11,1)=imapp( 5)
      idp1(3,2,12,1)=imapp( 5)
      idp1(3,2,13,1)=imapp( 2)
      idp1(3,2,14,1)=imapp( 2)
      idp1(3,2,15,1)=imapp( 2)
      idp1(3,2,16,1)=imapp( 2)
      idp1(3,2,17,1)=imapp( 3)
      idp1(3,2,18,1)=imapp( 3)
      idp1(3,2,19,1)=imapp( 3)
      idp1(3,2,20,1)=imapp( 3)
      idp1(3,2,21,1)=imapp( 5)
      idp1(3,2,22,1)=imapp( 5)
      idp1(3,2,23,1)=imapp( 5)
      idp1(3,2,24,1)=imapp( 5)
c
      idp2(3,2, 1,1)=imapp( 1)
      idp2(3,2, 2,1)=imapp( 4)
      idp2(3,2, 3,1)=-imapp(1)
      idp2(3,2, 4,1)=-imapp(4)
      idp2(3,2, 5,1)=imapp( 1)
      idp2(3,2, 6,1)=imapp( 4)
      idp2(3,2, 7,1)=-imapp(1)
      idp2(3,2, 8,1)=-imapp(4)
      idp2(3,2, 9,1)=imapp( 1)
      idp2(3,2,10,1)=imapp( 4)
      idp2(3,2,11,1)=-imapp(1)
      idp2(3,2,12,1)=-imapp(4)
      idp2(3,2,13,1)=imapp( 3)
      idp2(3,2,14,1)=imapp( 5)
      idp2(3,2,15,1)=-imapp(3)
      idp2(3,2,16,1)=-imapp(5)
      idp2(3,2,17,1)=imapp( 2)
      idp2(3,2,18,1)=imapp( 5)
      idp2(3,2,19,1)=-imapp(2)
      idp2(3,2,20,1)=-imapp(5)
      idp2(3,2,21,1)=imapp( 2)
      idp2(3,2,22,1)=imapp( 3)
      idp2(3,2,23,1)=-imapp(2)
      idp2(3,2,24,1)=-imapp(3)
c
      do itype=1,24
        idp3(3,2,itype,1)=idp2(3,2,itype,1)
      enddo
c
      do itype=1,24
        idp1(4,2,itype,1)=idp2(3,2,itype,1)
        idp2(4,2,itype,1)=idp1(3,2,itype,1)
        idp3(4,2,itype,1)=idp3(3,2,itype,1)
      enddo
c
      idp1(5,2,1,1)=imapp( 2)
      idp1(5,2,2,1)=imapp( 3)
      idp1(5,2,3,1)=imapp( 5)
c
      idp2(5,2,1,1)=-imapp(2)
      idp2(5,2,2,1)=-imapp(3)
      idp2(5,2,3,1)=-imapp(5)
c
      do itype=1,3
        idp3(5,2,itype,1)=idp2(5,2,itype,1)
      enddo
c
      do itype=1,24
        idp1(6,2,itype,1)=idp2(5,2,itype,1)
        idp2(6,2,itype,1)=idp1(5,2,itype,1)
        idp3(6,2,itype,1)=idp3(5,2,itype,1)
      enddo
c
      idp1(7,2,1,1)=imapp( 2)
      idp1(7,2,2,1)=imapp( 3)
      idp1(7,2,3,1)=imapp( 5)
c
      idp2(7,2,1,1)=imapp( 2)
      idp2(7,2,2,1)=imapp( 3)
      idp2(7,2,3,1)=imapp( 5)
c
      do itype=1,3
        idp3(7,2,itype,1)=idp2(7,2,itype,1)
      enddo
c jproc=3
      idp1(1,3,1,1)=imapp( 2)
      idp1(1,3,2,1)=imapp( 3)
      idp1(1,3,3,1)=imapp( 5)
c
      idp2(1,3,1,1)=imapp( 0)
      idp2(1,3,2,1)=imapp( 0)
      idp2(1,3,3,1)=imapp( 0)
c
      do itype=1,3
        idp3(1,3,itype,1)=imapp( 0)
      enddo
c
      do itype=1,3
        idp1(3,3,itype,1)=idp2(1,3,itype,1)
        idp2(3,3,itype,1)=idp1(1,3,itype,1)
        idp3(3,3,itype,1)=idp3(1,3,itype,1)
      enddo
c
c tbar production; charge conjugation of those relevant to t production
c
      do i=1,7
        do jproc=1,3
          do itype=1,24
            idp1(i,jproc,itype,2)=
     #        ichconj(idp1(i,jproc,itype,1))
            idp2(i,jproc,itype,2)=
     #        ichconj(idp2(i,jproc,itype,1))
            idp3(i,jproc,itype,2)=
     #        ichconj(idp3(i,jproc,itype,1))
          enddo
        enddo
      enddo
c
      if(ichkpid.eq.0)call parcheckpar()
      return
      end


      subroutine parcheckpar()
      implicit none
      integer iallzero,iz,i,jproc,itype,ich,itt,ihpro,i1,i2,i3,i4,i5
      parameter (iallzero=0)
      parameter (iz=0)
      integer ivbhpro(7,3,24,2)
      common/civbhpro/ivbhpro
      integer idp1(7,3,24,2),idp2(7,3,24,2),idp3(7,3,24,2)
      common/cidpart/idp1,idp2,idp3
c
      do i=1,7
        do jproc=1,3
          do itype=1,24
            do itt=1,2
              ihpro=ivbhpro(i,jproc,itype,itt)
              i1=idp1(i,jproc,itype,itt)
              i2=idp2(i,jproc,itype,itt)
              i3=idp3(i,jproc,itype,itt)
              i4=6*isign(1,1-2*(itt-1))
              i5=37*isign(1,-1+2*(itt-1))
              call parcheckfin(ihpro,i1,i2,i3,i4,i5,iallzero,iz,
     #                         i,jproc,itype,ich,itt)
            enddo
          enddo
        enddo
      enddo
      return
      end


      subroutine parcheckfin(ihpro,i1,i2,i3,i4,i5,iallzero,ic,
     #                       idr,jproc,itype,ich,itt)
      implicit none
      integer ihpro,i1,i2,i3,i4,i5,iallzero,ic,idr,jproc,itype,
     # ich,itt,isum,izero
      parameter (izero=0)
      real*8 tiny,chrg,chin,chout,chall,chprdct
      parameter (tiny=1.d-8)
      logical ferror,isalqrk
c
      ferror=.false.
      if(itt.eq.1)then
        chprdct=chrg(6)
      elseif(itt.eq.2)then
        chprdct=chrg(-6)
      else
        write(*,*)'Wrong itt in parcheckfin',itt
        stop
      endif
      if(chprdct.ne.chrg(i4))ferror=.true.
      chprdct=chprdct+sign(1.d0,dfloat(i5))
      isum=abs(i1)+abs(i2)+abs(i3)
      if(isum.ne.0)then
        chin=chrg(i1)+chrg(i2)
        chout=chrg(i3)
      endif
      if(iallzero.eq.0)then
c i1,i2,i3 must be either all nonzero, or all zero
        if( ( (i1.ne.0) .and. 
     #        (i2.eq.0.or.i3.eq.0) ) .or.
     #      ( (i2.ne.0) .and. 
     #        (i1.eq.0.or.i3.eq.0) ) .or.
     #      ( (i3.ne.0) .and. 
     #        (i1.eq.0.or.i2.eq.0) ) )ferror=.true.
      elseif(iallzero.eq.1)then
c all process parameters must be different from zero
        if(i1.eq.0.or.i2.eq.0.or.
     #     i3.eq.0.or.ihpro.eq.0)ferror=.true.
      else
        write(*,*)'parcheckfin called improperly'
        stop
      endif
      if(isum.ne.0)then
c charge must be conserved
        chall=chin-chout-chprdct
        if(abs(chall).gt.tiny)ferror=.true.
c 401 is qqbar
        if( ihpro.eq.401 .and.
     #      (i1.le.0 .or. i2.ge.0 .or. 
     #      (.not.isalqrk(i3,izero))) )ferror=.true.
c 402 is qg
        if( ihpro.eq.402 .and.
     #      (i1.le.0 .or. i2.ne.21 .or. i3.ne.21) )ferror=.true.
c 403 is qbarq
        if( ihpro.eq.403 .and.
     #      (i1.ge.0 .or. i2.le.0 .or. 
     #      (.not.isalqrk(i3,izero))) )ferror=.true.
c 404 is qbarg
        if( ihpro.eq.404 .and.
     #      (i1.ge.0 .or. i2.ne.21 .or. i3.ne.21) )ferror=.true.
c 405 is gq
        if( ihpro.eq.405 .and.
     #      (i1.ne.21 .or. i2.le.0 .or. i3.ne.21) )ferror=.true.
c 406 is gqbar
        if( ihpro.eq.406 .and.
     #      (i1.ne.21 .or. i2.ge.0 .or. i3.ne.21) )ferror=.true.
c 407 is gg
        if( ihpro.eq.407 .and.
     #      (i1.ne.21 .or. i2.ne.21 .or. 
     #      (.not.isalqrk(i3,izero))) )ferror=.true.
c 408 is qq
        if( ihpro.eq.408 .and.
     #      (i1.le.0 .or. i2.le.0 .or. 
     #      (.not.isalqrk(i3,izero))) )ferror=.true.
c 409 is qbarqbar
        if( ihpro.eq.409 .and.
     #      (i1.ge.0 .or. i2.ge.0 .or. 
     #      (.not.isalqrk(i3,izero))) )ferror=.true.
      endif
      if(ferror)then
        write(*,*)'Error in parcheckfin'
        write(*,*)'ihpro,i1,i2,i3,i4,i5:',ihpro,i1,i2,i3,i4,i5
        write(*,*)'idr,jproc,itype,ich,itt:',idr,jproc,itype,ich,itt
        write(*,*)'chin,chout,chprdct,chall:',chin,chout,chprdct,chall
        write(*,*)'crossing:',ic
        stop
      endif
      return
      end


      function chrg(id)
      implicit none
      real*8 chrg,tmp
      integer id,ia
      real*8 chup,chdn
      parameter (chup=2.d0/3.d0)
      parameter (chdn=-1.d0/3.d0)
c
      ia=abs(id)
      if(ia.eq.1.or.ia.eq.3.or.ia.eq.5)then
        tmp=chdn
      elseif(ia.eq.2.or.ia.eq.4.or.ia.eq.6)then
        tmp=chup
      elseif(ia.eq.21)then
        tmp=0.d0
      elseif(ia.eq.625.or.ia.eq.635.or.ia.eq.6235)then
        tmp=chdn
      else
        write(*,*)'Error in chrg: id=',id
        stop
      endif
      chrg=sign(1.d0,dfloat(id))*tmp
      return
      end


      function ichconj(ip)
c Charge conjugation
      implicit none
      integer ichconj,ip,itmp
c
      if(abs(ip).ne.21)then
        itmp=-ip
      else
        itmp=ip
      endif
      ichconj=itmp
      return
      end


      function isalqrk(ipart,iflag)
c Returns true is ipart is ID of a light quark
      implicit none
      logical isalqrk,ltmp
      integer ipart,iflag
c
      if(iflag.eq.0)then
        ltmp=( abs(ipart).le.5 .or. abs(ipart).eq.6235 .or.
     #         abs(ipart).eq.625 .or. abs(ipart).eq.635 ) .and. 
     #       ipart.ne.21
      else
        ltmp=abs(ipart).le.5 .and. ipart.ne.21
      endif
      isalqrk=ltmp
      return
      end


      subroutine getnloiproc(iprdct0hw)
c Converts the MC@NLO process codes for single top production into the codes
c used in the NLO computation. MC@NLO conventions are
c Process: iprdct=2040+IT   Ht channel      
c with
c                 IT=0  t+tbar production  
c                 IT=1  tbar production    
c                 IT=4  t production       
c This routine in meant to be called after setting 
c iprdct0hw=mod(-iprdct0hw,10000). Furthermore
c The NLO conventions are
c     ich=3    -> Ht-channel
c     ittbar=1 -> t production
c     ittbar=2 -> tbar production
c and this routines sets the ranges for ich and ittbar
c 
      implicit none
      integer iprdct0hw,itmp
      integer ichmin,ichmax
      common/cichrange/ichmin,ichmax
      integer ittmin,ittmax
      common/cittrange/ittmin,ittmax
c
      if(iprdct0hw.ge.2040.and.iprdct0hw.le.2049)then
        ichmin=3
        ichmax=3
        itmp=iprdct0hw-2040
      else
        write(*,*)'getnloiproc: wrong process number',iprdct0hw
        stop
      endif
      if(itmp.eq.0)then
        ittmin=1
        ittmax=2
      elseif(itmp.eq.1)then
        ittmin=2
        ittmax=2
      elseif(itmp.eq.4)then
        ittmin=1
        ittmax=1
      else
        write(*,*)'getnloiproc: wrong process number',iprdct0hw
        stop
      endif
      return
      end


      subroutine setpardec()
      implicit none
      include 'stpcblks.h'
      real * 8 pi,one,zero,xme,xmmu,xmtau,ze2_dec,xalfaem,
     # topdecw,xmt2,xia,xib,tmpmss(3)
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      parameter (zero=0.d0)
c Values from PDG 2003
      parameter (xme=0.510998902d-3)
      parameter (xmmu=105.6583568d-3)
      parameter (xmtau=1776.99d-3)
      real*8 xmw,gaw
      common/cwparam/xmw,gaw
      real*8 xmt,twidth
      common/ctparam/xmt,twidth
      real*8 xlep1mass(2),xlep2mass(2)
      common/clepmass/xlep1mass,xlep2mass
      real*8 pdglepmss(11:16)
      common/cpdglepmss/pdglepmss
      real*8 xm012,ga1,bw1delf,bw1fmmn,xm1low2,xm1upp2
      common/cbw1/xm012,ga1,bw1delf,bw1fmmn,xm1low2,xm1upp2
      real*8 brrtop1,brrv2msb
      common/brratios/brrtop1,brrv2msb
      real*8 xbrrtoplep,xbrrtophad
      common/xibrratios/xbrrtoplep,xbrrtophad
      real*8 xbrrwlep,xbrrwhad
      common/xwibrratios/xbrrwlep,xbrrwhad
      real*8 frac12,frac123
      common/cfracs123/frac12,frac123
      real*8 wfrac12,wfrac123
      common/cwfracs123/wfrac12,wfrac123
      real*8 cthw2,sthw2
      common/cweinan/sthw2,cthw2
      real*8 sthw20
      common/csthw20/sthw20
      real*8 ckm2(1:6,1:6)
      common/cckm2/ckm2
      real*8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      real*8 anorm,bnorm
      common/higgs/anorm,bnorm
      integer idec
      common/cidec/idec
      integer iwidth
      common/ciwidth/iwidth
      integer il1hw,il2hw
      common/cilhw/il1hw,il2hw
      integer inonbtop
      common/cinonbtop/inonbtop
c Parton identities: if top doesn't decay, 4 is top, 3 and 5 outgoing
c light parton and Higgs respectively. If top decays, (3,5)-->(3,4), and 
c top decay products are 5(l), 6(nu), and 7(b) for IL=1,2,3 and inonbtop=0; 
c when IL=0,4..7 and/or inonbtop=1, parton ids are trivially derived from 
c these. PDG conventions are used
      integer ip1s,ip2s,ip3s,ip4s,ip5s,ip6s,ip7s
      common/ci1parts/ip1s,ip2s,ip3s,ip4s,ip5s,ip6s,ip7s
c PDG codes for charged leptons and neutrinos for a given IL (NLO) code;
c the particle code (not the antiparticle) is entered here
c Charged lepton from W decay
      integer ichlw(0:6)
      data ichlw/0,11,13,15,0,0,0/
c Neutrino from W decay
      integer ineuw(0:6)
      data ineuw/0,12,14,16,0,0,0/
c
      if(idec.eq.1)then
        write(*,*)'Routine setpardec should not be called'
        stop
      endif
      if(sthw20.ne.sthw2)then
        write(*,*)'Inconsistency between setpar and reset_twdbr'
        stop
      endif
c Electron charge squared (computed at the W mass, used for decays)
      xmt2=xmt**2
      ze2_dec = 4*pi*xalfaem(xmw2)
c Lepton masses and identities: xlep#mass(i) is the mass of the decay 
c product # in the decay of the top (i=1, in which case #=1 ==> antiparticle,
c #=2 ==> particle) or of the W (i=2, in which case #=1 ==> particle,
c #=2 ==> antiparticle). If the top and W have a given leptonic decay, the
c masses are set here, but are in any case overwritten by the routine
c getpdecids(), in order to treat all the cases on equal footing.
c We assume here that only (top,W-) are produced, and change signs on 
c event-by-event basis if the event is (tbar,W+) (see getpdecids)
      pdglepmss(11)=xme
      pdglepmss(12)=0.d0
      pdglepmss(13)=xmmu
      pdglepmss(14)=0.d0
      pdglepmss(15)=xmtau
      pdglepmss(16)=0.d0
      tmpmss(1)=xme
      tmpmss(2)=xmmu
      tmpmss(3)=xmtau
      ip5s=-ichlw(il1hw)
      ip6s=ineuw(il1hw)
      ip7s=5
      if(il1hw.ge.1.and.il1hw.le.3)then
        xlep1mass(1)=tmpmss(il1hw)
        xlep2mass(1)=0.d0
      elseif(il1hw.eq.0.or.(il1hw.ge.4.and.il1hw.le.6))then
        xlep1mass(1)=-1.d0
        xlep2mass(1)=-1.d0
      else
        write(*,*)'Error #1 in setpardec: inconsistent entries',il1hw
        stop
      endif
c Insert here identities for Higgs decay if need be
c Fills MadEvent common blocks. Set positron charge and QCD coupling g 
c equal to one, and use the actual values in the main code
      xia=anorm
      xib=bnorm
      call setmeFRpar(xmw,gaw,zero,zero,
     #                xmt,twidth,xm2,zero,
     #                zero,sthw2,one,one,xia,xib)
c Compute reweight factors for event weights, that take into
c account the values of branching ratios.
c W from top decay
      if(iwidth.eq.1)then
c Correction factor for effective W mass range
        brrtop1=topdecw(xmt,xmw,gaw,xm1low2,xm1upp2,sthw2)/
     #          topdecw(xmt,xmw,gaw,zero,xmt2,sthw2)
      else
c W is on shell: the event weight is proportional to dBR_t/dQ_W^2 
c (at Q_W^2=MW^2); use dGamma_t/dQ_W^2 = Gamma_t/(pi*MW*Gamma_W)
        brrtop1=1.d0/(pi*xmw*gaw)
      endif
c Now include the branching ratios
      if(il1hw.eq.0)then
        brrtop1=(3*xbrrtoplep+2*xbrrtophad) * brrtop1
      elseif(il1hw.ge.1.and.il1hw.le.3)then
        brrtop1=xbrrtoplep * brrtop1
      elseif(il1hw.eq.4)then
        brrtop1=2*xbrrtoplep * brrtop1
      elseif(il1hw.eq.5)then
        brrtop1=2*xbrrtophad * brrtop1
      elseif(il1hw.eq.6)then
        brrtop1=(2*xbrrtoplep+2*xbrrtophad) * brrtop1
      endif
c Correct reweight factor for top decays if only t->Wb decays
c are selected
      if(inonbtop.eq.0)then
        brrtop1=ckm2(6,5)/rtckm * brrtop1
      elseif(inonbtop.ne.1)then
        write(*,*)'Unknown inonbtop in setpar',inonbtop
        stop
      endif
c frac12 is the fraction of decays W->e+mu/W->e+mu+all quarks
      frac12=2*xbrrtoplep/(2*xbrrtoplep+2*xbrrtophad)
c frac123 is the fraction of decays W->e+mu+tau/W->e+mu+tau+all quarks
      frac123=3*xbrrtoplep/(3*xbrrtoplep+2*xbrrtophad)
c Higgs: no decay performed here
      brrv2msb=1.d0
c
      if(brrtop1.eq.0.d0.or.brrv2msb.eq.0.d0)then
        write(*,*)
     #    'These decay channels will return a zero cross section'
        stop
      endif
c
      return
      end


      subroutine reset_twdbr2(xmt,twidth,xmw,gaw,
     #                        xbrrtoplep,xbrrtophad,
     #                        xbrrwlep,xbrrwhad,
     #                        gammax1,gammax2)
c If twidth>0 and gaw>0, use xbrrtoplep and xbrrtophad as branching ratios
c for top decays, and xbrrwlep,xbrrwhad as branching ratios for W decays.
c If twidth<0 and gaw<0, compute the total widths according to LO SM, and
c redefine xbrrtoplep, xbrrtophad, xbrrwlep and xbrrwhad accordingly.
c The width and branching ratios assume lepton and hadron universality;
c no CKM factors are included here
      implicit none
      real*8 xmt,twidth,xmw,gaw,xbrrtoplep,xbrrtophad,
     # xbrrwlep,xbrrwhad,gammax1,gammax2
      real*8 pi,zero,xmt2,xmw2,alfaem,xalfaem,ze2,wdtwmsb,brtop2,
     # topdecw,brtopw,brtop0
      parameter (pi=3.14159265358979312D0)
      parameter (zero=0.d0)
      real*8 sthw20
      common/csthw20/sthw20
c
c sthw20 must have the same value as sthw2 in setpar. If this is not
c the case, the code stops in setpar
      sthw20=0.23113d0
      if(twidth.gt.0.d0.and.gaw.gt.0.d0)then
        if( xbrrtoplep.lt.0.d0.or.xbrrtophad.lt.0.d0 .or.
     #      (xbrrtoplep.eq.0.d0.and.xbrrtophad.eq.0.d0) )then
          write(*,*)'Error in reset_twdbr2: negative or zero top BRs',
     #              xbrrtoplep,xbrrtophad
          stop
        endif
        if( xbrrwlep.lt.0.d0.or.xbrrwhad.lt.0.d0 .or.
     #      (xbrrwlep.eq.0.d0.and.xbrrwhad.eq.0.d0) )then
          write(*,*)'Error in reset_twdbr2: negative or zero W BRs',
     #              xbrrwlep,xbrrwhad
          stop
        endif
      elseif(twidth.lt.0.d0.and.gaw.lt.0.d0)then
        xmt2=xmt**2
        xmw2=xmw**2
        alfaem=xalfaem(xmw2)
        ze2=4*pi*alfaem
c
        xbrrtoplep=1/9.d0
        xbrrtophad=1/3.d0
        xbrrwlep=1/9.d0
        xbrrwhad=1/3.d0
c hard W first; use MSbar expression for partial width
        wdtwmsb=xmw*ze2/(4*pi*12*sthw20)
        gaw=wdtwmsb*9.d0
        if(gammax1.ne.0)then
c topdecw returns dGamma(t->blnu)/dq^2 integrated over q2 in the range
c (xm1low2,xm1upp2), up to a factor e^4*|Vtb|^2
          brtop2=topdecw(xmt,xmw,gaw,zero,xmt2,sthw20)
          twidth=brtop2*ze2**2*9.d0
        elseif(gammax1.eq.0)then
c W's are at the pole mass
c brtopw is dGamma(t->blnu)/dq^2, up to a factor gw^4*|Vtb|^2/2
          brtopw=(xmt2-xmw2)**2*(xmt2+2*xmw2)/(6144*pi**3*xmt**3)
c Multiply brtopw by the Breit Wigner for the W, in the narrow width
c approximation BR(W)->pi/(MW*GammaW). A factor 2/sin^4\theta_W is inserted
c in order to eliminate the factor 1/2 in brtopw, and convert gw^4 into e^4
          brtop0=brtopw* 2*pi/(sthw20**2 * xmw*gaw)
          twidth=brtop0*ze2**2*9.d0
        else
          write(*,*)'Error #1 in reset_twdbr2',gammax1,gammax2
          stop
        endif
      else
        write(*,*)'Error in reset_twdbr2: top or W width equal to zero'
        stop
      endif
      if( (3*xbrrtoplep+2*xbrrtophad).gt.1.0001d0)then
        write(*,*)'Error #2 in reset_twdbr2',xbrrtoplep,xbrrtophad,
     #            3*xbrrtoplep+2*xbrrtophad
        stop
      endif
      if( (3*xbrrwlep+2*xbrrwhad).gt.1.0001d0)then
        write(*,*)'Error #3 in reset_twdbr2',xbrrwlep,xbrrwhad,
     #            3*xbrrwlep+2*xbrrwhad
        stop
      endif
      return
      end


      function topdecw(xmt,xmw,wwidth,xmw2low,xmw2upp,sthw2)
c Returns top decay width integrated over W virtuality, as computed
c in topwidth.m. Insert a factor e^4*|Vtb|^2 for the correct normalization
      implicit none
      real*8 topdecw,xmt,xmw,wwidth,xmw2low,xmw2upp,sthw2
      real*8 pi,norm,tmp
      parameter (pi=3.1415926535897932d0)
c
      norm=1/(3072*pi**3*wwidth*xmw*sthw2**2*xmt**3)
      tmp=( xmt**6-6*wwidth**2*xmw**4+2*xmw**6+
     #      3*xmt**2*xmw**2*(wwidth**2-xmw**2) )*
     #    ( atan((xmw2upp-xmw**2)/(wwidth*xmw)) -
     #      atan((xmw2low-xmw**2)/(wwidth*xmw)) ) +
     #    wwidth*xmw* (
     #      (xmw2upp-xmw2low)*(4*xmw**2+xmw2low+xmw2upp-3*xmt**2)+
     #      xmw**2*(3*xmt**2+wwidth**2-3*xmw**2)*(
     #        log(wwidth**2*xmw**2 + (xmw2low - xmw**2)**2) - 
     #        log(wwidth**2*xmw**2 + (xmw2upp - xmw**2)**2) ) )
      topdecw=norm*tmp
      return
      end
c
c
c End initialization
c
c
c
c
c Begin of event file utilities
c
c
      subroutine whichone(iseed,itot,mx_of_evta,mx_of_evtb,iunit)
c Determines the type of event at random
      implicit none
      integer iseed,itot,mx_of_evta,mx_of_evtb,iunit,i0
      real*8 xpa,xpb,tiny,one,xsum,rnd,fk88random,prob
      parameter (tiny=1.d-4)
      logical flag
c
      if(itot.le.0)then
        write(6,*)'Fatal error #1 in whichone'
        stop
      endif
      xpa=dfloat(mx_of_evta)/dfloat(itot)
      xpb=dfloat(mx_of_evtb)/dfloat(itot)
      one=xpa+xpb
      if(abs(one-1.d0).gt.tiny)then
        write(6,*)'Error #1 in whichone: probability not normalized'
        stop
      endif
      i0=0
      flag=.true.
      xsum=0.d0
      rnd=fk88random(iseed)
      do while(flag)
        if(i0.gt.2)then
          write(6,*)'Fatal error #2 in whichone'
          stop
        endif
        i0=i0+1
        prob=xpa
        if(i0.gt.1)prob=xpb
        xsum=xsum+prob
        if(rnd.lt.xsum)then
          flag=.false.
          itot=itot-1
          if(i0.le.1)then
            mx_of_evta=mx_of_evta-1
          else
            mx_of_evtb=mx_of_evtb-1
          endif
          iunit=20+i0
        endif
      enddo
      return
      end


      subroutine crosscheck(itot,mx_of_evta,mx_of_evtb)
c Checks whether whichone did it right
      implicit none
      integer itot,mx_of_evta,mx_of_evtb
c
      if(itot.ne.0)then
        write(6,*)'Error: itot=',itot
        stop
      endif
      if(mx_of_evta.ne.0)then
        write(6,*)'Error: mx_of_evta=',mx_of_evta
        stop
      endif
      if(mx_of_evtb.ne.0)then
        write(6,*)'Error: mx_of_evtb=',mx_of_evtb
        stop
      endif
      return
      end


      subroutine retrieve_events(iunit,ii,dummy)
c Reads from disk the complete information on the events; see store_events
c for the conventions used
      implicit none
      integer iunit,ii,i,j
      real*8 dummy
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
      integer iccode
      common/ciccode/iccode
      integer idec
      common/cidec/idec
      integer np
      common/cnp/np
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 emsca
      common/cemsca/emsca
      real*8 xmom_lb(10,4)
      common/cxmomlb/xmom_lb
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
c
      read(iunit,901,end=997,err=998)i1hpro,iccode,np
      if(idec.eq.0)then
        read(iunit,902,end=997,err=998)ip1,ip2,ip3,
     #                                 ip4,ip5,ip6,ip7
        read(iunit,903,end=997,err=998)xevsign,emsca
        read(iunit,904,end=997,err=998)((xmom_lb(i,j),j=1,4),i=1,3),
     #                                 ((xmom_lb(i,j),j=1,4),i=5,8)
      elseif(idec.eq.1)then
        read(iunit,902,end=997,err=998)ip1,ip2,ip3,ip4,ip5
        read(iunit,903,end=997,err=998)xevsign,emsca
        read(iunit,904,end=997,err=998)((xmom_lb(i,j),j=1,4),i=1,5)
      endif
      read(iunit,905,end=997,err=998) ux1,ux2,uq2
      goto 999
 901  format(1x,i3,2(1x,i2))
 902  format(8(1x,i3))
 903  format(2(1x,d14.8))
 904  format(32(1x,d14.8))
 905  format(3(1x,d14.8))
 997  write(*,*)'unexpected end of file, iunit=',iunit
      stop
 998  write(*,*)'format error'
      write(77,*)'event #:',ii
      write(77,901)i1hpro,iccode,np
      write(77,902)ip1,ip2,ip3,ip4,ip5,ip6,ip7
      write(77,903)xevsign,emsca
      write(77,904)((xmom_lb(i,j),j=1,4),i=1,10)
      stop
 999  continue
      return
      end


      subroutine store_events(iunit,xpmone)
c Stores on disk the complete information on the events. Starting
c from version 3.1, each event has the following format:
c       IPR, IC, NP
c      (ID(I),I=1,NP)
c      ((P(J,I),J=1,4),I=1,NP)
c where IPR is the subprocess code (i1hpro), IC is the colour code
c (iccode, NON trivial here), NP is the number of partons entering the 
c reaction (thus, this includes the soft parton in the case of S events),
c ID(I) are the particle identities (ip1,...,ip9 here), and P(J,I) are 
c the particles four momenta in the lab frame (P(J,I)=xmom_lb(i,j) here).
c
c This routine is called with xpmone=1 when events are obtained from
c SPRING, and with xpmone=-1 after the events are read from the temporary
c files (via retrieve_events), to be stored in the final event file.
c When xpmone=1, one has xevsign=+1/-1, and the weight of the event is 
c xevsign*wgt[a,b]ev. When xpmone=-1, then xevsign is the weight of the event. 
c
c i1hpro has the following conventions:
c   i1hpro         process
c    401        q qbar    -> X
c    402        q g       -> X
c    403        qbar q    -> X
c    404        qbar g    -> X
c    405        g q       -> X
c    406        g qbar    -> X
c    407        g g       -> X
c    408        q q       -> X
c    409        qbar qbar -> X
c X being the top, the W and the light parton. Note that the light parton in
c the final state has typically a different flavour wrt the one in other 
c processes with the same i1hpro; i1hpro is however used here only in the
c context of the NLO computation.
c ipX is the parton code relevant to parton # X. PDG conventions are
c used: 1=d, 2=u, 3=s, 4=c, 5=b, 21=g
      implicit none
      integer iunit,i,j
      real*8 xpmone,xevwgt,xfact,brfact,dummy
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3,ip4,ip5,ip6,ip7
      common/ci1part/ip1,ip2,ip3,ip4,ip5,ip6,ip7
      integer iccode
      common/ciccode/iccode
      integer idec
      common/cidec/idec
      integer np
      common/cnp/np
      integer ievffmt
      common/cievffmt/ievffmt
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 emsca
      common/cemsca/emsca
c xmom_lb(i,j) is the j component of the four vector of the particle # i,
c given in the laboratory frame. j=4 is the energy for MC@NLO versions
c up to 2.31, the mass for version 3.1 onwards. i=1,2 are the incoming
c partons, 3 is the outgoing parton, 4 is t, 5 is H. When the top decays, 
c 6=l+, 7=nu, 8=b are the decay products of the top. Momentum conservation is 
c (1+2)-(3+4+5)=0 or (1+2)-(3+5+6+7+8)=0
c The order in which momenta and parton identities are stored is
c 1,2,3,4,5 and 1,2,3,5,6,7,8 for the undecayed and decayed cases respectively
      real*8 xmom_pass(9,4)
      integer IPS(9)
      real*8 xmom_lb(10,4)
      common/cxmomlb/xmom_lb
      integer iwgtnorm
      common/ciwgtnorm/iwgtnorm
      real*8 wgtaev,wgtbev
      common/cwgtev/wgtaev,wgtbev
      real*8 wgtmax
      common/cwgtmax/wgtmax
c Reweight factors that include branching ratios, to be inserted in 
c the case of decayed top and W
      real*8 brrtop1,brrv2msb
      common/brratios/brrtop1,brrv2msb
c PDF stuff
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
C--Les Houches Common Blocks
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
c
      if(xpmone.eq.-1)then
c Events are already stored in temporary files, and are passed to this
c routines through common blocks filled by retrieve_events
        xevwgt=xevsign
        xfact=1.d0
      elseif(xpmone.eq.1)then
c Events are obtained from SPRING, and are written to temporary files
c for the first time
        if(idec.eq.0)then
          np=7
          brfact=brrtop1*brrv2msb
        elseif(idec.eq.1)then
          np=5
          brfact=1.d0
        else
          write(6,*) 'Error in store_events: idec=',idec
          stop
        endif
        if(xmom_lb(3,4).eq.0.d0)then
          xevwgt=xevsign*wgtbev*brfact
        else
          xevwgt=xevsign*wgtaev*brfact
        endif
        if(abs(xevwgt).gt.wgtmax)wgtmax=abs(xevwgt)
        xfact=1.d0
      else
        write(*,*)'Fatal error in store_events: xpmone=',xpmone
        stop
      endif
      if(ievffmt.eq.0.or.xpmone.eq.1)then
        write(iunit,901)i1hpro,iccode,np
        if(idec.eq.0)then
          write(iunit,902)ip1,ip2,ip3,ip4,ip5,ip6,ip7
          write(iunit,903)xevwgt,emsca
          write(iunit,904)((xfact*xmom_lb(i,j),j=1,4),i=1,3),
     #                    ((xfact*xmom_lb(i,j),j=1,4),i=5,8)
        elseif(idec.eq.1)then
          write(iunit,902)ip1,ip2,ip3,ip4,ip5
          write(iunit,903)xevwgt,emsca
          write(iunit,904)((xmom_lb(i,j),j=1,4),i=1,5)
        endif
        write(iunit,905) ux1,ux2,uq2
      elseif(ievffmt.eq.1.and.xpmone.eq.-1)then
c The variable dummy can be replaced by xscale if need be
         dummy=0.d0
         if(idec.eq.0) then
            IPS(1)=ip1
            IPS(2)=ip2
            IPS(3)=ip3
            IPS(4)=ip4
            IPS(5)=ip5
            IPS(6)=ip6
            IPS(7)=ip7
            IPS(8)=0
            IPS(9)=0
            DO I=1,3
               DO J=1,4
                  xmom_pass(I,J)=xfact*xmom_lb(I,J)
               ENDDO
            ENDDO
            DO I=5,8
               DO J=1,4
                  xmom_pass(I-1,J)=xfact*xmom_lb(I,J)
               ENDDO
            ENDDO            
            DO J=1,4
               xmom_pass(8,J)=0
               xmom_pass(9,J)=0
            ENDDO
         elseif(idec.eq.1) then
            IPS(1)=ip1
            IPS(2)=ip2
            IPS(3)=ip3
            IPS(4)=ip4
            IPS(5)=ip5
            DO I=1,5
               DO J=1,4
                  xmom_pass(I,J)=xmom_lb(I,J)
               ENDDO
            ENDDO
            DO I=6,9
               IPS(I)=0
               DO J=1,4
                  xmom_pass(I,J)=0
               ENDDO
            ENDDO
         else
            write(*,*)'Error in store_events:',idec
            stop
         endif
         call write_lhef_event(iunit,idec,
     #        i1hpro,iccode,np,IPS,
     #        xevwgt,emsca,
     #        xmom_pass)
      else
         write(*,*)'Error in store_events: unknown file format',ievffmt
         stop
      endif
 901  format(1x,i3,2(1x,i2))
 902  format(8(1x,i3))
 903  format(2(1x,d14.8))
 904  format(32(1x,d14.8))
 905  format(3(1x,d14.8))
      return
      end


      subroutine getipdgpart(part,ipdg)
c Converts particle naming conventions, for Herwig to PDG
      implicit none
      character * 4 part
      integer ipdg
c
      ipdg=-100
      if(part.eq.'P   ')then
        ipdg=2212
      elseif(part.eq.'PBAR')then
        ipdg=-2212
      elseif(part.eq.'N   ')then
        ipdg=2112
      elseif(part.eq.'NBAR')then
        ipdg=-2112
      elseif(part.eq.'GAMA')then
        ipdg=22
      elseif(part.eq.'E-  ')then
        ipdg=11
      else
        write(*,*)'Error in getipdgpart'
        write(*,*)'No such particle in MLM:',part
        stop
      endif
      return
      end
c
c
c End of event file utilities
c
c
c
c
c Begin of MC subtraction terms
c
c
      subroutine xmcsubthtpp(jproc,z1,z2,xm12,xm22,s,x,yi,cth1,cth2,
     #                       gfactsf,gfactcl,flagmc,lzone,z,xmatmc)
c Computes the MC counterterms for single-top production, Ht channel. 
c The outputs of this routine are
c  gfactsf       -> the value of G_soft
c  gfactcl       -> dummy; kept for backward compatibility
c  flagmc        -> .true. if at least one lzone()=.true.
c  lzone(7,3,6)  -> .true. if the corresponding contribution to the MC
c                   counterterms is non zero, .false. otherwise
c  z(7,3,6)      -> Herwig shower variable z
c  xmatmc(7,3,6) -> the MC counterterm proper, defined as in the QQ paper
c For all the array outputs, (7,3,6) correspond to (idr,ileg,ie0sq)
      implicit none
      integer jproc
      real*8 z1,z2,xm12,xm22,s,x,yi,cth1,cth2,gfactsf,gfactcl,
     # z(7,3,6),xmatmc(7,3,6)
      logical flagmc,lzone(7,3,6)
      real*8 tiny,vcf,vca,zero,two,ff1,alpha
      parameter (tiny=1.d-4)
      parameter (vcf=4.d0/3.d0)
      parameter (vca=3.d0)
      parameter (zero=0.d0)
      parameter (two=2.d0)
      parameter (ff1=0.115d0)
      parameter (alpha=1.d0)
      integer jht,ithree
      parameter (jht=3)
      parameter (ithree=3)
      character*2 str
      parameter (str='p1')
      real*8 sx,si,tki,uki,q1qi,q2qi,so,tko,uko,q1qo,q2qo,ss,tk,uk,
     # q1q,q2q,e0sq,de0sqdx,de0sqdc,xz,zherppht,xxkp,xqherppht,
     # xalsfht,xbesfht,xalclht,xbeclht,xalazi,xbeazi,
     # gfunsoftht,gfactazi,gfunazi,sbar,tbar,tmp1,tmp,
     # ap_kern,ap_kern_qc,qin_kern,xfact,ap,
     # xjac,xjac_zqttoxyht,qk,dcosdcth1,
     # xmHt2,ptmin,ptmax,rrnd,fk88random,emscainv,emsca_bare,
     # pthw,emscafun,wsum,w1sum,
     # xinv(5),ae0sq(6),ade0sqdx(6),ade0sqdc(6),
     # x2to2(3,6),dx2to2dx(3,6),dx2to2dc(3,6),
     # xmcz(3,6),xmcxkp(3,6),xmce0sq(3,6),
     # xborn(7),xazi(7),ptresc(3,6),emscwgt(3,6),emscav(3,6)
      integer ilegmax,ileg,ie0sq,icllazi,i,icode,iborn,
     # iwsum,ileg0,ie0sq0
      logical flagxs(3),flxsec(3,6)
      real*8 xicut,deltai,deltao
      common/parsub/xicut,deltai,deltao
      real*8 alsfht(3),besfht(3)
      common/cgfunsfht/alsfht,besfht
      real*8 alclht(3),beclht(3)
      common/cgfunclht/alclht,beclht
      real*8 alazi(3),beazi(3)
      common/cgfunazi/alazi,beazi
      real*8 emsca_prc(3)
      common/cemscaprc/emsca_prc
      real*8 emwgt_prc(3)
      common/cemwgtprc/emwgt_prc
      real*8 sh
      common/shadr/sh
      integer ialwsplit(3,3,6)
      common/cialwsplit/ialwsplit
      integer idrmax(3,3)
      common/cidrmax/idrmax
      integer ie0ht(3,7,3,6)
      common/cie0ht/ie0ht
      integer imcmult(3,7,3)
      common/cimcmult/imcmult
      integer iapcp(3,7),iapcm(3,7)
      common/ciap/iapcp,iapcm
      integer idrlimcp(3,3,8),idrlimcm(3,3,8)
      common/cidrlims/idrlimcp,idrlimcm
      integer ifuntype
      common/cifuntype/ifuntype
      integer ifk88seed
      common/cifk88seed/ifk88seed
c
      sx=s*x
      flagxs(1)=.true.
      flagxs(2)=.true.
      flagxs(3)=.false.
      ilegmax=2
      if(jproc.eq.3)ilegmax=3
      do ileg=1,3
        do i=1,idrmax(jproc,jht)
          do ie0sq=1,6
            xmatmc(i,ileg,ie0sq)=0.d0
            lzone(i,ileg,ie0sq)=.false.
            z(i,ileg,ie0sq)=-1.d0
          enddo
        enddo
      enddo
c ptmin and ptmax define the veto region
      xmHt2=s*x/2.d0
      ptmin=max(ff1*sqrt(xmHt2),10.d0)
      ptmax=max(sqrt(xmHt2),ptmin+20.d0)
      rrnd=fk88random(ifk88seed)
      rrnd=emscainv(rrnd,alpha)
      emsca_bare=ptmin+rrnd*(ptmax-ptmin)
c Compute the invariants
      si=s
      call invar_in(xm12,xm22,s,x,yi,cth1,cth2,str,
     #              tki,uki,q1qi,q2qi,xinv)
      if(ifuntype.eq.1)then
        so=si
        tko=tki
        uko=uki
        q1qo=q1qi
        q2qo=q2qi
        if(jproc.eq.3)flagxs(3)=.true.
      elseif(ifuntype.eq.2.and.jproc.eq.3)then
        if((sx*x).gt.((sqrt(xm12)+sqrt(xm22))**2))then
          call invar_in(xm12,xm22,sx,x,yi,cth1,cth2,str,
     #                  tko,uko,q1qo,q2qo,xinv)
          so=sx
          flagxs(3)=.true.
        endif
      endif
c Generate the 2-->2 invariants; those for ileg=2 are identical to the
c corresponding ones for ileg=1, so we skip the generation in such a case
      do ileg=1,ilegmax
        if( ileg.eq.1 .or. (ileg.eq.3.and.jproc.eq.3) .and.
     #      flagxs(ileg) )then
          call xinvtoinv(ileg,
     #      si,tki,uki,q1qi,q2qi,
     #      so,tko,uko,q1qo,q2qo,
     #      ss,tk,uk,q1q,q2q)
          call gete0sq(jht,ileg,z1,z2,xm12,xm22,ss,x,yi,cth1,cth2,
     #                 two,zero,zero,tk,uk,q1q,q2q,
     #                 ae0sq,ade0sqdx,ade0sqdc)
          do ie0sq=1,6
            x2to2(ileg,ie0sq)=ae0sq(ie0sq)
            dx2to2dx(ileg,ie0sq)=ade0sqdx(ie0sq)
            dx2to2dc(ileg,ie0sq)=ade0sqdc(ie0sq)
          enddo
        endif
      enddo
      do ie0sq=1,6
        x2to2(2,ie0sq)=x2to2(1,ie0sq)
        dx2to2dx(2,ie0sq)=dx2to2dx(1,ie0sq)
        dx2to2dc(2,ie0sq)=dx2to2dc(1,ie0sq)
      enddo
c Now get z and xi for possible splittings
      flagmc=.false.
      do ileg=1,ilegmax
        call xinvtoinv(ileg,
     #    si,tki,uki,q1qi,q2qi,
     #    so,tko,uko,q1qo,q2qo,
     #    ss,tk,uk,q1q,q2q)
        do ie0sq=1,6
          if(ialwsplit(jproc,ileg,ie0sq).eq.1.and.flagxs(ileg))then
            e0sq=x2to2(ileg,ie0sq)
            de0sqdx=dx2to2dx(ileg,ie0sq)
            de0sqdc=dx2to2dc(ileg,ie0sq)
            xz=zherppht(ileg,e0sq,xm12,xm22,ss,x,yi,cth1,cth2,
     #                  tk,uk,q1q,q2q)
            xxkp=xqherppht(ileg,e0sq,xm12,xm22,ss,x,yi,cth1,cth2,
     #                    tk,uk,q1q,q2q)
            if(ileg.le.2)then
              if(ileg.eq.1.and.ie0sq.gt.3)then
                write(*,*)'Error #1 in xmcsubtpp',ie0sq
                stop
              endif
              if(ileg.eq.2.and.(ie0sq.eq.2.or.
     #            ie0sq.eq.3.or.ie0sq.eq.6))then
                write(*,*)'Error #2 in xmcsubtpp',ie0sq
                stop
              endif
              if(xz.ge.0.d0.and.xxkp.ge.0.d0.and. 
     #           xxkp.le.(2*e0sq))then
                flxsec(ileg,ie0sq)=.true.
                if(.not.flagmc)flagmc=.true.
              else
                flxsec(ileg,ie0sq)=.false.
              endif
            else
              if(ie0sq.eq.1.or.ie0sq.eq.3.or.
     #           ie0sq.eq.5.or.ie0sq.eq.6)then
                write(*,*)'Error #3 in xmcsubtpp',ie0sq
                stop
              endif
              if(xz.ge.0.d0.and.xxkp.ge.0.d0.and.
     #           xxkp.le.(2*e0sq+xm12))then
                flxsec(ileg,ie0sq)=.true.
                if(.not.flagmc)flagmc=.true.
              else
                flxsec(ileg,ie0sq)=.false.
              endif
            endif
            if(flxsec(ileg,ie0sq))then
              xmcz(ileg,ie0sq)=xz
              xmcxkp(ileg,ie0sq)=xxkp
              xmce0sq(ileg,ie0sq)=e0sq
              if(ileg.le.2)then
c Compute the rescaled relative pt for the branching
                pthw=(1-x)/2.d0*sqrt(s*(1-yi**2))
                ptresc(ileg,ie0sq)=(pthw-ptmin)/(ptmax-ptmin)
                if(ptresc(ileg,ie0sq).le.0.d0)then
                  emscwgt(ileg,ie0sq)=1.d0
                  emscav(ileg,ie0sq)=emsca_bare
                elseif(ptresc(ileg,ie0sq).lt.1.d0)then
                  emscwgt(ileg,ie0sq)=1-
     #              emscafun(ptresc(ileg,ie0sq),alpha)
                  emscav(ileg,ie0sq)=emsca_bare
                else
                  emscwgt(ileg,ie0sq)=0.d0
                  emscav(ileg,ie0sq)=ptmax
                endif
              else
                emscav(ileg,ie0sq)=sqrt(sh)
                emscwgt(ileg,ie0sq)=1.d0
              endif
            else
c Dead zone
              emscav(ileg,ie0sq)=sqrt(sh)
              emscwgt(ileg,ie0sq)=0.d0
            endif
          else
c The branching does not exist or it's outside kinematical accessible range
            flxsec(ileg,ie0sq)=.false.
            emscav(ileg,ie0sq)=sqrt(sh)
            emscwgt(ileg,ie0sq)=0.d0
          endif
        enddo
      enddo
c Even if flagmc=.false., evaluate the G functions, since they are used
c to compute the ME part of the MC subtraction term
      xalsfht=alsfht(jproc)
      xbesfht=besfht(jproc)
      xalclht=alclht(jproc)
      xbeclht=beclht(jproc)
      xalazi=alazi(jproc)
      xbeazi=beazi(jproc)
      gfactsf=gfunsoftht(x,ss,xm12,xm22,xalsfht,xbesfht)
      gfactcl=-1.d8
      gfactazi=gfunazi(yi,xalazi,xbeazi,deltai)
c Compute the cross sections if in the live zones. Also add the azimuthal
c correlation term (whose kernel is tmp1) when gfactazi#0; since this is
c of ME origin, the sum over ie0sq is not be performed. We (arbitrarily)
c associate it to the first non-zero ie0sq contribution
      if(flagmc)then
        wsum=0.d0
        do ileg=1,ilegmax
          call xinvtoinv(ileg,
     #      si,tki,uki,q1qi,q2qi,
     #      so,tko,uko,q1qo,q2qo,
     #      ss,tk,uk,q1q,q2q)
          sbar=2*x2to2(ileg,1)
          tbar=-2*x2to2(ileg,2)+xm12
          call fbornht(sbar,tbar,ithree,xborn)
          if( gfactazi.ne.0.d0 .and. ileg.le.2 .and.
     #        (jproc.eq.2.or.jproc.eq.3) )then
            icllazi=1-2*(ileg-1)
            call faziht(sbar,tbar,cth2,jproc,icllazi,xazi)
          endif
          do i=1,idrmax(jproc,jht)
            tmp1=0.d0
            do ie0sq=1,6
              if( flxsec(ileg,ie0sq) .and.
     #            ie0ht(jproc,i,ileg,ie0sq).eq.ie0sq )then
                lzone(i,ileg,ie0sq)=.true.
                z(i,ileg,ie0sq)=xmcz(ileg,ie0sq)
                e0sq=xmce0sq(ileg,ie0sq)
                tmp=0.d0
                if(1-x.lt.tiny.and.ileg.le.2)then
                  tmp=0.d0
                elseif(1-x.lt.tiny.and.ileg.eq.3.and.jproc.eq.3)then
                  tmp=0.d0
                elseif(1-yi.lt.tiny.and.ileg.eq.1)then
                  icode=iapcp(jproc,i)
                  tmp=8*ap_kern(x,icode)/ss
                  if( gfactazi.ne.0.d0.and.tmp1.eq.0.d0.and.
     #                (jproc.eq.2.or.jproc.eq.3) )
     #              tmp1=8*qin_kern(x,icode)/ss
                elseif(1+yi.lt.tiny.and.ileg.eq.2)then
                  icode=iapcm(jproc,i)
                  tmp=8*ap_kern(x,icode)/ss
                  if( gfactazi.ne.0.d0.and.tmp1.eq.0.d0.and.
     #                (jproc.eq.2.or.jproc.eq.3) )
     #              tmp1=8*qin_kern(x,icode)/ss
                else
                  xfact=4/s*(1-x)*(1-yi**2)
                  xz=xmcz(ileg,ie0sq)
                  xxkp=xmcxkp(ileg,ie0sq)
                  if(ileg.eq.1)then
                    icode=iapcp(jproc,i)
                  elseif(ileg.eq.2)then
                    icode=iapcm(jproc,i)
                  else
                    icode=4
                  endif
                  ap=ap_kern(xz,icode)/(1-xz)
                  de0sqdx=dx2to2dx(ileg,ie0sq)
                  de0sqdc=dx2to2dc(ileg,ie0sq)
                  xjac=xjac_zqttoxyht(ileg,e0sq,de0sqdx,de0sqdc,
     #                                xm12,xm22,ss,x,yi,cth1,cth2,
     #                                tk,uk,q1q,q2q)
                  tmp=xfact*xjac*ap/xxkp
                  if( gfactazi.ne.0.d0.and.tmp1.eq.0.d0.and.
     #                ileg.le.2.and.(jproc.eq.2.or.jproc.eq.3) )then
                    qk=qin_kern(xz,icode)/(1-xz)
                    tmp1=xfact*xjac*qk/xxkp
                  endif
                endif
                if(ileg.eq.1)then
                  iborn=idrlimcp(jht,jproc,i)
                elseif(ileg.eq.2)then
                  iborn=idrlimcm(jht,jproc,i)
                else
                  iborn=i
                endif
                xmatmc(i,ileg,ie0sq)=( tmp*xborn(iborn)+
     #                                 tmp1*xazi(i)*gfactazi )/
     #                               imcmult(jproc,i,ileg)
                xmatmc(i,ileg,ie0sq)=xmatmc(i,ileg,ie0sq)*gfactsf*
     #                               dcosdcth1(ileg,x,cth1)
              else
                xmatmc(i,ileg,ie0sq)=0.d0
                lzone(i,ileg,ie0sq)=.false.
                z(i,ileg,ie0sq)=-1.d0
              endif
              xmatmc(i,ileg,ie0sq)=xmatmc(i,ileg,ie0sq)*
     #                             emscwgt(ileg,ie0sq)
              wsum=wsum+abs(xmatmc(i,ileg,ie0sq))
            enddo
          enddo
        enddo
c Assign emsca on statistical basis
        if(wsum.gt.1.d-30)then
          rrnd=fk88random(ifk88seed)
          w1sum=0.d0
          iwsum=0
          ileg0=0
          ie0sq0=0
          do ileg=1,ilegmax
            do ie0sq=1,6
              if(flxsec(ileg,ie0sq).and.iwsum.eq.0)then
                do i=1,idrmax(jproc,jht)
                  w1sum=w1sum+abs(xmatmc(i,ileg,ie0sq))
                enddo
                if(w1sum.gt.rrnd*wsum)then
                  iwsum=1
                  ileg0=ileg
                  ie0sq0=ie0sq
                endif
              endif
            enddo
          enddo
          if(ileg0.eq.0.or.ie0sq0.eq.0)then
            write(*,*)'Fatal error #2 in xmcsubt:',ileg0,ie0sq0
            stop
          else
            emsca_prc(jproc)=emscav(ileg0,ie0sq0)
            emwgt_prc(jproc)=emscwgt(ileg0,ie0sq0)
          endif
        else
          emsca_prc(jproc)=sqrt(sh)
          emwgt_prc(jproc)=1.d0
        endif
      else
        emsca_prc(jproc)=sqrt(sh)
        emwgt_prc(jproc)=1.d0
      endif
      return
      end


      function dcosdcth1(ileg,x,cth1)
c Returns dcos(theta_in)/dcos(th1) or dcos(theta_out)/dcos(th1)
      implicit none
      real*8 dcosdcth1,x,cth1,tmp
      integer ileg
c
      tmp=1.d0
      if(ileg.eq.3)then
        tmp=4*x/(1+x-(1-x)*cth1)**2
      elseif(ileg.ne.1.and.ileg.ne.2)then
        write(*,*)'Error in dcosdcth1: unknown ileg',ileg
        stop
      endif
      dcosdcth1=tmp
      return
      end
    

      subroutine xinvtoinv(ileg,
     #  si,tki,uki,q1qi,q2qi,
     #  so,tko,uko,q1qo,q2qo,
     #  ss,tk,uk,q1q,q2q)
      implicit none
      real*8 si,tki,uki,q1qi,q2qi,so,tko,uko,q1qo,q2qo,
     #  ss,tk,uk,q1q,q2q
      integer ileg
c
      if(ileg.le.2)then
        ss=si
        tk=tki
        uk=uki
        q1q=q1qi
        q2q=q2qi
      elseif(ileg.eq.3)then
        ss=so
        tk=tko
        uk=uko
        q1q=q1qo
        q2q=q2qo
      else
        write(*,*)'Unknown leg in xinvtoinv:',ileg
        stop
      endif
      return
      end


      function gfunsoftht(x,s,xm12,xm22,alsf,besf)
c Gets smoothly to 0 in the soft limit. The same comments as for gfunsoft 
c apply, except for the dependence on the masses
      implicit none
      real*8 gfunsoftht,x,s,xm12,xm22,alsf,besf,xminsf,xgsoft,tt,tmp
      real*8 sh
      common/shadr/sh
      real*8 etacut
      common/cetacut/etacut
      integer isubttype
      common/cisubttype/isubttype
c
      if(besf.lt.0.d0)then
        xminsf=(sqrt(xm12)+sqrt(xm22))**2/sh
      else
        if(isubttype.eq.0)then
          xminsf=(sqrt(xm12)+sqrt(xm22))**2/s
        elseif(isubttype.eq.1)then
          xminsf=1-sqrt(etacut)
        else
          write(*,*)'Fatal error #1 in gfunsoftht',isubttype
          stop
        endif
      endif
      xgsoft=1.d0-(1-xminsf)*abs(besf)
      if(xgsoft.gt.0.99d0)xgsoft=0.99d0
      tt=(x-xgsoft)/(1.d0-xgsoft)
      if(tt.gt.1.d0)then
        write(6,*)'Fatal error #2 in gfunsoftht',x
        stop
      endif
      tmp=1.d0
      if(alsf.gt.0.d0)then
        if(tt.gt.0.d0.and.x.lt.0.99d0)
     #    tmp=(1-tt)**(2*alsf)/(tt**(2*alsf)+(1-tt)**(2*alsf))
        if(x.ge.0.99d0)tmp=0.d0
      endif
      gfunsoftht=tmp
      return
      end


      function gfuncoll(yy,alcl,becl,delta)
c Gets smoothly to 0 in the collinear limits; the function gfunsoft
c must be called before this function. The functional form is given
c in eq.(A.86) of FW, with alpha==alcl. tilde{x}_{DZ} is replaced here
c by ygcoll, and x_{DZ} by ymincl. The function is different from 1
c in the region ygcoll<|y|<1. Call with
c  becl<0  ==> ymincl=0
c  becl>0  ==> ymincl=Max(0,1-delta) for standard subtraction
c              ymincl=0 for zeta-subtraction
c  |becl|-->0 ==> ygcoll-->1
c  |becl|-->1 ==> ygcoll-->ymincl
c This function has been derived from the analogous function in the
c QQ code; alcl, becl and delta are now given in input rather than in common;
c the dependence on gacl has been eliminated (was only useful for testing
c purposes), and as a consequence the entry xx has been removed
      implicit none
      real * 8 gfuncoll,yy,alcl,becl,delta,y,ymincl,ygcoll,tt,tmp
      real * 8 etacut
      common/cetacut/etacut
      integer isubttype
      common/cisubttype/isubttype
c
      y=yy
      if(becl.lt.0.d0)then
        ymincl=0.d0
      else
        if(isubttype.eq.0)then
          ymincl=max(0.d0,1.d0-delta)
        elseif(isubttype.eq.1)then
          ymincl=0.d0
        else
          write(*,*)'Fatal error #1 in gfuncoll',isubttype
          stop
        endif
      endif
      ygcoll=1.d0-(1-ymincl)*abs(becl)
      if(ygcoll.gt.0.99d0)ygcoll=0.99d0
      tt=(abs(y)-ygcoll)/(1.d0-ygcoll)
      if(tt.gt.1.d0)then
        write(6,*)'Fatal error #2 in gfuncoll',tt
        stop
      endif
      tmp=1.d0
      if(alcl.gt.0.d0)then
        if(tt.gt.0.d0.and.abs(y).lt.0.99d0)
     #    tmp=(1-tt)**(2*alcl)/(tt**(2*alcl)+(1-tt)**(2*alcl))
        if(abs(y).ge.0.99d0)tmp=0.d0
      endif
      gfuncoll=tmp
      return
      end


      function gfunazi(y,alazi,beazi,delta)
c This function multiplies the azimuthal correlation term in the MC 
c subtraction kernel; it is not the same as in QQ code. We have
c   alazi<0  ==>  gfunazi=1-gfuncoll(|alazi|)
c   alazi>0  ==>  gfunazi=0
c ie in the testing phase (alazi<0) we include an azimuthal-dependent
c contribution in the MC subtraction terms
      implicit none
      real*8 gfunazi,y,alazi,beazi,delta,aalazi,tmp,gfuncoll
c
      tmp=0.d0
      if(alazi.lt.0.d0)then
        aalazi=abs(alazi)
        tmp=1.d0-gfuncoll(y,aalazi,beazi,delta)
      endif
      gfunazi=tmp
      return
      end


      function emscafun(x,alpha)
      implicit none
      real*8 emscafun,x,alpha
c
      if(x.lt.0.d0.or.x.gt.1.d0)then
        write(6,*)'Fatal error in emscafun'
        stop
      endif
      emscafun=x**(2*alpha)/(x**(2*alpha)+(1-x)**(2*alpha))
      return
      end


      function emscainv(r,alpha)
c Inverse of emscafun; implemented only for alpha=1 for the moment
      implicit none
      real*8 emscainv,r,alpha
c
      if(r.lt.0.d0.or.r.gt.1.d0.or.alpha.ne.1.d0)then
        write(6,*)'Fatal error in emscafun'
        stop
      endif
      emscainv=(r-sqrt(r-r**2))/(2*r-1)
      return
      end


      function ap_kern_qc(x,q2til,xm2,index)
      implicit none 
      real*8 ap_kern_qc,x,q2til,xm2
      integer index
      real*8 vcf,vtf,vca,tmp,ap_kern
      parameter (vcf=4.d0/3.d0)
      parameter (vtf=1.d0/2.d0)
      parameter (vca=3.d0)
c
      if(index.eq.1)then
        tmp=0.d0
      elseif(index.eq.2)then
        tmp=vtf*2*xm2/(x*q2til)
      elseif(index.eq.3)then
        tmp=-vcf*2*xm2/(x*q2til)
      elseif(index.eq.4)then
        tmp=-vcf*2*xm2/(x*q2til)
      else
        write(6,*)'Error in ap_kern_qc: wrong index value',index
        stop
      endif
      ap_kern_qc=ap_kern(x,index)+tmp
      if(ap_kern_qc.lt.0.d0)then
        write(*,*)'Error in ap_kern_qc',ap_kern_qc,x,q2til,xm2
        stop
      endif
      return
      end
c
c
c End of MC subtraction terms
c
c
c
c
c Begin of utility routines for q2til, z, and 2-->2 invariants
c
c
c The following routines are relevant to HERWIG shower variables, as
c computed for single-top production, Ht channel. Legs 1,2 are identical
c to those for s- and t-channel, and these routines are just wrappers;
c leg 3 has changed since it is treated here with an initial-state
c parametrization. As for the s- and t-channel, the invariants
c given in input to the routines are those computed by invar_in,
c and use FNR conventions (i.e., are defined as (p-k)^2). 
c Those used within these routines follow the notes, and therefore use
c MNR conventions (i.e., are defined as -2p.k). Using eq.(2.7) of FNR
c and the table of the draft we obtain
c
c  MNR   FNR
c  q1c = m12-s-tk-q1q
c  q2c = m22-s-uk-q2q
c  w1  = -q1q+q2q-tk
c  w2  = q1q-q2q-uk
c
      function zherppht(ileg,e0sq,xm12,xm22,s,x,yi,cth1,cth2,
     #                 xtk,xuk,xq1q,xq2q)
c Returns Herwig shower variable z; xm12 and xm22 are the top and W 
c masses squared
      implicit none
      real * 8 zherppht,e0sq,xm12,xm22,s,x,yi,cth1,cth2,
     # xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,dummy,zherpp,w1,w2,xm12red,xm22red,beta,sth1,
     # cthg,beta2,eps2,zeta1
      parameter (tiny=1.d-4)
      parameter (dummy=-1.d3)
c
c incoming partons #1 (left) and #2 (right)
      if(ileg.eq.1.or.ileg.eq.2)then
        zherppht=zherpp(ileg,e0sq,xm12,s,x,yi,dummy,xtk,xuk,xq1q,xq2q)
c outgoing parton #3 (top)
      elseif(ileg.eq.3)then
        w1=-xq1q+xq2q-xtk
        w2=-xq2q+xq1q-xuk
        xm12red=xm12/s
        xm22red=xm22/s
        beta=sqrt(1-2*(xm12+xm22)/s+(xm12-xm22)**2/s**2)
        if(1-x.lt.tiny)then
          sth1=sqrt(1-cth1**2)
          cthg=yi*cth1+sqrt(1-yi**2)*cth2*sth1
          zherppht=1-(1+cthg)*(1-x)/(1+xm12red-xm22red+beta)
        else
          beta2=sqrt((1-(xm12-xm22)/(s-w1))**2-4*s*xm22/(s-w1)**2)
          eps2=1-(xm12-xm22)/(s-w1)
          zeta1=( (2*s-(s-w1)*eps2)*w2+
     #            (s-w1)*((w1+w2)*beta2-eps2*w1) )/
     #          ( (s-w1)*beta2*(2*s-(s-w1)*eps2+(s-w1)*beta2) )
          zherppht=1-zeta1
        endif
      else
        write(6,*)'zherppht: unknown parton number'
        stop
      endif
      return
      end


      function xqherppht(ileg,e0sq,xm12,xm22,s,x,yi,cth1,cth2,
     #                  xtk,xuk,xq1q,xq2q)
c Returns Herwig shower variable ktilde; xm12 and xm22 are the top and W 
c masses squared
      implicit none
      real * 8 xqherppht,e0sq,xm12,xm22,s,x,yi,cth1,cth2,
     # xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,vtiny,dummy,xqherpp,w1,w2,xm12red,xm22red,beta,sth1,
     # cthg,z,zherppht
      parameter (tiny=1.d-3)
      parameter (vtiny=1.d-5)
      parameter (dummy=-1.d3)
c
c incoming partons #1 (left) and #2 (right)
      if(ileg.eq.1.or.ileg.eq.2)then
        xqherppht=xqherpp(ileg,e0sq,xm12,s,x,yi,dummy,xtk,xuk,xq1q,xq2q)
c outgoing parton #3 (top)
      elseif(ileg.eq.3)then
        w1=-xq1q+xq2q-xtk
        w2=-xq2q+xq1q-xuk
        xm12red=xm12/s
        xm22red=xm22/s
        if(1-x.lt.tiny)then
          beta=sqrt(1-2*(xm12+xm22)/s+(xm12-xm22)**2/s**2)
          sth1=sqrt(1-cth1**2)
          cthg=yi*cth1+sqrt(1-yi**2)*cth2*sth1
          xqherppht=(1+xm12red-xm22red+beta)*
     #              (1+xm12red-xm22red-beta*cthg)*s/
     #              (2.d0*(1+cthg))
        else
          z=zherppht(ileg,e0sq,xm12,xm22,s,x,yi,cth1,cth2,
     #               xtk,xuk,xq1q,xq2q)
          if(z.lt.1-vtiny.and.z.gt.vtiny)then
            xqherppht=w1/(z*(1-z))
          else
            xqherppht=-1.d0
          endif
        endif
      else
        write(6,*)'xqherppht: unknown parton number'
        stop
      endif
      return
      end


      function xjac_zqttoxyht(ileg,e0sq,de0sqdx,de0sqdc,xm12,xm22,s,
     #                        x,yi,cth1,cth2,xtk,xuk,xq1q,xq2q)
c Returns the jacobian d(z,ktil)/d(x,c), where z and xi are Herwig shower 
c variables, and x and c are FKS variables
      implicit none 
      real * 8 xjac_zqttoxyht,e0sq,de0sqdx,de0sqdc,xm12,xm22,s,
     # x,yi,cth1,cth2,xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,dummy,tmp,xjac_zqttoxy,w1,w2,xm12red,xm22red,sth1,
     # beta,cthg,zmo,dktil0sfdc,z,zherppht,xkp,xqherppht,cpsip,spsip,
     # beta2,betax,eps2,zeta1,dw1dx,dw2dx,dw1dc,dw2dc,dzeta1dw2
      parameter (tiny=1.d-4)
      parameter (dummy=-1.d3)
c
      tmp=0.d0
c incoming partons #1 (left) and #2 (right)
      if(ileg.eq.1.or.ileg.eq.2)then
        tmp=xjac_zqttoxy(ileg,e0sq,de0sqdx,de0sqdc,xm12,s,
     #                   x,yi,dummy,xtk,xuk,xq1q,xq2q)
c outgoing parton #3 (top)
      elseif(ileg.eq.3)then
        w1=-xq1q+xq2q-xtk
        w2=-xq2q+xq1q-xuk
        xm12red=xm12/s
        xm22red=xm22/s
        sth1=sqrt(1-cth1**2)
        beta=sqrt(1-2*(xm12+xm22)/s+(xm12-xm22)**2/s**2)
        if(1-x.lt.tiny)then
          cthg=yi*cth1+sqrt(1-yi**2)*cth2*sth1
          zmo=-(1+cthg)/(1+xm12red-xm22red+beta)
          dktil0sfdc=-beta*(1+xm12red-xm22red+beta)*s/(2.d0*(1+cthg))-
     #               (1+xm12red-xm22red+beta)*
     #               (1+xm12red-xm22red-beta*cthg)*s/(2.d0*(1+cthg)**2)
          tmp=-zmo*dktil0sfdc
        else
          cpsip=(1+yi-x*(1-yi))/(1+yi+x*(1-yi))
          spsip=sqrt(4*x*(1-yi**2))/(1+yi+x*(1-yi))
          cthg=cpsip*cth1+spsip*cth2*sth1
          betax=sqrt(1-2*(xm12+xm22)/(x*s)+(xm12-xm22)**2/(x*s)**2)
          beta2=sqrt((1-(xm12-xm22)/(s-w1))**2-4*s*xm22/(s-w1)**2)
          eps2=1-(xm12-xm22)/(s-w1)
          dw1dx=-s/2.d0*(1-betax*cthg)-
     #         s/2.d0*(1-x)*cthg*(s*x*(xm12+xm22)-(xm12-xm22)**2)/
     #         (betax*s**2*x**3)-(xm12-xm22)/(2*x**2)
          dw2dx=-s/2.d0*(1+betax*cthg)+
     #         s/2.d0*(1-x)*cthg*(s*x*(xm12+xm22)-(xm12-xm22)**2)/
     #         (betax*s**2*x**3)+(xm12-xm22)/(2*x**2)
          dw1dc=-s/2.d0*betax*(1-x)
          dw2dc=s/2.d0*betax*(1-x)
          zeta1=( (2*s-(s-w1)*eps2)*w2+
     #            (s-w1)*((w1+w2)*beta2-eps2*w1) )/
     #          ( (s-w1)*beta2*(2*s-(s-w1)*eps2+(s-w1)*beta2) )
          dzeta1dw2=1/((s-w1)*beta2)
          tmp=(dw1dx*dw2dc-dw1dc*dw2dx)*dzeta1dw2/(zeta1*(1-zeta1))
          tmp=tmp*4*x/(1+yi+x*(1-yi))**2
        endif
      else
        write(6,*)'xjac_zqttoxyht: unknown parton number'
        stop
      endif
      xjac_zqttoxyht=abs(tmp)
      return 
      end 


      subroutine gete0sq(ich,ileg,z1,z2,xm12,xm22,s,x,yi,cth1,cth2,
     #                   yj,phii,phij,xtk,xuk,xq1q,xq2q,
     #                   e0sq,de0sqdx,de0sqdc)
c Returns Herwig shower scale squared, and its derivatives, defined as
c follows for s- (ich=1) and t-channel (ich=2): de0sqdx=d(e0sq)/dx, 
c de0sqdc=F*d(e0sq)/dc, with F=(1-yi^2) for legs 1 and 2, F=1 for leg 3, 
c and F=(1-yj) for leg 4. In the case of initial-state emissions, we have 
c x=xii and c=yi; in the case of final-state emissions, we have x=xii and 
c c=yj. The scales and their derivatives are stored in arrays, with the 
c following conventions
c   e0sq(1) --> E0^2=|p1.p2|
c   e0sq(2) --> E0^2=|p1.k1|
c   e0sq(3) --> E0^2=|p1.k2|
c   e0sq(4) --> E0^2=|p2.k1|
c   e0sq(5) --> E0^2=|p2.k2|
c   e0sq(6) --> E0^2=|k1.k2|
c and analogously for de0sqdx and for de0sqdc.
c In the case of Ht channel (ich=3) the emission from leg 3 is treated
c similarly to what done in the case of QQbar production, since we use 
c an initial-state parametrization. In order to treat Ht channel, xm22
c is now an entry of this routine (rather than a parameter)
      implicit none
      real * 8 z1,z2,xm12,xm22,s,x,yi,cth1,cth2,yj,phii,phij,xtk,xuk,
     # xq1q,xq2q,e0sq(6),de0sqdx(6),de0sqdc(6)
      integer ich,ileg
      real * 8 dm12,sm12,sth1,betax,q1q,q2q,q1c,q2c,w1,w2,uk,tk,
     # y,sbar,dsbardx,dsbardc,xmn,xpl,galonred,betalon,cpsi,spsi,dtkdx,
     # dukdx,dq1qdx,dq2qdx,dxpldx,dxmndx,dtkdy,dukdy,dq1qdy,dq2qdy,
     # dxpldy,dxmndy,tbar,dtbardx,dtbardc,beta,beta1,beta2,si,sj,cphij,
     # sphij,dq1cdx,dw2dx,dq1cdy,dw2dy,ubar,dubardx,dubardc,tiny,
     # p2k1,p2k2,k1k2,cthg,sphigsthg,cpsip,spsip,ctho,dydcpsip,
     # dydx,dcpsipdx,dcpsipdcthg,dcth1dx,dw1dx,dtkdc,dw1dc,dq1qdc
      integer ia1ora2
      common/cia1ora2/ia1ora2
      parameter (tiny=1.d-4)
c
      if(ia1ora2.ne.1)then
        write(*,*)'gete0sq: unknown option',ia1ora2
        stop
      endif
      if(ich.lt.1.or.ich.gt.3)then
        write(*,*)'gete0sq: unknown channel',ich
        stop
      endif
      dm12=xm12-xm22
      sm12=xm12+xm22
      sth1=sqrt(1-cth1**2)
      beta=sqrt(1-2*(sm12)/s+dm12**2/s**2)
      betax=sqrt(1-2*(sm12)/(x*s)+dm12**2/(x*s)**2)
      q1q = xq1q-xm12
      q2q = xq2q-xm22
      q1c = xm12-s-xtk-xq1q
      q2c = xm22-s-xuk-xq2q
      w1  = -xq1q+xq2q-xtk
      w2  = xq1q-xq2q-xuk
      uk  = xuk
      tk  = xtk
      if(ileg.eq.1.or.ileg.eq.2)then
        y=yi
        sbar=x*s
        dsbardx=s
        dsbardc=0.d0
        if(1-x.lt.tiny)then
          tbar=-s/2.d0*(1+dm12/s-beta*cth1)
          dtbardx=-s/2.d0*( 1-(s-sm12)/(s*beta)*cth1-
     #            z1/(z1+z2)*beta*sqrt(1-y**2)*cth2*sth1 )
          dtbardc=0.d0
        elseif(1-y.lt.tiny)then
          tbar=-s*x/2.d0*(1+dm12/(s*x)-betax*cth1)
          dtbardx=-s/2.d0*( 1-(s*x-sm12)/(s*x*betax)*cth1 )
          dtbardc=0.d0
        elseif(1+y.lt.tiny)then
          tbar=-s*x/2.d0*(1+dm12/(s*x)-betax*cth1)
          dtbardx=-s/2.d0*( 1-(s*x-sm12)/(s*x*betax)*cth1 )
          dtbardc=0.d0
        else
          xmn=((s+uk)*z1/s-(s+tk)*z2/s)/2.d0
          xpl=((s+uk)*z1/s+(s+tk)*z2/s)/2.d0
          galonred=sqrt(xpl**2-z1*z2*tk*uk/s**2)
          betalon=-xmn/galonred
          cpsi=1-8*x/((1+y+x*(1-y))*(1-y+x*(1+y)))
          spsi=4*(1-x)*sqrt(x*(1-y**2))/
     #         ((1+y+x*(1-y))*(1-y+x*(1+y)))
          dtkdx=s*(1-y)/2.d0
          dukdx=s*(1+y)/2.d0
          dq1qdx=1/(4*betax*x**2)*(1+y+x*(1-y))*
     #      ( dm12*betax+(sm12-dm12**2/(s*x))*cth1 ) -
     #      s/4.d0*(1-y)*(1+dm12/(x*s)-betax*cth1)
          dq2qdx=-s/4.d0*(1+y)*( 1-dm12/(x*s)+
     #      betax*(cpsi*cth1+spsi*cth2*sth1) ) -
     #      (1-y+x*(1+y))*(sm12-dm12**2/(s*x))/
     #        (4*x**2*betax)*(cth2*sth1*spsi+cth1*cpsi) -
     #      (1-y+x*(1+y))*dm12/(4*x**2)+
     #      betax*s*( 2*(1-x**2)*(1-y**2)/
     #      ((1+y+x*(1-y))**2*(1-y+x*(1+y)))*cth1 -
     #      (1+x)*(1-y**2)*((1-y**2)*(1+x**2)-2*x*(3-y**2))/
     #      (2*(1+y+x*(1-y))**2*(1-y+x*(1+y))*sqrt(x*(1-y**2)))*
     #      cth2*sth1 )
          dxpldx=(z1*(1+y)+z2*(1-y))/4.d0
          dxmndx=(z1*(1+y)-z2*(1-y))/4.d0
          dtkdy=s*(1-x)/2.d0
          dukdy=-s*(1-x)/2.d0
          dq1qdy=-s/4.d0*(1-x)*(1+dm12/(x*s)-betax*cth1)
          dq2qdy=s/4.d0*(1-x)*( 1-dm12/(x*s)+
     #      betax*(cpsi*cth1+spsi*cth2*sth1) ) -
     #      betax*s*( -4*(1-x)**2*x*y/
     #                ((1+y+x*(1-y))**2*(1-y+x*(1+y)))*cth1 +
     #      (1-x)*x*y*((1-y**2)*(1+x**2)-2*x*(3-y**2))/
     #      ( (1+y+x*(1-y))**2*(1-y+x*(1+y))*sqrt(x*(1-y**2)) )*
     #      cth2*sth1 )
          dxpldy=-(1-x)*(z1-z2)/4.d0
          dxmndy=-(1-x)*(z1+z2)/4.d0
          tbar=-sbar/2.d0*( 1-(z2*(q1q-q1c)+z1*(q2q-q2c))/
     #                        (2*s*galonred) )
     #         -dm12/2.d0*(1-betalon) 
          dtbardx=-(s+tk+uk)/2.d0*(
     #      ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )*
     #      ( 2*dxpldx*xpl-dukdx*tk*z1*z2/s**2-dtkdx*uk*z1*z2/s**2 )/
     #      ( 4*s*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) -
     #      ( (2*dq2qdx+dukdx)*z1+(2*dq1qdx+dtkdx)*z2 )/
     #      ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) ) -
     #      (dtkdx+dukdx)/2.d0*(1 -
     #          ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )/
     #          ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) )-
     #      dm12/2.d0*( dxmndx/sqrt(xpl**2-z1*z2*tk*uk/s**2)-xmn*
     #      ( 2*dxpldx*xpl-dukdx*tk*z1*z2/s**2-dtkdx*uk*z1*z2/s**2 )/
     #      ( 2*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) )
          dtbardc=-(s+tk+uk)/2.d0*(
     #      ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )*
     #      ( 2*dxpldy*xpl-dukdy*tk*z1*z2/s**2-dtkdy*uk*z1*z2/s**2 )/
     #      ( 4*s*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) -
     #      ( (2*dq2qdy+dukdy)*z1+(2*dq1qdy+dtkdy)*z2 )/
     #      ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) ) -
     #      (dtkdy+dukdy)/2.d0*(1 -
     #          ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )/
     #          ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) )-
     #      dm12/2.d0*( dxmndy/sqrt(xpl**2-z1*z2*tk*uk/s**2)-xmn*
     #      ( 2*dxpldy*xpl-dukdy*tk*z1*z2/s**2-dtkdy*uk*z1*z2/s**2 )/
     #      ( 2*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) )
          dtbardc=(1-y**2)*dtbardc
        endif
      elseif(ileg.eq.3)then
        sbar=s
        dsbardx=0.d0
        dsbardc=0.d0
        if(ich.eq.1.or.ich.eq.2)then
          if(1-x.lt.tiny)then
            tbar=-s/2.d0*(1+xm12/s+beta*yi)
          else
            beta2=sqrt((1-dm12/(s-w1))**2-4*s*xm22/(s-w1)**2)
            tbar=-s/2.d0*(1-(q2q-q1c)/(s-w1)*beta/beta2)-dm12/2.d0
          endif
          dtbardx=0.d0
          dtbardc=0.d0
        else
          if(1-x.lt.tiny)then
            tbar=-s/2.d0*(1+dm12/s-beta*cth1)
            cthg=yi*cth1+sqrt(1-yi**2)*cth2*sth1
            sphigsthg=sqrt(1-yi**2)*cth1*cth2-yi*sth1
            dtbardx=s/4.d0*(cthg*cth1-yi+beta*(cthg-yi*cth1-sth1**2))-
     #              dm12/4.d0*sth1*sphigsthg
            dtbardc=0.d0
          else
            beta2=sqrt((1-dm12/(s-w1))**2-4*s*xm22/(s-w1)**2)
            tbar=-s/2.d0*(1-(q2q-q1c)/(s-w1)*beta/beta2)-dm12/2.d0
            cpsip=(1+yi-x*(1-yi))/(1+yi+x*(1-yi))
            spsip=sqrt(4*x*(1-yi**2))/(1+yi+x*(1-yi))
            cthg=cpsip*cth1+spsip*cth2*sth1
            ctho=-(1-x-(1+x)*cth1)/(1+x-(1-x)*cth1)
            dydcpsip=4*x/(1+x-cpsip*(1-x))**2
            dydx=2*(1-cpsip**2)/(1+x-cpsip*(1-x))**2
            sphigsthg=cth1*cth2*spsip-cpsip*sth1
            dcpsipdx=-sth1*(cthg*sth1+sphigsthg*cth1)/(2*x)
            dcpsipdcthg=cth1+sth1*sphigsthg*cthg/(1-cthg**2)
            dcth1dx=-2*(1-ctho**2)/(1+x+(1-x)*ctho)**2
            dtkdx=s*(1-yi)/2.d0+
     #            s*(1-x)/2.d0*(dydx+dydcpsip*dcpsipdx)
            dw1dx=-s/2.d0*(1-betax*cthg)-
     #             s/2.d0*(1-x)*cthg*(s*x*(xm12+xm22)-(xm12-xm22)**2)/
     #             (betax*s**2*x**3)-(xm12-xm22)/(2*x**2)
            dq1qdx=-dtkdx/2.d0*(1+(xm12-xm22)/(x*s)-betax*cth1)+
     #        s/4.d0*(1+yi+x*(1-yi))*( (s*x*(xm12+xm22)-
     #          (xm12-xm22)**2)/(betax*s**2*x**3)*cth1+
     #          betax*dcth1dx+(xm12-xm22)/(s*x**2) )
            dtkdc=s*(1-x)/2.d0*dydcpsip*dcpsipdcthg
            dw1dc=-s/2.d0*(1-x)*betax
            dq1qdc=-dtkdc/2.d0*(1+(xm12-xm22)/(x*s)-betax*cth1)
            dtbardx=s/2.d0*( 
     #        beta*(2*dq1qdx+2*dtkdx+dw1dx)/(beta2*(s-w1)) +
     #        beta*dw1dx*(4*xm22*s+(xm12-xm22)*(s-w1-xm12+xm22))*
     #           (2*q1q+s+2*tk+w1+xm12-xm22)/(beta2**3*(s-w1)**4) +
     #        beta*dw1dx*(2*q1q+s+2*tk+w1+xm12-xm22)/
     #                   (beta2*(s-w1)**2) )
            dtbardc=s/2.d0*( 
     #        beta*(2*dq1qdc+2*dtkdc+dw1dc)/(beta2*(s-w1)) +
     #        beta*dw1dc*(4*xm22*s+(xm12-xm22)*(s-w1-xm12+xm22))*
     #           (2*q1q+s+2*tk+w1+xm12-xm22)/(beta2**3*(s-w1)**4) +
     #        beta*dw1dc*(2*q1q+s+2*tk+w1+xm12-xm22)/
     #                   (beta2*(s-w1)**2) )
          endif
        endif
      elseif(ileg.eq.4)then
        if(ich.eq.3)then
          write(*,*)'gete0sq: Ht channel with leg 4'
          stop
        endif
        sbar=s
        dsbardx=0.d0
        dsbardc=0.d0
        if(1-x.lt.tiny)then
          tbar=-s/2.d0*(1+xm12/s+beta*yi)
          dtbardx=-s/2.d0*beta*sqrt(1-yi**2)*sqrt(1-yj**2)*
     #            cos(phij)/(1-xm12/s)
          dtbardc=0.d0
        elseif(1-yj.lt.tiny)then
          tbar=-s/2.d0*(1+xm12/s+beta*yi)
          dtbardx=0.d0
          dtbardc=0.d0
        elseif(1+yj.lt.tiny)then
c Insert this special case only to prevent a division by zero below (1/sj).
c The dtbardc is proportional to 1/sqrt(1+yj), but it is damped in the
c jacobian. Therefore, set it to zero here
          if(s*x**2-xm12.gt.0.d0)then
            tbar=-s/2.d0*(1+xm12/s+beta*yi)
          else
            tbar=-s/2.d0*(1+xm12/s-beta*yi)
          endif          
          dtbardx=0.d0
          dtbardc=0.d0
        else
          beta1=sqrt((1+dm12/(s-w2))**2-4*s*xm12/(s-w2)**2)
          si=sqrt(1-yi**2)
          sj=sqrt(1-yj**2)
          cphij=cos(phij)
          sphij=sin(phij)
          dtkdx=s/2.d0*(1-yi*yj+si*sj*cphij)
          dukdx=s/2.d0*(1+yi*yj-si*sj*cphij)
          dq1cdx=-(1-yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dq2qdx=-(1+yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dw2dx=(1-yj)*(s*(1+yj-x*(2*(1+yj)+x*(1-yj)))+2*xm12)/
     #          (1+yj+x*(1-yj))**2
          dtkdy=s/2.d0*(1-x)*(yi+yj*si*cphij/sj)
          dukdy=-s/2.d0*(1-x)*(yi+yj*si*cphij/sj)
          dq1cdy=(1-x)*(1-yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dq2qdy=(1-x)*(1+yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dw2dy=-2*(1-x)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          tbar=-s/2.d0*(1-(q1q-q2c)/(s-w2)*beta/beta1)-dm12/2.d0
          dtbardx=s*beta/(2*beta1*(s-w2))*( dukdx-dtkdx+dq2qdx-dq1cdx +
     #            dw2dx*(q1q-q2c)*(s-w2+xm12)/
     #            ((s-w2)**2-2*xm12*(s+w2)+xm12**2) )
          dtbardc=s*beta/(2*beta1*(s-w2))*( dukdy-dtkdy+dq2qdy-dq1cdy +
     #            dw2dy*(q1q-q2c)*(s-w2+xm12)/
     #            ((s-w2)**2-2*xm12*(s+w2)+xm12**2) )
          dtbardc=(1-yj)*dtbardc
        endif
      else
        write(*,*)'gete0sq: unknown parton number',ileg
        stop
      endif
c
      ubar=-sbar-tbar
      dubardx=-dsbardx-dtbardx
      dubardc=-dsbardc-dtbardc
c E0^2=|p1.p2|
      e0sq(1)=abs(sbar/2.d0)
      de0sqdx(1)=sign(1.d0,sbar)*dsbardx/2.d0
      de0sqdc(1)=sign(1.d0,sbar)*dsbardc/2.d0
c E0^2=|p1.k1|
      e0sq(2)=abs(tbar/2.d0)
      de0sqdx(2)=sign(1.d0,tbar)*dtbardx/2.d0
      de0sqdc(2)=sign(1.d0,tbar)*dtbardc/2.d0
c E0^2=|p1.k2|
      e0sq(3)=abs(ubar/2.d0)
      de0sqdx(3)=sign(1.d0,ubar)*dubardx/2.d0
      de0sqdc(3)=sign(1.d0,ubar)*dubardc/2.d0
c E0^2=|p2.k1|
      p2k1=-(ubar-xm12+xm22)/2.d0
      e0sq(4)=abs(p2k1)
      de0sqdx(4)=sign(1.d0,-p2k1)*dubardx/2.d0
      de0sqdc(4)=sign(1.d0,-p2k1)*dubardc/2.d0
c E0^2=|p2.k2|
      p2k2=-(tbar+xm12-xm22)/2.d0
      e0sq(5)=abs(p2k2)
      de0sqdx(5)=sign(1.d0,-p2k2)*dtbardx/2.d0
      de0sqdc(5)=sign(1.d0,-p2k2)*dtbardc/2.d0
c E0^2=|k1.k2|
      k1k2=(sbar-xm12-xm22)/2.d0
      e0sq(6)=abs(k1k2)
      de0sqdx(6)=sign(1.d0,k1k2)*dsbardx/2.d0
      de0sqdc(6)=sign(1.d0,k1k2)*dsbardc/2.d0
c
      return
      end
c
c
c Here include z,q2til,xjac taken from the s- and t-channel
c
c
      function zherpp(ileg,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      implicit none
      real * 8 zherpp,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,xw1,xw2,xmq2,xm2red,zeta1,eps1,beta1,zeta2num,
     # zeta2den,zeta2
      parameter (tiny=1.d-5)
c
c incoming parton (left)
      if(ileg.eq.1)then
        zherpp=(1-yi+x*(1+yi))/2.d0
c incoming parton (right)
      elseif(ileg.eq.2)then
        zherpp=(1+yi+x*(1-yi))/2.d0
c outgoing top
      elseif(ileg.eq.3)then
        xw1=-xq1q+xq2q-xtk
        xw2=xq1q-xq2q-xuk
        xmq2=xm12
        if(1-x.lt.tiny)then
          zherpp=1-(1-x)*(1-yj)/2.d0
        else
          zeta1=xw2/(s-xw1-xmq2)
          zherpp=1-zeta1
        endif
c outgoing light antiquark
      elseif(ileg.eq.4)then
        xw1=-xq1q+xq2q-xtk
        xw2=xq1q-xq2q-xuk
        xmq2=xm12
        xm2red=xmq2/s
        if(1-x.lt.tiny)then
          zherpp=1-(1-x)*(1+yj)/(2.d0*(1-xm2red))
        elseif(1-yj.lt.tiny)then
          zherpp=(x-xm2red)/(1-xm2red)+
     #           (x-xm2red)*(x-xm2red*(2-x))/
     #           (2*(1-xm2red)**3)*(1-x)*(1-yj)
        else
          eps1=1+xmq2/(s-xw2)
          beta1=sqrt(eps1**2-4.d0*s*xmq2/(s-xw2)**2)
          zeta2num=
     #       (2*s-(s-xw2)*eps1)*xw1+(s-xw2)*((xw1+xw2)*beta1-eps1*xw2)
          zeta2den=(s-xw2)*beta1*(2*s-(s-xw2)*(eps1-beta1))
          zeta2=zeta2num/zeta2den
          zherpp=1-zeta2
        endif
      else
        write(*,*)'zherpp: unknown parton number'
        stop
      endif
      if(zherpp.lt.0.d0.or.zherpp.gt.1.d0)then
        write(6,*)'zherpp: fatal error',zherpp
        stop
      endif
      return
      end


      function xqherpp(ileg,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
      implicit none
      real * 8 xqherpp,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,vtiny,xw1,xw2,xmq2,xm2red,z,zherpp
      parameter (tiny=1.d-5)
      parameter (vtiny=1.d-8)
c
c incoming parton (left)
      if(ileg.eq.1)then
        if(yi.gt.-1+tiny)then
          xqherpp=s*(1-yi)/(1+yi)
        else
          xqherpp=-1.d0
        endif
c incoming parton (right)
      elseif(ileg.eq.2)then
        if(yi.lt.1-tiny)then
          xqherpp=s*(1+yi)/(1-yi)
        else
          xqherpp=-1.d0
        endif
c outgoing top
      elseif(ileg.eq.3)then
        xw1=-xq1q+xq2q-xtk
        xmq2=xm12
        xm2red=xmq2/s
        if(1-x.lt.tiny)then
          xqherpp=s*(1+yj+xm2red*(1-yj))/(1-yj)
        else
          z=zherpp(ileg,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
          if(z.lt.1-vtiny.and.z.gt.vtiny)then
            xqherpp=xw1/(z*(1-z))
          else
            xqherpp=-1.d0
          endif
        endif
c outgoing light antiquark
      elseif(ileg.eq.4)then
        xw2=xq1q-xq2q-xuk
        xmq2=xm12
        xm2red=xmq2/s
        if(1-x.lt.tiny)then
          xqherpp=s*(1-xm2red)**2*(1-yj)/(1+yj)
        elseif(1-yj.lt.tiny)then
          xqherpp=s*(1-xm2red)**2*(1-yj)/2.d0
        else
          z=zherpp(ileg,e0sq,xm12,s,x,yi,yj,xtk,xuk,xq1q,xq2q)
          if(z.lt.1-vtiny.and.z.gt.vtiny)then
            xqherpp=xw2/(z*(1-z))
          else
            xqherpp=-1.d0
          endif
        endif
      else
        write(*,*)'xqherpp: unknown parton number'
        stop
      endif
      return
      end


      function xjac_zqttoxy(ileg,e0sq,de0sqdx,de0sqdc,xm12,s,x,yi,yj,
     #                      xtk,xuk,xq1q,xq2q)
      implicit none
      real * 8 xjac_zqttoxy,e0sq,de0sqdx,de0sqdc,xm12,s,x,yi,yj,
     # xtk,xuk,xq1q,xq2q
      integer ileg
      real * 8 tiny,tmp,xw1,xw2,xmq2,xm2red,zeta1,dzeta1dw2,dq1cdx,
     # dq2qdx,dw2dx,dw1dx,dq1cdy,dq2qdy,dw2dy,dw1dy,eps1,beta1,
     # zeta2num,zeta2den,zeta2,dzeta2dw1
      parameter (tiny=1.d-5)
c
c incoming parton (left)
      if(ileg.eq.1)then
        if(yi.gt.-1+tiny)then
          tmp=-s/(1+yi)
        else
          tmp=0.d0
        endif
c incoming parton (right)
      elseif(ileg.eq.2)then
        if(yi.lt.1-tiny)then
          tmp=s/(1-yi)
        else
          tmp=0.d0
        endif
c outgoing top
      elseif(ileg.eq.3)then
        xw1=-xq1q+xq2q-xtk
        xw2=xq1q-xq2q-xuk
        xmq2=xm12
        xm2red=xmq2/s
        if(1-x.lt.tiny)then
          tmp=s/(1-yj)
        else
          zeta1=xw2/(s-xw1-xmq2)
          dzeta1dw2=zeta1/xw2
          dq1cdx=-(1-yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dq2qdx=-(1+yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dw2dx=(1-yj)*(s*(1+yj-x*(2*(1+yj)+x*(1-yj)))+2*xm12)/
     #          (1+yj+x*(1-yj))**2
          dw1dx=dq1cdx+dq2qdx
          dq1cdy=(1-x)*(1-yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dq2qdy=(1-x)*(1+yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dw2dy=-2*(1-x)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dw1dy=dq1cdy+dq2qdy
          tmp=(dw1dx*dw2dy-dw1dy*dw2dx)*dzeta1dw2/(zeta1*(1-zeta1))
        endif
c outgoing light antiquark
      elseif(ileg.eq.4)then
        xw1=-xq1q+xq2q-xtk
        xw2=xq1q-xq2q-xuk
        xmq2=xm12
        xm2red=xmq2/s
        if(1-x.lt.tiny)then
          tmp=-s*(1-xm2red)/(1+yj)
        elseif(1-yj.lt.tiny)then
          tmp=-s*(1-xm2red)/2.d0
        else
          eps1=1+xmq2/(s-xw2)
          beta1=sqrt(eps1**2-4.d0*s*xmq2/(s-xw2)**2)
          zeta2num=
     #       (2*s-(s-xw2)*eps1)*xw1+(s-xw2)*((xw1+xw2)*beta1-eps1*xw2)
          zeta2den=(s-xw2)*beta1*(2*s-(s-xw2)*(eps1-beta1))
          zeta2=zeta2num/zeta2den
          dzeta2dw1=1.d0/((s-xw2)*beta1)
          dq1cdx=-(1-yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dq2qdx=-(1+yi)*(s*(1+yj)+xm12*(1-yj))/(1+yj+x*(1-yj))**2
          dw2dx=(1-yj)*(s*(1+yj-x*(2*(1+yj)+x*(1-yj)))+2*xm12)/
     #          (1+yj+x*(1-yj))**2
          dw1dx=dq1cdx+dq2qdx
          dq1cdy=(1-x)*(1-yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dq2qdy=(1-x)*(1+yi)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dw2dy=-2*(1-x)*(s*x-xm12)/(1+yj+x*(1-yj))**2
          dw1dy=dq1cdy+dq2qdy         
          tmp=-(dw1dx*dw2dy-dw1dy*dw2dx)*dzeta2dw1/(zeta2*(1-zeta2))
        endif
      else
        write(*,*)'xjac_zqttoxy: unknown parton number'
        stop
      endif
c
      xjac_zqttoxy=abs(tmp)
      return
      end 
c
c
c End of utility routines for q2til, z, and 2-->2 invariants
c
c
c
c
c Begin of utility routines for Bjorken x's
c
c
      function x1soft(xx1,xx2,xx,yy)
      implicit none
      real*8 x1soft,xx1,xx2,xx,yy,tiny,x1,x2,x,y,csi,rx,tmp,xa,xb
      parameter (tiny=1.d-5)
      integer iprespl
      common/ciprespl/iprespl
c
      x1=xx1
      x2=xx2
      x=xx
      y=yy
      if(iprespl.eq.0)then
        csi=sqrt( (2-(1-x)*(1+y))/(2-(1-x)*(1-y)) )
        rx=sqrt(x)
        tmp=x1*csi*rx
      elseif(iprespl.eq.1)then
        if(1-x.lt.tiny)then
          tmp=x1*(1-(1-x)*(1+y)/2.d0)
        elseif(1-y.lt.tiny)then
          tmp=x*x1*(1+(1-x)*(1-y)*(x1+x2)/(2.d0*(x*x1+x2)))
        elseif(1+y.lt.tiny)then
          tmp=x1*(1-(1-x)*(1+y)*(x1+x2)/(2.d0*(x1+x*x2)))
        else
          xa=x*x1*x2
          xb=0.5d0*((1-x)*y*(x1+x2)+(1+x)*(x2-x1))
          tmp=0.5d0*(sqrt(xb**2+4*xa)-xb)
        endif
      else
        write(*,*)'Error in x1soft',iprespl
        stop
      endif
c Do not use negative values in order to avoid crashes
      if(tmp.le.0.d0)tmp=2.d0
      x1soft=tmp
      return
      end


      function x2soft(xx1,xx2,xx,yy)
      implicit none
      real*8 x2soft,xx1,xx2,xx,yy,tiny,x1,x2,x,y,csi,rx,tmp,xa,xb
      parameter (tiny=1.d-5)
      integer iprespl
      common/ciprespl/iprespl
c
      x1=xx1
      x2=xx2
      x=xx
      y=yy
      if(iprespl.eq.0)then
        csi=sqrt( (2-(1-x)*(1+y))/(2-(1-x)*(1-y)) )
        rx=sqrt(x)
        tmp=x2*rx/csi
      elseif(iprespl.eq.1)then
        if(1-x.lt.tiny)then
          tmp=x2*(1-(1-x)*(1-y)/2.d0)
        elseif(1-y.lt.tiny)then
          tmp=x2*(1-(1-x)*(1-y)*(x1+x2)/(2.d0*(x*x1+x2)))
        elseif(1+y.lt.tiny)then
          tmp=x*x2*(1+(1-x)*(1+y)*(x1+x2)/(2.d0*(x1+x*x2)))
        else
          xa=x*x1*x2
          xb=0.5d0*((1-x)*y*(x1+x2)+(1+x)*(x2-x1))
          tmp=0.5d0*(sqrt(xb**2+4*xa)+xb)
        endif
      else
        write(*,*)'Error in x2soft',iprespl
        stop
      endif
c Do not use negative values in order to avoid crashes
      if(tmp.le.0.d0)tmp=2.d0
      x2soft=tmp
      return
      end


      function x1x2jac(xx1,xx2,xx,yy,iileg)
      implicit none
      real*8 x1x2jac,xx1,xx2,xx,yy,tiny,x1,x2,x,y,tmp,xa,xb
      parameter (tiny=1.d-5)
      integer iileg,ileg,iprespl
      common/ciprespl/iprespl
c
      x1=xx1
      x2=xx2
      x=xx
      y=yy
      ileg=iileg
      if(ileg.eq.1.or.ileg.eq.2)then
        if(iprespl.eq.0)then
          tmp=x
        elseif(iprespl.eq.1)then
          if(1-x.lt.tiny)then
            tmp=x+(1-x)**2*(1-y**2)*x1*x2/(2.d0*(x1+x2)**2)
          elseif(1-y.lt.tiny)then
            tmp=x+(1-x)**2*x*(1-y)*x1*x2/(x*x1+x2)**2
          elseif(1+y.lt.tiny)then
            tmp=x+(1-x)**2*x*(1+y)*x1*x2/(x1+x*x2)**2
          else
            xa=x*x1*x2
            xb=0.5d0*((1-x)*y*(x1+x2)+(1+x)*(x2-x1))
            tmp=x*((1-y+x*(1+y))*x1+(1+y+x*(1-y))*x2)/
     #          (2.d0*sqrt(xb**2+4*xa))
          endif
        else
          write(*,*)'Error # 1 in x1x2jac',iprespl
          stop
        endif
      elseif(ileg.eq.3.or.ileg.eq.4)then
        tmp=1.d0
      else
        write(*,*)'Error # 2 in x1x2jac',ileg
        stop
      endif
      x1x2jac=abs(tmp)
      return
      end
c
c
c Running couplings
c
c
      function zgmu2_nlo()
c Sets the scales for NLO subtraction terms
      implicit none
      real * 8 zgmu2_nlo
      real * 8 pi,ptv1,ptv2,ptvg
      common/perpen/ptv1(2),ptv2(2),ptvg(2)
      parameter (pi=3.14159265358979312D0)
      include 'stpcblks.h'
      real * 8 ptsum,xmu2,xmu2a,xmu2b,as,alfas
      integer inloscale
      common/cinloscale/inloscale
      real * 8 scptveto
      common/cscptveto/scptveto
c
      if(mod(inloscale,10).eq.1)then
        ptsum = ptv1(1)**2+ptv1(2)**2+
     #          ptv2(1)**2+ptv2(2)**2
        xmu2 = (ptsum+xm12+xm22)/2.d0
      elseif(mod(inloscale,10).eq.2)then
        xmu2 = xm12
      elseif(mod(inloscale,10).eq.3)then
        xmu2 = (xm12+xm22)/2.d0
      elseif(mod(inloscale,10).eq.4)then
        xmu2 = ( (sqrt(xm12)+sqrt(xm22))/2.d0 )**2
      else
        write(*,*)'Unknown option in zgmu2_nlo',inloscale
        stop
      endif
      if( inloscale.le.10 .or. 
     #    (inloscale.gt.20.and.inloscale.le.30) )then
        xmu2b = xmu2
      else
        xmu2b = scptveto**2
      endif
      if(inloscale.le.20)then
        xmu2a = xmu2
      else
        xmu2a = scptveto**2
      endif
c set the factorization scales for hadron 1 and 2, and the
c renormalization scale
      xmuf2h1 = xmu2b*xf2h1
      xmuf2h2 = xmu2b*xf2h2
      xmur2  = xmu2a*xren2
      as    = alfas(xmur2,xlam,nl)
      zgmu2_nlo = 4.d0*pi*as
      zg = sqrt(zgmu2_nlo)
      end


      function zgmu2_mc()
c Sets the scales for MC subtraction terms
      implicit none
      real * 8 zgmu2_mc
      real * 8 pi,ptv1,ptv2,ptvg
      common/perpen/ptv1(2),ptv2(2),ptvg(2)
      parameter (pi=3.14159265358979312D0)
      include 'stpcblks.h'
      real * 8 ptsum,xmu2,xmu2a,xmu2b,as,alfas
      integer imcscale
      common/cimcscale/imcscale
      real * 8 scptveto
      common/cscptveto/scptveto
c
      if(mod(imcscale,10).eq.1)then
        ptsum = ptv1(1)**2+ptv1(2)**2+
     #          ptv2(1)**2+ptv2(2)**2
        xmu2 = (ptsum+xm12+xm22)/2.d0
      elseif(mod(imcscale,10).eq.2)then
        xmu2 = xm12
      elseif(mod(imcscale,10).eq.3)then
        xmu2 = (xm12+xm22)/2.d0
      elseif(mod(imcscale,10).eq.4)then
        xmu2 = ( (sqrt(xm12)+sqrt(xm22))/2.d0 )**2
      else
        write(*,*)'Unknown option in zgmu2_mc',imcscale
        stop
      endif
      if( imcscale.le.10 .or. 
     #    (imcscale.gt.20.and.imcscale.le.30) )then
        xmu2b = xmu2
      else
        xmu2b = scptveto**2
      endif
      if(imcscale.le.20)then
        xmu2a = xmu2
      else
        xmu2a = scptveto**2
      endif
c set the factorization scales for hadron 1 and 2, and the
c renormalization scale
      xmumcf2h1 = xmu2b*xf2h1mc
      xmumcf2h2 = xmu2b*xf2h2mc
      xmumcr2  = xmu2a*xren2mc
      as    = alfas(xmumcr2,xlam,nl)
      zgmu2_mc = 4.d0*pi*as
      zg = sqrt(zgmu2_mc)
      end


c-------------------------------------------------------------------------
      function xalfaem(q2)
c Alpha_em(MSbar) at the scale q2 = q^2. 
c Uses alpha_Thomson below the electron mass, an interpolation between
c m_e and m_tau, and the evolution equation above m_b. This function is
c taken from the gamma*gamma* --> hadrons package
c-------------------------------------------------------------------------
      implicit real*8 (a-z)
      integer npoints,ideg
      parameter (npoints=3,ideg=3)
      real*4 ooa(npoints),xlogmu(npoints),divdif
c 1/alpha_em at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data ooa     / 137.036, 135.95, 133.513 /
c logs of sqrt(q2) at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data xlogmu  / -7.57914, -2.2481, 0.574927 /
      data zm/91.2d0/,ooaz/127.934d0/,pi/3.1415927d0/,nc/3/
c
      if(q2.lt.exp(2.*xlogmu(1))) then
         xalfaem = 1.d0/ooa(1)	 
      elseif(q2.lt.exp(2.*xlogmu(3))) then
         xlogq = log(q2)/2.d0
         xalfaem = 1.d0/divdif(ooa,xlogmu,npoints,sngl(xlogq),ideg)
      elseif(q2.lt.5.**2) then
         b = 3 + 2*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2) - 2.*xlogmu(3)
         xalfaem = 1d0/ooa(3)/(1.d0 - 1.d0/3.d0/pi/ooa(3)*b*xlq)
      else
         b = 3 + 3*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2/zm**2)
         xalfaem = 1d0/ooaz/(1.d0 - 1.d0/3.d0/pi/ooaz*b*xlq)
      endif
      return
      end
