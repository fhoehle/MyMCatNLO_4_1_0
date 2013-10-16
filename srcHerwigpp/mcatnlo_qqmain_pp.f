      PROGRAM MCATNLO_QQMAIN_PP
c Integrates heavy quark pair cross sections, and produces the event
c file which serves as input to Herwig. Quantities relevant to H and S
c events are indentified with letters a and b respectively
      implicit none
      real * 8 
     #  value(20),ecmlst(100),
     #  xmlst(100),fh1lst(100),fh2lst(100),renlst(100),
     #  fh1mclst(100),fh2mclst(100),renmclst(100),
     #  xlep1mass(2),xlep2mass(2),xmomshifts(4),vickm(1:6,1:6),
     #  alsf(3),besf(3),alcl(3),becl(3),alazi(3),beazi(3)
      real * 8 
     #  xm,xpdflam4,xpdflam5,tmp,xfh,xfh1,xfh2,ecm,xren,betfac,
     #  delta,deltas,deltac,dtot,avtot,ac1,ac2,xtotal,ytotal,av3a,
     #  av3nega,d3a,d3nega,av3b,av3negb,d3b,d3negb,ctime,pi,tmas,
     #  etacut,xalsf,xbesf,xalcl,xbecl,xalazi,xbeazi,wgtaev,wgtbev,
     #  dummy,xmone,xfhmc,xrenmc,evfrac,evprcfrac,xares,yares,xbres,
     #  ybres,xmw,gaw,xmt,twidth,xm012,ga1,bw1delf,bw1fmmn,xm022,ga2,
     #  bw2delf,bw2fmmn,xm1low2,xm1upp2,xm2low2,xm2upp2,brrtop1,
     #  brrtop2,xmw2,gammax1,xm1low,xm1upp,gammax2,xm2low,xm2upp,
     #  bw1mdpl,bw1mdmn,bw1fmpl,bw2mdpl,bw2mdmn,bw2fmpl,al_spcfun,
     #  be_spcfun,tga1,tga1mn,tga1pl,tbw1fmpl,tbw1fmmn,tbw1delf,
     #  ym1low2,ym1upp2,tga2,tga2mn,tga2pl,tbw2fmpl,tbw2fmmn,tbw2delf,
     #  ym2low2,ym2upp2,xmt2,gammay1,ym1low,ym1upp,gammay2,ym2low,
     #  ym2upp,tbw1mdpl,tbw1mdmn,tbw2mdpl,tbw2mdmn,wgtmax,xbrrtoplep,
     #  xbrrtophad,xsecpp,xerrpp
      integer 
     #  ih1,ih2,ndns1,ndns2,jloop,iseld,nlf,ncl3,jecm,
     #  loproc,maproc,iproc,iinput,iverbose,ichkmom,
     #  ibswrite,itmpih,itmpndns,idpdfset,ipdfih,ipdfgroup,ipdfndns,
     #  ifk88istrl,ifk88ih,ifk88ndns,maxevt,it1,it2,ifuntype,
     #  ndim,nwild,itd1,itd2,ibscall,iwgtnorm,iseed0,
     #  iseed,maxtrials,mode,lo,isubttype,iprespl,iwrong,iwrong1,
     #  ntotal,ndiff,nevts,ntrls,itot,ionshell,iunita,iunitb,
     #  ioutput,ii,iunit,i,j,itmpqq,itmpvv,mx_of_evta,mx_of_evtb,
     #  nlfp1sch,nlas,ia1ora2,iasmc,iassoft,ifk88seed,ip4,ip5,
     #  ip6,ip7,ip8,ip9,izero,ione,idec,iwidth,il1hw,il2hw,
     #  neventsuw,ifailuw,ncntuws,nmaxuw,nqmaxuw,nqeventsuw,
     #  nqcntuws,ideconsh,jwidth,mqeventsuw,inonbtop,ievffmt
      character * 2 scheme,xproc(3)
      character * 4 part1,part2
      character * 20 parm(20),gname
      character * 80 fname,fnamea,fnameb,fname1,fnamev
      character * 80 pref,prefn,prefev,prefnev
      character * 70 strin,strout,lhapdf
      logical iphflag,evgen
      external sig5a,sig5b
      parameter (pi=3.14159265358979312D0)
      parameter (xmone=-1.d0)
      parameter (izero=0)
      parameter (ione=1)
      include 'hvqcblks.h'
c
c common /strfun0/ is only in strfun:
c ndns = pdf type
c ih1,ih2 = beam type (0=(p+n)/2, 1=p, -1=pbar, 2=n, -2=nbar)
      common/strfun0/ih1,ih2,ndns1,ndns2
c quark and gluon masses, used by Herwig. PDF labeling convention
      real*8 xmass(-5:21)
      common/parmass/xmass
c CKM matrix elements entered by the user
      common/cvickm/vickm
c alsf and besf are the parameters that control gfunsoft
      common/cgfunsfp/alsf,besf
c alcl and becl are the parameters that control gfuncoll (unused)
      common/cgfunclp/alcl,becl
c alazi and beazi are the parameters that control gfunazi
      common/cgfunazi/alazi,beazi
c al_spcfun, be_spcfun are the parameters entering spcdamp
      common/cspcpar/al_spcfun,be_spcfun
c iwgtnorm=0 for weight=+1/-1, iwgtnorm=1 otherwise
      common/ciwgtnorm/iwgtnorm
c number of events generated
      common/cmaxevt/maxevt
c ievffmt=0 for MC@NLO event file format, ievffmt=1 for LHEF format
      common/cievffmt/ievffmt
c wgtaev and wgtbev are the norms of weights for H and S events respectively
      common/cwgtev/wgtaev,wgtbev
c wgtmax is the maximum absolute value of the weights
      common/cwgtmax/wgtmax
c iprespl=0 ==> preserves rapidity
c iprespl=1 ==> preserves longitudinal momentum
      common/ciprespl/iprespl
c ichkmom=0 --> enables checks on kinematics
      common/cichkmom/ichkmom
c----------------------------------------------------------
c Variables that control the integrations
c
      common/cisubttype/isubttype
      common/betfac/betfac,delta
      common/pmerge/deltas,deltac
c etacut is the maximum allowed for [2*kt(gluon)/sqrt(shat)]^2
      common/cetacut/etacut
      integer nsamp
      common/samp/nsamp
c----------------------------------------------------------
c Top decay variables
c Decay of the tops: idec=0    -->   tops decay
c                    idec=1    -->   tops don't decay, or b production
      common/cidec/idec
c Top mass ranges: jwidth=0    -->   top on shell
c                  jwidth=1    -->   top off shell
      common/cjwidth/jwidth
c Top mass ranges
      common/tbw1/tga1,tga1mn,tga1pl,tbw1fmpl,tbw1fmmn,tbw1delf,
     #            ym1low2,ym1upp2
      common/tbw2/tga2,tga2mn,tga2pl,tbw2fmpl,tbw2fmmn,tbw2delf,
     #            ym2low2,ym2upp2
c W mass ranges: iwidth=0    -->   W on shell
c                iwidth=1    -->   W off shell
      common/ciwidth/iwidth
c Type of W decays; il1hw and il2hw are entered with the following conventions,
c which control the decay of the W (top always decays into Wb)
c  IL=0     ==> all W decays (quark+leptons)
c  IL=1,2,3 ==> W -> e\nu_e, mu\nu_mu, tau\nu_tau
c  IL=4     ==> W -> e\nu_e + mu\nu_mu
c  IL=5     ==> W -> all quarks
c  IL=6     ==> W -> e\nu_e + mu\nu_mu + all quarks (ie all decays except tau)
c  IL=7     ==> the top does not decay
      common/cilhw/il1hw,il2hw
c W mass and width
      common/cwparam/xmw,gaw
c top mass and width; top mass squared is stored in fixvar; xmt must
c be used only in those parts of the code relevant to top decay
      common/ctparam/xmt,twidth
c W mass ranges
      common/cbw1/xm012,ga1,bw1delf,bw1fmmn
      common/cbw2/xm022,ga2,bw2delf,bw2fmmn
      common/bounds/xm1low2,xm1upp2,xm2low2,xm2upp2
c top branching ratios, for lepton and hadron decays; apart from testing 
c purposes, these should be about 0.111 and 0.333 respectively,
c ie for W->e nu_e and W->udbar+usbar+ubar
      common/xibrratios/xbrrtoplep,xbrrtophad
c reweight factors when tops decay, that include branching ratios,
c for top and antitop
      common/brratios/brrtop1,brrtop2
c mass of particles from W decays: not necessarily leptons!
      common/clepmass/xlep1mass,xlep2mass
c Decay of the tops: inonbtop=0    -->   t->Wb only
c                    inonbtop=1    -->   t->W+any down-type quark
      common/cinonbtop/inonbtop
c----------------------------------------------------------
c nlfp1sch=0 --> use nl light flavours, nlfp1sch=1 --> nl+1 scheme
      common/cnlfp1sch/nlfp1sch
c nlas is the number of light flavours used in the computation of alpha_S
      common/cnlas/nlas
c Identities of final-state particles, except for the light parton (included
c as ip3 in ci1part), according to MC particle numbering scheme; used 
c when writing the event file
      common/ci2part/ip4,ip5,ip6,ip7,ip8,ip9
c----------------------------------------------------------
c The following refer to the computation of MC subtraction terms
c ia1ora2=1 -> full invariants, ia1ora2=2 -> simplified invariants
      common/cia1ora2/ia1ora2
c iasmc=1 -> as_nlo**3, iasmc=2 -> as_mc**3, iasmc=3 -> as_mc*as_nlo**2
      common/ciasmc/iasmc
c iassoft=0 -> as_nlo(hard), iassoft=1 -> as_nlo(soft)
      common/ciassoft/iassoft
c----------------------------------------------------------
c Subprocesses: prc = 'gg', 'qq', 'qg', corresponding to jproc=1,2,3
c and equal to xproc(jproc)
      common/cxproc/xproc
c In the integration routines, loproc<=jproc<=maproc
      common/cwchproc/loproc,maproc
c Number of failures in flavour determination
      common/ciwrong/iwrong,iwrong1
c Common blocks for statistics relevant to secondary unweighting
      common/c1iunwgt/neventsuw,nqeventsuw,mqeventsuw,ifailuw
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
c Common blocks for general MC@NLO routines
c common block for internal rnd number generation, independent of bases
      common/cifk88seed/ifk88seed
c common block fk88ipdfs is filled by our interface to MLMPDF
      common/fk88ipdfs/ifk88ih,ifk88ndns
c common block w50511 and w50512 are filled by PDFLIB 
      common/w50511/ipdfih,ipdfgroup,ipdfndns,mode,nlf,lo,tmas
      common/w50512/xpdflam4,xpdflam5
c----------------------------------------------------------
c- list of subprocesses
      data xproc/'gg','qq','qg'/
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
      open(unit=11,file='hvqlog',status=newver)
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
c
c-----------------------------------------------------------------
c Kept for backward compatibility; when set to false there's 
c no photon around
      iphflag=.false.
c----------------------------------------------------------
c Parameters of the run
c
      if(iphflag) then
         write(*,*)'Only hadron collisions are implemented'
         stop
      else
         write(*,*)' '
         write(*,*)
     # 'Enter pair ECM(GeV),fren[NLO],ffact[NLO],fren[MC],ffact[MC]'
         write(*,*)' fren=mu_ren/mu0'
         write(*,*)' ffact=mu_fac/mu0'
         write(*,*)' mu_ren=renormalization scale'
         write(*,*)' mu_fac=factorization scale'
         write(*,*)' mu0=reference scale'
         jecm = 0
         ecm = 1.d-8
c Allow only one set of entries in this version; disallow negative
c entries, energies of the two beams will be entered separately eventually
         dowhile(jecm.lt.1.and.ecm.gt.0)
            read(*,*) ecm,xren,xfh,xrenmc,xfhmc
            write(11,'(5(1x,d10.4),1x,a)') ecm,xren,xfh,xrenmc,xfhmc,
     #    '! energy, fren, ffact, frenmc, ffactmc'
            jecm=jecm+1
            ecmlst(jecm)=ecm
            fh1lst(jecm)=xfh
            fh2lst(jecm)=xfh
            renlst(jecm)=xren
c Will use xfhmc and xrenmc in future versions
            fh1mclst(jecm)=xfh
            fh2mclst(jecm)=xfh
            renmclst(jecm)=xren
         enddo
         if(jecm.eq.100.and.xm.gt.0) then
            write(*,*) 'no more than 100 values'
            stop
         endif
      endif
c
c Process number (redundant with mass -- keep it for consistency with Herwig)
      write(*,*)' '
      write(*,*)'Enter -1705 for b-bbar production'
      write(*,*)'      -1706 for t-tbar'
      read(*,*) itmpvv
      if(itmpvv.ne.-1705 .and. itmpvv.ne.-1706 .and.
     #   itmpvv.ne.-11705 .and. itmpvv.ne.-11706 )then
         write(*,*) 'Error: wrong IPROC'
         stop
      endif
      write(11,'(1x,i6,27x,a)') itmpvv,'! -1705/1706=b/t'
c
      if(itmpvv.eq.-1705 .or. itmpvv.eq.-11705)then
        write(*,*)'b-bbar production not available with HW++'
        stop
      endif
c
c Heavy quark mass
      write(*,*)' '
      if(itmpvv.eq.-1705.or.itmpvv.eq.-11705)then
        write(*,*)'Enter the bottom mass (GeV)'
      elseif(itmpvv.eq.-1706.or.itmpvv.eq.-11706)then
        write(*,*)'Enter the top mass (GeV)'
      endif
      read(*,*)xm
      write(11,'(1x,d10.4,23x,a)') xm,'! M_Q'
      if(jecm.ne.1)then
        write(*,*)'Fatal error: multiple inputs',jecm
        stop
      endif
      xmlst(jecm)=xm                             
c
c Top decay parameters
      if(itmpvv.eq.-1706.or.itmpvv.eq.-11706)then
        write(*,*)' '
        write(*,*)'Enter IL=0..6 for t->W(->d1_IL d2_IL) b'
        write(*,*)'      IL=7 for undecayed tops'
        write(*,*)'for W+ and W- from top and tbar'
        read(*,*) il1hw,il2hw
        write(11,'(1x,i2,1x,i2,28x,a)') il1hw,il2hw,
     #    '! 0..6 -> t dec, 7 -> t undec'
        if( (il1hw.eq.7.and.il2hw.ne.7) .or.
     #      (il1hw.ne.7.and.il2hw.eq.7) )then
          write(*,*) 'Ws must both decay or being stable'
          stop
        elseif(il1hw.eq.7.and.il2hw.eq.7)then
          idec=1
        elseif( (il1hw.ge.0.and.il1hw.le.6) .and.
     #          (il2hw.ge.0.and.il2hw.le.6) )then
          idec=0
        else
          write(*,*) 'Unknown options:',il1hw,il2hw
          stop
        endif
        if(idec.eq.0)then
          xmt=xm
          xmt2=xmt**2
          write(*,*)' '
          write(*,*)'Enter top width'
          read(*,*)twidth
          write(11,'(1x,d10.4,23x,a)') twidth,'! top width'
c
          write(*,*)' '
          write(*,*)'Enter W mass and width (GeV)'
          read(*,*)xmw,gaw
          write(11,'(2(1x,d10.4),12x,a)') xmw,gaw,'! M_W, Gamma_W'
          xmw2=xmw**2
c
          write(*,*)' '
          write(*,*)'Enter GammaX, M_T(min), M_T(max) for top'
          write(*,*)
     #     '  If GammaX>0, the top mass is chosen in the range'
          write(*,*)'      M0-GammaX*width < M_T < M0+GammaX*width'
          write(*,*)'  and M_T(min), M_T(max) are ignored'
          write(*,*)
     #     '  If GammaX<0, the top mass is chosen in the range'
          write(*,*)'            M_T(min) < M_T < M_T(max)'
          write(*,*)
     #  '  If GammaX=0, the top mass is set equal to the pole mass'
          read(*,*)gammay1,ym1low,ym1upp
          write(11,'(3(1x,d10.4),1x,a)') gammay1,ym1low,ym1upp,
     #     '! GammaX, M_T(min), M_T(max)'
          if(gammay1.lt.0.and.ym1low.ge.ym1upp)then
            write(*,*)'Enter a non-zero range'
            stop
          endif
c
          write(*,*)' '
          write(*,*)'Enter GammaX, M_Tb(min), M_Tb(max) for tbar'
          write(*,*)
     #     '  If GammaX>0, the tbar mass is chosen in the range'
          write(*,*)'      M0-GammaX*width < M_Tb < M0+GammaX*width'
          write(*,*)'  and M_Tb(min), M_Tb(max) are ignored'
          write(*,*)
     #     '  If GammaX<0, the tbar mass is chosen in the range'
          write(*,*)'            M_Tb(min) < M_Tb < M_Tb(max)'
          write(*,*)
     #  '  If GammaX=0, the tbar mass is set equal to the pole mass'
          read(*,*)gammay2,ym2low,ym2upp
          write(11,'(3(1x,d10.4),1x,a)') gammay2,ym2low,ym2upp,
     #     '! GammaX, M_Tb(min), M_Tb(max)'
          if(gammay2.lt.0.and.ym2low.ge.ym2upp)then
            write(*,*)'Enter a non-zero range'
            stop
          endif
c
          write(*,*)' '
          write(*,*)'Enter GammaX, M_V1(min), M_V1(max) for W+'
          write(*,*)
     #     '  If GammaX>0, the boson mass is chosen in the range'
          write(*,*)'      M0-GammaX*width < M_W+ < M0+GammaX*width'
          write(*,*)'  and M_V1(min), M_V1(max) are ignored'
          write(*,*)
     #     '  If GammaX<0, the boson mass is chosen in the range'
          write(*,*)'            M_V1(min) < M_W+ < M_V1(max)'
          write(*,*)
     #  '  If GammaX=0, the boson mass is set equal to the pole mass'
          read(*,*)gammax1,xm1low,xm1upp
          write(11,'(3(1x,d10.4),1x,a)') gammax1,xm1low,xm1upp,
     #     '! GammaX, M_V1(min), M_V1(max)'
          if(gammax1.lt.0.and.xm1low.ge.xm1upp)then
            write(*,*)'Enter a non-zero range'
            stop
          endif
c
          write(*,*)' '
          write(*,*)'Enter GammaX, M_V2(min), M_V2(max) for W-'
          write(*,*)
     #     '  If GammaX>0, the boson mass is chosen in the range'
          write(*,*)'      M0-GammaX*width < M_W- < M0+GammaX*width'
          write(*,*)'  and M_V2(min), M_V2(max) are ignored'
          write(*,*)
     #     '  If GammaX<0, the boson mass is chosen in the range'
          write(*,*)'            M_V2(min) < M_W- < M_V2(max)'
          write(*,*)
     #  '  If GammaX=0, the boson mass is set equal to the pole mass'
          read(*,*)gammax2,xm2low,xm2upp
          write(11,'(3(1x,d10.4),1x,a)') gammax2,xm2low,xm2upp,
     #    '! GammaX, M_V2(min), M_V2(max)'
          if(gammax2.lt.0.and.xm2low.ge.xm2upp)then
            write(*,*)'Enter a non-zero range'
            stop
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
          write(11,'(3(1x,d10.4),1x,a)')
     #      vickm(1,2),vickm(1,3),vickm(1,5),
     #      '! |V_ud|,|V_us|,|V_ub|'
          write(*,*)'Enter |V_cd|, |V_cs|, |V_cb|'
          read(*,*)vickm(4,2),vickm(4,3),vickm(4,5)
          write(11,'(3(1x,d10.4),1x,a)')
     #      vickm(4,2),vickm(4,3),vickm(4,5),
     #      '! |V_cd|,|V_cs|,|V_cb|'
          write(*,*)'Enter |V_td|, |V_ts|, |V_tb|'
          read(*,*)vickm(6,2),vickm(6,3),vickm(6,5)
          write(11,'(3(1x,d10.4),1x,a)')
     #      vickm(6,2),vickm(6,3),vickm(6,5),
     #      '! |V_td|,|V_ts|,|V_tb|'
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
     #      '! t -> leptons branching ratio'
          write(*,*)' '
          write(*,*)'Enter top -> hadrons branching ratio'
          read(*,*)xbrrtophad
          write(11,'(1x,d10.4,23x,a)') xbrrtophad,
     #      '! t -> hadrons branching ratio'
c Redefine width or branching ratios if need be
          call reset_twdbr(xmt,twidth,xbrrtoplep,xbrrtophad,
     #                     gammax1,gammax2)
        else
          xmt=0.d0
          twidth=0.d0
          xmw=0.d0
          gaw=0.d0
          inonbtop=-1
        endif
      else
        idec=1
      endif
c
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
c Initialize parameters, such as labelling for parton processes
      call parsetpar()
c
      write(*,*)' '
      write(*,*)
     #  'Enter beam type for beam1 and beam2 (p, pbar, n, nbar):'
      read(*,*) part1,part2
      write(11,'(1x,a,2x,a,19x,a)') ''''//part1//'''',
     #  ''''//part2//'''','! colliding hadron types'
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
         write(*,*) 'Enter Lambda_QCD_5 (GeV), < 0 for default'
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
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)
        write(*,*)'Enter the number n of light flavours (4 or 5);'
        write(*,*)
     #    ' a negative entry will force the code to use the default:'
        write(*,*)
     #    ' n=3 for m<3, n=4 for 3<m<7, n=5 for m>7'
        read(*,*) nlf
        write(11,'(1x,i2,31x,a)') nlf,'! # of light flavours'
      else
        nlf=-1
      endif
c nl+1 or nl scheme; the former is only used for bottom production
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 for the nl scheme'
        write(*,*)'      1 for the nl+1 scheme'
        read(*,*) nlfp1sch
        write(11,'(1x,i2,31x,a)') nlfp1sch,
     #    '! 0 for nl, 1 for nl+1 scheme'
      else
        nlfp1sch=-1
      endif
c-----------------------------------------------------------------
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_soft[gg,qq]'
        write(*,*)' Defaults are: alpha=1, beta=-0.1'
        read(*,*) xalsf,xbesf
        write(11,'(2(2x,d10.4),10x,a)') xalsf,xbesf,
     #    '! alpha, beta [soft;gg,qq]'
      else
        xalsf=1.d0
        xbesf=-0.1d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_coll[gg,qq]'
        write(*,*)' Defaults are: alpha=0, beta=0'
        read(*,*) xalcl,xbecl
        write(11,'(2(2x,d10.4),10x,a)') xalcl,xbecl,
     #    '! alpha, beta [coll;gg,qq]'
      else
        xalcl=0.d0
        xbecl=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_azi[gg,qq]'
        write(*,*)' Defaults are: alpha=-1, beta=-0.1'
        read(*,*) xalazi,xbeazi
        write(11,'(2(2x,d10.4),10x,a)') xalazi,xbeazi,
     #    '! alpha, beta [azi;gg,qq]'
      else
        xalazi=-1.d0
        xbeazi=-0.1d0
      endif
      alsf(1)=xalsf
      besf(1)=xbesf
      alcl(1)=xalcl
      becl(1)=xbecl
      alazi(1)=xalazi
      beazi(1)=xbeazi
      alsf(2)=xalsf
      besf(2)=xbesf
      alcl(2)=xalcl
      becl(2)=xbecl
      alazi(2)=xalazi
      beazi(2)=xbeazi
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_soft[qg]'
        write(*,*)' Defaults are: alpha=-1, beta=0'
        read(*,*) xalsf,xbesf
        write(11,'(2(2x,d10.4),10x,a)') xalsf,xbesf,
     #    '! alpha, beta [soft;qg]'
      else
        xalsf=-1.d0
        xbesf=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_coll[qg]'
        write(*,*)' Defaults are: alpha=0, beta=0'
        read(*,*) xalcl,xbecl
        write(11,'(2(2x,d10.4),10x,a)') xalcl,xbecl,
     #    '! alpha, beta [coll;qg]'
      else
        xalcl=0.d0
        xbecl=0.d0
      endif
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha, beta for the function G_azi[qg]'
        write(*,*)' Defaults are: alpha=-1, beta=-0.1'
        read(*,*) xalazi,xbeazi
        write(11,'(2(2x,d10.4),10x,a)') xalazi,xbeazi,
     #    '! alpha, beta [azi;qg]'
      else
        xalazi=-1.d0
        xbeazi=-0.1d0
      endif
      alsf(3)=xalsf
      besf(3)=xbesf
      alcl(3)=xalcl
      becl(3)=xbecl
      alazi(3)=xalazi
      beazi(3)=xbeazi
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter alpha and beta for the function SPC_damp'
        write(*,*)' Defaults are: alpha=1, beta=0.5'
        write(*,*)' Allowed ranges: alpha>=1, 0<beta<=1'
        read(*,*) al_spcfun,be_spcfun
        write(11,'(2(1x,d10.4),12x,a)') al_spcfun,be_spcfun,
     #    '! alpha, beta (spin corr)'
      else
        al_spcfun=1.d0
        be_spcfun=0.5d0
      endif
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
        write(*,*)' '
        write(*,*)
     #    'For the computation of alpha_S in the MC subtraction terms'
        write(*,*)'Enter 1 to use as_nlo**3'
        write(*,*)'      2 to use as_mc**3'
        write(*,*)'      3 to use as_mc*as_nlo**2'
        write(*,*)' The default is 2'
        read(*,*) iasmc
        write(11,'(1x,i2,31x,a)') iasmc,
     #    '! 1->as_nlo^3, 2->as_mc^3, 3->as_mc*as_nlo^2'
        write(*,*)' '
        write(*,*)'When using alpha_nlo in the MC subtraction terms'
        write(*,*)'Enter 0 to use as_nlo(2->3 configuration)'
        write(*,*)'      1 to use as_nlo(soft 2->2 configuration)'
        write(*,*)' The default is 0'
        read(*,*) iassoft
        write(11,'(1x,i2,31x,a)') iassoft,
     #    '! 0->as_nlo(hard), 1->as_nlo(soft)'
      else
        ia1ora2=1
        iasmc=3
        iassoft=0
      endif
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
     #  '! 0 => wgt=+1/-1, 1 => wgt=+w/-w'
c iseed0 is the seed for the integration step, iseed is the seed
c for the event generation step
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
c-----------------------------------------------------------------
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to use standard subtraction'
        write(*,*)'      1 to use zeta subtraction'
        read(*,*)isubttype
        write(11,'(1(1x,i8),25x,a)') isubttype,
     #                               '! 0=std subt, 1=zeta subt'
      else
        isubttype=1
      endif
      if(isubttype.eq.0)then
        write(*,*)' '
        write(*,*)'Enter betfact and delta (defaults: 0.3, 0.2)'
        read(*,*)betfac,delta
        write(11,'(2(2x,d10.4),10x,a)') betfac,delta,'! betfac,delta'
      else
        write(*,*)' '
        write(*,*)'Enter zi ( [ 2*kt(gluon)/sqrt(shat) ]^2 < zi )'
        write(*,*)' Default is: zi=0.3'
        read(*,*) etacut
        write(11,'(1x,d10.4,23x,a)') etacut,'! zi'
        betfac = 1.d0
        delta = 1.d0
      endif
      deltas = 0
      deltac = 0
c 
c We should actually choose iprespl=0 for Herwig 6.5 and iprespl=1 for
c Herwig 6.4 or lower, but in any case versions different from 6.5 have
c Thomas precession not taken into account (a NNLO effect anyhow).
c We always recommend the most recent version
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*)'Enter 0 to preserve rapidity'
        write(*,*)'      1 to preserve longitudinal momentum'
        read(*,*)iprespl
        write(11,'(1(1x,i8),25x,a)') iprespl,'! 0=y, 1=k_3 preserved'
      else
        iprespl=0
      endif
c---------------------------------------------------------------
c Select subprocess
c
      if(iinput.eq.1)then
        write(*,*)' '
        write(*,*) 'Enter 1 for gg, 2 for qq, 3 for qg, 0 for all'
        write(*,*) 'to select the subprocess'
        read(*,*) iproc
        write(11,'(1x,i2,31x,a)') iproc,'! 1=gg, 2=qq, 3=qg, 0=all'
      else
        iproc=0
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
        write(*,*)'Enter 0 to leave the top decay products massless'
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
          if(ncl3.lt.0)ncl3=120000
          write(11,'(1x,i9,24x,a)')ncl3,'! # of calls for bases'
        endif
      else
        ncl3=120000
      endif
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
c When tops decay, compute the relevant parameters
      if(idec.eq.0)then
        if( (gammay1.ne.0.d0.and.twidth.eq.0.d0) .or.
     #      (gammay2.ne.0.d0.and.twidth.eq.0.d0) )then
          write(*,*)'Non-zero top mass range requires non-zero width'
          stop
        endif
        if(gammay1.eq.0.and.gammay2.eq.0)then
          jwidth=0
          ym1low2=-1.d0
          ym1upp2=-1.d0
          ym2low2=-1.d0
          ym2upp2=-1.d0
          tbw1delf=0.d0
          tbw2delf=0.d0
        elseif(gammay1.ne.0.and.gammay2.ne.0)then
          jwidth=1
          tga1=twidth
          tga2=twidth
          if(gammay1.ge.0)then
            ym1low2=(max( 10.d0,xlep1mass(1)+xlep2mass(1)+xmass(5),
     #                    xmt-gammay1*tga1 ))**2
            ym1upp2=(min(ecmlst(1)-10.d0,xmt+gammay1*tga1))**2
          else
            ym1low2=(max( 10.d0,xlep1mass(1)+xlep2mass(1)+xmass(5),
     #                    ym1low) )**2
            ym1upp2=(min(ecmlst(1)-10.d0,ym1upp))**2
          endif
          if(ym1low2.gt.ym1upp2)then
            write(*,*)'Error in top mass range #1'
            write(*,*)ym1low2,ym1upp2
            stop
          endif
c Parameters for the skewed Breit Wigner function
          tga1mn=tga1
          tga1pl=1.15d0*tga1
          tbw1mdpl=ym1upp2-xmt2
          tbw1mdmn=xmt2-ym1low2
          tbw1fmpl=tga1pl/tga1*atan(tbw1mdpl/(xmt*tga1pl))
          tbw1fmmn=tga1mn/tga1*atan(tbw1mdmn/(xmt*tga1mn))
          tbw1delf=(tbw1fmpl+tbw1fmmn)/pi
c
          if(gammay2.ge.0)then
            ym2low2=(max( 10.d0,xlep1mass(2)+xlep2mass(2)+xmass(5),
     #                    xmt-gammay2*tga2 ))**2
            ym2upp2=(min(ecmlst(1)-10.d0,xmt+gammay2*tga2))**2
          else
            ym2low2=(max( 10.d0,xlep1mass(2)+xlep2mass(2)+xmass(5),
     #                    ym2low) )**2
            ym2upp2=(min(ecmlst(1)-10.d0,ym2upp))**2
          endif
          if(ym2low2.gt.ym2upp2)then
            write(*,*)'Error in top mass range #2'
            write(*,*)ym2low2,ym2upp2
            stop
          endif
c Parameters for the skewed Breit Wigner function
          tga2mn=tga2
          tga2pl=1.15d0*tga2
          tbw2mdpl=ym2upp2-xmt2
          tbw2mdmn=xmt2-ym2low2
          tbw2fmpl=tga2pl/tga2*atan(tbw2mdpl/(xmt*tga2pl))
          tbw2fmmn=tga2mn/tga2*atan(tbw2mdmn/(xmt*tga2mn))
          tbw2delf=(tbw2fmpl+tbw2fmmn)/pi
        else
          write(*,*)'Both mass ranges must be non-zero'
          stop
        endif
c
        if( (gammax1.ne.0.d0.and.gaw.eq.0.d0) .or.
     #      (gammax2.ne.0.d0.and.gaw.eq.0.d0) )then
          write(*,*)'Non-zero W mass range requires non-zero width'
          stop
        endif
        xm012=xmw2
        xm022=xmw2
        if(gammax1.eq.0.and.gammax2.eq.0)then
          iwidth=0
          xm1low2=-1.d0
          xm1upp2=-1.d0
          xm2low2=-1.d0
          xm2upp2=-1.d0
          bw1delf=0.d0
          bw2delf=0.d0
        elseif(gammax1.ne.0.and.gammax2.ne.0)then
          iwidth=1
          ga1=gaw
          ga2=gaw
          if(gammax1.ge.0)then
            xm1low2=(max( 1.d-1,xlep1mass(1)+xlep2mass(1),
     #                    xmw-gammax1*ga1 ))**2
            xm1upp2=(min(xmt-1.d-1,xmw+gammax1*ga1))**2
          else
            xm1low2=(max( 1.d-1,xlep1mass(1)+xlep2mass(1),
     #                    xm1low ))**2
            xm1upp2=(min(xmt-1.d-1,xm1upp))**2
          endif
          if(jwidth.eq.0)then
            xm1upp2=min((xmt-xmass(5))**2-1.d-1,xm1upp2)
          else
            xm1upp2=min((sqrt(ym1upp2)-xmass(5))**2-1.d-1,xm1upp2)
          endif
          if(xm1low2.gt.xm1upp2)then
            write(*,*)'Error in pair mass range #1'
            write(*,*)xm1low2,xm1upp2
            stop
          endif
c Parameters for the Breit Wigner
          bw1mdpl=xm1upp2-xmw2
          bw1mdmn=xmw2-xm1low2
          bw1fmpl=atan(bw1mdpl/(xmw*ga1))
          bw1fmmn=atan(bw1mdmn/(xmw*ga1))
          bw1delf=(bw1fmpl+bw1fmmn)/pi
c
          if(gammax2.ge.0)then
            xm2low2=(max( 1.d-1,xlep1mass(2)+xlep2mass(2),
     #                    xmw-gammax2*ga2 ))**2
            xm2upp2=(min(xmt-1.d-1,xmw+gammax2*ga2))**2
          else
            xm2low2=(max( 1.d-1,xlep1mass(2)+xlep2mass(2),
     #                    xm2low ))**2
            xm2upp2=(min(xmt-1.d-1,xm2upp))**2
          endif
          if(jwidth.eq.0)then
            xm2upp2=min((xmt-xmass(5))**2-1.d-1,xm2upp2)
          else
            xm2upp2=min((sqrt(ym2upp2)-xmass(5))**2-1.d-1,xm2upp2)
          endif
          if(xm2low2.gt.xm2upp2)then
            write(*,*)'Error in pair mass range #2'
            write(*,*)xm2low2,xm2upp2
            stop
          endif
c Parameters for the Breit Wigner function
          bw2mdpl=xm2upp2-xmw2
          bw2mdmn=xmw2-xm2low2
          bw2fmpl=atan(bw2mdpl/(xmw*ga2))
          bw2fmmn=atan(bw2mdmn/(xmw*ga2))
          bw2delf=(bw2fmpl+bw2fmmn)/pi
        else
          write(*,*)'Both mass ranges must be non-zero'
          stop
        endif
c Initialize other parameters
        call setpar()
      endif
c
      do jloop=1,jecm
c main loop (over energies and scale factors); jecm>1, append
c loop number at prefix
         prefn = pref
         if(jecm.gt.1) call fk88strnum(prefn,jloop)
         prefnev = prefev
         if(jecm.gt.1) call fk88strnum(prefnev,jloop)
         ecm = ecmlst(jloop)
         sh = ecm**2
         xm = xmlst(jloop)
         xm2 = xm**2
         if(nlf.lt.0) then
c number of light flavours on the basis of the mass
           if(xm.lt.3) then
             nl = 3
           elseif(xm.lt.7) then
             nl = 4
           else
             nl = 5
           endif
         else
           nl=nlf
         endif
c nl used in alpha_S computation -- allow nl+1 only for charm or bottom
         nlas=nl
         if(nlfp1sch.ge.0)then
           if(xm.lt.7.and.nlfp1sch.eq.1)nlas=nl+1
           if(xm.gt.7.and.nlfp1sch.eq.1)then
             write(*,*)
     #    'nl+1 scheme should not be used for this mass value'
             stop
           endif
         else
           if(xm.lt.7)then
             nlas=nl+1
             nlfp1sch=1
           else
             nlfp1sch=0
           endif
         endif
c check inputs
         if(nl.le.4)then
           if(xmass(nl+1).ne.xm)then
             write(*,*)'Two different values have been assigned to M_Q'
             stop
           endif
         endif
c Herwig code, and consistency check
         itmpqq=-1701-nl
         if(itmpqq.ne.itmpvv.and.(-10000+itmpqq).ne.itmpvv)then
           write(*,*)'Error in process codes:',itmpqq,itmpvv,nl
           stop
         endif
c Heavy quarks identities according to MC particle numbering scheme;
c if tops decay, ip4 and ip5 are identities of decay products, set
c in setpar()
         if(idec.eq.1)then
           ip4 = nl+1
           ip5 = -ip4
         endif
c
         xfh1 = fh1lst(jloop)
         xfh2 = fh2lst(jloop)
         xren = renlst(jloop)
c- common block values for scale factor
         xren2 = xren**2
         xf2h1 = xfh1**2
         xf2h2 = xfh2**2
c- common block values for scale factor (MC terms)
         xren2mc = renmclst(jloop)**2
         xf2h1mc = fh1mclst(jloop)**2
         xf2h2mc = fh2mclst(jloop)**2
c tau generated according to a flat distribution in (1/tau)**nsamp
         nsamp = 1
c
         ndim=6
         nwild=5
         dtot=0.d0
c double differential
         if(iseld.eq.1)then
           xtotal=0.d0
           ytotal=0.d0
           xares=0.d0
           yares=0.d0
           xbres=0.d0
           ybres=0.d0
           mx_of_evta=0
           mx_of_evtb=0
           fname=prefn
c
           ifuntype=1
           call fk88strcat(fname,'_a',fnamea)
           call run_bases(sig5a,fnamea,ndim,nwild,ncl3,it1,it2,
     #       ac1,ac2,av3a,d3a,av3nega,d3nega,ctime,itd1,itd2,iseed0,
     #       ibswrite,ibscall)
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
           call fk88strcat(fname,'_b',fnameb)
           call run_bases(sig5b,fnameb,ndim,nwild,ncl3,it1,it2,
     #       ac1,ac2,av3b,d3b,av3negb,d3negb,ctime,itd1,itd2,iseed0,
     #       ibswrite,ibscall)
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
     #          form='formatted',status='unknown')
           write(21,240)xares
           write(21,240)xbres
           write(21,240)yares
           write(21,240)ybres
           close(21)
 240       format(1x,d14.8)
         endif
c Sanity check
         if(isubttype.eq.1.and.(betfac.ne.1.d0.or.delta.ne.1.d0))then
           write(*,*)'Fatal error: betfac, delta=',betfac,delta
           stop
         endif
         if(iseld.eq.0)then
c Read integrals from disk only if the integration step has been skipped
           call fk88strcat(prefn,'.integrals',fname)
           open(unit=21,file=fname,
     #          form='formatted',status='old')
           read(21,240)xares
           read(21,240)xbres
           read(21,240)yares
           read(21,240)ybres
           close(21)
         endif
c Generates events when evgen=.true.; if evgen=.false., maxevt=100000 in
c order to estimate the number of negative weights
         if(maxevt.ne.0)then
           ntotal=0
           xtotal=0.d0
           ytotal=0.d0
           xtotal=xtotal+xares+xbres
           ytotal=ytotal+yares+ybres
           avtot=ytotal
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
     #               (xares+yares)
           evprcfrac=evprcfrac/(1+evprcfrac)
           evfrac=evfrac+evprcfrac*mx_of_evta
           write(*,*)'Events[a]: w<0/all:',evprcfrac
           evprcfrac=(xbres-ybres)/
     #               (xbres+ybres)
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
     #          form='formatted',status='unknown')
           write(22,250)mx_of_evta
           close(22)
           call fk88strcat(fname,'_b.events',fname1)
           open(unit=22,file=fname1,
     #          form='formatted',status='unknown')
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
     #          form='formatted',status='old',access='append')
           call run_spring(sig5a,fnamea,mx_of_evta,maxtrials,
     #                     nevts,ntrls,ndim,nwild,iseed)
           close(22)
           if(iverbose.eq.1)then
             write(*,*)'   '
             write(*,*)'Events[a]'
             write(*,*)'Trials:',ntrls
             write(*,*)'Events generated:',nevts
             write(*,*)'Unlike sign events(1):',iwrong
             write(*,*)'Unlike sign events(2):',iwrong1
             write(*,*)'Unlike sign(1)/all events:',
     #                 iwrong/dfloat(nevts)
             write(*,*)'Unlike sign(2)/all events:',
     #                 iwrong1/dfloat(nevts)
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
     #                     nqeventsuw/dfloat(nqcntuws)
               endif
             endif
             write(*,*)'   '
             write(*,*)'Average momentum shifts due to masses'
             do i=1,4
               if(idec.eq.0)then
                 write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(11*nevts)
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
     #          form='formatted',status='old',access='append')
           call run_spring(sig5b,fnameb,mx_of_evtb,maxtrials,
     #                     nevts,ntrls,ndim,nwild,iseed)
           close(22)
           if(iverbose.eq.1)then
             write(*,*)'   '
             write(*,*)'Events[b]'
             write(*,*)'Trials:',ntrls
             write(*,*)'Events generated:',nevts
             write(*,*)'Unlike sign events(1):',iwrong
             write(*,*)'Unlike sign events(2):',iwrong1
             write(*,*)'Unlike sign(1)/all events:',
     #                 iwrong/dfloat(nevts)
             write(*,*)'Unlike sign(2)/all events:',
     #                 iwrong1/dfloat(nevts)
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
     #                     nqeventsuw/dfloat(nqcntuws)
               endif
             endif
             write(*,*)'   '
             write(*,*)'Average momentum shifts due to masses'
             do i=1,4
               if(idec.eq.0)then
                 write(*,*)'  ',i,': ',xmomshifts(i)/dfloat(10*nevts)
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
     #          status='unknown')
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
     #       ecmlst(jloop),renlst(jloop),fh1lst(jloop),
     #       renmclst(jloop),fh1mclst(jloop),
     #       ': CM energy, muR/mu0[NLO], muF/mu0[NLO], '//
     #       'muR/mu0[MC], muF/mu0[MC]'
           write(ioutput,802)abs(itmpvv),': 1705/1706=b/t'
           if(itmpvv.eq.-1705.or.itmpvv.eq.-11705)then
             write(ioutput,814)xmlst(jloop),': M_b'
           else
             write(ioutput,803)xmlst(jloop),twidth,
     #                         ': M_top, Gamma_top'
             write(ioutput,803)xmw,gaw,': M_W, Gamma_W'
             write(ioutput,810)il1hw,il2hw,': IL1, IL2 (0..7)'
           endif
           write(ioutput,804)xmass(1),xmass(2),
     #                       xmass(3),xmass(4),
     #                       xmass(5),xmass(21),
     #                       ': quark and gluon masses'
           write(ioutput,805)part1,part2,': colliding particles'
           write(ioutput,806)gname(1:8),idpdfset,
     #       ': PDF group and id number'
           write(ioutput,807)xlam,scheme,': Lambda_5, scheme'
           write(ioutput,811)'H++',': Herwig++ (v3.3 and higher)'
           write(ioutput,250)maxevt
           if(ievffmt.eq.1)then
c Close LHEF header
             if(idec.eq.0)then
               xsecpp=avtot*brrtop1*brrtop2
               xerrpp=dtot*brrtop1*brrtop2
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
             call write_lhef_init(ioutput,itmpvv,ecmlst(jloop),
     #                            xsecpp,xerrpp,part1,part2)
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
 111       continue
         endif
         if(idec.eq.0)then
           write(*,*) '   '
           write(*,*)'Branching ratios used in the computation:'
           write(*,*)' BR(t -> e nu [b+d+s])= ',xbrrtoplep
           write(*,*)' BR(t -> [udbar+usbar+ubbar] [b+d+s])= ',
     #               xbrrtophad
           write(*,*) '   '
           write(*,*)'Normalization factor due to decays:',
     #               brrtop1*brrtop2
         endif 
         write(*,*) '   '
         write(*,*) 'Total for fully inclusive'
         write(*,200)ih1,ih2,ndns1,ndns2,nl,xlam
         write(*,'(a)')
     #  ' ecm           mass      f1   f2   r    tot        err'
         write(*,300)
     #      abs(ecmlst(jloop)),xmlst(jloop),
     #      fh1lst(jloop),fh2lst(jloop),renlst(jloop),avtot,dtot
c end of the main loop
      enddo
 200  format(' had1=',i2,'  had2=',i2,'  strf1=',i6,'  strf2=',i6,
     #  '  nl=',i2,'  lambda5=',d10.4)
 300  format((1x,1pd9.3),4x,(1x,1pd9.3),3(1x,0pf4.2),
     # 2(1x,0pd10.4,1x,1pd6.0))
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


      subroutine strfun(x1,x2,sf)
c Return parton densities through the matrix
c  sf(idr,jproc,itype), with the following conventions:
c  idr=1 -> direct events
c  idr=2 -> charge-conjugated events
c  idr=3 -> reflected events
c  idr=4 -> charge-conjugated and reflected events
c  jproc=1,2,3 -> gg, qqbar, qg processes respectively
c  itype -> identifies the individual contribution to a given jproc
c
      implicit none
      real*4 fh1x1(-5:5),fh2x2(-5:5),fh1x2(-5:5),fh2x1(-5:5),
     #  smuf2h1,smuf2h2
      real*8 pi,x1,x2,sf(4,3,5)
      integer ih1,ih2,ndns1,ndns2,i,jproc,itype
      parameter(pi=3.14159265358979312D0)
      include 'hvqcblks.h'
      common/strfun0/ih1,ih2,ndns1,ndns2
      integer ipdfscale
      common/cipdfscale/ipdfscale
c ipdfscale=1 --> use NLO factorization scale
c ipdfscale=2 --> use MC factorization scale
c
      do i=1,4
        do jproc=1,3
          do itype=1,5
            sf(i,jproc,itype)=0.d0
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
      call mlmpdf(ndns1,ih1,smuf2h1,sngl(x2),fh1x2,5)
      call mlmpdf(ndns2,ih2,smuf2h2,sngl(x1),fh2x1,5)
c
      sf(1,1,1) = dble(fh1x1(0)*fh2x2(0))/4.d0
      sf(2,1,1) = sf(1,1,1)
      sf(3,1,1) = dble(fh2x1(0)*fh1x2(0))/4.d0
      sf(4,1,1) = sf(3,1,1)
c
      do itype=1,nl
        sf(1,2,itype) = dble(fh1x1( itype) * fh2x2(-itype))/2
        sf(2,2,itype) = dble(fh1x1(-itype) * fh2x2( itype))/2
        sf(3,2,itype) = dble(fh2x1( itype) * fh1x2(-itype))/2
        sf(4,2,itype) = dble(fh2x1(-itype) * fh1x2( itype))/2
      enddo
c
      do itype=1,nl
        sf(1,3,itype) = dble(fh1x1( itype) * fh2x2( 0))
        sf(2,3,itype) = dble(fh1x1(-itype) * fh2x2( 0))
        sf(3,3,itype) = dble(fh2x1( itype) * fh1x2( 0))
        sf(4,3,itype) = dble(fh2x1(-itype) * fh1x2( 0))
      enddo
c
      return
      end


      function sig5a(xx)
c Integrand function for H events
      implicit none
      real * 8 sig5a,xx
      real * 8 pi,tiny
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-5)
      dimension xx(6)
      include 'hvqcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      integer nsamp
      common/samp/nsamp
      integer iprespl
      common/ciprespl/iprespl
      real * 8 xjac,roh,zzz,ttt,th,th2,x,y,csi,cth1,rx,
     # ximax0,ximin0,ymax,ymin,s,ro,rox,rohx,tmp,tot5a,
     # taumax,xxa1,xxa2,xxc,xxymax,xxymin
c
c xx(1) --> tau, xx(2) --> ycm, xx(3) --> x, xx(4) --> y, xx(5) --> cth1,
c xx(6) --> th2
c
      xjac = 1
      roh   = 4*xm2/sh
c
c To improve convergence in the soft regions
      zzz = tiny+(1-tiny)*xx(3)**2
      xjac = xjac * xx(3) * 2
      x = 1 - zzz*(1-roh)
      xjac = xjac * (1-roh)
c
c To improve convergence in the collinear regions
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
c
      y    = cos(th)
      xjac = xjac * sin(th)
c
c Generation of tau and ycm values and computation of the integration
c limits:
c
      csi = sqrt((1-(1-x)*(1+y)/2)/(1-(1-x)*(1-y)/2))
      rx = sqrt(x)
      rohx = roh/x
      taumax = 1/x
      ximax0 = rohx**(-nsamp)
      ximin0 = taumax**(-nsamp)
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
      if(iprespl.eq.0)then
        ymax= -log(tau)/2 + log(1/(csi*rx))
        ymin=  log(tau)/2 - log(csi/rx)
      else
        xxa1 = (1+x-y*(1-x))/2.d0
        xxa2 = (1+x+y*(1-x))/2.d0
        xxc = (1-x*tau)/sqrt(tau)
        xxymax = (xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        xxymin = (-xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        ymax = max(log(xxymax),-log(tau)/2.d0)
        ymin = min(log(xxymin),log(tau)/2.d0)
      endif
      ycm = ymin + xx(2)*(ymax-ymin)
      xjac= xjac * (ymax-ymin)
c
      s = sh * tau
      ro = roh/tau
c
c Change variables from xx(5) to cth1, xjac--> xjac * d cth1/d xx(5)
c
      rox  = ro/x
      call zzchvar(xx(5),cth1,xjac,rox)
c
      th2 = xx(6) * pi
      xjac = xjac * pi
c
      sig5a = tot5a(s,x,y,cth1,th2,xjac)
      return
      end


      function tot5a(s,xx,xy,xcth1,xth2,xjac)
c Implements standard subtraction
      implicit none
      character*2 str
      parameter (str='p1')
      real*8 tot5a,s,xx,xy,xcth1,xth2,xjac
      real*8 pi,pi2,hc2
      parameter (pi=3.14159265358979312D0)
      parameter (pi2=pi*pi)
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
c TeV to picobarn conversion factor: sigma (pb) = hc2 * sigma (TeV^-2)
c sigma (pb) = 10^6 * sigma (mub)
      parameter (hc2=3.8937966d2)
      integer ione
      parameter (ione=1)
      character*2 xproc(3)
      common/cxproc/xproc
      real*8 betfac,delta,deltas,deltac
      common/betfac/betfac,delta
      common/pmerge/deltas,deltac
      include 'hvqcblks.h'
      real*8 ycm,tau
      common/x1x2/ycm,tau
      real*8 sf(4,3,5)
      integer ipdfscale
      common/cipdfscale/ipdfscale
      integer idec
      common/cidec/idec
      real*8 bsfsgn
      common/cbssgn/bsfsgn
      real*8 bsewgt
      common/cbswgt/bsewgt
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 ps,px,py,pcth1,pcth2
      common/cpsave/ps,px,py,pcth1,pcth2
      real*8 vv(4,3,5),vvs(4,3,5)
      common/cvv/vv
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
      real*8 x,y,cth1,th2,cth2,sx,rox,bx,ro,beta,
     # x1,x2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,zg2_nlo,
     # xnorm,f,xint,www,x1t,x2t,xtmp,ytmp,
     # zgmu2_nlo,zgmu6_mc,zg6_mc,xphsp,xphspcm,xphspcp,xfact,
     # x1soft,x2soft,x1x2j,x1x2jac,zhwfct,xsum,dummy,fpp,
     # gfactsf,gfactcl,xnormmc,xphspmc,betamc
      real*8 xmcxsec(1:4,1:3),xmce0sq(1:4,1:3),xmcz(1:4,1:3)
      real*8 xqrksc(1:3,1:4,1:3),xqbrsc(1:3,1:4,1:3)
      real*8 gfsf(3),gfcl(3)
      integer i,itype,jproc,loproc,maproc,ileg,ie0sq,i2b,
     # iret,itoosoftkin
      common/cjproc/jproc
      common/cwchproc/loproc,maproc
      logical flxsec(1:4,1:3),flagmc,fx1x2,gfflag
c
      x=xx
      y=xy
      cth1=xcth1
      th2=xth2
      cth2=cos(th2)
c
      sx=s*x
      ro=4*xm2/s
      beta=sqrt(1-ro)
      rox=4*xm2/sx
      bx=sqrt(1-rox)
c
      x1=sqrt(tau)*exp(ycm)
      x2=tau/x1
c
      do jproc=1,3
        do i=1,4
          do itype=1,nl
            vv(i,jproc,itype)=0.d0
            vvs(i,jproc,itype)=0.d0
          enddo
        enddo
      enddo
      xnorm=xjac / (s * 64*pi2 * 16*pi2)
      xphsp=bx/(1-x)*( 1/(1-y) + 1/(1+y) ) 
      xphspcp=bx/((1-x)*(1-y))
      xphspcm=bx/((1-x)*(1+y))
      xnormmc=xjac / (64*pi2 * 16*pi2)
      xphspmc=1/(1-x)*( 1/(1-y) + 1/(1+y) ) 
c
c Event
c
      if(x1.lt.1.and.x2.lt.1)then
        ipdfscale=1
        call invar(xm2,s,x,y,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        zg2_nlo=zgmu2_nlo()
        do jproc=loproc,maproc
          prc=xproc(jproc)
          f=fpp(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
          www=zg2_nlo**3*xnorm*xphsp*f
          call strfun(x1,x2,sf)
          do i=1,4
            do itype=1,nl
              vv(i,jproc,itype)=sf(i,jproc,itype)*www
            enddo
          enddo
        enddo
      endif
c
c MC subt term: pure MC
c
      if(x1.lt.1.and.x2.lt.1)then
        zg6_mc=zgmu6_mc(zg2_nlo)
      else
        zg6_mc=0.d0
      endif
      ipdfscale=2
      gfflag=.false.
      do jproc=loproc,maproc
        call xmcsubtpp(x1,x2,xm2,s,x,y,cth1,cth2,
     #    xmcxsec,xmce0sq,xmcz,xqrksc,xqbrsc,flxsec,flagmc,
     #    gfactsf,gfactcl)
        gfsf(jproc)=gfactsf
        gfcl(jproc)=gfactcl
        gfflag=gfflag.or.(gfactsf.lt.1.d0)
        if(flagmc)then
          do ileg=1,4
            do ie0sq=1,3
              if(flxsec(ileg,ie0sq))then
                if(ileg.eq.1)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)/zhwfct
                  x2t=x2soft(x1,x2,x,y)
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1
                elseif(ileg.eq.2)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)
                  x2t=x2soft(x1,x2,x,y)/zhwfct
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1
                else
                  zhwfct=1.d0
                  x1t=x1
                  x2t=x2
                  betamc=beta
                  fx1x2=x1t.lt.1.and.x2t.lt.1
                endif
                if(fx1x2)then
                  x1x2j=x1x2jac(x1,x2,x,y,ileg)/zhwfct
                  www=-zg6_mc*xnormmc*xphspmc*x1x2j*betamc*
     #                 xmcxsec(ileg,ie0sq)
                  call strfun(x1t,x2t,sf)
                  do i=1,4
                    do itype=1,nl
                      vv(i,jproc,itype)=vv(i,jproc,itype)+
     #                  sf(i,jproc,itype)*www
                    enddo
                  enddo
                endif
              endif
            enddo
          enddo
        endif
      enddo
c
c MC subt term: collinear ME
c
      if(gfflag)then
        if(y.gt.0.d0)then
          ytmp=1.d0
          x1t=x1soft(x1,x2,x,y)/x
          x2t=x2soft(x1,x2,x,y)
          xfact=xnorm*xphspcp
        else
          ytmp=-1.d0
          x1t=x1soft(x1,x2,x,y)
          x2t=x2soft(x1,x2,x,y)/x
          xfact=xnorm*xphspcm
        endif
        if(x1t.lt.1.and.x2t.lt.1.and.x1.lt.1.and.x2.lt.1)then
          ipdfscale=1
          x1x2j=x1x2jac(x1,x2,x,y,ione)/x
          call invar(xm2,s,x,ytmp,cth1,cth2,str,
     #               tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          zg2_nlo=zgmu2_nlo()
          do jproc=loproc,maproc
            if(gfsf(jproc).lt.1.d0)then
              prc=xproc(jproc)
              f=fpp(s,x,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              www=-zg2_nlo**3*xfact*x1x2j*f*(1-gfsf(jproc))
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            endif
          enddo
        endif
      endif
c
c MC subt term: soft ME
c
      if(gfflag)then
        xtmp=1.d0
        x1t=x1soft(x1,x2,x,y)
        x2t=x2soft(x1,x2,x,y)
        if(x1t.lt.1.and.x2t.lt.1.and.x1.lt.1.and.x2.lt.1)then
          ipdfscale=1
          x1x2j=x1x2jac(x1,x2,x,y,ione)
          call invar(xm2,sx,xtmp,y,cth1,cth2,str,
     #               tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          zg2_nlo=zgmu2_nlo()
          do jproc=loproc,maproc
            if(gfsf(jproc).lt.1.d0)then
              prc=xproc(jproc)
              f=fpp(sx,xtmp,y,xm2,q1q,q2q,w1h,w2h,cth2)
              www=-zg2_nlo**3*(xnorm/x)*xphsp*x1x2j*
     #            f*(1-gfsf(jproc))
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            endif
          enddo
        endif
      endif
c
c MC subt term: soft-collinear ME
c
      if(gfflag)then
        if(y.gt.0.d0)then
          ytmp=1.d0
          xfact=(xnorm/x)*xphspcp
        else
          ytmp=-1.d0
          xfact=(xnorm/x)*xphspcm
        endif
        xtmp=1.d0
        x1t=x1soft(x1,x2,x,y)
        x2t=x2soft(x1,x2,x,y)
        if(x1t.lt.1.and.x2t.lt.1.and.x1.lt.1.and.x2.lt.1)then
          ipdfscale=1
          x1x2j=x1x2jac(x1,x2,x,y,ione)
          call invar(xm2,sx,xtmp,ytmp,cth1,cth2,str,
     #               tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          zg2_nlo=zgmu2_nlo()
          do jproc=loproc,maproc
            if(gfsf(jproc).lt.1.d0)then
              prc=xproc(jproc)
              f=fpp(sx,xtmp,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              www=zg2_nlo**3*xfact*x1x2j*f*(1-gfsf(jproc))
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            endif
          enddo
        endif
      endif
c
      call checkvv(xsum,dummy,iret)
      if(iret.eq.1)then
        do jproc=loproc,maproc
          do i=1,4
            do itype=1,nl
              vvs(i,jproc,itype)=vv(i,jproc,itype)
            enddo
          enddo
        enddo
        call invar(xm2,s,x,y,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        if(idec.eq.0)then
          ps=s
          px=x
          py=y
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
          call invar(xm2,sx,xtmp,ytmp,cth1,cth2,str,
     #               tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          if(idec.eq.0)then
            ps=sx
            px=xtmp
            py=ytmp
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
      tot5a=abs(xint)
c
      return
      end


      function sig5b(xx)
c Integrand function for S events
      implicit none
      real * 8 sig5b,xx
      real * 8 pi,tiny
      parameter (pi=3.14159265358979312D0)
      parameter (tiny=1.d-5)
      dimension xx(6)
      include 'hvqcblks.h'
      real * 8 ycm,tau
      common/x1x2/ycm,tau
      integer nsamp
      common/samp/nsamp
      integer iprespl
      common/ciprespl/iprespl
      real * 8 xjac,roh,zzz,ttt,th,th2,x,y,csi,cth1,rx,
     # ximax0,ximin0,ymax,ymin,s,ro,rox,rohx,tmp,tot5b,
     # taumax,xxa1,xxa2,xxc,xxymax,xxymin
c
c xx(1) --> tau, xx(2) --> ycm, xx(3) --> x, xx(4) --> y, xx(5) --> cth1,
c xx(6) --> th2
c
      xjac = 1
      roh   = 4*xm2/sh
c
c To improve convergence in the soft regions
      zzz = tiny+(1-tiny)*xx(3)**2
      xjac = xjac * xx(3) * 2
      x = 1 - zzz*(1-roh)
      xjac = xjac * (1-roh)
c
c To improve convergence in the collinear regions
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
c
      y    = cos(th)
      xjac = xjac * sin(th)
c
c Generation of tau and ycm values and computation of the integration
c limits:
c
      csi = sqrt((1-(1-x)*(1+y)/2)/(1-(1-x)*(1-y)/2))
      rx = sqrt(x)
      rohx = roh/x
      taumax = 1/x**2
      ximax0 = rohx**(-nsamp)
      ximin0 = taumax**(-nsamp)
      tmp  = ximin0 + xx(1)*(ximax0-ximin0)
      tau = tmp**(-1/dfloat(nsamp))
      xjac= xjac/nsamp*tau**(nsamp+1)*(ximax0-ximin0)
      if(iprespl.eq.0)then
        ymax= -log(tau)/2 + log(1/(csi*rx))
        ymin=  log(tau)/2 - log(csi/rx)
      else
        xxa1 = (1+x-y*(1-x))/2.d0
        xxa2 = (1+x+y*(1-x))/2.d0
        xxc = (1-x*tau)/sqrt(tau)
        xxymax = (xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        xxymin = (-xxc+sqrt(xxc**2+4*xxa1*xxa2))/(2*xxa1)
        ymax = max(log(xxymax),-log(tau)/2.d0)
        ymin = min(log(xxymin),log(tau)/2.d0)
      endif
      ycm = ymin + xx(2)*(ymax-ymin)
      xjac= xjac * (ymax-ymin)
c
      s = sh * tau
      ro = roh/tau
c
c Change variables from xx(5) to cth1, xjac--> xjac * d cth1/d xx(5)
c
      rox  = ro/x
      call zzchvar(xx(5),cth1,xjac,rox)
c
      th2 = xx(6) * pi
      xjac = xjac * pi
c
      sig5b = tot5b(s,x,y,cth1,th2,xjac)
      return
      end


      function tot5b(s,xx,xy,xcth1,xth2,xjac)
      implicit none
      real * 8 tot5b,tot5bs,tot5bz,s,xx,xy,xcth1,xth2,
     #  x,y,cth1,th2,xjac,tmp
      integer isubttype
      common/cisubttype/isubttype
c
      x      = xx
      y      = xy
      cth1   = xcth1
      th2    = xth2
      if(isubttype.eq.0)then
        tmp=tot5bs(s,x,y,cth1,th2,xjac)
      elseif(isubttype.eq.1)then
        tmp=tot5bz(s,x,y,cth1,th2,xjac)
      else
        write(*,*)'Fatal error in tot5b:',isubttype
        stop
      endif
      tot5b=tmp
      return
      end


      function tot5bs(s,xx,xy,xcth1,xth2,xjac)
c Implements standard subtraction
      implicit none
      character*2 str
      parameter (str='p1')
      real*8 tot5bs,s,xx,xy,xcth1,xth2,xjac
      real*8 pi,pi2,hc2
      parameter (pi=3.14159265358979312D0)
      parameter (pi2=pi*pi)
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
c TeV to picobarn conversion factor: sigma (pb) = hc2 * sigma (TeV^-2)
c sigma (pb) = 10^6 * sigma (mub)
      parameter (hc2=3.8937966d2)
      integer ione
      parameter (ione=1)
      character*2 xproc(3)
      common/cxproc/xproc
      real*8 betfac,delta,deltas,deltac
      common/betfac/betfac,delta
      common/pmerge/deltas,deltac
      include 'hvqcblks.h'
      real*8 ycm,tau
      common/x1x2/ycm,tau
      real*8 sf(4,3,5)
      integer ipdfscale
      common/cipdfscale/ipdfscale
      integer idec
      common/cidec/idec
      real*8 bsfsgn
      common/cbssgn/bsfsgn
      real*8 bsewgt
      common/cbswgt/bsewgt
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 ps,px,py,pcth1,pcth2
      common/cpsave/ps,px,py,pcth1,pcth2
      real*8 vv(4,3,5),vvs(4,3,5)
      common/cvv/vv
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
      real*8 x,y,cth1,th2,cth2,sx,rox,bx,btildex,rotildx,
     # x1,x2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,zg2_nlo,
     # xnorm,f,xint,www,xlgomx,x1t,x2t,xtmp,ytmp,
     # xlmude,xnormc,xcplus,xcminus,xintcm,xintffs,
     # xnormsv,xintcps,xintcms,xintcp,xnormb,bbb,x1soft,
     # x2soft,x1x2j,x1x2jac,hvqborn,ppsv,fpp,ppcolp,ppcoll,
     # xphsp,xphspcp,xphspcm,zgmu6_mc,zg6_mc,gfactsf,gfactcl,zhwfct,
     # dfact,gcpfact,zgmu2_nlo,gcmfact,gsffact,dfact1,
     # xsum,dummy,xnormmc,xphspmc,ro,beta,betamc,xnlfscheme
      real*8 xcs(3),xsv(3),xints(3),xborn(3)
      real*8 xmcxsec(1:4,1:3),xmce0sq(1:4,1:3),xmcz(1:4,1:3)
      real*8 xqrksc(1:3,1:4,1:3),xqbrsc(1:3,1:4,1:3)
      real*8 gfsf(3),gfcl(3)
      integer i,itype,ileg,ie0sq,iret,iproc,iproclo,iprocma
      integer jproc
      common/cjproc/jproc
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      logical flxsec(1:4,1:3),flagmc,fx1x2,gfflag
c
      x=xx
      y=xy
      cth1=xcth1
      th2=xth2
      cth2=cos(th2)
c
      sx=s*x
      ro=4*xm2/s
      beta=sqrt(1-ro)
      rox=4*xm2/sx
      bx=sqrt(1-rox)
      btildex=bx*betfac
      rotildx=1-btildex**2
c
      x1=sqrt(tau)*exp(ycm)
      x2=tau/x1
c
      xlgomx=log(1-x)
      do jproc=1,3
        do i=1,4
          do itype=1,nl
            vv(i,jproc,itype)=0.d0
            vvs(i,jproc,itype)=0.d0
          enddo
        enddo
      enddo
      xnorm=xjac / (s * 64*pi2 * 16*pi2)
      xphsp=bx/(1-x)*( 1/(1-y) + 1/(1+y) ) 
      xphspcp=bx/((1-x)*(1-y))
      xphspcm=bx/((1-x)*(1+y))
      xnormmc=xjac / (64*pi2 * 16*pi2)
      xphspmc=1/(1-x)*( 1/(1-y) + 1/(1+y) ) 
c
c MC subt term: pure MC
c
      if(x1.lt.1.and.x2.lt.1.and.s.gt.4*xm2)then
        call invar(xm2,s,x,y,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        zg2_nlo=zgmu2_nlo()
        zg6_mc=zgmu6_mc(zg2_nlo)
      else
        zg6_mc=0.d0
      endif
      ipdfscale=2
      gfflag=.false.
      do jproc=loproc,maproc
        call xmcsubtpp(x1,x2,xm2,s,x,y,cth1,cth2,
     #    xmcxsec,xmce0sq,xmcz,xqrksc,xqbrsc,flxsec,flagmc,
     #    gfactsf,gfactcl)
        gfsf(jproc)=gfactsf
        gfcl(jproc)=gfactcl
        gfflag=gfflag.or.(gfactsf.lt.1.d0)
        if(flagmc)then
          do ileg=1,4
            do ie0sq=1,3
              if(flxsec(ileg,ie0sq))then
                if(ileg.eq.1)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)/zhwfct
                  x2t=x2soft(x1,x2,x,y)
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1.and.
     #                  (x*tau).lt.1
                elseif(ileg.eq.2)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)
                  x2t=x2soft(x1,x2,x,y)/zhwfct
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1.and.
     #                  (x*tau).lt.1
                else
                  zhwfct=1.d0
                  x1t=x1soft(x1,x2,x,y)
                  x2t=x2soft(x1,x2,x,y)
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  (x**2*tau).lt.1
                endif
                if(fx1x2)then
                  x1x2j=x1x2jac(x1,x2,x,y,ione)/zhwfct
                  www=zg6_mc*xnormmc*xphspmc*x1x2j*betamc*
     #                xmcxsec(ileg,ie0sq)
                  call strfun(x1t,x2t,sf)
                  do i=1,4
                    do itype=1,nl
                      vv(i,jproc,itype)=vv(i,jproc,itype)+
     #                  sf(i,jproc,itype)*www
                    enddo
                  enddo
                endif
              endif
            enddo
          enddo
        endif
      enddo
c
c Counter-event (x,y)=(x,1) and MC subt term, collinear ME part
c
      if(y.gt.1-delta.or.gfflag) then
         dfact=0.d0
         if(y.gt.1-delta)dfact=1.d0
         ytmp=1.d0
         x1t=x1soft(x1,x2,x,y)/x
         x2t=x2soft(x1,x2,x,y)
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)/x
            call invar(xm2,s,x,ytmp,cth1,cth2,'p1',
     #                 tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gcpfact=0.d0
              if(y.gt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #          gcpfact=1-gfsf(jproc)
              prc=xproc(jproc)
              xintcp=0.d0
              xcplus=0.d0
              f=fpp(s,x,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              xintcp=zg2_nlo**3*xnorm*xphspcp*x1x2j*f*(gcpfact-dfact)
              if(dfact.eq.1.d0)then
c Adding the collinear contribution
                xlmude=log(s/xmuf2h1)+log(delta/2)
                xnormc=zg2_nlo**3*x1x2j*xjac/(delta*16*pi2)
                xcplus=xnormc*bx/(1-x)*(
     #                        ppcolp(ytmp,s,q2q,x,xm2,xlmude)
     #                +xlgomx*ppcoll(ytmp,s,q2q,x,xm2) )
              endif
              www=xintcp + xcplus
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Counter-event (x,y)=(x,-1) and MC subt term, collinear ME part
c
      if(y.lt.-1+delta.or.gfflag) then
         dfact=0.d0
         if(y.lt.-1+delta)dfact=1.d0
         ytmp=-1.d0
         x1t=x1soft(x1,x2,x,y)
         x2t=x2soft(x1,x2,x,y)/x
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)/x
            call invar (xm2,s,x,ytmp,cth1,cth2,'p1',
     #            tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gcmfact=0.d0
              if(y.lt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #          gcmfact=1-gfsf(jproc)
              prc=xproc(jproc)
              xintcm=0.d0
              xcminus=0.d0
              f=fpp(s,x,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              xintcm=zg2_nlo**3*xnorm*xphspcm*x1x2j*f*(gcmfact-dfact)
              if(dfact.eq.1.d0)then
c Adding the collinear contribution
                xlmude=log(s/xmuf2h2)+log(delta/2)
                xnormc=zg2_nlo**3*x1x2j*xjac/(delta*16*pi2)
                xcminus=xnormc*bx/(1-x)*(
     #                        ppcolp(ytmp,s,q1q,x,xm2,xlmude)
     #                +xlgomx*ppcoll(ytmp,s,q1q,x,xm2) )
              endif
              www=xintcm + xcminus
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Soft counter-events, and MC subt term, soft and soft-collinear ME parts
c
      if(x.gt.rotildx.or.gfflag) then
         dfact=0.d0
         if(x.gt.rotildx)dfact=1.d0
         xtmp=1.d0
         x1t=x1soft(x1,x2,x,y)
         x2t=x2soft(x1,x2,x,y)
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)
            call invar(xm2,sx,xtmp,y,cth1,cth2,'p1',
     #           tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gsffact=0.d0
              if(x1.lt.1.and.x2.lt.1)gsffact=1-gfsf(jproc)
              prc=xproc(jproc)
              xints(jproc)=0.d0
              xborn(jproc)=0.d0
              xcs(jproc)=0.d0
              xsv(jproc)=0.d0
              f=fpp(sx,xtmp,y,xm2,q1q,q2q,w1h,w2h,cth2)
              xintffs=zg2_nlo**3*(xnorm/x)*xphsp*x1x2j*f*(gsffact-dfact)
              xints(jproc)=xints(jproc)+xintffs
              if(dfact.eq.1.d0)then
c Adding the soft-virtual contribution
                xnormsv=zg2_nlo**3*x1x2j*xjac/ 
     #                  (32*pi2*16*pi2*(1-rotildx))
                xsv(jproc)=xnormsv*bx*
     #                     ppsv(sx,q1q,xm2,xmur2,xmuf2h1,xmuf2h2)
c Adding the Born term
                xnormb=zg2_nlo**2*x1x2j*xjac/(16*pi*2*pi*(1-rotildx))
                bbb=hvqborn(sx,q1q,xm2,jproc)
                xborn(jproc)=xborn(jproc)+xnormb*bx*bbb*( 1.d0+
     #            xnlfscheme(xm2,xmur2,xmuf2h1,xmuf2h2,zg2_nlo,jproc) )
              endif
            enddo
            if(y.gt.1-delta.or.gfflag) then
              ipdfscale=1
              dfact1=0.d0
              if(y.gt.1-delta.and.dfact.eq.1.d0)dfact1=1.d0
              ytmp=1.d0
              call invar(xm2,sx,xtmp,ytmp,cth1,cth2,'p1',
     #          tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
              zg2_nlo=zgmu2_nlo()
              do jproc=loproc,maproc
                gcpfact=0.d0
                if(y.gt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #            gcpfact=1-gfsf(jproc)
                prc=xproc(jproc)
                f=fpp(sx,xtmp,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
                xintcps=-zg2_nlo**3*(xnorm/x)*xphspcp*x1x2j*f*
     #                  (gcpfact-dfact1)
                xints(jproc)=xints(jproc)+xintcps
                if(dfact1.eq.1.d0)then
c Adding the collinear contribution
                  xlmude=log(sx/xmuf2h1)+log(delta/2)
                  xnormc=zg2_nlo**3*x1x2j*xjac/(delta*16*pi2)
                  xcs(jproc)=xcs(jproc) - xnormc*bx/(1-x)*(
     #                       ppcolp(ytmp,sx,q2q,xtmp,xm2,xlmude)
     #               +xlgomx*ppcoll(ytmp,sx,q2q,xtmp,xm2) )
                endif
              enddo
            endif
            if(y.lt.-1+delta.or.gfflag) then
              ipdfscale=1
              dfact1=0.d0
              if(y.lt.-1+delta.and.dfact.eq.1.d0)dfact1=1.d0
              ytmp=-1.d0
              call invar(xm2,sx,xtmp,ytmp,cth1,cth2,'p1',
     #          tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
              zg2_nlo=zgmu2_nlo()
              do jproc=loproc,maproc
                gcmfact=0.d0
                if(y.lt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #            gcmfact=1-gfsf(jproc)
                prc=xproc(jproc)
                f=fpp(sx,xtmp,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
                xintcms=-zg2_nlo**3*(xnorm/x)*xphspcm*x1x2j*f*
     #                  (gcmfact-dfact1)
                xints(jproc)=xints(jproc) + xintcms
                if(dfact1.eq.1.d0)then
c Adding the collinear contribution
                  xlmude=log(sx/xmuf2h2)+log(delta/2)
                  xnormc=zg2_nlo**3*x1x2j*xjac/(delta*16*pi2)
                  xcs(jproc)=xcs(jproc) - xnormc*bx/(1-x)*(
     #                       ppcolp(ytmp,sx,q1q,xtmp,xm2,xlmude)
     #               +xlgomx*ppcoll(ytmp,sx,q1q,xtmp,xm2) )
                endif
              enddo
            endif
c
            call strfun(x1t,x2t,sf)
            do jproc=loproc,maproc
              www=xints(jproc)+xsv(jproc)+xborn(jproc)+xcs(jproc)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Compute Born ME times luminosity to get flavour assignment; the
c normalization is irrelevant
c
      xtmp=1.d0
      x1t=x1soft(x1,x2,x,y)
      x2t=x2soft(x1,x2,x,y)
      if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
        ipdfscale=1
        call invar(xm2,sx,xtmp,y,cth1,cth2,'p1',
     #      tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        zg2_nlo=zgmu2_nlo()
        call strfun(x1t,x2t,sf)
        if(loproc.eq.3.and.maproc.eq.3)then
          iproclo=1
          iprocma=2
        else
          iproclo=loproc
          iprocma=maproc
        endif
        do iproc=iproclo,iprocma
          prc=xproc(iproc)
          www=hvqborn(sx,q1q,xm2,iproc)
          do i=1,4
            do itype=1,nl
              vvs(i,iproc,itype)=sf(i,iproc,itype)*www
            enddo
          enddo
        enddo
      endif
c
      call checkvv(xsum,dummy,iret)
      if(iret.eq.1)then
        xtmp=1.d0
        ytmp=1.d0
        call invar(xm2,sx,xtmp,ytmp,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        x1t=x1soft(x1,x2,x,y)
        x2t=x2soft(x1,x2,x,y)
        ycm=0.5d0*log(x1t/x2t)
        tau=x*tau
        if(idec.eq.0)then
          ps=sx
          px=xtmp
          py=ytmp
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
      tot5bs=abs(xint)
c
      return
      end


      function tot5bz(s,xx,xy,xcth1,xth2,xjac)
c Implements zeta subtraction
      implicit none
      character*2 str
      parameter (str='p1')
      real*8 tot5bz,s,xx,xy,xcth1,xth2,xjac
      real*8 pi,pi2,hc2
      parameter (pi=3.14159265358979312D0)
      parameter (pi2=pi*pi)
c GeV to microbarn conversion factor: sigma (mub) = hc2 * sigma (GeV^-2)
c TeV to picobarn conversion factor: sigma (pb) = hc2 * sigma (TeV^-2)
c sigma (pb) = 10^6 * sigma (mub)
      parameter (hc2=3.8937966d2)
      integer ione
      parameter (ione=1)
      character*2 xproc(3)
      common/cxproc/xproc
      real*8 betfac,delta,deltas,deltac,etacut
      common/betfac/betfac,delta
      common/pmerge/deltas,deltac
      common/cetacut/etacut
      include 'hvqcblks.h'
      real*8 ycm,tau
      common/x1x2/ycm,tau
      real*8 sf(4,3,5)
      integer ipdfscale
      common/cipdfscale/ipdfscale
      integer idec
      common/cidec/idec
      real*8 bsfsgn
      common/cbssgn/bsfsgn
      real*8 bsewgt
      common/cbswgt/bsewgt
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 ps,px,py,pcth1,pcth2
      common/cpsave/ps,px,py,pcth1,pcth2
      real*8 vv(4,3,5),vvs(4,3,5)
      common/cvv/vv
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
      real*8 x,y,cth1,th2,cth2,sx,rox,bx,btildex,rotildx,
     # x1,x2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,zg2_nlo,
     # xnorm,f,xint,www,xlgomx,x1t,x2t,xtmp,ytmp,
     # xlmude,xnormc,xcplus,xcminus,xintcm,xintffs,
     # xnormsv,xintcps,xintcms,xintcp,xnormb,bbb,x1soft,
     # x2soft,x1x2j,x1x2jac,hvqborn,ppsv,fpp,ppcolp,ppcoll,
     # xphsp,xphspcp,xphspcm,zgmu6_mc,zg6_mc,gfactsf,gfactcl,zhwfct,
     # dfact,gcpfact,zgmu2_nlo,gcmfact,gsffact,dfact1,
     # xsum,dummy,xktrel,bdelta,svn,delppsv,bsub,xnormmc,
     # xphspmc,ro,beta,betamc,xnlfscheme
      real*8 xcs(3),xsv(3),xints(3),xborn(3)
      real*8 xmcxsec(1:4,1:3),xmce0sq(1:4,1:3),xmcz(1:4,1:3)
      real*8 xqrksc(1:3,1:4,1:3),xqbrsc(1:3,1:4,1:3)
      real*8 gfsf(3),gfcl(3)
      integer i,itype,ileg,ie0sq,iret,iproc,iproclo,iprocma
      integer jproc
      common/cjproc/jproc
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      logical flxsec(1:4,1:3),flagmc,fx1x2,gfflag
c
      x=xx
      y=xy
      cth1=xcth1
      th2=xth2
      cth2=cos(th2)
c
      sx=s*x
      ro=4*xm2/s
      beta=sqrt(1-ro)
      rox=4*xm2/sx
      bx=sqrt(1-rox)
      btildex=bx*betfac
      rotildx=1-btildex**2
      xktrel = (1-x)**2*(1-y**2)
c
      x1=sqrt(tau)*exp(ycm)
      x2=tau/x1
c
      xlgomx=log(1-x)
      do jproc=1,3
        do i=1,4
          do itype=1,nl
            vv(i,jproc,itype)=0.d0
            vvs(i,jproc,itype)=0.d0
          enddo
        enddo
      enddo
      xnorm=xjac / (s * 64*pi2 * 16*pi2)
      xphsp=bx/(1-x)*( 1/(1-y) + 1/(1+y) ) 
      xphspcp=bx/((1-x)*(1-y))
      xphspcm=bx/((1-x)*(1+y))
      xnormmc=xjac / (64*pi2 * 16*pi2)
      xphspmc=1/(1-x)*( 1/(1-y) + 1/(1+y) ) 
c
c MC subt term: pure MC
c
      if(x1.lt.1.and.x2.lt.1.and.s.gt.4*xm2)then
        call invar(xm2,s,x,y,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        zg2_nlo=zgmu2_nlo()
        zg6_mc=zgmu6_mc(zg2_nlo)
      else
        zg6_mc=0.d0
      endif
      ipdfscale=2
      gfflag=.false.
      do jproc=loproc,maproc
        call xmcsubtpp(x1,x2,xm2,s,x,y,cth1,cth2,
     #    xmcxsec,xmce0sq,xmcz,xqrksc,xqbrsc,flxsec,flagmc,
     #    gfactsf,gfactcl)
        gfsf(jproc)=gfactsf
        gfcl(jproc)=gfactcl
        gfflag=gfflag.or.(gfactsf.lt.1.d0)
        if(flagmc)then
          do ileg=1,4
            do ie0sq=1,3
              if(flxsec(ileg,ie0sq))then
                if(ileg.eq.1)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)/zhwfct
                  x2t=x2soft(x1,x2,x,y)
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1.and.
     #                  (x*tau).lt.1
                elseif(ileg.eq.2)then
                  zhwfct=xmcz(ileg,ie0sq)
                  x1t=x1soft(x1,x2,x,y)
                  x2t=x2soft(x1,x2,x,y)/zhwfct
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  x1.lt.1.and.x2.lt.1.and.
     #                  (x*tau).lt.1
                else
                  zhwfct=1.d0
                  x1t=x1soft(x1,x2,x,y)
                  x2t=x2soft(x1,x2,x,y)
                  betamc=bx
                  fx1x2=x1t.lt.1.and.x2t.lt.1.and.
     #                  (x**2*tau).lt.1
                endif
                if(fx1x2)then
                  x1x2j=x1x2jac(x1,x2,x,y,ione)/zhwfct
                  www=zg6_mc*xnormmc*xphspmc*x1x2j*betamc*
     #                xmcxsec(ileg,ie0sq)
                  call strfun(x1t,x2t,sf)
                  do i=1,4
                    do itype=1,nl
                      vv(i,jproc,itype)=vv(i,jproc,itype)+
     #                  sf(i,jproc,itype)*www
                    enddo
                  enddo
                endif
              endif
            enddo
          enddo
        endif
      enddo
c
c Counter-event (x,y)=(x,1) and MC subt term, collinear ME part
c
      if( (xktrel.lt.etacut.and.y.gt.1-delta) .or. gfflag) then
         dfact=0.d0
         if(xktrel.lt.etacut.and.y.gt.1-delta)dfact=1.d0
         ytmp=1.d0
         x1t=x1soft(x1,x2,x,y)/x
         x2t=x2soft(x1,x2,x,y)
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)/x
            call invar(xm2,s,x,ytmp,cth1,cth2,'p1',
     #                 tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gcpfact=0.d0
              if(y.gt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #          gcpfact=1-gfsf(jproc)
              prc=xproc(jproc)
              xintcp=0.d0
              xcplus=0.d0
              f=fpp(s,x,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              xintcp=zg2_nlo**3*xnorm*xphspcp*x1x2j*f*(gcpfact-dfact)
              if(dfact.eq.1.d0)then
c Adding the collinear contribution
                xlmude=log(s/xmuf2h1)+log(delta/2)+
     #                 log( (1-bdelta(x))/delta )
                xnormc=zg2_nlo**3*x1x2j*xjac/(16*pi2*(1-bdelta(x)))
                xcplus=xnormc*bx/(1-x)*(
     #                        ppcolp(ytmp,s,q2q,x,xm2,xlmude)
     #                +xlgomx*ppcoll(ytmp,s,q2q,x,xm2) )
              endif
              www=xintcp + xcplus
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Counter-event (x,y)=(x,-1) and MC subt term, collinear ME part
c
      if( (xktrel.lt.etacut.and.y.lt.-1+delta) .or. gfflag) then
         dfact=0.d0
         if(xktrel.lt.etacut.and.y.lt.-1+delta)dfact=1.d0
         ytmp=-1.d0
         x1t=x1soft(x1,x2,x,y)
         x2t=x2soft(x1,x2,x,y)/x
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)/x
            call invar (xm2,s,x,ytmp,cth1,cth2,'p1',
     #            tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gcmfact=0.d0
              if(y.lt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #          gcmfact=1-gfsf(jproc)
              prc=xproc(jproc)
              xintcm=0.d0
              xcminus=0.d0
              f=fpp(s,x,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
              xintcm=zg2_nlo**3*xnorm*xphspcm*x1x2j*f*(gcmfact-dfact)
              if(dfact.eq.1.d0)then
c Adding the collinear contribution
                xlmude=log(s/xmuf2h2)+log(delta/2)+
     #                 log( (1-bdelta(x))/delta )
                xnormc=zg2_nlo**3*x1x2j*xjac/(16*pi2*(1-bdelta(x)))
                xcminus=xnormc*bx/(1-x)*(
     #                        ppcolp(ytmp,s,q1q,x,xm2,xlmude)
     #                +xlgomx*ppcoll(ytmp,s,q1q,x,xm2) )
              endif
              www=xintcm + xcminus
              call strfun(x1t,x2t,sf)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Soft counter-events, and MC subt term, soft and soft-collinear ME parts
c
      if( (xktrel.lt.etacut.and.x.gt.rotildx) .or. gfflag) then
         dfact=0.d0
         if(xktrel.lt.etacut.and.x.gt.rotildx)dfact=1.d0
         xtmp=1.d0
         x1t=x1soft(x1,x2,x,y)
         x2t=x2soft(x1,x2,x,y)
         if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
            ipdfscale=1
            x1x2j=x1x2jac(x1,x2,x,y,ione)
            call invar(xm2,sx,xtmp,y,cth1,cth2,'p1',
     #           tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            zg2_nlo=zgmu2_nlo()
            do jproc=loproc,maproc
              gsffact=0.d0
              if(x1.lt.1.and.x2.lt.1)gsffact=1-gfsf(jproc)
              prc=xproc(jproc)
              xints(jproc)=0.d0
              xborn(jproc)=0.d0
              xcs(jproc)=0.d0
              xsv(jproc)=0.d0
              f=fpp(sx,xtmp,y,xm2,q1q,q2q,w1h,w2h,cth2)
              xintffs=zg2_nlo**3*xnorm*xphsp*(x1x2j/x)*f*(gsffact-dfact)
              xints(jproc)=xints(jproc)+xintffs
              if(dfact.eq.1.d0)then
c Adding the soft-virtual contribution
                xnormsv=zg2_nlo**3*x1x2j*xjac/ 
     #                  (32*pi2*16*pi2*(1-rotildx+svn(rotildx)))
                delppsv=-bsub(sx,xm2,cth1,etacut,jproc)/(4.d0*sx)
                xsv(jproc)=xnormsv * bx * ( delppsv +
     #                ppsv(sx,q1q,xm2,xmur2,xmuf2h1,xmuf2h2) )
c Adding the Born term
                xnormb=zg2_nlo**2*x1x2j*xjac/
     #                 (16*pi*2*pi*(1-rotildx+svn(rotildx)))
                bbb=hvqborn(sx,q1q,xm2,jproc)
                xborn(jproc)=xborn(jproc)+xnormb*bx*bbb*( 1.d0+
     #            xnlfscheme(xm2,xmur2,xmuf2h1,xmuf2h2,zg2_nlo,jproc) )
              endif
            enddo
            if( (xktrel.lt.etacut.and.y.gt.1-delta) .or. gfflag ) then
              ipdfscale=1
              dfact1=0.d0
              if(xktrel.lt.etacut.and.y.gt.1-delta.and.
     #           dfact.eq.1.d0)dfact1=1.d0
              ytmp=1.d0
              call invar(xm2,sx,xtmp,ytmp,cth1,cth2,'p1',
     #          tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
              zg2_nlo=zgmu2_nlo()
              do jproc=loproc,maproc
                gcpfact=0.d0
                if(y.gt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #            gcpfact=1-gfsf(jproc)
                prc=xproc(jproc)
                f=fpp(sx,xtmp,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
                xintcps=-zg2_nlo**3*xnorm*xphspcp*(x1x2j/x)*f*
     #                  (gcpfact-dfact1)
                xints(jproc)=xints(jproc)+xintcps
                if(dfact1.eq.1.d0)then
c Adding the collinear contribution
                  xlmude=log(sx/xmuf2h1)+log(delta/2)+
     #                   log( (1-bdelta(x))/delta )
                  xnormc=zg2_nlo**3*x1x2j*xjac/(16*pi2*(1-bdelta(x)))
                  xcs(jproc)=xcs(jproc) - xnormc*bx/(1-x)*(
     #                       ppcolp(ytmp,sx,q2q,xtmp,xm2,xlmude)
     #               +xlgomx*ppcoll(ytmp,sx,q2q,xtmp,xm2) )
                endif
              enddo
            endif
            if( (xktrel.lt.etacut.and.y.lt.-1+delta) .or. gfflag) then
              ipdfscale=1
              dfact1=0.d0
              if(xktrel.lt.etacut.and.y.lt.-1+delta.and.
     #           dfact.eq.1.d0)dfact1=1.d0
              ytmp=-1.d0
              call invar(xm2,sx,xtmp,ytmp,cth1,cth2,'p1',
     #          tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
              zg2_nlo=zgmu2_nlo()
              do jproc=loproc,maproc
                gcmfact=0.d0
                if(y.lt.0.d0.and.x1.lt.1.and.x2.lt.1)
     #            gcmfact=1-gfsf(jproc)
                prc=xproc(jproc)
                f=fpp(sx,xtmp,ytmp,xm2,q1q,q2q,w1h,w2h,cth2)
                xintcms=-zg2_nlo**3*xnorm*xphspcm*(x1x2j/x)*f*
     #                  (gcmfact-dfact1)
                xints(jproc)=xints(jproc) + xintcms
                if(dfact1.eq.1.d0)then
c Adding the collinear contribution
                  xlmude=log(sx/xmuf2h2)+log(delta/2)+
     #                   log( (1-bdelta(x))/delta )
                  xnormc=zg2_nlo**3*x1x2j*xjac/(16*pi2*(1-bdelta(x)))
                  xcs(jproc)=xcs(jproc) - xnormc*bx/(1-x)*(
     #                       ppcolp(ytmp,sx,q1q,xtmp,xm2,xlmude)
     #               +xlgomx*ppcoll(ytmp,sx,q1q,xtmp,xm2) )
                endif
              enddo
            endif
c
            call strfun(x1t,x2t,sf)
            do jproc=loproc,maproc
              www=xints(jproc)+xsv(jproc)+xborn(jproc)+xcs(jproc)
              do i=1,4
                do itype=1,nl
                  vv(i,jproc,itype)=vv(i,jproc,itype)+
     #              sf(i,jproc,itype)*www
                enddo
              enddo
            enddo
         endif
      endif
c
c Compute Born ME times luminosity to get flavour assignment; the
c normalization is irrelevant
c
      xtmp=1.d0
      x1t=x1soft(x1,x2,x,y)
      x2t=x2soft(x1,x2,x,y)
      if(x1t.lt.1.and.x2t.lt.1.and.(x*tau).lt.1)then
        ipdfscale=1
        call invar(xm2,sx,xtmp,y,cth1,cth2,'p1',
     #      tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        zg2_nlo=zgmu2_nlo()
        call strfun(x1t,x2t,sf)
        if(loproc.eq.3.and.maproc.eq.3)then
          iproclo=1
          iprocma=2
        else
          iproclo=loproc
          iprocma=maproc
        endif
        do iproc=iproclo,iprocma
          prc=xproc(iproc)
          www=hvqborn(sx,q1q,xm2,iproc)
          do i=1,4
            do itype=1,nl
              vvs(i,iproc,itype)=sf(i,iproc,itype)*www
            enddo
          enddo
        enddo
      endif
c
      call checkvv(xsum,dummy,iret)
      if(iret.eq.1)then
        xtmp=1.d0
        ytmp=1.d0
        call invar(xm2,sx,xtmp,ytmp,cth1,cth2,str,
     #             tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        x1t=x1soft(x1,x2,x,y)
        x2t=x2soft(x1,x2,x,y)
        ycm=0.5d0*log(x1t/x2t)
        tau=x*tau
        if(idec.eq.0)then
          ps=sx
          px=xtmp
          py=ytmp
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
      tot5bz=abs(xint)
c
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
      integer i0,jproc0,itype0
      common/cidproc/i0,jproc0,itype0
      integer idec
      common/cidec/idec
      integer iret
      real*8 ycm0
c
      call xout(iret)
      if(iret.eq.1)then
        if(idec.eq.0)call getspincorr(jproc0)
        if(i0.eq.1)then
          call labmom(ycm)
          ycm0=ycm
        elseif(i0.eq.2)then
          call conjug(ycm)
          ycm0=ycm
        elseif(i0.eq.3)then
          call reflex(ycm)
          ycm0=-ycm
        elseif(i0.eq.4)then
          call refcon(ycm)
          ycm0=-ycm
        else
          write(*,*)'Fatal error in sprfin'
          stop
        endif
        call getx1x2(tau,ycm0)
        call getmom(tau,ycm0,i0)
        call store_events(iunit,xone)
      endif
      return
      end


      subroutine labmom(y)
c boost CM momenta to the lab system
c
      implicit none
      real * 8 y
      real * 8 yq10,yq20,yp0,pq10,pq20,pp0
      common/ycmvar/yq10,yq20,yp0
      common/perpen/pq10(2),pq20(2),pp0(2)
      include 'hvqcblks.h'
      integer j
      yq1 = yq10 + y
      yq2 = yq20 + y
      yp  = yp0  + y
      do j=1,2
         pq1(j) = pq10(j)
         pq2(j) = pq20(j)
         pp(j)  = pp0(j)
      enddo
      return
      entry conjug(y)
      yq1 = yq20 + y
      yq2 = yq10 + y
      yp  = yp0  + y
      do j=1,2
         pq1(j) = pq20(j)
         pq2(j) = pq10(j)
         pp(j)  = pp0(j)
      enddo
      return
      entry reflex(y)
      yq1 = - yq10 - y
      yq2 = - yq20 - y
      yp  = - yp0  - y
      do j=1,2
         pq1(j) = - pq10(j)
         pq2(j) = - pq20(j)
         pp(j)  = - pp0(j)
      enddo
      return
      entry refcon(y)
      yq1 = - yq20 - y
      yq2 = - yq10 - y
      yp  = - yp0  - y
      do j=1,2
         pq1(j) = - pq20(j)
         pq2(j) = - pq10(j)
         pp(j)  = - pp0(j)
      enddo
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


      subroutine getmom(tau,ycm,i0)
      implicit none
      integer i0
      real*8 tau,ycm
      include 'hvqcblks.h'
      real*8 pi
      parameter (pi=3.14159265358979312D0)
      integer i,j,k,imax,itype
      real*8 xsign,xtmp,theta,cth,sth,fk88random,sqsh,ycmnew
      real*8 x1,x2
      common/cx1x2/x1,x2
      real*8 tq12,tq22
      common/ctvirt/tq12,tq22
      real*8 q12,q22
      common/cvirt/q12,q22
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_lb(11,4)
      common/cxmomlb/xmom_lb
      real*8 xmom_prime(11,4)
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
      if(i0.le.2)then
        xsign=1.d0
      elseif(i0.le.4)then
        xsign=-1.d0
      else
        write(*,*)'Fatal error in getmom'
        stop
      endif
      imax=5
      if(idec.eq.0)imax=11
      do j=1,3
        xmom_cm(3,j)=xsign*xmom_cm(3,j)
        if(i0.eq.1.or.i0.eq.3)then
          do i=4,imax
            xmom_cm(i,j)=xsign*xmom_cm(i,j)
          enddo
        else
          xtmp=xsign*xmom_cm(5,j)
          xmom_cm(5,j)=xsign*xmom_cm(4,j)
          xmom_cm(4,j)=xtmp
          if(idec.eq.0)then
            do k=1,3
              xtmp=xsign*xmom_cm(k+8,j)
              xmom_cm(k+8,j)=xsign*xmom_cm(k+5,j)
              xmom_cm(k+5,j)=xtmp
            enddo
          endif
        endif
      enddo
      if(i0.eq.2.or.i0.eq.4)then
        xtmp=xmom_cm(5,4)
        xmom_cm(5,4)=xmom_cm(4,4)
        xmom_cm(4,4)=xtmp
        if(idec.eq.0)then
          do k=1,3
            xtmp=xmom_cm(k+8,4)
            xmom_cm(k+8,4)=xmom_cm(k+5,4)
            xmom_cm(k+5,4)=xtmp
          enddo
        endif
c Exchange W and top invariant masses as well: used by put_on_shell
        xtmp=q22
        q22=q12
        q12=xtmp
        xtmp=tq22
        tq22=tq12
        tq12=xtmp
      endif
c perform a random rotation in the transverse plane
      theta=2*pi*fk88random(ifk88seed)
      cth=cos(theta)
      sth=sin(theta)
      call transrot(cth,sth,pq1(1),pq1(2))
      call transrot(cth,sth,pq2(1),pq2(2))
      call transrot(cth,sth,pp(1),pp(2))
      do i=3,imax
        call transrot(cth,sth,xmom_cm(i,1),xmom_cm(i,2))
      enddo
      if(ichkmom.eq.0)call checkmom(xmom_cm,sh,0.d0,3,2)
c determine colour connections
      call getcolconn()
c put partons on Herwig mass shell
      if(ionshell.eq.0.and.ideconsh.eq.0)then
c keep the parton massless; getspincorr may change shat, so boost all
        sqsh=sqrt(sh)
        ycmnew=ycm
        do i=1,imax
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
        itype=idec+1
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
      real*8 xmom_lb(11,4)
      common/cxmomlb/xmom_lb
      real*8 xmss(1:11)
      common/procmass/xmss
      integer i
c
      do i=1,11
        if(xmom_lb(i,4).ne.0.d0)xmom_lb(i,4)=xmss(i)
      enddo
      return
      end


      subroutine setxmss()
c Fills the common block xmss. Used only if put_on_shell is not called;
c thus, set all masses equal to zero, except those of the primary top and W
      implicit none
      include 'hvqcblks.h'
      integer i
      real*8 tq12,tq22
      common/ctvirt/tq12,tq22
      real*8 xmss(1:11)
      common/procmass/xmss
      integer idec
      common/cidec/idec
c
      do i=1,11
        xmss(i)=0.d0
      enddo
      if(idec.eq.0)then
        xmss(4) = sqrt(tq12)
        xmss(5) = sqrt(tq22)
      elseif(idec.eq.1)then
        xmss(4) = sqrt(xm2)
        xmss(5) = sqrt(xm2)
      endif
      return
      end


      subroutine getcolconn()
c Determines colour connections. Derived from Bryan's subroutine UPFLOW
      implicit none
      include 'hvqcblks.h'
      real*8 xm,t1r,t2r,trn,crnd,fk88random,dotprod,s,tk,uk,q1q,q2q,
     #  s2,q1c,q2c,w1,w2,t(6)
      integer i
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_cross(5,4)
      common/cxmomcross/xmom_cross
      integer i1hpro
      common/ci1hpro/i1hpro
      integer iccode
      common/ciccode/iccode
      integer ifk88seed
      common/cifk88seed/ifk88seed
c
      if(xmom_cm(3,4).eq.0.d0)then
c 2-body kinematics
        if(i1hpro.eq.401)then
          iccode=1
        elseif(i1hpro.eq.403)then
          iccode=2
        elseif(i1hpro.eq.407)then
          crnd=fk88random(ifk88seed)
          t1r=
     # dotprod(xmom_cm(1,1),xmom_cm(1,2),xmom_cm(1,3),xmom_cm(1,4),
     #         xmom_cm(5,1),xmom_cm(5,2),xmom_cm(5,3),xmom_cm(5,4))
          t2r=
     # dotprod(xmom_cm(1,1),xmom_cm(1,2),xmom_cm(1,3),xmom_cm(1,4),
     #         xmom_cm(4,1),xmom_cm(4,2),xmom_cm(4,3),xmom_cm(4,4))
          iccode=3
          if( (t1r**2).lt.(crnd*(t1r**2+t2r**2)) )iccode=4
        else
          write(*,*)'Fatal error #1 in getcolconn: i1hpro=',i1hpro
          stop
        endif
      else
c 3-body kinematics
        if(i1hpro.lt.401.or.i1hpro.gt.407)then
          write(*,*)'Fatal error #2 in getcolconn: i1hpro=',i1hpro
          stop
        endif
        crnd=fk88random(ifk88seed)
        call xcrossing(i1hpro)
        s=   2*dotprod(xmom_cross(1,1),xmom_cross(1,2),
     #                 xmom_cross(1,3),xmom_cross(1,4),
     #                 xmom_cross(2,1),xmom_cross(2,2),
     #                 xmom_cross(2,3),xmom_cross(2,4))
        tk= -2*dotprod(xmom_cross(1,1),xmom_cross(1,2),
     #                 xmom_cross(1,3),xmom_cross(1,4),
     #                 xmom_cross(3,1),xmom_cross(3,2),
     #                 xmom_cross(3,3),xmom_cross(3,4))
        uk= -2*dotprod(xmom_cross(2,1),xmom_cross(2,2),
     #                 xmom_cross(2,3),xmom_cross(2,4),
     #                 xmom_cross(3,1),xmom_cross(3,2),
     #                 xmom_cross(3,3),xmom_cross(3,4))
        q1q=-2*dotprod(xmom_cross(1,1),xmom_cross(1,2),
     #                 xmom_cross(1,3),xmom_cross(1,4),
     #                 xmom_cross(4,1),xmom_cross(4,2),
     #                 xmom_cross(4,3),xmom_cross(4,4))
        q2q=-2*dotprod(xmom_cross(2,1),xmom_cross(2,2),
     #                 xmom_cross(2,3),xmom_cross(2,4),
     #                 xmom_cross(5,1),xmom_cross(5,2),
     #                 xmom_cross(5,3),xmom_cross(5,4))
        s2=s+tk+uk 
        q1c=-s-tk-q1q 
        q2c=-s-uk-q2q
        w1=-q1q+q2q-tk 
        w2=q1q-q2q-uk 
        xm=sqrt(xm2)
        if(i1hpro.lt.407)then
          call qqbplanar(s,tk,uk,q1q,q2q,s2,q1c,q2c,w1,w2,xm,t1r,t2r)
          if (t1r.gt.(t1r+t2r)*crnd) then
            iccode=2*(i1hpro-400)-1
          else
            iccode=2*(i1hpro-400)
          endif
        else
          call ggplanar(s,tk,uk,q1q,q2q,s2,q1c,q2c,w1,w2,xm,
     #                  t(1),t(2),t(3),t(4),t(5),t(6))
          do i=2,6
            t(i)=t(i)+t(i-1)
          enddo
          trn=t(6)*crnd
          do i=1,5
            if (trn.lt.t(i)) goto 10
          enddo
          i=6
 10       iccode=i+12
        endif
      endif
      return
      end


      subroutine xcrossing(i1hpro)
c Crosses parton 4-momenta, in order to determine colour connections for
c 2->3 processes. Derived from Bryan's subroutine UPFLOW
      integer i1hpro
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_cross(5,4)
      common/cxmomcross/xmom_cross
      real*8 xsign
      integer ihpro,i,j,ip,icros(1:3,1:7)
      data icros/
     #  1, 2, 3, 
     #  1,-3,-2, 
     #  2, 1, 3,
     # -3, 1,-2, 
     #  2,-3,-1,
     # -3, 2,-1,
     #  1, 2, 3/
c
      ihpro=i1hpro-400
c heavy quark 4-momenta are not affected
      do i=4,5
        do j=1,4
          xmom_cross(i,j)=xmom_cm(i,j)
        enddo
      enddo
c cross parton 4-momenta
      do i=1,3
        ip=icros(i,ihpro)
        xsign=1.d0
        if(ip.lt.0)xsign=-1.d0
        do j=1,4
          xmom_cross(i,j)=xsign*xmom_cm(iabs(ip),j)
        enddo
      enddo
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
      include 'hvqcblks.h'
      integer i2b,i,j,it,il,in,ib,ii
      real*8 xmss(1:11),xtmp(1:4),xk1tmp(1:4),ytmp1(1:4),ytmp2(1:4),
     #  xavg3(1:3),wvec(1:4),wvec2(1:4)
      real*8 ycm,ycmnew,pi,one,delta_thrs,shat,xkp2prime_norm2,
     #  xkp2prime_norm,xkprime_0,xsign,xnorm_3,delta,gamma,xmprime,
     #  xk1prime_norm,fakemass,xk1tmp_norm,xkprime_norm,xavgnorm,
     #  xnormsq,xbwnorm,xlepnorm,tmplmass,qw2,qw
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      parameter (delta_thrs=0.5d-3)
      real*8 tq12,tq22
      common/ctvirt/tq12,tq22
      real*8 q12,q22
      common/cvirt/q12,q22
      common/procmass/xmss
      real*8 xmass(-5:21)
      common/parmass/xmass
c W mass and width
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
c partons, 3 is the outgoing parton, 4 is Q, 5 is Qbar. When the tops
c decay, 6=l+, 7=nu, 8=b are the decay products of the top, 9=l-, 10=nubar,
c 11=bbar are the decay products of the tbar. Momentum conservation is 
c (1+2)-(3+4+5)=0 or (1+2)-(3+6+7+8+9+10+11)=0
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
c new momenta (put on shell) are stored here
      real*8 xmom_prime(11,4)
      common/cxmomprime/xmom_prime
c ipX is the parton code relevant to parton # X. PDG conventions are
c used: 1=d, 2=u, 3=s, 4=c, 5=b, 21=g
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
      integer ip4,ip5,ip6,ip7,ip8,ip9
      common/ci2part/ip4,ip5,ip6,ip7,ip8,ip9
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
        xmss(5) = sqrt(tq22)
        if(ideconsh.eq.0)then
          do i=6,11
            xmss(i) = 0.d0
          enddo
        elseif(ideconsh.eq.2)then
          xmss(6) = xlep1mass(1)
          xmss(7) = xlep2mass(1)
          xmss(8) = xmass(ip6)
          xmss(9) = xlep1mass(2)
          xmss(10) = xlep2mass(2)
          xmss(11) = xmass(ip9)
        elseif(ideconsh.eq.12)then
          if( abs(ip4).ge.11.and.abs(ip4).le.16 )then
            xmss(6) = xlep1mass(1)
            xmss(7) = xlep2mass(1)
          else
            xmss(6) = 0.d0
            xmss(7) = 0.d0
          endif
          xmss(8) = 0.d0
          if( abs(ip7).ge.11.and.abs(ip7).le.16 )then
            xmss(9) = xlep1mass(2)
            xmss(10) = xlep2mass(2)
          else
            xmss(9) = 0.d0
            xmss(10) = 0.d0
          endif
          xmss(11) = 0.d0
        else
          write(*,*)'Error in put_on_shell: unknown ideconsh',ideconsh
          stop
        endif
      elseif(idec.eq.1)then
        xmss(4) = sqrt(xm2)
        xmss(5) = sqrt(xm2)
        do i=6,11
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
        write(*,*)'Fatal error in put_on_shell: unknown ionshell'
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
c of the masses of the heavy quarks
      delta=sqrt(xkprime_0**2-xkprime_norm**2)-xmss(4)-xmss(5)
      if(delta.lt.delta_thrs)then
c parton 3-momenta cannot be kept fixed: the total available energy
c is not sufficient; modify 3-momenta of the incoming partons
        gamma=sqrt( (xmss(4)+xmss(5)+delta_thrs)**2+xkprime_norm**2 )+
     #        xmom_prime(3,4)
        if(gamma.lt.(xmss(1)+xmss(2)))then
          write(6,*)'Fatal error #0 in put_on_shell'
          write(6,*)gamma,xmom_prime(3,4)
          stop
        endif
        xkp2prime_norm2=( gamma**2-2*(xmss(1)**2+xmss(2)**2)+
     #                    (xmss(1)**2-xmss(2)**2)**2/gamma**2 )/4.d0
        xkp2prime_norm=sqrt(xkp2prime_norm2)
        xmom_prime(1,3)=sign(1.d0,xmom_cm(1,3))*xkp2prime_norm
        xmom_prime(1,4)=sqrt(xkp2prime_norm2+xmss(1)**2)
        xmom_prime(2,3)=sign(1.d0,xmom_cm(2,3))*xkp2prime_norm
        xmom_prime(2,4)=sqrt(xkp2prime_norm2+xmss(2)**2)
        xkprime_0=xmom_prime(1,4)+xmom_prime(2,4)-xmom_prime(3,4)
        shat=(xmom_prime(1,4)+xmom_prime(2,4))**2 -
     #       (xmom_prime(1,3)+xmom_prime(2,3))**2
      endif
c now the parton 3-momenta have been defined in such a way
c that the momenta of the heavy quarks can be transformed.
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
          write(6,*)'Fatal error #1 in put_on_shell'
          write(6,*)i,xmss(i),fakemass
          stop
        endif
        xk1tmp_norm=xnorm_3(xk1tmp)
c xavg is the direction along which the Q1 and Q2 momenta are placed
c in the new QQ rest frame. It is arbitrarily defined by averaging 
c (hence the 1/2 in the definition) the directions of the original 
c Q1 and Q2 momenta. It may not have modulus 1, so normalize it
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
        xk1tmp(4)=sqrt(xk1prime_norm**2+xmss(i)**2)
        call xhwulob(xtmp,xmprime,
     #               xk1tmp,xmss(i),
     #               ytmp2,fakemass)
        if(abs(fakemass-xmss(i)).gt.1.d-4)then
          write(6,*)'Fatal error #2 in put_on_shell'
          write(6,*)i,xmss(i),fakemass
          stop
        endif
        call getvec(ytmp2,xmom_prime(i,1),xmom_prime(i,2),
     #                    xmom_prime(i,3),xmom_prime(i,4))
      enddo
      if(idec.eq.0)then
        do it=4,5
          il=it+2+2*(it-4)
          in=it+3+2*(it-4)
          ib=it+4+2*(it-4)
          call fillvec(xmom_prime(it,1),xmom_prime(it,2),
     #                 xmom_prime(it,3),xmom_prime(it,4),xtmp)
c First deal with the Wb pair; define W momentum, and compute W mass
c (when iwidth=1, W is off shell)
          call vecsum(xmom_cm(il,1),xmom_cm(il,2),
     #                xmom_cm(il,3),xmom_cm(il,4),one,
     #                xmom_cm(in,1),xmom_cm(in,2),
     #                xmom_cm(in,3),xmom_cm(in,4),one,wvec)
          qw2=xnormsq(wvec)
          qw=sqrt(qw2)
          if( ichkmom.eq.0 .and. 
     #        ( (iwidth.eq.0.and.abs(qw/xmw-1.d0).gt.1.d-4) .or.
     #          (iwidth.eq.1.and.it.eq.4.and.
     #           abs(qw/sqrt(q12)-1.d0).gt.1.d-4) .or.
     #          (iwidth.eq.1.and.it.eq.5.and.
     #           abs(qw/sqrt(q22)-1.d0).gt.1.d-4) ) )then
            write(6,*)'Error #3 in put_on_shell'
            write(6,*)qw,it,il,in
            stop
          endif
          if( ichkmom.eq.0 .and. iwidth.eq.1 .and.
     #        qw.gt.xmss(it) )then
            write(6,*)'Error #4 in put_on_shell'
            write(6,*)qw,it,il,in
            stop
          endif
          xbwnorm=xmss(it)**2-2*(xmss(ib)**2+qw2)+
     #            (xmss(ib)**2-qw2)**2/xmss(it)**2
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
     #                 xmom_cm(ib,3),xmom_cm(ib,4),ytmp1)
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
     #              (1+xsign*(qw2-xmss(ib)**2)/xmss(it)**2)
          call xhwulob(xtmp,xmss(it),xk1tmp,qw,wvec2,fakemass)
          xsign=-1.d0
          do j=1,3
            xk1tmp(j)=xsign*xbwnorm*xavg3(j)
          enddo
          xk1tmp(4)=xmss(it)/2.d0*
     #              (1+xsign*(qw2-xmss(ib)**2)/xmss(it)**2)
          call xhwulob(xtmp,xmss(it),xk1tmp,xmss(ib),ytmp2,fakemass)
          call getvec(ytmp2,xmom_prime(ib,1),xmom_prime(ib,2),
     #                      xmom_prime(ib,3),xmom_prime(ib,4))
c Next deal with the lepton pair; W has momentum wvec2
          xlepnorm=qw2-2*(xmss(il)**2+xmss(in)**2)+
     #             (xmss(il)**2-xmss(in)**2)**2/qw2
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
     #                   xmom_cm(ii,3),xmom_cm(ii,4),ytmp1)
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
     #        (1+xsign*(xmss(il)**2-xmss(in)**2)/qw2)
            call xhwulob(wvec2,qw,xk1tmp,tmplmass,ytmp2,fakemass)
            call getvec(ytmp2,xmom_prime(ii,1),xmom_prime(ii,2),
     #                        xmom_prime(ii,3),xmom_prime(ii,4))
          enddo
        enddo
      else
        do i=6,11
          do j=1,4
            xmom_prime(i,j)=0.d0
          enddo
        enddo
      endif
      if(ichkmom.eq.0)then
        if(idec.eq.0)then
          call checktdec2(xmom_prime,4,6,7,8)
          call checktdec2(xmom_prime,5,9,10,11)
          call checkmom(xmom_prime,shat,0.d0,4,1)
        else
          call checkmom(xmom_prime,shat,0.d0,4,2)
        endif
        if(xmass(1).eq.0.and.xmass(2).eq.0.and.xmass(3).eq.0.and.
     #     xmass(4).eq.0.and.xmass(5).eq.0.and.xmass(21).eq.0.and.
     #     xlep1mass(1).eq.0.and.xlep2mass(1).eq.0.and.
     #     xlep1mass(2).eq.0.and.xlep2mass(2).eq.0)then
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
      e1=sqrt(p13**2+xm1**2)
      e2=sqrt(p23**2+xm2**2)
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
      REAL*8 PF4,FN,PS(5),PI(4),PF(4)
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
      REAL*8 PF4,FN,PS(5),PI(4),PF(4)
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
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      real*8 xmom_prime(11,4)
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
      if(idec.eq.0)imax=11
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
      integer iret,iretvv,iretvvs,iproc,iproclo,iprocma,i,itype,
     #  iwh,iflag,ihpro,i1,i2,i3,i1hproo,ip1o,ip2o,ip3o
      real*8 wwx(4,3,5),xsum,xsumabs,xsumvvs,xsumabsvvs,xstsign,
     #  xg,wh,rmax,fk88random
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      integer ifuntype
      common/cifuntype/ifuntype
      real*8 vv(4,3,5)
      common/cvv/vv
      real*8 vvs(4,3,5)
      common/cvvs/vvs
      integer iwrong,iwrong1
      common/ciwrong/iwrong,iwrong1
      integer i0,jproc0,itype0
      common/cidproc/i0,jproc0,itype0
      integer ivbhpro(4,3,5)
      common/civbhpro/ivbhpro
      integer idp1(4,3,5),idp2(4,3,5),idp3(4,3,5)
      common/cidpart/idp1,idp2,idp3
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
      integer nl
      common/nl/nl
      integer ifk88seed
      common/cifk88seed/ifk88seed
      integer idec
      common/cidec/idec
c
      i0=0
      jproc0=0
      itype0=0
      iret=0
      call checkvv(xsum,xsumabs,iretvv)
      call checkvvs(xsumvvs,xsumabsvvs,iretvvs)
c vvs is computed but not meant to be used; set iretvvs=0 for safety.
c If used, reinstate cvvs common blocks in functions tot5#
      iretvvs=0
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
          if(loproc.eq.3.and.maproc.eq.3)then
            iproclo=1
            iprocma=2
          else
            iproclo=loproc
            iprocma=maproc
          endif
        else
          write(*,*)'Fatal error in xout: ifuntype=',ifuntype
          stop
        endif
        if(iretvvs.eq.1)then
          xsum=xsumvvs
          xsumabs=xsumabsvvs
          do iproc=iproclo,iprocma
            do i=1,4
              do itype=1,nl
                wwx(i,iproc,itype)=vvs(i,iproc,itype)
              enddo
            enddo
          enddo
        else
          do iproc=iproclo,iprocma
            do i=1,4
              do itype=1,nl
                wwx(i,iproc,itype)=vv(i,iproc,itype)
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
          do i=1,4
            do itype=1,nl
              if(iwh.eq.0)then
                wh=wh+abs(wwx(i,iproc,itype))/xsumabs
                if(wh.gt.xg)then
                  i0=i
                  jproc0=iproc
                  itype0=itype
                  iwh=1
                endif
              endif
              if(wwx(i,iproc,itype).ne.0.d0)then
                if(xstsign.ne.sign(1.d0,wwx(i,iproc,itype)))then
                  if(iflag.eq.0)then
                    iwrong=iwrong+1
                    iflag=1
                  endif
                  rmax=max(rmax,abs(wwx(i,iproc,itype)))
                endif
              endif
            enddo
          enddo
        enddo
        if(iflag.eq.1)then
          if(xsum.ne.0.d0)rmax=rmax/xsum
          if(rmax.gt.0.05d0.or.xsum.eq.0.d0)iwrong1=iwrong1+1
        endif
        if(i0.eq.0.or.jproc0.eq.0.or.itype0.eq.0)then
          write(*,*)'Fatal error in xout'
          stop
        endif
        ihpro=ivbhpro(i0,jproc0,itype0)
        i1=idp1(i0,jproc0,itype0)
        i2=idp2(i0,jproc0,itype0)
        i3=idp3(i0,jproc0,itype0)
        call parcrossing(jproc0,ihpro,i1,i2,i3,i1hproo,ip1o,ip2o,ip3o)
        i1hpro=i1hproo
        ip1=ip1o
        ip2=ip2o
        ip3=ip3o
c The tops decay. In such a case, determine the identities of the decay 
c products, and assign them the relevant masses (to be used by put_on_shell)
        if(idec.eq.0)call getpdecids()
      endif
      return
      end


      subroutine getpdecids()
c Determine the identities of top and antitop decay products
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
      integer ip4,ip5,ip6,ip7,ip8,ip9
      common/ci2part/ip4,ip5,ip6,ip7,ip8,ip9
      integer ip4s,ip5s,ip6s,ip7s,ip8s,ip9s
      common/ci2parts/ip4s,ip5s,ip6s,ip7s,ip8s,ip9s
c
      if(il1hw.eq.7.or.il2hw.eq.7)then
        write(*,*)'Error #1 in getpdecids()',il1hw,il2hw
        stop
      endif
      if(inonbtop.eq.0)then
c t->Wb
        ip6=ip6s
        ip9=ip9s
      elseif(inonbtop.eq.1)then
c t->W+any down-type quark
        call dectopuwgt(ipone,ip6)
        call dectopuwgt(imone,ip9)
      else
        write(*,*)'Unknown option in getpdecids()',inonbtop
        stop
      endif
c W+ decay
      if(il1hw.eq.0)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac123)then
          call declepuwgt(-iemutau,ip4,ip5)
        else
          call decqrkuwgt(ipone,ip4,ip5)
        endif
      elseif(il1hw.ge.1.and.il1hw.le.3)then
        ip4=ip4s
        ip5=ip5s
      elseif(il1hw.eq.4)then
        call declepuwgt(-iemu,ip4,ip5)
      elseif(il1hw.eq.5)then
        call decqrkuwgt(ipone,ip4,ip5)
      elseif(il1hw.eq.6)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac12)then
          call declepuwgt(-iemu,ip4,ip5)
        else
          call decqrkuwgt(ipone,ip4,ip5)
        endif
      else
        write(*,*)'Error #2 in getpdecids()',il1hw
        stop
      endif
      if(abs(ip4).le.5.and.abs(ip5).le.5)then
        xlep1mass(1)=xmass(ip4)
        xlep2mass(1)=xmass(ip5)
      elseif( abs(ip4).ge.11.and.abs(ip4).le.16 .and.
     #        abs(ip5).ge.11.and.abs(ip5).le.16 )then
        xlep1mass(1)=pdglepmss(abs(ip4))
        xlep2mass(1)=pdglepmss(abs(ip5))
      else
        write(*,*)'Error #4 in getpdecids()',ip4,ip5
        stop
      endif
c W- decay
      if(il2hw.eq.0)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac123)then
          call declepuwgt(iemutau,ip7,ip8)
        else
          call decqrkuwgt(imone,ip7,ip8)
        endif
      elseif(il2hw.ge.1.and.il2hw.le.3)then
        ip7=ip7s
        ip8=ip8s
      elseif(il2hw.eq.4)then
        call declepuwgt(iemu,ip7,ip8)
      elseif(il2hw.eq.5)then
        call decqrkuwgt(imone,ip7,ip8)
      elseif(il2hw.eq.6)then
        xg=fk88random(ifk88seed)
        if(xg.lt.frac12)then
          call declepuwgt(iemu,ip7,ip8)
        else
          call decqrkuwgt(imone,ip7,ip8)
        endif
      else
        write(*,*)'Error #3 in getpdecids()',il2hw
        stop
      endif
      if(abs(ip7).le.5.and.abs(ip8).le.5)then
        xlep1mass(2)=xmass(ip7)
        xlep2mass(2)=xmass(ip8)
      elseif( abs(ip7).ge.11.and.abs(ip7).le.16 .and.
     #        abs(ip8).ge.11.and.abs(ip8).le.16 )then
        xlep1mass(2)=pdglepmss(abs(ip7))
        xlep2mass(2)=pdglepmss(abs(ip8))
      else
        write(*,*)'Error #5 in getpdecids()',ip7,ip8
        stop
      endif
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


      subroutine parcrossing(jproc0,ihpro,i1,i2,i3,
     #                       i1hproo,ip1o,ip2o,ip3o)
      implicit none
      integer jproc0,ihpro,i1,i2,i3,i1hproo,ip1o,ip2o,ip3o,
     # iallzero,iahprotrans(402:406),ibhprotrans(402:406)
      parameter (iallzero=1)
      real*8 xg,fk88random
      integer ifuntype
      common/cifuntype/ifuntype
      integer ifk88seed
      common/cifk88seed/ifk88seed
      data iahprotrans/401,0,403,407,407/
      data ibhprotrans/407,0,407,403,401/
c
      if( (ifuntype.eq.1) .or. (ifuntype.eq.2.and.jproc0.ne.3) )then
        i1hproo=ihpro
        ip1o=i1
        ip2o=i2
        ip3o=i3
      elseif(ifuntype.eq.2.and.jproc0.eq.3)then
        if(ihpro.eq.401.or.ihpro.eq.403.or.ihpro.eq.407)then
          write(*,*)'Error #1 in parcrossing:',ihpro,i1,i2,i3
          stop
        endif
        xg=fk88random(ifk88seed)
        if(xg.lt.0.5d0)then
          i1hproo=iahprotrans(ihpro)
        else
          i1hproo=ibhprotrans(ihpro)
        endif
        if(i1hproo.ne.407)then
          if(i1.eq.21)then
            ip1o=-i3
            ip2o=i2
            ip3o=i1
          elseif(i2.eq.21)then
            ip1o=i1
            ip2o=-i3
            ip3o=i2
          endif
        else
          ip1o=21
          ip2o=21
          ip3o=21
        endif
        if(i1.ne.21.and.i2.ne.21)then
          write(*,*)'Error #2 in parcrossing:',ihpro,i1,i2,i3
          stop
        endif
      else
        write(*,*)'parcrossing: do not know what to do'
        write(*,*)ifuntype,jproc0
        stop
      endif
      call parcheckfin(i1hproo,ip1o,ip2o,ip3o,iallzero)
      return
      end


      subroutine checkvv(xsum,xsumabs,iret)
c iret=0 -> all vv entries are equal to zero
c iret=1 -> there is at least one entry which is not zero
c xsum is the sum of all the entries of vv
c xsumabs is the sum of the absolute value of all the entries of vv
      implicit none
      integer jproc,iret,i,itype
      integer nl
      common/nl/nl
      integer loproc,maproc
      common/cwchproc/loproc,maproc
      real * 8 vv(4,3,5)
      common/cvv/vv
      real * 8 xsum,xsumabs
c
      xsum=0.d0
      xsumabs=0.d0
      iret=0
      do jproc=loproc,maproc
        do i=1,4
          do itype=1,nl
            if(vv(i,jproc,itype).ne.0.d0)iret=1
            xsum=xsum+vv(i,jproc,itype)
            xsumabs=xsumabs+abs(vv(i,jproc,itype))
          enddo
        enddo
      enddo
      return
      end


      subroutine checkvvs(xsum,xsumabs,iret)
c identical to checkvv, except for the fact that works on vvs instead of vv,
c and jproc is not fixed
      implicit none
      integer jproc,iret,i,itype
      integer nl
      common/nl/nl
      real * 8 vvs(4,3,5)
      common/cvvs/vvs
      real * 8 xsum,xsumabs
c
      xsum=0.d0
      xsumabs=0.d0
      iret=0
      do jproc=1,2
        do i=1,4
          do itype=1,nl
            if(vvs(i,jproc,itype).ne.0.d0)iret=1
            xsum=xsum+vvs(i,jproc,itype)
            xsumabs=xsumabs+abs(vvs(i,jproc,itype))
          enddo
        enddo
      enddo
      return
      end


      subroutine getspincorr(jproc0)
c Determines the lepton momenta, by performing an unweighting using
c the exact real and Born lepton matrix elements. This is done assuming
c idr=1; therefore, this routine must be called before labmom...refcon
c and getmom()
      implicit none
      integer jproc0
      real*8 pi,tolerance,pdftiny,bdredfact
      parameter (pi=3.14159265358979312D0)
      parameter (tolerance=1.d-2)
      parameter (pdftiny=1.d-8)
c Divide the bound by bdredfact to compensate for peculiar behaviour of
c ttbar cross section. May become a process-dependent correction if need be
      parameter (bdredfact=2.d0)
      integer ione
      parameter (ione=1)
      character*2 str,prcsave
      parameter (str='p1')
      include 'hvqcblks.h'
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      real*8 ps,px,py,pcth1,pcth2
      common/cpsave/ps,px,py,pcth1,pcth2
      real*8 tq12,tq22
      common/ctvirt/tq12,tq22
      real*8 q12,q22
      common/cvirt/q12,q22
      real*8 ycm,tau
      common/x1x2/ycm,tau
      real*8 x1,x2
      common/cx1x2/x1,x2
      real*8 tga1,tga1mn,tga1pl,tbw1fmpl,tbw1fmmn,tbw1delf,
     #       ym1low2,ym1upp2
      common/tbw1/tga1,tga1mn,tga1pl,tbw1fmpl,tbw1fmmn,tbw1delf,
     #            ym1low2,ym1upp2
      real*8 tga2,tga2mn,tga2pl,tbw2fmpl,tbw2fmmn,tbw2delf,
     #       ym2low2,ym2upp2
      common/tbw2/tga2,tga2mn,tga2pl,tbw2fmpl,tbw2fmmn,tbw2delf,
     #            ym2low2,ym2upp2
      real*8 xm012,ga1,bw1delf,bw1fmmn
      common/cbw1/xm012,ga1,bw1delf,bw1fmmn
      real*8 xm022,ga2,bw2delf,bw2fmmn
      common/cbw2/xm022,ga2,bw2delf,bw2fmmn
      real*8 xm1low2,xm1upp2,xm2low2,xm2upp2
      common/bounds/xm1low2,xm1upp2,xm2low2,xm2upp2
      real*8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      real*8 xmt,twidth
      common/ctparam/xmt,twidth
      real*8 xmass(-5:21)
      common/parmass/xmass
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
      integer ifuntype
      common/cifuntype/ifuntype
      integer jwidth
      common/cjwidth/jwidth
      integer iwidth
      common/ciwidth/iwidth
      integer ichkmom
      common/cichkmom/ichkmom
      integer neventsuw,nqeventsuw,mqeventsuw,ifailuw
      common/c1iunwgt/neventsuw,nqeventsuw,mqeventsuw,ifailuw
      integer ncntuws,nqcntuws,nmaxuw,nqmaxuw
      common/c2iunwgt/ncntuws,nqcntuws,nmaxuw,nqmaxuw
      integer mqcntuws,mqmaxuw
      common/c3iunwgt/mqcntuws,mqmaxuw
      integer idec
      common/cidec/idec
      integer ifk88seed
      common/cifk88seed/ifk88seed
      character*2 xproc(3)
      common/cxproc/xproc
      real*8 xtmp,prob,spcdamp,rrnd,fk88random,e1,f1,g1,h1,e2,f2,
     # g2,h2,phitq1,cthtq1,phitq2,cthtq2,phiat1,cthat1,phiat2,cthat2,
     # o,p,xbwmass3,rat1,qphsp,q1,q2,tk,uk,q1q,q2q,q1c,q2c,w1,
     # w2,w1h,w2h,xdec,xmadevttb,unxdec,xttb,dmfactb1,dmfact1,
     # dmfactb2,dmfact2,phspfact1,phspfact2,xboundb,rat,q,r,
     # xbwmass3_mod,xp,xpsave,psp,staup,xlstaup,ycmp,x1p,x2p,
     # wm1upp2,wbw1mdpl,wbw1fmpl,wbw1delf,wm2upp2,wbw2mdpl,wbw2fmpl,
     # wbw2delf,tq1,tq2,xlumund,getlumspc,xlumdec,bwfunc_mod,
     # xtq(4),xbq(4),xel(4),xnu(4)
      integer iborn,iproj,icross,jjprc,icntuw,iqcntuw,kp1,kp2,kp3,
     # jqcntuw
c
      if(ichkmom.eq.0.and.jwidth.eq.0)call spccheck(1)
      if(jwidth.eq.1)call getx1x2(tau,ycm)
      if(ifuntype.eq.2)then
        if(px.ne.1.d0.or.xmom_cm(3,4).ne.0.d0)then
          write(*,*)'Error #1 in getspincorr'
          stop
        else
          iborn=0
          iproj=0
          if(jproc0.ne.3)then
            icross=0
          else
            icross=1
          endif
          xtmp=px
        endif
      endif
      if(ifuntype.eq.1)then
        prob=spcdamp(px,py)
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
c When iproj=0, the Born and real kinematics are used to perform unweighting
c for S and H events respectively. When iproj=1, the real kinematics is close 
c to the soft/collinear limits, and the Born is used to unweight. In the case 
c of the qg process, such Born is either gg (q||q <==> y->1) or qq 
c (g||q <==> y->-1). We choose gg (qq) if py>0 (py<0). This strategy,
c which serves to set here the local value of jproc (jjprc), cannot be
c be adopted in the case of S events due to the qg contribution, since
c for such event we always have py=1. For S events, we rather define jjprc
c using the value of i1hpro defined in xout(), where the conversion 
c qg->qqbar or qg->gg has already been made.
c When this routine is called, the parton identities have already been
c determined by xout (ip1,ip2,ip3). This includes reflection/conjugation,
c which must be undone since we assume here idr=1. In other words, the
c matrix elements called here compute the processes
c  gg->QQbarg; qqbar->QQbarg; qg->QQbarq
c (with only the first two at the Born level, excluding the final-state g)
c and therefore the local parton identities (kp1,kp2,kp3) must be
c consistent with such computations; the parton coming from the left (right)
c will have Bjorken x equal to x1 (x2).
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
      if(icross.eq.0)then
        jjprc=jproc0
      elseif(icross.eq.1)then
        if(jproc0.ne.3)then
          jjprc=jproc0
        else
          if(ifuntype.eq.1)then
            if(py.ge.0.d0)then
              jjprc=1
            else
              jjprc=2
            endif
          else
            if(i1hpro.eq.407)then
              jjprc=1
            elseif(i1hpro.eq.401.or.i1hpro.eq.403)then
              jjprc=2
            else
              write(*,*)'Error #10 in getspincorr'
              stop
            endif
          endif
        endif
      else
        write(*,*)'Error #2 in getspincorr'
        stop
      endif
      prcsave=prc
      prc=xproc(jjprc)
c Get local parton identities if top are off-shell
      if(jwidth.eq.1)then
        if(iproj.eq.0)then
          kp1=ip1
          kp2=ip2
          kp3=ip3
          call undorcspc(jjprc,kp1,kp2,kp3)
        elseif(iproj.eq.1)then
          if(ifuntype.ne.1)then
            write(*,*)'Error #5 in getspincorr'
            stop
          else
            if(jproc0.ne.3)then
              kp1=ip1
              kp2=ip2
              kp3=ip3
              call undorcspc(jjprc,kp1,kp2,kp3)
            else
              if(jjprc.eq.1)then
                kp1=21
                kp2=21
                kp3=21
              elseif(jjprc.eq.2)then
                kp1=abs(ip3)
                kp2=-abs(ip3)
                kp3=21
              else
                write(*,*)'Error #6 in getspincorr'
                stop
              endif
            endif
          endif
        else
          write(*,*)'Error #7 in getspincorr'
          stop
        endif
        call checkpidspc(jjprc,kp1,kp2,kp3,jproc0)
      endif
c
      neventsuw=neventsuw+1
      icntuw=0
 100  icntuw=icntuw+1
      e1=fk88random(ifk88seed)
      f1=fk88random(ifk88seed)
      g1=fk88random(ifk88seed)
      h1=fk88random(ifk88seed)
      e2=fk88random(ifk88seed)
      f2=fk88random(ifk88seed)
      g2=fk88random(ifk88seed)
      h2=fk88random(ifk88seed)
      phitq1=2*pi*e1
      cthtq1=-1.d0+2*f1
      phitq2=2*pi*g1
      cthtq2=-1.d0+2*h1
      phiat1=2*pi*e2
      cthat1=-1.d0+2*f2
      phiat2=2*pi*g2
      cthat2=-1.d0+2*h2
 300  continue
      if(jwidth.eq.1)then
        jqcntuw=0
        jqcntuw=jqcntuw+1
        q=fk88random(ifk88seed)
        r=fk88random(ifk88seed)
c Since the upper bound is q-dependent, distribute q according to the
c form of the bound, which is a skewed Breit Wigner
        tq12=xbwmass3_mod(q,xm2,tga1,tga1mn,tga1pl,tbw1delf,tbw1fmmn)
        tq22=xbwmass3_mod(r,xm2,tga2,tga2mn,tga2pl,tbw2delf,tbw2fmmn)
        mqcntuws=mqcntuws+jqcntuw
        if(jqcntuw.gt.mqmaxuw)mqmaxuw=jqcntuw
        mqeventsuw=mqeventsuw+1
        xp=xtmp
        xpsave=px
        psp=min( ((sqrt(tq12)+sqrt(tq22))/(2*xmt))**2*ps , sh )
        staup=sqrt(psp/sh)
        xlstaup=log(staup)
        ycmp=ycm
        if(ycmp.le.xlstaup)ycmp=0.999*xlstaup
        if(ycmp.ge.-xlstaup)ycmp=-0.999*xlstaup
        x1p=staup*exp(ycmp)
        x2p=staup/exp(ycmp)
c Define local upper bounds for W mass ranges
        wm1upp2=(sqrt(tq12)-xmass(5))**2-1.d-1
        if(xm1upp2.gt.wm1upp2)then
          if(xm1low2.gt.wm1upp2)then
            write(*,*)'Error in pair mass range #1: getspincorr'
            write(*,*)xm1low2,xm1upp2,tq12
            stop
          endif
          wbw1mdpl=wm1upp2-xm012
          wbw1fmpl=atan(wbw1mdpl/(sqrt(xm012)*ga1))
          wbw1delf=(wbw1fmpl+bw1fmmn)/pi
        else
          wbw1delf=bw1delf
        endif
        wm2upp2=(sqrt(tq22)-xmass(5))**2-1.d-1
        if(xm2upp2.gt.wm2upp2)then
          if(xm2low2.gt.wm2upp2)then
            write(*,*)'Error in pair mass range #2: getspincorr'
            write(*,*)xm2low2,xm2upp2,tq22
            stop
          endif
          wbw2mdpl=wm2upp2-xm022
          wbw2fmpl=atan(wbw2mdpl/(sqrt(xm022)*ga2))
          wbw2delf=(wbw2fmpl+bw2fmmn)/pi
        else
          wbw2delf=bw2delf
        endif
      else
        tq12=xm2
        tq22=xm2
        xp=xtmp
        xpsave=px
        psp=ps
        x1p=-1.d8
        x2p=-1.d8
        wbw1delf=bw1delf
        wbw2delf=bw2delf
      endif
      tq1=sqrt(tq12)
      tq2=sqrt(tq22)
      if(iwidth.eq.1)then
        iqcntuw=0
 200    iqcntuw=iqcntuw+1
        o=fk88random(ifk88seed)
        p=fk88random(ifk88seed)
c First distribute q's according to the matrix element upper bound,
c which can be done exactly the upper bound being a Breit Wigner
        q12=xbwmass3(o,xm012,ga1,wbw1delf,bw1fmmn)
        q22=xbwmass3(p,xm022,ga2,wbw2delf,bw2fmmn)
c Then reject some of the values generated according to the phase-space
c q-dependent factor. A 1->1+(1->2) phase-space decomposition has been used.
c Much better here than after computing matrix elements. The following
c form works since qphsp is a function decreasing with q2
        rat1=( qphsp(q12,tq12)/qphsp(xm1low2,tq12) )*
     #       ( qphsp(q22,tq22)/qphsp(xm2low2,tq22) )
        rrnd=fk88random(ifk88seed)
        if(rat1.lt.rrnd)goto 200
        nqcntuws=nqcntuws+iqcntuw
        if(iqcntuw.gt.nqmaxuw)nqmaxuw=iqcntuw
        nqeventsuw=nqeventsuw+1
        q1=sqrt(q12)
        q2=sqrt(q22)
      else
        q12=xm012
        q1=sqrt(q12)
        q22=xm022
        q2=sqrt(q22)
      endif
      if(q12.gt.tq12.or.q22.gt.tq22)then
        write(*,*)'Error #8 in getspincorr'
        write(*,*)q12,tq12,q22,tq22
        stop
      endif
      call invar(xm2,ps,xtmp,py,pcth1,pcth2,str,
     #           tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      unxdec=xttb(iborn,jjprc,ione,xm2,ps,xtmp,py,
     #            tk,uk,q1q,q2q,w1h,w2h,pcth2)
      if(jwidth.eq.0)then
        xlumund=1.d0
      else
        xlumund=getlumspc(x1,x2,kp1,kp2)
        call invar2(tq12,tq22,psp,xp,py,pcth1,pcth2,str,
     #              tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      endif
      call gentopdmom(tq1,q1,cthtq1,phitq1,cthtq2,phitq2,
     #                xtq,xbq,xel,xnu,1)
      if(ichkmom.eq.0)call checktdec1(tq1,xtq,xbq,xel,xnu,1)
      call gentopdmom(tq2,q2,cthat1,phiat1,cthat2,phiat2,
     #                xtq,xbq,xel,xnu,2)
      if(ichkmom.eq.0)call checktdec1(tq2,xtq,xbq,xel,xnu,2)
      if(ichkmom.eq.0)call checkmom(xmom_cm,psp,0.d0,10,1)
      xdec=xmadevttb(iborn,jjprc,ione,psp,tk,uk,xmom_cm)*
     #     psp/ps
      if(jwidth.eq.0)then
        xlumdec=1.d0
      else
        xlumdec=getlumspc(x1p,x2p,kp1,kp2)
      endif
      dmfactb1=256*xm2**2/16.d0
      dmfact1=1/(64.d0*sthw2**2)*
     #        1.d0/((q12-xm012)**2+xm012*ga1**2)
      dmfactb2=256*xm2**2/16.d0
      dmfact2=1/(64.d0*sthw2**2)*
     #        1.d0/((q22-xm022)**2+xm022*ga2**2)
      if(jwidth.eq.0)then
        phspfact1=1.d0/(xm2*twidth**2)
        phspfact2=1.d0/(xm2*twidth**2)
      else
        phspfact1=pi/(xmt*tga1)*bwfunc_mod(tq12,xm2,tga1,tga1mn,tga1pl)
        phspfact2=pi/(xmt*tga2)*bwfunc_mod(tq22,xm2,tga2,tga2mn,tga2pl)
      endif
      xboundb=dmfactb1*dmfact1*phspfact1*
     #        dmfactb2*dmfact2*phspfact2
      if(xlumdec.gt.pdftiny.and.xlumund.gt.pdftiny)then
        rat=xdec*xlumdec/((1+tolerance)*unxdec*xlumund*xboundb)
      else
        rat=xdec/((1+tolerance)*unxdec*xboundb)
      endif
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
c configuration)
      if(iproj.eq.0)then
        if(xp.ne.xpsave)then
          write(*,*)'Error #3 in getspincorr'
          stop
        endif
      elseif(iproj.eq.1)then
        if(xp.ne.1.d0)then
          write(*,*)'Error #9 in getspincorr'
          stop
        endif
        call invar2(tq12,tq22,psp,xpsave,py,pcth1,pcth2,str,
     #              tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        call gentopdmom(tq1,q1,cthtq1,phitq1,cthtq2,phitq2,
     #                  xtq,xbq,xel,xnu,1)
        call gentopdmom(tq2,q2,cthat1,phiat1,cthat2,phiat2,
     #                  xtq,xbq,xel,xnu,2)
        if(ichkmom.eq.0)call checkmom(xmom_cm,psp,0.d0,20,1)
      else
        write(*,*)'Error #4 in getspincorr'
        stop
      endif
      if(ichkmom.eq.0.and.jwidth.eq.0)call spccheck(2)
      prc=prcsave
      return
      end


      subroutine undorcspc(jjprc,kp1,kp2,kp3)
c Returns parton identities corresponding to idr=1
      implicit none
      integer jjprc,kp1,kp2,kp3
c
      if(jjprc.eq.1)then
        continue
      elseif(jjprc.eq.2)then
        if(kp1.lt.0)then
          kp1=-kp1
          kp2=-kp2
        endif
      elseif(jjprc.eq.3)then
        if(kp1.eq.21)then
          kp1=abs(kp2)
          kp3=kp1
          kp2=21
        else
          if(kp1.lt.0)then
            kp1=-kp1
            kp3=-kp3
          endif
        endif
      else
        write(*,*)'Error in undorcspc: unknown process',jjprc
        stop
      endif
      return
      end


      subroutine checkpidspc(jjprc,kp1,kp2,kp3,jproc0)
c Checks parton identities after manipulations performed in getspincorr
      implicit none
      integer jjprc,kp1,kp2,kp3,jproc0,ierr
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
c
      ierr=0
      if(jjprc.eq.1)then
        if(kp1.ne.21.or.kp2.ne.21.or.kp3.ne.21)ierr=1
      elseif(jjprc.eq.2)then
        if( kp1.lt.1.or.kp1.gt.5 .or.
     #      kp2.lt.-5.or.kp2.gt.-1 .or.
     #      kp3.ne.21 )ierr=1
      elseif(jjprc.eq.3)then
        if( kp1.lt.1.or.kp1.gt.5 .or.
     #      kp2.ne.21 .or.
     #      kp3.lt.1.or.kp3.gt.5 )ierr=1
      else
        write(*,*)'Error in checkpidspc: unknown process',jjprc
        stop
      endif
      if(ierr.eq.1)then
        write(*,*)'Error in checkpidspc'
        write(*,*)jjprc,kp1,kp2,kp3
        write(*,*)jproc0,ip1,ip2,ip3
        stop
      endif
      return
      end


      function getlumspc(x1,x2,kp1,kp2)
      implicit none
      include 'hvqcblks.h'
      real*8 getlumspc,x1,x2
      integer kp1,kp2
      integer ih1,ih2,ndns1,ndns2
      common/strfun0/ih1,ih2,ndns1,ndns2
      integer invimapp(-5:21)
      common/cinvimapp/invimapp
      real*8 zg2,zgmu2_nlo
      real*4 fh1x1(-5:5),fh2x2(-5:5),smuf2h1,smuf2h2
      integer ip1,ip2
c
      if( x1.lt.0.d0.or.x1.gt.1.d0 .or.
     #    x2.lt.0.d0.or.x2.gt.1.d0 )then
        write(*,*)'Error in getlumspc: Bjorken x',x1,x2
        stop
      endif
      zg2=zgmu2_nlo()
      smuf2h1=sngl(xmuf2h1)
      smuf2h2=sngl(xmuf2h2)
      call mlmpdf(ndns1,ih1,smuf2h1,sngl(x1),fh1x1,5)
      call mlmpdf(ndns2,ih2,smuf2h2,sngl(x2),fh2x2,5)
      ip1=invimapp(kp1)
      ip2=invimapp(kp2)
      if( ip1.lt.-5.or.ip1.gt.5 .or.
     #    ip2.lt.-5.or.ip2.gt.5 )then
        write(*,*)'Error in getlumspc: parton IDs',ip1,ip2,kp1,kp2
        stop
      endif
      getlumspc=dble(fh1x1(ip1)*fh2x2(ip2))
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
      real*8 xmom_cm(11,4)
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


      function xttb(iborn,jproc,idr,xmt2,s,x,y,
     #              tk,uk,q1q,q2q,w1h,w2h,cth2)
c Wrapper for the undecayed matrix elements of the original code.
c For Born matrix elements, q1q is t (consistently with the
c routine invar). This function is called by getspincorr, where only
c idr=1 is considered. Stop if other options are given in input
      implicit none
      real*8 xttb,xmt2,s,x,y,tk,uk,q1q,q2q,w1h,w2h,cth2
      integer iborn,jproc,idr
      real*8 s0,x0,y0,tk0,uk0,q1q0,q2q0,w1h0,w2h0,cth20,t0,
     # ggborn,qqborn,fpp
c
      if(idr.eq.1)then
        s0=s
        if(iborn.eq.0)then
          t0=q1q
        else
          x0=x
          y0=y
          tk0=tk
          uk0=uk
          q1q0=q1q
          q2q0=q2q
          w1h0=w1h
          w2h0=w2h
          cth20=cth2
        endif
      else
        write(*,*)'Error in xttb: use only direct events',idr
        stop
      endif
      if(iborn.eq.0.and.jproc.eq.1)then
        xttb=ggborn(s0,t0,xmt2)
      elseif(iborn.eq.0.and.jproc.eq.2)then
        xttb=qqborn(s0,t0,xmt2)
      else
        xttb=fpp(s0,x0,y0,xmt2,q1q0,q2q0,w1h0,w2h0,cth20)
      endif
      return
      end


      subroutine gentopdmom(xmt,xmw,cth1,phi1,cth2,phi2,
     #                      xtq,xbq,xel,xnu,iqrk)
c Generates the four-momenta of the decay products of the top (if iqrk=1)
c or of the tbar (if iqrk=2). These four-momenta are returned in the top/tbar 
c rest frame (xbq, xel, xnu; the trivial top/tbar momentum is returned as 
c well, xtq). The four-momenta are also boosted to the frame in which the 
c top/tbar has momentum xmom_cm(4,*)/xmom_cm(5,*), and the common block 
c xmomcm is filled according to the identifications
c   l+ --> xmom_cm(6,*), nu --> xmom_cm(7,*), b --> xmom_cm(8,*), 
c   l- --> xmom_cm(9,*), nub --> xmom_cm(10,*), bb --> xmom_cm(11,*), 
c consistently with the labelling conventions used in MC@NLO:
c   x(1)y(2) -> z(3)t(4)[->l+(6)nu(7)b(8)]tb(5)[->l-(9)nub(10)bb(11)]
c The inputs of the routine are cth1,phi1,cth2,phi2, which are cosines of
c polar angles and azimuthal angles, with
c   (cth1,phi1) --> direction of W in the top/tbar rest frame
c   (cth2,phi2) --> direction of l in the W rest frame
      implicit none
      real*8 xmt,xmw,cth1,phi1,cth2,phi2,xtq(4),xbq(4),xel(4),xnu(4)
      integer iqrk
      real*8 tiny,xmt2,xmw2,sth1,sth2,ew,eb,pwx,pwy,pwz,pbx,pby,pbz,
     # eel,enu,pex,pey,pez,pnx,pny,pnz,xmq,tmp(5),tmp1(4),tmp2(4)
      parameter (tiny=1.d-8)
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      integer itop,iel,inu,ib
c
      if(iqrk.eq.1)then
        itop=4
        iel=6
        inu=7
        ib=8
      elseif(iqrk.eq.2)then
        itop=5
        iel=9
        inu=10
        ib=11
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
      real*8 xmom_cm(11,4)
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
c
c
c End of event-generation routines
c
c
c
c
c Begin of phase-space routines
c
c
      subroutine invar(xm2,s,x,y,cth1,cth2,str,
     #     tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
c This routine upgrades the one used up to version 3.2, which was derived
c from MNR. Only str='p1' is a legal choice, which coincides with the 
c formulae of the MC@NLO paper. The dimensions of xmom_cm have been changed
c   xmom_cm(11,4) --> xmom_cm(5,4)
c to include the four-momenta of the leptons and b quarks emerging from
c top decays. The calls to the routine fillmom have been eliminated, to
c avoid potential (although very unlikely) numerical problems, in favour
c of the formulae adopted for VH production.
c This routine is meant to be called to generate the kinematics of 
c direct events, i.e. those for which
c   a(p1)+b(p2) --> Q(k1)+Qbar(k2)+c(k)
c The outputs are given in the common blocks perpen (pq1(2)==Q, 
c pq2(2)==Qbar, pp(2)==c transverse momenta), ycmvar (yq1==Q,
c yq2==Qbar ,yp==c rapidities in the parton c.m. frame), and stored
c as four-momenta in xmom_cm(ipart,icomp), with the conventions:
c   icomp=1 -> px, icomp=2 -> py, icomp=3 -> pz, icomp=4 -> E;
c   ipart=1 -> p1, ipart=2 -> p2, ipart=3 -> k, ipart=4 -> k1, ipart=5 -> k2
c The four-momenta of the decay products are computed elsewhere 
c (see the routine gentopdmom())
      implicit none
      real * 8 xm2,s,x,y,cth1,cth2,tk,uk,q1q,q2q,q1c,q2c,
     # w1,w2,w1h,w2h
      character * 2 str
      real * 8 pq1,pq2,pp,yq1,yq2,yp
      common/perpen/pq1(2),pq2(2),pp(2)
      common/ycmvar/yq1,yq2,yp
      real * 8 tiny,s2,drs2,p10,p20,k0,k10,k20,bx,sth1,cpsi,
     # spsi,cpsi2,spsi2,cpsi1,spsi1,xktsq,xkt1sq,xkt2sq,
     # xkt,xkt1,xkt2,tmp,sqs,e1lab,pl1lab,e2lab,pl2lab
      parameter (tiny=1.d-14)
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      integer ichkmom
      common/cichkmom/ichkmom
c
      tk  = -s/2*(1-x)*(1-y)
      uk  = -s/2*(1-x)*(1+y)
      s2  = tk+uk+s
      drs2 = 2*sqrt(s2)
      p10 = (s+tk)/drs2
      p20 = (s+uk)/drs2
      k0  = -(tk+uk)/drs2
      k10 = drs2/4
      k20 = drs2/4
      bx = sqrt(1-4*xm2/s2)
      sth1 = sqrt(1-cth1**2)
      if(str.eq.'p1') then
         cpsi2 = 1
         spsi2 = 0
         cpsi = 1-8*x/((1+y+x*(1-y))*(1-y+x*(1+y)))
         spsi = 4*(1-x)*sqrt(x*(1-y**2))/
     #          ((1+y+x*(1-y))*(1-y+x*(1+y)))
         cpsi1 = (1+y-x*(1-y))/(1+y+x*(1-y))
         spsi1 = sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
      else
         write(6,*) 'error in invar: str=',str
         stop
      endif
      q1q = - 2*p10*k10*(1-bx*(cth2*sth1*spsi2+cth1*cpsi2))
      q2q = - 2*p20*k20*(1+bx*(cth2*sth1*spsi +cth1*cpsi ))
      q1c = - s - tk - q1q
      q2c = - s - uk - q2q
      w1  = - q1q + q2q - tk
      w2  = - q2q + q1q - uk
      w1h = 1-bx*(cth2*sth1*spsi1+cth1*cpsi1)
      w2h = 1+bx*(cth2*sth1*spsi1+cth1*cpsi1)
c
      if(abs(q1q).lt.tiny) then
        yq1  = 1.d8
      elseif(abs(q2c).lt.tiny) then
        yq1  = -1.d8
      else
        yq1 = .5d0*log( q2c/q1q )
      endif
      if(abs(q1c).lt.tiny) then
        yq2  = 1.d8
      elseif(abs(q2q).lt.tiny) then
        yq2  = -1.d8
      else
        yq2 = .5d0*log( q2q/q1c )
      endif
      if(tk.eq.0) then
         yp  = 1.d8
      elseif(uk.eq.0) then
         yp  = -1.d8
      else
         yp  = .5d0*log( uk/tk )
      endif
c-----------------------------------------------------------------
c xktsq, xkt1sq e xkt2sq are the square of transverse momenta of g, Q, 
c and Qb respectively. The axis orientation is such that Q is always
c along the x direction
c
      xktsq  = uk*tk/s
      if(xktsq.eq.0) then
         pq1(1) = sqrt(x*s)/2.d0*bx*sth1
         pq1(2) = 0.d0
         pq2(1) = -pq1(1)
         pq2(2) = 0.d0
         pp(1) = 0.d0
         pp(2) = 0.d0
      else
         xkt1sq = q2c*q1q/s - xm2
         xkt2sq = q2q*q1c/s - xm2
         xkt = sqrt(xktsq)
         xkt1 = sqrt(xkt1sq)
         xkt2 = sqrt(xkt2sq)
         pq1(1) = xkt1
         pq1(2) = 0.d0
         pq2(1) = (xktsq-xkt1sq-xkt2sq)/(2.d0*xkt1)
         tmp = xkt2sq-pq2(1)**2
         if(tmp.gt.0.d0)then
            pq2(2) = sqrt(tmp)
         else
            pq2(2) = 0.d0
         endif
         pp(1) = (xkt2sq-xkt1sq-xktsq)/(2.d0*xkt1)
         tmp = xktsq-pp(1)**2
         if(tmp.gt.0.d0)then
            pp(2) = -sqrt(tmp)
         else
            pp(2) = 0.d0
         endif
      endif
c Incoming parton
      sqs=sqrt(s)
      xmom_cm(1,1)=0.d0
      xmom_cm(1,2)=0.d0
      xmom_cm(1,3)=sqs/2.d0
      xmom_cm(1,4)=sqs/2.d0
      xmom_cm(2,1)=0.d0
      xmom_cm(2,2)=0.d0
      xmom_cm(2,3)=-sqs/2.d0
      xmom_cm(2,4)=sqs/2.d0
c Outgoing light parton
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
        xmom_cm(3,1)=pp(1)
        xmom_cm(3,2)=pp(2)
        xmom_cm(3,3)=sqs/2.d0*(1-x)*y
        xmom_cm(3,4)=sqs/2.d0*(1-x)
      endif
c Heavy quark
      e1lab=-(q1q+q2c)/(2*sqs)
      pl1lab=(q1q-q2c)/(2*sqs)
      xmom_cm(4,1)=pq1(1)
      xmom_cm(4,2)=pq1(2)
      xmom_cm(4,3)=pl1lab
      xmom_cm(4,4)=e1lab
c Heavy antiquark
      e2lab=-(q1c+q2q)/(2*sqs)
      pl2lab=(q1c-q2q)/(2*sqs)
      xmom_cm(5,1)=pq2(1)
      xmom_cm(5,2)=pq2(2)
      xmom_cm(5,3)=pl2lab
      xmom_cm(5,4)=e2lab
c
      if(ichkmom.eq.0)call checkmom(xmom_cm,s,0.d0,1,2)
      return
      end


      subroutine invar2(xm12,xm22,s,x,y,cth1,cth2,str,
     #            tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
c This routine is derived from invar() of the VV package. It serves to
c fill xmom_cm in the case of top and tbar unequal masses, and must be
c called only by getspincorr(). It must not be used to compute MNR matrix
c elements. The common blocks have been redefined to be consistent with
c those of invar() of the QQ package. The return values of the invariants
c which are the entries of this routine will be set to zero to force the
c code to crash if an improper use is attempted; exception is made for
c tk and uk, since these are used to define the FKS damping factor
      implicit none
      real * 8 xm12,xm22,s,x,y,cth1,cth2,tk,uk,q1q,q2q,q1c,q2c,
     # w1,w2,w1h,w2h
      character * 2 str
      real * 8 pq1,pq2,pp,yq1,yq2,yp
      common/perpen/pq1(2),pq2(2),pp(2)
      common/ycmvar/yq1,yq2,yp
      real * 8 tiny,s2,drs2,p10,p20,k0,k10,k20,bx,sth1,cpsi2,spsi2,
     # cpsi,spsi,cpsi1,spsi1,xktsq,xkt1,xkt2,xkt1sq,xkt2sq,xkt,
     # tmp,sqs,e1lab,pl1lab,e2lab,pl2lab
      parameter (tiny=1.d-14)
      real*8 xmom_cm(11,4)
      common/cxmomcm/xmom_cm
      integer ichkmom
      common/cichkmom/ichkmom
c
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
c
      if(abs(q1q-xm12).lt.tiny) then
        yq1  = 1.d8
      elseif(abs(q2c-xm12).lt.tiny) then
        yq1  = -1.d8
      else
        yq1 = .5d0*log( (xm12-q2c)/(xm12-q1q) )
      endif
      if(abs(q1c-xm22).lt.tiny) then
        yq2  = 1.d8
      elseif(abs(q2q-xm22).lt.tiny) then
        yq2  = -1.d8
      else
        yq2 = .5d0*log( (xm22-q2q)/(xm22-q1c) )
      endif
      if(abs(tk).lt.tiny) then
        yp  = 1.d8
      elseif(abs(uk).lt.tiny) then
        yp  = -1.d8
      else
        yp  = .5d0*log( uk/tk )
      endif
c-----------------------------------------------------------------
c xktsq, xkt1sq e xkt2sq are the square of transverse momenta of g, Q, 
c and Qb respectively. The axis orientation is such that Q is always
c along the x direction. The component of p_T(Qb) along the y direction
c is always positive or zero
c
      xktsq = uk*tk/s
      if(xktsq.eq.0) then
         pq1(1) = bx*sth1
         pq1(2) = 0.d0
         pq2(1) = -pq1(1)
         pq2(2) = 0.d0
         pp(1) = 0.d0
         pp(2) = 0.d0
         xkt1 = pq1(1)
         xkt2 = xkt1
      else
         xkt1sq = (xm12-q2c)*(xm12-q1q)/s - xm12
         xkt2sq = (xm22-q2q)*(xm22-q1c)/s - xm22
         xkt = sqrt(xktsq)
         xkt1 = sqrt(xkt1sq)
         xkt2 = sqrt(xkt2sq)
         pq1(1) = xkt1
         pq1(2) = 0.d0
         pq2(1) = (xktsq-xkt1sq-xkt2sq)/(2.d0*xkt1)
         tmp = xkt2sq-pq2(1)**2
         if(tmp.gt.0.d0)then
            pq2(2) = sqrt(tmp)
         else
            pq2(2) = 0.d0
         endif
         pp(1) = (xkt2sq-xkt1sq-xktsq)/(2.d0*xkt1)
         tmp = xktsq-pp(1)**2
         if(tmp.gt.0.d0)then
            pp(2) = -sqrt(tmp)
         else
            pp(2) = 0.d0
         endif
      endif
      if(ichkmom.eq.0)call checkptcon(pq1,pq2,pp)
c
c Incoming partons
      sqs=sqrt(s)
      xmom_cm(1,1)=0.d0
      xmom_cm(1,2)=0.d0
      xmom_cm(1,3)=sqs/2.d0
      xmom_cm(1,4)=sqs/2.d0
      xmom_cm(2,1)=0.d0
      xmom_cm(2,2)=0.d0
      xmom_cm(2,3)=-sqs/2.d0
      xmom_cm(2,4)=sqs/2.d0
c Outgoing light parton
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
        xmom_cm(3,1)=pp(1)
        xmom_cm(3,2)=pp(2)
        xmom_cm(3,3)=sqs/2.d0*(1-x)*y
        xmom_cm(3,4)=sqs/2.d0*(1-x)
      endif
c Heavy quark
      e1lab=(2*xm12-q1q-q2c)/(2*sqs)
      pl1lab=(q1q-q2c)/(2*sqs)
      xmom_cm(4,1)=pq1(1)
      xmom_cm(4,2)=pq1(2)
      xmom_cm(4,3)=pl1lab
      xmom_cm(4,4)=e1lab
c Heavy antiquark
      e2lab=(2*xm22-q1c-q2q)/(2*sqs)
      pl2lab=(q1c-q2q)/(2*sqs)
      xmom_cm(5,1)=pq2(1)
      xmom_cm(5,2)=pq2(2)
      xmom_cm(5,3)=pl2lab
      xmom_cm(5,4)=e2lab
c
      if(ichkmom.eq.0)call checkmom(xmom_cm,s,0.d0,1,2)
c Set all return values to zero, since invariants must not be used 
c after calling this routine. This may be easily changed. However,
c the entries of this routine follow MNR conventions, while the local
c values have been computed using FNR conventions. Thus, the assignments
c below need be changed into
c  q1q=q1q-xm12;  q2q=q2q-xm22;  q1c=q1c-xm22;  q2c=q2c-xm12;
c  w1=w1-xm12;  w2=w2-xm22
c and w1h, w2h need be recomputed explicitly
      q1q=0.d0
      q2q=0.d0
      q1c=0.d0
      q2c=0.d0
      w1=0.d0
      w2=0.d0
      w1h=0.d0
      w2h=0.d0
c
      return
      end


      subroutine checkmom(xmom,smax,ybst,iflag,itype)
      implicit none
      real * 8 xmom(11,4)
      real * 8 smax,ybst,xpmax
      real*8 x1,x2
      common/cx1x2/x1,x2
      real * 8 tiny,vtiny,xsum(4),xsuma(4),xsign,xrat(4)
      parameter (tiny=5.d-3)
      parameter (vtiny=1.d-5)
      integer iflag,itype,i,j,jj,jflag,jeflag,jmax
c
      if(itype.eq.1)then
        jmax=11
      elseif(itype.eq.2)then
        jmax=5
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
          if((itype.eq.1.and.j.ne.4.and.j.ne.5).or.itype.eq.2)then
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
          do j=1,11
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
        do j=1,11
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
      real * 8 xmom(11,4)
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


      subroutine checkptcon(ptv1,ptv2,ptvg)
      implicit none
      real*8 ptv1(2),ptv2(2),ptvg(2),tiny,pt1,pt2,ptmax
      parameter (tiny=1.d-5)
      integer jj
c
      ptmax=max(abs(ptv1(1)),abs(ptv2(1)),abs(ptvg(1)),
     #          abs(ptv1(2)),abs(ptv2(2)),abs(ptvg(2)))
      pt1=ptv1(1)+ptv2(1)+ptvg(1)
      pt2=ptv1(2)+ptv2(2)+ptvg(2)
      if(pt1.gt.ptmax*tiny.or.pt2.gt.ptmax*tiny)then
        write(*,*)'Transverse momentum is not conserved'
        write(*,'(4(d14.8,1x))') (ptv1(jj),jj=1,2)
        write(*,'(4(d14.8,1x))') (ptv2(jj),jj=1,2)
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


      function bwfunc_mod(s,xm02,ga0,gamn,gapl)
c Returns a skewed Breit Wigner function, defined as follows
c   skBW(M)=BW(M,M0,G0,G_mn)  if  M<M0
c   skBW(M)=BW(M,M0,G0,G_pl)  if  M>M0
c where 
c   BW(M,M0,G0,G)=(M0*G^2)/(pi*G0) 1/((M^2-M0^2)^2+M0^2 G^2)
c whose normalization is therefore such that
c   BW(M0,M0,G0,G)=1/(pi*M0*G0)
      implicit none
      real*8 bwfunc_mod,s,xm02,ga0,gamn,gapl
      real*8 pi,xm0,tmp
      parameter (pi=3.1415926535897932d0)
c
      xm0=sqrt(xm02)
      if(s.le.xm02)then
        tmp=xm0*gamn**2/(ga0*pi*((s-xm02)**2+xm02*gamn**2))
      else
        tmp=xm0*gapl**2/(ga0*pi*((s-xm02)**2+xm02*gapl**2))
      endif
      bwfunc_mod=tmp
      return
      end


      function xbwmass3_mod(t,xm02,ga0,gamn,gapl,bwdelf,bwfmmn)
c Returns the boson mass squared, given 0<t<1, the nominal mass (xm0),
c the widths of the skewed Breit Wigner function, and the mass range 
c (implicit in bwdelf and bwfmmn). This function is the inverse of 
c F(M^2), where
c   F(M^2)=\int_{xmlow2}^{M^2} ds skBW(sqrt(s))
c   skBW(M)=BW(M,M0,G0,G_mn)  if  M<M0
c   skBW(M)=BW(M,M0,G0,G_pl)  if  M>M0
c with
c   BW(M,M0,G0,G)=(M0*G^2)/(pi*G0) 1/((M^2-M0^2)^2+M0^2 G^2)
c whose normalization is therefore such that
c   BW(M0,M0,G0,G)=1/(pi*M0*G0)
c and therefore eats up the skewed Breit-Wigner when changing integration 
c variable M^2 --> t
      implicit none
      real*8 xbwmass3_mod,t,xm02,ga0,gamn,gapl,bwdelf,bwfmmn
      real*8 pi,xm0,tmp
      parameter (pi=3.1415926535897932d0)
c
      xm0=sqrt(xm02)
      if(t.le.(bwfmmn/pi))then
        tmp=xm02+xm0*gamn*tan( ga0/gamn*(pi*bwdelf*t-bwfmmn) )
      else
        tmp=xm02+xm0*gapl*tan( ga0/gapl*(pi*bwdelf*t-bwfmmn) )
      endif
      xbwmass3_mod=tmp
      return
      end


      subroutine chvar(z,x,xjac,del)
c
      implicit none
      real * 8 z,x,xjac,del,odel
      real * 8 xld,x0,y0,yh,yl,ydif,y,t
      data odel/0.d0/
      if(odel.ne.del) then
         odel = del
         xld = log(del)
         x0  = 1 + del + 1/xld
         if(x0.le.0.or.x0.ge.1) then
            call hvqwarn('CHVAR')
            write(*,*) 'inappropriate delta'
            stop
         endif
         y0  = log(1-x0+del)
         yh = y0 - x0*xld
         yl = xld
         ydif = yh-yl
      endif
      y = yl + z*ydif
      xjac = xjac * ydif
      if(y.gt.y0) then
          x = - (yh - y)/xld
          xjac = xjac / (-xld)
      else
          t = exp(y)
          xjac = xjac * t
          x = 1 - t + del
      endif
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
      REAL*8 MINV,FK88RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      FK88RANDOM = SEED*MINV
      END


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
c Initialization
c
c
      subroutine setpar()
c Electroweak parameters for top decay; this routine is a modification
c of that of the VB code
      implicit none
      include 'hvqcblks.h'
      integer i,j
      real * 8 pi,zero,one,xme,xmmu,xmtau
      parameter (pi=3.14159265358979312D0)
      parameter (one=1.d0)
      parameter (zero=0.d0)
c Values from PDG 2003
      parameter (xme=0.510998902d-3)
      parameter (xmmu=105.6583568d-3)
      parameter (xmtau=1776.99d-3)
      real * 8 ckm(1:6,1:6),ckm2(1:6,1:6),vickm(1:6,1:6)
      common/cckm2/ckm2
      common/cvickm/vickm
      real * 8 ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      common/cckmfct/ruckm,rcckm,rtckm,rducckm,rsucckm,rbucckm
      real * 8 xmw,gaw
      common/cwparam/xmw,gaw
      real * 8 xlep1mass(2),xlep2mass(2)
      common/clepmass/xlep1mass,xlep2mass
      real * 8 pdglepmss(11:16)
      common/cpdglepmss/pdglepmss
      real * 8 xm1low2,xm1upp2,xm2low2,xm2upp2
      common/bounds/xm1low2,xm1upp2,xm2low2,xm2upp2
      real * 8 brrtop1,brrtop2
      common/brratios/brrtop1,brrtop2
      real * 8 xbrrtoplep,xbrrtophad
      common/xibrratios/xbrrtoplep,xbrrtophad
      real * 8 frac12,frac123
      common/cfracs123/frac12,frac123
      real * 8 sthw2,cthw2
      common/cweinan/sthw2,cthw2
      real*8 sthw20
      common/csthw20/sthw20
      real * 8 xmt,twidth
      common/ctparam/xmt,twidth
      integer idec
      common/cidec/idec
      integer iwidth
      common/ciwidth/iwidth
      integer inonbtop
      common/cinonbtop/inonbtop
c Type of V decays, with HERWIG conventions; see the beginning of this file
      integer il1hw,il2hw
      common/cilhw/il1hw,il2hw
c Identities of the vector bosons or leptons in the final state 
c (PDG conventions)
      integer ip4s,ip5s,ip6s,ip7s,ip8s,ip9s
      common/ci2parts/ip4s,ip5s,ip6s,ip7s,ip8s,ip9s
c PDG codes for charged leptons and neutrinos for a given IL (NLO) code;
c the particle code (not the antiparticle) is entered here
c Charged lepton from W decay
      integer ichlw(0:6)
      data ichlw/0,11,13,15,0,0,0/
c Neutrino from W decay
      integer ineuw(0:6)
      data ineuw/0,12,14,16,0,0,0/
      real * 8 alfaem,xalfaem,xmwme,gawme,topdecw,xmt2,xmw2,tmpmss(3)
c
c Electron charge squared (computed at the W mass)
      xmt2=xmt**2
      xmw2=xmw**2
      alfaem = xalfaem(xmw2)
      ze2 = 4*pi*alfaem
c sin and cos squared of theta_W; MSbar scheme, from PDG2003
      sthw2=0.23113d0
      cthw2=1-sthw2
      if(sthw20.ne.sthw2)then
        write(*,*)'Inconsistency between setpar and reset_twdbr'
        stop
      endif
c Lepton masses and identities: xlep#mass(i) is the mass of the decay 
c product # in the decay of top (i=1, in which case #=1 ==> antiparticle,
c #=2 ==> particle) or tbar (i=2, in which case #=1 ==> particle,
c #=2 ==> antiparticle). If the tops have a given leptonic decay, the
c masses are set here, but are in any case overwritten by the routine
c getpdecids(), in order to treat all the cases on equal footing
      pdglepmss(11)=xme
      pdglepmss(12)=0.d0
      pdglepmss(13)=xmmu
      pdglepmss(14)=0.d0
      pdglepmss(15)=xmtau
      pdglepmss(16)=0.d0
      tmpmss(1)=xme
      tmpmss(2)=xmmu
      tmpmss(3)=xmtau
      ip4s=-ichlw(il1hw)
      ip5s=ineuw(il1hw)
      ip6s=5
      ip7s=ichlw(il2hw)
      ip8s=-ineuw(il2hw)
      ip9s=-5
      if(il1hw.ge.1.and.il1hw.le.3)then
        xlep1mass(1)=tmpmss(il1hw)
        xlep2mass(1)=0.d0
      elseif(il1hw.eq.0.or.(il1hw.ge.4.and.il1hw.le.6))then
        xlep1mass(1)=-1.d0
        xlep2mass(1)=-1.d0
      else
        write(*,*)'Error in setpar: inconsistent entries'
        stop
      endif
      if(il2hw.ge.1.and.il2hw.le.3)then
        xlep1mass(2)=tmpmss(il2hw)
        xlep2mass(2)=0.d0
      elseif(il2hw.eq.0.or.(il2hw.ge.4.and.il2hw.le.6))then
        xlep1mass(2)=-1.d0
        xlep2mass(2)=-1.d0
      else
        write(*,*)'Error in setpar: inconsistent entries'
        stop
      endif
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
c Combinations used for unweighting
      ruckm=ckm2(1,2)+ckm2(1,3)+ckm2(1,5)
      rcckm=ckm2(4,2)+ckm2(4,3)+ckm2(4,5)
      rtckm=ckm2(6,2)+ckm2(6,3)+ckm2(6,5)
      rducckm=ckm2(1,2)+ckm2(4,2)
      rsucckm=ckm2(1,3)+ckm2(4,3)
      rbucckm=ckm2(1,5)+ckm2(4,5)
c Fills MadEvent common blocks. Set positron charge and QCD coupling g 
c equal to one, and use the actual values in the main code
      xmwme=xmw
      gawme=gaw
      call setmepar(xmwme,gawme,zero,zero,
     #              xmt,twidth,zero,sthw2,one,one)
c Compute reweight factors for event weights, that take into
c account the values of branching ratios
      if(iwidth.eq.1)then
c Correction factors for effective W mass ranges
        brrtop1=topdecw(xmt,xmw,gaw,xm1low2,xm1upp2,sthw2)/
     #          topdecw(xmt,xmw,gaw,zero,xmt2,sthw2)
        brrtop2=topdecw(xmt,xmw,gaw,xm2low2,xm2upp2,sthw2)/
     #          topdecw(xmt,xmw,gaw,zero,xmt2,sthw2)
      else
c W's are on shell: the event weight is proportional to dBR_t/dQ_W^2 
c (at Q_W^2=MW^2); use dGamma_t/dQ_W^2 = Gamma_t/(pi*MW*Gamma_W)
        brrtop1=1.d0/(pi*xmw*gaw)
        brrtop2=1.d0/(pi*xmw*gaw)
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
      if(il2hw.eq.0)then
        brrtop2=(3*xbrrtoplep+2*xbrrtophad) * brrtop2
      elseif(il2hw.ge.1.and.il2hw.le.3)then
        brrtop2=xbrrtoplep * brrtop2
      elseif(il2hw.eq.4)then
        brrtop2=2*xbrrtoplep * brrtop2
      elseif(il2hw.eq.5)then
        brrtop2=2*xbrrtophad * brrtop2
      elseif(il2hw.eq.6)then
        brrtop2=(2*xbrrtoplep+2*xbrrtophad) * brrtop2
      endif
c Correct reweight factors for top decays if only t->Wb decays
c are selected
      if(inonbtop.eq.0)then
        brrtop1=ckm2(6,5)/rtckm * brrtop1
        brrtop2=ckm2(6,5)/rtckm * brrtop2
      elseif(inonbtop.ne.1)then
        write(*,*)'Unknown inonbtop in setpar',inonbtop
        stop
      endif
c frac12 is the fraction of decays W->e+mu/W->e+mu+all quarks
      frac12=2*xbrrtoplep/(2*xbrrtoplep+2*xbrrtophad)
c frac123 is the fraction of decays W->e+mu+tau/W->e+mu+tau+all quarks
      frac123=3*xbrrtoplep/(3*xbrrtoplep+2*xbrrtophad)
c
      if(brrtop1.eq.0.d0.or.brrtop2.eq.0.d0)then
        write(*,*)
     #    'These decay channels will return a zero cross section'
        stop
      endif
c
      return
      end


      subroutine reset_twdbr(xmt,twidth,xbrrtoplep,xbrrtophad,
     #                       gammax1,gammax2)
c If twidth>0, use xbrrtoplep and xbrrtophad as branching ratios.
c If twidth<0, compute the total width according to LO SM, and
c redefine xbrrtoplep xbrrtophad accordingly.
c The width and branching ratios assume lepton and hadron universality;
c no CKM factors are included here
      implicit none
      real*8 xmt,twidth,xbrrtoplep,xbrrtophad,gammax1,gammax2
      real*8 pi,zero,xmt2,xmw2,alfaem,xalfaem,ze2,brtop2,topdecw,
     # brtopw,brtop0
      parameter (pi=3.14159265358979312D0)
      parameter (zero=0.d0)
      real*8 xmw,gaw
      common/cwparam/xmw,gaw
      real*8 sthw20
      common/csthw20/sthw20
c
c sthw20 must have the same value as sthw2 in setpar. If this is not
c the case, the code stops in setpar
      sthw20=0.23113d0
      if(twidth.gt.0.d0)then
        if( xbrrtoplep.lt.0.d0.or.xbrrtophad.lt.0.d0 .or.
     #      (xbrrtoplep.eq.0.d0.and.xbrrtophad.eq.0.d0) )then
          write(*,*)'Error in reset_twdbr: negative or zero BRs',
     #              xbrrtoplep,xbrrtophad
          stop
        endif
      elseif(twidth.lt.0.d0)then
        xmt2=xmt**2
        xmw2=xmw**2
        alfaem=xalfaem(xmw2)
        ze2=4*pi*alfaem
c
        xbrrtoplep=1/9.d0
        xbrrtophad=1/3.d0
        if(gammax1.ne.0.and.gammax2.ne.0)then
c topdecw returns dGamma(t->blnu)/dq^2 integrated over q2 in the range
c (xm1low2,xm1upp2), up to a factor e^4*|Vtb|^2
          brtop2=topdecw(xmt,xmw,gaw,zero,xmt2,sthw20)
          twidth=brtop2*ze2**2*9.d0
        elseif(gammax1.eq.0.and.gammax2.eq.0)then
c W's are at the pole mass
c brtopw is dGamma(t->blnu)/dq^2, up to a factor gw^4*|Vtb|^2/2
          brtopw=(xmt2-xmw2)**2*(xmt2+2*xmw2)/(6144*pi**3*xmt**3)
c Multiply brtopw by the Breit Wigner for the W, in the narrow width
c approximation BR(W)->pi/(MW*GammaW). A factor 2/sin^4\theta_W is inserted
c in order to eliminate the factor 1/2 in brtopw, and convert gw^4 into e^4
          brtop0=brtopw* 2*pi/(sthw20**2 * xmw*gaw)
          twidth=brtop0*ze2**2*9.d0
        else
          write(*,*)'Error #1 in reset_twdbr',gammax1,gammax2
          stop
        endif
      else
        write(*,*)'Error in reset_twdbr: top width equal to zero'
        stop
      endif
      if( (3*xbrrtoplep+2*xbrrtophad).gt.1.0001d0)then
        write(*,*)'Error #2 in reset_twdbr',xbrrtoplep,xbrrtophad,
     #            3*xbrrtoplep+2*xbrrtophad
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


      subroutine parsetpar()
      implicit none
      integer jproc,ileg,ie0sq,i,itype
      integer imapp(0:5),invimapp(-5:21)
      integer ialwsplit(1:3,1:4,1:3)
      integer icllborn(1:3,1:4,1:3),icllkern(1:3,1:4,1:3)
      integer icolconn(1:3,1:4,1:3)
      integer ivbhpro(4,3,5)
      integer idp1(4,3,5),idp2(4,3,5),idp3(4,3,5)
      character * 2 xproc(3)
      common/cimapp/imapp
      common/cinvimapp/invimapp
      common/cialwsplit/ialwsplit
      common/cicllsplit/icllborn,icllkern
      common/cicolconn/icolconn
      common/civbhpro/ivbhpro
      common/cidpart/idp1,idp2,idp3
      common/cxproc/xproc
c
      xproc(1)='gg'
      xproc(2)='qq'
      xproc(3)='qg'
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
c invimapp(i) is the inverse map of imapp(i)
      do i=-5,21
        invimapp(i)=-100000
      enddo
      invimapp(1)=2
      invimapp(2)=1
      invimapp(3)=3
      invimapp(4)=4
      invimapp(5)=5
      invimapp(21)=0
      do i=-5,-1
        invimapp(i)=-invimapp(-i)
      enddo
c
c ialwsplit(jproc,ileg,ie0sq) returns 0 if there's no splitting from leg ileg
c in the 2-->3 process whose initial state is identified by jproc (jproc=1,2,3
c ==> prc=gg,qq,qg respectively), and whose colour connection corresponds to
c the scale choice ie0sq (ie0sq=1,2,3 ==> e0sq=s/2,-t/2,-u/2)
      do jproc=1,3
        do ileg=1,4
          do ie0sq=1,3
            ialwsplit(jproc,ileg,ie0sq)=0
          enddo
        enddo
      enddo
c jproc=1 (gg)
      do ileg=1,4
        do ie0sq=1,3
          if(ileg.le.2.or.(ileg.gt.2.and.ie0sq.ne.1))
     #      ialwsplit(1,ileg,ie0sq)=1
        enddo
      enddo
c jproc=2 (qq)
      do ileg=1,4
        ialwsplit(2,ileg,2)=1
      enddo
c jproc=3 (qg)
      do ie0sq=1,3
        ialwsplit(3,1,ie0sq)=1
      enddo
      ialwsplit(3,2,2)=1
c
c icllborn(jproc,ileg,ie0sq) returns the code (1==gg, 2==qq) for the initial
c state of the Born that factorizes in the MC subtraction terms
      do jproc=1,3
        do ileg=1,4
          do ie0sq=1,3
            icllborn(jproc,ileg,ie0sq)=0
          enddo
        enddo
      enddo
c Final-state emissions
      do jproc=1,3
        do ileg=3,4
          do ie0sq=1,3
            if(ialwsplit(jproc,ileg,ie0sq).eq.1)
     #        icllborn(jproc,ileg,ie0sq)=jproc
          enddo
        enddo
      enddo
c jproc=1 (gg)
      do ileg=1,2
        do ie0sq=1,3
          if(ialwsplit(1,ileg,ie0sq).eq.1)
     #     icllborn(1,ileg,ie0sq)=1
        enddo
      enddo
c jproc=2 (qq)
      do ileg=1,2
        do ie0sq=1,3
          if(ialwsplit(2,ileg,ie0sq).eq.1)
     #     icllborn(2,ileg,ie0sq)=2
        enddo
      enddo
c jproc=3 (qg)
      do ie0sq=1,3
        icllborn(3,1,ie0sq)=1
      enddo
      icllborn(3,2,2)=2
c
c icllkern(jproc,ileg,ie0sq) returns the code (1,..,4, see the function
c ap_kern) of the Altarelli-Parisi splitting function that factorizes in 
c the MC subtraction terms
      do jproc=1,3
        do ileg=1,4
          do ie0sq=1,3
            icllkern(jproc,ileg,ie0sq)=0
          enddo
        enddo
      enddo
c Final-state emissions
      do jproc=1,3
        do ileg=3,4
          do ie0sq=1,3
            if(ialwsplit(jproc,ileg,ie0sq).eq.1)
     #        icllkern(jproc,ileg,ie0sq)=4
          enddo
        enddo
      enddo
c jproc=1 (gg)
      do ileg=1,2
        do ie0sq=1,3
          if(ialwsplit(1,ileg,ie0sq).eq.1)
     #     icllkern(1,ileg,ie0sq)=1
        enddo
      enddo
c jproc=2 (qq)
      do ileg=1,2
        do ie0sq=1,3
          if(ialwsplit(2,ileg,ie0sq).eq.1)
     #     icllkern(2,ileg,ie0sq)=4
        enddo
      enddo
c jproc=3 (qg)
      do ie0sq=1,3
        icllkern(3,1,ie0sq)=3
      enddo
      icllkern(3,2,2)=2
c
c If icolconn(jproc,ileg,ie0sq)<0, the corresponding 2-->2 reduced matrix
c element squared (d\bar{sigma}) is multiplied by 1/|icolconn|; if it's 
c positive, it is also multiplied by u**2/(u**2+t**2) or t**2/(u**2+t**2), 
c depending on the colour structure chosen. See hvqborncol
      do jproc=1,3
        do ileg=1,4
          do ie0sq=1,3
            icolconn(jproc,ileg,ie0sq)=0
          enddo
        enddo
      enddo
c jproc=1 (gg)
      icolconn(1,1,1)=-2
      icolconn(1,1,2)=2
      icolconn(1,1,3)=2
      icolconn(1,2,1)=-2
      icolconn(1,2,2)=2
      icolconn(1,2,3)=2
      do ileg=3,4
        do ie0sq=1,3
          if(ialwsplit(1,ileg,ie0sq).eq.1)
     #     icolconn(1,ileg,ie0sq)=1
        enddo
      enddo
c jproc=2 (qq)
      do ileg=1,4
        do ie0sq=1,3
          if(ialwsplit(2,ileg,ie0sq).eq.1)
     #     icolconn(2,ileg,ie0sq)=-1
        enddo
      enddo
c jproc=3 (qg)
      icolconn(3,1,1)=-2
      icolconn(3,1,2)=2
      icolconn(3,1,3)=2
      icolconn(3,2,2)=-1
c
c ivbhpro returns the process number associated to the entries; this is
c identical to i1hpro (see the routine store_events)
      do i=1,4
        do jproc=1,3
          do itype=1,5
            ivbhpro(i,jproc,itype)=0
          enddo
        enddo
      enddo
c 
      do i=1,4
        ivbhpro(i,1,1)=407
      enddo
c 
      do itype=1,5
        ivbhpro(1,2,itype)=401
        ivbhpro(2,2,itype)=403
        ivbhpro(3,2,itype)=403
        ivbhpro(4,2,itype)=401
        ivbhpro(1,3,itype)=402
        ivbhpro(2,3,itype)=404
        ivbhpro(3,3,itype)=405
        ivbhpro(4,3,itype)=406
      enddo
c
c idpX returns the flavour of parton number X (1=coming from the left,
c 2=coming from the right, 3=outgoing) in the process associated to the
c entries. The labelling scheme of PDG has been used
      do i=1,4
        do jproc=1,3
          do itype=1,5
            idp1(i,jproc,itype)=0
            idp2(i,jproc,itype)=0
            idp3(i,jproc,itype)=0
          enddo
        enddo
      enddo
c
      do i=1,4
        idp1(i,1,1)=21
        idp2(i,1,1)=21
        idp3(i,1,1)=21
      enddo
c
      do itype=1,5
        idp1(1,2,itype)=imapp(itype)
        idp1(2,2,itype)=-imapp(itype)
        idp1(3,2,itype)=-imapp(itype)
        idp1(4,2,itype)=imapp(itype)
c
        idp1(1,3,itype)=imapp(itype)
        idp1(2,3,itype)=-imapp(itype)
        idp1(3,3,itype)=21
        idp1(4,3,itype)=21
c
        idp2(1,2,itype)=-imapp(itype)
        idp2(2,2,itype)=imapp(itype)
        idp2(3,2,itype)=imapp(itype)
        idp2(4,2,itype)=-imapp(itype)
c
        idp2(1,3,itype)=21
        idp2(2,3,itype)=21
        idp2(3,3,itype)=imapp(itype)
        idp2(4,3,itype)=-imapp(itype)
c
        idp3(1,2,itype)=21
        idp3(2,2,itype)=21
        idp3(3,2,itype)=21
        idp3(4,2,itype)=21
c
        idp3(1,3,itype)=imapp(itype)
        idp3(2,3,itype)=-imapp(itype)
        idp3(3,3,itype)=imapp(itype)
        idp3(4,3,itype)=-imapp(itype)
      enddo
c
      call parcheckpar()
      return
      end


      subroutine parcheckpar()
      implicit none
      integer iallzero,i,jproc,itype,ihpro,i1,i2,i3
      parameter (iallzero=0)
      integer ivbhpro(4,3,5)
      common/civbhpro/ivbhpro
      integer idp1(4,3,5),idp2(4,3,5),idp3(4,3,5)
      common/cidpart/idp1,idp2,idp3
c
      call parcheckinit()
      do i=1,4
        do jproc=1,3
          do itype=1,5
            ihpro=ivbhpro(i,jproc,itype)
            i1=idp1(i,jproc,itype)
            i2=idp2(i,jproc,itype)
            i3=idp3(i,jproc,itype)
            call parcheckfin(ihpro,i1,i2,i3,iallzero)
          enddo
        enddo
      enddo
      return
      end


      subroutine parcheckfin(ihpro,i1,i2,i3,iallzero)
      implicit none
      integer ihpro,i1,i2,i3,iallzero,isum,chin,chout,chall
      real*8 tiny
      parameter (tiny=1.d-8)
      logical ferror
      real*8 chrg(-5:21),chprdct
      common/ccharges/chrg,chprdct
c
      ferror=.false.
      isum=abs(i1)+abs(i2)+abs(i3)
      chin=chrg(i1)+chrg(i2)
      chout=chrg(i3)
      chall=chin-chout-chprdct
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
        if(abs(chall).gt.tiny)ferror=.true.
c 401 is qqbar
        if( ihpro.eq.401 .and.
     #      (i1.le.0 .or. i2.ge.0 .or.
     #       i3.ne.21 .or. (i1+i2).ne.0) )ferror=.true.
c 402 is qg
        if( ihpro.eq.402 .and.
     #      (i1.le.0 .or. i2.ne.21 .or. 
     #       i3.le.0 .or. i1.ne.i3) )ferror=.true.
c 403 is qqbar
        if( ihpro.eq.403 .and.
     #      (i1.ge.0 .or. i2.le.0 .or. 
     #       i3.ne.21 .or. (i1+i2).ne.0) )ferror=.true.
c 404 is qbarg
        if( ihpro.eq.404 .and.
     #      (i1.ge.0 .or. i2.ne.21 .or. 
     #       i3.ge.0 .or. i1.ne.i3) )ferror=.true.
c 405 is gq
        if( ihpro.eq.405 .and.
     #      (i1.ne.21 .or. i2.le.0 .or. 
     #       i3.le.0 .or. i2.ne.i3) )ferror=.true.
c 406 is gqbar
        if( ihpro.eq.406 .and.
     #      (i1.ne.21 .or. i2.ge.0 .or. 
     #       i3.ge.0 .or. i2.ne.i3) )ferror=.true.
c 407 is gg
        if( ihpro.eq.407 .and.
     #      (i1.ne.21 .or. i2.ne.21 .or. 
     #       i3.ne.21) )ferror=.true.
      endif
      if(ferror)then
        write(*,*)'Error in parcheckfin'
        write(*,*)'ihpro,i1,i2,i3:',ihpro,i1,i2,i3
        stop
      endif
      return
      end


      subroutine parcheckinit()
      implicit none
      integer i
      real*8 chup,chdn
      parameter (chup=2.d0/3.d0)
      parameter (chdn=-1.d0/3.d0)
      real*8 chrg(-5:21),chprdct
      common/ccharges/chrg,chprdct
c
      do i=-5,21
        chrg(i)=1000.d0
      enddo
      chrg(1)=chdn
      chrg(2)=chup
      chrg(3)=chdn
      chrg(4)=chup
      chrg(5)=chdn
      chrg(21)=0.d0
      do i=1,5
        chrg(-i)=-chrg(i)
      enddo
      chprdct=0.d0
      return
      end
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
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
      integer ip4,ip5,ip6,ip7,ip8,ip9
      common/ci2part/ip4,ip5,ip6,ip7,ip8,ip9
      integer iccode
      common/ciccode/iccode
      integer idec
      common/cidec/idec
      integer np
      common/cnp/np
      real*8 xevsign
      common/cxevsign/xevsign
      real*8 xmom_lb(11,4)
      common/cxmomlb/xmom_lb
      real*8 ux1,ux2,uq2
      common/uxx/ux1,ux2,uq2
c
      read(iunit,901,end=997,err=998)i1hpro,iccode,np
      if(idec.eq.0)then
        read(iunit,902,end=997,err=998)ip1,ip2,ip3,
     #                                 ip4,ip5,ip6,ip7,ip8,ip9
        read(iunit,903,end=997,err=998)xevsign
        read(iunit,904,end=997,err=998)((xmom_lb(i,j),j=1,4),i=1,3),
     #                                 ((xmom_lb(i,j),j=1,4),i=6,11)
      elseif(idec.eq.1)then
        read(iunit,902,end=997,err=998)ip1,ip2,ip3,ip4,ip5
        read(iunit,903,end=997,err=998)xevsign
        read(iunit,904,end=997,err=998)((xmom_lb(i,j),j=1,4),i=1,5)
      endif
      read(iunit,905,end=997,err=998) ux1,ux2,uq2
      goto 999
 901  format(1x,i3,2(1x,i2))
 902  format(9(1x,i3))
 903  format(1x,d14.8)
 904  format(36(1x,d14.8))
 905  format(3(1x,d14.8))
 997  write(*,*)'unexpected end of file, iunit=',iunit
      stop
 998  write(*,*)'format error'
      write(77,*)'event #:',ii
      write(77,901)i1hpro,iccode,np
      write(77,902)ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9
      write(77,903)xevsign
      write(77,904)((xmom_lb(i,j),j=1,4),i=1,11)
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
c    401        q qbar -> g Q Qb
c    402        q g    -> q Q Qb
c    403        qbar q -> g Q Qb
c    404        qbar g -> qbar Q Qb
c    405        g q    -> q Q Qb
c    406        g qbar -> qbar Q Qb
c    407        g g    -> g Q Qb
c ipX is the parton code relevant to parton # X. PDG conventions are
c used: 1=d, 2=u, 3=s, 4=c, 5=b, 21=g. Note that, at variance with
c what happens for xmom_lb, in the case if ipX only the identities of
c final-state particles are kept. Thus, (ip4,ip5)=(Q,Qbar) identities
c when tops don't decay. But (ip4,ip5,ip6)=identities of top decay
c products, and (ip7,ip8,ip9)=identities of tbar decay products when
c tops decay.
c iccode is the (internal) code which identifies the colour connection
      implicit none
      integer iunit,i,j
      real*8 xpmone,xevwgt,xfact,brfact,dummy
      integer i1hpro
      common/ci1hpro/i1hpro
      integer ip1,ip2,ip3
      common/ci1part/ip1,ip2,ip3
      integer ip4,ip5,ip6,ip7,ip8,ip9
      common/ci2part/ip4,ip5,ip6,ip7,ip8,ip9
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
c xmom_lb(i,j) is the j component of the four vector of the particle # i,
c given in the laboratory frame. j=4 is the energy for MC@NLO versions
c up to 2.31, the mass for version 3.1 onwards. i=1,2 are the incoming
c partons, 3 is the outgoing parton, 4 is Q, 5 is Qbar. When the tops
c decay, 6=l+, 7=nu, 8=b are the decay products of the top, 9=l-, 10=nubar,
c 11=bbar are the decay products of the tbar. Momentum conservation is 
c (1+2)-(3+4+5)=0 or (1+2)-(3+6+7+8+9+10+11)=0
      real*8 xmom_pass(9,4)
      integer IPS(9)
      real*8 xmom_lb(11,4)
      common/cxmomlb/xmom_lb
      integer iwgtnorm
      common/ciwgtnorm/iwgtnorm
      real*8 wgtaev,wgtbev
      common/cwgtev/wgtaev,wgtbev
      real*8 wgtmax
      common/cwgtmax/wgtmax
c Reweight factors that include branching ratios, to be inserted in 
c the case of decayed tops
      real*8 brrtop1,brrtop2
      common/brratios/brrtop1,brrtop2
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
          np=9
          brfact=brrtop1*brrtop2
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
          write(iunit,902)ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9
          write(iunit,903)xevwgt
          write(iunit,904)((xfact*xmom_lb(i,j),j=1,4),i=1,3),
     #                    ((xfact*xmom_lb(i,j),j=1,4),i=6,11)
        elseif(idec.eq.1)then
          write(iunit,902)ip1,ip2,ip3,ip4,ip5
          write(iunit,903)xevwgt
          write(iunit,904)((xmom_lb(i,j),j=1,4),i=1,5)
        else
          write(6,*) 'Error in store_events: idec=',idec
          stop
        endif
        write(iunit,905) ux1,ux2,uq2
      elseif(ievffmt.eq.1.and.xpmone.eq.-1)then
c     The variable dummy can be replaced by xscale if need be
         dummy=0.d0
         if(idec.eq.0) then
            DO I=1,3
               DO J=1,4
                  xmom_pass(I,J)=xfact*xmom_lb(I,J)
               ENDDO
            ENDDO
            do I=6,11
               do J=1,4
                  xmom_pass(I-2,J)=xfact*xmom_lb(I,J)
               ENDDO
            ENDDO            
            IPS(1)=ip1
            IPS(2)=ip2
            IPS(3)=ip3
            IPS(4)=ip4
            IPS(5)=ip5
            IPS(6)=ip6
            IPS(7)=ip7
            IPS(8)=ip8
            IPS(9)=ip9
         else
            DO I=1,5
               DO J=1,4
                  xmom_pass(I,J)=xmom_lb(I,J)
               ENDDO
            ENDDO
            DO I=6,9
               IPS(I)=0
               DO J=1,4
                  xmom_pass(I,J)=0.
               ENDDO
            ENDDO
            IPS(1)=ip1
            IPS(2)=ip2
            IPS(3)=ip3
            IPS(4)=ip4
            IPS(5)=ip5
         ENDIF
         call write_lhef_event(iunit,idec,
     #        i1hpro,iccode,np,IPS,
     #        xevwgt,dummy,
     #        xmom_pass)
      else
         write(*,*)'Error in store_events: unknown file format',ievffmt
         stop
      endif
 901  format(1x,i3,2(1x,i2))
 902  format(9(1x,i3))
 903  format(1x,d14.8)
 904  format(36(1x,d14.8))
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
      subroutine xmcsubtpp(xz1,xz2,xxmq2,xs,xx,xy,xcth1,xcth2,
     #  xmcxsec,xmce0sq,xmcz,xqrksc,xqbrsc,flxsec,flagmc,
     #  gfactsf,gfactcl)
c xmcxsec is the main output of this routine, and is the analogue of the 
c G(x)-part in eqs.(A.83) and (A.84) of FW. It is defined as follows:
c  g^4 xmcxsec = P(z)/((1-x)ktil) (4 tk uk/s**2) 
c                4*M_born Jac(ktil,z,phi;x,y,th2)
c In the arrays xmcxsec, xmce0sq, xmcz, and flxsec the two arguments 
c (ileg,ie0sq) are such that
c   ileg=1,2,3,4 --> emitting leg
c   ie0sq=1,2,3 --> scale choice, corresponding to s/2, -t/2, -u/2.
c xmce0sq the value of the Herwig scale E_0.
c xmcz the value of the Herwig variable z.
c xqrksc(iinv,ileg,ie0sq) have the formal meaning of 2p1.p2, -2p1.k1,
c and -2p2.k1 for iinv=1,2,3, which are the invariants needed to compute 
c the transverse momentum of the quark; xqbrsc is the analogous quantity 
c for the antiquark; pi and ki need not coincide with the actual one, they
c are only used for setting the scales.
c flxsec=.false. when xmcxsec=0, and flagmc=.false. if all flxsec=.false. 
c gfactsf is the value of G_soft(); gfactcl used to be the value of G_coll(),
c is a dummy in this version
      implicit none
      character*2 str
      parameter (str='p1')
      real*8 tiny,vcf,vca
      parameter (tiny=1.d-5)
      parameter (vcf=4.d0/3.d0)
      parameter (vca=3.d0)
      real*8 xz1,xz2,xxmq2,xs,xx,xy,xcth1,xcth2,gfactsf,gfactcl
      real*8 xmcxsec(1:4,1:3),xmce0sq(1:4,1:3),xmcz(1:4,1:3)
      real*8 xqrksc(1:3,1:4,1:3),xqbrsc(1:3,1:4,1:3)
      logical flxsec(1:4,1:3),flagmc,flagxs(1:4)
      real*8 betfac,delta
      common/betfac/betfac,delta
      real*8 alsf(3),besf(3)
      common/cgfunsfp/alsf,besf
      real*8 alcl(3),becl(3)
      common/cgfunclp/alcl,becl
      real*8 alazi(3),beazi(3)
      common/cgfunazi/alazi,beazi
      real*8 z1,z2,xmq2,x,y,cth1,cth2,
     #  si,tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi,
     #  so,tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho,
     #  s,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #  e0sq,de0sqdx,de0sqdc,xz,zherpp,xkp,xqherpp,beta,
     #  gfunsoft,gfactazi,gfunazi,tmp,tmp1,ap,ap_kern,ap_kern_qc,
     #  qin_kern,xfact,qk,xborn,hvqborncol,xazicorr,
     #  hvqcorr,xjac,xjac_zqttoxy,sx,dcosdcth1,
     #  xalsf,xbesf,xalcl,xbecl,xalazi,xbeazi
      real*8 xmckp(1:4,1:3)
      real*8 x2to2(1:4,1:3),dx2to2dx(1:4,1:3),dx2to2dc(1:4,1:3)
      integer jproc0,ifuntype,ileg,ie0sq,index,iborn
      integer ialwsplit(1:3,1:4,1:3)
      integer icllborn(1:3,1:4,1:3),icllkern(1:3,1:4,1:3)
      integer icolconn(1:3,1:4,1:3)
      common/cialwsplit/ialwsplit
      common/cicllsplit/icllborn,icllkern
      common/cicolconn/icolconn
      common/cjproc/jproc0
      common/cifuntype/ifuntype
c
      z1=xz1
      z2=xz2
      xmq2=xxmq2
      s=xs
      x=xx
      y=xy
      cth1=xcth1
      cth2=xcth2
      sx=s*x
      flagxs(1)=.true.
      flagxs(2)=.true.
      flagxs(3)=.true.
      flagxs(4)=.true.
c Compute the invariants
      call invar(xmq2,s,x,y,cth1,cth2,str,
     #           tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi)
      si=s
      if(ifuntype.eq.1)then
        so=si
        tko=tki
        uko=uki
        q1qo=q1qi
        q2qo=q2qi
        q1co=q1ci
        q2co=q2ci
        w1o=w1i
        w2o=w2i
        w1ho=w1hi
        w2ho=w2hi
      elseif(ifuntype.eq.2.and.jproc0.ne.3)then
        if((sx*x).gt.(4*xmq2))then
          call invar(xmq2,sx,x,y,cth1,cth2,str,
     #               tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho)
          so=sx
        else
          flagxs(3)=.false.
          flagxs(4)=.false.
        endif
      endif
c Generate the 2-->2 invariants; those for ileg=2 are identical to the
c corresponding ones for ileg=1, so we skip the generation in such a case
      do ileg=1,4
        if( (ileg.eq.1 .or. (ileg.ge.3.and.jproc0.ne.3)) .and.
     #      flagxs(ileg) )then
          call xinvtoinv(ileg,
     #      si,tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi,
     #      so,tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho,
     #      s,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          call get2to2wr(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #      tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #      x2to2(ileg,1),x2to2(ileg,2),x2to2(ileg,3),
     #      dx2to2dx(ileg,1),dx2to2dx(ileg,2),dx2to2dx(ileg,3),
     #      dx2to2dc(ileg,1),dx2to2dc(ileg,2),dx2to2dc(ileg,3))
        endif
      enddo
      do ie0sq=1,3
        x2to2(2,ie0sq)=x2to2(1,ie0sq)
        dx2to2dx(2,ie0sq)=dx2to2dx(1,ie0sq)
        dx2to2dc(2,ie0sq)=dx2to2dc(1,ie0sq)
      enddo
c Define the invariants for the settings of the scales
      do ileg=1,2
        do ie0sq=1,3
          xqrksc(1,ileg,ie0sq)=si
          xqrksc(2,ileg,ie0sq)=q1qi
          xqrksc(3,ileg,ie0sq)=q2ci
          xqbrsc(1,ileg,ie0sq)=si
          xqbrsc(2,ileg,ie0sq)=q1ci
          xqbrsc(3,ileg,ie0sq)=q2qi
        enddo
      enddo
      do ileg=3,4
        do ie0sq=1,3
          xqrksc(1,ileg,ie0sq)=so
          xqrksc(2,ileg,ie0sq)=q1qo
          xqrksc(3,ileg,ie0sq)=q2co
          xqbrsc(1,ileg,ie0sq)=so
          xqbrsc(2,ileg,ie0sq)=q1co
          xqbrsc(3,ileg,ie0sq)=q2qo
        enddo
      enddo
c Now get z and ktilde for the allowed splittings
      flagmc=.false.
      do ileg=1,4
        call xinvtoinv(ileg,
     #    si,tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi,
     #    so,tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho,
     #    s,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
        do ie0sq=1,3
          if(ialwsplit(jproc0,ileg,ie0sq).eq.1.and.flagxs(ileg))then
            e0sq=abs(x2to2(ileg,ie0sq))/2.d0
            de0sqdx=dx2to2dx(ileg,ie0sq)/2.d0*
     #              sign(1.d0,x2to2(ileg,ie0sq))
            de0sqdc=dx2to2dc(ileg,ie0sq)/2.d0*
     #              sign(1.d0,x2to2(ileg,ie0sq))
            xz=zherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #           cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            xkp=xqherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #           cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
            if(ileg.le.2)then
              if( xz.ge.0.d0 .and. xkp.ge.0.d0. and.
     #            xkp.le.(2*e0sq) )then
                flxsec(ileg,ie0sq)=.true.
                if(.not.flagmc)flagmc=.true.
              else
                flxsec(ileg,ie0sq)=.false.
              endif
            else
              beta=sqrt(1-4*xmq2/s)
              if( xz.ge.0.d0 .and. xkp.ge.0.d0. and.
     #            ( (ie0sq.eq.1.and.xkp.le.(e0sq*(1+beta))) .or.
     #              (ie0sq.gt.1.and.xkp.le.(2*e0sq+xmq2)) ) )then
                flxsec(ileg,ie0sq)=.true.
                if(.not.flagmc)flagmc=.true.
              else
                flxsec(ileg,ie0sq)=.false.
              endif
            endif
            if(flxsec(ileg,ie0sq))then
              xmcz(ileg,ie0sq)=xz
              xmckp(ileg,ie0sq)=xkp
              xmce0sq(ileg,ie0sq)=e0sq
            endif
          else
            flxsec(ileg,ie0sq)=.false.
          endif
        enddo
      enddo
c Even if flagmc=.false., evaluate the G functions, since they are used
c to compute the ME part of the MC subtraction term
      xalsf=alsf(jproc0)
      xbesf=besf(jproc0)
      xalcl=alcl(jproc0)
      xbecl=becl(jproc0)
      xalazi=alazi(jproc0)
      xbeazi=beazi(jproc0)
      gfactsf=gfunsoft(x,s,xmq2,xalsf,xbesf)
      gfactcl=-1.d8
      gfactazi=gfunazi(y,xalazi,xbeazi,delta)
c Compute the cross sections if in the live zones. Also add the azimuthal
c correlation term (whose kernel is tmp1) when gfactazi#0; since this is
c of ME origin, the sum over ie0sq should not be performed. We (arbitrarily)
c associate it to the first non-zero ie0sq contribution
      if(flagmc)then
        do ileg=1,4
          call xinvtoinv(ileg,
     #      si,tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi,
     #      so,tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho,
     #      s,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          tmp1=0.d0
          do ie0sq=1,3
            if(flxsec(ileg,ie0sq))then
              e0sq=xmce0sq(ileg,ie0sq)
              if(1-x.lt.tiny)then
                tmp=0.d0
              elseif(1-y.lt.tiny.and.ileg.eq.1)then
                index=icllkern(jproc0,ileg,ie0sq)
                tmp=(1+y)*ap_kern(x,index)
                if(gfactazi.ne.0.d0.and.tmp1.eq.0.d0)
     #            tmp1=(1+y)*qin_kern(x,index)
              elseif(1+y.lt.tiny.and.ileg.eq.2)then
                index=icllkern(jproc0,ileg,ie0sq)
                tmp=(1-y)*ap_kern(x,index)
                if(gfactazi.ne.0.d0.and.tmp1.eq.0.d0)
     #            tmp1=(1-y)*qin_kern(x,index)
              else
                xfact=(1-x)*(1-y**2)
                xz=xmcz(ileg,ie0sq)
                xkp=xmckp(ileg,ie0sq)
                index=icllkern(jproc0,ileg,ie0sq)
                if(ileg.le.2)then
                  ap=ap_kern(xz,index)/(1-xz)
                else
                  ap=ap_kern_qc(xz,xkp,xmq2,index)/(1-xz)
                endif
                de0sqdx=dx2to2dx(ileg,ie0sq)/2.d0*
     #                  sign(1.d0,x2to2(ileg,ie0sq))
                de0sqdc=dx2to2dc(ileg,ie0sq)/2.d0*
     #                  sign(1.d0,x2to2(ileg,ie0sq))
                xjac=xjac_zqttoxy(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,
     #               y,cth1,cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
                tmp=xfact*xjac*ap/xkp
                if(ileg.le.2.and.gfactazi.ne.0.d0.and.tmp1.eq.0.d0)then
                  qk=qin_kern(xz,index)/(1-xz)
                  tmp1=xfact*xjac*qk/xkp
                endif
              endif
              iborn=icllborn(jproc0,ileg,ie0sq)
              xborn=hvqborncol(x2to2(ileg,1),x2to2(ileg,2),xmq2,
     #                         iborn,jproc0,ileg,ie0sq)
              if(ileg.le.2)xazicorr=hvqcorr(x2to2(ileg,1),
     #          x2to2(ileg,2),cth2,xmq2,iborn)
              xmcxsec(ileg,ie0sq)=4*tmp*xborn+
     #                            4*tmp1*xazicorr*gfactazi
              xmcxsec(ileg,ie0sq)=xmcxsec(ileg,ie0sq)*gfactsf*
     #                            dcosdcth1(ileg,x,cth1)
            else
              xmcxsec(ileg,ie0sq)=0.d0
            endif
          enddo
        enddo
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
      if(ileg.eq.3.or.ileg.eq.4)then
        tmp=4*x/(1+x-(1-x)*cth1)**2
      elseif(ileg.ne.1.and.ileg.ne.2)then
        write(*,*)'Error in dcosdcth1: unknown ileg',ileg
        stop
      endif
      dcosdcth1=tmp
      return
      end
    

      subroutine xinvtoinv(ileg,
     #  si,tki,uki,q1qi,q2qi,q1ci,q2ci,w1i,w2i,w1hi,w2hi,
     #  so,tko,uko,q1qo,q2qo,q1co,q2co,w1o,w2o,w1ho,w2ho,
     #  s,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      implicit real*8(a-h,o-z)
      integer ileg
c
      if(ileg.le.2)then
        s=si
        tk=tki
        uk=uki
        q1q=q1qi
        q2q=q2qi
        q1c=q1ci
        q2c=q2ci
        w1=w1i
        w2=w2i
        w1h=w1hi
        w2h=w2hi
      elseif(ileg.le.4)then
        s=so
        tk=tko
        uk=uko
        q1q=q1qo
        q2q=q2qo
        q1c=q1co
        q2c=q2co
        w1=w1o
        w2=w2o
        w1h=w1ho
        w2h=w2ho
      else
        write(*,*)'Unknown leg in xinvtoinv:',ileg
        stop
      endif
      return
      end


      function gfunsoft(xx,xs,xxmq2,alsf,besf)
c Gets smoothly to 0 in the soft limit. The functional form is given
c in eq.(A.86) of FW, with alpha==alsf. tilde{x}_{DZ} is replaced here
c by xgsoft, and x_{DZ} by xminsf. The function is different from 1
c in the region xgsoft<x<1. Call with
c  besf<0  ==> xminsf=4*m2/S_{hadr}
c  besf>0  ==> xminsf=tilde{rho} for standard subtraction
c              xminsf=1-sqrt{zeta} for zeta-subtraction
c When besf>0, besf-->0 ==> xgsoft-->1
c              besf-->1 ==> xgsoft-->xminsf
c If alsf<0, gfunsoft equals 1 everywhere. This option should be used
c for testing purposes only
      implicit none
      real*8 gfunsoft,xx,xs,xxmq2,alsf,besf,x,s,xmq2,xminsf,xgsoft,
     # tt,tmp
      real*8 sh
      common/shadr/sh
      real*8 betfac,delta
      common/betfac/betfac,delta
      real*8 etacut
      common/cetacut/etacut
      integer isubttype
      common/cisubttype/isubttype
c
      x=xx
      s=xs
      xmq2=xxmq2
      if(besf.lt.0.d0)then
        xminsf=4*xmq2/sh
      else
        if(isubttype.eq.0)then
c This is tilde{ro}; don't use tilde{rox}, G_soft gets too fast to zero
          xminsf=1.d0-(1-4*xmq2/s)*betfac**2
        elseif(isubttype.eq.1)then
          xminsf=1-sqrt(etacut)
        else
          write(*,*)'Fatal error #1 in gfunsoft',isubttype
          stop
        endif
      endif
      xgsoft=1.d0-(1-xminsf)*abs(besf)
      if(xgsoft.gt.0.99d0)xgsoft=0.99d0
      tt=(x-xgsoft)/(1.d0-xgsoft)
      if(tt.gt.1.d0)then
        write(6,*)'Fatal error #2 in gfunsoft',x
        stop
      endif
      tmp=1.d0
      if(alsf.gt.0.d0)then
        if(tt.gt.0.d0.and.x.lt.0.99d0)
     #    tmp=(1-tt)**(2*alsf)/(tt**(2*alsf)+(1-tt)**(2*alsf))
        if(x.ge.0.99d0)tmp=0.d0
      endif
      gfunsoft=tmp
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
c This function differs from the original one of the QQ code; 
c alcl, becl and delta are now given in input rather than in common;
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
c subtraction kernel; it is not the same as in the old QQ code. We have
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


c Taken from hdyjetdiff.for; the splitting is b(p)-->a(x*p)+c((1-x)*p)
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
        write(6,*)'Error in ap_kern: wrong index value',index
        stop
      endif
      return
      end


      function qin_kern(x,index)
c This function returns the quantity (1-x)*Q_{a*b}(x), where
c Q_{a*b} are the kernels relevant to the azimuthal correlations
c in the splittings of incoming partons; their explicit form can be
c found in eqs.(B.42)--(B.45) of FKS (NPB467(96)399). The splitting 
c partons {a*b} are defined with the following conventions
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
        write(6,*)'Error in qin_kern: wrong index value',index
        stop
      endif
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


      function hvqborn(xs,xt,xxmq2,ijproc)
c Returns the Born cross section for the partonic process whose initial 
c state is identified by ijproc. The normalization is such that
c dsigma_born(s,t) = g^4 hvqborn(s,t) dphi2(s)
c (see the function ggborn)
      implicit none
      real*8 hvqborn,xs,xt,xxmq2,s,t,xmq2,tmp,ggborn,qqborn
      integer ijproc,jproc
c
      s=xs
      t=xt
      xmq2=xxmq2
      jproc=ijproc
      if(jproc.eq.1)then
        tmp=ggborn(s,t,xmq2)
      elseif(jproc.eq.2)then
        tmp=qqborn(s,t,xmq2)
      elseif(jproc.eq.3)then
        tmp=0.d0
      else
        write(*,*)'Unknown process in hvqborn',jproc
        stop
      endif
      hvqborn=tmp
      return
      end


      function hvqborncol(xs,xt,xxmq2,iborn,jproc,ileg,ie0sq)
c Returns the 2-->2 reduced matrix element squared (d\bar{sigma}) which
c factorizes in the MC subtraction term; the definition is
c   hvqborncol = factor * hvqborn,   where
c   factor = 1/|icolconn|*[1,u**2/(u**2+t**2),t**2/(u**2+t**2)], 
c the term is square brackets depending on the colour structure; icolconn 
c is defined in parsetpar, iborn is the 2-->2 initial state, whereas jproc 
c identifies the 2-->3 one.
      implicit none
      real*8 hvqborncol,xs,xt,xxmq2,s,t,xmq2,u,xfact,xborn,hvqborn
      integer iborn,jproc,ileg,ie0sq
      integer icolconn(1:3,1:4,1:3)
      common/cicolconn/icolconn
c
      s=xs
      t=xt
      xmq2=xxmq2
      u=-s-t
      if(icolconn(jproc,ileg,ie0sq).eq.0)then
        xfact=0.d0
      else
        xborn=hvqborn(s,t,xmq2,iborn)
        xfact=1.d0/dfloat(abs(icolconn(jproc,ileg,ie0sq)))
        if(icolconn(jproc,ileg,ie0sq).gt.0)then
          if(ie0sq.eq.1)then
            write(*,*)
     #        'hvqborncol: no such configuration in this process'
            stop
          elseif(ie0sq.eq.2)then
            xfact=xfact*u**2/(t**2+u**2)
          elseif(ie0sq.eq.3)then
            xfact=xfact*t**2/(t**2+u**2)
          else
            write(*,*)'Fatal error in hvqborncol',ie0sq
            stop
          endif
        endif
      endif
      hvqborncol=xfact*xborn
      return
      end


      function hvqcorr(xs,xt,xcth2,xxmq2,iborn)
c Returns the azimuthal correlation term; it is only used in the case
c in which G_coll does not vanish in the collinear limit for any value
c of x
      implicit real*8(a-h,o-z)
      parameter (xnc=3.d0)
c
      s=xs
      t=xt
      cth2=xcth2
      xmq2=xxmq2
      u=-s-t
      if(iborn.eq.1)then
        tmp=-1.d0/(xnc**2-1)*(u/t+t/u-s**2/(xnc**2*t*u))*
     #       (xmq2/s-xmq2**2/(t*u))
        tmp=(2*cth2**2-1)*tmp
      else
        tmp=0.d0
      endif
      hvqcorr=tmp/(2*s)
      return
      end


      function xnlfscheme(xm2,xmur2,xmuf2h1,xmuf2h2,zg2,jproc)
c Returns the factor which multiplies sigma_Born with nl (nlfp1sch=0)
c or nl+1 (nlfp1sch=1) schemes; process code is 1,2,3 for gg,qq,qg
      implicit none
      real*8 xnlfscheme,xm2,xmur2,xmuf2h1,xmuf2h2,zg2,alfas,tmp
      real * 8 tf
      parameter (tf=0.5d0)
      real * 8 pi
      parameter (pi=3.14159265358979312D0)
      integer jproc,nlfp1sch
      common/cnlfp1sch/nlfp1sch
c
      if(nlfp1sch.eq.0)then
        tmp=0.d0
      elseif(nlfp1sch.eq.1)then
        alfas=zg2/(4.d0*pi)
        if(jproc.eq.1)then
          tmp=-alfas*tf/(3.d0*pi)*( log(xmur2/xmuf2h1)+
     #                              log(xmur2/xmuf2h2) )
        elseif(jproc.eq.2)then
          tmp=-alfas*2*tf/(3.d0*pi)*log(xmur2/xm2)
        elseif(jproc.eq.3)then
          tmp=0.d0
        else
          write(*,*)'Unknown process in xnlfscheme:',jproc
        endif
      else
        write(*,*)'Unknown scheme in xnlfscheme:',nlfp1sch
      endif
      xnlfscheme=tmp
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
      function zherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #                cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      implicit none
      real*8 zherpp,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,cth2,tk,uk,
     # q1q,q2q,q1c,q2c,w1,w2,w1h,w2h
      integer ileg
      real*8 tiny,xm2red,beta,sth1,cpsip,spsip,cthg,beta1,beta2,
     # zeta1,zeta2
      parameter (tiny=1.d-5)
c
      sth1=sqrt(1-cth1**2)
c incoming parton (left)
      if(ileg.eq.1)then
        zherpp=(1-y+x*(1+y))/2.d0
c incoming parton (right)
      elseif(ileg.eq.2)then
        zherpp=(1+y+x*(1-y))/2.d0
c outgoing heavy quark
      elseif(ileg.eq.3)then
        xm2red=xmq2/s
        beta=sqrt(1-4*xm2red)
        if(1-x.lt.tiny)then
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          zherpp=1-(1+cthg)*(1-x)/(1+beta)
        else
          beta2=sqrt(1-4*xmq2*s/(s-w1)**2)
          zeta1=( (s+w1)*w2+(s-w1)*((w1+w2)*beta2-w1) )/
     #          ( (s-w1)*beta2*(s+w1+(s-w1)*beta2) )
          zherpp=1-zeta1
        endif
c outgoing heavy antiquark
      elseif(ileg.eq.4)then
        xm2red=xmq2/s
        beta=sqrt(1-4*xm2red)
        if(1-x.lt.tiny)then
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          zherpp=1-(1-cthg)*(1-x)/(1+beta)
        else
          beta1=sqrt(1-4*xmq2*s/(s-w2)**2)
          zeta2=( (s+w2)*w1+(s-w2)*((w2+w1)*beta1-w2) )/
     #          ( (s-w2)*beta1*(s+w2+(s-w2)*beta1) )
          zherpp=1-zeta2
        endif
      else
        write(6,*)'zherpp: unknown parton number'
        stop
      endif
      if(zherpp.lt.0.d0.or.zherpp.gt.1.d0)then
        write(6,*)'zherpp: fatal error',zherpp
        stop
      endif
      return
      end


      function xqherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #                cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      implicit none
      real*8 xqherpp,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,cth2,tk,uk,
     # q1q,q2q,q1c,q2c,w1,w2,w1h,w2h
      integer ileg
      real*8 tiny,vtiny,z,zherpp,sth1,beta,cpsip,spsip,cthg
      parameter (tiny=1.d-5)
      parameter (vtiny=1.d-8)
c
      sth1=sqrt(1-cth1**2)
c incoming parton (left)
      if(ileg.eq.1)then
        if(y.gt.-1+tiny)then
          xqherpp=s*(1-y)/(1+y)
        else
          xqherpp=-1.d0
        endif
c incoming parton (right)
      elseif(ileg.eq.2)then
        if(y.lt.1-tiny)then
          xqherpp=s*(1+y)/(1-y)
        else
          xqherpp=-1.d0
        endif
c outgoing heavy quark
      elseif(ileg.eq.3)then
        if(1-x.lt.tiny)then
          beta=sqrt(1-4*xmq2/s)
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          xqherpp=(1+beta)*(1-beta*cthg)*s/(2*(1+cthg))
        else
          z=zherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #             cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          if(z.lt.1-vtiny.and.z.gt.vtiny)then
            xqherpp=w1/(z*(1-z))
          else
            xqherpp=-1.d0
          endif
        endif
c outgoing heavy antiquark
      elseif(ileg.eq.4)then
        if(1-x.lt.tiny)then
          beta=sqrt(1-4*xmq2/s)
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          xqherpp=(1+beta)*(1+beta*cthg)*s/(2*(1-cthg))
        else
          z=zherpp(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #             cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
          if(z.lt.1-vtiny.and.z.gt.vtiny)then
            xqherpp=w2/(z*(1-z))
          else
            xqherpp=-1.d0
          endif
        endif
      else
        write(6,*)'xqherpp: unknown parton number'
        stop
      endif
      return
      end


      function xjac_zqttoxy(ileg,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,
     #                      cth2,tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h)
      implicit none
      real*8 xjac_zqttoxy,e0sq,de0sqdx,de0sqdc,xmq2,s,x,y,cth1,cth2,
     # tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h
      integer ileg
      real*8 tiny,tmp,sth1,xm2red,beta,cpsip,spsip,cthg,zmo,dktil0sfdc,
     # betax,beta2,dw1dx,dw2dx,dw1dc,dw2dc,zeta1,dzeta1dw2,beta1,zeta2,
     # dzeta2dw1
      parameter (tiny=1.d-5)
c
      tmp=0.d0
      sth1=sqrt(1-cth1**2)
c incoming parton (left)
      if(ileg.eq.1)then
        if(y.gt.-1+tiny)then
          tmp=-s/(1+y)
        else
          tmp=0.d0
        endif
c incoming parton (right)
      elseif(ileg.eq.2)then
        if(y.lt.1-tiny)then
          tmp=s/(1-y)
        else
          tmp=0.d0
        endif
c outgoing heavy quark
      elseif(ileg.eq.3)then
        xm2red=xmq2/s
        beta=sqrt(1-4*xm2red)
        cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
        spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
        cthg=cpsip*cth1+spsip*cth2*sth1
        if(1-x.lt.tiny)then
          zmo=-(1+cthg)/(1+beta)
          dktil0sfdc=-beta*(1+beta)*s/(2*(1+cthg))-
     #               (1+beta)*(1-beta*cthg)*s/(2*(1+cthg)**2)
          tmp=-zmo*dktil0sfdc
        else
          betax=sqrt(1-4*xmq2/(s*x))
          beta2=sqrt(1-4*xmq2*s/(s-w1)**2)
          dw1dx=-s/2.d0*(1-betax*cthg)
     #          -cthg*xmq2*(1-x)/(betax*x**2)
          dw2dx=-s/2.d0*(1+betax*cthg)
     #          +cthg*xmq2*(1-x)/(betax*x**2)
          dw1dc=-s/2.d0*betax*(1-x)
          dw2dc=s/2.d0*betax*(1-x)
          zeta1=( (s+w1)*w2+(s-w1)*((w1+w2)*beta2-w1) )/
     #          ( (s-w1)*beta2*(s+w1+(s-w1)*beta2) )
          dzeta1dw2=1/((s-w1)*beta2)
          tmp=(dw1dx*dw2dc-dw1dc*dw2dx)*dzeta1dw2/(zeta1*(1-zeta1))
          tmp=tmp*4*x/(1+y+x*(1-y))**2
        endif
c outgoing heavy antiquark
      elseif(ileg.eq.4)then
        xm2red=xmq2/s
        beta=sqrt(1-4*xm2red)
        cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
        spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
        cthg=cpsip*cth1+spsip*cth2*sth1
        if(1-x.lt.tiny)then
          zmo=-(1-cthg)/(1+beta)
          dktil0sfdc=beta*(1+beta)*s/(2*(1-cthg))+
     #               (1+beta)*(1+beta*cthg)*s/(2*(1-cthg)**2)
          tmp=-zmo*dktil0sfdc
        else
          betax=sqrt(1-4*xmq2/(s*x))
          beta1=sqrt(1-4*xmq2*s/(s-w2)**2)
          dw1dx=-s/2.d0*(1-betax*cthg)
     #          -cthg*xmq2*(1-x)/(betax*x**2)
          dw2dx=-s/2.d0*(1+betax*cthg)
     #          +cthg*xmq2*(1-x)/(betax*x**2)
          dw1dc=-s/2.d0*betax*(1-x)
          dw2dc=s/2.d0*betax*(1-x)
          zeta2=( (s+w2)*w1+(s-w2)*((w2+w1)*beta1-w2) )/
     #          ( (s-w2)*beta1*(s+w2+(s-w2)*beta1) )
          dzeta2dw1=1/((s-w2)*beta1)
          tmp=-(dw1dx*dw2dc-dw1dc*dw2dx)*dzeta2dw1/(zeta2*(1-zeta2))
          tmp=tmp*4*x/(1+y+x*(1-y))**2
        endif
      else
        write(6,*)'xjac_zqttoxy: unknown parton number'
        stop
      endif
      xjac_zqttoxy=abs(tmp)
      return 
      end 


      subroutine get2to2wr(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #                     tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #                     sbar,tbar,ubar,dsbardx,dtbardx,dubardx,
     #                     dsbardc,dtbardc,dubardc)
c Wrapper for get2to2a1 and get2to2a2
      implicit real*8(a-h,o-z)
      common/cia1ora2/ia1ora2
c
      if(ia1ora2.eq.1)then
        call get2to2a1(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #                 tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #                 sbar,tbar,ubar,dsbardx,dtbardx,dubardx,
     #                 dsbardc,dtbardc,dubardc)
      elseif(ia1ora2.eq.2)then
        call get2to2a2(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #                 tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #                 sbar,tbar,ubar,dsbardx,dtbardx,dubardx,
     #                 dsbardc,dtbardc,dubardc)
      else
        write(6,*)'Fatal error in get2to2wr: unknown option'
        write(6,*)ia1ora2
        stop
      endif
      return
      end


      subroutine get2to2a1(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #                     tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #                     sbar,tbar,ubar,dsbardx,dtbardx,dubardx,
     #                     dsbardc,dtbardc,dubardc)
c Given the 2->3 kinematics, computes the 2->2 invariants sbar, tbar,
c and ubar. In the computation of the jacobian, the derivatives of
c these quantities are used as well. d{s,t,u}bardx is the derivative
c with respect to x, at y constant for legs 1 and 2, and at cthg and phig
c constant for legs 3 and 4. d{s,t,u}bardc is the derivative with respect
c to cthg for legs 3 and 4, whereas d{s,t,u}bardc=(1-y**2)*diff({s,t,b}bar,y)
c for legs 1 and 2
      implicit real*8(a-h,o-z)
      parameter (tiny=1.d-5)
c
      sth1=sqrt(1-cth1**2)
      beta=sqrt(1-4*xmq2/s)
      betax=sqrt(1-4*xmq2/(s*x))
      if(ileg.eq.1.or.ileg.eq.2)then
        sbar=x*s
        dsbardx=s
        dsbardc=0.d0
        if(1-x.lt.tiny)then
          tbar=-s/2.d0*(1-beta*cth1)
          dtbardx=-1/2.d0*( s - (s-2*xmq2)*cth1/beta -
     #          s*beta*sqrt(1-y**2)*cth2*sth1*z1/(z1+z2) )
          dtbardc=0.d0
        elseif(1-y.lt.tiny)then
          tbar=-s*x/2.d0*(1-betax*cth1)
          dtbardx=-1/2.d0*( s - (s*x-2*xmq2)*cth1/(x*betax) )
          dtbardc=0.d0
        elseif(1+y.lt.tiny)then
          tbar=-s*x/2.d0*(1-betax*cth1)
          dtbardx=-1/2.d0*( s - (s*x-2*xmq2)*cth1/(x*betax) )
          dtbardc=0.d0
        else
          cpsi=1-8*x/((1+y+x*(1-y))*(1-y+x*(1+y)))
          spsi=4*(1-x)*sqrt(x*(1-y**2))/
     #         ((1+y+x*(1-y))*(1-y+x*(1+y)))
          xpl=((s+uk)*z1/s+(s+tk)*z2/s)/2.d0
          dtkdx=s*(1-y)/2.d0
          dukdx=s*(1+y)/2.d0
          dq1qdx=xmq2*(1+y+x*(1-y))*cth1/(2*betax*x**2) - 
     #            s/4.d0*(1-y)*(1-betax*cth1)
          dq2qdx=-s/4.d0*(1+y)*(1+betax*(cpsi*cth1+
     #                          spsi*cth2*sth1)) -
     #            xmq2*(1-y+x*(1+y))*(cpsi*cth1+
     #            spsi*cth2*sth1)/(2*betax*x**2) +
     #            2*betax*s*(1-x**2)*(1-y**2)*cth1/
     #            ((1+y+x*(1-y))**2*(1-y+x*(1+y))) -
     #            betax*s*cth2*sth1*
     #            (1+x)*(1-y**2)*((1-y**2)*(1+x**2)-2*x*(3-y**2))/
     #            (2*(1+y+x*(1-y))**2*(1-y+x*(1+y))*sqrt(x*(1-y**2)))
          dxpldx=(z1*(1+y)+z2*(1-y))/4.d0
          dtkdy=s*(1-x)/2.d0
          dukdy=-s*(1-x)/2.d0
          dq1qdy=-s/4.d0*(1-x)*(1-betax*cth1)
          dq2qdy=s/4.d0*(1-x)*(1+betax*(cpsi*cth1+
     #                         spsi*cth2*sth1)) +
     #           4*betax*s*(1-x)**2*x*y*cth1/
     #           ((1+y+x*(1-y))**2*(1-y+x*(1+y))) -
     #           betax*s*(1-x)*x*y*((1-y**2)*(1+x**2)-
     #           2*x*(3-y**2))*cth2*sth1/
     #           ( (1+y+x*(1-y))**2*(1-y+x*(1+y))*sqrt(x*(1-y**2)) )
          dxpldy=-(1-x)*(z1-z2)/4.d0
          tbar=-(s+tk+uk)/2.d0*( 1-(z2*(q1q-q1c)+z1*(q2q-q2c))/
     #                    (2*s*sqrt(xpl**2-z1*z2*tk*uk/s**2)) )
          dtbardx=-(s+tk+uk)/2.d0*(
     #      ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )*
     #      ( 2*dxpldx*xpl-dukdx*tk*z1*z2/s**2-dtkdx*uk*z1*z2/s**2 )/
     #      ( 4*s*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) -
     #      ( (2*dq2qdx+dukdx)*z1+(2*dq1qdx+dtkdx)*z2 )/
     #      ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) ) -
     #      (dtkdx+dukdx)/2.d0*(1 -
     #          ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )/
     #          ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) )
          dtbardc=-(s+tk+uk)/2.d0*(
     #      ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )*
     #      ( 2*dxpldy*xpl-dukdy*tk*z1*z2/s**2-dtkdy*uk*z1*z2/s**2 )/
     #      ( 4*s*(xpl**2-tk*uk*z1*z2/s**2)**(1.5d0) ) -
     #      ( (2*dq2qdy+dukdy)*z1+(2*dq1qdy+dtkdy)*z2 )/
     #      ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) ) -
     #      (dtkdy+dukdy)/2.d0*(1 -
     #          ( (2*q2q+s+uk)*z1+(2*q1q+s+tk)*z2 )/
     #          ( 2*s*sqrt(xpl**2-tk*uk*z1*z2/s**2) ) )
          dtbardc=(1-y**2)*dtbardc
        endif
      elseif(ileg.eq.3)then
        sbar=s
        dsbardx=0.d0
        dsbardc=0.d0
        if(1-x.lt.tiny)then
          tbar=-s/2.d0*(1-beta*cth1)
          cthg=y*cth1+sqrt(1-y**2)*cth2*sth1
          dtbardx=s/4.d0*(cthg*cth1-y+beta*(cthg-y*cth1-sth1**2))
          dtbardc=0.d0
        else
          beta2=sqrt(1-4*xmq2*s/(s-w1)**2)
          tbar=-s/2.d0*(1-(q2q-q1c)/(s-w1)*beta/beta2)
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          ctho=-(1-x-(1+x)*cth1)/(1+x-(1-x)*cth1)
          dydcpsip=4*x/(1+x-cpsip*(1-x))**2
          dydx=2*(1-cpsip**2)/(1+x-cpsip*(1-x))**2
          dcpsipdx=-sth1*( cthg*sth1 + 
     #                     cth1*(cth1*cth2*spsip-cpsip*sth1) )/
     #             (2.d0*x)
          dcpsipdcthg=cth1+sth1*(cth2*cth1*spsip-cpsip*sth1)*
     #                cthg/(1-cthg**2)
          dcth1dx=-2*(1-ctho**2)/(1+x+(1-x)*ctho)**2
          dtkdx=s*(1-y)/2.d0+
     #          s*(1-x)/2.d0*(dydx+dydcpsip*dcpsipdx)
          dukdx=s*(1+y)/2.d0-
     #          s*(1-x)/2.d0*(dydx+dydcpsip*dcpsipdx)
          dw1dx=-s/2.d0*(1-betax*cthg)-xmq2*(1-x)*cthg/(betax*x**2)
          dq1qdx=-dtkdx/2.d0*(1-betax*cth1)+
     #            s/4.d0*(1+y+x*(1-y))*
     #            (2*xmq2/(s*betax*x**2)*cth1+betax*dcth1dx)
          dtkdc=s/2.d0*(1-x)*dydcpsip*dcpsipdcthg
          dukdc=-s/2.d0*(1-x)*dydcpsip*dcpsipdcthg
          dw1dc=-s/2.d0*(1-x)*betax
          dq1qdc=-dtkdc/2.d0*(1-betax*cth1)
          dtbardx=s/2.d0*( 
     #                  beta*(2*dq1qdx+2*dtkdx+dw1dx)/(beta2*(s-w1)) 
     #                  +4*beta*dw1dx*xmq2*s*(2*q1q+s+2*tk+w1)/
     #                   (beta2**3*(s-w1)**4) 
     #              +beta*dw1dx*(2*q1q+s+2*tk+w1)/(beta2*(s-w1)**2) )
          dtbardc=s/2.d0*( 
     #                  beta*(2*dq1qdc+2*dtkdc+dw1dc)/(beta2*(s-w1)) 
     #                  +4*beta*dw1dc*xmq2*s*(2*q1q+s+2*tk+w1)/
     #                   (beta2**3*(s-w1)**4) 
     #              +beta*dw1dc*(2*q1q+s+2*tk+w1)/(beta2*(s-w1)**2) )
        endif
      elseif(ileg.eq.4)then
        sbar=s
        dsbardx=0.d0
        dsbardc=0.d0
        if(1-x.lt.tiny)then
          tbar=-s/2.d0*(1-beta*cth1)
          cthg=y*cth1+sqrt(1-y**2)*cth2*sth1
          dtbardx=s/4.d0*(-cthg*cth1+y+beta*(cthg-y*cth1-sth1**2))
          dtbardc=0.d0
        else
          beta1=sqrt(1-4*xmq2*s/(s-w2)**2)
          tbar=-s/2.d0*(1-(q1q-q2c)/(s-w2)*beta/beta1)
          cpsip=(1+y-x*(1-y))/(1+y+x*(1-y))
          spsip=sqrt(4*x*(1-y**2))/(1+y+x*(1-y))
          cthg=cpsip*cth1+spsip*cth2*sth1
          ctho=-(1-x-(1+x)*cth1)/(1+x-(1-x)*cth1)
          dydcpsip=4*x/(1+x-cpsip*(1-x))**2
          dydx=2*(1-cpsip**2)/(1+x-cpsip*(1-x))**2
          dcpsipdx=-sth1*( cthg*sth1 + 
     #                     cth1*(cth1*cth2*spsip-cpsip*sth1) )/
     #             (2.d0*x)
          dcpsipdcthg=cth1+sth1*(cth2*cth1*spsip-cpsip*sth1)*
     #                cthg/(1-cthg**2)
          dcth1dx=-2*(1-ctho**2)/(1+x+(1-x)*ctho)**2
          dtkdx=s*(1-y)/2.d0+
     #          s*(1-x)/2.d0*(dydx+dydcpsip*dcpsipdx)
          dukdx=s*(1+y)/2.d0-
     #          s*(1-x)/2.d0*(dydx+dydcpsip*dcpsipdx)
          dw2dx=-s/2.d0*(1+betax*cthg)+xmq2*(1-x)*cthg/(betax*x**2)
          dq1qdx=-dtkdx/2.d0*(1-betax*cth1)+
     #            s/4.d0*(1+y+x*(1-y))*
     #            (2*xmq2/(s*betax*x**2)*cth1+betax*dcth1dx)
          dtkdc=s/2.d0*(1-x)*dydcpsip*dcpsipdcthg
          dukdc=-s/2.d0*(1-x)*dydcpsip*dcpsipdcthg
          dw2dc=s/2.d0*(1-x)*betax
          dq1qdc=-dtkdc/2.d0*(1-betax*cth1)
          dtbardx=s/2.d0*( 
     #                  beta*(2*dq1qdx-dw2dx)/(beta1*(s-w2))
     #                  +4*beta*dw2dx*xmq2*s*(2*q1q+s-w2)/
     #                  (beta1**3*(s-w2)**4) +
     #                  beta*dw2dx*(2*q1q+s-w2)/(beta1*(s-w2)**2) )
          dtbardc=s/2.d0*( 
     #                  beta*(2*dq1qdc-dw2dc)/(beta1*(s-w2))
     #                  +4*beta*dw2dc*xmq2*s*(2*q1q+s-w2)/
     #                  (beta1**3*(s-w2)**4) +
     #                  beta*dw2dc*(2*q1q+s-w2)/(beta1*(s-w2)**2) )
        endif
      else
        write(6,*)'Fatal error in get2to2a1: unknown leg'
        write(6,*)ileg
        stop
      endif
      ubar=-sbar-tbar
      dubardx=-dsbardx-dtbardx
      dubardc=-dsbardc-dtbardc
      return
      end


      subroutine get2to2a2(ileg,z1,z2,xmq2,s,x,y,cth1,cth2,
     #                     tk,uk,q1q,q2q,q1c,q2c,w1,w2,w1h,w2h,
     #                     sbar,tbar,ubar,dsbardx,dtbardx,dubardx,
     #                     dsbardc,dtbardc,dubardc)
c Given the 2->3 kinematics, computes the 2->2 invariants sbar, tbar,
c and ubar. In the computation of the jacobian, the derivatives of
c these quantities are used as well. d{s,t,u}bardx is the derivative
c with respect to x, at y constant for legs 1 and 2, and at cthg and phig
c constant for legs 3 and 4. d{s,t,u}bardc is the derivative with respect
c to cthg for legs 3 and 4, whereas d{s,t,u}bardc=(1-y**2)*diff({s,t,b}bar,y)
c for legs 1 and 2
      implicit real*8(a-h,o-z)
      parameter (tiny=1.d-5)
c
      sth1=sqrt(1-cth1**2)
      beta=sqrt(1-4*xmq2/s)
      betax=sqrt(1-4*xmq2/(s*x))
      if(ileg.eq.1.or.ileg.eq.2)then
        sbar=x*s
        dsbardx=s
        dsbardc=0.d0
        tbar=-s*x/2.d0*(1-betax*cth1)
        dtbardx=-s/2.d0*(1-betax*cth1)+xmq2*cth1/(betax*x)
        dtbardc=0.d0
        dtbardc=(1-y**2)*dtbardc
      elseif(ileg.eq.3.or.ileg.eq.4)then
        ctho=-(1-x-(1+x)*cth1)/(1+x-(1-x)*cth1)
        sbar=s
        dsbardx=0.d0
        dsbardc=0.d0
        tbar=-s/2.d0*(1-beta*ctho)
        dtbardx=0.d0
        dtbardc=0.d0
      else
        write(6,*)'Fatal error in get2to2a2: unknown leg'
        write(6,*)ileg
        stop
      endif
      ubar=-sbar-tbar
      dubardx=-dsbardx-dtbardx
      dubardc=-dsbardc-dtbardc
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
c End of utility routines for Bjorken x's
c
c
c
c
c Begin of zeta-subtraction routines
c
c
      function bsubint(y)
      implicit none
      real * 8 bsubint,y
      real * 8 s,xm2,cth1,zeta
      integer iproc
c =1 for gg, 2 for qq
      common/kinparms/s,xm2,cth1,zeta,iproc
      real * 8 rp1p2,rk1k2,rp1k1,rp1k2,rp2k1,rp2k2,rk1k1,rk2k2
      real * 8 sth1, rho,b,t,u,m,m1,m2,m12,mmm
      real * 8 tf,ca,cf,da,df,bf
      parameter (tf=0.5d0,ca=3d0,cf=4d0/3d0,da=8d0,df=3d0,bf=5d0/12d0)
c 
      sth1=sqrt(1-cth1**2)
      rho=4*xm2/s
      b=sqrt(1-rho)
c 2.3 in MNR
      t=-s/2*(1-b*cth1)
      u=-s/2*(1+b*cth1)
c Eikonal factor integrated in d th2/(2 pi)
      rp1p2 = -8/(s*(y-1)*(y+1))
      rk1k2 =(4/(s*(b*cth1*y+1)*SQRT(1-b**2*sth1**2*(1-y**2)/(b*cth1*y+1
     1   )**2))-4/(s*(b*cth1*y-1)*SQRT(1-b**2*sth1**2*(1-y**2)/(b*cth1*y
     2   -1)**2)))*(1-2*xm2/s)
      rp1k1 = -4*(b*cth1-1)/(s*(y-1)*(b*cth1*y-1)*SQRT(1-b**2*sth1**2*(1
     1   -y**2)/(b*cth1*y-1)**2))
      rp1k2 = -4*(b*cth1+1)/(s*(y-1)*(b*cth1*y+1)*SQRT(1-b**2*sth1**2*(1
     1   -y**2)/(b*cth1*y+1)**2))
      rp2k1 = -4*(b*cth1+1)/(s*(y+1)*(b*cth1*y-1)*SQRT(1-b**2*sth1**2*(1
     1   -y**2)/(b*cth1*y-1)**2))
      rp2k2 = -4*(b*cth1-1)/(s*(y+1)*(b*cth1*y+1)*SQRT(1-b**2*sth1**2*(1
     1   -y**2)/(b*cth1*y+1)**2))
      rk1k1 = 16*xm2*(1-b**2*sth1**2*(1-y**2)/(b*cth1*y-1)**2)**((-3d0)/
     1   2d0)/(s**2*(b*cth1*y-1)**2)
      rk2k2 = 16*xm2*(1-b**2*sth1**2*(1-y**2)/(b*cth1*y+1)**2)**((-3d0)/
     1   2d0)/(s**2*(b*cth1*y+1)**2)
c
      if(iproc.eq.1) then
c gg case;
c A.13 and A.14 in MNR
      m1=1/(8*s)*(8*t*(t**2+u**2)/(s**2*u)+8*rho*t/u-2*rho**2*s**2/u**2)
      m2=1/(8*s)*(8*u*(t**2+u**2)/(s**2*t)+8*rho*u/t-2*rho**2*s**2/t**2)
      m12=1/(8*s)*(16*(u**2+t**2)/s**2+16*rho-4*s**2*rho**2/(t*u))
c A.12 in MNR
      mmm=m1*tf/da*( (rp1k2+rp2k1)*ca*cf + (rp1k1+rp2k2)*ca*(cf-ca/2)
     # +rp1p2*ca**2/2-(rk2k2+rk1k1)*cf**2+2*rk1k2*(cf-ca/2)**2   )
     # +m2*tf/da*( (rp2k2+rp1k1)*ca*cf + (rp2k1+rp1k2)*ca*(cf-ca/2)
     # + rp1p2*ca**2/2-(rk2k2+rk1k1)*cf**2+2*rk1k2*(cf-ca/2)**2)
     # +m12*tf/da*(cf-ca/2)*( (rp1k2+rp2k2+rp1k1+rp2k1)*ca
     # -(rk2k2+rk1k1)*cf+2*rk1k2*(cf-ca))
      elseif(iproc.eq.2) then
c qq case
c A.22
      m=1/(2*s)*cf*tf/df*(-4*t*u/s**2+rho+2)
c A.21
      mmm=m*(-2*(rp1k2+rp2k1)*(-ca/4+bf)+2*(rp1k1+rp2k2)*(ca/4+bf)
     # -(rk1k1+rk2k2)*cf+2*(rp1p2+rk1k2)*(cf-ca/2))
      else
         write(*,*) ' bsubint: wrong iproc=',iproc
         write(*,*) ' must be 1 (gg) or 2 (qq)'
         stop
      endif
c I should multiply by 4 tk uk /(1-y^2) = s^2
      bsubint=mmm*(log(b**4/zeta)+log(1-y**2))*s**2
c the above expression, integrated in y between -ybar and ybar,
c corresponds to                
c     _
c   / y   /        /      4                  2  \
c  |      |       |  log b / zeta  + log (1-y )  |                 |
c  |  dy  | d th  |  --------------------------- | f  (x,y,th ,th )|
c  | _    |     2 |                2             |  gg/qq    1   2 |x=1
c / -y   /         \          1 - y             /                      
c      
c   up to a factor of: 2 pi g^6
c (2 pi missing from theta 2 integration)
c
      end


      function bsub(as,axm2,acth1,azeta,iiproc)
      implicit none
      real * 8 bsub,as,axm2,acth1,azeta
      integer iiproc
      real * 8 s,xm2,cth1,zeta
      integer iproc
      common/kinparms/s,xm2,cth1,zeta,iproc
      real * 8 yb,b
      real * 8 dgaussfwn,bsubint
      external bsubint
c
      if(iiproc.eq.3)then
        bsub=0.d0
      else
        s=as
        xm2=axm2
        cth1=acth1
        zeta=azeta
        iproc=iiproc
        b=sqrt(1-4*xm2/s)
        if(1-zeta/b**4.lt.0) then
          bsub=0
        else
          yb=sqrt(1-zeta/b**4)
          bsub=dgaussfwn(bsubint,-yb,yb,1.d-6)
        endif
      endif
      end


c This is dgauss; the name has been changed to avoid potential conflicts
      function dgaussfwn(f,a,b,eps)
c.----------------------------------------------------------------------
c.
c.    gauss integral of the function f in interval a,b
c.    last update: 10/04/88
c.
c.----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension w(12),x(12)
      external f
      data const/1.e-12/
      data w
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      data x
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
c--
c--   initialise
      delta=const*abs(a-b)
      dgaussfwn=0.
      aa=a
c--
c--   iteration loop
   10 y=b-aa
c--
c--   epsilon reached ??
      if (abs(y).le.delta) return
   20 bb=aa+y
      c1=0.5*(aa+bb)
      c2=c1-aa
      s8=0.
      s16=0.
      do 30 i=1,4
        u=x(i)*c2
   30 s8=s8+w(i)*(f(c1+u)+f(c1-u))
      do 40 i=5,12
        u=x(i)*c2
   40 s16=s16+w(i)*(f(c1+u)+f(c1-u))
      s8=s8*c2
      s16=s16*c2
      if (abs(s16-s8).gt.eps*(1.0+abs(s16))) goto 50
      dgaussfwn=dgaussfwn+s16
      aa=bb
      goto 10
   50 y=0.5*y
      if (abs(y).le.delta) write(6,9040)
      goto 20
9040  format(1H ,'**** DGAUSS: Too high Accuracy required !!     ****')
      end
c
c
c End of zeta-subtraction routines
c
c
c
c
c Begin of planar amplitudes squared
c
c
C----------------------------------------------------------------------
      SUBROUTINE qqbplanar(s,tk,uk,q1,q2,s2,q1c,q2c,w1,w2,m,t1r,t2r)
C----------------------------------------------------------------------
      IMPLICIT none
      REAL * 8 s,tk,uk,q1,q2,s2,q1c,q2c,w1,w2,m,t1r,t2r
c Planar invariant squared amplitudes, summed over spins and color 
c for the planar q qbar -> Q Qbar g process; misses a factor Ncolor^4 
c k1: quark momentum 
c k2: anti-quark momentum 
c p1: incoming light quark momentum 
c p2: incoming light antiquark momentum 
c k: radiated gluon momentum 
c 
c return values: 
c t1r: planar configuration p1,k,k1 
c t2r: planar configuration p2,k,k2 
c
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
c
c    d(p1,p1)=0 
c    d(p2,p2)=0 
c    d(k,k)=0 
c    d(k1,k1)=m2 
c    d(k2,k2)=m2 
c    d(p1,p2)=s/2 
c    d(p1,k)=-tk/2 
c    d(p2,k)=-uk/2 
c    d(p1,k1)=-q1/2 
c    d(p2,k2)=-q2/2 
c    d(k1,k2)=(s2-2*m2)/2 
c    d(p1,k2)=-q1c/2 
c    d(p2,k1)=-q2c/2 
c    d(k,k1)=w1/2 
c    d(k,k2)=w2/2 
c 
c  Relations among invariants 
c 
c          s2=s+tk+uk 
c         q1c=-s-tk-q1 
c         q2c=-s-uk-q2 
c         w1=-q1+q2-tk 
c         w2=q1-q2-uk 
      T1R = -(4*S*S2*W1**4+(8*S*S2*TK+(-4*S**2-12*Q2*S)*S2)*W1**3+(8*M**
     1   2*TK**3+((8*S-16*M**2)*S2+16*M**2*S)*TK**2+((8*M**2-4*S)*S2**2+
     2   ((-20*Q2-12*M**2)*S-4*S**2)*S2+12*M**2*S**2)*TK+2*S*S2**3+(4*Q2
     3   +4*M**2)*S*S2**2+(2*S**3+(8*Q2+4*M**2)*S**2+16*Q2**2*S)*S2+4*M*
     4   *2*S**3)*W1**2+(4*S*S2*TK**3+(((8*M**2-12*Q2)*S+16*M**2*Q2)*S2-
     5   4*S*S2**2)*TK**2+(2*S*S2**3+(8*Q2*S-16*M**2*Q2)*S2**2+(2*S**3+(
     6   4*Q2+12*M**2)*S**2+(16*Q2**2+16*M**2*Q2)*S)*S2)*TK-2*Q2*S*S2**3
     7   +(-4*Q2**2-4*M**2*Q2)*S*S2**2+(-2*Q2*S**3+(-4*Q2**2-4*M**2*Q2)*
     8   S**2-8*Q2**3*S)*S2)*W1+(4*M**2*S**2+(8*M**2*Q2+8*M**4)*S+8*M**2
     9   *Q2**2)*S2**2*TK)
     :   *tk*(uk*q1*q2*q1c*q2c*w2)**2
c     :   /(S**2*S2**2*TK*W1**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      T2R = -(4*S*S2*W2**4+(8*S*S2*UK+(-4*S**2-12*Q1*S)*S2)*W2**3+(8*M**
     1   2*UK**3+((8*S-16*M**2)*S2+16*M**2*S)*UK**2+((8*M**2-4*S)*S2**2+
     2   ((-20*Q1-12*M**2)*S-4*S**2)*S2+12*M**2*S**2)*UK+2*S*S2**3+(4*Q1
     3   +4*M**2)*S*S2**2+(2*S**3+(8*Q1+4*M**2)*S**2+16*Q1**2*S)*S2+4*M*
     4   *2*S**3)*W2**2+(4*S*S2*UK**3+(((8*M**2-12*Q1)*S+16*M**2*Q1)*S2-
     5   4*S*S2**2)*UK**2+(2*S*S2**3+(8*Q1*S-16*M**2*Q1)*S2**2+(2*S**3+(
     6   4*Q1+12*M**2)*S**2+(16*Q1**2+16*M**2*Q1)*S)*S2)*UK-2*Q1*S*S2**3
     7   +(-4*Q1**2-4*M**2*Q1)*S*S2**2+(-2*Q1*S**3+(-4*Q1**2-4*M**2*Q1)*
     8   S**2-8*Q1**3*S)*S2)*W2+(4*M**2*S**2+(8*M**2*Q1+8*M**4)*S+8*M**2
     9   *Q1**2)*S2**2*UK)
     :   *uk*(tk*q1*q2*q1c*q2c*w1)**2
c     :   /(S**2*S2**2*UK*W2**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      END
C----------------------------------------------------------------------
      SUBROUTINE ggplanar(s,tk,uk,q1,q2,s2,q1c,q2c,w1,w2,m,tr12,t1r2,t12
     1   r,tr21,t2r1,t21r)
      IMPLICIT none
      REAL * 8 s,tk,uk,q1,q2,s2,q1c,q2c,w1,w2,m,tr12,t1r2,t12r,tr21,t2r1
     1   ,t21r
c Planar invariant squared amplitudes, summed over spins and color   
c for the planar g g -> Q Qbar g process; misses a factor Ncolor^4   
c k1: quark moment   
c k2: anti-quark moment   
c p1: incoming gluon moment   
c p2: incoming gluon moment   
c k: radiated gluon moment   
c   
c Planar configurations are obtained with the gluon all on the same
c side of the fermion line; a given colour flow is specified by
c the ordering of the attachment of the gluons on the fermion
c line. For example: k1,k,p1,p2,k2 is the color structure
c                             
c    k1 -<----||--<--  ---<---  ---<--- k2
c             ||     ||       ||       
c             ^V     ^V       ^V       
c             ||     ||       ||       
c             ||     ||       ||       
c             k      p1       p2
c
c
c return values:   
c tr12: planar configuration k1,k,p1,p2,k2   
c t1r2:                      k1,p1,k,p2,k2   
c t12r:                      k1,p1,p2,k,k2   
c tr21:                      k1,k,p2,p1,k2   
c t2r1:                      k1,p2,k,p1,k2   
c t21r:                      k1,p2,p1,k,k2   
c txyz: xyz stand for incoming gluon 1, 2, and radiated gluon (r)   
c
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
c   
c    d(p1,p1)=0   
c    d(p2,p2)=0   
c    d(k,k)=0   
c    d(k1,k1)=m2   
c    d(k2,k2)=m2   
c    d(p1,p2)=s/2   
c    d(p1,k)=-tk/2   
c    d(p2,k)=-uk/2   
c    d(p1,k1)=-q1/2   
c    d(p2,k2)=-q2/2   
c    d(k1,k2)=(s2-2*m2)/2   
c    d(p1,k2)=-q1c/2   
c    d(p2,k1)=-q2c/2   
c    d(k,k1)=w1/2   
c    d(k,k2)=w2/2   
c   
c  Relations among invariants   
c   
c          s2=s+tk+uk   

c         q1c=-s-tk-q1   
c         q2c=-s-uk-q2   
c         w1=-q1+q2-tk   
c         w2=q1-q2-uk   
      TR12 = -(8*Q2*S*S2*TK*W1**5+(8*Q2*S*S2*TK**2+(8*Q2*S*S2**2+(-16*Q2
     1   *S**2-16*Q2**2*S)*S2)*TK-8*M**2*S**2*S2**2)*W1**4+(6*Q2*S*S2*TK
     2   **3+((8*M**2*Q2-24*Q2**2)*S-12*Q2*S**2)*S2*TK**2+(6*Q2*S*S2**3+
     3   ((-12*Q2-16*M**2)*S**2+8*M**2*Q2*S)*S2**2+(12*Q2*S**3+(24*Q2**2
     4   -16*M**2*Q2)*S**2+24*Q2**3*S)*S2)*TK+16*M**2*Q2*S**2*S2**2-16*M
     5   **2*Q2*S**3*S2)*W1**3+((2*Q2*S*S2-8*M**2*Q2**2)*TK**4+((-6*Q2*S
     6   **2+(8*M**2*Q2-12*Q2**2)*S+16*M**2*Q2**2)*S2-16*M**2*Q2**2*S)*T
     7   K**3+((-12*M**2*S**2-8*M**2*Q2**2)*S2**2+(6*Q2*S**3+(24*Q2**2-1
     8   6*M**2*Q2)*S**2+(24*Q2**3-8*M**2*Q2**2)*S)*S2-24*M**2*Q2**2*S**
     9   2)*TK**2+(2*Q2*S*S2**4+(8*M**2*Q2*S-6*Q2*S**2)*S2**3+(6*Q2*S**3
     :   +(-8*M**2*Q2-16*M**4)*S**2)*S2**2+(-4*Q2*S**4+(-12*Q2**2-16*M**
     ;   2*Q2)*S**3+(-24*Q2**3-8*M**2*Q2**2)*S**2-16*Q2**4*S)*S2-16*M**2
     <   *Q2**2*S**3)*TK-8*M**2*Q2**2*S**2*S2**2+16*M**2*Q2**2*S**3*S2-8
     =   *M**2*Q2**2*S**4)*W1**2+((-2*Q2*S**2-4*Q2**2*S)*S2*TK**4+((6*Q2
     >   **2*S-4*M**2*S**2)*S2**2+((6*Q2**2-4*M**2*Q2)*S**2+(12*Q2**3-16
     ?   *M**2*Q2**2+16*M**4*Q2)*S-16*M**2*Q2**3)*S2)*TK**3+(-6*Q2**2*S*
     @   S2**3+((-12*M**2*Q2-16*M**4)*S**2+(-12*Q2**3-8*M**2*Q2**2-16*M*
     1   *4*Q2)*S+16*M**2*Q2**3)*S2**2+(-2*Q2*S**4+(-6*Q2**2-4*M**2*Q2)*
     2   S**3+(-12*Q2**3-16*M**2*Q2**2+16*M**4*Q2)*S**2+(-16*Q2**4-16*M*
     3   *2*Q2**3)*S)*S2)*TK**2+(2*Q2**2*S*S2**4+(6*Q2**3+8*M**2*Q2**2+1
     4   6*M**4*Q2)*S*S2**3+((8*Q2**4+8*M**2*Q2**3)*S-16*M**4*Q2*S**2)*S
     5   2**2+(2*Q2**2*S**4+(6*Q2**3+8*M**2*Q2**2+16*M**4*Q2)*S**3+(8*Q2
     6   **4+8*M**2*Q2**3)*S**2+8*Q2**5*S)*S2)*TK)*W1+(-4*M**2*Q2*S**3+(
     7   -12*M**2*Q2**2-16*M**4*Q2-16*M**6)*S**2+(-16*M**2*Q2**3-16*M**4
     8   *Q2**2)*S-8*M**2*Q2**4)*S2**2*TK**2)
     9   *(uk*q1*q1c*q2c*w2)**2
c     9   /(Q2**2*S**2*S2**2*TK**2*W1**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      T1R2 = ((2*Q1*Q2*S2*TK**2+(4*Q1**2*Q2-2*Q1*Q2**2)*S2*TK+8*M**2*Q1*
     1   *2*Q2**2)*UK**4+((4*M**2*Q2*S2**2+(6*Q1*Q2**2+(4*M**2*Q1-6*Q1**
     2   2)*Q2)*S2)*TK**2+(-6*Q1**2*Q2*S2**2+(-6*Q1*Q2**3+(12*Q1**2-8*M*
     3   *2*Q1)*Q2**2+(-12*Q1**3+16*M**2*Q1**2-16*M**4*Q1)*Q2)*S2+16*M**
     4   2*Q1**2*Q2**2)*TK+(16*M**2*Q1**3*Q2-16*M**2*Q1**2*Q2**2)*S2)*UK
     5   **3+(2*Q1*Q2*S2*TK**4+(4*M**2*Q1*S2**2+((6*Q1**2+4*M**2*Q1)*Q2-
     6   6*Q1*Q2**2)*S2)*TK**3+((12*M**2*Q2**2+(12*M**2*Q1+16*M**4)*Q2+1
     7   2*M**2*Q1**2+16*M**4*Q1+16*M**6)*S2**2+(12*Q1*Q2**3+(16*M**2*Q1
     8   -24*Q1**2)*Q2**2+(12*Q1**3+16*M**2*Q1**2-16*M**4*Q1)*Q2)*S2+24*
     9   M**2*Q1**2*Q2**2)*TK**2+(6*Q1**2*Q2*S2**3+((12*Q1**3+8*M**2*Q1*
     :   *2+16*M**4*Q1)*Q2+16*M**2*Q1**3+16*M**4*Q1**2)*S2**2+(-8*Q1*Q2*
     ;   *4+(24*Q1**2-8*M**2*Q1)*Q2**3+(8*M**2*Q1**2-24*Q1**3)*Q2**2+(16
     <   *Q1**4+16*M**2*Q1**3)*Q2)*S2)*TK+(8*M**2*Q1**2*Q2**2-16*M**2*Q1
     =   **3*Q2+8*M**2*Q1**4)*S2**2)*UK**2+((4*Q1*Q2**2-2*Q1**2*Q2)*S2*T
     >   K**4+(-6*Q1*Q2**2*S2**2+(-12*Q1*Q2**3+(12*Q1**2+16*M**2*Q1)*Q2*
     ?   *2+(-6*Q1**3-8*M**2*Q1**2-16*M**4*Q1)*Q2)*S2+16*M**2*Q1**2*Q2**
     @   2)*TK**3+(6*Q1*Q2**2*S2**3+((12*Q1+16*M**2)*Q2**3+(8*M**2*Q1+16
     1   *M**4)*Q2**2+16*M**4*Q1*Q2)*S2**2+(16*Q1*Q2**4+(16*M**2*Q1-24*Q
     2   1**2)*Q2**3+(24*Q1**3+8*M**2*Q1**2)*Q2**2+(-8*Q1**4-8*M**2*Q1**
     3   3)*Q2)*S2)*TK**2+((-2*Q1*Q2**2-2*Q1**2*Q2)*S2**4+(-6*Q1*Q2**3-8
     4   *M**2*Q1*Q2**2+(-6*Q1**3-8*M**2*Q1**2-16*M**4*Q1)*Q2)*S2**3+(-8
     5   *Q1*Q2**4-8*M**2*Q1*Q2**3+(-8*Q1**4-8*M**2*Q1**3)*Q2)*S2**2+(-8
     6   *Q1*Q2**5+16*Q1**2*Q2**4-24*Q1**3*Q2**3+16*Q1**4*Q2**2-8*Q1**5*
     7   Q2)*S2)*TK)*UK+8*M**2*Q1**2*Q2**2*TK**4+(16*M**2*Q1*Q2**3-16*M*
     8   *2*Q1**2*Q2**2)*S2*TK**3+(8*M**2*Q2**4-16*M**2*Q1*Q2**3+8*M**2*
     9   Q1**2*Q2**2)*S2**2*TK**2)
     :   *(s*q1c*q2c*w1*w2)**2
c     :   /(Q1**2*Q2**2*S2**2*TK**2*UK**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      T12R = -(8*Q1*S*S2*UK*W2**5+(8*Q1*S*S2*UK**2+(8*Q1*S*S2**2+(-16*Q1
     1   *S**2-16*Q1**2*S)*S2)*UK-8*M**2*S**2*S2**2)*W2**4+(6*Q1*S*S2*UK
     2   **3+((8*M**2*Q1-24*Q1**2)*S-12*Q1*S**2)*S2*UK**2+(6*Q1*S*S2**3+
     3   ((-12*Q1-16*M**2)*S**2+8*M**2*Q1*S)*S2**2+(12*Q1*S**3+(24*Q1**2
     4   -16*M**2*Q1)*S**2+24*Q1**3*S)*S2)*UK+16*M**2*Q1*S**2*S2**2-16*M
     5   **2*Q1*S**3*S2)*W2**3+((2*Q1*S*S2-8*M**2*Q1**2)*UK**4+((-6*Q1*S
     6   **2+(8*M**2*Q1-12*Q1**2)*S+16*M**2*Q1**2)*S2-16*M**2*Q1**2*S)*U
     7   K**3+((-12*M**2*S**2-8*M**2*Q1**2)*S2**2+(6*Q1*S**3+(24*Q1**2-1
     8   6*M**2*Q1)*S**2+(24*Q1**3-8*M**2*Q1**2)*S)*S2-24*M**2*Q1**2*S**
     9   2)*UK**2+(2*Q1*S*S2**4+(8*M**2*Q1*S-6*Q1*S**2)*S2**3+(6*Q1*S**3
     :   +(-8*M**2*Q1-16*M**4)*S**2)*S2**2+(-4*Q1*S**4+(-12*Q1**2-16*M**
     ;   2*Q1)*S**3+(-24*Q1**3-8*M**2*Q1**2)*S**2-16*Q1**4*S)*S2-16*M**2
     <   *Q1**2*S**3)*UK-8*M**2*Q1**2*S**2*S2**2+16*M**2*Q1**2*S**3*S2-8
     =   *M**2*Q1**2*S**4)*W2**2+((-2*Q1*S**2-4*Q1**2*S)*S2*UK**4+((6*Q1
     >   **2*S-4*M**2*S**2)*S2**2+((6*Q1**2-4*M**2*Q1)*S**2+(12*Q1**3-16
     ?   *M**2*Q1**2+16*M**4*Q1)*S-16*M**2*Q1**3)*S2)*UK**3+(-6*Q1**2*S*
     @   S2**3+((-12*M**2*Q1-16*M**4)*S**2+(-12*Q1**3-8*M**2*Q1**2-16*M*
     1   *4*Q1)*S+16*M**2*Q1**3)*S2**2+(-2*Q1*S**4+(-6*Q1**2-4*M**2*Q1)*
     2   S**3+(-12*Q1**3-16*M**2*Q1**2+16*M**4*Q1)*S**2+(-16*Q1**4-16*M*
     3   *2*Q1**3)*S)*S2)*UK**2+(2*Q1**2*S*S2**4+(6*Q1**3+8*M**2*Q1**2+1
     4   6*M**4*Q1)*S*S2**3+((8*Q1**4+8*M**2*Q1**3)*S-16*M**4*Q1*S**2)*S
     5   2**2+(2*Q1**2*S**4+(6*Q1**3+8*M**2*Q1**2+16*M**4*Q1)*S**3+(8*Q1
     6   **4+8*M**2*Q1**3)*S**2+8*Q1**5*S)*S2)*UK)*W2+(-4*M**2*Q1*S**3+(
     7   -12*M**2*Q1**2-16*M**4*Q1-16*M**6)*S**2+(-16*M**2*Q1**3-16*M**4
     8   *Q1**2)*S-8*M**2*Q1**4)*S2**2*UK**2)
     9   *(tk*q2*q1c*q2c*w1)**2
c     9   /(Q1**2*S**2*S2**2*UK**2*W2**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      TR21 = -(8*Q1C*S*S2*UK*W1**5+(8*Q1C*S*S2*UK**2+(8*Q1C*S*S2**2+(-16
     1   *Q1C*S**2-16*Q1C**2*S)*S2)*UK-8*M**2*S**2*S2**2)*W1**4+(6*Q1C*S
     2   *S2*UK**3+((8*M**2*Q1C-24*Q1C**2)*S-12*Q1C*S**2)*S2*UK**2+(6*Q1
     3   C*S*S2**3+((-12*Q1C-16*M**2)*S**2+8*M**2*Q1C*S)*S2**2+(12*Q1C*S
     4   **3+(24*Q1C**2-16*M**2*Q1C)*S**2+24*Q1C**3*S)*S2)*UK+16*M**2*Q1
     5   C*S**2*S2**2-16*M**2*Q1C*S**3*S2)*W1**3+((2*Q1C*S*S2-8*M**2*Q1C
     6   **2)*UK**4+((-6*Q1C*S**2+(8*M**2*Q1C-12*Q1C**2)*S+16*M**2*Q1C**
     7   2)*S2-16*M**2*Q1C**2*S)*UK**3+((-12*M**2*S**2-8*M**2*Q1C**2)*S2
     8   **2+(6*Q1C*S**3+(24*Q1C**2-16*M**2*Q1C)*S**2+(24*Q1C**3-8*M**2*
     9   Q1C**2)*S)*S2-24*M**2*Q1C**2*S**2)*UK**2+(2*Q1C*S*S2**4+(8*M**2
     :   *Q1C*S-6*Q1C*S**2)*S2**3+(6*Q1C*S**3+(-8*M**2*Q1C-16*M**4)*S**2
     ;   )*S2**2+(-4*Q1C*S**4+(-12*Q1C**2-16*M**2*Q1C)*S**3+(-24*Q1C**3-
     <   8*M**2*Q1C**2)*S**2-16*Q1C**4*S)*S2-16*M**2*Q1C**2*S**3)*UK-8*M
     =   **2*Q1C**2*S**2*S2**2+16*M**2*Q1C**2*S**3*S2-8*M**2*Q1C**2*S**4
     >   )*W1**2+((-2*Q1C*S**2-4*Q1C**2*S)*S2*UK**4+((6*Q1C**2*S-4*M**2*
     ?   S**2)*S2**2+((6*Q1C**2-4*M**2*Q1C)*S**2+(12*Q1C**3-16*M**2*Q1C*
     @   *2+16*M**4*Q1C)*S-16*M**2*Q1C**3)*S2)*UK**3+(-6*Q1C**2*S*S2**3+
     1   ((-12*M**2*Q1C-16*M**4)*S**2+(-12*Q1C**3-8*M**2*Q1C**2-16*M**4*
     2   Q1C)*S+16*M**2*Q1C**3)*S2**2+(-2*Q1C*S**4+(-6*Q1C**2-4*M**2*Q1C
     3   )*S**3+(-12*Q1C**3-16*M**2*Q1C**2+16*M**4*Q1C)*S**2+(-16*Q1C**4
     4   -16*M**2*Q1C**3)*S)*S2)*UK**2+(2*Q1C**2*S*S2**4+(6*Q1C**3+8*M**
     5   2*Q1C**2+16*M**4*Q1C)*S*S2**3+((8*Q1C**4+8*M**2*Q1C**3)*S-16*M*
     6   *4*Q1C*S**2)*S2**2+(2*Q1C**2*S**4+(6*Q1C**3+8*M**2*Q1C**2+16*M*
     7   *4*Q1C)*S**3+(8*Q1C**4+8*M**2*Q1C**3)*S**2+8*Q1C**5*S)*S2)*UK)*
     8   W1+(-4*M**2*Q1C*S**3+(-12*M**2*Q1C**2-16*M**4*Q1C-16*M**6)*S**2
     9   +(-16*M**2*Q1C**3-16*M**4*Q1C**2)*S-8*M**2*Q1C**4)*S2**2*UK**2)
     :   *(tk*q1*q2*q2c*w2)**2
c     :   /(Q1C**2*S**2*S2**2*UK**2*W1**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      T2R1 = ((2*Q1C*Q2C*S2*TK**2+(4*Q1C**2*Q2C-2*Q1C*Q2C**2)*S2*TK+8*M*
     1   *2*Q1C**2*Q2C**2)*UK**4+((4*M**2*Q2C*S2**2+(6*Q1C*Q2C**2+(4*M**
     2   2*Q1C-6*Q1C**2)*Q2C)*S2)*TK**2+(-6*Q1C**2*Q2C*S2**2+(-6*Q1C*Q2C
     3   **3+(12*Q1C**2-8*M**2*Q1C)*Q2C**2+(-12*Q1C**3+16*M**2*Q1C**2-16
     4   *M**4*Q1C)*Q2C)*S2+16*M**2*Q1C**2*Q2C**2)*TK+(16*M**2*Q1C**3*Q2
     5   C-16*M**2*Q1C**2*Q2C**2)*S2)*UK**3+(2*Q1C*Q2C*S2*TK**4+(4*M**2*
     6   Q1C*S2**2+((6*Q1C**2+4*M**2*Q1C)*Q2C-6*Q1C*Q2C**2)*S2)*TK**3+((
     7   12*M**2*Q2C**2+(12*M**2*Q1C+16*M**4)*Q2C+12*M**2*Q1C**2+16*M**4
     8   *Q1C+16*M**6)*S2**2+(12*Q1C*Q2C**3+(16*M**2*Q1C-24*Q1C**2)*Q2C*
     9   *2+(12*Q1C**3+16*M**2*Q1C**2-16*M**4*Q1C)*Q2C)*S2+24*M**2*Q1C**
     :   2*Q2C**2)*TK**2+(6*Q1C**2*Q2C*S2**3+((12*Q1C**3+8*M**2*Q1C**2+1
     ;   6*M**4*Q1C)*Q2C+16*M**2*Q1C**3+16*M**4*Q1C**2)*S2**2+(-8*Q1C*Q2
     <   C**4+(24*Q1C**2-8*M**2*Q1C)*Q2C**3+(8*M**2*Q1C**2-24*Q1C**3)*Q2
     =   C**2+(16*Q1C**4+16*M**2*Q1C**3)*Q2C)*S2)*TK+(8*M**2*Q1C**2*Q2C*
     >   *2-16*M**2*Q1C**3*Q2C+8*M**2*Q1C**4)*S2**2)*UK**2+((4*Q1C*Q2C**
     ?   2-2*Q1C**2*Q2C)*S2*TK**4+(-6*Q1C*Q2C**2*S2**2+(-12*Q1C*Q2C**3+(
     @   12*Q1C**2+16*M**2*Q1C)*Q2C**2+(-6*Q1C**3-8*M**2*Q1C**2-16*M**4*
     1   Q1C)*Q2C)*S2+16*M**2*Q1C**2*Q2C**2)*TK**3+(6*Q1C*Q2C**2*S2**3+(
     2   (12*Q1C+16*M**2)*Q2C**3+(8*M**2*Q1C+16*M**4)*Q2C**2+16*M**4*Q1C
     3   *Q2C)*S2**2+(16*Q1C*Q2C**4+(16*M**2*Q1C-24*Q1C**2)*Q2C**3+(24*Q
     4   1C**3+8*M**2*Q1C**2)*Q2C**2+(-8*Q1C**4-8*M**2*Q1C**3)*Q2C)*S2)*
     5   TK**2+((-2*Q1C*Q2C**2-2*Q1C**2*Q2C)*S2**4+(-6*Q1C*Q2C**3-8*M**2
     6   *Q1C*Q2C**2+(-6*Q1C**3-8*M**2*Q1C**2-16*M**4*Q1C)*Q2C)*S2**3+(-
     7   8*Q1C*Q2C**4-8*M**2*Q1C*Q2C**3+(-8*Q1C**4-8*M**2*Q1C**3)*Q2C)*S
     8   2**2+(-8*Q1C*Q2C**5+16*Q1C**2*Q2C**4-24*Q1C**3*Q2C**3+16*Q1C**4
     9   *Q2C**2-8*Q1C**5*Q2C)*S2)*TK)*UK+8*M**2*Q1C**2*Q2C**2*TK**4+(16
     :   *M**2*Q1C*Q2C**3-16*M**2*Q1C**2*Q2C**2)*S2*TK**3+(8*M**2*Q2C**4
     ;   -16*M**2*Q1C*Q2C**3+8*M**2*Q1C**2*Q2C**2)*S2**2*TK**2)
     <    *(s*q1*q2*w1*w2)**2
c     <   /(Q1C**2*Q2C**2*S2**2*TK**2*UK**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      T21R = -(8*Q2C*S*S2*TK*W2**5+(8*Q2C*S*S2*TK**2+(8*Q2C*S*S2**2+(-16
     1   *Q2C*S**2-16*Q2C**2*S)*S2)*TK-8*M**2*S**2*S2**2)*W2**4+(6*Q2C*S
     2   *S2*TK**3+((8*M**2*Q2C-24*Q2C**2)*S-12*Q2C*S**2)*S2*TK**2+(6*Q2
     3   C*S*S2**3+((-12*Q2C-16*M**2)*S**2+8*M**2*Q2C*S)*S2**2+(12*Q2C*S
     4   **3+(24*Q2C**2-16*M**2*Q2C)*S**2+24*Q2C**3*S)*S2)*TK+16*M**2*Q2
     5   C*S**2*S2**2-16*M**2*Q2C*S**3*S2)*W2**3+((2*Q2C*S*S2-8*M**2*Q2C
     6   **2)*TK**4+((-6*Q2C*S**2+(8*M**2*Q2C-12*Q2C**2)*S+16*M**2*Q2C**
     7   2)*S2-16*M**2*Q2C**2*S)*TK**3+((-12*M**2*S**2-8*M**2*Q2C**2)*S2
     8   **2+(6*Q2C*S**3+(24*Q2C**2-16*M**2*Q2C)*S**2+(24*Q2C**3-8*M**2*
     9   Q2C**2)*S)*S2-24*M**2*Q2C**2*S**2)*TK**2+(2*Q2C*S*S2**4+(8*M**2
     :   *Q2C*S-6*Q2C*S**2)*S2**3+(6*Q2C*S**3+(-8*M**2*Q2C-16*M**4)*S**2
     ;   )*S2**2+(-4*Q2C*S**4+(-12*Q2C**2-16*M**2*Q2C)*S**3+(-24*Q2C**3-
     <   8*M**2*Q2C**2)*S**2-16*Q2C**4*S)*S2-16*M**2*Q2C**2*S**3)*TK-8*M
     =   **2*Q2C**2*S**2*S2**2+16*M**2*Q2C**2*S**3*S2-8*M**2*Q2C**2*S**4
     >   )*W2**2+((-2*Q2C*S**2-4*Q2C**2*S)*S2*TK**4+((6*Q2C**2*S-4*M**2*
     ?   S**2)*S2**2+((6*Q2C**2-4*M**2*Q2C)*S**2+(12*Q2C**3-16*M**2*Q2C*
     @   *2+16*M**4*Q2C)*S-16*M**2*Q2C**3)*S2)*TK**3+(-6*Q2C**2*S*S2**3+
     1   ((-12*M**2*Q2C-16*M**4)*S**2+(-12*Q2C**3-8*M**2*Q2C**2-16*M**4*
     2   Q2C)*S+16*M**2*Q2C**3)*S2**2+(-2*Q2C*S**4+(-6*Q2C**2-4*M**2*Q2C
     3   )*S**3+(-12*Q2C**3-16*M**2*Q2C**2+16*M**4*Q2C)*S**2+(-16*Q2C**4
     4   -16*M**2*Q2C**3)*S)*S2)*TK**2+(2*Q2C**2*S*S2**4+(6*Q2C**3+8*M**
     5   2*Q2C**2+16*M**4*Q2C)*S*S2**3+((8*Q2C**4+8*M**2*Q2C**3)*S-16*M*
     6   *4*Q2C*S**2)*S2**2+(2*Q2C**2*S**4+(6*Q2C**3+8*M**2*Q2C**2+16*M*
     7   *4*Q2C)*S**3+(8*Q2C**4+8*M**2*Q2C**3)*S**2+8*Q2C**5*S)*S2)*TK)*
     8   W2+(-4*M**2*Q2C*S**3+(-12*M**2*Q2C**2-16*M**4*Q2C-16*M**6)*S**2
     9   +(-16*M**2*Q2C**3-16*M**4*Q2C**2)*S-8*M**2*Q2C**4)*S2**2*TK**2)
     :   *(uk*q1*q2*q1c*w1)**2
c     :   /(Q2C**2*S**2*S2**2*TK**2*W2**2)
C--N.B. new definitions Ti-->Ti*(s*tk*uk*q1*q2*s2*q1c*q2c*w1*w2)**2
      END
c
c
c End of planar amplitudes squared
c
c
      function zgmu2_nlo()
c Sets the desired factorization scale and returns the strong coupling squared
c To be called is association with pure NLO terms
      implicit none
      real * 8 zgmu2_nlo
      real * 8 pi
      parameter (pi=3.14159265358979312D0)
      real * 8 pq10,pq20,pp0
      common/perpen/pq10(2),pq20(2),pp0(2)
      integer nlas
      common/cnlas/nlas
      include 'hvqcblks.h'
      real * 8 pt12,pt22,pt2,xmu2,as
      real * 8 alfas
c
      pt12= pq10(1)**2 + pq10(2)**2
      pt22= pq20(1)**2 + pq20(2)**2
      pt2 = (pt12+pt22)/2.d0
      xmu2 = pt2 + xm2
c set the factorization scales for hadron 1 and 2, and the
c renormalization scale
      xmuf2h1 = xmu2*xf2h1
      xmuf2h2 = xmu2*xf2h2
      xmur2  = xmu2*xren2
      as    = alfas(xmur2,xlam,nlas)
      zgmu2_nlo = 4.d0*pi*as
      zg = sqrt(zgmu2_nlo)
      end


      function zgmu2_mc()
c Sets the desired factorization scale and returns the strong coupling squared
c To be called is association with MC subtraction terms
      implicit none
      real * 8 zgmu2_mc
      real * 8 pi
      parameter (pi=3.14159265358979312D0)
      real * 8 pq10,pq20,pp0
      common/perpen/pq10(2),pq20(2),pp0(2)
      integer nlas
      common/cnlas/nlas
      include 'hvqcblks.h'
      real * 8 pt12,pt22,pt2,xmu2,as
      real * 8 alfas
c
      pt12= pq10(1)**2 + pq10(2)**2
      pt22= pq20(1)**2 + pq20(2)**2
      pt2 = (pt12+pt22)/2.d0
      xmu2 = pt2 + xm2
c set the factorization scales for hadron 1 and 2, and the
c renormalization scale
      xmumcf2h1 = xmu2*xf2h1mc
      xmumcf2h2 = xmu2*xf2h2mc
      xmumcr2  = xmu2*xren2mc
      as    = alfas(xmumcr2,xlam,nlas)
      zgmu2_mc = 4.d0*pi*as
      zg = sqrt(zgmu2_mc)
      end


      function zgmu6_mc(zg2_nlo)
c Computes the coupling constant for MC subtraction terms
      implicit none
      real * 8 zgmu6_mc,zg2_nlo,zgmu2_mc,zg2_mc,tmp
      integer iasmc
      common/ciasmc/iasmc
c
      zg2_mc=zgmu2_mc()
      if(iasmc.eq.1)then
        tmp=zg2_nlo**3
      elseif(iasmc.eq.2)then
        tmp=zg2_mc**3
      elseif(iasmc.eq.3)then
        tmp=zg2_mc*zg2_nlo**2
      else
        write(*,*)'Fatal error in zgmu6_mc: unknown iasmc',iasmc
        stop
      endif
      zgmu6_mc=tmp
      return
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


      subroutine HVQWARN(str)
      character *(*) str
      write(*,*) '********** WARNING **********'
      write(*,*) '*********  ',str,'  *********'
      end
