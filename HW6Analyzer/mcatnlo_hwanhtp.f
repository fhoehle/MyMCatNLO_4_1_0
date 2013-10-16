C----------------------------------------------------------------------
      SUBROUTINE RCLOS()
C     DUMMY IF HBOOK IS USED
C----------------------------------------------------------------------
      END


C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      real * 4 pi
      parameter (pi=3.14160E0)
      integer j,k,ivlep1,ivlep2,idec
      common/vvlin/ivlep1,ivlep2
      character*8 cc(2)
      data cc/' no veto',' veto=30'/
c
      if(ivlep2.ne.ivlep1)call hwwarn('HWABEG',500)
c
      if(ivlep1.eq.7)then
        idec=1
      else
        if(ivlep1.gt.3)call hwwarn('HWABEG',501)
        idec=0
      endif
c
      call inihist 
      if(idec.eq.0)then
c Spin correlations are included
        do j=1,2
          k=(j-1)*20
          call mbook(k+ 1,'l[t] pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 2,'nu[t] pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 3,'b[t] pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 4,'l[H] pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 5,'nu[H] pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 6,'l[t] y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+ 7,'nu[t] y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+ 8,'b[t] y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+ 9,'l[H] y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+10,'nu[H] y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+11,'ll pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+12,'ll dphi'//cc(j),pi/20.e0,0.e0,pi)
          call mbook(k+13,'ll m'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+14,'ll DelR'//cc(j),pi/20.e0,0.e0,3*pi)
        enddo
      elseif(idec.eq.1)then
c Spin correlations are not included
        do j=1,2
          k=(j-1)*20
          call mbook(k+ 1,'Ht pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 2,'Ht log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
          call mbook(k+ 3,'Ht inv m'//cc(j),10.e0,240.e0,800.e0)
          call mbook(k+ 4,'Ht azimt'//cc(j),pi/20.e0,0.e0,pi)
          call mbook(k+ 5,'Ht log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
          call mbook(k+ 6,'Ht del R'//cc(j),pi/20.e0,0.e0,3*pi)
          call mbook(k+ 7,'t pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 8,'H pt'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+ 9,'t y'//cc(j),0.2e0,-9.e0,9.e0)
          call mbook(k+10,'H y'//cc(j),0.2e0,-9.e0,9.e0)
        enddo
      else
        call hwwarn('HWABEG',502)
      endif
c     
      END


C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      REAL*8 XNORM
      INTEGER I,J,K,IVLEP1,IVLEP2,IDEC
      COMMON/VVLIN/IVLEP1,IVLEP2
C
      OPEN(UNIT=99,FILE='HERHT.TOP',STATUS='UNKNOWN')
C XNORM IS SUCH THAT THE CROSS SECTION PER BIN IS IN PB, SINCE THE HERWIG 
C WEIGHT IS IN NB, AND CORRESPONDS TO THE AVERAGE CROSS SECTION
      XNORM=1.D3/DFLOAT(NEVHEP)
      DO I=1,100              
 	CALL MFINAL3(I)             
        CALL MCOPY(I,I+100)
        CALL MOPERA(I+100,'F',I+100,I+100,SNGL(XNORM),0.E0)
 	CALL MFINAL3(I+100)             
      ENDDO                          
C
      if(ivlep1.eq.7)then
        idec=1
      else
        if(ivlep1.gt.3)call hwwarn('HWAEND',501)
        idec=0
      endif
c
      if(idec.eq.0)then
        do j=1,2
          k=(j-1)*20
          call multitop(100+k+ 1,99,3,2,'l[t] pt',' ','LOG')
          call multitop(100+k+ 2,99,3,2,'nu[t] pt',' ','LOG')
          call multitop(100+k+ 3,99,3,2,'b[t] pt',' ','LOG')
          call multitop(100+k+ 4,99,3,2,'l[H] pt',' ','LOG')
          call multitop(100+k+ 5,99,3,2,'nu[H] pt',' ','LOG')
          call multitop(100+k+ 6,99,3,2,'l[t] y',' ','LOG')
          call multitop(100+k+ 7,99,3,2,'nu[t] y',' ','LOG')
          call multitop(100+k+ 8,99,3,2,'b[t] y',' ','LOG')
          call multitop(100+k+ 9,99,3,2,'l[H] y',' ','LOG')
          call multitop(100+k+10,99,3,2,'nu[H] y',' ','LOG')
          call multitop(100+k+11,99,3,2,'ll pt',' ','LOG')
          call multitop(100+k+12,99,3,2,'ll dphi',' ','LOG')
          call multitop(100+k+13,99,3,2,'ll m',' ','LOG')
          call multitop(100+k+14,99,3,2,'ll DelR',' ','LOG')
        enddo
      elseif(idec.eq.1)then
        do j=1,2
          k=(j-1)*20
          call multitop(100+k+ 1,99,2,3,'Ht pt',' ','LOG')
          call multitop(100+k+ 2,99,2,3,'Ht log[pt]',' ','LOG')
          call multitop(100+k+ 3,99,2,3,'Ht inv m',' ','LOG')
          call multitop(100+k+ 4,99,2,3,'Ht azimt',' ','LOG')
          call multitop(100+k+ 5,99,2,3,'Ht log[pi-azimt]',' ','LOG')
          call multitop(100+k+ 6,99,2,3,'Ht del R',' ','LOG')
          call multitop(100+k+ 7,99,2,3,'t pt',' ','LOG')
          call multitop(100+k+ 8,99,2,3,'H pt',' ','LOG')
          call multitop(100+k+ 9,99,2,3,'t y',' ','LOG')
          call multitop(100+k+10,99,2,3,'H y',' ','LOG')
        enddo
      else
        call hwwarn('HWAEND',502)
      endif
c
      CLOSE(99)
      END

C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWVDOT,PSUM(4)
      INTEGER JPR,ICHSUM,ICHINI,IDEC,IT,IW,IB,IHEP,ID1,IHADR,I,IPT,
     # IPTN,I1,ISGNT,ISGNW,ILPT,INUT,ILPW,INUW,IBQ,ITMP,KK,JB(30),
     # JT(10),JW(10),JBMAX(2)
      LOGICAL DIDSOF
      DOUBLE PRECISION RCUT,PTB1MAX,PTB2MAX,PTB,PTCALC,ETAB,
     # GETPSEUDORAP,XPTLT,XYLT,GETRAPIDITY,XPTNUT,XYNUT,XPTBT,XYBT,
     # XPTLW,XYLW,XPTNUW,XYNUW,XPTLL,DPHILL,GETDELPHI,XMLL,GETINVM,
     # XLLDELR,GETDR,XPTT,XYT,XPTW,XYW,XPTWT,WTM,AZI,AZINORM,WTDELR,
     # PB(4,30),PPTP(4,10),PPW(4,10),PLPT(4),PNUT(4),PLPW(4),PNUW(4),
     # PBQT(4),YPPLL(4),YPPTPW(4)
      REAL*8 WWW0
      REAL*8 PI
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 BVETO
      PARAMETER (BVETO=30.d0)
      LOGICAL VETOCUT
      INTEGER IVLEP1,IVLEP2
      COMMON/VVLIN/IVLEP1,IVLEP2
c
      IF (IERROR.NE.0) RETURN
c
      IF(IPROC.GT.0)CALL HWWARN('HWANAL',500)
      JPR=MOD(ABS(IPROC),10000)/10
      IF(JPR.NE.204)CALL HWWARN('HWANAL',502)
C INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT'S A POWER-SUPPRESSED
C EFFECT, SO THROW THE EVENT AWAY
      IF(SIGN(1.D0,PHEP(3,4)).EQ.SIGN(1.D0,PHEP(3,5)))THEN
        CALL HWWARN('HWANAL',111)
        GOTO 999
      ENDIF
      WWW0=EVWGT
      CALL HWVSUM(4,PHEP(1,1),PHEP(1,2),PSUM)
      CALL HWVSCA(4,-1D0,PSUM,PSUM)
      ICHSUM=0
      ICHINI=ICHRG(IDHW(1))+ICHRG(IDHW(2))
      DIDSOF=.FALSE.
      IF(IVLEP1.EQ.7)THEN
        IDEC=1
      ELSE
        IDEC=0
      ENDIF
c pseudorapidity cut
      IF(PBEAM1.GT.2500)THEN
        RCUT=2.5D0
      ELSE
        RCUT=1.0D0
      ENDIF
c IT is the number of top quarks; their four momenta are PPTP(1-4,IT).
c Positions in event record are JT(IT)
      IT=0
c IW is the number of Ws; their four momenta are PPW(1-4,IW).
c Positions in event record are JW(IW)
      IW=0
c IB is the number of stable final-state b hadrons; their four momenta 
c are PB(1-4,IB). Positions in event record are JB(IB)
      IB=0
c
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
C---FIND FINAL STATE HADRONS AND PHOTONS
        IF (ISTHEP(IHEP).EQ.1 .AND.
     &       (ABS(IDHEP(IHEP)).GT.100.OR.IDHEP(IHEP).EQ.22)) THEN
          ID1=IHADR(IDHEP(IHEP))
          IF(ID1.EQ.5)THEN
C FOUND A B-FLAVOURED HADRON
            IB=IB+1
            JB(IB)=IHEP
            DO I=1,4
               PB(I,IB)=PHEP(I,IHEP)
            ENDDO
          ENDIF
        ENDIF
        IF(ISTHEP(IHEP).EQ.155.AND.ABS(IDHEP(IHEP)).EQ.6)THEN
C FOUND A TOP
          IT=IT+1
          JT(IT)=IHEP
          DO I=1,4
             PPTP(I,IT)=PHEP(I,IHEP)
          ENDDO
        ENDIF
C Intermediate undecayed Higgs is given status 155 (was 195 for W)
        IF( ((ISTHEP(IHEP).EQ.155.AND.IVLEP2.EQ.7).OR.
     #       (ISTHEP(IHEP).EQ.155.AND.IVLEP2.LT.7)) .AND.
     #      ABS(IDHEP(IHEP)).EQ.37 )THEN
C FOUND A H+/H-
          IW=IW+1
          JW(IW)=IHEP
          DO I=1,4
             PPW(I,IW)=PHEP(I,IHEP)
          ENDDO
        ENDIF
  100 CONTINUE
      IF(IT.EQ.0.OR.IW.EQ.0)THEN
        CALL HWUEPR
        CALL HWWARN('HWANAL',503)
      ENDIF
C CHECK MOMENTUM AND CHARGE CONSERVATION
      IF (HWVDOT(3,PSUM,PSUM).GT.1.E-4*PHEP(4,1)**2) THEN
         CALL HWUEPR
         CALL HWWARN('HWANAL',112)
         GOTO 999
      ENDIF
      IF (ICHSUM.NE.ICHINI) THEN
         CALL HWUEPR
         CALL HWWARN('HWANAL',113)
         GOTO 999
      ENDIF
C FIND THE PRIMARY TOP
      IF(IT.EQ.1)THEN
        IPT=1
      ELSE
        IPTN=0
        DO I=1,IT
          IHEP=JMOHEP(1,JT(I))
          I1=0
          DOWHILE(JMOHEP(1,IHEP).NE.6)
            IHEP=JMOHEP(1,IHEP)
            I1=I1+1
            IF(I1.GT.30)THEN
              CALL HWUEPR
              WRITE(*,*)'IHEP,I1=',IHEP,I1
              CALL HWWARN('HWANAL',515)
            ENDIF
          ENDDO
          IF(ABS(IDHEP(IHEP)).EQ.6)THEN
            IPT=I
            IPTN=IPTN+1
          ENDIF
        ENDDO
C MORE THAN ONE PRIMARY TOP FOUND: THIS MAY HAPPEN IN THE EXTREMELY 
C RARE CASE IN WHICH A GLUON WHICH BRANCHES INTO A TTBAR PAIR.
C DO NOT END RUN, JUST KILL EVENT
        IF(IPTN.NE.1)THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',116)
          GOTO 999
        ENDIF
      ENDIF
      ISGNT=SIGN(1,IDHEP(JT(IPT)))
      ISGNW=-ISGNT
C FIND THE DECAY PRODUCTS OF THE PRIMARY TOP AND W IF SPIN CORRELATIONS 
C ARE INCLUDED: ASSUME THE HARD W IS THE FIRST IN THE EVENT RECORD
      IF(IDEC.EQ.0)THEN
        ILPT=JDAHEP(1,JDAHEP(1,JDAHEP(1,JDAHEP(1,JDAHEP(1,JT(IPT))))))
        INUT=JDAHEP(1,JDAHEP(2,JDAHEP(1,JDAHEP(1,JDAHEP(1,JT(IPT))))))
        ILPW=JDAHEP(1,JDAHEP(1,JW(1)))
        INUW=JDAHEP(1,JDAHEP(2,JW(1)))
c IBQ is the position in the event record of the b quark, emerging from
c top decay, before hadronization. Consider only the primary top here
        IBQ=0
        IF(JDAHEP(2,JT(IPT))-JDAHEP(1,JT(IPT)).EQ.1)THEN
          DO I1=JDAHEP(1,JDAHEP(1,JDAHEP(2,JT(IPT)))),
     #          JDAHEP(2,JDAHEP(1,JDAHEP(2,JT(IPT))))
            IF(ISTHEP(I1).EQ.2.AND.ABS(IDHEP(I1)).EQ.5)IBQ=I1
          ENDDO
        ELSEIF(JDAHEP(2,JT(IPT))-JDAHEP(1,JT(IPT)).EQ.2)THEN
          DO I1=JDAHEP(1,JDAHEP(1,JDAHEP(2,JT(IPT))-1)),
     #          JDAHEP(2,JDAHEP(1,JDAHEP(2,JT(IPT))-1))
            IF(ISTHEP(I1).EQ.2.AND.ABS(IDHEP(I1)).EQ.5)IBQ=I1
          ENDDO
        ELSE
          CALL HWUEPR
          CALL HWWARN('HWANAL',505)
        ENDIF
C THE IDENTIFICATIONS ABOVE ARE CORRECT FOR TW- PRODUCTION; PARTICLES AND
C ANTIPARTICLES MAY BE EXCHANGED IN THE CASE OF TBARW+ PRODUCTION
        IF(ISGNT.LT.0)THEN
          IF(IDHEP(ILPT).LT.0)THEN
            ITMP=ILPT
            ILPT=INUT
            INUT=ITMP
          ENDIF
          IF(IDHEP(ILPW).GT.0)THEN
            ITMP=ILPW
            ILPW=INUW
            INUW=ITMP
          ENDIF
        ENDIF
C CHECK THAT THE DECAY PRODUCTS ARE WHAT THEY ARE SUPPOSED TO BE
        IF( ( (IDHEP(ILPT).NE.-11*ISGNT .AND.
     #         IDHEP(ILPT).NE.-13*ISGNT .AND.
     #         IDHEP(ILPT).NE.-15*ISGNT) .OR.
     #        (ISTHEP(ILPT).NE.1.AND.ISTHEP(ILPT).NE.195) ) .OR.
     #      ( (IDHEP(INUT).NE.12*ISGNT .AND.
     #         IDHEP(INUT).NE.14*ISGNT .AND.
     #         IDHEP(INUT).NE.16*ISGNT).OR.
     #        (ISTHEP(INUT).NE.1.AND.ISTHEP(INUT).NE.195) ) .OR.
     #      ( (IDHEP(ILPW).NE.-11*ISGNW .AND.
     #         IDHEP(ILPW).NE.-13*ISGNW .AND.
     #         IDHEP(ILPW).NE.-15*ISGNW).OR.
     #        (ISTHEP(ILPW).NE.1.AND.ISTHEP(ILPW).NE.195) ) .OR.
     #      ( (IDHEP(INUW).NE.12*ISGNW .AND.
     #         IDHEP(INUW).NE.14*ISGNW .AND.
     #         IDHEP(INUW).NE.16*ISGNW).OR.
     #        (ISTHEP(INUW).NE.1.AND.ISTHEP(INUW).NE.195) ) .OR.
     #        (IDHEP(IBQ).NE.5*ISGNT.OR.ISTHEP(IBQ).NE.2) )THEN
          CALL HWUEPR
          WRITE(*,*)'IT, IPT, IPTN',IT,IPT,IPTN
          DO I=1,IT
            WRITE(*,*)'JT',JT(I)
          ENDDO
          WRITE(*,*)ILPT,INUT,ILPW,INUW,IBQ
          CALL HWWARN('HWANAL',514)
        ENDIF
C FILL THE FOUR-MOMENTA OF THE DECAY PRODUCTS
        DO I=1,4
          PLPT(I)=PHEP(I,ILPT)
          PNUT(I)=PHEP(I,INUT)
          PBQT(I)=PHEP(I,IBQ)
          PLPW(I)=PHEP(I,ILPW)
          PNUW(I)=PHEP(I,INUW)
          YPPLL(I)=PLPT(I)+PLPW(I)
        ENDDO
      ELSE
C CONSTRUCT THE Ht PAIR MOMENTUM; PRIMARY TOP AND W ARE THE FIRST ON RECORD
        DO I=1,4
          YPPTPW(I)=PPTP(I,1)+PPW(I,1)
        ENDDO
      ENDIF
C ORDER STABLE B HADRONS IN PT; KEEP ONLY THOSE WITHIN |ETA|<RCUT==2.5 FOR LHC.
C USE THE SECOND HARDEST (IF PRESENT) TO APPLY THE VETO
      PTB1MAX=0.d0
      PTB2MAX=0.d0
      JBMAX(1)=0
      JBMAX(2)=0
      DO I=1,IB
        PTB=PTCALC(PB(1,I))
        ETAB=GETPSEUDORAP(PB(4,I),PB(1,I),PB(2,I),PB(3,I))
        IF(PTB.GT.PTB1MAX.AND.ABS(ETAB).LE.RCUT)THEN
          JBMAX(2)=JBMAX(1)
          PTB2MAX=PTB1MAX
          JBMAX(1)=I
          PTB1MAX=PTB
        ELSEIF(PTB.GT.PTB2MAX.AND.ABS(ETAB).LE.RCUT)THEN
          JBMAX(2)=I
          PTB2MAX=PTB
        ENDIF
      ENDDO
C If VETOCUT=.TRUE., accept the event
      VETOCUT=.FALSE.
      IF( JBMAX(2).EQ.0 .OR.
     #    (JBMAX(2).NE.0.AND.PTB2MAX.LE.BVETO) )VETOCUT=.TRUE.
C FILL THE HISTOS
      if(idec.eq.0)then
        xptlt=ptcalc(plpt(1))
        xylt=getrapidity(plpt(4),plpt(3))
        xptnut=ptcalc(pnut(1))
        xynut=getrapidity(pnut(4),pnut(3))
        xptbt=ptcalc(pbqt(1))
        xybt=getrapidity(pbqt(4),pbqt(3))
        xptlw=ptcalc(plpw(1))
        xylw=getrapidity(plpw(4),plpw(3))
        xptnuw=ptcalc(pnuw(1))
        xynuw=getrapidity(pnuw(4),pnuw(3))
        xptll=ptcalc(yppll)
        dphill=getdelphi(plpt(1),plpt(2),
     #                   plpw(1),plpw(2))
        xmll=getinvm(yppll(4),yppll(1),yppll(2),yppll(3))
        xlldelr=getdr(plpt(4),plpt(1),plpt(2),plpt(3),
     #                plpw(4),plpw(1),plpw(2),plpw(3))
c no veto
        kk=0
        call mfill(kk+ 1,sngl(xptlt),sngl(WWW0))
        call mfill(kk+ 2,sngl(xptnut),sngl(WWW0))
        call mfill(kk+ 3,sngl(xptbt),sngl(WWW0))
        call mfill(kk+ 4,sngl(xptlw),sngl(WWW0))
        call mfill(kk+ 5,sngl(xptnuw),sngl(WWW0))
        call mfill(kk+ 6,sngl(xylt),sngl(WWW0))
        call mfill(kk+ 7,sngl(xynut),sngl(WWW0))
        call mfill(kk+ 8,sngl(xybt),sngl(WWW0))
        call mfill(kk+ 9,sngl(xylw),sngl(WWW0))
        call mfill(kk+10,sngl(xynuw),sngl(WWW0))
        call mfill(kk+11,sngl(xptll),sngl(WWW0))
        call mfill(kk+12,sngl(dphill),sngl(WWW0))
        call mfill(kk+13,sngl(xmll),sngl(WWW0))
        call mfill(kk+14,sngl(xlldelr),sngl(WWW0))
c with veto
        kk=20
        if(vetocut)then
          call mfill(kk+ 1,sngl(xptlt),sngl(WWW0))
          call mfill(kk+ 2,sngl(xptnut),sngl(WWW0))
          call mfill(kk+ 3,sngl(xptbt),sngl(WWW0))
          call mfill(kk+ 4,sngl(xptlw),sngl(WWW0))
          call mfill(kk+ 5,sngl(xptnuw),sngl(WWW0))
          call mfill(kk+ 6,sngl(xylt),sngl(WWW0))
          call mfill(kk+ 7,sngl(xynut),sngl(WWW0))
          call mfill(kk+ 8,sngl(xybt),sngl(WWW0))
          call mfill(kk+ 9,sngl(xylw),sngl(WWW0))
          call mfill(kk+10,sngl(xynuw),sngl(WWW0))
          call mfill(kk+11,sngl(xptll),sngl(WWW0))
          call mfill(kk+12,sngl(dphill),sngl(WWW0))
          call mfill(kk+13,sngl(xmll),sngl(WWW0))
          call mfill(kk+14,sngl(xlldelr),sngl(WWW0))
        endif
      elseif(idec.eq.1)then
        xptt=ptcalc(pptp(1,1))
        xyt=getrapidity(pptp(4,1),pptp(3,1))
        xptw=ptcalc(ppw(1,1))
        xyw=getrapidity(ppw(4,1),ppw(3,1))
        xptwt=ptcalc(ypptpw)
        wtm=getinvm(ypptpw(4),ypptpw(1),ypptpw(2),ypptpw(3))
        azi=getdelphi(pptp(1,1),pptp(2,1),ppw(1,1),ppw(2,1))
        azinorm=(pi-azi)/pi
        wtdelr=getdr(pptp(4,1),pptp(1,1),pptp(2,1),pptp(3,1),
     #               ppw(4,1),ppw(1,1),ppw(2,1),ppw(3,1))
c no veto
        kk=0
        call mfill(kk+ 1,sngl(xptwt),sngl(WWW0))
        if(xptwt.gt.0.d0)
     #    call mfill(kk+ 2,sngl(log10(xptwt)),sngl(WWW0))
        call mfill(kk+ 3,sngl(wtm),sngl(WWW0))
        call mfill(kk+ 4,sngl(azi),sngl(WWW0))
        if(azinorm.gt.0.d0) 
     #    call mfill(kk+ 5,sngl(log10(azinorm)),sngl(WWW0))
        call mfill(kk+ 6,sngl(wtdelr),sngl(WWW0))
        call mfill(kk+ 7,sngl(xptt),sngl(WWW0))
        call mfill(kk+ 8,sngl(xptw),sngl(WWW0))
        call mfill(kk+ 9,sngl(xyt),sngl(WWW0))
        call mfill(kk+10,sngl(xyw),sngl(WWW0))
c with veto
        kk=20
        if(vetocut)then
          call mfill(kk+ 1,sngl(xptwt),sngl(WWW0))
          if(xptwt.gt.0.d0)
     #      call mfill(kk+ 2,sngl(log10(xptwt)),sngl(WWW0))
          call mfill(kk+ 3,sngl(wtm),sngl(WWW0))
          call mfill(kk+ 4,sngl(azi),sngl(WWW0))
          if(azinorm.gt.0.d0) 
     #      call mfill(kk+ 5,sngl(log10(azinorm)),sngl(WWW0))
          call mfill(kk+ 6,sngl(wtdelr),sngl(WWW0))
          call mfill(kk+ 7,sngl(xptt),sngl(WWW0))
          call mfill(kk+ 8,sngl(xptw),sngl(WWW0))
          call mfill(kk+ 9,sngl(xyt),sngl(WWW0))
          call mfill(kk+10,sngl(xyw),sngl(WWW0))
        endif
      else
        call hwwarn('HWANAL',555)
      endif
c
      RETURN
 999  END


      FUNCTION PTCALC(P)
      IMPLICIT NONE
      DOUBLE PRECISION PTCALC,P(4),PTSQ
      PTSQ=P(1)**2+P(2)**2
      IF (PTSQ.EQ.0D0) THEN
         PTCALC=0D0
      ELSE
         PTCALC=SQRT(PTSQ)
      ENDIF
      END


      FUNCTION ETCALC(P)
      IMPLICIT NONE
      DOUBLE PRECISION ETCALC,P(4),PTSQ
      PTSQ=P(1)**2+P(2)**2
      IF (PTSQ.EQ.0D0) THEN
         ETCALC=0D0
      ELSE
         ETCALC=P(4)*SQRT(PTSQ/(PTSQ+P(3)**2))
      ENDIF
      END


      FUNCTION YCALC(P)
      IMPLICIT NONE
      DOUBLE PRECISION YCALC,P(4),XXLOG
      IF (ABS(P(4)-ABS(P(3))).LT.1.D-8) THEN
         YCALC=SIGN(1.D8,P(3))
      ELSE
         YCALC=0.5D0*XXLOG((P(4)+P(3))/(P(4)-P(3)))
      ENDIF
      END


      FUNCTION XXLOG(X)
      IMPLICIT NONE
      DOUBLE PRECISION XXLOG,X,TMP
C
      IF(X.LT.1.D-8)THEN
        TMP=-1.D5
      ELSE
        TMP=LOG(X)
      ENDIF
      XXLOG=TMP
      RETURN
      END


      FUNCTION IHADR(ID)
c Returns the PDG code of the heavier quark in the hadron of PDG code ID
      IMPLICIT NONE
      INTEGER IHADR,ID,ID1
C
      IF(ID.NE.0)THEN
        ID1=ABS(ID)
        IF(ID1.GT.10000)ID1=ID1-1000*INT(ID1/1000)
        IHADR=ID1/(10**INT(LOG10(DFLOAT(ID1))))
      ELSE
        IHADR=0
      ENDIF
      RETURN
      END


      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-8)
c
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
        if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny )then
          y=0.5d0*log( xplus/xminus )
        else
          y=sign(1.d0,pl)*1.d8
        endif
      else
        y=sign(1.d0,pl)*1.d8
      endif
      getrapidity=y
      return
      end


      function getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-8)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end


      function getinvm(en,ptx,pty,pl)
      implicit none
      real*8 getinvm,en,ptx,pty,pl,tiny,tmp
      parameter (tiny=1.d-5)
c
      tmp=en**2-ptx**2-pty**2-pl**2
      if(tmp.gt.0.d0)then
        tmp=sqrt(tmp)
      elseif(tmp.gt.-tiny)then
        tmp=0.d0
      else
        write(*,*)'Attempt to compute a negative mass'
        stop
      endif
      getinvm=tmp
      return
      end


      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end


      function getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
      implicit none
      real*8 getdr,en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2,deta,dphi,
     # getpseudorap,getdelphi
c
      deta=getpseudorap(en1,ptx1,pty1,pl1)-
     #     getpseudorap(en2,ptx2,pty2,pl2)
      dphi=getdelphi(ptx1,pty1,ptx2,pty2)
      getdr=sqrt(dphi**2+deta**2)
      return
      end
