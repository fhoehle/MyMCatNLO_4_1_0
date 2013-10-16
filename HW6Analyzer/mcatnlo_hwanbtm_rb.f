C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      parameter (pi=3.14160E0)
      character*8 cc(4)
      data cc/' b      ',' b, cuts',' B      ',' B, cuts'/
c
      call rinit("HERQQ.root")
      do j=1,4
      k=(j-1)*15
      call rbook(k+ 1,'QQ pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 2,'QQ log[pt]'//cc(j),0.05e0,0.0e0,4.e0)
      call rbook(k+ 3,'QQ azimt'//cc(j),pi/20.e0,0.e0,pi)
      call rbook(k+ 4,'QQ log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
      call rbook(k+ 5,'QQ del R'//cc(j),pi/20.e0,0.e0,3*pi)
      call rbook(k+ 6,'Qb pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 7,'Q pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 8,'y_Qb'//cc(j),0.1e0,-4.e0,4.e0)
      call rbook(k+ 9,'y_Q'//cc(j),0.1e0,-4.e0,4.e0)
      call rbook(k+10,'# Q+Qb'//cc(j),1.e0,-0.5e0,40.5e0)
      call rbook(k+11,'# QQ pairs'//cc(j),1.e0,-0.5e0,40.5e0)
      enddo
      END



C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
c
      CALL RWRIT()
c
      END


C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWVDOT,PSUM(4)
      INTEGER ICHSUM,ICHINI,IHEP
      LOGICAL DIDSOF
      LOGICAL QQFLAG,SIQ1FLAG,SIQ2FLAG
      LOGICAL HHFLAG,SIH1FLAG,SIH2FLAG
      LOGICAL BHADRN,BMESON,BBARYON
      INTEGER ID,IST,IJ,IJ1,IJ2,IP
      INTEGER IQ,IQ1,IQ2,ICQ1,ICQ2,IQQPAIR,ICQQPAIR
      DOUBLE PRECISION PTQ1,PTQ2,ENQ1,ENQ2,PLQ1,PLQ2,YQ1,YQ2,PTQQ,
     # THQ1,THQ2,ETAQ1,ETAQ2,DQQ,CQQ,AZIQQ,AZQQNORM,DRQQ
      INTEGER IH,ICH,IHHPAIR,ICHHPAIR
      DOUBLE PRECISION PTH1,PTH2,ENH1,ENH2,PLH1,PLH2,YH1,YH2,PTHH,
     # THH1,THH2,ETAH1,ETAH2,DHH,CHH,AZIHH,AZHHNORM,DRHH
      REAL*8 PI
      PARAMETER (PI=3.14159265358979312D0)
      DOUBLE PRECISION TINY
      DOUBLE PRECISION PQ1(5,100),PQ2(5,100),PQ(5,100)
      DOUBLE PRECISION PH(5,100)
      REAL*8 WWW0
      DATA TINY/.1D-5/
c
      IF (IERROR.NE.0) RETURN
c
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
      IQ1=0
      IQ2=0
      IQ=0
      IH=0
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
        IST=ISTHEP(IHEP)      
        ID=IDHEP(IHEP)
        CALL BHAD(ID,BHADRN,BMESON,BBARYON)
        IF(ABS(IPROC).EQ.1705.OR.ABS(IPROC).EQ.11705)THEN
C LOOK FOR Q AND Qb; KEEP ALL OF THEM. THESE ARE THE QUARKS AFTER THE
C SHOWER -- THEY SHOULD BE FULLY PERTURBATIVE. IF THE GLUON SPLITTING-PHASE
C NEEDS TO BE INCLUDED, USE IST.EQ.158
          IF (ABS(ID).EQ.5.AND.IST.EQ.2) THEN
            IF(ID.EQ.5)THEN
C FOUND A Q
              IQ1=IQ1+1
              IQ=IQ+1
              DO IJ=1,5
                PQ1(IJ,IQ1)=PHEP(IJ,IHEP)
                PQ(IJ,IQ)=PHEP(IJ,IHEP)
              ENDDO
            ELSEIF(ID.EQ.-5)THEN
C FOUND A Qb
              IQ2=IQ2+1
              IQ=IQ+1
              DO IJ=1,5
                PQ2(IJ,IQ2)=PHEP(IJ,IHEP)
                PQ(IJ,IQ)=PHEP(IJ,IHEP)
              ENDDO
            ENDIF
          ELSEIF (BHADRN.AND.IST.GT.194.AND.IST.LT.200
     &           .AND.RLTIM(IDHW(IHEP)).GT.1D-13) THEN
C FOUND A WEAKLY-DECAYING Q-MESON
            IH=IH+1
            DO IJ=1,5
              PH(IJ,IH)=PHEP(IJ,IHEP)
            ENDDO
          ENDIF
        ELSE
          CALL HWWARN('HWANAL',500)
        ENDIF
  100 CONTINUE
      IF(IQ.NE.(IQ1+IQ2))CALL HWWARN('HWANAL',501)
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
C DEAL WITH THE QUARKS
      IP=0
      ICQ1=0
      ICQ2=0
C LOOP OVER Q'S
      DO IJ=1,IQ1
        PTQ1 = DSQRT(PQ1(1,IJ)**2+PQ1(2,IJ)**2)
        PLQ1 = PQ1(3,IJ)
        ENQ1 = PQ1(4,IJ)
        YQ1 = 0.5D0*LOG((ENQ1+PLQ1)/(ENQ1-PLQ1))
        CALL RFILL(IP+7,SNGL(PTQ1),SNGL(WWW0))
        CALL RFILL(IP+9,SNGL(YQ1),SNGL(WWW0))
        IF(ABS(YQ1).LT.1.D0)CALL RFILL(IP+15+7,SNGL(PTQ1),SNGL(WWW0))
        IF(PTQ1.GT.5.D0)CALL RFILL(IP+15+9,SNGL(YQ1),SNGL(WWW0))
        IF(ABS(YQ1).LT.1.D0.AND.PTQ1.GT.5.D0)ICQ1=ICQ1+1
      ENDDO
C LOOP OVER Qb'S
      DO IJ=1,IQ2
        PTQ2 = DSQRT(PQ2(1,IJ)**2+PQ2(2,IJ)**2)
        PLQ2 = PQ2(3,IJ)
        ENQ2 = PQ2(4,IJ)
        YQ2 = 0.5D0*LOG((ENQ2+PLQ2)/(ENQ2-PLQ2))
        CALL RFILL(IP+6,SNGL(PTQ2),SNGL(WWW0))
        CALL RFILL(IP+8,SNGL(YQ2),SNGL(WWW0))
        IF(ABS(YQ2).LT.1.D0)CALL RFILL(IP+15+6,SNGL(PTQ2),SNGL(WWW0))
        IF(PTQ2.GT.5.D0)CALL RFILL(IP+15+8,SNGL(YQ2),SNGL(WWW0))
        IF(ABS(YQ2).LT.1.D0.AND.PTQ2.GT.5.D0)ICQ2=ICQ2+1
      ENDDO
      CALL RFILL(IP+10,REAL(IQ1+IQ2),SNGL(WWW0))
      CALL RFILL(IP+15+10,REAL(ICQ1+ICQ2),SNGL(WWW0))
C LOOP OVER Q-Q, Q-Qb, Qb-Qb PAIRS
      IQQPAIR=0
      ICQQPAIR=0
      DO IJ1=1,IQ
        DO IJ2=IJ1+1,IQ
          IQQPAIR=IQQPAIR+1
          PTQ1 = DSQRT(PQ(1,IJ1)**2+PQ(2,IJ1)**2)
          PLQ1 = PQ(3,IJ1)
          ENQ1 = PQ(4,IJ1)
          PTQ2 = DSQRT(PQ(1,IJ2)**2+PQ(2,IJ2)**2)
          PLQ2 = PQ(3,IJ2)
          ENQ2 = PQ(4,IJ2)
          YQ1 = 0.5D0*LOG((ENQ1+PLQ1)/(ENQ1-PLQ1))
          YQ2 = 0.5D0*LOG((ENQ2+PLQ2)/(ENQ2-PLQ2))
          PTQQ = DSQRT( (PQ(1,IJ1)+PQ(1,IJ2))**2+
     #                  (PQ(2,IJ1)+PQ(2,IJ2))**2 )
          THQ1 = ATAN2(PTQ1+TINY,PLQ1)
          THQ2 = ATAN2(PTQ2+TINY,PLQ2)
          ETAQ1= -LOG(TAN(THQ1/2))
          ETAQ2= -LOG(TAN(THQ2/2))
          DQQ = PQ(1,IJ1)*PQ(1,IJ2)+PQ(2,IJ1)*PQ(2,IJ2)
          CQQ=0
          IF(PTQ1.NE.0.AND.PTQ2.NE.0) CQQ = DQQ * (1-TINY)/(PTQ1*PTQ2)
          IF(ABS(CQQ).GT.1) THEN
	    WRITE(*,*) ' COSINE = ',CQQ ,DQQ,PTQ1,PTQ2
            CQQ = - 1
          ENDIF
          AZIQQ = (1-TINY)*ACOS(CQQ)
          AZQQNORM = (PI-AZIQQ)/PI
          DRQQ  = DSQRT(AZIQQ**2+(ETAQ1-ETAQ2)**2)
          SIQ1FLAG=PTQ1.GT.5.D0.AND.ABS(YQ1).LT.1.D0
          SIQ2FLAG=PTQ2.GT.5.D0.AND.ABS(YQ2).LT.1.D0
          QQFLAG=SIQ1FLAG.AND.SIQ2FLAG
          CALL RFILL(IP+1,SNGL(PTQQ),SNGL(WWW0))
          IF(PTQQ.GT.0) CALL RFILL(IP+2,SNGL(LOG10(PTQQ)),SNGL(WWW0))
          CALL RFILL(IP+3,SNGL(AZIQQ),SNGL(WWW0))
          IF(AZQQNORM.GT.0) 
     #      CALL RFILL(IP+4,SNGL(LOG10(AZQQNORM)),SNGL(WWW0))
          CALL RFILL(IP+5,SNGL(DRQQ),SNGL(WWW0))
          IF(QQFLAG)THEN
            ICQQPAIR=ICQQPAIR+1
            CALL RFILL(IP+15+1,SNGL(PTQQ),SNGL(WWW0))
            IF(PTQQ.GT.0) 
     #        CALL RFILL(IP+15+2,SNGL(LOG10(PTQQ)),SNGL(WWW0))
            CALL RFILL(IP+15+3,SNGL(AZIQQ),SNGL(WWW0))
            IF(AZQQNORM.GT.0) 
     #        CALL RFILL(IP+15+4,SNGL(LOG10(AZQQNORM)),SNGL(WWW0))
            CALL RFILL(IP+15+5,SNGL(DRQQ),SNGL(WWW0))
          ENDIF
        ENDDO
      ENDDO
      CALL RFILL(IP+11,REAL(IQQPAIR),SNGL(WWW0))
      CALL RFILL(IP+15+11,REAL(ICQQPAIR),SNGL(WWW0))
C DEAL WITH THE HADRONS
      IP=30
      ICH=0
C LOOP OVER Q-HADRONS AND Qb-HADRONS
      DO IJ=1,IH
        PTH1 = DSQRT(PH(1,IJ)**2+PH(2,IJ)**2)
        PLH1 = PH(3,IJ)
        ENH1 = PH(4,IJ)
        YH1 = 0.5D0*LOG((ENH1+PLH1)/(ENH1-PLH1))
        CALL RFILL(IP+7,SNGL(PTH1),SNGL(WWW0))
        CALL RFILL(IP+9,SNGL(YH1),SNGL(WWW0))
        IF(ABS(YH1).LT.1.D0)CALL RFILL(IP+15+7,SNGL(PTH1),SNGL(WWW0))
        IF(PTH1.GT.5.D0)CALL RFILL(IP+15+9,SNGL(YH1),SNGL(WWW0))
        IF(ABS(YH1).LT.1.D0.AND.PTH1.GT.5.D0)ICH=ICH+1
      ENDDO
      CALL RFILL(IP+10,REAL(IH),SNGL(WWW0))
      CALL RFILL(IP+15+10,REAL(ICH),SNGL(WWW0))
C LOOP OVER Q-Q, Q-Qb, Qb-Qb HADRON PAIRS
      IHHPAIR=0
      ICHHPAIR=0
      DO IJ1=1,IH
        DO IJ2=IJ1+1,IH
          IHHPAIR=IHHPAIR+1
          PTH1 = DSQRT(PH(1,IJ1)**2+PH(2,IJ1)**2)
          PLH1 = PH(3,IJ1)
          ENH1 = PH(4,IJ1)
          PTH2 = DSQRT(PH(1,IJ2)**2+PH(2,IJ2)**2)
          PLH2 = PH(3,IJ2)
          ENH2 = PH(4,IJ2)
          YH1 = 0.5D0*LOG((ENH1+PLH1)/(ENH1-PLH1))
          YH2 = 0.5D0*LOG((ENH2+PLH2)/(ENH2-PLH2))
          PTHH = DSQRT( (PH(1,IJ1)+PH(1,IJ2))**2+
     #                  (PH(2,IJ1)+PH(2,IJ2))**2 )
          THH1 = ATAN2(PTH1+TINY,PLH1)
          THH2 = ATAN2(PTH2+TINY,PLH2)
          ETAH1= -LOG(TAN(THH1/2))
          ETAH2= -LOG(TAN(THH2/2))
          DHH = PH(1,IJ1)*PH(1,IJ2)+PH(2,IJ1)*PH(2,IJ2)
          CHH=0
          IF(PTH1.NE.0.AND.PTH2.NE.0) CHH = DHH * (1-TINY)/(PTH1*PTH2)
          IF(ABS(CHH).GT.1) THEN
	    WRITE(*,*) ' COSINE = ',CHH ,DHH,PTH1,PTH2
            CHH = - 1
          ENDIF
          AZIHH = (1-TINY)*ACOS(CHH)
          AZHHNORM = (PI-AZIHH)/PI
          DRHH  = DSQRT(AZIHH**2+(ETAH1-ETAH2)**2)
          SIH1FLAG=PTH1.GT.5.D0.AND.ABS(YH1).LT.1.D0
          SIH2FLAG=PTH2.GT.5.D0.AND.ABS(YH2).LT.1.D0
          HHFLAG=SIH1FLAG.AND.SIH2FLAG
          CALL RFILL(IP+1,SNGL(PTHH),SNGL(WWW0))
          IF(PTHH.GT.0) CALL RFILL(IP+2,SNGL(LOG10(PTHH)),SNGL(WWW0))
          CALL RFILL(IP+3,SNGL(AZIHH),SNGL(WWW0))
          IF(AZHHNORM.GT.0) 
     #      CALL RFILL(IP+4,SNGL(LOG10(AZHHNORM)),SNGL(WWW0))
          CALL RFILL(IP+5,SNGL(DRHH),SNGL(WWW0))
          IF(HHFLAG)THEN
            ICHHPAIR=ICHHPAIR+1
            CALL RFILL(IP+15+1,SNGL(PTHH),SNGL(WWW0))
            IF(PTHH.GT.0) 
     #        CALL RFILL(IP+15+2,SNGL(LOG10(PTHH)),SNGL(WWW0))
            CALL RFILL(IP+15+3,SNGL(AZIHH),SNGL(WWW0))
            IF(AZHHNORM.GT.0) 
     #        CALL RFILL(IP+15+4,SNGL(LOG10(AZHHNORM)),SNGL(WWW0))
            CALL RFILL(IP+15+5,SNGL(DRHH),SNGL(WWW0))
          ENDIF
        ENDDO
      ENDDO
      CALL RFILL(IP+11,REAL(IHHPAIR),SNGL(WWW0))
      CALL RFILL(IP+15+11,REAL(ICHHPAIR),SNGL(WWW0))
 999  END

C----------------------------------------------------------------------
      SUBROUTINE BHAD(IDPDG,BHADRN,BMESON,BBARYON)
C     TEST FOR A B-FLAVOURED HADRON
C----------------------------------------------------------------------
      INTEGER IDPDG,ID,IM,IB
      LOGICAL BHADRN,BMESON,BBARYON
C
      ID=ABS(IDPDG)
      IM=MOD(ID/100,100)
      IB=MOD(ID/1000,100)
      BMESON=IM.EQ.5
      BBARYON=IB.EQ.5
      IF(BMESON.AND.BBARYON)CALL HWWARN('BHAD  ',500)
      BHADRN=BMESON.OR.BBARYON
 999  END

