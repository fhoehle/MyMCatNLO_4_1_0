C Collects all of the Les Houches interface routines, plus utilities
C for colour codes
C Should be backward compatible with v3.1 and v3.2 (thanks to B. Kersevan)
C
C----------------------------------------------------------------------
      SUBROUTINE UPEVNT
C----------------------------------------------------------------------
C  Reads MC@NLO input files and fills Les Houches event common HEPEUP
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
C---Les Houches Event Common Block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP,
     & XMP2,XMA2,XMB2,BETA,VA,VB,SIGMA,DELTA,S2,XKA,XKB,PTF,E,PL,
     & XSCALE,XEPHO
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
      DOUBLE PRECISION UX1,UX2,UQ2
      COMMON/CPDFRWGT/UX1,UX2,UQ2
      DOUBLE PRECISION WGTACP(28)
      COMMON/CWGTACP/WGTACP
      DOUBLE PRECISION PCM(5),PTR,XMTR,HWVDOT,HWULDO,PDB(5)
      INTEGER I,J,IC,JPR,MQQ,NQQ,IUNIT,ISCALE,I1HPRO,IBOS,NP,IG,IP,
     & ILEP,ID,IA,IB,ICOL4(4,5),ICOL5(5,18),ICOP5(5,4),ICOST(4,20),
     & IVBF4(4,4),IVBF5(5,8,2),JJPROC,IVHVEC,IVHLEP,MUP,JC
      PARAMETER (IUNIT=61)
      LOGICAL BOPRO,NODEC,REMIT
      COMMON/NQQCOM/MQQ,NQQ
      COMMON/VHLIN/IVHVEC,IVHLEP
      COMMON/PHOCOM/XEPHO
C---Colour flows for heavy quark pair production
      DATA ICOL4/
     & 10,02,10,02,01,20,20,01,12,23,10,03,12,31,30,02,00,12,10,02/
      DATA ICOL5/
     & 10,02,13,30,02, 10,02,32,10,03,
     & 10,21,30,20,03, 10,23,20,10,03,
     & 01,20,23,30,01, 01,20,31,20,03,
     & 01,23,03,20,01, 01,12,03,30,02,
     & 12,20,30,10,03, 12,30,10,30,02,
     & 12,03,02,10,03, 12,01,03,30,02,
     & 12,23,14,40,03, 12,34,32,10,04,
     & 12,23,43,10,04, 12,31,34,40,02,
     & 12,34,14,30,02, 12,31,42,30,04/
C---Colour flows for heavy quark pair photoproduction
      DATA ICOP5/
     & 00,12,32,10,03, 00,12,13,30,02,
     & 00,10,20,10,02, 00,01,02,20,01/
C---Colour flows for single top production
      DATA ICOST/
     & 20,12,00,10, 12,20,00,10, 10,02,02,10, 02,10,02,10,
     & 20,10,20,10, 20,02,01,10, 10,20,20,10, 02,20,01,10,
     & 20,02,01,10, 10,02,02,10, 02,20,01,10, 02,10,02,10,
     & 20,10,20,10, 10,20,20,10, 12,23,03,10, 23,12,03,10,
     & 30,12,32,10, 30,23,21,10, 12,30,32,10, 23,30,21,10/
C---Colour flows for vector boson fusion
      DATA IVBF4/
     & 10,20,10,20, 10,02,10,02, 01,20,01,20, 01,02,01,02/
      DATA IVBF5/
     & 10,20,13,30,20, 10,02,13,30,02, 01,20,31,03,20, 01,02,31,03,02,
     & 12,30,10,02,30, 12,03,10,02,03, 10,23,20,10,03, 01,23,20,01,03,
     & 10,20,23,10,30, 10,02,32,10,03, 01,20,23,01,30, 01,02,32,01,03,
     & 12,30,02,10,30, 12,03,02,10,03, 10,23,03,10,20, 01,23,03,01,20/

      IF (IERROR.NE.0) RETURN
C---CHECK PROCESS CODE
      JJPROC=MOD(ABS(IPROC),10000)
      JPR=JJPROC/100
C---READ AN EVENT
      IF(NQQ.GE.MQQ)CALL HWWARN('UPEVNT',201)
      READ(IUNIT,901) I1HPRO,IC,NP
      READ(IUNIT,902) (IDUP(I),I=1,NP)
      IF (JPR.EQ.51) THEN
         READ(IUNIT,903) XWGTUP,XSCALE,XEPHO
      ELSE
         READ(IUNIT,903) XWGTUP,XSCALE
      ENDIF
C---Les Houches expects mean weight to be the cross section in pb
      XWGTUP= XWGTUP*MQQ
      READ(IUNIT,904) ((PUP(J,I),J=1,4),I=1,NP)
      READ(IUNIT,905) UX1,UX2,UQ2
C---The following are relevant only to WZ and WW production in v4.0
      IF(JJPROC.EQ.2850)THEN
        READ(IUNIT,907) (WGTACP(J),J=1,28)
      ELSEIF(JJPROC.EQ.2870.OR.JJPROC.EQ.2880)THEN
        READ(IUNIT,906) (WGTACP(J),J=1,10)
      ENDIF
      NQQ=NQQ+1
C---Input format is now (px,py,pz,m)
      DO I=1,NP
         E=SQRT(HWVDOT(4,PUP(1,I),PUP(1,I)))
         PUP(5,I)=PUP(4,I)
         PUP(4,I)=E
      ENDDO
      CALL HWVSUM(4,PUP(1,1),PUP(1,2),PCM)
      CALL HWUMAS(PCM)
C---REMIT MEANS A REAL PARTON EMISSION OCCURRED
      REMIT=PUP(4,3).NE.ZERO
C---NODEC MEANS DECAYS NOT YET DONE
      NODEC=NP.EQ.5
      NUP=NP
      BOPRO=JPR.EQ.13.OR.JPR.EQ.14.OR.JPR.EQ.16.OR.JPR.EQ.36
      IF (BOPRO) THEN
C----------------------------------------------------------------------
C   SINGLE GAUGE OR HIGGS BOSON PRODUCTION
C   B = Z/gamma, W or H (SM or any MSSM neutral Higgs)
C-----------------------------------------------------------------------
C I1HPRO IDENTIFIES THE PARTONIC SUBPROCESS, WITH THE FOLLOWING CONVENTIONS:
C   I1HPRO         PROCESS
C    401        q qbar -> g B
C    402        q g    -> q B
C    403        qbar q -> g B
C    404        qbar g -> qbar B
C    405        g q    -> q B
C    406        g qbar -> qbar B
C    407        g g    -> g B
C-----------------------------------------------------------------------
C---NODEC=.TRUE. FOR HIGGS AND UNDECAYED EW BOSON
         NODEC=NP.EQ.4
         IHPRO=I1HPRO-400
         ISCALE=0
         IF(JPR.EQ.16)ISCALE=2
      ELSEIF (JPR.EQ.17.OR.JPR.EQ.20.OR.JPR.EQ.51) THEN
C----------------------------------------------------------------------
C   HEAVY Q and/or QBAR PRODUCTION
C   IPROC=-1705,-1706 for Q=b,t
C   IPROC=-5105,-5106 for Q=b,t photoproduction
C   IPROC=-2000,-2030/1/4 for single top, Wt+Wtb/Wtb/Wt
C-----------------------------------------------------------------------
C I1HPRO IDENTIFIES THE PARTONIC SUBPROCESS, WITH THE FOLLOWING CONVENTIONS:
C   I1HPRO         PROCESS
C    401        q qbar -> g Q Qbar
C    402        q g    -> q Q Qbar
C    403        qbar q -> g Q Qbar
C    404        qbar g -> qbar Q Qbar
C    405        g q    -> q Q Qbar
C    406        g qbar -> qbar Q Qbar
C    407        g g    -> g Q Qbar
C    408        q q    -> g t q
C    409        qbar qbar -> g tbar qbar
C-----------------------------------------------------------------------
C    505        gamma q    -> q Q Qbar
C    506        gamma qbar -> qbar Q Qbar
C    507        gamma g    -> g Q Qbar
C-----------------------------------------------------------------------
C IC SPECIFIES THE COLOUR CONNECTION (NOW IN INPUT FILE)
C-----------------------------------------------------------------------
C---SET IHPRO AS FOR DIRECT PHOTON (IPROC=1800)
         IHPRO=I1HPRO-360
         ISCALE=0
         IF(ABS(IPROC).EQ.1705.OR.ABS(IPROC).EQ.11705)ISCALE=5
      ELSEIF (JPR.EQ.19) THEN
C----------------------------------------------------------------------
C   Higgs production by vector boson fusion
C   IPROC=-1900-ID where ID controls Higgs decay
C-----------------------------------------------------------------------
C I1HPRO IDENTIFIES THE PARTONIC SUBPROCESS, WITH THE FOLLOWING CONVENTIONS:
C   I1HPRO         PROCESS
C    401        q q    -> g q q h
C    402        q qbar -> g q qbar h
C    403        qbar q -> g qbar q h
C    404        qbar qbar -> g qbar qbar h
C    405        g q    -> q qbar q h OR qbar q q h 
C    406        g qbar -> q qbar qbar h OR qbar q qbar h 
C    407        q g    -> q q qbar h OR qbar q q h 
C    408        qbar g -> q qbar qbar h OR qbar qbar q h 
C-----------------------------------------------------------------------
C IC SPECIFIES THE COLOUR CONNECTION (NOW IN INPUT FILE)
C For processes 401-4, IC = 1,2 for gluon (3) connection to parton IC
C For processes 405-8, IC = 1,2 for 'emitted' parton (3) = q,qbar
C-----------------------------------------------------------------------
         IHPRO=I1HPRO-400
C--Scale??
         ISCALE=0
      ELSEIF (JPR.EQ.28) THEN
C----------------------------------------------------------------------
C   GAUGE BOSON PAIR PRODUCTION
C   VV=WW,ZZ,ZW+,ZW- FOR IPROC=-2850,-2860,-2870,-2880
C-----------------------------------------------------------------------
C I1HPRO IDENTIFIES THE PARTONIC SUBPROCESS, WITH THE FOLLOWING CONVENTIONS:
C   I1HPRO         PROCESS
C    401        q qbar -> g V V
C    402        q g    -> q V V
C    403        qbar q -> g V V
C    404        qbar g -> qbar V V
C    405        g q    -> q V V
C    406        g qbar -> qbar V V
C-----------------------------------------------------------------------
         IHPRO=I1HPRO-400
         ISCALE=0
      ELSEIF (JPR.EQ.26.OR.JPR.EQ.27) THEN
C----------------------------------------------------------------------
C   GAUGE BOSON PLUS HIGGS PRODUCTION
C   VH=WH,ZH FOR IPROC=-2600-ID,-2700-ID
C   WHERE ID CONTROLS HIGGS DECAY AS IN STANDARD HERWIG
C-----------------------------------------------------------------------
         IHPRO=I1HPRO-400
         ISCALE=0
      ELSEIF (JPR.EQ.1) THEN
C----------------------------------------------------------------------
C   E+E- ANNIHILATION
C----------------------------------------------------------------------
         ISCALE=1
      ELSE
         CALL HWWARN('UPEVNT',202)
      ENDIF
C---HARD SCALE
      SCALUP=PCM(5)
      IF (XSCALE.GT.0D0.AND.XSCALE.LT.PCM(5)) SCALUP=XSCALE
      IF (REMIT) THEN
         IF (ISCALE.EQ.0) THEN
            PTR=SQRT(PUP(1,3)**2+PUP(2,3)**2)
            SCALUP=PCM(5)-2.*PTR
         ELSEIF(ISCALE.EQ.1)THEN
            SCALUP=PCM(5)
         ELSEIF (ISCALE.EQ.2) THEN
            SCALUP=SQRT(PUP(1,3)**2+PUP(2,3)**2)
         ELSEIF (ISCALE.EQ.3.OR.ISCALE.EQ.4.OR.ISCALE.EQ.5) THEN
            PTR=SQRT(PUP(1,3)**2+PUP(2,3)**2)
            IA=4
            IB=5
            XMP2=PUP(5,3)**2
            XMA2=PUP(5,IA)**2
            XMB2=PUP(5,IB)**2
            S2=XMA2+XMB2+2*HWULDO(PUP(1,IA),PUP(1,IB))
            SIGMA=XMA2+XMB2
            DELTA=XMA2-XMB2
            BETA=SQRT(1-2*SIGMA/S2+(DELTA/S2)**2)
            VA=BETA/(1+DELTA/S2)
            VB=BETA/(1-DELTA/S2)
            XKA=HWULDO(PUP(1,3),PUP(1,IA))
            XKB=HWULDO(PUP(1,3),PUP(1,IB))
            E=(XKA+XKB)/SQRT(S2)
            PL=-2.0/((VA+VB)*BETA*SQRT(S2))*(VA*XKA-VB*XKB)
            PTF=E**2-PL**2-XMP2
            IF (PTF.LE.ZERO) THEN
               CALL HWWARN('UPEVNT',103)
               GOTO 999
            ENDIF
            PTF=SQRT(PTF)
            IF(ISCALE.EQ.3)THEN
              SCALUP=PCM(5)-2.*MIN(PTR,PTF)
            ELSEIF(ISCALE.EQ.4)THEN
              SCALUP=MIN(PTR,PTF)
            ELSE
              SCALUP=(MIN(PTR,PTF))**2+(XMA2+XMB2)/2.D0
              SCALUP=SQRT(SCALUP)
            ENDIF
            IF (SCALUP.LE.ZERO) THEN
               CALL HWWARN('UPEVNT',100)
               GOTO 999
            ENDIF
         ELSEIF (ISCALE.EQ.6) THEN
            XMTR=SQRT(PUP(5,4)**2+PUP(1,4)**2+PUP(2,4)**2)
            PTR=SQRT(PUP(1,3)**2+PUP(2,3)**2)
            SCALUP=PCM(5)-PTR-XMTR
            IF (SCALUP.LE.ZERO) THEN
               CALL HWWARN('UPEVNT',100)
               GOTO 999
            ENDIF
         ELSEIF (ISCALE.EQ.7) THEN
            SCALUP=SQRT(PUP(5,4)**2+PUP(1,4)**2+PUP(2,4)**2)
         ELSE
            CALL HWWARN('UPEVNT',501)
         ENDIF
      ELSE
         NUP=NUP-1
      ENDIF
C---FINAL STATE (INITIAL STATE DONE BELOW)
      DO I=3,NUP
         ISTUP(I)=1
         MOTHUP(1,I)=1
         MOTHUP(2,I)=2
      ENDDO
      IF (BOPRO.AND.NODEC) THEN
C---SINGLE BOSON (UNDECAYED)
         IF (REMIT) THEN
C---SET COLOUR CONNECTIONS
            DO I=1,3
               ICOLUP(1,I)=501
               ICOLUP(2,I)=502
            ENDDO
            IF (IHPRO.EQ.1) THEN
               ICOLUP(2,1)=0
               ICOLUP(1,2)=0
            ELSEIF (IHPRO.EQ.2) THEN
               ICOLUP(1,1)=502
               ICOLUP(2,1)=0
               ICOLUP(2,3)=0
            ELSEIF (IHPRO.EQ.3) THEN
               ICOLUP(1,1)=0
               ICOLUP(2,2)=0
            ELSEIF (IHPRO.EQ.4) THEN
               ICOLUP(1,1)=0
               ICOLUP(2,1)=501
               ICOLUP(1,3)=0
            ELSEIF (IHPRO.EQ.5) THEN
               ICOLUP(1,2)=502
               ICOLUP(2,2)=0
               ICOLUP(2,3)=0
            ELSEIF (IHPRO.EQ.6) THEN
               ICOLUP(1,2)=0
               ICOLUP(2,2)=501
               ICOLUP(1,3)=0
            ELSEIF (IHPRO.EQ.7) THEN
               ICOLUP(1,2)=502
               ICOLUP(2,2)=503
               ICOLUP(2,3)=503
            ELSE
               CALL HWWARN('UPEVNT',101)
               GOTO 999
            ENDIF
         ELSE
            CALL HWVEQU(5,PUP(1,4),PUP(1,3))
C---SET COLOUR CONNECTIONS
            DO I=1,2
               ICOLUP(1,I)=0
               ICOLUP(2,I)=0
            ENDDO
            IF (IDUP(1).GT.0) THEN
               ICOLUP(1,1)=501
               ICOLUP(2,2)=501
               IF (IDUP(1).GT.0) THEN
C---GG FUSION
                  ICOLUP(2,1)=502
                  ICOLUP(1,2)=502
               ENDIF
            ELSE
C---QBAR Q
               ICOLUP(2,1)=501
               ICOLUP(1,2)=501
            ENDIF
         ENDIF
         ICOLUP(1,NUP)=0
         ICOLUP(2,NUP)=0
C---LOAD BOSON ID
         IF (JPR.EQ.13) THEN
            IDUP(NUP)=23
         ELSEIF (JPR.EQ.16) THEN
            IDUP(NUP)=25
         ELSEIF (JPR.EQ.36) THEN
            IBOS=MOD(JJPROC,100)
            IF (IBOS.EQ.10) THEN
               IDUP(NUP)=26
            ELSEIF (IBOS.EQ.20) THEN
               IDUP(NUP)=35
            ELSEIF (IBOS.EQ.30) THEN
               IDUP(NUP)=36 
            ELSE
               CALL HWWARN('UPEVNT',502)
            ENDIF
         ELSEIF (JPR.EQ.14) THEN
            IBOS=0
            DO I=1,NUP-1
               ID=IDUP(I)
               IF (ID.EQ.21) THEN
                  IC=0
               ELSEIF (ID.GT.0) THEN
                  IC=ICHRG(ID)
               ELSE
                  IC=ICHRG(6-ID)
               ENDIF
               IBOS=IBOS+IC
            ENDDO
            IF (REMIT) IBOS=IBOS-2*IC
            IF (ABS(IBOS).NE.3) CALL  HWWARN('UPEVNT',503)
            IDUP(NUP)=8*IBOS
         ENDIF
      ELSEIF (JPR.EQ.51.AND.IDUP(1).EQ.22) THEN
C---HEAVY QUARK DIRECT PHOTOPRODUCTION
         IF (REMIT) THEN
C---3-BODY FINAL STATE
C---SET COLOUR CONNECTIONS
            IF (IC.GE.4.AND.IC.LE.7) THEN
               JC=IC-3
               DO I=1,5
                  CALL UPCODE(ICOP5(I,JC),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',105)
               GOTO 999
            ENDIF
         ELSE
C---2-BODY FINAL STATE
            DO IP=3,NUP
               IDUP(IP)=IDUP(IP+1)
               CALL HWVEQU(5,PUP(1,IP+1),PUP(1,IP))
            ENDDO
C---SET COLOUR CONNECTIONS
            IF (IC.GE.1.AND.IC.LE.3) THEN
               JC=IC
               IF (IC.EQ.3) JC=5
               DO I=1,4
                  CALL UPCODE(ICOL4(I,JC),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',104)
               GOTO 999
            ENDIF
         ENDIF
      ELSEIF (JPR.EQ.19) THEN
C---HIGGS VIA VBF
         IF (REMIT) THEN
C---4-BODY FINAL STATE
C---SET COLOUR CONNECTIONS
            IF (IHPRO.GE.5) THEN
               IC=1
               IF (IDUP(3).LT.0) IC=2
            ENDIF
            IF (IHPRO.LE.8.AND.IC.LE.2) THEN
               DO I=1,5
                  CALL UPCODE(IVBF5(I,IHPRO,IC),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',105)
               GOTO 999
            ENDIF
         ELSE
C---3-BODY FINAL STATE
            DO IP=3,NUP
               IDUP(IP)=IDUP(IP+1)
               CALL HWVEQU(5,PUP(1,IP+1),PUP(1,IP))
            ENDDO
C---SET COLOUR CONNECTIONS
            IF (IHPRO.LE.4) THEN
               DO I=1,4
                  CALL UPCODE(IVBF4(I,IHPRO),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',104)
               GOTO 999
            ENDIF
         ENDIF
         ICOLUP(1,NUP)=0
         ICOLUP(2,NUP)=0
      ELSEIF (JPR.EQ.17.OR.(JPR.EQ.51.AND.IDUP(1).NE.22)) THEN
C---HEAVY QUARKS
         IF (REMIT) THEN
C---3-BODY FINAL STATE
C---SET COLOUR CONNECTIONS
            IF (IC.LE.18) THEN
               DO I=1,5
                  CALL UPCODE(ICOL5(I,IC),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',105)
               GOTO 999
            ENDIF
         ELSE
C---2-BODY FINAL STATE
            DO IP=3,NUP
               IDUP(IP)=IDUP(IP+1)
               CALL HWVEQU(5,PUP(1,IP+1),PUP(1,IP))
            ENDDO
C---SET COLOUR CONNECTIONS
            IF (IC.LE.4) THEN
               DO I=1,4
                  CALL UPCODE(ICOL4(I,IC),ICOLUP(1,I))
               ENDDO
            ELSE
               CALL HWWARN('UPEVNT',104)
               GOTO 999
            ENDIF
         ENDIF
         IF (.NOT.NODEC) THEN
C---T-TBAR WITH DECAYS
C   LES HOUCHES COMMON IS FILLED AS FOLLOWS:
C
C  1 INCOMING PARTON
C  2 INCOMING PARTON
C  3 OUTGOING LIGHT PARTON (IF ANY, OTHERWISE 4-13 BECOME 3-12)
C  4 TOP
C  5 ANTITOP
C  6 QUARK FROM TOP
C  7 W+ FROM TOP
C  8 ANTIQUARK FROM ANTITOP
C  9 W- FROM ANTITOP
C 10 DECAY PRODUCT FROM W+ FROM TOP
C 11 DECAY PRODUCT FROM W+ FROM TOP
C 12 DECAY PRODUCT FROM W- FROM ANTITOP
C 13 DECAY PRODUCT FROM W- FROM ANTITOP
C
C---RECONSTRUCT TOP DECAYS
            IF (MOD(JJPROC,10).NE.6) THEN
               CALL HWWARN('UPEVNT',210)
               GOTO 999
            ENDIF
            NP=NUP
C--W DECAYS
            IDUP(NP+1)=IDUP(NP-5)
            IDUP(NP+2)=IDUP(NP-4)
            IDUP(NP+3)=IDUP(NP-2)
            IDUP(NP+4)=IDUP(NP-1)
            IDUP(NP-1)=IDUP(NP)
            CALL HWVEQU(5,PUP(1,NP-5),PUP(1,NP+1))
            CALL HWVEQU(5,PUP(1,NP-4),PUP(1,NP+2))
            CALL HWVEQU(5,PUP(1,NP-2),PUP(1,NP+3))
            CALL HWVEQU(5,PUP(1,NP-1),PUP(1,NP+4))
            CALL HWVEQU(5,PUP(1,NP  ),PUP(1,NP-1))
            CALL HWVSUM(4,PUP(1,NP+1),PUP(1,NP+2),PUP(1,NP-2))
            CALL HWVSUM(4,PUP(1,NP+3),PUP(1,NP+4),PUP(1,NP))
            CALL HWUMAS(PUP(1,NP-2))
            CALL HWUMAS(PUP(1,NP))
            IDUP(NP-2)=-24*(MOD(IDUP(NP+1),2)+MOD(IDUP(NP+2),2))
            IDUP(NP)=-IDUP(NP-2)
            DO IP=NP-3,NP+4
               ISTUP(IP)=1
            ENDDO
            ISTUP(NP-2)=2
            ISTUP(NP)=2
C--TOP DECAYS
            CALL HWVSUM(4,PUP(1,NP-3),PUP(1,NP-2),PUP(1,NP-5))
            CALL HWVSUM(4,PUP(1,NP-1),PUP(1,NP  ),PUP(1,NP-4))
            CALL HWUMAS(PUP(1,NP-5))
            CALL HWUMAS(PUP(1,NP-4))
            IDUP(NP-5)=IDUP(NP-2)/4
            IDUP(NP-4)=-IDUP(NP-5)
            ISTUP(NP-5)=2
            ISTUP(NP-4)=2
            MOTHUP(1,NP-5)=1
            MOTHUP(2,NP-5)=2
            MOTHUP(1,NP-4)=1
            MOTHUP(2,NP-4)=2
            MOTHUP(1,NP-3)=NP-5
            MOTHUP(1,NP-2)=NP-5
            MOTHUP(1,NP-1)=NP-4
            MOTHUP(1,NP  )=NP-4
            MOTHUP(1,NP+1)=NP-2
            MOTHUP(1,NP+2)=NP-2
            MOTHUP(1,NP+3)=NP
            MOTHUP(1,NP+4)=NP
            DO IP=NP-3,NP+4
               MOTHUP(2,IP)=MOTHUP(1,IP)
               ICOLUP(1,IP)=0
               ICOLUP(2,IP)=0
            ENDDO
            ICOLUP(1,NP-3)=ICOLUP(1,NP-5)
            ICOLUP(2,NP-3)=ICOLUP(2,NP-5)
            ICOLUP(1,NP-1)=ICOLUP(1,NP-4)
            ICOLUP(2,NP-1)=ICOLUP(2,NP-4)
            IF (IDUP(NP+1).GT.0) THEN
               ICOLUP(1,NP+1)=505
               ICOLUP(2,NP+2)=505
            ELSE
               ICOLUP(2,NP+1)=505
               ICOLUP(1,NP+2)=505
            ENDIF
            IF (IDUP(NP+3).GT.0) THEN
               ICOLUP(1,NP+3)=506
               ICOLUP(2,NP+4)=506
            ELSE
               ICOLUP(2,NP+3)=506
               ICOLUP(1,NP+4)=506
            ENDIF
            NUP=NP+4
         ENDIF
      ELSEIF (JPR.EQ.20) THEN
         IF (JJPROC.GE.2030) THEN
C---TOP+W OR TOP+H FINAL STATE
C---SET COLOUR CONNECTIONS
            IF (IC.GT.20) THEN
               CALL HWWARN('UPEVNT',111)
               GOTO 999
            ENDIF
            DO I=1,4
               CALL UPCODE(ICOST(I,IC),ICOLUP(1,I))
            ENDDO
            IF (JJPROC.LT.2040.OR.NODEC) THEN
               ICOLUP(1,5)=0
               ICOLUP(2,5)=0
            ELSE
               ICOLUP(1,5)=ICOLUP(1,4)
               ICOLUP(2,5)=ICOLUP(2,4)
               ICOLUP(1,4)=0
               ICOLUP(2,4)=0
            ENDIF
            IF (.NOT.REMIT) THEN
C---2-BODY FINAL STATE
               DO IP=3,NUP
                  IDUP(IP)=IDUP(IP+1)
                  ICOLUP(1,IP)=ICOLUP(1,IP+1)
                  ICOLUP(2,IP)=ICOLUP(2,IP+1)
                  CALL HWVEQU(5,PUP(1,IP+1),PUP(1,IP))
               ENDDO
            ENDIF
            IF (NODEC) THEN
               IP=NUP-1
            ELSEIF (JJPROC.LT.2040) THEN
               IP=NUP-2
C---T W (OR TBAR W) WITH DECAYS
C   LES HOUCHES COMMON IS FILLED AS FOLLOWS:
C
C  1 INCOMING PARTON
C  2 INCOMING PARTON
C  3 OUTGOING LIGHT PARTON (IF ANY, OTHERWISE 4-11 BECOME 3-10)
C  4 TOP
C  5 W FROM TOP
C  6 QUARK FROM TOP
C  7 DECAY PRODUCT FROM W FROM TOP
C  8 DECAY PRODUCT FROM W FROM TOP
C  9 W NOT FROM TOP
C 10 DECAY PRODUCT FROM W NOT FROM TOP
C 11 DECAY PRODUCT FROM W NOT FROM TOP
C
C---RECONSTRUCT W DECAY
               NP=NUP
               IDUP(NP+1)=-24*(MOD(IDUP(NP-1),2)+MOD(IDUP(NP),2))
               IDUP(NP+2)=IDUP(NP-1)
               IDUP(NP+3)=IDUP(NP)
               CALL HWVEQU(5,PUP(1,NP-1),PUP(1,NP+2))
               CALL HWVEQU(5,PUP(1,NP),PUP(1,NP+3))
               CALL HWVSUM(4,PUP(1,NP-1),PUP(1,NP),PUP(1,NP+1))
               CALL HWUMAS(PUP(1,NP+1))
               ISTUP(NP+1)=2
               ISTUP(NP+2)=1
               ISTUP(NP+3)=1
               MOTHUP(1,NP+1)=1
               MOTHUP(2,NP+1)=2
               MOTHUP(1,NP+2)=NP+1
               MOTHUP(2,NP+2)=NP+1
               MOTHUP(1,NP+3)=NP+1
               MOTHUP(2,NP+3)=NP+1
               ICOLUP(1,NP+1)=0
               ICOLUP(2,NP+1)=0
               ICOLUP(1,NP+2)=505
               ICOLUP(2,NP+2)=0
               ICOLUP(1,NP+3)=0
               ICOLUP(2,NP+3)=505
               NUP=NP+3
            ELSE
               IP=NUP
            ENDIF                  
            IF (IDUP(IP).LT.0) THEN
C---CHARGE CONJUGATE PROCESS
               DO I=1,NUP
                  IA=ICOLUP(1,I)
                  ICOLUP(1,I)=ICOLUP(2,I)
                  ICOLUP(2,I)=IA
               ENDDO
            ENDIF
         ELSE
C---SINGLE TOP: IA,IB ARE THE QUARKS THAT ARE COLOUR CONNECTED
C   I.E. (FOR H EVENTS) THOSE THAT ARE NOT CONNECTED TO GLUON
            IA=IC/10
            IB=IC-10*IA
            IF (IA.LT.1.OR.IA.GT.5) THEN
               CALL HWWARN('UPEVNT',108)
            ELSEIF (IB.LT.1.OR.IB.GT.5) THEN
               CALL HWWARN('UPEVNT',109)
            ELSEIF (IA.EQ.IB) THEN
               CALL HWWARN('UPEVNT',110)
            ENDIF
            IF (IERROR.NE.0) GOTO 999
            IF (.NOT.NODEC) THEN
               ID=IDUP(5)
               IDUP(5)=IDUP(7)
            ENDIF
            DO I=1,5
               IF (I.EQ.IA.OR.I.EQ.IB) THEN
                  IF (IDUP(I).GT.0) THEN
                     ICOLUP(1,I)=501
                     ICOLUP(2,I)=0
                  ELSE
                     ICOLUP(1,I)=0
                     ICOLUP(2,I)=501
                  ENDIF
               ELSEIF (IDUP(I).EQ.21) THEN
                  IG=I
                  ICOLUP(1,I)=502
                  ICOLUP(2,I)=503
               ELSEIF (IDUP(I).GT.0) THEN
                  ICOLUP(1,I)=502
                  ICOLUP(2,I)=0
               ELSE
                  ICOLUP(1,I)=0
                  ICOLUP(2,I)=502
               ENDIF
            ENDDO
            IF (.NOT.NODEC) IDUP(5)=ID
            IF (REMIT) THEN
C---3-BODY FINAL STATE
C---COMPLETE GLUON COLOUR CONNECTIONS
               IF (.NOT.NODEC) THEN
                  ID=IDUP(5)
                  IDUP(5)=IDUP(7)
               ENDIF
               DO I=1,5
                  IF (I.NE.IA.AND.I.NE.IB.AND.I.NE.IG) THEN
                     IF (IDUP(I).GT.0) THEN
                        IF((I.LT.3.AND.IG.LT.3)
     &                       .OR.(I.GT.2.AND.IG.GT.2)) ICOLUP(1,I)=503
                     ELSE
                        IF((I.LT.3.AND.IG.GT.2)
     &                       .OR.(I.GT.2.AND.IG.LT.3)) ICOLUP(2,I)=503
                     ENDIF
                  ENDIF
               ENDDO
               IF (.NOT.NODEC) IDUP(5)=ID
            ELSE
C---2-BODY FINAL STATE
               DO IP=3,NUP
                  IDUP(IP)=IDUP(IP+1)
                  CALL HWVEQU(5,PUP(1,IP+1),PUP(1,IP))
               ENDDO
C---SET COLOUR CONNECTIONS
               ICOLUP(1,3)=ICOLUP(1,4)
               ICOLUP(2,3)=ICOLUP(2,4)
               ICOLUP(1,4)=ICOLUP(1,5)
               ICOLUP(2,4)=ICOLUP(2,5)
            ENDIF
         ENDIF
         IF (.NOT.NODEC) THEN
            IF (JJPROC.GE.2030.AND.JJPROC.LT.2040) THEN
               NP=NUP-5
            ELSE
               NP=NUP
C---SINGLE TOP WITH DECAYS
C   LES HOUCHES COMMON IS FILLED AS FOLLOWS:
C
C  1 INCOMING PARTON
C  2 INCOMING PARTON
C  3 OUTGOING LIGHT PARTON (IF ANY, OTHERWISE 4-9 BECOME 3-8)
C  4 OUTGOING BBAR-TYPE QUARK
C  5 TOP
C  6 W FROM TOP
C  7 QUARK FROM TOP
C  8 DECAY PRODUCT FROM W FROM TOP
C  9 DECAY PRODUCT FROM W FROM TOP
C
C---RECONSTRUCT TOP DECAY
            ENDIF
            IDUP(NP+1)=IDUP(NP-2)
            IDUP(NP+2)=IDUP(NP-1)
            IDUP(NP-1)=-24*(MOD(IDUP(NP+1),2)+MOD(IDUP(NP+2),2))
            IDUP(NP-2)=IDUP(NP-1)/4
            CALL HWVEQU(5,PUP(1,NP-2),PUP(1,NP+1))
            CALL HWVEQU(5,PUP(1,NP-1),PUP(1,NP+2))
            CALL HWVSUM(4,PUP(1,NP+1),PUP(1,NP+2),PUP(1,NP-1))
            CALL HWVSUM(4,PUP(1,NP-1),PUP(1,NP  ),PUP(1,NP-2))
            CALL HWUMAS(PUP(1,NP-1))
            CALL HWUMAS(PUP(1,NP-2))
            DO IP=NP-3,NP+2
               ISTUP(IP)=1
            ENDDO
            ISTUP(2)=-1
            ISTUP(NP-1)=2
            ISTUP(NP-2)=2
            MOTHUP(1,NP-3)=1
            MOTHUP(2,NP-3)=2
            MOTHUP(1,NP-2)=1
            MOTHUP(2,NP-2)=2
            MOTHUP(1,NP-1)=NP-2
            MOTHUP(1,NP  )=NP-2
            MOTHUP(1,NP+1)=NP-1
            MOTHUP(1,NP+2)=NP-1
            DO IP=NP-1,NP+2
               MOTHUP(2,IP)=MOTHUP(1,IP)
               ICOLUP(1,IP)=0
               ICOLUP(2,IP)=0
            ENDDO
            ICOLUP(1,NP)=ICOLUP(1,NP-2)
            ICOLUP(2,NP)=ICOLUP(2,NP-2)
            IF (IDUP(NP+1).GT.0) THEN
               ICOLUP(1,NP+1)=504
               ICOLUP(2,NP+2)=504
            ELSE
               ICOLUP(2,NP+1)=504
               ICOLUP(1,NP+2)=504
            ENDIF
            IF (JJPROC.LT.2030.OR.JJPROC.GT.2039) NUP=NP+2
         ENDIF
       ELSEIF (JPR.EQ.1) THEN
C---E+E- ANNIHILATION
         DO I=1,NUP
            ICOLUP(1,I)=0
            ICOLUP(2,I)=0
         ENDDO
C---RESCALE 3-MOMENTA TO PUT PARTONS ON-SHELL
         PUP(5,1)=RMASS(121)
         PUP(5,2)=PUP(5,1)
         CALL HWURSC(2,PUP)
         PUP(5,3)=RMASS(13)
         PUP(5,4)=RMASS(ABS(IDUP(4)))
         PUP(5,5)=PUP(5,4)
         IF (REMIT) THEN
            CALL HWURSC(3,PUP(1,3))
            ICOLUP(1,3)=501
            ICOLUP(2,3)=502
            IF (IDUP(4).GT.0) THEN
               ICOLUP(1,4)=502
               ICOLUP(2,5)=501
            ELSE
               ICOLUP(2,4)=501
               ICOLUP(1,5)=502
            ENDIF
         ELSE
            CALL HWURSC(2,PUP(1,4))
            DO I=3,4
               CALL HWVEQU(5,PUP(1,I+1),PUP(1,I))
               IDUP(I)=IDUP(I+1)
               ISTUP(I)=1
            ENDDO
            IF (IDUP(3).GT.0) THEN
               ICOLUP(1,3)=501
               ICOLUP(2,4)=501
            ELSE
               ICOLUP(2,3)=501
               ICOLUP(1,4)=501
            ENDIF
         ENDIF
      ELSE
C---BOSON PAIR OR LEPTON PAIR
         IF (BOPRO.OR.NODEC) THEN
            NUP=6
            DO I=6,5,-1
               CALL HWVEQU(5,PUP(1,I-1),PUP(1,I))
               IDUP(I)=IDUP(I-1)
               ISTUP(I)=1
            ENDDO
         ELSE
C---BOSON PAIR: ONE OR BOTH DECAYED
C---ADD BOSON(S) TO EVENT RECORD
            IF (ABS(IDUP(6)).LT.20) THEN
               NUP=8
               I=2
               IF (ABS(IDUP(4)).LT.20) THEN
                  NUP=10
                  I=3
               ENDIF
               MUP=NUP-1
               CALL HWVEQU(5,PUP(1,MUP-I),PUP(1,MUP))
               CALL HWVEQU(5,PUP(1,NUP-I),PUP(1,NUP))
               CALL HWVSUM(4,PUP(1,MUP),PUP(1,NUP),PUP(1,6))
               CALL HWUMAS(PUP(1,6))
               IDUP(MUP)=IDUP(MUP-I)
               IDUP(NUP)=IDUP(NUP-I)
               ISTUP(MUP)=1
               ISTUP(NUP)=1
               MOTHUP(1,MUP)=6
               MOTHUP(2,MUP)=6
               MOTHUP(1,NUP)=6
               MOTHUP(2,NUP)=6
               ISTUP(6)=2
               ID=MOD(IDUP(MUP),2)+MOD(IDUP(NUP),2)
               IF (ID.EQ.0) THEN
                  IDUP(6)=23
               ELSEIF (ABS(ID).EQ.1) THEN
                  IDUP(6)=-24*ID
               ELSE
                  CALL HWWARN('UPEVNT',106)
                  GOTO 999            
               ENDIF
            ENDIF
            IF (ABS(IDUP(4)).LT.20) THEN
               CALL HWVZRO(4,PDB)
               DO I=8,7,-1
                  CALL HWVEQU(5,PUP(1,I-3),PUP(1,I))
                  CALL HWVSUM(4,PUP(1,I),PDB,PDB)
                  IDUP(I)=IDUP(I-3)
                  ISTUP(I)=1
                  MOTHUP(1,I)=5
                  MOTHUP(2,I)=5
               ENDDO
               CALL HWUMAS(PDB)
               CALL HWVEQU(5,PDB,PUP(1,5))
               ISTUP(5)=2
               ID=MOD(IDUP(7),2)+MOD(IDUP(8),2)
               IF (ID.EQ.0) THEN
                  IDUP(5)=23
               ELSEIF (ABS(ID).EQ.1) THEN
                  IDUP(5)=-24*ID
               ELSE
                  CALL HWWARN('UPEVNT',107)
                  GOTO 999            
               ENDIF
            ELSE
               CALL HWVEQU(5,PUP(1,4),PUP(1,5))
               IDUP(5)=IDUP(4)
               ISTUP(5)=1
               MOTHUP(1,5)=4
               MOTHUP(2,5)=4
            ENDIF
         ENDIF
C---ADD DIBOSON OR DILEPTON TO EVENT RECORD (TO FIX ITS MASS)
         CALL HWVZRO(4,PDB)
         DO I=6,5,-1
            CALL HWVSUM(4,PUP(1,I),PDB,PDB)
            MOTHUP(1,I)=4
            MOTHUP(2,I)=4
         ENDDO
         CALL HWUMAS(PDB)
         CALL HWVEQU(5,PDB,PUP(1,4))
         ISTUP(4)=2
         IDUP(4)=0
         IF (REMIT) THEN
C---SET COLOUR CONNECTIONS
            DO I=1,3
               ICOLUP(1,I)=501
               ICOLUP(2,I)=502
            ENDDO
            IF (IHPRO.EQ.1) THEN
               ICOLUP(2,1)=0
               ICOLUP(1,2)=0
            ELSEIF (IHPRO.EQ.2) THEN
               ICOLUP(1,1)=502
               ICOLUP(2,1)=0
               ICOLUP(2,3)=0
            ELSEIF (IHPRO.EQ.3) THEN
               ICOLUP(1,1)=0
               ICOLUP(2,2)=0
            ELSEIF (IHPRO.EQ.4) THEN
               ICOLUP(1,1)=0
               ICOLUP(2,1)=501
               ICOLUP(1,3)=0
            ELSEIF (IHPRO.EQ.5) THEN
               ICOLUP(1,2)=502
               ICOLUP(2,2)=0
               ICOLUP(2,3)=0
            ELSEIF (IHPRO.EQ.6) THEN
               ICOLUP(1,2)=0
               ICOLUP(2,2)=501
               ICOLUP(1,3)=0
            ELSE
               CALL HWWARN('UPEVNT',102)
               GOTO 999
            ENDIF
            DO I=4,NUP
               ICOLUP(1,I)=0
               ICOLUP(2,I)=0
            ENDDO
         ELSE
            DO I=5,NUP
               CALL HWVEQU(5,PUP(1,I),PUP(1,I-2))
               IDUP(I-2)=IDUP(I)
               ISTUP(I-2)=ISTUP(I)
               MOTHUP(1,I-2)=MOTHUP(1,I)-2
               MOTHUP(2,I-2)=MOTHUP(1,I)-2
            ENDDO
            MOTHUP(1,3)=1
            MOTHUP(1,4)=1
            NUP=NUP-2
C---SET COLOUR CONNECTIONS
            DO I=1,NUP
               ICOLUP(1,I)=0
               ICOLUP(2,I)=0
            ENDDO
            IF (IDUP(1).GT.0) THEN
               ICOLUP(1,1)=501
               ICOLUP(2,2)=501
            ELSE
               ICOLUP(2,1)=501
               ICOLUP(1,2)=501
            ENDIF
         ENDIF
         IF (BOPRO) THEN
C---DILEPTON PRODUCTION
            IBOS=MOD(JJPROC,100)
            ILEP=MOD(JJPROC,10)
            IBOS=IBOS-ILEP
C---LOAD LEPTON AND BOSON ID
            I=NUP-1
            J=NUP
            IF ( IBOS.EQ.50 .OR.
     #          (IBOS.EQ.60.AND.JPR.EQ.13) .OR.
     #          (IBOS.EQ.70.AND.JPR.EQ.13) ) THEN
               IDUP(I)=-ILEP-10
               IDUP(J)=-IDUP(I)
               IF (REMIT) IDUP(4)=23
            ELSEIF (IBOS.EQ.60.AND.JPR.EQ.14) THEN
               IDUP(I)=-9-2*ILEP
               IDUP(J)=1-IDUP(I)
               IF (REMIT) IDUP(4)=24
            ELSEIF (IBOS.EQ.70.AND.JPR.EQ.14) THEN
               IDUP(I)=-10-2*ILEP
               IDUP(J)=-1-IDUP(I)
               IF (REMIT) IDUP(4)=-24
            ELSE
               CALL HWWARN('UPEVNT',504)
            ENDIF
         ENDIF
      ENDIF
C---INITIAL STATE
      DO I=1,2
         ISTUP(I)=-1
         MOTHUP(1,I)=0
         MOTHUP(2,I)=0
      ENDDO
 999  CONTINUE
      IF(IERROR.LT.100) RETURN
      PRINT *
      DO I=1,NUP
         PRINT '(4I4,5F8.2)',IDUP(I),ISTUP(I),(ICOLUP(J,I),J=1,2),
     &        (PUP(J,I),J=1,5)
      ENDDO
c       IPR, IC, NP
 901  FORMAT(1X,I3,2(1X,I2))
c      (ID(I),I=1,NP)
 902  FORMAT(9(1X,I3))
c       XEVWGT,EMSCA,XEPHO
 903  FORMAT(3(1X,D14.8))
c      ((P(J,I),J=1,4),I=1,NP)
 904  FORMAT(36(1X,D14.8))
c      X1, X2, Q2
 905  FORMAT(3(1X,D14.8))
c      WGTACP(1)..WGTACP(10)
 906  FORMAT(10(1X,D14.8))
c      WGTACP(1)..WGTACP(28)
 907  FORMAT(28(1X,D14.8))
      END
C----------------------------------------------------------------------
      SUBROUTINE UPCODE(ICODE,ICOL)
C--DECODES COLOUR CONNECTIONS
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ICODE,ICOL(2)
      ICOL(1)=ICODE/10
      IF (ICOL(1).NE.0) ICOL(1)=ICOL(1)+500
      ICOL(2)=MOD(ICODE,10)
      IF (ICOL(2).NE.0) ICOL(2)=ICOL(2)+500
      END
C----------------------------------------------------------------------
      SUBROUTINE UPINIT
C----------------------------------------------------------------------
C  Reads MC@NLO input headers and fills Les Houches run common HEPRUP
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
C--Les Houches Common Blocks
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &                IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &                XMAXUP(MAXPUP),LPRUP(MAXPUP)
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
      DOUBLE PRECISION XCKECM,XCKPB1,XTMP1,XTMP2,XTMP3,XTMP4,XMT,XMW,
     & XMZ,XMH,XMV,XM1,XM2,XM3,XM4,XM5,XM21,XLAM,GAH,GAT,GAW,XMB,TINY
      DOUBLE PRECISION XMV1,GAV1,GAMAX1,XMV2,GAV2,GAMAX2
      INTEGER IVVCODE,IFAIL,MQQ,NQQ,IHW,I,NDNS1,NDNS2,JPR,JPR0,IH,IHQ,
     & IVHVEC,IVHLEP,IVLEP1,IVLEP2
      CHARACTER*60 TMPSTR,TMPSTR2
      CHARACTER*4 STRP1,STRP2
      CHARACTER*8 STRGRP1,STRGRP2
      CHARACTER*2 STRSCH
      CHARACTER*3 STRFMT
      CHARACTER*50 QQIN
      LOGICAL FK88STRNOEQ,OLDFORM
      DATA TINY/1.D-3/
      COMMON/NQQCOM/MQQ,NQQ
      COMMON/VVJIN/QQIN
      COMMON/VHLIN/IVHVEC,IVHLEP
      COMMON/VVLIN/IVLEP1,IVLEP2
C
      IF (IERROR.NE.0) RETURN
      OLDFORM=.FALSE.
C--SET UP INPUT FILES
      OPEN(UNIT=61,FILE=QQIN,STATUS='UNKNOWN')
C--READ HEADERS OF EVENT FILE
      READ(61,801)XCKECM,XTMP1,XTMP2,XTMP3,XTMP4,TMPSTR
      READ(61,802)IVVCODE,TMPSTR
      IVVCODE=MOD(IVVCODE,10000)
C---CHECK PROCESS CODE
      JPR0=MOD(ABS(IPROC),10000)
      JPR=JPR0/100
      IF (JPR.NE.IVVCODE/100) CALL HWWARN('UPINIT',500)
      IF ((JPR.EQ.17.OR.JPR.EQ.28.OR.JPR.EQ.36).AND.
     & IVVCODE.NE.MOD(ABS(IPROC),10000)) CALL HWWARN('UPINIT',501)
      IF (JPR.EQ.13.OR.JPR.EQ.14) THEN
         IF(JPR0.EQ.1396.OR.JPR0.EQ.1371.OR.
     #      JPR0.EQ.1372.OR.JPR0.EQ.1373)THEN
           READ(61,808)EMMIN,EMMAX,TMPSTR
         ELSE
           READ(61,809)XMV,GAH,GAMMAX,TMPSTR
         ENDIF
C-- CHECK VECTOR BOSON MASS
         IF( (IVVCODE.EQ.1397.AND.ABS(1-XMV/RMASS(200)).GT.TINY) .OR.
     #       (IVVCODE.EQ.1497.AND.ABS(1-XMV/RMASS(198)).GT.TINY) .OR.
     #       (IVVCODE.EQ.1498.AND.ABS(1-XMV/RMASS(199)).GT.TINY) )
     #      CALL HWWARN('UPINIT',502)
      ELSEIF (JPR.EQ.26.OR.JPR.EQ.27.OR.JPR.EQ.19) THEN
         READ(61,810)IVHVEC,IVHLEP,TMPSTR
         READ(61,809)XMV,GAH,GAMMAX,TMPSTR
         READ(61,809)XMH,GAH,GAMMAX,TMPSTR
         IF( (JPR.EQ.26.AND.ABS(1-XMV/RMASS(199)).GT.TINY) .OR.
     #       (JPR.EQ.27.AND.ABS(1-XMV/RMASS(200)).GT.TINY) )
     #      CALL HWWARN('UPINIT',508)
         IF(ABS(1-XMH/RMASS(201)).GT.TINY) CALL HWWARN('UPINIT',509)
      ELSEIF (JPR.EQ.28) THEN
         READ(61,808)XMW,XMZ,TMPSTR
C-- CHECK VECTOR BOSON MASSES
         IF(ABS(1-XMW/RMASS(198)).GT.TINY .OR.
     #      ABS(1-XMZ/RMASS(200)).GT.TINY) CALL HWWARN('UPINIT',502)
         READ(61,810)IVLEP1,IVLEP2,TMPSTR
         READ(61,809)XMV1,GAV1,GAMAX1,TMPSTR
         READ(61,809)XMV2,GAV2,GAMAX2,TMPSTR
      ELSEIF (JPR.EQ.16.OR.JPR.EQ.36) THEN
         READ(61,813)XMH,GAH,XMT,XMB,TMPSTR
C-- CHECK HIGGS AND TOP MASSES
         IH=201
         IF (JPR.EQ.36) IH=IVVCODE/10-158
         IF(ABS(1-XMH/RMASS(IH)).GT.TINY) CALL HWWARN('UPINIT',503)
         IF(ABS(1-XMT/RMASS(6)) .GT.TINY) CALL HWWARN('UPINIT',504)
      ELSEIF (JPR.EQ.17.OR.JPR.EQ.51) THEN
         IHQ=MOD(JPR0,10)
         IF (IHQ.EQ.6) THEN
            READ(61,'(A)') TMPSTR2
            IF (TMPSTR2(17:19).EQ.'M_Q') THEN
               OLDFORM=.TRUE.
               READ(TMPSTR2,803)XMT,TMPSTR
            ELSE
               READ(TMPSTR2,808)XMT,GAT,TMPSTR
            ENDIF
         ELSE
            READ(61,803)XMT,TMPSTR
         ENDIF
C-- CHECK HEAVY QUARK MASS
         IF(ABS(1-XMT/RMASS(IHQ)).GT.TINY) CALL HWWARN('UPINIT',505)
         IF (IHQ.EQ.6) THEN
            IF(.NOT.OLDFORM)THEN
               READ(61,808)XMW,GAW,TMPSTR
               READ(61,810)IVLEP1,IVLEP2,TMPSTR
C-- CHECK W BOSON MASS WHEN TOPS DECAY
               IF( IVLEP1.NE.7.AND.IVLEP2.NE.7 .AND.
     #             ABS(1-XMW/RMASS(198)).GT.TINY ) 
     #            CALL HWWARN('UPINIT',502)
            ELSE
               XMW=0.D0
               GAW=0.D0
               IVLEP1=7
               IVLEP2=7
            ENDIF
         ENDIF
      ELSEIF (JPR.EQ.20) THEN
         READ(61,'(A)') TMPSTR2
         IF (TMPSTR2(28:43).EQ.'M_top, Gamma_top') THEN
            READ(TMPSTR2,808)XMT,GAT,TMPSTR
         ELSE
            OLDFORM=.TRUE.
            READ(TMPSTR2,803)XMT,TMPSTR
         ENDIF
C-- CHECK TOP QUARK MASS
         IF(ABS(1-XMT/RMASS(6)).GT.TINY) CALL HWWARN('UPINIT',511)
         IF(OLDFORM)GOTO 444
         READ(61,808)XMW,GAW,TMPSTR
         IF(JPR0.LT.2030.OR.JPR0.GT.2039)THEN
            READ(61,812)IVLEP1,TMPSTR
C-- CHECK W BOSON MASS WHEN TOPS DECAY
            IF(JPR0.LT.2030 .AND. IVLEP1.NE.7 .AND.
     #        ABS(1-XMW/RMASS(198)).GT.TINY ) CALL HWWARN('UPINIT',502)
            IF(JPR0.GT.2039) THEN
C--CHARGED HIGGS: IVLEP2 IS DUMMY
              IVLEP2=IVLEP1
C--XMW HERE MEANS CHARGED HIGGS MASS
              IF(ABS(1-XMW/RMASS(206)).GT.TINY ) 
     #          CALL HWWARN('UPINIT',502)
           ENDIF
         ELSE
            READ(61,810)IVLEP1,IVLEP2,TMPSTR
C-- CHECK W BOSON MASS
            IF(ABS(1-XMW/RMASS(198)).GT.TINY) CALL HWWARN('UPINIT',502)
         ENDIF
 444     CONTINUE
      ELSEIF (JPR.NE.1) THEN
         CALL HWWARN('UPINIT',506)
      ENDIF
      READ(61,804)XM1,XM2,XM3,XM4,XM5,XM21,TMPSTR
      IF (JPR.NE.1) THEN
         READ(61,805)STRP1,STRP2,TMPSTR
         READ(61,806)STRGRP1,NDNS1,TMPSTR
         IF (JPR.EQ.51) THEN
            READ(61,806)STRGRP2,NDNS2,TMPSTR
            READ(61,803)XCKPB1,TMPSTR
            IF(ABS(XCKPB1-PBEAM1).GT.TINY) CALL HWWARN('UPINIT',512)
         ELSE
            STRGRP2=STRGRP1
            NDNS2=NDNS1
         ENDIF
      ENDIF
      READ(61,807)XLAM,STRSCH,TMPSTR
C--CHECK THAT EVENT FILE HAS BEEN GENERATED CONSISTENTLY WITH 
C--HERWIG PARAMETERS ADOPTED HERE
      IFAIL=0
C-- CM ENERGY
      IF( ABS(XCKECM-2D0*SQRT(PBEAM1*PBEAM2)).GT.TINY .OR.
C-- QUARK AND GLUON MASSES
     #     ABS(XM1-RMASS(1)).GT.TINY .OR.
     #     ABS(XM2-RMASS(2)).GT.TINY .OR.
     #     ABS(XM3-RMASS(3)).GT.TINY .OR.
     #     ABS(XM4-RMASS(4)).GT.TINY .OR.
     #     ABS(XM5-RMASS(5)).GT.TINY .OR.
     #     ABS(XM21-RMASS(13)).GT.TINY) IFAIL=1
C-- LAMBDA_QCD: NOW REMOVED TO ALLOW MORE FLEXIBILITY (NNLO EFFECT ANYHOW)
C     #     ABS(XLAM-QCDLAM).GT.TINY .OR.
C-- REPLACE THE FOLLOWING WITH A CONDITION ON STRSCH, IF CONSISTENT 
C-- INFORMATION ON PDF SCHEME WILL BE AVAILABLE FROM PDF LIBRARIES AND HERWIG
C-- COLLIDING PARTICLE TYPE
      IF (JPR.NE.1.AND.IFAIL.EQ.0) THEN
         IF(
     #        FK88STRNOEQ(STRP1,PART1) .OR.
     #        FK88STRNOEQ(STRP2,PART2) )IFAIL=1
C--IF PDF LIBRARY IS USED, CHECK PDF CONSISTENCY
         IF( IFAIL.EQ.0 .AND. MODPDF(1).NE.-1)THEN
            IF( 
     #          (FK88STRNOEQ(STRGRP1,AUTPDF(1)) .AND.
     #          (STRGRP1.EQ.'LHAPDF'.AND.AUTPDF(1).NE.'HWLHAPDF')) .OR.
     #          (FK88STRNOEQ(STRGRP2,AUTPDF(2)) .AND.
     #          (STRGRP2.EQ.'LHAPDF'.AND.AUTPDF(2).NE.'HWLHAPDF')) .OR.
     #          ABS(NDNS1-MODPDF(1)).GT.TINY .OR.
     #          ABS(NDNS2-MODPDF(2)).GT.TINY )IFAIL=1
C--WHEN LHAPDF IS LINKED, AUTPDF() IS A MC@NLO-DEFINED STRING
            IF(AUTPDF(1).EQ.'LHAPDF'.OR.AUTPDF(1).EQ.'LHAEXT'.OR.
     #         AUTPDF(1).EQ.'HWLHAPDF')THEN
               AUTPDF(1)='DEFAULT'
               AUTPDF(2)='DEFAULT'
            ENDIF
         ENDIF
      ENDIF
      IF(IFAIL.EQ.1) CALL HWWARN('UPINIT',507)
      CALL HWUIDT(3,IDBMUP(1),IHW,PART1)
      CALL HWUIDT(3,IDBMUP(2),IHW,PART2)
      EBMUP(1)=PBEAM1
      EBMUP(2)=PBEAM2
      DO I=1,2
         PDFGUP(I)=-1
         PDFSUP(I)=-1
      ENDDO
      IDWTUP=-4
      NPRUP=1
      LPRUP(1)=IVVCODE
C-- TEST FOR NEW FORMAT INPUT MOMENTA: (PX,PY,PZ,M)
      READ(61,811) STRFMT,TMPSTR
      IF (STRFMT.NE.'P,M') CALL HWWARN('UPINIT',510)
      READ(61,900) MQQ
      NQQ=0
C-- LARGEST EXPECTED NUMBER OF LEGS
      NUP=10
      AQEDUP=ZERO
      AQCDUP=ZERO
      DO I=1,NUP
         VTIMUP(I)=ZERO
         SPINUP(I)=9.
      ENDDO
 801  FORMAT(5(1X,D10.4),1X,A)
 802  FORMAT(1X,I6,1X,A)
 803  FORMAT(1X,D10.4,1X,A)
 804  FORMAT(6(1X,D10.4),1X,A)
 805  FORMAT(2(1X,A4),1X,A)
 806  FORMAT(1X,A8,1X,I6,1X,A)
 807  FORMAT(1X,D10.4,1X,A2,1X,A)
 808  FORMAT(2(1X,D10.4),1X,A)
 809  FORMAT(3(1X,D10.4),1X,A)
 810  FORMAT(2(1X,I2),1X,A)
 811  FORMAT(1X,A3,1X,A)
 812  FORMAT(1X,I2,1X,A)
 813  FORMAT(4(1X,D10.4),1X,A)
 900  FORMAT(I9)
 999  END


C----------------------------------------------------------------------
      SUBROUTINE HWURSC(NP,PP)
C  RESCALES A SET OF NP (<21) 3-MOMENTA PP(1-3,*) IN
C  THEIR CMF TO PUT PP ON MASS-SHELL AT MASSES PP(5,*) 
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NP,IP,IT,NT
      DOUBLE PRECISION PP(5,*),P(5,20),P2(20),M2(20),SP(5),
     & TINY,FAC,ECM,DCM,EP,STEP,FRT,HWUSQR
      DATA TINY,NT/1D-9,20/
      IF (NP.GT.20) CALL HWWARN('HWURSC',300+NP)
C--COMPUTE CM MOMENTUM
      CALL HWVZRO(4,SP)
      DO IP=1,NP
         CALL HWVSUM(4,PP(1,IP),SP,SP)
      ENDDO
      CALL HWUMAS(SP)
C--BOOST TO CMF
      DO IP=1,NP
         CALL HWULOF(SP,PP(1,IP),P(1,IP))
         P2(IP)=P(1,IP)**2+P(2,IP)**2+P(3,IP)**2
         M2(IP)=P(5,IP)**2
      ENDDO
C--ITERATE RESCALING OF 3-MOMENTA
      FAC=1D0
      DO IT=1,NT
         ECM=0D0
         DCM=0D0
         DO IP=1,NP
            EP=HWUSQR(M2(IP)+FAC*P2(IP))
            IF (EP.GT.0D0) THEN
               ECM=ECM+EP
               DCM=DCM+P2(IP)/EP
            ENDIF
         ENDDO
         IF (DCM.EQ.0D0) CALL HWWARN('HWURSC',390)
         STEP=2D0*(ECM-SP(5))/DCM
         FAC=FAC-STEP
         IF (ABS(STEP).LT.TINY) GOTO 100
      ENDDO
C--FAILED TO CONVERGE
      CALL HWWARN('HWURSC',1)
C--CONVERGED: RESCALE 3-MOMENTA AND BOOST BACK 
 100  IF (FAC.LT.0D0) CALL HWWARN('HWURSC',391)
      FRT=SQRT(FAC)
      DO IP=1,NP
         CALL HWVSCA(3,FRT,P(1,IP),P(1,IP))
         P(4,IP)=SQRT(M2(IP)+FAC*P2(IP))
         CALL HWULOB(SP,P(1,IP),PP(1,IP))
      ENDDO
      END
