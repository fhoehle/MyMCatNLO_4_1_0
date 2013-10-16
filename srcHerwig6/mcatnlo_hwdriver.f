C Driver, to be modified by the user (keeping the read statements).
C Should be backward compatible with v3.1 and v3.2 (thanks to B. Kersevan)
      PROGRAM HWIGPR
      INCLUDE 'HERWIG65.INC'
      INTEGER N,NSTEP,I,JPR0,JPR
C QQIN IS THE EVENT FILE
      CHARACTER*50 QQIN
      COMMON/VVJIN/QQIN
      CHARACTER*60 TMPSTR2
      REAL*8 TMPLAM,GAMW0,GAMZ0,VGAMMAX,V2GAMMAX
      INTEGER IPDF
      CHARACTER * 70 LHAPDF
      LOGICAL LHACRTL,OLDFORM
      PARAMETER (LHACRTL=.TRUE.)
c
      OLDFORM=.FALSE.
      WRITE(*,*)'Enter filename for events'
      READ(*,*)QQIN
      WRITE(*,*)'Enter maximum number of events to generate'
      WRITE(*,*)'MUST coincide with the number of events on tape'
      READ(*,*)MAXEV
      WRITE(*,*)'Enter 0 to use Herwig default PDFs'
      WRITE(*,*)'      1 to use PDFs from library'
      READ(*,*)IPDF
C OUTPUT THE RESULTS AFTER ANY NSTEP EVENTS
      NSTEP=20000
C---BEAM PARTICLES
      WRITE(*,*)'Enter colliding particles (PART1, PART2)'
      READ(*,*)PART1,PART2
C---BEAM MOMENTA
      WRITE(*,*)'Enter beam momenta (PBEAM1, PBEAM2)'
      READ(*,*)PBEAM1,PBEAM2
C---PROCESS
      WRITE(*,*)'Enter process number (IPROC) (IPROC<0 for Les Houches)'
      READ(*,*)IPROC
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
C
C************************************************************************
C---UNCOMMENT THE ASSIGNMENT PRESPL=.FALSE. WHEN USING HERWIG VERSION 6.52 
C---OR HIGHER (SEE MC@NLO MANUAL, APPENDIX A.8, PAGE 25)
      PRESPL=.FALSE.
C************************************************************************
C LHSOFT=.FALSE. SWITCHES UNDERLYING EVENT OFF
      LHSOFT=ABS(IPROC).LE.10000
      WRITE(*,*)
      IF (LHSOFT) THEN
         WRITE(*,*)'Underlying event WILL be generated'
      ELSE
         WRITE(*,*)'Underlying event WILL NOT be generated'
      ENDIF
      WRITE(*,*)
C
C Select W/Z boson decay modes
      MODBOS(1)=2
      MODBOS(2)=2
C
      IF(IPDF.EQ.1)THEN
         DO I=1,2
            WRITE(*,*)'   Incoming particle # ',I
            WRITE(*,*)'Enter PDF group name (AUTPDF)'
            READ(*,*)AUTPDF(I)
            WRITE(*,*)'Enter PDF set number (MODPDF)'
            READ(*,*)MODPDF(I)
         ENDDO
C---SET LHACRTL TO FALSE IF LHAPDF DEFAULTS FOR 16, 18, AND 19 ARE OK
         IF(LHACRTL.AND.
     #      (AUTPDF(1).EQ.'LHAPDF'.OR.AUTPDF(1).EQ.'LHAEXT'))THEN
            LHAPDF='FREEZE'
            IF(AUTPDF(1).EQ.'LHAEXT')LHAPDF='EXTRAPOLATE'
            CALL SETLHACBLK(LHAPDF)
C---MODERN VERSIONS OF LHAPDF REQUIRE THE FOLLOWING SETTING
            DO I=1,2
               AUTPDF(I)='HWLHAPDF'
            ENDDO
         ENDIF
      ENDIF
      WRITE(*,*)'Enter Lambda_QCD, <0 for Herwig default'
      READ(*,*)TMPLAM
      IF(TMPLAM.GE.0.D0)QCDLAM=TMPLAM
C---CHECK PROCESS CODE
      JPR0=MOD(ABS(IPROC),10000)
      JPR=JPR0/100
      IF (JPR.EQ.13.OR.JPR.EQ.14) THEN
         IF(JPR0.EQ.1396)THEN
           WRITE(*,*)'Enter M_GAMMA*(INF), M_GAMMA*(SUP)'
           READ(*,*)EMMIN,EMMAX
         ELSEIF(JPR0.EQ.1397)THEN
           WRITE(*,*)'Enter Z0 mass, width, and GammaX'
           READ(*,*)RMASS(200),GAMZ0,GAMMAX
           IF(GAMZ0.NE.0.D0)GAMZ=GAMZ0
         ELSEIF(JPR0.EQ.1497.OR.JPR0.EQ.1498)THEN
           WRITE(*,*)'Enter W mass, width, and GammaX'
           READ(*,*)RMASS(198),GAMW0,GAMMAX
           RMASS(199)=RMASS(198)
           IF(GAMW0.NE.0.D0)GAMW=GAMW0
         ELSEIF( (JPR0.GE.1350.AND.JPR0.LE.1356) .OR.
     #           (JPR0.GE.1361.AND.JPR0.LE.1366) )THEN
           WRITE(*,*)'Enter Z0 mass, width'
           READ(*,*)RMASS(200),GAMZ
           WRITE(*,*)'Enter GAMMAX, M_Z*(INF), M_Z*(SUP)'
           READ(*,*)GAMMAX,EMMIN,EMMAX
           IF(IPROC.GT.0)THEN
             EMMIN=RMASS(200)-GAMZ*GAMMAX
             EMMAX=RMASS(200)+GAMZ*GAMMAX
           ENDIF
         ELSEIF(JPR0.GE.1371.AND.JPR0.LE.1373)THEN
           WRITE(*,*)'Enter M_LL(INF), M_LL(SUP)'
           READ(*,*)EMMIN,EMMAX
         ELSEIF( (JPR0.GE.1450.AND.JPR0.LE.1453) .OR.
     #           (JPR0.GE.1461.AND.JPR0.LE.1463) .OR.
     #           (JPR0.GE.1471.AND.JPR0.LE.1473) )THEN
           WRITE(*,*)'Enter W mass, width'
           READ(*,*)RMASS(198),GAMW
           RMASS(199)=RMASS(198)
           WRITE(*,*)'Enter GAMMAX, M_W*(INF), M_W*(SUP)'
           READ(*,*)GAMMAX,EMMIN,EMMAX
         ENDIF
      ELSEIF (JPR.EQ.28) THEN
         WRITE(*,*)'Enter W mass, width'
         READ(*,*)RMASS(198),GAMW0
         RMASS(199)=RMASS(198)
         WRITE(*,*)'Enter Z mass, width'
         READ(*,*)RMASS(200),GAMZ0
         WRITE(*,*)'Enter VGAMMAX, V2GAMMAX'
         READ(*,*)VGAMMAX,V2GAMMAX
         IF(GAMW0.NE.0.D0)GAMW=GAMW0
         IF(GAMZ0.NE.0.D0)GAMZ=GAMZ0
         GAMMAX=MAX(VGAMMAX,V2GAMMAX)
      ELSEIF (JPR.EQ.16) THEN
         WRITE(*,*)'Enter Higgs boson and top masses'
         READ(*,*)RMASS(201),RMASS(6)
      ELSEIF (JPR.EQ.17) THEN
         IF(ABS(IPROC).EQ.1705.OR.ABS(IPROC).EQ.11705)THEN
           WRITE(*,*)'Enter bottom mass'
           READ(*,*)RMASS(5)
         ELSEIF(ABS(IPROC).EQ.1706.OR.ABS(IPROC).EQ.11706)THEN
           WRITE(*,*)'Enter top mass, W mass'
           READ(*,'(A)') TMPSTR2
           READ(TMPSTR2,*,ERR=616)RMASS(6),RMASS(198)
           RMASS(199)=RMASS(198)
           GOTO 617
 616       OLDFORM=.TRUE.
           READ(TMPSTR2,*) RMASS(6)
 617       CONTINUE
         ENDIF
      ELSEIF (JPR.EQ.26) THEN
         WRITE(*,*)'Enter W mass, width'
         READ(*,*)RMASS(198),GAMW0
         RMASS(199)=RMASS(198)
         WRITE(*,*)'Enter Higgs boson mass'
         READ(*,*)RMASS(201)
         IF(GAMW0.NE.0.D0)GAMW=GAMW0
      ELSEIF (JPR.EQ.27) THEN
         WRITE(*,*)'Enter Z mass, width'
         READ(*,*)RMASS(200),GAMZ0
         WRITE(*,*)'Enter Higgs boson mass'
         READ(*,*)RMASS(201)
         IF(GAMZ0.NE.0.D0)GAMZ=GAMZ0
      ELSEIF (JPR.EQ.20) THEN
         IF(JPR0.LE.2039)THEN
           WRITE(*,*)'Enter top mass, W mass'
           READ(*,'(A)') TMPSTR2
           READ(TMPSTR2,*,ERR=618)RMASS(6),RMASS(198)
           RMASS(199)=RMASS(198)
         ELSEIF(JPR0.GT.2039.AND.JPR0.LE.2044)THEN
           WRITE(*,*)'Enter top mass, H+ mass, W mass'
           READ(*,'(A)') TMPSTR2
           READ(TMPSTR2,*,ERR=618)RMASS(6),RMASS(206),RMASS(198)
           RMASS(207)=RMASS(206)
c Default lifetime for charged Higgs will prevent its decay
           RLTIM(206)=0D0
           RLTIM(207)=0D0
         ELSE
           WRITE(*,*)'Unknow single-top channel: IPROC =',IPROC
           CALL HWWARN('HWIGPR',501)
         ENDIF
         GOTO 619
 618     OLDFORM=.TRUE.
         READ(TMPSTR2,*) RMASS(6)
 619     CONTINUE
      ELSE
         WRITE(*,*)'   Bad process code IPROC =',IPROC
         CALL HWWARN('HWIGPR',502)
      ENDIF
      WRITE(*,*)'Enter quark (d,u,s,c,b) and gluon masses'
      READ(*,*)RMASS(1),RMASS(2),RMASS(3),
     #         RMASS(4),RMASS(5),RMASS(13)
      DO I=1,5
         RMASS(I+6)=RMASS(I)
      ENDDO
C NO SOFT ME CORRECTIONS (ALREADY INCLUDED IN MC@NLO)
      IF(IPROC.LT.0)THEN
        SOFTME=.FALSE.
      ELSE
        SOFTME=.TRUE.
      ENDIF
C HERWIG RUNS IN WEIGHTED MODE EXCEPT IN THE CASE OF SINGLE TOP IN
C ORDER TO AVOID TOO MANY EVENTS KILLED BY SHOWER
      NOWGT=.FALSE.
      IF(JPR.EQ.20.AND.IPROC.GT.0)NOWGT=.TRUE.
C NEGATIVE WEIGHTS ALLOWED
      NEGWTS=.TRUE.
      MAXPR=0
      MAXER=MAXEV/100
      LRSUD=0
      LWSUD=77
C IN THE CASE HERWIG PDFS ARE USED, ADOPT MRST
      NSTRU=8
      PRVTX=.FALSE.
      PTMIN=0.5
      NRN(1)=1973774260
      NRN(2)=1099242306
C THE FOLLOWING SHOULD BE USED ONLY IN WEIGHTED MODE
      IF(.NOT.NOWGT)THEN
        WGTMAX=1.000001D0
        AVABW=1.000001D0
      ENDIF
C FOR TOP PRODUCTION (HARMLESS ELSEWHERE)
      RLTIM(6)=1.D-23
      RLTIM(12)=1.D-23
C---B FRAGMENTATION PARAMETERS (FOR B PRODUCTION ONLY)
      IF(ABS(IPROC).EQ.1705.OR.ABS(IPROC).EQ.11705)THEN
        PSPLT(2)=0.5
      ENDIF
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
C---USE THE FOLLOWING FOR SINGLE TOP -- AVOIDS TROUBLES WITH ISR
      IF(JPR.EQ.20)THEN
        CALL HWUSTA('B+      ')
        CALL HWUSTA('B-      ')
        CALL HWUSTA('B_D0    ')
        CALL HWUSTA('B_DBAR0 ')
        CALL HWUSTA('B_S0    ')
        CALL HWUSTA('B_SBAR0 ')
        CALL HWUSTA('SIGMA_B+')
        CALL HWUSTA('LMBDA_B0')
        CALL HWUSTA('SIGMA_B-')
        CALL HWUSTA('XI_B0   ')
        CALL HWUSTA('XI_B-   ')
        CALL HWUSTA('OMEGA_B-')
        CALL HWUSTA('B_C-    ')
        CALL HWUSTA('UPSLON1S')
        CALL HWUSTA('SGM_BBR-')
        CALL HWUSTA('LMD_BBR0')
        CALL HWUSTA('SGM_BBR+')
        CALL HWUSTA('XI_BBAR0')
        CALL HWUSTA('XI_B+   ')
        CALL HWUSTA('OMG_BBR+')
        CALL HWUSTA('B_C+    ')
        IF(JPR0.GT.2039.AND.JPR0.LE.2044)THEN
C---INSERT CHARGED HIGGS DECAYS HERE
          CALL HWMODK(37,ONE,0,-13,14,0,0,0)
          CALL HWMODK(-37,ONE,0,13,-14,0,0,0)
        ENDIF
      ENDIF
C---USER'S INITIAL CALCULATIONS
      CALL HWABEG
C---INITIALISE ELEMENTARY PROCESS
      CALL HWEINI
C---LOOP OVER EVENTS
      DO 100 N=1,MAXEV
C---INITIALISE EVENT
         CALL HWUINE
C---GENERATE HARD SUBPROCESS
         CALL HWEPRO
C---GENERATE PARTON CASCADES
         CALL HWBGEN
C---DO HEAVY OBJECT DECAYS
         CALL HWDHOB
C---DO CLUSTER FORMATION
         CALL HWCFOR
C---DO CLUSTER DECAYS
         CALL HWCDEC
C---DO UNSTABLE PARTICLE DECAYS
         CALL HWDHAD
C---DO HEAVY FLAVOUR HADRON DECAYS
         CALL HWDHVY
C---ADD SOFT UNDERLYING EVENT IF NEEDED
         CALL HWMEVT
C---FINISH EVENT
         CALL HWUFNE
C---USER'S EVENT ANALYSIS
         CALL HWANAL
         IF(MOD(NEVHEP,NSTEP).EQ.0) THEN
            WRITE(*,*)'# of events processed=',NEVHEP
            CALL HWAEND
         ENDIF
  100 CONTINUE
C---TERMINATE ELEMENTARY PROCESS
      CALL HWEFIN
C---USER'S TERMINAL CALCULATIONS
      WRITE(*,*)'# of events processed=',NEVHEP
      CALL HWAEND
C---CLEAN EXIT IF USING ROOT; DUMMY OTHERWISE
      CALL RCLOS()
 999  STOP
      END
