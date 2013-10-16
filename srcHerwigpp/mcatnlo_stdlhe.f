      subroutine write_lhef_init(ifile,ivvcode,xckecm,xsec,xsecerr,
     #                           part1,part2)
c Writes compulsory initialization information for LHEF on file #ifile;
c this corresponds to the HEPRUP common block (see hep-ph/0609017).
c This code is derived where possible from mcatnlo_hwlhin.f
      implicit none
      integer ifile,ivvcode,i
      real*8 xckecm,xsec,xsecerr,zero
      character * 4 part1,part2
      parameter (zero=0.d0)
      real*8 wgtmax
      common/cwgtmax/wgtmax
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
c
      call getipdgpart(part1,IDBMUP(1))
      call getipdgpart(part2,IDBMUP(2))
      DO I=1,2
         EBMUP(I)=XCKECM/2.d0
c The current setup ignores PDFGUP(*) -- the proper PDFs are determined
c by the scripts. If need be, set PDFGUP(*)=0, PDFSUP(*)=LHAGLUE_NUMBER
c in the case external PDF are used
         PDFGUP(I)=-1
         PDFSUP(I)=-1
      ENDDO
      IDWTUP=-4
      NPRUP=1
      write(ifile,'(a)')
     # '  <init>'
      write(ifile,501)IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     #                PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),
     #                IDWTUP,NPRUP
      XSECUP(1)=xsec
      XERRUP(1)=xsecerr
      if(XERRUP(1).eq.0.d0)XERRUP(1)=XSECUP(1)*5.d-4
      XMAXUP(1)=wgtmax
      IVVCODE=MOD(ABS(IVVCODE),10000)
      LPRUP(1)=IVVCODE
      write(ifile,502)XSECUP(1),XERRUP(1),XMAXUP(1),LPRUP(1)
      write(ifile,'(a)')
     # '  </init>'
c This are event-specific parameters, to be used in store_events. Some
c are set here as done in UPINIT in mcatnlo_hwlhin.f
      NUP=10
      AQEDUP=ZERO
      AQCDUP=ZERO
      IDPRUP=LPRUP(1)
      DO I=1,NUP
         VTIMUP(I)=ZERO
         SPINUP(I)=9.
      ENDDO
 501  format(2(1x,i6),2(1x,e14.8),2(1x,i2),2(1x,i6),1x,i2,1x,i3)
 502  format(3(1x,e14.8),1x,i6)
c
      return
      end


      subroutine write_lhef_event(ifile,idec,
     #     i1hpro,ic,np,ips,
     #     xevwgt,xscale,
     #     xmom_lb)
c Writes event information in LHEF format on file #ifile;
c this corresponds to the HEPEUP common block (see hep-ph/0609017).
c This code is derived where possible from mcatnlo_hwlhin.f, in particular
c from the routine UPEVNT() [of version 4.0X] (up to processes that will
c never be implemented in MC@NLO/HW++)
      implicit none
      integer ifile,idec,i1hpro,ic,np,ips(9)
      real*8 zero,xevwgt,xscale,xmom_lb(9,4)
      integer maxevt
      common/cmaxevt/maxevt
      parameter (zero=0.d0)
C--Les Houches Common Blocks
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &     IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &     XMAXUP(MAXPUP),LPRUP(MAXPUP)
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
      DOUBLE PRECISION UX1,UX2,UQ2
      COMMON/UXX/UX1,UX2,UQ2
      DOUBLE PRECISION WGTACP(28)
      COMMON/CWGTACP/WGTACP
      DOUBLE PRECISION E,PCM(5),PTR,XMP2,XMA2,XMB2,S2,SIGMA,DELTA,
     #     BETA,VA,VB,XKA,XKB,PL,PTF,XMTR,PDB(5),HWULDO123
      INTEGER I,J,MQQ,ISCALE,JJPROC,JPR,IHPRO,IA,IB,IBOS,ILEP,IPROC

      INTEGER ICOL4(4,5),ICOL5(5,18),ICOST(4,20),IP,IG,ID,MUP

      LOGICAL BOPRO,NODEC,REMIT
c     

      integer iwgtnorm
      common/ciwgtnorm/iwgtnorm

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
C---Colour flows for single top production
      DATA ICOST/
     & 20,12,00,10, 12,20,00,10, 10,02,02,10, 02,10,02,10,
     & 20,10,20,10, 20,02,01,10, 10,20,20,10, 02,20,01,10,
     & 20,02,01,10, 10,02,02,10, 02,20,01,10, 02,10,02,10,
     & 20,10,20,10, 10,20,20,10, 12,23,03,10, 23,12,03,10,
     & 30,12,32,10, 30,23,21,10, 12,30,32,10, 23,30,21,10/

c     
      MQQ=maxevt
      DO I=1,NP
         IDUP(I)=ips(I)
         DO J=1,4
            PUP(J,I)=xmom_lb(I,J)
         ENDDO
      ENDDO
C---  Les Houches expects mean weight to be the cross section in pb
      XWGTUP=xevwgt
      if(iwgtnorm.gt.0) XWGTUP= XWGTUP*MQQ
      
C---  Input format is now (px,py,pz,m)
      DO I=1,NP
         E=PUP(1,I)**2+PUP(2,I)**2+PUP(3,I)**2+PUP(4,I)**2
         E=SQRT(E)
         PUP(5,I)=PUP(4,I)
         PUP(4,I)=E
      ENDDO
c     Define CM four-momentum
      do j=1,4
         PCM(J)=PUP(J,1)+PUP(J,2)
      enddo
      PCM(5)=PCM(4)**2-PCM(1)**2-PCM(2)**2-PCM(3)**2
      if(pcm(5).le.0.d0)then
         write(*,*)'Error #2 in write_lhef_event:',pcm(5)
         stop
      endif
      PCM(5)=SQRT(PCM(5))
C---REMIT MEANS A REAL PARTON EMISSION OCCURRED
      REMIT=PUP(4,3).NE.ZERO
C---NODEC MEANS DECAYS NOT YET DONE
      NODEC=NP.EQ.5
      NUP=NP
C---CHECK PROCESS CODE
      IPROC=LPRUP(1)
      JJPROC=MOD(ABS(LPRUP(1)),10000)
      JPR=JJPROC/100
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
C         IF(ABS(IPROC).EQ.1705.OR.ABS(IPROC).EQ.11705)ISCALE=5
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
      ELSE         
         write(*,*)'Error #999 in write_lhef_event: Process',
     .        IPROC,'not supported'
         stop
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
            S2=XMA2+XMB2+2*HWULDO123(PUP(1,IA),PUP(1,IB))
            SIGMA=XMA2+XMB2
            DELTA=XMA2-XMB2
            BETA=SQRT(1-2*SIGMA/S2+(DELTA/S2)**2)
            VA=BETA/(1+DELTA/S2)
            VB=BETA/(1-DELTA/S2)
            XKA=HWULDO123(PUP(1,3),PUP(1,IA))
            XKB=HWULDO123(PUP(1,3),PUP(1,IB))
            E=(XKA+XKB)/SQRT(S2)
            PL=-2.0/((VA+VB)*BETA*SQRT(S2))*(VA*XKA-VB*XKB)
            PTF=E**2-PL**2-XMP2
            IF (PTF.LE.ZERO) THEN
               write(*,*)'Error #103 in write_lhef_event:',ptf
               PTF=0.d0
CCC TODO: SURE TO CONTINUE? NOT STOP?     
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
               write(*,*)'Error #100 in write_lhef_event:',SCALUP
               SCALUP=0.d0
CCC TODO: SURE TO CONTINUE? NOT STOP?              
            ENDIF
         ELSEIF (ISCALE.EQ.6) THEN
            XMTR=SQRT(PUP(5,4)**2+PUP(1,4)**2+PUP(2,4)**2)
            PTR=SQRT(PUP(1,3)**2+PUP(2,3)**2)
            SCALUP=PCM(5)-PTR-XMTR
            IF (SCALUP.LE.ZERO) THEN
               write(*,*)'Error #100 in write_lhef_event:',SCALUP
               SCALUP=0.d0
CCC TODO: SURE TO CONTINUE? NOT STOP?
            ENDIF
         ELSEIF (ISCALE.EQ.7) THEN
            SCALUP=SQRT(PUP(5,4)**2+PUP(1,4)**2+PUP(2,4)**2)
         ELSE
            write(*,*)'Error #501 in write_lhef_event:',ISCALE
            stop
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
               write(*,*)'Error #103 in write_lhef_event:',IHPRO
               stop
            ENDIF
         ELSE
            CALL HWVEQU123(5,PUP(1,4),PUP(1,3))
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
         IF (JPR.EQ.16) THEN
            IDUP(NUP)=25
         ELSE
           write(*,*)'Error #4 in write_lhef_event:',JPR
           stop
         ENDIF
      ELSEIF (JPR.EQ.51.AND.IDUP(1).EQ.22) THEN
C---HEAVY QUARK DIRECT PHOTOPRODUCTION
         write(*,*) 'Error 997 in write_lhef_event: Process',
     .        IPROC,'not supported.'
         stop
      ELSEIF (JPR.EQ.19) THEN
C---HIGGS VIA VBF
         write(*,*) 'Error 996 in write_lhef_event: Process',
     .        IPROC,'not supported.'
         stop
      ELSEIF (JPR.EQ.17.OR.(JPR.EQ.51.AND.IDUP(1).NE.22)) THEN
C--- DOUBLE X CHECK
         IF(LPRUP(1).NE.1705.AND.LPRUP(1).NE.1706)THEN
            write(*,*) 'Error 995 in write_lhef_event: Process',
     .           IPROC,'not supported.'
            stop
         ENDIF         
C---  HEAVY QUARKS
         IF (REMIT) THEN
C---3-BODY FINAL STATE
C---SET COLOUR CONNECTIONS
            IF (IC.LE.18) THEN
               DO I=1,5
                  CALL UPCODE(ICOL5(I,IC),ICOLUP(1,I))
               ENDDO
            ELSE
               write(*,*)'Error #105 in write_lhef_event:',IC
               stop
            ENDIF
         ELSE
C---2-BODY FINAL STATE
            DO IP=3,NUP
               IDUP(IP)=IDUP(IP+1)
               CALL HWVEQU123(5,PUP(1,IP+1),PUP(1,IP))
            ENDDO
C---SET COLOUR CONNECTIONS
            IF (IC.LE.4) THEN
               DO I=1,4
                  CALL UPCODE(ICOL4(I,IC),ICOLUP(1,I))
               ENDDO
            ELSE
               write(*,*)'Error #104 in write_lhef_event:',IC
               stop
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
            IF (LPRUP(1).NE.1706) THEN
               write(*,*)'Error #210 in write_lhef_event:',LPRUP(1)
               stop
            ENDIF
            NP=NUP
C--W DECAYS
            IDUP(NP+1)=IDUP(NP-5)
            IDUP(NP+2)=IDUP(NP-4)
            IDUP(NP+3)=IDUP(NP-2)
            IDUP(NP+4)=IDUP(NP-1)
            IDUP(NP-1)=IDUP(NP)
            CALL HWVEQU123(5,PUP(1,NP-5),PUP(1,NP+1))
            CALL HWVEQU123(5,PUP(1,NP-4),PUP(1,NP+2))
            CALL HWVEQU123(5,PUP(1,NP-2),PUP(1,NP+3))
            CALL HWVEQU123(5,PUP(1,NP-1),PUP(1,NP+4))
            CALL HWVEQU123(5,PUP(1,NP  ),PUP(1,NP-1))
            CALL HWVSUM123(4,PUP(1,NP+1),PUP(1,NP+2),PUP(1,NP-2))
            CALL HWVSUM123(4,PUP(1,NP+3),PUP(1,NP+4),PUP(1,NP))
            CALL HWUMAS123(PUP(1,NP-2))
            CALL HWUMAS123(PUP(1,NP))
            IDUP(NP-2)=-24*(MOD(IDUP(NP+1),2)+MOD(IDUP(NP+2),2))
            IDUP(NP)=-IDUP(NP-2)
            DO IP=NP-3,NP+4
               ISTUP(IP)=1
            ENDDO
            ISTUP(NP-2)=2
            ISTUP(NP)=2
C--TOP DECAYS
            CALL HWVSUM123(4,PUP(1,NP-3),PUP(1,NP-2),PUP(1,NP-5))
            CALL HWVSUM123(4,PUP(1,NP-1),PUP(1,NP  ),PUP(1,NP-4))
            CALL HWUMAS123(PUP(1,NP-5))
            CALL HWUMAS123(PUP(1,NP-4))
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
               write(*,*)'Error #111 in write_lhef_event:',IC
               stop
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
                  CALL HWVEQU123(5,PUP(1,IP+1),PUP(1,IP))
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
               CALL HWVEQU123(5,PUP(1,NP-1),PUP(1,NP+2))
               CALL HWVEQU123(5,PUP(1,NP),PUP(1,NP+3))
               CALL HWVSUM123(4,PUP(1,NP-1),PUP(1,NP),PUP(1,NP+1))
               CALL HWUMAS123(PUP(1,NP+1))
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
               write(*,*)'Error #108 in write_lhef_event:',IA
               stop
            ELSEIF (IB.LT.1.OR.IB.GT.5) THEN
               write(*,*)'Error #109 in write_lhef_event:',IB
               stop
            ELSEIF (IA.EQ.IB) THEN
               write(*,*)'Error #110 in write_lhef_event:',IA,IB
               stop
            ENDIF
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
                  CALL HWVEQU123(5,PUP(1,IP+1),PUP(1,IP))
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
            CALL HWVEQU123(5,PUP(1,NP-2),PUP(1,NP+1))
            CALL HWVEQU123(5,PUP(1,NP-1),PUP(1,NP+2))
            CALL HWVSUM123(4,PUP(1,NP+1),PUP(1,NP+2),PUP(1,NP-1))
            CALL HWVSUM123(4,PUP(1,NP-1),PUP(1,NP  ),PUP(1,NP-2))
            CALL HWUMAS123(PUP(1,NP-1))
            CALL HWUMAS123(PUP(1,NP-2))
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
         write(*,*)'Error #5 in write_lhef_event: Process',
     .        IPROC,'not implemented.'
         stop
      ELSE
C---BOSON PAIR OR LEPTON PAIR
         IF (BOPRO.OR.NODEC) THEN
            NUP=6
            DO I=6,5,-1
               CALL HWVEQU123(5,PUP(1,I-1),PUP(1,I))
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
               CALL HWVEQU123(5,PUP(1,MUP-I),PUP(1,MUP))
               CALL HWVEQU123(5,PUP(1,NUP-I),PUP(1,NUP))
               CALL HWVSUM123(4,PUP(1,MUP),PUP(1,NUP),PUP(1,6))
               CALL HWUMAS123(PUP(1,6))
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
                  write(*,*) 'Error #106 in write_lhef_event'
                  stop
               ENDIF
            ENDIF
            IF (ABS(IDUP(4)).LT.20) THEN
               CALL HWVZRO123(4,PDB)
               DO I=8,7,-1
                  CALL HWVEQU123(5,PUP(1,I-3),PUP(1,I))
                  CALL HWVSUM123(4,PUP(1,I),PDB,PDB)
                  IDUP(I)=IDUP(I-3)
                  ISTUP(I)=1
                  MOTHUP(1,I)=5
                  MOTHUP(2,I)=5
               ENDDO
               CALL HWUMAS123(PDB)
               CALL HWVEQU123(5,PDB,PUP(1,5))
               ISTUP(5)=2
               ID=MOD(IDUP(7),2)+MOD(IDUP(8),2)
               IF (ID.EQ.0) THEN
                  IDUP(5)=23
               ELSEIF (ABS(ID).EQ.1) THEN
                  IDUP(5)=-24*ID
               ELSE
                  write(*,*) 'Error #107 in write_lhef_event'
                  stop
               ENDIF
            ELSE
               CALL HWVEQU123(5,PUP(1,4),PUP(1,5))
               IDUP(5)=IDUP(4)
               ISTUP(5)=1
               MOTHUP(1,5)=4
               MOTHUP(2,5)=4
            ENDIF
         ENDIF
C---ADD DIBOSON OR DILEPTON TO EVENT RECORD (TO FIX ITS MASS)
         CALL HWVZRO123(4,PDB)
         DO I=6,5,-1
            CALL HWVSUM123(4,PUP(1,I),PDB,PDB)
            MOTHUP(1,I)=4
            MOTHUP(2,I)=4
         ENDDO
         CALL HWUMAS123(PDB)
         CALL HWVEQU123(5,PDB,PUP(1,4))
         ISTUP(4)=2
         IDUP(4)=0
         IF (REMIT) THEN
            IF(.NOT.BOPRO.AND..TRUE.) THEN
C--- REMOVE SPURIOUS 0 PARTICLE FROM EVENT LISTING, HW++ DOES NO LIKE THEM
               DO I=5,NUP
                  CALL HWVEQU123(5,PUP(1,I),PUP(1,I-1))
                  IDUP(I-1)=IDUP(I)
                  ISTUP(I-1)=ISTUP(I)
                  IF(MOTHUP(1,I).EQ.4.AND.MOTHUP(2,I).EQ.4) THEN
                     MOTHUP(1,I-1)=1
                     MOTHUP(2,I-1)=2
                  ELSEIF(MOTHUP(1,I).EQ.4.OR.MOTHUP(2,I).EQ.4) THEN
                     write(*,*) 'Error in removing 0 particle'
                     stop
                  ELSE
                     MOTHUP(1,I-1)=MOTHUP(1,I)-1
                     MOTHUP(2,I-1)=MOTHUP(1,I)-1
                  ENDIF
               ENDDO
               NUP=NUP-1
            ENDIF
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
               write(*,*)'Error #102 in write_lhef_event:',IHPRO
               stop
            ENDIF
            DO I=4,NUP
               ICOLUP(1,I)=0
               ICOLUP(2,I)=0
            ENDDO
         ELSE
            DO I=5,NUP
               CALL HWVEQU123(5,PUP(1,I),PUP(1,I-2))
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
               write(*,*)'Error #504 in write_lhef_event:',IBOS,JPR
               stop
            ENDIF
c Best to include the intermediate vector boson also in the case of 2->2 
c events. Do it here rather than modifying the code above
            IF(.NOT.REMIT) THEN
c Move leptons one place up in the event record
               NUP=5
               DO I=5,4,-1
                  CALL HWVEQU123(5,PUP(1,I-1),PUP(1,I))
                  IDUP(I)=IDUP(I-1)
                  ISTUP(I)=1
               ENDDO
c Create vector boson entry
               CALL HWVZRO123(4,PDB)
               DO I=5,4,-1
                  CALL HWVSUM123(4,PUP(1,I),PDB,PDB)
                  MOTHUP(1,I)=3
                  MOTHUP(2,I)=3
               ENDDO
               CALL HWUMAS123(PDB)
               CALL HWVEQU123(5,PDB,PUP(1,3))
               ISTUP(3)=2
               MOTHUP(1,3)=1
               MOTHUP(2,3)=2
               IF ( IBOS.EQ.50 .OR.
     #                 (IBOS.EQ.60.AND.JPR.EQ.13) .OR.
     #                 (IBOS.EQ.70.AND.JPR.EQ.13) ) THEN
                  IDUP(3)=23
               ELSEIF (IBOS.EQ.60.AND.JPR.EQ.14) THEN
                  IDUP(3)=24
               ELSEIF (IBOS.EQ.70.AND.JPR.EQ.14) THEN
                  IDUP(3)=-24
               ELSE
                  write(*,*)'Error #506 in write_lhef_event:',IBOS,JPR
                  stop
               ENDIF               
            ENDIF
         ENDIF
      ENDIF
C---INITIAL STATE
      DO I=1,2
         ISTUP(I)=-1
         MOTHUP(1,I)=0
         MOTHUP(2,I)=0
      ENDDO
      do i=1,nup
         id=abs(idup(i))
         if (id.gt.10.and.id.lt.17) then
            icolup(1,i)=0
            icolup(2,i)=0
         endif
      enddo
c Write all on event file
      write(ifile,'(a)')
     # '  <event>'
      write(ifile,503)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      do i=1,nup
        write(ifile,504)IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),
     #                  ICOLUP(1,I),ICOLUP(2,I),
     #                  PUP(1,I),PUP(2,I),PUP(3,I),PUP(4,I),PUP(5,I),
     #                  VTIMUP(I),SPINUP(I)
      enddo
      write(ifile,505) '#pdf',ux1,ux2,uq2
      if(jjproc.eq.2850)then
        write(ifile,507) '#an_cpl_wgt',(wgtacp(j),j=1,28)
      elseif(jjproc.eq.2870.or.jjproc.eq.2880)then
        write(ifile,506) '#an_cpl_wgt',
     #                    wgtacp(1), wgtacp(2), wgtacp(3),
     #                    wgtacp(4), wgtacp(8), wgtacp(9),
     #                    wgtacp(10),wgtacp(14),wgtacp(15),
     #                    wgtacp(19)
      endif
      write(ifile,'(a)')
     # '  </event>'
 503  format(1x,i2,1x,i6,4(1x,e14.8))
 504  format(1x,i6,1x,i2,4(1x,i4),5(1x,e14.8),2(1x,e10.4))
 505  format(1x,a4,3(1x,e14.8))
 506  format(1x,a11,10(1x,e14.8))
 507  format(1x,a11,28(1x,e14.8))
c
      return
      end


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
      FUNCTION HWULDO123(P,Q)
C----------------------------------------------------------------------
C   LORENTZ 4-VECTOR DOT PRODUCT
C SF, 2/6/08: renamed to avoid conflicts
C----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION HWULDO123,P(4),Q(4)
      HWULDO123=P(4)*Q(4)-(P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3))
      END


C-----------------------------------------------------------------------
      SUBROUTINE HWVEQU123(N,P,Q)
C-----------------------------------------------------------------------
C     VECTOR EQUALITY
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION P(N),Q(N)
      DO 10 I=1,N
   10 Q(I)=P(I)
      END


C-----------------------------------------------------------------------
      SUBROUTINE HWVSUM123(N,P,Q,R)
C-----------------------------------------------------------------------
C    VECTOR SUM
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION P(N),Q(N),R(N)
      DO 10 I=1,N
   10 R(I)=P(I)+Q(I)
      END


C-----------------------------------------------------------------------
      SUBROUTINE HWUMAS123(P)
C-----------------------------------------------------------------------
C     PUTS INVARIANT MASS IN 5TH COMPONENT OF VECTOR
C     (NEGATIVE SIGN IF SPACELIKE)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION HWUSQR123,P(5)
      EXTERNAL HWUSQR123
      P(5)=HWUSQR123((P(4)+P(3))*(P(4)-P(3))-P(1)**2-P(2)**2)
      END


C-----------------------------------------------------------------------
      FUNCTION HWUSQR123(X)
C-----------------------------------------------------------------------
C     SQUARE ROOT WITH SIGN RETENTION
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION HWUSQR123,X
      HWUSQR123=SIGN(SQRT(ABS(X)),X)
      END


C-----------------------------------------------------------------------
      SUBROUTINE HWVZRO123(N,P)
C-----------------------------------------------------------------------
C     ZERO VECTOR
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION P(N)
      DO 10 I=1,N
   10 P(I)=0D0
      END
