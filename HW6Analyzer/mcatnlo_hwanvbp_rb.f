C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
c See the documentation for the labels: (V1,V2) are (W+,W-), (Z,Z),
c (W+,Z), and (W-,Z) for IPROC=-2850, -2860, -2870, and 2880 respectively.
c In the case of IPROC=-2850, l1 and l2 are the charged leptons resulting
c from the decay of W+ (==V1) and of W- (==V2)
      INCLUDE 'HERWIG65.INC'
      real * 8 xm01,xm02,bin1,bin2,ga1,ga2,xm1low,xm1upp,xm2low,xm2upp
      real * 4 pi,xm1i,xm1s,xm2i,xm2s
      parameter (pi=3.14160E0)
      integer jpr,j,k,ivlep1,ivlep2,idec
      common/vvlin/ivlep1,ivlep2
      character*5 cc(2)
      data cc/'     ',' cuts'/
c
      if(ivlep1.gt.6)then
        if(ivlep2.le.6)call hwwarn('HWABEG',501)
        idec=1
      else
        if(ivlep2.gt.6)call hwwarn('HWABEG',501)
        idec=0
      endif
c
      jpr=mod(abs(iproc),10000)/10
      if(jpr.eq.285)then
        xm01=rmass(198)
        xm02=rmass(199)
        ga1=gamw
        ga2=gamw
      elseif(jpr.eq.286)then
        xm01=rmass(200)
        xm02=rmass(200)
        ga1=gamz
        ga2=gamz
      elseif(jpr.eq.287)then
        xm01=rmass(198)
        xm02=rmass(200)
        ga1=gamw
        ga2=gamz
      elseif(jpr.eq.288)then
        xm01=rmass(199)
        xm02=rmass(200)
        ga1=gamw
        ga2=gamz
      else
        call hwwarn('HWABEG',500)
      endif
c
      if(idec.eq.0)then
        xm1low=max(0.d0,xm01-gammax*ga1)
        xm1upp=xm01+gammax*ga1
        bin1=(xm1upp-xm1low)/100.d0
        xm1i=sngl(xm01-(49*bin1+bin1/2.d0))
        xm1s=sngl(xm01+(50*bin1+bin1/2.d0))
        xm2low=max(0.d0,xm02-gammax*ga2)
        xm2upp=xm02+gammax*ga2
        bin2=(xm2upp-xm2low)/100.d0
        xm2i=sngl(xm02-(49*bin2+bin2/2.d0))
        xm2s=sngl(xm02+(50*bin2+bin2/2.d0))
      else
        bin1=0.5d0
        xm1i=sngl(xm01-24.75d0)
        xm1s=sngl(xm01+25.25d0)
        bin2=0.5d0
        xm2i=sngl(xm02-24.75d0)
        xm2s=sngl(xm02+25.25d0)
      endif
c
      call rinit("HERVB.root")
      do j=1,2
      k=(j-1)*50
      call rbook(k+ 1,'V1 pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 2,'V1 log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+ 3,'V2 pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 4,'V2 log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+ 5,'y_V1'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 6,'y_V2'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 7,'y_V1-y_V2'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 8,'eta_V1-eta_V2'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 9,'VV inv m'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+10,'y_VV'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+11,'VV pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+12,'VV log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+13,'VV azimt'//cc(j),pi/20.e0,0.e0,pi)
      call rbook(k+14,'VV log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
      call rbook(k+15,'M_V1'//cc(j),sngl(bin1),xm1i,xm1s)
      call rbook(k+16,'M_V2'//cc(j),sngl(bin2),xm2i,xm2s)
c
      if(idec.eq.0)then
      k=(j-1)*50+16
      call rbook(k+ 1,'l1 pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 2,'l1 log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+ 3,'l2 pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 4,'l2 log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+ 5,'y_l1'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 6,'y_l2'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 7,'y_l1-y_l2'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+ 8,'ll inv m'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+ 9,'y_ll'//cc(j),0.2e0,-4.e0,4.e0)
      call rbook(k+10,'ll pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+11,'ll log[pt]'//cc(j),0.05e0,-0.5e0,3.e0)
      call rbook(k+12,'ll azimt'//cc(j),pi/20.e0,0.e0,pi)
      endif
c
      enddo
c
 999  END


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
      INTEGER ICHSUM,ICHINI,IHEP,IST,ID,ID1,IJ,JJV1,JJV2,IV1,IV2,ILL1,
     # ILLB1,ILL2,ILLB2,JPR,IDENT1,IDENT2
      LOGICAL DIDSOF
      REAL * 8 PPV1(4),PPV2(4),PPE1(4),PPNU1(4),PPE2(4),PPNU2(4),
     # PPVV(4),PPEE(4)
      REAL * 8 PTV1,YV1,GETRAPIDITY,ETAV1,GETPSEUDORAP,PTV2,YV2,ETAV2,
     # XMVV,GETINVM,YVV,PTE1,YE1,ETAE1,PTE2,YE2,ETAE2,PTEE,XMEE,DPHIEE,
     # GETDELPHI,YEE,PTVV,DPHIVV,XLDPHIVV,XMV1,XMV2
      INTEGER I,KK,IVLEP1,IVLEP2,IDEC
      COMMON/VVLIN/IVLEP1,IVLEP2
      REAL*8 PI
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 WWW0
c
      IF (IERROR.NE.0) RETURN
c
      WWW0=EVWGT
      CALL HWVSUM(4,PHEP(1,1),PHEP(1,2),PSUM)
      CALL HWVSCA(4,-1D0,PSUM,PSUM)
      ICHSUM=0
      ICHINI=ICHRG(IDHW(1))+ICHRG(IDHW(2))
      DIDSOF=.FALSE.
      JJV1=0
      JJV2=0
      JPR=MOD(ABS(IPROC),10000)/10
      IF(JPR.EQ.285)THEN
        IDENT1=24
        IDENT2=-24
      ELSEIF(JPR.EQ.286)THEN
        IDENT1=23
        IDENT2=23
      ELSEIF(JPR.EQ.287)THEN
        IDENT1=24
        IDENT2=23
      ELSEIF(JPR.EQ.288)THEN
        IDENT1=-24
        IDENT2=23
      ELSE
        CALL HWWARN('HWANAL',500)
      ENDIF
      IF(IVLEP1.GT.6)THEN
        IDEC=1
      ELSE
        IDEC=0
      ENDIF
C
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
        IST=ISTHEP(IHEP)      
        ID=IDHW(IHEP)
        ID1=IDHEP(IHEP)
C LOOK FOR V1 AND V2
        IF(JPR.NE.286)THEN
C W+W-, W+Z, AND W-Z
          IF( (IST.EQ.155.AND.IDEC.EQ.0) .OR.
     #        (IST.EQ.195.AND.IDEC.EQ.1) )THEN
            IF(ID1.EQ.IDENT1)THEN
C FOUND THE V1
              JJV1=JJV1+1
              IV1=IHEP
              DO IJ=1,4
	        PPV1(IJ)=PHEP(IJ,IHEP)
	      ENDDO
            ELSEIF(ID1.EQ.IDENT2)THEN
C FOUND THE V2
              JJV2=JJV2+1
              IV2=IHEP
              DO IJ=1,4
	        PPV2(IJ)=PHEP(IJ,IHEP)
	      ENDDO
            ENDIF
          ENDIF
        ELSE
C ZZ
          IF( (IST.EQ.155.AND.IDEC.EQ.0) .OR.
     #        (IST.EQ.195.AND.IDEC.EQ.1) )THEN
            IF(ID1.EQ.IDENT1)THEN
C FOUND A Z
              IF(JJV1.EQ.0)THEN
                JJV1=JJV1+1
                DO IJ=1,5
	          PPV1(IJ)=PHEP(IJ,IHEP)
	        ENDDO
              ELSE
                JJV2=JJV2+1
                DO IJ=1,5
	          PPV2(IJ)=PHEP(IJ,IHEP)
	        ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF
  100 CONTINUE
      IF( (JJV1+JJV2).LT.2 .AND. IERROR.EQ.0)
     #  CALL HWWARN('HWANAL',501)
      DO IJ=1,4
        PPVV(IJ)=PPV1(IJ)+PPV2(IJ)
      ENDDO
C FIND THE LEPTONS
      IF(IDEC.EQ.0)THEN
        IF( IDHEP(JDAHEP(1,JDAHEP(1,IV1))).GT.0 .AND.
     #      IDHEP(JDAHEP(1,JDAHEP(2,IV1))).LT.0 )THEN
          ILL1=JDAHEP(1,JDAHEP(1,IV1))
          ILLB1=JDAHEP(1,JDAHEP(2,IV1))
        ELSEIF( IDHEP(JDAHEP(1,JDAHEP(2,IV1))).GT.0 .AND.
     #          IDHEP(JDAHEP(1,JDAHEP(1,IV1))).LT.0 )THEN
          ILL1=JDAHEP(1,JDAHEP(2,IV1))
          ILLB1=JDAHEP(1,JDAHEP(1,IV1))
        ELSE
          CALL HWUEPR
          CALL HWWARN('HWANAL',504)
        ENDIF
        IF( IDHEP(JDAHEP(1,JDAHEP(1,IV2))).GT.0 .AND.
     #      IDHEP(JDAHEP(1,JDAHEP(2,IV2))).LT.0 )THEN
          ILL2=JDAHEP(1,JDAHEP(1,IV2))
          ILLB2=JDAHEP(1,JDAHEP(2,IV2))
        ELSEIF( IDHEP(JDAHEP(1,JDAHEP(2,IV2))).GT.0 .AND.
     #          IDHEP(JDAHEP(1,JDAHEP(1,IV2))).LT.0 )THEN
          ILL2=JDAHEP(1,JDAHEP(2,IV2))
          ILLB2=JDAHEP(1,JDAHEP(1,IV2))
        ELSE
          CALL HWUEPR
          CALL HWWARN('HWANAL',504)
        ENDIF
        IF(JPR.EQ.285.OR.JPR.EQ.286)THEN
C W+W- AND ZZ
          DO IJ=1,4
            PPNU1(IJ)=PHEP(IJ,ILL1)
            PPE1(IJ)=PHEP(IJ,ILLB1)
            PPE2(IJ)=PHEP(IJ,ILL2)
            PPNU2(IJ)=PHEP(IJ,ILLB2)
          ENDDO
        ELSEIF(JPR.EQ.287)THEN
C W+Z
          DO IJ=1,4
            PPNU1(IJ)=PHEP(IJ,ILL1)
            PPE1(IJ)=PHEP(IJ,ILLB1)
            PPE2(IJ)=PHEP(IJ,ILL2)
            PPNU2(IJ)=PHEP(IJ,ILLB2)
          ENDDO
        ELSEIF(JPR.EQ.288)THEN
C W-Z
          DO IJ=1,4
            PPNU1(IJ)=PHEP(IJ,ILLB1)
            PPE1(IJ)=PHEP(IJ,ILL1)
            PPE2(IJ)=PHEP(IJ,ILLB2)
            PPNU2(IJ)=PHEP(IJ,ILL2)
          ENDDO
        ENDIF
C CHECK THAT THE LEPTONS ARE FINAL-STATE LEPTONS
        IF( ABS(IDHEP(ILL1)).LT.11 .OR. ABS(IDHEP(ILL1)).GT.16 .OR.
     #      ABS(IDHEP(ILLB1)).LT.11 .OR. ABS(IDHEP(ILLB1)).GT.16 .OR.
     #      ABS(IDHEP(ILL2)).LT.11 .OR. ABS(IDHEP(ILL2)).GT.16 .OR.
     #      ABS(IDHEP(ILLB2)).LT.11 .OR. ABS(IDHEP(ILLB2)).GT.16 )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',505)
        ENDIF
        IF( (ISTHEP(ILL1).NE.1.AND.ISTHEP(ILL1).NE.195) .OR.
     #      (ISTHEP(ILLB1).NE.1.AND.ISTHEP(ILLB1).NE.195) .OR.
     #      (ISTHEP(ILL2).NE.1.AND.ISTHEP(ILL2).NE.195) .OR.
     #      (ISTHEP(ILLB2).NE.1.AND.ISTHEP(ILLB2).NE.195) )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',506)
        ENDIF
C
        DO IJ=1,4
          PPEE(IJ)=PPE1(IJ)+PPE2(IJ)
        ENDDO
      ENDIF
C FILL THE HISTOS
      ptv1=sqrt(ppv1(1)**2+ppv1(2)**2)
      yv1=getrapidity(ppv1(4),ppv1(3))
      etav1=getpseudorap(ppv1(4),ppv1(1),ppv1(2),ppv1(3))
      xmv1=getinvm(ppv1(4),ppv1(1),ppv1(2),ppv1(3))
c
      ptv2=sqrt(ppv2(1)**2+ppv2(2)**2)
      yv2=getrapidity(ppv2(4),ppv2(3))
      etav2=getpseudorap(ppv2(4),ppv2(1),ppv2(2),ppv2(3))
      xmv2=getinvm(ppv2(4),ppv2(1),ppv2(2),ppv2(3))
c
      xmvv=getinvm(ppvv(4),ppvv(1),ppvv(2),ppvv(3))
      yvv=getrapidity(ppvv(4),ppvv(3))
      ptvv=sqrt(ppvv(1)**2+ppvv(2)**2)
      dphivv=getdelphi(ppv1(1),ppv1(2),ppv2(1),ppv2(2))
      xldphivv=(pi-dphivv)/pi
c
      if(idec.eq.0)then
        pte1=sqrt(ppe1(1)**2+ppe1(2)**2)
        ye1=getrapidity(ppe1(4),ppe1(3))
        etae1=getpseudorap(ppe1(4),ppe1(1),ppe1(2),ppe1(3))
c
        pte2=sqrt(ppe2(1)**2+ppe2(2)**2)
        ye2=getrapidity(ppe2(4),ppe2(3))
        etae2=getpseudorap(ppe2(4),ppe2(1),ppe2(2),ppe2(3))
c
        ptee=sqrt(ppee(1)**2+ppee(2)**2)
        xmee=getinvm(ppee(4),ppee(1),ppee(2),ppee(3))
        dphiee=getdelphi(ppe1(1),ppe1(2),ppe2(1),ppe2(2))
        yee=getrapidity(ppee(4),ppee(3))
      endif
c
      kk=0
      call rfill(kk+1,sngl(ptv1),sngl(WWW0))
      if(ptv1.gt.0.d0)call rfill(kk+2,sngl(log10(ptv1)),sngl(WWW0))
      call rfill(kk+3,sngl(ptv2),sngl(WWW0))
      if(ptv2.gt.0.d0)call rfill(kk+4,sngl(log10(ptv2)),sngl(WWW0))
      call rfill(kk+5,sngl(yv1),sngl(WWW0))
      call rfill(kk+6,sngl(yv2),sngl(WWW0))
      call rfill(kk+7,sngl(yv1-yv2),sngl(WWW0))
      call rfill(kk+8,sngl(etav1-etav2),sngl(WWW0))
      call rfill(kk+9,sngl(xmvv),sngl(WWW0))
      call rfill(kk+10,sngl(yvv),sngl(WWW0))
      call rfill(kk+11,sngl(ptvv),sngl(WWW0))
      if(ptvv.gt.0.d0)call rfill(kk+12,sngl(log10(ptvv)),sngl(WWW0))
      call rfill(kk+13,sngl(dphivv),sngl(WWW0))
      if(xldphivv.gt.0.d0)
     #  call rfill(kk+14,sngl(log10(xldphivv)),sngl(WWW0))
      call rfill(kk+15,sngl(xmv1),sngl(WWW0))
      call rfill(kk+16,sngl(xmv2),sngl(WWW0))
c
      if(idec.eq.0)then
        kk=16
        call rfill(kk+1,sngl(pte1),sngl(WWW0))
        if(pte1.gt.0.d0)call rfill(kk+2,sngl(log10(pte1)),sngl(WWW0))
        call rfill(kk+3,sngl(pte2),sngl(WWW0))
        if(pte2.gt.0.d0)call rfill(kk+4,sngl(log10(pte2)),sngl(WWW0))
        call rfill(kk+5,sngl(ye1),sngl(WWW0))
        call rfill(kk+6,sngl(ye2),sngl(WWW0))
        call rfill(kk+7,sngl(ye1-ye2),sngl(WWW0))
        call rfill(kk+8,sngl(xmee),sngl(WWW0))
        call rfill(kk+9,sngl(yee),sngl(WWW0))
        call rfill(kk+10,sngl(ptee),sngl(WWW0))
        if(ptee.gt.0.d0)call rfill(kk+11,sngl(log10(ptee)),sngl(WWW0))
        call rfill(kk+12,sngl(dphiee),sngl(WWW0))
      endif
c
      kk=50
      if(abs(yv1).lt.2.d0)then
        call rfill(kk+1,sngl(ptv1),sngl(WWW0))
        if(ptv1.gt.0.d0)call rfill(kk+2,sngl(log10(ptv1)),sngl(WWW0))
        if(ptv1.gt.20.d0)call rfill(kk+15,sngl(xmv1),sngl(WWW0))
      endif
      if(abs(yv2).lt.2.d0)then
        call rfill(kk+3,sngl(ptv2),sngl(WWW0))
        if(ptv2.gt.0.d0)call rfill(kk+4,sngl(log10(ptv2)),sngl(WWW0))
        if(ptv2.gt.20.d0)call rfill(kk+16,sngl(xmv2),sngl(WWW0))
      endif
      if(ptv1.gt.20.d0)call rfill(kk+5,sngl(yv1),sngl(WWW0))
      if(ptv2.gt.20.d0)call rfill(kk+6,sngl(yv2),sngl(WWW0))
      if(ptv1.gt.20.d0.and.ptv2.gt.20.d0.and.
     #   abs(yv1).lt.2.d0.and.abs(yv2).lt.2.d0)then
        call rfill(kk+7,sngl(yv1-yv2),sngl(WWW0))
        call rfill(kk+8,sngl(etav1-etav2),sngl(WWW0))
        call rfill(kk+9,sngl(xmvv),sngl(WWW0))
        call rfill(kk+10,sngl(yvv),sngl(WWW0))
        call rfill(kk+11,sngl(ptvv),sngl(WWW0))
        if(ptvv.gt.0.d0)call rfill(kk+12,sngl(log10(ptvv)),sngl(WWW0))
        call rfill(kk+13,sngl(dphivv),sngl(WWW0))
        if(xldphivv.gt.0.d0)
     #    call rfill(kk+14,sngl(log10(xldphivv)),sngl(WWW0))
      endif
c
      if(idec.eq.0)then
        kk=66
        if(abs(ye1).lt.2.d0)then
          call rfill(kk+1,sngl(pte1),sngl(WWW0))
          if(pte1.gt.0.d0)call rfill(kk+2,sngl(log10(pte1)),sngl(WWW0))
        endif
        if(abs(ye2).lt.2.d0)then
          call rfill(kk+3,sngl(pte2),sngl(WWW0))
          if(pte2.gt.0.d0)call rfill(kk+4,sngl(log10(pte2)),sngl(WWW0))
        endif
        if(pte1.gt.20.d0)call rfill(kk+5,sngl(ye1),sngl(WWW0))
        if(pte2.gt.20.d0)call rfill(kk+6,sngl(ye2),sngl(WWW0))
        if(pte1.gt.20.d0.and.pte2.gt.20.d0.and.
     #     abs(ye1).lt.2.d0.and.abs(ye2).lt.2.d0)then
          call rfill(kk+7,sngl(ye1-ye2),sngl(WWW0))
          call rfill(kk+8,sngl(xmee),sngl(WWW0))
          call rfill(kk+9,sngl(yee),sngl(WWW0))
          call rfill(kk+10,sngl(ptee),sngl(WWW0))
          if(ptee.gt.0.d0)call rfill(kk+11,sngl(log10(ptee)),sngl(WWW0))
          call rfill(kk+12,sngl(dphiee),sngl(WWW0))
        endif
      endif
c
C CHECK MOMENTUM CONSERVATION
      IF (HWVDOT(3,PSUM,PSUM).GT.1.E-4*PHEP(4,1)**2) THEN
         CALL HWWARN('HWANAL',1)
  101     PRINT 102, PSUM
  102     FORMAT(' Bad Momentum Sum = ',4F9.3/)
         CALL HWUEPR
         STOP
      ENDIF
      IF (ICHSUM.NE.ICHINI) THEN
        CALL HWWARN('HWANAL',2)
  201     PRINT 202, ICHSUM
  202     FORMAT(' Bad Charge Sum = ',I5/)
        CALL HWUEPR
        STOP
      ENDIF
 999  END


      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-5)
c
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
        if( (xplus/xminus).gt.tiny )then
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
      parameter (tiny=1.d-5)
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
