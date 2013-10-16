C----------------------------------------------------------------------
      SUBROUTINE RCLOS()
C     DUMMY IF HBOOK IS USED
C----------------------------------------------------------------------
      END


C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      real * 8 xmh0
      real * 4 xmhi,xmhs
      integer j,k
      character*5 cc(2)
      data cc/'     ',' cuts'/
c
      xmh0=RMASS(201)
      xmhi=sngl(xmh0-24.75d0)
      xmhs=sngl(xmh0+25.25d0)
      call inihist
      do j=1,2
      k=(j-1)*50
      call mbook(k+ 1,'Higgs pt'//cc(j),2.e0,0.e0,200.e0)
      call mbook(k+ 5,'Higgs log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call mbook(k+ 9,'Higgs y'//cc(j),0.2e0,-4.e0,4.e0)
      call mbook(k+13,'jet y'//cc(j),0.2e0,-4.e0,4.e0)
      call mbook(k+17,'mH'//cc(j),0.5e0,xmhi,xmhs)
      call mbook(k+21,'Higgs pt'//cc(j),2.e0,200.e0,400.e0)
      call mbook(k+25,'Higgs pt'//cc(j),10.e0,0.e0,500.e0)
      enddo
      END


C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      REAL*8 XNORM
      INTEGER I,J,K
      OPEN(UNIT=99,FILE='HERHG.TOP',STATUS='UNKNOWN')
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
      do j=1,2
      k=(j-1)*50
      call multitop(100+k+ 1,99,2,2,'Higgs pt',' ','LOG')
      call multitop(100+k+ 5,99,2,2,'Higgs log[pt]',' ','LOG')
      call multitop(100+k+ 9,99,2,2,'Higgs y',' ','LOG')
      call multitop(100+k+13,99,2,2,'jet y',' ','LOG')
      enddo
      do j=1,2
      k=(j-1)*50
      call multitop(100+k+17,99,2,2,'mH',' ','LOG')
      call multitop(100+k+21,99,2,2,'Higgs pt',' ','LOG')
      call multitop(100+k+25,99,2,2,'Higgs pt',' ','LOG')
      enddo
c
      CLOSE(99)
      END

C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWVDOT,PSUM(4),PPH(5),XMH,PTH,YH,PTP,YP
      INTEGER ICHSUM,ICHINI,IHEP,IFH,IST,ID,IJ,ID1
      LOGICAL DIDSOF
      REAL*8 WWW0,TINY
      INTEGER KK
      DATA TINY/.1D-5/
c
      IF (IERROR.NE.0) RETURN
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
      IFH=0
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
        IST=ISTHEP(IHEP)      
        ID=IDHW(IHEP)
        ID1=IDHEP(IHEP)
        IF(IST.EQ.195)THEN
          IF(ID1.EQ.25)THEN
            IFH=IFH+1
            DO IJ=1,5
	      PPH(IJ)=PHEP(IJ,IHEP)
	    ENDDO
          ENDIF
        ENDIF
  100 CONTINUE
      IF(IFH.NE.1.AND.IERROR.EQ.0)THEN
         CALL HWUEPR
         CALL HWWARN('HWANAL',501)
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
C FILL THE HISTOS
      xmh=pph(5)
      pth=sqrt(pph(1)**2+pph(2)**2)
      if(abs(pph(4)-abs(pph(3))).gt.tiny)then
        yh=0.5d0*log( (pph(4)+pph(3))/
     #                (pph(4)-pph(3)) )
      else
        yh=sign(1.d0,pph(3))*1.d8
      endif
C May define proper jet variables eventually
      ptp=pth
      yp=yh
c
      kk=0
      call mfill(kk+1,sngl(pth),sngl(WWW0))
      call mfill(kk+21,sngl(pth),sngl(WWW0))
      call mfill(kk+25,sngl(pth),sngl(WWW0))
      if(pth.gt.0.d0)call mfill(kk+5,sngl(log10(pth)),sngl(WWW0))
      call mfill(kk+9,sngl(yh),sngl(WWW0))
      call mfill(kk+13,sngl(yh),sngl(WWW0))
      call mfill(kk+17,sngl(xmh),sngl(WWW0))
      kk=50
      if(abs(yh).lt.2.d0)then
        call mfill(kk+1,sngl(pth),sngl(WWW0))
        call mfill(kk+21,sngl(pth),sngl(WWW0))
        call mfill(kk+25,sngl(pth),sngl(WWW0))
        if(pth.gt.0.d0)call mfill(kk+5,sngl(log10(pth)),sngl(WWW0))
      endif
      if(pth.gt.20.d0)then
        call mfill(kk+9,sngl(yh),sngl(WWW0))
        call mfill(kk+13,sngl(yp),sngl(WWW0))
      endif
      if(abs(yh).lt.2.d0.and.pth.gt.20.d0)
     #  call mfill(kk+17,sngl(xmh),sngl(WWW0))
 999  END
