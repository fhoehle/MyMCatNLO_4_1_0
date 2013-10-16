C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      real * 8 xm0,xmlow,xmupp,bin
      real * 4 xmi,xms,pi
      parameter (pi=3.14160E0)
      integer j,k,jpr
      character*5 cc(2)
      data cc/'     ',' cuts'/
c
      jpr=mod(abs(iproc),10000)
      if(jpr.eq.1396) then
        xm0=(emmin+emmax)/2.d0
        xmupp=emmax
        xmlow=emmin
      elseif(jpr.eq.1397) then
        xm0=rmass(200)
        xmupp=xm0+gammax*gamz
        xmlow=xm0-gammax*gamz
      elseif(jpr.eq.1497) then
        xm0=rmass(198)
        xmupp=xm0+gammax*gamw
        xmlow=xm0-gammax*gamw
      elseif(jpr.eq.1498) then
        xm0=rmass(199)
        xmupp=xm0+gammax*gamw
        xmlow=xm0-gammax*gamw
      else
        call hwwarn('HWABEG',501)
      endif
      if(abs(xmlow-xmupp).lt.1.d-3)then
        bin=0.5d0
        xmi=sngl(xm0-24.75d0)
        xms=sngl(xm0+25.25d0)
      else
        bin=(xmupp-xmlow)/100.d0
        xmi=sngl(xm0-(49*bin+bin/2))
        xms=sngl(xm0+(50*bin+bin/2))
      endif
      call rinit("HERSB.root")
      do j=1,2
      k=(j-1)*50
c
      call rbook(k+ 1,'V pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 2,'V pt'//cc(j),100.e0,0.e0,6000.e0)
      call rbook(k+ 3,'V log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+ 4,'V y'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+ 5,'V eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+ 6,'mV'//cc(j),sngl(bin),xmi,xms)
c
      call rbook(k+ 7,'l pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 8,'l pt'//cc(j),100.e0,0.e0,6000.e0)
      call rbook(k+ 9,'l log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+10,'l eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+11,'lb pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+12,'lb pt'//cc(j),100.e0,0.e0,6000.e0)
      call rbook(k+13,'lb log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+14,'lb eta'//cc(j),0.2e0,-9.e0,9.e0)
c
      call rbook(k+15,'llb delta eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+16,'llb azimt'//cc(j),pi/20.e0,0.e0,pi)
      call rbook(k+17,'llb log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
      call rbook(k+18,'llb inv m'//cc(j),sngl(bin),xmi,xms)
      call rbook(k+19,'llb pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+20,'llb log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
c
      enddo
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
      DOUBLE PRECISION HWVDOT,PSUM(4),PPV(5),YCUT,XMV,PTV,YV,THV,ETAV,
     #  PPL(5),PPLB(5),PTL,YL,THL,ETAL,PLL,ENL,PTLB,YLB,THLB,ETALB,
     #  PLLB,ENLB,PTPAIR,DLL,CLL,AZI,AZINORM,XMLL,DETALLB
      INTEGER ICHSUM,ICHINI,IHEP,IV,IFV,IST,ID,IJ,ID1,JPR,IDENT,
     #  ILL,ILLB
      LOGICAL DIDSOF,TEST1,TEST2
      REAL*8 PI
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 WWW0,TINY
      INTEGER KK
      DATA TINY/.1D-5/
c
      IF (IERROR.NE.0) RETURN
c
      JPR=MOD(ABS(IPROC),10000)
      IF(JPR.EQ.1396) THEN
        IDENT=23
      ELSEIF(JPR.EQ.1397) THEN
        IDENT=23
      ELSEIF(JPR.EQ.1497) THEN
        IDENT=24
      ELSEIF(JPR.EQ.1498) THEN
        IDENT=-24
      ELSE
        CALL HWWARN('HWANAL',502)
      ENDIF
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
      IFV=0
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
        IST=ISTHEP(IHEP)      
        ID=IDHW(IHEP)
        ID1=IDHEP(IHEP)
        TEST1=IST.EQ.195
        TEST2=ID1.EQ.IDENT
        IF(TEST1.AND.TEST2)THEN
          IV=IHEP
          IFV=IFV+1
          DO IJ=1,5
	    PPV(IJ)=PHEP(IJ,IHEP)
	  ENDDO
        ENDIF
  100 CONTINUE
      IF(IFV.EQ.0.AND.IERROR.EQ.0) THEN
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
      IF(IFV.GT.1.AND.IERROR.EQ.0) THEN
         CALL HWUEPR
         CALL HWWARN('HWANAL',55)
         GOTO 999
      ENDIF
C FIND THE LEPTONS
      IF( IDHEP(JDAHEP(1,IV)).GT.0 .AND.
     #    IDHEP(JDAHEP(2,IV)).LT.0 )THEN
        ILL=JDAHEP(1,IV)
        ILLB=JDAHEP(2,IV)
      ELSEIF( IDHEP(JDAHEP(2,IV)).GT.0 .AND.
     #        IDHEP(JDAHEP(1,IV)).LT.0 )THEN
        ILL=JDAHEP(2,IV)
        ILLB=JDAHEP(1,IV)
      ELSE
        CALL HWUEPR
        CALL HWWARN('HWANAL',504)
      ENDIF
      DO IJ=1,5
        PPL(IJ)=PHEP(IJ,ILL)
        PPLB(IJ)=PHEP(IJ,ILLB)
      ENDDO
C CHECK THAT THE LEPTONS ARE FINAL-STATE LEPTONS
      IF( ABS(IDHEP(ILL)).LT.11 .OR. ABS(IDHEP(ILL)).GT.16 .OR.
     #    ABS(IDHEP(ILLB)).LT.11 .OR. ABS(IDHEP(ILLB)).GT.16 )THEN
        CALL HWUEPR
        CALL HWWARN('HWANAL',505)
      ENDIF
      IF( (ISTHEP(ILL).NE.1.AND.ISTHEP(ILL).NE.195) .OR.
     #    (ISTHEP(ILLB).NE.1.AND.ISTHEP(ILLB).NE.195) )THEN
        CALL HWUEPR
        CALL HWWARN('HWANAL',506)
      ENDIF
C FILL THE HISTOS
      IF(PBEAM1.GT.2500)THEN
        YCUT=2.5D0
      ELSE
        YCUT=1.0D0
      ENDIF
C Variables of the vector boson
      xmv=ppv(5)
      ptv=sqrt(ppv(1)**2+ppv(2)**2)
      if(abs(ppv(4)-abs(ppv(3))).gt.tiny)then
        yv=0.5d0*log( (ppv(4)+ppv(3))/
     #                (ppv(4)-ppv(3)) )
      else
        yv=sign(1.d0,ppv(3))*1.d8
      endif
      thv = atan2(ptv+tiny,ppv(3))
      etav= -log(tan(thv/2))
C Variables of the leptons
      ptl=sqrt(ppl(1)**2+ppl(2)**2)
      if(abs(ppl(4)-abs(ppl(3))).gt.tiny)then
        yl=0.5d0*log( (ppl(4)+ppl(3))/
     #                (ppl(4)-ppl(3)) )
      else
        yl=sign(1.d0,ppl(3))*1.d8
      endif
      thl = atan2(ptl+tiny,ppl(3))
      etal= -log(tan(thl/2))
      pll = ppl(3)
      enl = ppl(4)
c
      ptlb=sqrt(pplb(1)**2+pplb(2)**2)
      if(abs(pplb(4)-abs(pplb(3))).gt.tiny)then
        ylb=0.5d0*log( (pplb(4)+pplb(3))/
     #                 (pplb(4)-pplb(3)) )
      else
        ylb=sign(1.d0,pplb(3))*1.d8
      endif
      thlb = atan2(ptlb+tiny,pplb(3))
      etalb= -log(tan(thlb/2))
      pllb = pplb(3)
      enlb = pplb(4)
c
      ptpair = dsqrt((ppl(1)+pplb(1))**2+(ppl(2)+pplb(2))**2)
      dll = ppl(1)*pplb(1)+ppl(2)*pplb(2)
      cll=0
      if(ptl.ne.0.and.ptlb.ne.0) cll = dll * (1-tiny)/(ptl*ptlb)
      if(abs(cll).gt.1) then
	write(*,*) ' cosine = ',cll ,dll,ptl,ptlb
	cll = - 1
      endif
      azi = (1-tiny)*acos(cll)
      azinorm = (pi-azi)/pi
      xmll = dsqrt( ppl(5)**2 + ppl(5)**2 + 
     #              2*(enl*enlb - pll*pllb - dll) )
      detallb = etal-etalb
c
      kk=0
      call rfill(kk+1,sngl(ptv),sngl(WWW0))
      call rfill(kk+2,sngl(ptv),sngl(WWW0))
      if(ptv.gt.0.d0)call rfill(kk+3,sngl(log10(ptv)),sngl(WWW0))
      call rfill(kk+4,sngl(yv),sngl(WWW0))
      call rfill(kk+5,sngl(etav),sngl(WWW0))
      call rfill(kk+6,sngl(xmv),sngl(WWW0))
c
      call rfill(kk+7,sngl(ptl),sngl(WWW0))
      call rfill(kk+8,sngl(ptl),sngl(WWW0))
      if(ptl.gt.0.d0)call rfill(kk+9,sngl(log10(ptl)),sngl(WWW0))
      call rfill(kk+10,sngl(etal),sngl(WWW0))
      call rfill(kk+11,sngl(ptlb),sngl(WWW0))
      call rfill(kk+12,sngl(ptlb),sngl(WWW0))
      if(ptlb.gt.0.d0)call rfill(kk+13,sngl(log10(ptlb)),sngl(WWW0))
      call rfill(kk+14,sngl(etalb),sngl(WWW0))
c
      call rfill(kk+15,sngl(detallb),sngl(WWW0))
      call rfill(kk+16,sngl(azi),sngl(WWW0))
      if(azinorm.gt.0.d0)
     #  call rfill(kk+17,sngl(log10(azinorm)),sngl(WWW0))
      call rfill(kk+18,sngl(xmll),sngl(WWW0))
      call rfill(kk+19,sngl(ptpair),sngl(WWW0))
      if(ptpair.gt.0)call rfill(kk+20,sngl(log10(ptpair)),sngl(WWW0))
c
      kk=50
      if(abs(etav).lt.ycut)then
        call rfill(kk+1,sngl(ptv),sngl(WWW0))
        call rfill(kk+2,sngl(ptv),sngl(WWW0))
        if(ptv.gt.0.d0)call rfill(kk+3,sngl(log10(ptv)),sngl(WWW0))
      endif
      if(ptv.gt.20.d0)then
        call rfill(kk+4,sngl(yv),sngl(WWW0))
        call rfill(kk+5,sngl(etav),sngl(WWW0))
      endif
      if(abs(etav).lt.ycut.and.ptv.gt.20.d0)
     #  call rfill(kk+6,sngl(xmv),sngl(WWW0))
c
      if(abs(etal).lt.ycut)then
        call rfill(kk+7,sngl(ptl),sngl(WWW0))
        call rfill(kk+8,sngl(ptl),sngl(WWW0))
        if(ptl.gt.0.d0)call rfill(kk+9,sngl(log10(ptl)),sngl(WWW0))
      endif
      if(ptl.gt.20.d0)call rfill(kk+10,sngl(etal),sngl(WWW0))
      if(abs(etalb).lt.ycut)then
        call rfill(kk+11,sngl(ptlb),sngl(WWW0))
        call rfill(kk+12,sngl(ptlb),sngl(WWW0))
        if(ptlb.gt.0.d0)call rfill(kk+13,sngl(log10(ptlb)),sngl(WWW0))
      endif
      if(ptlb.gt.20.d0)call rfill(kk+14,sngl(etalb),sngl(WWW0))
c
      if( abs(etal).lt.ycut.and.abs(etalb).lt.ycut .and.
     #    ptl.gt.20.d0.and.ptlb.gt.20.d0)then
        call rfill(kk+15,sngl(detallb),sngl(WWW0))
        call rfill(kk+16,sngl(azi),sngl(WWW0))
        if(azinorm.gt.0.d0)
     #    call rfill(kk+17,sngl(log10(azinorm)),sngl(WWW0))
        call rfill(kk+18,sngl(xmll),sngl(WWW0))
        call rfill(kk+19,sngl(ptpair),sngl(WWW0))
        if(ptpair.gt.0) 
     #    call rfill(kk+20,sngl(log10(ptpair)),sngl(WWW0))
      endif
 999  END
