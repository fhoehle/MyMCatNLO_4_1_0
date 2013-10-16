C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      real * 8 xm0v,xm0h,binv,xmvlow,xmvupp,binh,xmhlow,xmhupp,bini,
     #  vgammax,hgammax,gav,gah
      real * 4 pi,xmvi,xmvs,xmhi,xmhs,xmii,xmis
      parameter (pi=3.14160E0)
      parameter (vgammax=30.d0)
      parameter (hgammax=30.d0)
      parameter (gah=0.5d0)
      integer j,k,jpr,idec,ivhvec,ivhlep
      common/vhlin/ivhvec,ivhlep
      character*5 cc(2)
      data cc/'     ',' cuts'/
c
      jpr=mod(abs(iproc),10000)/100
      if(jpr.eq.27)then
        xm0v=rmass(200)
        gav=gamz
      elseif(jpr.eq.26)then
        xm0v=rmass(198)
        gav=gamw
      else
        call hwwarn('HWABEG',500)
      endif
      xm0h=rmass(201)
c      gah=gamh
      if(ivhlep.gt.6)then
        idec=1
      else
        idec=0
      endif
c
      if(gav.eq.0)then
        binv=0.5d0
        xmvi=sngl(xm0v-24.75d0)
        xmvs=sngl(xm0v+25.25d0)
      else
        xmvlow=max(0.d0,xm0v-vgammax*gav)
        xmvupp=xm0v+vgammax*gav
        binv=(xmvupp-xmvlow)/100.d0
        xmvi=sngl(xm0v-(49*binv+binv/2))
        xmvs=sngl(xm0v+(50*binv+binv/2))
      endif
      if(gah.eq.0)then
        binh=0.5d0
        xmhi=sngl(xm0h-24.75d0)
        xmhs=sngl(xm0h+25.25d0)
      else
        xmhlow=max(0.d0,xm0h-hgammax*gah)
        xmhupp=xm0h+hgammax*gah
        binh=(xmhupp-xmhlow)/100.d0
        xmhi=sngl(xm0h-(49*binh+binh/2))
        xmhs=sngl(xm0h+(50*binh+binh/2))
      endif
      bini=1.d0
      xmii=sngl(xm0v+xm0h-4.25d0)
      xmis=sngl(xm0v+xm0h+95.75d0)
c
      call rinit("HERVH.root")
      do j=1,2
      k=(j-1)*50
c
      call rbook(k+ 1,'V pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 2,'V pt'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+ 3,'V log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+ 4,'V y'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+ 5,'V eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+ 6,'mV'//cc(j),sngl(binv),xmvi,xmvs)
c
      call rbook(k+ 7,'H pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+ 8,'H pt'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+ 9,'H log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+10,'H y'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+11,'H eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+12,'mH'//cc(j),sngl(binh),xmhi,xmhs)
c
      call rbook(k+13,'VH pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+14,'VH pt'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+15,'VH log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+16,'VH y'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+17,'VH eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+18,'mVH'//cc(j),sngl(bini),xmii,xmis)
      call rbook(k+19,'VH azimt'//cc(j),pi/20.e0,0.e0,pi)
      call rbook(k+20,'VH log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
      call rbook(k+21,'VH delta eta'//cc(j),0.2e0,-9.e0,9.e0)
c
      enddo
c
      if(idec.eq.0)then
      do j=1,2
      k=(j-1)*50
      call rbook(k+25+ 1,'l pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+25+ 2,'l pt'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+25+ 3,'l log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+25+ 4,'l eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+25+ 5,'lb pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+25+ 6,'lb pt'//cc(j),10.e0,0.e0,1000.e0)
      call rbook(k+25+ 7,'lb log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
      call rbook(k+25+ 8,'lb eta'//cc(j),0.2e0,-9.e0,9.e0)
c                       
      call rbook(k+25+ 9,'llb delta eta'//cc(j),0.2e0,-9.e0,9.e0)
      call rbook(k+25+10,'llb azimt'//cc(j),pi/20.e0,0.e0,pi)
      call rbook(k+25+11,'llb log[pi-azimt]'//cc(j),0.05e0,-4.e0,0.1e0)
      call rbook(k+25+12,'llb inv m'//cc(j),sngl(binv),xmvi,xmvs)
      call rbook(k+25+13,'llb pt'//cc(j),2.e0,0.e0,200.e0)
      call rbook(k+25+14,'llb log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
c
      enddo
      endif
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
      DOUBLE PRECISION HWVDOT,PSUM(4),YCUT,PPL(5),PPLB(5),
     #  PPV(5),PPH(5),PPVH(5),PTL,YL,THL,ETAL,PLL,ENL,PTLB,YLB,
     #  THLB,ETALB,PLLB,ENLB,PTLL,DLL,CLL,AZILL,AZILLNORM,XMLL,
     #  DETALLB,PTV,YV,THV,ETAV,PTH,YH,THH,ETAH,PTVH,YVH,TVHH,ETAVH,
     #  XMVH,DVH,CVH,AZIVH,AZIVHNORM,DETAVH
      INTEGER ICHSUM,ICHINI,IHEP,JPR,IDENT,IDEC,IST,ID,ID1,IH,
     #  ILL,ILLB,IJ,NIH,NIL,NILB,IVHVEC,IVHLEP,IV,NIV
      LOGICAL DIDSOF
      COMMON/VHLIN/IVHVEC,IVHLEP
      REAL*8 PI,GETRAPIDITY
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 WWW0,TINY
      INTEGER KK
      DATA TINY/.1D-5/
c
      IF (IERROR.NE.0) RETURN
c
      IF(IPROC.GT.0)CALL HWWARN('HWANAL',500)
      JPR=MOD(ABS(IPROC),10000)/100
      IF(JPR.EQ.27) THEN
        IDENT=23
      ELSEIF(JPR.EQ.26.AND.IVHVEC.EQ.1) THEN
        IDENT=24
      ELSEIF(JPR.EQ.26.AND.IVHVEC.EQ.-1) THEN
        IDENT=-24
      ELSE
        CALL HWWARN('HWANAL',502)
      ENDIF
      IF(IVHLEP.GT.6)THEN
        IDEC=1
      ELSE
        IDEC=0
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
      NIH=0
      NIV=0
      NIL=0
      NILB=0
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
  100 CONTINUE
      IF(IDEC.EQ.0.AND.ISTHEP(7).EQ.123.AND.IDHEP(7).EQ.25)THEN
         IH=JDAHEP(1,7)
         ILL=JDAHEP(1,JDAHEP(1,JDAHEP(1,JDAHEP(1,8))))
         ILLB=JDAHEP(1,JDAHEP(2,JDAHEP(1,JDAHEP(1,8))))
         NIH=1
         NIL=1
         NILB=1
      ELSEIF(IDEC.EQ.1.AND.ISTHEP(8).EQ.124.AND.IDHEP(8).EQ.25)THEN
         IH=JDAHEP(1,8)
         IV=JDAHEP(1,7)
         NIH=1
         NIV=1
      ELSE
        DO 200 IHEP=1,NHEP
          IST=ISTHEP(IHEP)      
          ID=IDHW(IHEP)
          ID1=IDHEP(IHEP)
          IF(IST.EQ.195.AND.ID1.EQ.25)THEN
            IH=IHEP
            NIH=NIH+1
          ENDIF
          IF(IDEC.EQ.0)THEN
            IF( IST.EQ.1 .AND. (ID1.GE.11.AND.ID1.LE.16) .AND.
     #          IDHEP(JMOHEP(1,JMOHEP(1,IHEP))).EQ.IDENT .AND. 
     #          ISTHEP(JMOHEP(1,JMOHEP(1,IHEP))).EQ.155 )THEN
              ILL=IHEP
              NIL=NIL+1
            ENDIF
            IF( IST.EQ.1 .AND. (ID1.GE.-16.AND.ID1.LE.-11).AND.
     #          IDHEP(JMOHEP(1,JMOHEP(1,IHEP))).EQ.IDENT .AND. 
     #          ISTHEP(JMOHEP(1,JMOHEP(1,IHEP))).EQ.155 )THEN
              ILLB=IHEP
              NILB=NILB+1
            ENDIF
          ELSE
            IF(IST.EQ.195.AND.ID1.EQ.IDENT)THEN
              IV=IHEP
              NIV=NIV+1
            ENDIF
          ENDIF
  200   CONTINUE
      ENDIF
      IF(IDEC.EQ.0)THEN
        IF(NIH.NE.1.OR.NIL.NE.1.OR.NILB.NE.1)THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',114)
          GOTO 999
        ENDIF
        IF( ABS(IDHEP(ILL)).LT.11 .OR. ABS(IDHEP(ILL)).GT.16 .OR.
     #      ABS(IDHEP(ILLB)).LT.11 .OR. ABS(IDHEP(ILLB)).GT.16 )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',505)
        ENDIF
        IF( (ISTHEP(ILL).NE.1.AND.ISTHEP(ILL).NE.195) .OR.
     #      (ISTHEP(ILLB).NE.1.AND.ISTHEP(ILLB).NE.195) )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',506)
        ENDIF
        IF( ISTHEP(IH).NE.195 )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',507)
        ENDIF
        DO IJ=1,5
          PPL(IJ)=PHEP(IJ,ILL)
          PPLB(IJ)=PHEP(IJ,ILLB)
          PPH(IJ)=PHEP(IJ,IH)
          PPV(IJ)=PPL(IJ)+PPLB(IJ)
          PPVH(IJ)=PPV(IJ)+PPH(IJ)
        ENDDO
        PPV(5)=SQRT(PPV(4)**2-PPV(1)**2-PPV(2)**2-PPV(3)**2)
      ELSE
        IF(NIH.NE.1.OR.NIV.NE.1)THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',114)
          GOTO 999
        ENDIF
        IF( ISTHEP(IH).NE.195.OR.ISTHEP(IV).NE.195 )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',506)
        ENDIF
        DO IJ=1,5
          PPH(IJ)=PHEP(IJ,IH)
          PPV(IJ)=PHEP(IJ,IV)
          PPVH(IJ)=PPV(IJ)+PPH(IJ)
        ENDDO
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
      IF(PBEAM1.GT.2500)THEN
        YCUT=2.5D0
      ELSE
        YCUT=1.0D0
      ENDIF
c
      if(idec.eq.0)then
        ptl=sqrt(ppl(1)**2+ppl(2)**2)
        yl=getrapidity(ppl(4),ppl(3))
        thl = atan2(ptl+tiny,ppl(3))
        etal= -log(tan(thl/2))
        pll = ppl(3)
        enl = ppl(4)
c
        ptlb=sqrt(pplb(1)**2+pplb(2)**2)
        ylb=getrapidity(pplb(4),pplb(3))
        thlb = atan2(ptlb+tiny,pplb(3))
        etalb= -log(tan(thlb/2))
        pllb = pplb(3)
        enlb = pplb(4)
c
        ptll = dsqrt((ppl(1)+pplb(1))**2+(ppl(2)+pplb(2))**2)
        dll = ppl(1)*pplb(1)+ppl(2)*pplb(2)
        cll=0
        if(ptl.ne.0.and.ptlb.ne.0) cll = dll * (1-tiny)/(ptl*ptlb)
        if(abs(cll).gt.1) then
	  write(*,*) ' cos[ll] = ',cll,dll,ptl,ptlb
          cll = - 1
        endif
        azill = (1-tiny)*acos(cll)
        azillnorm = (pi-azill)/pi
        xmll = dsqrt( ppl(5)**2 + pplb(5)**2 + 
     #                2*(enl*enlb - pll*pllb - dll) )
        detallb = etal-etalb
      endif
c
      ptv=sqrt(ppv(1)**2+ppv(2)**2)
      yv=getrapidity(ppv(4),ppv(3))
      thv = atan2(ptv+tiny,ppv(3))
      etav= -log(tan(thv/2))
c
      pth=sqrt(pph(1)**2+pph(2)**2)
      yh=getrapidity(pph(4),pph(3))
      thh = atan2(pth+tiny,pph(3))
      etah= -log(tan(thh/2))
c
      ptvh=sqrt(ppvh(1)**2+ppvh(2)**2)
      yvh=getrapidity(ppvh(4),ppvh(3))
      tvhh = atan2(ptvh+tiny,ppvh(3))
      etavh= -log(tan(tvhh/2))
      xmvh = dsqrt( ppv(5)**2 + pph(5)**2 + 
     #              2*( ppv(4)*pph(4)-ppv(1)*pph(1)-
     #                  ppv(2)*pph(2)-ppv(3)*pph(3) ) )
      dvh = ppv(1)*pph(1)+ppv(2)*pph(2)
      cvh=0
      if(ptv.ne.0.and.pth.ne.0) cvh = dvh * (1-tiny)/(ptv*pth)
      if(abs(cvh).gt.1) then
        write(*,*) ' cos[vh] = ',cvh,dvh,ptv,pth
        cvh = - 1
      endif
      azivh = (1-tiny)*acos(cvh)
      azivhnorm = (pi-azivh)/pi
      detavh = etav-etah
c
      kk=0
      call rfill(kk+1,sngl(ptv),sngl(WWW0))
      call rfill(kk+2,sngl(ptv),sngl(WWW0))
      if(ptv.gt.0.d0)call rfill(kk+3,sngl(log10(ptv)),sngl(WWW0))
      call rfill(kk+4,sngl(yv),sngl(WWW0))
      call rfill(kk+5,sngl(etav),sngl(WWW0))
      call rfill(kk+6,sngl(ppv(5)),sngl(WWW0))
c
      call rfill(kk+7,sngl(pth),sngl(WWW0))
      call rfill(kk+8,sngl(pth),sngl(WWW0))
      if(pth.gt.0.d0)call rfill(kk+9,sngl(log10(pth)),sngl(WWW0))
      call rfill(kk+10,sngl(yh),sngl(WWW0))
      call rfill(kk+11,sngl(etah),sngl(WWW0))
      call rfill(kk+12,sngl(pph(5)),sngl(WWW0))
c
      call rfill(kk+13,sngl(ptvh),sngl(WWW0))
      call rfill(kk+14,sngl(ptvh),sngl(WWW0))
      if(ptvh.gt.0.d0)call rfill(kk+15,sngl(log10(ptvh)),sngl(WWW0))
      call rfill(kk+16,sngl(yvh),sngl(WWW0))
      call rfill(kk+17,sngl(etavh),sngl(WWW0))
      call rfill(kk+18,sngl(xmvh),sngl(WWW0))
      call rfill(kk+19,sngl(azivh),sngl(WWW0))
      if(azivhnorm.gt.0.d0)
     #  call rfill(kk+20,sngl(log10(azivhnorm)),sngl(WWW0))
      call rfill(kk+21,sngl(detavh),sngl(WWW0))
c
      if(idec.eq.0)then
        kk=25
        call rfill(kk+1,sngl(ptl),sngl(WWW0))
        call rfill(kk+2,sngl(ptl),sngl(WWW0))
        if(ptl.gt.0.d0)call rfill(kk+3,sngl(log10(ptl)),sngl(WWW0))
        call rfill(kk+4,sngl(etal),sngl(WWW0))
        call rfill(kk+5,sngl(ptlb),sngl(WWW0))
        call rfill(kk+6,sngl(ptlb),sngl(WWW0))
        if(ptlb.gt.0.d0)call rfill(kk+7,sngl(log10(ptlb)),sngl(WWW0))
        call rfill(kk+8,sngl(etalb),sngl(WWW0))
c
        call rfill(kk+9,sngl(detallb),sngl(WWW0))
        call rfill(kk+10,sngl(azill),sngl(WWW0))
        if(azillnorm.gt.0.d0)
     #    call rfill(kk+11,sngl(log10(azillnorm)),sngl(WWW0))
        call rfill(kk+12,sngl(xmll),sngl(WWW0))
        call rfill(kk+13,sngl(ptll),sngl(WWW0))
        if(ptll.gt.0)call rfill(kk+14,sngl(log10(ptll)),sngl(WWW0))
      endif
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
     #  call rfill(kk+6,sngl(ppv(5)),sngl(WWW0))
c
      if(abs(etah).lt.ycut)then
        call rfill(kk+7,sngl(pth),sngl(WWW0))
        call rfill(kk+8,sngl(pth),sngl(WWW0))
        if(pth.gt.0.d0)call rfill(kk+9,sngl(log10(pth)),sngl(WWW0))
      endif
      if(pth.gt.20.d0)then
        call rfill(kk+10,sngl(yh),sngl(WWW0))
        call rfill(kk+11,sngl(etah),sngl(WWW0))
      endif
      if(abs(etah).lt.ycut.and.pth.gt.20.d0)
     #  call rfill(kk+12,sngl(pph(5)),sngl(WWW0))
c
      if( abs(etav).lt.ycut.and.abs(etah).lt.ycut .and.
     #    ptv.gt.20.d0.and.pth.gt.20.d0)then
        call rfill(kk+13,sngl(ptvh),sngl(WWW0))
        call rfill(kk+14,sngl(ptvh),sngl(WWW0))
        if(ptvh.gt.0.d0)call rfill(kk+15,sngl(log10(ptvh)),sngl(WWW0))
        call rfill(kk+16,sngl(yvh),sngl(WWW0))
        call rfill(kk+17,sngl(etavh),sngl(WWW0))
        call rfill(kk+18,sngl(xmvh),sngl(WWW0))
        call rfill(kk+19,sngl(azivh),sngl(WWW0))
        if(azivhnorm.gt.0.d0)
     #    call rfill(kk+20,sngl(log10(azivhnorm)),sngl(WWW0))
        call rfill(kk+21,sngl(detavh),sngl(WWW0))
      endif
c
      if(idec.eq.0)then
        kk=75
        if(abs(etal).lt.ycut)then
          call rfill(kk+1,sngl(ptl),sngl(WWW0))
          call rfill(kk+2,sngl(ptl),sngl(WWW0))
          if(ptl.gt.0.d0)call rfill(kk+3,sngl(log10(ptl)),sngl(WWW0))
        endif
        if(ptl.gt.20.d0)call rfill(kk+4,sngl(etal),sngl(WWW0))
        if(abs(etalb).lt.ycut)then
          call rfill(kk+5,sngl(ptlb),sngl(WWW0))
          call rfill(kk+6,sngl(ptlb),sngl(WWW0))
          if(ptlb.gt.0.d0)call rfill(kk+7,sngl(log10(ptlb)),sngl(WWW0))
        endif
        if(ptlb.gt.20.d0)call rfill(kk+8,sngl(etalb),sngl(WWW0))
c
        if( abs(etal).lt.ycut.and.abs(etalb).lt.ycut .and.
     #      ptl.gt.20.d0.and.ptlb.gt.20.d0)then
          call rfill(kk+9,sngl(detallb),sngl(WWW0))
          call rfill(kk+10,sngl(azill),sngl(WWW0))
          if(azillnorm.gt.0.d0)
     #      call rfill(kk+11,sngl(log10(azillnorm)),sngl(WWW0))
          call rfill(kk+12,sngl(xmll),sngl(WWW0))
          call rfill(kk+13,sngl(ptll),sngl(WWW0))
          if(ptll.gt.0) 
     #      call rfill(kk+14,sngl(log10(ptll)),sngl(WWW0))
        endif
      endif
c
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
