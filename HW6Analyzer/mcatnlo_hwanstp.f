C----------------------------------------------------------------------
      SUBROUTINE RCLOS()
C     DUMMY IF HBOOK IS USED
C----------------------------------------------------------------------
      END


C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      REAL*4 pi
      parameter (pi=3.14160E0)
      integer j,k,ivlep1,ivlep2,idec
      common/vvlin/ivlep1,ivlep2
      character*5 cc(2)
      data cc/'     ',' cuts'/
c
      if(ivlep1.eq.7)then
        idec=1
      else
        if(ivlep1.gt.3)call hwwarn('HWABEG',501)
        idec=0
      endif
c
      call inihist
      do j=1,2
        k=(j-1)*50
        call mbook(k+ 1,'t pt'//cc(j),4.e0,0.e0,400.e0)
        call mbook(k+ 2,'t eta'//cc(j),0.2e0,-9.e0,9.e0)
        call mbook(k+ 3,'tb pt'//cc(j),4.e0,0.e0,400.e0)
        call mbook(k+ 4,'tb eta'//cc(j),0.2e0,-9.e0,9.e0)
        call mbook(k+ 5,'t log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
        call mbook(k+ 6,'t y'//cc(j),0.2e0,-9.e0,9.e0)
        call mbook(k+ 7,'tb log[pt]'//cc(j),0.05e0,0.1e0,5.e0)
        call mbook(k+ 8,'tb y'//cc(j),0.2e0,-9.e0,9.e0)
        if(idec.eq.0)then
          call mbook(k+11,'l+ pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+12,'l+ y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+13,'nu pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+14,'nu y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+15,'b pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+16,'b y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+19,'M[W+]'//cc(j),2.e0,0.e0,200.e0)
          call mbook(k+21,'l- pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+22,'l- y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+23,'nubar pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+24,'nubar y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+25,'bbar pt'//cc(j),4.e0,0.e0,400.e0)
          call mbook(k+26,'bbar y'//cc(j),0.2e0,-4.e0,4.e0)
          call mbook(k+29,'M[W-]'//cc(j),2.e0,0.e0,200.e0)
        endif
      enddo
      END


C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      REAL*8 XNORM
      INTEGER I,J,K,IVLEP1,IVLEP2,IDEC
      COMMON/VVLIN/IVLEP1,IVLEP2
      OPEN(UNIT=99,FILE='HERST.TOP',STATUS='UNKNOWN')
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
      do j=1,2
        k=(j-1)*50
        call multitop(100+k+ 1,99,3,2,'t pt',' ','LOG')
        call multitop(100+k+ 2,99,3,2,'t eta',' ','LOG')
        call multitop(100+k+ 3,99,3,2,'tb pt',' ','LOG')
        call multitop(100+k+ 4,99,3,2,'tb eta',' ','LOG')
        call multitop(100+k+ 5,99,3,2,'t log[pt]',' ','LOG')
        call multitop(100+k+ 6,99,3,2,'t y',' ','LOG')
        call multitop(100+k+ 7,99,3,2,'tb log[pt]',' ','LOG')
        call multitop(100+k+ 8,99,3,2,'tb y',' ','LOG')
        if(idec.eq.0)then
          call multitop(100+k+11,99,3,2,'l+ pt',' ','LOG')
          call multitop(100+k+12,99,3,2,'l+ y',' ','LOG')
          call multitop(100+k+13,99,3,2,'nu pt',' ','LOG')
          call multitop(100+k+14,99,3,2,'nu y',' ','LOG')
          call multitop(100+k+15,99,3,2,'b pt',' ','LOG')
          call multitop(100+k+16,99,3,2,'b y',' ','LOG')
          call multitop(100+k+19,99,3,2,'M[W+]',' ','LOG')
          call multitop(100+k+21,99,3,2,'l- pt',' ','LOG')
          call multitop(100+k+22,99,3,2,'l- y',' ','LOG')
          call multitop(100+k+23,99,3,2,'nubar pt',' ','LOG')
          call multitop(100+k+24,99,3,2,'nubar y',' ','LOG')
          call multitop(100+k+25,99,3,2,'bbar pt',' ','LOG')
          call multitop(100+k+26,99,3,2,'bbar y',' ','LOG')
          call multitop(100+k+29,99,3,2,'M[W-]',' ','LOG')
        endif
      enddo
c
      CLOSE(99)
      END


C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION HWVDOT,PSUM(4)
      INTEGER ICHSUM,ICHINI,IHEP
      LOGICAL DIDSOF,sicuts,flcuts
      INTEGER ID,ID1,IST,IQ,IT1,ILP,INU,IBQ,IJ
      DOUBLE PRECISION YCUT,PTCUT,pt1,eta1,getpseudorap,yt1,
     # getrapidity,ptlp,ylp,ptnu,ynu,ptbq,ybq,xmw1,getinvm
      DOUBLE PRECISION XPTQ(5),XPLP(5),XPNU(5),XPBQ(5),YPW1(5)
      REAL*8 PI
      PARAMETER (PI=3.14159265358979312D0)
      REAL*8 WWW0
      INTEGER KK,JJ,IVLEP1,IVLEP2,IDEC
      COMMON/VVLIN/IVLEP1,IVLEP2
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
      IQ=0
      IF(IVLEP1.EQ.7)THEN
        IDEC=1
      ELSE
        IDEC=0
      ENDIF
      DO 100 IHEP=1,NHEP
        IF (IDHW(IHEP).EQ.16) DIDSOF=.TRUE.
        IF (ISTHEP(IHEP).EQ.1) THEN
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))
        ENDIF
        IST=ISTHEP(IHEP)      
        ID=IDHW(IHEP)
        ID1=IDHEP(IHEP)
        IF(IST.EQ.155.AND.ID1.EQ.6)THEN
C FOUND A TOP; KEEP ONLY THE FIRST ON RECORD
          IQ=IQ+1
          IF(IQ.EQ.1)IT1=IHEP
        ELSEIF(IST.EQ.155.AND.ID1.EQ.-6)THEN
C FOUND AN ANTITOP; KEEP ONLY THE FIRST ON RECORD
          IQ=IQ+1
          IF(IQ.EQ.1)IT1=IHEP
        ENDIF
  100 CONTINUE
      IF(IQ.EQ.0.AND.IERROR.EQ.0)CALL HWWARN('HWANAL',501)
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
C FIND THE DECAY PRODUCTS IF SPIN CORRELATIONS ARE INCLUDED
      IF(IDEC.EQ.0)THEN
        ILP=JDAHEP(1,JDAHEP(1,JDAHEP(1,JDAHEP(1,JDAHEP(1,IT1)))))
        INU=JDAHEP(1,JDAHEP(2,JDAHEP(1,JDAHEP(1,JDAHEP(1,IT1)))))
        IBQ=0
        IF(JDAHEP(2,IT1)-JDAHEP(1,IT1).EQ.1)THEN
          DO IJ=JDAHEP(1,JDAHEP(1,JDAHEP(2,IT1))),
     #          JDAHEP(2,JDAHEP(1,JDAHEP(2,IT1)))
            IF(ISTHEP(IJ).EQ.2.AND.
     #         IDHEP(IJ).EQ.5*SIGN(1,IDHEP(IT1)))IBQ=IJ
          ENDDO
        ELSEIF(JDAHEP(2,IT1)-JDAHEP(1,IT1).EQ.2)THEN
          DO IJ=JDAHEP(1,JDAHEP(1,JDAHEP(2,IT1)-1)),
     #          JDAHEP(2,JDAHEP(1,JDAHEP(2,IT1)-1))
            IF(ISTHEP(IJ).EQ.2.AND.
     #         IDHEP(IJ).EQ.5*SIGN(1,IDHEP(IT1)))IBQ=IJ
          ENDDO
        ELSE
          CALL HWUEPR
          CALL HWWARN('HWANAL',502)
        ENDIF
C CHECK THAT THE DECAY PRODUCTS ARE WHAT THEY ARE SUPPOSED TO BE
        IF( ( (ABS(IDHEP(ILP)).NE.11.AND.ABS(IDHEP(ILP)).NE.13.AND.
     #         ABS(IDHEP(ILP)).NE.15).OR.
     #      (ISTHEP(ILP).NE.1.AND.ISTHEP(ILP).NE.195) ) .OR.
     #      ( (ABS(IDHEP(INU)).NE.12.AND.ABS(IDHEP(INU)).NE.14.AND.
     #         ABS(IDHEP(INU)).NE.16).OR.
     #        (ISTHEP(INU).NE.1.AND.ISTHEP(INU).NE.195) ) .OR.
     #      ( ABS(IDHEP(IBQ)).NE.5 .OR. ISTHEP(IBQ).NE.2 ) .OR.
     #        SIGN(1,IDHEP(IT1)).NE.-SIGN(1,IDHEP(ILP)) )THEN
          CALL HWUEPR
          CALL HWWARN('HWANAL',504)
        ENDIF
      ENDIF
C FILL THE FOUR-MOMENTA
      DO IJ=1,5
        XPTQ(IJ)=PHEP(IJ,IT1)
        IF(IDEC.EQ.0)THEN
          XPLP(IJ)=PHEP(IJ,ILP)
          XPNU(IJ)=PHEP(IJ,INU)
          XPBQ(IJ)=PHEP(IJ,IBQ)
        ENDIF
      ENDDO
      IF(IDEC.EQ.0)THEN
        DO IJ=1,4
          YPW1(IJ)=XPLP(IJ)+XPNU(IJ)
        ENDDO
      ENDIF
C FILL THE HISTOS
      IF(PBEAM1.GT.2500)THEN
        YCUT=2.5D0
        PTCUT=30.D0
      ELSE
        YCUT=1.0D0
        PTCUT=15.D0
      ENDIF
      JJ=0
      IF(SIGN(1,IDHEP(IT1)).EQ.-1)JJ=2
c
      pt1=sqrt(xptq(1)**2+xptq(2)**2)
      eta1=getpseudorap(xptq(4),xptq(1),xptq(2),xptq(3))
      yt1=getrapidity(xptq(4),xptq(3))
      sicuts=abs(yt1).le.ycut.and.pt1.ge.ptcut
      if(idec.eq.0)then
        ptlp=sqrt(xplp(1)**2+xplp(2)**2)
        ylp=getrapidity(xplp(4),xplp(3))
        ptnu=sqrt(xpnu(1)**2+xpnu(2)**2)
        ynu=getrapidity(xpnu(4),xpnu(3))
        ptbq=sqrt(xpbq(1)**2+xpbq(2)**2)
        ybq=getrapidity(xpbq(4),xpbq(3))
        xmw1=getinvm(ypw1(4),ypw1(1),ypw1(2),ypw1(3))
        flcuts=abs(ybq).le.ycut.and.ptbq.ge.ptcut .and.
     #         abs(ylp).le.ycut.and.ptlp.ge.ptcut
      endif
C
C WITHOUT CUTS
C
      kk=0
      call mfill(kk+jj+1,sngl(pt1),sngl(WWW0))
      call mfill(kk+jj+2,sngl(eta1),sngl(WWW0))
      if(pt1.gt.0.d0)call mfill(kk+jj+5,sngl(log10(pt1)),sngl(WWW0))
      call mfill(kk+jj+6,sngl(yt1),sngl(WWW0))
      if(idec.eq.0)then
        call mfill(kk+5*jj+11,sngl(ptlp),sngl(WWW0))
        call mfill(kk+5*jj+12,sngl(ylp),sngl(WWW0))
        call mfill(kk+5*jj+13,sngl(ptnu),sngl(WWW0))
        call mfill(kk+5*jj+14,sngl(ynu),sngl(WWW0))
        call mfill(kk+5*jj+15,sngl(ptbq),sngl(WWW0))
        call mfill(kk+5*jj+16,sngl(ybq),sngl(WWW0))
        call mfill(kk+5*jj+19,sngl(xmw1),sngl(WWW0))
      endif
C
C WITH CUTS
C
      kk=50
      if(sicuts)then
        call mfill(kk+jj+1,sngl(pt1),sngl(WWW0))
        call mfill(kk+jj+2,sngl(eta1),sngl(WWW0))
        if(pt1.gt.0.d0)call mfill(kk+jj+5,sngl(log10(pt1)),sngl(WWW0))
        call mfill(kk+jj+6,sngl(yt1),sngl(WWW0))
      endif
      if(idec.eq.0.and.flcuts)then
        call mfill(kk+5*jj+11,sngl(ptlp),sngl(WWW0))
        call mfill(kk+5*jj+12,sngl(ylp),sngl(WWW0))
        call mfill(kk+5*jj+13,sngl(ptnu),sngl(WWW0))
        call mfill(kk+5*jj+14,sngl(ynu),sngl(WWW0))
        call mfill(kk+5*jj+15,sngl(ptbq),sngl(WWW0))
        call mfill(kk+5*jj+16,sngl(ybq),sngl(WWW0))
        call mfill(kk+5*jj+19,sngl(xmw1),sngl(WWW0))
      endif
 999  RETURN
      END


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
