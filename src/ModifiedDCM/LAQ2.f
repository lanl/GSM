******************************************************************
*          For  LAQGSM-MARS interface      KKG  04/09/07
******************************************************************
* aa2k7g.f,  gengamn7.f, gadd7.f,    ->     laqgsm2007_1.f
* qgsm7.f                            ->     laqgsm2002_2.f
* coales07.f, gemaa07.f, preco07g.f, spegem7g.f, fermi07.f,
* gambetm7.f                         ->     laqgsm2004_3.f
*
c
******************************************************************
c          laqgsm2007_2.f   subroutines :
******************************************************************
c
c last changes by KKG for interface with MARS code       23.03.07
c     last correction in CLUSLE, REP/KKG, Oct. 2004
c     sigma=0.51             G. 25.05.01
C   * * * * * * * * *  11-16-94 03:27PM  MV <==> 5999 * * * * * *
      SUBROUTINE HEINEN(PIN,IIN,PN,IPN,MV,NP,NIN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL *8 PIN,PN
      CHARACTER*8 DTYP(11),TDIAG
      COMMON/MEMORYLAQ/ PME(9,5999),IME(5,5999)
      DIMENSION  PIN(9),IIN(5),PN(9),IPN(5)
      COMMON /IACT/ IACT
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON /LOWMIS/ LOWMIS
      COMMON /CSLID/ CLIDER(5999)
      COMMON/COMENB/ ENBOU
      COMMON/IORHEI/ IORHEI
      COMMON/ITHEA/ITHEA(11)
      COMMON /HELE/ IHELE
      COMMON /IDPME/ IDPME(5999)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      DIMENSION PSU(3)
      DATA DTYP/'DIFTRI','PLANAR','UNCYLI','ELAST  ','ANNIH','DIFSMA',
     *          'CHAINS','BINAR ','??????','REGTRI ','DOUBDI'/
      sigmadf=SIGMA
c      IF(NCAS.GE.NCPRI)
c     *write(16,601) PIN,IIN,PN,IPN
  601 FORMAT(1X,'HEINEN:',9(1X,E10.3),4I2,I10/
     *       1X,'      +',9(1X,E10.3),4I2,I10)
      ES0=SQRT((PIN(8)+PIN(9)+PN(8)+PN(9))**2-(PIN(4)+PN(4))**2-
     -(PIN(5)+PN(5))**2-(PIN(6)+PN(6))**2)
      CS0=DBLE(IIN(1)+IPN(1))
      SS0=DBLE(IIN(3)+IPN(3))
      BS0=DBLE(IIN(4)+IPN(4))
      ESU=0.
      CSU=0.
      SSU=0.
      BSU=0.
      LOWMIS=0
      DO  10  K=1,3
   10 PSU(K)=0.
c
      CALL  ACTNAM(PIN,PN,MV,NP,NIN)
c
      DO  19  I=1,11
   19 IF(ITHEA(I).NE.0) TDIAG=DTYP(I)
      if(SIGMA.ne.sigmadf)  then
        write( *,2001) sigmadf,SIGMA,TDIAG
 2001   format(' Default parameter SIGMA(',F7.3,
     &         ' ) was changed(=',F7.3,' in ',A8)
      endif
      IF(NIN.EQ.1)  RETURN
      DO  14   K=1,NP
      M=MV+K
      CSU=CSU+DBLE(IME(1,M))
      SSU=SSU+DBLE(IME(3,M))
      BSU=BSU+DBLE(IME(4,M))
      PSU(1)=PSU(1)+PME(4,M)
      PSU(2)=PSU(2)+PME(5,M)
      PSU(3)=PSU(3)+PME(6,M)
      ETM=SQRT(PME(4,M)**2+PME(5,M)**2+PME(6,M)**2+PME(9,M)**2)
      ESU=ESU+ETM
   14 CONTINUE
      ESU=SQRT(ESU**2-PSU(1)**2-PSU(2)**2-PSU(3)**2)
      IF(ABS(ESU-ES0).GT.0.100.OR.ABS(CSU-CS0).GT.0.1.
     *OR.ABS(SSU-SS0).GT.0.1.OR.ABS(BSU-BS0).GT.0.1)
     *LOWMIS=1
c      IF(LOWMIS.EQ.1)
c     *write(16,1000) ES0,CS0,SS0,BS0,ESU,CSU,SSU,BSU
 1000 FORMAT(1X,'ES0,CS0,SS0,BS0,ESU,CSU,SSU,BSU=',
     *F9.3,3(2X,F3.0)/33X,F9.3,3(2X,F3.0))
      IF(ABS(PSU(1)).GT.0.001.OR.ABS(PSU(2)).GT.0.001.OR.
     *ABS(PSU(3)).GT.0.001)
     *LOWMIS=1
c      IF(LOWMIS.EQ.1)
c     *write(16,1001) PSU,TDIAG
 1001 FORMAT(1X,'PSU=',3(F9.3,2X),5X,A8)
C     write(16,*) 'NP,IORHEI=', NP,IORHEI
      IF(IORHEI.NE.1)                                  GO  TO  112
      IF(NP.EQ.1.OR.(IABS(IIN(4))+IABS(IPN(4))).LE.0)  GO  TO  112
C
      NBA=0
      MB=MV+1
      DO  102  K=1,NP
      M=MV+K
      IF(IME(4,M).EQ.0)    GO  TO  102
      NBA=NBA+1
      DO  100  L=1,9
      TEMP=PME(L,MB)
      PME(L,MB)=PME(L,M)
  100 PME(L,M)=TEMP
      DO  101  L=1,5
      ITEMP=IME(L,MB)
      IME(L,MB)=IME(L,M)
  101 IME(L,M)=ITEMP
      TEMC=CLIDER(MB)
      CLIDER(MB)=CLIDER(M)
      CLIDER(M)=TEMC
      IDTE=IDPME(MB)
      IDPME(MB)=IDPME(M)
      IDPME(M)=IDTE
C      write(16,*) 'M,NBA,MB,IDPME(MB) =', M,NBA,MB,IDPME(MB)
      MB=MV+3
  102 CONTINUE
      IF(NP.NE.2)          GO  TO  105
      IF(NBA.EQ.2)         GO  TO  107
      DO 103 L=1,9
  103 PME(L,MV+3)=PME(L,MV+2)
      DO 104 L=1,5
  104 IME(L,MV+3)=IME(L,MV+2)
      CLIDER(MV+3)=CLIDER(MV+2)
      IDPME(MV+3)=IDPME(MV+2)
                        GO  TO  112
  105 IF(NBA.EQ.2)      GO  TO  107
      MLI=MV+1
      YL=YMEM(MLI)
      DO  106  N=2,NP
      YN=YMEM(MV+N)
      IF(YL.GT.YN)      GO  TO  106
      YL=YN
      MLI=MV+N
  106 CONTINUE
                        GO  TO  109
  107 Y3=YMEM(MV+3)
      Y1=YMEM(MV+1)
      IF(Y3.GT.Y1)      GO  TO  112
      MLI=MV+1
  109 DO  110  L=1,9
      TEMP=PME(L,MV+3)
      PME(L,MV+3)=PME(L,MLI)
  110 PME(L,MLI)=TEMP
      DO  111  L=1,5
      ITEMP=IME(L,MV+3)
      IME(L,MV+3)=IME(L,MLI)
  111 IME(L,MLI)=ITEMP
      TEMC=CLIDER(MV+3)
      CLIDER(MV+3)=CLIDER(MLI)
      CLIDER(MLI)=TEMC
      IDTE=IDPME(MV+3)
      IDPME(MV+3)=IDPME(MLI)
      IDPME(MLI)=IDTE
  112 CONTINUE
      IF(NCAS.LT.NCPRI)  RETURN
      write(16,599) NCAS,TDIAG,NP
  599 FORMAT(1X,'NCAS=',I6,' RESULTS FROM ',A8,' NP=',I5)
      DO  600  L=1,NP
      M=MV+L
      IF(NP.EQ.2.AND.L.EQ.2.AND.(IABS(IIN(4))+IABS(IPN(4))).GT.0)
     *M=M+1
      write(16,602) M,(PME(I,M),I=4,6),PME(9,M),(IME(J,M),J=1,5),
     * CLIDER(M),IDPME(M)
  602 FORMAT(1X,I3,4(1X,E11.4),4I2,I10,1X,F5.3,I5)
  600 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION YMEM(M)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      COMMON /MEMORYLAQ/ PME(9,5999),IME(5,5999)
      E=SQRT(PME(4,M)**2+PME(5,M)**2+PME(6,M)**2+PME(9,M)**2)
      YMEM=0.5*LOG((E+PME(6,M))/(E-PME(6,M)))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =============================================================
! LOWPEC removed by CMJ (XCP-3, LANL) on 09/08/2017
!    LOWPEC unused within LAQGSM code
! Purpose: UNKNOWN (no comments)
! =============================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      BLOCK DATA  BPNAME
C                -----
      CHARACTER*8  PNAME
      integer*4 identp
      COMMON /PNAME/PNAME(65),IDENTP(65)
      DATA PNAME/
     *  'PI+  ','PI-  ','K+   ','K-   ','K0   ','AK0  ','PI0  ',
     *  'ETA  ','ETAP ','RHO+ ','RHO- ','K*+  ','K*-  ','K*0  ',
     *  'AK*0 ','RHO0 ','OMEG ','PHI  ',
     *  'AP   ','AN   ','AS+  ','AS-  ','AS0  ','AXI- ','AXI0 ',
     *  'AL   ','ADL++','ADL+ ','ADL- ','ADL0 ','AS*+ ','AS*- ',
     *  'AS*0 ','AXI*-','AXI*0','AOM- ',
     *  'P    ','N    ','S+   ','S-   ','S0   ','XI-  ','XI0  ',
     *  'L    ','DL++ ','DL+  ','DL-  ','DL0  ','S*+  ','S*-  ',
     *  'S*0  ','XI*- ','XI*0 ','OM-  ','KS   ','KL   ',
     *  'GM   ','E-   ','E+   ','MU+  ','MU-  ','NUE  ','NUM  ',
     *  'ANUE ','ANUM '/
C
      DATA IDENTP/
     *     120 ,  -120 ,   130 ,  -130 ,   230 ,  -230 ,   110 ,
     *     220 ,   330 ,   121 ,  -121 ,   131 ,  -131 ,   231 ,
     *    -231 ,   111 ,   221 ,   331 ,
     *   -1120 , -1220 , -1130 , -2230 , -1230 , -2330 , -1330 ,
     *   -2130 , -1111 , -1121 , -2221 , -1221 , -1131 , -2231 ,
     *   -1231 , -2331 , -1331 , -3331 ,
     *    1120 ,  1220 ,  1130 ,  2230 ,  1230 ,  2330 ,  1330 ,
     *    2130 ,  1111 ,  1121 ,  2221 ,  1221 ,  1131 ,  2231 ,
     *    1231 ,  2331 ,  1331 ,  3331 ,20 ,   -20 ,
     *      10 ,    12 ,   -12 ,   -14 ,    14 ,    11 ,    13 ,
     *     -11 ,   -13 /
C     -----------------------------------------------------------------
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  ACTNAM(PIN,PN,MV,NP,NIN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL *8 PIN,PN,MA
      CHARACTER*8 PNA1,PNA2,PNAJ
      LOGICAL GH1H2,YESELA
      COMMON /INDDEC/INDDEC
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON /ITHEA/ ITHEA(11)
      COMMON /H1H2/ GH1H2(11)
      COMMON /YESELA/ YESELA
      COMMON /COMLID/ PLIDER(499)
      COMMON /IDN12/ ID1,ID2
      COMMON /IDN120/ ID10,ID20
      DIMENSION  PIN(9),PN(9),RS(4),PS(3)
C
   1  CONTINUE
      ID1=ID10
      ID2=ID20
      NPTCL=0
      DO  10  I=1,11
      ITHEA(I)=0
   10 GH1H2(I)=.TRUE.
      GH1H2(4)=YESELA
C
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      PX2= PN(4)
      PY2= PN(5)
      PZ2= PN(6)
      AM1=PIN(9)
      AM2= PN(9)
      ECM=SQRT((PIN(8)+PIN(9)+PN(8)+PN(9))**2-(PIN(4)+PN(4))**2-
     -(PIN(5)+PN(5))**2-(PIN(6)+PN(6))**2)
      CALL  PANUID(ID1,IK1,PNA1)
      CALL  PANUID(ID2,IK2,PNA2)
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      CALL  CROSEC(3,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIAN,0)
c     IF(SITO.LE.0..OR.SIEL.LE.0.)
c    *write(16,601) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2,
c    *             SITO,SIEL,SIEX,SIAN
 601  FORMAT(1X,'ACTNAM1: (',A8,')',4(F8.4,1X)/
     *8X,'+(',A8,')',4(F8.4,1X)/
     *8X,'SITO,SIEL,SIEX,SIAN=',4(E11.4,1X))
      IF(SITO.LE.0.)  SITO=1.
      IF(SIEL.LE.0.)  SIEL=SITO
C
c      IF(NCAS.GE.NCPRI)
c     *write(16,600) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2,
c     *             SITO,SIEL,SIEX,SIAN
c      IF(NCAS.GE.NCPRI)
c     *write( *,600) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2,
c     *             SITO,SIEL,SIEX,SIAN
 600  FORMAT(1X,'ACTNAM: (',A8,')',4(F8.4,1X)/
     *8X,'+(',A8,')',4(F8.4,1X)/
     *8X,'SITO,SIEL,SIEX,SIAN=',4(E11.4,1X))
      CALL  COLLHH(ID1,AM1,PX1,PY1,PZ1,ID2,AM2,PX2,PY2,PZ2,
     *SITO,SIAN,SIEL)
      NP=0
      NIN=0
C
      IF(INDDEC.EQ.0)  GO  TO 116
C
      IF(NPTCL.LE.0)   RETURN
      DO 115  ND=1,5
      IF(NPTCL.LE.0)   GO  TO  1
      NDEC=0
      NPTCL1=NPTCL
      DO 114  I=1,NPTCL1
      IDI=IDENT(I)
      IF(IDI.EQ.110.OR.IDI.EQ.230.OR.IDI.EQ.-230)
     *                 GO  TO  114
      IF(IDI.EQ.1230.OR.IDI.EQ.-1230)
     *                 GO  TO  114
      CALL DECAYQ(I,IPOINT)
      IF(IPOINT.LT.0)  GO  TO  114
      NDEC=NDEC+1
      DO  113  J=1,9
      PPTCL(J,I)=PPTCL(J,NPTCL)
  113 CONTINUE
      IDENT(I)=IDENT(NPTCL)
      IORIG(I)=IORIG(NPTCL)
      IDCAY(I)=IDCAY(NPTCL)
      PLIDER(I)=PLIDER(NPTCL)
      NPTCL=NPTCL-1
  114 CONTINUE
      IF(NDEC.EQ.0)     GO  TO  116
  115 CONTINUE
  116 IF(NPTCL.LE.0)  GO  TO 1
      DO  13  J=1,NPTCL
      MA=PPTCL(5,J)
      PS(1)=PPTCL(1,J)
      PS(2)=PPTCL(2,J)
      PS(3)=PPTCL(3,J)
      RS(1)=PPTCL(6,J)
      RS(2)=PPTCL(7,J)
      RS(3)=PPTCL(8,J)
      RS(4)=PPTCL(9,J)
      IDJ=IDENT(J)
      CALL  PANUID(IDJ,JP,PNAJ)
      CALL  RASPAN(IDJ,JP,MA,PS,RS,NP,MV)
   13 CONTINUE
      DO  14  I=1,11
      IF(GH1H2(I))   ITHEA(I)=1
   14 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  RASPAN(IDJ,JP,MA,P,R,NP,MV)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL *8 P,MA,R
      COMMON/IACT/ IACT
      COMMON /MEMORYLAQ/ PME(9,5999),IME(5,5999)
      COMMON /COMCHA/ ICHAM(18),ICHAB(18)
      COMMON /COMSTR/ ISTR(36)
      COMMON /COMLID/ PLIDER(499)
      COMMON /CSLID/ CLIDER(5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION  P(3),R(4)
      IF(JP.LT.1.OR.JP.GT.65)         GO  TO  12
      M=MV+NP+1
      IF(M.GT.5999)                   GO  TO  14
      CHARx=CHARGE(IDJ)*1.001d0
      IME(1,M)=INT(CHARx)
c      write(16,*) 'from RASPAN: IDJ,CHARx,IME(1,M)=', IDJ,CHARx,IME(1,M)
      IME(2,M)=0
      IF(JP.GE.55.AND.JP.LE.65) IME(2,M)=IDJ
      IME(3,M)=IS(IDJ)
      IME(4,M)=IB(IDJ)
      IF(JP.LE. 7)               IME(5,M)=0
      IF(JP.GE. 8.AND.JP.LE.18)  IME(5,M)=INTG(1000.  *TAUN(JP))
      IF(JP.GE.19.AND.JP.LE.26)  IME(5,M)=0
      IF(JP.GE.27.AND.JP.LE.35)  IME(5,M)=INTG(1000.  *TAUN(JP))
      IF(JP.GE.36.AND.JP.LE.44)  IME(5,M)=0
      IF(JP.GE.45.AND.JP.LE.53)  IME(5,M)=INTG(1000.  *TAUN(JP))
      IF(JP.GE.54)               IME(5,M)=0
   99 CONTINUE
      IF(IACT.GT.2)  GO  TO  98
      PME(1,M)=0.
      PME(2,M)=0.
      PME(3,M)=0.
      PME(7,M)=0.
      GO  TO  97
   98 PME(1,M)=R(1)
      PME(2,M)=R(2)
      PME(3,M)=R(3)
      PME(7,M)=R(4)
   97 PME(4,M)=P(1)
      PME(5,M)=P(2)
      PME(6,M)=P(3)
      PME(9,M)=MA
      NP=NP+1
      CLIDER(M)=PLIDER(NP)
      IDPME(M)=IDJ
      RETURN
   12 write(16,13) JP
   13 FORMAT(5X,'RASPAN: JP= ',I5)
      RETURN
   14 write(16,15)
   15 FORMAT(5X,'RASPAN: MV+NP+1>5999')
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE   PANUID(ID,N,PN)
      IMPLICIT INTEGER (I-N)
      CHARACTER*8 PN,PNAME
      COMMON /PNAME/PNAME(65),IDENTP(65)
      N=0
      DO  10  I=1,65
      IF(IDENTP(I).EQ.ID)  N=I
      IF(IDENTP(I).EQ.ID)  PN=PNAME(I)
   10 CONTINUE
      IF(N.EQ.0)   write(16,991) ID
      IF(N.EQ.0)   write( *,991) ID
  991 FORMAT(2X,'PANUID:  ID=',I5)
      RETURN
      ENTRY IDPANU(ID,N,PN)
      IF(N.LT.1.OR.N.GT.65)  write(16,992) N
  992 FORMAT(2X,'IDPANU:   N=',I5)
      ID=IDENTP(N)
      PN=PNAME(N)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PANUN(P,IP,N)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL *8 P,M
      CHARACTER*8 PN
      DIMENSION P(9),IP(5)
      N=0
      IC=IP(1)
      IL=IP(2)
      IS=IP(3)
      IB=IP(4)
      IR=IP(5)
      M =P(9)
      IF(IL.NE.0)  CALL  PANUID(IL,N,PN)
      IF(IP(2).NE.0)  RETURN
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0)     N=1
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0)     N=2
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.EQ.0)     N=3
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.EQ.0)     N=4
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.EQ.0)     N=5
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.EQ.0)     N=6
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND.
     *ABS(M-.140).LT.0.01)                                 N=7
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.
     *ABS(M-.549).LT.0.01)                                 N=8
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.
     *ABS(M-.958).LT.0.01)                                 N=9
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0)     N=10
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0)     N=11
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.NE.0)     N=12
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.NE.0)     N=13
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.NE.0)     N=14
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.NE.0)     N=15
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.
     *ABS(M-.770).LT.0.01)                                 N=16
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.
     *ABS(M-.783).LT.0.01)                                 N=17
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.
     *ABS(M-1.020).LT.0.01)                                N=18
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.EQ.0)    N=19
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.EQ.0)    N=20
      IF(IC.EQ.-1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=21
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=22
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=23
      IF(IC.EQ. 1.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.EQ.0)    N=24
      IF(IC.EQ. 0.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.EQ.0)    N=25
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=26
      IF(IC.EQ.-2.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=27
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=28
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=29
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=30
      IF(IC.EQ.-1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=31
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=32
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=33
      IF(IC.EQ. 1.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.NE.0)    N=34
      IF(IC.EQ. 0.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.NE.0)    N=35
      IF(IC.EQ. 1.AND.IS.EQ. 3.AND.IB.EQ.-1.AND.IR.EQ.0)    N=36
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.EQ.0)     N=37
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.EQ.0)     N=38
      IF(IC.EQ. 1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=39
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=40
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=41
      IF(IC.EQ.-1.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.EQ.0)     N=42
      IF(IC.EQ. 0.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.EQ.0)     N=43
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=44
      IF(IC.EQ. 2.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=45
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=46
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=47
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=48
      IF(IC.EQ. 1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=49
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=50
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=51
      IF(IC.EQ.-1.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.NE.0)     N=52
      IF(IC.EQ. 0.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.NE.0)     N=53
      IF(IC.EQ.-1.AND.IS.EQ.-3.AND.IB.EQ.1.AND.IR.EQ.0)     N=54
C
C  VARIABLE MASS NONSRTANGE RESONANSES ===>RO0
c      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.N.EQ.0)
c     *write(16,1001) M,IC
1001  FORMAT(1X, 'PANUN: GETTING RHO WITH MASS=',F7.3,', CHARGE=',I3)
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.N.EQ.0)
     *                                                      N=16
C
C  VARIABLE MASS NONSTRANGE PARTICLE    ===>PI0
c      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND.N.EQ.0)
c     *write(16,1002) M,IC
1002  FORMAT(1X, 'PANUN: GETTING PI0 WITH MASS=',F7.3,', CHARGE=',I3)
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND.N.EQ.0)
     *                                                      N= 7
C
      IF(N.EQ.0)   write(16,1000) IP,M
 1000 FORMAT(2X,'PANUN: IP=',5I5,2X,'M=',E13.6)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION TAUN(IK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18),
     *PGAMB(18),MESO(9,2)
      IF(IK.LE.18)   GO  TO  2
      IF(IK.GE.37)   GO  TO  1
      IF(PGAMB(IK-18).EQ.0.)  GO  TO  3
      TAU0=0.197/PGAMB(IK-18)
      GO  TO  4
    1 IF(PGAMB(IK-36).EQ.0.)  GO  TO  3
      TAU0=0.197/PGAMB(IK-36)
      GO  TO  4
    2 IF(PGAMM(IK).EQ.0.)     GO  TO  3
      TAU0=0.197/PGAMM(IK)
      GO  TO  4
    3 TAUN=10.0E+9
      RETURN
    4 DRND=RNDM(-1.)
      TAUN=-TAU0*LOG(DRND)
      IF(TAUN.LE.0.001)     TAUN=0.0011
      IF(TAUN.GT.1.0D+3)    TAUN=1.0D+3
      RETURN
      END
c************ last correction 11.10.22 16:57(elasle)******

! =============================================================
! HHQGSE removed by CMJ (XCP-3, LANL) on 09/08/2017
!    HHQGSE unused within LAQGSM code
! Purpose: UNKNOWN (no comments)
! =============================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE TYPRE(IK1,AM1,IK2,AM2,
     *SQS,IEL,IDIF,IUNCY,IPLAN,IBINA,SIGTOT,SIGEL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C     COMPUTE REACTION CHANNEL AT LOW ENERGY HADRON COLLISION
C
C
C     PRODUCED BY DR.N.S.AMELIN FROM LCTA OF JOINT INSTITUTE FOR NUCLEAR
C     RESEARCH.  1986 -VERSION (QUARK-GLUON STRINGS MODEL USED)
C     FOR HADRON HADRON COLLISIONS SIMULATION BY MONTE-CARLO METHOD
C
C               WILL BE ONLY,IF
C     IEL=1-ELASTIC SCATTERING
C     IDIF=1-DIFFRACTIVE SCATTERING
C     IUNCY=1-UNDEVELOPED CYLINDER SCATTERING
C     IPLAN=1-PLANAR SCATTERING
C     IBINA=1-TWO PARTICLE REACTION
C     IBINA=IEL=IDIF=IUNCY=IPLAN=0-CYLINDER SCATTERING.
C
C         PROJECTILE PARAMETERS
C     PX1-X MOMENTUM COMPONENT(GEV/C)
C     PY1-Y MOMENTUM COMPONENT(GEV/C)
C     PZ1-Z MOMENTUM COMPONENT(GEV/C)
C     IK1-PARTICLE TABLE NUMBER
C
C         TARGET PARAMETERS
C     PX2-X...
C     PY2-Y...
C     PZ2-Z...
C     IK2- PARTICLE...
C
C   FOR EXAMPLE PP-COLLISIONS AT 5GEV/C IN LABORATORY FRAME
C     IK1=37,IK2=37,PX1=PX2=PY1=PY2=PZ2=0.,PZ1=5.
C
C  PARTICLE TABLE NUMBER IK            PARTICLE LABEL
C           1                                  PI+
C           2                                  PI-
C           3                                  K+
C           4                                  K-
C           5                                  K0
C           6                                  AK0
C           7                                  PI0
C           8                                  ETA
C           9                                  ETAP
C          10                                  RHO+
C          11                                  RHO-
C          12                                  K*+
C          13                                  K*-
C          14                                  K*0
C          15                                  AK*0
C          16                                  RH0
C          17                                  OMEGA
C          18                                  PHI0
C          37                                  P
C          38                                  N
C          39                                  S+
C          40                                  S-
C          41                                  S0
C          42                                  KSI-
C          43                                  KSI0
C          44                                  L
C          45                                  DL++
C          46                                  DL+
C          47                                  DL-
C          48                                  DL0
C          49                                  S*+
C          50                                  S*-
C          51                                  S*0
C          52                                  KSI*-
C          53                                  KSI*0
C          54                                  OM-
C
C    ESSENTIAL PARAMETERS FOR STRING OR CLUSTER DECAY
C     PUD-STRANGE PARTICLE PRODUCTION PROBABILITY
C     PS1-VECTOR MESON PRODUCTION PROBABILITY
C     SIGMA-CONSTANT TO COMPUTE TRANSFERSE MOMENTUM OF PARTICLES
C     PARM,PARB,SWMAX ARE CUTS FOR STRING DECAY
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /DATA2/ PUD,PS1,SIGMA,CX2
      COMMON /SIGDIA/ CROSD(5),DST
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON /PRINTS/ IPRINT
      COMMON /CPRSIG/ ISGCOL
      COMMON /P0LAB1/ P0,DSM,TKIN
      COMMON /YESELA/YESELA
      COMMON/STREXC/STREXC
      LOGICAL IPRINT
      LOGICAL YESELA
C     REAL *8 LAB1,LAB2
      CHARACTER*8 LAB1,LAB2
      CALL LABEL(LAB1,IDIN(1))
      CALL LABEL(LAB2,IDIN(2))
      IEL=0
      IDIF=0
      IUNCY=0
      IPLAN=0
      IBINA=0
16    SIGIN=SIGTOT-SIGEL
      IB1=IBLE(IK1)
      IB2=IBLE(IK2)
      SBIN   = CROSD(1)
      SIGDIF = CROSD(2)
      SUNC   = CROSD(3)
      SIPLAN = CROSD(4)
      SCYLIN = CROSD(5)
      S=SQS**2
      PS1=0.75
      IF(.NOT.(IB1.EQ.0.AND.IB2.EQ.1)) GO TO 7655
      IF(SQS-AM1-AM2.GT.0.30) GO TO 7655
C  FOR M-N REACTION
      SIGDIF=0.
      SCYLIN=0.
      SUNC=0.
* ---V.T.  05.11.94 renormalization for strangeness exchange
      IS1=ISLE(IK1)
      IS2=ISLE(IK2)
      ICHA=IQLE(IK1)+IQLE(IK2)
      if(IS1 == 0.and.IS2 == 0.and.
     &  (ICHA > -1.and.ICHA < 3) ) then
        SIGTOT=SIGTOT-SBIN
        SIGIN=SIGIN-SBIN
        SBIN=SBIN-STREXC
        if(SBIN < 0.d0) SBIN=0.
        SIGTOT=SIGTOT+SBIN
        SIGIN=SIGIN+SBIN
      endif
* -----
      SIGEL=SIGTOT-SBIN-SIPLAN
      GO TO 7654
C      DELTA+NUCLEON
7655  CONTINUE
c     IF(SQS-AM1-AM2.GT.0.15) GO TO 7654
      IF(SQS.GT.3.00) GO TO 7654
      IF(.NOT.(IK2.EQ.37.OR.IK2.EQ.38)) GO TO 7654
      IF(.NOT.(IK1.EQ.46.OR.IK1.EQ.45.OR.IK1.EQ.47.OR.
     *IK1.EQ.48)) GO TO 7654
      ICHA=IQLE(IK1)+IQLE(IK2)
      IF(ICHA.EQ.3.OR.ICHA.EQ.-1) GO TO 7654
      SBIN=CRNDNN(S,ICHA)
      SIGDIF=0.
      SCYLIN=0.
      SUNC=0.
      SIPLAN=0.
      SIGEL=SIGTOT-SBIN
7654  CONTINUE
       ISGCOL=0
      IF((ISGCOL.NE.0).OR.(.NOT.IPRINT)) GO TO 525
      WRITE(ITLIS,1000) LAB1,LAB2,S,P0
1000  FORMAT(/15X,A8,A8,18H COLLISION AT SCM=,F10.4,7H GEV**2,
     *'( PLAB OF ',F10.4)
      IF(SIGEL.GT.0) WRITE(ITLIS,1001)  SIGEL
      IF(SBIN.GT.0) WRITE(ITLIS,1002)   SBIN
      IF(SIGDIF.GT.0) WRITE(ITLIS,1003) SIGDIF
      IF(SUNC  .GT.0) WRITE(ITLIS,1004) SUNC
      IF(SIPLAN.GT.0) WRITE(ITLIS,1005) SIPLAN
      IF(SCYLIN.GT.0) WRITE(ITLIS,1006) SCYLIN
1001  FORMAT(3X,'  ELASTIC  REACTION CROSS SECTION =',F10.4,'MB')
1002  FORMAT(3X,'    BINAR  DIAGRAM  CROSS SECTION =',F10.4,'MB')
1003  FORMAT(3X,'DIFRACTION DIAGRAM  CROSS SECTION =',F10.4,'MB')
1004  FORMAT(3X,'UNDEV.CYL. DIAGRAM  CROSS SECTION =',F10.4,'MB')
1005  FORMAT(3X,' PLANAR    DIAGRAM  CROSS SECTION =',F10.4,'MB')
1006  FORMAT(3X,' CYLINDR   DIAGRAM  CROSS SECTION =',F10.4,'MB')
 525  RND=RNDM(-1.)
      IF(.NOT.YESELA)         GO TO 526
      IF(RND.GT.SIGEL/SIGTOT) GO TO 526
      IEL=1
      GO TO 10
  526 IS1=ISLE(IK1)
      IS2=ISLE(IK2)
      IF (IS1.NE.0.OR.IS2.NE.0)        GO  TO  6
      IF (IK1.LE.36.OR.IK2.LE.36)      PS1=0.50
    6 IF (SIGIN.GT.0.001)              GO  TO  7
      IEL=1
                                       GO  TO  10
    7 RND=RNDM(-1.)
      IF(RND.GT.SBIN/SIGIN)            GO  TO  77
      IBINA=1
                                       GO  TO  10
   77 IF(RND.GT.(SIPLAN+SBIN)/SIGIN)   GO  TO  8
      IPLAN=1
                                       GO  TO  10
    8 IF(RND.GT.(SIPLAN+SBIN+SIGDIF)/SIGIN) GO  TO  9
      IDIF=1
                                       GO  TO  10
    9 IF(RND.GT.(SIPLAN+SIGDIF+SBIN+SUNC)/SIGIN)     GO  TO  10
      IUNCY=1
10    CONTINUE
      ISGCOL=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SLOPE(B,BFOR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   SLOPE EVALUATES THE ELASTIC SLOPES B AND BFOR, WHERE THE
C CROSS-SECTION VARIES AS EXP(BFOR*TFORWARD)*EXP(B*(T-TFORWARD))
C   *** FOR THE CHOICE OF ISL, ISLFOR
C
      COMMON /CUTOF2/BA,BB,BINF,POW,ISL,ISLFOR
      COMMON /COMKI1/ HLA2,HLB2,W,INUMA
      COMMON /COMKI2/ELA,ELB,PLALB
      COMMON /CALC/HA,HB,HA2,HB2
      DIMENSION  COPE(3),SLFOR(2)
      S=W**2
      HLA=SQRT(HLA2)
      HLB=SQRT(HLB2)
      FACT=ALAMB(1.0D0,HLA2/S,HLB2/S)
      COPE(1)=0.5*FACT*(BA*(HLA/HA)**POW + BB*(HLB/HB)**POW)
      COPE(2)=0.5*FACT*(BA+BB)
      DIFFA=(BA-BINF)*(HLA/HA)**POW
      DIFFB=(BB-BINF)*(HLB/HB)**POW
      COPE(3)=0.5*FACT*(BINF+DIFFA+BINF+DIFFB)
      B0=COPE(ISL)
      SLFOR(1)=0.0
      SLFOR(2)=B0
      BFOR=SLFOR(ISLFOR)
      IF(INUMA.EQ.0.OR.INUMA.EQ.2) B0=B0*0.8
C----  FOR MESON-MESON COLLISION
      TLAB=SQRT(PLALB**2+HLA2)-HLA
      IF(INUMA.EQ.2.AND.TLAB.LT.2.4)  B0=2.0
      B=B0+0.7*LOG(PLALB)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ANG(TFOR,TBACK,T,Z,PHI)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   ANG CALCULATES (RANDOMLY) THE POLAR AND AZIMUTHAL SCATTERING ANGLES
C (Z=COS(THETA) AND PHI) FOR THE TWO-BODY PROCESS, USING THE ELASTIC
C   SLOPES B AND BFOR
C
C   TFOR=MOMENTUM TRANSFER FOR FORWARD SCATTERING INTO A AND B
C   TBACK=MOMENTUM TRANSFER FOR BACKWARD SCATTERING INTO A AND B
C
      COMMON /COMKI1/ HLA2,HLB2,W,INUMA
      COMMON /COMKI2/ELA,ELB,PLALB
      COMMON /CALC/HA,HB,HA2,HB2
      COMMON /BEES/B,BFOR
      S=W**2
      EA=(S+HA2-HB2)/(2.0*W)
      EB=(S+HB2-HA2)/(2.0*W)
      PAB2=EA**2-HA2
      PAB=SQRT(PAB2)
      TFOR=HA2+HLA2-2.0*(EA*ELA-PAB*PLALB)
      TBACK=HA2+HLA2-2.0*(EA*ELA+PAB*PLALB)
      TDIFF=TBACK-TFOR
      TB=B*TDIFF
      IF(TB.GT.-30.) GO TO 5
      TBB=0.
      GO TO 6
    5 TBB=EXP(TB)
    6 R3=RNDM(-1.)
      ZOT=TBB*(1.0-R3) + R3
      TPRIME=LOG(ZOT)
      TPRIME=TPRIME/B
      T=TPRIME+TFOR
      TRAT=TPRIME/TDIFF
      Z=1.0-2.0*TRAT
      IF(Z.LT.-1.0) Z=-1.0
      IF(Z.GT.1.0) Z=1.0
      R4=RNDM(-1.)
      PHI=twpi*R4
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ELZPLE(IK1,IK2,TKIN,Z,PHI,IEXE)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      PHI=twpi*RNDM(-1.)
      IF(IEXE.EQ.0) GO TO 1
       Z=COSP(TKIN,12)
       RETURN
 1    CALL MARLE(IK1,IK2,KS)
      IBP=IBLE(IK1)
      Z=COSAM(IBP,TKIN,KS)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE MARLE(IK01,IK02,KS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      IK1=IK01
      IK2=IK02
      IB1=IBLE(IK1)
      IF(IK2.EQ.38) GO TO 112
      IK2=37
112   IB2=IBLE(IK2)
      IQ1=IQLE(IK1)
      IQ2=IQLE(IK2)
      MQ=IQ1+IQ2
      IF(IB1+IB2.LE.1) GO TO 1
C   NUCLEON-NUCLEON COLLISION
      IF(MQ-1) 3,4,3
C  MESON-NUCLEON COLLISION
 1    IF(MQ.EQ.2.OR.MQ.EQ.-1) GO TO 3
      IF(MQ.EQ.0) GO TO 2
      IF(IQ1-1) 5,4,5
 2    IF(IQ1+1) 5,4,5
 3    KS=1
      RETURN
 4    KS=2
      RETURN
 5    KS=3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSAM(IB1,TKIN,KS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      IF(IB1.EQ.1) GO TO 4
C  MESON-NUCLEON SCATTERING
      IF(KS-2) 1,2,3
C  POSETIVE MESON+PROTON
 1    COSAM=COSP(TKIN,4)
      RETURN
C  NEGATIVE MESON+PROTON
 2    COSAM=COSP(TKIN,8)
      RETURN
C   NEUTRAL MESON+PROTON
 3    IF(RNDM(-1.)-0.5) 1,1,2
C  NUCLEON-NUCLEON SCATTERING
 4    IF(KS.EQ.1) GO TO 5
C  NEUTRON-PROTON SCATTERING
      IF(TKIN.GT.0.97) GO TO 6
      COSAM=CMJ(TKIN,3)
      RETURN
C  PROTON-PROTON SCATTERING
 5    IF(TKIN.GT.0.46) GO TO 6
      COSAM=1.-2.0*RNDM(-1.)
      RETURN
C  NEUTRON-PROTON AND PROTON-PROTON SCATTERING
 6    IF(TKIN.GT.2.8) GO TO 7
      COSAM=(1.0+CMJ(TKIN,1))/4.0
      RETURN
 7    COSAM=(3.0+CMJ(TKIN,2))/4.0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSP(T,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF SCATTERED ANGLE
      JJ=J                       ! 10.02.2002
      IF(T.LE.0.08) GO TO 1
      IF(T.LE.0.3) GO TO 2
      IF(T.LE.1.) GO TO 3
      IF(J.EQ.12) JJ=8
      COSP=CMJ(T,JJ+3)
      RETURN
 1    COSP=CMJ(T,J)
      RETURN
 2    COSP=CMJ(T,J+1)
      RETURN
 3    IF(J.EQ.12) JJ=8
      COSP=CMJ(T,JJ+2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION CMJ(T,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF COSINE BY MEANS OF THE TABLE COEFFICIENTS
      COMMON /COEF3/ ANKJ(4,4,13)
      DIMENSION C(4,4)
      S1=0.
      S2=0.
      R=RNDM(-1.)
      R1=SQRT(R)
      DO 2 K=1,4
      DO 2 N=1,4
      C(N,K)=ANKJ(N,K,J)
 2    CONTINUE
      DO 3 N=1,4
      DO 3 K=1,4
      S1=S1+C(N,K)*(R**(N-1))*(T**(K-1))
 3    S2=S2+C(N,K)*(T**(K-1))
      IF(S2.GE.1.)  S2=0.999999
      D=2.0*R1*(S1+(1.-S2)*R**4)-1.
      IF(ABS(D).GT.1.) GO TO 4
      CMJ=D
      RETURN
 4    CMJ=SIGN(1.D0,D)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IQLE(IK)
      IMPLICIT  INTEGER (I-N)
C  COMPUTE PARTICLE CHARGE
      COMMON /COMCHA/ ICHAM(18),ICHAB(18)
      IF(IK.GT.36) GO TO 1
C   MESON CHARGE
      IQLE=ICHAM(IK)
      RETURN
C  BARION CHARGE
 1    IQLE=ICHAB(IK-36)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION GAM(IK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18),
     *PGAMB(18),MESO(9,2)
      IF(IK-36) 1,1,2
 1    GAM=PGAMM(IK)
      RETURN
 2    GAM=PGAMB(IK-36)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION XDIST(XMIN,IB1,IS1,IFL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON /COMABM/ ALFAM,BETAM
      COMMON /COMABB/ ALFAB,BETAB
      IF(IABS(IB1).EQ.1) GO TO 3
      IF(IS1.NE.0) GO TO 33
 2    CALL SBETA(X,ALFAM,BETAM)
      IF(X.LE.XMIN) GO TO 2
      IF(RNDM(-1.).GE.0.5) X=1.-X
      XDIST=X
      RETURN
 33   BETAN=1.
      GO TO 1
 3    BETAN=BETAB
      IF(IABS(IFL).GT.1) BETAN=BETAB+1
 1    CALL SBETA(X,ALFAB,BETAN)
      IF(X.LE.XMIN) GO TO 1
      XDIST=X
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SBETA(X,ALFA,BETA)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  SIMULATION OF BETA DISTRIBUTION U(X)=C*X**(ALFA-1)*(1-X)**(BETA-1)
C  IONK,S METHOD
 1    RAN1=RNDM(-1.)
      RAN2=RNDM(-1.)
      R1A=RAN1**(1./ALFA)
      R2B=RAN2**(1./BETA)
      R12=R1A+R2B
      IF(R12.GE.1.) GO TO 1
      X=R1A/R12
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION WMX(KS,IF1,IF2,PARM,PARB)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON /DATA4/QMAS(9)
      IF(KS)1,1,2
 1    WMX=PARM+QMAS(IABS(IF1))+QMAS(IABS(IF2))
      RETURN
 2    WMX=PARB+QMAS(IABS(IF1))+QMAS(IABS(IF2))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION AMASF(IK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  COMPUTE FIXED PARTICLE MASS
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18),
     *PGAMB(18),MESO(9,2)
      IF(IK-36) 1,1,2
 1    AMASF=PMASM(IK)
      RETURN
 2    AMASF=PMASB(IK-36)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION AMAS(IK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  COMPUTE PARTICLE MASS
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18),
     *PGAMB(18),MESO(9,2)
      COMMON /ISOB3/ISOB3
      COMMON /INTTYP/ITYP
      COMMON /ITHEA/ITHEA(11)
      if ((ik - 36) <= 0) then
         AMAS=PMASM(IK)
         if (( IK /= 10 .AND. IK /= 11 .AND. IK /= 16) .OR.
     *(isob3 /= 1)) then
            return
         end if
         CALL MRHO(AMR,GD)
         AMAS=AMR
         RETURN
      else
         AMAS=PMASB(IK-36)
         IF(IK < 45 .OR. 48 < IK) RETURN
         IF(ISOB3.NE.1)  RETURN
      end if
      CALL  MDELTA(AMD,GD)
      AMAS=AMD
      RETURN
      END
C   ******************************************************************
        SUBROUTINE MRHO(AMRHO,GD)
C *** SIMULATION OF RHO_MESON MASS DISTRIBUTION  ***
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
        DATA AMPI/0.140/,AMR0/0.770/
        DATA AMRMIN/0.281/,AMRMAX/1.200/
   10 CONTINUE
        AMR=AMRMIN+RNDM(-1.)*(AMRMAX-AMRMIN)
        Q=0.5*SQRT(AMR*AMR-4.*AMPI*AMPI)
        GD=0.095 *q*((q/AMPI)/(1.+(q/AMR0)**2))**2
*ti       GD=0.150 *
*ti  *  AMR/AMR0*((1.-(2.*AMPI/AMR)**2)/(1.-(2.*AMPI/AMR0)**2))**1.5
*ti     F=    (GD*AMR0)**2/((AMR**2 -AMR0**2)**2+(AMR0*GD)**2)
        F=0.25*GD**2/((AMR-AMR0)**2+0.25*GD**2)
      IF(RNDM(-1.).GT.F)  GO  TO  10
c       write( 6,*) 'MR,GD,q=', AMR,GD,Q
        AMRHO=AMR
        RETURN
        END
C   ******************************************************************
        SUBROUTINE MDELT1(AMDEL,GD)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/IDN12/ ID1,ID2
      COMMON/SCOLLD/ SCOLLD
        DATA AMN/0.940/,AMPI/0.140/,AMD0/1.232/
        DATA AMDMIN/1.081/,AMDMAX/1.700/
        AMD02=AMD0*AMD0
        IBB=IB(ID1)+IB(ID2)
        IF(IBB.EQ.0)  write( 6,*) 'MDELT1: IBB', IBB
        AMX=AMN
        IF(IBB.EQ.1)  AMX=AMPI
        S=SCOLLD+0.03
        SSX=SQRT(S)-AMX
        AMDMAXX=AMDMAX
        IF(SSX.LT.AMDMAX)  AMDMAXX=SSX
        IF(AMDMAXX.LT.AMDMIN) WRITE( 6,*) 'SSX,ID1,ID2=', SSX,ID1,ID2
        IF(AMDMAXX.LT.AMDMIN) GO TO 13
        L=0
   10 CONTINUE
        LI=0
   11 CONTINUE
        AMD=AMDMIN+RNDM(-1.)*(AMDMAXX-AMDMIN)
        GD=GDM(AMD)
        GD2=GD*GD
        X=SQRT(1.-1.25*GD2/AMD02)
        XM=AMD0*SQRT(0.6+0.4*X)
        IF(XM.GT.AMDMAXX)  XM=AMDMAXX
        F1=((AMD02-XM*XM)**2+AMD02*GD2)*AMD /
     1    (((AMD02-AMD*AMD)**2+AMD02*GD2)*XM )
      IF(RNDM(-1.).GT.F1)           GO  TO  11
        AX=(AMD+AMX)**2
        IF(AX.GE.S)                 LI=LI+1
        IF(AX.GE.S.AND.LI.LE.99)    GO  TO  11
        IF(LI.GT.99)                GO  TO  12
        XS=((S-AX)*(S-(AMD-AMX)**2)) /
     1   ((S-(AMDMIN+AMX)**2)*(S-(AMDMIN-AMX)**2) )
        F2=SQRT(XS)
        L=L+1
        IF(L.GT.50) WRITE(16,*) 'MDELT1: L,S =', L,S
        IF(L.GT.50) GO TO 2
      IF(RNDM(-1.).GT.F2)  GO  TO  10
    2   AMDEL=AMD
        RETURN
   12 write( 6,*) 'MDELT1:S,AMD,SSX,ID1,ID2=', S,AMD,SSX, ID1,ID2
        AMDEL=AMD
        RETURN
   13   AMDEL=AMDMIN
        RETURN
        END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE KSPIN(IFL1,KS1,IK1,IB1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C  CALCULATION OF KS1
C  KS1=0 FOR MESON,KS1=1 FOR DIQUARK WITH S=0,
C  KS1=2 FOR DIQUARK WITH S=1
      COMMON /DATA6/POD810(3,6),PODSA(8,3),KBAR(18,2)
      IF(IB1.EQ.0) KS1=0
      IF(IB1.EQ.1.AND.IK1.GE.45) KS1=2
      IF(IB1.EQ.1.AND.IK1.LT.45)
     *KS1=1+INT(PODSA(IK1-36,IFL1)+RNDM(-1.))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FLAVD(IFLQ,IFLJ,IFLL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C  COMPUTE QUARK FLAVOUR IN DIQUARK
c
      KD=IFLQ-3
      GO TO (61,62,63,64,65,66),KD
 61   IFLJ=1
      IFLL=1
      RETURN
 62   IFLJ=2
      IFLL=2
      RETURN
 63   IFLJ=1
      IFLL=2
      IF(RNDM(-1.).GT.0.5D0) RETURN
      IFLJ=2
      IFLL=1
      RETURN
 64   IFLJ=1
      IFLL=3
      IF(RNDM(-1.).GT.0.5D0) RETURN
      IFLJ=3
      IFLL=1
      RETURN
 65   IFLJ=2
      IFLL=3
      IF(RNDM(-1.).GT.0.5D0) RETURN
      IFLJ=3
      IFLL=2
      RETURN
 66   IFLJ=3
      IFLL=3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION KI2(IFL01,IFL02,KSA,IR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C   TABLE NUMBER FOR MESON OR BARION DETERMINATION
      COMMON /ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /DATA6/POD810(3,6),PODSA(8,3),KBAR(18,2)
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18),
     *PGAMB(18),MESO(9,2)
      IFL1=IFL01
      IFL2=IFL02
      IF(KSA.NE.0) GO TO 2
C  TABLE NUMBER FOR MESON
      IFLSGN=(10-IFL1)/5
      IFL1=IABS(IFL1)
      IFL2=IABS(IFL2)
      IFL12=3*(IFL1-1)+IFL2
      K1=MESO(IFL12,IFLSGN)
      ISPIN=INT(PS1+RNDM(-1.))
      IF(IR.EQ.1) ISPIN=1
      IF(IR.EQ.2) ISPIN=0
      K2=9*ISPIN+K1
      IF(K1.LE.6) GO TO 1
      TMIX=RNDM(-1.)
      KM=K1-6+3*ISPIN
      K2=7+9*ISPIN+INT(TMIX+CMIX(KM,1))+INT(TMIX+CMIX(KM,2))
 1    KI2=K2
      RETURN
C  TABLE NUMBER FOR BARION
 2    IF(KSA.EQ.1)K=1
      IFLQ=IFL1
      IFQQ=IFL2
      IF(IFL2.GT.3) GO TO 21
      IFLQ=IFL2
      IFQQ=IFL1
 21   CONTINUE
      IF(KSA.EQ.2) K=1+INT(POD810(IFLQ,IFQQ-3)+RNDM(-1.))
      IF(IR.EQ.1) K=2
      IF(IR.EQ.2) K=1
 10   CONTINUE
      IFQ12=3*(IFQQ-4)+IFLQ
      K2=KBAR(IFQ12,K)
      IF(K2.NE.0) GO TO 11
      K=2
      GO TO 10
 11   CONTINUE
      IF(K2.EQ.41.OR.K2.EQ.44)GO TO 3
      KI2=K2
      RETURN
 3    IF(KSA.EQ.1) GO TO 6
      IF(IFQQ.EQ.6)GO TO 4
      K2=41+3*INT(0.75+RNDM(-1.))
      GO TO 5
 4    K2=41
 5    KI2=K2
      RETURN
 6    IF(IFQQ.EQ.6)GO TO 7
      K2=41+3*INT(0.25+RNDM(-1.))
      GO TO 8
 7    K2=44
 8    KI2=K2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! PTDIST removed by CMJ (XCP-3, LANL) on 09/07/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: COMPUTE NEW PRODUCTED QUARK TRANSFERSE MOMENTUM SIMULATION
!          OF PT DISTRIBUTION FROM U(PT)=SQRT(1./PI*SIGMA**2)
!          *EXP(-PT**2/(SIGMA**2))
! =====================================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTAM(CT,ST,CFI,SFI,PX1,PX2,J)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  ROTATE OF VECTOR PX1
      DIMENSION ROT(3,3),PX1(3),PX2(3)
      ROT(1,1)=CT*CFI
      ROT(1,2)=-SFI
      ROT(1,3)=ST*CFI
      ROT(2,1)=CT*SFI
      ROT(2,2)=CFI
      ROT(2,3)=ST*SFI
      ROT(3,1)=-ST
      ROT(3,2)=0.
      ROT(3,3)=CT
      IF(J.EQ.0) GO TO 2
      DO 1 I=1,3
 1    PX2(I)=ROT(I,1)*PX1(1)+ROT(I,2)*PX1(2)+ROT(I,3)*PX1(3)
      RETURN
 2    DO 3 I=1,3
 3    PX2(I)=ROT(1,I)*PX1(1)+ROT(2,I)*PX1(2)+ROT(3,I)*PX1(3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ANGLE(PX,CT,ST,CFI,SFI)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  COMPUTE ROTOR PARAMETERS
C  AND COSINE,SINE OF THETA AND PHI
      DIMENSION PX(3)
      P=SQRT(SPQ(PX,PX))
      CT=PX(3)/P
	if(CT*CT > 1.D0)  CT=SIGN(1.D0,CT)   ! 23.12.01
      ST=SQRT(1.D0-CT*CT)
      PM=P*ST
      IF(PM.EQ.0.) GO TO 1
      CFI=PX(1)/PM
      SFI=PX(2)/PM
      RETURN
 1    CFI=1.
      SFI=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORLLE(V,PX,E,L)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  LORENTZ TRANSFORMATION OF PX MOMENTUM COMPONENTS
      DIMENSION V(3),PX(3)
      VV=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
      GA=1.D0/SQRT(ABS(1.D0-VV))
      BEP=SPQ(V,PX)
      GABEP=GA*(GA*BEP/(1.+GA)-L*E)
      DO 1 I=1,3
 1    PX(I)=PX(I)+GABEP*V(I)
        RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SPQ(A,B)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     SCALAR PRODUCT OF THREE DIMENSIONAL VEKTORS
      DIMENSION A(3),B(3)
      SPQ=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE KINEM(PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,V,S,P1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF CM VELOCITY V(I),CM TOTAL ENERGY S
C  AND MOMENTUM COMPONENTS IN CM FRAME
      DIMENSION V(3),P1(3),P2(3)
      P1(1)=PX1
      P1(2)=PY1
      P1(3)=PZ1
      P2(1)=PX2
      P2(2)=PY2
      P2(3)=PZ2
      E1=SQRT(AM1**2+SPQ(P1,P1))
      E2=SQRT(AM2**2+SPQ(P2,P2))
      E=E1+E2
      S=AM1**2+AM2**2+2.*E1*E2-2.*SPQ(P1,P2)
      V(1)=(P1(1)+P2(1))/E
      V(2)=(P1(2)+P2(2))/E
      V(3)=(P1(3)+P2(3))/E
      CALL LORLLE(V,P1,E1,1)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FLAVO(IB1,IK1,IB2,IK2,IFL1,IFL2,IFL3,IFL4)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      COMMON /DATA5/IFLM(18,2),IFLB(18,3),DQQ(3,3)
C  FLAVOUR OF INTERACTING QUARK FOR BARYON AND MESON
      IF(IB1) 2,2,3
C  FLAVOURS OF MESON QUARKS
 2    INR1=INT(1.+2.*RNDM(-1.))
      IFL1=IFLM(IK1,INR1)
      I1=1
      IF(INR1.EQ.1)I1=2
      IFL2=IFLM(IK1,I1)
      IF(IK1.EQ.7.OR.IK1.EQ.16.OR.IK1.EQ.17) GO TO 13
      IF(IK1.EQ.8.OR.IK1.EQ.9)GOTO 12
      GO TO 4
 12   IFL1=2+INT(0.25+RNDM(-1.))
      IF(IFL1.EQ.2) IFL1=1+INT(0.5+RNDM(-1.))
      IFL1=IFL1*(-1.)**(INR1-1)
      IFL2=-IFL1
      GO TO 4
 13   IFL1=(1+INT(0.5+RNDM(-1.)))*(-1.)**(INR1-1)
      IFL2=-IFL1
      GO TO 4
C  FLAVOURS OF BARION QUARKS
 3    INR1=INT(1.+3.*RNDM(-1.))
      IFL1=IFLB(IK1-36,INR1)
      IF(INR1-2)200,201,202
 200  I1=2
      I2=3
      GO TO 205
 201  I1=1
      I2=3
      GO TO 205
 202  I1=1
      I2=2
 205  IF1=IFLB(IK1-36,I1)
      IF2=IFLB(IK1-36,I2)
      IFL2=DQQ(IF1,IF2)
C  FLAVOURS OF TARGET HADRON QUARKS
 4    IF(IB2) 20,20,30
 20   INR2=INT(1.+2.*RNDM(-1.))
      IFL3=IFLM(IK2,INR2)
      I2=1
      IF(INR2.EQ.1) I2=2
      IFL4=IFLM(IK2,I2)
      IF(IK2.EQ.7.OR.IK2.EQ.16.OR.IK2.EQ.17) GO TO 33
      IF(IK2.EQ.8.OR.IK2.EQ.9) GO TO 32
      RETURN
 32   IFL3=2+INT(0.25+RNDM(-1.))
      IF(IFL3.EQ.2) IFL3=1+INT(0.5+RNDM(-1.))
      IFL3=IFL3*(-1.)**(INR2-1)
      IFL4=-IFL3
      RETURN
 33   IFL3=(1+INT(0.5+RNDM(-1.)))*(-1.)**(INR2-1)
      IFL4=-IFL3
      RETURN
 30   INR2=INT(1.+3.*RNDM(-1.))
      IFL3=IFLB(IK2-36,INR2)
      IF(INR2-2) 500,501,502
 500  I1=2
      I2=3
      GO TO 505
 501  I1=1
      I2=3
      GO TO 505
 502  I1=1
      I2=2
 505  IF1=IFLB(IK2-36,I1)
      IF2=IFLB(IK2-36,I2)
      IFL4=DQQ(IF1,IF2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IBLE(IK)
      IMPLICIT  INTEGER (I-N)
C   CALCULATION OF PARTICLE BARION NUMBER
      IF(IK-36) 1,1,2
 1    IBLE=0
      RETURN
 2    IBLE=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORPLE(BV,NIN,NFIN,L)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  LORENTZ TRANSFORMATION OF MOMENTUM COMPONENTS PR(I,J)
C  AND ENERGY PR(4,J) (L=1)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      DIMENSION BV(3)
      BVV=BV(1)*BV(1)+BV(2)*BV(2)+BV(3)*BV(3)
      GA=1.D0/SQRT(ABS(1.D0-BVV))
      DO 1 J=NIN,NFIN
      BEP=BV(1)*PR(1,J)+BV(2)*PR(2,J)+BV(3)*PR(3,J)
      GABEP=GA*(GA*BEP/(1.+GA)-L*PR(4,J))
      PR(1,J)=PR(1,J)+GABEP*BV(1)
      PR(2,J)=PR(2,J)+GABEP*BV(2)
      PR(3,J)=PR(3,J)+GABEP*BV(3)
 1    PR(4,J)=GA*(PR(4,J)-L*BEP)
        RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORCLE(BV,NIN,NFIN,L)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  LORENTZ TRANSFORMATION OF COORDINATES (5,6,7)
C  AND TIME (8)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      DIMENSION BV(3)
      BVV=BV(1)*BV(1)+BV(2)*BV(2)+BV(3)*BV(3)
      GA=1.D0/SQRT(ABS(1.D0-BVV))
      DO 1 J=NIN,NFIN
      BEV=BV(1)*PR(5,J)+BV(2)*PR(6,J)+BV(3)*PR(7,J)
      GABEP=GA*(BEV*GA/(GA+1.)-DBLE(L)*PR(8,J))
      PR(5,J)=PR(5,J)+GABEP*BV(1)
      PR(6,J)=PR(6,J)+GABEP*BV(2)
      PR(7,J)=PR(7,J)+GABEP*BV(3)
      PR(8,J)=GA*(PR(8,J)-DBLE(L)*BEV)
 1    CONTINUE
        RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION ISLE(K)
      IMPLICIT  INTEGER (I-N)
C  COMPUTE PARTICLE STRANGE
      COMMON /COMSTR/ ISTR(36)
      IF(K.GT.36) GO TO 1
      ISLE=ISTR(K)
      RETURN
 1    ISLE=ISTR(K-18)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION ZDIST(IFL1,IFL2,ZMIN,ZMAX,IBP)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  COMPUTE FRACTION Z FROM E+PZ OF QUARK OR ANTIQUARK
C  SIMULATION OF Z DISTRIBUTION FROM U(Z)=1.-CX2+3.*CX2*(1.-Z)**2
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      IF(IFL1.LE.3.AND.IBP.EQ.1) GO TO 200
      IF(IFL1.LE.3.AND.IBP.EQ.0) GO TO 100
 201  ZDIST=RNDM(-1.)*(ZMAX-ZMIN)+ZMIN
      YF=0.8+0.2*EXP(-25.*(1.-ZDIST))/(1.-EXP(-25.D0))
      YP=1.
      IF(RNDM(-1.)*YP.LE.YF) RETURN
      GO TO 201
  100 IF(IABS(IFL2).EQ.3)        GO TO 300
  150 ZDIST=RNDM(-1.)*(ZMAX-ZMIN)+ZMIN
      YF=2.*CX2*(1.-ZDIST)+1.-CX2
      YP=2.
      IF(RNDM(-1.)*YP.LE.YF)RETURN
      GO TO 150
 200  ZDIST=RNDM(-1.)*(ZMAX-ZMIN)+ZMIN
      YF=3.*(1.-ZDIST)**2
      YP=3.
      IF(RNDM(-1.)*YP.LE.YF) RETURN
      GO TO 200
  300 ZDIST=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZDIST)**0.5
      IF(RNDM(-1.).GT.YF)    GO  TO  300
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CLUSLE(IFL1,IFL2,KSD,AMCTR)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  HADRONS PRODUCTION BY MEANS CLUSTER BREAKING
C  WITH QUARK AND ANTIQUARK OR QUARK AND DIQUARK IFL1 AND IFL2
C  ON ENDS
C  AMCTR IS MASS OF CLUSTER
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /PRODMA/ PPMAS(50)
      DIMENSION IFL(2),U(3)
      DOUBLE PRECISION PCM,A,B,C
      COMMON /COLRET/LRET
      LOGICAL LRET
      LRET = .FALSE.
      NREP=0
      NFIX=NP
 100  I=NFIX
      NREPA=0          ! REP 21OCT2004
      NREP = NREP + 1
      IF(NREP.LE.NTRIES) GOTO 102
      LRET = .TRUE.
      IF(IPRINT) WRITE(ITLIS,1200) NREP,IFL1,IFL2,AMCTR
1200  FORMAT(3X,' IN CLUSLE NREP GT ',3I8,' AMCTR=',F12.4)
      RETURN
 102  CONTINUE
      KSPIN=0
      IFL(1)=IFL1
      IFL(2)=IFL2
      I=I+2
C  CHOOSE SIDE OF BREAK
      JSIDE=1
C  IF ANY IFL IS A DIQUARK
      IF(KSD.NE.0) GO TO 150
C  IFL(1) AND IFL(2) ARE QUARKS
C  Q,QBAR PAIR
      IFLN=ISIGN(INT(RNDM(-1.)/PUD)+1,-IFL(JSIDE))
      KSD1=KSD
      KSD2=KSD
      GO TO 200
C  IFL(1) OR IFL(2) IS DIQUARK
C Q,QBAR PAIR
150   IPSIGN=IFL(JSIDE)
      IF(IFL(JSIDE).GT.3) GO TO 130
      IPSIGN=-IFL(JSIDE)
      KSD1=0
      KSD2=KSD
      GO TO 135
130   KSD1=KSD
      KSD2=0
135   IFLN=ISIGN(INT(RNDM(-1.)/PUD)+1,IPSIGN)
C+++++++++ SIVOKL 11.11.92
c     NREPA=0     ! REP 21OCT2004
9100  CONTINUE
      NREPA=NREPA+1
      IF(NREPA.GT.10)   GO TO 100
C+++++++++++++++++++++++++
C  IDENTS AND MASSES OF PARTICLES
 200  IPR(I-1)=KI2(IFL(JSIDE),IFLN,KSD1,0)
      IPR(I)=KI2(IFL(3-JSIDE),-IFLN,KSD2,0)
      AM1=AMAS(IPR(I-1))
      AM2=AMAS(IPR(I))
C  IF TOO LOW MASS,START ALL OVER
C+++++++++ SIVOKL 11.11.92
c     IF(AMCTR.LE.AM1+AM2) GO TO 100
      IF(AMCTR.LE.AM1+AM2) GO TO 9100
C+++++++++++++++++++++++++
      A=AMCTR
      B=AM1
      C=AM2
      PCM=SQRT((A*A-B*B-C*C)**2-(2.D0*B*C)**2)/(2.D0*A)
      PA=PCM
C     PROB=2.*PA/AMCTR
C   PROB IS TWO-BODY PHASE SPACE FACTOR
C     IF(RNDM(-1.).GT.PROB) GO TO 100
      U(3)=2.*RNDM(-1.)-1.
      PHI=twpi*RNDM(-1.)
      ST=SQRT(1.-U(3)**2)
      U(1)=ST*COS(PHI)
      U(2)=ST*SIN(PHI)
      PR(1,I-1)=PA*U(1)
      PR(1,I)=-PA*U(1)
      PR(2,I-1)=PA*U(2)
      PR(2,I)=-PA*U(2)
      PR(3,I-1)=PA*U(3)
      PR(3,I)=-PA*U(3)
      PA2=PA**2
      PR(4,I-1)=SQRT(PA2+AM1**2)
      PR(4,I)=SQRT(PA2+AM2**2)
      PPMAS(I-1)=AM1
      PPMAS(I)=AM2
      NP=I
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION ALAMB(X,Y,Z)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C    COMPUTE KINEMATIC FUNCTION
C
      ALAMB=(X-Y-Z)*(X-Y-Z) - 4.*Y*Z
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GETPT(PT0,SIGMA)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   GENERATE DISTRIBUTION WITH 1/(1+B*PT**2)**4
      DATA CON1/1.697652726/,CON2/-.3333333333/
      PT0=CON1*SIGMA*SQRT(RNDM(-1.)**CON2-1.)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PTDGET(PX,PY,SIGMA)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   GENERATE DISTRIBUTION WITH 1/(1+B*PT**2)**4
      DATA CON1/1.697652726/,CON2/-.3333333333/
      PT0=CON1*SIGMA*SQRT(RNDM(-1.)**CON2-1.)
      PHI=twpi*RNDM(-1.)
      PX=PT0*COS(PHI)
      PY=PT0*SIN(PHI)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE STRILE(IFL1,IFL2,KSD0,ECS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  HADRONS PRODUCTION BY MEANS STRING BREAK
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /COMLD/PLDER(50)
      COMMON /CINSID/ INSIDE
      COMMON /KAPPA/ XAP
      COMMON /DATA4/QMAS(9)
      COMMON /DATA5/IFLM(18,2),IFLB(18,3),DQQ(3,3)
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON /COLRET/ LRET
      COMMON /COMQSE/ QSEE,QVSEE
      COMMON /PRODMA/ PPMAS(50)
      LOGICAL QSEE,QVSEE
      DIMENSION W(2),PX1(2),PY1(2),IFL(2),PMTS(2)
      DIMENSION PX1L(2),PY1L(2),V(3),NIN(2),NFIN(2)
      DIMENSION PRR(4,50),PRL(4,50),IPRR(50),IPRL(50)
      DIMENSION P7(50),P8(50),PPMAL(50),PPMAR(50)
      LOGICAL LRET
      LRET = .FALSE.
      NREP = 0
      NFIX=NP
 100  I=NFIX
      NPR=0
      NPL=0
      NP=NFIX
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 102
      LRET=.TRUE.
      RETURN
102   CONTINUE
      IFL(1)=IFL1
      IFL(2)=IFL2
      KSD=KSD0
      IBP=0
C     PROPER MOMENTUM OF STRING
      SIGML=0.000001
C
      IF(KSD.EQ.0) GO TO 13
      CALL PTDGET(PX1L(1),PY1L(1),SIGML)
      PX1L(2)=-PX1L(1)
      PY1L(2)=-PY1L(1)
      IFFL=IFL2
      IF(IFL2.LE.3)IFFL=IFL1
      CALL FLAVD(IFFL,IFLJ,IFLL)
 13   CONTINUE
      DO 1 JT=1,2
      W(JT)=ECS
 1    CONTINUE
      CALL PTDGET(PX1(1),PY1(1),SIGML)
      IF(KSD.EQ.0) GO TO 14
      PX1(1)=PX1L(1)
      PY1(1)=PY1L(1)
14    PX1(2)=-PX1(1)
      PY1(2)=-PY1(1)
C  IF ENERGY IS LOW,ONLY ONE BREAK
      WMAX=WMX(KSD,IFL1,IFL2,PARM,PARB)
      IF(W(1)*W(2).LE.WMAX**2) GO TO 12
C  CHOOSE SIDE. GENERATE A QUARK-ANTIQUARK PAIR FORM HADRON
 2    I=I+1
      JT=1.+2.*RNDM(-1.)
      IF(JT.EQ.1) NPR=NPR+1
      IF(JT.EQ.2) NPL=NPL+1
      KSDN=KSD
      IF(IFL(JT).LE.3) GO TO 5
      IF(RNDM(-1.).GT.0.1) GO TO 5
      IBP=1
      CALL PTDGET(PXLL,PYLL,SIGMA)
      PX1(JT)=-PXLL
      PY1(JT)=-PYLL
      KSDN=0
      IFL(JT)=IFLJ
 5    IFS=-IFL(JT)
      IF(IFL(JT).GT.3) IFS=IFL(JT)
      IFLN=ISIGN(1+INT(RNDM(-1.)/PUD),IFS)
      IFQ=IFL(JT)
      IFQN=IFLN
      IF(IFQ.LT.4) KSDN=0
      IPR(I)=KI2(IFQ,IFQN,KSDN,0)
      IF(KSDN.NE.0) KSD=0
      CALL PTDGET(PX2,PY2,SIGMA)
      PR(1,I)=PX1(JT)+PX2
      PR(2,I)=PY1(JT)+PY2
      PMAS=AMAS(IPR(I))
      PPMAS(I)=PMAS
C   GENERATE Z
      PMTS(3-JT)=0.6
      PMTS(JT)=PMAS**2+PR(1,I)**2+PR(2,I)**2
      IF(PMTS(1)+PMTS(2).GE.0.9*W(1)*W(2)) GO TO 100
      ZMIN=PMTS(JT)/(W(1)*W(2))
      ZMAX=1.-PMTS(3-JT)/(W(1)*W(2))
      IF(ZMIN.GE.ZMAX) GO TO 100
      Z=ZDIST(IFL(JT),IFQN,ZMIN,ZMAX,IBP)
      PR(3,I)=0.5*(Z*W(JT)-PMTS(JT)/(Z*W(JT)))*(-1.)**(JT+1)
      PR(4,I)=0.5*(Z*W(JT)+PMTS(JT)/(Z*W(JT)))
      IF(.NOT.(JT.EQ.1)) GO TO 282
      IPRR(NPR)=IPR(I)
      PRR(1,NPR)=PR(1,I)
      PRR(2,NPR)=PR(2,I)
      PRR(3,NPR)=PR(3,I)
      PRR(4,NPR)=PR(4,I)
      PPMAR(NPR)=PPMAS(I)
282   IF(.NOT.(JT.EQ.2)) GO TO 283
      IPRL(NPL)=IPR(I)
      PRL(1,NPL)=PR(1,I)
      PRL(2,NPL)=PR(2,I)
      PRL(3,NPL)=PR(3,I)
      PRL(4,NPL)=PR(4,I)
      PPMAL(NPL)=PPMAS(I)
283   IF(IBP.EQ.0) GO TO 10
      IFL(JT)=DQQ(IABS(IFLN),IFLL)
      IFLJ=IABS(IFLN)
      PX1L(JT)=PX1L(JT)+PXLL-PX2
      PY1L(JT)=PY1L(JT)+PYLL-PY2
      PX1(JT)=PX1L(JT)
      PY1(JT)=PY1L(JT)
      IF(IFL(JT).EQ.4.OR.IFL(JT).EQ.5.OR.IFL(JT).EQ.9) KSD=2
      IBP=0
      GO TO 11
 10   CONTINUE
      IFL(JT)=-IFLN
      PX1(JT)=-PX2
      PY1(JT)=-PY2
 11   CONTINUE
      W(1)=W(1)-PR(4,I)-PR(3,I)
      W(2)=W(2)-PR(4,I)+PR(3,I)
C   IF ENOUCH ENERGY LEFT,CONTINUE GENERATE
 12   IKB=KI2(IFL(1),IFL(2),KSD,2)
      PARC=0.2
      IF(IABS(IFL(1)).EQ.3.OR.IABS(IFL(2)).EQ.3) PARC=0.5
      AMB=AMAS(IKB)+PARC
      P1X=PX1(1)+PX1(2)
      P1Y=PY1(1)+PY1(2)
      PT12=P1X**2+P1Y**2
      W12=W(1)*W(2)
      AMS2=W12-PT12
      IF(AMS2.LT.AMB**2) GO TO 100
      WMAX=WMX(KSD,IFL(1),IFL(2),PARM,PARB)
      IF(W12.GE.WMAX**2) GO TO 2
C  GIVEN FINAL TWO HADRON
 3    NP=I
      AMC=SQRT(AMS2)
      EC=(W(1)+W(2))/2.0
      V(1)=P1X/EC
      V(2)=P1Y/EC
      V(3)=(W(1)-W(2))/(2.0*EC)
      NIN(1)=NP+1
      CALL CLUSLE(IFL(1),IFL(2),KSD,AMC)
      IF(LRET) GO TO 100
      NFIN(1)=NP
      CALL LORPLE(V,NIN(1),NFIN(1),-1)
      NI=NIN(1)
      NF=NFIN(1)
      NPR=NPR+1
      NPL=NPL+1
      IPRR(NPR)=IPR(NI)
      PRR(1,NPR)=PR(1,NI)
      PRR(2,NPR)=PR(2,NI)
      PRR(3,NPR)=PR(3,NI)
      PRR(4,NPR)=PR(4,NI)
      PPMAR(NPR)=PPMAS(NI)
      IPRL(NPL)=IPR(NF)
      PRL(1,NPL)=PR(1,NF)
      PRL(2,NPL)=PR(2,NF)
      PRL(3,NPL)=PR(3,NF)
      PRL(4,NPL)=PR(4,NF)
      PPMAL(NPL)=PPMAS(NF)
      JJ=NFIX
      DO 284 J=1,NPR
      JJ=JJ+1
      IPR(JJ)=IPRR(J)
      PR(1,JJ)=PRR(1,J)
      PR(2,JJ)=PRR(2,J)
      PR(3,JJ)=PRR(3,J)
      PR(4,JJ)=PRR(4,J)
      PPMAS(JJ)=PPMAR(J)
284   CONTINUE
      JJ=NFIX+NPR
      DO 285 J=1,NPL
      JJ=JJ+1
      K=NPL-J+1
      IPR(JJ)=IPRL(K)
      PR(1,JJ)=PRL(1,K)
      PR(2,JJ)=PRL(2,K)
      PR(3,JJ)=PRL(3,K)
      PR(4,JJ)=PRL(4,K)
      PPMAS(JJ)=PPMAL(K)
285   CONTINUE
      N1=NFIX+1
      N2=NFIX+NPR+NPL-1
      IF(INSIDE.NE.0) GO TO 1252
C------------------------------------------------------C
C-----  CONSTITUENT     TIME           ----------------C
C------------------------------------------------------C
      DO 286 J=N1,N2
      P3S=0.
      ES=0.
      DO 287 L=N1,J
      P3S=P3S+PR(3,L)
 287  ES=ES+PR(4,L)
      TI=(ECS-2.*P3S)/(2.*XAP)
      ZI=(ECS-2.*ES)/(2.*XAP)
      IF(J.NE.N2) GO TO 288
      TII=TI
      ZII=ZI
 288  PR(5,J)=0.
      PR(6,J)=0.
      PR(7,J)=ZI
      PR(8,J)=TI
 286  CONTINUE
      PR(5,N2+1)=0.
      PR(6,N2+1)=0.
      PR(7,N2+1)=ZII
      PR(8,N2+1)=TII
C]]]]]]]]]]]]
      IF(N2.LE.1) GO TO 1253
      LN11=N1+1
      DO 1389 L=LN11,N2
      P7(L)=0.5*(PR(7,L-1)+PR(7,L))
      P8(L)=0.5*(PR(8,L-1)+PR(8,L))
1389  CONTINUE
      DO 1489 L=LN11,N2
      PR(7,L)=P7(L)
      PR(8,L)=P8(L)
1489  CONTINUE
C]]]]]]]]]]]]
      GO TO 1253
C------------------------------------------------------C
C-----  INSIDE-OUTSIDE  TIME           ----------------C
C------------------------------------------------------C
1252  CONTINUE
      DO 1286 J=N1,NP
      P3S=0.
      ES=0.
      NJ=J-1
      IF(NJ.EQ.0) GO TO 1289
      DO 1287 L=N1,NJ
      P3S=P3S+PR(3,L)
1287  ES=ES+PR(4,L)
1289  TI=(ECS-2.*P3S+PR(4,J)-PR(3,J))/(2.*XAP)
      ZI=(ECS-2.*ES-PR(4,J)+PR(3,J))/(2.*XAP)
      PR(5,J)=0.
      PR(6,J)=0.
      PR(7,J)=ZI
      PR(8,J)=TI
1286  CONTINUE
1253  CONTINUE
C-------------------------------------------------------------
      DO 386 J=N1,NP
386   PLDER(J)=0.
      IB1=IBLE(IPR(N1))
      IB2=IBLE(IPR(NP))
      PLDER(N1)=.667
      PLDER(NP)=.667
      IF(IB1.EQ.0) PLDER(N1)=.5
      IF(IB2.EQ.0) PLDER(NP)=.5
      IF(.NOT.QVSEE) RETURN
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 387
      IF(IB1.EQ.0) PLDER(N1)=0.
      IF(IB2.EQ.0) PLDER(NP)=0.
      RETURN
387   RM=RNDM(-1.)
      IF(RM.GT.0.5) PLDER(N1)=0.
      IF(RM.LE.0.5) PLDER(NP)=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PLANLE(IK1,IB1,IK2,IB2,S,IEL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C  CALCULATION OF PLANAR GRAPH
C  ONLY FOR MESON BARYON COLLISION
C
      COMMON /COMLD/ PLDER(50)
      COMMON /DATA5/IFLM(18,2),IFLB(18,3),DQQ(3,3)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /PRODT/ IORD(50)
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON /COLRET/ LRET
      LOGICAL LRET
      DIMENSION NIN(2),NFIN(2)
      LRET = .FALSE.
      IF(IB1.EQ.1) GO TO 1003
      IF1P=IABS(IFLM(IK1,2))
      IF1PN=IFLM(IK1,1)
      IF1=IFLB(IK2-36,1)
      IF2=IFLB(IK2-36,2)
      IF3=IFLB(IK2-36,3)
      IF(IF1P.EQ.IF1) IF2P=DQQ(IF2,IF3)
      IF(IF1P.EQ.IF2) IF2P=DQQ(IF1,IF3)
      IF(IF1P.EQ.IF3) IF2P=DQQ(IF1,IF2)
      IF(IF1P.NE.IF1.AND.IF1P.NE.IF2.AND.IF1P.NE.IF3) GO TO 100
      CALL KSPIN(IF1P,KSP,IK2,IB2)
      GO TO 1004
 1003 IF1P=IABS(IFLM(IK2,2))
      IF1PN=IFLM(IK2,1)
      IF1=IFLB(IK1-36,1)
      IF2=IFLB(IK1-36,2)
      IF3=IFLB(IK1-36,3)
      IF(IF1P.EQ.IF1) IF2P=DQQ(IF2,IF3)
      IF(IF1P.EQ.IF2) IF2P=DQQ(IF1,IF3)
      IF(IF1P.EQ.IF3) IF2P=DQQ(IF1,IF2)
      IF(IF1P.NE.IF1.AND.IF1P.NE.IF2.AND.IF1P.NE.IF3) GO TO 100
      CALL KSPIN(IF1P,KSP,IK1,IB1)
 1004 IKN0=KI2(IF2P,IF1PN,KSP,2)
      IKN=KI2(IF2P,IF1PN,KSP,1)
      PAM0=AMAS(IKN0)
      PAM=AMAS(IKN)
      NIN(1)=NP+1
      AMS=SQRT(S)
      PARC=0.
      IF(AMS.GT.PAM0+PARC) GO TO 1007
      IEL=1
      GO TO 100
1007  IF(AMS.GE.PAM+SWMAX) GO TO 1008
      CALL CLUSLE(IF1PN,IF2P,KSP,AMS)
      IF(LRET) RETURN
      CALL TIFILE(NIN(1),NP,AMS)
      GO TO 1014
1008  CALL STRILE(IF1PN,IF2P,KSP,AMS)
      NIN1=NIN(1)
      IF(LRET) RETURN
 1014 NFIN(2)=NP
      NIN1=NIN(1)
      NFIN2=NFIN(2)
      DO 1006 JO=NIN1,NFIN2
 1006 IORD(JO)=0
      RETURN
 100  LRET=.TRUE.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CYLLE(IK1,IB1,AM1,IK2,IB2,AM2,P1,IBINA)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  CALCULATION OF CYLINDRICAL GRAPH
C
      COMMON /COMLD/ PLDER(50)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /PRODT/ IORD(50)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON /COLRET/ LRET
      DIMENSION P1(3),NIN(2),NFIN(2),VS1(3),VS2(3)
      DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      LOGICAL LRET
C ______   PRIMORDIAL MOMENTUM OF PARTONS  _____________________
*      DATA SIGMAI/0.4/
      DATA SIGMAI/0.5/    ! 15.06.98
C ______________________________________________________________
      LRET = .FALSE.
      NREP = 0
      NPOLD=NP
      AZ12=AM1**2
      AZ22=AM2**2
      PZER2=P1(3)**2
      P0I=SQRT(SPQ(P1,P1))
      EI1=SQRT(AM1**2+P0I**2)
      EI2=SQRT(AM2**2+P0I**2)
      IF(EI1+EI2-AM1-AM2.GT.0.3) GO TO 100
      LRET = .TRUE.
      IBINA=1
      RETURN
100   CONTINUE
      NP = NPOLD
      NREP=NREP+1
      IF(NREP.LT.  NTRIES) GO TO 200
C     WRITE(ITLIS,1001)IK1,IK2,AM1,AM2,PZER2
1001  FORMAT(1X,'CYLLE:NREP > NTRIES, IK1,IK2,AM1,AM2,PZER2=',2I5,
     *3E10.4)
      IBINA=1
      LRET = .TRUE.
      RETURN
200   CONTINUE
      CALL FLAVO(IB1,IK1,IB2,IK2,IFL1,IFL2,IFL3,IFL4)
      CALL KSPIN(IFL3,KS2,IK2,IB2)
      CALL KSPIN(IFL1,KS1,IK1,IB1)
C  HADRON GENERATION BY MEANS FORMING AND BREAKING STRINGS
      IF(IFL1) 1,1,2
 1    IF11=IFL2
      IF22=IFL1
      GO TO 3
 2    IF11=IFL1
      IF22=IFL2
 3    CONTINUE
      IF(IB2.EQ.1) GO TO 102
      IF(IFL3) 101,101,102
 101  IF44=IFL3
      IF33=IFL4
      GO TO 103
 102  IF33=IFL3
      IF44=IFL4
 103  CONTINUE
      IKN1=KI2(IF44,IF11,KS2,1)
      IKN2=KI2(IF22,IF33,KS1,1)
      PAM1=AMAS(IKN1)
      PAM2=AMAS(IKN2)
C  MOMENTUM OF QUARKS
      X1MIN=0.
      X3MIN=0.
      IS1=ISLE(IK1)
      IF11S=IF11
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IF11).EQ.3) IF11S=IF22
      X1=XDIST(X1MIN,IB1,IS1,IF11S)
      IS2=ISLE(IK2)
      IF33S=IF33
      IF(IS2.NE.0.AND.IB2.EQ.0.AND.IABS(IF33).EQ.3) IF33S=IF44
      X3=XDIST(X3MIN,IB2,IS2,IF33S)
      X2=1.-X1
      X4=1.-X3
      PZ11=P1(3)*X1
      PZ22=P1(3)*X2
      PZ33=-P1(3)*X3
      PZ44=-P1(3)*X4
C    COMPUTE PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
  160 CALL GETPT(PT1,SIGMAI)
      AMQ21=AZ12*(AZ12+4.*X1*X2*PZER2)/(4.*(PZER2+AZ12))-PT1**2
      PX11=PT1*COS(PHI)
      PY11=PT1*SIN(PHI)
      PX22=-PX11
      PY22=-PY11
      PHI=twpi*RNDM(-1.)
  170 CALL GETPT(PT3,SIGMAI)
      AMQ22=AZ22*(AZ22+4.*X3*X4*PZER2)/(4.*(PZER2+AZ22))-PT3**2
      PX33=PT3*COS(PHI)
      PY33=PT3*SIN(PHI)
      PX44=-PX33
      PY44=-PY33
C  COMPUTE OF STABLE PARTICLE MASSES
      IKN01=KI2(IF44,IF11,KS2,2)
      IKN02=KI2(IF22,IF33,KS1,2)
      PAM01=AMAS(IKN01)
      PAM02=AMAS(IKN02)
C   WILL BE START ALL OVER
      PARC1=0.2
      PARC2=0.2
C  STRINGS OR CLUSTER DECAY
      E11=PX11**2+PY11**2+PZ11**2+AMQ21
      IF(E11.LT.0.) GO TO 100
      E22=PX22**2+PY22**2+PZ22**2+AMQ21
      IF(E22.LT.0.) GO TO 100
      E33=PX33**2+PY33**2+PZ33**2+AMQ22
      IF(E33.LT.0.) GO TO 100
      E44=PX44**2+PY44**2+PZ44**2+AMQ22
      IF(E44.LT.0.) GO TO 100
      E11=SQRT(E11)
      E22=SQRT(E22)
      E33=SQRT(E33)
      E44=SQRT(E44)
      E1=E11+E44
      E2=E22+E33
      AMS2=E2**2-(PX22+PX33)**2-(PY22+PY33)**2-(PZ22+PZ33)**2
      IF(AMS2.LT.(PAM02+PARC2)**2) GO TO 100
      AMS1=E1**2-(PX11+PX44)**2-(PY11+PY44)**2-(PZ11+PZ44)**2
      IF(AMS1.LT.(PAM01+PARC1)**2) GO TO 100
      AMS1=SQRT(AMS1)
      AMS2=SQRT(AMS2)
C  VELOCITIES OF CM STRINGS
      VS1(1)=(PX11+PX44)/E1
      VS1(2)=(PY11+PY44)/E1
      VS1(3)=(PZ11+PZ44)/E1
      VS2(1)=(PX22+PX33)/E2
      VS2(2)=(PY22+PY33)/E2
      VS2(3)=(PZ22+PZ33)/E2
C  BREAK OF STRING WITH QUARKS OF NUMBER 1 AND 4
      NIN(1)=NP+1
      IF(AMS1.LE.PAM1+SWMAX) GO TO 714
      CALL STRILE(IF11,IF44,KS2,AMS1)
      IF(LRET) GO TO 100
      NFIN(1)=NP
      GO TO 715
714   CALL CLUSLE(IF11,IF44,KS2,AMS1)
      IF(LRET) GO TO 100
      NFIN(1)=NP
      CALL TIFILE(NIN(1),NFIN(1),AMS1)
715   NIN1=NIN(1)
      NFIN1=NFIN(1)
      L=1
      PPX1(1)=PX11
      PPX1(2)=PY11
      PPX1(3)=PZ11
      CALL LORLLE(VS1,PPX1,E11,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PR(1,J)
      PPX1(2)=PR(2,J)
      PPX1(3)=PR(3,J)
      CALL ROTAM(CT,ST,CFI,SFI,PPX1,PPX2,L)
      PR(1,J)=PPX2(1)
      PR(2,J)=PPX2(2)
      PR(3,J)=PPX2(3)
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
510   CONTINUE
      NIN(2)=NP+1
      IF(AMS2.LE.PAM2+SWMAX)  GO TO 914
C  BREAK OF STRING WITH QUARKS OF NUMBER 3 AND 2 ON ENDS
      CALL STRILE(IF22,IF33,KS1,AMS2)
      IF(LRET) GO TO 100
      NFIN(2)=NP
      GO TO 915
914   CALL CLUSLE(IF22,IF33,KS1,AMS2)
      IF(LRET) GO TO 100
      NFIN(2)=NP
      CALL TIFILE(NIN(2),NFIN(2),AMS2)
915   NIN2=NIN(2)
      NFIN2=NFIN(2)
      L=1
      PPX1(1)=PX22
      PPX1(2)=PY22
      PPX1(3)=PZ22
      CALL LORLLE(VS2,PPX1,E22,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN2,NFIN2
      PPX1(1)=PR(1,J)
      PPX1(2)=PR(2,J)
      PPX1(3)=PR(3,J)
      CALL ROTAM(CT,ST,CFI,SFI,PPX1,PPX2,L)
      PR(1,J)=PPX2(1)
      PR(2,J)=PPX2(2)
      PR(3,J)=PPX2(3)
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
610   CONTINUE
C  RETURN IN OVERALL CM FRAME
      CALL LORPLE(VS1,NIN(1),NFIN(1),-1)
      CALL LORCLE(VS1,NIN(1),NFIN(1),-1)
      NIN1=NIN(1)
      NFIN1=NFIN(1)
      NIN2=NIN(2)
      NFIN2=NFIN(2)
      DO 813 JO=NIN1,NFIN1
 813  IORD(JO)=0
      DO 814 JO=NIN2,NFIN2
 814  IORD(JO)=0
      CALL LORPLE(VS2,NIN(2),NFIN(2),-1)
      CALL LORCLE(VS2,NIN(2),NFIN(2),-1)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE UNCYLE(IK1,IB1,AM1,IK2,IB2,AM2,P1,IBINA)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C      COMPUTE UNDEVELOPED CYLINDER DIAGRAM
C
      COMMON /COMLD/ PLDER(50)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /PRODT/ IORD(50)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON /COLRET/ LRET
      COMMON /PRODMA/ PPMAS(50)
      DIMENSION V(3),P1(3)
      DIMENSION NIN(2),NFIN(2)
      DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      LOGICAL LRET
C    INITIALIZE
      LRET = .FALSE.
C ______   PRIMORDIAL MOMENTUM OF PARTONS  _____________________
*      DATA SIGMAI/0.4/
      DATA SIGMAI/0.5/     ! 15.06.98
C ______________________________________________________________
C    INITIALIZE
      LRET = .FALSE.
      NREP=0
      P0=SQRT(SPQ(P1,P1))
      AZ12=AM1**2
      AZ22=AM2**2
      PZER2=P0**2
      E1=SQRT(AM1**2+P0**2)
      E2=SQRT(AM2**2+P0**2)
      ECM=E1+E2
      SCM=ECM**2
      IF(ECM-AM1-AM2.GT.0.3) GO TO 150
      LRET = .TRUE.
      IBINA=1
      RETURN
150   PSIGN=-1.
      NREP = NREP + 1
      IF(NREP.LE.  NTRIES) GO TO 200
C     WRITE(ITLIS,1001)IK1,IK2,AMA,AMB,ZER2B
1001   FORMAT(1X,'UCYLLE:NREP > NTRIES,IK1,IK2,AM1,AM2,PZER2',
     *2I4,3(1X,F7.3))
      IBINA=1
      LRET = .TRUE.
      RETURN
200   CONTINUE
      CALL FLAVO(IB1,IK1,IB2,IK2,IFL11,IFL22,IFL33,IFL44)
      CALL KSPIN(IFL11,KS11,IK1,IB1)
      CALL KSPIN(IFL33,KS22,IK2,IB2)
C    COMPUTE X VALUES FOR PARTONS
      XMIN=0.
      RND=RNDM(-1.)
      IS1=ISLE(IK1)
      IFL11S=IFL11
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL11).EQ.3) IFL11S=IFL22
      X11=XDIST(XMIN,IB1,IS1,IFL11S)
      X22=1.-X11
      IS2=ISLE(IK2)
      IFL33S=IFL33
      IF(IS2.NE.0.AND.IB2.EQ.0.AND.IABS(IFL33).EQ.3) IFL33S=IFL44
      X33=XDIST(XMIN,IB2,IS2,IFL33S)
      X44=1.-X33
C    COMPUTE PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
  160 CALL GETPT(PT11,SIGMAI)
      AMQ21=AZ12*(AZ12+4.*X11*X22*PZER2)/(4.*(PZER2+AZ12))-PT11**2
      PX11=PT11*COS(PHI)
      PY11=PT11*SIN(PHI)
      PX22=-PX11
      PY22=-PY11
      PHI=twpi*RNDM(-1.)
  170 CALL GETPT(PT33,SIGMAI)
      AMQ22=AZ22*(AZ22+4.*X33*X44*PZER2)/(4.*(PZER2+AZ22))-PT33**2
      PX33=PT33*COS(PHI)
      PY33=PT33*SIN(PHI)
      PX44=-PX33
      PY44=-PY33
      IF(IFL11.GT.0) GO TO 130
      IFL1=IFL22
      PX1=PX22
      PY1=PY22
      IFL2=IFL11
      PX2=PX11
      PY2=PY11
      KS1=KS11
      X1=X22
      X2=X11
      GO TO 140
130   IFL1=IFL11
      PX1=PX11
      PY1=PY11
      IFL2=IFL22
      PX2=PX22
      PY2=PY22
      KS1=KS11
      X1=X11
      X2=X22
140   IF(IB2.EQ.1) GO TO 102
      IF(IFL33.GT.0) GO TO 102
      IFL4=IFL33
      PX4=PX33
      PY4=PY33
      IFL3=IFL44
      PX3=PX44
      PY3=PY44
      KS2=KS22
      X3=X44
      X4=X33
      GO TO 103
102   IFL3=IFL33
      PX3=PX33
      PY3=PY33
      IFL4=IFL44
      PX4=PX44
      KS2=KS22
      PY4=PY44
      X3=X33
      X4=X44
103   CONTINUE
      PXH=PX1+PX4
      PYH=PY1+PY4
      X01=X2
      IF(X01.EQ.0.) X01=X1
      X02=X3
      IF(X02.EQ.0.) X02=X4
      PX01=PX2
      PY01=PY2
      PX02=PX3
      PY02=PY3
      KS02=KS1
      IDH=KI2(IFL1,IFL4,KS2,0)
      AMH=AMAS(IDH)
      IFL01=IFL2
      IFL02=IFL3
      IF(RND.GE.0.5) GO TO 100
      PXH=PX2+PX3
      PYH=PY2+PY3
      X01=X1
      IF(X01.EQ.0.) X01=X2
      X02=X4
      IF(X02.EQ.0.) X02=X3
      IFL01=IFL1
      IFL02=IFL4
      IDH=KI2(IFL2,IFL3,KS1,0)
      AMH=AMAS(IDH)
      PX01=PX1
      PY01=PY1
      PX02=PX4
      PY02=PY4
      KS02=KS2
 100  P01=X01*P0
      P02=X02*P0*PSIGN
      E01=AMQ21+P01**2+PX01**2+PY01**2
      IF(E01.LT.0.) GO TO 150
      E02=AMQ22+P02**2+PX02**2+PY02**2
      IF(E02.LT.0.) GO TO 150
      E01=SQRT(E01)
      E02=SQRT(E02)
      AMDTR=(E01+E02)**2-(P01+P02)**2-(PX01+PX02)**2-
     *(PY01+PY02)**2
      IDH1=KI2(IFL01,IFL02,KS02,2)
      PARC=0.001
      AMHB=AMAS(IDH1)+PARC
      IF(AMDTR.LE.AMHB**2) GO TO 150
      AMD=SQRT(AMDTR)
      IF(ECM.LE.AMD+AMH) GO TO 150
      ALA=ALAMB(SCM,AMDTR,AMH**2)
      P0H=SQRT(ALA)/(2.*ECM)
      PTHX=-(PX01+PX02)
      PTHY=-(PY01+PY02)
      DTRM=P0H**2-PTHX**2-PTHY**2
      IF(DTRM.LT.0.) GO TO 150
      PZH0=SQRT(DTRM)
      PZH=SIGN(PZH0,-(P01+P02))
      ED=SQRT(AMDTR+P0H**2)
      EH=SQRT(AMH**2+P0H**2)
      PSIGN=SIGN(1.D0,-PZH)
      V(1)=-PTHX/ED
      V(2)=-PTHY/ED
      V(3)=PSIGN*PZH0/ED
      IDHR=KI2(IFL01,IFL02,KS02,1)
      AMHS=AMAS(IDHR)+SWMAX
      IF(AMD.GT.AMHS) GO TO 300
      NFIX=NP
      NIN(1)=NP+1
      CALL CLUSLE(IFL01,IFL02,KS02,AMD)
      IF(LRET) GO TO 150
      NFIN(1)=NP
      CALL TIFILE(NIN(1),NFIN(1),AMD)
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      NIN1=NIN(1)
      NFIN1=NFIN(1)
      L=1
      CALL LORLLE(V,PPX1,E01,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
610   CONTINUE
      CALL LORPLE(V,NIN(1),NFIN(1),-1)
      CALL LORCLE(V,NIN(1),NFIN(1),-1)
      NPRODS=NP-NFIX
      GO TO 400
300   NFIX=NP
      NIN(1)=NP+1
      CALL STRILE(IFL01,IFL02,KS02,AMD)
      IF(LRET) GO TO 150
      NFIN(1)=NP
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      NIN1=NIN(1)
      NFIN1=NFIN(1)
      L=1
      CALL LORLLE(V,PPX1,E01,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PR(1,J)
      PPX1(2)=PR(2,J)
      PPX1(3)=PR(3,J)
      CALL ROTAM(CT,ST,CFI,SFI,PPX1,PPX2,L)
      PR(1,J)=PPX2(1)
      PR(2,J)=PPX2(2)
      PR(3,J)=PPX2(3)
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
510   CONTINUE
      CALL LORPLE(V,NIN(1),NFIN(1),-1)
      CALL LORCLE(V,NIN(1),NFIN(1),-1)
      NPRODS=NP-NFIX
400   NIN1=NIN(1)
      NFIN1=NFIN(1)
      DO 350 J=NIN1,NFIN1
350   IORD(J)=0
      NP=NP+1
      IORD(NP)=0
      PR(1,NP)=PTHX
      PR(2,NP)=PTHY
      PR(3,NP)=PZH
      PR(4,NP)=EH
      PPMAS(NP)=AMH
      IPR(NP)=IDH
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DIFSCA(IFL01,IFL02,KS01,IK1,AM1,
     *IFL03,IFL04,KS02,IK2,AM2,P1,IBINA)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C     COMPUTE LOW MASS DIFFRACTION
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /PRODT/ IORD(50)
      COMMON /ORDER/ IRD1,IRD2
      COMMON /PRODMA/ PPMAS(50)
      COMMON /COMLD/ PLDER(50)
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /COMCUT/ PARM,PARB,SWMAX
      DIMENSION V(3),P1(3)
      DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      DIMENSION GAMA(3),AMR(3)
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DATA SIGMA1/0.23/,SIGMA2/0.30/
      DATA GAMA/0.032,0.032,0.032/
      DATA AMR /1.47,1.10,1.30/
      LRET = .FALSE.
      NREP=0
      XMIN=0.
      P0=SQRT(SPQ(P1,P1))
      E1=SQRT(P0**2+AM1**2)
      E2=SQRT(P0**2+AM2**2)
      ECM=E1+E2
      SCM=ECM**2
      DS=ECM-AM1-AM2
      IF(DS.GT.0.3) GO TO 150
      IBINA=1
      RETURN
150   CONTINUE
      AMQ1=0.
      AMQ2=0.
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 151
C     WRITE(ITLIS,1001)IK1,IK2,AM1,AM2,P0
1001   FORMAT(1X,'DIFSCA:NREP > NTRIES,IK1,IK2,AM1,AM2,P0=',
     *2I4,3(1X,F7.3))
      LRET = .TRUE.
      RETURN
151   CONTINUE
      IFL1=IFL01
      IFL2=IFL02
      IF(IFL2.GT.3) AMQ2=0
      KS1=KS01
      IKA=IK1
      IKB=IK2
      AMA=AM1
      AMB=AM2
      IRDA=IRD1
      IRDB=IRD2
      PSIGN=-1.
      IF(RNDM(-1.).GT.0.5) GO TO 100
      PSIGN=1.
      IFL1=IFL03
      IFL2=IFL04
      IF(IFL2.GT.3) AMQ2=0.
      KS1=KS02
      IKA=IK2
      IKB=IK1
      AMB=AM1
      AMA=AM2
      IRDA=IRD2
      IRDB=IRD1
C    COMPUTE X VALUES FOR PARTONS
100   IBA=IBLE(IKA)
      ISA=ISLE(IKA)
      IFL1S=IFL1
      IF(ISA.NE.0.AND.IBA.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      X1=XDIST(XMIN,IBA,ISA,IFL1S)
      X2=1.-X1
C     COMPUTE PT VALUE FOR HADRON
140   CONTINUE
      CALL GETPT(PT1,SIGMA2)
      PHI=twpi*RNDM(-1.)
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      CALL GETEXP(PT,SIGMA1)
      PHI=twpi*RNDM(-1.)
      PTX=PT*COS(PHI)
      PTY=PT*SIN(PHI)
      PX2=-PX1+PTX
      PY2=-PY1+PTY
      PT2=PX2**2+PY2**2
      AMH=AMB
      E1=SQRT(AMQ1**2+(P0*X1)**2+PX1**2+PY1**2)
      E2=SQRT(AMQ2**2+(P0*X2)**2+PX2**2+PY2**2)
      AMDTR=(E1+E2)**2-P0**2-PT**2
      IDH=KI2(IFL1,IFL2,KS1,2)
      IDHR=KI2(IFL1,IFL2,KS1,1)
      AMHR=AMAS(IDHR)
      PARC1=0.001
      AMHB=AMAS(IDH)+PARC1
      IF(AMDTR.GE.AMHB**2)       GO TO 200
                                 GO TO 140
 200  AMD=SQRT(AMDTR)
      IF(ECM.LE.AMD+AMH)         GO TO 140
      IBR=IBLE(IDHR)
      ISR=ISLE(IDHR)
      GAMRES=GAMA(3)
      AMRES=AMR(3)
      IF(IBR.NE.0)               AMRES=AMR(1)
      IF(IBR.EQ.0.AND.ISR.NE.0)  AMRES=AMR(2)
      IF(AMD.GT.AMRES)           GO  TO  162
      IF(DS.LT.1.0.AND.IBA.EQ.0) GO  TO  162
      ARGWG=-(AMD-AMRES)**2/GAMRES
      IF(ARGWG.LE.-30.) ARGWG=-30.
      WG=EXP(ARGWG)
      IF(RNDM(-1.).GT.WG)         GO  TO  140
 162  ALA=ALAMB(SCM,AMDTR,AMH**2)
      P0H=SQRT(ALA)/(2.*ECM)
      DTRM=P0H**2-PT**2
      IF(DTRM.LT.0.)             GO TO 140
      PZH=SIGN(SQRT(DTRM),-PSIGN)
      ED=SQRT(AMD**2+P0H**2)
      V(1)=PTX/ED
      V(2)=PTY/ED
      V(3)=PZH/ED
      IF(AMD.GT.AMHR+SWMAX) GO TO 300
      NFIX=NP
      NIN1=NP+1
      CALL CLUSLE(IFL1,IFL2,KS1,AMD)
      IF(LRET) GO TO 150
      NFIN1=NP
      CALL TIFILE(NIN1,NFIN1,AMD)
      L=1
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=-PSIGN*P0*X1
      CALL LORLLE(V,PPX1,E1,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
610   CONTINUE
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
      GO TO 400
300   NFIX=NP
      NIN1=NP+1
      CALL STRILE(IFL1,IFL2,KS1,AMD)
      IF(LRET) GO TO 150
      NFIN1=NP
      L=1
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=-PSIGN*P0*X1
      CALL LORLLE(V,PPX1,E1,L)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PR(1,J)
      PPX1(2)=PR(2,J)
      PPX1(3)=PR(3,J)
      CALL ROTAM(CT,ST,CFI,SFI,PPX1,PPX2,L)
      PR(1,J)=PPX2(1)
      PR(2,J)=PPX2(2)
      PR(3,J)=PPX2(3)
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,L)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
510   CONTINUE
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
400   CONTINUE
      DO 350 J=NIN1,NFIN1
C     PPMAS(J)=AMAS(IPR(J))
C     PPMAS(J)=SQRT(PR(4,J)**2-PR(1,J)**2-PR(2,J)**2-PR(3,J)**2)
350   IORD(J)=IRDA
      NP=NP+1
      IORD(NP)=IRDB
      IPR(NP)=IKB
      PR(1,NP)=-PTX
      PR(2,NP)=-PTY
      PR(3,NP)=SIGN(PZH,-PZH)
      PR(4,NP)=SQRT(AMH**2+PR(1,NP)**2+
     *PR(2,NP)**2+PR(3,NP)**2)
      PPMAS(NP)=AMH
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GETEXP(PT,SIGMAD)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      DRND=RNDM(-1.)
      PT=SIGMAD*SQRT(-LOG(DRND))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ELASLE(P,AMA,AMB,IK01,IK02)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  MONTE CARLO SIMULATION OF ELASTIC HADRONIC COLLISION
C--------------------------------------------------------
C   HADRON-NUCLEON COLLISION ONLY
C--------------------------------------------------------
      COMMON /COMKI1/ HLA2,HLB2,W,INUMA
      COMMON /COMKI2/ELA,ELB,PLALB
      COMMON /CALC/HA,HB,HA2,HB2
      COMMON /BEES/B,BFOR
      DIMENSION P(3),PA(3)
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /PRODT/ IORD(50)
      COMMON /ORDER/ IRD1,IRD2
      COMMON /DATA3/ POPB(10)
      COMMON /PRODMA/PPMAS(50)
      COMMON /COMELX/ SIGEL
      COMMON /COMLD/PLDER(50)
      COMMON /COMCRO/ SIGTOT
      DIMENSION P1(3),P2(3)
      IEXE =0
      HLA=AMA
      HLB=AMB
      HLA2=HLA*HLA
      HLB2=HLB*HLB
      E1=SQRT(HLA2+SPQ(P,P))
      E2=SQRT(HLB2+SPQ(P,P))
      S=(E1+E2)**2
C   S=(PA+PB)**2
C   W= CENTRE OF MASS (C.M.) ENERGY
      W=SQRT(S)
C  TKIN=KINETIC ENERGY OF PROJECTILE IN TARGET REST FRAME
      TKIN=(S-HLA2-HLB2)/(2.0*HLB)-HLA
C   PLALB=CM MOMENTUM OF A OR B IN ELASTIC EVENTS
      PLALB=SQRT(ALAMB(S,HLA2,HLB2))/(2.0*W)
C   ELA=CM ENERGY OF A IN ELASTIC EVENT *** ELB=SAME FOR B
      ELA=(S+HLA2-HLB2)/(2.0*W)
      ELB=(S+HLB2-HLA2)/(2.0*W)
      IK1=IK01
      IK2=IK02
      AMN1=HLA
      AMN2=HLB
      HA=AMN1
      HB=AMN2
      INUMA=1
      IF(IK02.LE.36.OR.IK01.LE.36) INUMA=0
      IF(IK02.LE.36.AND.IK01.LE.36) INUMA=2
      HA2=HA*HA
      HB2=HB*HB
      TOBR=10.0
      IF(IK01.GE.37) GOTO 71
        TOBR=2.4
      IF(AMAS(IK01).GT.0.14) GO TO 71
      IF(IK02.LE.36.AND.AMAS(IK02).GT.0.14) GO TO 71
      IF(IK02.GT.38) GO TO 71
      IQSUM = IQLE(IK01)+IQLE(IK02)
      IF(IQSUM.EQ.-2.OR.IQSUM.EQ.3) GO TO 71
      IF(IK01.LE.36.AND.IK02.LE.36.AND.
     *IQSUM.EQ.2) GO TO 71
C     IF(RNDM(-1.).GT.SIGEX/(SIGEL+SIGEX)) GO TO 79
C        IK1=7
C        IK2=38
C     IF(IK01.EQ.2.AND.IK02.EQ.37) GO TO 75
C        IK1=1
C        IK2=38
C     IF(IK01.EQ.7.AND.IK02.EQ.37) GO TO 75
C        IK1= 2
C        IK2=37
C     IF(IK01.EQ.7.AND.IK02.EQ.38) GO TO 75
C        IK1= 2
C        IK2=37
C     IF(IK01.EQ.1.AND.IK02.EQ.38) GO TO 75
C        IK1=IK01
C        IK2=IK02
C75   IEXE=1
 79   IF(TKIN.GT.2.5) GO TO 71
      ISOB=0
      P1(1)=P(1)
      P1(2)=P(2)
      P1(3)=P(3)
      P2(1)=-P(1)
      P2(2)=-P(2)
      P2(3)=-P(3)
      IF(IK01.GT.36.OR.IK02.GT.36) GO TO 60
      if(W > 0.281)
     & CALL FOROM(IK1,P1,AMN1,IK2,P2,AMN2,SIGEL,
     *IKD,PXD,PYD,PZD,DMAS,ISOB)
      GO TO 61
60    CALL FISOB(IK1,P1,AMN1,IK2,P2,AMN2,SIGEL,
     *IKD,PXD,PYD,PZD,DMAS,ISOB)
61    IF(ISOB.EQ.0) GO TO 71
      NP=1
      IORD(NP)=0
      IPR(NP)=IKD
      PPMAS(NP)=DMAS
      PR(1,NP)=PXD
      PR(2,NP)=PYD
      PR(3,NP)=PZD
      PR(4,NP)=SQRT(DMAS**2+PXD**2+PYD**2+PZD**2)
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      RETURN
 71   IF(TKIN-TOBR) 72,72,73
 72   CALL ELZPLE(IK01,IK02,TKIN,Z,PHI,IEXE)
C !!!!!!!!!!!!!!!! corrected 11.10.00 !!!!!!!!!!!
      IF( (IK01.EQ.10.OR.IK01.EQ.11.OR.IK01.EQ.16.OR.
     &     (IK01.GT.45.and.IK01.LE.48)).AND.
     &    (IK02.EQ.37.OR.IK02.EQ.38))  then
         AM1NEW=AMAS(IK01)
         if(AM1NEW.LT.(W-HB))  AMN1=AM1NEW
         HA=AMN1
         HA2=HA*HA
      ENDIF
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      GO TO 74
 73   CONTINUE
C   SLOPE CALCULATES THE ELASTIC SLOPES FOR THE CHOSEN MASSES
      CALL SLOPE(B,BFOR)
C   ANG CALCULATES THE TWO-BODY SCATTERING ANGLES (AZIMUTHAL ANGLE PHI
C   AND POLAR ANGLE THETA,WHERE Z=COS(THETA)
      CALL ANG(TFOR,TBACK,T,Z,PHI)
 74   CONTINUE
c
      if(ALAMB(S,HA2,HB2) < 0.0.or.(1.-Z**2) < 0.0)  then
        write( *,*) 'IK01,IK02,AMA,AMB,HA,HB,W,Z=',
     &                  IK01,IK02,AMA,AMB,HA,HB,W,Z
        stop
      endif
c
      PAMOD=SQRT(ALAMB(S,HA2,HB2))/(2.0*W)
      PAN=PAMOD*SQRT(1.-Z**2)
      PA(1)=PAN*COS(PHI)
      PA(2)=PAN*SIN(PHI)
      PA(3)=PAMOD*Z
      EA=SQRT(PAMOD**2+HA**2)
      EB=SQRT(PAMOD**2+HB**2)
      NP=NP+1
      IPR(NP)=IK1
      PR(1,NP)=PA(1)
      PR(2,NP)=PA(2)
      PR(3,NP)=PA(3)
      PR(4,NP)=EA
      IORD(NP)=IRD1
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      PPMAS(NP)=AMN1
      NP=NP+1
      IPR(NP)=IK2
      PR(1,NP)=-PA(1)
      PR(2,NP)=-PA(2)
      PR(3,NP)=-PA(3)
      PR(4,NP)=EB
      IORD(NP)=IRD2
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PPMAS(NP)=AMN2
      PLDER(NP)=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FOROM(IK01,P01,AM01,IK02,P02,AM02,SIGEL,
     *IKD,PXD,PYD,PZD,DMAS,ISOB)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  FORMATION OF  RHO,OMEGA,PHI AND K* MESONS
      DIMENSION P01(3),P02(3),P1(3),P2(3)
      ISOB=0
      E1=SQRT( SPQ(P01,P01)+AM01**2)
      E2=SQRT( SPQ(P02,P02)+AM02**2)
      S=AM01**2+AM02**2+2.*E1*E2-2.* SPQ(P01,P02)
      PXC=SQRT(ALAMB(S,AM01**2,AM02**2))/(2.*SQRT(S))
      PT=5.067*PXC
      DM=SQRT(S)
      IK1=IK01
      AM1=AM01
      IK2=IK02
      AM2=AM02
      DO 1 J=1,3
      P1(J)=P01(J)
      P2(J)=P02(J)
 1    CONTINUE
      Q1=CHARGE(IK1)
      Q2=CHARGE(IK2)
      IS1=IS(IK1)
      IS2=IS(IK2)
      QS=Q1+Q2
      IF(IS1.NE.0.OR.IS2.NE.0) GO TO 10
      IF(Q1.EQ.0..AND.Q2.EQ.0.) RETURN
C  RHO OR OMEGA MESONS FROM PIONS
      IF(QS) 3,4,5
 3    IKD=-121
      GO TO 7
 4    IKRHO=111
      SIGRHO=SGRO(IKRHO,SIGEL,DM,PT)
      IKOME=221
      SIGOME=SGRO(IKOME,SIGEL,DM,PT)
      IKD=IKRHO
      IF(SIGRHO.LT.SIGOME) IKD=IKOME
      GO TO 7
 5    IKD=121
      GO TO 7
10    IF(IS1.NE.0.AND.IS2.NE.0) GO TO 20
C   K* MESONS FROM PI AND K
C PI0
      IF(IS1.EQ.0.AND.Q1.EQ.0.) GO TO 21
      IF(IS2.EQ.0.AND.Q2.EQ.0.) GO TO 22
C PI+
      IF(IS1.EQ.0.AND.Q1.EQ.1.) GO TO 23
      IF(IS2.EQ.0.AND.Q2.EQ.1.) GO TO 24
C PI-
      IF(IS1.EQ.0.AND.Q1.EQ.-1.) GO TO 25
      IF(IS2.EQ.0.AND.Q2.EQ.-1.) GO TO 26
      RETURN
21    IKD=IABS(IK2)+1
      IKD=ISIGN(IKD,IK2)
      GO TO 7
22    IKD=IABS(IK1)+1
      IKD=ISIGN(IKD,IK1)
      GO TO 7
23    IF(IK2.EQ.-230) RETURN
      IKD=131
      IF(IK2.EQ.230) GO TO 7
      IKD=-231
      GO TO 7
24    IF(IK1.EQ.-230) RETURN
      IKD=131
      IF(IK1.EQ.230) GO TO 7
      IKD=-231
      GO TO 7
25    IF(IK2.EQ.230) RETURN
      IKD=231
      IF(IK2.EQ.130) GO TO 7
      IKD=-131
      GO TO 7
26    IF(IK1.EQ.230) RETURN
      IKD=231
      IF(IK1.EQ.130) GO TO 7
      IKD=-131
      GO TO 7
C  ******* SIVOKL'S CHANGES ******************************
20    IF(IK1+IK2.EQ.0) GO TO 27
C  *******************************************************
C20    IF(IABS(IK1).EQ.130.AND.IABS(IK2).NE.130) GO TO 27
C      IF(IABS(IK1).EQ.230.AND.IABS(IK2).NE.230) GO TO 27
      RETURN
C PHI FROM K+ AND K- OR K0 AND ANTK0
27    IKD=331
C  COMPUTE RESONANCE CROSS SECTION
 7    SIGR=SGRO(IKD,SIGEL,DM,PT)
      PR=SIGR/SIGEL
      IF(RNDM(-1.).GE.PR) RETURN
C  RESONANCE-MESON PARAMETERS
      PXD=P1(1)+P2(1)
      PYD=P1(2)+P2(2)
      PZD=P1(3)+P2(3)
      DMAS=DM
      ISOB=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SGRO(IKR,SIGEL,DM,PT)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF RESONANCE CROSS SECTION
      PT2=PT**2
      DM0=AMASS(IKR)
      GM=GAMHE(IKR)
      DMM0=(DM**2-DM0**2)**2
      DMG=(DM0*GM)**2
      ANORM=SIGEL*PT2
      SGRO=ANORM*DMG/(PT2*(DMG+DMM0))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION WIDTLE(GAM)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE WIDTH OF PARTICLE
C
      SIGMA=GAM
100   DRND=RNDM(-1.)
      IF(DRND.LT.1.0D-10)      GO TO 100
      GT=SIGMA*SQRT(-LOG(DRND))
      PHI=twpi*RNDM(-1.)
      WIDTLE=GT*COS(PHI)
      IF(ABS(WIDTLE).GT.GAM) GO TO 100
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SLOPEB(IB1,IB2,PL,B)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C  COMPUTE SLOPE IN TWO-BODY REACTION
C
      COMMON /CALC/HA,HB,HA2,HB2
      DATA BM/ 3.0/,BB/ 3.0/
      IF(IB1.EQ.0.OR.IB2.EQ.0) GO TO 100
      B0=BB
      GO TO 200
100   B0=BM
200   B=B0+0.7*LOG(PL)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        SUBROUTINE XDIST2(X1,X2)
      use modifiedDCMParams, only: pi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C         U(X)=1./SQRT(X1*X2)*DELTA(1.-X1-X2)
          X1=0.5+0.5*SIN(PI*(RNDM(-1.)-0.5))
          X2=1.-X1
          RETURN
           END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE XCORLE(IFL1,IFL2,KS1,PX1,PY1,PX2,PY2,X1,X2,
     *PSIGN,NPRODS,RETU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   CORRECT X-FRACTION AND DECAY STRING OR CLUSTER
C
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /COMLD/ PLDER(50)
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
       DIMENSION V(3)
       DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON/KAPPA/ XAP
      COMMON/PRIMP0/ P0
      LOGICAL RETU
      COMMON/COLRET/LRET
      COMMON /PRODMA/ PPMAS(50)
      LOGICAL LRET
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL  QSEE,QVSEE
      LRET = .FALSE.
C     INITIALIZE
      NPRODS=0
      NRET=0
      AMQ21=0.
      AMQ22=0.
      NFIX=NP
      RETU=.FALSE.
      PTX=PX1+PX2
      PTY=PY1+PY2
      PT12=PX1**2+PY1**2
      PT22=PX2**2+PY2**2
      P1=X1*P0
      P2=X2*P0*PSIGN
      E12=P1**2+PT12+AMQ21
      IF(E12.GE.0.) GO TO 200
      RETU=.TRUE.
      RETURN
200   E22=P2**2+PT22+AMQ22
      IF(E22.GE.0.) GO TO 210
      RETU=.TRUE.
      RETURN
210   E1=SQRT(E12)
      E2=SQRT(E22)
      AMSS12=(E1+E2)**2-(P1+P2)**2-PTX**2-PTY**2
       IKHR1=KI2(IFL1,IFL2,KS1,2)
      PARBE=0.2
      IF(IABS(IFL1).EQ.3.OR.IABS(IFL2).EQ.3) PARBE=0.3
       AMHR =AMAS(IKHR1)
       AMHRB=AMHR+PARBE
      IKHR=KI2(IFL1,IFL2,KS1,1)
      AMHR1=AMAS(IKHR)
       IF(AMSS12.GE.AMHRB**2) GO TO 400
      IF(NRET.EQ.1) GO TO 420
      NP=NP+1
      NPRODS=1
      PR(1,NP)=PTX
      PR(2,NP)=PTY
      PR(3,NP)=P1+P2
      PR(4,NP)=E1+E2
      PPMAS(NP)=AMHR
      IPR(NP)=IKHR1
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=PR(4,NP)/XAP
      AMT=SQRT(PTX**2+PTY**2+AMHR**2)
      PR(8,NP)=SQRT(2.D0)*AMT/XAP*PR(4,NP)/AMHR
      PLDER(NP)=1.
      IF(QVSEE) PLDER(NP)=0.
      IF(AMSS12.GE.AMHR**2)  GO  TO  419
      PR(4,NP)=SQRT(AMHR**2+PR(1,NP)**2+PR(2,NP)**2+PR(3,NP)**2)
419   RETURN
420    RETU=.TRUE.
       RETURN
400   AMSS1=SQRT(AMSS12)
      PZ1=P1+P2
      ESS1=E1+E2
      V(1)=PTX/ESS1
      V(2)=PTY/ESS1
      V(3)=PZ1/ESS1
      NIN1=NP+1
      IF(AMSS1.GE.AMHR1+SWMAX) GO TO 600
      CALL CLUSLE(IFL1,IFL2,KS1,AMSS1)
      IF(LRET) GO TO 800
      NFIN1=NP
      CALL TIFILE(NIN1,NFIN1,AMSS1)
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=P1
      CALL LORLLE(V,PPX1,E1,1)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 710 J=NIN1,NFIN1
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,1)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
710   CONTINUE
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
      NPRODS=NP-NFIX
      RETURN
600    CALL STRILE(IFL1,IFL2,KS1,AMSS1)
      IF(LRET) GO TO 800
       NFIN1=NP
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=P1
      CALL LORLLE(V,PPX1,E1,1)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 700 J=NIN1,NFIN1
      PPX1(1)=PR(1,J)
      PPX1(2)=PR(2,J)
      PPX1(3)=PR(3,J)
      CALL ROTAM(CT,ST,CFI,SFI,PPX1,PPX2,1)
      PR(1,J)=PPX2(1)
      PR(2,J)=PPX2(2)
      PR(3,J)=PPX2(3)
      PRX1(1)=PR(5,J)
      PRX1(2)=PR(6,J)
      PRX1(3)=PR(7,J)
      CALL ROTAM(CT,ST,CFI,SFI,PRX1,PRX2,1)
      PR(5,J)=PRX2(1)
      PR(6,J)=PRX2(2)
      PR(7,J)=PRX2(3)
700   CONTINUE
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
      NPRODS=NP-NFIX
       RETURN
800   CONTINUE
      RETU = .TRUE.
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION XSEE(XMIN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C   SIMULATE U(X)=1/X
C
      XSEE=XMIN**RNDM(-1.)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DIFTLE(IFL01,IFL02,KS01,IK1,AM1,
     *                  IFL03,IFL04,KS02,IK2,AM2,P1,IBINA)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE TRIPLE POMERON VERTEX DIFFRACTION
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      LOGICAL RETU
      COMMON /ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /DATA2/PUD,PS1,SIGMA,CX2
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /PRODT/ IORD(50)
      COMMON /ORDER/ IRD1,IRD2
      COMMON /PRODMA/ PPMAS(50)
      COMMON /COMCUT/ PARM,PARB,SWMAX
      COMMON /COLRET/ LRET
      COMMON /COMLD/PLDER(50)
      COMMON/PRIMP0/ P0
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/NPTCLZ/ NPTCLZ
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL  QSEE,QVSEE
      LOGICAL LRET
      DIMENSION V(3),P1(3),PSUM(5)
*      DATA SIGMAD/0.45/   ! IT IS NOT USED ?? 11.16.94  V.T.
      ZER=0.0
      LRET = .FALSE.
      QVSEE = .TRUE.
      NREP=0
C     INITIALIZE
C----- DON'T CHANGE ! (?)
      PARBE=0.15
C     NPI=NP
      NPI=0
      P0=SQRT(SPQ(P1,P1))
      E1=SQRT(P0**2+AM1**2)
      E2=SQRT(P0**2+AM2**2)
      ECM=E1+E2
      SCM=ECM**2
      IF(ECM-AM1-AM2.GT.0.3) GO TO 856
      LRET = .TRUE.
      IBINA=1
      RETURN
856   CONTINUE
      DO  96 I=1,3
   96 PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      XMIN=0.
      P0OLD=P0
100   CONTINUE
C      NP=NPI
      PSIGN=-1.
      PSOR=-1.
      IKA=IK1
      IKB=IK2
      IRDH=IRD2
      AMH=AM2
      AMB=AM1
      IFL1=IFL01
      IFL2=IFL02
      KS1=KS01
      EP=E1
      IF(RNDM(-1.).GT.0.5) GO TO 150
      IKA=IK2
      IKB=IK1
      IRDH=IRD1
      AMH=AM1
      AMB=AM2
      IFL1=IFL03
      IFL2=IFL04
      KS1=KS02
      EP=E2
      PSIGN=1.
      PSOR=1.
150   XMINS=(PARBE+AMB)**2/SCM
      P0=P0OLD
      NREP = NREP + 1
      IF(NREP.LT.NTRIES) GO TO 102
C     WRITE(ITLIS,101)IK1,IK2,AM1,AM2,P0
101   FORMAT(1X,'DIFTLE:NREP > NTRIES,IK1,IK2,AM1,AM2,P0=',
     *2I4,3(1X,F7.3))
      IBINA=1
      LRET = .TRUE.
      RETURN
102   CONTINUE
C   COMPUTE X VALUE FOR SEE QUARKS
      XS=XSEE(XMINS)
C    COMPUTE PT VALUE FOR HADRON
      CALL GETPT(PTH,SIGMA)
      PHI=twpi*RNDM(-1.)
      PTHX=PTH*COS(PHI)
      PTHY=PTH*SIN(PHI)
      PS=XS*P0
      ES=PS
      AMD2=(EP+ES)**2-(P0-PS)**2
      AMD=SQRT(AMD2)
      IF(ECM.LE.AMD+AMH)         GO TO 150
      ALA=ALAMB(SCM,AMD2,AMH**2)
      P0H=SQRT(ALA)/(2.0*ECM)
      DTRM=P0H**2-PTH**2
      IF(DTRM.LT.0.)             GO TO 150
      PZH=SQRT(DTRM)*PSIGN
      EH=SQRT(AMH**2+P0H**2)
      ED=SQRT(AMD2+P0H**2)
      V(1)=-PTHX/ED
      V(2)=-PTHY/ED
      V(3)=-PZH/ED
C    COMPUTE X VALUES FOR PARTONS
170   IFLS1=1+INT(RNDM(-1.)/PUD)
      IFLS2=IFLS1
      IBA=IBLE(IKA)
      ISA=ISLE(IKA)
      IFL1S=IFL1
      IF(ISA.NE.0.AND.IBA.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      XQVAL=XDIST(XMIN,IBA,ISA,IFL1S)
      XQQVA=1.-XQVAL
C   COMPUTE X VALUE FOR PION
      CALL XDIST2(XPI,X2PI)
C     NIN1=NP+1
      NIN1=1
      IF(IFL2.GT.0) GO TO 160
      IFL1T=IFL1
      IFLS1T=-IFLS1
      XQVALT=XQVAL
      XPIT=XPI
      IFL2T=IFL2
      IFLS2T=IFLS2
      XQQVAT=XQQVA
      XPIT1=1.-XPI
      PSI=1.
      IF(PSIGN.LT.0.) GO TO 400
      IFL1T=-IFLS1
      IFLS1T=IFL1
      XQVALT=XPI
      XPIT=XQVAL
      IFL2T=IFLS2
      IFLS2T=IFL2
      XQQVAT=1.-XPI
      XPIT1=XQQVA
      PSI=-1.
400   PSIGN=PSIGN*PSI
      P0=AMD/2.0
      KS1T=0
      IF(IABS(IFL1T).GT.3.OR.IABS(IFLS1T).GT.3) KS1T=KS1
      CALL XCORLE(IFL1T,IFLS1T,KS1T,ZER,ZER,ZER,ZER,XQVALT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 170
      KS2T=0
      IF(IABS(IFL2T).GT.3.OR.IABS(IFLS2T).GT.3) KS2T=0
      CALL XCORLE(IFL2T,IFLS2T,KS2T,ZER,ZER,ZER,ZER,XQQVAT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NP=NP-NPRD
      GO TO 170
130   NFIN1=NP
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
      GO TO 300
160   IFL1T=IFL1
      IFLS2T=ISIGN(IFLS2,-IFL1)
      XQVALT=XQVAL
      XPIT1=1.-XPI
      IFL2T=IFL2
      IFLS1T=IFLS1
      IF(IFL2.LE.3)  IFLS1T=ISIGN(IFLS1,-IFL2)
      XQQVAT=XQQVA
      XPIT=XPI
      PSI=1.
      IF(PSIGN.LT.0.) GO TO 450
      IFL1T =ISIGN(IFLS2,-IFL1)
      IFLS2T=IFL1
      XQVALT=1.-XPI
      XPIT1=XQVAL
      IFLS1T=IFL2
      IFL2T=IFLS2
      IF(IFL2.LE.3)  IFL2T =ISIGN(IFLS2,-IFL2)
      XQQVAT=XPI
      XPIT=XQQVA
      PSI=-1.
450   PSIGN=PSIGN*PSI
      P0=AMD/2.0
      KS1T=0
      IF(IABS(IFL1T).GT.3.OR.IABS(IFLS2T).GT.3) KS1T=KS1
      CALL XCORLE(IFL1T,IFLS2T,KS1T,ZER,ZER,ZER,ZER,XQVALT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 170
      KS2T=0
      IF(IABS(IFL2T).GT.3.OR.IABS(IFLS1T).GT.3) KS2T=KS1
      CALL XCORLE(IFL2T,IFLS1T,KS2T,ZER,ZER,ZER,ZER,XQQVAT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 140
      NP=NP-NPRD
      GO TO 170
140   NFIN1=NP
      CALL LORPLE(V,NIN1,NFIN1,-1)
      CALL LORCLE(V,NIN1,NFIN1,-1)
300   P0=P0OLD
      DO 500 I=NIN1,NFIN1
      IORD(I)=0
C     PPMAS(I)=AMAS(IPR(I))
500   CONTINUE
      NP=NP+1
      NFIN1=NP
      NPRD=NP-NPI
      PR(1,NP)=PTHX
      PR(2,NP)=PTHY
      PR(3,NP)=PZH
      PR(4,NP)=EH
      IPR(NP)=IKB
      IORD(NP)=IRDH
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      PPMAS(NP)=AMH
      NPTCL=0
      DO  502  J=NIN1,NFIN1
      NPTCL=NPTCL+1
      PPTCL(1,J)=PR(1,J)
      PPTCL(2,J)=PR(2,J)
      PPTCL(3,J)=PR(3,J)
      PPTCL(4,J)=PR(4,J)
      PPTCL(5,J)=PPMAS(J)
  502 CONTINUE
      CALL  RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0)   GO  TO  501
      NP=NP-NPRD
      GO  TO  100
  501 CONTINUE
      DO  503  J=NIN1,NFIN1
      PR(1,J)=PPTCL(1,J)
      PR(2,J)=PPTCL(2,J)
      PR(3,J)=PPTCL(3,J)
      PR(4,J)=PPTCL(4,J)
      PPMAS(J)=PPTCL(5,J)
  503 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BINALE(P,AMA,AMB,IK1,IK2,IEL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  SIMULATION TWO PARTICLE REACTION
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /COMKI1/ HLA2,HLB2,W,INUMA
      COMMON /COMKI2/ELA,ELB,PLALB
      COMMON /CALC/HA,HB,HA2,HB2
      COMMON /BEES/B,BFOR
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON /PRODT/ IORD(50)
      COMMON /PRODMA/PPMAS(50)
      COMMON /COMLD/PLDER(50)
      COMMON /COLRET/ LRET
      COMMON /ISOB3/ISOB3
      DIMENSION P(3),PA(3)
      LOGICAL LRET
      LRET = .FALSE.
C===== IF DL- + DL- , DL++ + DL++OR FI0 + OM- =====C
      IF((IK1.EQ.47.AND.IK2.EQ.47).OR.
     *   (IK1.EQ.18.AND.IK2.EQ.18).OR.
     *   (IK1.EQ.54.AND.IK2.EQ.54).OR.
     *   (IK1.EQ.45.AND.IK2.EQ.45).OR.
     *   (IK1.EQ.18.AND.IK2.EQ.54)) IEL=1
      IF(IEL.EQ.1) RETURN
C  INITIALIZE
      INUMA=0
      IREP2 =0
      E1=SQRT(AMA**2+SPQ(P,P))
      E2=SQRT(AMB**2+SPQ(P,P))
      S=(E1+E2)**2
      W=SQRT(S)
      IB1=IBLE(IK1)
      IB2=IBLE(IK2)
      HLA=AMA
      HLB=AMB
      HLA2=HLA*HLA
      HLB2=HLB*HLB
      PLALB=SQRT(ALAMB(S,HLA2,HLB2))/(2.0*W)
      ELA=(S+HLA2-HLB2)/(2.0*W)
      ELB=(S+HLB2-HLA2)/(2.0*W)
C
c !!!      20.10.92T     DL+N->N+N'
      IF(W.GT.3.00.OR.(IB1+IB2).NE.2)     GO TO 99
      IF(AMA.GT.AMB)  THEN
        IF(.NOT.(IK2.EQ.37.OR.IK2.EQ.38))   GO TO 99
        IF(.NOT.(IK1.EQ.46.OR.IK1.EQ.45.OR.IK1.EQ.47.OR.
     *  IK1.EQ.48))                         GO TO 99
      ELSE
        IF(.NOT.(IK1.EQ.37.OR.IK1.EQ.38))   GO TO 99
        IF(.NOT.(IK2.EQ.46.OR.IK2.EQ.45.OR.IK2.EQ.47.OR.
     * IK2.EQ.48))                         GO TO 99
      END IF
      ICHA1=IQLE(IK1)+IQLE(IK2)+1
      IF(ICHA1.EQ.4.OR.ICHA1.EQ.0)        GO TO 99
      GOTO (95,96,97), ICHA1
  95  IKH1=38
      IKH2=38
      GO  TO  98
  96  IKH1=38
      IKH2=37
      GO  TO  98
  97  IKH1=37
      IKH2=37
  98  AMP1=AMASF(IKH1)
      AMP2=AMASF(IKH2)
      GO  TO  109
  99  CONTINUE
c !!!
C   CHOOSE INTERACTIVE QUARKS
105   CALL FLAVO(IB1,IK1,IB2,IK2,IFL1,IFL2,IFL3,IFL4)
      IREP2 = IREP2 + 1
      IF(IREP2.LE.NTRIES)   GOTO 305
C     WRITE(ITLIS,1200)IK1,IK2,PLALB
1200  FORMAT(1X,'BINALE:NREP > NTRIES,IK1,IK2,PLALB=',
     *2I4,1X,F7.3)
      IEL=1
C     LRET=.TRUE.
      RETURN
305   CONTINUE
C   CHOOSE DIQUARK SPIN
      CALL KSPIN(IFL3,KS2,IK2,IB2)
      CALL KSPIN(IFL1,KS1,IK1,IB1)
      IREP1=0
C  HADRONS GENERATE BY MEANS QUARKS EXCHANGE
      IF(IFL1)1,1,2
1     IF11=IFL2
      IF22=IFL1
      GO TO 3
2     IF11=IFL1
      IF22=IFL2
3     CONTINUE
      IF(IB2.EQ.1) GO TO 102
      IF(IFL3) 101,101,102
101   IF44=IFL3
      IF33=IFL4
      GO TO 103
102   IF33=IFL3
      IF44=IFL4
103   CONTINUE
104   IKH2=KI2(IF44,IF11,KS2,0)
      IKH1=KI2(IF22,IF33,KS1,0)
      IREP1=IREP1+1
      IF(IREP1.GT.NTRIES) GO TO 105
C  SELECT ELASTIC COLLISION
      IF(IKH1.EQ.IK1.AND.IKH2.EQ.IK2.AND.IREP1.LE.NTRIES)
     *GO  TO  104
C  CHOOSE TABLE MASSES AND TABLE WIDTH OF HADRONS
      AMH1=AMASF(IKH1)
      AMH2=AMASF(IKH2)
      GAMH1=GAM(IKH1)
      GAMH2=GAM(IKH2)
C  COMPUTE MASSES OF PARTICLES
      IREP3=0
205   CONTINUE
      IREP3=IREP3+1
      IF(IREP3.LE.NTRIES) GO TO 108
C     WRITE(ITLIS,1200)IK1,IK2,PLALB
      IEL=1
      RETURN
108   GAM1=WIDTLE(GAMH1)
      GAM2=WIDTLE(GAMH2)
      AMP1=AMH1+GAM1
      AMP2=AMH2+GAM2
      IF(ISOB3 == 1) THEN
         IF(IKH1.GE.45.AND.IKH1.LE.48)               AMP1=AMAS(IKH1)
         IF(IKH2.GE.45.AND.IKH2.LE.48)               AMP2=AMAS(IKH2)
         IF(IKH1.EQ.10.OR.IKH1.EQ.11.OR.IKH1.EQ.16)  AMP1=AMAS(IKH1)
         IF(IKH2.EQ.10.OR.IKH2.EQ.11.OR.IKH2.EQ.16)  AMP2=AMAS(IKH2)
      end if
C  CHECK ENERGY THRESHOLD
107   IF(W.LT.(AMP1+AMP2)) GO TO 205
109   HA=AMP1
      HA2=AMP1**2
      HB=AMP2
      HB2=AMP2**2
C  COMPUTE SCATTERING ANGLE
      CALL SLOPEB(IB1,IB2,PLALB,B)
      CALL ANG(TFOR,TBACK,T,Z,PHI)
      PAMOD=SQRT(ALAMB(S,HA2,HB2))/(2.0*W)
      PAN=PAMOD*SQRT(1.-Z**2)
      PA(1)=PAN*COS(PHI)
      PA(2)=PAN*SIN(PHI)
      PA(3)=PAMOD*Z
      NP=NP+1
      IPR(NP)=IKH1
      PR(1,NP)=PA(1)
      PR(2,NP)=PA(2)
      PR(3,NP)=PA(3)
      PR(4,NP)=SQRT(PAMOD**2+HA2)
      IORD(NP)=0
      PPMAS(NP)=AMP1
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      NP=NP+1
      IPR(NP)=IKH2
      PR(1,NP)=-PA(1)
      PR(2,NP)=-PA(2)
      PR(3,NP)=-PA(3)
      PR(4,NP)=SQRT(PAMOD**2+HB2)
      IORD(NP)=0
      PPMAS(NP)=AMP2
      PR(5,NP)=0.
      PR(6,NP)=0.
      PR(7,NP)=0.
      PR(8,NP)=0.
      PLDER(NP)=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RESCAL(N1,N2,PSUM,IFAIL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C          RESCALE MOMENTA OF PARTICLES N1...N2 TO GIVE TOTAL
C          FOUR-MOMENTUM PSUM.
C          RETURN IFAIL=0 IF OK, IFAIL=1 IF NO GOOD.
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      DIMENSION PSUM(5),PADD(5),BETA(3)
      DATA ERRLIM/.00001/
C          ORIGIONAL MOMENTUM IS PADD.
      IFAIL=1
      IF(N1.GE.N2) RETURN
      DO 100 K=1,5
100   PADD(K)=0.
      DO 110 IP=N1,N2
      DO 110 K=1,5
      PADD(K)=PADD(K)+PPTCL(K,IP)
110   CONTINUE
      IF(PADD(5).GE.PSUM(5)) RETURN
      PADD(5)=PADD(4)**2-PADD(1)**2-PADD(2)**2-PADD(3)**2
      IF(PADD(5).LE.0) RETURN
      PADD(5)=SQRT(PADD(5))
      DO 120 K=1,3
120   BETA(K)=-PADD(K)/PADD(5)
      GAMMA=PADD(4)/PADD(5)
C          BOOST PARTICLES TO REST.
200   CONTINUE
      DO 210 IP=N1,N2
      BP=0.
      DO 220 K=1,3
220   BP=BP+PPTCL(K,IP)*BETA(K)
      DO 230 K=1,3
230   PPTCL(K,IP)=PPTCL(K,IP)+BETA(K)*PPTCL(4,IP)
     $+BETA(K)*BP/(GAMMA+1.)
      PPTCL(4,IP)=GAMMA*PPTCL(4,IP)+BP
210   CONTINUE
      IF(IFAIL.EQ.0) RETURN
C          RESCALE MOMENTA IN REST FRAME.
      SCAL=1.
      DO 301 IPASS=1,500
      SUM=0.
      DO 310 IP=N1,N2
      DO 320 K=1,3
320   PPTCL(K,IP)=SCAL*PPTCL(K,IP)
      PPTCL(4,IP)=SQRT(PPTCL(1,IP)**2+PPTCL(2,IP)**2+PPTCL(3,IP)**2
     $+PPTCL(5,IP)**2)
      SUM=SUM+PPTCL(4,IP)
310   CONTINUE
      SCAL=PSUM(5)/SUM
      IF(ABS(SCAL-1.).LE.ERRLIM) GO TO 300
301   CONTINUE
300   CONTINUE
C          BOOST BACK WITH PSUM.
      BMAG=0.
      DO 400 K=1,3
      BETA(K)=PSUM(K)/PSUM(5)
      BMAG=BMAG+ABS(BETA(K))
400   CONTINUE
      GAMMA=PSUM(4)/PSUM(5)
      IFAIL=0
      IF(BMAG.EQ.0.) RETURN
      GO TO 200
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE COLLHH(ID01,AM01,PX01,PY01,PZ01,ID02,AM02,
     *PX02,PY02,PZ02,SIGTO0,SIGAN0,SIGEL0)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      real*8 :: psum(5)
C
C          MAIN SUBROUTINE FOR HHQGS, A MONTE CARLO EVENT GENERATOR
C          FOR  H+H COLLISION AT LOW AND HIGH ENERGY
C          JTDKY = +/- UNIT NUMBER FOR DECAY TABLE FILE.
C                      IF IT IS NEGATIVE, DECAY TABLE IS NOT PRINTED.
C          JTEVT =     UNIT NUMBER FOR OUTPUT EVENT FILE.
C          JTCOM =     UNIT NUMBER FOR COMMAND FILE.
C          JTLIS =     UNIT NUMBER FOR LISTING.
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMELX/ SIGEL
      COMMON/COMCRO/ SIGTOT
      COMMON/CSIGA/ SIGAN
      COMMON/ORDER/ IRD1,IRD2
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL  QSEE,QVSEE
      COMMON/PROD/PR(8,50),IPR(50),NP
      COMMON/PRIMP0/ P0
      COMMON/PRIMPL/ PL
      COMMON/COMFR/ ICMS
      COMMON/PRODT/ IORD(50)
      COMMON/PARORD/ IORDP(499)
      COMMON/PRODMA/ PPMAS(50)
      COMMON/COMLID/ PLIDER(499)
      COMMON/COMLD/ PLDER(50)
      COMMON/COMASS/ AM1N,AM2N
      LOGICAL GH1H2
      COMMON/H1H2/ GH1H2(11)
      COMMON/COMMUL/ MULTP
      LOGICAL MULTP
      DIMENSION P1(3),P2(3),V(3),P1R(3),P2R(3)
      COMMON/COMENB/ ENBOU
      COMMON/COLRET/ LRET
      COMMON /SIGDIA/ CROSD(5),DST
      COMMON/CPRSIG/ ISGCOL
      COMMON/PRINTS/IPRINT
      COMMON/P0LAB1/ P00,DSM,TKIN
!      COMMON/KEYEL/IEL
      LOGICAL IPRINT
      COMMON/KEYHH/KEYHH
      LOGICAL KEYHH
      COMMON/KEYPLA/KEYPLA
      LOGICAL KEYPLA
!      COMMON/COMCOR/ PSUM(5)
      COMMON/VALON/ VALON
      COMMON/HELE/ IHELE
      COMMON/NCASCA/NCAS,NCPRI
      LOGICAL VALON
      DIMENSION IK(5),PO(3),IRL(250),ILL(250)
C
      CHARACTER*8 LAB1,LAB2
      LOGICAL BACK
      LOGICAL LRET
      LRET=.FALSE.
      NREP=0
C          ENTRY.
      if(NCAS >= NCPRI) PRINT *, ' in COLLHH'
      DO  110 I=1,11
  110 GH1H2(I)=.FALSE.
      IB1=IB(ID01)
      IB2=IB(ID02)
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 12
      IF((ID01.EQ.1120.OR.ID01.EQ.1220) .AND.
     *    ID02.NE.1120.AND.ID02.NE.1220) GO TO 13
      IF(IB1.GT.0.AND.IB2.LE.0) GO TO 13
      IF(IB1.LT.0.AND.IB2.EQ.0) GO TO 13
C     IF INCIDENT PARTICLE IS MESON OR
C     IF TARGET PARTICLE IS NUCLEON OR
C     IF INCIDENT AND TARGET PARTICLES ARE NOT NUCLEONS
C
12    ID1=ID01
      PX1=PX01
      PY1=PY01
      PZ1=PZ01
      AM1=AM01
      IRD1=1
      ID2=ID02
      PX2=PX02
      PY2=PY02
      PZ2=PZ02
      AM2=AM02
      IRD2=2
      GO TO 14
C       ELSE DO EXCHANGE:
13    CONTINUE
      ID1=ID02
      PX1=PX02
      PY1=PY02
      PZ1=PZ02
      AM1=AM02
      IRD1=2
      ID2=ID01
      PX2=PX01
      PY2=PY01
      PZ2=PZ01
      AM2=AM01
      IRD2=1
C     DETERMINE THE COMMON/COMELX/ AND /CSIGA/ ELEMENTS
14    SIGTOT=SIGTO0
      SIGAN=SIGAN0
      SIGEL=SIGEL0
C
      BACK=.TRUE.
1912  NP=0
      IEL=0
      NREP = NREP + 1
      IF(NREP.LE.3*NTRIES) GO TO 1913
      LRET=.TRUE.
c     WRITE(ITLIS,1917)ID01,ID02,PX01,PY01,PZ01,PX02,PY02,PZ02
c    *,S,P00,SIGTOT,SIGEL,SIGAN,CROSD
1917  FORMAT(5X,' STOP IN COLLHH :ID01,ID02,PX,Y,Z,1-2=',
     *2I6,6F8.3/,' S=',E10.4,' P0=',E10.4,' SIGTOT=',E10.4,
     *' SIGEL',E10.4,' SIGAN=',E10.4,/,' CROSD(5)=',5E10.4)
      STOP                          ! 11.06.2001 G.
C   COMPUTE OF CM MOMENTUM COMPONENTS
 1913 CALL KINEM(PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,
     *AM2,V,S,P1)
C     IN CM SYSTEM:
      DO  9 I=1,3
  9      PO(I)=P1(I)
C     IN LAB SYSTEM:
      PL=SQRT(ALAMB(S,AM1**2,AM2**2))/(2.0*AM2)
C     ROTATE CM SYSTEM
      CALL ANGLE(P1,CT,ST,CFI,SFI)
      CALL ROTR(CT,ST,CFI,SFI,P1,P2,BACK)
      IB1=IB(ID1)
      IB2=IB(ID2)
      SQS=SQRT(S)
 120  IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 210
      ITOT=0
      IKS=0
      IF(SQS.GT.ENBOU) GO TO 210
      CALL PARCRO(ITOT,IK,IKS,ID1,ID2,PX1,PY1,
     *PZ1,AM1,PX2,PY2,PZ2,AM2)
      IF(SQS.GT.ENBOU.OR.(IB1.LT.0.OR.IB2.LT.0.)) GO TO 210
C
C  COMPUTE LOW ENERGY HADRON-BARION COLLISION =============
C       (QUARK-GLUON STRING MODEL)
C
      IHELE=1
      CALL BACKID(ID1,ID1N)
      CALL BACKID(ID2,ID2N)
      IDIN(1)=ID01
      IDIN(2)=ID02
220   CALL TYPRE(ID1N,AM1,ID2N,AM2,
     *SQS,IEL,IDIF,IUNCY,IPLAN,IBINA,SIGTOT,SIGEL)
      DO 10 I=1,3
 10      P1(I)=P2(I)
      QVSEE=.FALSE.
C
      IF(LRET) IEL=1
C
      IF(IEL.EQ.1) GO TO 133
      IF(IBINA.EQ.1) GO TO 135
      IF(IDIF.EQ.1) GO TO 132
      IF(IPLAN.EQ.0) GO TO 100
C
C  CALCULATION OF PLANAR GRAPH -------------------------
C     NIN1=NP+1
      NIN1=1
      if(NCAS >= NCPRI) PRINT *, ' TO PLANLE'
      CALL PLANLE(ID1N,IB1,ID2N,IB2,S,IEL)
      if(NCAS >= NCPRI) PRINT *, 'from PLANLE'
      IF(IEL.EQ.1) GO TO 133
      IF(NP.EQ.2.AND.((IPR(1).EQ.ID1N.AND.IPR(2).EQ.ID2N)
     *.OR.(IPR(2).EQ.ID1N.AND.IPR(1).EQ.ID2N))) GO TO 1912
      IF(LRET) GO TO 1912
      GH1H2(2)=.TRUE.
      NFIN2=NP
      GO TO 300
100   CONTINUE
      IF(IUNCY.EQ.0) GO TO 104
C     NIN1=NP+1
      NIN1=1
C
C   UNDEVELOPED CYLINDER DIAGRAM -------------------
      if(NCAS >= NCPRI) PRINT *, ' TO UNCYLE'
      CALL UNCYLE(ID1N,IB1,AM1,ID2N,IB2,AM2,P1,IBINA)
      if(NCAS >= NCPRI) PRINT *, 'from UNCYLE'
      IF(IBINA.EQ.1) GO TO 135
      IF(LRET) GO TO 1912
      GH1H2(3)=.TRUE.
      NFIN2=NP
      GO TO 300
C
C   CYLINDER DIAGRAM -------------------------------
104   CONTINUE
C     NIN1=NP+1
      NIN1=1
      if(NCAS >= NCPRI) PRINT *, ' TO CYLLE'
      CALL CYLLE(ID1N,IB1,AM1,ID2N,IB2,AM2,P1,IBINA)
      if(NCAS >= NCPRI) PRINT *, 'from CYLLE'
      IF(IBINA.EQ.1) GO TO 135
      IF(LRET) GO TO 1912
      GH1H2(7)=.TRUE.
      NFIN2=NP
      GO TO 300
C
C  DIFRACTIVE SCATTERING CALCULATION ----------------
 132  CONTINUE
      CALL FLAVO(IB1,ID1N,IB2,ID2N,IFL1,IFL2,IFL3,IFL4)
      CALL KSPIN(IFL3,KS2,ID2N,IB2)
      CALL KSPIN(IFL1,KS1,ID1N,IB1)
C     NIN1=NP+1
      NIN1=1
      DRND=RNDM(-1.)
      IF(DRND.GT.DST)  GO  TO  222
      if(NCAS >= NCPRI) PRINT *, ' TO DIFTLE'
      CALL DIFTLE(IFL1,IFL2,KS1,ID1N,AM1,IFL3,IFL4,KS2,ID2N,AM2,P1,
     *IBINA)
      if(NCAS >= NCPRI) PRINT *, 'from DIFTLE'
      IF(IBINA.EQ.1) GO TO 135
      IF(LRET) GO TO 1912
      GH1H2(1)=.TRUE.
      GO  TO  223
 222  CONTINUE
      if(NCAS >= NCPRI) PRINT *, 'TO DIFSCA'
      CALL DIFSCA(IFL1,IFL2,KS1,ID1N,AM1,IFL3,IFL4,KS2,ID2N,AM2,P1,
     *IBINA)
      if(NCAS >= NCPRI) PRINT *, 'from DIFSCA'
      IF(IBINA.EQ.1) GO TO 135
      IF(LRET) GO TO 1912
      GH1H2(6)=.TRUE.
 223  CONTINUE
      NFIN2=NP
C     GO  TO  302
 300  CONTINUE
      DO 301 I=NIN1,NFIN2
C     PPMAS(I)=AMAS(IPR(I))
C     PPMAS(I)=SQRT(PR(4,I)**2-PR(1,I)**2-PR(2,I)**2-PR(3,I)**2)
 301  CONTINUE
      GO TO 302
C
C   TWO-PARTICLE REACTION --------------------------------
135   CONTINUE
      DO 1350 I=1,3
1350  P1(I)=P2(I)
C     NIN1=NP+1
      NIN1=1
      if(NCAS >= NCPRI) PRINT *, ' TO BINALE'
      CALL BINALE(P1,AM1,AM2,ID1N,ID2N,IEL)
      if(NCAS >= NCPRI) PRINT *, 'FROM BINALE'
      IF(LRET) GO TO 1912
      IF(IEL.EQ.1) GO TO 133
      GH1H2(8)=.TRUE.
      NFIN2=NP
      GO TO 302
C
C   ELASTIC SCATTERING ----------------------------------
133   CONTINUE
      DO 1330 I=1,3
1330  P1(I)=P2(I)
C     NIN1=NP+1
      NIN1=1
      if(NCAS >= NCPRI) PRINT *, ' TO ELASLE'
      CALL ELASLE(P1,AM1,AM2,ID1N,ID2N)
      if(NCAS >= NCPRI) PRINT *, 'FROM ELASLE'
      GH1H2(4)=.TRUE.
      NFIN2=NP
 302  CONTINUE
      GO TO 500
C
C   AT HIGH ENERGY=====================================
C
210   CONTINUE
      IHELE=2
      IRET=0
      NREP = NREP + 1
      IF(NREP.LE.3*NTRIES) GO TO 2913
      LRET=.TRUE.
c     WRITE(ITLIS,1905)ID01,ID02,PX01,PY01,PZ01,PX02,PY02,PZ02
c    *,S,PL,SIGTOT,SIGEL,SIGAN
1905  FORMAT(5X,' STOP IN COLLHH (HIGH EN.):ID01,ID02,PX,Y,Z,1-2=',
     *2I6,6F8.3/,' S=',E10.4,' P0=',E10.4,' SIGTOT=',E10.4,
     *' SIGEL',E10.4,' SIGAN=',E10.4)
      RETURN
2913  NPTCL=0
      NIN1=NPTCL+1
C FILL  /PRIMAR/
      IDIN(1)=ID1
      IDIN(2)=ID2
      AM1N=AM1
      AM2N=AM2
      SCM=S
      ECM=SQS
      P0=SQRT(SPQ(P1,P1))
      HALFE=ECM/2.0
      PSUM(1)=0.
      PSUM(2)=0.
      PSUM(3)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      CALL LABEL(LAB1,ID01)
      CALL LABEL(LAB2,ID02)
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,1000) LAB1,LAB2,SCM
1000  FORMAT(//15X,47HI SELECT THE NEXT REACTIONS WITH CROSS SECTIONS
     *,/20X,3HFOR,1X,A8,A8,18H COLLISION AT SCM=,E10.4,7H GEV**2/)
C
C       HH AND AHH COLLISIONS AT HIGH ENERGY AND ANTIBARION
C       BARION COLLISION OR MESON-MESON COLISION
C         (QUARK-GLUON STRINGS MODEL)
C-------
      IF(.not.KEYHH) then
         MULTP=.TRUE.
         QSEE=.FALSE.
      end if
C-------
      DO 399 JJ=1,11
 399  GH1H2(JJ)=.TRUE.
      CALL SIGIN
      CALL REACTL
      KEYPLA=.FALSE.
      IF(.NOT.GH1H2(1)) GO TO 600
      if(NCAS >= NCPRI) PRINT *, ' TO DIFTRI'
      CALL DIFTRI(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM DIFTRI'
      GO TO 707
600   IF(.NOT.GH1H2(8)) GO TO 700
      if(NCAS >= NCPRI) PRINT *, ' TO BINAR'
      CALL BINAR(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM BINAR'
      GO TO 707
700   IF(.NOT.GH1H2(10)) GO TO 701
      if(NCAS >= NCPRI) PRINT *, ' TO REGTRI'
      CALL REGTRI(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM REGTRI'
      GO TO 707
701   IF(.NOT.GH1H2(2)) GO TO 702
      if(NCAS >= NCPRI) PRINT *, ' TO PLANAR'
      CALL PLANAR(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM PLANAR'
      GO TO 707
702   IF(.NOT.GH1H2(3)) GO TO 703
      if(NCAS >= NCPRI) PRINT *, ' TO UNCYLY'
      CALL UNCYLI(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM UNCYLY'
      GO TO 707
703   IF(.NOT.GH1H2(4)) GO TO 704
      IEL=1
      if(NCAS >= NCPRI) PRINT *, ' TO ELAST'
      CALL ELASTQ(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM ELAST'
      GO TO 707
704   IF(.NOT.GH1H2(5)) GO TO 705
      if(NCAS >= NCPRI) PRINT *, ' TO ANNIH'
      CALL ANNIH(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM ANNIH'
      GO TO 707
705   IF(.NOT.GH1H2(6)) GO TO 706
      if(NCAS >= NCPRI) PRINT *, ' TO DIFSMA'
      CALL DIFSMA(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM DIFSMA'
      GO TO 707
706   IF(.NOT.GH1H2(7)) GO TO 709
      if(MULTP.and.NCAS >= NCPRI) PRINT *, ' TO CHAINS'
      IF(MULTP) CALL CHAINS(IRET)
      if(MULTP.and.NCAS >= NCPRI) PRINT *, 'FROM CHAINS'
      if(.NOT.MULTP.and.NCAS >= NCPRI) PRINT *, ' TO CYLIN'
      IF(.NOT.MULTP) CALL CYLIN(IRET)
      if(.NOT.MULTP.and.NCAS >= NCPRI) PRINT *, 'FROM CYLIN'
709   IF(.NOT.GH1H2(11)) GO TO 707
      if(NCAS >= NCPRI) PRINT *, ' TO DOUBDI'
      CALL  DOUBDI(IRET)
      if(NCAS >= NCPRI) PRINT *, 'FROM DOUBDI'
707   IF(IRET.EQ.1) GO TO 210
      GO TO 520
500   CONTINUE
      NPTCL=0
      DO 510 I=NIN1,NFIN2
      NPTCL=NPTCL+1
      CALL FORID(IPR(I),IDENT(NPTCL))
      PPTCL(1,NPTCL)=PR(1,I)
      PPTCL(2,NPTCL)=PR(2,I)
      PPTCL(3,NPTCL)=PR(3,I)
      PPTCL(4,NPTCL)=PR(4,I)
      PPTCL(5,NPTCL)=PPMAS(I)
      PPTCL(6,NPTCL)=PR(5,I)
      PPTCL(7,NPTCL)=PR(6,I)
      PPTCL(8,NPTCL)=PR(7,I)
      PPTCL(9,NPTCL)=PR(8,I)
      IDCAY(NPTCL)=0
      IORIG(NPTCL)=0
      IORDP(NPTCL)=IORD(I)
      PLIDER(NPTCL)=PLDER(I)
510   CONTINUE
C   BACKWARD ROTATE
520   CONTINUE
      NFIN2=NPTCL
      IF(NFIN2.EQ.1) GO TO 99
      NRL=0
      NLL=0
      DO 530 I=NIN1,NFIN2
      DO 540 J=1,3
      P2R(J)=PPTCL(5+J,I)
540   P2(J)=PPTCL(J,I)
C   BACKWARD ROTATE
      BACK=.FALSE.
      CALL ROTR(CT,ST,CFI,SFI,P2,P1,BACK)
      CALL ROTR(CT,ST,CFI,SFI,P2R,P1R,BACK)
      DO 550 J=1,3
      PPTCL(J+5,I)=P1R(J)
550   PPTCL(J,I)=P1(J)
      PPTCL(4,I)=SQRT(PPTCL(1,I)**2+PPTCL(2,I)**2+
     *PPTCL(3,I)**2+PPTCL(5,I)**2)
      SCAL=PO(1)*PPTCL(1,I)+PO(2)*PPTCL(2,I)+PO(3)*PPTCL(3,I)
      IF(VALON) GO TO 779
      IF(SCAL.GT.0.) NRL=NRL+1
      IF(SCAL.LE.0.) NLL=NLL+1
      IF(SCAL.GT.0.) IRL(NRL)=I
      IF(SCAL.LE.0.) ILL(NLL)=I
779   IF(.NOT.(IORDP(I).EQ.0.AND.PLIDER(I).GT..3)) GO TO 530
      IF(SCAL.GE.0.) IORDP(I)=IRD1
      IF(SCAL.LT.0.) IORDP(I)=IRD2
530   CONTINUE
      IF(NRL.EQ.0.OR.NLL.EQ.0) GO TO 99
      IF(VALON) GO TO 99
      IF(NRL.LE.1) GO TO 1537
      DO 537 I=1,NRL
      IF(PLIDER(IRL(I)).LT.0.1.OR.PLIDER(IRL(I)).GT.0.9) GO TO 537
      PMAX=PPTCL(4,IRL(I))
      IMAX=I
      PLI=PLIDER(IRL(I))
      GO TO 1538
537   CONTINUE
      GO TO 1537
1538  CONTINUE
      DO 533 J=1,NRL
      IF(PLIDER(IRL(J)).GT..9.OR.PLIDER(IRL(J)).LT.0.1) GO TO 533
      PLIDER(IRL(J))=0.
      IF(PPTCL(4,IRL(J)).LT.PMAX) GO TO 533
      PLIDER(IRL(IMAX))=0.
      IMAX=J
      PMAX=PPTCL(4,IRL(J))
      PLIDER(IRL(J))=PLI
533   CONTINUE
1537  CONTINUE
      IF(NLL.LE.1) GO TO 99
      DO 637 I=1,NLL
      IF(PLIDER(ILL(I)).LT.0.1.OR.PLIDER(ILL(I)).GT.0.9) GO TO 637
      PMAX=PPTCL(4,ILL(I))
      IMAX=I
      PLI=PLIDER(ILL(I))
      GO TO 638
637   CONTINUE
      GO TO 99
638   CONTINUE
      DO 633 J=1,NLL
      IF(PLIDER(ILL(J)).GT..9.OR.PLIDER(ILL(J)).LT.0.1) GO TO 633
      PLIDER(ILL(J))=0.
      IF(PPTCL(4,ILL(J)).LT.PMAX) GO TO 633
      PLIDER(ILL(IMAX))=0.
      IMAX=J
      PMAX=PPTCL(4,ILL(J))
      PLIDER(ILL(J))=PLI
633   CONTINUE
C  LORENTZ TRANSFORMATION FROM CMS TO LAB. FRAME
 99   IF(ICMS.EQ.1) GO TO 9700
      BACK=.TRUE.
C LORENTZ BOOST OF ENERGIES & MOMENTA
      CALL LORTR(V,NIN1,NFIN2,BACK)
C LORENTZ BOOST OF COORDINATES & TIMES
      CALL LORCO(V,NIN1,NFIN2,BACK)
9700   CONTINUE
      if(NCAS >= NCPRI) PRINT *, ' FROM COLLHH'
      RETURN
      END
c   ********************************************************
      SUBROUTINE DIFSMA(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE LOW MASS DIFFRACTION
C
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMLID/PLIDER(499)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON /ORDER/ IRD1,IRD2
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/PRIMP0/ P0
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PARCUT/ SWMAX
      COMMON/COMASS/AM1,AM2
      DIMENSION V(3),PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      LOGICAL BACK
      LOGICAL SPINT
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION GAM(3),AMR(3)
      DATA SIGMAD/0.23/
      DATA GAM/0.03,0.03,0.03/
      DATA AMR/1.40,1.10,1.30/
C     INITIALIZE
      NREP=0
      IRET=0
      MXPTCL=499
      XMINO=XMIN
      SIGMAO=SIGMA
      IPACK=1000
      DS=ECM-AM1-AM2
      PARBE=0.21
      IF(ECM.GT.AM1+AM2+PARBE) GO TO 150
      IRET=1
      go  to  2001
150   BACK=.TRUE.
      SIGMA=SIGMAO
      XMIN=XMINO
      SPINT=.TRUE.
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 101
999   IRET=1
C     WRITE(ITLIS,1200) NREP
1200  FORMAT(3X,' IN DIFSMA NREP GT ',I8)
      go  to  2001
101   CONTINUE
      SIGMA=0.30
      XMIN=0.
      IRET=0
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      IRDA=IRD1
      IRDB=IRD2
      PSIGN=-1.
      IF(RNDM(-1.).GT.0.5) GO TO 100
      PSIGN=1.
      IKA=IDIN(2)
      IKB=IDIN(1)
      AMA=AM2
      AMB=AM1
      IRDA=IRD2
      IRDB=IRD1
100   CALL FLAVOB(IKA,IFL1,IFL2)
C    COMPUTE X VALUES FOR PARTONS
      IB1=IB(IKA)
      IS1=IS(IKA)
      IFL1S=IFL1
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      X1=XDIST(XMIN,IB1,IS1,IFL1S)
      X2=1.-X1
C     COMPUTE PT VALUE FOR HADRON
160   CALL GAUSPT(PT2,SIGMA)
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 102
C   *********** 18.09.91 *************
      XMIN=XMINO
      SIGMA=SIGMAO
c   **********************************
      IRET=1
      go  to  2001
102   CONTINUE
      PHI=twpi*RNDM(-1.)
      PX2=PT2*COS(PHI)
      PY2=PT2*SIN(PHI)
      CALL GAUSPT(PT,SIGMAD)
      PHI=twpi*RNDM(-1.)
      PTX=PT*COS(PHI)
      PTY=PT*SIN(PHI)
      PX1=-PX2+PTX
      PY1=-PY2+PTY
      PT1=PX1**2+PY1**2
      AMH=AMB
      E1=SQRT((P0*X1)**2+PX1**2+PY1**2)
      E2=SQRT((P0*X2)**2+PX2**2+PY2**2)
      AMDTR=(E1+E2)**2-P0**2-PT**2
C     IF(AMDTR.GE.3.3) GO TO 160
      IDH=IDPARS(IFL1,IFL2,SPINT,2)
      IDHR=IDPARS(IFL1,IFL2,SPINT,1)
      AMHR=AMASS(IDHR)
      PARBE=0.2
      IF(IABS(IFL1).EQ.3.OR.IABS(IFL2).EQ.3) PARBE=0.3
      AMHB=AMASS(IDH)+PARBE
      IF(AMDTR.GE.AMHB**2) GO TO 200
      GO TO 150
200   AMD=SQRT(AMDTR)
      IF(ECM.LE.AMD+AMH) GO TO 160
      IBR=IB(IDHR)
      ISR=IS(IDHR)
      GAMRES=GAM(3)
      AMRES=AMR(3)
      IF(IBR.NE.0) AMRES=AMR(1)
      IF(IBR.EQ.0.AND.ISR.EQ.0) AMRES=AMR(2)
      IF(AMD.GE.AMRES) GO TO 162
      ARGWG=-(AMD-AMRES)**2/GAMRES
      IF(ARGWG.LE.-30.) ARGWG=-30.
      WG=EXP(ARGWG)
      IF(DS.LE.1.0.AND.IB1.EQ.0) GO TO 162
      DRND=RNDM(-1.)
      IF(DRND.GT.WG) GO TO 160
162   ALA=ALAMB(SCM,AMDTR,AMH**2)
      P0H=SQRT(ALA)/(2.*ECM)
      DTRM=P0H**2-PT**2
      IF(DTRM.LT.0.) GO TO 160
      PZH=SIGN(SQRT(DTRM),-PSIGN)
      ED=SQRT(AMD**2+P0H**2)
      V(1)=PTX/ED
      V(2)=PTY/ED
      V(3)=PZH/ED
      IF(AMD.GT.AMHR+SWMAX) GO TO 300
      NFIX=NPTCL
      NIN1=NPTCL+1
      CALL CLUSTR(IFL1,IFL2,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      CALL TIFILL(NIN1,NFIN1,AMD,IFL1,IFL2)
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=-PSIGN*P0*X1
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E1,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
610   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      GO TO 400
300   NFIX=NPTCL
      NIN1=NPTCL+1
      CALL STRING(IFL1,IFL2,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      NPRODS=NPTCL-NFIX
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=-PSIGN*P0*X1
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E1,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PPTCL(1,J)
      PPX1(2)=PPTCL(2,J)
      PPX1(3)=PPTCL(3,J)
      CALL ROTR(CT,ST,CFI,SFI,PPX1,PPX2,BACK)
      PPTCL(1,J)=PPX2(1)
      PPTCL(2,J)=PPX2(2)
      PPTCL(3,J)=PPX2(3)
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
510   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
400   CONTINUE
      DO 500 J=NIN1,NFIN1
      IORDP(J)=IRDA
500   IORIG(J)=6
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      PPTCL(5,NPTCL)=AMH
      PPTCL(1,NPTCL)=-PTX
      PPTCL(2,NPTCL)=-PTY
      PPTCL(3,NPTCL)=SIGN(PZH,-PZH)
      PPTCL(4,NPTCL)=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+
     *PPTCL(2,NPTCL)**2+PPTCL(3,NPTCL)**2)
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IKB
      IORIG(NPTCL)=6
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=IRDB
2001  XMIN=XMINO
      SIGMA=SIGMAO
      RETURN
9999  WRITE(ITLIS,1000) SCM,NPTCL
1000  FORMAT(//10X,38H....STOP IN DIFSMA..ENRGY TOO LOW SCM=,E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DIFTRI(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE TRIPLE POMERON VERTEX DIFFRACTION
C
      COMMON/COMLID/PLIDER(499)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON /ORDER/ IRD1,IRD2
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PRIMP0/ P0
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON /NEEDR/ NRET
      COMMON/COMASS/ AM1,AM2
      COMMON /COMQSE/ QSEE,QVSEE
      LOGICAL  QSEE,QVSEE
      LOGICAL RETU
      LOGICAL BACK
      DIMENSION V(3),PSUM(5)
      DATA SIGMAD/0.45/
C     INITIALIZE
      NREP=0
      ZER=0.0
      PARBE=0.30
      IRET=0
      IF(ECM.GT.(AM1+AM2+PARBE)) GO TO 999
      IRET=1
      RETURN
999   MXPTCL=499
      IPACK=1000
      NPTCLI=NPTCL
      DO 96 I=1,3
96    PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      NRET=0
      XMINO=XMIN
      XMIN=0.
      P0OLD=P0
      AMQ21=0.0
      AMQ22=0.0
100   IRET=0
      P0=P0OLD
      RETU=.FALSE.
      BACK=.TRUE.
      PSIGN=-1.
      PSOR=-1.
      EP=SQRT(P0**2+AM1**2)
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      IRDB=IRD2
      IF(RNDM(-1.).GT.0.5) GO TO 151
      IKA=IDIN(2)
      IKB=IDIN(1)
      AMA=AM2
      AMB=AM1
      IRDB=IRD1
      PSIGN=1.
      PSOR=1.
      EP=SQRT(P0**2+AM2**2)
151   B0=1/SIGMAD**2
      B=B0+0.6*LOG(SCM/2.0)
C     IF(IB(IDIN(1)).EQ.-1.AND.IB(IDIN(2)).NE.-1) B=12.0
      SIGMAN=SQRT(1./B)
C   DON'T CHANGE !!!
      PARBE=0.15
150   XMINS=(PARBE+AMA)**2/SCM
      NREP=NREP+1
C     WRITE(17,1993) NREP,NTRIES
1993  FORMAT(1X,'150C: ',2(1X,I10))
      IF(NREP.LT.NTRIES)  GO  TO 1994
      P0=P0OLD
      IRET=1
      RETURN
1994  CONTINUE
C   COMPUTE X VALUE FOR SEE QUARKS
C      PRINT  *, 'TO XSDIS'
C     WRITE(17,1995) XMINS,XMAX
1995  FORMAT(1X,'XMINS,XMAX=',2(F9.5,1X))
      CALL XSDIS(XS,XMINS,XMAX)
C      PRINT  *, 'FROM XSDIS'
C    COMPUTE PT VALUE FOR HADRON
C      PRINT *, 'TO GAUSPT'
      CALL GAUSPT(PTH,SIGMAN)
C      PRINT  *, 'FROM GAUSPT'
      PHI=twpi*RNDM(-1.)
      PTHX=PTH*COS(PHI)
      PTHY=PTH*SIN(PHI)
      PS=XS*P0
      ES=PS
      AMH=AMB
      AMD2=(EP+ES)**2-(P0-PS)**2
      AMD=SQRT(AMD2)
      IF(ECM.LE.AMD+AMH) GO TO 150
      ALA=ALAMB(SCM,AMD2,AMH**2)
      P0H=SQRT(ALA)/(2.0*ECM)
      DTRM=P0H**2-PTH**2
      IF(DTRM.LT.0.) GO TO 150
      PZH=SQRT(DTRM)*PSIGN
      EH=SQRT(AMH**2+P0H**2)
      ED=SQRT(AMD2+P0H**2)
      V(1)=-PTHX/ED
      V(2)=-PTHY/ED
      V(3)=-PZH/ED
      NREP=0
170   IFLS1=1+INT(RNDM(-1.)/PUD)
      NREP=NREP+1
C     WRITE(17,1991) NREP,NTRIES
1991  FORMAT(1X,'170C: ',2(1X,I10))
      IF(NREP.LT.NTRIES)  GO  TO 1992
      P0=P0OLD
      IRET=1
      RETURN
1992  CONTINUE
      IFLS2=IFLS1
C    COMPUTE X VALUES FOR PARTONS
      CALL FLAVOB(IKA,IFL1,IFL2)
      IB1=IB(IKA)
      IS1=IS(IKA)
      IFL1S=IFL1
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      XQVAL=XDIST(XMIN,IB1,IS1,IFL1S)
      XQQVA=1.-XQVAL
C   COMPUTE X VALUE FOR PION
      CALL XDIST2(XPI,X2PI)
      NIN1=NPTCL+1
C   IS THERE ANTIDIQUARK
      IF(MOD(IFL2,100).EQ.0.AND.IFL2.GT.0) GO TO 160
      IF(MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 160
      IFL1T=IFL1
      IFLS1T=IFLS1
      XQVALT=XQVAL
      XPIT=XPI
      IFL2T=IFL2
      IFLS2T=-IFLS2
      XQQVAT=XQQVA
      XPIT1=1.-XPI
      PSI=1.
      IF(PSIGN.LT.0.) GO TO 400
      IFL1T=IFLS1
      IFLS1T=IFL1
      XQVALT=XPI
      XPIT=XQVAL
      IFL2T=-IFLS2
      IFLS2T=IFL2
      XQQVAT=1.-XPI
      XPIT1=XQQVA
      PSI=-1.
400   PSIGN=PSIGN*PSI
      P0=AMD/2.0
      CALL XCORR(IFL1T,IFLS1T,ZER,ZER,ZER,ZER,XQVALT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 170
      CALL XCORR(IFL2T,IFLS2T,ZER,ZER,ZER,ZER,XQQVAT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      GO TO 170
130   NFIN1=NPTCL
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      GO TO 300
160   IFL1T=IFL1
      IFLS2T=-IFLS2
      XQVALT=XQVAL
      XPIT1=1.-XPI
      IFL2T=IFL2
      IFLS1T=IFLS1
      XQQVAT=XQQVA
      XPIT=XPI
      PSI=1.
      IF(PSIGN.LT.0.) GO TO 450
      IFL1T=-IFLS2
      IFLS2T=IFL1
      XQVALT=1.-XPI
      XPIT1=XQVAL
      IFL2T=IFLS2
      IFLS1T=IFL2
      XQQVAT=XPI
      XPIT=XQQVA
      PSI=-1.
450   PSIGN=PSIGN*PSI
      P0=AMD/2.0
      CALL XCORR(IFL1T,IFLS2T,ZER,ZER,ZER,ZER,XQVALT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 170
      CALL XCORR(IFL2T,IFLS1T,ZER,ZER,ZER,ZER,XQQVAT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 140
      NPTCL=NPTCL-NPRD
      GO TO 170
140   NFIN1=NPTCL
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
300   P0=P0OLD
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
500   IORIG(I)=1
      NPTCL=NPTCL+1
      NFIN1=NFIN1+1
      NPRD=NPTCL-NPTCLI
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      PPTCL(1,NPTCL)=PTHX
      PPTCL(2,NPTCL)=PTHY
      PPTCL(3,NPTCL)=PZH
      PPTCL(4,NPTCL)=EH
      PPTCL(5,NPTCL)=AMH
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IKB
      IORIG(NPTCL)=1
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=IRDB
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 501
      NPTCL=NPTCL-NPRD
      GO TO 100
501   XMIN=XMINO
      RETURN
9999  WRITE(ITLIS,1000) SCM,NPTCL
1000  FORMAT(//10X,38H...STOP IN DIFTRI..ENERGY TOO LOW SCM=,E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ANNIH(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C      COMPUTE ANNIHILATION DIAGRAM
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON/COMANN/ DIQAN
      LOGICAL DIQAN
      COMMON/CANPOM/ POAGEN(15)
      COMMON/COMCOL/ NAC(100,4),NBC(100,4),NCOL
      COMMON/COMPLI/ LIMP
      LOGICAL MULTP
      COMMON/COMMUL/MULTP
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMENB/ ENBOU
      COMMON/CSIGA/ SIGAN
      COMMON/SIGDIA/ CROSD(5),DST
      DIMENSION IF1(3),IF2(3),IFD1(3),IFD2(3)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      DATA COEFPL/ 300./,COEFTW/100./
      DATA ICALL/0/
C   INITIALIZE
      IRET=0
      ICALL=ICALL+1
      DIQAN=.TRUE.
      CALL FLAVOR(IDIN(1),IF1(1),IF1(2),IF1(3),JSPIN1,INDEX1)
      CALL FLAVOR(IDIN(2),IF2(1),IF2(2),IF2(3),JSPIN2,INDEX2)
C    CAN BE ONE SHEET OR TWO SHEETS ANNIHILATION
      DO 3 I=1,3
      IF1(I)=IABS(IF1(I))
      DO 3 J=1,3
      IF2(J)=IABS(IF2(J))
      IF(IF1(I).NE.IF2(J)) GO TO 3
      GO TO 4
3     CONTINUE
      SIGPL=0.
      SIGTWO=0.
      GO TO 8
C   CAN BE ONE SHEET ANNIHILATION
4     IFD1(1)=IABS(IF1(1))*1000+IABS(IF1(2))*100
      IFD1(2)=IABS(IF1(2))*1000+IABS(IF1(3))*100
      IFD1(3)=IABS(IF1(1))*1000+IABS(IF1(3))*100
      IFD2(1)=IABS(IF2(1))*1000+IABS(IF2(2))*100
      IFD2(2)=IABS(IF2(2))*1000+IABS(IF2(3))*100
      IFD2(3)=IABS(IF2(1))*1000+IABS(IF2(3))*100
      DO 5 I=1,3
      DO 5 J=1,3
      IF(IFD1(I).NE.IFD2(J)) GO TO 5
      GO TO 6
5     CONTINUE
      SIGPL=0.
      GO TO 7
 6     SIGPL=COEFPL*SCM**(-1.5)
      IF(ECM.LT.ENBOU) SIGPL=CROSD(4)
 7     SIGTWO=COEFTW/SCM
      IF(ECM.LT.ENBOU) SIGTWO=CROSD(5)
 8     SIGTH=SIGAN-SIGPL-SIGTWO
      IF(ICALL.EQ.1.AND.IPRINT)WRITE(ITLIS,9540) SIGPL,SIGTWO,SIGTH
 9540 FORMAT(//10X,' ANNIHILATION CROSS SECTION CONSISTS OF'/
     *   10X,' ONE SHEET DIAGRAMM WITH CR.SEC.=',E13.3,' MB'/
     *   10X,' TWO SHEET DIAGRAMM WITH CR.SEC.=',E13.3,' MB'/
     *   9X,'THREE SHEET DIAGRAMM WITH CR.SEC.=',E13.3,' MB'//)
      TRY=RNDM(-1.)
      IF(TRY.GT.SIGPL/SIGAN) GO TO 500
C   ONE SHEET ANNIHILATION
      CALL PLANAR(IRET)
C     SIGMA=SIGMA1
      RETURN
500   IF(TRY.GT.(SIGPL+SIGTWO)/SIGAN) GO TO 600
C    TWO SHEETS ANNIHILATION
      CALL TWOSHE(IRET)
C     SIGMA=SIGMA1
      RETURN
C   THREE SHEETS ANNIHILATION
600   IF(.NOT.MULTP) GO TO 700
C  SELECT NUMBER OF POMERONS
      TRY=RNDM(-1.)
      DO 710 NPOM=1,LIMP
      NC=NPOM
      IF(POAGEN(NPOM).GT.TRY) GO TO 800
710   CONTINUE
800   CONTINUE
      NCOL=NC
      DO 900 J=1,NCOL
      NAC(J,4)=0
      NBC(J,4)=0
      NAC(J,3)=0
      NBC(J,3)=0
      NAC(J,1)=1
      NBC(J,1)=1
      NAC(J,2)=IDIN(1)
      NBC(J,2)=IDIN(2)
900   CONTINUE
      CALL CHAINS(IRET)
C     SIGMA=SIGMA1
      RETURN
700   CALL THREES(IRET)
C     SIGMA=SIGMA1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DOUBDI(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C-----------------------------------------------------------------------
C  COMPUTE DOUBLE DIFFRACTION DISSOCIATION
C-----------------------------------------------------------------------
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMECB/ ECMB
      COMMON/NCASCA/NCAS,NCPRI
      IF(ECM.GT.ECMB) GO TO 95
      if(NCAS >= NCPRI) PRINT *, ' TO DOUBSM'
      CALL DOUBSM(IRET)
      if(NCAS >= NCPRI) PRINT *, ' from DOUBSM'
      RETURN
95    IF(RNDM(-1.).GT.0.5) GO TO 96
      if(NCAS >= NCPRI) PRINT *, ' TO DOUBLO ECM,ECMB=',ECM,ECMB
      CALL DOUBLO(IRET)
      if(NCAS >= NCPRI) PRINT *, ' from DOUBLO'
      RETURN
96    continue
      if(NCAS >= NCPRI) PRINT *, ' TO DOUBY'
      CALL DOUBY(IRET)
      if(NCAS >= NCPRI) PRINT *, ' from DOUBY'
      RETURN
      END
C***********************************************************************
      SUBROUTINE DOUBLO(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE ENHANCEMENT (LOOP)-POMERONS DIFFRACTION
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMASS/ AM1,AM2
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMECB/ ECMB
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PRIMP0/ P0
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/REACOE/ COEF(11),COEF1(11)
      COMMON /COMQSE/ QSEE,QVSEE
      COMMON /NEEDR/ NRET
      COMMON/NCASCA/NCAS,NCPRI
      LOGICAL QSEE,QVSEE
      LOGICAL RETU
      LOGICAL BACK
      DIMENSION V1(3),V2(3),PSUM(5)
      DATA SIGMAD/0.45/
C     INITIALIZE
      if(NCAS >= NCPRI) PRINT *, ' in DOUBLO NRTIES=',NTRIES
      NREP=0
      ZE=0.
      PARBE=0.21
      QVSEE=.TRUE.
      IRET=0
      IF(ECM.GT.AM1+AM2+3.*PARBE) GO TO 999
      IRET=1
      GO  TO  9999
999   MXPTCL=499
      IPACK=1000
      NPTCLI=NPTCL
      DO 96 I=1,3
96    PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      NRET=0
      XMINO=XMIN
      XMIN=0.
      P0OLD=P0
      AMQ21=0.0
      AMQ22=0.0
100   IRET=0
      P0=P0OLD
      RETU=.FALSE.
      BACK=.TRUE.
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      B0=1/SIGMAD**2
      B=B0+0.6*LOG(SCM/2.0)
C     IF(IB(IKA).EQ.-1.AND.IB(IKB).NE.-1) B=12.0
      SIGMAN=SQRT(1./B)
      PARBE=0.3
150   XMINS1=(PARBE+AMA)**2/SCM
      NREP=NREP+1
      IF(NREP.LT.NTRIES)  GO  TO  1994
      P0=P0OLD
      XMIN=XMINO                         ! 17.05.2002
      IRET=1
      RETURN
1994  CONTINUE
C   COMPUTE X VALUE FOR SEE QUARKS
      CALL XSDIS(XS1,XMINS1,XMAX)
      XS3=1.-XS1
C    COMPUTE PT VALUE
      CALL GAUSPT(PTH,SIGMAN)
      PHI=twpi*RNDM(-1.)
      PTHX=PTH*COS(PHI)
      PTHY=PTH*SIN(PHI)
      XMINS2=(PARBE+AMB)**2/SCM
C   COMPUTE X VALUE FOR SEE QUARKS
      CALL XSDIS(XS2,XMINS2,XMAX)
      XS4=1.-XS2
      P1=P0*XS1
      P2=PSIGN*P0*XS2
      P3=P0*XS3
      P4=PSIGN*P0*XS4
      E1=P1
      E2=ABS(P2)
      E3=P3
      E4=ABS(P4)
      E14=E1+E4
      P14=P1+P4
      AMD14=E14**2-P14**2
      E23=E2+E3
      P23=P2+P3
      AMD23=E23**2-P23**2
      AMD1=SQRT(AMD14)
      AMD2=SQRT(AMD23)
      IF(ECM.LE.AMD1+AMD2) GO TO 150
      ALA=ALAMB(SCM,AMD14,AMD23)
      P0H=SQRT(ALA)/(2.0*ECM)
      DTRM=P0H**2-PTH**2
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO1: NREP,PTH,P0H',NREP,PTH,P0H
      IF(DTRM.LT.0.) GO TO 150
      PZH14=SQRT(DTRM)*PSIGN
      PZH23=-PZH14
      V1(1)=PTHX/E14
      V1(2)=PTHY/E14
      V1(3)=PZH14/E14
      V2(1)=-PTHX/E23
      V2(2)=-PTHY/E23
      V2(3)=PZH23/E23
      NIN1=NPTCL+1
      NREP=0
170   IFLS1=1+INT(RNDM(-1.)/PUD)
      NREP=NREP+1
      IF(NREP.LT.NTRIES)  GO  TO  1992
2002  P0=P0OLD                                ! 17.05.2002
      XMIN=XMINO                              ! 17.05.2002
      IRET=1
      RETURN
1992  CONTINUE
      CALL FLAVOB(IKA,IFL1,IFL2)
      IB1=IB(IKA)
      IS1=IS(IKA)
      IFL1S=IFL1
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      XQVAL=XDIST(XMIN,IB1,IS1,IFL1S)
      XQQVA=1.-XQVAL
C   COMPUTE X VALUE FOR PION
      CALL XDIST2(XPI,X2PI)
C   IS THERE ANTIDIQUARK
      IF(MOD(IFL2,100).EQ.0.AND.IFL2.GT.0) GO TO 160
      IF(MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 160
      IFL1T=IFL1
      IFLS1T=IFLS1
      XQVALT=XQVAL
      XPIT=XPI
      IFL2T=IFL2
      IFLS2T=-IFLS1
      XQQVAT=XQQVA
      XPIT1=1.-XPI
      P0=AMD1/2.0
      CALL XCORR(IFL1T,IFLS1T,ZE,ZE,ZE,ZE,XQVALT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO2: NREP,IFL1T,IFLS1T',
     &                         	NREP,IFL1T,IFLS1T
      IF(RETU) GO TO 170
      CALL XCORR(IFL2T,IFLS2T,ZE,ZE,ZE,ZE,XQQVAT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO3: NREP,IFL2T,IFLS2T',
     &                         	NREP,IFL2T,IFLS2T
      GO TO 170
130   NFIN1=NPTCL
      CALL LORTR(V1,NIN1,NFIN1,BACK)
      CALL LORCO(V1,NIN1,NFIN1,BACK)
      GO TO 180
160   IFL1T=IFL1
      IFLS2T=-IFLS1
      XQVALT=XQVAL
      XPIT1=1.-XPI
      IFL2T=IFL2
      IFLS1T=IFLS1
      XQQVAT=XQQVA
      XPIT=XPI
      P0=AMD1/2.0
      CALL XCORR(IFL1T,IFLS2T,ZE,ZE,ZE,ZE,XQVALT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO4: NREP,IFL1T,IFLS2T',
     &                         	NREP,IFL1T,IFLS2T
      IF(RETU) GO TO 170
      CALL XCORR(IFL2T,IFLS1T,ZE,ZE,ZE,ZE,XQQVAT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 140
      NPTCL=NPTCL-NPRD
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO5: NREP,IFL2T,IFLS1T',
     &                         	NREP,IFL2T,IFLS1T
      GO TO 170
140   NFIN1=NPTCL
      CALL LORTR(V1,NIN1,NFIN1,BACK)
      CALL LORCO(V1,NIN1,NFIN1,BACK)
180   NPRDAL=NPRD
      NIN1=NPTCL+1
181   IFLS1=1+INT(RNDM(-1.)/PUD)
      NREP=NREP+1
	if(NREP > NTRIES)  go  to  2002       ! 17.05.2002
      CALL FLAVOB(IKB,IFL1,IFL2)
      IB1=IB(IKB)
      IS1=IS(IKB)
      IFL1S=IFL1
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      XQVAL=XDIST(XMIN,IB1,IS1,IFL1S)
      XQQVA=1.-XQVAL
C   COMPUTE X VALUE FOR PION
      CALL XDIST2(XPI,X2PI)
C   IS THERE ANTIDIQUARK
      IF(MOD(IFL2,100).EQ.0.AND.IFL2.GT.0) GO TO 190
      IF(MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 190
      IFL1T=IFLS1
      IFLS1T=IFL1
      XQVALT=XPI
      XPIT=XQVAL
      IFL2T=-IFLS1
      IFLS2T=IFL2
      XQQVAT=1.-XPI
      XPIT1=XQQVA
      P0=AMD2/2.0
      CALL XCORR(IFL1T,IFLS1T,ZE,ZE,ZE,ZE,XQVALT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO6: NREP,IFL1T,IFLS1T',
     &                         	NREP,IFL1T,IFLS1T
      IF(RETU) GO TO 181
      CALL XCORR(IFL2T,IFLS2T,ZE,ZE,ZE,ZE,XQQVAT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 230
      NPTCL=NPTCL-NPRD
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO7: NREP,IFL2T,IFLS2T',
     &                         	NREP,IFL2T,IFLS2T
      GO TO 181
230   NFIN1=NPTCL
      CALL LORTR(V2,NIN1,NFIN1,BACK)
      CALL LORCO(V2,NIN1,NFIN1,BACK)
      GO TO 250
190   IFL1T=-IFLS1
      IFLS2T=IFL1
      XQVALT=1.-XPI
      XPIT1=XQVAL
      IFL2T=IFLS1
      IFLS1T=IFL2
      XQQVAT=XPI
      XPIT=XQQVA
      P0=AMD2/2.0
      CALL XCORR(IFL1T,IFLS2T,ZE,ZE,ZE,ZE,XQVALT,XPIT1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO8: NREP,IFL1T,IFLS2T',
     &                         	NREP,IFL1T,IFLS2T
      IF(RETU) GO TO 181
      CALL XCORR(IFL2T,IFLS1T,ZE,ZE,ZE,ZE,XQQVAT,XPIT,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 240
      NPTCL=NPTCL-NPRD
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO9: NREP,IFL2T,IFLS1T',
     &                         	NREP,IFL2T,IFLS1T
      GO TO 181
240   NFIN1=NPTCL
      CALL LORTR(V2,NIN1,NFIN1,BACK)
      CALL LORCO(V2,NIN1,NFIN1,BACK)
250   P0=P0OLD
      NIN1=NPTCLI+1
      NPRDAL=NPRDAL+NPRD
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
500   IORIG(I)=11
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 501
      NPTCL=NPTCL-NPRDAL
      if(NCAS >= NCPRI) PRINT *, ' DOUBLO10: NREP,NIN1,NFIN1=',
     &                         	NREP,NIN1,NFIN1
      GO TO 100
501   XMIN=XMINO
      RETURN
9999  CONTINUE
C     WRITE(ITLIS,1000) SCM,NPTCL
1000  FORMAT(//10X,38H...STOP IN DOUBLO..ENERGY TOO LOW SCM=,E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DOUBY(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE DOUBLE ENHANCEMENT (Y)-POMERONS DIFFRACTION
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMLID/PLIDER(499)
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMASS/ AM1,AM2
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PRIMP0/ P0
      COMMON /ORDER/ IRD1,IRD2
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/REACOE/ COEF(11),COEF1(11)
      COMMON /COMQSE/ QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      COMMON /NEEDR/ NRET
      LOGICAL RETU
      LOGICAL BACK
      DIMENSION V1(3),V2(3),PSUM(5)
      DATA SIGMAD/0.45/
C     INITIALIZE
      NREP=0
      PARBE=0.21
      QSEE=.TRUE.
      QVSEE=.FALSE.
      IRET=0
      IF(ECM.GT.AM1+AM2+3.*PARBE) GO TO 999
      IRET=1
      RETURN
999   MXPTCL=499
      IPACK=1000
      NPTCLI=NPTCL
      DO 96 I=1,3
96    PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      NRET=0
      XMINO=XMIN
      XMIN=0.
      AMQ21=0.0
      AMQ22=0.0
      IRET=0
      RETU=.FALSE.
      BACK=.TRUE.
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      IRDA=IRD1
      IRDB=IRD2
      B0=1/SIGMAD**2
      B=B0+0.6*LOG(SCM/2.0)
C     IF(IB(IKA).EQ.-1.AND.IB(IKB).NE.-1) B=12.0
      SIGMAN=SQRT(1./B)
      PARBE=0.3
150   XMINS1=(PARBE+AMA)**2/SCM
      NREP=NREP+1
      IF(NREP.LT.NTRIES)  GO  TO  1992
      IRET=1
      RETURN
1992  CONTINUE
      NPRDAL=0
C   COMPUTE X VALUE FOR SEE QUARKS
      CALL XSDIS(XS1,XMINS1,XMAX)
      XS3=1.-XS1
      XS11=XMINS1+(XS1-XMINS1)*RNDM(-1.)
      XS12=XS1-XS11
C    COMPUTE PT VALUE
      CALL GAUSPT(PTH,SIGMAN)
      PHI=twpi*RNDM(-1.)
      PTHX1=PTH*COS(PHI)
      PTHY1=PTH*SIN(PHI)
      PTHX11=PTHX1*RNDM(-1.)
      PTHX12=PTHX1-PTHX11
      PTHY11=PTHY1*RNDM(-1.)
      PTHY12=PTHY1-PTHY11
      PTHX3=-PTHX1
      PTHY3=-PTHY1
      XMINS2=(PARBE+AMB)**2/SCM
C   COMPUTE X VALUE FOR SEE QUARKS
      CALL XSDIS(XS2,XMINS2,XMAX)
      XS4=1.-XS2
      XS21=XMINS2+(XS2-XMINS2)*RNDM(-1.)
      XS22=XS2-XS21
C    COMPUTE PT VALUE
      CALL GAUSPT(PTH,SIGMAN)
      PHI=twpi*RNDM(-1.)
      PTHX2=PTH*COS(PHI)
      PTHY2=PTH*SIN(PHI)
      PTHX21=PTHX2*RNDM(-1.)
      PTHX22=PTHX2-PTHX21
      PTHY21=PTHY2*RNDM(-1.)
      PTHY22=PTHY2-PTHY21
      PTHX4=-PTHX2
      PTHY4=-PTHY2
      P11=P0*XS11
      P12=P0*XS12
      P21=PSIGN*P0*XS21
      P22=PSIGN*P0*XS22
      P3=P0*XS3
      P4=PSIGN*P0*XS4
      E11=SQRT(P11**2+PTHX11**2+PTHY11**2)
      E12=SQRT(P12**2+PTHX12**2+PTHY12**2)
      E21=SQRT(P21**2+PTHX21**2+PTHY21**2)
      E22=SQRT(P22**2+PTHX22**2+PTHY22**2)
      E3=SQRT(P3**2+PTHX3**2+PTHY3**2+AMA**2)
      E4=SQRT(P4**2+PTHX4**2+PTHY4**2+AMB**2)
C  Q11+Q22
      E1=E11+E22
      P1=P11+P22
      PX1=PTHX11+PTHX22
      PY1=PTHY11+PTHY22
      AMD1=E1**2-P1**2-PX1**2-PY1**2
C  Q12+Q21
      E2=E12+E21
      P2=P12+P21
      PX2=PTHX12+PTHX21
      PY2=PTHY12+PTHY21
      AMD2=E2**2-P2**2-PX2**2-PY2**2
C
      AMD1=SQRT(AMD1)
      AMD2=SQRT(AMD2)
      IF(ECM.LE.AMD1+AMD2+2.0) GO TO 150
      V1(1)=PX1/E1
      V1(2)=PY1/E1
      V1(3)=P1/E1
      V2(1)=PX2/E2
      V2(2)=PY2/E2
      V2(3)=P2/E2
      NIN1=NPTCL+1
      IFL11=1+INT(RNDM(-1.)/PUD)
      IFL12=-IFL11
      IFL21=1+INT(RNDM(-1.)/PUD)
      IFL22=-IFL21
C
      CALL XCORR(IFL11,IFL22,PTHX11,PTHY11,PTHX22,PTHY22,XS11,XS22,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 150
      NFIN1=NPTCL
      CALL LORTR(V1,NIN1,NFIN1,BACK)
      CALL LORCO(V1,NIN1,NFIN1,BACK)
      NIN1=NPTCL+1
      CALL XCORR(IFL12,IFL21,PTHX12,PTHY12,PTHX21,PTHY21,XS12,XS21,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      GO TO 150
130   NFIN1=NPTCL
      CALL LORTR(V2,NIN1,NFIN1,BACK)
      CALL LORCO(V2,NIN1,NFIN1,BACK)
C
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      PPTCL(1,NPTCL)=PTHX3
      PPTCL(2,NPTCL)=PTHY3
      PPTCL(3,NPTCL)=P3
      PPTCL(4,NPTCL)=E3
      PPTCL(5,NPTCL)=AMA
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IKA
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=IRDA
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      PPTCL(1,NPTCL)=PTHX4
      PPTCL(2,NPTCL)=PTHY4
      PPTCL(3,NPTCL)=P4
      PPTCL(4,NPTCL)=E4
      PPTCL(5,NPTCL)=AMB
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IKB
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=IRDB
      NIN1=NPTCLI+1
      NPRDAL=NPRDAL+NPRD+2
      NFIN1=NPTCL
      DO 500 I=NIN1,NFIN1
      IF(I.GE.NPTCL-1) GO TO 500
      IORDP(I)=0
500   IORIG(I)=11
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 501
      NPTCL=NPTCL-NPRDAL
      GO TO 150
501   XMIN=XMINO
      RETURN
9999  WRITE(ITLIS,1000) SCM,NPTCL
1000  FORMAT(//10X,38H...STOP IN DOUBY ..ENERGY TOO LOW SCM=,E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C***********************************************************************
      SUBROUTINE DOUBSM(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE DOUBLE SMALL MASS DIFFRACTION
C
      COMMON/PRIMAR/ SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PRIMP0/ P0
      COMMON/COMASS/ AM1,AM2
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON /NEEDR/ NRET
      LOGICAL RETU
      LOGICAL BACK
      LOGICAL SPINT
      DIMENSION PSUM(5),GAM(3),AMR(3)
      DATA SIGMAN/0.23/
      DATA GAM/0.03,0.03,0.03/
      DATA AMR/1.40,1.10,1.30/
C     INITIALIZE
      XMINO=XMIN
      SIGMA0=SIGMA
      NREP=0
      PARBE=0.23
      SPINT=.TRUE.
      IRET=0
      IF(ECM.GT.AM1+AM2+2.*PARBE) GO TO 999
      IRET=1
      GO  TO  9999
 999  DS=ECM-AM1-AM2
      MXPTCL=499
      IPACK=1000
      DO 96 I=1,3
96    PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
      NRET=1
      XMIN=0.
      SIGMA=0.3
      AMQ21=0.0
      AMQ22=0.0
100   IRET=0
      NPRDAL=0
      RETU=.FALSE.
      BACK=.TRUE.
      PSIGN=1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      PARBE=0.22
150   CALL FLAVOB(IKA,IFL1,IFL2)
C   COMPUTE X VALUE FOR VALENCE QUARK
      IB1=IB(IKA)
      IS1=IS(IKA)
      IFL1S=IFL1
      IF(IS1.NE.0.AND.IB1.EQ.0.AND.IABS(IFL1).EQ.3) IFL1S=IFL2
      X1=XDIST(XMIN,IB1,IS1,IFL1S)
      X2=1.-X1
      CALL FLAVOB(IKB,IFL3,IFL4)
C   COMPUTE X VALUE FOR VALENCE QUARK
      IB2=IB(IKB)
      IS2=IS(IKB)
      IFL3S=IFL3
      IF(IS2.NE.0.AND.IB2.EQ.0.AND.IABS(IFL3).EQ.3) IFL3S=IFL4
      X03=XDIST(XMIN,IB2,IS2,IFL3S)
      X04=1.-X03
      X3=-X03
      X4=-X04
C COMPUTE PT VALUE FOR HADRON
160   CALL GAUSPT(PTH,SIGMAN)
      NREP=NREP+1
      IF(NREP.LT.NTRIES)  GO  TO  1994
      IRET=1
      go  to  2001
1994  CONTINUE
      PHI=twpi*RNDM(-1.)
      PTXH=PTH*COS(PHI)
      PTYH=PTH*SIN(PHI)
C    COMPUTE PT VALUE
      CALL GAUSPT(PT1,SIGMA)
      PHI=twpi*RNDM(-1.)
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      CALL GAUSPT(PT2,SIGMA)
      PHI=twpi*RNDM(-1.)
      PX2=PT2*COS(PHI)
      PY2=PT2*SIN(PHI)
      PX1H=PTXH-PX1
      PY1H=PTYH-PY1
      PX2H=-PTXH-PX2
      PY2H=-PTYH-PY2
      P1=P0*X1
      P2=P0*X2
      P3=P0*X3
      P4=P0*X4
      E11=SQRT(P1**2+PX1**2+PY1**2)
      E12=SQRT(P2**2+PX1H**2+PY1H**2)
      E23=SQRT(P3**2+PX2**2+PY2**2)
      E24=SQRT(P4**2+PX2H**2+PY2H**2)
      E1=E11+E12
      E2=E23+E24
      AMD1=E1**2-P0**2-PTXH**2-PTYH**2
      IDH1=IDPARS(IFL1,IFL2,SPINT,2)
      AMHB1=AMASS(IDH1)+PARBE
      IF(AMD1.LT.AMHB1**2) GO TO 160
      IDHR1=IDPARS(IFL1,IFL2,SPINT,1)
      IBR1=IB(IDHR1)
      ISR1=IS(IDHR1)
      GAMRES=GAM(3)
      AMRES=AMR(3)
      IF(IBR1.NE.0) AMRES=AMR(1)
      IF(IBR1.EQ.0.AND.ISR1.EQ.0) AMRES=AMR(2)
      AMD10=SQRT(AMD1)
      IF(AMD10.GE.AMRES) GO TO 162
      ARGWG1=-(AMD10-AMRES)**2/GAMRES
      IF(ARGWG1.LE.-30.) ARGWG1=-30.
      WG1=EXP(ARGWG1)
      IF(DS.LE.1.0.AND.IB1.EQ.0) GO TO 162
      IF(RNDM(-1.).GT.WG1) GO TO 160
162   AMD2=E2**2-P0**2-PTXH**2-PTYH**2
      IDH2=IDPARS(IFL3,IFL4,SPINT,2)
      AMHB2=AMASS(IDH2)+PARBE
      IF(AMD2.LT.AMHB2**2) GO TO 160
      IDHR2=IDPARS(IFL3,IFL4,SPINT,1)
      IBR2=IB(IDHR2)
      ISR2=IS(IDHR2)
      GAMRES=GAM(3)
      AMRES=AMR(3)
      IF(IBR2.NE.0) AMRES=AMR(1)
      IF(IBR2.EQ.0.AND.ISR2.EQ.0) AMRES=AMR(2)
      AMD20=SQRT(AMD2)
      IF(AMD20.GE.AMRES) GO TO 163
      ARGWG2=-(AMD20-AMRES)**2/GAMRES
      IF(ARGWG2.LE.-30.) ARGWG2=-30.
      WG2=EXP(ARGWG2)
      IF(DS.LE.1.0.AND.IB2.EQ.0) GO TO 163
      IF(RNDM(-1.).GT.WG2) GO TO 160
163   CONTINUE
      IF(ECM.LE.AMD10+AMD20) GO TO 160
      ALA=ALAMB(SCM,AMD1,AMD2)
      P0H=SQRT(ALA)/(2.0*ECM)
      DTRM=P0H**2-PTH**2
      IF(DTRM.LT.0.) GO TO 160
      NIN1=NPTCL+1
      CALL XCORR(IFL1,IFL2,PX1,PY1,PX1H,PY1H,X1,X2,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 160
      CALL XCORR(IFL3,IFL4,PX2,PY2,PX2H,PY2H,X3,X4,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      GO TO 160
130   NFIN2=NPTCL
      NPRDAL=NPRDAL+NPRD
      IF(NPRDAL.NE.2) GO TO 131
      IF(IKA.NE.IDENT(NPTCL-1).OR.IKB.NE.IDENT(NPTCL)) GO TO 131
      NPTCL=NPTCL-NPRDAL
131   CALL RESCAL(NIN1,NFIN2,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 501
      NPTCL=NPTCL-NPRDAL
      GO TO 100
501   continue
      DO 500 I=NIN1,NFIN2
      IORDP(I)=0
500   IORIG(I)=11
2001  XMIN=XMINO
      SIGMA=SIGMA0
      RETURN
9999  CONTINUE
C     WRITE(ITLIS,1000) SCM,NPTCL
1000  FORMAT(//10X,38H...STOP IN DOUBSM..ENERGY TOO LOW SCM=,E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE CHAINS(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C   FORM COLOUR NEUTRAL CHAINS FROM QUARKS
C
C
      DIMENSION PSUM(5)
       COMMON/COMFLA/MNASEA(12),MNBSEA(12),IFLAS(12),IFLBS(12),
     * NUAVAL,
     * NUBVAL,IFLQA1,IFLQB1,IFLQA2,IFLQB2,IFAQQ,IFBQQ
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMLID/PLIDER(499)
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      COMMON/PRIMP0/ P0
      COMMON/COMCOL/ NAC(100,4),NBC(100,4),NCOL
       COMMON/COMVA/ XAVAL1,XAVAL2,XAQQ,XASEA1(12),
     * XASEA2(12),NPOMA
       COMMON/COMVB/ XBVAL1,XBVAL2,XBQQ,XBSEA1(12),
     * XBSEA2(12),NPOMB
       COMMON/COMPYA/ PYAV1,PYAV2,PYAQQ,
     *PYAS1(12),PYAS2(12)
       COMMON/COMPXA/ PXAV1,PXAV2,PXAQQ,
     *PXAS1(12),PXAS2(12)
       COMMON/COMPXB/ PXBV1,PXBV2,PXBQQ,
     *PXBS1(12),PXBS2(12)
       COMMON/COMPYB/ PYBV1,PYBV2,PYBQQ,
     *PYBS1(12),PYBS2(12)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/COMQMA/ AMQUA1,AMQUA2,AMQQA,
     *AMQAS1(12),AMQAS2(12)
      COMMON/COMQMB/ AMQUB1,AMQUB2,AMQQB,
     *AMQBS1(12),AMQBS2(12)
      COMMON/MASQUA/ AMQ21,AMQ22
      DIMENSION JSA0(12),JSB0(12)
      COMMON/COMXM/ XMIN,XMAX
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMDIF/ NDIFA,NDIFB
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/NEEDR/ NRET
      COMMON /UNCYS/ NUNCY
!      COMMON /UNCYS1/ NUNC
      LOGICAL RETU
      LOGICAL DIQAN
      COMMON/COMANN/ DIQAN
      COMMON/COMASS/ AM1,AM2
      COMMON/KAPPA/ XAP
      COMMON/NPTCLZ/NPTCLZ
C   COMPUTE QUARK PARAMETERS
C     INITIALIZE
      MXPTCL=499
      NUNC=0
      PARBE=0.3
      IRET=0
      IF(ECM.GT.AM1+AM2+3.*PARBE) GO TO 999
      IRET=1
      RETURN
999   IREP=0
      IPACK=1000
      NRET=0
      DO 151 I=1,3
151   PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
1111  NPTCL=0
      QSEE=.FALSE.
      NIN=NPTCL+1
      IREP=IREP+1
      IF(IREP.LE.NTRIES) GO TO 1112
      IRET=1
      RETURN
1112  PSIGN=-1.
       NVSA=0
       NVSB=0
       NPOMA=1
       NPOMB=1
      NAVAL=0
      NBVAL=0
       DO 114 J=1,12
       JSA0(J)=0
       JSB0(J)=0
  114  CONTINUE
       NASEA=0
       NBSEA=0
       DO 206 JC=1,NCOL
       IPAR=NAC(JC,1)
       JPAR=NBC(JC,1)
       IKA=NAC(JC,2)
       IBA=IB(IKA)
       IKB=NBC(JC,2)
       IBB=IB(IKB)
      IDIF=NAC(JC,3)
      JDIF=NBC(JC,3)
      IF(NAVAL.EQ.0) GO TO 102
       IF(NUAVAL.NE.IPAR) GO TO 101
       NPOMA=NPOMA+1
      NASEA=NASEA+1
      MNASEA(NASEA)=1
      GO TO 103
  101  CONTINUE
102   NAVAL=NAVAL+1
       CALL FLAVOB(IKA,IFL1,IFL2)
       IFLQA1=IFL1
       IFLQA2=0
       IFAQQ=IFL2
       IFAHQ=IKA
      IF(IBA.NE.0) GO TO 123
      IF(IABS(IFL1).NE.3) GO TO 123
         IFAQQ=IFL1
         IFLQA1=IFL2
123    NUAVAL=IPAR
      NDIFA=IDIF
      IF(.NOT.DIQAN) GO TO 103
      CALL FLAVOR(IFL2,IFLQA2,IFAQQ,IFK ,JSPIN,INDEX)
      IF(RNDM(-1.).GE.0.5) GO TO 103
      ISWAP=IFAQQ
      IFAQQ=IFLQA2
      IFLQA2=ISWAP
103   IF(NBVAL.EQ.0) GO TO 105
       IF(NUBVAL.NE.JPAR) GO TO 104
       NPOMB=NPOMB+1
       NBSEA=NBSEA+1
       MNBSEA(NBSEA)=1
       GO TO 206
  104  CONTINUE
105   NBVAL=NBVAL+1
       CALL FLAVOB(IKB,IFL1,IFL2)
       IFLQB1=IFL1
       IFLQB2=0
       IFBQQ=IFL2
       IFBHQ=IKB
       NUBVAL=JPAR
       NDIFB=JDIF
      IF(.NOT.DIQAN) GO TO 206
      CALL FLAVOR(IFL2,IFLQB2,IFBQQ,IFK ,JSPIN,INDEX)
      IF(RNDM(-1.).GE.0.5) GO TO 206
      ISWAP=IFBQQ
      IFBQQ=IFLQB2
      IFLQB2=ISWAP
  206  CONTINUE
C  COMPUTE X AND PT VALUES FOR PARTONS
       XMIE=0.015/ECM
       CALL XQUARK(0,XMIE,XMAX,ALFA,BETA,IBA)
       CALL XQUARK(1,XMIE,XMAX,ALFA,BETA,IBB)
       CALL PTQUAR(0)
       CALL PTQUAR(1)
       IF(NASEA.EQ.0) GO TO 506
       JSD=1
       DO 505 JSA=1,NASEA
       JSA0(JSA)=JSA0(JSA)+JSD
       JSD=JSD+1
  505  CONTINUE
  506  CONTINUE
      IF(NBSEA.EQ.0) GO TO 634
       JSD=1
       DO 632 JSB=1,NBSEA
       JSB0(JSB)=JSB0(JSB)+JSD
       JSD=JSD+1
  632  CONTINUE
 634  CONTINUE
C
C    SELECT SEA-SEA CHAINS
C
      NAVAL=0
      NBVAL=0
       NASEA=0
       NBSEA=0
       DO 6 JC=1,NCOL
       IPAR=NAC(JC,1)
       JPAR=NBC(JC,1)
       IAA=0
       IBB=0
      IF(NAVAL.EQ.0) GO TO 2
       IF(NUAVAL.NE.IPAR) GO TO 1
       IAA=1
       NASEA=NASEA+1
       JS1=NASEA
       IFLAS(NASEA)=1+INT(RNDM(-1.)/PUD)
C      IFLAS(NASEA)=INT(3.*RNDM(-1.)+1.)
       GO TO 3
  1    CONTINUE
2     NAVAL=NAVAL+1
3     IF(NBVAL.EQ.0) GO TO 5
       IF(NUBVAL.NE.JPAR) GO TO 4
       NBSEA=NBSEA+1
       IBB=1
       JS2=NBSEA
       IFLBS(NBSEA)=1+INT(RNDM(-1.)/PUD)
C      IFLBS(NBSEA)=INT(3.*RNDM(-1.)+1.)
       GO TO 7
  4    CONTINUE
5     NBVAL=NBVAL+1
  7    CONTINUE
       IF(IAA.EQ.0.OR.IBB.EQ.0) GO TO 6
C   DECAY OF SEA-SEA CHAINS
       QSEE = .TRUE.
       QVSEE =.FALSE.
       JVA=MNASEA(JS1)
       JVB=MNBSEA(JS2)
       JS10=JSA0(JS1)
       JS20=JSB0(JS2)
      NIN1=NPTCL+1
      AMQ21=AMQAS1(JS10)
      AMQ22=AMQBS2(JS20)
      CALL XCORR(IFLAS(JS1),-IFLBS(JS2),PXAS1(JS10),
     1PYAS1(JS10),PXBS2(JS20),PYBS2(JS20),
     2XASEA1(JS10),XBSEA2(JS20),PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQAS2(JS10)
      AMQ22=AMQBS1(JS20)
      CALL XCORR(-IFLAS(JS1),IFLBS(JS2),PXAS2(JS10),
     1PYAS2(JS10),PXBS1(JS20),PYBS1(JS20),
     2XASEA2(JS10),XBSEA1(JS20),PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      NFIN1=NPTCL
      DO 500 J=NIN1,NFIN1
      IORIG(J)=7+NCOL*10
      IORDP(J)=0
      IF(DIQAN) IORIG(J)=5+10*(NCOL+2)
500   CONTINUE
      IF(NUNCY.EQ.1) NUNC=NUNC+1
6     CONTINUE
C
      IF(NPOMA.NE.1.OR.NDIFA.NE.1) GO TO 177
      IFLQA2=1+INT(RNDM(-1.)/PUD)
      IFLQA1=-IFLQA2
      IFAQQ=IFAHQ
177   CONTINUE
      IF(NPOMB.NE.1.OR.NDIFB.NE.1) GO TO 178
      IFLQB2=1+INT(RNDM(-1.)/PUD)
      IFLQB1=-IFLQB2
      IFBQQ=IFBHQ
178   CONTINUE
C
C        SELECT VALENCE-VALENCE CHAINS
C
      QSEE=.FALSE.
       NAVAL=0
       NBVAL=0
       NASEA=0
       NBSEA=0
       DO 86 JC=1,NCOL
      IDIFA=0
      IDIFB=0
       IPAR=NAC(JC,1)
       JPAR=NBC(JC,1)
       IAA=0
       IBB=0
      IF(NAVAL.EQ.0) GO TO 92
       IF(NUAVAL.NE.IPAR) GO TO 91
       IAA=1
       GO TO 93
  91   CONTINUE
92    NAVAL=NAVAL+1
  93   IF(NBVAL.EQ.0) GO TO 95
       IF(NUBVAL.NE.JPAR) GO TO 94
       IBB=1
       GO TO 97
  94   CONTINUE
95    NBVAL=NBVAL+1
  97   IF(IAA.EQ.1.OR.IBB.EQ.1) GO TO 86
C   DECAY OF VALENCE-VALENCE CHAINS
      NIN1=NPTCL+1
      IF(DIQAN) GO TO 914
      IF(IFLQA1.GT.0) GO TO 815
      IF(NPOMA.EQ.1.AND.NDIFA.EQ.1) GO TO 951
      IF(NPOMB.EQ.1.AND.NDIFB.EQ.1) GO TO 952
C     IS THERE ANTIDIQUARK
      IF(IFLQB1.GT.0) GO TO 941
      AMQ21=AMQUA1
      AMQ22=AMQQB
      CALL XCORR(IFLQA1,IFBQQ,PXAV1,PYAV1,
     *PXBQQ,PYBQQ,XAVAL1,1.-XBQQ,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQUB1
      CALL XCORR(IFAQQ,IFLQB1,PXAQQ,PYAQQ,
     *PXBV1,PYBV1,1.-XAQQ,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      GO TO 953
941   AMQ21=AMQUA1
      AMQ22=AMQUB1
      CALL XCORR(IFLQA1,IFLQB1,PXAV1,PYAV1,
     *PXBV1,PYBV1,XAVAL1,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQQB
      CALL XCORR(IFAQQ,IFBQQ,PXAQQ,PYAQQ,
     *PXBQQ,PYBQQ,1.-XAQQ,1.-XBQQ,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      GO TO 953
951   CONTINUE
      AMQ21=AMQUA1
      AMQ22=AMQUB1
      CALL XCORR(IFLQA1,IFLQB1,PXAV1,PYAV1,
     *PXBV1,PYBV1,XAVAL1,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQUA2
      AMQ22=AMQQB
      CALL XCORR(IFLQA2,IFBQQ,PXAV2,PYAV2,
     *PXBQQ,PYBQQ,XAVAL2,1.-XBQQ,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      IDENT(NPTCL)=IFAHQ
      IDCAY(NPTCL)=0
      PPTCL(1,NPTCL)=PXAQQ
      PPTCL(2,NPTCL)=PYAQQ
      PPTCL(3,NPTCL)=(1.-XAQQ)*P0
      PPTCL(5,NPTCL)=AMASS(IDENT(NPTCL))
      PPTCL(4,NPTCL)=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+
     *PPTCL(2,NPTCL)**2+PPTCL(3,NPTCL)**2)
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=PPTCL(4,NPTCL)/XAP
      AMT=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+PPTCL(2,NPTCL)**2)
      PPTCL(9,NPTCL)=SQRT(2.D0)*AMT/XAP*PPTCL(4,NPTCL)/PPTCL(5,NPTCL)
      PLIDER(NPTCL)=0.
      GO TO 953
952   CONTINUE
      AMQ21=AMQUA1
      AMQ22=AMQUB2
      CALL XCORR(IFLQA1,IFLQB2,PXAV1,PYAV1,
     *PXBV2,PYBV2,XAVAL1,XBVAL2,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQUB1
      CALL XCORR(IFAQQ,IFLQB1,PXAQQ,PYAQQ,
     *PXBV1,PYBV1,1.-XAQQ,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      IDENT(NPTCL)=IFBHQ
      IDCAY(NPTCL)=0
      PPTCL(1,NPTCL)=PXBQQ
      PPTCL(2,NPTCL)=PYBQQ
      PPTCL(3,NPTCL)=(1.-XBQQ)*P0*PSIGN
      PPTCL(5,NPTCL)=AMASS(IDENT(NPTCL))
      PPTCL(4,NPTCL)=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+
     *PPTCL(2,NPTCL)**2+PPTCL(3,NPTCL)**2)
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=PPTCL(4,NPTCL)/XAP
      AMT=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+PPTCL(2,NPTCL)**2)
      PPTCL(9,NPTCL)=SQRT(2.D0)*AMT/XAP*PPTCL(4,NPTCL)/PPTCL(5,NPTCL)
      PLIDER(NPTCL)=0.
953   NFIN1=NPTCL
      DO 750 J=NIN1,NFIN1
      IORIG(J)=7+10*NCOL
750   CONTINUE
      GO TO 86
815   IF(NPOMB.EQ.1.AND.NDIFB.EQ.1) GO TO 955
      IF(IFLQB1.GT.0) GO TO 943
      AMQ21=AMQUA1
      AMQ22=AMQUB1
      CALL XCORR(IFLQA1,IFLQB1,PXAV1,PYAV1,
     *PXBV1,PYBV1,XAVAL1,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQQB
      CALL XCORR(IFAQQ,IFBQQ,PXAQQ,PYAQQ,
     *PXBQQ,PYBQQ,1.-XAQQ,1.-XBQQ,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      GO TO 956
943   AMQ21=AMQUA1
      AMQ22=AMQQB
      CALL XCORR(IFLQA1,IFBQQ,PXAV1,PYAV1,
     1PXBQQ,PYBQQ,XAVAL1,1.-XBQQ,
     2PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQUB1
      CALL XCORR(IFAQQ,IFLQB1,PXAQQ,
     1PYAQQ,PXBV1,PYBV1,1.-XAQQ,XBVAL1,
     2PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      GO TO 956
955   CONTINUE
      AMQ21=AMQUA1
      AMQ22=AMQUB1
      CALL XCORR(IFLQA1,IFLQB1,PXAV1,PYAV1,
     *PXBV1,PYBV1,XAVAL1,XBVAL1,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQUB2
      CALL XCORR(IFAQQ,IFLQB2,PXAQQ,PYAQQ,
     *PXBV2,PYBV2,1.-XAQQ,XBVAL2,
     *PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      IDENT(NPTCL)=IFBHQ
      IDCAY(NPTCL)=0
      PPTCL(1,NPTCL)=PXBQQ
      PPTCL(2,NPTCL)=PYBQQ
      PPTCL(3,NPTCL)=(1.-XBQQ)*P0*PSIGN
      PPTCL(5,NPTCL)=AMASS(IDENT(NPTCL))
      PPTCL(4,NPTCL)=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+
     *PPTCL(2,NPTCL)**2+PPTCL(3,NPTCL)**2)
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=PPTCL(4,NPTCL)/XAP
      AMT=SQRT(PPTCL(5,NPTCL)**2+PPTCL(1,NPTCL)**2+PPTCL(2,NPTCL)**2)
      PPTCL(9,NPTCL)=SQRT(2.D0)*AMT/XAP*PPTCL(4,NPTCL)/PPTCL(5,NPTCL)
      PLIDER(NPTCL)=0.
956   NFIN1=NPTCL
      DO 700 J=NIN1,NFIN1
      IORIG(J)=7+10*NCOL
700   CONTINUE
      GO TO 86
C     DECAY OF VALENCE-VALENCE CHAINS IN
C         ANNIHILATION CASE:
 914  AMQ21=AMQUA1
      AMQ22=AMQUB1
      CALL XCORR(IFLQA1,IFLQB1,PXAV1,PYAV1,PXBV1,PYBV1,
     *   XAVAL1,XBVAL1,PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQUA2
      AMQ22=AMQUB2
      CALL XCORR(IFLQA2,IFLQB2,PXAV2,PYAV2,PXBV2,PYBV2,
     * XAVAL2,XBVAL2,PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      AMQ21=AMQQA
      AMQ22=AMQQB
      CALL XCORR(IFAQQ,IFBQQ,PXAQQ,PYAQQ,PXBQQ,PYBQQ,
     * XAQQ,XBQQ,PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF(NUNCY.EQ.1) NUNC=NUNC+1
      NFIN1=NPTCL
      DO 950 J=NIN1,NFIN1
      IORIG(J)=7+10*NCOL
      IF(DIQAN) IORIG(J)=5+10*(NCOL+2)
 950  CONTINUE
86    CONTINUE
      NFIN=NPTCL
      DO 861 J=NIN,NFIN
861   IORDP(J)=0
      CALL RESCAL(NIN,NFIN,PSUM,IFAIL)
      IF(IFAIL.EQ.0) RETURN
      GO TO 1111
9999  WRITE(ITLIS,9998) SCM,NPTCL
9998  FORMAT(//10X,'...STOP IN CHAINS..ENERGY TOO LOW SCM=',E10.4/
     *10X,26H..OR NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE REGTRI(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C      COMPUTE TRIPLE REGGEON DIAGRAM
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMLID/PLIDER(499)
      COMMON /PARCUT/ SWMAX
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/COMASS/ AM1,AM2
      COMMON/PRIMP0/ P0
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION V(3)
      DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      LOGICAL BACK
      LOGICAL RETU
      LOGICAL SPINT
      LOGICAL DIQQ
C    INITIALIZE
      NREP=0
      IPACK=1000
      PARBE=0.3
      IRET=0
      IF(ECM.GT.AM1+AM2+   PARBE)
     *GO TO 150
      IRET=1
      GO  TO  9999
150   SPINT=.FALSE.
      DIQQ=.FALSE.
      RETU=.FALSE.
      BACK=.TRUE.
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 102
      IRET=1
C     WRITE(ITLIS,1200) NREP
1200  FORMAT(3X,' IN REGTRI NREP GT ',I8)
      RETURN
102   CONTINUE
      IRET=0
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      CALL FLAVOB(IKA,IFL01,IFL02)
      CALL FLAVOB(IKB,IFL03,IFL04)
      IFL1=IFL01
      IFL2=IFL02
      IFL3=IFL03
      IFL4=IFL04
      RND=RNDM(-1.)
C   COMPUTE X VALUES FOR PARTONS
      IBB=IB(IKB)
      ISB=IS(IKB)
      X3=XDIST(XMIN,IBB,ISB,IFL03)
      X4=1.-X3
      IBA=IB(IKA)
      ISA=IS(IKA)
      IFL01S=IFL01
      IF(ISA.NE.0.AND.IBA.EQ.0.AND.IABS(IFL01).EQ.3) IFL01S=IFL2
      X1=XDIST(XMIN,IBA,ISA,IFL01S)
      X2=1.-X1
C    COMPUTE PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT1,SIGMA)
      AMZER2=AMA**2
      PZER2=P0**2
      AMQ21=AMZER2*(AMZER2+4.*X1*X2*PZER2)/(4.*(AMZER2+PZER2))-PT1**2
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      PX2=-PX1
      PY2=-PY1
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT3,SIGMA)
      AMZER2=AMB**2
      AMQ22=AMZER2*(AMZER2+4.*X3*X4*PZER2)/(4.*(AMZER2+PZER2))-PT3**2
      PX3=PT3*COS(PHI)
      PY3=PT3*SIN(PHI)
      PX4=-PX3
      PY4=-PY3
C    IS THERE ANTIDIQUARK
      IF(IBB.GT.0.AND.MOD(IFL2,100).EQ.0.AND.IFL2.LT.0) GO TO 140
      IF(IBB.GT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.GT.0) GO TO 140
      IF(IBB.LT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 140
      PXH=PX1+PX4
      PYH=PY1+PY4
      X01=1.
      X02=X3
      PX01=PX2
      PY01=PY2
      PX02=PX3
      PY02=PY3
      IDH=IDPARS(IFL1,IFL4,SPINT,0)
      AMH=AMASS(IDH)
      IFL01=IFL2
      IFL02=IFL3
      IF(RND.GE.0.5) GO TO 100
      PXH=PX2+PX3
      PYH=PY2+PY3
      X01=X1
      X02=1.
      IFL01=IFL1
      IFL02=IFL4
      IDH=IDPARS(IFL2,IFL3,SPINT,0)
      AMH=AMASS(IDH)
      PX01=PX1
      PY01=PY1
      PX02=PX4
      PY02=PY4
      GO TO 100
140   PXH=PX1+PX3
      PYH=PY1+PY3
      X01=X2
      X02=1.
      IF(MOD(IABS(IFL2),100).NE.0) GO TO 141
      IF(RNDM(-1.).GT.0.5) GO TO 141
      X01=1.
      X02=X4
141   IDH=IDPARS(IFL1,IFL3,SPINT,0)
      AMH=AMASS(IDH)
      IFL01=IFL2
      IFL02=IFL4
      PX01=PX2
      PY01=PY2
      PX02=PX4
      PY02=PY4
      IF(MOD(IABS(IFL2),100).EQ.0) GO TO 101
      IF(RND.GE.0.5) GO TO 100
      PXH=PX2+PX4
      PYH=PY2+PY4
      X01=1.
      X02=X3
      IDH=IDPARS(IFL2,IFL4,SPINT,0)
      IFL01=IFL1
      IFL02=IFL3
      AMH=AMASS(IDH)
      PX01=PX1
      PY01=PY1
      PX02=PX3
      PY02=PY3
      GO TO 100
101   DIQQ=.TRUE.
      IFCN=1
      IF(RNDM(-1.).GT.0.5) IFCN=2
      IFLC1=-IFCN
      IF(IFL01.GT.0) IFLC1=IFCN
      IFLC2=-IFLC1
      IKH1=IDPARS(IFL01,IFLC1,SPINT,2)
      IKHR1=IDPARS(IFL01,IFLC1,SPINT,1)
      IKH2=IDPARS(IFL02,IFLC2,SPINT,2)
      IKHR2=IDPARS(IFL02,IFLC2,SPINT,1)
      AMHB=AMASS(IKH1)+AMASS(IKH2)+SWMAX
      AMHS=AMASS(IKHR1)+AMASS(IKHR2)+SWMAX
100   P01=X01*P0
      AMQ1=0.
      AMQ2=0.
      P02=X02*P0*PSIGN
      E01=SQRT(AMQ1**2+P01**2+PX01**2+PY01**2)
      E02=SQRT(AMQ2**2+P02**2+PX02**2+PY02**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-(PX01+PX02)**2-
     *(PY01+PY02)**2
      SPINT=.TRUE.
      IF(DIQQ) GO TO 160
      IDH1=IDPARS(IFL01,IFL02,SPINT,2)
      PARBE=0.2
      IF(IABS(IFL01).EQ.3.OR.IABS(IFL02).EQ.3) PARBE=0.3
      AMHB=AMASS(IDH1)+PARBE
160   IF(AMDTR.LE.AMHB**2) GO TO 150
      AMD=SQRT(AMDTR)
      IF(ECM.LE.AMD+AMH) GO TO 150
      ALA=ALAMB(SCM,AMDTR,AMH**2)
      P0H=SQRT(ALA)/(2.*ECM)
      PTHX=-(PX01+PX02)
      PTHY=-(PY01+PY02)
      DTRM=P0H**2-PTHX**2-PTHY**2
      IF(DTRM.LT.0.) GO TO 150
      PZH0=SQRT(DTRM)
      PZH=SIGN(PZH0,-(P01+P02))
      ED=SQRT(AMDTR+P0H**2)
      EH=SQRT(AMH**2+P0H**2)
      PSIGN=SIGN(1.D0,-PZH)
      V(1)=-PTHX/ED
      V(2)=-PTHY/ED
      V(3)=PSIGN*PZH0/ED
      IF(DIQQ) GO TO 170
      IDHR=IDPARS(IFL01,IFL02,SPINT,1)
      AMHS=AMASS(IDHR)+SWMAX
170   IF(AMD.GT.AMHS) GO TO 300
      NFIX=NPTCL
      NIN1=NPTCL+1
      CALL CLUSTR(IFL01,IFL02,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      CALL TIFILL(NIN1,NFIN1,AMD,IFL01,IFL02)
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E01,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
 610  CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      GO TO 400
300   NFIX=NPTCL
      NIN1=NPTCL+1
      CALL STRING(IFL01,IFL02,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      NPRODS=NPTCL-NFIX
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E01,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PPTCL(1,J)
      PPX1(2)=PPTCL(2,J)
      PPX1(3)=PPTCL(3,J)
      CALL ROTR(CT,ST,CFI,SFI,PPX1,PPX2,BACK)
      PPTCL(1,J)=PPX2(1)
      PPTCL(2,J)=PPX2(2)
      PPTCL(3,J)=PPX2(3)
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
510   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
400   CONTINUE
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
500   IORIG(I)=10
      NPTCL=NPTCL+1
      PPTCL(1,NPTCL)=PTHX
      PPTCL(2,NPTCL)=PTHY
      PPTCL(5,NPTCL)=AMH
      PPTCL(3,NPTCL)=PZH
      PPTCL(4,NPTCL)=EH
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IDH
      IORIG(NPTCL)=10
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=0
      RETURN
9999  CONTINUE
C     WRITE(ITLIS,1000) SCM
1000  FORMAT(//10X,39H...STOP IN REGTRI...ENERGY TOO LOW SCM=,E10.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE UNCYLI(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C      COMPUTE UNDEVELOPED CYLINDER DIAGRAM
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMLID/PLIDER(499)
      COMMON /PARCUT/ SWMAX
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/PRIMP0/ P0
      COMMON/COMASS/ AM1,AM2
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION V(3)
      DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      LOGICAL BACK
      LOGICAL RETU
      LOGICAL SPINT
      LOGICAL DIQQ
C    INITIALIZE
      NREP=0
      IPACK=1000
      PARBE=0.3
      IRET=0
      IF(ECM.GT.AM1+AM2+PARBE)
     *GO TO 150
      IRET=1
      GO  TO  9999
150   SPINT=.FALSE.
      DIQQ=.FALSE.
      RETU=.FALSE.
      BACK=.TRUE.
      PARBE=0.2
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 102
      IRET=1
C     WRITE(ITLIS,1200) NREP
1200  FORMAT(3X,' IN UNCYLI NREP GT ',I8)
      RETURN
102   CONTINUE
      IRET=0
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      CALL FLAVOB(IKA,IFL01,IFL02)
      CALL FLAVOB(IKB,IFL03,IFL04)
      IFL1=IFL01
      IFL2=IFL02
      IFL3=IFL03
      IFL4=IFL04
      RND=RNDM(-1.)
C   COMPUTE X VALUES FOR PARTONS
      IBB=IB(IKB)
      ISB=IS(IKB)
      X3=XDIST(XMIN,IBB,ISB,IFL03)
      X4=1.-X3
      IBA=IB(IKA)
      ISA=IS(IKA)
      IFL01S=IFL01
      IF(ISA.NE.0.AND.IBA.EQ.0.AND.IABS(IFL01).EQ.3) IFL01S=IFL2
      X1=XDIST(XMIN,IBA,ISA,IFL01S)
      X2=1.-X1
C    COMPUTE PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT1,SIGMA)
      AMZER2=AMA**2
      PZER2=P0**2
      AMQ21=AMZER2*(AMZER2+4.*X1*X2*PZER2)/(4.*(AMZER2+PZER2))-PT1**2
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      PX2=-PX1
      PY2=-PY1
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT3,SIGMA)
      AMZER2=AMB**2
      AMQ22=AMZER2*(AMZER2+4.*X3*X4*PZER2)/(4.*(AMZER2+PZER2))-PT3**2
      PX3=PT3*COS(PHI)
      PY3=PT3*SIN(PHI)
      PX4=-PX3
      PY4=-PY3
C    IS THERE ANTIDIQUARK
      IF(IBB.GT.0.AND.MOD(IFL2,100).EQ.0.AND.IFL2.LT.0) GO TO 140
      IF(IBB.GT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.GT.0) GO TO 140
      IF(IBB.LT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 140
      PXH=PX1+PX4
      PYH=PY1+PY4
      X01=X2
      X02=X3
      PX01=PX2
      PY01=PY2
      PX02=PX3
      PY02=PY3
      IDH=IDPARS(IFL1,IFL4,SPINT,0)
      AMH=AMASS(IDH)
      IFL01=IFL2
      IFL02=IFL3
      IDH1=IDPARS(IFL1,IFL4,SPINT,1)
      AMH1=AMASS(IDH1)
      P01=X1*P0
      P02=X4*P0*PSIGN
      E01=SQRT(P01**2+PX1**2+PY1**2)
      E02=SQRT(P02**2+PX4**2+PY4**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-PXH**2-PYH**2
      AMHB=AMH1+PARBE
      IF(AMDTR.LE.AMHB**2) GO TO 100
      PXH=PX2+PX3
      PYH=PY2+PY3
      X01=X1
      X02=X4
      IFL01=IFL1
      IFL02=IFL4
      IDH=IDPARS(IFL2,IFL3,SPINT,0)
      AMH=AMASS(IDH)
      PX01=PX1
      PY01=PY1
      PX02=PX4
      PY02=PY4
      IDH1=IDPARS(IFL2,IFL3,SPINT,1)
      AMH1=AMASS(IDH1)
      AMHB=AMH1+PARBE
      P01=X2*P0
      P02=X3*P0*PSIGN
      E01=SQRT(P01**2+PX2**2+PY2**2)
      E02=SQRT(P02**2+PX3**2+PY3**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-PXH**2-PYH**2
      IF(AMDTR.LE.AMHB**2) GO TO 100
      GO TO 150
140   PXH=PX1+PX3
      PYH=PY1+PY3
      X01=X2
      X02=X4
      IDH=IDPARS(IFL1,IFL3,SPINT,0)
      AMH=AMASS(IDH)
      IFL01=IFL2
      IFL02=IFL4
      PX01=PX2
      PY01=PY2
      PX02=PX4
      PY02=PY4
      IDH1=IDPARS(IFL1,IFL3,SPINT,1)
      AMH1=AMASS(IDH1)
      AMHB=AMH1+PARBE
      P01=X1*P0
      P02=X3*P0*PSIGN
      E01=SQRT(P01**2+PX1**2+PY1**2)
      E02=SQRT(P02**2+PX3**2+PY3**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-PXH**2-PYH**2
      IF(AMDTR.LE.AMHB**2.AND.MOD(IABS(IFL2),100).EQ.0) GO TO 101
      IF(AMDTR.LE.AMHB**2) GO TO 100
      IF(MOD(IFL2,100).EQ.0.AND.MOD(IFL4,100).EQ.0) GO TO 150
      PXH=PX2+PX4
      PYH=PY2+PY4
      X01=X1
      X02=X3
      IDH=IDPARS(IFL2,IFL4,SPINT,0)
      IFL01=IFL1
      IFL02=IFL3
      AMH=AMASS(IDH)
      PX01=PX1
      PY01=PY1
      PX02=PX3
      PY02=PY3
      IDH1=IDPARS(IFL2,IFL4,SPINT,1)
      AMH1=AMASS(IDH1)
      AMHB=AMH1+PARBE
      P01=X2*P0
      P02=X4*P0*PSIGN
      E01=SQRT(P01**2+PX2**2+PY2**2)
      E02=SQRT(P02**2+PX4**2+PY4**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-PXH**2-PYH**2
      IF(AMDTR.LE.AMHB**2) GO TO 100
      GO TO 150
101   DIQQ=.TRUE.
      IFCN=1
      IF(RNDM(-1.).GT.0.5) IFCN=2
      IFLC1=-IFCN
      IF(IFL01.GT.0) IFLC1=IFCN
      IFLC2=-IFLC1
      IKH1=IDPARS(IFL01,IFLC1,SPINT,2)
      IKHR1=IDPARS(IFL01,IFLC1,SPINT,1)
      IKH2=IDPARS(IFL02,IFLC2,SPINT,2)
      IKHR2=IDPARS(IFL02,IFLC2,SPINT,1)
      AMHB=AMASS(IKH1)+AMASS(IKH2)+SWMAX
      AMHS=AMASS(IKHR1)+AMASS(IKHR2)+SWMAX
100   P01=X01*P0
      AMQ1=0.
      AMQ2=0.
      P02=X02*P0*PSIGN
      E01=SQRT(AMQ1**2+P01**2+PX01**2+PY01**2)
      E02=SQRT(AMQ2**2+P02**2+PX02**2+PY02**2)
      AMDTR=(E01+E02)**2-(P01+P02)**2-(PX01+PX02)**2-
     *(PY01+PY02)**2
      SPINT=.TRUE.
      IF(DIQQ) GO TO 160
      IDH1=IDPARS(IFL01,IFL02,SPINT,2)
      PARBE=0.2
      IF(IABS(IFL01).EQ.3.OR.IABS(IFL02).EQ.3) PARBE=0.3
      AMHB=AMASS(IDH1)+PARBE
160   IF(AMDTR.LE.AMHB**2) GO TO 150
      AMD=SQRT(AMDTR)
      IF(ECM.LE.AMD+AMH) GO TO 150
      ALA=ALAMB(SCM,AMDTR,AMH**2)
      P0H=SQRT(ALA)/(2.*ECM)
      PTHX=-(PX01+PX02)
      PTHY=-(PY01+PY02)
      DTRM=P0H**2-PTHX**2-PTHY**2
      IF(DTRM.LT.0.) GO TO 150
      PZH0=SQRT(DTRM)
      PZH=SIGN(PZH0,-(P01+P02))
      ED=SQRT(AMDTR+P0H**2)
      EH=SQRT(AMH**2+P0H**2)
      PSIGN=SIGN(1.D0,-PZH)
      V(1)=-PTHX/ED
      V(2)=-PTHY/ED
      V(3)=PSIGN*PZH0/ED
      IF(DIQQ) GO TO 170
      IDHR=IDPARS(IFL01,IFL02,SPINT,1)
      AMHS=AMASS(IDHR)+SWMAX
170   IF(AMD.GT.AMHS) GO TO 300
      NFIX=NPTCL
      NIN1=NPTCL+1
      CALL CLUSTR(IFL01,IFL02,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      CALL TIFILL(NIN1,NFIN1,AMD,IFL01,IFL02)
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E01,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
610   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      GO TO 400
300   NFIX=NPTCL
      NIN1=NPTCL+1
      CALL STRING(IFL01,IFL02,AMD)
      IF(LRET) GO TO 150
      NFIN1=NPTCL
      NPRODS=NPTCL-NFIX
      PPX1(1)=PX01
      PPX1(2)=PY01
      PPX1(3)=P01
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E01,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PPTCL(1,J)
      PPX1(2)=PPTCL(2,J)
      PPX1(3)=PPTCL(3,J)
      CALL ROTR(CT,ST,CFI,SFI,PPX1,PPX2,BACK)
      PPTCL(1,J)=PPX2(1)
      PPTCL(2,J)=PPX2(2)
      PPTCL(3,J)=PPX2(3)
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
510   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
400   CONTINUE
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
      IORIG(I)=3
500   CONTINUE
      NPTCL=NPTCL+1
      PPTCL(1,NPTCL)=PTHX
      PPTCL(2,NPTCL)=PTHY
      PPTCL(5,NPTCL)=AMH
      PPTCL(3,NPTCL)=PZH
      PPTCL(4,NPTCL)=EH
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDENT(NPTCL)=IDH
      IORIG(NPTCL)=3
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=0
      RETURN
9999  CONTINUE
C     WRITE(ITLIS,1000) SCM
1000  FORMAT(//10X,39H...STOP IN UNCYLI...ENERGY TOO LOW SCM=,E10.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CYLIN(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE CYLINDER TYPE DIAGRAM
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /PARCUT/ SWMAX
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMXM/ XMIN,XMAX
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/NEEDR/ NRET
      COMMON/COMASS/ AM1,AM2
      COMMON/PRIMP0/ P0
      DIMENSION PSUM(5)
      LOGICAL RETU
C   INITIALIZE
      IPACK=1000
      NRET=1
      NREP=0
      PARBE=0.3
      IF(ECM.GT.AM1+AM2+   PARBE) GO TO 999
      IRET=1
      GO  TO  9999
999   CONTINUE
      DO 151 I=1,3
 151  PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
100   IRET=0
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 101
      IRET=1
      RETURN
101   CONTINUE
      RETU=.FALSE.
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      AMQ21=0.
      AMQ22=0.
      CALL FLAVOB(IKA,IFL01,IFL02)
      CALL FLAVOB(IKB,IFL03,IFL04)
      IFL1=IFL01
      IFL2=IFL02
      IFL3=IFL03
      IFL4=IFL04
C   COMPUTE X VALUES FOR PARTONS
      IBA=IB(IKA)
      ISA=IS(IKA)
      IFL01S=IFL01
      IF(ISA.NE.0.AND.IBA.EQ.0.AND.IABS(IFL01).EQ.3) IFL01S=IFL2
      X1=XDIST(XMIN,IBA,ISA,IFL01S)
      X2=1.-X1
      IBB=IB(IKB)
      ISB=IS(IKB)
      X3=XDIST(XMIN,IBB,ISB,IFL03)
      X4=1.-X3
      IF(IBA.EQ.0.AND.IBB.EQ.0) NRET=0
C   COMPUTE PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT1,SIGMA)
      AMZER2=AMA**2
      PZER2=P0**2
      AMQ21=AMZER2*(AMZER2+4.*X1*X2*PZER2)/(4.*(AMZER2+PZER2))-PT1**2
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      PX2=-PX1
      PY2=-PY1
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT3,SIGMA)
      AMZER2=AMB**2
      AMQ22=AMZER2*(AMZER2+4.*X3*X4*PZER2)/(4.*(AMZER2+PZER2))-PT3**2
      PX3=PT3*COS(PHI)
      PY3=PT3*SIN(PHI)
      PX4=-PX3
      PY4=-PY3
C    IS THERE ANTIDIQUARK
      NIN1=NPTCL+1
      IF(IBB.GT.0.AND.MOD(IFL2,100).EQ.0.AND.IFL2.GT.0) GO TO 150
      IF(IBB.GT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.LT.0) GO TO 150
      IF(IBB.LT.0.AND.IBA.LT.0) GO TO 150
      IF(IBB.LT.0.AND.MOD(IFL2,100).NE.0.AND.IFL2.GT.0) GO TO 150
      IF(.NOT.(IBA.EQ.0.AND.IBB.EQ.0)) GO TO 160
      IF(IFL1*IFL3.GT.0) GO TO 150
160   CALL XCORR(IFL1,IFL3,PX1,PY1,PX3,PY3,X1,X3,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 100
      CALL XCORR(IFL2,IFL4,PX2,PY2,PX4,PY4,X2,X4,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      GO TO 100
130   NFIN1=NPTCL
      DO 550 I=NIN1,NFIN1
550   IORIG(I)=7
C@@@@@@@@ 14.08.91   SIVOKL
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) RETURN
      NPTCL=NPTCL-NPRD
      GO TO 100
C*@@@ RETURN
150   CALL XCORR(IFL1,IFL4,PX1,PY1,PX4,PY4,X1,X4,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 100
      CALL XCORR(IFL2,IFL3,PX2,PY2,PX3,PY3,X2,X3,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 140
      NPTCL=NPTCL-NPRD
      GO TO 100
140   NFIN1=NPTCL
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
500   IORIG(I)=7
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) RETURN
      NPTCL=NPTCL-NPRD
      GO TO 100
9999  CONTINUE
C     WRITE(ITLIS,1000) SCM
1000  FORMAT(//10X,'...STOP IN CYLIN..ENERGY TOO LOW SCM=',E10.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE TIFILL(N1,N2,AMS,IFL1,IFL2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C----  COMPUTE ZI Z-COORD. & TI TIME OF HADRONS AFTER CLUSTER DECAY
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/KAPPA/XAP
      COMMON /CINSID/ INSIDE
      COMMON/COMLID/PLIDER(499)
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      IF(INSIDE.NE.0)   GO  TO  1
C-----  CONSTITUENT   TIME -----------------C
      TI=(AMS-2.*PPTCL(3,N1))/(2.*XAP)
      ZI=(AMS-2.*PPTCL(4,N1))/(2.*XAP)
      PPTCL(6,N1)=0.
      PPTCL(7,N1)=0.
      PPTCL(8,N1)=ZI
      PPTCL(9,N1)=TI
      PPTCL(6,N2)=0.
      PPTCL(7,N2)=0.
      PPTCL(8,N2)=ZI
      PPTCL(9,N2)=TI
      GO  TO  2
   1  CONTINUE
C-----  INSIDE - OUTSIDE TIME -----------------C
      T1=(AMS+PPTCL(4,N1)-PPTCL(3,N1))/(2.*XAP)
      T2=(AMS-2.*PPTCL(3,N1)+PPTCL(4,N2)-PPTCL(3,N2))/(2.*XAP)
      Z1=(AMS-PPTCL(4,N1)+PPTCL(3,N1))/(2.*XAP)
      Z2=(AMS-2.*PPTCL(4,N1)-PPTCL(4,N2)+PPTCL(3,N2))/(2.*XAP)
      PPTCL(6,N1)=0.
      PPTCL(7,N1)=0.
      PPTCL(8,N1)=Z1
      PPTCL(9,N1)=T1
      PPTCL(6,N2)=0.
      PPTCL(7,N2)=0.
      PPTCL(8,N2)=Z2
      PPTCL(9,N2)=T2
   2  CONTINUE
C-------------------------------------------C
      PLIDER(N1)=0.
      PLIDER(N2)=0.
      IF(QSEE) RETURN
      IB1=IB(IDENT(N1))
      IB2=IB(IDENT(N2))
      PLIDER(N1)=.667
      PLIDER(N2)=.667
      IF(IB1.EQ.0) PLIDER(N1)=.5
      IF(IB2.EQ.0) PLIDER(N2)=.5
      IF(PLIDER(N1).GT.0.6.AND.MOD(IFL1,100).NE.0) PLIDER(N1)=0.333
      IF(PLIDER(N2).GT.0.6.AND.MOD(IFL2,100).NE.0) PLIDER(N2)=0.333
      IF(.NOT.QVSEE) RETURN
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 387
      IF(IB1.EQ.0) PLIDER(N1)=0.
      IF(IB2.EQ.0) PLIDER(N2)=0.
      IF(PLIDER(N1).GT.0.6.AND.MOD(IFL1,100).NE.0) PLIDER(N1)=0.333
      IF(PLIDER(N2).GT.0.6.AND.MOD(IFL2,100).NE.0) PLIDER(N2)=0.333
      RETURN
387   RM=RNDM(-1.)
      IF(RM.GT.0.5) PLIDER(N1)=0.
      IF(RM.LE.0.5) PLIDER(N2)=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE TIFILE(N1,N2,AMS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C----  COMPUTE ZI Z-COORD. & TI TIME OF HADRONS AFTER CLUSTER DECAY
      COMMON /PROD/ PR(8,50),IPR(50),NP
      COMMON/KAPPA/XAP
      COMMON /CINSID/ INSIDE
      COMMON/COMLD/PLDER(50)
      COMMON/COMQSE/QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      IF(INSIDE.NE.0)  GO  TO  1
C-----  CONSTITUENT      TIME -----------------C
      TI=(AMS-2.*PR(3,N1))/(2.*XAP)
      ZI=(AMS-2.*PR(4,N1))/(2.*XAP)
      PR(5,N1)=0.
      PR(6,N1)=0.
      PR(7,N1)=ZI
      PR(8,N1)=TI
      PR(5,N2)=0.
      PR(6,N2)=0.
      PR(7,N2)=ZI
      PR(8,N2)=TI
      GO  TO  2
C-----  INSIDE - OUTSIDE TIME -----------------C
   1  CONTINUE
      T1=(AMS+PR(4,N1)-PR(3,N1))/(2.*XAP)
      T2=(AMS-2.*PR(3,N1)+PR(4,N2)-PR(3,N2))/(2.*XAP)
      Z1=(AMS-PR(4,N1)+PR(3,N1))/(2.*XAP)
      Z2=(AMS-2.*PR(4,N1)-PR(4,N2)+PR(3,N2))/(2.*XAP)
      PR(5,N1)=0.
      PR(6,N1)=0.
      PR(7,N1)=Z1
      PR(8,N1)=T1
      PR(5,N2)=0.
      PR(6,N2)=0.
      PR(7,N2)=Z2
      PR(8,N2)=T2
   2  CONTINUE
C-----------------------------------------------C
      IB1=IBLE(IPR(N1))
      IB2=IBLE(IPR(N2))
      PLDER(N1)=.667
      PLDER(N2)=.667
      IF(IB1.EQ.0) PLDER(N1)=.5
      IF(IB2.EQ.0) PLDER(N2)=.5
      IF(.NOT.QVSEE) RETURN
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 387
      IF(IB1.EQ.0) PLDER(N1)=0.
      IF(IB2.EQ.0) PLDER(N2)=0.
      RETURN
387   RM=RNDM(-1.)
      IF(RM.GT.0.5) PLDER(N1)=0.
      IF(RM.LE.0.5) PLDER(N2)=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c************** last correction 12-21-95 05:40pm*************

! =====================================================================
! READHH removed by CMJ (XCP-3, LANL) on 09/08/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: UNKNOWN (purpose not commented)
! =====================================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE   FORID(IDOLD,IDNEW)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      DIMENSION IDENT(54)
      DATA IDENT/
     *     120 ,  -120 ,   130 ,  -130 ,   230 ,  -230 ,   110 ,
     *     220 ,   330 ,   121 ,  -121 ,   131 ,  -131 ,   231 ,
     *    -231 ,   111 ,   221 ,   331 ,
     *   -1120 , -1220 , -1130 , -2230 , -1230 , -2330 , -1330 ,
     *   -2130 , -1111 , -1121 , -2221 , -1221 , -1131 , -2231 ,
     *   -1231 , -2331 , -1331 , -3331 ,
     *    1120 ,  1220 ,  1130 ,  2230 ,  1230 ,  2330 ,  1330 ,
     *    2130 ,  1111 ,  1121 ,  2221 ,  1221 ,  1131 ,  2231 ,
     *    1231 ,  2331 ,  1331 ,  3331 /
      N=0
      DO  10  I=1,54
      IF(      I .NE.IDOLD)  GO TO 10
      N=I
      GO TO 20
   10 CONTINUE
   20 IF(N.EQ.0) WRITE(ITLIS,991)IDOLD
  991 FORMAT(2X,'FORID:  IDOLD =',I5)
      IDNEW=IDENT(N)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE   BACKID(IDOLD,IDNEW)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
C
      CHARACTER*8 PNAME,LAB
      DIMENSION PNAME(54)
      DATA PNAME/
     *  'PI+  ','PI-  ','K+   ','K-   ','K0   ','AK0  ','PI0  ',
     *  'ETA  ','ETAP ','RHO+ ','RHO- ','K*+  ','K*-  ','K*0  ',
     *  'AK*0 ','RHO0 ','OMEG ','PHI  ',
     *  'AP   ','AN   ','AS+  ','AS-  ','AS0  ','AXI- ','AXI0 ',
     *  'AL   ','ADL++','ADL+ ','ADL- ','ADL0 ','AS*+ ','AS*- ',
     *  'AS*0 ','AXI*-','AXI*0','AOM- ',
     *  'P    ','N    ','S+   ','S-   ','S0   ','XI-  ','XI0  ',
     *  'L    ','DL++ ','DL+  ','DL-  ','DL0  ','S*+  ','S*-  ',
     *  'S*0  ','XI*- ','XI*0 ','OM-  '/
C
      CALL LABEL(LAB,IDOLD)
      N=0
      DO 100 I=1,54
      IF(LAB.NE.PNAME(I)) GO TO 100
      N=I
100   CONTINUE
      IDNEW=N
      IF(N.EQ.0) WRITE(ITLIS,991) IDOLD
991   FORMAT(10X,'BACKID;IDOLD=',I6)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE MARK(IK01,IK02,KS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      IK1=IK01
      IK2=IK02
      IB1=IB(IK1)
      IB2=IB(IK2)
      IF(IB1.GE.0.OR.IB2.GE.0) GO TO 110
C     AB-AB -> B-B
      IK1=IABS(IK1)
      IK2=IABS(IK2)
      IF(IK2.EQ.1120.OR.IK2.EQ.1220) GO TO 112
      IF(IK1.NE.1120.AND.IK1.NE.1220) GO TO 111
      IK11=IK1
      IK1=IK2
      IK2=IK11
      GO TO 112
110   IF(IB1.NE.0.OR.IB2.GE.0) GO TO 111
C     M-AB -> AM-B
      IK2=IABS(IK2)
      IF(IK1/100.NE.MOD(IK1,100)/10) IK1=-IK1
111   IF(IK2.NE.1220) IK2=1120
112   IB1=IB(IK1)
      IB2=IB(IK2)
      IQ1=INT(CHARGE(IK1)*1.001)
      IQ2=INT(CHARGE(IK2)*1.001)
      MQ=IQ1+IQ2
      IF(IB1+IB2.LE.1) GO TO 1
C   NUCLEON-NUCLEON COLLISION
      IF(MQ-1) 3,4,3
C  MESON-NUCLEON COLLISION
 1    IF(MQ.EQ.2.OR.MQ.EQ.-1) GO TO 3
      IF(MQ.EQ.0) GO TO 2
      IF(IQ1-1) 5,4,5
 2    IF(IQ1+1) 5,4,5
 3    KS=1
      RETURN
 4    KS=2
      RETURN
 5    KS=3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! comp removed by CMJ (XCP-3, LANL) on 09/07/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: Calculates Px/y/z from Ptot and cos(theta), and (phi)
! =====================================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSDD(I)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
 100  R1=RNDM(-1.)
      CT=2.*R1-1.
      IF(I.EQ.0) GO TO 2
      R2=RNDM(-1.)
      W=0.25+0.75*CT**2
      IF(W.LT.R2) GO TO 100
 2    COSDD=CT
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION DBLPCM(A,B,C)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      VAL=(A**2-B**2-C**2)**2-(2.D0*B*C)**2
      DBLPCM=0.
      IF(VAL.GT.0.D0)
     *DBLPCM=SQRT(VAL)/(2.*A)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION AMASSF(ID)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C          THIS FUNCTION RETURNS THE MASS OF THE PARTICLE WITH
C          IDENT CODE ID.
C          QUARK-BASED IDENT CODE
      DIMENSION AMMES0(10),AMMES1(10),AMBAR0(30),AMBAR1(30)
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
C          0- MESON MASS TABLE
      DATA AMMES0/.13496,.13957,.5488,.49367,.49767,.9576,1.8633
     1,1.8683,2.030,2.976/
C          1- MESON MASS TABLE
      DATA AMMES1/.776,.776,.7826,.8881,.8922,1.0196,2.006,2.0086
     1,2.140,3.097/
C          1/2+ BARYON MASS TABLE
      DATA AMBAR0/-1.,.93828,.93957,2*-1.,1.1894,1.1925,1.1974
     1,1.1156,1.3149,1.3213,3*-1.,2.43,2.43,2.43,2.26
     2,2.50,2.50,2.60,2.40,2.40,3.55,3.55,3.70,4*-1./
C          3/2+ BARYON MASS TABLE
      DATA AMBAR1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     1,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     2,2.63,2.63,2.63,-1.,2.70,2.70,2.80,2*-1.,3.75,3.75
     3,3.90,4.80,3*-1./
C          ENTRY
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 400
      IF(IABS(IFL1).GT.4.OR.IABS(IFL2).GT.4.OR.IABS(IFL3).GT.4)
     1GO TO 300
      IF(IFL2.EQ.0) GO TO 200
      IF(IFL1.EQ.0) GO TO 100
C          BARYONS
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-11
      AMASSF=(1-JSPIN)*AMBAR0(INDEX)+JSPIN*AMBAR1(INDEX)
      RETURN
C          MESONS
100   CONTINUE
      INDEX=INDEX-36*JSPIN-NQLEP
      INDEX=INDEX-11
      AMASSF=(1-JSPIN)*AMMES0(INDEX)+JSPIN*AMMES1(INDEX)
      RETURN
C          QUARKS AND LEPTONS
200   CONTINUE
      AMASSF=AMLEP(INDEX)
      RETURN
C          B AND T PARTICLES
300   CONTINUE
      AMASSF=AMLEP(IABS(IFL2))+AMLEP(IABS(IFL3))-.03+.04*JSPIN
      IF(IFL1.NE.0) AMASSF=AMASSF+AMLEP(IABS(IFL1))
      RETURN
C          DIQUARKS
400   AMASSF=AMLEP(IABS(IFL1))+AMLEP(IABS(IFL2))+0.5
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION AMASS(ID)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C          THIS FUNCTION RETURNS THE MASS OF THE PARTICLE WITH
C          IDENT CODE ID.
C          QUARK-BASED IDENT CODE
      DIMENSION AMMES0(10),AMMES1(10),AMBAR0(30),AMBAR1(30)
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
      COMMON/ISOB3/ISOB3
      COMMON/INTTYP/ITYP
      COMMON/ITHEA/ITHEA(11)
C          0- MESON MASS TABLE
      DATA AMMES0/.13496,.13957,.5488,.49367,.49767,.9576,1.8633
     1,1.8683,2.030,2.976/
C          1- MESON MASS TABLE
      DATA AMMES1/.776,.776,.7826,.8881,.8922,1.0196,2.006,2.0086
     1,2.140,3.097/
C          1/2+ BARYON MASS TABLE
      DATA AMBAR0/-1.,.93828,.93957,2*-1.,1.1894,1.1925,1.1974
     1,1.1156,1.3149,1.3213,3*-1.,2.43,2.43,2.43,2.26
     2,2.50,2.50,2.60,2.40,2.40,3.55,3.55,3.70,4*-1./
C          3/2+ BARYON MASS TABLE
      DATA AMBAR1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     1,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     2,2.63,2.63,2.63,-1.,2.70,2.70,2.80,2*-1.,3.75,3.75
     3,3.90,4.80,3*-1./
C          ENTRY
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 400
      IF(IABS(IFL1).GT.4.OR.IABS(IFL2).GT.4.OR.IABS(IFL3).GT.4)
     1GO TO 300
      IF(IFL2.EQ.0) GO TO 200
      IF(IFL1.EQ.0) GO TO 100
C          BARYONS
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-11
      AMASS=(1-JSPIN)*AMBAR0(INDEX)+JSPIN*AMBAR1(INDEX)
      IF(ISOB3 /= 1)  RETURN
      IF(ID.EQ.1111.OR.ID.EQ.1121.OR.ID.EQ.2221.OR.ID.EQ.1221)
     *GO  TO  1991
      RETURN
1991  CONTINUE
C     IF((ITHEA(8).EQ.1).OR.(ITYP.EQ.3))  THEN
C           CALL  MDELT1(AMD,GD)
C     ELSE
            CALL  MDELTA(AMD,GD)
C     ENDIF
      AMASS=AMD
      RETURN
C          MESONS
100   CONTINUE
      INDEX=INDEX-36*JSPIN-NQLEP
      INDEX=INDEX-11
      AMASS=(1-JSPIN)*AMMES0(INDEX)+JSPIN*AMMES1(INDEX)
      IF(ISOB3.NE.1)  RETURN
      IF(ID.EQ.10.OR.ID.EQ.11.OR.ID.EQ.16)
     *GO  TO  1992
      RETURN
 1992       CALL  MRHO(AMRHO,GD)
      AMASS=AMRHO
      RETURN
C          QUARKS AND LEPTONS
200   CONTINUE
      AMASS=AMLEP(INDEX)
      RETURN
C          B AND T PARTICLES
300   CONTINUE
      AMASS=AMLEP(IABS(IFL2))+AMLEP(IABS(IFL3))-.03+.04*JSPIN
      IF(IFL1.NE.0) AMASS=AMASS+AMLEP(IABS(IFL1))
      RETURN
C          DIQUARKS
400   AMASS=AMLEP(IABS(IFL1))+AMLEP(IABS(IFL2))+0.5
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        DOUBLE PRECISION FUNCTION GDM(X)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
c
        DATA AMN/0.940/,AMPI/0.140/,B/0.300/,GD0/0.110/,AMD0/1.232/
        AMD=X
        EN0=(AMD0**2+AMN**2-AMPI**2)/(2.*AMD0)
        Q0=SQRT(EN0**2-AMN**2)
        EN =(AMD**2+AMN**2-AMPI**2)/(2.*AMD )
        QX=      EN**2-AMN**2
        IF(QX.LT.0.)  WRITE(ITLIS,*) 'GDM: QX,AMD=', QX,AMD
        Q =SQRT(ABS(QX))
        V0=B**2/(B**2+Q0**2)
        V =B**2/(B**2+Q**2)
        GDM=(Q/Q0)**3*(AMD0/AMD)*(V/V0)**2*GD0
        RETURN
        END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        SUBROUTINE MDELTA(AMDEL,GD)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
        DATA AMD0/1.232/
        DATA AMDMIN/1.081/,AMDMAX/1.700/
   10 CONTINUE
        AMD=AMDMIN+RNDM(-1.)*(AMDMAX-AMDMIN)
        GD=GDM(AMD)
        F=0.25*GD**2/((AMD-AMD0)**2+0.25*GD**2)
      IF(RNDM(-1.).GT.F)  GO  TO  10
        AMDEL=AMD
        RETURN
        END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
           SUBROUTINE LABEL(LABEL1,ID)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C          RETURN THE LABEL FOR THE PARTICLE ID.
C          QUARK-BASED IDENT CODE.
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
C
      CHARACTER*8   LABEL1, LLEP(104)
     *,LMES0(64),LMES1(64)
     *,LBAR0(109),LABAR0(109),LBAR1(109),LABAR1(109)
     *,LQQ(21),LAQQ(21)
C          DIQUARK LABELS
      DATA LQQ/
     1'UU0. ','UD0. ','DD0. ','US0. ','DS0. ','SS0. ','UC0. ','DC0. ',
     2'SC0. ','CC0. ','UB0. ','DB0. ','SB0. ','CB0. ','BB0. ','UT0. ',
     3'DT0. ','ST0. ','CT0. ','BT0. ','TT0. '/
      DATA LAQQ/
     1'AUU0.','AUD0.','ADD0.','AUS0.','ADS0.','ASS0.','AUC0.','ADC0.',
     2'ASC0.','ACC0.','AUB0.','ADB0.','ASB0.','ACB0.','ABB0.','AUT0.',
     3'ADT0.','AST0.','ACT0.','ABT0.','ATT0.'/
C          QUARK AND LEPTON LABELS
      DATA LLEP/
     $'     ','UP   ','UB   ','DN   ','DB   ','ST   ','SB   ','CH   ',
     $'CB   ','BT   ','BB   ','TP   ','TB   ','Y    ','YB   ','X    ',
     $'XB   ','GL   ','ERR  ','GM   ','ERR  ','NUE  ','ANUE ','E-   ',
     $'E+   ','NUM  ','ANUM ','MU-  ','MU+  ','NUT  ','ANUT ','TAU- ',
     $'TAU+ ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','KS   ',
     $'ERR  ','ERR  ','KL   ',
     $'UPSS ','UBSS ','DNSS ','DBSS ','STSS ','SBSS ','CHSS ','CBSS ',
     $'BTSS ','BBSS ','TPSS ','TBSS ','ERR  ','ERR  ','ERR  ','ERR  ',
     $'GLSS ','ERR  ','GMSS ','ERR  ','NESS ','ANESS','E-SS ','E+SS ',
     $'NMSS ','ANMSS','MU-SS','MU+SS','NTSS ','ANTSS','T-SS ','T+SS ',
     $'ERR  ','ERR  ','ERR  ','ERR  ','W+SS ','W-SS ','Z0SS ','ERR  ',
     $'W+   ','W-   ','H10  ','AH10 ','H20  ','AH20 ','H30  ','AH30 ',
     $'H4+  ','H4-  ','H5+  ','H5-  ','H6+  ','H6-  ','H7++ ','H7-- ',
     $'H8++ ','H8-- ','H9++ ','H9-- ','Z0   '/
C          0- MESON LABELS
      DATA LMES0/
     1'PI0  ','PI+  ','ETA  ','PI-  ','K+   ','K0   ','ETAP ','AK0  ',
     2'K-   ','AD0  ','D-   ','F-   ','ETAC ','F+   ','D+   ','D0   ',
     2'UB.  ','DB.  ','SB.  ','CB.  ','BB.  ','BC.  ','BS.  ','BD.  ',
     3'BU.  ','UT.  ','DT.  ','ST.  ','CT.  ','BT.  ','TT.  ','TB.  ',
     4'TC.  ','TS.  ','TD.  ','TU.  ','UY.  ','DY.  ','SY.  ','CY.  ',
     5'BY.  ','TY.  ','YY.  ','YT.  ','YB.  ','YC.  ','YS.  ','YD.  ',
     6'YU.  ','UX.  ','DX.  ','SX.  ','CX.  ','BX.  ','TX.  ','YX.  ',
     7'XX.  ','XY.  ','XT.  ','XB.  ','XC.  ','XS.  ','XD.  ','XU.  '/
C          1- MESON LABELS
      DATA LMES1/
     1'RHO0 ','RHO+ ','OMEG ','RHO- ','K*+  ','K*0  ','PHI  ','AK*0 ',
     2'K*-  ','AD*0 ','D*-  ','F*-  ','JPSI ','F*+  ','D*+  ','D*0  ',
     3'UB*  ','DB*  ','SB*  ','CB*  ','UPSL ','BC*  ','BS*  ','BD*  ',
     4'BU*  ','UT*  ','DT*  ','ST*  ','CT*  ','BT*  ','TT*  ','TB*  ',
     5'TC*  ','TS*  ','TD*  ','TU*  ','UY*  ','DY*  ','SY*  ','CY*  ',
     6'BY*  ','TY*  ','YY*  ','YT*  ','YB*  ','YC*  ','YS*  ','YD*  ',
     7'YU*  ','UX*  ','DX*  ','SX*  ','CX*  ','BX*  ','TX*  ','YX*  ',
     8'XX*  ','XY*  ','XT*  ','XB*  ','XC*  ','XS*  ','XD*  ','XU*  '/
C          1/2+ BARYON LABELS
      DATA LBAR0/
     1'ERR  ','P    ','N    ','ERR  ','ERR  ','S+   ','S0   ','S-   ',
     2'L    ','XI0  ','XI-  ','ERR  ','ERR  ','ERR  ','SC++ ','SC+  ',
     3'SC0  ','LC+  ','USC. ','DSC. ','SSC. ','SDC. ','SUC. ','UCC. ',
     4'DCC. ','SCC. ','ERR  ','ERR  ','ERR  ','ERR  ','UUB. ','UDB. ',
     5'DDB. ','DUB. ','USB. ','DSB. ','SSB. ','SDB. ','SUB. ','UCB. ',
     6'DCB. ','SCB. ','CCB. ','CSB. ','CDB. ','CUB. ','UBB. ','DBB. ',
     7'SBB. ','CBB. ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','UTT. ',
     8'UDT. ','DDT. ','DUT. ','UST. ','DST. ','SST. ','SDT. ','SUT. ',
     9'UCT. ','DCT. ','SCT. ','CCT. ','CST. ','CDT. ','CUT. ','UBT. ',
     1'DBT. ','SBT. ','CBT. ','BBT. ','BCT. ','BST. ','BDT. ','BUT. ',
     2'UTT. ','DTT. ','STT. ','CTT. ','BTT. ','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY. ','UDY. ','DDY. ','DUY. ','USY. ',
     4'DSY. ','SSY. ','SDY. ','SUY. ','UUX. ','UDX. ','DDX. ','DUX. ',
     5'USX. ','DSX. ','SSX. ','SDX. ','SUX. '/
      DATA LABAR0/
     1'ERR  ','AP   ','AN   ','ERR  ','ERR  ','AS-  ','AS0  ','AS+  ',
     2'AL   ','AXI0 ','AXI+ ','ERR  ','ERR  ','ERR  ','ASC--','ASC- ',
     3'ASC0 ','ALC- ','AUSC.','ADSC.','ASSC.','ASDC.','ASUC.','AUCC.',
     4'ADCC.','ASCC.','ERR  ','ERR  ','ERR  ','ERR  ','AUUB.','AUDB.',
     5'ADDB.','ADUB.','AUSB.','ADSB.','ASSB.','ASDB.','ASUB.','AUCB.',
     6'ADCB.','ASCB.','ACCB.','ACSB.','ACDB.','ACUB.','AUBB.','ADBB.',
     7'ASBB.','ACBB.','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','AUTT.',
     8'AUDT.','ADDT.','ADUT.','AUST.','ADST.','ASST.','ASDT.','ASUT.',
     9'AUCT.','ADCT.','ASCT.','ACCT.','ACST.','ACDT.','ACUT.','AUBT.',
     1'ADBT.','ASBT.','ACBT.','ABBT.','ABCT.','ABST.','ABDT.','ABUT.',
     2'AUTT.','ADTT.','ASTT.','ACTT.','ABTT.','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY.','AUDY.','ADDY.','ADUY.','AUSY.',
     4'ADSY.','ASSY.','ASDY.','ASUY.','AUUX.','AUDX.','ADDX.','ADUX.',
     5'AUSX.','ADSX.','ASSX.','ASDX.','ASUX.'/
C          3/2+ BARYON LABELS
      DATA LBAR1/
     1'DL++ ','DL+  ','DL0  ','DL-  ','ERR  ','S*+  ','S*0  ','S*-  ',
     2'ERR  ','XI*0 ','XI*- ','OM-  ','ERR  ','ERR  ','UUC* ','UDC* ',
     3'DDC* ','ERR  ','USC* ','DSC* ','SSC* ','ERR  ','ERR  ','UCC* ',
     4'DCC* ','SCC* ','CCC* ','ERR  ','ERR  ','ERR  ','UUB* ','UDB* ',
     5'DDB* ','ERR  ','USB* ','DSB* ','SSB* ','ERR  ','ERR  ','UCB* ',
     6'DCB* ','SCB* ','CCB* ','ERR  ','ERR  ','ERR  ','UBB* ','DBB* ',
     7'SBB* ','CBB* ','BBB* ','ERR  ','ERR  ','ERR  ','ERR  ','UTT* ',
     8'UDT* ','DDT* ','ERR  ','UST* ','DST* ','SST* ','ERR  ','ERR  ',
     9'UCT* ','DCT* ','SCT* ','CCT* ','ERR  ','ERR  ','ERR  ','UBT* ',
     1'DBT* ','SBT* ','CBT* ','BBT* ','ERR  ','ERR  ','ERR  ','ERR  ',
     2'UTT* ','DTT* ','STT* ','CTT* ','BTT* ','TTT* ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY* ','UDY* ','DDY* ','ERR  ','USY* ',
     4'DSY* ','SSY* ','ERR  ','ERR  ','UUX* ','UDX* ','DDX* ','ERR  ',
     5'USX* ','DSX* ','SSX* ','ERR  ','ERR  '/
      DATA LABAR1/
     1'ADL--','ADL- ','ADL0 ','ADL+ ','ERR  ','AS*- ','AS*0 ','AS*+ ',
     2'ERR  ','AXI*0','AXI*+','AOM+ ','ERR  ','ERR  ','AUUC*','AUDC*',
     3'ADDC*','ERR  ','AUSC*','ADSC*','ASSC*','ERR  ','ERR  ','AUCC*',
     4'ADCC*','ASCC*','ACCC*','ERR  ','ERR  ','ERR  ','AUUB*','AUDB*',
     5'ADDB*','ERR  ','AUSB*','ADSB*','ASSB*','ERR  ','ERR  ','AUCB*',
     6'ADCB*','ASCB*','ACCB*','ERR  ','ERR  ','ERR  ','AUBB*','ADBB*',
     7'ASBB*','ACBB*','ABBB*','ERR  ','ERR  ','ERR  ','ERR  ','AUTT*',
     8'AUDT*','ADDT*','ERR  ','AUST*','ADST*','ASST*','ERR  ','ERR  ',
     9'AUCT*','ADCT*','ASCT*','ACCT*','ERR  ','ERR  ','ERR  ','AUBT*',
     1'ADBT*','ASBT*','ACBT*','ABBT*','ERR  ','ERR  ','ERR  ','ERR  ',
     2'AUTT*','ADTT*','ASTT*','ACTT*','ABTT*','ATTT*','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY*','AUDY*','ADDY*','ERR  ','AUSY*',
     4'ADSY*','ASSY*','ERR  ','ERR  ','AUUX*','AUDX*','ADDX*','ERR  ',
     5'AUSX*','ADSX*','ASSX*','ERR  ','ERR  '/
C          ENTRY
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
C@@@@@@@@@@ SIVOCL @@@@@
      IF(ID.EQ.110.OR.ID.EQ.111.OR.ID.EQ.221.
     *OR.ID.EQ.220.OR.ID.EQ.330) THEN
      IDABS=IABS(ID)
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      IFL1=0
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,-ID)
      ENDIF
C@@@@@@@@@@ SIVOCL @@@@@
      IF(IABS(ID).LT.100) GO TO 200
      IF(IABS(ID).LT.1000) GO TO 100
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 300
C          BARYONS
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-11
      IF(JSPIN.EQ.0.AND.ID.GT.0) LABEL1=LBAR0(INDEX)
      IF(JSPIN.EQ.0.AND.ID.LT.0) LABEL1=LABAR0(INDEX)
      IF(JSPIN.EQ.1.AND.ID.GT.0) LABEL1=LBAR1(INDEX)
      IF(JSPIN.EQ.1.AND.ID.LT.0) LABEL1=LABAR1(INDEX)
      RETURN
C          MESONS
100   CONTINUE
      I=MAX0(IFL2,IFL3)
      J=-MIN0(IFL2,IFL3)
      INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)
      IF(JSPIN.EQ.0) LABEL1=LMES0(INDEX)
      IF(JSPIN.EQ.1) LABEL1=LMES1(INDEX)
      RETURN
C          QUARKS, LEPTONS, ETC.
200   CONTINUE
      INDEX=2*INDEX
      IF(ID.LE.0) INDEX=INDEX+1
      LABEL1=LLEP(INDEX)
      RETURN
300   I=IABS(IFL1)
      J=IABS(IFL2)
      INDEX=I+J*(J-1)/2
      IF(ID.GT.0) LABEL1=LQQ(INDEX)
      IF(ID.LT.0) LABEL1=LAQQ(INDEX)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION CHARGE(ID)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C          COMPUTE CHARGE OF PARTICLE WITH IDENT CODE ID
C          ICHRG MUST BE DIMENSIONED NQLEP+12
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      DIMENSION ICHRG(53),IFL(3)
      DATA ICHRG/0,2,-1,-1,2,-1,2,-1,2,0,0,0,-3,0,-3,0,-3,0,-3
     $,0,0,0,2,-1,-1,2,-1,2,-1,2,0,0,0,-3,0,-3,0,-3,0,-3,3,0
     $,3,0,0,0,3,3,3,6,6,6,0/
      IDABS=IABS(ID)
      CALL FLAVOR(ID,IFL(1),IFL(2),IFL(3),JSPIN,INDEX)
c      write(16,*) 'IFL(3)=', IFL
      if(IDABS >= 100) then
        ISUM=0
        DO I=1,3
          IS1=1
          if(IFL(I).ne.0)  IS1=ISIGN(1,IFL(I))
          ISUM=ISUM+ICHRG(IABS(IFL(I))+1)*IS1
c         ISUM=ISUM+ICHRG(IABS(IFL(I))+1)*ISIGN(1,IFL(I))
c          write(16,*) 'I,ISUM,IS1,IFL(I)=',I,ISUM,IS1,IFL(I)
           ISUM=ISUM    !!!! Lena
        ENDDO
      else
        ISUM=ICHRG(INDEX+1)*ISIGN(1,ID)
      endif
      CHARGE=ISUM/3.d0
c      write(16,*) 'ID,ISUM,CHARGE=', ID,ISUM,CHARGE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DECAYQ(IP,IPOINT)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      integer*4 :: NJSET,JORIG(100),JTYPE(100),
     &     JDCAY(100)
      real*8 :: pjset(5, 100)
C          THIS SUBROUTINE DECAYS PARTICLE IP FROM /PARTCL/ USING THE
C          BRANCHING RATIOS FROM /DKYTAB/ AND ADDS THE DECAY PRODUCTS
C          TO /PARTCL/ WITH IORIG=IP.
C          QUARK-BASED IDENT CODE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/COMLID/PLIDER(499)
      COMMON/PARORD/IORDP(499)
C     LOOK MUST BE DIMENSIONED TO THE MAXIMUM VALUE OF INDEX
      COMMON/DKYTAB/LOOK(400),CBR(600),MODE(5,600)
      LOGICAL NODCAY,NOETA,NOPI0,NOKA0
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NOKA0
      COMMON/JWORK/ZZC(100),P1CM(4),E1CM,E2CM,E3CM,E4CM,E5CM,
     &     J1,J2,J3,J4,J5,JMATCH(100),TNEW
      LOGICAL TNEW
      DIMENSION JJ(5),EE(5)
      EQUIVALENCE (J1,JJ(1)),(E1CM,EE(1))
      DIMENSION PGEN(5,5),RND(5),U(3),BETA(3),ROT(3,3),PSAVE(3),
     &     REDUCE(5), PSUM(5)
      DATA REDUCE/1.,1.,2.,5.,15./
      DATA PSUM/5*0./
      DATA TWOME/1.022006E-3/
C          FUNCTION DEFINITIONS
C          PCM AND DOT MUST BE CALCULATED IN DOUBLE PRECISION
C          ON 32-BIT MACHINES.
CC    PCM(A,B,C)=SQRT((A**2-B**2-C**2)**2-(2.*B*C)**2)/(2.*A)
C          SELECT DECAY MODE
 92   MXPTCL=499
      IPACK=1000
      IDPACK=10000
      MXJSET=100
      JPACK=1000
      IDLV1=IDENT(IP)
      CALL FLAVOR(IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX)
C     IF(IABS(IFL1).GT.4.OR.IABS(IFL2).GT.4.OR.IABS(IFL3).GT.4)GOTO101
      IF(NOPI0.AND.IDLV1.EQ.110) GO TO 101
      IF(NOETA.AND.IDLV1.EQ.220) GO TO 101
      IF(NOKA0.AND.(IDLV1.EQ.230.OR.IDLV1.EQ.-230)) GO TO 101
      GO TO 102
101   IPOINT=-1
      RETURN
102   CONTINUE
1     CONTINUE
      IPOINT=LOOK(INDEX)-1
      IF(IPOINT.LT.0) RETURN
      IF(IPRINT) WRITE(ITLIS,2000) IDENT(IP),(PPTCL(J,IP),J=1,4)
2000  FORMAT(' DECAY OF ',I6,' WITH MOM.',3F10.3,'& ENERGY',F10.3)
      TRY=RNDM(-1.)
100   IPOINT=IPOINT+1
      IF(TRY.GT.CBR(IPOINT)) GO TO 100
      NADD=0
      SUM=0.
      NSTART=NPTCL+1
      DO 110 I=1,5
      IF(MODE(I,IPOINT).EQ.0) GO TO 110
      IF(NPTCL+NADD+1.GT.MXPTCL) GO TO 9999
      NADD=NADD+1
      NEW=NPTCL+NADD
      IDENT(NEW)=MODE(I,IPOINT)
      IDLV1=IDENT(NEW)
      PPTCL(5,NEW)=AMASSF(IDLV1)
      SUM=SUM+PPTCL(5,NEW)
110   CONTINUE
       IF(SUM.GT.PPTCL(5,IP)) GO TO 1991
      NADD1=NADD-1
      DO 120 J=1,5
      PGEN(J,1)=PPTCL(J,IP)
120   CONTINUE
      PGEN(5,NADD)=PPTCL(5,NPTCL+NADD)
      IF(NADD.EQ.1) GO TO 700
      IF(NADD.EQ.2) GO TO 400
C          USE KROLL-WADA DISTRIBUTION FOR DALITZ DECAYS.
      IF(.NOT.((IDENT(IP).EQ.110.OR.IDENT(IP).EQ.220).AND.
     1IABS(IDENT(NPTCL+2)).EQ.12)) GO TO 130
125   AMEE=TWOME*(PPTCL(5,IP)/TWOME)**RNDM(-1.)
      REE=(TWOME/AMEE)**2
      WTEE=(1.-(AMEE/PPTCL(5,IP))**2)**3*SQRT(1.-REE)*(1.+.5*REE)
      DRND=RNDM(-1.)
      IF(WTEE.LT.DRND) GO TO 125
      PGEN(5,2)=AMEE
      GO TO 400
130   CONTINUE
C          CALCULATE MAXIMUM PHASE-SPACE WEIGHT
      WTMAX=1./REDUCE(NADD)
      SUM1=PGEN(5,1)
      SUM2=SUM-PPTCL(5,NPTCL+1)
      DO 200 I=1,NADD1
      WTMAX=WTMAX*DBLPCM(SUM1,SUM2,PPTCL(5,NPTCL+I))
      SUM1=SUM1-PPTCL(5,NPTCL+I)
      SUM2=SUM2-PPTCL(5,NPTCL+I+1)
200   CONTINUE
C          GENERATE UNIFORM NADD-BODY PHASE SPACE
300   CONTINUE
      RND(1)=1.
      DO 310 I=2,NADD1
      RNEW=RNDM(-1.)
      I1=I-1
      DO 320 JJ1=1,I1
      J=I-JJ1
      JSAVE=J+1
      IF(RNEW.LE.RND(J)) GO TO 310
      RND(JSAVE)=RND(J)
320   CONTINUE
310   RND(JSAVE)=RNEW
      RND(NADD)=0.
      WT=1.
      SUM1=SUM
      DO 330 I=2,NADD
      SUM1=SUM1-PPTCL(5,NPTCL+I-1)
      PGEN(5,I)=SUM1+RND(I)*(PGEN(5,1)-SUM)
      WT=WT*DBLPCM(PGEN(5,I-1),PGEN(5,I),PPTCL(5,NPTCL+I-1))
330   CONTINUE
      IF(WT.LT.RNDM(-1.)*WTMAX) GO TO 300
C          CARRY OUT TWO-BODY DECAYS IN PGEN FRAMES
400   CONTINUE
      DO 410 I=1,NADD1
      QCM=DBLPCM(PGEN(5,I),PGEN(5,I+1),PPTCL(5,NPTCL+I))
      U(3)=2.*RNDM(-1.)-1.
      PHI=twpi*RNDM(-1.)
      U(1)=SQRT(1.-U(3)**2)*COS(PHI)
      U(2)=SQRT(1.-U(3)**2)*SIN(PHI)
      DO 420 J=1,3
      PPTCL(J,NPTCL+I)=QCM*U(J)
      PGEN(J,I+1)=-PPTCL(J,NPTCL+I)
420   CONTINUE
      PPTCL(4,NPTCL+I)=SQRT(QCM**2+PPTCL(5,NPTCL+I)**2)
      PGEN(4,I+1)=SQRT(QCM**2+PGEN(5,I+1)**2)
410   CONTINUE
      DO 430 J=1,4
      PPTCL(J,NPTCL+NADD)=PGEN(J,NADD)
430   CONTINUE
C          BOOST PGEN FRAMES TO LAB FRAME
      DO 500 II=1,NADD1
      I=NADD-II
      DO 510 J=1,3
      BETA(J)=PGEN(J,I)/PGEN(4,I)
510   CONTINUE
      GAMMA=PGEN(4,I)/PGEN(5,I)
      DO 520 K=I,NADD
      K1=NPTCL+K
      BP=BETA(1)*PPTCL(1,K1)+BETA(2)*PPTCL(2,K1)+BETA(3)*PPTCL(3,K1)
      DO 530 J=1,3
      PPTCL(J,K1)=PPTCL(J,K1)+GAMMA*BETA(J)*(PPTCL(4,K1)
     1+BP*GAMMA/(GAMMA+1.))
530   CONTINUE
      PPTCL(4,K1)=GAMMA*(PPTCL(4,K1)+BP)
520   CONTINUE
500   CONTINUE
C          MATRIX ELEMENTS
      IF(NADD.EQ.3.AND.(IDENT(IP).EQ.221.OR.IDENT(IP).EQ.331)) GO TO 610
      IF(NADD.EQ.3.AND.IABS(IDENT(NPTCL+1)).LT.20.AND.
     1IDENT(NPTCL+1).NE.10) GO TO 620
      GO TO 800
C          OMEG AND PHI DECAY
610   WT=(PPTCL(5,NPTCL+1)*PPTCL(5,NPTCL+2)*PPTCL(5,NPTCL+3))**2
     1-(PPTCL(5,NPTCL+1)*DOT(NPTCL+2,NPTCL+3))**2
     2-(PPTCL(5,NPTCL+2)*DOT(NPTCL+1,NPTCL+3))**2
     3-(PPTCL(5,NPTCL+3)*DOT(NPTCL+1,NPTCL+2))**2
     4+2.*DOT(NPTCL+1,NPTCL+2)*DOT(NPTCL+2,NPTCL+3)*DOT(NPTCL+1,NPTCL+3)
      IF(WT.LT.RNDM(-1.)*PPTCL(5,IP)**6/108.) GO TO 300
      GO TO 800
C          SEMILEPTONIC AND QUARK DECAYS
C          INCLUDE W PROPAGATOR
620   WT=DOT(IP,NPTCL+2)*DOT(NPTCL+1,NPTCL+3)
      S12=PPTCL(5,NPTCL+1)**2+PPTCL(5,NPTCL+2)**2
     1+2.*DOT(NPTCL+1,NPTCL+2)
      S12MAX=PPTCL(5,IP)**2
      WT=WT*WPROP(S12MAX)/WPROP(S12)
      IF(WT.LT.RNDM(-1.)*PPTCL(5,IP)**4/16.) GO TO 300
      GO TO 800
C          ONE-PARTICLE DECAYS
700   CONTINUE
      DO 710 J=1,5
      PPTCL(J,NPTCL+1)=PPTCL(J,IP)
710   CONTINUE
C          SWAP PARTICLES AND ANTIPARTICLES IF IDENT(IP)<0
800   CONTINUE
      IF(IDENT(IP).GE.0.OR.IABS(IDENT(IP)).EQ.20) GO TO 900
      DO 810 I=1,NADD
      IDABS=IABS(IDENT(NPTCL+I))
      IFL1=IDABS/1000
      IFL2=MOD(IDABS/100,10)
      IFL3=MOD(IDABS/10,10)
      IF(IFL1.EQ.0.AND.IFL2.NE.0.AND.IFL2.EQ.IFL3) GO TO 810
      IF(IDABS.EQ.9.OR.IDABS.EQ.10.OR.IDABS.EQ.20) GO TO 810
      IF(IDABS.EQ.29.OR.IDABS.EQ.30.OR.IDABS.EQ.40) GO TO 810
      IDENT(NPTCL+I)=-IDENT(NPTCL+I)
810   CONTINUE
C          REMOVE QUARKS FROM /PARTCL/ AND TRANSFORM BACK TO REST FRAME
900   CONTINUE
      NPTCL=NPTCL+NADD
      NQK=0
      IF(IABS(IDENT(NPTCL)).GE.10.AND.MOD(IDENT(NPTCL),100).NE.0)
     1GO TO 1000
      NOFF=NPTCL-NSTART+1
      DO 910 II=1,NOFF
      I=NPTCL+1-II
      IF(IABS(IDENT(NPTCL)).LT.10.OR.MOD(IDENT(NPTCL),100).NE.0)
     1NQK=NQK+1
      BP=BETA(1)*PPTCL(1,I)+BETA(2)*PPTCL(2,I)+BETA(3)*PPTCL(3,I)
      DO 911 J=1,3
      PPTCL(J,I)=PPTCL(J,I)-GAMMA*BETA(J)*(PPTCL(4,I)
     1-BP*GAMMA/(GAMMA+1.))
911   CONTINUE
      PPTCL(4,I)=GAMMA*(PPTCL(4,I)-BP)
910   CONTINUE
C          COPY DECAY PRODUCTS INTO /JETSET/
      IF(NJSET+NADD.GT.MXJSET) GO TO 9998
      NJSAVE=NJSET
      NPTCL=NPTCL-NADD
      DO 920 I=1,NADD
      NJSET=NJSET+1
      DO 921 K=1,5
921   PJSET(K,NJSET)=PPTCL(K,NPTCL+I)
      JORIG(NJSET)=0
      JTYPE(NJSET)=IDENT(NPTCL+I)
      JDCAY(NJSET)=0
      JMATCH(NJSET)=JPACK*(NJSAVE+1)+NJSAVE+NADD
C          QCD EVOLUTION STARTS FROM PARENT MASS
C          BUT USE NADD*ENERGY TO PRESERVE TP --> W+ BT
      IF(IABS(JTYPE(NJSET)).GE.10.AND.MOD(JTYPE(NJSET),100).NE.0)
     1GO TO 920
      JDCAY(NJSET)=-1
      PJSET(5,NJSET)=DMIN1(PPTCL(5,IP),NADD*PJSET(4,NJSET))
920   CONTINUE
C          PERFORM QCD JET EVOLUTION
      PRINT*,'DECAY, IDENT(IP)= ',IDENT(IP)
c     CALL QCDFRG(NJSAVE+1)
C          DECAY QUARKS AND ROTATE TO PROPER ANGLES
C          HADRONIZE JETS
      NJ1=NJSAVE+1
      DO 931 I=NJ1,NJSET
      IF(JDCAY(I).NE.0) GO TO 931
      IF(IABS(JTYPE(I)).GE.10.AND.MOD(JTYPE(I),100).NE.0)
     1GO TO 935
      NEXT=NPTCL+1
      PJET=SQRT(PJSET(1,I)**2+PJSET(2,I)**2+PJSET(3,I)**2)
      CTHQK=PJSET(3,I)/PJET
      STHQK=SQRT(1.-CTHQK**2)
      CPHIQK=PJSET(1,I)/(PJET*STHQK)
      SPHIQK=PJSET(2,I)/(PJET*STHQK)
c      CALL JETGEN(I)
      IF(NEXT.GT.NPTCL) GO TO 931
      ROT(1,1)=CPHIQK*CTHQK
      ROT(2,1)=SPHIQK*CTHQK
      ROT(3,1)=-STHQK
      ROT(1,2)=-SPHIQK
      ROT(2,2)=CPHIQK
      ROT(3,2)=0.
      ROT(1,3)=CPHIQK*STHQK
      ROT(2,3)=SPHIQK*STHQK
      ROT(3,3)=CTHQK
      DO 932 K=NEXT,NPTCL
      DO 933 J=1,3
      PSAVE(J)=PPTCL(J,K)
      PPTCL(J,K)=0.
933   CONTINUE
      DO 932 J=1,3
      DO 932 JJ1=1,3
      PPTCL(J,K)=PPTCL(J,K)+ROT(J,JJ1)*PSAVE(JJ1)
932   CONTINUE
      GOTO 931
C          ADD LEPTON TO /PARTCL/
  935 NPTCL=NPTCL+1
      DO 936 K=1,5
  936 PPTCL(K,NPTCL)=PJSET(K,I)
      IDENT(NPTCL)=JTYPE(I)
931   CONTINUE
C          RESET NJSET SO DECAY JETS DO NOT APPEAR IN /JETSET/
      NJSET=NJSAVE
C          CHECK FOR AT LEAST TWO PARTICLES
      IF(NPTCL.GT.NSTART) GO TO 939
      NPTCL=NSTART-1
      GO TO 1
939   CONTINUE
C          CONSERVE CHARGE
      SUMQ=0.
      DO 960 I=NSTART,NPTCL
      IDLV1=IDENT(I)
      SUMQ=SUMQ+CHARGE(IDLV1)
960   CONTINUE
      IDLV1=IDENT(IP)
      SUMQ=SUMQ-CHARGE(IDLV1)
      IF(SUMQ.EQ.0.) GO TO 970
      DO 961 I=NSTART,NPTCL
      ID1=IDENT(I)
      IF(IABS(ID1).GT.1000) GO TO 961
      I1=MOD(IABS(ID1)/100,10)
      I2=MOD(IABS(ID1)/10,10)
      I3=MOD(IABS(ID1),10)
      IF(I1.EQ.1.AND.I2.GT.2.AND.SUMQ*ID1.GT.0.) GO TO 962
      IF(I1.EQ.2.AND.I2.GT.2.AND.SUMQ*ID1.LT.0.) GO TO 963
      IF(I1.EQ.1.AND.I2.EQ.2.AND.SUMQ*ID1.GT.0.) GO TO 964
      IF(I1.EQ.1.AND.I2.EQ.1) GO TO 965
      GO TO 961
962   IDENT(I)=ISIGN(200+10*I2+I3,ID1)
      GO TO 969
963   IDENT(I)=ISIGN(100+10*I2+I3,ID1)
      GO TO 969
964   IDENT(I)=110+I3
      GO TO 969
965   IDENT(I)=(120+I3)*(-SIGN(1.D0,SUMQ))
969   SUMQ=SIGN(ABS(SUMQ)-1.,SUMQ)
      IDLV1=IDENT(I)
      PPTCL(5,I)=AMASSF(IDLV1)
      PPTCL(4,I)=SQRT(PPTCL(1,I)**2+PPTCL(2,I)**2+PPTCL(3,I)**2
     1+PPTCL(5,I)**2)
      IF(SUMQ.EQ.0.) GO TO 970
961   CONTINUE
      NPTCL=NSTART-1
      GO TO 1
C          RESCALE MOMENTA FOR CORRECT MASS
970   CONTINUE
      PSUM(4)=PPTCL(5,IP)
      PSUM(5)=PSUM(4)
      NPTLV1=NPTCL
      CALL RESCAL(NSTART,NPTLV1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 940
      NPTCL=NSTART-1
      GO TO 1
940   CONTINUE
C          BOOST BACK TO LAB FRAME
      DO 950 I=NSTART,NPTCL
      BP=BETA(1)*PPTCL(1,I)+BETA(2)*PPTCL(2,I)+BETA(3)*PPTCL(3,I)
      DO 951 J=1,3
      PPTCL(J,I)=PPTCL(J,I)+GAMMA*BETA(J)*(PPTCL(4,I)
     1+BP*GAMMA/(GAMMA+1.))
951   CONTINUE
      PPTCL(4,I)=GAMMA*(PPTCL(4,I)+BP)
950   CONTINUE
C          SET IORIG AND IDCAY
1000  CONTINUE
      IDCAY(IP)=IDPACK*NSTART+NPTCL
      IF(IPRINT) WRITE(ITLIS,2001)
2001  FORMAT(' PRODUCTS OF DECAY:   ID     .    PX    .  PY    ',
     *'    PZ    .   E       ')
      DO 1010 I=NSTART,NPTCL
      IF(IPRINT) WRITE(ITLIS,2002) IDENT(I),(PPTCL(J,I),J=1,4)
2002  FORMAT(18X,I6,4F10.3)
C      IORIG(I)=IP + IORIG(IP)*IPACK
       IORIG(I)=IP
      IDCAY(I)=0
      PLIDER(I)=PLIDER(IP)
      IORDP(I)=IORDP(IP)
      PPTCL(6,I)=PPTCL(6,IP)
      PPTCL(7,I)=PPTCL(7,IP)
      PPTCL(8,I)=PPTCL(8,IP)
      PPTCL(9,I)=PPTCL(9,IP)
1010  CONTINUE
      RETURN
C
9999  WRITE(ITLIS,10) NPTCL
10    FORMAT(//5X,25HERROR IN DECAY...NPTCL > ,I5)
      RETURN
C
1991  CONTINUE                          ! no print: 02.04.01
c      WRITE(ITLIS,1992)  SUM,PPTCL(5,IP),IDENT(IP),NADD
c      WRITE(17,1992)  SUM,PPTCL(5,IP),IDENT(IP),NADD
1992  FORMAT(1X,'ERROR IN DECAY:SUMM>MD ',2(F7.3,1X),2I5)
      PPTCL(5,IP)=SUM
      GO  TO  92
C
9998  WRITE(ITLIS,20) NJSET
20    FORMAT(//5X,25HERROR IN DECAY...NJSET > ,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION DOT(I1,I2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      DOT=PPTCL(4,I1)*PPTCL(4,I2)-PPTCL(1,I1)*PPTCL(1,I2)
     *   -PPTCL(2,I1)*PPTCL(2,I2)-PPTCL(3,I1)*PPTCL(3,I2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     SUBROUTINE  JETGEN(I)
c     K=I
c     RETURN
c     END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION WPROP(Z)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP
C          CHARGED W PROPAGATOR.
      WPROP=(Z-WMASS(2)**2)**2+(WMASS(2)*WGAM(2))**2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     SUBROUTINE  QCDFRG(I)
c     K=I
c     RETURN
c     END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IB(ID)
      IMPLICIT  INTEGER (I-N)
C
C    COMPUTE BARYON NUMBER
C
      IB=0
      IF(IABS(ID)/1000.NE.0) IB=ISIGN(1,ID)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IS(ID)
      IMPLICIT  INTEGER (I-N)
C
C   COMPUTE STRANGENESS
C
      CALL FLAVOR(ID,IF1,IF2,IF3,JSPIN,INDEX)
      IS=0
      IF(IABS(IF1).EQ.3) IS=IS+ISIGN(1,-IF1)
      IF(IABS(IF2).EQ.3) IS=IS+ISIGN(1,-IF2)
      IF(IABS(IF3).EQ.3) IS=IS+ISIGN(1,-IF3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CLUSTR(IFL1,IFL2,AMCTR)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  HADRONS PRODUCTION BY MEANS CLUSTER BREAKING
C  WITH QUARK AND ANTIQUARK OR QUARK AND DIQUARK OR DIQUARK AND
C  ANTIDIQUARK IFL1 AND IFL2 ON ENDS
C  AMCTR IS MASS OF CLUSTER
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/ITAPES/ ITVKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/FRGCPA/ PUDC,PUDCC,PSPINC,PJSPNC,PMIX1C(3,2),PMIX2C(3,2),
     *PBARC
      COMMON/COMNPT/ NOPTC
      LOGICAL NOPTC
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION IFL(2),U(3)
      LOGICAL SPINT
      COMMON/COMWTI/ WTIME
      LOGICAL WTIME
      MXPTCL=499
C     IF(.NOT.WTIME) GO TO 99
C     NOPTC=.FALSE.
C     CALL CLUSTT(IFL1,IFL2,AMCTR)
C     RETURN
99    NFIX=NPTCL
      NREP=0
      LRET=.FALSE.
 100  I=NFIX
      IF(NREP.LT.NTRIES) GO TO 101
      LRET=.TRUE.
      IF(IPRINT) WRITE(ITLIS,1200) NREP,IFL1,IFL2,AMCTR
1200  FORMAT(3X,' IN CLUSTR NREP GT ',3I8,' AMCTR=',F12.4)
      RETURN
101   CONTINUE
      KSPIN=0
      IFL(1)=IFL1
      IFL(2)=IFL2
      SPINT=.FALSE.
      I=I+2
      IF(I.GT.MXPTCL) GO TO 9999
C  CHOOSE SIDE OF BREAK
      JSIDE=1
C  IF ANY IFL IS A DIQUARK
      IF(MOD(IFL(1),100).EQ.0.OR.MOD(IFL(2),100).EQ.0) GO TO 150
C  IFL(1) AND IFL(2) ARE QUARKS
C  SELECT Q,QBARPAIR OR QQ,QQBAR PAIR
      DRND=RNDM(-1.)
      IF(DRND.LT.PBARC) GO TO 140
C  Q,QBAR PAIR
      IFLN=ISIGN(INT(RNDM(-1.)/PUDC)+1,-IFL(JSIDE))
      IF(IABS(IFLN).EQ.3)
     *IFLN=ISIGN(INT(RNDM(-1.)/PUDCC)+3,-IFL(JSIDE))
      GO TO 200
C  QQ,QQBAR PAIR
140   IQ1=INT(RNDM(-1.)/PUDC)+1
      IF(IQ1.EQ.3) IQ1=INT(RNDM(-1.)/PUDCC)+3
      IQ2=INT(RNDM(-1.)/PUDC)+1
      IF(IQ2.EQ.3) IQ2=INT(RNDM(-1.)/PUDCC)+3
      IF(IQ1.LE.IQ2) GO TO 145
      ISWAP=IQ1
      IQ1=IQ2
      IQ2=ISWAP
145   IFQQ=1000*IQ1+100*IQ2
      IFLN=ISIGN(IFQQ,IFL(JSIDE))
      GO TO 200
C  IFL(1) OR IFL(2) IS DIQUARK
C Q,QBAR PAIR
150   IPSIGN=IFL(JSIDE)
      IF(MOD(IFL(JSIDE),100).EQ.0) GO TO 130
      IPSIGN=-IFL(JSIDE)
130   IFLN=ISIGN(INT(RNDM(-1.)/PUDC)+1,IPSIGN)
      IF(IABS(IFLN).EQ.3)
     *IFLN=ISIGN(INT(RNDM(-1.)/PUDCC)+3,IPSIGN)
C  IDENTS AND MASSES OF PARTICLES
 200  IDENT(I-1)=IDPARC(IFL(JSIDE),IFLN,SPINT,KSPIN)
      IDENT(I)=IDPARC(IFL(3-JSIDE),-IFLN,SPINT,KSPIN)
      PPTCL(5,I-1)=AMASS(IDENT(I-1))
      PPTCL(5,I)=AMASS(IDENT(I))
C  IF TOO LOW MASS,START ALL OVER
      IF(AMCTR.GT.PPTCL(5,I-1)+PPTCL(5,I)) GO TO 102
      NREP=NREP+1
      GO TO 100
102   CONTINUE
      PA=DBLPCM(AMCTR,PPTCL(5,I-1),PPTCL(5,I))
C      PROB=2.*PA/AMCTR
C   PROB IS TWO-BODY PHASE SPACE FACTOR
C      DRND=RNDM(-1.)
C      IF(DRND.GT.PROB) GO TO 100
      U(3)=COSDD(0)
      PHI=twpi*RNDM(-1.)
      ST=SQRT(1.-U(3)**2)
      U(1)=ST*COS(PHI)
      U(2)=ST*SIN(PHI)
      PPTCL(1,I-1)=PA*U(1)
      PPTCL(1,I)=-PA*U(1)
      PPTCL(2,I-1)=PA*U(2)
      PPTCL(2,I)=-PA*U(2)
      PPTCL(3,I-1)=PA*U(3)
      PPTCL(3,I)=-PA*U(3)
      PA2=PA**2
      PPTCL(4,I-1)=SQRT(PA2+PPTCL(5,I-1)**2)
      PPTCL(4,I)=SQRT(PA2+PPTCL(5,I)**2)
      IDCAY(I-1)=0
      IDCAY(I)=0
      NPTCL=I
      RETURN
9999  WRITE(ITLIS,9998) I
9998  FORMAT(//10X,40H...STOP IN CLUSTR..NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IDPARC(IFL01,IFL02,SPINT,IR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C   CONSTRUCT MESON FROM QUARK AND ANTIQUARK WITH FLAVORS IFL01,IFL02
C   OR CONSTRUCT BARION FROM DIQUARK AND QUARK OR ANTIDIQUARK AND
C  ANTIQUARK WITH FLAVORS IFL01,IFL02
      COMMON/FRGCPA/ PUDC,PUDCC,PSPINC,PJSPNC,PMIX1C(3,2),PMIX2C(3,2),
     *PBARC
      LOGICAL SPINT
      IFL1=IFL01
      IFL2=IFL02
C  CONSTRUCT MESON WITH ACCOUNT FLAVOR MIXING
      IF(MOD(IFL1,100).EQ.0) GO TO 420
      IF(MOD(IFL2,100).EQ.0) GO TO 425
      JSPIN=INT(PSPINC+RNDM(-1.))
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
      ID1=IFL1
      ID2=IFL2
      IF(ID1+ID2.NE.0) GO TO 400
      RND=RNDM(-1.)
      ID1=IABS(ID1)
      ID1=INT(PMIX1C(ID1,JSPIN+1)+RND)+
     +INT(PMIX2C(ID1,JSPIN+1)+RND)+1
      ID2=-ID1
 400  IF(IABS(ID1).LE.IABS(ID2)) GO TO 410
      ISAVE=ID1
      ID1=ID2
      ID2=ISAVE
 410  IDHAD=ISIGN(100*IABS(ID1)+10*IABS(ID2)+JSPIN,ID1)
      GO TO 470
C  CONSTRUCT BARYON IDENT
 420  ID3=ISIGN(MOD(IFL1/100,10),IFL1)
      ID2=IFL1/1000
      ID1=IFL2
      GO TO 430
 425  ID3=ISIGN(MOD(IFL2/100,10),IFL2)
      ID2=IFL2/1000
      ID1=IFL1
 430  IF(IABS(ID1).LE.IABS(ID2)) GO TO 431
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 431  IF(IABS(ID2).LE.IABS(ID3)) GO TO 432
      ISWAP=ID2
      ID2=ID3
      ID3=ISWAP
 432  IF(IABS(ID1).LE.IABS(ID2)) GO TO 440
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 440  JSPIN=1
      IF(ID1.EQ.ID2.AND.ID2.EQ.ID3) GO TO 450
      JSPIN=INT(RNDM(-1.)+PJSPNC)
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
 450  IF(JSPIN.EQ.1.OR.ID1.EQ.ID2.OR.ID2.EQ.ID3) GO TO 460
      DRND=RNDM(-1.)
      IF(DRND.GT.PJSPNC) GO TO 460
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 460  IDHAD=1000*IABS(ID1)+100*IABS(ID2)+10*IABS(ID3)+JSPIN
      IDHAD=ISIGN(IDHAD,IFL1)
 470  IDPARC=IDHAD
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE STRING(IFL1,IFL2,AMSTR)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  HADRONS PRODUCTION BY MEANS STRING BREAKING
C  WITH QUARK AND ANTIQUARK OR QUARK AND DIQUARK OR DIQUARK AND
C  ANTIDIQUARK IFL1 AND IFL2 ON ENDS
C  AMSTR IS MASS OF STRING
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/COMLID/ PLIDER(499)
      COMMON/COMQSE/ QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      COMMON/FRGSPA/ PUDS,PUDSC,PSPINS,PJSPNS,PMIX1S(3,2),PMIX2S(3,2),
     *SIGQTS,WENDM,WENDB,PBARS,PRIQS(9),PARDBS,PARQLS,PARRS,PUNDS
      COMMON /PARCUT/ SWMAX
      DIMENSION PX1L(2),PY1L(2),PX1(2),PY1(2),PMTS(2),W(2),IFL(2)
      LOGICAL DIQBR,SPINT,BACK
      LOGICAL MESON
      DIMENSION VC(3)
      DIMENSION PPTL(5,200),PPTR(5,200),IDENTL(200),IDENTR(200)
      COMMON/KAPPA/ XAP
      COMMON /CINSID/ INSIDE
      COMMON/COMWTI/ WTIME
      LOGICAL WTIME
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON/COLRET/ LRET
      LOGICAL LRET
      DIMENSION P8(500),P9(500)
      LRET=.FALSE.
      NREP = 0
      MXPTCL=499
C     IF(.NOT.WTIME) GO TO 99
C     CALL STRINT(IFL1,IFL2,AMSTR)
C     RETURN
99    DIQBR=.TRUE.
      NFIX=NPTCL
      BACK=.TRUE.
 100  I=NFIX
      NPR=0
      NPL=0
      NPTCL=NFIX
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 102
      LRET=.TRUE.
      IF(IPRINT) WRITE(ITLIS,1200) NREP
1200  FORMAT(3X,' IN STRING NREP GT ',I8)
      RETURN
102   CONTINUE
      IFL(1)=IFL1
      IFL(2)=IFL2
      DO 110 J=1,2
 110  W(J)=AMSTR
      DO 120 J=1,2
      PX1L(J)=0.
      PY1L(J)=0.
      PX1(J)=0.
 120  PY1(J)=0.
C  WILL BE ONLY ONE BREAK OR NOT
      SPINT=.TRUE.
      KSPIN=1
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GO TO 131
      IDR=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      PARC1=0.35
      IF(IABS(IFL(1)).EQ.3.OR.IABS(IFL(2)).EQ.3) PARC1=0.5
      IF(IABS(IFL(1)).GE.4.AND.IABS(IFL(1)).LT.6) PARC1=0.7
      IF(IABS(IFL(2)).GE.4.AND.IABS(IFL(2)).LT.6) PARC1=0.7
      WEND=(AMASS(IDR)+PARC1)**2
      GO TO 151
131   IFCN=1
      IF(RNDM(-1.).GT.0.5) IFCN=2
      IFLC1=IFCN
      IF(IFL(1).LT.0) IFLC1=-IFCN
      IDR1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IDR2=IDPARS(IFL(2),-IFLC1,SPINT,KSPIN)
      WEND=(AMASS(IDR1)+AMASS(IDR2)+SWMAX)**2
151   SPINT=.FALSE.
      KSPIN=0
      IF(W(1)*W(2).LE.WEND) GO TO 225
 130  I=I+1
      IF(I.GT.MXPTCL) GO TO 9999
      MESON=.FALSE.
C  CHOOSE SIDE OF BREAK
      JSIDE=INT(1.+2.*RNDM(-1.))
      IF(JSIDE.EQ.1) NPR=NPR+1
      IF(JSIDE.EQ.2) NPL=NPL+1
      IF(NPR.GT.200.OR.NPL.GT.200) GO TO 9999
C  IF IFL(JSIDE) A QUARK OR A DIQUARK
      IF(MOD(IFL(JSIDE),100).EQ.0) GO TO 150
C  IFL(JSIDE) IS A QUARK
C SELECT Q,QBAR PAIR OR QQ,QQBAR PAIR
      DIQBR=.FALSE.
      DRND=RNDM(-1.)
      IF(DRND.LT.PBARS) GO TO 140
C  Q,QBAR PAIR
      IFLN=ISIGN(INT(RNDM(-1.)/PUDS)+1,-IFL(JSIDE))
      IF(IABS(IFLN).EQ.3)
     *IFLN=ISIGN(INT(RNDM(-1.)/PUDSC)+3,-IFL(JSIDE))
      GO TO 200
C  QQ,QQBAR PAIR
140   IQ1=INT(RNDM(-1.)/PUDS)+1
      IF(IQ1.EQ.3) IQ1=INT(RNDM(-1.)/PUDSC)+3
      IQ2=INT(RNDM(-1.)/PUDS)+1
      IF(IQ2.EQ.3) IQ2=INT(RNDM(-1.)/PUDSC)+3
      IF(IQ1.LE.IQ2) GO TO 145
      ISWAP=IQ1
      IQ1=IQ2
      IQ2=ISWAP
145   IFQQ=1000*IQ1+100*IQ2
      IFLN=ISIGN(IFQQ,IFL(JSIDE))
      GO TO 200
C  IFL(JSIDE) IS A DIQUARK
C  CAN DIQUARK BREAK OR NOT
 150  DRND=RNDM(-1.)
      IF(DRND.LE. PARDBS) GO TO 190
C DIQUARK BREAK
      MESON=.TRUE.
      CALL FLAVOR(IFL(JSIDE),IFLD1,IFLD2,IFLD3,JSPIN,INDEX)
      IFLL=IFLD1
      IFL(JSIDE)=IFLD2
      DRND=RNDM(-1.)
      IF(DRND.GE.PARQLS) GO TO 160
      IFLL=IFLD2
      IFL(JSIDE)=IFLD1
 160  DIQBR=.TRUE.
C  LEADING QUARK TRANSFERSE MOMENTUM
C     CALL PTDGET(PXL,PYL,SIGQTS)
      CALL GAUSPT(PTL0,SIGQTS)
      PHI=twpi*RNDM(-1.)
      PXL=PTL0*COS(PHI)
      PYL=PTL0*SIN(PHI)
      PX1L(JSIDE)=PX1(JSIDE)
      PY1L(JSIDE)=PY1(JSIDE)
      PX1(JSIDE)=-PXL
      PY1(JSIDE)=-PYL
C Q,QBAR PAIR
      IFLN=ISIGN(INT(RNDM(-1.)/PUDS)+1,-IFL(JSIDE))
      IF(IABS(IFLN).EQ.3)
     *IFLN=ISIGN(INT(RNDM(-1.)/PUDSC)+3,-IFL(JSIDE))
      GO TO 200
C DIQUARK DOES NOT BREAK
C Q,QBAR PAIR
 190  IFLN=ISIGN(INT(RNDM(-1.)/PUDS)+1,IFL(JSIDE))
      IF(IABS(IFLN).EQ.3)
     *IFLN=ISIGN(INT(RNDM(-1.)/PUDSC)+3,IFL(JSIDE))
      DIQBR=.FALSE.
C  IDENT,MASS AND TRANSFERSE MOMENTUM OF PARTICLE
 200  IDENT(I)=IDPARS(IFL(JSIDE),IFLN,SPINT,KSPIN)
      PPTCL(5,I)=AMASS(IDENT(I))
C     CALL PTDGET(PX2,PY2,SIGQTS)
      CALL GAUSPT(PT2,SIGQTS)
      PHI=twpi*RNDM(-1.)
      PX2=PT2*COS(PHI)
      PY2=PT2*SIN(PHI)
      PPTCL(1,I)=PX1(JSIDE)+PX2
      PPTCL(2,I)=PY1(JSIDE)+PY2
C GENERATE Z
      PMTS(3-JSIDE)=AMASS(IABS(IFL(3-JSIDE)))**2
      PTS=PPTCL(1,I)**2+PPTCL(2,I)**2
      PMTS(JSIDE)=PPTCL(5,I)**2+PTS
      IF(PMTS(JSIDE)+PMTS(3-JSIDE).GE.PARRS*W(1)*W(2)) GO TO 100
      ZMIN=PMTS(JSIDE)/(W(1)*W(2))
      ZMAX=1.-PMTS(3-JSIDE)/(W(1)*W(2))
      IF(ZMIN.GE.ZMAX) GO TO 100
C  WARNING. VERY IMPORTANT THE ORDER OF IFL AND IFLN IN ZFRAGS
ca     Z=ZFRAG0(IFL(JSIDE),IFLN,MESON,PTS,ZMIN,ZMAX)
      Z=ZFRAGS(IFL(JSIDE),IFLN,PTS,ZMIN,ZMAX)
      PPTCL(3,I)=0.5*(Z*W(JSIDE)-PMTS(JSIDE)/(Z*W(JSIDE)))
     **(-1.)**(JSIDE+1)
      PPTCL(4,I)=0.5*(Z*W(JSIDE)+PMTS(JSIDE)/(Z*W(JSIDE)))
      IDCAY(I)=0
      IF(.NOT.(JSIDE.EQ.1)) GO TO 282
      IDENTR(NPR)=IDENT(I)
      PPTR(1,NPR)=PPTCL(1,I)
      PPTR(2,NPR)=PPTCL(2,I)
      PPTR(3,NPR)=PPTCL(3,I)
      PPTR(4,NPR)=PPTCL(4,I)
      PPTR(5,NPR)=PPTCL(5,I)
 282  IF(.NOT.(JSIDE.EQ.2)) GO TO 283
      IDENTL(NPL)=IDENT(I)
      PPTL(1,NPL)=PPTCL(1,I)
      PPTL(2,NPL)=PPTCL(2,I)
      PPTL(3,NPL)=PPTCL(3,I)
      PPTL(4,NPL)=PPTCL(4,I)
      PPTL(5,NPL)=PPTCL(5,I)
 283  IF(DIQBR) GO TO 210
      IFL(JSIDE)=-IFLN
      PX1(JSIDE)=-PX2
      PY1(JSIDE)=-PY2
      GO TO 220
C  NEW DIQUARK CREATION
210   ID1=IABS(IFLL)
      ID2=IABS(IFLN)
      IF(ID1.LE.ID2) GO TO 215
      ISWAP=ID1
      ID1=ID2
      ID1=ISWAP
215   IFL(JSIDE)=ISIGN(1000*ID1+100*ID2,IFLL)
      PX1L(JSIDE)=PX1L(JSIDE)+PXL-PX2
      PY1L(JSIDE)=PY1L(JSIDE)+PYL-PY2
      PX1(JSIDE)=PX1L(JSIDE)
      PY1(JSIDE)=PY1L(JSIDE)
 220  W(1)=W(1)-PPTCL(4,I)-PPTCL(3,I)
      W(2)=W(2)-PPTCL(4,I)+PPTCL(3,I)
      SPINT=.TRUE.
      KSPIN=2
      PARC=0.2
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GO TO 240
      IDB=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      IF(IABS(IFL(1)).GE.4.AND.IABS(IFL(1)).LT.6) PARC=0.7
      IF(IABS(IFL(2)).GE.4.AND.IABS(IFL(2)).LT.6) PARC=0.7
      IF(IABS(IFL(1)).EQ.3.OR.IABS(IFL(2)).EQ.3) PARC=0.5
      AMB=AMASS(IDB)+PARC
      GO TO 211
240   IFCN=1
      IF(RNDM(-1.).GT.0.5) IFCN=2
      IFLC1=-IFCN
      IF(IFL(1).GT.0) IFLC1=IFCN
      IFLC2=-IFLC1
      IKH1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IKH2=IDPARS(IFL(2),IFLC2,SPINT,KSPIN)
      AMB=AMASS(IKH1)+AMASS(IKH2)+PARC
211   P1X=PX1(1)+PX1(2)
      P1Y=PY1(1)+PY1(2)
      PT12=P1X**2+P1Y**2
      W12=W(1)*W(2)
      AMS2=W12-PT12
      IF(AMS2.LT.AMB**2) GO TO 100
      SPINT=.TRUE.
      KSPIN=1
      PARC1=0.2
      IF(MOD(IFL(1),100).EQ.0.AND.MOD(IFL(2),100).EQ.0) GO TO 231
      IF(IABS(IFL(1)).EQ.3.OR.IABS(IFL(2)).EQ.3) PARC1=0.5
      IF(IABS(IFL(1)).GE.4.AND.IABS(IFL(1)).LT.6) PARC1=0.7
      IF(IABS(IFL(2)).GE.4.AND.IABS(IFL(2)).LT.6) PARC1=0.7
      IDR=IDPARS(IFL(1),IFL(2),SPINT,KSPIN)
      WEND=(AMASS(IDR)+PARC1)**2
      GO TO 232
231   IKHR1=IDPARS(IFL(1),IFLC1,SPINT,KSPIN)
      IKHR2=IDPARS(IFL(2),IFLC2,SPINT,KSPIN)
      WEND=(AMASS(IKHR1)+AMASS(IKHR2)+PARC1)**2
232   SPINT=.FALSE.
      KSPIN=0
      IF(W(1)*W(2).GE.WEND) GO TO 130
         GO TO 230
225   P1X=PX1(1)+PX1(2)
      P1Y=PY1(1)+PY1(2)
      PT12=P1X**2+P1Y**2
      W12=W(1)*W(2)
      AMS2=W12-PT12
C  LAST BREAK OF STRING
 230  NPTCL=I
      AMC=SQRT(AMS2)
      EC=(W(1)+W(2))/2.0
      VC(1)=P1X/EC
      VC(2)=P1Y/EC
      VC(3)=(W(1)-W(2))/(2.0*EC)
      NIN1=NPTCL+1
      CALL CLUSTR(IFL(1),IFL(2),AMC)
      IF(LRET) GO TO 100
      NFIN1=NPTCL
      CALL LORTR(VC,NIN1,NFIN1,BACK)
      NPR=NPR+1
      NPL=NPL+1
      IF(NPR.GT.200.OR.NPL.GT.200) GO TO 9999
      IDENTL(NPL)=IDENT(NFIN1)
      PPTL(1,NPL)=PPTCL(1,NFIN1)
      PPTL(2,NPL)=PPTCL(2,NFIN1)
      PPTL(3,NPL)=PPTCL(3,NFIN1)
      PPTL(4,NPL)=PPTCL(4,NFIN1)
      PPTL(5,NPL)=PPTCL(5,NFIN1)
      IDENTR(NPR)=IDENT(NIN1)
      PPTR(1,NPR)=PPTCL(1,NIN1)
      PPTR(2,NPR)=PPTCL(2,NIN1)
      PPTR(3,NPR)=PPTCL(3,NIN1)
      PPTR(4,NPR)=PPTCL(4,NIN1)
      PPTR(5,NPR)=PPTCL(5,NIN1)
      JJ=NFIX
      DO 284 J=1,NPR
      JJ=JJ+1
      IDENT(JJ)=IDENTR(J)
      PPTCL(1,JJ)=PPTR(1,J)
      PPTCL(2,JJ)=PPTR(2,J)
      PPTCL(3,JJ)=PPTR(3,J)
      PPTCL(4,JJ)=PPTR(4,J)
      PPTCL(5,JJ)=PPTR(5,J)
284   CONTINUE
      JJ=NFIX+NPR
      DO 285 J=1,NPL
      JJ=JJ+1
      K=NPL-J+1
      IDENT(JJ)=IDENTL(K)
      PPTCL(1,JJ)=PPTL(1,K)
      PPTCL(2,JJ)=PPTL(2,K)
      PPTCL(3,JJ)=PPTL(3,K)
      PPTCL(4,JJ)=PPTL(4,K)
      PPTCL(5,JJ)=PPTL(5,K)
285   CONTINUE
      N1=NFIX+1
      N2=NFIX+NPR+NPL-1
      IF(INSIDE.NE.0)   GO  TO  385
C------------------------------------------------------C
C----- CONSTITUENT      TIME           ----------------C
C------------------------------------------------------C
      DO 286 J=N1,N2
      P3S=0.
      ES=0.
      DO 287 L=N1,J
      P3S=P3S+PPTCL(3,L)
 287  ES=ES+PPTCL(4,L)
      TI=(AMSTR-2.*P3S)/(2.*XAP)
      ZI=(AMSTR-2.*ES)/(2.*XAP)
      IF(J.NE.N2) GO TO 288
      TII=TI
      ZII=ZI
 288  PPTCL(6,J)=0.
      PPTCL(7,J)=0.
      PPTCL(8,J)=ZI
      PPTCL(9,J)=TI
 286  CONTINUE
      PPTCL(6,N2+1)=0.
      PPTCL(7,N2+1)=0.
      PPTCL(8,N2+1)=ZII
      PPTCL(9,N2+1)=TII
C]]]]]]]]]]]]
      IF(N2.LE.1) GO TO 485
      LN11=N1+1
      DO 1289 L=LN11,N2
      P8(L)=0.5*(PPTCL(8,L-1)+PPTCL(8,L))
      P9(L)=0.5*(PPTCL(9,L-1)+PPTCL(9,L))
1289  CONTINUE
      DO 1389 L=LN11,N2
      PPTCL(8,L)=P8(L)
      PPTCL(9,L)=P9(L)
1389  CONTINUE
C]]]]]]]]]]]]
      GO  TO  485
C------------------------------------------------------C
C-----  INSIDE-OUTSIDE  TIME           ----------------C
C------------------------------------------------------C
 385  CONTINUE
      DO 386 J=N1,NPTCL
      P3S=0.
      ES=0.
      NJ=J-1
      IF(NJ.EQ.0) GO TO 389
      DO 387 L=N1,NJ
      P3S=P3S+PPTCL(3,L)
 387  ES=ES+PPTCL(4,L)
 389  TI=(AMSTR-2.*P3S+PPTCL(4,J)-PPTCL(3,J))/(2.*XAP)
      ZI=(AMSTR-2.*ES-PPTCL(4,J)+PPTCL(3,J))/(2.*XAP)
      PPTCL(6,J)=0.
      PPTCL(7,J)=0.
      PPTCL(8,J)=ZI
      PPTCL(9,J)=TI
 386  CONTINUE
 485  CONTINUE
C------------------------------------------------------C
      DO 486 J=N1,NPTCL
 486  PLIDER(J)=0.
      IF(QSEE) RETURN
      IB1=IB(IDENT(N1))
      IB2=IB(IDENT(NPTCL))
      PLIDER(N1)=.667
      PLIDER(NPTCL)=.667
      IF(IB1.EQ.0) PLIDER(N1)=.5
      IF(IB2.EQ.0) PLIDER(NPTCL)=.5
      IF(PLIDER(N1).GT.0.6.AND.MOD(IFL1,100).NE.0) PLIDER(N1)=0.333
      IF(PLIDER(NPTCL).GT..6.AND.MOD(IFL2,100).NE.0)PLIDER(NPTCL)=.333
      IF(.NOT.QVSEE) RETURN
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 1387
      IF(IB1.EQ.0) PLIDER(N1)=0.
      IF(IB2.EQ.0) PLIDER(NPTCL)=0.
      IF(PLIDER(N1).GT.0.6.AND.MOD(IFL1,100).NE.0) PLIDER(N1)=0.333
      IF(PLIDER(NPTCL).GT..6.AND.MOD(IFL2,100).NE.0)PLIDER(NPTCL)=.333
      RETURN
1387  RM=RNDM(-1.)
      IF(RM.GT.0.5) PLIDER(N1)=0.
      IF(RM.LE.0.5) PLIDER(NPTCL)=0.
      RETURN
9999  WRITE(ITLIS,9998) I
9998  FORMAT(//10X,40H...STOP IN STRING..NPTCL TOO HIGH NPTCL=,I5)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER FUNCTION IDPARS(IFL01,IFL02,SPINT,IR)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C   CONSTRUCT MESON FROM QUARK AND ANTIQUARK WITH FLAVORS IFL01,IFL02
C   OR CONSTRUCT BARION FROM DIQUARK AND QUARK OR ANTIDIQUARK AND
C  ANTIQUARK WITH FLAVORS IFL01,IFL02
      COMMON/FRGSPA/ PUDS,PUDSC,PSPINS,PJSPNS,PMIX1S(3,2),PMIX2S(3,2),
     *SIGQTS,WENDM,WENDB,PBARS,PRIQS(9),PARDBS,PARQLS,PARRS,PUNDS
      LOGICAL SPINT
      IFL1=IFL01
      IFL2=IFL02
C  CONSTRUCT MESON WITH ACCOUNT FLAVOR MIXING
      IF(MOD(IFL1,100).EQ.0) GO TO 420
      IF(MOD(IFL2,100).EQ.0) GO TO 425
      JSPIN=INT(PSPINS+RNDM(-1.))
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
      ID1=IFL1
      ID2=IFL2
      IF(ID1+ID2.NE.0) GO TO 400
      RND=RNDM(-1.)
      ID1=IABS(ID1)
      ID1=INT(PMIX1S(ID1,JSPIN+1)+RND)+
     +INT(PMIX2S(ID1,JSPIN+1)+RND)+1
      ID2=-ID1
 400  IF(IABS(ID1).LE.IABS(ID2)) GO TO 410
      ISAVE=ID1
      ID1=ID2
      ID2=ISAVE
 410  IDHAD=ISIGN(100*IABS(ID1)+10*IABS(ID2)+JSPIN,ID1)
      GO TO 470
C  CONSTRUCT BARYON IDENT
 420  ID3=ISIGN(MOD(IFL1/100,10),IFL1)
      ID2=IFL1/1000
      ID1=IFL2
      GO TO 430
 425  ID3=ISIGN(MOD(IFL2/100,10),IFL2)
      ID2=IFL2/1000
      ID1=IFL1
 430  IF(IABS(ID1).LE.IABS(ID2)) GO TO 431
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 431  IF(IABS(ID2).LE.IABS(ID3)) GO TO 432
      ISWAP=ID2
      ID2=ID3
      ID3=ISWAP
 432  IF(IABS(ID1).LE.IABS(ID2)) GO TO 440
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
  440 JSPIN=1
      IF(ID1.EQ.ID2.AND.ID2.EQ.ID3) GO TO 450
      JSPIN=INT(RNDM(-1.)+PJSPNS)
      IF(SPINT.AND.IR.EQ.2) JSPIN=0
      IF(SPINT.AND.IR.EQ.1) JSPIN=1
  450 IF(JSPIN.EQ.1.OR.ID1.EQ.ID2.OR.ID2.EQ.ID3) GO TO 460
      DRND=RNDM(-1.)
      IF(DRND.GT.PJSPNS) GO TO 460
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 460  IDHAD=1000*IABS(ID1)+100*IABS(ID2)+10*IABS(ID3)+JSPIN
      IDHAD=ISIGN(IDHAD,IFL1)
 470  IDPARS=IDHAD
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORTR(V,NIN,NFIN,BACK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   LORENTZ BOOST ENERGY AND MOMENTUM OF HADRON
C
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      DIMENSION V(3)
      LOGICAL BACK
      L=1
      IF(BACK) L=-1
      VV=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
      GA=1.D0/SQRT(ABS(1.D0-VV))
      DO 100 J=NIN,NFIN
      VP=V(1)*PPTCL(1,J)+V(2)*PPTCL(2,J)+V(3)*PPTCL(3,J)
      GAVP=GA*(GA*VP/(1.+GA)-DBLE(L)*PPTCL(4,J))
      PPTCL(1,J)=PPTCL(1,J)+GAVP*V(1)
      PPTCL(2,J)=PPTCL(2,J)+GAVP*V(2)
      PPTCL(3,J)=PPTCL(3,J)+GAVP*V(3)
      PMAS=PPTCL(5,J)
      PPTCL(4,J)=SQRT(PPTCL(1,J)**2+PPTCL(2,J)**2+PPTCL(3,J)**2+
     +PMAS**2)
100   CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORCO(V,NIN,NFIN,BACK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   LORENTZ BOOST TIME AND COORDINATES OF HADRON
C
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      DIMENSION V(3)
      LOGICAL BACK
      L=1
      IF(BACK) L=-1
      VV=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
      GA=1.D0/SQRT(ABS(1.D0-VV))
      DO 100 J=NIN,NFIN
      VR=V(1)*PPTCL(6,J)+V(2)*PPTCL(7,J)+V(3)*PPTCL(8,J)
      PPTCL(6,J)=PPTCL(6,J)+GA*V(1)*(VR*GA/(GA+1.)-DBLE(L)*PPTCL(9,J))
      PPTCL(7,J)=PPTCL(7,J)+GA*V(2)*(VR*GA/(GA+1.)-DBLE(L)*PPTCL(9,J))
      PPTCL(8,J)=PPTCL(8,J)+GA*V(3)*(VR*GA/(GA+1.)-DBLE(L)*PPTCL(9,J))
      PPTCL(9,J)=GA*(PPTCL(9,J)-DBLE(L)*VR)
100   CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C          THIS SUBROUTINE UNPACKS THE IDENT CODE ID=+/-IJKL
C          I.E., GIVEN THE PARTICLE "ID", LOOKUP CHARACTERISTICS AND
C          RETURN
C
C          MESONS--
C          I=0, J<=K, +/- IS SIGN FOR J
C          ID=110 FOR PI0, ID=220 FOR ETA, ETC.
C
C          BARYONS--
C          I<=J<=K IN GENERAL
C          J<I<K FOR SECOND STATE ANTISYMMETRIC IN (I,J), EG. L = 2130
C
C          OTHER--
C          ID=1,...,6 FOR QUARKS
C          ID=9 FOR GLUON
C          ID=10 FOR PHOTON
C          ID=11,...,16 FOR LEPTONS
C          ID=20 FOR KS, ID=-20 FOR KL
C
C          I=21...26 FOR SCALAR QUARKS
C          I=29 FOR GLUINO
C          I=30 FOR PHOTINO
C          I=31...36 FOR SCALAR LEPTONS
C          I=39 FOR WINO
C          I=40 FOR ZINO
C
C          ID=80 FOR W+
C          ID=81,...,89 FOR HIGGS MESONS
C          ID=90 FOR Z0
C
C          DIQUARKS--
C          ID=+/-IJ00, I<J FOR DIQUARK COMPOSED OF I,J.
C
C          INDEX IS A SEQUENCE NUMBER USED INTERNALLY
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
*@@@@@@@@@@   SIVOKL  @@@@@@@@@@
      COMMON/FLACOM/NFLA,NFL1,NFL2,NFL3,NSPIN,NNDEX
*@@@@@@@@@@@@@@@@@@@@
      IDABS=IABS(ID)
      I=IDABS/1000
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      JSPIN=MOD(IDABS,10)
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 300
      IF(J.EQ.0) GO TO 200
      IF(I.EQ.0) GO TO 100
C          BARYONS
C         ONLY X,Y BARYONS ARE QQX, QQY, Q=U,D,S.
      IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,ID)
      IF(.NOT.(K.LE.6)) GO TO 1
        INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)+(K-1)*K*(2*K-1)/6
     1  +109*JSPIN+36*NMES+NQLEP+11
      GO TO 2
1       CONTINUE
        INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)+9*(K-7)+91
     1  +109*JSPIN+36*NMES+NQLEP+11
2     CONTINUE
      RETURN
C          MESONS
100   CONTINUE
      IFL1=0
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,-ID)
      INDEX=J+K*(K-1)/2+36*JSPIN+NQLEP
      INDEX=INDEX+11
*@@@@@@@@@@@@@ SIVOKL - TONEEV @@@@@
      IF(ID.EQ.110.OR.ID.EQ.111.OR.ID.EQ.221) GO TO 13
      IF(ID.EQ.220.OR.ID.EQ.330)GOTO 12
      RETURN
 12   IFL2=2+INT(0.25+RNDM(-1.))
      IF(IFL2.EQ.2) IFL2=1+INT(0.5+RNDM(-1.))
      IFL2=ISIGN(IFL2,ID)
      IFL3=ISIGN(IFL2,-ID)

      IF(NFLA.EQ.-1) THEN
      NFL1=IFL1
      NFL2=IFL2
      NFL3=IFL3
      NSPIN=JSPIN
      NNDEX=INDEX
      NFLA=ID
      ENDIF
      RETURN
 13   IFL2= 1+INT(0.5+RNDM(-1.))
      IFL2=ISIGN(IFL2,ID)
      IFL3=ISIGN(IFL2,-ID)

      IF(NFLA.EQ.-1) THEN
      NFL1=IFL1
      NFL2=IFL2
      NFL3=IFL3
      NSPIN=JSPIN
      NNDEX=INDEX
      NFLA=ID
c      RETURN
c in the next line NFLB was not defined 11.16.94 V.T. !
c      ELSE IF(NFLB.EQ.-1) THEN
c      MFL1=IFL1
c      MFL2=IFL2
c      MFL3=IFL3
c      MSPIN=JSPIN
c      MNDEX=INDEX
c      NFLB=ID
      ENDIF
*@@@@@@@@@@@@@@@@@@@@
      RETURN
200   CONTINUE
      IFL1=0
      IFL2=0
      IFL3=0
      JSPIN=0
      INDEX=IDABS
      IF(IDABS.LT.20) RETURN
C          DEFINE INDEX=20 FOR KS, INDEX=21 FOR KL
      INDEX=IDABS+1
      IF(ID.EQ.20) INDEX=20
C          INDEX=NQLEP+1,...,NQLEP+11 FOR W+, HIGGS, Z0
      IF(IDABS.LT.80) RETURN
      INDEX=NQLEP+IDABS-79
      RETURN
300   IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=0
      JSPIN=0
      INDEX=0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GAUSPT(PT0,SIGQT)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  GENERATE PT WITH 1/SQRT(PI*SIGQT**2)*EXP(-PT**2/SIGQT**2)
C  DISTRIBUTION
      RND=RNDM(-1.)
c                                 02.09.2000
      if(RND <= 0.0)  then
        PT0=1.0D-6
      else
        PT0=SIGQT*SQRT(-LOG(RND))
      endif
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FLAVOB(IDB,IFLQ,IFLQQ)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  COMPUTE FLAVOUR OF INTERACTING QUARK
C
      CALL FLAVOR(IDB,IFL01,IFL02,IFL03,JSPIN,INDEX)
      IF(IFL01.NE.0) GO TO 150
      IFLQ=IFL02
      IFLQQ=IFL03
      RETURN
150   INR=INT(1.+3.*RNDM(-1.))
      IF(INR.EQ.1) GO TO 100
      IF(INR.EQ.2) GO TO 200
      IFLQ=IFL03
      ID1=IABS(IFL01)
      ID2=IABS(IFL02)
      ISWAP=ID1
      GO TO 300
100   IFLQ=IFL01
      ID1=IABS(IFL02)
      ID2=IABS(IFL03)
      ISWAP=ID1
      GO TO 300
200   IFLQ=IFL02
      ID1=IABS(IFL01)
      ID2=IABS(IFL03)
      ISWAP=ID1
300   IF(ID1.LE.ID2) GO TO 400
      ID1=ID2
      ID2=ISWAP
400   IFLQQ=ISIGN((1000*ID1+100*ID2),IFLQ)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! XDISP2 removed by CMJ (XCP-3, LANL) on 09/07/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: Sets bin X1 and X2 bounds
! =====================================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION XDISTP(XMIN,IB1,IS1,IFL1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/COMABM/ ALFAM,BETAM
      COMMON/COMABB/ ALFAB,BETAB
      IFL=IABS(IFL1)
      IF(IABS(IB1).EQ.1) GO TO 3
      IF(IS1.NE.0) GO TO 33
 2    CALL SBETA(X,ALFAM,BETAM)
      IF(RNDM(-1.).GE.0.5) X=1.-X
      IF(X.LE.XMIN.OR.X.GT.1.-XMIN ) GO TO 2
      XDISTP=X
      RETURN
 33   BETAN=1.
      GO TO 1
 3    BETAN=BETAB
      IF(IFL.GT.1) BETAN=BETAB+1
 1    CALL SBETA(X,ALFAB,BETAN)
      IF(X.LE.XMIN.OR.X.GT.1.-XMIN ) GO TO 1
      XDISTP=X
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! XDIST3 removed by CMJ (XCP-3, LANL) on 09/07/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: Calculates X1, X2, and X3
! =====================================================================

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION GAMHE(ID)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C    THIS FUNCTION RETURNS THE WIDTH OF THE PARTICLE
C    WITH IDENT CODE ID
C    QUARK-BASED IDENT OF PARTICLE
C
      COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      DIMENSION GAMES0(10),GAMES1(10),GAMBR0(30),GAMBR1(30)
C    0-MESON WIDTH TABLE
C-------------- GAM(ETA) MUST BE 1.05 KEV (0.000001 GEV) ------C
      DATA GAMES0/ 2*0.,.83D-6,2*0.,0.00029,4*0./
C    1-MESON WIDTH TABLE
      DATA GAMES1/2*0.152,0.01,0.049,0.049,0.004,4*0./
C    1/2-BARION MASS TABLE
      DATA GAMBR0/-1.,2*0.,2*-1.,6*0.,3*-1.,12*0.,4*-1./
C    3/2-BARION MASS TABLE
      DATA GAMBR1/0.122,0.122,0.122,0.122,-1.,0.035,0.034,
     *0.042,-1.,0.0106,0.0091,0.,2*-1.,3*0.,-1.,3*0.,2*-1.,
     *4*0.,3*-1./
C    ENTRY
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF(IFL1.EQ.0) GO TO 100
C     BARYONS
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-11
      GAMHE=(1-JSPIN)*GAMBR0(INDEX)+JSPIN*GAMBR1(INDEX)
      RETURN
C    MESONS
100   CONTINUE
      INDEX=INDEX-36*JSPIN-NQLEP
      INDEX=INDEX-11
      GAMHE=(1-JSPIN)*GAMES0(INDEX)+JSPIN*GAMES1(INDEX)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE XSDIS(X,XMIN,XMAX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE XSEE FROM 1./XSEE*(1.-XSEE)**(BETA-1) DISTRIBUTION
C
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      BETAO=BETA
      BETA=3.
      IF(XMAX.GT.0.99999999) XMAX=0.99999999
      PE1=(1.-XMIN)**(BETA+1.)/(BETA+1.)-(1.-XMAX)**(BETA+1.)/
     *(BETA+1.)
      PE=PE1+LOG(XMAX/XMIN)
      PE1=PE1/PE
100   RND=RNDM(-1.)
      IF(RNDM(-1.).GT.PE1) GO TO 200
      X=1.-((1.-XMIN)**(BETA+1.)*(1.-RND)+(1.-XMAX)**(BETA+1.)*
     *RND)**(1./(BETA+1.))
      GO TO 300
200   X=XMIN*(XMAX/XMIN)**RND
300   PPE1=(1.-X)**BETA
      PPE2=1./X
      PPE1=PPE1+PPE2
      PPE2=PPE1*PPE2
      IF(PPE1*RNDM(-1.).GT.PPE2) GO TO 100
      BETA=BETAO
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE X2DIST(X1,X2,IFL1,IFL2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/COMABB/ ALFAB,BETAB
      BETA=BETAB-1.
      ALFA=BETAB-1.
      IF(IABS(IFL1).GT.1) ALFA=BETAB+1.
      IF(IABS(IFL2).GT.1) BETA=BETAB+1.
      XMAX=(BETA-.5)/(ALFA+BETA-1.)
      FMAX=XMAX**(BETA-.5)*(1.-XMAX)**(ALFA-.5)
 100  X1=RNDM(-1.)
      FX1=X1**(BETA-.5)*(1.-X1)**(ALFA-.5)
      IF(FMAX*RNDM(-1.).GT.FX1) GO TO 100
      X2=1.-X1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE X3DIST(X1,X2,X3,IFL1,IFL2,IFL3)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/COMABB/ ALFAB,BETAB
      BETA1=BETAB
      IF(IABS(IFL1).GT.1) BETA1=BETAB+1.
      CALL SBETA(X1,ALFAB,BETA1)
      CALL X2DIST(X2,X3,IFL2,IFL3)
      XS=1.-X1
      X2=XS*X2
      X3=XS*X3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE XCORR(IFL1,IFL2,PX1,PY1,PX2,PY2,X1,X2,
     *PSIGN,NPRODS,RETU)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   DECAY STRING OR CLUSTER
C
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/COMLID/PLIDER(499)
      COMMON/PRIMP0/ P0
       DIMENSION V(3)
       DIMENSION PPX1(3),PPX2(3),PRX1(3),PRX2(3)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON /PARCUT/ SWMAX
      COMMON /NEEDR/ NRET
      COMMON /UNCYS/ NUNCY
      COMMON/COLRET/ LRET
      LOGICAL LRET
      COMMON/COMQSE/ QSEE,QVSEE
      LOGICAL QSEE,QVSEE
      COMMON/KAPPA/ XAP
      COMMON/KEYPLA/ KEYPLA
      LOGICAL KEYPLA
      LOGICAL BACK
      LOGICAL SPINT
      LOGICAL RETU
C     INITIALIZE
      NREP=0
      SWMAX0=SWMAX
1     CONTINUE
      SWMAX=SWMAX0
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 101
      RETU=.TRUE.
      SWMAX=SWMAX0
c     WRITE(ITLIS,1200) NREP,NTRIES
1200  FORMAT(3X,' IN XCORR  NREP(=',I8,') GT NTRIES(=',I8,')')
      RETURN
101   CONTINUE
      NUNCY=0
      MXPTCL=499
      NPRODS=0
      NFIX=NPTCL
      BACK=.TRUE.
      SPINT=.TRUE.
      RETU=.FALSE.
      PTX=PX1+PX2
      PTY=PY1+PY2
      PT12=PX1**2+PY1**2
      PT22=PX2**2+PY2**2
      P1=X1*P0
      P2=X2*P0*PSIGN
      E12=P1**2+PT12+AMQ21
      IF(E12.GE.0.) GO TO 200
      RETU=.TRUE.
      RETURN
200   E22=P2**2+PT22+AMQ22
      IF(E22.GE.0.) GO TO 210
      RETU=.TRUE.
      RETURN
210   E1=SQRT(E12)
      E2=SQRT(E22)
      AMSS12=(E1+E2)**2-(P1+P2)**2-PTX**2-PTY**2
C   ARE THERE ANTIDIQUARK AND DIQUARK
      IF(MOD(IFL1,100).EQ.0.AND.MOD(IFL2,100).EQ.0) GO TO 100
       IKHR1=IDPARS(IFL1,IFL2,SPINT,2)
      PARBE=0.31
      IF(KEYPLA) PARBE=0.01
C     IF(IABS(IFL1).EQ.3.OR.IABS(IFL2).EQ.3) PARBE=0.5
      AMHR=AMASS(IKHR1)
       AMHRB=AMHR+PARBE
      IKHR=IDPARS(IFL1,IFL2,SPINT,1)
      AMHR1=AMASS(IKHR)
      IF(AMSS12.GT.AMHR**2) GO   TO  220
      RETU=.TRUE.
      RETURN
220   CONTINUE
       IF(AMSS12.GE.AMHRB**2) GO TO 400
      IF(NRET.EQ.1) GO TO 420
      NUNCY=1
      NPTCL=NPTCL+1
      IF(NPTCL.GT.MXPTCL) GO TO 9999
      IDENT(NPTCL)=IKHR1
      PPTCL(1,NPTCL)=PTX
      PPTCL(2,NPTCL)=PTY
      PPTCL(3,NPTCL)=P1+P2
      PPTCL(4,NPTCL)=E1+E2
C   !!!!!!!!!!!!!!!!!!!
C     AMHR=SQRT(AMSS12)     14.08.91 !!!
C   !!!!!!!!!!!!!!!!!!!
      PPTCL(5,NPTCL)=AMHR
      IDCAY(NPTCL)=0
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=PPTCL(4,NPTCL)/XAP
      AMT=SQRT(PTX**2+PTY**2+AMHR**2)
      PPTCL(9,NPTCL)=SQRT(2.D0)*AMT/XAP*PPTCL(4,NPTCL)/PPTCL(5,NPTCL)
      PLIDER(NPTCL)=1.
      IF(QSEE.OR.QVSEE) PLIDER(NPTCL)=0.
      NPRODS=1
      IF(SQRT(AMSS12).GE.AMHR) GO TO 419
      PPTCL(4,NPTCL)=SQRT(AMHR**2+PPTCL(1,NPTCL)**2+
     + PPTCL(2,NPTCL)**2+PPTCL(3,NPTCL)**2)
 419    RETURN
 420   RETU=.TRUE.
       RETURN
100   IFCN=1
      IF(RNDM(-1.).GT.0.5) IFCN=2
      IFLC1=-IFCN
      IF(IFL1.GT.0) IFLC1=IFCN
      IFLC2=-IFLC1
      IF(KEYPLA) SWMAX=0.
      IKH1=IDPARS(IFL1,IFLC1,SPINT,2)
      IKH2=IDPARS(IFL2,IFLC2,SPINT,2)
      AMB2=(AMASS(IKH1)+AMASS(IKH2)+SWMAX)**2
      IF(AMSS12.GE.AMB2) GO TO 350
      IF(NRET.EQ.1) GO TO 430
      NUNCY=1
      AMSS12=AMB2
      E1PE2=AMB2+(P1+P2)**2+PTX**2+PTY**2
      E2   =SQRT(E1PE2)-E1
      GO TO 350
430    RETU=.TRUE.
      SWMAX=SWMAX0
       RETURN
350   IKHR1=IDPARS(IFL1,IFLC1,SPINT,1)
      IKHR2=IDPARS(IFL2,IFLC2,SPINT,1)
      IF(AMSS12.GT.(AMASS(IKHR1)+AMASS(IKHR2)+SWMAX)**2) GO TO 500
      AMSS1=SQRT(AMSS12)
      PZ1=P1+P2
      ESS1=E1+E2
      V(1)=PTX/ESS1
      V(2)=PTY/ESS1
      V(3)=PZ1/ESS1
      NIN1=NPTCL+1
      CALL CLUSTR(IFL1,IFL2,AMSS1)
      IF(LRET) GO TO 1
      NFIN1=NPTCL
      CALL TIFILL(NIN1,NFIN1,AMSS1,IFL1,IFL2)
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=P1
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E1,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 610 J=NIN1,NFIN1
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
610   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      SWMAX=SWMAX0
      RETURN
400   AMSS1=SQRT(AMSS12)
      PZ1=P1+P2
      ESS1=E1+E2
      V(1)=PTX/ESS1
      V(2)=PTY/ESS1
      V(3)=PZ1/ESS1
      NIN1=NPTCL+1
      IF(AMSS1.GE.AMHR1+SWMAX) GO TO 600
      CALL CLUSTR(IFL1,IFL2,AMSS1)
      IF(LRET) GO TO 1
      NFIN1=NPTCL
      CALL TIFILL(NIN1,NFIN1,AMSS1,IFL1,IFL2)
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=P1
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E1,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 710 J=NIN1,NFIN1
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
710   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      SWMAX=SWMAX0
      RETURN
500   AMSS1=SQRT(AMSS12)
      PZ1=P1+P2
      ESS1=E1+E2
      V(1)=PTX/ESS1
      V(2)=PTY/ESS1
      V(3)=PZ1/ESS1
      NFIX=NPTCL
      NIN1=NPTCL+1
600    CALL STRING(IFL1,IFL2,AMSS1)
      IF(LRET) GO TO 1
       NFIN1=NPTCL
      PPX1(1)=PX1
      PPX1(2)=PY1
      PPX1(3)=P1
      BACK=.FALSE.
      CALL LORLC(V,PPX1,E1,BACK)
      CALL ANGLE(PPX1,CT,ST,CFI,SFI)
      DO 510 J=NIN1,NFIN1
      PPX1(1)=PPTCL(1,J)
      PPX1(2)=PPTCL(2,J)
      PPX1(3)=PPTCL(3,J)
      CALL ROTR(CT,ST,CFI,SFI,PPX1,PPX2,BACK)
      PPTCL(1,J)=PPX2(1)
      PPTCL(2,J)=PPX2(2)
      PPTCL(3,J)=PPX2(3)
      PRX1(1)=PPTCL(6,J)
      PRX1(2)=PPTCL(7,J)
      PRX1(3)=PPTCL(8,J)
      CALL ROTR(CT,ST,CFI,SFI,PRX1,PRX2,BACK)
      PPTCL(6,J)=PRX2(1)
      PPTCL(7,J)=PRX2(2)
      PPTCL(8,J)=PRX2(3)
510   CONTINUE
      BACK=.TRUE.
      CALL LORTR(V,NIN1,NFIN1,BACK)
      CALL LORCO(V,NIN1,NFIN1,BACK)
      NPRODS=NPTCL-NFIX
      SWMAX=SWMAX0
       RETURN
9999  WRITE(ITLIS,9998) NPTCL
9998  FORMAT(//10X,39H...STOP IN XCORR..NPTCL TOO HIGH NPTCL=,I5)
      RETURN
       END
      DOUBLE PRECISION FUNCTION CRINT(T,I)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     CALCULATION OF THE CROSS SECTION BY MEANS OF INTERPOLATION
      COMMON/TABLE/ SIGMA(30,21),ENERGY(30,4)
      COMMON/KSI/ KSI(4)
      COMMON/INTERP/ F(6)
      S=2.*T
      IF(I .GT. 6) GO TO 1
C     COMPUTE THE ENERGY TABLE NUMBER
      J=2
            GO TO 3
 1    IF(I .GT. 10) GO TO 2
      J=1
            GO TO 3
 2    IF(I.GT.19) GO TO 14
      J=3
      GO TO 3
14    J=4
C     IF  30  < T< 220 CRINT = CONST <=================
3     IF(.NOT.(T.GE.ENERGY(30,J).AND.T.LT.220.)) GO TO 33
      CRINT=SIGMA(30,I)
      RETURN
C     COMPUTE THREE POINTS FOR THE INTERPOLATION
33    L=1
 4    IF(T-ENERGY(L,J))6,5,11
 5    CRINT=SIGMA(L,I)
                         RETURN
C     (INTERPOLATION IS UNNECESSARY)
 6    IF(L .GT. 1) GO TO 7
      CRINT=SIGMA(L,I)
                 RETURN
 7    IF(L .GE. 29) GO TO 8
      L1=L-2
               GO TO 9
 8    L1=27
 9    DO 10 K=1,3
      LL1=L1+K
                 F(K)=SIGMA(LL1,I)
                                      F(K+3)=ENERGY(LL1,J)
 10   CONTINUE
                 GO TO 12
 11   L=L+1
C     IF(L.EQ.30)GOTO 13
      IF(L.EQ.30)GOTO 15
      GOTO 4
 12   CRINT=SINTER(T)
      RETURN
C13   IF(I.LE.19) GO TO 15
C     CRINT=SIGMA(30,I)
C     RETURN
15    IF(KSI(1).EQ.5) GO TO 21
      IF(IABS(KSI(2)).EQ.2) GO TO 22
      COEF=SIGMA(30,I)/(1.+0.0023*LOG(1.16*ENERGY(L,J))**2)
      CRINT=COEF*(1.+0.0023*LOG(1.16*T)**2)
      RETURN
21    IB=-1.
      IF(I.GT.20) GO TO 23
C     ANTI-BARYON-BARYON CROSS SECTION
C     TOTAL CROSS SECTION
      CRINT=PPCRSE(S,IB)
      RETURN
C     ELASTIC CROSS SECTION
23    IF(I.GT.21) GO TO 24
      CRINT=PPELSE(S,IB)
      RETURN
24    CRINT=0.
      RETURN
C     BARYON-BARYON CROSS SECTION
22    IB=1.
      IF(.NOT.(I.EQ.7.OR.I.EQ.9)) GO TO 26
C     TOTAL CROSS SECTION
      CRINT=PPCRSE(S,IB)
      RETURN
26    IF(.NOT.(I.EQ.8.OR.I.EQ.10)) GO TO 27
C     ELASTIC CROSS SECTION
      CRINT=PPELSE(S,IB)
      RETURN
27    CRINT=0.
      RETURN
               END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION PPCRSE(S,IB)
      use modifiedDCMParams, only: pi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE PP-TOTAL CROSS SECTION
C
C   M.M. BLOCK AND R.N. CAHN REV.MOD.PHYS.V.57,N2(1985) (SET 2)
C
      COMMON/RHOPP/ RHOPP
      COMMON/CREGGE/ CREGGE
      DATA A/41.30/,D/-40.51/,BETA/0.62/,ALFA/0.47/,S0/293.46/
      DATA C/8.40/,AMU/0.5/
C
      SIGN=1.
      IF(IB.EQ.-1.) SIGN=-1.
      CREGGE=D*COS(PI*ALFA/2.)*S**(ALFA-1.)
      PPCRSE=A+BETA*(LOG(S/S0)**2-PI**2/4.)+
     +C*SIN(PI*AMU/2.0)*S**(AMU-1.)+SIGN*CREGGE
      RHOPP=(PI*BETA*LOG(S/S0)-C*COS(PI*AMU/2.)*S**(AMU-1.)+SIGN*
     *D*SIN(PI*ALFA/2.)*S**(ALFA-1.))/PPCRSE
      RETURN
      END
C**********************************************************************
      DOUBLE PRECISION FUNCTION PPELSE(S,IB)
      use modifiedDCMParams, only: pi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE PP-ELASTIC CROSS SECTION
C
C   M.M. BLOCK AND R.N. CAHN REV.MOD.PHYS.V.57,N2(1985) (SET 1)
C
      COMMON/CPPTOT/ CPPTOT
      COMMON/RHOPP/ RHOPP
      DATA CP/10.90/,DP/-0.08/,EP/0.043/,CM/23.27/,DM/0.93/
C
C  COMPUTE SLOPE BPP
      BPPLUS=CP+DP*LOG(S)+EP*LOG(S)**2
      BPMIN=CM+DM*LOG(S)
C  COMPUTE ELASTIC CROSS SECTION
      IBB=-1
      CPPTPL=PPCRSE(S,IBB)
      IF(.NOT.(IB.EQ.IBB)) GO TO 100
      RHO=RHOPP
      CPPTOT=CPPTPL
100   IBB=1
      CPPTMI=PPCRSE(S,IBB)
      IF(.NOT.(IB.EQ.IBB)) GO TO 200
      RHO=RHOPP
      CPPTOT=CPPTPL
200   SIGMAP=CPPTPL+CPPTMI
      SIGMAM=CPPTPL-CPPTMI
      BPP=BPPLUS+(SIGMAM/SIGMAP)*(BPMIN-BPPLUS)
      PPELSE=CPPTOT**2*(1.+RHO**2)/(16.*PI*BPP*0.389)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SMARK(IK01,IK02)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C  COMPUTE CROSS SECTION TYPE KSI(1),SUMMARY BARION
C  NUMBER KSI(2),SUMMARY CHARGE KSI(3),SUMMARY STRANGE KSI(4)
C  ----------IK2-BARYON & MESON ------
      COMMON/KSI/ KSI(4)
      IK1=IK01
      IK2=IK02
      IB1=IB(IK1)
      IB2=IB(IK2)
      IF(IB1.GE.0.OR.IB2.GE.0) GO TO 110
C     AB-AB -> B-B
      IK1=IABS(IK1)
      IK2=IABS(IK2)
      IF(IK2.EQ.1120.OR.IK2.EQ.1220) GO TO 112
      IF(IK1.NE.1120.AND.IK1.NE.1220) GO TO 111
      IK11=IK1
      IK1=IK2
      IK2=IK11
      GO TO 112
110   IF(IB1.NE.0.OR.IB2.GE.0) GO TO 111
C     M-AB -> AM-B
      IK2=IABS(IK2)
      IF(IK1/100.NE.MOD(IK1,100)/10) IK1=-IK1
111   IF(IK2.NE.1220) IK2=1120
112   IB1=IB(IK1)
      IB2=IB(IK2)
      IQ1=INT(CHARGE(IK1)*1.01)
      IQ2=INT(CHARGE(IK2)*1.01)
      IS1=IS(IK1)
      IS2=IS(IK2)
      KSI(2)=IB1+IB2
      KSI(3)=IQ1+IQ2
      KSI(4)=IS1+IS2
      IF(KSI(4).EQ.0) GO TO 10
      IF(KSI(2).LE.1) GO TO 7
C  STRANGE BARYON NUCLEON COLLISION
      IF(KSI(4)+2)1,2,3
C  OMEGA NUCLEON COLLISION
 1    KSI(1)=1
                 RETURN
C  KSI NUCLEON COLLISION
 2    KSI(1)=1
      IF(KSI(3).EQ.-1)KSI(1)=2
      RETURN
C  SIGMA OR LAMBDA NUCLEON COLLISION
 3    IF(KSI(3).EQ.0.OR.KSI(3).EQ.1) GO TO 4
      KSI(1)=1
                 RETURN
 4    IF(IQ1)6,5,6
 5    KSI(1)=3
                 RETURN
 6    KSI(1)=2
                 RETURN
7     IF(IB1.EQ.0.OR.IB2.EQ.0) GO TO 17
C    ANTIBARYON NUCLEON COLLISION
      KSI(1)=5
      RETURN
C  STRANGE MESON NUCLEON COLLISION
 17   IF(KSI(4)-1)9,8,9
 8    KSI(1)=1
      IF(IQ1.EQ.1) RETURN
      KSI(1)=3
                 RETURN
 9    KSI(1)=2
      IF(KSI(3).EQ.0)RETURN
      KSI(1)=3
      K=IK1
      NF=0
      IF(K.EQ.-230.OR.K.EQ.-231) NF=1
      IF(KSI(3).EQ.-1.OR.NF.EQ.1)KSI(1)=4
      RETURN
C  NONSTRANGE PARTICLE COLLISION
 10   IF(KSI(2).EQ.2) GO TO 13
      IF(IB1.EQ.0.OR.IB2.EQ.0) GO TO 18
C   ANTIBARYON NUCLEON COLLISION
      KSI(1)=5
      RETURN
C  MESON NUCLEON COLLISION
18    IK=KSI(3)+2
      GOTO(11,12,12,11),IK
 11   KSI(1)=1
                 RETURN
 12   KSI(1)=2
      IF(IQ1.EQ.0) KSI(1)=3
      RETURN
C  BARION NUCLEON COLLISION
 13   K=IABS(IK1)
      IF(K.EQ.1120.OR.K.EQ.1220) GO TO 14
      IF(KSI(3).EQ.3.OR.KSI(3).EQ.-1) GO TO 16
      KSI(1)=2
      IF(KSI(3).NE.1) RETURN
      KSI(1)=3
                 RETURN
 14   KSI(1)=1
      IF(KSI(3).EQ.1)KSI(1)=2
      RETURN
 16   KSI(1)=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SFICRI(T,IKS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF CROSS SECTION FOR FI0 NUCLEON INTERACTION
      COMMON/DATA10/ PAROM,PARKSI,PARSIG,PARF0
      GO TO (1,2,3,3,4,4,5),IKS
 1    SIF=(CRINT(T,1)+CRINT(T,3))/2.
                                      GOTO 6
 2    SIF=(CRINT(T,2)+CRINT(T,4)-CRINT(T,5))/2.
                                                  GO TO 6
 3    SIF=CRINT(T,IKS+2)
                           GO TO 6
 4    SIF=(CRINT(T,IKS+10)+CRINT(T,IKS+12))/2.
                                                 GOTO 6
 5    SIF=(CRINT(T,15)+CRINT(T,16)+CRINT(T,17)+
     * CRINT(T,18)+CRINT(T,19))/2.
 6    SFICRI=SIF*PARF0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SBCRI(T,IKS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF CROSS SECTION FOR STRANGE BARYON NUCLEON INTERACTION
      COMMON/DATA10/ PAROM,PARKSI,PARSIG,PARF0
      COMMON/KSI/ KSI(4)
      IF(KSI(4)+2)3,2,1
 1    CS=PARSIG
                  GO TO 4
 2    CS=PARKSI
                  GO TO 4
 3    CS=PAROM
 4    IF(KSI(1)-2)5,9,12
C  POSITIVE SIGMA PROTON REACTION
 5    GO TO(6,6,7,7,6,6,8),IKS
 6    SIG=CRINT(T,IKS+6)
                           GO TO 14
 7    SIG=0.
               GO TO 14
 8    SIG=CRINT(T,11)+CRINT(T,12)
                                    GO TO 14
C  NEGATIVE SIGMA PROTON REACTION
 9    GOTO(10,10,7,7,10,10,11),IKS
 10   SIG=CRINT(T,IKS+8)
                           GO TO 14
 11   SIG=2.*CRINT(T,13)+CRINT(T,14)
                                       GO TO 14
C  NEUTRAL SIGMA PROTON REACTION
 12   GOTO(13,13,7,7,13,13,11),IKS
 13   SIG=(CRINT(T,IKS+8)+CRINT(T,IKS+6))/2.
 14   SBCRI=SIG*CS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SCRIK(T1,IKS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF CROSS SECTION FOR STRANGE KAON NUCLEON REACTION
      COMMON/KSI/ KSI(4)
      T=T1
      P0=SQRT(T*T+2.*0.497*T)
      IF(T1.GE.20.0.AND.IKS.EQ.3) GO TO 4
      IF(T.GE.49.5) T=49.5
      IK=KSI(1)
      GO TO(1,5,9,13),IK
C  POSITIVE KAON PROTON REACTION
 1    GO TO(2,3,4,4,4,4,4),IKS
 2    IF(P0.GE.10.) SCRIK=17.042+0.0188*T+0.000093*T**2
      IF(P0.LT.0.8) SCRIK = 13.0
      IF(P0.GE.0.8.AND.P0.LT.10.)
     * SCRIK=23.4/(P0+1.)+15.-7.5/(P0-0.3)
               GO TO 16
 3    IF(P0.GE.10.) SCRIK=3.77-0.047*T+0.00035*T**2
      IF(P0.LT.0.8) SCRIK = 13.0
      IF(P0.GE.0.8.AND.P0.LT.10.) SCRIK=23.4/(P0+1.)
               GO TO 16
 4    SCRIK=0.
               GO TO 16
C  NEGATIVE KAON PROTON REACTION CROSS SECTION
 5    GO TO(6,7,8,4,4,4,4),IKS
 6    IF(P0.GE.4.) SCRIK=23.97-0.179*T+0.0022*T**2
      IF(P0.LT.0.25) SCRIK=60.
      IF(P0.GE.0.25.AND.P0.LT.0.8) SCRIK=15./P0 + 13./P0
      IF(P0.GE.0.8.AND.P0.LT.4.)
     *SCRIK=15./P0+9.-1.2/(P0-.6)+13./SQRT(P0)
               GO TO 16
 7    IF(P0.GE.4.) SCRIK=4.45-0.103*T+0.0012*T**2
      IF(P0.LT.0.25) SCRIK=60.
      IF(P0.GE.0.25.AND.P0.LT.4.) SCRIK=15./P0
               GO TO 16
 8    SCRIK=0.23-0.022*T+0.00053*T**2
               GO TO 16
C  POSITIVE KAON NEUTRON REACTION CROSS SECTION
 9    GO TO(10,11,12,4,4,4,4),IKS
 10   IF(P0.GE.10.) SCRIK=17.49+0.0192*T+0.000039*T**2
      IF(P0.LT.0.8) SCRIK = 13.0
      IF(P0.GE.0.8.AND.P0.LT.10.)
     * SCRIK=23.4/(P0+1.)+15.-7.5/(P0-0.3)
               GO TO 16
 11   IF(P0.GE.10.) SCRIK=3.77-0.047*T+0.00035*T**2
      IF(P0.LT.0.8) SCRIK = 13.0
      IF(P0.GE.0.8.AND.P0.LT.10.) SCRIK=23.4/(P0+1.)
               GO TO 16
 12   SCRIK=0.9487-0.176*T+0.0084*T**2
               GO TO 16
C  NEGATIVE MESON NEUTRON REACTION
 13   GO TO(14,15,4,4,4,4,4),IKS
 14   IF(P0.GE.4.) SCRIK=21.53-0.114*T+0.0014*T**2
      IF(P0.LT.0.25) SCRIK=60.
      IF(P0.GE.0.25.AND.P0.LT.0.8) SCRIK=15./P0
      IF(P0.GE.0.8.AND.P0.LT.4.)
     *SCRIK=15./P0+9.-1.2/(P0-.6)+13./SQRT(P0)
               GO TO 16
 15   IF(P0.GE.4.) SCRIK=19.51-10.2*T+1.52*T**2
      IF(P0.LT.0.25) SCRIK=60.
      IF(P0.GE.0.25.AND.P0.LT.4.) SCRIK=15./P0
 16   IF(T1.LE.49.5)  RETURN
      COEF=SCRIK/(1.+0.0023*LOG(1.16*T)**2)
      T=T1
      SCRIK=COEF*(1.+0.0023*LOG(1.16*T)**2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SIGMATQ(T,IKS,IK1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     CALCULATION OF THE CROSS-SECTION
      COMMON/KSI/ KSI(4)
      IF(KSI(1).EQ.5) GO TO 26
      IS1=IS(IK1)
C  SEPARATION STRANGE PARTICLES
      IF(IS1.NE.0) GOTO 24
C     SEPARATION OF THE MESON AND NUCLEON REACTION
      IF(KSI(2) .GT. 1) GO TO 17
C     MESON-NUCLEON REACTIONS
      IF(KSI(1)-2)1,7,11
C     POSITIVE MESON+PROTON
 1    GO TO (2,2,3,4,5,5,6),IKS
 2    SIGMATQ=CRINT(T,IKS)
                            RETURN
 3    SIGMATQ=0.
                  RETURN
 4    SIGMATQ=CRINT(T,6)
                          RETURN
 5    SIGMATQ=CRINT(T,IKS+10)
                               RETURN
 6    SIGMATQ=CRINT(T,15)+CRINT(T,16)
                                       RETURN
C     NEGATIVE MESON+PROTON
 7    GO TO (8,8,8,8,9,9,10),IKS
 8    SIGMATQ=CRINT(T,IKS+2)
                              RETURN
 9    SIGMATQ=CRINT(T,IKS+12)
                               RETURN
 10   SIGMATQ=CRINT(T,17)+CRINT(T,18)+CRINT(T,19)
                                                   RETURN
C     NEUTRAL MESON+PROTON
 11   KN=IABS(IK1)
      IF(KN.NE.331) GO TO 36
C  NEUTRAL FI-MESON NUCLEON REACTION
      SIGMATQ=SFICRI(T,IKS)
                             RETURN
 36   GOTO(12,13,14,14,15,15,16),IKS
 12   SIGMATQ=(CRINT(T,1)+CRINT(T,3))/2.
                                          RETURN
 13   SIGMATQ=(CRINT(T,2)+CRINT(T,4)-CRINT(T,5))/2.
                                                     RETURN
 14   SIGMATQ=CRINT(T,IKS+2)
                              RETURN
 15   SIGMATQ=(CRINT(T,IKS+10)+CRINT(T,IKS+12))/2.
                                                    RETURN
 16   SIGMATQ=(CRINT(T,15)+CRINT(T,16)
     *+CRINT(T,17)+CRINT(T,18)+CRINT(T,19))/2.
                                                 RETURN
C     NUCLEON-NUCLEON REACTIONS
C
 17   IF(KSI(1)-2)18,21,21
C     PROTON+PROTON
 18   GO TO (19,19,3,3,19,19,20),IKS
 19   SIGMATQ=CRINT(T,IKS+6)
                              RETURN
 20   SIGMATQ=CRINT(T,11)+CRINT(T,12)
                                       RETURN
C     NEUTRON +PROTON
 21   GO TO (22,22,3,3,22,22,23),IKS
 22   SIGMATQ=CRINT(T,IKS+8)
                              RETURN
 23   SIGMATQ=2.*CRINT(T,13)+CRINT(T,14)
                                          RETURN
C  SEPARATION OF MESON AND BARYON REACTION
 24   IF(KSI(2).GT.1)GO TO 25
C  STRANGE MESON NUCLEON REACTION
      SIGMATQ=SCRIK(T,IKS)
                            RETURN
C  STRANGE BARYON NUCLEON REACTION
 25   SIGMATQ=SBCRI(T,IKS)
                            RETURN
C    ANTIBARYON NUCLEON REACTION
26    SIGMATQ=APPCRI(T,IKS)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SINTER(X)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C   INTERPOLATION BY MEANS OF LAGRANGE POLINOMS
      COMMON/INTERP/ F(6)
      SINTER=0.
      DO 2 I=1,3
      R=1.0
      K=3+I
      DO 1 J=1,3
      IF(I.EQ.J) GO TO 1
      L=3+J
      R=R*(X-F(L))/(F(K)-F(L))
 1    CONTINUE
 2    SINTER=SINTER+R*F(I)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION  CRNNDN(S)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C  NN--NDELTA CROSS SECTION
C
      SQS=SQRT(S)
      IF(SQS.GT.2.015) GO TO 100
      CRNNDN=0.
      RETURN
 100  CRNNDN=20.*(SQS-2.015)**2/(0.015+(SQS-2.015)**2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION CRNDNN(S,IE)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C       X2=Mnucl+Mpion
      DATA X2 /1.080/
C
C  N+DELTA->N+N' CROSS SECTION (NEW) with D-B correction 04.12.92
C
      T=(S-(2.*.940)**2)/1.88
      P0=SQRT(T*(T+1.88))
      PIN=(S+0.63)**2/(4.*S)-1.51
      PF=S/4.-0.88
c------------- degeneration factor ------------------
      G=4.D0
      IF(IE.EQ.1)    G=2.D0
c-------------Danielewitz-Bertsch correction (appr.) -----------
      X1=SQRT(S)-.940
      GD1=GDM(X1)
      ARG1=(X1**2-1.232**2)/(1.232*GD1)
      ARG2=(X2**2-1.232**2)/(1.232*0.118)
      COR=1./3.1415927*(ATAN(ARG1)-ATAN(ARG2))
c---------------------------------------------------------------
c     CRNDNN=PF*CRNNDN(S)/(8.*PIN)    Amel
      CRNDNN=PF*CROSS(P0,1,5)/(G*PIN*COR)
      if(CRNDNN.LE.0.) CRNDNN=0.
      RETURN
      END
C******************************************************************
      SUBROUTINE CROSEC(ITIN,IK01,IK02,PX01,PY01,PZ01,AM01,
     *PX02,PY02,PZ02,AM02,SIGMA,INEL)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C  CALCULATION OF TOTAL(ITIN=1) OR ELASTIC(ITIN=0)
C  OR CHARGE-EXCHANGE (ITIN=2) OR ANNIHILATION(ITIN=3)
C  HADRON-HADRON CROSS SECTION
C   IF INEL=1 ==> SIGTOT=SIGTOT-SIGEL
C   IF INEL=2 ==> SIGTOT AND SIGEL ARE CONSTANTS
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/ IPRINT
      LOGICAL IPRINT
      ITOT=ITIN
      IB1=IB(IK01)
      IB2=IB(IK02)
      IB12=IB1*IB2
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 10
      IF((IK01.EQ.1120.OR.IK01.EQ.1220) .AND.
     *    IK02.NE.1120.AND.IK02.NE.1220) GO TO 13
      IF(IB1.GT.0.AND.IB2.LE.0) GO TO 13
      IF(IB1.LT.0.AND.IB2.EQ.0) GO TO 13
 10   IK1=IK01
      IK2=IK02
      PX1=PX01
      PY1=PY01
      PZ1=PZ01
      AM1=AM01
      PX2=PX02
      PY2=PY02
      PZ2=PZ02
      AM2=AM02
      GO TO 20
13    IK1=IK02
      IK2=IK01
      PX1=PX02
      PY1=PY02
      PZ1=PZ02
      AM1=AM02
      PX2=PX01
      PY2=PY01
      PZ2=PZ01
      AM2=AM01
 20   CONTINUE
      SIGMA=0.
      IF(ITOT.EQ.3.AND.IB12.GE.0) RETURN
      IB1=IB(IK1)
      IB2=IB(IK2)
      IS1=IABS(IS(IK1))
      IS2=IABS(IS(IK2))
C
      AM1N=0.139
      IF(IS1.NE.0) AM1N=0.497
      IF(IB1.NE.0) AM1N=0.939
      IF(IB1.NE.0.AND.IS1.NE.0) AM1N=1.1156
C
      COEFMS=1.
C--- IF MESON-MESON COLL. SIGMA=(2/3)*SIGMA
      IF(IB2.EQ.0) COEFMS=0.667
      IF(IS2.NE.0.AND.IB12.GE.0) COEFMS=COEFMS*0.81**IS2
C--- IF STRANGE HADRONS   SIGMA=0.81**(IS1+IS2)*SIGMA
      IF(IB12.LT.0) COEFMS=COEFMS*0.81**(IS1+IS2)
      E1=SQRT(AM1**2+PX1**2+PY1**2+PZ1**2)
      E2=SQRT(AM2**2+PX2**2+PY2**2+PZ2**2)
      S=AM1**2+AM2**2+2.*E1*E2-2.*(PX1*PX2+PY1*PY2+PZ1*PZ2)
      SQR=SQRT(S)
c !!!
      TKINM=(S-AM1**2-AM2**2)/(2.*AM2)-AM1
c !!!
      DSR=0.20
C     IF(IS1.EQ.0.AND.IS2.EQ.0.AND.IB1.EQ.1.AND.IB2.EQ.1) DSR=0.17
C
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 21
      AM2N=0.93900000
      TKIN=(S-AM1N**2-AM2N**2)/(2.*AM2N)-AM1N
C
      IF(TKIN.LE.0.0001.AND.IB12.GE.0) RETURN
C
C  M--N AND N--N
      IF(SQR-AM1N-AM2N.LT.DSR.AND.IB12.GE.0) ITOT=0
C
      IKS=4
      GO TO 22
C PI---PI
21    IF(TKINM.LE.0.0001) RETURN
      IF(SQR-AM1-AM2.LT.DSR) ITOT=0
C
22    IF(ITOT.EQ.0) IKS=2
      IF(ITOT.EQ.1) IKS=1
      IF(ITOT.EQ.2) IKS=3
      IF(ITOT.EQ.0) SIGCON=7.
      IF(ITOT.EQ.1) SIGCON=39.
      IF(IB1.EQ.0) SIGCON=0.6667*SIGCON
      IF(ITOT.EQ.3) GO TO 35
      IF(IB12.GE.0) GO TO 40
C
      IF(TKIN.LT.0.0001.AND.(IKS.EQ.2.OR.IKS.EQ.3)) RETURN
C]]]]]]]]]PBAR+P]]]]]]]]]]]]]]
35    IK1=-1120
      IK2=1120
40    SIGMA=SIGCON
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 500
      IF(IB1.GE.0.AND.IB2.EQ.1.AND.AM2.GT.0.96) GO TO 500
      CALL SMARK(IK1,IK2)
      IF(IK1.EQ.-1120.AND.IK2.EQ.1120)  GO TO 1825
      IF(INEL.NE.2) GO TO 400
      SIGMA=SIGCON
      GO TO 500
400   SIGMA=SIGMATQ(TKIN,IKS,IK1)
C
c !!!         20.10.92T
      IF(ITOT.NE.1.OR.SQR.GT.3.0)       GO TO 500
      IF(IK1.EQ.1111.OR.IK1.EQ.1121.OR.IK1.EQ.2221.OR.IK1.EQ.1221) THEN
        IE1=INT(CHARGE(IK1)*1.001)
        IE2=INT(CHARGE(IK2)*1.001)
        IE12=IE1+IE2
        IF(IE12.EQ.3.OR.IE12.EQ.-1)      GO TO 500
        P0=SQRT(TKIN*(TKIN+2.*AM2))
        SBIN=CROSS(P0,1,5)
        SBIND=CRNDNN(S,IE12)
        SIGMA=SIGMA-SBIN+SBIND
      END IF
c !!!
500   SIGMA=SIGMA*COEFMS
C   FOR LOW-ENERGY PION-PION COLLISIONS (A LA WOLF)
c     IF(AM1.LT.0.141.AND.AM2.LT.0.141) THEN
c       IF(SQR-AM1N-AM2N.LT.0.800.AND.ITIN.LE.1) SIGMA=30.
c       GO TO  1827
c     ENDIF
C   FOR LOW-ENERGY PION-PION COLLISIONS (A LA Bao-An Li)
        IF(AM1.LT.0.141.AND.AM2.LT.0.141) THEN
          IF(SQR-AM1N-AM2N.LT.1.000.AND.ITIN.LE.1)
     &   SIGMA=SIGMA+PIPICS(IK1,IK2,SQR,ITOT)
          GO TO  1827
        ENDIF
C
      GO TO 1826
1825  SIGMA=SIGMATQ(TKIN,IKS,IK1)*COEFMS
1826  CONTINUE
      IF(IB12.LT.0.AND.TKIN.LT.0.02) TKIN=0.02
C]]]]]]]]]]]PBAR+P ]]]]]]]]]]]]
      IF(INEL.EQ.1.AND.ITOT.EQ.1)
     *SIGMA=SIGMA-SIGMATQ(TKIN,2,IK1)*COEFMS
1827  CONTINUE
      IF(IPRINT.AND.ITOT.EQ.1) WRITE(ITLIS,2001) SIGMA
     *,IK01,IK02,PX01,PY01,PZ01,PX02,PY02,PZ02
      IF(IPRINT.AND.ITOT.EQ.0) WRITE(ITLIS,2000) SIGMA
     *,IK01,IK02,PX01,PY01,PZ01,PX02,PY02,PZ02
      IF(IPRINT.AND.ITOT.EQ.3) WRITE(ITLIS,2002) SIGMA
     *,IK01,IK02,PX01,PY01,PZ01,PX02,PY02,PZ02
2000  FORMAT('  SIGMA ELASTIC =',F10.4,2I6,6F12.4)
2001  FORMAT('  SIGMA TOTAL   =',F10.4,2I6,6F12.4)
2002  FORMAT('  SIGMA ANNIH.  =',F10.4,2I6,6F12.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION PIPICS(ID1,ID2,SS,ICS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      REAL*8 MPI,MR,MS
      DATA MPI/0.140/,PI/3.141592/
C  to transforme 1/GeV**2 -> mb
	MR=0.770
	MS=5.8*MPI
      C=10./(5.06**2)
      PIPICS=0.
      Q2=(SS/2.)**2-MPI**2
      IF(Q2.LE.0.)  RETURN
      Q=SQRT(Q2)
      GR=0.095*Q*((Q/MPI)/(1.+(Q/MR)**2))**2
	GS=2.060*Q
	S2D11=(GR/2.)**2/((MR-SS)**2+(GR/2.)**2)
	S2D00=(GS/2.)**2/((MS-SS)**2+(GS/2.)**2)
	S2D20=SIN(-0.12*Q/MPI)**2
      S00=8.*C*PI/(Q**2)*1.*S2D00
      S11=8.*C*PI/(Q**2)*3.*S2D11
      S20=8.*C*PI/(Q**2)*1.*S2D20
      PIPICS=0.
      IF(ICS.EQ.-1)    THEN
C         ISOSPIN AVERAGED CROSS SECTION
	  PIPICS=1./9.*S00+1./3.*S11+5./9.*S20
	  RETURN
      ELSEIF(ICS.EQ.1.OR.ICS.EQ.0) THEN
C         TOTAL=ELASTIC CROSS SECTION
C        PI+ PI-  OR PI- PI-
	  IF((ID1.EQ.120.AND.ID2.EQ.120).OR.
     &       (ID1.EQ.-120.AND.ID2.EQ.-120))
     &    PIPICS=S20
C        PI+ PI0  OR PI- PI0 OR PI0 PI- OR PI0 PI+
	  IF((IABS(ID1).EQ.120.AND.ID2.EQ.110).OR.
     &	     (ID1.EQ.110.AND.IABS(ID1).EQ.120))
     &    PIPICS=S20/2.+S11/2.
C        PI+ PI-  OR PI- PI+
	  IF((ID1.EQ.120.AND.ID2.EQ.-120).OR.
     &       (ID1.EQ.-120.AND.ID2.EQ.120))
     &    PIPICS=S20/6.+S11/2.+S00/3.
C        PI0 PI0
	  IF(ID1.EQ.110.AND.ID2.EQ.110)
     &    PIPICS=S00/3.+2.*S20/3.
      ELSEIF(ICS.EQ.2) THEN
C CHARGE-EXCHANGE PI+PI- <-> PI0PI0
	  IF((ID1.EQ.120.AND.ID2.EQ.-120).OR.
     &       (ID1.EQ.-120.AND.ID2.EQ.120).OR.
     &       (ID1.EQ.110.AND.ID2.EQ.110))
     &    PIPICS=S00/3.+S20/3.
      ELSEIF(ICS.EQ.3) THEN
C         RHO FORMATION CROSS SECTION
	  PIPICS=S11/2.
      ENDIF
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SIGFIH(SIG,I1,I2,I3,I4)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C      FILL /HADSIG/
C
C   I1 IS ABSOLUTE VALUE OF PROJECTILE IDENT
C   I2 IS ABSOLUTE VALUE OF TARGET IDENT
C   I3 MARK REACTION TYPE
C   I4 IS SIGN OF PROJECTILE IDENT
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      COMMON/HADSIG/SIGS(100),SIGEVT,NSIGS,INOUT(2,100)
      COMMON/COMMUL/ MULTP
      LOGICAL MULTP
      COMMON/CSIGSU/ SIGSUM
      COMMON/CPRSIG/ ISGCOL
      CHARACTER*8 AGH1H2(5,11)
      DATA AGH1H2 /
     * 'TRIPLE P','OMERON V','ERTEX DI','FRACTION',' DIAGRAM',
     * 'PLANAR (','ONE REGG','EON) DIA','GRAM    ','        ',
     * 'UNDEVELO','PED CYLI','NDER DIA','GRAM    ','        ',
     * 'ELASTIC ','SCATTERI','NG DIAGR','AM      ','        ',
     * 'ANNIHILA','TION DIA','GRAM    ','        ','        ',
     * 'SMALL MA','SS DIFRA','CTION DI','AGRAM   ','        ',
     * 'MULTI PO','MERON SC','ATTER.(C','YLINDER)',' DIAGRAM',
     * 'TWO PART','ICLE REA','CTION   ','        ','        ',
     * 'MULTI PO','MERON SC','ATTER.(C','HAINS  )',' DIAGRAM',
     * 'TRIPLE R','EGGEON V','ERTEX DI','AGRAM   ','        ',
     * 'DOUBLE D','IFRACTIO','N DIAGRA','M       ','        '/
      IOPAK=10000
      NSIGS=NSIGS+1
      SIGS(NSIGS)=SIG
      SIGSUM=SIGSUM+SIG
C   @@@@@@@@@@@@@@@@@@@@
C     INOUT1(1,NSIGS)=I1
C     INOUT2(1,NSIGS)=I2
C     INOUT1(2,NSIGS)=I3
C     INOUT2(2,NSIGS)=I4
      INOUT(1,NSIGS)=I1+IOPAK*I2
      INOUT(2,NSIGS)=I3+IOPAK*I4
C   @@@@@@@@@@@@@@@@@@@@
      K3=I3
      IF(I3.EQ.7.AND.MULTP) K3=9
      IF(ISGCOL.EQ.0.AND.IPRINT)
     *WRITE(ITLIS,1000) I3,SIG,(AGH1H2(J1,K3),J1=1,5)
1000  FORMAT(10X,17HREACTION TYPE I3=,I2,25H  REACTION CROSS SECTION=,
     *E10.4,5HMB - ,5A8)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION APPCRI(T,IKS)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C    COMPUTE ANTIPROTON-PROTON CROSS SECTION
C
      GO TO (100,200,300,400),IKS
C     TOTAL CROSS SECTION
100   INUM=20
      APPCRI=CRINT(T,INUM)
      RETURN
C   ELASTIC CROSS SECTION
200   INUM=21
      APPCRI=CRINT(T,INUM)
      RETURN
C    CHARGE-EXCHANGE CROSS SECTION
300   APPCRI=0.
      RETURN
C   ANNIHILATION CROSS SECTION
400   IF(T.GE.0.9) GO TO 500
      IF(T.GE.0.02) GO TO 502
      APPCRI=162.0
      RETURN
502   PMOD=SQRT(T*(T+1.88))
      APPCRI=61.6/PMOD**0.6
      RETURN
500   INUM=20
      APPCRI=CRINT(T,INUM)
      INUM=7
      PPCRI=CRINT(T,INUM)
      APPCRI=APPCRI-PPCRI
      IF(APPCRI.LE.0.) APPCRI=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PLANAR(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE PLANAR DIAGRAM
C
C     DECAY Q-QBAR OR QQ-QQBAR STRING
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/PRIMP0/ P0
      COMMON /NEEDR/ NRET
      LOGICAL RETU
      COMMON/COMANN/ DIQAN
      LOGICAL DIQAN
      COMMON/KEYPLA/ KEYPLA
      LOGICAL KEYPLA
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/COMASS/ AM1,AM2
      DIMENSION IFORD(3),IFORDP(3)
      DIMENSION IMQ1(2),IMQ2(2)
      COMMON/FLACOM/NFLA,NFL1,NFL2,NFL3,NSPIN,NNDEX
      ZER=0.0
      NREP=0
      KEYPLA=.TRUE.
      NPOLD=NPTCL
      P0OLD=P0
1111  CONTINUE
      P0=P0OLD
      NPTCL=NPOLD
      NREP=NREP+1
      IF(NREP.LE.NTRIES) GO TO 201
C     WRITE(ITLIS,1200) IDIN(1),IDIN(2),P0
1200  FORMAT(1X,'IN PLANAR:NREP > NTRIES,IK1,IK2,P0=',
     *2I5,1X,F7.3)
      IRET=1
      RETURN
201   CONTINUE
      IPACK=1000
      NRET=1
C    INITIALIZE
 110  RETU=.FALSE.
      IRET=0
      PSIGN=-1.
      P0=HALFE
      IKA=IDIN(1)
      IKB=IDIN(2)
      IF(IB(IKA).EQ.0.AND.IB(IKB).EQ.0) NRET=0
      I1=IABS(IDIN(1))
      IFL=I1/1000
      IF(IFL.NE.0) GO TO 260
C
C    INITIAL MESON
      IQ=0
      JFL=MOD(I1/100,10)
      KFL=MOD(I1/10,10)
      J=ISIGN(JFL,IDIN(1))
      K=ISIGN(KFL,-IDIN(1))
*@@@@@@@@@@@@@@@@@ SIVOKL @@@@@@
      IF(IKA.EQ.110.OR.IKA.EQ.111.OR.IKA.EQ.221.
     *OR.IKA.EQ.220.OR.IKA.EQ.330) THEN
      J=NFL2
      K=NFL3
      ENDIF
*@@@@@@@@@@@@@@@@@
      IFLU=J
      IFLQM=K
      IF(K.GT.0) GO TO 150
      IFLQM=J
      IFLU=K
150   IAFLU=IABS(IFLU)
*@@@@@@@@@@@@@@@@@ SIVOKL @@@@@@
      IF(IKB.EQ.110.OR.IKB.EQ.111.OR.IKB.EQ.221.
     *OR.IKB.EQ.220.OR.IKB.EQ.330) THEN
      IFORD(2)=NFL2
      IFORD(3)=NFL3
      ELSE
      CALL FLAVOR(IKB,IFORD(1),IFORD(2),IFORD(3),JSPIN,INDEX)
      ENDIF
*@@@@@@@@@@@@@@@@@
C     CALL FLAVOR(IKB,IFORD(1),IFORD(2),IFORD(3),JSPIN,INDEX)
      IF(IB(IKB).EQ.0) GO TO 101
      ISIG=ISIGN(1,IKB)
      IF(ISIG.LT.0) GO TO 255
      DO 155 I=1,3
      IF(IAFLU.NE.IFORD(I)) GO TO 155
      IQ=I
      GO TO 156
155   CONTINUE
      GO TO 156
255   CONTINUE
      DO 256 I=1,3
      IF(IFLQM.NE.IABS(IFORD(I))) GO TO 256
      IFLQM0=IFLQM
      IFLQM =IFLU
      IFLU  =IFLQM0
       IQ=I
       GO TO 156
256   CONTINUE
156   IF(IQ-2) 10,20,30
10    ID1=IFORD(2)
      ID2=IFORD(3)
      GO TO 40
20    ID1=IFORD(1)
      ID2=IFORD(3)
      GO TO 40
30    ID1=IFORD(1)
      ID2=IFORD(2)
40    IF(IABS(ID1).LE.IABS(ID2)) GO TO 50
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
 50     IFLQ=ID1*1000+ID2*100
        IFLQ=ISIGN(IFLQ,IKB)
      GO TO 100
C
C    INITIAL BARION
260   IF(DIQAN) GO TO 400
      CALL FLAVOR(IKB,IFORD(1),IFORD(2),IFORD(3),JSPIN,INDEX)
      CALL FLAVOR(IKA,IFORDP(1),IFORDP(2),IFORDP(3),JSPIN1,INDEX1)
170   INR=INT(1.+3.*RNDM(-1.))
      IFLU=IFORDP(INR)
      IAFLU=IABS(IFLU)
      IQ=0
      DO 160 I=1,3
      IF(IAFLU.NE.IFORD(I)) GO TO 160
      IQ=I
      IFL=IFORD(I)
      GO TO 165
160   CONTINUE
      IF(IQ.EQ.0) GO TO 170
165   IF(INR-2) 175,176,177
175   ID1=IABS(IFORDP(2))
      ID2=IABS(IFORDP(3))
      GO TO 180
176   ID1=IABS(IFORDP(1))
      ID2=IABS(IFORDP(3))
      GO TO 180
177   ID1=IABS(IFORDP(1))
      ID2=IABS(IFORDP(2))
180   IF(ID1.LE.ID2) GO TO 185
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
185   IFQM=ID1*1000+ID2*100
      IFLQM=ISIGN(IFQM,IKA)
      IF(IQ-2) 195,196,197
195   ID1=IFORD(2)
      ID2=IFORD(3)
      GO TO 190
196   ID1=IFORD(1)
      ID2=IFORD(3)
      GO TO 190
197   ID1=IFORD(1)
      ID2=IFORD(2)
190   IF(ID1.LE.ID2) GO TO 194
      ISWAP=ID1
      ID1=ID2
      ID2=ISWAP
194   IFLQ=ID1*1000+ID2*100
      GO TO 100
C     COMPUTE ANNIHILATION DIAGRAM
400   CALL FLAVOR(IKA,IFORDP(1),IFORDP(2),IFORDP(3),JSPIN,INDEX)
      IFT1=IFORDP(1)*1000+IFORDP(2)*100
      IFT2=IFORDP(2)*1000+IFORDP(3)*100
      IFT3=IFORDP(1)*1000+IFORDP(3)*100
410   CALL FLAVOB(IKB,IFL1,IFL2)
      IFLQM=IFORDP(3)
      IF(IABS(IFT1).EQ.IABS(IFL2)) GO TO 420
      IFLQM=IFORDP(1)
      IF(IABS(IFT2).EQ.IABS(IFL2)) GO TO 420
      IFLQM=IFORDP(2)
      IF(IABS(IFT3).EQ.IABS(IFL2)) GO TO 420
      GO TO 410
420   IFLQ=IFL1
      GO TO 100
C---  MESON-MESON COLLISION ---------------------------
101   CONTINUE
      IAN=0
      IF(.NOT.(IFORD(2).GT.0.AND.IFORD(2).EQ.IAFLU)) GO TO 102
      IAN=IAN+1
      IMQ1(IAN)=IFORD(3)
      IMQ2(IAN)=IFLQM
102   IF(.NOT.(IFORD(2).LT.0.AND.IABS(IFORD(2)).EQ.IFLQM))GO TO 103
      IAN=IAN+1
      IMQ2(IAN)=IFLU
      IMQ1(IAN)=IFORD(3)
103   IF(.NOT.(IFORD(3).GT.0.AND.IFORD(3).EQ.IAFLU)) GO TO 104
      IAN=IAN+1
C     IMQ1(IAN)=IFORD(3)
C     IMQ2(IAN)=IFLU
      IMQ1(IAN)=IFORD(2)
      IMQ2(IAN)=IFLQM
104   IF(.NOT.(IFORD(3).LT.0.AND.IABS(IFORD(3)).EQ.IFLQM))GO TO 105
      IAN=IAN+1
C     IMQ2(IAN)=IFLQM
C     IMQ1(IAN)=IFORD(3)
      IMQ2(IAN)=IFLU
      IMQ1(IAN)=IFORD(2)
105   IFLQM=IMQ2(1)
      IFLQ =IMQ1(1)
      IF(IAN.EQ.1) GO TO 100
      IF(RNDM(-1.).GT.0.5) GO TO 100
      IFLQM=IMQ2(2)
      IFLQ =IMQ1(2)
100   NIN=NPTCL+1
C  ******
C     IF(IB(IKA).EQ.0.AND.IB(IKB).EQ.0)
C    *       WRITE(16,*) 'IKA,IKB,IFLQM,IFLQ=',IKA,IKB,IFLQM,IFLQ
C  ******
      AMQ21=0.
      AMQ22=0.
      CALL XCORR(IFLQM,IFLQ,ZER,ZER,ZER,ZER,
     *1.D0,1.D0,PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 1111
      IF((NPTCL-NPOLD).EQ.2.AND.((IDENT(NPTCL).EQ.IDIN(1).AND.
     *IDENT(NPTCL-1).EQ.IDIN(2)).OR.(IDENT(NPTCL).EQ.IDIN(2).AND.
     *IDENT(NPTCL-1).EQ.IDIN(1)))) GO TO 1111
      NFIN=NPTCL
      IRD=0
      DO 300 I=NIN,NFIN
      IORDP(I)=IRD
      IORIG(I)=2
      IF(.NOT.DIQAN) GO TO 300
      IORIG(I)=15
300   CONTINUE
      P0=P0OLD
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE REACTL
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  SELECT REACTION TYPE AT LOW ENERGY IN H-N AND HBAR-N COLLISIONS
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      LOGICAL GH1H2
      COMMON/H1H2/ GH1H2(11)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/HADSIG/SIGS(100),SIGEVT,NSIGS,INOUT(2,100)
      DIMENSION INITYP(2),IREATY(2)
      COMMON/COMANN/ DIQAN
      COMMON/PRIMP0/ P0
      LOGICAL DIQAN
      COMMON/COMMUL/ MULTP
      LOGICAL MULTP
      COMMON /YESELA/YESELA
      LOGICAL YESELA
      COMMON/COMCOL/ NAC(100,4),NBC(100,4),NCOL
      COMMON /PRINTS/ IPRINT
      COMMON /CPRSIG/ ISGCOL
      LOGICAL IPRINT
C     INITIALIZE
  1   SIGEVT=0.
      NCOL=0
      DIQAN=.FALSE.
      IOPAK=10000
      IF(.NOT.MULTP) GO TO 109
      CALL REACTH
      RETURN
109   SIGTOT=0.
      DO 110 I=1,11
110   GH1H2(I)=.FALSE.
      DO 100 ISIGN=1,NSIGS
100   SIGTOT=SIGTOT+SIGS(ISIGN)
      IF(SIGTOT.EQ.0.) GO TO 9999
C    FIND REACTION
      TRY=RNDM(-1.)
      SUM=0.
      DO 200 I=1,NSIGS
      ISIGS=I
      SUM=SUM+SIGS(I)/SIGTOT
      IF(SUM.GT.TRY) GO TO 300
200   CONTINUE
300   SIGEVT=SIGS(ISIGS)
      IF(IPRINT)
     *  WRITE(ITLIS,1001) NSIGS,(SIGS(K),K=1,NSIGS)
1001  FORMAT(1X,'RL:',1X,I3,11(1X,E10.3))
      IF(IPRINT) WRITE(ITLIS,1002) ISIGS,SUM,TRY
1002  FORMAT(1X,'RL: ITLIS,SUM,TRY=',I3,2(1X,E13.6))
C   @@@@@@@@@@@@@@@@@@
C     INITYP(1)=INOUT1(1,ISIGS)
C     INITYP(2)=INOUT2(1,ISIGS)
C     IREATY(1)=INOUT1(2,ISIGS)
C     IREATY(2)=INOUT2(2,ISIGS)
      I1=INOUT(1,ISIGS)
      DO 400 K=1,2
      INITYP(K)=MOD(I1,IOPAK)
  400 I1=I1/IOPAK
      I2=INOUT(2,ISIGS)
      DO 500 K=1,2
      IREATY(K)=MOD(I2,IOPAK)
  500 I2=I2/IOPAK
      IF(IREATY(2).EQ.2.OR.IREATY(2).EQ.4) INITYP(1)=-INITYP(1)
      IF(IREATY(2).EQ.3.OR.IREATY(2).EQ.4) INITYP(2)=-INITYP(2)
C   @@@@@@@@@@@@@@@@@@
      GH1H2(IREATY(1))=.TRUE.
      IF(IPRINT) WRITE(ITLIS,1003) INITYP,IREATY
1003  FORMAT(1X,'RL: INITYP,IREATY=',4(I6,1X))
      IF(YESELA)  RETURN
      IF(GH1H2(4).AND.NSIGS.GT.1)  GO  TO  1
      RETURN
9999  WRITE(ITLIS,1000) SIGEVT
1000  FORMAT(/10X,38H..STOP IN REACTL.. CHECK YOUR INPUT...,
     *7HSIGEVT=,E10.4/)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE TWOSHE(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C       COMPUTE TWO SHEETS ANNIHILATION DIAGRAM
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/PRIMP0/ P0
      COMMON/NEEDR/ NRET
      COMMON/COMASS/ AM1,AM2
      DIMENSION IFORD(3),IFORDP(3)
      DIMENSION PSUM(5)
      LOGICAL RETU
C      INITIALIZE
      IPACK=1000
      NRET=0
      DO 151 I=1,3
151   PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
100   IRET=0
      RETU=.FALSE.
      PSIGN=-1.
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      CALL FLAVOR(IKA,IFORDP(1),IFORDP(2),IFORDP(3),JSPIN,INDEX)
      CALL FLAVOR(IKB,IFORD(1),IFORD(2),IFORD(3),JSPIN,INDEX)
110   INR=INT(1.+3.*RNDM(-1.))
      IFLN=IFORDP(INR)
      IFL1=IABS(IFLN)
      IQ=0
      DO 155 I=1,3
      IF(IFL1.NE.IFORD(I)) GO TO 155
      IQ=I
      GO TO 156
155   CONTINUE
      IF(IQ.EQ.0) GO TO 110
156   IF(IQ-2) 10,20,30
10    IFL01=IFORD(2)
      IFL02=IFORD(3)
      IF(RNDM(-1.).GE.0.5) GO TO 40
      IFL01=IFORD(3)
      IFL02=IFORD(2)
      GO TO 40
20    IFL01=IFORD(1)
      IFL02=IFORD(3)
      IF(RNDM(-1.).GE.0.5) GO TO 40
      IFL01=IFORD(3)
      IFL02=IFORD(1)
      GO TO 40
30    IFL01=IFORD(1)
      IFL02=IFORD(2)
      IF(RNDM(-1.).GE.0.5) GO TO 40
      IFL01=IFORD(2)
      IFL02=IFORD(1)
40    CONTINUE
      IQ=0
      DO 50 I=1,3
      IF(IFL1.NE.IABS(IFORDP(I))) GO TO 50
      IQ=I
      GO TO 60
50    CONTINUE
60    IF(IQ-2) 70,80,90
70    IFL03=IFORDP(2)
      IFL04=IFORDP(3)
      IF(RNDM(-1.).GE.0.5) GO TO 95
      IFL03=IFORDP(3)
      IFL04=IFORDP(2)
      GO TO 95
80    IFL03=IFORDP(1)
      IFL04=IFORDP(3)
      IF(RNDM(-1.).GE.0.5) GO TO 95
      IFL03=IFORDP(3)
      IFL04=IFORDP(1)
      GO TO 95
90    IFL03=IFORDP(1)
      IFL04=IFORDP(2)
      IF(RNDM(-1.).GE.0.5) GO TO 95
      IFL03=IFORDP(2)
      IFL04=IFORDP(1)
95    CONTINUE
C    COMPUTE X VALUES FOR PARTONS
      CALL X2DIST(X1,X2,IFL01,IFL02)
      CALL X2DIST(X3,X4,IFL03,IFL04)
C     PT VALUES FOR PARTONS
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT1,SIGMA)
      AMZER2=AMB**2
      PZER2=P0**2
      AMQ21=AMZER2*(AMZER2+4.*X1*X2*PZER2)/(4.*(AMZER2+PZER2))-PT1**2
      PX1=PT1*COS(PHI)
      PY1=PT1*SIN(PHI)
      PX2=-PX1
      PY2=-PY1
      PHI=twpi*RNDM(-1.)
      CALL GAUSPT(PT3,SIGMA)
      AMZER2=AMA**2
      AMQ22=AMZER2*(AMZER2+4.*X3*X4*PZER2)/(4.*(AMZER2+PZER2))-PT3**2
      PX3=PT3*COS(PHI)
      PY3=PT3*SIN(PHI)
      PX4=-PX3
      PY4=-PY3
      NIN=NPTCL+1
      CALL XCORR(IFL03,IFL01,PX3,PY3,PX1,PY1,X3,X1,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRODS
      IF(RETU) GO TO 100
      CALL XCORR(IFL04,IFL02,PX4,PY4,PX2,PY2,X4,X2,
     *PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 140
      NPTCL=NPTCL-NPRD
      GO TO 100
140   NFIN=NPTCL
      IRD=0
      DO 200 I=NIN,NFIN
      IORDP(I)=IRD
200   IORIG(I)=25
      CALL RESCAL(NIN,NFIN,PSUM,IFAIL)
      IF(IFAIL.EQ.0) RETURN
      NPTCL=NPTCL-NPRD
      GO TO 100
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE THREES(IRET)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     COMPUTE THREE SHEETS ANNIHILATION DIAGRAM
C
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PRIMP0/ P0
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/NEEDR/ NRET
      COMMON/MASQUA/ AMQ21,AMQ22
      COMMON/COAUX1/ AMA2,XI2(20),INUM
      COMMON/COMASS/AM1,AM2
      DIMENSION IFORD(3),IFORD1(3)
      DIMENSION PSUM(5)
      LOGICAL RETU
      EXTERNAL FMQ2
      DATA AM/0./,BM/15.0/,MAXFUN/200/
      DATA EPSI/0.001/
C       INITIALIZE
      IPACK=1000
      NRET=0
      DO 151 I=1,3
151   PSUM(I)=0.
      PSUM(4)=ECM
      PSUM(5)=ECM
150   IRET=0
      RETU=.FALSE.
      PSIGN=-1.
      P02=P0**2
      INUM=0
      IKA=IDIN(1)
      IKB=IDIN(2)
      AMA=AM1
      AMB=AM2
      CALL FLAVOR(IKB,IFORD1(1),IFORD1(2),IFORD1(3),JSPIN,INDEX)
      CALL FLAVOR(IKA,IFORD(1),IFORD(2),IFORD(3),JSPIN,INDEX)
      INR=INT(1.+3.*RNDM(-1.))
      IFL11=IFORD(INR)
      IF(INR-2) 10,20,30
10    IFL22=IFORD(2)
      IFL33=IFORD(3)
      IF(RNDM(-1.).GT.0.5) GO TO 160
      IFL22=IFORD(3)
      IFL33=IFORD(2)
      GO TO 160
20    IFL22=IFORD(1)
      IFL33=IFORD(3)
      IF(RNDM(-1.).GT.0.5) GO TO 160
      IFL22=IFORD(3)
      IFL33=IFORD(1)
      GO TO 160
30    IFL22=IFORD(1)
      IFL33=IFORD(2)
      IF(RNDM(-1.).GT.0.5) GO TO 160
      IFL22=IFORD(2)
      IFL33=IFORD(1)
160   INR=INT(1.+3.*RNDM(-1.))
      IFT11=IFORD1(INR)
      IF(INR-2) 40,50,60
40    IFT22=IFORD1(2)
      IFT33=IFORD1(3)
      IF(RNDM(-1.).GE.0.5) GO TO 260
      IFT22=IFORD1(3)
      IFT33=IFORD1(2)
      GO TO 260
50    IFT22=IFORD1(1)
      IFT33=IFORD1(3)
      IF(RNDM(-1.).GE.0.5) GO TO 260
      IFT22=IFORD1(3)
      IFT33=IFORD1(1)
      GO TO 260
60    IFT22=IFORD1(1)
      IFT33=IFORD1(2)
      IF(RNDM(-1.).GE.0.5) GO TO 260
      IFT22=IFORD1(2)
      IFT33=IFORD1(1)
C    COMPUTE X VALUES FOR PARTONS
260   CALL X3DIST(X11,X12,X13,IFL11,IFL22,IFL33)
      AMA2=AMA**2
      INUM=3
      XI2(1)=X11**2
      XI2(2)=X12**2
      XI2(3)=X13**2
      CALL RZERO(AM,BM,AMQ1P,RES,EPSI,MAXFUN,FMQ2)
      IF(RES.LT.0.) GO TO 260
C  COMPUTE PT VALUES FOR PARTONS
      CALL GAUSPT(PT11,SIGMA)
      PHI1=twpi*RNDM(-1.)
      CALL GAUSPT(PT12,SIGMA)
      PHI2=twpi*RNDM(-1.)
      PT11X=PT11*COS(PHI1)
      PT11Y=PT11*SIN(PHI1)
      PT12X=PT12*COS(PHI2)
      PT12Y=PT12*SIN(PHI2)
      PT13X=-(PT11X+PT12X)
      PT13Y=-(PT11Y+PT12Y)
C  COMPUTE X VALUES FOR PARTONS
261   CALL X3DIST(X21,X22,X23,IFT11,IFT22,IFT33)
      AMA2=AMB**2
      INUM=3
      XI2(1)=X21**2
      XI2(2)=X22**2
      XI2(3)=X23**2
      CALL RZERO(AM,BM,AMQ2P,RES,EPSI,MAXFUN,FMQ2)
      IF(RES.LT.0.) GO TO 261
C  COMPUTE PT VALUES FOR PARTONS
      CALL GAUSPT(PT21,SIGMA)
      PHI1=twpi*RNDM(-1.)
      CALL GAUSPT(PT22,SIGMA)
      PHI2=twpi*RNDM(-1.)
      PT21X=PT21*COS(PHI1)
      PT21Y=PT21*SIN(PHI1)
      PT22X=PT22*COS(PHI2)
      PT22Y=PT22*SIN(PHI2)
      PT23X=-(PT21X+PT22X)
      PT23Y=-(PT21Y+PT22Y)
      NIN1=NPTCL+1
      AMQ21=AMQ1P-PT11**2
      AMQ22=AMQ2P-PT21**2
      CALL XCORR(IFL11,IFT11,PT11X,PT11Y,PT21X,PT21Y,
     *X11,X21,PSIGN,NPRODS,RETU)
      IF(RETU) GO TO 150
      NPRD=NPRODS
      AMQ21=AMQ1P-PT12**2
      AMQ22=AMQ2P-PT22**2
      CALL XCORR(IFL22,IFT22,PT12X,PT12Y,PT22X,PT22Y,
     *X12,X22,PSIGN,NPRODS,RETU)
      NPRD=NPRD+NPRODS
      IF(.NOT.RETU) GO TO 130
      NPTCL=NPTCL-NPRD
      GO TO 150
130   AMQ21=AMQ1P-(PT13X**2+PT13Y**2)
      AMQ22=AMQ2P-(PT23X**2+PT23Y**2)
      CALL XCORR(IFL33,IFT33,PT13X,PT13Y,PT23X,PT23Y,
     *X13,X23,PSIGN,NPRODS,RETU)
      IF(.NOT.RETU) GO TO 140
      NPTCL=NPTCL-NPRD
      GO TO 150
140   NFIN1=NPTCL
      DO 500 I=NIN1,NFIN1
      IORDP(I)=0
500   IORIG(I)=35
      NPRD=NPRD+NPRODS
      CALL RESCAL(NIN1,NFIN1,PSUM,IFAIL)
      IF(IFAIL.EQ.0) GO TO 501
      NPTCL=NPTCL-NPRD
      GO TO 150
501   RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BINAR(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  SIMULATION TWO PARTICLE REACTION
C
      COMMON /COMKI1/ HLA2,HLB2,W,INUMA
      COMMON/COMLID/PLIDER(499)
      COMMON /COMKI2/ELA,ELB,PLALB
      COMMON /CALC/HA,HB,HA2,HB2
      COMMON /BEES/B,BFOR
      COMMON/COMASS/ AM1,AM2
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      DIMENSION PA(3)
      LOGICAL SPINT
C  INITIALIZE
      INUMA=1
      NREP =0
      IK1=IDIN(1)
      IK2=IDIN(2)
      SPINT=.FALSE.
      IRET=0
      IB1=IB(IK1)
      IB2=IB(IK2)
       W= ECM
      HLA=AM1
      HLB=AM2
      PARBE=0.2
      IF(ECM.GT.HLA+HLB+PARBE) GO TO 999
      IRET=1
      RETURN
999   HLA2=HLA*HLA
      HLB2=HLB*HLB
      PLALB=SQRT(ALAMB(SCM,HLA2,HLB2))/(2.0*ECM)
      ELA=(SCM+HLA2-HLB2)/(2.0*ECM)
      ELB=(SCM+HLB2-HLA2)/(2.0*ECM)
C   SELECT INTERACTIVE QUARKS
105   CALL FLAVOB(IK1,IFL1,IFL2)
      CALL FLAVOB(IK2,IFL3,IFL4)
      NREP=NREP+1
      IF(NREP.LT.NTRIES) GO TO 106
C     WRITE(ITLIS,1200) IDIN(1),IDIN(2),PLAB
1200  FORMAT(1X,'IN BINAR:NREP > NTRIES,IK1,IK2,PLAB=',
     *2I4,1X,F7.3)
      IRET=1
      RETURN
106   CONTINUE
      IREP1=0
C  HADRONS GENERATE BY MEANS QUARKS EXCHANGE
      IF(IFL1)1,1,2
1     IF11=IFL2
      IF22=IFL1
      GO TO 3
2     IF11=IFL1
      IF22=IFL2
3     CONTINUE
      IF(IB2.EQ.1) GO TO 102
      IF(IFL3) 101,101,102
101   IF44=IFL3
      IF33=IFL4
      GO TO 103
102   IF33=IFL3
      IF44=IFL4
103   CONTINUE
104   IKH2=IDPARS(IF11,IF44,SPINT,0)
      IKH1=IDPARS(IF33,IF22,SPINT,0)
      IREP1=IREP1+1
      IF(IREP1.GT.NTRIES) GO TO 105
C  SELECT ELASTIC COLLISION
C     IF(IKH1.EQ.IK1.AND.IKH2.EQ.IK2.AND.IREP1.LE.NTRIES) GO TO 104
      IF(IKH1.EQ.IK1.AND.IKH2.EQ.IK2.AND.IREP1.LE.NTRIES) GO TO 105
C  SELECT TABLE MASSES AND TABLE WIDTH OF HADRONS
      AMH1=AMASSF(IKH1)
      AMH2=AMASSF(IKH2)
      GAMH1=GAMHE(IKH1)
      GAMH2=GAMHE(IKH2)
C  COMPUTE MASSES OF PARTICLES
      IREP3=0
205   CONTINUE
      IREP3=IREP3+1
C     IF(IREP3.GT.NTRIES) GO TO 104
      IF(IREP3.GT.NTRIES) GO TO 105
      GAM1=WIDTH(GAMH1)
      GAM2=WIDTH(GAMH2)
      AMP1=AMH1+GAM1
      AMP2=AMH2+GAM2
C  CHECK ENERGY THRESHOLD
      IF(W.LT.AMP1+AMP2) GO TO 205
      HA=AMP1
      HA2=AMP1**2
      HB=AMP2
      HB2=AMP2**2
C  COMPUTE SCATTERING ANGLE
      IF(IB1.EQ.0.OR.IB2.EQ.0) INUMA=0
      IF(IB1.EQ.0.AND.IB2.EQ.0) INUMA=2
      CALL SLOPEB(IB1,IB2,PLALB,B)
      CALL ANG(TFOR,TBACK,T,Z,PHI)
      PAMOD=SQRT(ALAMB(SCM,HA2,HB2))/(2.0*ECM)
      PAN=PAMOD*SQRT(1.-Z**2)
      PA(1)=PAN*COS(PHI)
      PA(2)=PAN*SIN(PHI)
      PA(3)=PAMOD*Z
      NPTCL=NPTCL+1
      IDENT(NPTCL)=IKH1
      PPTCL(1,NPTCL)=PA(1)
      PPTCL(2,NPTCL)=PA(2)
      PPTCL(3,NPTCL)=PA(3)
      PPTCL(4,NPTCL)=SQRT(PAMOD**2+HA2)
      PPTCL(5,NPTCL)=AMP1
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IDCAY(NPTCL)=0
      IORIG(NPTCL)=8
      IORDP(NPTCL)=0
      NPTCL=NPTCL+1
      IDENT(NPTCL)=IKH2
      PPTCL(1,NPTCL)=-PA(1)
      PPTCL(2,NPTCL)=-PA(2)
      PPTCL(3,NPTCL)=-PA(3)
      PPTCL(4,NPTCL)=SQRT(PAMOD**2+HB2)
      PPTCL(5,NPTCL)=AMP2
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      IORIG(NPTCL)=8
      IDCAY(NPTCL)=0
      IORDP(NPTCL)=0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION WIDTH(GAM)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE WIDTH OF PARTICLE
C
100   DRND=RNDM(-1.)
      GT=GAM*SQRT(-LOG(DRND))
      PHI=twpi*RNDM(-1.)
      WIDTH=GT*COS(PHI)
      IF(ABS(WIDTH).GT.GAM) GO TO 100
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION COSAN(T)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C    COMPUTE ANGULAR DISTRIBUTION
C
C    FOR AP P COLLISION
C
      COMMON/TABLE1/ B12(30,3),ENER(30)
      LOGICAL SWANG
      RAN=RNDM(-1.)
      IF(T.LE.0.084) GO TO 200
      AMP=0.939
      PP=APP(T,AMP)
      WCOS=0.0105+0.0165*PP
      IF(T.GE.59.0) WCOS=1.
      SWANG=.FALSE.
      IF(RNDM(-1.).LE.WCOS) SWANG=.TRUE.
      IF(T.LE.0.226) GO TO 300
      IF(SWANG) GO TO 100
      B1=ANGINT(T,1)
      CANG0=ANGINT(T,3)
      COSAN=LOG(RAN*EXP(B1)+(1.-RAN)*EXP(B1*CANG0))/B1
      RETURN
100   B2=ANGINT(T,2)
      CANG0=ANGINT(T,3)
      COSAN=LOG(RAN*EXP(B2*CANG0)+(1.-RAN)*EXP(-B2))/B2
      RETURN
200   B1=ANGINT(T,1)
      COSAN=1.+LOG(RAN*(1.-EXP(-2.*B1))+EXP(-2.*B1))/B1
      RETURN
300   IF(SWANG) GO TO 400
      B1=ANGINT(T,1)
      CANG0=ANGINT(T,3)
      COSAN=LOG(RAN*EXP(B1)+(1.-RAN)*EXP(B1*CANG0))/B1
      RETURN
400   CANG0=ANGINT(T,3)
      COSAN=RAN*(CANG0+1.)-1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION APP(T,W)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C     CALCULATION OF THE MOMENTUM
      APP=SQRT(T*(T+2.*W))
      RETURN
               END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION ANGINT(T,I)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C    COMPUTE CONSTANTS BY INTERPOLATION
C    FOR AP P ELASTIC SCATTERING
C
      COMMON/TABLE1/ B12(30,3),ENER(30)
      COMMON/INTERP/ F(6)
C   COMPUTE OF THREE POINTS FOR INTERPOLATION
      L=1
4     IF(T-ENER(L)) 6,5,11
5     ANGINT=B12(L,I)
      RETURN
6     IF(L.GT.1) GO TO 7
      ANGINT=B12(L,I)
      RETURN
7     IF(L.GE.29) GO TO 8
      L1=L-2
      GO TO 9
8     L1=27
9     DO 10 K=1,3
      LL1=L1+K
      F(K)=B12(LL1,I)
      F(K+3)=ENER(LL1)
10    CONTINUE
      GO TO 12
11    L=L+1
      IF(L.EQ.30) GO TO 13
      GO TO 4
12    ANGINT=SINTER(T)
      RETURN
13    ANGINT=B12(30,I)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ELASTQ(IRET)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C     MONTE CARLO SIMULATION ELASTIC HADRON NUCLEON COLLISION
C
      COMMON/COMKI1/ HLA2,HLB2,W,INUMA
      COMMON/COMKI2/ ELA,ELB,PLALB
      COMMON/CALC/ HA,HB,HA2,HB2
      COMMON/BEES/ B,BFOR
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON/PARORD/ IORDP(499)
      COMMON /ORDER/ IRD1,IRD2
      COMMON/COMLID/ PLIDER(499)
      COMMON/COMELX/ SIGEL
      COMMON/COMCRO/ SIGTOT
      COMMON/COMASS/ AM1,AM2
      COMMON/NPTCLZ/ NPTCLZ
      COMMON/PRIMP0/ P0
      DIMENSION PA(3),P1(3),P2(3)
      IRET=0
      SIGEL0=SIGEL
      P02=P0**2
      IEXE=0
      IK01=IDIN(1)
      IK02=IDIN(2)
      HLA=AM1
      HLB=AM2
      HLA2=HLA*HLA
      HLB2=HLB*HLB
C   W= CENTRE OF MASS (C.M.) ENERGY
      W=ECM
C  TKIN=KINETIC ENERGY OF PROJECTILE IN TARGET REST FRAME
      TKIN=(SCM-HLA2-HLB2)/(2.0*HLB)-HLA
C   PLALB=CM MOMENTUM OF A OR B IN ELASTIC EVENTS
      PLALB=SQRT(ALAMB(SCM,HLA2,HLB2))/(2.0*W)
C   ELA=CM ENERGY OF A IN ELASTIC EVENT *** ELB=SAME FOR B
      ELA=(SCM+HLA2-HLB2)/(2.0*W)
      ELB=(SCM+HLB2-HLA2)/(2.0*W)
      IK1=IK01
      IK2=IK02
      IB1=IB(IK1)
      IB2=IB(IK2)
      HA=HLA
      HB=HLB
      INUMA=1
      IF(IB1.EQ.0.OR.IB2.EQ.0) INUMA=0
      IF(IB1.EQ.0.AND.IB2.EQ.0) INUMA=2
      HA2=HA*HA
      HB2=HB*HB
      TOBR=10.0
      IF(IB1.NE.0) GOTO 71
        TOBR=2.4
      IF(IB2.NE.0) GO TO 71
      TOBR=0.
      IF(TKIN.GT.2.5) GO TO 71
      IF(AM1.GT.0.50.OR.AM2.GT.0.50) GO TO 71
C  COMPUTE OF RESONANCE PARAMETERS
C   PI+PI--RHO,OMEGA,  PI+K--K*,   K+K--PHI
C
      QSUM = CHARGE(IK01)+CHARGE(IK02)
      IF(ABS(QSUM).GE.2.) GO TO 71
      ISOB=0
      P1(1)=0.
      P1(2)=0.
      P1(3)=P0
      P2(1)=0.
      P2(2)=0.
      P2(3)=-P0
      CALL FOROM(IK1,P1,AM1,IK2,P2,AM2,SIGEL0,
     *IKD,PXD,PYD,PZD,DMAS,ISOB)
      IF(ISOB.EQ.0) GO TO 71
      NPTCL=1
      IORDP(NPTCL)=0
      IORIG(NPTCL)=0
      IDCAY(NPTCL)=0
      IDENT(NPTCL)=IKD
      PPTCL(5,NPTCL)=DMAS
      PPTCL(1,NPTCL)=PXD
      PPTCL(2,NPTCL)=PYD
      PPTCL(3,NPTCL)=PZD
      PPTCL(4,NPTCL)=SQRT(DMAS**2+PXD**2+PYD**2+PZD**2)
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      RETURN
 71   IF(TKIN-TOBR) 72,72,73
 72   CALL ELZPHI(IK01,IK02,TKIN,Z,PHI,IEXE)
      GO TO 74
 73   CONTINUE
C   SLOPE CALCULATES THE DIFFRACTIVE SLOPES FOR THE CHOSEN MASSES
      CALL SLOPE(B,BFOR)
C   ANG CALCULATES THE TWO-BODY SCATTERING ANGLES (AZIMUTHAL ANGLE PHI
C   AND POLAR ANGLE THETA,WHERE Z=COS(THETA)
      IB1=IB(IK1)
      IB2=IB(IK2)
C     IF(IB1.EQ.-1.AND.IB2.NE.-1) B=11.0
      CALL ANG(TFOR,TBACK,T,Z,PHI)
74    IF(IEXE.EQ.0) GO TO 76
      HA=AMASS(IK1)
      HB=AMASS(IK2)
      HA2=HA*HA
      HB2=HB*HB
 76   PAMOD=SQRT(ALAMB(SCM,HA2,HB2))/(2.0*W)
      PAN=PAMOD*SQRT(ABS(1.-Z**2))
      PA(1)=PAN*COS(PHI)
      PA(2)=PAN*SIN(PHI)
      PA(3)=PAMOD*Z
      EA=SQRT(PAMOD**2+HA2)
      EB=SQRT(PAMOD**2+HB2)
      NPTCL=NPTCL+1
      IDCAY(NPTCL)=0
      IORIG(NPTCL)=4
      IORDP(NPTCL)=IRD1
      IDENT(NPTCL)=IK1
      PPTCL(1,NPTCL)=PA(1)
      PPTCL(2,NPTCL)=PA(2)
      PPTCL(3,NPTCL)=PA(3)
      PPTCL(4,NPTCL)=EA
      PPTCL(5,NPTCL)=HA
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      NPTCL=NPTCL+1
      IDCAY(NPTCL)=0
      IORIG(NPTCL)=4
      IORDP(NPTCL)=IRD2
      IDENT(NPTCL)=IK2
      PPTCL(1,NPTCL)=-PA(1)
      PPTCL(2,NPTCL)=-PA(2)
      PPTCL(3,NPTCL)=-PA(3)
      PPTCL(4,NPTCL)=EB
      PPTCL(5,NPTCL)=HB
      PPTCL(6,NPTCL)=0.
      PPTCL(7,NPTCL)=0.
      PPTCL(8,NPTCL)=0.
      PPTCL(9,NPTCL)=0.
      PLIDER(NPTCL)=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        SUBROUTINE ELZPHI(IK1,IK2,TKIN,Z,PHI,IEXE)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      PHI=twpi*RNDM(-1.)
      IF(IEXE.EQ.0) GO TO 1
       Z=COSP(TKIN,12)
       RETURN
1     IB1=IB(IK1)
      IB2=IB(IK2)
      IF(IB1.LT.0.AND.IB2.GT.0) GO TO 2
      CALL MARK(IK1,IK2,KS)
      IBP=IB(IK1)
      Z=COSAM(IBP,TKIN,KS)
      RETURN
2     Z=COSAN(TKIN)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION FMQ2(AMQ2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/PRIMP0/ P0
      COMMON/COAUX1/ AMA2,XI2(20),INUM
      P02=P0**2
      ENER=SQRT(AMA2+P02)
      SUME=0.
      DO 100 I=1,INUM
      SUME=SUME+SQRT(AMQ2+XI2(I)*P02)
100   CONTINUE
      FMQ2=ENER-SUME
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE REACTH
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  SELECT REACTION TYPE AT HIGH ENERGY IN H-N AND HBAR-N COLLISIONS
C
      LOGICAL GH1H2
      COMMON/H1H2/ GH1H2(11)
      COMMON/ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/HADSIG/SIGS(100),SIGEVT,NSIGS,INOUT(2,100)
      COMMON/COMPOM/ POMGEN(15)
      COMMON/COMCOL/ NAC(100,4),NBC(100,4),NCOL
      DIMENSION INITYP(2),IREATY(2)
      COMMON/COMPLI/ LIMP
      COMMON/YESELA/YESELA
      COMMON /PRINTS/ IPRINT
      COMMON /CPRSIG/ ISGCOL
      LOGICAL YESELA
      LOGICAL IPRINT
C      INITIALIZE
  1   SIGEVT=0.
      IOPAK=10000
      SIGTOT=0.
      DO 110 I=1,11
110   GH1H2(I)=.FALSE.
      DO 100 ISIGN=1,NSIGS
100   SIGTOT=SIGTOT+SIGS(ISIGN)
      IF(SIGTOT.EQ.0.) GO TO 9999
C   FIND REACTION
      TRY=RNDM(-1.)
      SUM=0.
      DO 200 I=1,NSIGS
      ISIGS=I
      SUM=SUM+SIGS(I)/SIGTOT
      IF(SUM.GT.TRY) GO TO 300
200   CONTINUE
300   SIGEVT=SIGS(ISIGS)
      IF(IPRINT)
     *  WRITE(ITLIS,1001) NSIGS,(SIGS(K),K=1,NSIGS)
1001  FORMAT(1X,'RH:',1X,I3,11(1X,E10.3))
      IF(IPRINT) WRITE(ITLIS,1002) ISIGS,SUM,TRY
1002  FORMAT(1X,'RH: ITLIS,SUM,TRY=',I3,2(1X,E13.6))
C   @@@@@@@@@@@@@@@@@
C     INITYP(1)=INOUT1(1,ISIGS)
C     INITYP(2)=INOUT2(1,ISIGS)
C     IREATY(1)=INOUT1(2,ISIGS)
C     IREATY(2)=INOUT2(2,ISIGS)
      I1=INOUT(1,ISIGS)
      DO 400 K=1,2
      INITYP(K)=MOD(I1,IOPAK)
  400 I1=I1/IOPAK
      I2=INOUT(2,ISIGS)
      DO 500 K=1,2
      IREATY(K)=MOD(I2,IOPAK)
  500 I2=I2/IOPAK
      IF(IREATY(2).EQ.2.OR.IREATY(2).EQ.4) INITYP(1)=-INITYP(1)
      IF(IREATY(2).EQ.3.OR.IREATY(2).EQ.4) INITYP(2)=-INITYP(2)
C   @@@@@@@@@@@@@@@@@
      GH1H2(IREATY(1))=.TRUE.
      IF(IPRINT) WRITE(ITLIS,1003) INITYP,IREATY
1003  FORMAT(1X,'RH: INITYP,IREATY=',4(I6,1X))
      IF(IREATY(1).NE.7) GO TO 910
C  SELECT NUMBER OF POMERONS
      TRY=RNDM(-1.)
      DO 700 NPOM=1,LIMP
      NC=NPOM
      IF(POMGEN(NPOM).GT.TRY) GO TO 800
700   CONTINUE
800   CONTINUE
      NCOL=NC
      DO 900 J=1,NCOL
      NAC(J,3)=0
      NBC(J,3)=0
      NAC(J,1)=1
      NBC(J,1)=1
      NAC(J,4)=0
      NBC(J,4)=0
      NAC(J,2)=INITYP(1)
      NBC(J,2)=INITYP(2)
900   CONTINUE
910   CONTINUE
      IF(YESELA)  RETURN
      IF(GH1H2(4))  GO  TO  1
      RETURN
9999  WRITE(ITLIS,1000) SIGEVT
1000  FORMAT(//10X,28H...CHECK YOUR INPUT..SIGEVT=,E10.4)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION FP02(PNEW2)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      COMMON/PRIMP0/ P0
      COMMON/COAUX1/ AMA2,XI2(20),INUM
      COMMON/COAUX2/ X12,X22,X32,PT12,PT22,PT32
      P02=P0**2
      ENER=SQRT(P02+AMA2)
      SUME=SQRT(AMA2+X32*PNEW2+PT32)+SQRT(X22*PNEW2+PT22)+
     *SQRT(X12*PNEW2+PT12)
      FP02=ENER-SUME
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION RNFAC(K)
C
C      RETURN N FACTORIAL
C
      use, intrinsic:: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), intent(in) :: k

      integer(int32) :: j
      IF(K.GT.1) GO TO 1
      RNFAC=1.
      RETURN
1     RNFAC=DBLE(K)
      DO 2 J=2,K
2     RNFAC=RNFAC*DBLE(K-J+1)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION FUNIT(Z)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C    RETURN UNITARIZATION FUNCTION FOR POMERON
C
      DATA LIMP2/12/
      FUNIT=1.
      DO 100 J=2,LIMP2
100   FUNIT=FUNIT+(-Z)**(J-1)/(DBLE(J)*RNFAC(J))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SIGIN
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  COMPUTE DIAGRAM WEIGHTS FOR INELASTIC
C  HADRON BARYON COLLISIONS
C
C  SIGMA    =CROSS SETION SUMMED OVER TYPES ALLOWED BY GH1H2
C  SIGS(I)  =PARTIAL CROSS SECTION FOR DUAL TYPE DIAGRAM
C
C
      CHARACTER*8 LAB1,LAB2
      LOGICAL GH1H2
      LOGICAL GH
      LOGICAL IPRINT
      LOGICAL MULTP
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/HADSIG/SIGS(100),SIGEVT,NSIGS,INOUT(2,100)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/REACOE/ COEF(11),COEF1(11)
      COMMON/COMPLI/ LIMP
      COMMON/COMECB/ ECMB
      COMMON/H1H2/ GH1H2(11)
      COMMON/COMPOM/ POMGEN(15)
      COMMON/CANPOM/ POAGEN(15)
      COMMON/CSIGSU/ SIGSUM
      COMMON/PRIMPL/ PL
      COMMON/CSIGA/ SIGAN
      COMMON/COMELX/ SIGEL
      COMMON/COMCRO/ SIGTOT
      COMMON/COMASS/ AM1,AM2
      COMMON/CPRSIG/ ISGCOL
      COMMON/PRINTS/IPRINT
      COMMON/COMENB/ ENBOU
      COMMON /SIGDIA/ CROSD(5),DST
      COMMON/FRGCPA/ PUDC,PUDCC,PSPINC,PJSPNC,PMIX1C(3,2),PMIX2C(3,2),
     *PBARC
      COMMON/FRGSPA/ PUDS,PUDSC,PSPINS,PJSPNS,PMIX1S(3,2),PMIX2S(3,2),
     *SIGQTS,WENDM,WENDB,PBARS,PRIQS(9),PARDBS,PARQLS,PARRS,PUNDS
      COMMON/COMMUL/ MULTP
*@@@@@@@@@@
      COMMON/FLACOM/NFLA,NFL1,NFL2,NFL3,NSPIN,NNDEX
*@@@@@@@@@@
      DIMENSION GH(11)
      DIMENSION COEF01(11)
      DIMENSION IFL1(3),IFL2(3)
C
      DATA ALFR/0.5/
      DATA BPOM/2.05/,BPOMA/2.0/,APOM/3.5/,APOMA/3.5/
C
      DO 1 I=1,11
 1    COEF01(I)=COEF1(I)
C
C   INITIALIZE CROSS SECTION
C             SCALE=0.389
      ISGCOL=0
      SIGSUM=0.
      NSIGS=0
      SIGM =0.
      SIGPS=0.
      SIG=0.
      IPL=0
       PJSC=PJSPNC
       PSPC=PSPINC
       PJSS=PJSPNS
       PSPS=PSPINS
      DO 99 I=1,11
99    GH(I)=.TRUE.
*@@@@@@@@@@
      NFLA=-1
*@@@@@@@@@@
      CALL FLAVOR(IDIN(1),IFL1(1),IFL1(2),IFL1(3),JSPIN1,INDEX1)
      CALL FLAVOR(IDIN(2),IFL2(1),IFL2(2),IFL2(3),JSPIN2,INDEX2)
*@@@@@@@@@@
      IF(NFLA.EQ.-1)  NFLA=-2
*@@@@@@@@@@
      IB1=IB(IDIN(1))
      IB2=IB(IDIN(2))
      IF(IB1.GT.0.AND.IB2.GT.0) GO TO 5
      IF(IB1.LT.0.AND.IB2.LT.0) GO TO 5
      L1=1
      IF(IB1.EQ.0.AND.IB2.EQ.0) L1=2
      DO 3 I=L1,3
      I1=IFL1(I)
      DO 3 J=L1,3
      I2=-IFL2(J)
      IF(I1.NE.I2) GO TO 39
      IPL=1
       GO TO 4
 39   CONTINUE
 3    CONTINUE
 4    IF(IPL.EQ.1) GO TO 77
 5    GH(2) = .FALSE.
      GH1H2(2)=.FALSE.
77    IF(IB1.LT.0.AND.IB2.GT.0) GO TO 6
      GH(5)=.FALSE.
      GH1H2(5)=.FALSE.
      GO TO 7
  6   GH(8)=.FALSE.
      GH1H2(8)=.FALSE.
 7    CONTINUE
      IF(IB1.EQ.0.AND.IB2.EQ.0) CALL  GHGH00(GH)
 8    PARBE=0.38
      DSM=ECM-AM1-AM2
      IF(ECM.GE.ENBOU) GO TO 33
      PJSPNC=.75
      PSPINC=.75
      PJSPNS=.75
      PSPINS=.75
      GH1H2(3)=.FALSE.
      GH(3)=.FALSE.
      GH1H2(11)=.FALSE.
      GH(11)=.FALSE.
C ONLY PLAN,DIFSMA,DIFTRI,ANNIH ALLOWED FOR ECM < ENBOW
      MULTP=.FALSE.
      IF(IB1.EQ.0.AND.IB2.LT.0) GO TO 30
      IF(IB1.LE.0.AND.IB2.LE.0) GO TO 30
      GH1H2(7)=.FALSE.
      GH1H2(8)=.FALSE.
      GH1H2(10)=.FALSE.
      GH(7)=.FALSE.
      GH(8)=.FALSE.
      GH(10)=.FALSE.
      GO TO 33
 30   IF(DSM.GE.PARBE) GO TO 33
      GH1H2(1)=.FALSE.
C    **** PLANAR DIAG. SHOULD NOT BE FORBIDDEN  ***  SIVOKL.01.08.91
      GH1H2(2)=.FALSE.
      GH1H2(6)=.FALSE.
      GH1H2(7)=.FALSE.
C    **** BINAR  DIAG. SHOULD NOT BE FORBIDDEN  ***  SIVOKL.01.08.91
      GH1H2(8)=.FALSE.
      GH1H2(10)=.FALSE.
      GH(1)=.FALSE.
      GH(2)=.FALSE.
      GH(6)=.FALSE.
      GH(7)=.FALSE.
      GH(8)=.FALSE.
      GH(10)=.FALSE.
33    I1=IABS(IDIN(1))
      I2=IABS(IDIN(2))
      I4=1
C   @@@@@@@@@@@@@@@
      IF(IDIN(1).LT.0) I4=2
      IF(IDIN(2).LT.0) I4=3
      IF(IDIN(1).LT.0.AND.IDIN(2).LT.0) I4=4
C   @@@@@@@@@@@@@@@
      IF(IB1.EQ.0.AND.IB2.EQ.0) MULTP=.FALSE.
      CALL LABEL(LAB1,IDIN(1))
      CALL LABEL(LAB2,IDIN(2))
C     IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,1000) LAB1,LAB2,SCM
1000  FORMAT(//15X,47HI SELECT THE NEXT REACTIONS WITH CROSS SECTIONS
     *,/20X,3HFOR,1X,A8,A8,18H COLLISION AT SCM=,E10.4,7H GEV**2/)
      IFL=IDIN(1)/1000
      IF(IFL.NE.0) GO TO 107
      JFL=MOD(I1/100,10)
      KFL=MOD(I1/10,10)
      IF(JFL.NE.3.AND.KFL.NE.3) GO TO 106
       COEF01(2)= COEF1(2)*0.4
       COEF01(3)= COEF1(3)*0.33
       COEF01(8)= COEF1(8)*0.33
       COEF01(10)= COEF1(10)*0.33
      GO TO 107
106   CONTINUE
       COEF01(2)= COEF1(2)*0.4
       COEF01(3)= COEF1(3)*0.33
       COEF01(8)= COEF1(8)*0.33
       COEF01(10)= COEF1(10)*0.33
107   CONTINUE
      IF(ECM.GE.ECMB) GO TO 400
C  TWO PARTICLE REACTION CROSS SECTION
C
      SIG     = COEF01(8)/SCM
      IF(ECM.LT.ENBOU) SIG=CROSD(1)
      IF(GH1H2(8)) GO TO 155
      IF(.NOT.GH(8)) GO TO 200
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,80) SIG
80    FORMAT(10X,34HI DO NOT USE BINAR REACTION.SIGMA=,E10.4)
       GO TO 200
155   CALL SIGFIH(SIG,I1,I2,8,I4)
200   CONTINUE
C    PLANAR TYPE DIAGRAM (Q,QBAR-ANNIHILATION)
      SIG     = COEF01(2)*SCM**(ALFR-1.)
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 292
      IF(ECM.LT.ENBOU.AND.IB1*IB2.EQ.-1) SIG = CROSD(3)
      IF(ECM.LT.ENBOU.AND.IB1*IB2.NE.-1) SIG = CROSD(4)
      IF(ECM.LT.ENBOU.AND..NOT.GH1H2(2)) GH(2)=.TRUE.
292   IF(GH1H2(2)) GO TO 291
      IF(.NOT.GH(2)) GO TO 400
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,82) SIG
 82   FORMAT(10X,35HI DO NOT USE PLANAR REACTION.SIGMA=,E10.4)
      GO TO 400
291   CALL SIGFIH(SIG,I1,I2,2,I4)
C      ELASTIC SCATTERING DIAGRAM
400   SIG=SIGEL
      IF(GH1H2(4)) GO TO 491
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,84) SIG
 84   FORMAT(10X,36HI DO NOT USE ELASTIC REACTION.SIGMA=,E10.4)
      GO TO 500
491   CALL SIGFIH(SIG,I1,I2,4,I4)
500   CONTINUE
      IF(ECM.GT.ECMB) GO TO 600
C     ANNIHILATION DIAGRAM
      SIG=SIGAN
      IF(GH1H2(5)) GO TO 591
      IF(.NOT.GH(5)) GO TO 600
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,85) SIG
 85   FORMAT(10X,41HI DO NOT USE ANNIHILATION REACTION.SIGMA=,E10.4)
      GO TO 600
591   CALL SIGFIH(SIG,I1,I2,5,I4)
      BPOMA=2.35-0.25*LOG(ECM)
      SIGPA= 0.
      DO 605 I=1,LIMP
      POAGEN(I)=(1.+APOMA*I**2)*EXP(-BPOMA*I)
 605  SIGPA=SIGPA+POAGEN(I)
      PSUM=0.
      DO 606 I=1,LIMP
      PSUM=PSUM+POAGEN(I)/SIGPA
      POAGEN(I)=PSUM
 606  CONTINUE
C   HIGH MASS DIFFRACTION
600   SIG = COEF01(1)*SIGEL
      IF(ECM.LT.ENBOU.AND.IB1*IB2.EQ.-1) SIG = CROSD(2)
      IF(ECM.LT.ENBOU.AND.IB1*IB2.NE.-1) SIG = CROSD(2)*DST
      IF(GH1H2(1)) GO TO 160
      IF(.NOT.GH(1)) GO TO 802
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,81) SIG
81    FORMAT(10X,37HI DO NOT USE DIFR.OF HIGH MASS.SIGMA=,E10.4)
      GO TO 802
160   CALL SIGFIH(SIG,I1,I2,1,I4)
C   TRIPLE REGGEON DIAGRAM
802   SIG = COEF01(10)*SCM**(ALFR-1.)
      IF(ECM.LT.ENBOU) SIG = CROSD(3)
      IF(GH1H2(10)) GO TO 162
      IF(.NOT.GH(10)) GO TO 800
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,181) SIG
181   FORMAT(10X,39HI DO NOT USE TRIPLE REGGEON DIAG.SIGMA=,E10.4)
      GO TO 800
162   CALL SIGFIH(SIG,I1,I2,10,I4)
C   DOUBLE DIFRACTION DIAGRAM
C
800    CONTINUE
      SIG =SIGEL* COEF01(11)
      IF(GH1H2(11)) GO TO 172
      IF(.NOT.GH(11)) GO TO 890
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,191) SIG
191   FORMAT(10X,42HI DO NOT USE DOUBLE DIFRACTION DIAG.SIGMA=,E10.4)
      GO TO 890
172   CALL SIGFIH(SIG,I1,I2,11,I4)
C    SMALL MASS DIFFRACTION
890    CONTINUE
      SIG =SIGEL* COEF01(6)
      IF(IB1.EQ.0.AND.IB2.EQ.0) GO TO 692
      IF(ECM.LT.ENBOU.AND.IB1*IB2.EQ.-1) SIG = CROSD(1)
      IF(ECM.LT.ENBOU.AND.IB1*IB2.NE.-1) SIG = CROSD(2)*(1.-DST)
      GO TO 693
692   SIG=SIG*1.9
693   IF(GH1H2(6)) GO TO 691
      IF(.NOT.GH(6)) GO TO 300
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,86) SIG
 86   FORMAT(10X,39HI DO NOT USE DIFR. OF SMALL MASS.SIGMA=,E10.4)
      GO TO 300
691   CALL SIGFIH(SIG,I1,I2,6,I4)
300   IF(ECM.GT.ECMB) GO TO 700
C           UNDEVELOPED CYLINDER TYPE DIAGRAM
      SIG      = COEF01(3)*SCM**(ALFR-1.)
      IF(GH1H2(3)) GO TO 391
      IF(.NOT.GH(3)) GO TO 700
      SIGM=SIGM+SIG
      IF(ISGCOL.EQ.0.AND.IPRINT) WRITE(ITLIS,83) SIG
 83   FORMAT(10X,39HI DO NOT USE UNCYLINDER REACTION.SIGMA=,E10.4)
      GO TO 700
391   CALL SIGFIH(SIG,I1,I2,3,I4)
C      MULTI POMERON SCATTERING DIAGRAM
 700  CONTINUE
      IF(ECM.GE.ECMB) BPOM=2.840-0.215*LOG(ECM)
      IF(ECM.LT.ECMB) BPOM=2.731-0.4500*LOG(ECM)
       DO 805 I=1,LIMP
      POMGEN(I)=(1.+APOM*I**2)*EXP(-BPOM*I)
805   SIGPS=SIGPS+POMGEN(I)
      PSUM=0.
      DO 810 I=1,LIMP
      PSUM=PSUM+POMGEN(I)/SIGPS
      POMGEN(I)=PSUM
810   CONTINUE
      IF(.NOT.GH1H2(7)) GO TO 900
      SIG=SIGTOT-SIGSUM-SIGM
      IF(SIG.LT.0.)  SIG=0.
      CALL SIGFIH(SIG,I1,I2,7,I4)
900   CONTINUE
      ISGCOL=1
      PJSPNC=PJSC
      PSPINC=PSPC
      PJSPNS=PJSS
      PSPINS=PSPS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GHGH00(GH)
      LOGICAL GH1H2
      LOGICAL GH
      COMMON/H1H2/ GH1H2(11)
      DIMENSION GH(11)
      GH(1)=.FALSE.
      GH(3)=.FALSE.
      GH(5)=.FALSE.
      GH(8)=.FALSE.
      GH(10)=.FALSE.
      GH(11)=.FALSE.
      GH1H2(1)=.FALSE.
      GH1H2(3)=.FALSE.
      GH1H2(5)=.FALSE.
      GH1H2(8)=.FALSE.
      GH1H2(10)=.FALSE.
      GH1H2(11)=.FALSE.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RZERO(A,B,X,R,ETA,MAXFUN,FCN)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
      EXTERNAL FCN
C
C     DATA TETA/1.E-12/
      DATA TETA/1.D-6/
C     TETA,MACHINE DEPENDENT,IS THE COMPUTER PRECISION
      EPSI=ETA
      IF(EPSI.LE.TETA) EPSI=TETA
      FLOW=1.D30
      E=1.
      MC=0
      XA=DMIN1(A,B)
      XB=DMAX1(A,B)
      I=1
c      FA=FCN(XA,I)
      FA=FCN(XA)
      MC=MC+1
      I=2
c      FB=FCN(XB,I)
      FB=FCN(XB)
      IF(FA*FB.GT.0.) GO TO 16
      MC=MC+1
C
    4 X=0.5*(XA+XB)
      R=X-XA
      EE=ABS(X)+E
      IF(R.LE.EE*EPSI) GO TO 18
      F1=FA
      X1=XA
      F2=FB
      X2=XB
    1 CONTINUE
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
c      FX=FCN(X,I)
      FX=FCN(X)
C
      IF(FX*FA.GT.0) GO TO 2
      FB=FX
      XB=X
      GO TO 3
    2 XA=X
      FA=FX
    3 CONTINUE
C
C     PARABOLA ITERATION
C
      F3=FX
      X3=X
      IF(ABS(F1-F2).GE.FLOW*ABS(X1-X2)) GO TO 4
      U1=(F1-F2)/(X1-X2)
      IF(ABS(F2-FX).GE.FLOW*ABS(X2-X)) GO TO 4
      U2=(F2-FX)/(X2-X)
      CA=U1-U2
      CB=(X1+X2)*U2-(X2+X)*U1
      CC=(X1-X)*F1-X1*(CA*X1+CB)
      IF(ABS(CB).GE.FLOW*ABS(CA)) GO TO 8
      U3=0.5*CB/CA
      IF(ABS(CC).GE.FLOW*ABS(CA)) GO TO 4
      U4=U3**2-CC/CA
      IF(U4.LT.0.) GO TO4
      U5=SQRT(U4)
      IF(X.GE.-U3) GO TO 10
      X=-U3-U5
      GO TO9
   10 X=-U3+U5
      GO TO 9
    8 IF(ABS(CC).GE.FLOW*ABS(CB)) GO TO 4
      X=-CC/CB
    9 CONTINUE
      IF(X.LT.XA) GO TO 4
      IF(X.GT.XB) GO TO 4
C
C     TEST FOR OUTPUT
C
      R=ABS(X-X3)
      R1=ABS(X-X2)
      IF(R.GT.R1) R=R1
      EE=ABS(X)+E
      IF(R/EE.GT.EPSI) GO TO 5
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
c      FX=FCN(X,I)
      FX=FCN(X)
      IF(FX.EQ.0.) GO TO 18
      IF(FX*FA.LT.0.) GO TO 7
      XX=X+EPSI*EE
      IF(XX.GE.XB) GO TO 18
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
c      FF=FCN(XX,I)
      FF=FCN(XX)
      FA=FF
      XA=XX
      GO TO 6
    7 XX=X-EPSI*EE
      IF(XX.LE.XA) GO TO 18
      MC=MC+1
c      FF=FCN(XX,I)
      FF=FCN(XX)
      FB=FF
      XB=XX
    6 IF(FX*FF.GT.0.) GO TO 14
   18 CONTINUE
      R=EPSI*EE
      I=3
c      FF=FCN(X,I)
      FF=FCN(X)
      RETURN
   14 F1=F3
      X1=X3
      F2=FX
      X2=X
      X=XX
      FX=FF
      GO TO 3
C
    5 CONTINUE
      F1=F2
      X1=X2
      F2=F3
      X2=X3
      GO TO 1
C
   16 WRITE(16,301)
  301 FORMAT(5X,'RZERO    FCN(A)  AND FCN(B)  HAVE THE SAME SIGN ')
      R=-2.*(XB-XA)
      X=0.
      RETURN
C
   17 WRITE(16,300) MC
  300 FORMAT(10X,'RZERO  ',I5,' CALLS OF THE FUNCTION'/10X,
     1'CALL LIMIT EXCEEDED'///)
      R=-0.5*ABS(XB-XA)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION ZFRAGS(IFL,IFLN,PT2,ZMIN,ZMAX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  RETURN FRACTION Z FROM E+PZ FOR QUARK OR DIQUARK
C  SIMULATION OF U(Z) DISTRIBUTION FROM A.B.KAIDALOV
C         ITEP-116
C
      PARAMETER(ALFT=0.5,ARHO=0.5,APHI=0.,APSI=-2.)
      PARAMETER(AN=-0.5,ALA=-0.75,ALAC=-1.75)
      PARAMETER(AKSI=-1.0,AUSC=-2.0,AUCC=-2.0)
C
      ID1=IABS(IFL)
      ID2=IABS(IFLN)
      IF(MOD(ID2,100).EQ.0) GO TO 15
      GO TO(1,2,3,4),ID2
C  UU-TRAJECTORY
1     ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.0-ZFRAGS)**(ALFT-ARHO)
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 1
C DD-TRAJECTORY
2     ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-ARHO)
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 2
C SS-TRAJECTORY
3     ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-APHI)
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 3
C CC-TRAJECTORY
4     ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-APSI)
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 4
C
15    CALL FLAVOR(ID2,IFL2,IFL3,IFL1,ISPIN,INDEX)
      IF(.NOT.(IFL2.EQ.1.AND.IFL3.EQ.1)) GO TO 16
      IF(.NOT.(IFL2.EQ.1.AND.IFL3.EQ.2)) GO TO 17
      IF(.NOT.(IFL2.EQ.1.AND.IFL3.EQ.3)) GO TO 18
      IF(.NOT.(IFL2.EQ.1.AND.IFL3.EQ.4)) GO TO 19
      IF(.NOT.(IFL2.EQ.2.AND.IFL3.EQ.2)) GO TO 20
      IF(.NOT.(IFL2.EQ.2.AND.IFL3.EQ.3)) GO TO 21
      IF(.NOT.(IFL2.EQ.2.AND.IFL3.EQ.4)) GO TO 22
      IF(.NOT.(IFL2.EQ.3.AND.IFL3.EQ.3)) GO TO 23
      IF(.NOT.(IFL2.EQ.3.AND.IFL3.EQ.4)) GO TO 24
      IF(.NOT.(IFL2.EQ.4.AND.IFL3.EQ.4)) GO TO 25
C UUUU-TRAJECTORY
16    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 16
C UDUD-TRAJECTORY
17    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 17
C USUS-TRAJECTORY
18    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALA-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 18
C UCUC-TRAJECTORY
19    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALAC-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 19
C DDDD-TRAJECTORY
20    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AN-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 16
C DSDS-TRAJECTORY
21    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALA-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 21
C DCDC-TRAJECTORY
22    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*ALAC-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 22
C SSSS-TRAJECTORY
23    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AKSI-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 23
C SCSC-TRAJECTORY
24    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AUSC-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 24
C CCCC-BARYON
25    ZFRAGS=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAGS)**(ALFT-(2.*AUCC-ARHO))
      IF(RNDM(-1.).LE.YF) RETURN
      GO TO 25
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION ZFRAG0(IFL,IFLN,MESON,PT2,ZMIN,ZMAX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  RETURN FRACTION Z FROM E+PZ FOR QUARK OR DIQUARK
C  SIMULATION OF U(Z) DISTRIBUTION FROM A.B.KAIDALOV
C         ITEP-116
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      LOGICAL MESON
      DATA BSLOPE/1.5/
C
      ID1=IABS(IFL)
      ID2=IABS(IFLN)
      IF(MOD(ID1,100).EQ.0) GO TO 400
      IF(MESON) GO TO 601
      IF(MOD(ID2,100).EQ.0) GO TO 600
C     IFL QUARK FRAGMENTS INTO MESON:
      IF(ID1.EQ.3) GO TO 200
      IF(ID2.EQ.3) GO TO 100
      IF(ID1.GE.4.OR.ID2.GE.4) GO TO 700
C     NONSTRANGE QUARK AND NONSTRANGE ANTIQUARK:
C       U(Z)=(1.-Z)**(1.6*PT2-0.5)
      ALFT=1.6*PT2
      IF(ALFT.GE.2.5) ALFT=2.5
50    ZFRAG0=1.-(RNDM(-1.)*(SQRT(1.-ZMAX)-SQRT(1.-ZMIN))+
     *SQRT(1.-ZMIN))**2
         YF=(1.-ZFRAG0)**(ALFT-0.5)
         YP=1./SQRT(1.-ZFRAG0)
       IF(YP*RNDM(-1.).LE.YF) RETURN
         GO TO 50
C
C     NONSTRANGE QUARK AND STRANGE ANTIQUARK:
C       U(Z)=(1.-Z)**ALFT
100   ALFT=1.6*PT2 - 0.45
      IF(ALFT.GE.2.5) ALFT=2.5
      ZFRAG0=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAG0)**ALFT
       IF(RNDM(-1.).LE.YF) RETURN
        GO TO 100
 200  IF(ID2.EQ.3) GO TO 300
C     STRANGE QUARK AND NONSTRANGE ANTIQUARK:
C       U(Z)=(1.-Z)**(1.6*PT2-0.5)
      ALFT=1.6*PT2
      IF(ALFT.GE.2.5) ALFT=2.5
51    ZFRAG0=1.-(RNDM(-1.)*(SQRT(1.-ZMAX)-SQRT(1.-ZMIN))+
     *SQRT(1.-ZMIN))**2
         YF=(1.-ZFRAG0)**(ALFT-0.5)
         YP=1./SQRT(1.-ZFRAG0)
       IF(YP*RNDM(-1.).LE.YF) RETURN
         GO TO 51
C
C        STRANGE QUARK AND STRANGE ANTIQUARK:
C       U(Z)=(1.-Z)**ALFT
300   ALFT=1.6*PT2 - 0.45
      IF(ALFT.GE.2.5) ALFT=2.5
      ZFRAG0=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=(1.-ZFRAG0)**ALFT
       IF(RNDM(-1.).LE.YF) RETURN
        GO TO 300
C     IFL DIQUARK FRAGMENTS INTO BARION%
 400  IF(ID2.EQ.3) GO TO 510
C       DIQUARK AND NONSTRANGE QUARK:
C        U(Z)=Z**1.5
610   ZFRAG0=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=ZFRAG0**1.5
      IF(RNDM(-1.).LE.YF) RETURN
        GO TO 610
C       DIQUARK AND    STRANGE QUARK:
C        U(Z)=Z**2*SQRT(1.-Z)
510   ZFRAG0=ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=ZFRAG0**2*SQRT(1.-ZFRAG0)
      IF(RNDM(-1.).LE.YF) RETURN
        GO TO 510
C     IFL   QUARK FRAGMENTS INTO BARION (IN CASE OF DIQUARK SPLITTING)
C        U(Z)=(1.-Z)**0
 601  ZFRAG0= ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
        RETURN
C     IFL   QUARK FRAGMENTS INTO BARION!
C        U(Z)=(1.-Z)**2
 600  ZFRAG0= ZMIN+RNDM(-1.)*(ZMAX-ZMIN)
      YF=3.*(1.-ZFRAG0)**2
      IF(3.*RNDM(-1.).LE.YF) RETURN
        GO TO 600
C      HEAVY QUARK OR ANTIQUARK
 700   GAMMA=BSLOPE*PT2
       ZM=GAMMA
       IF(ZM.LE.ZMIN) ZM=ZMIN
       IF(ZM.GE.ZMAX) ZM=ZMAX
        UMAX=EXP(-GAMMA/ZM)/ZM
 710    ZFRAG0=RNDM(-1.)*(ZMAX-ZMIN)+ZMIN
        UF=EXP(-GAMMA/ZFRAG0)/ZFRAG0
        IF(RNDM(-1.)*UMAX.GT.UF) GO TO 710
        RETURN
        END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTR(CT,ST,CFI,SFI,PX1,PX2,BACK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  ROTATE OF VECTOR PX1
      DIMENSION ROT(3,3),PX1(3),PX2(3)
      LOGICAL BACK
      ROT(1,1)=CT*CFI
      ROT(1,2)=-SFI
      ROT(1,3)=ST*CFI
      ROT(2,1)=CT*SFI
      ROT(2,2)=CFI
      ROT(2,3)=ST*SFI
      ROT(3,1)=-ST
      ROT(3,2)=0.
      ROT(3,3)=CT
      IF(BACK) GO TO 2
      DO 1 I=1,3
 1    PX2(I)=ROT(I,1)*PX1(1)+ROT(I,2)*PX1(2)+ROT(I,3)*PX1(3)
      RETURN
 2    DO 3 I=1,3
 3    PX2(I)=ROT(1,I)*PX1(1)+ROT(2,I)*PX1(2)+ROT(3,I)*PX1(3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE LORLC(V,PX,E,BACK)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  LORENTZ TRANSFORMATION OF PX MOMENTUM COMPONENTS
      DIMENSION V(3),PX(3)
      LOGICAL BACK
      REAL L
      L=1.
      IF(BACK) L=-1.
      VV=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
      GA=1.D0/SQRT(ABS(1.D0-VV))
      BEP=SPQ(V,PX)
      GABEP=GA*(GA*BEP/(1.+GA)-L*E)
      DO 1 I=1,3
 1    PX(I)=PX(I)+GABEP*V(I)
        RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE FISOB(IK01,P01,AM01,IK02,P02,AM02,SIGEL,
     *IKD,PXD,PYD,PZD,DMAS,ISOB)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C  FORM ISOBAR FROM PION AND NUCLEON
      DIMENSION P01(3),P02(3),P1(3),P2(3)
      ISOB=0
      E1=SQRT( SPQ(P01,P01)+AM01**2)
      E2=SQRT( SPQ(P02,P02)+AM02**2)
      S=AM01**2+AM02**2+2.*E1*E2-2.* SPQ(P01,P02)
      PXC=SQRT(ALAMB(S,AM01**2,AM02**2))/(2.*SQRT(S))
      PT=5.067*PXC
      DM=SQRT(S)
      IK1=IK01
      AM1=AM01
      IK2=IK02
      AM2=AM02
      DO 1 J=1,3
      P1(J)=P01(J)
      P2(J)=P02(J)
 1    CONTINUE
      IF(IBLE(IK2).NE.0) GO TO 10
      IK1=IK02
      AM1=AM02
      IK2=IK01
      AM2=AM01
      DO 2 J=1,3
      P1(J)=P02(J)
      P2(J)=P01(J)
 2    CONTINUE
 10   IQ1=IQLE(IK1)
      IQ2=IQLE(IK2)
      IKS=IQ1+IQ2+2
      GO TO (3,4,5,6),IKS
 3    AK=1.
      IKD=47
      GO TO 7
 4    AK=0.3334
      IF(IQ1.EQ.0) AK=0.6667
      IKD=48
      GO TO 7
 5    AK=0.3334
      IF(IQ1.EQ.0) AK=0.6667
      IKD=46
      GO TO 7
 6    AK=1.
      IKD=45
C  COMPUTE RESONANCE CROSS SECTION
 7    SIGR=SGR(DM,SIGEL,PT)
C     PR=AK*SIGR/SIGEL
      PR=SIGR/SIGEL
      IF(RNDM(-1.).GE.PR) RETURN
C  ISOBAR PARAMETERS
      PXD=P1(1)+P2(1)
      PYD=P1(2)+P2(2)
      PZD=P1(3)+P2(3)
      DMAS=DM
      ISOB=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DOUBLE PRECISION FUNCTION SGR(DM,SIGEL,PT)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C  CALCULATION OF RESONANCE CROSS SECTION
      DM0=1.23000000
      GM=0.12700000
      DMM0=(DM**2-DM0**2)**2
      DMG=(DM0*GM)**2
      PT2=PT**2
      ANORM=SIGEL*PT2
      SGR=ANORM*DMG/(PT2*(DMG+DMM0))
C     SGR=(251.32740000*DMG)/((PT**2)*(DMG+DMM0))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE XQUARK(IAB,XMIN,XMAX,ALFA,BETA,IB1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C   COMPUTE PART OF HADRON MOMENTUM FOR PARTONS
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/COMFLA/MNASEA(12),MNBSEA(12),IFLAS(12),IFLBS(12)
     * ,NUAVAL,
     * NUBVAL,IFLQA1,IFLQB1,IFLQA2,IFLQB2,IFAQQ,IFBQQ
       COMMON/COMVA/ XAVAL1,XAVAL2,XAQQ,XASEA1(12),
     * XASEA2(12),NPOMA
       COMMON/COMVB/ XBVAL1,XBVAL2,XBQQ,XBSEA1(12),
     * XBSEA2(12),NPOMB
      COMMON/COMDIF/ NDIFA,NDIFB
      LOGICAL DIQAN
      COMMON/COMANN/ DIQAN
       IF(IAB.EQ.1) GO TO 10
       NPA=NPOMA
      IF(NPA.EQ.1.AND.NDIFA.EQ.1) GO TO 500
  2    JS=1
       NPAJ=2*NPA-1-JS
       BETA1=BETA
       IF(IABS(IFLQA1).GT.1) BETA1=BETA+1.
C ONLY FOR STRANGE QUARK
       IF(IABS(IFAQQ).EQ.3)  BETA1=0.
       ALFA1=-ALFA
       BETAV=DBLE(NPAJ)*(1.-ALFA)+BETA1
  3    RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 3
       XAVAL1=RNDAV/RNDV
       IF(XAVAL1.LT.XMIN.OR.XAVAL1.GT.XMAX) GO TO 3
       XAQQ=XAVAL1
      XAVAL2=0.
      GO TO 501
500   CALL XSDIS(XAVAL1,XMIN,XMAX)
      CALL XSDIS(XAVAL2,XMIN,XMAX)
      XAQQ=XAVAL1+XAVAL2
      IF(XAQQ.GE.1.0) GO TO 500
501    IF(NPOMA.EQ.1) GO TO 1
       NSA=NPOMA-1
      IF(NSA.GE.13) WRITE(ITLIS,999) NSA
999   FORMAT(/10X,'NUMBER OF COLLISIONS TOO HIGH,NSA=',I5)
      IF(XMIN*DBLE(2*NSA+1).LT.1.0) GO TO 20
      WRITE(ITLIS,995) XMIN,NSA
995   FORMAT(/10X,'...STOP IN XQUARK..XMIN=',E10.4,'NSA=',I3)
      RETURN
20    CONTINUE
       DO 6 JS1=1,NSA
       JS=JS+1
       NPAJ=2*NPA-1-JS
       BETAV=DBLE(NPAJ)*(1.-ALFA)+BETA1
       XB=1.-XAQQ
       IF(XB.LE.XMIN) GO TO 2
  4    RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 4
       XASEA1(JS1)=RNDAV/RNDV*(XB-XMIN)+XMIN
       IF(XASEA1(JS1).GT.XMAX) GO TO 4
       XAQQ=XAQQ+XASEA1(JS1)
       JS=JS+1
       NPAJ=2*NPA-1-JS
       BETAV=DBLE(NPAJ)*(1.-ALFA)+BETA1
       XB=1.-XAQQ
       IF(XB.LE.XMIN) GO TO 2
  5    RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 5
       XASEA2(JS1)=RNDAV/RNDV*(XB-XMIN)+XMIN
       IF(XASEA2(JS1).GT.XMAX) GO TO 5
       XAQQ=XAQQ+XASEA2(JS1)
  6    CONTINUE
  1    CONTINUE
       IF(IB1.EQ.0.AND.NDIFA.EQ.0) GO TO 502
      IF(.NOT.DIQAN) RETURN
      IFL1=IFLQA2
      IFL2=IFAQQ
      CALL X2DIST(X1,X2,IFL1,IFL2)
      XAVAL2=(1.-XAQQ)*X1
      XAQQ  =(1.-XAQQ)*X2
      IF(XAVAL2.LT.XMIN.OR.XAQQ.LT.XMIN) GO TO 2
       RETURN
 502   IF(IABS(IFAQQ).EQ.3) RETURN
       IF(RNDM(-1.).GE.0.5) RETURN
       XSWAP=XAVAL1
       XAVAL1=1.-XAQQ
      XAQQ=1.-XSWAP
       RETURN
  10   CONTINUE
       NPB=NPOMB
      IF(NPB.EQ.1.AND.NDIFB.EQ.1) GO TO 600
  12   JS=1
       NPBJ=2*NPB-1-JS
       BETA1=BETA
       ALFA1=-ALFA
       IF(IABS(IFLQB1).GT.1) BETA1=BETA+1
       BETAV=DBLE(NPBJ)*(1.-ALFA)+BETA1
  13   RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 13
       XBVAL1=RNDAV/RNDV
       IF(XBVAL1.LT.XMIN.OR.XBVAL1.GT.XMAX) GO TO 13
       XBQQ=XBVAL1
      XBVAL2=0.
      GO TO 601
600   CALL XSDIS(XBVAL1,XMIN,XMAX)
      CALL XSDIS(XBVAL2,XMIN,XMAX)
      XBQQ=XBVAL1+XBVAL2
      IF(XBQQ.GE.1.0) GO TO 600
601    IF(NPOMB.EQ.1) GO TO 11
       NSB=NPOMB-1
      IF(NSB.GE.13) WRITE(ITLIS,998) NSB
998   FORMAT(/10X,'NUMBER OF COLLISIONS TOO HIGH,NSB=',I5)
      IF(XMIN*DBLE(2*NSB+1).LT.1.0) GO TO 30
      WRITE(ITLIS,996) XMIN,NSB
996   FORMAT(/10X,'..STOP IN XQUARK...XMIN=',E10.4,'NSB=',I3)
      RETURN
30    CONTINUE
       DO 16 JS2=1,NSB
       JS=JS+1
       NPBJ=2*NPB-1-JS
       BETAV=DBLE(NPBJ)*(1.-ALFA)+BETA1
       XB=1.-XBQQ
       IF(XB.LE.XMIN) GO TO 12
  14   RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 14
       XBSEA1(JS2)=RNDAV/RNDV*(XB-XMIN)+XMIN
       IF(XBSEA1(JS2).GT.XMAX) GO TO 14
       XBQQ=XBQQ+XBSEA1(JS2)
       JS=JS+1
       NPBJ=2*NPB-1-JS
       BETAV=DBLE(NPBJ)*(1.-ALFA)+BETA1
       XB=1.-XBQQ
       IF(XB.LE.XMIN) GO TO 12
  15   RNDAV=RNDM(-1.)**(1./(1.+ALFA1))
       RNDBV=RNDM(-1.)**(1./(1.+BETAV))
       RNDV=RNDAV+RNDBV
       IF(RNDV.GE.1.) GO TO 15
       XBSEA2(JS2)=RNDAV/RNDV*(XB-XMIN)+XMIN
       IF(XBSEA2(JS2).GT.XMAX) GO TO 15
       XBQQ=XBQQ+XBSEA2(JS2)
  16   CONTINUE
  11   CONTINUE
      IF(.NOT.DIQAN) RETURN
      IFL1=IFLQB2
      IFL2=IFBQQ
      CALL X2DIST(X1,X2,IFL1,IFL2)
      XBVAL2=(1.-XBQQ)*X1
      XBQQ  =(1.-XBQQ)*X2
      IF(XBVAL2.LT.XMIN.OR.XBQQ.LT.XMIN) GO TO 10
       RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE PTQUAR(KEY)
      use modifiedDCMParams, only: twpi
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C    COMPUTE PARTON TRANSFERSE MOMENTA
C
C        TRANSFERSE MOMENTUM OF EACH DIQUARK IS EQUAL SUM
C        MOMENTA VALENCE AND SEA QUARKS
       COMMON/COMVA/ XAVAL1,XAVAL2,XAQQ,XASEA1(12),
     * XASEA2(12),NPOMA
       COMMON/COMVB/ XBVAL1,XBVAL2,XBQQ,XBSEA1(12),
     * XBSEA2(12),NPOMB
       COMMON/COMPXA/ PXAV1,PXAV2,PXAQQ,
     *PXAS1(12),PXAS2(12)
       COMMON/COMPYA/ PYAV1,PYAV2,PYAQQ,
     *PYAS1(12),PYAS2(12)
       COMMON/COMPXB/ PXBV1,PXBV2,PXBQQ,
     *PXBS1(12),PXBS2(12)
       COMMON/COMPYB/ PYBV1,PYBV2,PYBQQ,
     *PYBS1(12),PYBS2(12)
      COMMON/PRIMAR/SCM,HALFE,ECM,NJET,IDIN(2),NEVENT,NTRIES
      COMMON/COMIND/ PUD,SIGMA,ALFA,BETA
      COMMON/COMQMA/ AMQUA1,AMQUA2,AMQQA,
     *AMQAS1(12),AMQAS2(12)
      COMMON/COMQMB/ AMQUB1,AMQUB2,AMQQB,
     *AMQBS1(12),AMQBS2(12)
      COMMON/PRIMP0/ P0
      COMMON/COMDIF/ NDIFA,NDIFB
      COMMON/COAUX1/ AMA2,XI2(20),INUM
      COMMON/COAUX2/ X12,X22,X32,PT12,PT22,PT32
      COMMON/COMASS/ AM1,AM2
      LOGICAL DIQAN
      COMMON/COMANN/ DIQAN
      EXTERNAL FMQ2
      EXTERNAL FP02
      DATA AM/0./,BM/15./,MAXFUN/200/,EPSI/0.001/
      DATA PAM2/0./
      DATA SIGMDQ/0.45/
      P02=P0**2
      AMA=AM1
      AMB=AM2
26    IF(KEY.EQ.1) GO TO 20
      AMA2=AMA**2
      XQQ=1.-XAQQ
       NPA=NPOMA
      IF(NPA.EQ.1.AND.NDIFA.EQ.1) GO TO 500
      IF(.NOT.DIQAN) GO TO 27
      XQQ=XAQQ
27     NS=NPA-1
      PTQXS=0.
      PTQYS=0.
      PTQXS0=0.
      PTQYS0=0.
      INUM=0
      IF(NPA.EQ.1) GO TO 5
      DO 4 JS=1,NS
      CALL GAUSPT(PTS1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      CALL GAUSPT(PTS2,SIGMA)
      PHI2=twpi*RNDM(-1.)
       PXAS1(JS)=PTS1*COS(PHI1)
       PYAS1(JS)=PTS1*SIN(PHI1)
       PXAS2(JS)=PTS2*COS(PHI2)
       PYAS2(JS)=PTS2*SIN(PHI2)
      PTQXS=PTQXS+PXAS1(JS)+PXAS2(JS)
      PTQYS=PTQYS+PYAS1(JS)+PYAS2(JS)
      INUM=INUM+2
      XI2(INUM-1)=XASEA1(JS)**2
      XI2(INUM)=XASEA2(JS)**2
4     CONTINUE
      PTQXS0=PTQXS
      PTQYS0=PTQYS
5     CALL GAUSPT(PTV1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      PXAV1=PTV1*COS(PHI1)
      PYAV1=PTV1*SIN(PHI1)
      PXAV2=0.
      PYAV2=0.
      IF(.NOT.DIQAN) GO TO 15
      CALL GAUSPT(PTV2,SIGMDQ)
      PHI2=twpi*RNDM(-1.)
      PXAV2=PTV2*COS(PHI2)
      PYAV2=PTV2*SIN(PHI2)
 15   PTQXS=PTQXS0
      PTQYS=PTQYS0
      PTQXS=PTQXS+PXAV1+PXAV2
      PTQYS=PTQYS+PYAV1+PYAV2
      PXAQQ=-PTQXS
      PYAQQ=-PTQYS
      IF(DIQAN) GO TO 28
      INUM=INUM+2
      XI2(INUM-1)=XAVAL1**2
      XI2(INUM)=XQQ**2
      GO TO 29
 28   INUM=INUM+3
      XI2(INUM-2)=XAVAL1**2
      XI2(INUM-1)=XAVAL2**2
      XI2(INUM  )=XQQ**2
29    FA=FMQ2(AM)
      FB=FMQ2(BM)
      IF(FA*FB.GE.0.) GO TO 26
      CALL RZERO(AM,BM,AMQ2,RES,EPSI,MAXFUN,FMQ2)
      IF(RES.LT.0.) GO TO 26
      IF(NPA.EQ.1) GO TO 3
      DO 2 JS=1,NS
      AMQAS1(JS)=AMQ2-PTS1**2
      AMQAS2(JS)=AMQ2-PTS2**2
2     CONTINUE
3     AMQUA1=AMQ2-PTV1**2
      AMQUA2=0.
      AMQQA=AMQ2-PXAQQ**2-PYAQQ**2
      IF(.NOT.DIQAN) GO TO 1
      AMQUA2=AMQ2-PTV2**2
      GO TO 1
500   CONTINUE
      PBM2=P02+25.
      CALL GAUSPT(PTV1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      PXAV1=PTV1*COS(PHI1)
      PYAV1=PTV1*SIN(PHI1)
      CALL GAUSPT(PTV2,SIGMA)
      PHI2=twpi*RNDM(-1.)
      PXAV2=PTV2*COS(PHI2)
      PYAV2=PTV2*SIN(PHI2)
      PXAQQ=-(PXAV1+PXAV2)
      PYAQQ=-(PYAV1+PYAV2)
      X12=XAVAL1**2
      X22=XAVAL2**2
      X32=XQQ**2
      PT12=PTV1**2
      PT22=PTV2**2
      PT32=PXAQQ**2+PYAQQ**2
      ENER=SQRT(P02+AMA2)
      ESUM=SQRT(AMA2+PT32)+PTV1+PTV2
      IF(ESUM.GE.ENER) GO TO 500
      FA=FP02(PAM2)
      FB=FP02(PBM2)
      IF(FA*FB.GE.0.) GO TO 500
      CALL RZERO(PAM2,PBM2,PNEW2,RES,EPSI,MAXFUN,FP02)
      IF(RES.LT.0.) GO TO 500
      AMQUA1=0.
      AMQUA2=0.
      AMQQA=AMA**2
      PNR=SQRT(PNEW2)/P0
      XAVAL1=XAVAL1*PNR
      XAVAL2=XAVAL2*PNR
      XAQQ=XAVAL1+XAVAL2
1     CONTINUE
      RETURN
  20   CONTINUE
      AMA2=AMB**2
      INUM=0
       XQQ=1.-XBQQ
       NPB=NPOMB
      IF(NPB.EQ.1.AND.NDIFB.EQ.1) GO TO 600
      IF(.NOT.DIQAN) GO TO 37
      XQQ=XBQQ
 37    NS=NPB-1
      PTQXS=0.
      PTQYS=0.
      PTQXS0=0.
      PTQYS0=0.
      IF(NPB.EQ.1) GO TO 7
       DO 6 JS=1,NS
      INUM=INUM+2
      XI2(INUM-1)=XBSEA1(JS)**2
      XI2(INUM)=XBSEA2(JS)**2
      CALL GAUSPT(PTS1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      CALL GAUSPT(PTS2,SIGMA)
      PHI2=twpi*RNDM(-1.)
       PXBS1(JS)=PTS1*COS(PHI1)
       PYBS1(JS)=PTS1*SIN(PHI1)
       PXBS2(JS)=PTS2*COS(PHI2)
       PYBS2(JS)=PTS2*SIN(PHI2)
       PTQXS=PTQXS+PXBS1(JS)+PXBS2(JS)
      PTQYS=PTQYS+PYBS1(JS)+PYBS2(JS)
 6    CONTINUE
      PTQXS0=PTQXS
      PTQYS0=PTQYS
 7    CALL GAUSPT(PTV1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      PXBV1=PTV1*COS(PHI1)
      PYBV1=PTV1*SIN(PHI1)
      PXBV2=0.
      PYBV2=0.
      IF(.NOT.DIQAN) GO TO 35
      CALL GAUSPT(PTV2,SIGMDQ)
      PHI2=twpi*RNDM(-1.)
      PXBV2=PTV2*COS(PHI2)
      PYBV2=PTV2*SIN(PHI2)
 35   PTQXS=PTQXS0
      PTQYS=PTQYS0
      PTQXS=PTQXS+PXBV1+PXBV2
      PTQYS=PTQYS+PYBV1+PYBV2
      PXBQQ=-PTQXS
      PYBQQ=-PTQYS
      IF(DIQAN) GO TO 38
      INUM=INUM+2
      XI2(INUM-1)=XBVAL1**2
      XI2(INUM)=XQQ**2
      GO TO 39
 38   INUM=INUM+3
      XI2(INUM-2)=XBVAL1**2
      XI2(INUM-1)=XBVAL2**2
      XI2(INUM  )=XQQ**2
39    FA=FMQ2(AM)
      FB=FMQ2(BM)
      IF(FA*FB.GE.0.) GO TO 20
      CALL RZERO(AM,BM,AMQ2,RES,EPSI,MAXFUN,FMQ2)
      IF(RES.LT.0.) GO TO 20
      IF(NPB.EQ.1) GO TO 13
      DO 12 JS=1,NS
      AMQBS1(JS)=AMQ2-PTS1**2
      AMQBS2(JS)=AMQ2-PTS2**2
12    CONTINUE
13    AMQUB1=AMQ2-PTV1**2
      AMQUB2=0.
      AMQQB=AMQ2-PXBQQ**2-PYBQQ**2
      IF(.NOT.DIQAN) GO TO 9
      AMQUB2=AMQ2-PTV2**2
      GO TO 9
600   CONTINUE
      CALL GAUSPT(PTV1,SIGMA)
      PHI1=twpi*RNDM(-1.)
      PXBV1=PTV1*COS(PHI1)
      PYBV1=PTV1*SIN(PHI1)
      CALL GAUSPT(PTV2,SIGMA)
      PHI2=twpi*RNDM(-1.)
      PXBV2=PTV2*COS(PHI2)
      PYBV2=PTV2*SIN(PHI2)
      PXBQQ=-(PXBV1+PXBV2)
      PYBQQ=-(PYBV1+PYBV2)
      AMQUB1=0.
      AMQUB2=0.
      AMQQB=AMB**2
      X12=XBVAL1**2
      X22=XBVAL2**2
      X32=XQQ**2
      PT12=PTV1**2
      PT22=PTV2**2
      PT32=PXBQQ**2+PYBQQ**2
      PBM2=P02+25.
      ENER=SQRT(P02+AMA2)
      ESUM=SQRT(AMA2+PT32)+PTV1+PTV2
      IF(ESUM.GE.ENER) GO TO 600
      FA=FP02(PAM2)
      FB=FP02(PBM2)
      IF(FA*FB.GE.0.) GO TO 600
      CALL RZERO(PAM2,PBM2,PNEW2,RES,EPSI,MAXFUN,FP02)
      IF(RES.LT.0.) GO TO 600
      PNR=SQRT(PNEW2)/P0
      XBVAL1=XBVAL1*PNR
      XBVAL2=XBVAL2*PNR
      XBQQ=XBVAL1+XBVAL2
9     CONTINUE
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! PRTEVT removed by CMJ (XCP-3, LANL) on 09/08/2017, it is not called
!    anywhere within LAQGSM (or GSM)
! Purpose: if IPRT is selected by NEVPRT and NJUMP...?
! =====================================================================

C-----------------------------------------------------------------------
      SUBROUTINE PARCRO(ITOT,IK,IKS,IK01,IK02,PX01,PY01,
     *PZ01,AM01,PX02,PY02,PZ02,AM02)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C  CALCULATION OF THE DIAGRAM CROSS SECTION
C  ATTENTION: IK2 - IS BARYON ONLY
C  ITOT=0  COMPUTE ALL POSSIBLE DIAGRAMS
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/P0LAB1/P0,DSM,TKIN
      COMMON/SIGDIA/ CROSD(5),DST
      COMMON/PRINTS/IPRINT
      LOGICAL IPRINT
      DIMENSION IK(5)
      IB1=IB(IK01)
      IB2=IB(IK02)
      IF(IB1.NE.0.OR.IB2.NE.0) GO TO 1
      WRITE(ITLIS, 30) IK01,IK02
 30   FORMAT(10X,'YOU CAN HADRON-NUCLEON CROSS SECTION TREAT'/
     *,' COLLIDING PARTICLES ARE ',I6,' AND ',I6)
      RETURN
 1    IF((IK01.EQ.1120.OR.IK01.EQ.1220) .AND.
     *    IK02.NE.1120.AND.IK02.NE.1220) GO TO 13
      IF(IB1.GT.0.AND.IB2.LE.0.OR.(IB1.LT.0.AND.IB2.EQ.0))GOTO13
C     IF INCIDENT PARTICLE IS MESON OR
C     IF TARGET PARTICLE IS NUCLEON OR
C     IF INCIDENT AND TARGET PARTICLES ARE NOT NUCLEONS
C
      IK1=IK01
      PX1=PX01
      PY1=PY01
      PZ1=PZ01
      AM1=AM01
      IK2=IK02
      PX2=PX02
      PY2=PY02
      PZ2=PZ02
      AM2=AM02
      GO TO 20
C      ELSE  EXCHANGE 1->2  2->1
13    IK1=IK02
      PX1=PX02
      PY1=PY02
      PZ1=PZ02
      AM1=AM02
      IK2=IK01
      PX2=PX01
      PY2=PY01
      PZ2=PZ01
      AM2=AM01
 20   CALL SMARK(IK1,IK2)
      E1=SQRT(AM1**2+PX1**2+PY1**2+PZ1**2)
      E2=SQRT(AM2**2+PX2**2+PY2**2+PZ2**2)
      S=AM1**2+AM2**2+2.*E1*E2-2.*(PX1*PX2+PY1*PY2+PZ1*PZ2)
      DSM=SQRT(S)-AM1-AM2
      TKIN=(S-AM1**2-AM2**2)/(2.*AM2)-AM1
      IF(TKIN.LE.0.02) TKIN=0.02
      P0 = SQRT(TKIN*(TKIN+2*AM1))
C]]]]]]]]
C TO TAKE INTO ACCOUNT NON-MESON AND NON-NUCLEON INELASTIC
C REACTION MORE CAREFULLY I CHANGE MASSES
      IB1=IB(IK1)
      IB2=IB(IK2)
      IS1=IABS(IS(IK1))
      IS2=IABS(IS(IK2))
      AM1N=0.139
      AM2N=0.939
      IF(IS1.NE.0) AM1N=0.497
      IF(IB1.NE.0) AM1N=0.939
      IF(IB1.NE.0.AND.IS1.NE.0) AM1N=1.1156
      DSM=SQRT(ABS(S))-AM1N-AM2N
      TKIN=(S-AM1N**2-AM2N**2)/(2.*AM2N)-AM1N
      IF(TKIN.LE.0.02) TKIN=0.02
      P0 = SQRT(TKIN*(TKIN+2*AM1N))
C]]]]]]]
      DO 99 I = 1,5
99    IK(I) = 0
      IF(ITOT.NE.0) GO TO 40
      IK(1)=1
      IK(2)=2
      IK(3)=3
      IK(4)=4
      IK(5)=5
      CALL CRODIA(P0,ITOT,IK,IK1)
      IF(IPRINT) WRITE(ITLIS,333) ITOT,IK01,IK02,P0,CROSD
 333  FORMAT(4X,' ##### PARCRO: ITOT=',I4,' IK01,IK02==>',2I6,' P0='
     *,E10.4,' CROSD(5)=',5E10.4)
      RETURN
40    GOTO(50,60,70,80,90),IKS
50    IK(1)=1
      GO TO 100
60    IK(2)=2
      GO TO 100
70    IK(3)=3
      GO TO 100
80    IK(4)=4
      GO TO 100
90    IK(5)=5
100   CALL CRODIA(P0,ITOT,IK,IK1)
      RETURN
      END
C-------------------------------------------------------------C
      SUBROUTINE CRODIA(P0,ITOT,IK,IK1)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C     CALCULATION OF THE DIAGRAM CROSS-SECTION
C
      COMMON/KSI/ KSI(4)
      COMMON/COMELX/ SIGEL
      COMMON/COMCRO/ SIGTOT
      COMMON/CSIGA/ SIGAN
      COMMON/DATA10/ PAROM,PARKSI,PARSIG,PARF0
      COMMON/SIGDIA/ CROSD(5),DST
      COMMON/P0LAB1/P00,DSM,TKIN
      DIMENSION IK(5)
      SINEL = SIGTOT-SIGEL
      DO 80 I=1,5
80    CROSD(I)=0.
      DST=0.5
      IF(P0.GE.4.0) DST=0.75
      IF(P0.LT.1.5) DST=0.
      IF(SINEL.LT.0.0001.AND.SIGAN.LT.0.0001) RETURN
      IF(KSI(1).EQ.5) GO TO 50
      IS1=IS(IK1)
C  SEPARATION STRANGE PARTICLES
      IF(IS1.NE.0) GOTO 26
C     SEPARATION OF THE MESON AND NUCLEON REACTION
      IF(KSI(2) .GT. 1) GO TO 19
C     MESON-NUCLEON REACTIONS
      IF(KSI(1)-2)1,7,13
C     POSITIVE MESON+PROTON
 1    IF(IK(1).EQ.0) GO TO 2
      CROSD(1) = CROSS(P0,1,1)
 2    IF(IK(2).EQ.0) GO TO 3
      IF(DSM.LT.0.3) GO TO 3
      CROSD(2) = CROSS(P0,2,1)
 3    IF(IK(3).EQ.0) GO TO 4
      CROSD(3) = CROSS(P0,3,1)
 4    IF(IK(4).EQ.0) GO TO 5
      CROSD(4) = CROSS(P0,4,1)
 5    IF(IK(5).EQ.0) GO TO 6
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 6      GO TO 500
C     NEGATIVE MESON+PROTON
 7    IF(IK(1).EQ.0) GO TO 8
      CROSD(1) = CROSS(P0,1,2)
 8    IF(IK(2).EQ.0) GO TO 9
      IF(DSM.LT.0.3) GO TO 9
      CROSD(2) = CROSS(P0,2,2)
 9    IF(IK(3).EQ.0) GO TO 10
      CROSD(3) = CROSS(P0,3,2)
 10   IF(IK(4).EQ.0) GO TO 11
      CROSD(4) = CROSS(P0,4,2)
 11   IF(IK(5).EQ.0) GO TO 12
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 12       GO TO 500
C     NEUTRAL MESON+PROTON
 13   KN=IABS(IK1)
      COEF = 0.5
      IF(KN.NE.331) GO TO 36
C  NEUTRAL FI-MESON NUCLEON REACTION
      COEF = PARF0/2.
 36   CONTINUE
 14   IF(IK(1).EQ.0) GO TO 15
      CROSD(1) =(CROSS(P0,1,1)+CROSS(P0,1,2))*COEF
 15   IF(IK(2).EQ.0) GO TO 16
      IF(DSM.LT.0.3) GO TO 16
      CROSD(2) =(CROSS(P0,2,1)+CROSS(P0,2,2))*COEF
 16   IF(IK(3).EQ.0) GO TO 17
      CROSD(3) =(CROSS(P0,3,1)+CROSS(P0,3,2))*COEF
 17   IF(IK(4).EQ.0) GO TO 18
      CROSD(4) =(CROSS(P0,4,1)+CROSS(P0,4,2))*COEF
 18   IF(IK(5).EQ.0) GO TO 181
      IF(DSM.LT.0.3) GO TO 181
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 181       GO TO 500
C     NUCLEON-NUCLEON REACTIONS
C
 19   CONTINUE
      IF(KSI(1)-2) 20,200,200
C  PROTON - PROTON REACTION
 20   IF(IK(1).EQ.0) GO TO 21
      CROSD(1) = CROSS(P0,1,5)
 21   IF(IK(2).EQ.0) GO TO 22
      IF(DSM.LT.0.3) GO TO 22
      CROSD(2) = CROSS(P0,2,5)
 22   IF(IK(3).EQ.0) GO TO 23
      CROSD(3) = CROSS(P0,3,5)
 23   IF(IK(4).EQ.0) GO TO 24
      CROSD(4) = CROSS(P0,4,5)
 24   IF(IK(5).EQ.0) GO TO 25
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 25     GO TO 500
C  NEUTRON - PROTON REACTION
 200  IF(IK(1).EQ.0) GO TO 210
      CROSD(1) = CROSS(P0,1,5)
 210  IF(IK(2).EQ.0) GO TO 220
      IF(DSM.LT.0.3) GO TO 220
      CROSD(2) = CROSS(P0,2,5)
 220  IF(IK(3).EQ.0) GO TO 230
      CROSD(3) = CROSS(P0,3,5)
 230  IF(IK(4).EQ.0) GO TO 240
      CROSD(4) = CROSS(P0,4,5)
 240  IF(IK(5).EQ.0) GO TO 250
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 250    GO TO 500
C  SEPARATION OF MESON AND BARYON REACTION
 26   IF(KSI(2).GT.1)GO TO 40
C  STRANGE MESON NUCLEON REACTION
      KSIGO=KSI(1)
      GOTO (27,32,27,32),KSIGO
C  STRANGE MESON NUCLEON ---> K0 N, K+ P, K+ N
 27   IF(IK(1).EQ.0) GO TO 28
      CROSD(1) = CROSS(P0,1,3)
 28   IF(IK(2).EQ.0) GO TO 29
      IF(DSM.LT.0.3) GO TO 29
      CROSD(2) = CROSS(P0,2,3)
 29   IF(IK(3).EQ.0) GO TO 30
      CROSD(3) = CROSS(P0,3,3)
 30   IF(IK(4).EQ.0) GO TO 31
      CROSD(4) = CROSS(P0,4,3)
 31   IF(IK(5).EQ.0) GO TO 500
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
       GO TO 500
C  STRANGE MESON NUCLEON --->AK0 N, K- P, K- N
 32   IF(IK(1).EQ.0) GO TO 33
      CROSD(1) = CROSS(P0,1,4)
 33   IF(IK(2).EQ.0) GO TO 34
      IF(DSM.LT.0.3) GO TO 34
      CROSD(2) = CROSS(P0,2,4)
 34   IF(IK(3).EQ.0) GO TO 35
      CROSD(3) = CROSS(P0,3,4)
 35   IF(IK(4).EQ.0) GO TO 37
      CROSD(4) = CROSS(P0,4,4)
 37   IF(IK(5).EQ.0) GO TO 500
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
       GO TO 500
C  STRANGE BARYON NUCLEON REACTION
 40   IF(2 - IABS(KSI(4))) 41,42,43
 41   COEF = PARSIG
      GO TO 44
 42   COEF = PARKSI
      GO TO 44
 43   COEF = PAROM
 44   IF(IK(1).EQ.0) GO TO 45
      CROSD(1) = CROSS(P0,1,5)*COEF
 45   IF(IK(2).EQ.0) GO TO 46
      IF(DSM.LT.0.3) GO TO 46
      CROSD(2) = CROSS(P0,2,5)*COEF
 46   IF(IK(3).EQ.0) GO TO 47
      CROSD(3) = CROSS(P0,3,5)*COEF
 47   IF(IK(4).EQ.0) GO TO 48
      CROSD(4) = CROSS(P0,4,5)*COEF
 48   IF(IK(5).EQ.0) GO TO 49
      CROSD(5)=SINEL-CROSD(1)-CROSD(2)-CROSD(3)-CROSD(4)
 49      GO TO 500
C    ANTIBARYON NUCLEON REACTION
 50   IF(IK(1).EQ.0) GO TO 51
      IF(TKIN.LT.0.0001) GO TO 53
      CROSD(1) = CROSS(P0,1,6)
 51   IF(IK(2).EQ.0) GO TO 52
      IF(TKIN.LT.0.0001) GO TO 53
      CROSD(2) = CROSS(P0,2,6)
 52   IF(IK(3).EQ.0) GO TO 53
      IF(TKIN.LT.0.0001) GO TO 53
      CROSD(3)=SINEL-SIGAN-CROSD(1)-CROSD(2)
 53   IF(IK(4).EQ.0) GO TO 54
      CROSD(4) = CROSS(P0,3,6)
 54   IF(IK(5).EQ.0) GO TO 500
      CROSD(5) = CROSS(P0,4,6)
 500  CONTINUE
      IF(ITOT.NE.0) RETURN
      IF(CROSD(5).LT.0) CROSD(5)=0.
      SUM=CROSD(1)+CROSD(2)+CROSD(3)+CROSD(4)+CROSD(5)
      IF(SUM.LE.0.00001)ACOE=0.
      IF(SUM.GT.0.00001) ACOE=SINEL/SUM
      DO 501 I=1,5
 501  CROSD(I)=CROSD(I)*ACOE
      RETURN
      END
C-------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION CROSS(P01,NDIAGR,ITYP)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
c
C
C     CALCULATION OF THE DIAGRAM CROSS SECTION
C         BY MEANS OF INTERPOLATION
C
      COMMON /SIGDGR/ SDGR(50,24)
      COMMON /INTERP/ F(6)
      I = NDIAGR+(ITYP-1)*4
      P0 = P01
C     COMPUTE THREE POINTS FOR THE INTERPOLATION
 3    L=1
 4    IF(P0-0.2*L+0.1)6,5,11
 5    CROSS=SDGR(L,I)
      RETURN
C     (INTERPOLATION IS UNNECCESSARY)
 6    IF(L .GT. 1) GO TO 7
      CROSS=SDGR(L,I)
      RETURN
 7    IF(L .GE. 49) GO TO 8
      L1=L-2
      GO TO 9
 8    L1=47
 9    DO 10 K=1,3
      LL1=L1+K
      F(K)=SDGR(LL1,I)
      F(K+3)= 0.2*(LL1-1)
 10   CONTINUE
                 GO TO 12
 11   L=L+1
      IF(L.EQ.50)GOTO 13
      GOTO 4
 12   CROSS=SINTER(P0)
      RETURN
 13   CROSS=SDGR(50,I)
      RETURN
      END
C-------------------------------------------------------------C
