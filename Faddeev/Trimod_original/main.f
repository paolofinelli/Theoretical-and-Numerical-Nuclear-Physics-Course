      PROGRAM TRIMOD
C======================================================================C
C     Output is printed in file output.dat
C     The output file is automatically generated by the program
C
C     PARAMETERS:
C     NQ = NUMBER OF Q POINTS IN THE INTERVAL[0,QMAX] WITH
C          Q0 BEEING THE MIDDLE
C     NP1 = NUMBER OF P POINTS IN THE INTERVAL [0,PMAX] WITH
C           P0 BEEING THE MIDDLE
C     NP2 = NUMBER OF P POINTS IN THE INTERVAL [PMAX,PCUT] WITH
C           PM BEEING THE MIDDLE
C     NPTOT = TOTAL NUMBER OF P POINTS AS NEEDED FOR THE SOLUTION
C             OF THE LIPPMANN SCHWINGER EQUATION IN TMAT
C     NP  = NUMBER OF P POINTS FOR THE SPLINE INTERPOLATION
C     NX  = NUMBER OF X POINTS IN THE INTERVAL [-1.,1.]
C     M   = NUCLEON MASS IN MEV
C     MF  = NUCLEON MASS IN FM-1
C     E0  = STARTING ENERGY IN MEV
C     DE  = ENERGY STEP IN MEV
C     IPF,IQF  INDICES FOR ARBITRARILY CHOSEN P AND Q POINTS FOR THE
C              VECTOR ITERATION PROCEDURE
C
C
C     PRINCIPAL GLOBAL VARIABLES
C
C     Q   ARRAY OF Q POINTS
C     P   ARRAY OF P POINTS
C     S1,S2   SPLINE ELEMENTS
C     PSI0  STARTING VECTOR AT EACH ENERGY
C     X1,X2  VECTORS,X2 RESULTS FROM APPLICATION OF THE FADDEEV KERNEL
C                  TO X1
C     ROLD,RNEW   RATIOS OF X2 AND X1
C     E    ENERGY VALUES DURING THE ENERGY SEARCH
C     ETA  EIGENVALUE OF THE FADDEEV KERNEL
C     T    THE TWO-BODY OFF-SHELL T-MATRIX
C     RATMAX,RATMIN   THE MAXIMAL AND MINIMUM RATIOS OF X2 AND X1
C                     WITH RESPECT TO P AND Q VARIATIONS
C
C======================================================================C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL M,MF
C
      PARAMETER(NQ=10,Q0=1.,QMAX=4.)
C
      PARAMETER(NP1=16,P0=1.,PMAX=1.5*QMAX+0.3,
     1          NP2=18,PM=10.,PCUT=50.,
     2          NPTOT=1+NP1+NP2,NP=NP1+1)
C
      PARAMETER(NX=10)
C
      PARAMETER(M=938.903,HBARC=197.3289,MF=M/HBARC)
C
      PARAMETER(E0=-7.4,DE=-0.1)
C
      PARAMETER(IPF=8,IQF=5)
C
      common/qpoint/Q(NQ),Q2(NQ),WQ(NQ)
C
      common/ppoint/P(NPTOT),P2(NPTOT),WP(NPTOT)
C
      COMMON/SPLIN/WX(NX),X(NX),S1(NX,NP,NQ,NQ),S2(NX,NP,NQ,NQ)
C
      COMMON/TM/TH(NPTOT,NPTOT)
C
      DIMENSION PSI0(NP,NQ),X1(NP,NQ),X2(NP,NQ),ROLD(NP,NQ),RNEW(NP,NQ)
C
      DIMENSION E(10),ETA(10),EIG(10)
C
      DIMENSION H(NP),T(NPTOT,NPTOT,NQ)
C
C
      open(1,file='output.dat',status='replace')

C
      WRITE(1,600)NP1,NP2,NPTOT
      WRITE(1,610)NQ,QMAX,NX
      WRITE(1,620)P0,PMAX,PM,PCUT
  600 FORMAT(/5X,' NP1 =',I3,5X,' NP2 =',I3,5X,' NPTOT =',I3/)
  610 FORMAT(/5X,' NQ =',I3,5X,' QMAX =',E9.3,5X,' NX =',I3/)
  620 FORMAT(/5X,' P0 =',E9.3,5X,' PMAX =',E9.3,5X,' PM =',
     1E9.3,5X,' PCUT =',E9.3/)
C
C
C     GAUSSIAN QUADRATURE POINTS IN Q AND P ARE GENERATED:
C
      CALL GAUSS(0.d0,Q0,QMAX,NQ,Q,Q2,WQ)
      CALL GAUSS(0.d0,P0,PMAX,NP1,P(2),P2(2),WP(2))
      CALL GAUSS(PMAX,PM,PCUT,NP2,P(NP1+2),P2(NP1+2),WP(NP1+2))
      P(1)=0.
      P2(1)=0.
      WP(1)=0.
C
C
C     THE SPLINE ELEMENTS ARE GENERATED:
C
      CALL SPLINE
C
C
C     THE POTENTIAL IN MOMENTUM SPACE IS GENERATED:
C
      CALL POT
C
C
C     THE STARTING ENERGY AND THE STARTING VECTOR ARE GIVEN:
C
      E(1)=(E0/HBARC)*MF
      DO2 IQ=1,NQ
      DO2 IP=1,NP
    2 PSI0(IP,IQ)=1.
C
      IE=0
C
   17 IE=IE+1
      IF(IE.GT.10) STOP
C
      EMEV=E(IE)*HBARC/MF
      WRITE(1,500)EMEV
  500 FORMAT(/'     EMEV =',E12.5/)
C
C     THE TWO-BODY OFF-SHELL T-MATRIX FOR ALL Q-VALUES
C     IS GENERATED:
C
      DO 25 IQ=1,NQ
      E2=E(IE)-0.75*Q2(IQ)
      CALL TMAT(E2)
      DO 26 IP=1,NPTOT
      DO 26 IPP=1,NPTOT
   26 T(IP,IPP,IQ)=TH(IP,IPP)
   25 CONTINUE
C
C     INITIALISATION FOR THE VECTOR ITERATION:
C
      ITER=0
      DO3 IQ=1,NQ
      DO3 IP=1,NP
      X1(IP,IQ)=PSI0(IP,IQ)
    3 ROLD(IP,IQ)=1.
C
   18 ITER=ITER+1
C
      DO4 IQ=1,NQ
      DO5 IP1=1,NP
      SHH=0.
      DO6 IQP=1,NQ
      SH=0.
      DO 7 IP2=1,NP
      S=0.
      DO 8 IX=1,NX
    8 S=S+WX(IX)*S1(IX,IP1,IQP,IQ)*S2(IX,IP2,IQP,IQ)
    7 SH=SH+S*X1(IP2,IQP)
    6 SHH=SHH+SH*WQ(IQP)
    5 H(IP1)=SHH
      DO 9 IP=1,NP
      SHHH=0.
      DO 10 IP1=1,NP
   10 SHHH=SHHH+T(IP,IP1,IQ)*H(IP1)
    9 X2(IP,IQ)=SHHH/(E(IE)-0.75*Q2(IQ)-P2(IP))
    4 CONTINUE
C
C     THE NEW VECTOR X2 HAS BEEN FOUND AFTER APPLICATION OF THE
C     FADDEEV KERNEL TO THE OLD VECTOR X1
C
      DO 11 IP=1,NP
      DO 11 IQ=1,NQ
   11 RNEW(IP,IQ)=X2(IP,IQ)/X1(IP,IQ)
      RATMAX=0.5
      RATMIN=2.
      DO 12 IP=1,NP
      DO 12 IQ=1,NQ
      IF(RNEW(IP,IQ).GT.RATMAX)RATMAX=RNEW(IP,IQ)
      IF(RNEW(IP,IQ).LT.RATMIN)RATMIN=RNEW(IP,IQ)
   12 CONTINUE
      IF(ITER.GT.1)WRITE(1,100)ITER,RATMIN,RNEW(IPF,IQF),RATMAX
  100 FORMAT(I10,3E15.6)
      DF=(RNEW(IPF,IQF)-ROLD(IPF,IQF))/ROLD(IPF,IQF)
C
C     DF IS TAKEN AS A MEASURE FOR THE CONVERGENCE OF THE
C     VECTOR ITERATION AT EACH ENERGY
C
      IF(ABS(DF).LT.1.E-7) GOTO 13
      DO 14 IP=1,NP
      DO 14 IQ=1,NQ
      X1(IP,IQ)=X2(IP,IQ)
   14 ROLD(IP,IQ)=RNEW(IP,IQ)
      GOTO 18
C
   13 ETA(IE)=RNEW(IPF,IQF)
C
C     ETA IS THE APPROXIMATE EIGENVALUE OF THE FADDEEV KERNEL
C
      WRITE(1,200) IE,E(IE),ETA(IE)
  200 FORMAT(/5X,' IE =',I3,5X,' E(IE) =',E12.5,5X,' ETA(IE) =',E12.5/)
C
C     EIG(IE) VANISHES AT THE BINDING ENERGY
C
      EIG(IE)=ETA(IE)-1.
      WRITE(1,300)EIG(IE)
  300 FORMAT(/10X,'  EIG(IE)  =',E12.5/)
      IF(ABS(EIG(IE)).LT.1.E-6) GOTO 15
C
C     ENERGY SEARCH PART:
C
      DO 16 IP=1,NP
      DO 16 IQ=1,NQ
   16 PSI0(IP,IQ)=X2(IP,IQ)
      IF(IE.EQ.1)THEN
      E(IE+1)=E(IE)+(DE/HBARC)*MF
      GOTO 17
      ELSE
      CALL SEARCH(EIG(IE-1),EIG(IE),E(IE-1),E(IE),E(IE+1))
      GOTO 17
      END IF
C
   15 CONTINUE
      WRITE(1,700) EMEV
  700 FORMAT(//10X,'  3-BODY BINDING ENERGY =',E12.5)
      close(1)
      END