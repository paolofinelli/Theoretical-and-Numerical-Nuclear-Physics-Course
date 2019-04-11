       SUBROUTINE POT
C     ===================
C
C     THIS SUBROUTINE DETERMINES A MALFLIET-TJION POTENTIAL
C     IN MOMENTUM SPACE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL M,MF
C
      PARAMETER(HBARC=197.3289,M=938.903,MF=M/HBARC)
C
      PARAMETER(NP1=16,NP2=18,NPTOT=1+NP1+NP2)
C
      common/ppoint/P(NPTOT),P2(NPTOT),WP(NPTOT)
C
      COMMON/VPOT/V(NPTOT,NPTOT)
C
      DIMENSION V0(2),A(2),V0M(2)
C
      DATA V0/-570.316,1438.4812/,A/1.55,3.11/
C
      PI=3.1415926
      DO 2 I=1,2
    2 V0M(I)=(V0(I)/HBARC)*MF
C
      DO 4 IP=1,NPTOT
      DO 4 IPP=1,IP
      X1=P(IP)
      X2=P(IPP)
      VV=0.
      DO 1 I=1,2
      IF(X1*X2.EQ.0.)THEN
      QV=2./(X1**2+X2**2+A(I)**2)
      ELSE
      Z=(A(I)**2+X1**2+X2**2)/(2.*X1*X2)
      IF(Z.LT.1.E2)THEN
      QV=(LOG(Z+1.)-LOG(Z-1.))/(2.*X1*X2)
      ELSE
      X=1./Z
      QV=(X+X**3/3.+X**5/5.+X**7/7.)*2./(2.*X1*X2)
      END IF
      END IF
      VV=VV+QV*V0M(I)/PI
    1 CONTINUE
      V(IP,IPP)=VV
    4 V(IPP,IP)=V(IP,IPP)
      RETURN
      END