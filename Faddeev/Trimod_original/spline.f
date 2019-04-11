      SUBROUTINE SPLINE
C     ===================
C
C     THIS SUBROUTINE DETERMINES THE SPLINE ELEMENTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NP1=16,NP2=18,NPTOT=1+NP1+NP2,NP=NP1+1,NQ=10,NX=10)
C
      DIMENSION SPL1(NP),SPL2(NP)
C
      common/qpoint/Q(NQ),Q2(NQ),WQ(NQ)
C
      common/ppoint/P(NPTOT),P2(NPTOT),WP(NPTOT)
C
      COMMON/SPLIN/WX(NX),X(NX),S1(NX,NP,NQ,NQ),S2(NX,NP,NQ,NQ)
C
C

      CALL SPREP(P,NP)
c      IF=0
c      CALL D01BCF(0,-1.d0,1.d0,0.d0,0.d0,NX,WX,X,IF)
      call gauleg(-1.0d0,1.0d0,X,WX,NX)

      DO 1 IQ=1,NQ
      DO 1 IQP=1,NQ
      C1=0.25*Q2(IQ)+Q2(IQP)
      C2=Q2(IQ)+0.25*Q2(IQP)
      C3=Q(IQ)*Q(IQP)
      DO 1 IX=1,NX
      C4=C3*X(IX)
      PI1=SQRT(C1+C4)
      PI2=SQRT(C2+C4)
      CALL SELEM(P,NP,PI1,SPL1)
      CALL SELEM(P,NP,PI2,SPL2)
      DO 2 IP=1,NP
      S1(IX,IP,IQP,IQ)=SPL1(IP)
    2 S2(IX,IP,IQP,IQ)=SPL2(IP)
    1 CONTINUE
      END
C======================================================================C
      SUBROUTINE SPREP (X,N)
C     =======================
C     THIS SUBROUTINE PREPARES COEFFICIENTS FOR THE SPLINE ELEMENTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),HI(40),U(40)
C
      COMMON /FAKTOR/ FAK1(40,40),FAK2(40,40),FAK3(40,40)
     1,Q(40,40),C(40,40)
C
      IF(N.GT.40)STOP
      U(1)=0.
      HI(2)=X(2)-X(1)
      DO 5 I=1,N
    5  Q(1,I)=0.
      DO 10 I=2,N-1
       AX=X(I+1)-X(I)
       HI(I+1)=AX
       BX=X(I+1)-X(I-1)
       CX=X(I)-X(I-1)
       AL=AX/BX
       AM=1.-AL
       PI=1./(2.-AM*U(I-1))
       U(I)=AL*PI
       DO 15 J=1,N
   15   Q(I,J)=-PI*AM*Q(I-1,J)
       Q(I,I-1)=Q(I,I-1)+PI/(CX*BX)
       Q(I,I)=Q(I,I)-PI/(CX*AX)
   10  Q(I,I+1)=Q(I,I+1)+PI/(AX*BX)
      DO 20 J=1,N
       C(N,J)=0.
       FAK1(N,J)=0.
       FAK2(N,J)=0.
   20  FAK3(N,J)=0.
      DO 25 I=N-1,1,-1
       H1=1./HI(I+1)
       DO 30 J=1,N
        C(I,J)=Q(I,J)-C(I+1,J)*U(I)
   30   FAK1(I,J)=-HI(I+1)*(2.*C(I,J)+C(I+1,J))
       FAK1(I,I)=FAK1(I,I)-H1
       FAK1(I,I+1)=FAK1(I,I+1)+H1
       DO 25 J=1,N
        FAK2(I,J)=3*C(I,J)
   25   FAK3(I,J)=(C(I+1,J)-C(I,J))*H1
      END
C======================================================================C
      SUBROUTINE SELEM(X,N,XA,SPL)
C     ==============================
C
C     THIS SUBROUTINE DETERMINES THE SPLINE ELEMENTS SPL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),SPL(N)
C
      COMMON /FAKTOR/ FAK1(40,40),FAK2(40,40),FAK3(40,40)
     1,QDUM(40,40),CDUM(40,40)
C
      IF (XA .LT.X(1).OR. XA.GT.X(N))STOP
      I=-1
   10  I=I+1
       IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10
      IF (I .EQ. 0) I=1
      DX=XA-X(I)
      DO 20 J=1,N
   20  SPL(J)=((FAK3(I,J)*DX+FAK2(I,J))*DX+FAK1(I,J))*DX
      SPL(I)=SPL(I)+1.
      END