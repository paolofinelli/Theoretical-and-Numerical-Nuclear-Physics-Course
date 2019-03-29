This file has a fortran program followed by an inpurt data file.
        PROGRAM RPA3
C   Written by G. Bertsch, based on a program by Shalom Shlomo.
C   Originally described in Nuclear Physics A243 (1975), p. 507
C   Parameters:
C     NRMAX  - maximum number of points in spatial grid
C     NHMAX  - maximum number of hole wave functions
C     MATMAX - maximum dimensionality of matrix
C
        Parameter (NRMAX=50,NHMAX=50,MATMAX=40)
C
C   Principal Variables:
C     DEL         - step size in fm
C     NGRID       - number of points in wave function
C     N=NGRID/NX  - number of spatial points in response matrix
C     V,VS,VC     - central,spin-orbit, and Coulomb potential
C     VNN,VNP     - nn and np residual interactions
C     UH          - hole wavefunction
C     LH,JH,NQ,EH - quantum numbers,charge, and energy of hole
C     L,EX        - multipolarity and frequency of response
C     G           - independent particle response matrix
C     B           - rpa polarizability matrix (1+VG)**(-1)
C  
        COMPLEX A(MATMAX,MATMAX),G(MATMAX,MATMAX),
     &  B(MATMAX,MATMAX),UP(NRMAX),WP(NRMAX),UP1(NRMAX),
     &  WP1(NRMAX),DRHOF(MATMAX),DRHORPA(MATMAX),SFREE,SRPA
        DIMENSION VNN(MATMAX),VNP(MATMAX),VEXT(MATMAX)
        COMMON/PARAM/DEL,DEL2,N,NX,N2
        COMMON/WFN/LH(NHMAX),JH(NHMAX),EH(NHMAX),
     &             NQ(NHMAX),UH(NHMAX,NRMAX),NH
        COMMON/POTEN/V(NRMAX),VS(NRMAX),VC(NRMAX),EL,NGRID,H
C  Number of mesh points and approximate spacing for response
c  Initialize clebsch routine
        CALL FACTOR
C  Set up wave functions for ground state
        CALL STATIC(AMASS,Z)
        PRINT*,'AMASS=',AMASS,'Z=',Z
        READ(5,*) N,DEL2
        WRITE(*,'(I5,'' POINTS IN RESPONSE MESH, WITH MESH SIZE='',
     &            F4.2)') N,DEL2
        NX=DEL2/DEL
        DEL2=NX*DEL
C  Residual interaction
        CALL VRESID(VNN,VNP,AMASS)
c  Multipolarity,starting energy, ending energy, and energy step
c  for response function
        READ(5,*) L,EX,EXM,DEX,GAM
        LL=2*L
        PRINT*,' L= ',L,' GAM= ', GAM  
C  Make multipole external field
        CALL EXTF(VEXT,L)
C  Loop on energy
        WRITE(*,'(''    Energy         Free response'',
     &            ''          RPA Response'')')
   11   IF(EX.GT.EXM) GO TO 40
C  Construct independent particle response
        N2=2*N
        DO 13 M1=1,N2
          DO 12 M2=1,N2
12        G(M1,M2)=0.0
13      CONTINUE
C       LOOP OVER HOLES
        DO 15 I2=1,NH
C  Loop over particle angular momentum (twice j)
          JMIN=IABS(LL-JH(I2))
          JMAX=LL+JH(I2)
C  Angular momentum triangle condition
          DO 14 JP=JMIN,JMAX,2
C  Determine orbital angular momentum from parity 
            IF(MOD((JP+1)/2+LH(I2)+L,2).NE.0) THEN
              LP=(JP-1)/2 
            ELSE 
              LP=(JP+1)/2
            END IF
            CALL CLEBSCH(JH(I2),JP,LL,1,-1,0,TEMP)
C  Angular momentum factor
            TEMP=(JP+1)*(JH(I2)+1)*TEMP**2/4./3.1416/FLOAT(LL+1) 
C  Particle Green's function
            CALL GREEN(UP,WP,EX+EH(I2),GAM,NQ(I2),LP,JP)
            CALL GREEN(UP1,WP1,-EX+EH(I2),-GAM,NQ(I2),LP,JP)
C  Independent particle response
            DO 17 M1=1,N
              MA=M1*NX
              DO 16 M2 = M1,N
                MB=M2*NX
C  Protons and neutrons are in separate blocks
                IF (NQ(I2).EQ.1) THEN  
                  MX=0
                ELSE 
                  MX=N
                END IF
   16         G(M1+MX,M2+MX)=G(M1+MX,M2+MX)+TEMP*UH(I2,MA)
     &            *UH(I2,MB)*(UP(MA)*WP(MB)+UP1(MA)*WP1(MB))
   17       CONTINUE
   14     CONTINUE
   15   CONTINUE 
C  Multiply G by the residual interaction V
        DO 22 M1=1,N2
          DO 21 M2=1,N2
            G(M2,M1)=G(M1,M2)
            IF(M1.LE.N) THEN
               IF(M2.LE.N) THEN 
                  A(M1,M2)=G(M1,M2)*VNN(M2)
               ELSE 
                  A(M1,M2)=G(M1,M2-N)*VNP(M2-N)
               END IF
            ELSE 
              IF(M2.LE.N) THEN 
                A(M1,M2)=G(M1,M2+N)*VNP(M2)
              ELSE 
                A(M1,M2)=G(M1,M2)*VNN(M2-N)
              END IF
            END IF
            IF(M1.EQ.M2) A(M1,M2)=1.0+A(M1,M2)
   21     CONTINUE
   22   CONTINUE
C  Invert   1+VG
        CALL MATR(N2,A,B)
C  Evaluate response to field VEXT
        SFREE=0.0
        DO 35 M1=1,N2
          DRHOF(M1)=0.0
          DO 34 M2=1,N2
   34     DRHOF(M1)=DRHOF(M1)+G(M1,M2)*VEXT(M2)*DEL2
C  Free response to field VEXT
   35   SFREE=SFREE+VEXT(M1)*DRHOF(M1)*DEL2
        SRPA=0.0
        DO 36 M1=1,N2
          DRHORPA(M1)=0.0
          DO 37 M2 = 1,N2
   37     DRHORPA(M1)=DRHORPA(M1)+B(M1,M2)*DRHOF(M2)
C  RPA response to field VEXT
   36   SRPA=SRPA+VEXT(M1)*DRHORPA(M1)*DEL2
        WRITE(*,'(F11.5,4E12.4)') EX,SFREE,SRPA
        SUMF=SUMF+AIMAG(SFREE)*EX*DEX/3.1416
        SUMRPA=SUMRPA+AIMAG(SRPA)*EX*DEX/3.1416
        EX=EX+DEX
        GO TO 11
   40   CALL SUMRULE(VEXT,VALUE,L)
        PRINT*,' TOTAL STRENGTH IN FREE RESPONSE, RPA RESPONSE, 
     &  AND SUM RULE'
        WRITE(*,'(21X,E12.4,E14.4,E15.4)') SUMF,SUMRPA,VALUE
        END
        SUBROUTINE STATIC(A,Z)
C   A static potential and occupied single particle levels
C   are created here, using a Woods-Saxon potential.  The user
C   should replace this routine with corresponding Hartree-Fock
C   quantities.  Input: DEL, NH, A, Z, and quantum numbers of
C   occupied orbits.  Output:  wavefunctions UH and energies E.
        PARAMETER (NRMAX=50,NHMAX=50)
        DIMENSION U(NRMAX)
        COMMON/DENSITY/RHOP(NRMAX),RHON(NRMAX)
        COMMON/PARAM/DEL,DEL2,N,NX,N2
        COMMON/WFN/LH(NHMAX),JH(NHMAX),EH(NHMAX),
     &             NQ(NHMAX),UH(NHMAX,NRMAX),NH
        COMMON/POTEN/V(NRMAX),VS(NRMAX),VC(NRMAX),EL,NGRID,H
        DATA VR,VT,VSO,A0,R0/-53.0,20.0,-15.5,0.65,1.25/
C  Usual Woods-Saxon parameters
        READ(5,*) DEL,NGRID
        PRINT*,'DEL=',DEL,'NGRID=',NGRID
        READ(5,*) A,Z
        PRINT*,'A=',A,'Z=',Z
C  Mass and charge of nucleus
        RR=(A-1)**(0.3333)*R0
        H=20.75
        DD=DEL**2/H
        IF(Z.GT.0) Z=Z-1
C  Set up the Woods-Saxon Potential with central,spin-orbit, and
C  Coulomb
        DO 6 I = 1, NGRID
          R=DEL*I
          TEMP=(R-RR)/A0
          IF(TEMP.GT.30) TEMP=30
          EX=EXP(TEMP)
          V(I)=(VR+VT*(A-2*Z-1)/A)/(1+EX)
          VS(I)=VSO*EX/((1.+EX)**2*R*A0)
          VC(I)=1.44*Z/R
          IF(R.LT.RR) VC(I)=1.44*Z/RR*(1.5-0.5*(R/RR)**2)
          VC(I)=VC(I)-2*VT*(A-2*Z-1)/A/(1+EX)
   6    CONTINUE
C  Determine occupied wave functions
        A=0.
        Z=0.
        NH=1
        DO 3 I=1,NGRID
          RHOP(I)=0.0
   3    RHON(I)=0.0
   1    READ(5,*) LH(NH),JH(NH),NQ(NH),NODE
        IF (LH(NH).LT.0) THEN
          NH=NH-1
          RETURN
        END IF
        ALL = LH(NH)*(LH(NH)+1)*H
        IF (2*LH(NH).LT.JH(NH)) THEN 
          FL=LH(NH)
        ELSE 
          FL=-(LH(NH)+1)
        END IF     
C  Energy of state is found by bisection
        EMIN=-50.0
        EMAX=0.0
        DO 7 K = 1,25
C  Integrate schroedinger equation
          ETRIAL=(EMIN+EMAX)/2.0
          U(1)=(0.1)**(lh(nh)+1)
          U(2)=(0.2)**(lh(nh)+1)
          ND=0
          DO 8 I =2,NGRID-1
            R=DEL*I
            VV=V(I)+NQ(NH)*VC(I)+FL*VS(I)+ALL/R**2
            SX=DD*(VV-ETRIAL)
            U(I+1)=(2+SX)*U(I)-U(I-1)  
            IF(U(I+1)*U(I).LT.0)ND=ND+1
   8      CONTINUE
          IF(ND.GT.NODE) THEN
            EMAX=ETRIAL
          ELSE
            EMIN=ETRIAL
          END IF
   7    CONTINUE
        SUM=0.0
        DO 9 I=1,NGRID
          SUM=SUM+U(I)**2*DEL
   9    CONTINUE 
        SUM=SQRT(SUM)
        DO 10 I=1,NGRID
          R=DEL*I
          UH(NH,I)=U(I)/SUM
          IF(NQ(NH).EQ.1) THEN
           RHOP(I)=RHOP(I)+(JH(NH)+1)*(UH(NH,I)/R)**2/4./3.1416
          ELSE 
           RHON(I)=RHON(I)+(JH(NH)+1)*(UH(NH,I)/R)**2/4./3.1416
          END IF
  10    CONTINUE
        EH(NH)=ETRIAL
        A=A+JH(NH)+1
        IF(NQ(NH).EQ.1)Z=Z+JH(NH)+1
        WRITE(*,'('' L,2J,Q,NODE,E='',4I3,F10.5)')
     &               LH(NH),JH(NH),NQ(NH),NODE,ETRIAL
        NH=NH+1
        GOTO 1            
        END
C
        SUBROUTINE VRESID(VNN,VNP,A)
        PARAMETER (NRMAX=50)
        DIMENSION VNN(1),VNP(1)
        COMMON/PARAM/DEL,DEL2,N,NX,N2
        DATA RHO,R0,A0/0.16,1.20,0.50/
        READ(5,*) T0,T3,X,VSCAL
        WRITE(*,'('' T0='',F10.2,''   T3='',F8.0,''   X='',F4.2,  
     &            ''   VSCAL='',F4.2)') T0,T3,X,VSCAL
        DO 1 I=1,N
          R=I*DEL2
          TEMP=0.5*T3*RHO/(1.0+EXP((R-R0*A**0.33333)/A0))
C  PN interaction
          VNP(I)=DEL2*VSCAL*(T0*(1+X/2)+TEMP)/(I*DEL2)**2
C  PP interaction
          VNN(I)=DEL2*VSCAL*(T0*(1-X)+TEMP)/(I*DEL2)**2/2.
1       CONTINUE
        RETURN
        END
C
        SUBROUTINE EXTF(F,L)
C  pure multipole field
        DIMENSION F(1)
        COMMON/PARAM/DEL,DEL2,N,NX,N2
        READ(5,*) I
        PRINT*,' type of field: ',I
        IF(L.EQ.0) THEN
          K=2
        ELSE
          K=L
        END IF
        DO 1 M1=1,N
          R=DEL2*M1
          F(M1)=R**K
          IF(I.EQ.0) THEN
            F(M1+N)=F(M1)
          ELSE IF (I.EQ.1) THEN
            F(M1+N)=-F(M1)
          ELSE
            F(M1+N)=0.0
          ENDIF
    1   CONTINUE
        END
      SUBROUTINE FACTOR
      COMMON/BL1/FACLOG(130)
      FACLOG(1)=0
      FACLOG(2)=0
      FN=1
      DO 95 M=3,130
        FN=FN+1
  95  FACLOG(M)=FACLOG(M-1)+ALOG(FN)
      RETURN
      END
      SUBROUTINE CLEBSCH(J1,J2,J3,M1,M2,M3,ANS)
C  calculates Clebsch-Gordon coefficients
C  Definition found in: "Elementary Theory of Angular Momentum", M.E. Ros
C  Equation 3.19 (John Wiley)
      COMMON/BL1/FACLOG(130)
      IF(M1+M2.NE.M3)GO TO 130
      IF(J3.LE.J1+J2.AND.J3.GE.ABS (J1-J2))GO TO 110
      GO TO 130
  110 IA=(J1+J2-J3+2)/2.0
      IB=(J3+J1-J2+2)/2.0
      IC=(J3+J2-J1+2)/2.0
      ID=(J1+J2+J3+4)/2.0
      IE=(J1+M1+2)/2.0
      IF=(J1-M1+2)/2.0
      IG=(J2+M2+2)/2.0
      IH=(J2-M2+2)/2.0
      II=(J3+M3+2)/2.0
      IJ=(J3-M3+2)/2.0
      FIRST=0.5*(ALOG(J3+1.0)+FACLOG(IA)+FACLOG(IB)+FACLOG(IC)
     &   +FACLOG(IE)+FACLOG(IF)+FACLOG(IG)+FACLOG(IH)+FACLOG(II)
     &   +FACLOG(IJ)-FACLOG(ID))
      PART1=EXP (FIRST)
      NUMIN1= ABS(AMIN1(((J3-J2+M1)/2.0),((J3-J1-M2)/2.0),0.0))+1.0
      NUMAX1=MIN1(((J1+J2-J3)/2.0),((J1-M1)/2.0),((J2+M2)/2.0))+1
      KB=(J3-J2+M1)/2.0+1.0
      KC=(J3-J1-M2)/2.0+1.0
      FF=0
      DO 120  NUM1=NUMIN1,NUMAX1
        NU=NUM1-1
        SECOND=-(FACLOG(NUM1)+FACLOG(IA-NU)+FACLOG(IF-NU)
     &          +FACLOG(IG-NU)+FACLOG(KB+NU)+FACLOG(KC+NU))
        CON=((-1)**NU)*EXP (SECOND)
        FF=FF+CON
  120 CONTINUE
      ANS=PART1*FF
      GO TO 140
  130 ANS=0
  140 CONTINUE
      RETURN
      END
      SUBROUTINE GREEN(UP,WP,E,GAM,NT,L,J)
      PARAMETER (NRMAX=50)
      COMPLEX UP(1),WP(1),ONE,EI,WW,ECOMP,S(NRMAX)
      COMMON/PARAM/DEL,DEL2,N,NX,N2
      COMMON/POTEN/V(NRMAX),VS(NRMAX),VC(NRMAX),EL,NGRID,H
      DD=DEL**2/H
      ONE=(1.,0.)
      EI=(0.,1.)
      ECOMP=E*ONE+GAM*EI/2.
      ALL=L*(L+1)*H
      FL=-(L+1)
      IF(2*L.LT.J) FL=L
      NMAX=NGRID
      NMAX1=NMAX+1
      NL=NMAX-1
C  Initialize regular solution
      UP(1)=(0.1)**(L+1)
      UP(2)=(0.2)**(L+1)
      DO 10 I=2,NMAX1
        K=I-1
        R=DEL*(I-1)
        S(I-1)=DD*(V(K)+NT*VC(K)+FL*VS(K)+ALL/R**2-ECOMP)
        IF(I.LT.3.OR.I.GT.NMAX) GO TO 10
        UP(I)=(2.+S(I-1))*UP(I-1)-UP(I-2)
   10 CONTINUE
      Z=-S(NL)
      IF(Z.GT.0.0) THEN
C  Initialize irregular solution to outgoing wave
        PK=SQRT(Z)
        WP(NL)=(1.0,0.0)
        WP(NMAX)=ONE*(1.-PK**2./2.)+EI*(PK)
      ELSE
        WP(NL)=(0.0,0.0)
        WP(NMAX)=0.001
      ENDIF
      DO 30 K=3,NMAX
        I=NMAX1-K
        WP(I)=(2.+S(I+1) )*WP(I+1)-WP(I+2)
   30 CONTINUE
C  Wronskian
      NQ=NMAX/2
      WW=-(UP(NQ)*WP(NQ+1)-UP(NQ+1)*WP(NQ))/DEL
      DO 50 I=1,NMAX
        WP(I)=WP(I)/H/WW
   50 CONTINUE
      RETURN
      END
        SUBROUTINE SUMRULE(F,VALUE,J)
C  calculates energy-weighted sum rule
        PARAMETER (NRMAX=50,MATMAX=40)
        DIMENSION F(1),FP(MATMAX)
        COMMON/DENSITY/RHOP(NRMAX),RHON(NRMAX)
        COMMON/PARAM/DEL,DEL2,N,NX,N2
        DO  3 I=1,N2
          IF(I.EQ.1.OR.I.EQ.N+1) THEN 
            FL=0.0
          ELSE 
            FL=F(I-1)
          ENDIF
          IF(I.EQ.N.OR.I.EQ.2*N) THEN 
            FG=(2*F(I)-F(I-1))
          ELSE 
            FG=F(I+1)
          ENDIF
          FP(I)=(FG-FL)/2.0/DEL2
    3    CONTINUE
        DO 5 I=1,N
          R=I*DEL2
          SS=SS+((FP(I)**2+(F(I)/R)**2*J*(J+1))*RHOP(I*NX)+
     &    (FP(I+N)**2+(F(I+N)/R)**2*J*(J+1))*RHON(I*NX))*
     &    R**2*DEL2*4.*3.1416
    5   CONTINUE
        VALUE=20.75/4./3.1416*SS
        RETURN
        END
      SUBROUTINE MATR(NMAX,C,D)
      PARAMETER (MATMAX=40)
      IMPLICIT COMPLEX (A-H,O-Z)
      DIMENSION U(MATMAX),V(MATMAX),C(MATMAX,MATMAX),D(MATMAX,MATMAX)
      DET=1.0
      DO 1 N=1,NMAX
        DO 1 M=1,NMAX
          D(N,M)=0.0
          IF(N.EQ.M) D(N,M)=1.0
    1 CONTINUE
      DO 10 N=1,NMAX
        T=C(N,N)
        IF(CABS(T).LE.1.E-10) GO TO 3
        GO TO 6
    3   J=N
    4   IF(J.GT.NMAX) GO TO 11
        J=J+1
        T=C(N,J)
        IF(CABS(T).LE.1.E-10) GO TO 4
        DO 5 K=1,NMAX
          U(K)=C(N,K)
          V(K)=D(N,K)
          C(N,K)=C(J,K)
          D(N,K)=D(J,K)
          C(J,K)=U(K)
    5     D(J,K)=V(K)
    6   DO 8 K=1,NMAX
          IF(K.EQ.N) GO TO 8
          A=C(K,N)/C(N,N)
          DO 7 L=1,NMAX
            C(K,L)=C(K,L)-A*C(N,L)
            D(K,L)=D(K,L)-A*D(N,L)
    7       CONTINUE
    8     CONTINUE
        B=C(N,N)
        DET=DET*B
        DO 10 M=1,NMAX
          C(N,M)=C(N,M)/B
          D(N,M)=D(N,M)/B
   10 CONTINUE
      RETURN
   11 WRITE (6,'('' Matrix not invertable '')')
      RETURN
      END
 0.25 50
16 8
0 1 0 0
0 1 1 0
1 3 0 0 
1 3 1 0
1 1 0 0
1 1 1 0
-1 0 0 0
10 0.5
-1100.  15000. 0.5  0.93
1 0.0   40.5     1.0  0.0 
-1

