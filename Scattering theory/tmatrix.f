C     *******************************************************
C         Example program used to evaluate the 
C         T-matrix following Kowalski's method (eqs V88 & V89
C         in Brown and Jackson)  
C         for positive energies only
C         The program is set up for S-waves only
C         Coded by : Morten Hjorth-Jensen
C         Language : FORTRAN 90
C     *******************************************************



C               ******************************
C                 Def of global variables
C               ******************************


      MODULE constants
         DOUBLE PRECISION , PUBLIC :: p_mass, hbarc
         PARAMETER (p_mass =938.926D0, hbarc = 197.327D0)
      END MODULE constants

      MODULE mesh_variables           
         INTEGER, PUBLIC :: n_rel
         PARAMETER(n_rel=48)
         DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: ra(:), wra(:)
      END MODULE  mesh_variables


C               ******************************
C                   Begin of main program
C               ******************************


      PROGRAM t_matrix
      USE mesh_variables
      IMPLICIT NONE
      INTEGER istat

      ALLOCATE( ra (n_rel), wra (n_rel),  STAT=istat )
      CALL rel_mesh               ! rel mesh & weights
      CALL t_channel              ! calculate the T-matrix
      DEALLOCATE( ra,wra, STAT=istat )

      END PROGRAM t_matrix

C     *********************************************************
C                    obtain the t-mtx
C                    vkk is the box potential
C                    f_mtx is equation V88 og Brown & Jackson
C     *********************************************************

      SUBROUTINE t_channel
      USE mesh_variables
      IMPLICIT NONE
      INTEGER istat, i,j
      DIMENSION vkk(:,:),f_mtx(:),t_mtx(:)
      DOUBLE PRECISION, ALLOCATABLE :: vkk,t_mtx,f_mtx
      DOUBLE PRECISION t_shell

      ALLOCATE(vkk (n_rel,n_rel), STAT=istat)
      CALL v_pot_yukawa(vkk) ! set up the box potential in routine vpot
      ALLOCATE(t_mtx (n_rel), STAT=istat) ! allocate space in heap for T
      ALLOCATE(f_mtx (n_rel), STAT=istat) ! allocate space for f
      DO i=1,n_rel               ! loop over positive energies e=k^2
         CALL f_mtx_eq(f_mtx,vkk,i)  ! solve eq. V88
         CALL principal_value(vkk,f_mtx,i,t_shell) ! solve Eq. V89 
         DO j=1,n_rel    ! the t-matrix
            t_mtx(j)=f_mtx(j)*t_shell
            IF(j == i) WRITE(6,*) ra(i) ,t_mtx(i)
c     &                            DATAN(-ra(i)*t_mtx(i))
         ENDDO
      ENDDO
      DEALLOCATE(vkk ,  STAT=istat)
      DEALLOCATE(t_mtx, f_mtx,  STAT=istat)
 1000 FORMAT( I3, 2F12.6) 

      END SUBROUTINE t_channel

C     ***********************************************************
C          The analytical expression for the box potential
C          of exercise 1 and 12
C          vkk is in units of fm^-2 (14 MeV/41.47Mevfm^2, where 
C          41.47= \hbarc^2/mass_nucleon),  
C          ra are mesh points in rel coordinates, units of fm^-1
C     ***********************************************************

      SUBROUTINE v_pot_box(vkk)
      USE mesh_variables
      USE constants
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION  vkk, r_0, v_0, a, b, fac
      PARAMETER(r_0=2.7d0,v_0=0.33759d0)  !r_0 in fm, v_0 in fm^-2 
      DIMENSION vkk(n_rel,n_rel)

      DO i=1,n_rel     ! set up the free potential
         DO j=1,i-1 
            a=ra(i)+ra(j)
            b=ra(i)-ra(j)
            fac=v_0/(2.d0*ra(i)*ra(j))
            vkk(j,i)=fac*(DSIN(a*r_0)/a-DSIN(b*r_0)/b)
            vkk(i,j)=vkk(j,i)
         ENDDO
         fac=v_0/(2.d0*(ra(i)**2))
         vkk(i,i)=fac*(DSIN(2.d0*ra(i)*r_0)/(2.d0*ra(i))-r_0)
      ENDDO

      END  SUBROUTINE v_pot_box

C     ***********************************************************
C          The analytical expression for  a Yukawa potential
C          in the l=0 channel
C          vkk is in units of fm^-2,  
C          ra are mesh points in rel coordinates, units of fm^-1
C          The parameters here are those of the Reid-Soft core
C          potential, see Brown and Jackson eq. A(4)
C     ***********************************************************

      SUBROUTINE v_pot_yukawa(vkk)
      USE mesh_variables
      USE constants
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION  vkk, mu1, mu2, mu3, v_1, v_2, v_3, a, b, fac
      PARAMETER(mu1=0.49d0,v_1=-0.252d0) 
      PARAMETER(mu2=7.84d0,v_2=-39.802d0) 
      PARAMETER(mu3=24.01d0,v_3=156.359d0) 
      DIMENSION vkk(n_rel,n_rel)

      DO i=1,n_rel     ! set up the free potential
         DO j=1,i 
            a=(ra(j)+ra(i))**2
            b=(ra(j)-ra(i))**2
            fac=1./(4.d0*ra(i)*ra(j))
            vkk(j,i)=v_1*fac*DLOG((a+mu1)/(b+mu1))+
     &               v_2*fac*DLOG((a+mu2)/(b+mu2))+
     &               v_3*fac*DLOG((a+mu3)/(b+mu3))
            vkk(i,j)=vkk(j,i)
         ENDDO
      ENDDO

      END  SUBROUTINE v_pot_yukawa
      
C     **************************************************
C         Solves eq. V88 
C         and returns < p | f_mtx | n_pole =k>
C     **************************************************

      SUBROUTINE f_mtx_eq(f_mtx,vkk,n_pole)
      USE mesh_variables
      USE constants
      IMPLICIT NONE
      INTEGER i, j, int, istat, n_pole
      DOUBLE PRECISION f_mtx,vkk,dp,deriv,pih,xsum
      DIMENSION dp(1),deriv(1)
      DIMENSION f_mtx(n_rel),vkk(n_rel,n_rel),a(:,:),fu(:),q(:),au(:)
      DOUBLE PRECISION, ALLOCATABLE :: fu, q, au, a

      pih=2.D0/ACOS(-1.D0)
      ALLOCATE( a (n_rel,n_rel), STAT=istat)
      DO i=1,n_rel
         ALLOCATE(fu(n_rel), q(n_rel), au(n_rel), STAT=istat)
         DO j=1,n_rel
            fu(j)=vkk(i,j)-vkk(i,n_pole)*vkk(n_pole,j)/
     &                 vkk(n_pole,n_pole)
         ENDDO
         DO j=1,n_rel
            IF(j /= n_pole ) THEN     ! regular part
               a(j,i)=pih*fu(j)*wra(j)*(ra(j)**2)/
     &                (ra(j)**2-ra(n_pole)**2)
            ELSEIF(j == n_pole) THEN  ! use l'Hopitals rule to get pole term
               dp(1)=ra(j)             
               CALL spls3(ra,fu,n_rel,dp,deriv(1),1,q,au,2,0) 
               a(j,i)=pih*wra(j)*ra(j)/2.d0*deriv(1)
            ENDIF
         ENDDO
         DEALLOCATE(fu, q, au, STAT=istat)   ! free space in heap 
         a(i,i)=a(i,i)+1.D0
      ENDDO
      CALL matinv(a, n_rel)      ! Invert the matrix a
      DO j=1,n_rel               ! multiply inverted matrix a with dim less pot
         xsum=0.D0
         DO i=1,n_rel
            xsum=xsum+a(i,j)*vkk(i,n_pole)/vkk(n_pole,n_pole)  ! gives f-matrix in V88
         ENDDO
         f_mtx(j)=xsum
      ENDDO
      DEALLOCATE (a, STAT=istat)

      END SUBROUTINE f_mtx_eq

C     **************************************************
C         Solves the principal value integral of V89
C         returns the t-matrix for k=k, t_shell
C     **************************************************

      SUBROUTINE principal_value(vkk,f_mtx,n_pole,t_shell) 
      USE mesh_variables
      IMPLICIT NONE
      DOUBLE PRECISION vkk, f_mtx, t_shell, sum, pih, deriv, term
      DIMENSION deriv(1)
      DIMENSION vkk(n_rel, n_rel), f_mtx(n_rel),fu(:), q(:), au(:)
      DOUBLE PRECISION, ALLOCATABLE :: fu, q, au
      INTEGER n_pole, i, istat

      ALLOCATE(fu(n_rel), q(n_rel), au(n_rel), STAT=istat)
      sum=0.D0
      pih=2.D0/ACOS(-1.D0)
      DO i=1,n_rel
         fu(i)=vkk(n_pole,i)*f_mtx(i)
      ENDDO
      DO i=1,n_pole-1  ! integrate up to the pole - 1 mesh
         term=fu(i)*(ra(i)**2)-fu(n_pole)*(ra(n_pole)**2)
         sum=sum+pih*wra(i)*term/(ra(i)**2-ra(n_pole)**2)
      ENDDO       ! here comes the pole part
      CALL spls3(ra,fu,n_rel,ra(n_pole),deriv,1,au,q,2,0)      
      sum=sum+pih*wra(n_pole)*(fu(n_pole)+ra(n_pole)*deriv(1)/2.d0) 
      DO i=n_pole+1,n_rel  ! integrate from pole + 1mesh pt to infty
         term=fu(i)*(ra(i)**2)-fu(n_pole)*(ra(n_pole)**2)
         sum=sum+pih*wra(i)*term/(ra(i)**2-ra(n_pole)**2)
      ENDDO
      t_shell=vkk(n_pole,n_pole)/(1.d0+sum)
      DEALLOCATE (fu, q, au, STAT=istat)

      END SUBROUTINE principal_value

C         ***********************************************
C             Set up of relative mesh and weights
C         ***********************************************

      SUBROUTINE rel_mesh
      USE mesh_variables
      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION pih,u,s,xx,c,h_max
      PARAMETER (c=0.75)
      DIMENSION u(n_rel), s(n_rel)

      pih=ACOS(-1.D0)/2.D0
      CALL gausslegendret (0.D0,1.d0,n_rel,u,s)
      DO i=1,n_rel
         xx=pih*u(i)
         ra(i)=DTAN(xx)*c
         wra(i)=pih*c/DCOS(xx)**2*s(i)
      ENDDO

      END SUBROUTINE rel_mesh

C     *********************************************************
C            Routines to do mtx inversion, from Numerical
C            Recepies, Teukolsky et al. Routines included
C            below are MATINV, LUDCMP and LUBKSB. See chap 2
C            of Numerical Recepies for further details
C            Recoded in FORTRAN 90 by M. Hjorth-Jensen
C     *********************************************************

      SUBROUTINE matinv(a,n)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION a(n,n)
      INTEGER istat
      DOUBLE PRECISION, ALLOCATABLE :: y(:,:)
      INTEGER, ALLOCATABLE :: indx(:)

      ALLOCATE (y( n, n), STAT =istat)
      ALLOCATE ( indx (n), STAT =istat)
      DO i=1,n
         DO j=1,n
            y(i,j)=0.
         ENDDO
         y(i,i)=1.
      ENDDO
      CALL  ludcmp(a,n,indx,d)
      DO j=1,n
         call lubksb(a,n,indx,y(1,j))
      ENDDO
      DO i=1,n
         DO j=1,n
            a(i,j)=y(i,j)
         ENDDO
      ENDDO
      DEALLOCATE ( y, STAT=istat)
      DEALLOCATE ( indx, STAT=istat)

      END SUBROUTINE matinv
 

      SUBROUTINE LUDCMP(A,N,INDX,D)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (TINY=1.0E-20)
      DIMENSION A(N,N),INDX(N)
      INTEGER istat
      DOUBLE PRECISION, ALLOCATABLE :: vv(:)

      ALLOCATE ( vv(n), STAT = istat)
      D=1.
      DO I=1,N
         AAMAX=0.
         DO J=1,N
            IF (ABS(A(I,J)) > AAMAX) AAMAX=ABS(A(I,J))
         ENDDO
         IF (AAMAX == 0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      ENDDO
      DO J=1,N
         IF (J > 1) THEN
            DO I=1,J-1
               SUM=A(I,J)
               IF (I > 1)THEN
                  DO K=1,I-1
                     SUM=SUM-A(I,K)*A(K,J)
                  ENDDO
                  A(I,J)=SUM
               ENDIF
            ENDDO
         ENDIF
         AAMAX=0.
         DO I=J,N
            SUM=A(I,J)
            IF (J > 1)THEN
               DO K=1,J-1
                  SUM=SUM-A(I,K)*A(K,J)
               ENDDO
               A(I,J)=SUM
            ENDIF
            DUM=VV(I)*ABS(SUM)
            IF (DUM >= AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
         ENDDO
         IF (J /= IMAX)THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            ENDDO
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF(J /= N)THEN
            IF(A(J,J) == 0.) A(J,J)=TINY
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            ENDDO
         ENDIF
      ENDDO
      IF(A(N,N) == 0.)  A(N,N)=TINY
      DEALLOCATE ( vv, STAT = istat)

      END SUBROUTINE LUDCMP
 
      SUBROUTINE LUBKSB(A,N,INDX,B)
      implicit real*8(a-h,o-z)
      DIMENSION A(N,N),INDX(N),B(N)
      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II /= 0)THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            ENDDO
         ELSE IF (SUM /= 0.) THEN
            II=I
         ENDIF
         B(I)=SUM
      ENDDO
      DO I=N,1,-1
         SUM=B(I)
         IF (I < N)THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            ENDDO
         ENDIF
         B(I)=SUM/A(I,I)
      ENDDO

      END SUBROUTINE lubksb