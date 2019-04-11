!================================================================
!================================================================
!     Evaluates the reduced radial wave function for the 
!     deuteron using the variational algorithm
!     Returns the binding energy and writes in util/fort.18
!     the variational coefficients.
!     It also prints the radial wave function, the partial wave
!     with l=0 (l=2) is saved to out/deuterons.dat (out/deuterond.dat)
!================================================================
!================================================================


      subroutine deuteron(iiv, binding)
      implicit real*8 (a-h,o-z)
      character*1 jobz,uplo, segno
      

      parameter(nne=200,nnn=2*nne,nnr=50000)

      data htc /197.326968/

      dimension a(nnr),xx(nnr),yy(nnr),fun(nnr)
      dimension lla(2)
      dimension icont(2,nne),ax(nnn,nnn),am(nnn,nnn),
     x          atc(nnn,nnn)
      dimension vv(2,2,nnr),wr(nnn),aux(4*nnn)
      dimension axx(nnn,nnn),amm(nnn,nnn)
      dimension u0(0:nne,nnr),u1(0:nne,nnr),u2(0:nne,nnr)
      dimension v0(2,nne,nnr),v1(2,nne,nnr),v2(2,nne,nnr)
	dimension aondaa(nnr), vvv(nnr)
	dimension v2c(2,2)

	dimension ysol(2,nnr), u0r(nnr), u2r(nnr)
      dimension u0r1(nnr), u0r2(nnr), u2r1(nnr),u2r2(nnr)

      namelist /name/
     x               af,rang,h,
     x               nnl,alf0,alf2,htm, neq, jtot

      common /energia/e0,htm


      af=1.01d0
      rang=90.d0
      h=0.01d0
      nnl=80
      alf0=4.d0
      alf2=4.d0

      am2=939.656
      am4=938.565
      amr=am2*am4/(am2+am4)
      htc = 197.326968d0
      htm=htc**2/(2*amr)
      write(*,*) htm,2*amr/htc**2,2*amr


      jtot=1
      neq=2

      write(*,*) '---------------------------------------'
 12     write(*,*) '---------------------------------------'
      write(6,*) 'htm: ',htm
 1    write(*,name)
      write(*,*) '---------------------------------------'
	


      if(jtot.eq.1.and.neq.eq.1) then
	lla(1)=0
	endif
      if(jtot.eq.1.and.neq.eq.2)then
         lla(1)=0
         lla(2)=2
      endif
      if(jtot.eq.3.and.neq.eq.1) then 
	lla(1)=2
	endif
      if(jtot.eq.3.and.neq.eq.2)then
         lla(1)=2
         lla(2)=4
      endif

	write(*,*) lla(1), lla(2)
      

c     Creazione griglia
      h5=h/22.5d0
      range=rang
      af1 =af-1.d0
      nx=1+dlog(1.d0+range*af1/h)/dlog(af)

      write(*,*) '---------------------------------------'
      if(nx.gt.nnr)then
      write(6,*)'Range troppo grande!'
      write(6,*) range,nx,nnr
      write(*,*) '---------------------------------------'
      go to 12
      endif

      np=nx+1
      range=(af**nx-1.d0)*h/af1


      do i=1,nx
      xx(i)=h/af1*(af**i-1.d0)
      a(i)=dlog(af)/af1*af**i
      enddo

      write(6,1020) nx,range
 1020 format(' Numero punti: ',i5,/,' Range effettivo: ',f15.5)

      nmx=nnl-1


c	Onda s
      do i=1,nx
         yy(i)=alf0*xx(i)
      enddo

      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)

      do il=0,nmx
      anj=dsqrt( alf0**3*dgamma(il+1.d0)/dgamma(il+3.d0) )

      do kl=1,nx
      ud0=anj*dexp(-yy(kl)/2.d0)   
      v0(1,il+1,kl)=ud0*u0(il,kl)
      v1(1,il+1,kl)=alf0*(-0.5d0*v0(1,il+1,kl)+ud0*u1(il,kl) )
      v2(1,il+1,kl)=alf0*(-0.5d0*v1(1,il+1,kl)
     x                     +alf0*ud0*(-0.5d0*u1(il,kl)+u2(il,kl)) )
      enddo

      enddo

      if(neq.eq.2)then



c	Onda d
      do i=1,nx
         yy(i)=alf2*xx(i)
      enddo

      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)

      do il=0,nmx
      anj=dsqrt( alf2**3*dgamma(il+1.d0)/dgamma(il+3.d0) )

      do kl=1,nx
      ud0=anj*dexp(-yy(kl)/2.d0)   
      v0(2,il+1,kl)=ud0*u0(il,kl)
      v1(2,il+1,kl)=alf2*(-0.5d0*v0(2,il+1,kl)+ud0*u1(il,kl) )
      v2(2,il+1,kl)=alf2*(-0.5d0*v1(2,il+1,kl)
     x                     +alf2*ud0*(-0.5d0*u1(il,kl)+u2(il,kl)) )
      enddo

      enddo
	endif


      ii=0
      do i=1,neq
      do k=1,nnl
         ii=ii+1
         icont(i,k)=ii
      enddo
      enddo


      ndim=ii
      write(*,*)'Dimensione matrice', ndim

      do i=1,ndim
      do j=1,ndim
         ax(i,j)=0.d0
         am(i,j)=0.d0
         axx(i,j)=0.d0
         amm(i,j)=0.d0
      enddo
      enddo


c	Norma e energia cinetica

      do il=1,nnl
      do ir=1,nnl

      do isv=1,neq
      do iqv=1,neq

      if(isv.eq.iqv) then

            iss=icont(isv,il)   
            iqq=icont(iqv,ir)   
     
            ax(iss,iqq)=0.d0
            if(il.eq.ir) then
            ax(iss,iqq)=1.d0
            endif
     
c     calcolo numerico dell'energia cinetica
            fun(1)=0.d0
            agr=lla(iqv)*(lla(iqv)+1)
            do kl=1,nx
               rho=xx(kl)
               ww =a(kl)
               fun(kl+1)=v0(isv,il,kl)*(     rho**2*v2(iqv,ir,kl)
     x                                + 2.d0*rho   *v1(iqv,ir,kl)      
     x                                         -agr*v0(iqv,ir,kl) )*ww   
            enddo
            wwtl=b5(1,nx,0,0,h5,0.d0,0.d0,fun,1)

            am(iss,iqq)=-htm*wwtl
            atc(iss,iqq)=-htm*wwtl

      endif

      enddo ! ir
      enddo ! il
      enddo ! iqv
      enddo ! isv


c	Energia potenziale
      

      iii=(iqv+iqm)+30
      do i=1,nx
            r = xx(i)
            call av18pw(1,lla(1),1,jtot,0,1,-1,r,v2c)
            do iqv = 1,2
            do iqm = 1,2
                  vv(iqv,iqm,i)=v2c(iqv,iqm)
            enddo
            enddo
	enddo

c      stop
          
      do isv=1,neq
      do iqv=1,neq

         do il=1,nnl
         iss=icont(isv,il)
         do ir=1,nnl
         iqq=icont(iqv,ir)

         fun(1)=0.d0
         do i=1,nx
         rr=xx(i)
         ww=rr*rr*a(i)
         fun(i+1)=vv(isv,iqv,i)*v0(isv,il,i)*v0(iqv,ir,i)*ww
         enddo
         gun=b5(1,nx,0,0,h5,0.d0,0.d0,fun,1)
         am(iss,iqq)=am(iss,iqq)+gun

         enddo
         enddo

      enddo
      enddo



c     Simmetrizzazione

      do i=1,ndim
      do j=1,ndim

      amm(i,j)=am(i,j)+am(j,i)
      axx(i,j)=ax(i,j)+ax(j,i)

      enddo
      enddo

      jobz='v'
      uplo='u'
      lwork=4*nnn
      call dsygv(1,jobz,uplo,ndim,amm,nnn,axx,nnn,wr,aux,
     x           lwork,info)

      write(*,*)'info',info



      do i=1,ndim
         if(jtot.eq.1)then
            if(wr(i).lt.0.d0) then
               write(*,*)i,wr(i)
               nsol=i
            endif
         else
            write(*,*)i,wr(i)
         endif
      enddo

	indice=2
      binding = -wr(indice)

c!	Stampa della funzione d'onda in onda.dat
c	open(unit=20,file='out/onda.dat',status='unknown'
c     x ,form='formatted')
	do i=1, nx
	aondaa(i)=0.d0
	do j=1, nnl
	aondaa(i) = aondaa(i) + amm(j,2)*v0(1,j,i)
	enddo
	enddo
      close(unit=20)












c**************************************************************************************************************************************
      e0=wr(nsol)
      do i=1,neq
      do il=1,nnl
         j=icont(i,il)
         ysol(i,il)=amm(j,nsol)
      enddo
      enddo
c      
      ann=0.d0
      pws=0.d0
      pwd=0.d0
      tcin=0.d0
      ene=0.d0
      do i=1,neq
      do j=1,neq
c
      do il=1,nnl
      do ir=1,nnl
         iss=icont(i,il)
         iqq=icont(j,ir)
         y12=ysol(i,il)*ysol(j,ir)
         ann=ann+ax(iss,iqq)*y12 ! norma
         if(i.eq.1.and.j.eq.1)
     x        pws=pws+ax(iss,iqq)*y12
            if(i.eq.2.and.j.eq.2)
     x         pwd=pwd+ax(iss,iqq)*y12
         tcin=tcin+y12*(atc(iss,iqq)+atc(iqq,iss))
         ene=ene+y12*(am(iss,iqq)+am(iqq,iss))
      enddo
      enddo
c
      enddo
      enddo
c
      pws=pws/ann
      pwd=pwd/ann
      ene=ene
c      
c
      nmx=nnl-1
c
c     onda s
      do i=1,nx
         yy(i)=alf0*xx(i)
      enddo
      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)
c
      do i=1,nx
      rr=xx(i)
      u0r(i)=0.d0
      dx=dexp(-0.5d0*yy(i))
      do il=0,nmx
         anj=dx*dsqrt( alf0**3*dgamma(il+1.d0)/dgamma(il+3.d0) )
         uu0=anj*u0(il,i)
         u0r(i)=u0r(i)+ysol(1,il+1)*uu0
      enddo
      enddo
c
c     onda d
      do i=1,nx
         yy(i)=alf2*xx(i)
      enddo
      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)
c
      do i=1,nx
      rr=xx(i)
      u2r(i)=0.d0
      dx=dexp(-0.5d0*yy(i))
      do il=0,nmx
         anj=dx*dsqrt( alf2**3*dgamma(il+1.d0)/dgamma(il+3.d0) )
         uu0=anj*u0(il,i)
         u2r(i)=u2r(i)+ysol(2,il+1)*uu0
      enddo
      enddo
c
      do i=1,nx
      rr=xx(i)
      enddo
 
      write(*,1030)100.d0*pws,100.d0*pwd,tcin,ene,e0
      
 1030 format(//,' pws(%), pwd(%), <t>(MeV), <h>(MeV), E(MeV)  =',/,
     x          5d16.8/)
c
c
      do i=1,nx
         yy(i)=alf0*xx(i)
      enddo
      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)
c
      do i=1,nx
      rr=xx(i)
      u0r(i)=0.d0
      u0r1(i)=0.d0
      u0r2(i)=0.d0
      dx=dexp(-0.5d0*yy(i))
      do il=0,nmx
         anj=dx*dsqrt( alf0**3*dgamma(il+1.d0)/dgamma(il+3.d0) )
         uu0=anj*u0(il,i)
         u0r(i)=u0r(i)+ysol(1,il+1)*uu0
         uu1=alf0*(-0.5d0*uu0+anj*u1(il,i))
         u0r1(i)=u0r1(i)+ysol(1,il+1)*uu1
         uu2=alf0*(-0.5d0*uu1
     x              +alf0*anj*(-0.5d0*u1(il,i)+u2(il,i)) )
         u0r2(i)=u0r2(i)+ysol(1,il+1)*uu2
      enddo
      enddo
c
c
      do i=1,nx
         yy(i)=alf2*xx(i)
      enddo
      call plag2(yy,nx,2.d0,nmx,u0,u1,u2)
c
      do i=1,nx
      rr=xx(i)
      u2r(i)=0.d0
      u2r1(i)=0.d0
      u2r2(i)=0.d0
      dx=dexp(-0.5d0*yy(i))
      do il=0,nmx
         anj=dx*dsqrt( alf2**3*dgamma(il+1.d0)/dgamma(il+3.d0) )
         uu0=anj*u0(il,i)
         u2r(i)=u2r(i)+ysol(2,il+1)*uu0
         uu1=alf2*(-0.5d0*uu0+anj*u1(il,i))
         u2r1(i)=u2r1(i)+ysol(2,il+1)*uu1
         uu2=alf2*(-0.5d0*uu1
     x              +alf2*anj*(-0.5d0*u1(il,i)+u2(il,i)) )
         u2r2(i)=u2r2(i)+ysol(2,il+1)*uu2
      enddo
      enddo
c
      fun(1)=0.d0
      do i=1,nx
      rr=xx(i)   
      wrr=rr*rr*a(i)
      fun(i+1)=wrr*(u0r(i)**2+u2r(i)**2)
      enddo
      unu=b5(1,nx,0,0,h5,0.d0,0.d0,fun,1)
      anorma=1./dsqrt(unu)

      open(unit=20,file='out/deuterons.dat',
     x      status='unknown',form='formatted')
      open(unit=21,file='out/deuterond.dat',
     x      status='unknown',form='formatted')

      do i=1,nx
      rr=xx(i)
      u0r(i)=u0r(i)*anorma
      u2r(i)=u2r(i)*anorma
      u0r1(i)=u0r1(i)*anorma
      u2r1(i)=u2r1(i)*anorma
      u0r2(i)=u0r2(i)*anorma
      u2r2(i)=u2r2(i)*anorma

      write(20,*) rr,u0r(i)*rr,u2r(i)*rr
      write(21,*) rr,u2r(i)*rr
      enddo
      close(20)
      close(21)
c
      fun(1)=0.d0
      do i=1,nx
      rr=xx(i)
      wrr=rr*rr*a(i)
      fun(i+1)=wrr*( u0r(i)*(u0r2(i)+2.d0*u0r1(i)/rr)
     x             +u2r(i)*(u2r2(i)+2.d0*u2r1(i)/rr-6.d0*u2r(i)/rr/rr) )
      enddo
      tcin=-htm*b5(1,nx,0,0,h5,0.d0,0.d0,fun,1)
c
      fun(1)=0.d0
      do i=1,nx
      rr=xx(i)
      wrr=rr*rr*a(i)
      fun(i+1)=wrr*( 
     x          u0r(i)*vv(1,1,i)*u0r(i)+u2r(i)*vv(1,2,i)*u0r(i)
     x         +u2r(i)*vv(2,1,i)*u0r(i)+u2r(i)*vv(2,2,i)*u2r(i) )
      enddo
      vpot=b5(1,nx,0,0,h5,0.d0,0.d0,fun,1)
c
      ene=(tcin+vpot)!/unu

c     creo fort.18
      open(unit=18,file='util/fort.18',
     x      status='unknown',form='formatted')
      write(18,*)neq,alf0,nnl
      do i=1,neq
         if(i.eq.1)l2i=0
         if(i.eq.2)l2i=2
         write(18,'(4i2)')l2i,1,1
         do il=1,nnl
         ysol(i,il)=   ysol(i,il)*anorma  
         write(18,'(d20.10)')   ysol(i,il)    
         enddo
      enddo
      write(18,*)e0,tcin,e0-tcin,pws,pwd
      close(unit=18)
c


c
      write(*,*)'valor medio di t,v,e h (MeV) = ',tcin,vpot,ene
	write(*,*) '--------------------------------------------'

	return
	end















	subroutine plag2(xx,nm,apf,n1,u,u1,u2)
c*calcula el polinomio de laguerre en las absisas xx(1,...nm) con numeros
c*cuanticos apf,n con n=0,....n1
      implicit real*8(a-h,o-z)
      parameter(nnr=500,nne=200)
      dimension xx(nnr),u(0:nne,nnr),u1(0:nne,nnr),u2(0:nne,nnr)
      nmax=n1-1
      if(nmax.gt.nne-1)then
         write(*,*)'n1 too large',n1,nne
         stop
      endif
      do 10 i=1,nm
      x=xx(i)
      u(0,i)=1.d0
      u(1,i)=apf+1.d0-x
      u1(0,i)=0.d0
      u1(1,i)=-1.d0
      u2(0,i)=0.d0
      u2(1,i)=0.d0
 10   continue
      do 20 n=1,nmax
      n1=n+1
      d=1.d0/(n+1.d0)
      a1=2.d0*n+apf+1.d0
      a0=n+apf
      a2=n1+apf
      do 11 i=1,nm
      x=xx(i)
      u (n1,i)=((a1-x)*u(n,i)-a0*u(n-1,i))*d
      u1(n1,i)=(n1*u(n1,i)-a2*u(n,i))/x
      u2(n1,i)=(-u(1,i)*u1(n1,i)-n1*u(n1,i))/x
 11   continue
 20   continue

      return
      end
