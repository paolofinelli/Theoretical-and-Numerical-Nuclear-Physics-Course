C======================== varmc.f ==============================

      program nuclei
c **********************************************************************
c main program for ground-state energy and density of 3h - 3he - 4he
c using variational monte carlo
c
c switches: lpot   1=reid
c           ltbp   0=no 3-body potential 1=urbana vii
c           lprt   0=no 1=print out wave function tables
c           lrw    0=no 1=write 2=read random walk
c
c parameters: npart = # of particles
c             nprot = # of protons
c             ns = # of spin states
c             nt = # of ispin states
c             nprs = # of pairs
c             npts = # of grid points in tabulated functions
c             h = grid spacing in tabulated functions
c             nin = i/o unit for input
c             nout = i/o unit for output
c             nrw = i/o unit for random walk file
c             ilenwk = length of temporary working space (integers)
c             iflnwk = length of temporary working space (floating pt)
c
c conventions: complex variables begin with 'c'
c              do loops run over indices: i for spin state #
c                                         j for ispin state #
c                                         k for pair #
c                                         m for space (x,y,z)
c                                         n for particle #
c
c for different calculations set the parameters and formats:
c 3h:  npart=3, nprot=1, nt=3, nfact=6,  write(nout,3...)
c 3he: npart=3, nprot=2, nt=3, nfact=6,  write(nout,3...)
c 4he: npart=4, nprot=2, nt=6, nfact=24, write(nout,4...)
c ----------------------------------------------------------------------
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nin=5,nout=6,nrw=9)
      parameter (ilenwk=2*npart+2*nfact+ns+npart*nt)
      parameter (iflnwk=16*ns*nt+1)
c ----------------------------------------------------------------------
c                 common block variables
c
c     /const/   h2m=hbar**2/m   pi=3.14159...
c     /corr/    tables used to interpolate correlations
c               fc          central correlation
c               us,ut,...   noncentral correlations: s(spin), t(tensor),
c                      ct(isospin), st(spin-isospin), tt(tensor-isospin)
c               rr          r grid
c     /dens/    bins for one-body number densities
c               rhoa        central density bin
c               rhot        isospin density bin
c     /ispin/   isospin matrix operators
c               ntexch(j,k) isospin state resulting from interchanging
c                             isospins of pair k in isospin state j
c               mtval(j,n)  third component of isospin for particle n
c                             in isospin state j
c     /lbp/     pair-particle lookup tables
c               n1p(k) = particle # for first particle of pair k
c               n2p(k) = particle # for second particle of pair k
c               k12p(n1,n2) = pair # for particles n1 & n2 (0 if n1=n2)
c     /nabla/   gradient components and laplacian of wave function
c               cdelr(i,j,m,n) derivative of wave function
c                                (spin state i, isospin state j,
c                                 space component m, particle n)
c               cdlsq(i,j)     del**2 acting on wave function
c                                (summed over all particles)
c     /order/   iordl,iordr    order of pair corrrelation operators
c                                on right and left hand sides
c     /pairs/
c               drh(m,k)       m component unit vector for pair k
c               drm(k)         magnitude of pair k separation
c               fcv,us,..      interpolated value of corr. operators
c               vcv,vsv,...    interpolated value of potential terms
c     /pass/
c               steps          step size used in metropolis walk
c               psi2(2)        new Õpsi2(itrial)å and last accepted
c                                Õpsi2(igood)å values for square of
c                                wave function
c               r(m,n,2)       particle positions, last index is either
c                                igood or itrial
c               fr0,frpl,frmn  weights used for interpolating function
c               ils            grid pt to use for interpolations
c               igood,itrial   pointers to last accepted and trial
c                                values of psi2 and r
c               nmoves         number of moves (metropolis steps) between
c                                energy calculations
c               mok            counter for number of accepted moves
c               iwarn(1)       number of steps rejected for
c                                being within 3 grid pts of origin
c               iwarn(2)       number of steps rejected for
c                                being beyond the grid
c               iwarn(3)       number of times real part of overlap < 0
c                                possible because of sampling of pair
c                                operator orders
c               iop            set to 1 if a move must be rejected for
c                                going outside tabulation range
c     /poten/   tables used to interpolate potentials
c               vc,vs,...   labels c,s,t,ct,st,tt as in /corr/, plus
c                           ls(spin-orbit), lst(spin-orbit-isospin),
c               vypi,vtpi   Yukawa and tensor functions for Vijk
c     /speed/
c               mflops      counter for number of floating pt operations
c     /spin/
c               nsexch(i,k) spin state resulting from interchanging
c                             spins of pair k in spin state i
c               msval(i,n)  third component of spin for particle n
c                             in spin state i
c               nsflip(i,n) spin state resulting from flipping
c                             spin of particle n in spin state i
c     /tbodyc/  three-body non-central correlation parameters
c               t1s,nt2s     coefficient and power for spin operators
c               t1t,nt2t     coefficient and power for tensor operators
c               t3s         exponent (used for both spin and tensor)
c     /wavf/
c               rphi(i,j)   initial uncorrelated wave function ]phi>
c               cpsl(i,j)   correlated wave function on left <psi]
c               cpsr(i,j)   correlated wave function on right ]psi>
c                             (spin state i, isospin state j)
c     /fwksp/   temporary working space
c     /iwksp/   temporary working space
c ----------------------------------------------------------------------
      common /const/ h2m,pi
      common /corr/ fc(npts),us(npts),ut(npts),uct(npts),ust(npts)
     &,utt(npts),rr(npts)
      common /dens/ rhoa(100),rhot(100)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /nabla/ cdelr(ns,nt,3,npart),cdlsq(ns,nt)
      common /order/ iordl(nprs),iordr(nprs)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /poten/ vc(npts),vs(npts),vt(npts),vls(npts)
     &,vct(npts),vst(npts),vtt(npts),vlst(npts),vtpi(npts),vypi(npts)
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /tbodyc/ t1s,t1t,nt2s,nt2t,t3s,t3t
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /fwksp/ dum(iflnwk)
      common /iwksp/ idum(ilenwk)
      real*8 dtot(2,11),dmod(11),etot(2,11),emod(11)
     &,en(11),ed(11),er(11),vn(11),vd(11),vr(11),und(11),unr(11)
     &,udr(11),enx(11),edx(11),erx(11),vnx(11),vdx(11),vrx(11)
     &,undx(11),unrx(11),udrx(11),vv(10)
     &,xrho(100),rhop(100),rhon(100),barp(100),barn(100)
      real*8 esep(4),eta(2),ac(6),aa(6),ar(6),alpha(6),gamma(6)
      character*2 nucnam(2),ncnams
      character*8 nrwfil
      character*20 potnam(1),timdat,tmdats
      character*40 sysdat
      data nucnam/'h ','he'/
      data potnam/'reid soft-core v8'/
      safert(z)=sign(sqrt(abs(z)),z)
c --------------------------------------------------
c start timing with function timer
c get time, date, system info from subroutine header
c --------------------------------------------------
ctss  call link("unit5=(i,open),unit6=(o,create,text)//")
      zero=0.
      time0=timer(zero)
      mflops=0
      call header(sysdat,timdat)
      write(nout,1001) npart,nucnam(nprot),sysdat,timdat
 1001 format(/,80('*')//3x,i1,a2,6x,a40,5x,a20)
c ------------------------------------------------
c read in and set hamiltonian
c lpot = two-body potential (1=reid v8)
c ltbp = three-body potential (0=no, 1=urbana vii)
c lprt = wave function print (0=no, 1=yes)
c subroutine pot generates potential functions
c ------------------------------------------------
      read(nin,*) lpot,ltbp,lprt
      write(nout,1002) potnam(lpot)
 1002 format(/3x,'hamiltonian h = ti + vij + uijk'//3x,'vij : ',a20)
      if (ltbp.eq.1) write(nout,1003)
 1003 format(/3x,'uijk : urbana model vii v2pi + v3nr')
      do 10 l=1,npts
        rr(l)=h*l
        call pot(lpot,rr(l),vv)
        vc(l)=vv(1)
        vs(l)=vv(3)
        vt(l)=vv(5)
        vls(l)=vv(7)
        vct(l)=vv(2)
        vst(l)=vv(4)
        vtt(l)=vv(6)
        vlst(l)=vv(8)
        vtpi(l)=vv(9)
        vypi(l)=vv(10)
   10 continue
c -------------------------------------------------------
c generate spin-isopin vector and two-body correlations
c
c subroutine setspn sets one-body wave function ]phi>
c
c two-body correlation parameters read include:
c   esep(4)  separation energy for 1s0, 1p1, 3s1, 3p2
c   eta(2)   asymptotic tensor/central ratio for 3d1, 3f2
c   ac(6),aa(6),ar(6) parameters in lambda(r)
c   alpha(6) potential quencher
c   gamma(6) lagrange multipliers (initial guess)
c
c subroutine f6cor generates correlations fc,us,...
c
c three-body correlation parameters read include:
c   t1x,nt2x,t3x where x = s(spin-ispin), or t(tensor)
c
c subroutine wfprt prints spin matrices and correlations
c -------------------------------------------------------
      call setspn
      read(nin,*) esep,eta
      read(nin,*) ac
      read(nin,*) aa
      read(nin,*) ar
      read(nin,*) alpha
      read(nin,*) gamma
      write(nout,1005) esep,eta
 1005 format(/3x,'esep(1s0)',4x,'esep(1p1)',4x,'esep(3s1)',4x
     &,'esep(3pj)',5x,'eta(3d1)',5x,'eta(3fj)',/3x,6(f9.3,4x))
      write(nout,1010) ac,aa,ar
 1010 format(3x,'ac',/3x,6(f9.3,4x),/3x,'aa',/3x,6(f9.3,4x)
     &,/3x,'ar',/3x,6(f9.3,4x))
      call f6cor(npart-1,lpot,esep,eta,ac,aa,ar,alpha,gamma
     &          ,fc,us,ut,uct,ust,utt,rr)
      write(nout,1015) alpha
 1015 format (3x,'alpha',/3x,6(f9.3,4x))
      write(nout,1020) gamma
 1020 format (3x,'gamma',/3x,6(f9.3,4x))
      read(nin,*) t1s,t1t,nt2s,nt2t,t3s,t3t
      write(nout,1025) t1s,t1t,nt2s,nt2t,t3s,t3t
 1025 format(3x,'t1s',10x,'t1t',10x,'nt2s',9x,'nt2t',9x,'t3s'
     &,10x,'t3t'/3x,2(f9.3,4x),2(i5,8x),2(f9.3,4x))
      if (lprt.ge.1) call wfprt
c --------------------------------------------------------------
c read random walk parameters
c   dseed = seed for random number generator
c   initmo = number of initial randomizing moves
c   normmo = number of moves between energy samples
c   steps = random move step size
c   ncycls = total number of energies to calculate
c   nprt = number of moves between printing local sums
c   lrw = random walk file switch
c   nrwfil = name for optional random walk file
c
c initialize random number generator, or read data from old walk
c start file if saving new walk
c --------------------------------------------------------------
      read(nin,*) dseed,initmo,normmo,steps,ncycls,nprt
      read(nin,*) lrw,nrwfil
      if (lrw.le.1) then
        write(nout,1050) dseed
 1050   format(/3x,'random number seed ',f10.0)
        call setrnd(dseed)
      else if (lrw.eq.2) then
        open(unit=nrw,file=nrwfil,status='unknown',form='unformatted'
     &      ,err=2001)
        read(nrw) nparts,ncnams,nsamps,tmdats
        write(nout,1055) nparts,ncnams,nrwfil,nsamps,tmdats
 1055   format(/3x,'reading ',i1,a2,' random walk file ',a8,' of',i6,
     &         ' points'/3x,'originally written ',a20)
      end if
      if (lrw.eq.1) then
        open(unit=nrw,file=nrwfil,status='unknown',form='unformatted'
     &      ,err=2002)
        write(nrw) npart,nucnam(nprot),ncycls,timdat
        write(nout,1060) nrwfil
 1060   format(/3x,'random walk file ',a8,' created')
      end if
      timein=timer(time0)
      write(nout,1065) timein
 1065 format(/3x,'wave function set-up time =',f8.3,' seconds')
c -----------------------
c initialize r,igood,psi2
c -----------------------
      one=1.
      pi=acos(-one)
      iwarn(1)=0
      iwarn(2)=0
      iwarn(3)=0
      igood=1
      itrial=2
      do 20 k=1,nprs
        iordl(k)=k
        iordr(k)=k
   20 continue
      do 25 mn=1,3*npart
        r(mn,1,1)=0.
   25 continue
      psi2(1)=0.
      ieee=0
      iee=0
      iprt=nprt
c ---------------------------
c initialize energy sums
c ea,va,uab are local sums
c eax,vax,uabx are total sums
c where a(b)=n,d,r
c ---------------------------
      do 30 l=1,11
        dmod(l)=0.
        en(l)=0.
        ed(l)=0.
        er(l)=0.
        vn(l)=0.
        vd(l)=0.
        vr(l)=0.
        und(l)=0.
        unr(l)=0.
        udr(l)=0.
        enx(l)=0.
        edx(l)=0.
        erx(l)=0.
        vnx(l)=0.
        vdx(l)=0.
        vrx(l)=0.
        undx(l)=0.
        unrx(l)=0.
        udrx(l)=0.
   30 continue
      v3a=0.
      v3c=0.
      v3u=0.
c ----------------------------------------------
c move initmo # of steps & zero density counters
c ----------------------------------------------
      timemv=0.
      if (lrw.le.1) then
        nmoves=initmo
        mok=0
        call moveem
        write(nout,1070) steps,initmo,mok
 1070   format(/3x,'random walk step size = ',f7.4
     &        ,/3x,'of',i5,' initial moves',i4,' were accepted')
        timemv=timer(time0)-timein
        write(nout,1073) timemv
 1073   format(/3x,'initial move time =',f8.3,' seconds')
        nmoves=normmo
        write(nout,1075) normmo
 1075   format(/3x,'metropolis walk attempts',i4,' steps per sample')
      else if (lrw.eq.2) then
        nmoves=1
      end if
      mok=0
      do 35 i=1,100
        rhoa(i)=0.
        rhot(i)=0.
   35 continue
c ----------------------------------------------------
c main loop for monte carlo sampling
c 1) attempt normmo # of moves with subroutine moveem
c or follow old walk with subroutine remove
c to obtain sample configuration with weight <psi]psi>
c 2) calculate expectation values of interest
c repeat ncycls # of times
c ----------------------------------------------------
      do 100 is=1,ncycls
        ieee=ieee+1
        iee=iee+1
        if (lrw.le.1) then
          call moveem
        else if (lrw.eq.2) then
          read(nrw,end=100) igood,iordr,iordl,psi2,r,dmod
          call remove
        end if
c -----------------------------------------------------
c subroutine rhor samples density
c subroutine ekin computes one-body kinetic energy ti
c subroutine vpot computes two-body potential energy
c   v6 = c,s,t,ct,st,tt sum
c   vso = ls,lst sum
c   vcf = coulomb part
c subroutine tbpot computes three-body potential energy
c   v3a = anticommutator part of v2pi3n
c   v3c = commutator part of v2pi3n
c   v3u = short-range repulsive part
c
c emod (dmod) = current (old) energy counter
c fac = 1 (psi2(new)/psi2(old)) for new (old) walk
c
c write data if saving random walk
c -----------------------------------------------------
        call rhor
        call ekin(ti)
        call vpot(v6,vso,vcf)
        if (ltbp.ge.1) call tbpot(v3a,v3c,v3u)
        vij=v6+vso+vcf
        uijk=v3a+v3c+v3u
        fac=1
        if (lrw.eq.2) fac=psi2(igood)/psi2(itrial)
        emod(1)=ti+vij
        emod(2)=emod(1)+uijk
        emod(3)=ti
        emod(4)=vij
        emod(5)=v6
        emod(6)=vso
        emod(7)=vcf
        emod(8)=uijk
        emod(9)=v3a
        emod(10)=v3c
        emod(11)=v3u
        if (lrw.eq.1) write (nrw) igood,iordr,iordl,psi2,r,emod
c ----------------------------------------------------
c summation counters:
c ea = expectation value sum
c va = variance sum
c uab = covariance sum
c a(b) = n (numerator), d (denominator), r (reference)
c ----------------------------------------------------
        do 50 l=1,11
          en(l)=en(l)+emod(l)*fac
          ed(l)=ed(l)+fac
          er(l)=er(l)+dmod(l)
          vn(l)=vn(l)+(emod(l)*fac)**2
          vd(l)=vd(l)+fac**2
          vr(l)=vr(l)+dmod(l)**2
          und(l)=und(l)+emod(l)*fac**2
          unr(l)=unr(l)+emod(l)*fac*dmod(l)
          udr(l)=udr(l)+fac*dmod(l)
   50   continue
        nflops=8+11*19
        if (ieee.lt.iprt) go to 99
c ------------------------------------------
c every nprt # of energies, save total sums,
c compute and print local energies
c etot = energy ( = num/den )
c dtot = energy difference ( = etot-ref )
c ------------------------------------------
        ze1=real(iee)
        ze2=real(ieee)
        do 60 l=1,11
          enx(l)=enx(l)+en(l)
          edx(l)=edx(l)+ed(l)
          erx(l)=erx(l)+er(l)
          vnx(l)=vnx(l)+vn(l)
          vdx(l)=vdx(l)+vd(l)
          vrx(l)=vrx(l)+vr(l)
          undx(l)=undx(l)+und(l)
          unrx(l)=unrx(l)+unr(l)
          udrx(l)=udrx(l)+udr(l)
          en(l)=en(l)/ze1
          ed(l)=ed(l)/ze1
          er(l)=er(l)/ze1
          vn(l)=vn(l)/ze1-en(l)**2
          vd(l)=vd(l)/ze1-ed(l)**2
          vr(l)=vr(l)/ze1-er(l)**2
          und(l)=und(l)/ze1-en(l)*ed(l)
          unr(l)=unr(l)/ze1-en(l)*er(l)
          udr(l)=udr(l)/ze1-ed(l)*er(l)
          etot(1,l)=en(l)/ed(l)
          etot(2,l)=sqrt((vn(l)/ed(l)**2+vd(l)*(en(l)/ed(l)**2)**2
     &                  -2*und(l)*en(l)/ed(l)**3)/ze1)
          dtot(1,l)=en(l)/ed(l)-er(l)
          dtot(2,l)=safert((vn(l)/ed(l)**2+vd(l)*(en(l)/ed(l)**2)**2
     &                    +vr(l)-2*und(l)*en(l)/ed(l)**3-2*unr(l)/ed(l)
     &                    +2*udr(l)*en(l)/ed(l)**2)/ze1)
   60   continue
        nflops=nflops+11*72
c ----------------
c print local sums
c ----------------
        write(nout,1080)
 1080   format(/,80('*'))
        iee=ieee-iee
        write(nout,1100) iee+1,ieee
 1100   format(/10x,'local average for points',i6,' - ',i6)
        write(nout,1140) etot(1,1),etot(2,1)
 1140   format(/5x,'ti + vij = ',f10.5,' +/- ',f8.5)
        if (ltbp.ge.1) write(nout,1145) etot(1,2),etot(2,2)
 1145   format(/5x,'+ uijk   = ',f10.5,' +/- ',f8.5)
        write(nout,1150) (etot(1,j),etot(2,j),j=3,7)
 1150   format(/5x,'   ti=',f10.5,' +/- ',f8.5
     &         ,5x,'  vij=',f10.5,' +/- ',f8.5
     &         /5x,'   v6=',f10.5,' +/- ',f8.5
     &         ,5x,'  vso=',f10.5,' +/- ',f8.5
     &         /5x,'  vcf=',f10.5,' +/- ',f8.5)
        if (ltbp.ge.1) write(nout,1160) (etot(1,j),etot(2,j),j=8,11)
 1160   format(/5x,' uijk=',f10.5,' +/- ',f8.5
     &         ,5x,'  v3a=',f10.5,' +/- ',f8.5
     &         /5x,'  v3c=',f10.5,' +/- ',f8.5
     &         ,5x,'  v3u=',f10.5,' +/- ',f8.5)
        if (lrw.eq.2) then
          write(nout,1165) iee+1,ieee
 1165     format(/10x,'local differences for points',i6,' - ',i6)
          write(nout,1166) dtot(1,1),dtot(2,1)
 1166     format(/5x,'delta( ti + vij ) = ',f10.5,' +/- ',f8.5)
          if (ltbp.ge.1) write(nout,1167) dtot(1,2),dtot(2,2)
 1167     format(/5x,'delta(  + uijk  ) = ',f10.5,' +/- ',f8.5)
        end if
c ----------------------
c compute total energies
c ----------------------
        do 70 l=1,11
          en(l)=enx(l)/ze2
          ed(l)=edx(l)/ze2
          er(l)=erx(l)/ze2
          vn(l)=vnx(l)/ze2-en(l)**2
          vd(l)=vdx(l)/ze2-ed(l)**2
          vr(l)=vrx(l)/ze2-er(l)**2
          und(l)=undx(l)/ze2-en(l)*ed(l)
          unr(l)=unrx(l)/ze2-en(l)*er(l)
          udr(l)=udrx(l)/ze2-ed(l)*er(l)
          etot(1,l)=en(l)/ed(l)
          etot(2,l)=sqrt((vn(l)/ed(l)**2+vd(l)*(en(l)/ed(l)**2)**2
     &                  -2*und(l)*en(l)/ed(l)**3)/ze2)
          dtot(1,l)=en(l)/ed(l)-er(l)
          dtot(2,l)=safert((vn(l)/ed(l)**2+vd(l)*(en(l)/ed(l)**2)**2
     &                    +vr(l)-2*und(l)*en(l)/ed(l)**3-2*unr(l)/ed(l)
     &                    +2*udr(l)*en(l)/ed(l)**2)/ze2)
   70   continue
        nflops=nflops+11*63
        sqrtn=1./sqrt(ze2)
c ----------------
c print total sums
c ----------------
        write(nout,1080)
        write(nout,1200) ieee
 1200   format(/10x,'total average after',i6,' points')
        write(nout,1140) etot(1,1),etot(2,1)
        if (ltbp.ge.1) write(nout,1145) etot(1,2),etot(2,2)
        write(nout,1150) (etot(1,j),etot(2,j),j=3,7)
        if (ltbp.ge.1) write(nout,1160) (etot(1,j),etot(2,j),j=8,11)
        if (lrw.eq.2) then
          write(nout,1202) ieee
 1202     format(/10x,'total differences after',i6,' points')
          write(nout,1166) dtot(1,1),dtot(2,1)
          if (ltbp.ge.1) write(nout,1167) dtot(1,2),dtot(2,2)
        end if
        write(nout,1205) mok,iwarn
 1205   format(/10x,'mok=',i8,10x,'iwarn=',4i8)
        write(nout,1080)
c --------------------
c reset local counters
c --------------------
        do 80 l=1,11
          en(l)=0.
          ed(l)=0.
          er(l)=0.
          vn(l)=0.
          vd(l)=0.
          vr(l)=0.
          und(l)=0.
          unr(l)=0.
          udr(l)=0.
   80   continue
        iee=0
        iprt=iprt+nprt
   99   mflops=mflops+nflops
  100 continue
c --------------------------------
c sum & print densities, rms radii
c pl denotes proton densities
c mn denotes neutron densities
c --------------------------------
      r4pl=0.
      r4mn=0.
      r2pl=0.
      r2mn=0.
      r0pl=0.
      r0mn=0.
      xnconf=1./(ncycls*nmoves)
      xl=0.
      do 120 i=1,100
        xrho(i)=xl+.05
        xu=xl+.1
        x2=xrho(i)*xrho(i)
        x4=x2*x2
        xvol=3./(4.*pi*(xu**3-xl**3))
        rhopl=.5*(rhoa(i)+rhot(i))
        rhomn=.5*(rhoa(i)-rhot(i))
        r4pl=r4pl+rhopl*x4
        r4mn=r4mn+rhomn*x4
        r2pl=r2pl+rhopl*x2
        r2mn=r2mn+rhomn*x2
        r0pl=r0pl+rhopl
        r0mn=r0mn+rhomn
        rhop(i)=rhopl*xnconf*xvol
        rhon(i)=rhomn*xnconf*xvol
        barp(i)=safert(rhop(i)*xnconf*xvol)
        barn(i)=safert(rhon(i)*xnconf*xvol)
        xl=xl+.1
  120 continue
      if (lrw.le.1) then
        write(nout,1300)
 1300   format(5x,'r',14x,'rho p',18x,'rho n'/)
        write(nout,1310) (xrho(i),rhop(i),barp(i)
     &                   ,rhon(i),barn(i),i=1,100)
 1310   format(3x,f5.2,6x,f6.4,' +/- ',f6.4,6x,f6.4,' +/- ',f6.4)
      end if
      r2pl=r2pl/r0pl
      r2mn=r2mn/r0mn
      var2pl=r4pl/r0pl-r2pl**2
      var2mn=r4mn/r0mn-r2mn**2
      rmspl=sqrt(r2pl)
      rmsmn=sqrt(r2mn)
      barmsp=.5*sqrt(var2pl*xnconf)/rmspl
      barmsm=.5*sqrt(var2mn*xnconf)/rmsmn
      write(nout,1350) rmspl,barmsp,rmsmn,barmsm
 1350 format(/10x,'rmspl=',f10.4,' +/- ',f6.4,
     &       /10x,'rmsmn=',f10.4,' +/- ',f6.4)
c -------------------------
c close files, stop timing,
c print performance data
c -------------------------
      if (lrw.ge.1) close(unit=nrw,status='keep')
      timeto=timer(time0)
      timeen=timeto-timemv-timein
      flops=mflops/(1e6*timeen)
      write(nout,2000) timeen,mflops,flops,timeto
 2000 format(/3x,'energy calculation time =',f8.3,' seconds'
     &,//3x,'total floating point operations =',i12
     &,//3x,'energy calculation speed =',f8.3,' mflops'
     &,//3x,'total job time =',f8.3,' seconds')
      call exit(0)
 2001 stop 2001
 2002 stop 2002
      end
c id* timer ************************************************************
c function for elapsed cpu time
c cray denotes cray-specific functions
c cibm denotes ibm-specific functions
c csun denotes sun-specific functions
c cvax denotes vax-specific functions
c ----------------------------------------------------------------------
      function timer(arg)
      implicit real*8 (a-h,o-z)
csun  real*4 sec(2)
      if (arg.eq.0.) then
        timer=1.e-6
cibm    call cputime(time,ircode)
csun    call etime(sec)
cvax    call lib$init_timer
        call lib$init_timer
      else
cray    timer=second(0)-arg
cibm    call cputime(time,ircode)
cibm    timer=1.e-6*time-arg
csun    call etime(sec)
csun    timer=sec(1)-arg
cvax    call lib$stat_timer(2,itime)
cvax    timer=.01*float(itime)
        call lib$stat_timer(2,itime)
        timer=.01*float(itime)
      end if
      return
      end
c id* header ***********************************************************
c subroutine for time, date, and system information
c cray denotes cray-specific functions
c cibm denotes ibm-specific functions
c csun denotes sun-specific functions
c cvax denotes vax-specific functions
c ----------------------------------------------------------------------
      subroutine header(sysdat,timdat)
      implicit real*8 (a-h,o-z)
      character*20 fordat,macdat,timdat
      character*40 sysdat
      character*9 thedat,thetim
cibm  integer*4 idtx(14)
csun  integer*4 idtx(3)
cray  data macdat/'cray x-mp/48'/,fordat/'cft77 3.1.1.5'/
cray  data macdat/'cray-2s/4-128'/,fordat/'cft77 3.1'/
cibm  data macdat/'ibm 3033'/,fordat/'vs fortran 2.3.0'/
cibm  data macdat/'ibm 3090/VF'/,fordat/'vs fortran 2.3.0'/
csun  data macdat/'sun-3/160'/,fordat/'f77'/
csun  data macdat/'sun-4/280'/,fordat/'f77'/
csun  data macdat/'sun sparcstation1'/,fordat/'f77'/
cvax  data macdat/'vax 11/780'/,fordat/'vax fortran v5.3-50'/
cvax  data macdat/'vaxstation 3100'/,fordat/'vax fortran v5.2-33'/
      data macdat/'vax 8700'/,fordat/'vax fortran v5.2-33'/
      sysdat=macdat//fordat
cray  write(thetim,1) clock()
cray  write(thedat,1) date()
cray1 format(a8)
cibm  call datimx(idtx)
cibm  write(thetim,1) idtx(5),idtx(4),idtx(3)
cibm1 format(i2,':',i2,':',i2)
cibm  write(thedat,2) idtx(7),idtx(6),idtx(14)
cibm2 format(i2,'/',i2,'/',i2)
csun  call itime(idtx)
csun  write(thetim,1) idtx
csun1 format(i2,':',i2,':',i2)
csun  call idate(idtx)
csun  idtx(3)=idtx(3)-1900
csun  write(thedat,2) idtx
csun2 format(i2,'/',i2,'/',i2)
      call time(thetim)
      call date(thedat)
      timdat=thetim//'  '//thedat
      return
      end
c id* pot **************************************************************
c subroutine for potential
c lpot: 1=reid v8
c sign(lpot): + return v(r)
c             - return r*v(r)
c ----------------------------------------------------------------------
      subroutine pot(lpot,rr,vv)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension vv(10)
      u=.7
      x=u*rr
      y1=exp(-x)/u
      y2=exp(-2.*x)/u
      y3=exp(-3.*x)/u
      y4=exp(-4.*x)/u
      y6=exp(-6.*x)/u
      y7=exp(-7.*x)/u
      if (rr.eq.0) then
        yr=23.5
      else
        yr=(1.+3./x+3./x**2)*y1-(12./x+3./x**2)*y4
      end if
      hr=10.463
      vv(1)=-19.874*y2+135.21*y3-1432.3*y4+4196.4*y6+1215.8*y7
      vv(2)=19.874*y2-135.21*y3+319.52*y4-1082.3*y6+405.3*y7
      vv(3)=46.241*y2-135.21*y3-64.78*y4+1398.8*y6-1215.8*y7
      vv(4)=(hr/3.)*y1-46.241*y2+135.21*y3+244.06*y4-360.76*y6-405.3*y7
      vv(5)=-26.194*y3+87.943*y4-418.38*y6
      vv(6)=(hr/3.)*yr-8.731*y3-87.943*y4+418.38*y6
      vv(7)=177.23*y4-2233.9*y6
      vv(8)=-177.23*y4+159.75*y6
c ---------------
      if (rr.eq.0.) then
        ypi=0.
        tpi=0.
      else if (rr.gt.4.) then
        ypi=exp(-x)/u
        tpi=(1.+3./x+3./x**2)*ypi
      else
        rcut=1.-exp(-2.*rr*rr)
        ypi=exp(-x)*rcut/u
        tpi=(1.+3./x+3./x**2)*ypi*rcut
      end if
      vv(9)=tpi
      vv(10)=ypi
c ---------------
      if (lpot.lt.0) return
      do 10 l=1,10
        vv(l)=vv(l)/rr
   10 continue
      return
      end
c id* setspn ***********************************************************
c subroutine for construction of particle, spin, ispin lookup tables
c & initial wave function ]phi>
c n1p(k) = particle # for first particle of pair k
c n2p(k) = particle # for second particle of pair k
c k12p(n1,n2) = pair # for particles n1 & n2 (0 if n1=n2)
c msval(i,n) = spin (+1 or -1) of particle n in spin state i
c mtval(j,n) = ispin (+1 or -1) of particle n in ispin state j
c nsexch(i,k) = spin state # on spin exchange of pair k in spin state i
c ntexch(j,k) = ispin state # on ispin exchange of pair k in ispin state
c nsflip(i,n) = spin state # on spin flip of particle n in spin state i
c rphi(i,j) = uncorrelated wave function ]phi> in spin-ispin state i,j
c ----------------------------------------------------------------------
      subroutine setspn
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /iwksp/ nm1p2(npart),istr(nfact),jstr(nfact)
     &,mapbe(ns),mapeb(npart,nt),lpoint(npart)
c --------------------------
c generate n1p, n2p and k12p
c --------------------------
      do 5 n=1,npart
        lpoint(n)=1
        k12p(n,n)=0
    5 continue
      k=1
      do 10 n1=1,npart-1
      do 10 n2=n1+1,npart
        n1p(k)=n1
        n2p(k)=n2
        k12p(n1,n2)=k
        k12p(n2,n1)=-k
        k=k+1
   10 continue
      do 15 n=1,npart
        nm1p2(n)=2**(n-1)
   15 continue
c ----------------------
c generate msval & mtval
c ----------------------
      do 18 n=1,npart
        msval(1,n)=-1
   18 continue
      do 30 i=1,ns-1
        mdsum=0
        intrm=i
        do 20 n=npart,1,-1
          mds=mod(intrm/nm1p2(n),2)
          mdsum=mds+mdsum
          intrm=intrm-mds*nm1p2(n)
          msval(i+1,n)=2*mds-1
   20   continue
        mapeb(mdsum,lpoint(mdsum))=i+1
        mapbe(i+1)=lpoint(mdsum)
        lpoint(mdsum)=lpoint(mdsum)+1
   30 continue
      do 32 n=1,npart
      do 32 j=1,nt
        mtval(j,n)=msval(mapeb(nprot,j),n)
   32 continue
c ----------------------------
c generate the exchange tables
c ----------------------------
      do 39 k=1,nprs
      do 39 i=1,ns
        call exch(i,nm1p2(n1p(k)),nm1p2(n2p(k)),ne)
        nsexch(i,k)=ne
   39 continue
      do 40 k=1,nprs
      do 40 j=1,nt
        call exch(mapeb(nprot,j),nm1p2(n1p(k)),nm1p2(n2p(k)),ne)
        ntexch(j,k)=mapbe(ne)
   40 continue
      do 45 n=1,npart
      do 45 i=1,ns
        n2=mod(((i-1)/nm1p2(n)),2)
        nsflip(i,n)=i+(-2*n2+1)*nm1p2(n)
   45 continue
c -------------------
c generate rphi table
c -------------------
      do 50 ij=1,ns*nt
        rphi(ij,1)=0.
   50 continue
c --------------------
c set first rphi state
c --------------------
      istr(1)=2*nprot+2
      jstr(1)=1
      rphi(istr(1),jstr(1))=1.
c --------------------------------
c generate rphi table by exchanges
c --------------------------------
      np=1
      do 80 n2=2,npart
        nnp=np
        do 70 n1=1,n2-1
          k=k12p(n1,n2)
          do 60 m=1,np
            i=istr(m)
            j=jstr(m)
            ie=nsexch(i,k)
            je=ntexch(j,k)
            rphi(ie,je)=-rphi(i,j)
            istr(m+nnp)=ie
            jstr(m+nnp)=je
   60     continue
          nnp=nnp+np
   70   continue
        np=nnp
   80 continue
      return
      end
c id* exch *************************************************************
c subroutine for finding exchange indices
c ----------------------------------------------------------------------
      subroutine exch(n,i,j,ne)
      implicit integer*4 (i-n)
      ni=mod((n-1)/i,2)
      nj=mod((n-1)/j,2)
      nd=(ni-nj)**2
      mi=nd*(-2*ni+1)
      mj=nd*(-2*nj+1)
      ne=n+mi*i+mj*j
      return
      end
c id* f6cor ************************************************************
c subroutine for finding few-body correlation functions for v8 problem
c
c differential equations solved iteratively in spin-isospin channels
c from inside out to point match using numerov method, then extended
c to point npts by asymptotic form, and projected to operator form
c
c nclstr = # of particles for wave function separation condition
c lpot = potential switch (1=reid v8)
c esep = separation energy for 1s0, 1p1, 3s1, 3p2
c eta = asymptotic tensor ratio for 3d1, 3f2
c ac,aa,ar: cutoffs for long- & short-range lambda(r)
c   lr=1-exp(-(r/ac)**2),   sr=1/(1+exp((r-ar)/aa))
c alpha = noncentral potential quencher
c gamma = lagrange multiplier for solving schroedinger equation
c   (initial guess required for convergence)
c fc = central correlation
c ux = noncentral correlations, where x = s(spin),t(tensor),
c                 ct(ispin),st(spin-ispin),tt(tensor-ispin)
c rr = r grid
c ----------------------------------------------------------------------
      subroutine f6cor(nclstr,lpot,esep,eta,ac,aa,ar,alpha,gamma
     &                ,fc,us,ut,uct,ust,utt,rr)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npts=750,h=.02)
      parameter (match=500,matchp=match+1,matchm=match-1)
      common /const/ h2m,pi
      dimension esep(4),eta(2),ac(6),aa(6),ar(6),alpha(6),gamma(6)
     &,fc(npts),us(npts),ut(npts),uct(npts),ust(npts),utt(npts)
     &,rr(npts)
      real*8 chi(matchp),kappa(4),kk(4)
     &,laml(matchp,6),lams(matchp,6),lam1,lam2,mh2
     &,phi0(6),phim(6),phip(6),psi(npts,6),ri(npts),rsi(matchp)
     &,rv(6),vv(10),vbar(matchp,10)
      h2m=41.47
      mh2=1/h2m
      dhs=h*h/12
      al=1/real(nclstr)
      do 10 l=1,4
        kk(l)=sqrt(2.*mh2*real(nclstr)*esep(l)/real(nclstr+1))
        kappa(l)=kk(l)*al
   10 continue
c ---------------------------
c set up functions
c   psi = fst*r**(l+1)
c   lams = short-range lambda
c   laml = long-range lambda
c   vbar = quenched potential
c   phix = match values for psi
c   rv = r*vbar at r=0
c ---------------------------
      do 30 k=1,2
        l=k+2
        m=k+4
        do 20 i=1,matchp
          ri(i)=1/rr(i)
          rsi(i)=ri(i)**2
          psi(i,k)=exp(-kappa(k)*rr(i))*rr(i)**k*ri(i)**al
          lams(i,k)=1/(1+exp((rr(i)-ar(k))/aa(k)))
          laml(i,k)=h2m*(kappa(k)**2+2*kappa(k)*(al-k)*ri(i)
     &             +al*(al-2*k+1)*rsi(i))*(1-exp(-(rr(i)/ac(k))**2))
          psi(i,l)=exp(-kappa(l)*rr(i))*rr(i)**k*ri(i)**al
          lams(i,l)=1/(1+exp((rr(i)-ar(l))/aa(l)))
          tt=1+3/(kk(l)*rr(i))+3/(kk(l)*rr(i))**2
          ttp=-3*(1/(kk(l)*rr(i))+2/(kk(l)*rr(i))**2)*ri(i)
          ttdp=6*(1/(kk(l)*rr(i))+3/(kk(l)*rr(i))**2)*rsi(i)
          cut=1-exp(-(rr(i)/ac(m))**2)
          cutp=2*rr(i)*(1-cut)/ac(m)**2
          cutdp=(1-2*(rr(i)/ac(m))**2)*cutp*ri(i)
          ten=tt*cut
          tenp=ttp*cut+tt*cutp
          tendp=ttdp*cut+2*ttp*cutp+tt*cutdp
          lam1=kappa(l)**2+2*kappa(l)*(al-k)*ri(i)
     &        +al*(al-2*k+1)*rsi(i)
          lam2=tendp-2*tenp*(kappa(l)+(al-k)*ri(i))+ten*lam1
          den=1-2*eta(k)*ten-8*(eta(k)*ten)**2
          laml(i,l)=h2m*((1-2*eta(k)*ten)*lam1-8*eta(k)**2*ten
     &             *(lam2-6*ten*rsi(i)))*(1-exp(-(rr(i)/ac(l))**2))/den
          psi(i,m)=eta(k)*ten*psi(i,l)
          lams(i,m)=1/(1+exp((rr(i)-ar(m))/aa(m)))
          laml(i,m)=h2m*eta(k)*(lam2-ten*(lam1+6*rsi(i)))
     &             *(1-exp(-(rr(i)/ac(l))**2))/den
   20   continue
   30 continue
      do 40 i=1,matchp
        xx=rr(i)
        call pot(lpot,xx,vv)
        vbar(i,1)=alpha(1)*vv(1)+alpha(2)*(vv(2)-3*(vv(3)+vv(4)))
        vbar(i,2)=alpha(1)*vv(1)-alpha(2)*(3*(vv(2)+vv(3)-3*vv(4)))
        vbar(i,3)=alpha(1)*vv(1)+alpha(2)*(vv(3)-3*(vv(2)+vv(4)))
        vbar(i,4)=alpha(1)*vv(1)+alpha(2)*(vv(2)+vv(3)+vv(4))
        vbar(i,5)=alpha(3)*(vv(5)-3*vv(6))
        vbar(i,6)=alpha(3)*(vv(5)+vv(6))
        vbar(i,7)=alpha(4)*(vv(7)-3*vv(8))
        vbar(i,8)=alpha(4)*(vv(7)+vv(8))
   40 continue
      do 60 k=1,2
        l=k+2
        m=k+4
        do 50 i=matchp,npts
          ri(i)=1/rr(i)
          psi(i,k)=exp(-kappa(k)*rr(i))*rr(i)**k*ri(i)**al
          psi(i,l)=exp(-kappa(l)*rr(i))*rr(i)**k*ri(i)**al
          tt=1+3/(kk(l)*rr(i))+3/(kk(l)*rr(i))**2
          cut=1-exp(-(rr(i)/ac(m))**2)
          ten=tt*cut
          psi(i,m)=eta(k)*ten*psi(i,l)
   50   continue
   60 continue
      do 70 k=1,6
        phim(k)=psi(matchm,k)
        phi0(k)=psi(match,k)
        phip(k)=psi(matchp,k)
   70 continue
      xx=0
      call pot(-lpot,xx,vv)
      rv(1)=vv(1)+vv(2)-3*(vv(3)+vv(4))
      rv(3)=vv(1)+vv(3)-3*(vv(2)+vv(4))
      rv(5)=vv(5)-3*vv(6)
      small=1.e-10
c --------------------------------------
c solve differential equations for psi
c k = 1 (2) solves 1s0 & 3s1 (1p1 & 3p2)
c --------------------------------------
      do 500 k=1,2
        lsq=2*(k-1)
c ---------------------------
c single-channel psi equation
c ---------------------------
        do 140 loop=1,20
          psim=0
          if (k.eq.1) then
            gpsim=dhs*mh2*rv(k)
          else if (k.eq.2) then
            gpsim=2*dhs
          end if
          psi0=psi(1,k)
          gpsi0=dhs*(mh2*(vbar(1,k)+gamma(k)*lams(1,k)+laml(1,k))
     &              +lsq*rsi(1))*psi0
          do 120 j=2,matchp
            gp=dhs*(mh2*(vbar(j,k)+gamma(k)*lams(j,k)+laml(j,k))
     &        +lsq*rsi(j))
            psi(j,k)=(2*psi0+10*gpsi0-psim+gpsim)/(1-gp)
            psim=psi0
            gpsim=gpsi0
            psi0=psi(j,k)
            gpsi0=gp*psi0
  120     continue
          dldif=(psi(matchp,k)-psi(matchm,k))/psi(match,k)
     &         -(phip(k)-phim(k))/phi0(k)
          if (loop.eq.1) then
            dldifo=dldif
            gammao=gamma(k)
            gamma(k)=gammao+1
          else
            if (abs(dldif).le.small) go to 150
            gamman=(dldifo*gamma(k)-dldif*gammao)/(dldifo-dldif)
            dldifo=dldif
            gammao=gamma(k)
            gamma(k)=gamman
          end if
  140   continue
  150   fac=phi0(k)/psi(match,k)
        do 190 j=1,matchp
          psi(j,k)=fac*psi(j,k)
  190   continue
c -----------------------------
c coupled-channel psi equations
c -----------------------------
        l=k+2
        lt=k+4
        lb=k+6
        do 400 iloop=1,20
c ---------------
c central channel
c ---------------
          do 240 loop=1,20
            psim=0.
            chim=0.
            fcm=3*psi(1,l)*ri(1)**k-3*psi(2,l)*ri(2)**k
     &           +psi(3,l)*ri(3)**k
            if (l.eq.3) then
              gpsim=dhs*mh2*rv(l)*fcm
            else if (l.eq.4) then
              gpsim=2*dhs*fcm
            end if
            gchim=gpsim
            hm=0.
            psi0=psi(1,l)
            chi0=psi(1,l)
            chi(1)=chi0
            gpsi0=dhs*(mh2*(vbar(1,l)+gamma(l)*lams(1,l)+laml(1,l))
     &           +lsq*rsi(1))*psi0
            gchi0=dhs*(mh2*(vbar(1,l)+gamma(l)*lams(1,l)+laml(1,l))
     &           +lsq*rsi(1))*chi0
            h0=dhs*mh2*(8*(vbar(1,lt)+gamma(lt)*lams(1,lt)
     &                 +laml(1,lt)))*psi(1,lt)
            do 220 j=2,matchp
              gp=dhs*(mh2*(vbar(j,l)+gamma(l)*lams(j,l)+laml(j,l))
     &               +lsq*rsi(j))
              hp=dhs*mh2*(8*(vbar(j,lt)+gamma(lt)*lams(j,lt)
     &                   +laml(j,lt)))*psi(j,lt)
              psi(j,l)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(j,l)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
  220       continue
            fac=(phi0(l)-psi(match,l))/chi(match)
            do 230 j=1,matchp
              psi(j,l)=psi(j,l)+fac*chi(j)
  230       continue
            dldif=(psi(matchp,l)-psi(matchm,l))/psi(match,l)
     &           -(phip(l)-phim(l))/phi0(l)
            if (loop.eq.1) then
              dldifo=dldif
              gammao=gamma(l)
              gamma(l)=gammao+1
            else
              if (abs(dldif).le.small) go to 300
              gamman=(dldifo*gamma(l)-dldif*gammao)/(dldifo-dldif)
              dldifo=dldif
              gammao=gamma(l)
              gamma(l)=gamman
            end if
  240     continue
c --------------
c tensor channel
c --------------
  300     do 340 loop=1,20
            psim=0.
            chim=0.
            gpsim=0.
            gchim=0.
            fcm=3*psi(1,l)/rr(1)**k-3*psi(2,l)/rr(2)**k
     &           +psi(3,l)/rr(3)**k
            if (l.eq.3) then
              hm=dhs*mh2*rv(lt)*fcm
            else if (l.eq.4) then
              hm=0
            end if
            psi0=psi(1,lt)
            chi0=psi(1,lt)
            chi(1)=chi0
            gpsi0=dhs*((6+lsq)*rsi(1)+mh2*(vbar(1,l)+gamma(l)*lams(1,l)
     &           +laml(1,l)-2*(vbar(1,lt)+gamma(lt)*lams(1,lt)
     &           +laml(1,lt))-3*vbar(1,lb)))*psi0
            gchi0=dhs*((6+lsq)*rsi(1)+mh2*(vbar(1,l)+gamma(l)*lams(1,l)
     &           +laml(1,l)-2*(vbar(1,lt)+gamma(lt)*lams(1,lt)
     &           +laml(1,lt))-3*vbar(1,lb)))*chi0
            h0=dhs*mh2*(vbar(1,lt)+gamma(lt)*lams(1,lt)+laml(1,lt))
     &        *psi(1,l)
            do 320 j=2,matchp
              gp=dhs*((6+lsq)*rsi(j)+mh2*(vbar(j,l)+gamma(l)*lams(j,l)
     &          +laml(j,l)-2*(vbar(j,lt)+gamma(lt)*lams(j,lt)
     &          +laml(j,lt))-3*vbar(j,lb)))
              hp=dhs*mh2*(vbar(j,lt)+gamma(lt)*lams(j,lt)+laml(j,lt))
     &          *psi(j,l)
              psi(j,lt)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(j,lt)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
  320       continue
            fac=(phi0(lt)-psi(match,lt))/chi(match)
            do 330 j=1,matchp
              psi(j,lt)=psi(j,lt)+fac*chi(j)
  330       continue
            dldif=(psi(matchp,lt)-psi(matchm,lt))/psi(match,lt)
     &           -(phip(lt)-phim(lt))/phi0(lt)
            if (loop.eq.1) then
              dldifo=dldif
              gammao=gamma(lt)
              gamma(lt)=gammao-1
            else
              if (abs(dldif).le.small) go to 400
              gamman=(dldifo*gamma(lt)-dldif*gammao)/(dldifo-dldif)
              dldifo=dldif
              gammao=gamma(lt)
              gamma(lt)=gamman
            end if
  340     continue
  400   continue
  500 continue
c ------------------------
c project out fc, us, etc.
c ------------------------
      fcmax=0.
      do 600 i=1,npts
        fc(i)=.0625*( 3*psi(i,1)+psi(i,2)*ri(i)+3*psi(i,3)
     &               +9*psi(i,4)*ri(i))*ri(i)
        uct(i)=.0625*(  psi(i,1)-psi(i,2)*ri(i)-3*psi(i,3)
     &               +3*psi(i,4)*ri(i))*ri(i)
        us(i)=.0625*(-3*psi(i,1)-psi(i,2)*ri(i)+  psi(i,3)
     &               +3*psi(i,4)*ri(i))*ri(i)
        ust(i)=.0625*( -psi(i,1)+psi(i,2)*ri(i)-  psi(i,3)
     &               +  psi(i,4)*ri(i))*ri(i)
        ut(i)= .25*( psi(i,5)+3*psi(i,6)*ri(i))*ri(i)
        utt(i)=.25*(-psi(i,5)+  psi(i,6)*ri(i))*ri(i)
        fcmax=max(fcmax,fc(i))
  600 continue
      do 620 i=1,npts
        us(i)=us(i)/fc(i)
        ut(i)=ut(i)/fc(i)
        uct(i)=uct(i)/fc(i)
        ust(i)=ust(i)/fc(i)
        utt(i)=utt(i)/fc(i)
        fc(i)=fc(i)/fcmax
  620 continue
      return
      end
c id* wfprt ************************************************************
c optional print of wave function info
c for 3h-3he use write(nout,3...)
c for    4he use write(nout,4...)
c ----------------------------------------------------------------------
      subroutine wfprt
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nin=5,nout=6)
      common /corr/ fc(npts),uct(npts),us(npts),ust(npts),ut(npts)
     &,utt(npts),rr(npts)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      write(nout,100)
  100 format(/1x,80('*'),/1x)
      write(nout,3010) (k,k=1,nprs)
 3010 format(/' nprs ',3i5/6x,3('_____'))
 4010 format(/' nprs ',6i5/6x,6('_____'))
      write(nout,3020) n1p
 3020 format(' n1p ]',3i5)
 4020 format(' n1p ]',6i5)
      write(nout,3030) n2p
 3030 format(' n2p ]',3i5)
 4030 format(' n2p ]',6i5)
      write(nout,3040) (n,n=1,npart)
 3040 format(/' k12p'/' n1]n2',3i5/6x,3('_____'))
 4040 format(/' k12p'/' n1]n2',4i5/6x,4('_____'))
      write(nout,3050) (n,(k12p(l,n),l=1,npart),n=1,npart)
 3050 format(i4,' ]',3i5)
 4050 format(i4,' ]',4i5)
      write(nout,3060) (k,k=1,nprs)
 3060 format(/' ntexch'/' nt]np',3i5/6x,3('_____'))
 4060 format(/' ntexch'/' nt]np',6i5/6x,6('_____'))
      write(nout,3070) (j,(ntexch(j,k),k=1,nprs),j=1,nt)
 3070 format(i4,' ]',3i5)
 4070 format(i4,' ]',6i5)
      write(nout,3080) (n,n=1,npart)
 3080 format(/' mtval'/' nt]nn',3i5/6x,3('_____'))
 4080 format(/' mtval'/' nt]nn',4i5/6x,4('_____'))
      write(nout,3050) (j,(mtval(j,l),l=1,npart),j=1,nt)
      write(nout,3090) (k,k=1,nprs)
 3090 format(/' nsexch'/' ns]np',3i5/6x,3('_____'))
 4090 format(/' nsexch'/' ns]np',6i5/6x,6('_____'))
      write(nout,3070) (i,(nsexch(i,k),k=1,nprs),i=1,ns)
      write(nout,3100) (n,n=1,npart)
 3100 format(/' msval'/' ns]nn',3i5/6x,3('_____'))
 4100 format(/' msval'/' ns]nn',4i5/6x,4('_____'))
      write(nout,3050) (i,(msval(i,l),l=1,npart),i=1,ns)
      write(nout,3110) (n,n=1,npart)
 3110 format(/' nsflip'/' ns]nn',3i5/6x,3('_____'))
 4110 format(/' nsflip'/' ns]nn',4i5/6x,4('_____'))
      write(nout,3050) (i,(nsflip(i,l),l=1,npart),i=1,ns)
      write(nout,3120) (j,j=1,nt)
 3120 format(/' rphi'/' ns]nt',3i5/6x,3('_____'))
 4120 format(/' rphi'/' ns]nt',6i5/6x,6('_____'))
      write(nout,3130) (i,(rphi(i,j),j=1,nt),i=1,ns)
 3130 format(i4,' ]',3f5.0)
 4130 format(i4,' ]',6f5.0)
      write(nout,140)
  140 format(/5x,'r',10x,'fc',9x,'us',9x,'ut',9x,'uct',8x,'ust'
     &      ,8x,'utt')
      write(nout,150) (rr(i),fc(i),us(i),ut(i),uct(i),ust(i),utt(i)
     &                ,i=5,npts,5)
  150 format(7f11.5)
      write(nout,100)
      return
      end
c id* setrnd rndnmb ****************************************************
c pseudo-random number generator
c setrnd is initialization call
c from seminumerical algorithms by knuth vol.2 p26 (1981)
c dseed and vmultp should be < vmultp
c nconst should be relatively prime to vmdul1
c vmultp should be odd and congruent to 1 mod 4 and not small
c rnmb(n) = set of n random numbers
c ----------------------------------------------------------------------
      subroutine setrnd(dseed)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension rnmb(1),y(55)
      save y,jang,kang,vmodul
      vmodul=2.**47
      vmdul1=2.**23
      vmultp=16805.
      nconst=0
      do 10 i=1,55
        dseed=mod(vmultp*dseed+nconst,vmdul1)
        pa1=dseed*vmdul1
        dseed=mod(vmultp*dseed+nconst,vmdul1)
        y(i)=(pa1+dseed)
   10 continue
      jang=24
      kang=55
      return
c ----------------------------------------------------------------------
      entry rndnmb(n,rnmb)
      do 100 i=1,n
        y(kang)=mod(y(kang)+y(jang),vmodul)
        rnmb(i)=y(kang)/vmodul
        jang=(jang+54)-55*int((jang+54)/56)
        kang=(kang+54)-55*int((kang+54)/56)
  100 continue
      return
      end
c id* moveem ***********************************************************
c subroutine for selecting points from monte carlo weight <psi]psi>
c trial function determined by specifying location of particles and
c order of noncentral pair operators in wave function
c r(3,npart,2) = x,y,z coordinates for each particle in last good and
c                new trial configurations: last index is specified by
c                igood = 1 (or 2) and itrial = 2 (or 1)
c iordr(nprs) = order of operators on right-hand side
c iordl(nprs) = order of operators on left-hand side
c psi2(2) = <psi]psi> for last good and current trial functions
c ----------------------------------------------------------------------
      subroutine moveem
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nrand=3*npart+2*nprs-1)
      parameter (nflops=9*npart+2*nprs-1)
      common /order/ iordl(nprs),iordr(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /wavtmp/ ctril(ns,nt),ctrir(ns,nt)
      dimension rnmb(nrand)
c ----------------------------------------------------------
c nmoves # of attempts are made to move to new configuration
c ----------------------------------------------------------
      do 100 move=1,nmoves
        call rndnmb(nrand,rnmb)
        nrn=1
c -----------------------------------------------------
c generate trial configuration
c size of attempted step controlled by steps
c iop flags if position cannot be interpolated reliably
c in which case move is rejected
c -----------------------------------------------------
        do 5 mn=1,3*npart
          r(mn,1,itrial)=r(mn,1,igood)+steps*(rnmb(nrn)-.5)
          nrn=nrn+1
    5   continue
        call interp(.true.)
        if (iop.ne.0) go to 100
c --------------------------------------
c generate trial order of pair operators
c --------------------------------------
        do 10 i=1,nprs-1
          i1=nprs+1-i
          ir=int(i1*rnmb(nrn))+1
          nrn=nrn+1
          mr=iordr(i1)
          iordr(i1)=iordr(ir)
          iordr(ir)=mr
          il=int(i1*rnmb(nrn))+1
          nrn=nrn+1
          ml=iordl(i1)
          iordl(i1)=iordl(il)
          iordl(il)=ml
   10   continue
c --------------------------------------------
c construct full trial function on rhs and lhs
c take dot product to obtain <psi]psi>
c compare to last good position and accept or
c reject new postition accordingly
c --------------------------------------------
        call wav(ctrir,iordr)
        call wav(ctril,iordl)
        psi2(itrial)=rcdot(ctril,ctrir)
        if (psi2(itrial).lt.0.) iwarn(3)=iwarn(3)+1
        if (psi2(itrial).gt.rnmb(nrn)*psi2(igood)) then
c ------------------------------------------
c if move is accepted save new wave function
c and interchange igood, itrial labels
c ------------------------------------------
          do 40 ij=1,ns*nt
            cpsl(ij,1)=ctril(ij,1)
            cpsr(ij,1)=ctrir(ij,1)
   40     continue
          mok=mok+1
          igood=itrial
          itrial=3-igood
          call saveit
        endif
c ---------------------------------------------
c now, whether itrial was accepted or rejected,
c the variable igood refers to the correct
c set of coordinates
c ---------------------------------------------
        mflops=mflops+nflops
  100 continue
      call restor
      return
      end
c id* remove ***********************************************************
c subroutine for constructing wave function when old random walk
c is being followed (called instead of moveem)
c ----------------------------------------------------------------------
      subroutine remove
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      common /order/ iordl(nprs),iordr(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
c ---------------------------------------------------------
c save reference <psi(old)]psi(old)> under itrial label
c set trial positions to reference values for interpolation
c compute new wave function and store under igood label
c ---------------------------------------------------------
      itrial=3-igood
      psi2(itrial)=psi2(igood)
      do 10 mn=1,3*npart
        r(mn,1,itrial)=r(mn,1,igood)
   10 continue
      call interp(.false.)
      call wav(cpsr,iordr)
      call wav(cpsl,iordl)
      psi2(igood)=rcdot(cpsl,cpsr)
      call saveit
      return
      end
c id* interp ***********************************************************
c subrtoutine for interpolating correlations
c fcv, uxv interpolated from fc, ux (x=s,t,ct,st,tt)
c drh(m,k) are components of unit pair vectors
c drm(k) are pair distances
c ----------------------------------------------------------------------
      subroutine interp(test)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nflops=69*nprs)
      logical test
      common /corr/ fc(npts),us(npts),ut(npts),uct(npts),ust(npts)
     &,utt(npts),rr(npts)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      iop=0
c --------------------------------------
c calculate pair distances, unit vectors
c --------------------------------------
      do 10 k=1,nprs
        n1=n1p(k)
        n2=n2p(k)
        drh(1,k)=r(1,n1,itrial)-r(1,n2,itrial)
        drh(2,k)=r(2,n1,itrial)-r(2,n2,itrial)
        drh(3,k)=r(3,n1,itrial)-r(3,n2,itrial)
        drm(k)=sqrt(drh(1,k)**2+drh(2,k)**2+drh(3,k)**2)
        drh(1,k)=drh(1,k)/drm(k)
        drh(2,k)=drh(2,k)/drm(k)
        drh(3,k)=drh(3,k)/drm(k)
   10 continue
c -------------------------
c find interpolation points
c -------------------------
      do 20 k=1,nprs
        frac=drm(k)/h
        il=nint(frac)
c ------------------------------------------------------------
c when called from moveem (test=true) following test is made:
c if pair distance is too small or too large so that
c interpolation is questionable, increment warning counter
c and return for new positions, rejecting this move
c when called from remove or ekin (test=false) no test is made
c ------------------------------------------------------------
        if (.not.test) go to 15
        if (il.lt.3) then
          iop=1
          iwarn(1)=iwarn(1)+1
          go to 100
        else if (il.gt.(npts-3)) then
          iop=1
          iwarn(2)=iwarn(2)+1
          go to 100
        end if
   15   ils(k)=il
        frac=frac-real(il)
        fras=frac**2
        fr0(k)=1.-fras
        frpl(k)=.5*(fras+frac)
        frmn(k)=.5*(fras-frac)
   20 continue
c ---------------------
c interpolate functions
c ---------------------
      do 30 k=1,nprs
        il=ils(k)
        ipl=il+1
        imn=il-1
        x0=fr0(k)
        xpl=frpl(k)
        xmn=frmn(k)
        fcv(k)=x0*fc(il)+xpl*fc(ipl)+xmn*fc(imn)
        usv(k)=x0*us(il)+xpl*us(ipl)+xmn*us(imn)
        utv(k)=x0*ut(il)+xpl*ut(ipl)+xmn*ut(imn)
        uctv(k)=x0*uct(il)+xpl*uct(ipl)+xmn*uct(imn)
        ustv(k)=x0*ust(il)+xpl*ust(ipl)+xmn*ust(imn)
        uttv(k)=x0*utt(il)+xpl*utt(ipl)+xmn*utt(imn)
   30 continue
      mflops=mflops+nflops
  100 return
      end
c id* wav **************************************************************
c subroutine for generating wave function
c starting from uncorrelated spin-ispin wave function rphi,
c multiplies in product of central pair correlations, then
c operates with successive sums of noncentral pair operators
c (including three-body correlation) in order specified by iord
c ----------------------------------------------------------------------
      subroutine wav(ctri,iord)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (ntrps=(nprs*(npart-2))/3)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nflops=ns*nt+6*nprs+90*ntrps)
      dimension ctri(ns,nt),iord(nprs),cdum(ns,nt)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /tbodyc/ t1s,t1t,nt2s,nt2t,t3s,t3t
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
c ------------------------------------------
c construct central correlated wave function
c ------------------------------------------
      ff=fcv(1)
      do 10 k=2,nprs
        ff=ff*fcv(k)
   10 continue
      do 20 ij=1,ns*nt
        cdum(ij,1)=ff*rphi(ij,1)
   20 continue
c ----------------------------------
c perform noncentral pair operations
c cdum & ctri are alternately used
c for input and output
c ----------------------------------
      iop=2
      do 40 k=1,nprs
        iop=3-iop
        k12=iord(k)
c ----------------------
c three-body correlation
c ----------------------
        f3bs=1.
        f3bt=1.
        do 30 l=1,npart
          k13=abs( k12p(n1p(k12),l) )
          k23=abs( k12p(n2p(k12),l) )
          if (k13.eq.0 .or. k23.eq.0) go to 30
          rsum=drm(k12)+drm(k13)+drm(k23)
          f3bt=f3bt*(1.-t1t*(drm(k12)/rsum)**nt2t*exp(-t3t*rsum))
          f3bs=f3bs*(1.-t1s*(drm(k12)/rsum)**nt2s*exp(-t3s*rsum))
   30   continue
c --------------------------------
c noncentral correlation functions
c --------------------------------
        wc=1.
        ws=f3bs*usv(k12)
        wt=f3bt*utv(k12)
        wct=f3bs*uctv(k12)
        wst=f3bs*ustv(k12)
        wtt=f3bt*uttv(k12)
c ------------------------
c call operator subroutine
c ------------------------
        if (iop.eq.1) then
          call f6op(k12,ctri,cdum,drh,wc,wct,ws,wst,wt,wtt)
        else
          call f6op(k12,cdum,ctri,drh,wc,wct,ws,wst,wt,wtt)
        end if
   40 continue
c ---------------------------
c ensure final output is ctri
c ---------------------------
      if (iop.eq.2) then
        do 50 ij=1,ns*nt
          ctri(ij,1)=cdum(ij,1)
   50   continue
      end if
      mflops=mflops+nflops
      return
      end
c id* f6op *************************************************************
c subroutine for pair operator
c
c this subroutine explicitly written in real arithmetic:
c cxxr = real component of cxx, cxxi = imaginary component of cxx
c
c kp = pair number
c cl = output wave function
c cr = input wave function
c drh = unit vectors for pair distances
c wx = correlation function values
c  x = c (central) ct (ispin) s (spin) st (spin-ispin)
c      t (tensor)  tt (tensor-ispin)
c ----------------------------------------------------------------------
      subroutine f6op(kp,cl,cr,drh,wc,wct,ws,wst,wt,wtt)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=10+(11+52*nt)*ns)
      dimension cr(2,ns,nt),cl(2,ns,nt),drh(3,nprs)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      dimension csmr(ns,nt),csmi(ns,nt),ctmr(ns,nt),ctmi(ns,nt)
     &,ctr(2),cti(2),et(2)
      n1=n1p(kp)
      n2=n2p(kp)
      et(1)=drh(3,kp)
      et(2)=-et(1)
      ctr(1)=drh(1,kp)
      ctr(2)=ctr(1)
      cti(1)=drh(2,kp)
      cti(2)=-cti(1)
      wcd=wc-wct
      wce=2.*wct
      wsd=ws-wst-wt+wtt
      wse=2.*(wst-wtt)
      wtd=3.*(wt-wtt)
      wte=6.*wtt
c --------------
c ispin exchange
c --------------
      do 20 j=1,nt
        je=ntexch(j,kp)
        do 10 i=1,ns
          cl(1,i,j)=wcd*cr(1,i,j)+wce*cr(1,i,je)
          cl(2,i,j)=wcd*cr(2,i,j)+wce*cr(2,i,je)
          csmr(i,j)=wsd*cr(1,i,j)+wse*cr(1,i,je)
          csmi(i,j)=wsd*cr(2,i,j)+wse*cr(2,i,je)
          ctmr(i,j)=wtd*cr(1,i,j)+wte*cr(1,i,je)
          ctmi(i,j)=wtd*cr(2,i,j)+wte*cr(2,i,je)
   10   continue
   20 continue
c ---------------------------
c spin exchange and spin flip
c ---------------------------
      do 40 i=1,ns
        ie=nsexch(i,kp)
        if1=nsflip(i,n1)
        if2=nsflip(i,n2)
        ifb=if1+if2-i
        m1=(3-msval(i,n1))/2
        m2=(3-msval(i,n2))/2
        fn=et(m1)*et(m2)
        cf1r=ctr(m1)*et(m2)
        cf1i=cti(m1)*et(m2)
        cf2r=et(m1)*ctr(m2)
        cf2i=et(m1)*cti(m2)
        cfbr=ctr(m1)*ctr(m2)-cti(m1)*cti(m2)
        cfbi=ctr(m1)*cti(m2)+cti(m1)*ctr(m2)
        do 30 j=1,nt
          cl(1,i,j)=cl(1,i,j)+fn*ctmr(i,j)-csmr(i,j)+2.*csmr(ie,j)
          cl(2,i,j)=cl(2,i,j)+fn*ctmi(i,j)-csmi(i,j)+2.*csmi(ie,j)
          cl(1,if1,j)=cl(1,if1,j)+cf1r*ctmr(i,j)-cf1i*ctmi(i,j)
          cl(2,if1,j)=cl(2,if1,j)+cf1r*ctmi(i,j)+cf1i*ctmr(i,j)
          cl(1,if2,j)=cl(1,if2,j)+cf2r*ctmr(i,j)-cf2i*ctmi(i,j)
          cl(2,if2,j)=cl(2,if2,j)+cf2r*ctmi(i,j)+cf2i*ctmr(i,j)
          cl(1,ifb,j)=cl(1,ifb,j)+cfbr*ctmr(i,j)-cfbi*ctmi(i,j)
          cl(2,ifb,j)=cl(2,ifb,j)+cfbr*ctmi(i,j)+cfbi*ctmr(i,j)
   30   continue
   40 continue
      mflops=mflops+nflops
      return
      end
c id* rcdot ************************************************************
c function for dot product <cl]cr>
c note: complex conjugation of cl done here
c this function explicitly real arithemtic
c ----------------------------------------------------------------------
      function rcdot(cl,cr)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=4*ns*nt)
      common /speed/ mflops
      dimension cl(2,ns,nt),cr(2,ns,nt)
      x1=0.
      do 10 ij=1,ns*nt
        x1=x1+cl(1,ij,1)*cr(1,ij,1)+cl(2,ij,1)*cr(2,ij,1)
   10 continue
      mflops=mflops+nflops
      rcdot=x1
      return
      end
c id* saveit restor ****************************************************
c subroutine keeps track of interpolation info for last good move
c ----------------------------------------------------------------------
      subroutine saveit
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      common /order/ iordl(nprs),iordr(nprs)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      dimension drhs(3,nprs),drms(nprs),iordls(nprs),iordrs(nprs)
     &,fr0s(nprs),frpls(nprs),frmns(nprs),ilss(nprs)
      save drhs,drms,iordls,iordrs,fr0s,frpls,frmns,ilss
      do 10 j=1,nprs
        iordrs(j)=iordr(j)
        iordls(j)=iordl(j)
        drhs(1,j)=drh(1,j)
        drhs(2,j)=drh(2,j)
        drhs(3,j)=drh(3,j)
        drms(j)=drm(j)
        fr0s(j)=fr0(j)
        frpls(j)=frpl(j)
        frmns(j)=frmn(j)
   10 ilss(j)=ils(j)
      return
c ----------------------------------------------------------------------
      entry restor
      do 20 j=1,nprs
        iordr(j)=iordrs(j)
        iordl(j)=iordls(j)
        drh(1,j)=drhs(1,j)
        drh(2,j)=drhs(2,j)
        drh(3,j)=drhs(3,j)
        drm(j)=drms(j)
        fr0(j)=fr0s(j)
        frpl(j)=frpls(j)
        frmn(j)=frmns(j)
   20 ils(j)=ilss(j)
      return
      end
c id* rhor *************************************************************
c subroutine for binning densities
c rhoa = central density
c rhot = ispin density (to differentiate between p & n densities)
c ----------------------------------------------------------------------
      subroutine rhor
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=3+(35+2*ns*nt)*npart)
      parameter (iflnwk=16*ns*nt+1,ilftwk=iflnwk-2*(ns*nt))
      common /dens/ rhoa(100),rhot(100)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /fwksp/ ctzn(ns,nt),dum(ilftwk)
c -----------------------
c find center of mass rcm
c -----------------------
      do 20 m=1,3
        rcm(m)=0.
        do 10 n=1,npart
          rcm(m)=rcm(m)+r(m,n,igood)
   10   continue
        rcm(m)=rcm(m)/npart
   20 continue
c ---------------------
c bin central densities
c ---------------------
      do 40 n=1,npart
        rpl=sqrt((r(1,n,igood)-rcm(1))**2
     &   +(r(2,n,igood)-rcm(2))**2+(r(3,n,igood)-rcm(3))**2)
        kx=int(10*rpl)+1
        if (kx.gt.100) kx=100
        rhoa(kx)=rhoa(kx)+1.
c --------------------------------
c compute <psi]tau(z)]psi> and bin
c --------------------------------
        do 30 j=1,nt
        do 30 i=1,ns
          ctzn(i,j)=mtval(j,n)*cpsr(i,j)
   30   continue
        rhot(kx)=rhot(kx)+rcdot(cpsl,ctzn)/psi2(igood)
   40 continue
      mflops=mflops+nflops
      return
      end
c id* ekin *************************************************************
c subroutine for calculating kinetic energy ti
c cdelr(i,j,m,n) = m component of grad]psi> wrt particle n
c                 (generated for vso calculation in subroutine vpot)
c cdlsq(i,j) = sum over all particles laplacian]psi>
c ----------------------------------------------------------------------
      subroutine ekin(ti)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nflops=3+(12+42*ns*nt)*npart)
      parameter (iflnwk=16*ns*nt+1,ilftwk=iflnwk-2*(2*ns*nt))
      common /const/ h2m,pi
      common /nabla/ cdelr(ns,nt,3,npart),cdlsq(ns,nt)
      common /order/ iordl(nprs),iordr(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /fwksp/ cplus(ns,nt),cmins(ns,nt),dum(ilftwk)
c -------------------------------------------------
c set trial configuration for use in interpolation
c -------------------------------------------------
      do 10 ij=1,ns*nt
        cdlsq(ij,1)=(0.,0.)
   10 continue
      do 20 mn=1,3*npart
        r(mn,1,itrial)=r(mn,1,igood)
   20 continue
c -------------------------------------------------
c move each particle in 3 dimensions, one at a time
c interpoalte and generate ]psi>+ and ]psi>-
c -------------------------------------------------
      do 50 n=1,npart
        do 40 m=1,3
          r(m,n,itrial)=r(m,n,igood)+h
          call interp(.false.)
          call wav(cplus,iordr)
          r(m,n,itrial)=r(m,n,igood)-h
          call interp(.false.)
          call wav(cmins,iordr)
          r(m,n,itrial)=r(m,n,igood)
c ---------------------------------------------------
c compute grad]psi> and laplacian]psi> by differences
c ---------------------------------------------------
          do 30 ij=1,ns*nt
            cdelr(ij,1,m,n)=(cplus(ij,1)-cmins(ij,1))/(2.*h)
            cdlsq(ij,1)=cdlsq(ij,1)+(cplus(ij,1)+cmins(ij,1)
     &       -2.*cpsr(ij,1))/h**2
   30     continue
   40   continue
   50 continue
c ------------------
c sum kinetic energy
c ------------------
      ti=-.5*h2m*rcdot(cpsl,cdlsq)/psi2(igood)
c ------------------------------------------------
c restore interp variables from last good position
c ------------------------------------------------
      call restor
      mflops=mflops+nflops
      return
      end
c id* vpot ************************************************************
c subroutine for computing potential expectation value
c v6 = six-operator contribution (c,s,t,ct,st,tt)
c vso = spin-orbit contribution (ls,lst)
c vcf = coulomb contribution (with form factor)
c ----------------------------------------------------------------------
      subroutine vpot(v6,vso,vcf)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=(34+5*ns+4*nt+64*ns*nt)*nprs)
      parameter (iflnwk=16*ns*nt+1,ilftwk=iflnwk-2*(7*ns*nt))
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /nabla/ cdelr(ns,nt,3,npart),cdlsq(ns,nt)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),spiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      common /fwksp/ c6psr(ns,nt),cbpsr(ns,nt),ccpsr(ns,nt)
     &,cls(ns,nt),cl(ns,nt,3),dum(ilftwk)
      ci=cmplx(0.,1.)
      v6=0.
      vso=0.
      vcf=0.
      call intrpv
c ----------------------------------------------
c six-operator contribution uses f6op subroutine
c ----------------------------------------------
      do 100 k=1,nprs
        call f6op(k,c6psr,cpsr,drh,vcv(k),vctv(k),vsv(k),vstv(k)
     &           ,vtv(k),vttv(k))
        v6=v6+rcdot(cpsl,c6psr)/psi2(igood)
c ---------------------------------------------------
c spin-orbit contribution
c uses grad]psi> (cdelr) generated by subroutine ekin
c to construct l]psi> (cl)
c then constructs l.s]psi> (cls) by
c using spin raising and lowering operations
c ---------------------------------------------------
        do 10 ij=1,ns*nt
          cls(ij,1)=(0.,0.)
   10   continue
        n1=n1p(k)
        n2=n2p(k)
        rx=drm(k)*drh(1,k)
        ry=drm(k)*drh(2,k)
        rz=drm(k)*drh(3,k)
        do 20 ij=1,ns*nt
          cl(ij,1,1)=-.5*ci*(ry*(cdelr(ij,1,3,n1)-cdelr(ij,1,3,n2))
     &                      -rz*(cdelr(ij,1,2,n1)-cdelr(ij,1,2,n2)))
          cl(ij,1,2)=-.5*ci*(rz*(cdelr(ij,1,1,n1)-cdelr(ij,1,1,n2))
     &                      -rx*(cdelr(ij,1,3,n1)-cdelr(ij,1,3,n2)))
          cl(ij,1,3)=-.5*ci*(rx*(cdelr(ij,1,2,n1)-cdelr(ij,1,2,n2))
     &                      -ry*(cdelr(ij,1,1,n1)-cdelr(ij,1,1,n2)))
   20   continue
        do 40 i=1,ns
          if1=nsflip(i,n1)
          if2=nsflip(i,n2)
          spl1=.25*(1+msval(i,n1))
          spl2=.25*(1+msval(i,n2))
          smn1=.25*(1-msval(i,n1))
          smn2=.25*(1-msval(i,n2))
          sz12=.5*(msval(i,n1)+msval(i,n2))
          do 30 j=1,nt
            cpl=cl(i,j,1)+ci*cl(i,j,2)
            cmn=cl(i,j,1)-ci*cl(i,j,2)
            cls(i,j)=cls(i,j)+sz12*cl(i,j,3)
            cls(if1,j)=cls(if1,j)+spl1*cpl+smn1*cmn
            cls(if2,j)=cls(if2,j)+spl2*cpl+smn2*cmn
   30     continue
   40   continue
c ----------------------------------------------------
c potential contributions now evaluated using l.s]psi>
c ----------------------------------------------------
        vlsd=vlsv(k)-vlstv(k)
        vlse=2.*vlstv(k)
        do 60 j=1,nt
          je=ntexch(j,k)
          do 50 i=1,ns
            cbpsr(i,j)=vlsd*cls(i,j)+vlse*cls(i,je)
   50     continue
   60   continue
        vso=vso+rcdot(cpsl,cbpsr)/psi2(igood)
c -----------------------------------------------------
c coulomb contribution <psi] ff*e**2/r (1+tau(z)) ]psi>
c where ff is form factor
c -----------------------------------------------------
        x=4.27*drm(k)
        ff=1.-exp(-x)*(1.+x*(33.+x*(9.+x))/48.)
        ffr=ff/drm(k)
        do 80 j=1,nt
          q1=.5*(1+mtval(j,n1))
          q2=.5*(1+mtval(j,n2))
          qsqr=q1*q2*ffr
          do 70 i=1,ns
            ccpsr(i,j)=qsqr*cpsr(i,j)
   70     continue
   80   continue
        vcf=vcf+1.44*rcdot(cpsl,ccpsr)/psi2(igood)
  100 continue
      mflops=mflops+nflops
      return
      end
c id* intrpv ***********************************************************
c subroutine for interpolating potential functions
c interpolation points saved from last good call to interp
c ----------------------------------------------------------------------
      subroutine intrpv
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (nflops=50*nprs)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /poten/ vc(npts),vs(npts),vt(npts),vls(npts)
     &,vct(npts),vst(npts),vtt(npts),vlst(npts),vtpi(npts),vypi(npts)
      common /speed/ mflops
      do 10 k=1,nprs
        il=ils(k)
        ipl=il+1
        imn=il-1
        x0=fr0(k)
        xpl=frpl(k)
        xmn=frmn(k)
        vcv(k)=vc(il)*x0+vc(ipl)*xpl+vc(imn)*xmn
        vctv(k)=vct(il)*x0+vct(ipl)*xpl+vct(imn)*xmn
        vsv(k)=vs(il)*x0+vs(ipl)*xpl+vs(imn)*xmn
        vstv(k)=vst(il)*x0+vst(ipl)*xpl+vst(imn)*xmn
        vtv(k)=vt(il)*x0+vt(ipl)*xpl+vt(imn)*xmn
        vttv(k)=vtt(il)*x0+vtt(ipl)*xpl+vtt(imn)*xmn
        vlsv(k)=vls(il)*x0+vls(ipl)*xpl+vls(imn)*xmn
        vlstv(k)=vlst(il)*x0+vlst(ipl)*xpl+vlst(imn)*xmn
        vypiv(k)=vypi(il)*x0+vypi(ipl)*xpl+vypi(imn)*xmn
        vtpiv(k)=vtpi(il)*x0+vtpi(ipl)*xpl+vtpi(imn)*xmn
   10 continue
      mflops=mflops+nflops
      return
      end
c id* tbpot ************************************************************
c subroutine for three-nucleon potential (urbana model vii)
c v3a = anticommutator part of two-pion exchange term
c v3c = commutator part of two-pion exchange term
c v3u = repulsive short-range part of urbana potential
c strengths set by parameters u2pi & u3nr
c ----------------------------------------------------------------------
      subroutine tbpot(v3a,v3c,v3u)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (ntrps=(nprs*(npart-2))/3)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (npts=750,h=.02)
      parameter (u2pi=-0.0333,u3nr=0.0038)
      parameter (nflops=33*ntrps)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /pairs/ drh(3,nprs),drm(nprs)
     &,fcv(nprs),usv(nprs),utv(nprs),uctv(nprs),ustv(nprs),uttv(nprs)
     &,vcv(nprs),vsv(nprs),vtv(nprs),vlsv(nprs),vctv(nprs),vstv(nprs)
     &,vttv(nprs),vlstv(nprs),vtpiv(nprs),vypiv(nprs)
      common /pass/ steps,psi2(2),r(3,npart,2),rcm(3)
     &,fr0(nprs),frpl(nprs),frmn(nprs),ils(nprs)
     &,igood,itrial,nmoves,mok,iwarn(3),iop
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /wavf/ rphi(ns,nt),cpsl(ns,nt),cpsr(ns,nt)
      dimension capsr(ns,nt),ccpsr(ns,nt),csav(ns,nt,nprs)
      v3a=0.
      v3c=0.
      v3u=0.
c ------------------------------------------------------
c calculate spin.spin and tensor operators for all pairs
c acting on ]psi> and store in csav
c ------------------------------------------------------
      call stone(csav,cpsr,drh,vypiv,vtpiv)
c --------------------
c sum over all triples
c --------------------
      do 20 n=1,npart
        do 10 k=1,nprs
          kp1=abs( k12p(n1p(k),n) )
          kp2=abs( k12p(n2p(k),n) )
          if (kp1.eq.0 .or. kp2.eq.0) go to 10
c ----------------------------------------------
c generate anticommutator & commutator operators
c capsr = {xp1,xp2}]psi>, ccpsr = Õxp1,xp2å]psi>
c sum potential contributions
c ----------------------------------------------
          call stac(kp1,kp2,capsr,ccpsr,csav,drh,vypiv,vtpiv)
          v3a=v3a+u2pi*rcdot(cpsl,capsr)/psi2(igood)
          v3c=v3c+.25*u2pi*rcdot(cpsl,ccpsr)/psi2(igood)
          v3u=v3u+u3nr*(vtpiv(kp1)*vtpiv(kp2))**2
   10   continue
   20 continue
      mflops=mflops+nflops
      return
      end
c id* stone ************************************************************
c subroutine for generating xp]psi> for all pairs
c where xp = ( ws(r)*spin.spin + wt(r)*tensor ) is pair operator
c results saved in array cs
c logic similar to subroutine f6op, but complex arithmetic used
c ----------------------------------------------------------------------
      subroutine stone(cs,cr,drh,ws,wt)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=(4+(18+36*nt)*ns)*nprs)
      dimension cs(ns,nt,nprs),cr(ns,nt),drh(3,nprs),ws(nprs),wt(nprs)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      dimension et(2),ctf(2)
      do 40 k=1,nprs
        n1=n1p(k)
        n2=n2p(k)
        et(1)=drh(3,k)
        et(2)=-drh(3,k)
        ctf(1)=cmplx(drh(1,k),drh(2,k))
        ctf(2)=cmplx(drh(1,k),-drh(2,k))
        wsd=-ws(k)+wt(k)
        wse=2*(ws(k)-wt(k))
        wtd=3*wt(k)
        do 10 ij=1,ns*nt
          cs(ij,1,k)=(0.,0.)
   10   continue
c --------------------------------------
c spin exchange and spin flip operations
c --------------------------------------
        do 30 i=1,ns
          ie=nsexch(i,k)
          if1=nsflip(i,n1)
          if2=nsflip(i,n2)
          ifb=if1+if2-i
          m1=(3-msval(i,n1))/2
          m2=(3-msval(i,n2))/2
          fn=et(m1)*et(m2)
          cf1=ctf(m1)*et(m2)
          cf2=et(m1)*ctf(m2)
          cfb=ctf(m1)*ctf(m2)
          do 20 j=1,nt
            cs(i,j,k)=cs(i,j,k)+wsd*cr(i,j)+wse*cr(ie,j)+wtd*fn*cr(i,j)
            cs(if1,j,k)=cs(if1,j,k)+wtd*cf1*cr(i,j)
            cs(if2,j,k)=cs(if2,j,k)+wtd*cf2*cr(i,j)
            cs(ifb,j,k)=cs(ifb,j,k)+wtd*cfb*cr(i,j)
   20     continue
   30   continue
   40 continue
      mflops=mflops+nflops
      return
      end
c id* stac *************************************************************
c subroutine to generate anticommutator & commutator operators
c ca = {tp1,tp2}{xp1,xp2}]psi>
c cc = Õtp1,tp2åÕxp1,xp2å]psi>
c where tpk = tau.tau for pair k
c expects stored values of xp]psi> for all pairs in cs
c logic similar to subroutines f6op and stone
c ----------------------------------------------------------------------
      subroutine stac(k1,k2,ca,cc,cs,drh,ws,wt)
      implicit real*8 (a-b,d-h,o-z)
      implicit integer*4 (i-n)
      implicit complex*16 (c)
      parameter (npart=3,ns=2**npart,nprs=(npart*(npart-1))/2)
      parameter (nprot=1,nt=3,nfact=6)
      parameter (nflops=8+(36+112*nt)*ns)
      parameter (iflnwk=16*ns*nt+1,ilftwk=iflnwk-2*(8*ns*nt))
      dimension ca(ns,nt),cc(ns,nt),cs(ns,nt,nprs)
     &,drh(3,nprs),ws(nprs),wt(nprs)
      common /ispin/ ntexch(nt,nprs),mtval(nt,npart)
      common /lbp/ n1p(nprs),n2p(nprs),k12p(npart,npart)
      common /speed/ mflops
      common /spin/ nsexch(ns,nprs),msval(ns,npart),nsflip(ns,npart)
      common /fwksp/ c12(ns,nt),c21(ns,nt),cpl(ns,nt),cmn(ns,nt)
     &,c1p(ns,nt),c1m(ns,nt),c2p(ns,nt),c2m(ns,nt),dum(ilftwk)
      dimension etk1(2),etk2(2),ctfk1(2),ctfk2(2)
c ---------------------------------
c xxk1 denotes xx value for pair k1
c xxk2 denotes xx value for pair k2
c ---------------------------------
      n1k1=n1p(k1)
      n2k1=n2p(k1)
      etk1(1)=drh(3,k1)
      etk1(2)=-drh(3,k1)
      ctfk1(1)=cmplx(drh(1,k1),drh(2,k1))
      ctfk1(2)=cmplx(drh(1,k1),-drh(2,k1))
      wsdk1=-ws(k1)+wt(k1)
      wsek1=2.*(ws(k1)-wt(k1))
      wtdk1=3.*wt(k1)
      n1k2=n1p(k2)
      n2k2=n2p(k2)
      etk2(1)=drh(3,k2)
      etk2(2)=-drh(3,k2)
      ctfk2(1)=cmplx(drh(1,k2),drh(2,k2))
      ctfk2(2)=cmplx(drh(1,k2),-drh(2,k2))
      wsdk2=-ws(k2)+wt(k2)
      wsek2=2.*(ws(k2)-wt(k2))
      wtdk2=3.*wt(k2)
      do 10 ij=1,ns*nt
        c12(ij,1)=(0.,0.)
        c21(ij,1)=(0.,0.)
   10 continue
c --------------------------------------
c spin exchange and spin flip operations
c --------------------------------------
      do 30 i=1,ns
        iek1=nsexch(i,k1)
        if1k1=nsflip(i,n1k1)
        if2k1=nsflip(i,n2k1)
        ifbk1=if1k1+if2k1-i
        m1k1=(3-msval(i,n1k1))/2
        m2k1=(3-msval(i,n2k1))/2
        fnk1=etk1(m1k1)*etk1(m2k1)
        cf1k1=ctfk1(m1k1)*etk1(m2k1)
        cf2k1=etk1(m1k1)*ctfk1(m2k1)
        cfbk1=ctfk1(m1k1)*ctfk1(m2k1)
        iek2=nsexch(i,k2)
        if1k2=nsflip(i,n1k2)
        if2k2=nsflip(i,n2k2)
        ifbk2=if1k2+if2k2-i
        m1k2=(3-msval(i,n1k2))/2
        m2k2=(3-msval(i,n2k2))/2
        fnk2=etk2(m1k2)*etk2(m2k2)
        cf1k2=ctfk2(m1k2)*etk2(m2k2)
        cf2k2=etk2(m1k2)*ctfk2(m2k2)
        cfbk2=ctfk2(m1k2)*ctfk2(m2k2)
c ------------------
c c12 = xp1*xp2]psi>
c c21 = xp2*xp1]psi>
c ------------------
        do 20 j=1,nt
          c12(i,j)=c12(i,j)+wsdk1*cs(i,j,k2)+wsek1*cs(iek1,j,k2)
     &                     +wtdk1*fnk1*cs(i,j,k2)
          c12(if1k1,j)=c12(if1k1,j)+wtdk1*cf1k1*cs(i,j,k2)
          c12(if2k1,j)=c12(if2k1,j)+wtdk1*cf2k1*cs(i,j,k2)
          c12(ifbk1,j)=c12(ifbk1,j)+wtdk1*cfbk1*cs(i,j,k2)
          c21(i,j)=c21(i,j)+wsdk2*cs(i,j,k1)+wsek2*cs(iek2,j,k1)
     &                     +wtdk2*fnk2*cs(i,j,k1)
          c21(if1k2,j)=c21(if1k2,j)+wtdk2*cf1k2*cs(i,j,k1)
          c21(if2k2,j)=c21(if2k2,j)+wtdk2*cf2k2*cs(i,j,k1)
          c21(ifbk2,j)=c21(ifbk2,j)+wtdk2*cfbk2*cs(i,j,k1)
   20   continue
   30 continue
      do 40 ij=1,ns*nt
        cpl(ij,1)=c12(ij,1)+c21(ij,1)
        cmn(ij,1)=c12(ij,1)-c21(ij,1)
   40 continue
c -----------------------------------
c ispin commutators & anticommutators
c -----------------------------------
      do 60 j=1,nt
        jek1=ntexch(j,k1)
        jek2=ntexch(j,k2)
        do 50 i=1,ns
          c1p(i,j)=-cpl(i,j)+2.*cpl(i,jek1)
          c2p(i,j)=-cpl(i,j)+2.*cpl(i,jek2)
          c1m(i,j)=-cmn(i,j)+2.*cmn(i,jek1)
          c2m(i,j)=-cmn(i,j)+2.*cmn(i,jek2)
   50   continue
   60 continue
      do 80 j=1,nt
        jek1=ntexch(j,k1)
        jek2=ntexch(j,k2)
        do 70 i=1,ns
          ca(i,j)=-c2p(i,j)+2.*c2p(i,jek1)-c1p(i,j)+2.*c1p(i,jek2)
          cc(i,j)=-c2m(i,j)+2.*c2m(i,jek1)+c1m(i,j)-2.*c1m(i,jek2)
   70   continue
   80 continue
      mflops=mflops+nflops
      return
      end
