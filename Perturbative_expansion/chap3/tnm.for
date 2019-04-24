      program tnm
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      logical  exfile,first,hbo,test
      character*6 chanid,headfile,inid*3,outid*3,testid*3,memo*64
      dimension chanid(nmx),ptab(ktot),etab(ktot),propag(lpmx,mmx,mmx)
     &,         grecor(ip,ip)
      common /grpoly/gref(2*nqmx,mmx,mmx)
      common /info/memo,efm,gmasq,pfm,pav,pkmx,emev,mdgr,kup,inid,outid
      common /kmesh/mxk,pk(ktot)
      common /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)
      data    first,lpup,njup,scale/.true.,lpmx,2,4.147d1/
10    data    headfile,ngrplo,ngrpup,hbo,testid,test/
     &        'HEAD',1,4,.true.,'INV',.true./

47    format(' tnm 12 or 15: Error in INQUIRE about file ',a6,' .')
48    format(' tnm 13: Header file ',a6,' missing.')
49    format(' tnm 14: Problems with header file ',a6,' .')
50    format(' tnm 12: Input file ',a6,' missing.')

      call greetings
*
* read input from header file & and scale spectrum
*
      inquire(file=headfile,iostat=io,err=12,exist=exfile)
13    if(.not.exfile)then
         write(*,47) headfile
         stop
      endif
      open(unit=11,iostat=io,err=14,file=headfile,status='unknown')
      read(11,*,iostat=io,err=14) memo
      read(11,*,iostat=io,err=14) efm,gmasq,pfm,pav,pkmx,emev,mup,kup
     &,                           inid,outid
      read(11,*,iostat=io,err=14) (ptab(k),etab(k),k=1,kup)
      close(11,iostat=io,err=14)
      mdgr=mup
      do 100 k=1,kup
         etab(k)=etab(k)/scale
100   continue
      estart=emev/scale
*
* label and group TR components and generate in/out filenames & display
*
      call connect(inid,chanid,njup)
      do 200 i=abs(link(ngrplo,ndim(ngrplo),ndim(ngrplo)))
     &,        abs(link(ngrpup,ndim(ngrpup),ndim(ngrpup)))
         inquire(file=chanid(i),iostat=io,err=15,exist=exfile)
15       if(iostat.ne.0) write(*,47) chanid(i)
         if(.not.exfile) write(*,50) chanid(i)
200   continue
      if(hbo)call show(ngrplo,ngrpup,chanid)
*
* difference in propagators
*
      call polyprop(mup,lpup,pkmx,pav,pfm,kup,ptab,etab,estart,efm
     &,             gmasq,propag)
*
* solve BETHE-GOLDSTONE equation for each spin/parity group separately
*
      do 300 ngroup=ngrplo,ngrpup
         call grtrafo(test,testid,ngroup,chanid,mup,pkmx,efm,gmasq)
         call lgnuc(ngroup,mup,propag,grecor)
         call lgback(ngroup,chanid,mup,grecor)
300   continue

12    if(io.gt.0) write(*,48) headfile
14    if(io.gt.0) write(*,49) headfile
      end


************************************************************from MAIN***
      subroutine polyprop(mup,lpup,pkmx,pav,pfm,kup,ptab,etab,estart
     &,                   efm,gmasq,propag)
*
*  Obtains difference between reference and nuclear matter propagators
*  in polynomial representation. The k-space integration for the Legind-
*  gaard transformation is done by a 'ngs'-point GAUSS-LEGENDRE inte-
*  gration. The k-range is divided into subintervals such that limits of
*  integration coincide with the discontinuities of the Pauli operator.
*  Subroutine 'tmpgrid' makes a k-space mesh with the appropriate abs-
*  cissa values for each subinterval. Subrt. 'angprop' does the angle
*  averaging on the meshpoints. After the polynomials are obtained, the
*  subrt. 'legprop' performs the transforming integration. The propaga-
*  tor is then accumulated over the various subintervals.The polynomials
*  'pledin()' are evaluated by subrt. 'polygen'.
*  OUTPUT    propag() - Propagator difference for each Pauli coupling
*              in polynomial representation
*  Variables for the GAUSS integration are
*     ngs   - half the order of the method (here: ngs=5)
*     gweigh(), gabsc() - relative weights and abscissa values
*     grang - range of subintervals
*  COMMON 'kmesh' is used here as dummy block for the subintervals.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(zero=1.0d-14,ngs=5)
      dimension ptab(ktot),etab(ktot),propag(lpmx,mmx,mmx),prorad(ktot)
      common / kmesh/ mxk,pk(ktot)
      common / gauss/ gweigh(ngs),gabsc(ngs)
      lp=1
      do 200 m=1,mup
         do 300 n=1,mup
            propag(lp,m,n)=0.0d0
300      continue
200   continue
*
* find integration regions for Q=1, 1>Q>0, Q=0; 0='zero'
*
      if(pav.ge.pfm)then
         plim1=zero
      else
         plim1=sqrt(pfm**2-pav**2)
         if(plim1.ge.pkmx)plim1=pkmx
         call tmpgrid(zero,plim1,grang)
         call angprop(1,1,pav,pfm,efm,gmasq,kup,ptab,etab,estart,prorad)
         call legprop(mup,lp,pkmx,grang,prorad,propag)
      endif
      if(plim1.ge.pkmx)then
         goto 400
      else
         plim2=pav+pfm
         if(plim1.lt.plim2)then
            if(plim2.ge.pkmx)plim2=pkmx
            call tmpgrid(plim1,plim2,grang)
            call angprop(2,1,pav,pfm,efm,gmasq,kup,ptab,etab
     &,                                                 estart,prorad)
             call legprop(mup,lp,pkmx,grang,prorad,propag)
         endif
      endif
      if(plim2.lt.pkmx)then
         call tmpgrid(plim2,pkmx,grang)
         call angprop(3,1,pav,pfm,efm,gmasq,kup,ptab,etab,estart,prorad)
         call legprop(mup,lp,pkmx,grang,prorad,propag)
      endif
400   continue

      return
      end


********************************************************from polyprop***
      subroutine tmpgrid(plow,phigh,grang)
*
*  Makes a momentum mesh and uses COMMON 'kmesh' to communicate it.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(ngs=5)
      common / kmesh/ mxk,pk(ktot)
      common / gauss/ gweigh(ngs),gabsc(ngs)

      gcentr=.50d0*(plow+phigh)
      grang=.50d0*(phigh-plow)
      mxk=2*ngs
      do 100 i=1,ngs
         ginc=grang*gabsc(ngs+1-i)
         pk(i)=gcentr-ginc
         pk(mxk+1-i)=gcentr+ginc
100   continue

      return
      end


********************************************************from polyprop***
      subroutine angprop(n,lp,pav,pfm,efm,gmasq,kup,ptab,etab,estart
     &,                  prorad)
*
*  Obtains the angle averaged (lp=1) propagator difference 'prorad(k)'
*  for each mesh point 'pk(k)'. The nuclear matter spectrum is inter-
*  polated from the table 'etab()' by calling the function 'energy'.
*  'n' serves as the index for the integration regions.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(imx=16)
      dimension ptab(ktot),etab(ktot),prorad(ktot),psing(2),xcos(0:imx)
     &,         proang(0:imx)
      common / kmesh/ mxk,pk(ktot)
      external  energy
*
* integrand for SIMPSON rule, 'refpro' is the reference propagator.
*
      do 100 k=1,mxk
         refpro=1.0d0/(pk(k)**2/efm+gmasq)
         if(n.eq.1)then
            prorad(k)=refpro
         else
            psqsum=pav**2+pk(k)**2
            pprod=pav*pk(k)*2.0d0
            if(n.eq.2)then
               xinc=(psqsum-pfm**2)/pprod/imx
            else if(n.eq.3)then
               xinc=1.0d0/imx
            endif
            do 200 j=0,imx
               xcos(j)=j*xinc
               psing(1)=sqrt(psqsum+pprod*xcos(j))
               psing(2)=sqrt(psqsum-pprod*xcos(j))
               proang(j)=1.0d0/(energy(ktot,kup,ptab,etab,psing)-estart)
200         continue
            call simpson(0,imx,1,imx,xcos,proang,proang(0),prorad(k))
            prorad(k)=refpro-prorad(k)
         endif
100   continue

      return
      end


********************************************************from polyprop***
      subroutine legprop(mup,lp,pkmx,grang,prorad,propag)
*
*  Transformation to Legindgaard basis in the current subinterval, given
*  in COMMON 'kmesh'. Output is 'propag()', the propagator accumulated
*  over the previous subintervals.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      dimension prorad(ktot),propag(lpmx,mmx,mmx),func(ktot)
      common / polys/ plegin(mmx,ktot)
      common / kmesh/ mxk,pk(ktot)

      call polygen(mup,pkmx)
      do 100 mi=1,mup
         do 200 mf=1,mi
            do 300 k=1,mxk
               func(k)=prorad(k)*plegin(mi,k)*plegin(mf,k)*pk(k)**2
300         continue
            call gaussint(grang,func,proleg)
            propag(lp,mf,mi)=propag(lp,mf,mi)+proleg
            propag(lp,mi,mf)=propag(lp,mf,mi)
200      continue
100   continue

      return
      end


********************************************************from polyprop***
      double precision function energy(kmx,kup,ptab,etab,prel)
*
*  Interpolates function 'energy(x)' at x='pki' from the table 'etab()'
*  over abcissa 'ptab()'. If 'pki' is not element of the table 'ptab(',
*  we find the closest entry index 'kclose' {call to function 'ksearch'}
*  and interpolate from the values 'etab(k)', k='kclose'-'idown'...
*  'kclose'+'kup' by polynomial intrapolation {call to 'approxp'}. The
*  order is 'imx'-1 = 'idown'+'iup'.
*  INPUT     ptab() - abscissa table (single particle momenta)
*            etab() - function table (single particle spectrum)
*            kup - number of table entries
*            pki - abscissa value to interpolate to
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(idown=1,iup=1,imx=1+idown+iup)
      dimension prel(2),ptab(kmx),etab(kmx),p(imx),e(imx),c(imx),d(imx)
      external  approxp,ksearch
      save kguess
      data kmin,kguess/1,1/

      energy=0.0d0
      do 100 j=1,2
         pki=prel(j)
         kclose=ksearch(kmx,kup,ptab,pki,kguess)
         if(ptab(kclose).eq.pki)then
            energy=energy+etab(kclose)
         else
            if(kclose+iup.gt.kup)then
               ilo=kup-imx-1
               iclose=imx-(kup-kclose)
            else if(kclose-idown.lt.kmin)then
               ilo=kmin-1
               iclose=kmin+(idown-kclose)
            else
               ilo=kclose-idown
               iclose=idown
            endif
            do 200 i=1,imx
               p(i)=ptab(i+ilo)
               e(i)=etab(i+ilo)
               c(i)=e(i)
               d(i)=e(i)
200         continue
            energy=energy+approxp(imx,iclose,p,pki,e,c,d)
         endif
      kguess=kclose
100   continue

      return
      end


**********************************************************from energy***
      integer function ksearch(kmx,kup,ptab,pkval,kguess)
*
*     Searches table ptab() and finds index 'ksearch' of the entry
*     closest to pkval. 'kguess' is the starting value for the
*     bisection routine used. ptab() must be an ordered table.
************************************************************************
      implicit double precision (a-h,o-z)
      dimension ptab(kmx)

      klow=kguess-1
      khigh=kup+1
100   if((khigh-klow).gt.1)then
         kp=(khigh+klow)/2
         if((pkval.gt.ptab(kp)).eqv.(ptab(kup).gt.ptab(1)))then
            klow=kp
         else if(pkval.eq.ptab(kp))then
            khigh=kp
            klow=kp
         else
            khigh=kp
         endif
      goto 100
      endif
      ksearch=klow

      return
      end


************************************************************************
      double precision function approxp(imx,istart,r,rex,f,g,h)
*
*  'approxp(rex)' is the value of the function 'f(r=rex)', interpolated
*  from i=1...imx values 'u(r(i))' as a polynomial of order 'imx'-1.
*  INPUT     imx - order of interpolating polynomial + 1
*            istart - r(istart) is closest to rex
*            r() - abscissa values
*            rex - abscissa value to interpolate
*            f() - functional values f(r()) - g(),d() = f()
************************************************************************
      implicit  double precision (a-h,o-z)
      dimension r(imx),f(imx),g(imx),h(imx)

      fex=f(istart)
      istart=istart-1
      do 100 k=1,imx-1
         do 200 i=1,imx-k
            rad1=r(i)-rex
            rad2=r(i+k)-rex
            diff=g(i+1)-h(i)
            drad=rad1-rad2
            if(drad.eq.0.) stop 'approxp: Extrapolation failed.'
            drad=diff/drad
            h(i)=rad2*drad
            g(i)=rad1*drad
200      continue
         if(2*istart.lt.imx-k)then
            df=g(istart+1)
         else
            df=h(istart)
            istart=istart-1
         endif
         fex=fex+df
100   continue
      approxp=fex

      return
      end


*********************************************************from legprop***
      subroutine gaussint(dx,func,fintgrl)
*
*  GAUSS LEGENDRE integral 'fintgrl' of function 'func'.
*  Order is 'ngs'*2, weights & abscissas are stored in COMMON 'gauss'
*  the integration range is 'dx' on the underlying mesh stored
*  in COMMON 'kmesh'.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(ngs=5)
      dimension func(ktot)
      common / kmesh/ mxk,pk(ktot)
      common / gauss/ gweigh(ngs),gabsc(ngs)

      fintgrl=0.0d0
      idumy=ngs+1
      do 100 i=1,ngs
         fintgrl=fintgrl+gweigh(idumy-i)*(func(i)+func(idumy+ngs-i))
100   continue
      fintgrl=fintgrl*dx

      return
      end


************************************************************from MAIN***
      subroutine grtrafo(test,testid,ngroup,chanid,mup,pkmx
     &,                  efm,gmasq)
*
*   Transforms TR matrix elements from momentum to polynomial represent-
*   ation. 'lkmesh' (subroutine) reads the input file for each channel,
*   'polygen' calculates the Legindgaard polynomials on the the (input)
*   k-mesh, 'grlegin' performs the transformation.
*   INPUT   ngroup - transform spin/parity group 'ngroup'
*           mup - highest polynomial degree
*           test - .true. if inverse transformtion of TR is checked
*   OUTPUT  gref(no,mi,mf) - transformed TR in channel 'no' and
*             in polynomial 'mi' & 'mf', stored in COMMON 'grpoly'
*           plegin(m,pk) - Legindgaard polynomial 'm' at momentum
*             'pk'; the k-space range is up to pmax='pk(mxk)'. The
*             polynomials are stored in common 'polys'.
*   REMARK  The input files read by 'lkmesh' are checked for consistent
*   maximal momenta, which should be the same for all. If the momentum
*   stepsizes differ for different input files the flag 'kmix' is set to
*   'true'. The polynomials are then recalculated in each call, which is
*   time consuming (as is the need to do so again when we transform back
*   to the momentum representation later).
*   One can do a quality check by calling 'grtest', which inverts the
*   transformation and obtains a L2 norm w/ respect to the input TR.
*   It outputs the diagonal matrix elements (k-space) from the inversion
*   as series in m=mtest...mtest-minver and the L2 norms 'gsigma(m)'.
*   The output goes for each channel to a file w/ a prefix 'testid'.
*   'minver' is set in the data block w/ label 10.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      logical   kmix,knew,test
      character*6 chanid,testid*3
      dimension chanid(nmx),refg(ktot,ktot),grepol(mmx,mmx)
      common   /grpoly/gref(2*nqmx,mmx,mmx)
      common   /kmesh/ mxk,pk(ktot)
      common   /polys/ plegin(mmx,ktot)
      common   /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)
      data      knew/.true./
10    data      minver/2/
      mtest=mup

      do 110 i=1,ndim(ngroup)
        do 120 j=1,i
          n=abs(link(ngroup,i,j))
          if(n.ne.0)then

* read TR in channel 'n' from file 'chanid(n)', i.e. 'refg()'
*
            call lkmesh(n,chanid(n),knew,kmix,pkmx,efm,gmasq,refg)
*
* if 1st run or momentum mesh changed calculate Legindgaard polynomials
*
            if(knew)call polygen(mup,pk(mxk))
*
* transform TR to polynomial basis, 'chanid(n)(6:6)' identifies crossed
* channels
*
            sign=n/dble(link(ngroup,i,j))
            call grlegin(mup,sign,refg,grepol)
            no=n-abs(link(ngroup,1,1))+1
            do 200 mi=1,mup
              do 300 mf=1,mup
                gref(no,mf,mi)=grepol(mf,mi)
300           continue
200         continue
*
* inversion: convergence test w/ respect to polynomial degree
*
            if(test)then
              call grtest(testid,chanid(n),sign,mup,mtest,minver,refg
     &,                   grepol)
            endif
          knew=.false.
          endif
120     continue
110   continue

      return
      end


*********************************************************from grtrafo***
      subroutine grlegin(mup,sign,refg,transg)
*
*  Legindgaard transformation of TR(k,k'). The character 'qualify' was
*  the last letter in the input file name and serves as indicator for
*  teh crossed coupled channels. The default here is 'qulidf'='C'.
*  INPUT   mup - maximum polynomial degree
*          sign - +/-1 for un-/coupled channels (|L-L'| = 2)
*          refg() - TR in momentum space representation
*  OUTPUT  transg() - TR in polynomial representation
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(gnorm=5.0660591821169D-02)
      dimension refg(ktot,ktot),transg(mmx,mmx),funcf(ktot),funci(ktot)
      common   /kmesh/ mxk,pk(ktot)
      common   /polys/ plegin(mmx,ktot)
*
* sign = -1 for crossed channels with |L-L'| = 2
*
        phase=gnorm*sign
*
* 2-d integration on kmesh via an incomplete SIMPSON rule, which does
* not include the origin pk=0.
*
      do 100 mi=1,mup
         do 200 mf=1,mup
            do 300 ki=1,mxk
               do 400 kf=1,mxk
                  funcf(kf)=pk(kf)*pk(kf)*plegin(mf,kf)*refg(kf,ki)
400            continue
               call simpson(1,ktot,1,mxk,pk,funcf,0.0d0,funco)
               funci(ki)=pk(ki)*pk(ki)*plegin(mi,ki)*funco
300         continue
            call simpson(1,ktot,1,mxk,pk,funci,0.0d0,funco)
            transg(mf,mi)=phase*funco
200      continue
100   continue

      return
      end


*********************************************************from grtrafo***
      subroutine lkmesh(n,chname,knew,kmix,pkmx,efm,gmasq,refg)
*
*  Reads parameters, momentum mesh and TR matrixelements from the
*  file for channel 'n', named 'chname'. The format is described
*  in the text section.
*  The program stops if any problems occur during reading.
*  INPUT   n,chname - channel number and file to read in
*          knew - =.true. for the first entry (see OUTPUT)
*          pkmx - subspace cut-off momentum
*          efm,gmasq - input reference spectrum parameters
*  OUTPUT  refg - TR matrixelements
*          pk() - momentum mesh, stored in COMMON 'kmesh' along
*            w/ upper and lower limits 'mxk     ','mxk'
*          knew - statusflag, =.false. if no new momentum mesh
*            appears in two subsequent calls
*          kmix - statusflag, =.true. if at least two different
*            meshes are recognized
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      logical   knew,kmix
      character *6 chname
      dimension refg(ktot,ktot),pkold(ktot)
      common  / kmesh/ mxk,pk(ktot)
      save  pkold
50    format(' lkmesh 15: Problems with file ',a6,' . Continue...')
51    format(' lkmesh 16: Data of file ',a6,' exceed array dimensions.')
52    format(' lkmesh 17: Inconsistent input. Current file: ',a6,' .')
*
* reading from file 'chname'
*
      open(unit=22,iostat=io,err=15,file=chname)
      read(22,*,iostat=io,err=15)efmin,gmain,mxk
      read(22,*,iostat=io,err=15)(pk(i),i=1,mxk)
      read(22,*,iostat=io,err=15)((refg(j,i),refg(i,j),i=1,j),j=1,mxk)
      close(22)
15    if(io.ne.0)then
         write(*,50) chname
         goto 200
      endif
*
* check the input data for consistency
*
16    if(mxk.gt.ktot) write(*,51) chname
17    if((pk(mxk).ne.pkmx).or.(efmin.ne.efm).or.(gmain.ne.gmasq))then
         write(*,52) chname
      endif
      do 100 i=1,mxk
         if((.not.knew).and.(pkold(i).ne.pk(i)))then
            knew=.true.
            kmix=.true.
         endif
         if(knew)pkold(i)=pk(i)
100   continue

200   return
      end


**********************************************from grtrafo & polyprop***
      subroutine polygen(mup,pkmx)
*
*  Legindgaard polynomials up to degree 'mup' w/ maximal momentum
*  'pkmx'='pk(mxk)' on mesh (COMMON kmesh). The values are stored
*  as 'plegin()' in COMMON 'polys'.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      common   /kmesh/ mxk,pk(ktot)
      common   /polys/ plegin(mmx,ktot)
      external polynom

      do 100 m=1,mup
         cnorm=sqrt((2.0d0*m+1.0d0)/pkmx**3)*(-1.0d0)**(m+1)
         do 200 k=1,mxk
            plegin(m,k)=polynom(m,pk(k)/pkmx,cnorm)
200      continue
100   continue

      return
      end


*********************************************************from polygen***
      double precision function polynom(m,parg,cnorm)
*
*  Legindgaard polynomial of order 'm'-1, scaled argument 'parg'=
*  'pk(k)'/'pmx', the norm 'cnorm' is also input.
*     REMARK   The polynomial involves factorials up to (2m)! . They are
*  predefined in a data block up to 2m=20. The program stops if this
*  limit is exceeded.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(mmax=10)
      dimension fac(0:20)
      data fac /1.0d0,1.0d0,2.0d0,6.0d0,2.4d1,1.2d2,7.2d2,5.04d3,4.032d4
     &,3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9,8.71782912d10
     &,1.307674368d12,2.0922789888d13,3.55687428096d14,6.402373705728d15
     &,1.21645100408832d17,2.432902008176640d18/

      if((m.gt.mmax).or.(parg.gt.1.0d0).or.(m.lt.1))then
15       stop 'polynom 15: Invalid arguments.'
      endif

      polynom=0.0d0
      do 100 n=0,m-1
         polynom=(-1.0d0)**(n+2)/fac(n)*fac(n+m+1)/fac(m-1-n)
     &           *parg**n/fac(n+2)+polynom
100   continue
      polynom=polynom*cnorm

      return
      end


***********************************************************************
      subroutine simpson(imin,imx,ilo,iup,x,f,df,fintgrl)
*
*  SIMPSON Rule for function 'f(x)' over abscissa array 'x()'.
*  Output is fintgrl=Integral of f(x) between 'x(ilo)' and 'x(iup)',
*  for 'df'=0, between 'x(ilo-1)' and 'x(iup)' for 'df'='f(mlo-1)'.
***********************************************************************
      implicit double precision (a-h,o-z)
      dimension x(imin:imx),f(imin:imx)
*
      oddsum=0.0d0
      evesum=0.0d0
      step=x(ilo+1)-x(ilo)
      do 100 i=ilo,iup-1,2
         oddsum=oddsum+f(i)
100   continue
      do 200 i=ilo+1,iup-2,2
         evesum=evesum+f(i)
200   continue
      oddsum=oddsum*4.0d0
      evesum=evesum*2.0d0
      fintgrl=(df+oddsum+evesum+f(iup))*step/3.0d0

      return
      end


************************************************************from MAIN***
      subroutine lgnuc(ngroup,mup,propag,grecor)
*
*  Calculate correction matrix Õ(1-TR(GR-GE))**-1 - 1þTR.
*  INPUT   ngroup - spin&parity group number
*          propag - propagator difference (from 'polyprop')
*  OUTPUT  grecor - correction term; as matrix consisting of blocks
*                   w/ the polynomial representation of each channel
************************************************************************
      implicit double precision (a-h, o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      dimension propag(lpmx,mmx,mmx),grecor(ip,ip),agref(ip,ip)
     &,         grcinv(ip,ip),indx(mmx*nqmx),hh(ip,ip)
      common / qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)
      common / grpoly/gref(2*nqmx,mmx,mmx)
*
* realize full correction matrix in polynomial representation
*
* ia= : What actual size does the correction matrix take ?
*
      ia=mup*ndim(ngroup)
      nchlo=abs(link(ngroup,1,1))
*
* calculate matrices TR & (GR-GE), 'agref()' & 'grcinv',
* determine elements in lower triangles, link() picks the correct
* angular momentum component TR
*
      do 100 nqi=1,ndim(ngroup)
         do 200 mi=1,mup
            i=(nqi-1)*mup+mi
            do 300 nqf=1,ndim(ngroup)
               ldum=link(ngroup,nqf,nqi)
               nchno=abs(ldum)-nchlo+1
               do 400 mf=1,mup
                  j=(nqf-1)*mup+mf
                  if(j.ge.i)then
                     if(ldum.ne.0)then
                        agref(j,i)=gref(nchno,mf,mi)
                     else
                        agref(j,i)=0.0d0
                     endif
                     if(lpauli(ngroup,1,nqf,nqi).eq.1)then
                        grcinv(j,i)=-propag(1,mf,mi)
                     else
                        grcinv(j,i)=0.0d0
                     endif
                  endif
                  hh(j,i)=-grcinv(j,i)
400            continue
300         continue
200      continue
100   continue
*
* calculate M=Õ1-TR(GR-GE)þ,'fillmat2' fills the upper triangle ensuring
* hermiticity
*
      call fillmat2(ia,ip,hh)
      call fillmat2(ia,ip,grcinv)
      call fillmat2(ia,ip,agref)
      call multmat(ia,ip,agref,grcinv,grecor)
      do 550 j=1,ia
         grecor(j,j)=grecor(j,j)+1.0d0
550   continue
*
* invert M, M**-1 is 'grcinv'
*
      do 700 i=1,ia
         do 800 j=1,ia
            grcinv(j,i)=0.0d0
800      continue
         grcinv(i,i)=1.0d0
700   continue
      call ludecom(grecor,ia,ip,indx,det)
      do 900 j=1,ia
         call lulineq(grecor,ia,ip,indx,grcinv(1,j))
900   continue
*
* calculate output 'grecor', which is ÕM**-1 - 1þTR
*
      do 1000 j=1,ia
         grcinv(j,j)=grcinv(j,j)-1.0d0
1000  continue
      call multmat(ia,ip,grcinv,agref,grecor)

      return
      end


************************************************************from MAIN***
      subroutine lgback(ngroup,chanid,mup,grecor)
*
*  Obtains correction matrix in momentum space representation and adds
*  TR to get the full nuclear matter T-matrix.
*  INPUT    ngroup - spin&parity group to be transformed back
*           mup - highest polynomial degree
*           grecor() - correction matrix for present spin&parity group
*  OUTPUT   T(k,k'), 'gkern()', written to the input file w/ TR but re-
*             named w/ prefix 'outid'(COMMON /info/).
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(gnorm=19.739208802179d0)
      logical  knew,kmix
      character *6 chanid,outname,outid *3,inid *3,memo *64
      dimension chanid(nmx),grecor(ip,ip),gkern(ktot,ktot)
      common /polys/ plegin(mmx,ktot)
      common /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)
      common /info/memo,efm,gmasq,pfm,pav,pkmx,emev,mdgr,kup,inid,outid
      common /kmesh/mxk,pk(ktot)
      kmix=.false.
50          format(a64,1x)
51          format(6(f8.4,1x),2(i3,1x),a3,1x,a3)
52          format(2(f8.4,1x),i3,1x)
53          format(d16.8,1x)
54          format(2(d16.8,1x))
55          format(' lgback: Could not open/close output file ',a10,'.')

      do 100 k=1,ndim(ngroup)
        do 200 l=1,k
          n=abs(link(ngroup,k,l))
          if(n.ne.0)then
*
* calculate correction matrix C=Õ(1-TR(GR-GE))**-1 - 1þTR in K-space
* locate the blocks for the angular momentum component 'n'
*
            ibegin=(l-1)*mup+1
            jbegin=(k-1)*mup+1
            phase=gnorm*dble(n)/dble(link(ngroup,k,l))
*
* read in again the reference matrix TR(n), 'gkern'; if the momentum
* meshes were different, we need to redo the Legindgaard polynomials.
*
            call lkmesh(n,chanid(n),knew,kmix,pkmx,efm,gmasq,gkern)
            if(kmix)call polygen(mup,pkmx)
*
* do inverse Legindgaard transformation, add correction and reference
* matrix elements to obtain nuclear matter reaction matrix 'gkern'
*
            do 400 ki=1,mxk
              do 500 kf=1,mxk
                ginv=0.0d0
                do 600 mi=1,mup
                  i=ibegin+mi-1
                  do 700 mf=1,mup
                    j=jbegin+mf-1
                    ginv=ginv+grecor(j,i)*plegin(mf,kf)*plegin(mi,ki)
700               continue
600             continue
                gkern(kf,ki)=phase*ginv+gkern(kf,ki)
500           continue
400         continue
*
* output in channel 'n'
*
            outname=outid(1:3)//chanid(n)(4:6)
            open(unit=33,iostat=io,err=15,file=outname)
15          continue
            if(io.eq.0)then
              write(33,50)memo
              write(33,51)efm,gmasq,pfm,pav,pkmx,emev,mup,kup,inid,outid
              write(33,52)efm,gmasq,mxk
              write(33,53)(pk(i),i=1,mxk)
              write(33,54)((gkern(kf,ki),gkern(ki,kf),ki=1,kf)
     &,                                               kf=1,mxk)
              close(33,iostat=io,err=15)
            else
              write(*,55) outname
            endif
          endif
200     continue
100   continue

      return
      end


***********************************************************from lgnuc***
      subroutine fillmat2(ia,ip,aa)
*  Input quadratic array 'aa()' w/ elements in lower triangle.
*  Output is thecomplete symmetyric array.
*     ip - physical dimension  ia - actual dimension
***********************************************************************
      implicit double precision (a-h,o-z)
      dimension aa(ip,ip)

      do 100 i=1,ia
         do 200 j=1,i-1
            aa(j,i)=aa(i,j)
200      continue
100   continue
      return
      end


***********************************************************from lgnuc***
      subroutine multmat(ia,ip,aa,bb,cc)
*  Input quadratic matrices aa & bb. Output is matrixproduct cc=aa*bb.
*     ip - physical dimension  ia - actual dimension
************************************************************************
      implicit double precision (a-h,o-z)
      dimension aa(ip,ip),bb(ip,ip),cc(ip,ip)

      do 100 i=1,ia
         do 200 j=1,ia
            dum=0.0d0
            do 300 k=1,ia
               dum=dum+aa(j,k)*bb(k,i)
300         continue
            cc(j,i)=dum
200      continue
100   continue
      return
      end


***********************************************************from lgnuc***
      subroutine ludecom(aa,na,np,indx,det)
*
* LU-decomposition of a square matrix Õc.f.Numerical Recipes,Ch.2þ.
* INPUT  aa - matrix to decompose
*        na - actual dimension of aa
*        np - physical dimension of aa array (np <= nparameter 'nmax')
* OUTPUT aa - decomposed input matrix
*        indx - permutation pointer for pivoting
*        det  - even/odd indicator for row changes
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(nmax=50)
      parameter(delta=1.0d-20)
      dimension aa(np,np),indx(na),vv(nmax)
*
* scale the columns
*
      det=1.0d0
      do 100 i=1,na
         amax=0.0d0
            do 200 j=1,na
               if(abs(aa(i,j)).gt.amax)amax=abs(aa(i,j))
200      continue
         if(amax.eq.0.0d0)stop 'Error ludecom: Matrix singular'
         vv(i)=1.0d0/amax
100   continue
*
* decomposition algorithm
*
      do 300 j=1,na
         do 400 i=1,j-1
            sum=aa(i,j)
            do 500 k=1,i-1
               sum=sum-aa(i,k)*aa(k,j)
500         continue
            aa(i,j)=sum
400      continue
         amax=0.0d0
         do 600 i=j,na
            sum=aa(i,j)
            do 700 k=1,j-1
               sum=sum-aa(i,k)*aa(k,j)
700         continue
            aa(i,j)=sum
            dum=vv(i)*abs(sum)
            if(dum.ge.amax)then
               imax=i
               amax=dum
            endif
600      continue
         if(j.ne.imax)then
            do 800 k=1,na
               dum=aa(imax,k)
               aa(imax,k)=aa(j,k)
               aa(j,k)=dum
800         continue
            det=-det
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(aa(j,j).eq.0.0d0)aa(j,j)=delta
         if(j.ne.na)then
            dum=1.0d0/aa(j,j)
            do 900 i=j+1,na
               aa(i,j)=aa(i,j)*dum
900         continue
         endif
300   continue

      return
      end


***********************************************************from lgnuc***
      subroutine lulineq(aa,na,np,indx,b)
*
*  Solves linear system  aa * x = b. Õc.f. Numerical Recipes, Ch. 2þ
*  INPUT     aa - the lu-decomposed of the original matrix aa.
*            na - actual dimension of system
*            np - physical dimension of aa array
*            indx - permutation pointer for pivoting (see ludecom.f)
*            b - the inhomogeneous vector
*  OUTPUT    b - solution vector x
*  REMARK    To be used after the lu decomposition of aa (ludecom.f).
************************************************************************
      implicit double precision (a-h,o-z)
      dimension aa(np,np),indx(na),b(na)
*
* search for non zero components of b & do forward substitution
*
      ii=0
      do 100 i=1, na
         ll=indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do 200 j=ii, i-1
               sum=sum-aa(i,j)*b(j)
200         continue
         else if (sum .ne. 0) then
            ii=i
         endif
         b(i)=sum
100   continue
*
* do back substitution & store solution b(i)
*
      do 300 i=na, 1, -1
         sum=b(i)
         if(i.lt.na) then
            do 400 j=i+1, na
               sum=sum-aa(i,j)*b(j)
400         continue
         endif
         b(i)=sum/aa(i,i)
300   continue

      return
      end


************************************************************from MAIN***
      subroutine connect(inid,chanid,njup)
*
*   Selects allowed matrix elements, generates filenames and internal
*   channel labels.
*   If a matrix element <jls|TR|JLS> is allowed,'link()' takes the label
*   of the in/out file for the proper TR-component as value, else it is
*   0. The internal label is taken < 0 for the tensor coupled channels.
*   The filenames are filename="identifying prefix" + "JLS", where the
*   spin quantum number S is replaced by 'crossid' for the tensor coup-
*   led channels. For an example see terminal output by subrt. 'show'.
*   For the angle averaged Pauli operator, 'lpauli()' takes the value 1
*   for diagonal matrix elements in L, else it is 0.
*   4 distinct groups for S = +/-1 and parity -1**L = +/-1 are created.
*   The selection rules are in subrt.'conmap','map' generates the names.
*   INPUT     njup - highest J quantum number
*             lpmx - number of Pauli couplings (global parameter)
*             inid - prefix of input files
*   OUTPUT    chanid() - input filenames w/ correct labelling
*             link() - label maps (to COMMON 'qnmap')
*             lpauli() - maps for Pauli operator (to COMMON 'qnmap')
*             ndim() - actual dimension of label map (to COMMON 'qnmap')
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      character*6 chanid,inid*3,crossid*1
      dimension chanid(nmx)
      common /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)
      data crossid/'C'/
*
* maps for TR, 4 distinct maps, 'ns' is spin quantum number
*
      ns=0
      call map(1,njup,ns,-1,inid,crossid,chanid)
      call map(2,njup,ns,1,inid,crossid,chanid)
      ns=1
      call map(3,njup,ns,-1,inid,crossid,chanid)
      call map(4,njup,ns,1,inid,crossid,chanid)
*
* map for Pauli couplings
*
      do 100 i=1,4
         do 200 j=1,ndim(i)
            do 300 k=1,j
               if(k.eq.j)lpauli(i,1,j,k)=1
300         continue
200      continue
100   continue
      end


*********************************************************from connect***
      subroutine map(k,njup,ns,npar,inid,crossid,chanid)
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      character*1 chaqn*3,crossid,inid*3,chanid*6
      dimension chanid(nmx),nqj(nqmx),nql(nqmx),linki(nqmx,nqmx)
      common /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)

      i=0
*
* construct angular momenta 'nqj'(J), 'nql'(L) from parity 'npar' & spin
*
      do 100 mj=0,njup
         do 200 ml=abs(mj-ns),abs(mj+ns),1
            if(((-1)**ml).eq.npar)then
               i=i+1
               nqj(i)=mj
               nql(i)=ml
            endif
200      continue
100   continue
      ndim(k)=0
      do 300 lp=1,lpmx
         call conmap(lp,ns,nqj,nql,i,linki)
         do 400 m=1,nqmx
            do 500 n=1,m
              link(k,m,n)=linki(n,m)
              if(linki(n,m).ne.0)then
                if(lp.eq.1)then
                    nocha=abs(linki(n,m))
                    write(chaqn,'(3i1)') nqj(m),nql(m),ns
                    if(linki(n,m).gt.0)then
                       ndim(k)=ndim(k)+1
                    else
                       write(chaqn,'(3i1)') nqj(m),nqj(m)-1,ns
                       chaqn(3:3)=crossid
                    endif
                   chanid(nocha)=inid(1:3)//chaqn
                 else
                   stop 'Angle averaged Pauli operator only.'
                 endif
               endif
500         continue
400      continue
300   continue
      return
      end


************************************************************from map****
      subroutine conmap(lp,ns,nqj,nql,nqup,linki)
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      dimension nqj(nqmx),nql(nqmx),linki(nqmx,nqmx)
      save i

      k=0
      lpdif=2*(lp-1)
      do 200 m=1,nqup
         do 300 n=1,m
            k=k+1
            linki(n,m)=0
            if(ns.eq.0)then
               if(lpdif.eq.0)then
                  if(nqj(m).eq.nqj(n))then
                     i=i+1
                     linki(n,m)=i
                  endif
               else if(lpdif.eq.abs(nql(m)-nql(n)))then
                  stop 'Angle averaged Pauli operator only.'
               endif
            else if(ns.eq.1)then
               if(lpdif.eq.0)then
                  if(nqj(m).eq.nqj(n))then
                     if(mod(nql(m)-nql(n),2).eq.0)then
                        i=i+1
                        if(nql(m).eq.nql(n))then
                           linki(n,m)=i
                        else
* negative number identifies coupled channels for lp=1
                           linki(n,m)=-i
                        endif
                     endif
                  endif
               else if(lpdif.eq.abs(nql(m)-nql(n)))then
                  stop 'Angle averaged Pauli operator only.'
               endif
            endif
300      continue
200   continue
      do 400 m=nqup+1,nqmx
         do 500 n=1,m
            linki(n,m)=0
500      continue
400   continue

      return
      end


************************************************************************
      subroutine show(ngrplo,ngrpup,chanid)
*
*  Display some information about parameters, input & output to standard
*  output.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      character *6 chanid(nmx),inid*3,outid*3,memo*64
      common /info/memo,efm,gmasq,pfm,pav,pkmx,emev,mdgr,kup,inid,outid
      common /qnmap/link(4,nqmx,nqmx),lpauli(4,lpmx,nqmx,nqmx),ndim(4)

40    format(//,' ARRAY DIMENSIONS :',/)
41    format(' Momentum Mesh Size KTOT',t28,': ',i4,t39,'Channel Numbers
     & NMX',t60,': ',i3)
42    format(' Polynomial Order MMX',t28,': ',i4,t39,'Pauli Coupling LPM
     &X',t60,': ',i3)
43    format(' Maximum Matrix Size IP',t28,': ',i4)
44    format(//,' PRESENT VALUES :',/)
45    format(' Fermi Momentum PFM',t28,': ',f8.4,' 1/fm .')
46    format(' Average Momentum PAV',t28,': ',f8.4,' 1/fm .')
47    format(' Cut-Off Momentum PKMX',t28,': ',f8.4,' 1/fm .')
48    format(' Starting Energy ESTART',t28,': ',f8.4,' * 41.47 MeV .')
49    format(' Healing Parameter GMASQ',t28,': ',f8.4,' 1/fmª2 .')
50    format(' Effective Mass EFM',t28,': ',f8.4,' .')
51    format(' Highest Polynomial MUP-1',t28,': ',i3,' .')
52    format(' Nuclear Matter Spectrum on ',i3,' points.')
53    format(//,' FILENAMES for Angular Momentum Channels :',//
     &,t8,'Label',7x,'Inputfile',3x,'Outputfile',4x,
     &'j  l  s',/)
54    format(8x,i3,6x,a10,3x,a10,5x,3(a,2x))
55    format(/,4x,'...',/)
      write(*,40)
      write(*,41) ktot,nmx
      write(*,42) mmx,lpmx
      write(*,43) ip
      write(*,44)
      write(*,45) pfm
      write(*,46) pav
      write(*,47) pkmx
      write(*,48) emev/4.147d1
      write(*,49) gmasq
      write(*,50) efm
      write(*,51) mdgr-1
      write(*,52) kup
      write(*,53)
      imax=abs(link(ngrpup,ndim(ngrpup),ndim(ngrpup)))
      do 100 i=abs(link(ngrplo,1,1)),imax
         write(*,54)i,chanid(i),outid(1:3)//chanid(i)(4:6)
     &,       chanid(i)(4:4),chanid(i)(5:5),chanid(i)(6:6)
100   continue
      write(*,55)

      return
      end


*********************************************************from grlegin***
      subroutine grtest(testid,name,sign,mup,mtry,minv,refg,grepol)
*
*  Inverted Legindgaard transformation (see subroutine 'grtrafo').
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=64,nmx=12,mmx=8,njmx=2,lpmx=1,nqmx=5,ip=mmx*nqmx)
      parameter(gnorm=19.739208802179d0)
      character *6 name,naminv,testid*3
      dimension refg(ktot,ktot),grepol(mmx,mmx),ginv(mmx,ktot,ktot)
     &,         sdiffe(ktot),sumdif(ktot),ssumme(ktot),sumsum(ktot)
     &,         gsigma(mmx),gsig(2)
      common   /kmesh/ mxk, pk(ktot)
      common   /polys/ plegin(mmx,ktot)

      do 100 m=mtry-minv,mtry
*
* inverse transformation w/ m-1 as highest polynomial degree
*
         do 200 i=1,mxk
            do 300 j=1,mxk
               gsum=0.0d0
               do 400 mi=1,m
                  do 500 mf=1,m
                     gfi=plegin(mf,j)*plegin(mi,i)
                     gsum=gsum+gfi*grepol(mf,mi)
500               continue
400            continue
               ginv(m,j,i)=gnorm*sign*gsum
300         continue
200      continue
*
* calculate normalized L2-norm
*
         do 600 j=1,mxk
            do 700 i=1,mxk
               sdiffe(i)=((ginv(m,i,j)-refg(i,j))*pk(i)*pk(j))**2
               ssumme(i)=(refg(i,j)*pk(i)*pk(j))**2
700         continue
            call simpson(1,ktot,1,mxk,pk,sdiffe,0.0d0,sumdif(j))
            call simpson(1,ktot,1,mxk,pk,ssumme,0.0d0,sumsum(j))
600      continue
         call simpson(1,ktot,1,mxk,pk,sumdif,0.0d0,gsig(1))
         call simpson(1,ktot,1,mxk,pk,sumsum,0.0d0,gsig(2))
         gsigma(m)=gsig(1)/gsig(2)
100   continue
*
* output specifics
*
51    format(i2,1x,f8.4,1x,7(f16.8,1x,:))
54    format(1x,a6,10x,'Inverse Legindgaard (mup=',i1,'): Diagonal TR.'
     &//,5x,'k',10x,'TR(original)',7x,7('TR(m=',i1,')':,10x)/)
55    format(/22x,'gsigma',1x,6(f16.8,1x,:))
57    format(' grtest 15: Could not open/close file: ',a6,'.')
      naminv=testid(1:3)//name(4:6)
      open(unit=11,iostat=io,err=15,file=naminv)
      write(11,54)naminv,mup,(m,m=mtry,mtry-minv,-1)
      do 800 k=1,mxk
         write(11,51)k,pk(k),refg(k,k),(ginv(m,k,k),m=mtry,mtry-minv,-1)
800   continue
      write(11,55) (gsigma(m),m=mtry,mtry-minv,-1)
      close(11,iostat=io,err=15)
15    if(io.ne.0) write(*,57) naminv

      return
      end


      subroutine greetings

50    format(//,t20,'***** T N M *****',//,' Nuclear Matter Effective In
     &teraction: The Reaction Matrix.',//,t10,'Martin Fuchs  &  Philip J
     &. Siemens',/,t10,'Department of Physics',/,t10,'Oregon State Unive
     &rsity',/,t10,'Corvallis, OR 97331',/,t10,'U.S.A.',//,t10,'June 15,
     & 1991',////)
      write(*,50)
      pause

      return
      end
      block data
      implicit double precision (a-h,o-z)
      parameter(zero=1.0d-14,ngs=5)
      common / gauss/ gweigh(ngs),gabsc(ngs)
      data (gabsc(i),i=1,ngs)/.1488743389d0,.4333953941d0,.6794095682d0
     &,                       .8650633666d0,.9739065285d0/
      data (gweigh(i),i=1,ngs)/.2955242247d0,.2692667193d0,.2190863625d0
     &,                        .1494513491d0,.0666713443d0/
      end


