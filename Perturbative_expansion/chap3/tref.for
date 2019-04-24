      program tref
      implicit double precision (a-h,o-z)
      parameter(ntot=500,ktot=50,nchmx=12)
      logical   display,hcheck
      character *6 chanid,chanid1,corefix*3,prefix*3
      dimension chanid(nchmx),chanid1(nchmx)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
10    data efm,gmasq,nchlo,nchup,prefix,(chanid1(i),i=1,nchmx)/
     &    1.0d0,1.4d0,1,12,'TRF','000','011','110','111','220'
     &,   '221','211','231','21C','101','121','10C'/
11    data display,hcheck,corefix /.true.,.true.,'TCR'/

      call greetings
*
* make filenames 'chanid(i)' for the i-th angular momentum channel
*
      do 100 i=1, nchup
         chanid(i)=prefix//chanid1(i)(1:3)
100   continue
*
* radial meshes in position and momentum space
*
      call rgrid
      call kgrid
*
* print information on screen
*
      if(display)then
         call show(efm,gmasq,nchlo,nchup,chanid)
      endif
*
* perform calculation loops in specified channels 'nchlo'...'nchup'
*
      call chanlop(efm,gmasq,nchlo,nchup,hcheck,corefix,chanid)

      end


**********************************************************************
      subroutine rgrid
*
*  Defines the radial mesh in position space.
*  INPUT noint   -  number of intervals w/ different stepsizes
*         numdf() -  number of points including the ()-th interval
*         minr - lower mesh boundary (should be 0)
*         maxr - upper mesh boundary
*         matr - matching point (for single channels)
*         mhco - position of the hard core
*         raddf() -  upper boundary of ()-th interval
*  OUTPUT rad(), along w/ 'minr,maxr,matr,mhco' in 'common /rmesh/'
*
*  The actual range of 'rad()' should not exceed 'rad(0:ntot-1)'.
**********************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)

      parameter(noint=1)
      dimension numdf(0:noint), raddf(0:noint)
      data (raddf(i),i=0,noint)/.1d-13,.1d2/
      data (numdf(i),i=0,noint)/0,200/
      mhco=1

      minr=numdf(0)
      maxr=numdf(noint)
      matr=maxr
      rad(minr)=raddf(0)
      do 100 i=0,noint-1
         rstep=(raddf(i+1)-raddf(i))/dble(numdf(i+1)-numdf(i))
         do 200 n=numdf(i),numdf(i+1),1
            rad(n)=raddf(i)+(n-numdf(i))*rstep
200      continue
100   continue

      return
      end


************************************************************from MAIN***
      subroutine kgrid
*
*     Defines the radial mesh in momentum space.
*     INPUT analogous to sbrt. 'rgrid'
*     OUTPUT momentum mesh 'kmesh', along w/ boundary indices
*            'mink', 'maxk' in COMMON /kmesh/ .
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ktot=50)
      common / kmesh/ mink,maxk,pk(ktot)

      parameter(noint=1)
      dimension numdf(0:noint),pkdf(0:noint)
      data (pkdf(i),i=0,noint)/.07d0,2.8d0/
      data (numdf(i),i=0,noint)/1,40/

      mink=numdf(0)
      maxk=numdf(noint)
      pk(mink)=pkdf(0)
      do 100 i=0,noint-1
         pstep=(pkdf(i+1)-pkdf(i))/dble(numdf(i+1)-numdf(i))
         do 200 n=numdf(i),numdf(i+1),1
            pk(n)=pkdf(i)+(n-numdf(i))*pstep
200      continue
100   continue

      return
      end


************************************************************from MAIN***
      subroutine chanlop(efm,gmasq,nchlo,nchup,hcheck,corefix,chanid)
*
*     Control routine to determine allowed matrix elements. We do the
*     uncoupled ones first and then the remaining coupled cases. The
*     internal labelling of the channels follows the scheme in MAIN !
*     The actual calculation and output follows 'klop1(2)' downwards.
*     INPUT    efm, gmasq - spectral parameters
*              nchlo - label of starting channel
*              nchup - label of last channel to calculate
*              chanid - set of filenames (where the output goes)
*     OUTPUT   j,l,s - channel quantum numbers (to subroutines)
*     REMARK   All 3 channels in a coupled case will be calculated
*              together. The label sequence is 'nchlo','nchlo+1',
*              '...+2'. A starting (last) value 'nchlo' ('nchup') MUST
*              be the first (last) one in this sequence.
**********************************************************************
      implicit double precision (a-h,o-z)
      parameter(nchmx=12)
      parameter(ktot=50,ntot=500)
      logical  hcheck
      character*6 chanid(nchmx),copid(3),corefix*3
      common / kmesh/ mink,maxk,pk(ktot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)

50    format(' Current channels: ',3(i4,1x,a6))
51    format(' Now calculating...'/)
53    format(i1)
      write(*,51)
      i=nchlo
*
* recover quantum numbers j,l,s from filenames
*
100   if(i.le.nchup)then
      read(chanid(i)(4:4),53)j
      read(chanid(i)(5:5),53)l
      if(chanid(i)(6:6).eq.'0')then
         ns=0
      else if(chanid(i)(6:6).eq.'1')then
         ns=1
      else
         ns=-1
      endif
*
* single channels
*
      if((ns.eq.0).or.((ns.eq.1).and.((j.eq.0).or.(j.eq.l))))then
         write(*,50) i, chanid(i)
      lold=l
         call klop1(hcheck,corefix,chanid(i),j,l,ns,efm,gmasq)
         i=i+1
      else if(((ns.eq.1).or.(ns.eq.-1)).and.(iabs(j-l).eq.1))then
*
* coupled channels
*
         copid(1)=chanid(i)
         copid(2)=chanid(i+1)
         copid(3)=chanid(i+2)
         write(*,50) i,chanid(i),i+1,chanid(i+1),i+2,chanid(i+2)
         call klop2(hcheck,corefix,copid,i,j,efm,gmasq)
         i=i+3
         endif
      goto 100
      endif

      return
      end


*********************************************************from chanlop***
      subroutine klop1(hcheck,corefix,name,j,l,ns,efm,gmasq)
*
*  Control routine to obtain a set of reaction matrix elements on the
*  2-d momentum mesh (common 'kmesh'). The full array is computed in
*  subrt. 'gref1'  and checked for hermiticity w/ respect to momentum.
*  The function 'trapcorr()' corrects the offset from zero-momentum for
*  momentum space integrations with teh subrt. 'trapez'.
*  INPUT   hcheck - .true. if hermiticity check to be performed
*          corefix - prefix for files w/ core contributions
*          name - name of output file
*          j,l,ns - channel quantum numbers
*          efm,gmasq - spectral parameters
*  OUTPUT  gsigma - measure of hermiticity w/ respect to momentum
*          gref(kf,ki) - the ME's <jls,kf|TR|jls,ki>, to file 'name'
************************************************************************
      implicit double precision (a-h,o-z)
      logical  hcheck
      character *6 name,corname,corefix*3
      parameter(ktot=50,ntot=500)
      dimension gref(ktot,ktot),gcore(ktot,ktot),gsigma(2)
     &,         funcf(ktot),funci(ktot)
      common / kmesh/ mink,maxk,pk(ktot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      trapcorr(x,y)=.5d0*y*x

50    format(32x,' Could not open file ',2(a10,1x))
51    format(32x,' Hermiticity:  gsigma = ',e16.8)
52    format(2(f8.4,1x),i4)
53    format(d16.8)
54    format(2(d16.8,1x))

      call gref1(j,l,ns,efm,gmasq,gcore,gref)
*
* Check hermiticity pk(ki) <-> pk(kf).
*
      call gref1(j,l,ns,efm,gmasq,gcore,gref)
      if(hcheck)then
         do 300 i=1,2
            dumy=(-1)**i
            do 400 ki=mink,maxk
               do 500 kf=mink,ki
                  funcf(kf)=((gref(kf,ki)+dumy*gref(ki,kf))
     &                      *pk(ki)*pk(kf))**2
500            continue
               call trapez(mink,maxk,mink,ki,pk,funcf,funco)
               funci(ki)=trapcorr(pk(mink),funcf(mink))+funco
400         continue
            call trapez(mink,maxk,mink,maxk,pk,funci,funco)
            gsigma(i)=trapcorr(pk(mink),funci(mink))+funco
300      continue
         gsigma(1)=gsigma(1)/gsigma(2)
         write(*,51) gsigma(1)
      endif

      open(unit=11,iostat=io,file=name)
      corname=corefix//name(4:6)
      open(unit=22,iostat=io,file=corname)
      if(io.ne.0)then
         write(*,50) name,corname
         goto 600
      endif
      write(11,52) efm,gmasq,maxk-mink+1
      write(11,53) (pk(i),i=mink,maxk)
      write(11,54) ((gref(kf,ki),gref(ki,kf),ki=mink,kf),kf=mink,maxk)
      write(22,52) efm,gmasq,maxk
      write(22,53) (pk(i),i=mink,maxk)
      write(22,54) ((gcore(kf,ki),gcore(ki,kf),ki=mink,kf),kf=mink,maxk)
      close(11)
      close(22)

600   lold=l
      return
      end


*********************************************************from chanlop***
      subroutine klop2(hcheck,corefix,name,n,j,efm,gmasq)
*
*  Control routine to obtain a set of reaction matrix elements for
*  coupled channels. See subroutine 'klop1'.
*  REMARK   'm' is used as label for matrixelements
*            m=1  <j-1|TR|j-1>    3  <j+1|TR|j-1>
*              2  <j+1|TR|j+1>    4  <j-1|TR|j+1>
*            gsigma(1)  Hermiticity check ki<->kf for m=1
*            gsigma(2)  ...for m=2
*            gsigma(3)  ...for m=3 & 4 (crossed channels)
************************************************************************
      implicit double precision (a-h,o-z)
      logical   hcheck
      character*6 name(3),chname,corename,corefix*3
      parameter(ktot=50)
      dimension gref(4,ktot,ktot),gcore(2,ktot,ktot)
     &,         gsigma(2),funcf(ktot),funci(ktot)
      common / kmesh/ mink,maxk,pk(ktot)
      trapcorr(x,y)=.5d0*x*y

50    format(32x,' Could not open file ',a10)
51    format(27x,' Hermiticity:  gsigma(',a3,') = ',e16.8)
52    format(2(f8.4,1x),i4)
53    format(d16.8)
54    format(2(d16.8,1x))
*
* Find matrix elements and check hermiticity for TR(l=j-1)&TR(l=j+1)
*
      call gref2(j,efm,gmasq,gcore,gref)
      if(hcheck)then
         do 400 m=1,3
            do 500 i=1,2
               dumy=(-1)**i
               do 600 ki=mink,maxk
                  if(m.lt.3)then
                     do 700 kf=mink,ki
                        funcf(kf)=((gref(m,kf,ki)+dumy*gref(m,ki,kf))
     &                             *pk(ki)*pk(kf))**2
700                  continue
                  else
                     do 710 kf=mink,ki
                        funcf(kf)=((gref(3,kf,ki)+dumy*gref(4,ki,kf))
     &                             *pk(ki)*pk(kf))**2
710                  continue
                 endif
                 call trapez(mink,maxk,mink,ki,pk,funcf,funco)
                 funci(ki)=trapcorr(pk(mink),funcf(mink))+funco
600            continue
               call trapez(mink,maxk,mink,maxk,pk,funci,funco)
               gsigma(i)=trapcorr(pk(mink),funci(mink))+funco
500         continue
            gsigma(1)=gsigma(1)/gsigma(2)
            write(*,51) name(m)(4:6),gsigma(1)
400      continue
      endif
*
* output
*
      do 800 m=1,3
         chname=name(m)
         corename=corefix//chname(4:6)
         no=m+22
         np=m+33
         open(unit=no,file=chname,iostat=io)
         open(unit=np,file=corename,iostat=io)
         if(io.ne.0)then
            write(*,50)chname
            goto 1000
         endif
         write(no,52) efm,gmasq,maxk-mink+1
         write(np,52) efm,gmasq,maxk-mink+1
         write(no,53) (pk(i),i=mink,maxk)
         write(np,53) (pk(i),i=mink,maxk)
         write(no,54) ((gref(m,kf,ki),gref(m,ki,kf),ki=mink,kf)
     &,                                           kf=mink,maxk)
         if(m.lt.3)write(np,54)((gcore(m,kf,ki),gcore(m,ki,kf)
     &,                          ki=mink,kf),kf=mink,maxk)
         close(no)
         close(np)
800   continue

      do 900 m=4,4
         chname=name(3)(1:4)//'CC'
         no=m+22
         open(unit=no,file=chname,iostat=io)
         if(io.ne.0)then
            write(*,50)chname
            goto 1000
         endif
         write(no,52) efm,gmasq,maxk-mink+1
         write(no,53) (pk(i),i=mink,maxk)
         write(no,54) ((gref(m,kf,ki),gref(m,ki,kf),ki=mink,kf)
     &,                                           kf=mink,maxk)
         close(no)
900   continue

1000  return
      end


***********************************************************from klop1***
      subroutine gref1(j,l,ns,efm,gmasq,gcore,gref)
*
*  Calculates the reaction matrix  <jls,kf|TR|jls,ki> for uncoupled
*  channels. The core scatterd waves 'fchan()' and the core contribution
*  'gcore()' to TR are obtained. The call to 'funcref1' returns the wave
*  functions 'urad()' and 'trapez' does the integration on the radial
*  mesh 'rmesh'. For the same 'l', 'gcore' and 'fchan' remain the same
*  and are saved.
*  OUTPUT  gcore() - core contribution
*          gref()  - outer contribution to <jls|TR|jls>
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500,ktot=50)
      parameter(fourpi=.125663706143592d2)
      dimension gcore(ktot,ktot),gref(ktot,ktot),fchan(0:ntot,ktot)
     &,         urad(0:ntot),pot(0:ntot),funci(0:ntot)
      dimension bessl(0:ntot,ktot),bhank(0:ntot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
      external bessj,bessh,dirpot1
      save     pot,fchan,nsold,lold,jold,bessl,bhank
      data     nsold,lold,jold,kiold/-1,-1,-1,-1/
*
* obtain Bessel & Hankel functions as arrays 'bessl()' & 'bhank()'
*
      gma=sqrt(gmasq)
      if(l.ne.lold)then
         do 100 i=minr,maxr
            bhank(i)=bessh(l,gma*rad(i))
            do 200 k=mink,maxk
               bessl(i,k)=bessj(l,pk(k)*rad(i))
200         continue
100      continue
*
* final state 'fchan()' and core volume & surface terms 'ghco' & 'gedg'
*
         yc=gma*rad(mhco)
         do 300 ki=mink,maxk
            hco=bessl(mhco,ki)/bhank(mhco)
            do 400 i=mhco,maxr
               fchan(i,ki)=bessl(i,ki)-hco*bhank(i)
400         continue
            pki=pk(ki)
            hco=gma*bessl(mhco,ki)/bhank(mhco)
            gmakk=gmasq+pki**2
            pc=pki*rad(mhco)
            do 500 kf=mink,maxk
               do 600 i=minr,mhco
                  funci(i)=bessl(i,ki)*bessl(i,kf)
600            continue
               call trapez(minr,maxr,minr,mhco,rad,funci,ghco)
               ghco=gmakk*ghco
               dbessh=hco*((l+1)*bhank(mhco)/yc-bessh(l+1,yc))
               dbessj=pki*((l+1)*bessl(mhco,ki)/pc-bessj(l+1,pc))
               gedg=bessl(mhco,kf)*(dbessj-dbessh)
               gnorm=fourpi/pki/pk(kf)
               gcore(kf,ki)=gnorm*(gedg+ghco)
500         continue
300      continue
      endif
*
* interaction: local potential 'pot()'
*
      if((l.ne.lold).or.(j.ne.jold).or.(ns.ne.nsold))then
         do 700 i=mhco,maxr
            pot(i)=dirpot1(j,l,ns,rad(i))
700      continue
      endif
*
* scattered reference wave function 'urad()' & outer term 'refg'
*
      do 800 ki=mink,maxk
         do 900 kf=mink,maxk
            gnorm=fourpi/pk(ki)/pk(kf)
            if(ki.ne.kiold)then
               call funcref1(j,l,ns,ki,efm,gmasq,pot,bessl,urad)
            endif
            do 1000 i=mhco,maxr
               funci(i)=fchan(i,kf)*pot(i)*urad(i)
1000        continue
*
* integrate and normalize, 'trapez' = trapezoidal rule
*
            call trapez(minr,maxr,mhco,maxr,rad,funci,refg)
            gref(kf,ki)=gnorm*refg+gcore(kf,ki)
         kiold=ki
900      continue
800   continue

      nsold=ns
      lold=l
      jold=j
      return
      end


***********************************************************from klop2***
      subroutine gref2(j,efm,gmasq,gcore,gref)
*
*  Reaction Matrix 'gref()'for coupled channels. The contributions
*  from the core and the outer region are computed.'funcref2' returns
*  the wavefunctions in both channels. To check for hermiticity, both
*  off-diagonal ME's (w/ respect to l) are evaluated.
*  OUTPUT gref(1,kf,ki) <j j-1|TR|j j-1>  gref(3,kf,ki) <j j+1|TR|j j-1>
*           ...2...     <j j+1|TR|j j+1>    ...4...     <j j-1|TR|j j+1>
*         gcore() - similar to 'gref()', core contributions only
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500,ktot=50)
      parameter(fourpi=.125663706143592d2)
      dimension gcore(2,ktot,ktot),gref(4,ktot,ktot),wrad(2,2,0:ntot)
     &,         pot(2,0:ntot),fchan(2,0:ntot),funci(0:ntot),ghco(2)
     &,         gedg(2),hco(2),bessl(2,0:ntot,ktot),bhank(2,0:ntot)
     &,         hou(2),cpot(0:ntot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
      external bessj,bessh,dirpot2,tenpot
      data     kiold/-1/
*
* interactions : 'pot()' - direct, 'cpot()' - coupling
* Bessel & Hankel functions as arrays 'bessl()' & 'bhank()'
*
      gma=sqrt(gmasq)
      lm=j-1
      lp=j+1
      do 100 i=minr, maxr
         x=gma*rad(i)
         pot(1,i)=dirpot2(j,lm,rad(i))
         pot(2,i)=dirpot2(j,lp,rad(i))
         cpot(i) =-tenpot(j,rad(i))
         bhank(1,i)=bessh(lm,x)
         bhank(2,i)=bessh(lp,x)
         do 200 k=mink,maxk
            x=pk(k)*rad(i)
            bessl(1,i,k)=bessj(lm,x)
            bessl(2,i,k)=bessj(lp,x)
200      continue
100   continue

      yc=rad(mhco)*gma
      do 300 ki=mink,maxk
         pki=pk(ki)
         hco(1)=gma*bessl(1,mhco,ki)/bhank(1,mhco)
         hco(2)=gma*bessl(2,mhco,ki)/bhank(2,mhco)
         pc=rad(mhco)*pki
         gmakk=gmasq+pki**2
         do 400 kf=mink,maxk
            pkf=pk(kf)
*
* core volume & core edge corrections 'ghco()' & 'gedg()'
*
            do 500 lc=1,2
              l=j+(-1)**lc
              do 600 i=minr,mhco
                 funci(i)=bessl(lc,i,ki)*bessl(lc,i,kf)
600           continue
              call trapez(minr,maxr,minr,mhco,rad,funci,refg)
              ghco(lc)=gmakk*refg
              dbessj=pki*((l+1)*bessl(lc,mhco,ki)/pc-bessj(l+1,pc))
              dbessh=hco(lc)*((l+1)*bhank(lc,mhco)/yc-bessh(l+1,yc))
              gedg(lc)=bessl(lc,mhco,kf)*(dbessj-dbessh)
*
* core scattered waves 'fchan()'
*
              hou(lc)=bessl(lc,mhco,kf)/bhank(lc,mhco)
              do 700 i=mhco,maxr
                 fchan(lc,i)=bessl(lc,i,kf)-hou(lc)*bhank(lc,i)
700           continue
500         continue
*
* reference wave functions: 'wrad(1..)' - l=j-1 is entrance channel
*                           'wrad(2..)' - l=j+1 is entrance channel
*
            if(ki.ne.kiold)then
               call funcref2(j,efm,gmasq,ki,pot,cpot,bessl,wrad)
            endif

            gnorm=fourpi/pki/pkf
            do 800 lc=1,2
               ll=3-lc
*
* diagonal matrix elements 'refg', with core correction
*
               do 900 i=mhco,maxr
                  funci(i) = fchan(lc,i)*(pot(lc,i)*wrad(lc,lc,i)
     &                      -cpot(i)*wrad(lc,ll,i))
900            continue
               call trapez(minr,maxr,mhco,maxr,rad,funci,refg)
               gref(lc,kf,ki)=(ghco(lc)+gedg(lc)+refg)*gnorm
               gcore(lc,kf,ki)=(gedg(lc)+ghco(lc))*gnorm
*
* off-diagonal matrix elements
*
               do 1000 i=mhco,maxr
                  funci(i) = fchan(ll,i)*(pot(ll,i)*wrad(lc,ll,i)
     &                      -cpot(i)*wrad(lc,lc,i))
1000           continue
               call trapez(minr,maxr,mhco,maxr,rad,funci,refg)
               gref(lc+2,kf,ki)=refg*gnorm
800         continue
         kiold = ki
400      continue
300   continue

      return
      end


***********************************************************from gref1***
      subroutine funcref1(j,l,ns,ki,efm,gmasq,pot,bessl,urad)
*
*  Finds the radial reference wavefunction in uncoupled channel (l,j).
*  Here we provide the full potential 'potfn' and the driving term
*  'ginfn' (i.e. w/ a Besselfunction). The boundary value problem is
*  initialized and the call 'usolve1' then solves for the wavefunction.
*  OUTPUT  urad() - radial wavefunction for momentum 'pk(ki)'.
*     ------ For NUMEROV w/ a 1-parameter NEWTON ITERATION. ------
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (ktot=50,ntot=500)
      parameter (ndim=1)
      dimension pot(0:ntot),urad(0:ntot),potfn(0:ntot),ginfn(0:ntot)
     &,         du(1),bessl(0:ntot,ktot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
      save      potfn,ginfn,lold,jold
      data      lold,jold,nsold/-1,-1,-1/
*
* set up the full potential: interaction + angular momentum term
*
      if((l.ne.lold).or.(j.ne.jold).or.(ns.ne.nsold))then
         dumy=-dble(l*(l+1))
         do 100 i=mhco,maxr
            potfn(i)=dumy/rad(i)**2-gmasq-efm*pot(i)
100      continue
      endif
*
* the driving term
*
         dumy=-(gmasq+pk(ki)**2)
         do 200 i=mhco,maxr
            ginfn(i)=dumy*bessl(i,ki)
200      continue
*
* initialize wavefunction (for a NUMEROV algorithm) and solve, 'du' is
* the trial starting value at the lower boundary
*
      i=l+1
      du(1)=1.0d-4
      urad(mhco)=0.0d0
      urad(mhco+1)=du(1)
      urad(maxr)=bessl(maxr,ki)
      call usolve1(ndim,du,mhco,potfn,ginfn,urad)

      lold=l
      jold=j
      nsold=ns
      return
      end


***********************************************************from gref2***
      subroutine funcref2(j,efm,gmasq,ki,pot,cpot,bessl,wrad)
*
*  Finds radial wavefunctions 'wrad()' for coupled channels. The
*  boundary value problem is initialized wnd solved by trial-integration
*  to a fitting  point ('matr',here: ='maxr') and linear improvement on
*  the initial conditions 'du()'.
*  INPUT    ki - index of entrance channel momentum
*  OUTPUT   wrad - scattered radial wavefunctions, see subrt. 'gref2'
*     ---- for NUMEROV w/ a 2-parameter NEWTON iteration -----
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (ntot=500,ktot=50)
      parameter (ndim=2)
      dimension pot(2,0:ntot),cpot(0:ntot),wrad(2,2,0:ntot)
     &,         potfn(2,0:ntot),ginfn(2,0:ntot),urad(2,0:ntot)
     &,         bessl(2,0:ntot,ktot),cpotfn(0:ntot),du(ndim)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
      save      potfn,cpotfn,ginfn,jold
      data      jold/-1/
*
* full potential in channels: interaction + angular momentum term
*
      if(j.ne.jold)then
         dumy1=-dble(j*(j-1))
         dumy2=-dble((j+1)*(j+2))
         do 100 i=mhco, maxr
            potfn(1,i)=dumy1/rad(i)**2-gmasq-efm*pot(1,i)
            potfn(2,i)=dumy2/rad(i)**2-gmasq-efm*pot(2,i)
            cpotfn(i)=efm*cpot(i)
100      continue
      endif
*
* driving terms, loop over the l=j+1 & j-1 entrance channels
*
      dumy1=-(gmasq+pk(ki)**2)
      do 200 lc=1,2
         l=j+(-1)**lc
         ll=3-lc
         do 300 i=mhco,maxr
              ginfn(lc,i)=dumy1*bessl(lc,i,ki)
           ginfn(ll,i)=0.0d0
300      continue
*
* initialize wavefunctions 'urad()' (NUMEROV algorithm) & solve
*
         do 400 i=minr,mhco
            urad(1,i)=0.0d0
            urad(2,i)=0.0d0
400      continue
         urad(1,mhco+1)=1.0d-3
         urad(2,mhco+1)=1.0d-3
         urad(lc,maxr)=bessl(lc,maxr,ki)
         urad(ll,maxr)=0.0d0
         du(1)=+.5d-3
         du(2)=-.5d-3
         call usolve2(ndim,du,mhco,potfn,cpotfn,ginfn,urad)
         do 500 i=minr,maxr
            wrad(lc,1,i)=urad(1,i)
            wrad(lc,2,i)=urad(2,i)
500      continue
200   continue

      jold = j
      return
      end


********************************************************from funcref1***
      subroutine usolve1(ndim,du,mpos,potfn,ginfn,urad)
*
*     Like subroutine usolve2, but for single channels. DIFFERENCES:
*     Input/ output potfn()...urad() are 1dimensional arrays.
*     Integration by call to nrv1.
*     INPUT (from funcref1)
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500)
      parameter(newtdm=1,dutiny=.5d-4,uacc=.5d-4)
      dimension potfn(0:ntot),ginfn(0:ntot),urad(0:ntot),udev(newtdm)
     &,         udevn(newtdm),dudev(newtdm,newtdm),uin(newtdm)
     &,         du(ndim),uou(newtdm)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      logical   encore
      encore=.true.
*
* Treat like 2 initial value problems at upper and lower boundaries.
* The outwards integrated solution is obtained by integrating to
* the fitting point ('matr'), where  only the function itself is
* required to be continuous w/ respect to the inward solution.
*
* evaluation at upper boundary
*
      uin(1)=urad(matr)
*
* trial outward integration from lower boundary
*
      mdir=1
      mend=matr
      call nrv1(mpos,mdir,mend,potfn,ginfn,urad)
      uou(1)=urad(matr)
      call check(newtdm,uin,uou,udev)
*
* now improving w/ NEWTON's method, newtdm=1 conditions.
*
100   continue
      do 110 i=1,newtdm
         uold=urad(mpos+1)
         urad(mpos+1)=urad(mpos+1)+du(i)
         call nrv1(mpos,mdir,mend,potfn,ginfn,urad)
         uou(1)=urad(matr)
         call check(newtdm,uin,uou,udevn)
         do 120 j=1,newtdm
            dudev(j,i)=(udevn(j)-udev(j))/du(i)
120      continue
         urad(mpos+1)=uold
110   continue
      dum=dudev(1,1)
      if(dum.eq.0.0)then
         du(1)=du(1)+dutiny
         goto 100
      endif
      du(1)=-udev(1)/dum
      urad(mpos+1)=urad(mpos+1)+du(1)
*
* improved trial; since we have a linear problem, this should yield
* the solution. If not, we dare to try NEWTON one more time.
*
      call nrv1(mpos,mdir,mend,potfn,ginfn,urad)
      uou(1)=urad(matr)
      call check(newtdm,uin,uou,udev)
50    format(' 50 usolve1: Matching failure:  udev(1) = ',e16.8)
      if(abs(udev(1)).gt.uacc)then
         if(encore)then
            encore=.false.
            goto 100
         else
            write(*,50) udev(1)
         endif
      endif

      return
      end


********************************************************from funcref2***
      subroutine usolve2(ndim,du,mpos,potfn,cpotfn,ginfn,urad)
*
*  Determines reference wave function urad(). The initialization for
*  each entrance channel is done in 'funcref2'. To integrate, a Numerov
*  algorithm is used {call to 'nrv2'}. See comments for method.
*  INPUT  du() - trial initial values
*         mpos - starting index on radial mesh common/rmesh/
*  OUTPUT urad() - reference wave function
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500)
      parameter(newtdm=2,dutiny=.5d-4,uacc=.5d-4)
      dimension potfn(2,0:ntot),cpotfn(0:ntot),ginfn(2,0:ntot),du(ndim)
     &,         urad(2,0:ntot),udev(newtdm),udevn(newtdm)
     &,         dudev(newtdm,newtdm),uin(newtdm),uou(newtdm)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      logical   encore
      encore=.true.
*
* Treat like 2 initial value problems at upper and lower boundaries.
* The outwards integrated solution is obtained by integrating to
* the fitting point ('matr'), where  only the function itself is
* required to be continuous w/ respect to the inward solution.
*
* evaluation at upper boundary
*
      uin(1)=urad(1,matr)
      uin(2)=urad(2,matr)
*
* trial outward integration from lower boundary
*
      mdir=1
      mend=matr
      call nrv2(mpos,mdir,mend,potfn,cpotfn,ginfn,urad)
      uou(1)=urad(1,matr)
      uou(2)=urad(2,matr)
      call check(newtdm,uin,uou,udev)
*
* now improving w/ NEWTON's method, newtdm=2 conditions.
*
100   continue
      do 110 i=1,newtdm
         uold=urad(i,mpos+1)
         urad(i,mpos+1)=urad(i,mpos+1)+du(i)
         call nrv2(mpos,mdir,mend,potfn,cpotfn,ginfn,urad)
         uou(1)=urad(1,matr)
         uou(2)=urad(2,matr)
         call check(newtdm,uin,uou,udevn)
         do 120 j=1,newtdm
            dudev(j,i)=(udevn(j)-udev(j))/du(i)
120      continue
         urad(i,mpos+1)=uold
110   continue
      dum=dudev(1,1)*dudev(2,2)-dudev(1,2)*dudev(2,1)
      if(dum.eq.0.0d0)then
         do 130 i=1,newtdm
            du(i)=du(i)+dutiny
130      continue
         goto 100
      endif
      du(1)=(-udev(1)*dudev(2,2)+udev(2)*dudev(1,2))/dum
      du(2)=(-udev(2)*dudev(1,1)+udev(1)*dudev(2,1))/dum
      urad(1,mpos+1)=urad(1,mpos+1)+du(1)
      urad(2,mpos+1)=urad(2,mpos+1)+du(2)
*
* improved trial; since we have a linear problem, this should yield
* the solution. If not, we dare to try NEWTON one more time.
*
      call nrv2(mpos,mdir,mend,potfn,cpotfn,ginfn,urad)
      uou(1)=urad(1,matr)
      uou(2)=urad(2,matr)
      call check(newtdm,uin,uou,udev)
50    format(' 50 usolve2: Matching failure:  udev(1) = ',e16.8,/,t33,'u
     &dev(2) = ',e16.8)
      if((abs(udev(1)).gt.uacc).or.(abs(udev(2)).gt.uacc))then
         if(encore)then
            encore=.false.
            goto 100
         else

            write(*,50) udev(1),udev(2)
         endif
      endif

      return
      end


*********************************************************from usolve1***
      subroutine nrv1(mpos,mdir,mend,potfn,ginfn,urad)
*
*     Numerov algorithm to solve ODE like
*            urad()'' = poterm()*urad() + giterm()
*     The Numerov recursion is performed on the abscissa mesh stored
*     in common/ rmesh/. Changing stepsizes are allowed and if a
*     change occurs, the algoritm is reinitialized by extrapolation
*     {function call to uex1(), currently via 6th order polynomial}.
*     INPUT     mpos - starting index on mesh
*               mdir - =1/-1 for upward/ downward integration
*               mend - ending index on mesh
*               potfn() - coefficient function
*               ginfn() - inhomogeneity function
*     OUTPUT    urad - solution to ODE on mesh sites u(x) between
*                      'mpos' & 'mend'.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500)
      dimension poterm(3),giterm(3),utrue(3)
     &,         potfn(0:ntot),ginfn(0:ntot),urad(0:ntot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      parameter(stcrit=1.0d-10)
      external uex1

      poterm(1) = potfn(mpos)
      poterm(2) = potfn(mpos+mdir)
      giterm(1) = ginfn(mpos)
      giterm(2) = ginfn(mpos+mdir)
      utrue(1) = urad(mpos)
      utrue(2) = urad(mpos+mdir)

      mcon=mdir*mend
      step = abs(rad(mpos)-rad(mpos+mdir))
      stepsq = step**2/12.0d0
      ir = mpos+2*mdir
100   if(mdir*ir.le.mcon)then
         poterm(3)=potfn(ir)
         giterm(3)=ginfn(ir)
         utrue(3)=(2.0d0*(1.0d0-5.0d0*stepsq*poterm(2))*utrue(2)
     &              -(1.0d0+poterm(1)*stepsq)*utrue(1)
     &              +(giterm(3)+1.0d1*giterm(2)+giterm(1))*stepsq)
     &            /(1.0d0+poterm(3)*stepsq)
         urad(ir) = utrue(3)
*
* next step: if stepsize is altered, reinitialize
*
         stdif = abs(rad(ir)-rad(ir+mdir))
         if(abs(stdif-step).ge.stcrit)then
            utrue(1)=utrue(3)
            utrue(2)=uex1(mdir,ir,urad)
            ir=ir+mdir
            step=stdif
            stepsq=step**2/12.0d0
            urad(ir)=utrue(2)
            poterm(1)=poterm(3)
            poterm(2)=potfn(ir)
            giterm(1)=giterm(3)
            giterm(2)=ginfn(ir)
            urad(ir)=utrue(2)
         else
            do 200 j = 1, 2
               utrue(j) = utrue(j+1)
               poterm(j) = poterm(j+1)
               giterm(j) = giterm(j+1)
200         continue
         endif
      ir = ir+mdir
      goto 100
      endif

      return
      end


************************************************************from nrv1***
      double precision function uex1(nrdir,nr,urad)
*
*          uex1 = urad(rad(nr+nrdir))
*     Comment at subroutine nrv1. 'nrdir'=+1/-1 for upward/ downward
*     extrapolation {call to extrapolate}, here for order 'idim'-1
*     polynomial.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500)
      dimension urad(0:ntot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      parameter(idim=7)
      dimension r(idim),u(idim),c(idim),d(idim)
      external  extrapol

      rex = rad(nr+nrdir)
      do 100 i=1,idim
         r(i)=rad(nr-nrdir*(idim-i))
         u(i)=urad(nr-nrdir*(idim-i))
         c(i)=u(i)
         d(i)=u(i)
100   continue
      ic=idim
      uex1=extrapol(idim,ic,r,rex,u,c,d)

      return
      end


*********************************************************from usolve2***
      subroutine nrv2(mpos,mdir,mend,potfn,cpotfn,ginfn,urad)
*
*  Numerov algorithm to solve system of 2 coupled ODEs like
*               urad()'' = pott()*urad() + gitt() .
*  pott() is the 2x2 coefficient matrix, gitt() the inhomogeneity.
*  For the reference wave function problem pott() is symmetric.
*  The Numerov recursion is performed on the abscissa mesh stored
*  in common/rmesh/ . Changing stepsizes are allowed and if a
*  change occcurs, the algorithm is reinitialized by extrapolation
*  {function call to uex2(), currently via 6th order polynomial}.
*  INPUT  mpos - starting index on rmesh
*         mdir - +1/-1 for upward/downward integration
*         mend - ending index on rmesh
*         potfn() - diagonal coefficient function
*         cpot() - off-diagonal coefficient function
*                  {currently 1dimensional due to symmetry}
*         ginhfn() - inhomogenity function
*  OUTPUT urad() - solution of the system
************************************************************************
      implicit  double precision (a-h,o-z)
      parameter (ntot=500)
      parameter (stcrit=1.0d-10)
      dimension potfn(2,0:ntot),cpotfn(0:ntot),ginfn(2,0:ntot)
     &,         urad(2,0:ntot),pott(2,2,3),gitt(2,3),utru(2,3)
     &,         potinv(2,2),potmid(2,2)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
*
* initialize
*
      step=abs(rad(mpos)-rad(mpos+mdir))
      stsq=step**2/12.0d0
      do 100 k=1,2
         m=mpos+(k-1)*mdir
         pott(1,1,k)=1.0d0+stsq*potfn(1,m)
         pott(1,2,k)=stsq*cpotfn(m)
         pott(2,1,k)=stsq*cpotfn(m)
         pott(2,2,k)=1.0d0+stsq*potfn(2,m)
         gitt(1,k)=ginfn(1,m)
         gitt(2,k)=ginfn(2,m)
         utru(1,k)=urad(1,m)
         utru(2,k)=urad(2,m)
100   continue

      mcon=mdir*mend
      m=mpos+2*mdir
200   if(mdir*m.le.mcon)then
         pott(1,1,3)=1.0d0+stsq*potfn(1,m)
         pott(2,2,3)=1.0d0+stsq*potfn(2,m)
         pott(1,2,3)=stsq*cpotfn(m)
         pott(2,1,3)=stsq*cpotfn(m)
         det=pott(1,1,3)*pott(2,2,3)-pott(1,2,3)*pott(2,1,3)
         potinv(1,1)=pott(2,2,3)/det
         potinv(2,2)=pott(1,1,3)/det
         potinv(1,2)=-pott(1,2,3)/det
         potinv(2,1)=-pott(2,1,3)/det
         potmid(1,1)=2.0d0*(-5.0d0*pott(1,1,2)+6.0d0)
         potmid(2,2)=2.0d0*(-5.0d0*pott(2,2,2)+6.0d0)
         potmid(1,2)=2.0d0*(-5.0d0*pott(1,2,2))
         potmid(2,1)=2.0d0*(-5.0d0*pott(2,1,2))
         gitt(1,3)=ginfn(1,m)
         gitt(2,3)=ginfn(2,m)
         do 210 j=1,2
            dum=0.0d0
            do 220 l=1,2
               dum=potmid(j,l)*utru(l,2)-pott(j,l,1)*utru(l,1)+dum
220         continue
            dum=stsq*(gitt(j,3)+1.0d1*gitt(j,2)+gitt(j,1))+dum
            utru(j,3)=dum
210      continue
         do 230 j=1,2
            dum=0.0d0
            do 240 l=1,2
               dum=potinv(j,l)*utru(l,3)+dum
240         continue
            utru(j,3)=dum
            urad(j,m)=dum
230      continue
*
* next step: if stepsize is altered, reinitialize
*
         stdif=abs(rad(m)-rad(m+mdir))
         if(abs(stdif-step).ge.stcrit)then
            utru(1,1)=utru(1,3)
            utru(2,1)=utru(2,3)
            utru(1,2)=uex2(1,mdir,m,urad)
            utru(2,2)=uex2(2,mdir,m,urad)
            m=m+mdir
            step=stdif
            stsq=stdif**2/12.0d0
            urad(1,m)=utru(1,2)
            urad(2,m)=utru(2,2)
            pott(1,1,1)=pott(1,1,3)
            pott(2,2,1)=pott(2,2,3)
            pott(1,2,1)=pott(1,2,3)
            pott(2,1,1)=pott(2,1,3)
            pott(1,1,2)=1.0d0+stsq*potfn(1,m)
            pott(2,2,2)=1.0d0+stsq*potfn(2,m)
            pott(1,2,2)=stsq*cpotfn(m)
            pott(2,1,2)=stsq*cpotfn(m)
            gitt(1,1)=gitt(1,3)
            gitt(2,1)=gitt(2,3)
            gitt(1,2)=ginfn(1,m)
            gitt(2,2)=ginfn(2,m)
         else
            do 250 k=1,2
               utru(1,k)=utru(1,k+1)
               utru(2,k)=utru(2,k+1)
               gitt(1,k)=gitt(1,k+1)
               gitt(2,k)=gitt(2,k+1)
               pott(1,1,k)=pott(1,1,k+1)
               pott(2,2,k)=pott(2,2,k+1)
               pott(1,2,k)=pott(1,2,k+1)
               pott(2,1,k)=pott(2,1,k+1)
250         continue
         endif
      m=m+mdir
      goto 200
      endif

      return
      end


************************************************************from nrv2***
      double precision function uex2(nc,nrdir,nr,urad)
*
*               uex2 = urad(nc,rad(nr+nrdir)
*     Comment at suroutine nrv2. 'nrdir'=+1/-1 for upward/downward
*     extrapolation {call to extrapolate}, here for order 'idim'-1
*     polynomial.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (ntot=500)
      dimension urad(2,0:ntot)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      parameter (idim=7)
      dimension r(idim),u(idim),c(idim),d(idim)
      external extrapol

      rex=rad(nr+nrdir)
      do 100 i=1,idim
         r(i)=rad(nr-nrdir*(idim-i))
         u(i)=urad(nc,nr-nrdir*(idim-i))
         c(i)=u(i)
         d(i)=u(i)
100   continue
      ic=idim
      uex2=extrapol(idim,ic,r,rex,u,c,d)

      return
      end


*****************************************************from usolve1 & 2***
      subroutine check(ndim,xx,yy,zz)
*
*     Pertains to subroutines usolve1&2. Output is zz()=xx()-yy().
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xx(ndim),yy(ndim),zz(ndim)

      do 100 i=1,ndim
         zz(i)=xx(i)-yy(i)
100   continue
      return
      end


********************************************************from uex1 & 2***
      double precision function extrapol(idim,ic,r,rex,u,c,d)
*
*     extrapol(rex) is the value of the function u(r=rex), extra-
*     polated from i=1...idim values u(r(i)) as a polynomial of order
*     idim-1.
*     INPUT idim - order of extrapolating polynomial + 1
*           ic - r(ic) is closest to rex
*           r() - abscissa values
*           rex - abscissa value to extrapolate to
*           u() - functional values u(r())
*           c(),d() = u(), altered in calculation
************************************************************************
      implicit  double precision (a-h,o-z)
      dimension r(idim),u(idim),c(idim),d(idim)

      uex=u(ic)
      ic=ic-1
      do 100 k=1,idim-1
         do 200 i=1,idim-k
            rad1=r(i)-rex
            rad2=r(i+k)-rex
            w=c(i+1)-d(i)
            drad=rad1-rad2
            if(drad.eq.0.0d0) stop 'extrapol: Extrapolation failed.'
            drad=w/drad
            d(i)=rad2*drad
            c(i)=rad1*drad
200      continue
         if(2*ic.lt.idim-k)then
            du=c(ic+1)
         else
            du=d(ic)
            ic=ic-1
         endif
         uex=uex+du
100   continue
      extrapol=uex

      return
      end


************************************************************************
      subroutine trapez(min,max,mlo,mup,x,func,fintgrl)
*
*     Trapezoidal Rule for function func(x) over abscissa array x().
*     Output is fintgrl=Integral of func(x) between x(mlo) and x(mup).
*     The stepsize in x may change, mlo < mup.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (delta=1.0d-10)
      dimension x(min:max),func(min:max)

      fintgrl=0.0d0
      sum=.5d0*func(mlo)
      step=x(mlo+1)-x(mlo)

      i=mlo+1
100   if(i.le.mup-1)then
         if(abs(x(i+1)-x(i)-step).gt.delta)then
            sum=sum+.5d0*func(i)
            fintgrl=fintgrl+sum*step
            step=x(i+1)-x(i)
            sum=.5d0*func(i)
         else
            sum=sum+func(i)
         endif
      i=i+1
      goto 100
      endif
      fintgrl=fintgrl+(sum+.5d0*func(mup))*step

      return
      end


************************************************************************
      double precision function bessj(l,x)
*
*     Calculates bessj(l,x) = x*j(l,x), where j(l,x) is the Spherical
*     Bessel Function of integer order l.
*     l=0,2 : The analytic expressions in sin(x) & cos(x) are used.
*     l>2   : The functional value is obtained by downward recursion,
*             if x<l, by upward recursion if x>l, calling 'besrec()'.
*     x<xmin,l>0: bessj()=x**(l+1)/(2*l-1)!!
*     If the argument is excessively small, 'bestiny' sets the output
*     equal to the lowest order term in the powerseries expansion.
************************************************************************
      implicit  double precision (a-h,o-z)
      parameter (xmin=1.0d-10)
      external  besrec,bestiny
      intrinsic sin,cos

      if(l.eq.0)then
         bessj=sin(x)
      else if(x.lt.xmin)then
         bessj=x*bestiny(l,x)
      else if(l.eq.1)then
         bessj=sin(x)/x-cos(x)
      else if(l.eq.2)then
         bessj=(3.0d0/x/x-1.0d0)*sin(x)-3.0d0*cos(x)/x
      else
         bessj=x*besrec(l,x)
      endif

      return
      end


************************************************************************
      double precision function besrec(l,x)
*
*     Recursion for Bessel functions. 'acclog' is about the negative
*     log of the desired accuracy.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(acclog=1.0d1)
      intrinsic sin,cos,dble

      if(x.gt.dble(l))then
*
* upward recursion
*
         downj=sin(x)/x
         tempj=(downj-cos(x))/x
         do 100 m=2,l
            upj=(2.0*m-1.0)/x*tempj-downj
              downj=tempj
              tempj=upj
100      continue
         besrec=upj
      else
*
* downward recursion & renormalization
*
         mstart=2*l+int(acclog*sqrt(dble(l)))
200      upj=0.0d0
         tempj=1.0d0
         do 300 m=mstart-1,0,-1
            downj=(2.0*m+3.0)/x*tempj-upj
            if(m.eq.l)besrec=downj
            upj=tempj
            tempj=downj
300      continue
         besrec=besrec/downj*sin(x)/x
      endif

400   return
      end


************************************************************************
      double precision function bestiny(l,x)
*
*     Lowest order powerseries term for spherical Bessel function.
************************************************************************
      implicit double precision (a-h,o-z)

      kdiv=1.0d0
      do 100 k=1,2*l-1,2
         kdiv=kdiv*k
100   continue
      bestiny=x**l/kdiv

      return
      end


************************************************************************
      double precision function bessh(l,x)
*
*     Calculates bessh(l,x)=2/pi*x*k(l,x), where k() is a Modified
*     Spherical Bessel function of integer order l.
*     l=0,1 : The analytic expressions are used.
*     l>1   : The functional value is obtained by upward recursion.
************************************************************************
      implicit  double precision (a-h,o-z)
      intrinsic exp

      if(l.eq.0)then
         bessh=1.0d0/exp(x)
      else if(l.eq.1)then
         bessh=(1.0d0+1.0d0/x)/exp(x)
      else
         downh=1.0d0/exp(x)
         temph=downh*(1.0d0+1.0d0/x)
         do 100 m=2,l
            uph=(2.0d0*m-1.0d0)/x*temph+downh
            downh=temph
            temph=uph
100      continue
         bessh=uph
      endif

      return
      end


************************************************************************
      double precision function dirpot1(j,l,ns,xx)
*
*     Local potentials for single channels.
*     INPUT   xx - radial coordinate
*             j,l,ns - quantum numbers of entrance channel
*     OUTPUT  dirpot2 - actual potential in fermi**-2
*     They are the REID SOFT CORE POTENTIALS.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (h=1.0463d1)
      parameter (scale=.241129d-1)
      parameter (xmu=.7d0)

      x=xmu*xx
      ex1=exp(x)
      ex2=ex1**2
      ex4=ex2*ex2
*
* cascade for jls combinations
*
      if(ns.eq.0)then
         if(l.eq.0)then
            rscp=-h/ex1-(1.6506d3-6.4842d3/ex2/ex1)/ex4
         elseif(l.eq.1)then
            rscp=3.0d0*h/ex1-(6.3439d2-2.1634d3/ex1)/ex2
         else if(l.eq.2)then
            rscp=-h/ex1-1.2322d1/ex2-(1.1126d3-6.4642d3/ex2/ex1)/ex4
         else
            goto 100
         endif
      else if(ns.eq.1)then
         if((l.eq.1).and.(j.eq.0))then
            rscp=-h*((1.0d0+4.0d0/x+4.0d0/x**2)/ex1 -(16.0d0/x
     &             +4.0d0/x**2)/ex4)+2.7133d1/ex2-7.9074d2/ex4
     &           +2.0662d4/ex4/ex2/ex1
         else if((l.eq.1).and.(j.eq.1))then
            rscp=+h*((1.0d0+2.0d0/x+2.0d0/x**2)/ex1-(8.0d0/x
     &             +2.0d0/x**2)/ex4)-(1.3525d2-4.7281d2/ex1)/ex2
         else if((l.eq.2).and.(j.eq.2))then
            rscp=-3.0d0*h*((1.0d0+2.0d0/x+2.0d0/x**2)/ex1 -(8.0d0/x
     &             +2.0d0/x**2)/ex4) -(2.2012d2-8.71d2/ex1)/ex2
         else
            goto 100
         endif
      else
         goto 100
      endif
      dirpot1=scale*rscp/x

      return
100   write(*,*) 'No potential for jls ', j,l,ns
      stop
      end


************************************************************************
      double precision function dirpot2(j,l,x)
*
*     Local potentials for direct terms in coupled cases.
*     INPUT   x - radial coordinate
*             j, l - quantum numbers of entrance channel
*     OUTPUT  dirpot2 - actual potential in fermi**-2
*     REMARK  The various terms are included WITH the eigenvalues
*        of the corresponding operators, whose functional
*        form is defined in the subsequent functions
*            tenpot    tensor interaction (Note REMARK there)
*               olspot    spin orbit interaction
*              cenpot    central interaction
*             They are the REID SOFT CORE POTENTIALS, Õ12þ.
************************************************************************
      implicit double precision (a-h,o-z)
      external  cenpot, tenpot, olspot

      if(j.eq.1)then
         tens=sqrt(8.0d0)
         if(l.eq.0)then
           rscp=-cenpot(1,x)
         else
           rscp=-cenpot(1,x)+2.0d0/tens*tenpot(1,x)+3.0d0*olspot(1,x)
         endif
      else if(j.eq.2)then
         tens=1.2d0*sqrt(6.0d0)
         if(l.eq.1)then
            rscp=-cenpot(2,x)+.4d0/tens*tenpot(2,x)-olspot(2,x)
         else
            rscp=-cenpot(2,x)+1.6d0/tens*tenpot(2,x)+4.0d0*olspot(2,x)
         endif
      else
         stop 'dirpot2: j out of range.'
      endif
      dirpot2=-rscp
      return
      end

      double precision function tenpot(j,xx)
**
**    REMARK   The varible 'tens' denotes the eigenvalue for the
**             off-diagonal matrixelement <j j-1 |TR|j j+1>.
**
      implicit double precision (a-h,o-z)
      parameter(h=1.0463d1)
      parameter(xmu=.7d0)
      parameter(scale=2.411388d-2)

      x=xmu*xx
      ex2=exp(x)**2
      ex4=ex2*ex2
      if(j.eq.1)then
         tens=sqrt(8.0d0)
         tenpot=-h*((1.0d0+3.0d0/x+3.0d0/x**2)/exp(x)
     &          -(12.0d0/x**1+3.0d0/x**2)/ex4)
     &          +3.5177d2/ex4-1.6735d3/ex4/ex2
      else if(j.eq.2)then
         tens=1.2d0*sqrt(6.0d0)
         tenpot=h*((1.0d0/3.0d0 +1.0d0/x +1.0d0/x**2)/exp(x)
     &          -(4.0d0+1.0d0/x)/x/ex4)
     &          -3.4925d1/exp(x)/ex2
      else
         stop 'dirpot2: Error in tenpot.'
      endif
         tenpot=tens*scale*tenpot/x
      return
      end


      double precision function olspot(j, xx)

      implicit double precision (a-h,o-z)
      parameter(xmu=.7d0)
      parameter(scale=2.411388d-2)

      x=xmu*xx
      ex2=exp(x)**2
      ex4=ex2*ex2
      if(j.eq.1)then
         olspot=(7.0891d2-2.7131d3/ex2)/ex4
      else if(j.eq.2)then
         olspot=-2.0741d3/ex2/ex4
      else
         stop 'dirpot2: Error in olspot.'
      endif
         olspot=scale*olspot/x
      return
      end

      double precision function cenpot(j,xx)

      implicit double precision (a-h,o-z)
      parameter(h=1.0463d1)
      parameter(xmu=.7d0)
      parameter(scale=2.411388d-2)

      x=xmu*xx
      ex2=exp(x)**2
      ex4=ex2*ex2
      if(j.eq.1)then
         cenpot=-h/exp(x)+1.05468d2/ex2-3.1878d3/ex4+9.9243d3/ex4/ex2
      else if(j.eq.2)then
         cenpot=h/3.0d0/exp(x)-9.3348d2/ex4+4.1521d3/ex4/ex2
      else
         stop 'dirpot2: Error in cenpot.'
      endif
         cenpot=scale*cenpot/x
      return
      end


************************************************************from MAIN***
      subroutine show(efm,gmasq,nchlow,nchup,chanid)
*
*  Terminal output on grids and parameters.
************************************************************************
      implicit double precision (a-h,o-z)
      parameter(ntot=500, ktot=50, nchmx=12)
      parameter(delta=1.d-10)
      character *6 chanid(nchmx)
      common / rmesh/ minr,maxr,matr,mhco,rad(0:ntot)
      common / kmesh/ mink,maxk,pk(ktot)
      data     rold,pold/0.,0./

48    format(//' Healing Parameter GMASQ     : ', f8.4,' 1/fmª2')
49    format(' Effective Mass Parameter EFM: ', f8.4)
50    format(//' FILENAMES for Angular Momentum Channels:'//
     &        ' Internal Number',2x,'Filename',4x,'j  l  s',
     &        4x,'(C = Crossterm for Coupled Channels)'/)
51    format(9x,i3,6x,a10,2x,3(a1,2x))
52    format(//' POSITION SPACE MESH - At Changing Stepsizes, in Õfmþ.'/
     &/'  i',5x,'r(i)',6x,'old step',2x,'i+1',3x,'r(i+1)',2x,'new step'/
     &)
53    format(/' Physical Size: ',i4,'   Used points    : ', i4)
54    format(' Minimum      : ',i4,'   Maximum        : ', i4)
55    format(' Matching at  : ',i4,'   Matching Radius: ', f8.4,' fm')
56    format(' Hard Core at : ',i4,'   Core Radius    : ',f8.4,' fm')
57    format(1x,2(i4,2x,f8.4,2x,f8.4,2x))
58    format(1x,i4,2x,2(f8.4,2x),'----',2x,'--------',2x,'--------')
62    format(//' MOMENTUM SPACE MESH - At Changing Stepsizes, in Õ1/fmþ.
     &'//'  i',5x,'p(i)',6x,'old step',2x,'i+1',3x,'p(i+1)',2x,'new step
     &'/)

      write(*,52)
      do 200 i=0,maxr
         rnew=rad(i+1)-rad(i)
         if(abs(rnew-rold).gt.delta)then
            if(i.eq.maxr)then
               write(*,58) i,rad(i),rold
            else
               write(*,57) i,rad(i),rold,i+1,rad(i+1),rnew
            endif
         endif
         rold=rnew
200   continue
      write(*,53)ntot+1,maxr+1
      write(*,54)minr,maxr
      write(*,55)matr,rad(matr)
      write(*,56)mhco,rad(mhco)

      write(*,62)
      do 300 i=mink,maxk
         pnew=pk(i+1)-pk(i)
         if(abs(pnew-pold).gt.delta)then
            if(i.eq.maxk)then
               write(*,58) i,pk(i),pold
            else
               write(*,57) i,pk(i),pold,i+1,pk(i+1),pnew
            endif
         endif
         pold=pnew
300   continue
      write(*,53)ktot,maxk
      write(*,54)mink,maxk

      write(*,48) gmasq
      write(*,49) efm
      write(*,50)
      do 100 i=nchlow, nchup
         write(*,51)i,chanid(i),chanid(i)(4:4),chanid(i)(5:5)
     &,               chanid(i)(6:6)
100   continue
      write(*,53) nchmx,nchup
      write(*,'(//)')

      return
      end


      subroutine greetings

50    format(//,t20,'***** T R E F *****',//,' Nuclear Matter Effective
     &Interaction: The Reference Reaction Matrix with Reid''s SHCP.',//,
     & t10,'Martin Fuchs  &  Philip J. Siemens',/,t10,'Department of Phy
     &sics',/,t10,'Oregon State University',/,t10,'Corvallis, OR 97331',
     &/,t10,'U.S.A',//,t10,'June 15, 1991',////)
      write(*,50)
      pause
      return
      end
