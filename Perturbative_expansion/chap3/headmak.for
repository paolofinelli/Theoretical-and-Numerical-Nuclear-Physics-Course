************************************************************************
      program headmak
*
*  Generates sample input file for program TNM.
*
*         VARIABLES
*  1    memo       memnonic (character*64)
*         REFERENCE SPECTRUM:
*  2    efm        effective mass parameter
*  3    gmasq      squared healing parameter Õ1/fmª2þ
*         NUCLEAR MATTER:
*  4    pfm        Fermi momentum Õ1/fmþ
*  5    pav        two particle average momentum Õ1/fmþ
*  6    pkmx       subspace cut-off Õ1/fmþ
*  7    emev       starting energy ÕMeVþ
*  8    etab()     nuclear matter single particle spectrum Õ1/fmþ,
*                  the function 'spectrum' gives the functional form
*  9    pk()       absissa values for spectrum
*  10   kup        number of points on which the spectrum is given
*         SPECIFICS OF "TNM"
*  11   mup        mup-1 = highest order of the Legindgaard polynomials
*  12   inid       prefix of input files  (character*3)
*  13   outid      prefix of output files (   - " -   )
*
************************************************************************
      implicit real*4 (a-h,o-z)
      character*5 inid,outid,memo*64
      parameter(kmx=50,scale=41.47)
      dimension pk(kmx),etab(kmx)
      data  memo/ '''Sample for TNM control file HEAD.'''/
      data  efm,gmasq,pfm,pav,pkmx,emev,mup,kup,inid,outid/
     &      1.0,1.4,1.3,0.7,2.8,162.3,5,20,'''TRF''','''SAM'''/

      spectrum(x)=scale*.5*x**2+100.0
*
* output formats: 50  --  line 1                memnonic
*                 51  --  lines 1-6 & 10-13     parameters
*                 52  --  lines 9,8             spectrum
50    format(a64,1x)
51    format(5(f8.4,1x),f10.4,2(1x,i3),2(1x,a5))
52    format(f8.4,1x,e16.8)

      pstep=pkmx/(kup-1)
      pk(1)=0.
      etab(1)=spectrum(0.)
      do 100 i=2,kup
         p=(i-1)*pstep
         pk(i)=p
         etab(i)=spectrum(p)
100   continue

      open(11,file='HEAD',form='formatted')
      write(11,50)memo
      write(11,51)efm,gmasq,pfm,pav,pkmx,emev,mup,kup,inid,outid
      write(11,52)(pk(i),etab(i),i=1,kup)
      endfile(11)
      close(11)

      end
