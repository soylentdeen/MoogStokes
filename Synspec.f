
      subroutine synspec (phi, psi)
c******************************************************************************
c     This routine does synthetic spectra                                
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      real*8 dd(5000)
c      real*8 Stokes(5)
c      real*8 NOK, NBAD, H_MIN, H_GUESS, EPS, B
c      real*8 factor, right
c      PARAMETER (NDGL=5,NRD=0)
c      PARAMETER (LWORK=11*NDGL+8*NRD+21, LIWORK=NRD+21)
c      DIMENSION STOKES(NDGL),WORK(LWORK),IWORK(LIWORK)
c      DIMENSION WORK(LWORK),IWORK(LIWORK)
c      real*8 STOKES(NDGL), rpar, atol, rtol
c      real*8 work(76), iwork(21)
c      integer ipar(1), lwork, liwork
c      real*8 rpar(1)
c      integer ipar
      real*8 phi, psi
      logical direction, prev_step
      real*8 Stokes(4)
c      EXTERNAL derivs, SOLOUT

c      NEQS = 5
c      RPAR=0.0
c      ipar = 0
c      counter = 0

c*****initialize the synthesis
      direction = .TRUE.
      prev_step = .FALSE.
      write (nf1out,1101)
      write (nf2out,1002) moditle(1:73)
      if (iunits .eq. 1) then
         write (nf2out,1103) oldstart,oldstop,oldstep,olddelta
      else
         write (nf2out,1102) start,sstop,step,delta
      endif
      if (iraf .eq. 1) then
         npoint = (sstop-start)/step
         write (nf4out,1104) npoint,wave,wave,step,step
         write (nf4out,1105)
         write (nf4out,1106) moditle
         do j=1,93
            if (pec(j) .gt. 0 ) then
               dummy1(j) = dlog10(xabund(j)) + 12.0
               write (nf4out,1107) names(j),dummy1(j)
            endif
         enddo
         write (nf4out,1108) vturb(1)
         write (nf4out,1109)
      endif
      n = 1           
      num = 0
      nsteps = 1
      lim1line = 0


c*****calculate continuum quantities at the spectrum wavelength
      wave = start
      left = wave
      right = wave
      wavl = 0.
      nf11out=11
c      nf12out=12
c      nf13out=13
      open(unit=nf11out, file=f11out)
c      open(unit=nf12out, file=f12out)
c      open(unit=nf13out, file=f13out)
30    if (dabs(wave-wavl)/wave .ge. 0.001) then
         wavl = wave   
         call opacit (2,wave)    
      endif


c*****find the appropriate set of lines for this wavelength, reading 
c     in a new set if needed
      if (mode .eq. 3) then
20       call linlimit
         if (lim2line .lt. 0) then
            call inlines (2)
            call nearly (1)
            go to 20
         endif
         lim1 = lim1line
         lim2 = lim2line
      endif


c*****compute a spectrum depth at this point
      call calcopacities
      call taukap(phi, psi, Stokes)
c      d(n) = Stokes(1)
      write (nf11out,12345) wave,Stokes
      write (*,*) wave,Stokes
c      write (nf12out,*) ' ' 
c      write (nf13out,*) ' '
c      write (*,*) wave, 1.0-d(n),idid
      if (mod(n,10) .eq. 0) then
         if (iraf .eq. 1) then
            do j=1,10
               dd(num+j) = 1. - d(num+j)
            enddo
            write (nf4out,1110) (dd(num+j),j=1,10)
         endif
         if (iunits .eq. 1) then
            wave3 = 1.d-4*(wave - 9.0*step)
            write (nf1out,1112) wave3,(d(num+j),j=1,10)
         else
            wave3 = wave - 9.0*step
            write (nf1out,1111) wave3,(d(num+j),j=1,10)
         endif
         if (nf2out .gt. 0) write (nf2out,1110) (d(num+j),j=1,10)
         num = num + 10
      endif


c*****step in wavelength and try again 
      if (d(n).gt.0.05) then
          stepsize = dopp(nstrong,50)*wave/2.997929e10/10.0
c            First step into a region with a line.  Need to reverse direction
          if (.not.prev_step) THEN
              direction = .FALSE.
              prev_step = .TRUE.
              right = wave
          ENDIF
          if (direction) THEN
              factor = 1.0
          else
              factor = -1.0
          endif
c          write (*,*) "line!  - ", d(n), wave, wave*dopp(nstrong,20)
c     .           /2.997929e10
          wave = wave + factor*stepsize
      else
          stepsize = dopp(nstrong, 50)*wave/2.997929e11
          prev_step = .FALSE.
          if (.not.direction) THEN
              direction = .TRUE.
              prev_step = .TRUE.
              stepsize = stepsize/50.0
              wave = right
          ENDIF
c          write (*,*) "continuum! - ", d(n), wave, wave*dopp(nstrong,20)
c     .           /2.997929e10
          wave = wave + stepsize
      endif
c      wave = oldstart + step*nsteps
      if (wave .le. sstop) then
         n = n + 1        
         nsteps = nsteps + 1
         if (n .gt. 5000) then
            n = 1                                      
            num = 0
         endif
         go to 30                   


c*****finish the synthesis
      else
         nn = mod(n,10)
         if (nn .ne. 0) then
            if (iraf .eq. 1) then
               do j=1,nn
                  dd(num+j) = 1. - d(num+j)
               enddo
               write (nf4out,1110) (dd(num+j),j=1,nn)
            endif
            if (iunits .eq. 1) then
               wave3 = 1.d-4*(wave - 9.0*step)
               write (nf1out,1112) wave3,(d(num+j),j=1,nn)
            else
               wave3 = wave - 9.0*step
               write (nf1out,1111) wave3,(d(num+j),j=1,nn)
            endif
            if (nf2out .gt. 0) write (nf2out,1110) (d(num+j),j=1,nn)
         endif
         if (iunits .eq. 1) then
            write (nf1out,1113) 1.d-4*wave
         else
            write (nf1out,1114) wave
         endif
         close(nf11out)
c         close(nf12out)
c         close(nf13out)
         return 
      endif


c*****format statements
1001  format ('  kaplam from 1 to ntau at wavelength',f10.2/
     .        (6(1pd12.4)))
1002  format ('MODEL: ',a73)
1003  format ('AT WAVELENGTH/FREQUENCY =',f11.7,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1004  format ('AT WAVELENGTH/FREQUENCY =',f11.3,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1101  format (/'SPECTRUM DEPTHS')
1102  format (4f11.3)
1103  format (4f10.7)
1104  format ('SIMPLE  =    t'/'NAXIS   =     1'/'NAXIS1  = ',i10,/
     .        'W0      =',f10.4/'CRVAL1  =',f10.4/'WPC     =',f10.4/
     .        'CDELT1  =',f10.4)
1105  format (16HORIGIN  = 'moog'/21HDATA-TYP= 'synthetic'/
     .        18HCTYPE1  = 'lambda'/21HCUNIT1  = 'angstroms')
1106  format (11HTITLE   = ',A65,1H')
1107  format ('ATOM    = ',1H',7x,a2,1H',/,'ABUND   = ',f10.2)
1108  format ('VTURB   = ',d10.4,'     /  cm/sec  ')
1109  format ('END')
1110  format (10f7.4)
1111  format (f10.3,': depths=',10f6.3)
1112  format (f10.7,': depths=',10f6.3)
1113  format ('FINAL WAVELENGTH/FREQUENCY =',f10.7/)
1114  format ('FINAL WAVELENGTH/FREQUENCY =',f10.3/)
12345 format (f10.4,5e15.5)

      end                                


c      subroutine Planck(temperature, B)
c      
c      implicit real*8 (a-h,o-z)
c      include 'Atmos.com'
c      include 'Linex.com'
c      include 'Factor.com'
c      include 'Pstuff.com'
c      include 'Dummy.com'
c      B = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
c     .      (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
c      end
c
c      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
c     &                     RPAR,IPAR,IRTRN,XOUT)
c      DIMENSION CON(8*ND),ICOMP(ND)
c      real*8 X, EI, EQ, EV, ZQ, ZV
c      real*8 Y(N)
c      include 'Atmos.com'
c      include 'Linex.com'
c      CALL LINTERPOLATE(ETA_I, 10.0**X, EI)
c      CALL LINTERPOLATE(ETA_Q, 10.0**X, EQ)
c      CALL LINTERPOLATE(ETA_V, 10.0**X, EV)
c      CALL LINTERPOLATE(ZET_Q, 10.0**X, ZQ)
c      CALL LINTERPOLATE(ZET_V, 10.0**X, ZV)
c      write (nf12out,322) X, Y
c      write (nf13out,322) X, EI, EQ, EV, ZQ, ZV
c322   format (e11.5, 5e13.5)
c      end 


c      subroutine dump_taus(value)
c
c      implicit real*8 (a-h,o-z)
c      include 'Atmos.com'
c      include 'Linex.com'
c
c      do 56 i=1, ntau
c56      write (nf3out,321) tauref(i), eta_I(i),eta_Q(i),eta_V(i),
c     .        zet_Q(i),zet_V(i)
c      close(nf3out)
c
c321   format (f11.3, 5e11.3)
c      end
