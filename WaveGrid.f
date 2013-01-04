
      subroutine wavegrid
c******************************************************************************
c     This routine goes through an initial line list and culls from it
c     absorption lines that are not substantial contributors.  This is
c     done in a simple fashion by eliminating lines weaker than X, where
c     X = kapnu/kaplam at the approximate line wavelength, calculated
c     at a continuue optical depth of 0.5.  The user will be prompted for
c     the desired value of X.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      real*8 xratio


c*****examine the parameter file

c      array = 'GIVE THE MINIMUM LINE/CONTINUUM OPACITY RATIO TO KEEP: '
      xratio = 0.0001
      

c*****compute the line opacities
      call inlines (1)
      call eqlib
1     call nearly (1)


c*****calculate continuum quantities at the line list wavelength middle
      wave = (wave1(1)+wave1(nlines))/2.
      call opacit (2,wave)
      
c*****   Generate original wavelength grid

      deltawave = 0.25
      nwave = 1

      wave = start
30    wavelength(nwave) = wave
      wave = wave + deltawave
      nwave = nwave + 1
      if (wave .le. sstop) then
          goto 30
      endif
      wave = start

c*****divide the lines into keepers and discards
      do j=1,nlines+nstrong
         if (strength(j)/kaplam(jtau5) .ge. xratio) then
             do k=-20,20
                 wavelength(nwave) = wave1(j)+0.25*asin(real(k/20.0))
                 nwave = nwave+1
             enddo
             wavelength(nwave) = wave1(j)+ 0.001
             nwave = nwave+1
             wavelength(nwave) = wave1(j)- 0.001
             nwave = nwave+1
         else
c             write (*,*) "Line at :", wave1(j), " is not strong enough!"
c             write (*,*) "Strength :", strength(j)/kaplam(jtau5)
         endif
      enddo
      if (nlines +nstrong .eq. 2500) then
         call inlines (6)
         go to 1
      endif
      nwave = nwave -1

      call sortwave

      end

      subroutine sortwave

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      INTEGER I,J
      REAL*8 X, temp(nwave)
      DO 30 I=2,nwave
      X=wavelength(I)
      J=I
   10 J=J-1
      IF(J.EQ.0 .OR. wavelength(J).LT.X) GO TO 20
      wavelength(J+1)=wavelength(J)
      GO TO 10
   20 wavelength(J+1)=X
   30 CONTINUE

c****   Now, go through and purge duplicates and very close together points
      temp(1) = wavelength(1)
      I=1
      do J=2,nwave
         if(abs(wavelength(J)-temp(I)) .lt. 0.001) THEN
            temp(I) = (wavelength(J)+temp(I))/2.0
            wavelength(J) = 0.0
         else
            I = I+1
            temp(I) = wavelength(J)
            wavelength(J) = 0.0
         endif
      enddo
      nwave = I
      do J=1, nwave
         wavelength(J)=temp(J)
      enddo
      END
