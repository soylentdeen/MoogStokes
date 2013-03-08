
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
      real*8 strongratio, weakratio


c*****examine the parameter file

c      array = 'GIVE THE MINIMUM LINE/CONTINUUM OPACITY RATIO TO KEEP: '
      strongratio = 0.1
      weakratio = 0.001
      ns_lines = 0
      nw_lines = 0
      

c*****compute the line opacities
      write (*,*) "Checkpoint A"
      call inlines (1)
      write (*,*) "Checkpoint B"
      call eqlib
1     call nearly (1)


c*****calculate continuum quantities at the line list wavelength middle
      wave = (wave1(1)+wave1(nlines))/2.
      call opacit (2,wave)
      
c*****divide the lines into keepers and discards
      do j=1,nlines+nstrong
         if (strength(j)/kaplam(jtau5) .ge. strongratio) then
             ns_lines = ns_lines+1
             strong(ns_lines) = wave1(j)
             write (*,*) "Strong line at ", wave1(j)
         elseif (strength(j)/kaplam(jtau5) .ge. weakratio) then
             nw_lines = nw_lines+1
             weak(nw_lines) = wave1(j)
             write (*,*) "Weak line at ", wave1(j)
         endif
      enddo
      if (nlines +nstrong .eq. 2500) then
         write (*,*) "Checkpoint C"
         call inlines (6)
         write (*,*) "Checkpoint D"
         go to 1
      endif

      end

c      subroutine sortwave
c
c      implicit real*8 (a-h,o-z)
c      include 'Atmos.com'
c      include 'Linex.com'
c      include 'Pstuff.com'
c      INTEGER I,J
c      REAL*8 X, temp(nwave)
c      DO 30 I=2,nwave
c      X=wavelength(I)
c      J=I
c   10 J=J-1
c      IF(J.EQ.0 .OR. wavelength(J).LT.X) GO TO 20
c      wavelength(J+1)=wavelength(J)
c      GO TO 10
c   20 wavelength(J+1)=X
c   30 CONTINUE
c
cc****   Now, go through and purge duplicates and very close together points
c      temp(1) = wavelength(1)
c      I=1
c      do J=2,nwave
c         if(abs(wavelength(J)-temp(I)) .lt. 0.001) THEN
c            temp(I) = (wavelength(J)+temp(I))/2.0
c            wavelength(J) = 0.0
c         else
c            I = I+1
c            temp(I) = wavelength(J)
c            wavelength(J) = 0.0
c         endif
c      enddo
c      nwave = I
c      do J=1, nwave
c         wavelength(J)=temp(J)
c      enddo
c      END
