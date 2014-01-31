
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
      include 'Stokes.com'
      real*8 strongratio, weakratio


c*****examine the parameter file

c      array = 'GIVE THE MINIMUM LINE/CONTINUUM OPACITY RATIO TO KEEP: '
      strongratio = 0.1
      weakratio = 0.01
      ns_lines = 0
      nw_lines = 0
      

c*****compute the line opacities
      call inlines (1)
      call eqlib
1     call nearly (1)


c*****calculate continuum quantities at the line list wavelength middle
      mode = 3
      wave = (wave1(1)+wave1(nlines))/2.
      call opacit (2,wave)
      call linlimit
      
c*****divide the lines into keepers and discards
      do j=1,nlines
         if (strength(j)/kaplam(jtau5) .ge. strongratio) then
             ns_lines = ns_lines+1
             strong(ns_lines) = wave1(j)
         elseif (strength(j)/kaplam(jtau5) .ge. weakratio) then
             nw_lines = nw_lines+1
             weak(nw_lines) = wave1(j)
         endif
      enddo
      if (nlines +nstrong .eq. 2500) then
         wave = wave1(nlines)
         call linlimit
         if (lim2line .lt. 0) then
             call inlines(2)
             call nearly(1)
         endif
         go to 1
      else
         do j=nlines+1, nlines+nstrong
            if (strength(j)/kaplam(jtau5) .ge. strongratio) then
               ns_lines = ns_lines+1
               strong(ns_lines) = wave1(j)
            elseif (strength(j)/kaplam(jtau5) .ge. weakratio) then
               nw_lines = nw_lines+1
               weak(nw_lines) = wave1(j)
            endif
         enddo
      endif
      call sortwave
      start = oldstart

      end

      subroutine sortwave

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Stokes.com'
      INTEGER I,J
      REAL*8 X
      DO 30 I=2,ns_lines
      X=strong(I)
      J=I
   10 J=J-1
      IF(J.EQ.0 .OR. strong(J).LT.X) GO TO 20
      strong(J+1)=strong(J)
      GO TO 10
   20 strong(J+1)=X
   30 CONTINUE

      do 60 I=2, nw_lines
      X=weak(I)
      J=I
   40 J=J-1
      if(J.eq.0 .or. weak(j).LT.X) GO TO 50
      weak(j+1)=weak(J)
      GO TO 40
   50 weak(j+1)=X
   60 CONTINUE

      END
