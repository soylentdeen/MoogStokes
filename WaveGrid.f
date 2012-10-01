
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
      real*4 shortnum


c*****examine the parameter file

c      array = 'GIVE THE MINIMUM LINE/CONTINUUM OPACITY RATIO TO KEEP: '
      xratio = 0.1
      

c*****compute the line opacities
      call inlines (1)
      call eqlib
1     call nearly (1)


c*****calculate continuum quantities at the line list wavelength middle
      wave = (wave1(1)+wave1(nlines))/2.
      call opacit (2,wave)
      
c*****   Generate original wavelength grid

      deltawave = 0.5
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
             do k=-10,10
                 wavelength(nwave) = wave1(j)+0.5*asin(real(k/10.0))
                 nwave = nwave+1
             enddo
             wavelength(nwave) = wave1(j)+ 0.001
             nwave = nwave+1
             wavelength(nwave) = wave1(j)- 0.001
             nwave = nwave+1
         endif
      enddo
      if (nlines +nstrong .eq. 2500) then
         call inlines (6)
         go to 1
      endif
      nwave = nwave -1

      do i=1, nwave
          write (*,*) wavelength(i)
      enddo
      read (*,*)
      call sortwave
      do i=1, nwave
          write (*,*) wavelength(i)
      enddo
      read (*,*)

c*****format statements
1001  format (/'DESIRED LINE-TO-CONTINUUM MINIMUM OPACITY RATIO: ', 
     .        1pe10.2)
1002  format ('THIS IS THE KEEPER LINE LIST')
1003  format ('THIS IS THE DISCARDED LINE LIST')
1004  format (f10.4, f10.1, f10.3, f10.3, 30x, f9.1)
1005  format (f10.4, f10.1, f10.3, f10.3, 10x, f10.3, 10x, f9.1)
1006  format ('  kaplam from 1 to ntau at wavelength',f10.2/
     .        (6(1pd12.4)))
1007  format (f10.4, f10.4, f10.3, f10.3, 30x, f9.1)
1008  format (f10.4, f10.5, f10.3, f10.3, 10x, f10.3, 10x, f9.1)


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
         if(abs(wavelength(J)-temp(I)) .gt. 0.005) THEN
            I = I+1
            temp(I) = (wavelength(J)+temp(I-1))/2.0
            wavelength(J) = 0.0
         endif
      enddo
      nwave = I
      do J=1, nwave
         wavelength(J)=temp(J)
      enddo
      END
