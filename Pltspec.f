

      subroutine pltspec (line,ncall)
c******************************************************************************
c     This subroutine controls the plotting of the synthesized and 
c     observed data
c******************************************************************************
 
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      real*8 xnum
      real*4 yyadd,xxadd,vveladd,yymult
      integer ncall



c*****initialize some variables, or re-set them to old values
3     if (ncall .eq. 1) then
         if (iscale .eq. 0) then
            choice = 's'
         else
            choice = '1'
         endif
         call plotremember (0)
         call plotremember (1)
      else
         call plotremember (4)
      endif


c*****make a special "choice" for grid syntheses, which go directly to
c*****postscript files if desired
      if (control.eq.'gridsyn' .or. control.eq.'gridend' .or.
     .    control.eq.'gridplo') then
         choice = 'g'
      endif


c*****begin by reading in an observed spectrum
      if (specfileopt.ge.1 .and. plotopt.ge.2 .and. ncall.eq.1) then
         call readobs (line)
         if (lount .eq. -1) then
            array = 'OBSERVED SPECTRUM FILE PROBLEM!  I QUIT.'
            istat = ivwrite (line+4,3,array,40)
            stop
         endif
      endif
      

c*****go through the option list; the routine may be exited at this point
      istat = ivcleof(4,1)
1     if (choice .eq. 'q') return


c*****or a default plot may be made upon entering the routine
      if (choice.eq.'1' .or. choice.eq.'g') then
c*****make cross correlation to line up synthetic and observed spectra 
c*****in velocity (wavelength) space.  Comment this line out if it is
c*****not working for you
         call smooth (-1,ncall)
c----------------------------------------------------------------------------
         if (plotopt.eq.2 .and. ncall.eq.1) then
            wmiddle = (start + sstop)/2.
            do i=1,lount
               xobs(i) = veladd*wmiddle/2.99795e+5 + xobs(i)
            enddo
c            call correl 
            do i=1,lount
               xobs(i) = deltavel*wmiddle/2.99795e+5 + xobs(i)
            enddo
         endif
c----------------------------------------------------------------------------
         if (ncall .eq. 1) then
            wmiddle = (start + sstop)/2.
            if (iunits .eq. 1) wmiddle = 1.d-4*wmiddle
            if (ymult .eq. 0.0) ymult = 1.0
            do i=1,lount
               yobs(i) = ymult*yobs(i)
               yobs(i) = yadd+yobs(i)
            enddo
            if (xlo .eq. xhi) then
               xlo = start
               xhi = sstop
               if (iunits .eq. 1) then
                  xlo = 1.d-4*xlo
                  xhi = 1.d-4*xhi
               endif
            endif
            if (ylo .eq. yhi) then
               ylo = 0.
               yhi = 1.1
            endif
            call plotremember (1)
         endif
         call boxit
      endif
            
         
c*****or the synthetic spectra may be resmoothed; if a problem develops in
c     a user-specified parameter (like a Gaussian FWHM that is too large),
c     then output a warning and let user decide what to do next
      if (choice .eq. 's') then
         call smooth (line+2,1)
2        if (smtype .eq. 'e') then
            array = 'REDO THE SMOOTHING (y/n)? '
            nchars = 26
            call getasci (nchars,line+9)
            smtype = chinfo(1:1)
            if (smtype .eq. 'n') then
               go to 100
            else
               istat = ivcleof (10,1)
               go to 2
            endif
         endif
         if (xlo .eq. 0.0 .and. xhi .eq. 0.0) then
            xlo = start
            xhi = sstop
            if (iunits .eq. 1) then
               xlo = 1.d-4*xlo
               xhi = 1.d-4*xhi
            endif
         endif
      endif


c*****or the observations may be rescaled
      if (choice .eq. 'r') then
         array = 'MULTIPLY THE OBSERVED POINTS BY WHAT FACTOR? '
         nchars = 45
         call getnum (nchars,13,xnum,yymult)
         ymult = ymult + yymult
         do i=1,lount
            yobs(i) = yymult*yobs(i)
         enddo
      endif


c*****or the observations may be shifted by an additive constant
      if (choice .eq. 'a') then
         array = 'ADD WHAT NUMBER TO THE OBSERVED POINTS? '
         nchars = 40
         call getnum (nchars,13,xnum,yyadd)
         yadd = yadd + yyadd
         do i=1,lount
            yobs(i) = yyadd + yobs(i)
         enddo
      endif


c*****or the observations may be shifted by a constant wavelength
      if (choice .eq. 'w') then
         array = 'SHIFT THE OBSERVED POINTS BY WHAT WAVELENGTH? '
         nchars = 46
         call getnum (nchars,13,xnum,xxadd)
         xadd = xadd + xxadd
         do i=1,lount
            xobs(i) = xxadd + xobs(i)
         enddo
      endif


c*****or the observations may be shifted by a constant velocity
      if (choice .eq. 'v') then
         array = 'SHIFT THE OBSERVED POINTS BY WHAT VELOCITY (KM/S)? '
         nchars = 51
         call getnum (nchars,13,xnum,vveladd)
         veladd = veladd * vveladd
         wmiddle = (start + sstop)/2.
         if (iunits .eq. 1) wmiddle = 1.d-4*wmiddle
         do i=1,lount
            xobs(i) = vveladd*wmiddle/2.99795e+5 + xobs(i)
         enddo
      endif


c*****or the plot boundaries may be changed
      if (choice .eq. 'c') then
         write (array,1001) xlo
         nchars = 29
         call getnum (nchars,15,xnum,xlo)
         write (array,1002) xhi
         nchars = 30
         call getnum (nchars,15,xnum,xhi)
         write (array,1003) ylo
         nchars = 34
         call getnum (nchars,15,xnum,ylo)
         write (array,1004) yhi
         nchars = 31
         call getnum (nchars,15,xnum,yhi)
         call boxit
      endif


c*****or the cross hairs can be used to zoom in on a part of the plot
C      if (choice .eq. 'z') then
C         array = 'MARK THE LOWER LEFT HAND CORNER WITH THE CURSOR'
C220      istat = ivcleof(13,1)
C         istat = ivwrite(13,3,array,47)
C         call pointcurs
C         xlo = xplotpos
C         ylo = yplotpos
C         array = 'MARK THE UPPER RIGHT HAND CORNER WITH THE CURSOR'
C212      istat = ivcleof(14,1)
C         istat = ivwrite (14,1,array,48)
C         call pointcurs
C         xhi = xplotpos
C         yhi = yplotpos
C         call boxit
C         if (iunits .eq. 1) then
C            xlo = 1.d-4*xlo
C            xhi = 1.d-4*xhi
C         endif
C         whichwin = '1of1'
C      endif



c*****or cursor position can be returned
C      if (choice .eq. 'p') then
C         array = 'MARK THE POSITION WITH THE CURSOR'
C 213     istat=ivcleof(21,1)
C         istat=ivwrite(13,3,array,34)
C         call drawcurs
C         go to 100
C      endif


c*****or the title of the model can be changed
C      if (choice .eq. 't') then
C         array = 'ENTER THE NEW TITLE'
C         istat = ivcleof(21,1)
C         istat = ivwrite (13,3,array,19)
C         read (*,33) moditle
C33       format(a80)
C      endif
   

c*****or the spectra can be replotted, with a separate plot showing the 
c     observed/synthtic spectrum differences
C      if (choice .eq. 'd') then
C         deviations = 1
C         whichwin = '2of2'
C      endif


c*****or the plot boundaries may be reset to the original values;
c     this is a basic starting over plot
C13    if (choice .eq. 'o') then
C         xlo = start
C         xhi = sstop
C         if (iunits .eq. 1) then
C            xlo = 1.d-4*xlo
C            xhi = 1.d-4*xhi
C         endif
C         xlo = origxlo
C         xhi = origxhi
C         ylo = origylo
C         yhi = origyhi
C         lim1obs = 1
C         lim2obs = lount
C         deviations = 0
C         call boxit
C      endif


c*****or the plot can simply be redone
C      if (choice .eq. 'm') then
C         go to 90
C      endif


c*****now either here make a hardcopy plot
C      if (choice .eq. 'h') then
C         if (control .eq. 'binary ') then
C            plotroutine = 'hard_land_bin '
C         else
C            plotroutine = 'hard_land_spec'
C         endif
C         if (deviations .eq. 0) then
C            whichwin = '1of1'
C         else
C            whichwin = '2of2'
C         endif
C         lscreen = 12
C         call makeplot (lscreen)
C         go to 100
C      endif


c*****or write the plot to a postscript file
C      if (choice.eq.'f' .or. choice.eq.'g') then
C         if (control .eq. 'binary ') then
C            plotroutine = 'file_land_bin '
C         else
C            plotroutine = 'file_land_spec'
C         endif
C         if (deviations .eq. 0) then
C            whichwin = '1of1'
C         else
C            whichwin = '2of2'
C         endif
C         lscreen = 12
C         call makeplot (lscreen)
C         if (choice .eq. 'g') then
C            return
C         else
C            go to 90
C         endif
C      endif

      
c*****or return to the calling routine in order to change abundances
C      if (choice .eq. 'n') then
C         call plotremember (3)
C         return
C      endif


c*****or add an additional uniform amount of flux, expressed in terms of
c     the current continuum flux; this only is approximately physically
c     correct if spectrograph smoothing is negligible compared to other
c     smoothing
      if (choice .eq. 'l') then
         write (array,*)
     .      'WHAT IS THE ADDITIONAL FLUX IN TERMS OF CONTINUUM [0.0]? '
         nchars = 57
         call getnum (nchars,15,xnum,addflux)
         call smooth (-1,ncall)
      endif
 
 
c*****or if total confusion has happened in the plotting, reset all parameters
c     to their original values, and replot
      if (choice .eq. 'u') then
         call plotremember (2)
         ncall = 1
         go to 3
      endif
 
 
c*****or plot on the terminal
90    if (control.eq.'gridsyn' .or. control.eq.'gridend' .or.
     .    control .eq. 'gridplo') return
      if (control .eq. 'binary ') then
         plotroutine = 'term_land_bin '
      else
         plotroutine = 'term_land_spec'
      endif
      lscreen = 12
      if (choice.eq.'f' .or. choice.eq.'h') choice = 'm'
      if (deviations .eq. 0) then
         whichwin = '1of1'
      else
         whichwin = '2of2'
      endif
      call makeplot (lscreen)


c*****finally, print the option table
100   istat = ivcleof (5,1)
      array = 'OPTIONS?    s=new smoothing     r=rescale obs.     '
      istat = ivwrite (5,3,array,48)
      array = '            a=add # to obs.     h=hardcopy         '
      istat = ivwrite (6,3,array,48)
      array = '            c=change bounds     q=quit             '
      istat = ivwrite (7,3,array,48)
      array = '            m=redo same plot    o=orig. plot bounds'
      istat = ivwrite (8,3,array,51)
      array = '            v=velocity shift    w=wavelength shift '
      istat = ivwrite (9,3,array,51)
      array = '            z=zoom in on plot   p=cursor position  '
      istat = ivwrite (10,3,array,51)
      array = '            t=change title      f=postscript file  '
      istat = ivwrite (11,3,array,51)
      array = '            n=new abundances    d=obs/syn deviation'
      istat = ivwrite (12,3,array,51)
      array = '            l=add veiling       u=undo all; replot '
      istat = ivwrite (13,3,array,51)
      array = 'What is your choice? '
      nchars = 21
      call getasci (nchars,15)
      choice = chinfo(1:1)


c*****reprint the option table if the choice is not understood
c     or take action on the choice
      if (choice.eq.'s' .or. choice.eq.'r' .or.
     .    choice.eq.'a' .or. choice.eq.'h' .or.
     .    choice.eq.'c' .or. choice.eq.'q' .or.
     .    choice.eq.'m' .or. choice.eq.'o' .or.
     .    choice.eq.'v' .or. choice.eq.'w' .or.
     .    choice.eq.'z' .or. choice.eq.'p' .or.
     .    choice.eq.'t' .or. choice.eq.'f' .or.
     .    choice.eq.'n' .or. choice.eq.'d' .or.
     .    choice.eq.'l' .or. choice.eq.'u') then
         go to 1
      else 
         go to 100
      endif


c*****format statements
1001     format ('LEFT WAVELENGTH (',f9.2,'): ')
1002     format ('RIGHT WAVELENGTH (',f9.2,'): ')
1003     format ('BOTTOM RELATIVE FLUX (',f9.2,'): ')
1004     format ('TOP RELATIVE FLUX (',f9.2,'): ')


      end
      




