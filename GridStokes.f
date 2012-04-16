
      subroutine gridstokes
c******************************************************************************
c     This program can synthesize the full stokes parameters for multiple
c     sections of spectra for multiple input model atmospheres
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'


c*****examine the parameter file
1     call params
      linprintopt = linprintalt


c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra
      nf1out = 20     
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21               
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)


c*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****open the line list file and the strong line list file
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
         lscreen = lscreen + 2
         array = 'THE STRONG LINE LIST'
         nchars = 20
         call infile ('input  ',nfslines,'formatted  ',0,nchars,
     .                 fslines,lscreen)
      endif
      

c*****Read in the line list
      call inlines (1)
      call eqlib
      call nearly (1)
c      call genwave
      do i=1,nwave
         wave = wavelength(i)
         call calcopacities
         do j=1,nrings
            chi_start = (i-1)*3.14159262/nrings
            chi_stop = i*3.14159262/nrings
            chi_angle = (chi_start + chi_stop)/2.0
            dchi = -cos(chi_stop) + cos(chi_start)
            ring_area = 3.14159262 * dchi ! total sterrad in semicircle
            ncells = nint(ring_area/cell_area) ! # of cells in ring
            cell_a = ring_area/float(ncells)
            dphi = 3.14159262/ncells
            do k=1,ncells
               phi_angle = -3.14159262/2.0+(j-0.5)*dphi
               call delo(phi_angle, chi_angle, Stokes)
               call appendStokes(phi_angle, chi_angle, cell_a, Stokes)
            enddo
         enddo
      enddo


c*****do the syntheses
      if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isorun = 1
         nlines = 0
         mode = 3
         call inlines (1)
         call eqlib
         call nearly (1)
         call synspec (dble(0.7853), dble(1.0))
c         call synspec (dble(0.7853), dble(-0.7853))
      else
         do n=1,numatomsyn
            isynth = n
            isorun = 1
            start = oldstart
            sstop = oldstop
            mode = 3
            call inlines (1)
            call eqlib
            call nearly (1)
c            call synspec (dble(0.0), dble(0.0))
            call synspec (dble(1.3853), dble(0.0))
c            call synspec (dble(0.7853), dble(-0.7853))
            linprintopt = 0
         enddo
      endif
         
c*****now plot the spectrum
c      if (plotopt.eq.2 .and. specfileopt.gt.0) then
c         nfobs = 33               
c         lscreen = lscreen + 2
c         array = 'THE OBSERVED SPECTRUM'
c         nchars = 21
c         if (plotopt.eq.1 .or. specfileopt.eq.3) then
c            call infile ('input  ',nfobs,'unformatted',2880,nchars,
c     .                   fobs,lscreen)
c         else
c            call infile ('input  ',nfobs,'formatted  ',0,nchars,
c     .                   fobs,lscreen)
c         endif
c      endif
c      if (plotopt .ne. 0) then
c         ncall = 1
c         call pltspec (lscreen,ncall)
c      endif


c*****finish
      if (control .ne. 'gridend') then
         call finish (1)
         go to 1
      else
         call finish (0)
      endif
      return

      end 






