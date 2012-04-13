
      subroutine gridsyn
c******************************************************************************
c     This program can synthesize multiple sections of spectra for multiple
c     input model atmospheres
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      integer status
      integer system
      integer ntot_cells, nrings, ncells
      real*8 chi_start, chi_stop, dchi, cell_area, cell_a, ring_area
      real*8 dphi, phi, chi_angle


c*****examine the parameter file
1     call params
      write (*,*) " This is a test.  This is only a test "
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
c      if (plotopt .gt. 0) then
c         nf3out = 22               
c         lscreen = lscreen + 2
c         array = 'SMOOTHED SYNTHESES OUTPUT'
c         nchars = 25
c         call infile ('output ',nf3out,'formatted  ',0,nchars,
c     .                f3out,lscreen)
c         nf5out = 26
c         lscreen = lscreen + 2
c         array = 'POSTSCRIPT PLOT OUTPUT'
c         nchars = 22
c         call infile ('output ',nf5out,'formatted  ',0,nchars,
c     .                f5out,lscreen)
c      endif
c      if (plotopt .gt. 1) then
c         nf6out = 27
c         lscreen = lscreen + 2
c         array = 'SPECTRUM COMPARISON OUTPUT'
c         nchars = 27
c         call infile ('output ',nf6out,'formatted  ',0,nchars,
c     .                f6out,lscreen)
c      endif
c      if (iraf .ne. 0) then
c         nf4out = 23               
c         lscreen = lscreen + 2
c         array = 'IRAF ("rtext") OUTPUT'
c         nchars = 24
c         call infile ('output ',nf4out,'formatted  ',0,nchars,
c     .                f4out,lscreen)
c      endif


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

      nrings = 20
      ntot_cells = 350
      cell_area = 4.0*3.14159262/ntot_cells
      

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
            call inlines (1)
            call eqlib
            call nearly (1)
            do i=1,nrings
               chi_start = (i-1)*3.14159262/nrings
               chi_stop = i*3.14159262/nrings
               dchi = -cos(chi_stop) + cos(chi_start)
               ring_area = 3.14159262 * dchi ! total sterrad in semicircle
               ncells = nint(ring_area/cell_area) ! number of cells in ring
               cell_a = ring_area/float(ncells)
               dphi = 3.14159262/ncells
               do j=1,ncells
                  phi_angle = -3.14159262/2.0+(j-0.5)*dphi
                  chi_angle = (chi_start + chi_stop)/2.0
                  isynth = n
                  isorun = 1
                  start = oldstart
                  sstop = oldstop
                  mode = 3
                  call synspec (dble(phi_angle), dble(chi_angle))
                  write (*,*) phi_angle*180.0/3.14159,chi_angle*180.0/3.14159
c                  status = system('gnuplot gnuplot.commands')
                  linprintopt = 0
              enddo
            enddo
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






