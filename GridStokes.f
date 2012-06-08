
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
      include 'Stokes.com'
      include 'Angles.com'
      integer n_cells
      real*8 chi_start, chi_stop, dchi
      real*8 ring_area, cell_area, cell_a


c*****examine the parameter file
1     call params
      linprintopt = linprintalt

      zeros(1,1) = 0.0
      zeros(1,2) = 0.0
      zeros(1,3) = 0.0
      zeros(2,1) = 0.0
      zeros(2,2) = 0.0
      zeros(2,3) = 0.0
      zeros(3,1) = 0.0
      zeros(3,2) = 0.0
      zeros(3,3) = 0.0


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

      B_sph(1) = 1.0
      B_sph(2) = 0.0
      B_sph(3) = 0.0
      
      B_xyz(1) = 0.0
      B_xyz(2) = 0.0
      B_xyz(3) = 0.0

      inclination = 3.1415926/2.0
      open(unit=nf11out, file=f11out)

c*****Read in the line list and calculate the equilibria
      call inlines (1)
      call eqlib
      call nearly (1)
c*****Perform the Synthesis
      wave = start
      wavl = 0.
      mode = 3
30    if (dabs(wave-wavl)/wave .ge. 0.001) then
         wavl = wave
         call opacit (2,wave)
      endif
20    call linlimit
      if (lim2line .lt. 0) then
          call inlines(2)
          call nearly (1)
          go to 20
      endif
      Stokes_I = 0.0
      Stokes_Q = 0.0
      Stokes_U = 0.0
      Stokes_V = 0.0
      total_weight = 0.0
      lim1 = lim1line
      lim2 = lim2line
      call calcopacities
      nrings = 23
      ncells = 700
      cell_area = 4.0*3.14159262/ncells
c      do i=1,nrings
c         chi_start = (i-1)*3.14159262/nrings
c         chi_stop = i*3.14159262/nrings
c         azimuth = (chi_start + chi_stop)/2.0
c         dchi = -cos(chi_stop) + cos(chi_start)
c         ring_area = 3.14159262 * dchi ! total sterrad in semicircle
c         n_cells = nint(ring_area/cell_area) ! # of cells in ring
c         cell_a = ring_area/float(n_cells)
c         dphi = 3.14159262/n_cells
c         do j=1,n_cells
c            longitude = -3.14159262/2.0+(j-0.5)*dphi
c            call computeRotations
c            call delo
c            call appendStokes(cell_a)
c         enddo
c      enddo
      azimuth = 3.14159262/2.0
      longitude = 0.0
      call computeRotations
      chi_angle = dble(0.0)
c      write (*,*) phi_angle, chi_angle, viewing_angle
      call delo
c     call rungeKutta
c      Stokes_I = Stokes_I/total_weight
c      Stokes_Q = Stokes_Q/total_weight
c      Stokes_U = Stokes_U/total_weight
c      Stokes_V = Stokes_V/total_weight
      Stokes_I = Stokes(1)!/continuum
      Stokes_Q = Stokes(2)!/continuum
      Stokes_U = Stokes(3)!/continuum
      Stokes_V = Stokes(4)!/continuum

c      write (*,*) wave, Stokes_I, Stokes_Q, Stokes_U, Stokes_V
      write (nf11out,12345) wave, Stokes_I, Stokes_Q, Stokes_U,Stokes_V,
     .      continuum
c      stepsize = dopp(nstrong, 50)*wave/2.997929e10/2.0
      stepsize = dopp(nstrong, 50)*wave/2.997929e11
      wave = wave + stepsize
      if (wave .le. sstop) then
          go to 30
      endif

c      control = 'gridend'


c*****finish
      if (control .ne. 'gridend') then
         call finish (1)
         go to 1
      else
         call finish (0)
      endif
      close(nf11out)
      return

12345 format (f10.4,5e15.5)
      end 


      subroutine computeRotations

      implicit real*8 (a-h,o-z)
      include "Angles.com"
      real*8 temp_A(3,3), Bmag

      T_rho(1,1) = 0.0
      T_rho(1,2) = 0.0
      T_rho(1,3) = 1.0
      T_rho(2,1) = -cos(longitude)
      T_rho(2,2) = sin(longitude)
      T_rho(2,3) = 0.0
      T_rho(3,1) = sin(longitude)
      T_rho(3,2) = cos(longitude)
      T_rho(3,3) = 0.0

      T_eta(1,1) = cos(azimuth+position_angle)
      T_eta(1,2) = -sin(azimuth+position_angle)
      T_eta(1,3) = 0.0
      T_eta(2,1) = sin(azimuth+position_angle)
      T_eta(2,2) = cos(azimuth+position_angle)
      T_eta(2,3) = 0.0
      T_eta(3,1) = 0.0
      T_eta(3,2) = 0.0
      T_eta(3,3) = 1.0

      T_i(1,1) = 1.0
      T_i(1,2) = 0.0
      T_i(1,3) = 0.0
      T_i(2,1) = 0.0
      T_i(2,2) = cos(inclination)
      T_i(2,3) = sin(inclination)
      T_i(3,1) = 0.0
      T_i(3,2) = -sin(inclination)
      T_i(3,3) = cos(inclination)

      call dcopy(9,zeros,1,temp_A,1)
      call dcopy(9,zeros,1,rotation_matrix,1)
      call dgemm('T','N',3,3,3,dble(1.0),T_rho,3,T_eta,
     .           3,dble(0.0),temp_A,3)
      call dgemm('N','N',3,3,3,dble(1.0),T_i,3,temp_A,
     .           3,dble(0.0),rotation_matrix,3)

      B_xyz(1) = 0.0
      B_xyz(2) = 0.0
      B_xyz(3) = 0.0
      call dgemm('N','N',3,3,3,dble(1.0),rotation_matrix,3,
     .           B_sph,3,dble(1.0),B_xyz,3)

      Bmag = sqrt(B_xyz(1)**2.0+B_xyz(2)**2.0+B_xyz(3)**2.0)
      phi_angle = acos(-B_xyz(1)/Bmag)
      chi_angle = atan(B_xyz(2)/B_xyz(3))
c      phi_angle = acos(B_xyz(3)/Bmag)
c      chi_angle = atan(B_xyz(2)/B_xyz(1))
      viewing_angle = acos(-B_xyz(1))
c      write (*,*) B_xyz, chi_angle, phi_angle, viewing_angle
c      read (*,*)
      return
      end
