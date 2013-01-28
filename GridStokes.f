
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
      integer n_cells, icell, testflag, wavecounter
      real*8 az_start, az_stop, daz, az, long, dlong
      real*8 ring_area, cell_area, cell_a, phi_ang, chi_ang, mu

      testflag = 1

      zeros(1,1) = 0.0
      zeros(1,2) = 0.0
      zeros(1,3) = 0.0
      zeros(2,1) = 0.0
      zeros(2,2) = 0.0
      zeros(2,3) = 0.0
      zeros(3,1) = 0.0
      zeros(3,2) = 0.0
      zeros(3,3) = 0.0

c*****examine the parameter file
1     call params
      write (*,*) fmodel
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
c      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)

c*****Open the output files for the Stokes Output
      nfAngles = 60               
c      lscreen = lscreen + 2
      array = 'ANGLES OUTPUT'
      nchars = 20
      call infile ('output ',nfAngles,'formatted  ',0,nchars,
     .             fAngles,lscreen)
      nfStokesI = 61               
c      lscreen = lscreen + 2
      array = 'Stokes I OUTPUT'
      nchars = 20
      call infile ('output ',nfStokesI,'formatted  ',0,nchars,
     .             fStokesI,lscreen)
      nfStokesQ = 62               
c      lscreen = lscreen + 2
      array = 'Stokes Q OUTPUT'
      nchars = 20
      call infile ('output ',nfStokesQ,'formatted  ',0,nchars,
     .             fStokesQ,lscreen)
      nfStokesU = 63               
c      lscreen = lscreen + 2
      array = 'Stokes U OUTPUT'
      nchars = 20
      call infile ('output ',nfStokesU,'formatted  ',0,nchars,
     .             fStokesU,lscreen)
      nfStokesV = 64               
c      lscreen = lscreen + 2
      array = 'Stokes V OUTPUT'
      nchars = 20
      call infile ('output ',nfStokesV,'formatted  ',0,nchars,
     .             fStokesV,lscreen)
      nfContinuum = 65               
c      lscreen = lscreen + 2
      array = 'Continuum OUTPUT'
      nchars = 20
      call infile ('output ',nfContinuum,'formatted  ',0,nchars,
     .             fContinuum,lscreen)

c*****open and read the model atmosphere file
      nfmodel = 30
c      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             trim(AtmosDir)//fmodel,lscreen)
      call inmodel


c*****open the line list file and the strong line list file
      nflines = 31
c      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
c         lscreen = lscreen + 2
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
      position_angle = 0.0

      write (nfAngles, 12347) ncells, nrings, inclination,
     .                        position_angle
      cell_area = 4.0*3.14159262/ncells
      radtodeg = 180.0/3.1459262
      icell = 1
      do i=1,nrings
         az_start = (i-1)*3.14159262/nrings - 3.14159262/2.0
         az_stop = i*3.14159262/nrings - 3.14159262/2.0
         az = (az_start + az_stop)/2.0
         daz= sin(az_stop) - sin(az_start)
         ring_area = 2.0*3.14159262 * daz ! total sterrad in circular ring
         n_cells = nint(ring_area/cell_area) ! # of cells in ring
         cell_a = ring_area/float(n_cells)
         dlong = 2.0*3.14159262/n_cells
         do j=1,n_cells
            long = -3.14159262+(j-0.5)*dlong
            call computeRotations(az, long, phi_ang, chi_ang, mu)
            if (mu .ge. 0.001) THEN
               phi_angle(icell) = phi_ang
               chi_angle(icell) = chi_ang
               azimuth(icell) = az
               longitude(icell) = long
               mus(icell) = mu
               write (nfAngles, 12346)icell, az, az_start, az_stop,long,
     .                dlong, phi_ang, chi_ang, mu
               icell = icell + 1
            endif
         enddo
      enddo

      icell = icell-1

      call wavegrid
c*****Read in the line list and calculate the equilibria
      call inlines (1)
      call eqlib
      call nearly (1)
      wavecounter = 1
c*****Perform the Synthesis
      wavl = 0.
      mode = 3
30    wave = wavelength(wavecounter)
c      wave = 4923.627
c      wave = 4923.927
      if (dabs(wave-wavl)/wave .ge. 0.001) then
         wavl = wave
         call opacit (2,wave)
      endif
20    call linlimit
      if (lim2line .lt. 0) then
          call inlines(2)
          call nearly (1)
          go to 20
      endif
      lim1 = lim1line
      lim2 = lim2line
      call calcopacities
      write (*,*) wave
      write (nfStokesI, 6520, advance='no') wave
      write (nfStokesQ, 6520, advance='no') wave
      write (nfStokesU, 6520, advance='no') wave
      write (nfStokesV, 6520, advance='no') wave
      write (nfContinuum, 6520, advance='no') wave
      if (testflag .eq. 1) then
c         call traceStokes(dble(0.0), dble(0.0), dble(1.0))
c         call traceStokes(dble(1.50707), dble(0.0), dble(1.0))
c         call traceStokes(dble(0.95532), dble(3.14159), dble(0.57735))
c         call traceStokes(dble(0.95532), dble(0.0), dble(0.57735))
c         call traceStokes(dble(1.5025), dble(4.712), dble(0.0682))

         call traceStokes(dble(0.698131), dble(0.0), dble(1.0))
c         call traceStokes(dble(0.27064), dble(0.0), dble(0.9636))
c         call traceStokes(dble(0.481286), dble(0.0), dble(0.88634))
c         call traceStokes(dble(0.640495), dble(0.0), dble(0.8018))
c         call traceStokes(dble(0.785389), dble(0.0), dble(0.7071))
c         call traceStokes(dble(0.929793), dble(0.0), dble(0.5979998))
c         call traceStokes(dble(1.089532), dble(0.0), dble(0.4629))
c         call traceStokes(dble(1.300206), dble(0.0), dble(0.2673))

         write (nfStokesI, 6521, advance='no') Stokes(1)/continuum
         write (nfStokesQ, 6521, advance='no') Stokes(2)/continuum
         write (nfStokesU, 6521, advance='no') Stokes(3)/continuum
         write (nfStokesV, 6521, advance='no') Stokes(4)/continuum
         write (nfContinuum, 6521, advance='no') continuum
      else
         do i = 1, icell
            call traceStokes(phi_angle(i), chi_angle(i), mus(i))
            write (nfStokesI, 6521, advance='no') Stokes(1)
            write (nfStokesQ, 6521, advance='no') Stokes(2)
            write (nfStokesU, 6521, advance='no') Stokes(3)
            write (nfStokesV, 6521, advance='no') Stokes(4)
            write (nfContinuum, 6521, advance='no') continuum
         enddo
      endif
      write (nfStokesI, *) ''
      write (nfStokesQ, *) ''
      write (nfStokesU, *) ''
      write (nfStokesV, *) ''
      write (nfContinuum, *) ''
      
      stepsize = dopp(nstrong, 50)*wave/2.997929e11
      wavecounter = wavecounter + 1
      if (wavecounter .le. nwave) then
          go to 30
      endif

c*****finish
      if (control .ne. 'gridend') then
         call finish (1)
         go to 1
      else
         call finish (0)
      endif
      return

c1001  format (a80)
c12345 format (f10.4,5e15.5)
12346 format (i5,8e16.5)
12347 format (2i5,2e16.5)
6520  format (f10.4)
6521  format (e16.8)
      end 


      subroutine computeRotations (az, long, phi_ang, chi_ang, mu)

      implicit real*8 (a-h,o-z)
      include "Angles.com"
      real*8 temp_A(3,3), Bmag, az, long, phi_ang, chi_ang, mu
      
      T_rho(1,1) = 0.0
      T_rho(1,2) = 0.0
      T_rho(1,3) = 1.0
      T_rho(2,1) = -cos(az)
      T_rho(2,2) = sin(az)
      T_rho(2,3) = 0.0
      T_rho(3,1) = sin(az)
      T_rho(3,2) = cos(az)
      T_rho(3,3) = 0.0

      T_eta(1,1) = cos(long+position_angle)
      T_eta(1,2) = -sin(long+position_angle)
      T_eta(1,3) = 0.0
      T_eta(2,1) = sin(long+position_angle)
      T_eta(2,2) = cos(long+position_angle)
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
      call dgemm('N','N',3,3,3,dble(1.0),T_eta,3,T_rho,
     .           3,dble(0.0),temp_A,3)

      call dgemm('N','N',3,3,3,dble(1.0),T_i,3,temp_A,
     .           3,dble(0.0),rotation_matrix,3)

      B_xyz(1) = 0.0
      B_xyz(2) = 0.0
      B_xyz(3) = 0.0
      call dgemm('N','N',3,3,3,dble(1.0),rotation_matrix,3,
     .           B_sph,3,dble(1.0),B_xyz,3)

      Bmag = sqrt(B_xyz(1)**2.0+B_xyz(2)**2.0+B_xyz(3)**2.0)
      phi_ang = acos(B_xyz(3)/Bmag)
      if (B_xyz(1) .gt. 0.0) then
          if (B_xyz(2) .gt. 0.0) then
              chi_ang = atan(B_xyz(2)/B_xyz(1))
          else
              chi_ang = 2.0*3.1415926+atan(B_xyz(2)/B_xyz(1))
          endif
      else
          chi_ang = 3.1415926+atan(B_xyz(2)/B_xyz(1))
      endif
      mu = B_xyz(3)
c      radtodeg = 180.0/3.14159
c      write (*,'(5f10.4)') az*radtodeg, long*radtodeg,
c     .           phi_ang*radtodeg, chi_ang*radtodeg, mu
c      write (*,'(3f10.4)') B_xyz
c      read (*,*)
      return
      end
