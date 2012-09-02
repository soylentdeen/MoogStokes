
      subroutine diskint
c******************************************************************************
c     This program synthesizes a section of spectrum and compares it
c     to an observation file.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      include 'Angles.com'


c*****examine the parameter file
      call params

      B_sph(1) = 1.0
      B_sph(2) = 0.0
      B_sph(3) = 0.0

      B_xyz(1) = 0.0
      B_xyz(2) = 0.0
      B_xyz(3) = 0.0

      inclination = 3.1415926/2.0

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
      if (plotopt .ne. 0) then
         nf3out = 22               
         lscreen = lscreen + 2
         array = 'SMOOTHED SYNTHESES OUTPUT'
         nchars = 25
         call infile ('output ',nf3out,'formatted  ',0,nchars,
     .                f3out,lscreen)
         if (f5out .ne. 'optional_output_file') then
            nf5out = 26
            lscreen = lscreen + 2
            array = 'POSTSCRIPT PLOT OUTPUT'
            nchars = 22
            call infile ('output ',nf5out,'formatted  ',0,nchars,
     .                   f5out,lscreen)
         endif
      endif
      if (iraf .ne. 0) then
         nf4out = 23               
         lscreen = lscreen + 2
         array = 'IRAF ("rtext") OUTPUT'
         nchars = 24
         call infile ('output ',nf4out,'formatted  ',0,nchars,
     .                f4out,lscreen)
      endif


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
      

c*****do the syntheses
10    isynth = 1
      isorun = 1
      nlines = 0
      mode = 3
      call inlines (1)
      call eqlib
      call nearly (1)
      nrings = 23
      ncells = 695
      cell_area = 4.0*3.14159262/ncells
      radtodeg = 180.0/3.1459262
      do i=1,nrings
         chi_start = (i-1)*3.14159262/nrings
         chi_stop = i*3.14159262/nrings
         azimuth = (chi_start + chi_stop)/2.0
         dchi = -cos(chi_stop) + cos(chi_start)
         ring_area = 2.0*3.14159262 * dchi ! total sterrad in semicircle
         n_cells = nint(ring_area/cell_area) ! # of cells in ring
         cell_a = ring_area/float(n_cells)
         dphi = 2.0*3.14159262/n_cells
c         write (*,*) wave, azimuth*radtodeg, n_cells
         do j=1,n_cells
            longitude = -3.14159262+(j-0.5)*dphi
            call computeRotations
c            write (*,*) phi_angle*radtodeg,chi_angle*radtodeg,
c     .                  longitude*radtodeg,viewing_angle*radtodeg
            if (viewing_angle*radtodeg .le. 90.0) THEN
               call synspectile
               call appendStokes(cell_a)
            endif
         enddo
      enddo

      do i = 1, len(spectrum)
         spectrum(i) = spectrum(i)/total_weight
      enddo

c*****now plot the spectrum
20    if (plotopt.eq.2 .and. specfileopt.gt.0) then
         nfobs = 33               
         lscreen = lscreen + 2
         array = 'THE OBSERVED SPECTRUM'
         nchars = 21
         if (specfileopt.eq.1 .or. specfileopt.eq.3) then
            call infile ('input  ',nfobs,'unformatted',2880,nchars,
     .                   fobs,lscreen)
         else
            call infile ('input  ',nfobs,'formatted  ',0,nchars,
     .                   fobs,lscreen)
         endif
      endif
      if (plotopt .ne. 0) then
c         call pltspec (lscreen,ncall)
      endif


c*****if the syntheses need to be redone: first rewind the output files,
c     then close/reopen line list(s), then rewrite model atmosphere output
      if (choice .eq. 'n') then
         call chabund
         if (choice .eq. 'x') go to 20
         rewind nf1out
         rewind nf2out
         if (nflines .ne. 0) then
            close (unit=nflines)
            open (unit=nflines,file=flines,access='sequential',
     .            form='formatted',blank='null',status='old',
     .            iostat=jstat,err=10)
         endif
         if (nfslines .ne. 0) then
            close (unit=nfslines)
            open (unit=nfslines,file=fslines,access='sequential',
     .            form='formatted',blank='null',status='old',
     .            iostat=jstat,err=10)
         endif
         if (plotopt .ne. 0) then
            rewind nf3out
         endif
         write (nf1out,1002) modtype
         if (modprintopt .ge. 1) then
            if (modtype .eq. 'begn      ' .or.
     .          modtype .eq. 'BEGN      ') write (nf1out,1003)
            write (nf1out,1102) moditle
            do i=1,ntau
               dummy1(i) = dlog10(pgas(i))
               dummy2(i) = dlog10(ne(i)*1.38054d-16*t(i))
            enddo
            write (nf1out,1103) wavref,(i,xref(i),tauref(i),t(i),
     .                          dummy1(i), pgas(i),dummy2(i),ne(i),
     .                          vturb(i),i=1,ntau)
            write (nf1out,1104)
            do i=1,95
               dummy1(i) = dlog10(xabund(i)) + 12.0
            enddo
            write (nf1out,1105) (names(i),i,dummy1(i),i=1,95)
            write (nf1out,1106) modprintopt, molopt, linprintopt, 
     .                          fluxintopt
            write (nf1out,1107) (kapref(i),i=1,ntau)
         endif
         linprintopt = linprintalt
         ncall = 2
         choice = '1'
         go to 10


c*****otherwise end the code gracefully
      else
         call finish (0)
      endif


c*****format statements
1001  format (/'KAPREF ARRAY:'/(6(1pd12.4)))
1002  format (13('-'),'MOOG OUTPUT FILE',10('-'),
     .        '(MOOG version from 23/04/07)',13('-')//
     .        'THE MODEL TYPE: ',a10)
1003  format ('   The Rosseland opacities and optical depths have ',
     .        'been read in')
1102  format (/'MODEL ATMOSPHERE HEADER:'/a80/)
1103  format ('INPUT ATMOSPHERE QUANTITIES',10x,
     .        '(reference wavelength =',f10.2,')'/3x,'i',2x,'xref',3x,
     .        'tauref',7x,'T',6x,'logPg',4x,'Pgas',6x,'logPe',
     .        5x,'Ne',9x,'Vturb'/
     .        (i4,0pf6.2,1pd11.4,0pf9.1,f8.3,1pd11.4,0pf8.3,
     .        1pd11.4,d11.2))
1104  format (/'INPUT ABUNDANCES: (log10 number densities, log H=12)'/
     .       '      Default solar abundances: Anders and Grevesse 1989')
1105  format (5(3x,a2,'(',i2,')=',f5.2))
1106  format (/'OPTIONS: atmosphere = ',i1,5x,'molecules  = ',i1/
     .        '         lines      = ',i1,5x,'flux/int   = ',i1)
1107  format (/'KAPREF ARRAY:'/(6(1pd12.4)))

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
      viewing_angle = acos(-B_xyz(1))
      mu = cos(viewing_angle)
      return
      end


      subroutine AppendStokes(cell_a)
c*****************************************************************
c     Adds an emergent spectrum (calculated at angles phi, chi) to
c     the resultant spectrum.  Calculated at a single wavelength.
c*****************************************************************

      implicit real*8 (a-h,o-z)
      include "Atmos.com"
      include "Linex.com"
      include "Angles.com"
      real*8 area, proj_area, limb_darkening
      real*8 surface_area, weight

      if ((1.0/(wave/10000.0)).lt.2.4) then
          alpha = -0.023 + 0.292/(wave/10000.0)
      else
          alpha = -0.507 + 0.441/(wave/10000.0)
      endif

      mu = cos(viewing_angle)

      limb_darkening = (1.0-(1.0-mu**alpha))

c*****   Calculate the projected area
      dotproduct = cos(viewing_angle)
c      dotproduct = sin(phi_angle)*sin(chi_angle)**2.0
      projected_area = cell_a*dotproduct/(4.0*3.14159262)

      weight = limb_darkening*projected_area
c      write (*,*) viewing_angle*180.0/3.14159262, projected_area
c      write (*,*) viewing_angle*180.0/3.14159262, mu, limb_darkening,
c     .            projected_area, weight

      do i = 1, len(d) 
          spectrum(i) = spectrum(i) + (1.0-d(i))*weight
      enddo


      total_weight = total_weight + weight

      end

