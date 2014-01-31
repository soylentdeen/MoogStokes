
      subroutine params 
c****************************************************************************** 
c     This subroutine reads in the commands from the parameter file 
c****************************************************************************** 

      implicit real*8 (a-h,o-z) 
      include 'Atmos.com' 
      include 'Linex.com' 
      include 'Factor.com' 
      include 'Pstuff.com' 
      include 'Obspars.com'
      include 'Mol.com'
      include 'Multistar.com'
      include 'Angles.com'
      include 'Stokes.com'
      real*8 deltalogab(5), temp
      character keyword*20, sandbox*80
      character arrayz*80, stokes_base*80
      integer kk
      data newcount, linecount /0, 0/

      if (linecount .eq. 0) oldcount = 0


c  IF DOING MULTIPLE RUNS: if this is not the first reading of the
c  parameter file, then read down to correct place in parameter file, 
c  using "linecount", and then skip the re-initialization of the 
c  various variables
      rewind nfparam
      read (nfparam,1001,end=100) arrayz
      if (linecount .ne. 0) then
         do i=1,linecount
            read (nfparam,1001,end=100) arrayz
         enddo
         go to 4
      else
         outflag = 0
         atmosflag = 0
      endif


c  INITIALIZE SOME VARIABLES: output file names and output file numbers
c  f1out, nf1out:    verbose standard output
c  f2out, nf2out:    raw synthetic spectra, or summary abundances, or
c                    summary curves-of-growth
c  f3out, nf3out:    smoothed synthetic spectra
c  f4out, nf4out:    IRAF-text-style synthetic spectra output
c  f5out, nf5out:    postscript plot output
c  f6out, nf6out:    synthetic/observed comparison text output (gridsyn)
c  f7out, nf7out:    raw synthetic spectrum of primary of a binary, 
c                    output summary run table, or
c                    to be used for some other purpose as needed; at
c                    present not given keyword input
c  f8out, nf8out:    raw synthetic spectrum of secondary of a binary,
c                    or the kept lines in a line weedout, or for
c                    some other purpose as needed
c  f9out, nf9out:    combined smoothed synthetic spectrum of a binary star,
c                    or the discarded lines in a line weedout, or the
c                    mean raw synthesis for a population of stars, or
c                    some other purpose as needed
c                    to be used for some other purpose as needed`
c  f10out, nf10out:  combined smoothed synthetic spectrum of a binary star,
c                    to be used for some other purpose as needed`
      f1out =   'no_filename_given'
      f2out =   'no_filename_given'
      f3out =   'no_filename_given'
      f4out =   'no_filename_given'
      f5out =   'optional_output_file'
      f6out =   'optional_output_file'
      f7out =   'no_filename_given'
      f8out =   'no_filename_given'
      f9out =   'no_filename_given'
      f10out =  'no_filename_given'
      fAngles = 'Stokes.angles'
      fStokesI = 'Stokes.spectrum_I'
      fStokesQ = 'Stokes.spectrum_Q'
      fStokesU = 'Stokes.spectrum_U'
      fStokesV = 'Stokes.spectrum_V'
      fContinuum = 'Stokes.spectrum_continuum'
      nf1out =   0
      nf2out =   0
      nf3out =   0
      nf4out =   0
      nf5out =   0
      nf6out =   0
      nf7out =   0
      nf8out =   0
      nf9out =   0
      nf10out =  0
      modelnum = 0
      nfAngles = 0
      nfStokesI = 0
      nfStokesQ = 0
      nfStokesU = 0
      nfStokesV = 0
      nfContinuum = 0
      modelnum = 0


c  INITIALIZE SOME VARIABLES: input file names and input file numbers
      fmodel =  'no_filename_given'
      flines =  'no_filename_given'
      fslines = 'no_filename_given'
      fobs =    'no_filename_given'
      ftable =  'no_filename_given' 
      nfmodel =  0 
      nflines =  0
      nfslines = 0
      nfobs =    0
      nftable =  0


c  INITIALIZE SOME VARIABLES: terminal type;
c  set smterm = ' ' for the plotting package 'sm'
      smterm = 'x11'
      

c  INITIALIZE SOME VARIABLES: 
c  atmosphere data printing:              modprintopt  [old ipr(1)]
c  molecular equilibrium:                 molopt       [old ipr(2)]
c  line data:                             linprintopt  [old ipr(3)] 
c  flux/intensity                         fluxintopt   [old ipr(4)]
c  synthesis/single line computations:    [deleted]    [old ipr(5)] 
c  line formation level                   [deleted]    [old ipr(6)]
c  plotting:                              plotopt      [old ipr(7)] 
c  damping                                dampingopt   [old ipr(8)]
c  observed spectrum file type            specfileopt  [old ipr(9)]
c  format of linelist (formatted/not)     linfileopt
c  source function with scat+abs          scatopt
      modprintopt  = 1
      molopt       = 1
      linprintopt  = 1
      fluxintopt   = 0
      plotopt      = 0
      dampingopt   = 0
      specfileopt  = 0
      linfileopt   = 0
      iunits       = 0
      itru         = 0
      iscale       = 0
      nlines       = 0
      iraf         = 0
      histoyes     = 0
      byteswap     = 0
      deviations   = 0
      scatopt      = 0
      gfstyle      = 0
      maxshift     = 0
      dostrong     = 0
 

c  INITIALIZE SOME VARIABLES:
c   if fudge is less than or equal to 0 then it does not scale it...make it
c   -1.0 just to be sure we dont have floating point problems with 0.0
      fudge = -1.0


c  INITIALIZE SOME VARIABLES: spectrum run parameters
      start = 0.
      sstop = 0.
      step = 0.
      delta = 0.
      cogatom = 0.


c  INITIALIZE SOME VARIABLES: line limit parameters
      ncurve = 0
      lim1line = 0
      lim2line = 0
      lim1obs = 0
      lim2obs = 0
      lim1 = 0
      lim2 = 0


c  INITIALIZE SOME VARIABLES: spectroscopic binary parameters
      deltaradvel = 0.
      lumratio = 1.

c  INITIALIZE SOME VARIABLES: MOOGStokes specific stellar geometry
      inclination = 3.1415926/2.0
      position_angle = 0.0
      diskflag = 1             ! 1 => annuli, 0 => disco-ball
      ncells = 695
      nrings = 23

c  INITIALIZE SOME VARIABLES: MOOGStokes test flag
      testflag = 0

c  INITIALIZE SOME VARIABLE: WaveGrid parameters
      beta_strong = 5.0
      beta_weak = 7.0
      R_strong = 0.5
      R_weak = 0.25

c  INITIALIZE SOME VARIABLES: elements with special abundance data; for
c  practical reasons, keyword "abundances" and "isotopes" must be specified
c  in each RUN; they must be reset each time.
4     neq = 0
      numpecatom = 0
      numatomsyn = 0
      newnumpecatom = 0
      newnumatomsyn = 0
      ninetynineflag = 0
      do i=1,95
         pec(i) = 0
         newpec(i) = 0
         do j=1,5
            abfactor(j) = 0.
            pecabund(i,j) = 0.
            newpecabund(i,j) = 0.
         enddo
      enddo
      numiso = 0
      numisosyn = 0
      newnumiso = 0
      newnumisosyn = 0
      do i=1,20
         isotope(i) = 0.
         newisotope(i) = 0.
         do j=1,5
            isoabund(i,j) = 0.
            newisoabund(i,j) = 0.
         enddo
      enddo


c  read a line of the parameter file into "arrayz"; decode it to get the 
c  key word, which is always within the first 20 characters;  dump the 
c  rest of arrayz into "array"
5     write (array,1007)
      read (nfparam,1001,end=98) arrayz
      linecount = linecount + 1
      i=index(arrayz,' ')
      keyword = arrayz(1:i-1)
      array = arrayz(i:)


c  keyword 'RUN' signals that there are either multiple syntheses being
c  done or multiple comparisons with observed spectra
      if (keyword .eq. 'RUN') then
         read (array,*) newcount
         if (newcount .gt. oldcount+1) then
            linecount = linecount - 1
            oldcount = syncount
            go to 100
         else
            syncount = newcount
            go to 5
         endif
      endif


c  keyword 'freeform' indicates whether or not the linelist will be read
c  in under format control (7e10.3) or will be free-form.  If freeform = 0,
c  the default value, then the old-style formatted input will be used;
c  If freeform = 1, unformatted read will be used, BUT the user must then
c  give values for all quantities (that is, explicit zeros will need to
c  be put instead of blank spaces.
c     MOOGStokes Modification - If freeform = 2, then Inlines will
c      expect 100 character lines in the linelist (delta_mj,
c      radiative_damping_coeff, C4_coefficient)
      if     (keyword .eq. 'freeform') then
         read (array,*) linfileopt
 
 
c  keyword 'standard_out' controls the name of the verbose standard output
      elseif (keyword .eq. 'standard_out') then
         read (array,*) f1out


c  keyword 'summary_out' controls the name of either the EW summary or
c  the raw synthesis output
      elseif (keyword .eq. 'summary_out') then
         read (array,*) f2out


c  keyword 'hardpost_out' controls the name of a postscript plot output
      elseif (keyword .eq. 'hardpost_out') then
         read (array,*) f5out


c  keyword 'speccomp_out' controls the name of a text file containing the
c  comparisons (wavelength shifts, sigmas, etc.) between observed and
c  synthetic spectra
      elseif (keyword .eq. 'speccomp_out') then
         read (array,*) f6out


c  keyword 'bin_raw_out' controls the name of a file containing the
c  raw synthesis of a spectroscopic binary, with an appropriate velocity 
c  difference and luminosity ratio dialed in
      elseif (keyword .eq. 'bin_raw_out') then
         read (array,*) f9out


c  keyword 'bin_smo_out' controls the name of a file containing the
c  smoothed synthesis of a spectroscopic binary
      elseif (keyword .eq. 'bin_smo_out') then
         read (array,*) f10out

 
c  keyword 'summary_in' controls the name of the raw synthesis file,
c  created previously, that will be read in for plotting purposes
      elseif (keyword .eq. 'summary_in') then
         read (array,*) f2out


c  keyword 'smoothed_out' controls the name of the smoothed synthesis output
      elseif (keyword .eq. 'smoothed_out') then
         read (array,*) f3out
c         sandbox = Jtrim(OutDir)//Jtrim(stokes_base)
c         f3out = Jtrim(sandbox)//'.moog'

c  keyword 'atmos_dir' specifies the location of the model atmosphere'
      elseif (keyword .eq. 'atmos_dir') then
         read (array,*) sandbox
         l = len(sandbox)
         AtmosDir = sandbox(1:l)
         atmosflag = 1

c  keyword 'out_dir' specifies the location of the model atmosphere'
      elseif (keyword .eq. 'out_dir') then
         read (array,*) sandbox
         OutDir = sandbox(1:len_trim(sandbox))
         outflag = 1

c  keyword 'stokes_out' controls the base name of all the Stokes-related output
      elseif (keyword .eq. 'stokes_out') then
         read (array,*) stokes_base
         if (outflag .eq. 1) then
             sandbox =OutDir(1:len_trim(Outdir))//
     ,               stokes_base(1:len_trim(stokes_base))
         else
             sandbox = stokes_base(1:len_trim(stokes_base))
         endif
         fAngles = sandbox(1:len_trim(sandbox))//'.angles'
         fStokesI = sandbox(1:len_trim(sandbox))//'.spectrum_I'
         fStokesQ = sandbox(1:len_trim(sandbox))//'.spectrum_Q'
         fStokesU = sandbox(1:len_trim(sandbox))//'.spectrum_U'
         fStokesV = sandbox(1:len_trim(sandbox))//'.spectrum_V'
         fContinuum = sandbox(1:len_trim(sandbox))//'.continuum'

c  keyword 'nrings' controls the number of latitude belts to synthesize in the
c        discoball (diskflag = 0)
      elseif (keyword .eq. 'nrings') then
         read (array,*) nrings

c  keyword 'ncells' controls the number of cells to synthesize in the
c        discoball (diskflag = 0)
      elseif (keyword .eq. 'ncells') then
         read (array,*) ncells

c  keyword 'inclination' controls the inclination of the stellar
c          rotation axis
      elseif (keyword .eq. 'inclination') then
         read (array,*) inclination

c  keyword 'pos_ang' controls the position angle of the stellar rotation
      elseif (keyword .eq. 'pos_ang') then
         read (array,*) position_angle

c  keyword 'diskflag' controls whether MOOGStokes calculates annuli (1) or a
c         disco-ball (0)
      elseif (keyword .eq. 'diskflag') then
         read (array,*) diskflag

c  keyword 'testflag' controls whether MOOGStokes calculates full disk
c         spectrum (0) or just a single pencil beam (1)
      elseif (keyword .eq. 'testflag') then
         read (array,*) testflag

c  keyword 'beta_strong' controls wavelength stepsize around strong lines
      elseif (keyword .eq. 'beta_strong') then
         read (array,*) beta_strong

c  keyword 'R_strong' controls wavelength stepsize around strong lines
      elseif (keyword .eq. 'R_strong') then
         read (array,*) R_strong

c  keyword 'beta_weak' controls wavelength stepsize around weak lines
      elseif (keyword .eq. 'beta_weak') then
         read (array,*) beta_weak

c  keyword 'R_weak' controls wavelength stepsize around weak lines
      elseif (keyword .eq. 'R_weak') then
         read (array,*) R_weak

c  keyword 'viewang' controls the viewing angle (measured in Degrees)
      elseif (keyword .eq. 'viewang') then
         read (array,*) temp 
         viewang = dble(3.14159262/180.0)*temp

c  keyword 'keeplines_out' controls the name of the list of kept lines
c  for future synthetic spectrum runs
      elseif (keyword .eq. 'keeplines_out') then
         read (array,*) f8out


c  keyword 'tosslines_out' controls the name of the list of discarded lines
c  that are too weak to keep in future synthetic spectrum runs
      elseif (keyword .eq. 'tosslines_out') then
         read (array,*) f9out


c  keyword 'iraf_out' controls the name of the optional IRAF output
      elseif (keyword .eq. 'iraf_out') then
         read (array,*) f4out


c  keyword 'model_in' controls the name of input model atmosphere file
      elseif (keyword .eq. 'model_in') then
         read (array,*) fmodel
         if (atmosflag .eq. 1) then
             sandbox = AtmosDir(1:len_trim(AtmosDir))//
     .          fmodel(1:len_trim(fmodel))
         else
             sandbox = fmodel(1:len_trim(fmodel))
         endif
         fmodel = sandbox


c  keyword 'lines_in' controls the name of the input line list
      elseif (keyword .eq. 'lines_in') then
         read (array,*) flines


c  keyword 'stronglines_in' controls the name of the input strong line list
      elseif (keyword .eq. 'stronglines_in') then
         read (array,*) fslines


c  keyword 'observed_in' controls the name of the input observed spectrum
      elseif (keyword .eq. 'observed_in') then
         read (array,*) fobs


c  keyword 'table_in' controls the name of the extra input instruction file
      elseif (keyword .eq. 'table_in   ') then
         read (array,*) ftable


c  keyword 'table_out' controls the name of the extra input instruction file
      elseif (keyword .eq. 'table_out  ') then
         read (array,*) f7out


c  keyword 'popsyn_out' controls the name of the extra input instruction file
      elseif (keyword .eq. 'popsyn_out ') then
         read (array,*) f9out


c  keyword 'rawbin_out ' controls the name of the input observed spectrum
      elseif (keyword .eq. 'rawbin_out ') then
         read (array,*) f9out


c  keyword 'smoobin_out' controls the name of the input observed spectrum
      elseif (keyword .eq. 'smoobin_out') then
         read (array,*) f10out


c  keyword 'atmosphere' controls the output of atmosphere quantities
c           0 = do not print out the atmosphere
c           1 = print out the standard things about an atmsophere
c           2 = print standard things and additional stuff like continuous
c                    opacities, etc.
      elseif (keyword .eq. 'atmosphere') then
         read (array,*) modprintopt


c  keyword 'molecules ' controls the molecular equilibrium calculations
c           0 = do not do molecular equilibrium
c           1 = do molecular equilibrium but do not print results
c           2 = do molecular equilibrium and print results
      elseif (keyword .eq. 'molecules') then
         read (array,*) molopt
         if (molopt .eq. 0) then
            nchars = 46
            array = 'WARNING: molecular eq. always done; OK (y/n)? '
            call getasci (nchars,l0)
            if (chinfo(1:1) .eq. 'n') then
               stop
            else
               molopt = 1
            endif
         endif


c  keyword 'deviations' controls whether, for synthetic spectrum computations,
c  an 'obs-comp' plot will be made in addition to the normal spectrum plot
c           0 = do not plot the obs-comp plot
c           1 = plot the obs-comp plot
      elseif (keyword .eq. 'deviations') then
         read (array,*) deviations


c  keyword 'lines     ' controls the output of line data
c           0 = print out nothing about the input lines
c           1 = print out standard information about the input line list
c           2 = gory line data print (usually for diagnostic purposes)
      elseif (keyword .eq. 'lines') then
         read (array,*) linprintopt
         linprintalt = linprintopt


c  keyword 'loggf     ' controls the output of line data
c           0 = straight gf values
c           1 = base-10 logarithms of the gf values
      elseif (keyword .eq. 'gfstyle') then
         read (array,*) gfstyle


c  keyword 'plotpars  ' allows you to set all of the plotting 
c     parameters if you know them in advance
c     0 = user supplied in plotting routine (default)
c     1 = given in following lines as follows 
c xlow         xhi         ylo       yhi
c vshift       lamshift    obsadd    obsmult
c smooth-type  FWHM-Gauss  vsini     L.D.C.    FWHM-Macro     FWHM-Loren
      elseif (keyword .eq. 'plotpars') then
         read (array,*) iscale
         read (nfparam,*) xlo, xhi, ylo, yhi
         linecount = linecount + 1
         read (nfparam,*) veladd, xadd, yadd, ymult
         if (xadd .ne. 0.) then
            veladd = 3.0d5*xadd/((xlo+xhi)/2.)
            xadd = 0.
         endif
         linecount = linecount + 1
         read (nfparam,*) smtype, fwhmgauss, vsini, limbdark, vmac,
     .                   fwhmloren
         linecount = linecount + 1
         call plotremember (1)


c  keyword 'pscale    ' is an old version of 'plotpars  '; if you use
c     it you will just get a warning message to switch to 'plotpars  '
      elseif (keyword .eq. 'pscale') then
         write(*,*) 'Warning: keyword changed to *plotpars*'
         write(*,*) 'for synthesis plotting parameters.'
         write(*,*) 'Here is the proper format:'
         write(*,*)
         write(*,*) 'plotpars 1'
         write(*,*) '<xlow> <xhi> <ylo> <yhi>'
         write(*,*) '<vshift> <lamshift> <obsadd> <obsmult>'
         write(*,*) '<smtype> <FWHMgauss> <vsini> <L.D.C.> <Vmacro> ',
     .     '<FWHMloren>'
         stop


c keyword 'trudamp     '  should moog use the detailed line damping for
c                         those transitions that have information stored in
c                         subroutine trudamp? (Default is *no*)
      elseif (keyword .eq. 'trudamp') then
         read (array,*) itru


c keyword 'veladjust   '  shoud moog try to do a cross-correlation between
c                         observed and synthetic spectra and use that to
c                         align the spectra better in wavelength
c                         (Default is *no*)
      elseif (keyword .eq. 'veladjust') then
         read (array,*) maxshift


c keyword 'units      ' controls the units in which moog 
c          outputs the final spectrum
c            0 = angs
c            1 = microns
c            2 = 1/cm
      elseif (keyword .eq. 'units') then
         read (array,*) iunits


c keyword 'iraf       ' allows the user to output a raw spectrum in 
c          a form suitable for IRAF's rtext input command
c            0 = don't do this, make output the normal way.
c            1 = make an IRAF-compatible output
      elseif (keyword .eq. 'iraf') then
         read (array,*) iraf


c  keyword 'scat       'allows the user to employ a source function
c          which has both scattering and absorption components     
c          0 = NO scattering
c          1 = scattering
      elseif (keyword .eq. 'scat') then
         read (array,*) scatopt


c  keyword 'flux/int  ' choses integrated flux or central intensity
c           0 = integrated flux calculations
c           1 = central intensity calculations
      elseif (keyword .eq. 'flux/int') then
         read (array,*) fluxintopt


c           0 = use the Unsold approximation, except multiply possibly
c                 by a factor read in for an individual line
c           1 = use the Unsold approximation multiplied by 6.3
c           2 = use the Unsold approximation multiplied by the
c                 factor recommended by the Blackwell group
      elseif (keyword .eq. 'damping') then
         read (array,*) dampingopt


c  keyword 'obspectrum' controls the file type of the observed spectrum
c           0 = no observed spectrum is to be input
c           1 = read a true FITS file with internal read statements
c          -1 = as if obspectrum = 1, but on a byte-swapping machine
c           2 = (not implemented yet)
c           3 = read a true Fits file with the FITSIO package
c           4 = (not implemented yet)
c           5 = read a special MONGO style (wavelength, flux pair) file
      elseif (keyword .eq. 'obspectrum') then
         read (array,*) specfileopt
         if (specfileopt .lt. 0) then
            byteswap = 1
            specfileopt = iabs(specfileopt)
         endif


c   keyword 'histogram' makes histogram plots of observed spectra if
c   histoyes = 1
      elseif (keyword .eq. 'histogram') then
         read (array,*) histoyes


c  keyword 'terminal  ' gives the sm plotting window type
c           smterm = a character string of the sm window type (see the 
c           appropriate sm manual for a list)
      elseif (keyword .eq. 'terminal') then
         read (array,*) smterm


c  keyword 'plot      ' decides whether or not to make a plot of results
c           0 = do not make a plot
c           For syntheses: 1 = plot only synthetic spectra
c                          2 = plot synthetic and observed spectra
c           For line analyses: # = the minimum number of lines of a 
c                                  species necessary to trigger a plot
c           For curves-of-growth: 1 = make plots
c           For flux curves: 1 = make plots
      elseif (keyword .eq. 'plot') then
         read (array,*) plotopt


c  keyword 'abundances' gives the changes to be applied to the abundances
c           # = the number of different syntheses to run
c               (the next line gives the different abundance factors
c               to use)
c  minimum error check:  numatomsyn must equal numisosyn or code will stop
      elseif (keyword .eq. 'abundances') then
         neq = 0
         numpecatom = 0
         numatomsyn = 0
         newnumpecatom = 0
         newnumatomsyn = 0
         ninetynineflag = 0
         do i=1,95
            pec(i) = 0
            newpec(i) = 0
            do j=1,5
               abfactor(j) = 0.
               pecabund(i,j) = 0.
               newpecabund(i,j) = 0.
            enddo
         enddo
         read (array,*) numpecatom,numatomsyn
         if (numisosyn .ne. 0) then
            if (numatomsyn .ne. numisosyn) then
               write (array,1002) numatomsyn, numisosyn
               call putasci (77,6)
               stop
            endif
         endif
         do l=1,numpecatom
            read (nfparam,*) jatom,(deltalogab(kk),kk=1,numatomsyn)
            linecount = linecount + 1
            if (jatom .eq. 99) then
               do kk=1,numatomsyn
                  abfactor (kk) = deltalogab(kk)
               enddo
            else
               do kk=1,numatomsyn
                  pecabund(jatom,kk) = deltalogab(kk)
               enddo
               pec(jatom) = 1
            endif
         enddo
         if (numpecatom.eq.1 .and. jatom.eq.99) ninetynineflag = 1


c keyword 'isotopes   ' gives the isotopes used in the line list and their
c                       abundance relative to the parent spiecies
c  minimum error check:  numatomsyn must equal numisosyn or code will stop
      elseif (keyword .eq. 'isotopes') then
         numiso = 0
         numisosyn = 0
         newnumiso = 0
         newnumisosyn = 0
         do i=1,20
            isotope(i) = 0.
            newisotope(i) = 0.
            do j=1,5
               isoabund(i,j) = 0.
               newisoabund(i,j) = 0.
            enddo
         enddo
         read (array,*) numiso,numisosyn
         if (numatomsyn .ne. 0) then
            if (numatomsyn .ne. numisosyn) then
               write (array,1002) numatomsyn, numisosyn
               call putasci (77,6)
               stop
            endif
         endif
         do  j=1,numiso
            read (nfparam,*) isotope(j),(isoabund(j,kk),kk=1,numisosyn)
            linecount = linecount + 1
         enddo


c  keyword 'lumratio' gives the ratio of the luminosity of two stars at a
c                     specific wavelength in a binary star system (used 
c                     only with driver "binary")
      elseif (keyword .eq. 'lumratio') then
         read (array,*) lumratio


c  keyword 'deltaradvel' gives the velocity difference between the stars
c                        binary star system (used only with driver "binary")
      elseif (keyword .eq. 'deltaradvel') then
         read (array,*) deltaradvel


c  keyword 'synlimits ' gives the wavelength parameters for syntheses;
c                       start and sstop are beginning and ending
c                       wavelengths, step is the step size in the
c                       syntheses, and delta is the wavelength range
c                       to either side of a synthesis point to consider
c                       for line opacity calculations
      elseif (keyword .eq. 'synlimits') then
         read (nfparam,*) start, sstop, step, delta
         linecount = linecount + 1


c  keyword 'fluxlimits' gives the wavelength parameters for flux curves;
c                       start and sstop are beginning and ending
c                       wavelengths, and step is the step size in the
c                        flux curve
      elseif (keyword .eq. 'fluxlimits') then
         read (nfparam,*) start, sstop, step
         linecount = linecount + 1


c  keyword 'blenlimits' gives the parameters for blended line abundance
c                       matches.  delwave is the wavelength offset 
c                       to the blue of first and to the red of the 
c                       last line in the blend to extend the syntheses; 
c                       step is the wavelength step size in the 
c                       computations; cogatom is the name of the
c                       element whose abundance should be varied
c                       to achieve an EW match with observations.
      elseif (keyword .eq. 'blenlimits') then
         read (nfparam,*) delwave, step, cogatom
         linecount = linecount + 1


c  keyword 'coglimits ' gives the log(W/lambda) limits for curves-of-growth
c                       rwlow and rwhigh are the beginning
c                       and ending points of the log(red.width) values,
c                       rwstep is the step in log(red.width),
c                       cogatom is the declaration of which element
c                       will have its abundance varied (necessary only
c                       for spectrum synthesis curves-of-growth,
c                       and wavestep is a forced (if desired) step size
c                       in wavelength along the line (this applies to
c                       single line computations only
      elseif (keyword .eq. 'coglimits') then
         read (nfparam,*) rwlow, rwhigh, rwstep, wavestep, cogatom
         linecount = linecount + 1


c  keyword 'limits    ' old limits format...tell the user to change the
c                       keyword and quit.

      elseif (keyword .eq. 'limits') then
      write(*,*) 'Warning: keyword changed to *synlimits*, *coglimits*'
      write(*,*) 'for Syntesis and COG calculations.'
      write(*,*) 'Here are the proper formats:'
      write(*,*) 
      write(*,*) 'synlimits <start> <stop> <step> <delta>'
      write(*,*) 'coglimits <rwlow> <rwhigh> <rwstep> <wavestep> ',
     .           '<cogatom>'
         stop


c   keyword of strong for lines which are to be considered for all of the 
c   synthesis
      elseif (keyword .eq. 'strong') then
         read (array,*) dostrong
     

c  keyword word of opacit which takes the continuus opacity and scales it
c  with the form of kaplam(i)= kaplam(i)*((factor*10000)/t(i))
c  in Opacit.f after it calulates the normal kaplam
c  if value is <= 0 then it does not do it
      elseif (keyword .eq. 'opacit') then
         read (array,*) fudge
     
      
c  any other keyword causes great grudge
      else
         write (array,1006) keyword
         call prinfo (5)
         stop

      endif
 

c  loop back to get another parameter
      go to 5
 

c  wrap things up with a few assignments
98    if (control .eq. 'gridsyn' .or. control .eq. 'gridplo' .or.
     .    control .eq. 'binary ' .or. control .eq. 'abandy ' .or.
     .    control .eq. 'gridsto') then
         control = 'gridend'
      endif


c  assign plotting window type; if no type has been given in the
c  parameter file, then ask for it
100   if (smterm .eq. ' ') then
         array = 'GIVE THE SM TERMINAL NAME : '
         nchar = 28
         call getasci (nchar,12)
         smterm = chinfo(1:nchar)
         ivstat = ivcleof(12,1)
      endif
      if (smterm.eq.'x11' .or. smterm.eq.'X11') then
         if    (control .eq. 'synth  ' .or.
     .          control .eq. 'synplot' .or.
     .          control .eq. 'isoplot' .or.
     .          control .eq. 'gridsyn' .or.
     .          control .eq. 'gridplo' .or.
     .          control .eq. 'doflux ' .or.
     .          control .eq. 'cogsyn ' .or.
     .          control .eq. 'cog    ' .or.
     .          control .eq. 'isotop ' .or.
     .          control .eq. 'binary ') then
             smterm = smt1
         else
             smterm = smt2
         endif
      endif


c  at the end of the parameter file, return after storing the starting
c  values of the synthesis parameters
c  wavlength units to Angstroms
c  exit normally
      if (iunits .eq. 1) then
         start = 1.d+4*start
         sstop  = 1.d+4*sstop
         step  = 1.d+4*step
         delta = 1.d+4*delta
      endif
      oldstart = start
      oldstop  = sstop
      oldstep  = step
      olddelta = delta
      return


c*****format statements
1001  format (a80)
1002  format ('# OF ABUNDANCE (',i1,') AND ISOTOPIC (',i1,')',
     .        ' SYNTHESES DO NOT AGREE!   I QUIT!       ')
1006  format ('THIS OPTION IS UNKNOWN TO MOOG: ', a10, ' I QUIT!')
1007  format (79(' '))


      end

c      character*(*) function JTRIM(str)
c      
c      character*(*) str
c      JTRIM = str(1:len_trim(str))
c      return
c      end

c      subroutine trim_spaces(str)
c****************************************************************************** 
c     This subroutine reads in the commands from the parameter file 
c****************************************************************************** 

c      implicit real*8 (a-h,o-z) 
c      character*80 str
c
c      i=1
c      do i=1,80
c          if str(i).eq.' ' then goto 10
c      enddo
c10    
c
c      end