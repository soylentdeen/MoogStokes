
c******************************************************************************
c     this common block has variables related to the lines.  Most
c     input quantities typically have single dimensions, while the 
c     things that are computed for each line at each atmosphere level
c     have double dimensions.  The variables "a", "dopp", and 
c     "kapnu0" are often over-written with plotting data,
c     so leave them alone or suffer unspeakable programming tortures.
c******************************************************************************

      real*8       deltamj(2500),phi_opacity(100, 3),psi_opacity(100,3),
     .             Stokes(4), continuum, wavelength, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V

      common/stokesparams/deltamj, phi_opacity, psi_opacity,
     .             Stokes, continuum, wavelength, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V

