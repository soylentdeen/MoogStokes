
      procedure AppendStokes(phi_angle, chi_angle, cell_a, Stokes)
c*****************************************************************
c     Adds an emergent spectrum (calculated at angles phi, chi) to
c     the resultant spectrum.  Calculated at a single wavelength.
c*****************************************************************

      implicit (a-h,o-z)
      include "Atmos.com"
      real*8 area, proj_area, limb_darkening
      real*8 surface_area, mu

      if ((1.0/wave).lt.2.4) then
          alpha = -0.023 + 0.292/wave
      else
          alpha = -0.507 + 0.441/wave
      endif

      mu = cos(phi_angle)

      limb_darkening = (1.0-(1.0-mu**alpha))
