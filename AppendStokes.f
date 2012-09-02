
      subroutine AppendStokes(cell_a)
c*****************************************************************
c     Adds an emergent spectrum (calculated at angles phi, chi) to
c     the resultant spectrum.  Calculated at a single wavelength.
c*****************************************************************

      implicit real*8 (a-h,o-z)
      include "Atmos.com"
      include "Stokes.com"
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
      
      projected_area = cell_a*dotproduct/(4.0*3.14159262)

      weight = limb_darkening*projected_area

      Stokes_I = Stokes_I + Stokes(1)/continuum*weight
      Stokes_Q = Stokes_Q + Stokes(2)/continuum*weight
      Stokes_U = Stokes_U + Stokes(3)/continuum*weight
      Stokes_V = Stokes_V + Stokes(4)/continuum*weight
      

      total_weight = total_weight + weight

      end
