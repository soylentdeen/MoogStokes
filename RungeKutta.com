
c******************************************************************************
c     this common block has variables related to stellar and observed geometry
c     
c     Stellar Rotation Axis Angles:
c==================================
c     inclination = angle of inclination (measured from pole-on)
c     theta_clock = clocking angle measured clockwise from North
c     
c     Observation Angles:
c==================================
c     polar_r = Not an angle, but distance from stellar disk center in polar
c               coordinates.
c     polar_phi = polar angle, measured from stellar disk equator
c     
c     Magnetic Field Angles: - depend on observation angles, stellar rotation
c     axis and the geometry of the magnetic field
c==================================
c     B_phi
c     B_chi
c
c     Geometrical Angles:
c==================================
c     Projected Area
c     Limb Darkening Angle
c******************************************************************************

      real*8 viewing_angle

      common/rungekutta/viewing_angle
