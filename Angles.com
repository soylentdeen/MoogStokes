
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

      real*8  inclination, clocking, position_angle, azimuth(1000),
     .      longitude(1000), B_sph(3), B_xyz(3), T_rho(3,3), T_eta(3,3),
     .      T_i(3,3), rotation_matrix(3,3), zeros(3,3),
     .      chi_angle(1000), phi_angle(1000), mus(1000)
      integer nrings, ncells

      common/angles/inclination, clocking, position_angle, azimuth,
     .              longitude, B_sph, B_xyz, T_rho, T_eta,
     .              T_i, rotation_matrix, zeros,
     .              chi_angle, phi_angle, mus,
     .              nrings, ncells
