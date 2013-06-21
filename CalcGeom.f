
      subroutine calcGeom (az, long, phi_ang, chi_ang, mu)

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
      return
      end
