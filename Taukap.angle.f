
      subroutine taukap (phi_angle, chi_angle)
c******************************************************************************
c     This routine calculates the line absorption coefficient and the line  
c     opacity at wavelength *wave* for all lines in the spectrum            
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      real*8 kappa(4, 4, 100), emission(4,4,100), kaptot
c      real*8 kapnu_I(100), kapnu_Q(100), kapnu_V(100), kapnu_U(100),
c     .       new_voigt, new_fv, voigt_x, voigt_y, gam_L, gam_D,
c     .       zetnu_Q(100), zetnu_V(100), zetnu_U(100), sqrtpi, old_voigt
c      complex Z, hui, w4

c*****compute the total line opacity at each depth    
      matX(1,1) = 1.0
      matX(1,2) = 0.0
      matX(1,3) = 0.0
      matX(1,4) = 0.0
      matX(2,1) = 0.0
      matX(2,2) = 1.0
      matX(2,3) = 0.0
      matX(2,4) = 0.0
      matX(3,1) = 0.0
      matX(3,2) = 0.0
      matX(3,3) = 1.0
      matX(3,4) = 0.0
      matX(4,1) = 0.0
      matX(4,2) = 0.0
      matX(4,3) = 0.0
      matX(4,4) = 1.0
      
      do i=1,ntau
         phi_I=(phi_opacity(i,2)*sin(phi_angle)**2.0+
     .                (phi_opacity(i,1)+phi_opacity(i,3))*(1.0+
     .                cos(phi_angle)**2.0)/2.0)/2.0
         phi_Q=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(1,3))
     .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
         phi_U=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(1,3))
     .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
         phi_V=(phi_opacity(i,1)-phi_opacity(i,3))*cos(phi_angle)/2.0
         psi_Q=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(1,3))
     .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
         psi_U=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(1,3))
     .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
         psi_V=(psi_opacity(i,1)-psi_opacity(i,3))*cos(phi_angle)/2.0

         kaptot = kaplam(i) + phi_I
         
         kappa(1,1,i)=0.0
         kappa(2,1,i)=phi_Q/kaptot
         kappa(3,1,i)=phi_U/kaptot
         kappa(4,1,i)=phi_V/kaptot
         kappa(1,2,i)=phi_Q/kaptot
         kappa(2,2,i)=0.0
         kappa(3,2,i)=psi_V/kaptot
         kappa(4,2,i)=(-1.0*psi_U)/kaptot
         kappa(1,3,i)=phi_U/kaptot
         kappa(2,3,i)=(-1.0*psi_V)/kaptot
         kappa(3,3,i)=0.0
         kappa(4,3,i)=psi_Q/kaptot
         kappa(1,4,i)=phi_V/kaptot
         kappa(2,4,i)=psi_U/kaptot
         kappa(3,4,i)=(-1.0*psi_Q)/kaptot
         kappa(4,4,i)=0.0

         emission(1,i)=source(i)
         emission(2,i)=source(i)*phi_Q/kaptot
         emission(3,i)=source(i)*phi_U/kaptot
         emission(4,i)=source(i)*phi_V/kaptot

c         eta_I(i) = kapnu_I(i)/(kaplam(i))
c         eta_Q(i) = kapnu_Q(i)/(kaplam(i))
c         eta_V(i) = kapnu_V(i)/(kaplam(i))
c         zet_Q(i) = zetnu_Q(i)/(kaplam(i))
c         zet_V(i) = zetnu_V(i)/(kaplam(i))
      enddo
      do i=1,ntau-2
         dtau = tau(i+1) - tau(i)
         etau = 2.71828183**(-dtau)

         alph = 1.0-etau
         bet = (1.0-(1.0+dtau)*etau)/dtau
         call daxpy(16,(alph-bet),kappa(i),1,matX,1)
         call daxpy(16,(alph-bet),kappa(i),1,matY,1)

      enddo
         
      matY = (etau*ones -bet*kappa(i+1)


      return                                              
321   format (f11.3, e11.3, f11.1, 5e11.3)
c322   format (f11.3, 4e11.3)
      end                                                


