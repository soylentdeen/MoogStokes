
      subroutine taukap (phi_angle, chi_angle)
c******************************************************************************
c     This routine calculates the line absorption coefficient and the line  
c     opacity at wavelength *wave* for all lines in the spectrum            
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      real*8 kappa(4, 4, 100), emission(4,4,100), kaptot, ones(4,4)
c      real*8 kapnu_I(100), kapnu_Q(100), kapnu_V(100), kapnu_U(100),
c     .       new_voigt, new_fv, voigt_x, voigt_y, gam_L, gam_D,
c     .       zetnu_Q(100), zetnu_V(100), zetnu_U(100), sqrtpi, old_voigt
c      complex Z, hui, w4
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V, dtau, etau
      real*8 matX(4,4), matY(4,4), matZ(4), old_I(4), matS1(4), matS2(4)
      real*8 alph, bet, gam, x, y, z, IPIV(4), INFO

c*****compute the total line opacity at each depth    
      ones(1,1) = 1.0
      ones(1,2) = 0.0
      ones(1,3) = 0.0
      ones(1,4) = 0.0
      ones(2,1) = 0.0
      ones(2,2) = 1.0
      ones(2,3) = 0.0
      ones(2,4) = 0.0
      ones(3,1) = 0.0
      ones(3,2) = 0.0
      ones(3,3) = 1.0
      ones(3,4) = 0.0
      ones(4,1) = 0.0
      ones(4,2) = 0.0
      ones(4,3) = 0.0
      ones(4,4) = 1.0
      
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

         source = Planck(t(i))

         emission(1,i)=source
         emission(2,i)=source*phi_Q/kaptot
         emission(3,i)=source*phi_U/kaptot
         emission(4,i)=source*phi_V/kaptot

c         eta_I(i) = kapnu_I(i)/(kaplam(i))
c         eta_Q(i) = kapnu_Q(i)/(kaplam(i))
c         eta_V(i) = kapnu_V(i)/(kaplam(i))
c         zet_Q(i) = zetnu_Q(i)/(kaplam(i))
c         zet_V(i) = zetnu_V(i)/(kaplam(i))
      enddo
      old_I(1) = Planck(t(ntau))
      old_I(2) = 0.0
      old_I(3) = 0.0
      old_I(4) = 0.0
      do i=ntau-2,1,-1
         dtau = tau(i+1) - tau(i)
         etau = 2.71828183**(-dtau)

         alph = 1.0-etau
         bet = (1.0-(1.0+dtau)*etau)/dtau

         call dcopy(16,ones, 1, matX, 1)
         call dcopy(16,ones, 1, matY, 1)
         call daxpy(16,(alph-bet),kappa((i-1)*16+1,i*16+1),1,matX,1)
         call dscal(16,etau, matY,1)
         call daxpy(16,(-1.0*bet),kappa(i*16+1,(i+1)*16+1),1,matY,1)
         
         x = 1 - etau
         y = dtau - x
         z = dtau**2.0 - 2 * y
         dtau_i = tau(i+2) - tau(i+1)
         alph = (z - dtau_i*y)/((dtau + dtau_i)*dtau)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau + 2*dtau_i)*y)/(dtau_i*(dtau+dtau_i))

         call dcopy(4, emission((i-1)*4+1,i*4+1), 1, matS1, 1)
         call dcopy(4, emission(i*4+1,(i+1)*4+1), 1, matS2, 1)
         call dcopy(4, emission((i+1)*4+1,(i+2)*4+1), 1, matZ, 1)
         call dscal(4, alph, matS1, 1)
         call dscal(4, bet, matS2, 1)
         call dscal(4, gam, matZ, 1)
         call daxpy(4, 1.0, matS1, 1, matS2, 1)
         call daxpy(4, 1.0, matS2, 1, matZ, 1)

c     Calculate the right hand side of the equation.  Store in matZ
         call dgemv('N',4,4,1.0,matY,4,old_I,1,1.0,matZ,1)

c     Solve the system of differential equations.
         call dgesv(4,1,matX,4,IPIV, matZ,4,INFO)

         write (*,*) matZ
         call dcopy(4, matZ, 1, old_I, 1)

      enddo
         
      matY = (etau*ones -bet*kappa(i+1)

      return
321   format (f11.3, e11.3, f11.1, 5e11.3)
c322   format (f11.3, 4e11.3)
      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      Planck =((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end

