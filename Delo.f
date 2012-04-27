
      subroutine delo
c**********************************************************************
c     This routine performs the DELO integration routine through the 
c     photosphere
c**********************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      include 'Stokes.com'
      include 'Angles.com'
      real*8 kappa(4,4,100), emission(4,100), kaptot(100), ones(4,4)
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V, dtau, etau
      real*8 matX(4,4), matY(4,4), matZ(4), matS1(4), matS2(4)
      real*8 alph, bet, gam, x, y, z, IPIV(4), INFO, dtau_i, dz

c*****Sets up constants

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

c***** For each layer in the atmosphere, calculate each element of the
c      opacity matrix and emission vector for the DELO algorithm
      do i=1,ntau
         phi_I=(phi_opacity(i,2)*sin(phi_angle)**2.0+
     .                (phi_opacity(i,1)+phi_opacity(i,3))*(1.0+
     .                cos(phi_angle)**2.0)/2.0)/2.0
         phi_Q=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
         phi_U=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
         phi_V=(phi_opacity(i,1)-phi_opacity(i,3))*cos(phi_angle)/2.0
         psi_Q=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
         psi_U=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
         psi_V=(psi_opacity(i,1)-psi_opacity(i,3))*cos(phi_angle)/2.0

c*****  The total opacity (line+continuum)
         kaptot(i) = kaplam(i) + phi_I

c*****  Assemble the Opacity matrix (K')
         kappa(1,1,i)=0.0
         kappa(1,2,i)=phi_Q/kaptot(i)
         kappa(1,3,i)=phi_U/kaptot(i)
         kappa(1,4,i)=phi_V/kaptot(i)
         kappa(2,1,i)=phi_Q/kaptot(i)
         kappa(2,2,i)=0.0
         kappa(2,3,i)=psi_V/kaptot(i)
         kappa(2,4,i)=(-1.0*psi_U)/kaptot(i)
         kappa(3,1,i)=phi_U/kaptot(i)
         kappa(3,2,i)=(-1.0*psi_V)/kaptot(i)
         kappa(3,3,i)=0.0
         kappa(3,4,i)=psi_Q/kaptot(i)
         kappa(4,1,i)=phi_V/kaptot(i)
         kappa(4,2,i)=psi_U/kaptot(i)
         kappa(4,3,i)=(-1.0*psi_Q)/kaptot(i)
         kappa(4,4,i)=0.0

c*****  Assumes LTE for the Source Function
         source = Planck(t(i))

c*****  Assembles the Emission matrix (J')
         emission(1,i)=source
         emission(2,i)=source*phi_Q/kaptot(i)
         emission(3,i)=source*phi_U/kaptot(i)
         emission(4,i)=source*phi_V/kaptot(i)

      enddo

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm

      Stokes(1) = source
      Stokes(2) = 0.0
      Stokes(3) = 0.0
      Stokes(4) = 0.0
      continuum = source
      do i=ntau-2,1,-1
         dz = -(tauref(i+1)-tauref(i))/kapref(i)
         dtau = -dz*kaptot(i)*cos(viewing_angle)
         etau = 2.71828183**(-dtau)

         alph = 1.0-etau
         bet =(1.0-(1.0+dtau)*etau)/dtau

         call dcopy(16,ones, 1, matX, 1)
         call dcopy(16,ones, 1, matY, 1)
         call daxpy(16,dble(alph-bet),kappa(:,:,i),1,matX,1)
         call dscal(16,etau, matY,1)
         call daxpy(16,dble(-1.0*bet),kappa(:,:,i+1),1,matY,1)

         x = 1 - etau
         y = dtau - x
         z = dtau**2.0 - 2 * y
         dz = -(tauref(i+2)-tauref(i))/kapref(i+1)
         dtau_i = -dz*kaptot(i+1)*cos(viewing_angle)
         alph = (z -dtau_i*y)/((dtau + dtau_i)*dtau)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau + 2*dtau_i)*y)/(dtau_i*(dtau+dtau_i))

         call dcopy(4, emission(:,i), 1, matS1, 1)
         call dcopy(4, emission(:,i+1), 1, matS2, 1)
         call dcopy(4, emission(:,i+2), 1, matZ, 1)
         call dscal(4, alph, matS1, 1)
         call dscal(4, bet, matS2, 1)
         call dscal(4, gam, matZ, 1)
         call daxpy(4, dble(1.0), matS1, 1, matS2, 1)
         call daxpy(4, dble(1.0), matS2, 1, matZ, 1)

c****     Calculate the right hand side of the equation.  Store in matZ
         call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)

c****     Solve the system of differential equations.
         call dgesv(4,1,matX,4,IPIV, matZ,4,INFO)

c****     Now do the same thing for the continuum
         dz = -(tauref(i+1)-tauref(i))/kapref(i)
         dtau = -dz*kaplam(i)*cos(viewing_angle)
         etau = 2.71828183**(-dtau)
         x = 1 - etau
         y = dtau - x
         z = dtau**2.0 - 2 * y
         dz = -(tauref(i+2)-tauref(i))/kapref(i+1)
         dtau_i = -dz*kaplam(i+1)*cos(viewing_angle)
         alph = (z -dtau_i*y)/((dtau + dtau_i)*dtau)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau + 2*dtau_i)*y)/(dtau_i*(dtau+dtau_i))
         continuum=etau*continuum+alph*emission(1,i)
     .             +bet*emission(1,i+1)
     .             +gam*emission(1,i+2)

         call dcopy(4, matZ, 1, Stokes, 1)
    
      enddo

c      write (*,*) continuum, Stokes

      return
      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
