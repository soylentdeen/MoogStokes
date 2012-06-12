
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
      real*8 matX(4,4), matY(4,4), ones(4,4), matZ(4)
      real*8 emiss_interp(4,3), kappa_interp(4,4,2)
      real*8 tau_interp(3), tau_interp_c(3), logtau
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V
      real*8 h1, h2, dx
      real*8 matS1(4), matS2(4)
      integer emiss_order(3), kappa_order(2)
      integer INFO, IPIV(4), LDA, LWORK
      parameter (LDA=4, LWORK=64*LDA)
      double precision WORK(LWORK)

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

c*****  zdepth is the physical depth scale
      do i=1,ntau
         if (i.eq.1) then
             zdepth(i) = 0.0
         else
             h1 = 1.0/(kapref(i-1))
             h2 = 1.0/(kapref(i))
             dtau = tauref(i)-tauref(i-1)
             zdepth(i) = zdepth(i-1)+(h1+h2)/2.0*dtau
         endif
      enddo

      phi_angle = dble(0.0)
      chi_angle = dble(0.0)
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

      Stokes(1) = emission(1,ntau)
      Stokes(2) = 0.0
      Stokes(3) = 0.0
      Stokes(4) = 0.0
      continuum = Stokes(1)

      delta_tau = -0.05
c      write (*,*) Stokes(1), emission(1, ntau)
      call dcopy(4, emission(:,ntau), 1, emiss_interp(:,1), 1)
c      tau_interp(1) = zdepth(ntau)*kaptot(ntau)
c      tau_interp_c(1) = zdepth(ntau)*kaplam(ntau)
      tau_interp(1) = tauref(ntau)*kaptot(ntau)/kapref(ntau)
      tau_interp_c(1) = tauref(ntau)*kaplam(ntau)/kapref(ntau)

c      write (*,*) emiss_interp
      call interp_opacities(log10(tauref(ntau)),delta_tau,kappa_interp,
     .        1, emiss_interp, 2, tau_interp, tau_interp_c)
      kappa_order(1) = 1
      kappa_order(2) = 2
      emiss_order(1) = 1
      emiss_order(2) = 2
      emiss_order(3) = 3
      do logtau=log10(tauref(ntau))+delta_tau,
     .              log10(tauref(1))-delta_tau,delta_tau
         call interp_opacities(logtau, delta_tau, kappa_interp,
     .        kappa_order(2), emiss_interp, emiss_order(3), tau_interp,
     .        tau_interp_c)
         
         dtau = (tau_interp(emiss_order(1))-tau_interp(emiss_order(2)))
     .              *cos(viewing_angle)
         etau = 2.71828183**(-dtau)

         alph = 1.0-etau
         bet =(1.0-(1.0+dtau)*etau)/dtau

         call dcopy(16,ones, 1, matX, 1)
         call dcopy(16,ones, 1, matY, 1)
         call daxpy(16,dble(alph-bet),kappa_interp(:,:,kappa_order(2)),
     .              1,matX,1)
         call dscal(16,etau, matY,1)
         call daxpy(16,dble(-1.0*bet),kappa_interp(:,:,kappa_order(1)),
     .              1,matY,1)

         x = 1.0 - etau
         y = dtau - x
         z = dtau**2.0 - 2.0 * y
         dtau_i=(tau_interp(emiss_order(2))-tau_interp(emiss_order(3)))
     .            *cos(viewing_angle)
         alph = (z -dtau*y)/((dtau + dtau_i)*dtau_i)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))

         call dcopy(4, emiss_interp(:,emiss_order(3)), 1, matS1, 1)
         call dcopy(4, emiss_interp(:,emiss_order(2)), 1, matS2, 1)
         call dcopy(4, emiss_interp(:,emiss_order(1)), 1, matZ, 1)
         call dscal(4, alph, matS1, 1)
         call dscal(4, bet, matS2, 1)
         call dscal(4, gam, matZ, 1)
         call daxpy(4, dble(1.0), matS1, 1, matS2, 1)
         call daxpy(4, dble(1.0), matS2, 1, matZ, 1)

c****      calculate the RHS of the equation.  Store in matZ
         call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)

c****      Solve the system of differential equations
         call dgesv(4,1,matX,4,IPIV,matZ,4,INFO)

         call dcopy(4, matZ, 1, Stokes, 1)

c****     Now do the same thing for the continuum
         dtau=(tau_interp_c(emiss_order(1))-
     .         tau_interp_c(emiss_order(2)))*cos(viewing_angle)
         etau = 2.71828183**(-dtau)
         x = 1.0 - etau
         y = dtau - x
         z = dtau**2.0 - 2.0 * y
         dtau_i = (tau_interp_c(emiss_order(2))-
     .             tau_interp_c(emiss_order(3)))*cos(viewing_angle)
         alph = (z -dtau*y)/((dtau + dtau_i)*dtau_i)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))
         continuum=etau*continuum+alph*emiss_interp(1,emiss_order(3))
     .             +bet*emiss_interp(1,emiss_order(2))
     .             +gam*emiss_interp(1,emiss_order(1))

         if (kappa_order(1).eq.1)then
             kappa_order(1) = 2
             kappa_order(2) = 1
         else
             kappa_order(1) = 1
             kappa_order(2) = 2
         endif
         if (emiss_order(1).eq.1) then
             emiss_order(1) = 2
             emiss_order(2) = 3
             emiss_order(3) = 1
         elseif (emiss_order(1).eq.2) then
             emiss_order(1) = 3
             emiss_order(2) = 1
             emiss_order(3) = 2
         else
             emiss_order(1) = 1
             emiss_order(2) = 2
             emiss_order(3) = 3
         endif

c         continuum=etau*continuum+alph+bet+gam
      enddo

c      write (*,*) Stokes
      return
      end

      subroutine interp_opacities(logtau, delta_tau, k_interp,
     .       k_ord, e_interp, e_ord, tau_interp, tau_interp_c)
c**********************************************************************
c     interp_opacities interpolates the following quantities relevant
c        to the DELO method:
c      kapp_interp(kap_ord) - kappa matrix interpolated at
c                      tau = 10.0**logtau
c      emiss_interp(emiss_ord) - emission matrix interpolated at
c                      tau = 10.0**(logtau+dtau)
c      tau_interp(emiss_ord) - total line opacity calculated at:
c                      tau = 10.0**(logtau+dtau) (reference tau)
c      tau_interp_c(emiss_ord) - continuum opacity calculated at:
c                      tau = 10.0**(logtau+dtau) (reference tau)
c***********************************************************************
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Stokes.com'
      real*8 slope, k_interp(4,4,2), e_interp(4,3)
      real*8 tau_interp(3), tau_interp_c(3), logtau, delta_tau
      integer k_ord, e_ord

      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.10.0**logtau) then
             goto 10
         endif
      enddo
10    do j=1,4
         do k=1,4
            slope=(kappa(j,k,i+1)-kappa(j,k,i))/(tauref(i+1)
     .             -tauref(i))
            k_interp(j,k,k_ord)=kappa(j,k,i)+slope*(10.0**logtau
     .             -tauref(i))
         enddo
      enddo

      line_kappa = 0.0
      cont_kappa = 0.0
      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.10.0**(logtau+delta_tau)) then
             goto 20
         endif
         dz = zdepth(i+1)-zdepth(i)
         line_kappa = line_kappa + (kaptot(i)+kaptot(i+1))/2.0*dz
         cont_kappa = cont_kappa + (kaplam(i)+kaplam(i+1))/2.0*dz
      enddo
20    do j=1,4
         slope=(emission(j,i+1)-emission(j,i))/(tauref(i+1)-tauref(i))
         e_interp(j,e_ord) = emission(j,i)+slope*(10.0**(logtau+
     .       delta_tau)-tauref(i))
      enddo
c      dz = zdepth(i+1)-zdepth(i)
c      slope = (kaptot(i+1)+kaptot(i))/2.0*dz/(tauref(i+1)-tauref(i))
c      write (*,*) 'slope: Zdepth method : ', i, slope
      slope = (tauref(i+1)*kaptot(i+1)/kapref(i+1)-
     .         tauref(i)*kaptot(i)/kapref(i))/(tauref(i+1)-tauref(i))
c      write (*,*) tauref(i+1)/kapref(i+1)
c      write (*,*) 'slope: tauref method : ', slope_a, slope/slope_a
c      tau_interp(e_ord) = line_kappa+
      tau_interp(e_ord) = tauref(i)*kaptot(i)/kapref(i)+
     .        slope*(10.0**(logtau+delta_tau)-tauref(i))

c      read (*,*)

      slope = (tauref(i+1)*kaplam(i+1)/kapref(i+1)-
     .         tauref(i)*kaplam(i)/kapref(i))/(tauref(i+1)-tauref(i))
c      slope = (zdepth(i+1)*kaplam(i+1)-
c     .         zdepth(i)*kaplam(i))/(tauref(i)-tauref(i))
c      tau_interp_c(e_ord) = cont_kappa+
      tau_interp_c(e_ord) = tauref(i)*kaplam(i)/kapref(i)+
     .        slope*(10.0**(logtau+delta_tau)-tauref(i))
      
      return
      end

c      do i=ntau-1,2,-1
c         dz = (tauref(i+1)-tauref(i))/kapref(i)
c         dtau = dz*kaptot(i)*cos(viewing_angle)
c         etau = 2.71828183**(-dtau)
c
c         alph = 1.0-etau
c         bet =(1.0-(1.0+dtau)*etau)/dtau
c
c         call dcopy(16,ones, 1, matX, 1)
c         call dcopy(16,ones, 1, matY, 1)
c         call daxpy(16,dble(alph-bet),kappa(:,:,i),1,matX,1)
c         call dscal(16,etau, matY,1)
c         call daxpy(16,dble(-1.0*bet),kappa(:,:,i+1),1,matY,1)
c
c         x = 1.0 - etau
c         y = dtau - x
c         z = dtau**2.0 - 2 * y
c         dz = (tauref(i)-tauref(i-1))/kapref(i-1)
c         dtau_i = dz*kaptot(i-1)*cos(viewing_angle)
c         alph = (z -dtau*y)/((dtau + dtau_i)*dtau_i)
c         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
c         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))
c
c         call dcopy(4, emission(:,i-1), 1, matS1, 1)
c         call dcopy(4, emission(:,i), 1, matS2, 1)
c         call dcopy(4, emission(:,i+1), 1, matZ, 1)
c         call dscal(4, alph, matS1, 1)
c         call dscal(4, bet, matS2, 1)
c         call dscal(4, gam, matZ, 1)
c         call daxpy(4, dble(1.0), matS1, 1, matS2, 1)
c         call daxpy(4, dble(1.0), matS2, 1, matZ, 1)
c
cc****     Calculate the right hand side of the equation.  Store in matZ
c         call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)
c
cc****     Solve the system of differential equations.
c         call dgesv(4,1,matX,4,IPIV, matZ,4,INFO)
c
cc****     Now do the same thing for the continuum
c         dz = (tauref(i+1)-tauref(i))/kapref(i)
c         dtau = dz*kaplam(i)*cos(viewing_angle)
c         etau = 2.71828183**(-dtau)
c         x = 1.0 - etau
c         y = dtau - x
c         z = dtau**2.0 - 2 * y
c         dz = (tauref(i)-tauref(i-1))/kapref(i-1)
c         dtau_i = dz*kaplam(i-1)*cos(viewing_angle)
c         alph = (z -dtau*y)/((dtau + dtau_i)*dtau_i)
c         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
c         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))
c         continuum=etau*continuum+alph*emission(1,i-1)
c     .             +bet*emission(1,i)
c     .             +gam*emission(1,i+1)
c
c         call dcopy(4, matZ, 1, Stokes, 1)
c    
c      enddo
c
c      write (*,*) continuum, Stokes
c
c      return
c      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
