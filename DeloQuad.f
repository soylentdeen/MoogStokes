
      subroutine traceStokes (phi_ang, chi_ang, mu)
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
      real*8 matX(4,4), matY(4,4), ones(4,4), matZ(4), bk(4,4)
      real*8 emiss_interp(4,3), kappa_interp(4,4,2)
      real*8 tau_interp(3), tau_interp_c(3), logtau
      real*8 phi_I, phi_U, phi_V, psi_Q, psi_U, psi_V
      real*8 h1, h2, dtau, etau, alph, bet, gam
      real*8 matS1(4), matS2(4), bgfl
      real*8 blah, k1, k2, z1, z2, deltaz(100), qmax
      real*8 phi_ang, chi_ang, mu, midpoint
      integer emiss_order(3), kappa_order(2)
      integer INFO, IPIV(4), LDA, LWORK
      parameter (LDA=4, LWORK=64*LDA)
      double precision WORK(LWORK)

c*****Sets up constants

      bgfl = dble(2e30)

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
      call spline(xref, kapref, ntau, bgfl, bgfl, dkref)
      do i=1,ntau
         if (i.eq.1) then
             zdepth(i) = 0.0
         else
             h1 = 1.0/kapref(i-1)
             midpoint = log10((tauref(i)+tauref(i-1))/2.0)
             call splint(xref, kapref, dkref, ntau, midpoint, k1)
             h2 = 4.0/k1
             h3 = 1.0/kapref(i)
             dtau = (tauref(i)-tauref(i-1))/(6.0*mu)
             zdepth(i) = zdepth(i-1)+(h1+h2+h3)*dtau
         endif
      enddo


c***** For each layer in the atmosphere, calculate each element of the
c      opacity matrix and emission vector for the DELO algorithm
      do i=1,ntau
         phi_I=(phi_opacity(i,2)*sin(phi_ang)**2.0+
     .                (phi_opacity(i,1)+phi_opacity(i,3))*(1.0+
     .                cos(phi_ang)**2.0)/2.0)/2.0
         phi_Q=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*
     .                      cos(2.0*chi_ang)/2.0
         phi_U=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*
     .                      sin(2.0*chi_ang)/2.0
         phi_V=(phi_opacity(i,1)-phi_opacity(i,3))*
     .                cos(phi_ang)/2.0
         psi_Q=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*
     .                cos(2.0*chi_ang)/2.0
         psi_U=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*
     .                sin(2.0*chi_ang)/2.0
         psi_V=(psi_opacity(i,1)-psi_opacity(i,3))*
     .                cos(phi_ang)/2.0

c*****  The total opacity (line+continuum)
         kaptot(i) = (kaplam(i) + phi_I)

C*****   Assemble the elements of the opacity matrix (K')
         phiQ(i) = phi_Q/kaptot(i)
         phiU(i) = phi_U/kaptot(i)
         phiV(i) = phi_V/kaptot(i)
         psiQ(i) = psi_Q/kaptot(i)
         psiU(i) = psi_U/kaptot(i)
         psiV(i) = psi_V/kaptot(i)

c*****  Assumes LTE for the Source Function
         source = Planck(t(i))

c*****  Assembles the Emission matrix (J')
         emission(1,i)=source
         emission(2,i)=source*phi_Q/kaptot(i)
         emission(3,i)=source*phi_U/kaptot(i)
         emission(4,i)=source*phi_V/kaptot(i)
      enddo

      call spline(xref, kaplam, ntau, bgfl, bgfl, dklam)
      call spline(xref, kaptot, ntau, bgfl, bgfl, dktot)
      call spline(xref, zdepth, ntau, bgfl, bgfl, deltaz)

      do i=1,ntau
         if (i.eq.1) then
            tlam(i) = kaplam(i)/kapref(i)*tauref(i)
            tautot(i) = kaptot(i)/kapref(i)*tauref(i)
         else
            dz = (zdepth(i)-zdepth(i-1))/6.0
            h1 = kaptot(i-1)
            midpoint = log10((tauref(i)+tauref(i-1))/2.0)
            call splint(xref, kaptot, dktot, ntau, midpoint, h2)
            h3 = kaptot(i)
            dtautot = dz*(h1+4.0*h2+h3)
            h1 = kaplam(i-1)
            call splint(xref, kaplam, dklam, ntau, midpoint, h2)
            h3 = kaplam(i)
            dtaulam = dz*(h1+4.0*h2+h3)
            tautot(i) = tautot(i-1)+dtautot
            tlam(i) = tlam(i-1)+dtaulam
         endif
      enddo

      call spline(xref, phiQ, ntau, bgfl, bgfl, dphiQ)
      call spline(xref, phiU, ntau, bgfl, bgfl, dphiU)
      call spline(xref, phiV, ntau, bgfl, bgfl, dphiV)
      call spline(xref, psiQ, ntau, bgfl, bgfl, dpsiQ)
      call spline(xref, psiU, ntau, bgfl, bgfl, dpsiU)
      call spline(xref, psiV, ntau, bgfl, bgfl, dpsiV)

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm

      dbdz = (Planck(t(1)) - Planck(t(2)) )/(zdepth(1)-zdepth(2))
      bk(1,1) = kaptot(1)
      bk(1,2) = phiQ(1)
      bk(1,3) = phiU(1)
      bk(1,4) = phiV(1)
      bk(2,1) = phiQ(1)
      bk(2,2) = kaptot(1)
      bk(2,3) = psiV(1)
      bk(2,4) = -psiU(1)
      bk(3,1) = phiU(1)
      bk(3,2) = -psiV(1)
      bk(3,3) = kaptot(1)
      bk(3,4) = psiQ(1)
      bk(4,1) = phiV(1)
      bk(4,2) = psiU(1)
      bk(4,3) = -psiQ(1)
      bk(4,4) = kaptot(1)

      CALL DGETRF( 4, 4, bk, 4, IPIV, INFO )
      CALL DGETRI(4, bk, 4, IPIV, WORK, LWORK, INFO)

      Stokes(1) = Planck(t(1)) - dbdz*bk(1,1)
      Stokes(2) = -dbdz*bk(2,1)
      Stokes(3) = -dbdz*bk(3,1)
      Stokes(4) = -dbdz*bk(4,1)
      continuum = Stokes(1)

      delta_tau = -0.05
      call dcopy(4, emission(:,ntau), 1, emiss_interp(:,1), 1)
      tau_interp(1) = tautot(ntau)
      tau_interp_c(1) = tlam(ntau)

      call interp_opacities(log10(tauref(ntau)),
     .   kappa_interp, 1, emiss_interp, 2, tau_interp,tau_interp_c,
     .   mu, delta_tau)
      kappa_order(1) = 1
      kappa_order(2) = 2
      emiss_order(1) = 1
      emiss_order(2) = 2
      emiss_order(3) = 3
      do logtau=log10(tauref(ntau))+2.0*delta_tau,
     .              log10(tauref(1)),delta_tau
         call interp_opacities(logtau, kappa_interp,
     .        kappa_order(2), emiss_interp, emiss_order(3), tau_interp,
     .        tau_interp_c, mu, delta_tau)
         dtau = (tau_interp(emiss_order(1))-tau_interp(emiss_order(2)))
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


         x = 1.0-etau
         y = dtau - x
         z = dtau**2.0 - 2.0 * y
         dtau_i =(tau_interp(emiss_order(2))-tau_interp(emiss_order(3)))
c         alph = (z - dtau*y)/((dtau + dtau_i)*dtau_i)
c         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
c         gam = x+(z-(dtau_i+2*dtau)*y)/(dtau*(dtau+dtau_i))
         alph = (z - dtau_i*y)/((dtau + dtau_i)*dtau)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(2*dtau_i+dtau)*y)/(dtau_i*(dtau+dtau_i))

         call dcopy(4, emiss_interp(:,emiss_order(3)), 1, matS1, 1)
         call dcopy(4, emiss_interp(:,emiss_order(2)), 1, matS2, 1)
         call dcopy(4, emiss_interp(:,emiss_order(1)), 1, matZ, 1)
         call dscal(4, alph, matS1, 1)
         call dscal(4, bet, matS2, 1)
         call dscal(4, gam, matZ, 1)
         call daxpy(4, dble(1.0), matS1, 1, matS2, 1)
         call daxpy(4, dble(1.0), matS2, 1, matZ, 1)

         call splint(xref, kaptot, dktot, ntau, logtau-delta_tau, k1)
         call splint(xref, kaptot, dktot, ntau, logtau, k2)
         call splint(xref, zdepth, deltaz, ntau, logtau-delta_tau, z1)
         call splint(xref, zdepth, deltaz, ntau, logtau, z2)

         do k=1, 4
             qmax = 0.5*(emiss_interp(k,emiss_order(1))*k1 +
     .               emiss_interp(k,emiss_order(2))*k2)*(z1-z2)
             if (qmax .gt. 0.0) then
                 matZ(k) = max(min(qmax, matZ(k)), dble(0.0))
             else
                 matZ(k) = min(max(qmax, matZ(k)), dble(0.0))
             endif

         enddo
c****      calculate the RHS of the equation.  Store in matZ
         call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)

c****      Solve the system of differential equations
         call dgesv(4,1,matX,4,IPIV,matZ,4,INFO)

         call dcopy(4, matZ, 1, Stokes, 1)

c****     Now do the same thing for the continuum
         dtau=(tau_interp_c(emiss_order(1))-
     .         tau_interp_c(emiss_order(2)))
         etau = 2.71828183**(-dtau)
         x = 1.0 - etau
         y = dtau - x
         z = dtau**2.0 - 2.0*y
         dtau_i = (tau_interp_c(emiss_order(2))-
     .             tau_interp_c(emiss_order(3)))
         alph = (z-dtau*y)/((dtau+dtau_i)*dtau_i)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))
c         alph = (z-dtau_i*y)/((dtau+dtau_i)*dtau)
c         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
c         gam = x+(z-(2.0*dtau_i + dtau)*y)/(dtau_i*(dtau+dtau_i))
c         continuum=etau*continuum+alph*emiss_interp(1,emiss_order(3))+
c     .              bet*emiss_interp(1,emiss_order(2))+
c     .              gam*emiss_interp(1,emiss_order(1))

         blah = alph*emiss_interp(1,emiss_order(3))+
     .          bet*emiss_interp(1,emiss_order(2))+
     .          gam*emiss_interp(1,emiss_order(1))
         call splint(xref, kaplam, dklam, ntau, logtau-delta_tau, k1)
         call splint(xref, kaplam, dklam, ntau, logtau, k2)

         qmax = 0.5*(emiss_interp(1,emiss_order(1))*k1 +
     .               emiss_interp(1,emiss_order(2))*k2)*(z1-z2)
         continuum = etau*continuum+max(min(blah, qmax), dble(0.0))

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

      enddo

      return
      end

      subroutine interp_opacities(logtau, k_interp,
     .       k_ord, e_interp, e_ord, tau_interp, tau_interp_c, mu, dtau)
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
      real*8 k_interp(4,4,2), e_interp(4,3)
      real*8 tau_interp(3), tau_interp_c(3), logtau
      real*8 t_tot, t_lam
      real*8 mu
      real*8 e_1, e_2, e_3, e_4
      real*8 phQ, phU, phV, psQ, psU, psV
      integer k_ord, e_ord

      call splint(xref, tlam, dtlam, ntau, logtau+dtau, t_lam)
      call splint(xref, tautot, dttot, ntau, logtau+dtau, t_tot)

      call splint(xref, phiQ, dphiQ, ntau, logtau, phQ)
      call splint(xref, phiU, dphiU, ntau, logtau, phU)
      call splint(xref, phiV, dphiV, ntau, logtau, phV)
      call splint(xref, psiQ, dpsiQ, ntau, logtau, psQ)
      call splint(xref, psiU, dpsiU, ntau, logtau, psU)
      call splint(xref, psiV, dpsiV, ntau, logtau, psV)

      call splint(xref, emission(1,:), de1, ntau, logtau+dtau, e_1)
      call splint(xref, emission(2,:), de2, ntau, logtau+dtau, e_2)
      call splint(xref, emission(3,:), de3, ntau, logtau+dtau, e_3)
      call splint(xref, emission(4,:), de4, ntau, logtau+dtau, e_4)

      k_interp(1,1,k_ord)=0.0
      k_interp(1,2,k_ord)=phQ
      k_interp(1,3,k_ord)=phU
      k_interp(1,4,k_ord)=phV
      k_interp(2,1,k_ord)=phQ
      k_interp(2,2,k_ord)=0.0
      k_interp(2,3,k_ord)=psV
      k_interp(2,4,k_ord)=-psU
      k_interp(3,1,k_ord)=phU
      k_interp(3,2,k_ord)=-psV
      k_interp(3,3,k_ord)=0.0
      k_interp(3,4,k_ord)=psQ
      k_interp(4,1,k_ord)=phV
      k_interp(4,2,k_ord)=psU
      k_interp(4,3,k_ord)=-psQ
      k_interp(4,4,k_ord)=0.0

      e_interp(1,e_ord) = e_1/mu
      e_interp(2,e_ord) = e_2/mu
      e_interp(3,e_ord) = e_3/mu
      e_interp(4,e_ord) = e_4/mu

      tau_interp(e_ord) = t_tot
      tau_interp_c(e_ord) =t_lam

      return
      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
