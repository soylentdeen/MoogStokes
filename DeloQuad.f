
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
      real*8 matS1(4), matS2(4), bgfl, ztau
      real*8 k1, k2, k3, z1, z2, deltaz(100), qmax
      real*8 phi_ang, chi_ang, mu, midpoint, delz, delta_tau
c      real*8 zdepth(20000), z_knots(20000), z_coeffs(20000),taus(20000)
      real*8 ktot, klam, kref
      integer emiss_order(3), kappa_order(2)!, n_z_knots
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
         emission(i)=source
      enddo

      call spl_def(ntau, xref, kaplam, klam_knots, n_klam_knots,
     .           klam_coeffs)
      call spl_def(ntau, xref, kaptot, ktot_knots, n_ktot_knots,
     .           ktot_coeffs)
      call spl_def(ntau, xref, emission, e_knots, n_e_knots,
     .           e_coeffs)

      do i=1, nz
         if (i.eq.1) then
            klam=spl_ev(klam_knots,n_klam_knots,klam_coeffs,taus(i))
            ktot=spl_ev(ktot_knots,n_ktot_knots,ktot_coeffs,taus(i))
            kref=spl_ev(kref_knots,n_kref_knots,kref_coeffs,taus(i))
            tlam(i) = klam/kref*10.0**(taus(i))
            ttot(i) = ktot/kref*10.0**(taus(i))
         else
            delz = (zdepth(i)-zdepth(i-1))/mu
            h1=spl_ev(ktot_knots,n_ktot_knots,ktot_coeffs,taus(i))
            h2=spl_ev(ktot_knots,n_ktot_knots,ktot_coeffs,taus(i-1))
            dtautot = delz*(h1+h2)/2.0
            ttot(i) = ttot(i-1)+dtautot
            h1=spl_ev(klam_knots,n_klam_knots,klam_coeffs,taus(i))
            h2=spl_ev(klam_knots,n_klam_knots,klam_coeffs,taus(i-1))
            dtaulam = delz*(h1+h2)/2.0
            tlam(i) = tlam(i-1)+dtaulam
         endif
c         write (*,*) i, taus(i), ttot(i)
      enddo
c      read (*,*)
c      h1=kaptot(ntau)
c      h2=spl_ev(ktot_knots,n_ktot_knots,ktot_coeffs,taus(i))
c      dtautot = delz*(h1+h2)/2.0
c      ttot(nz) = ttot(nz-1)+dtautot
c      h1=kaplam(ntau)
c      h2=spl_ev(klam_knots,n_klam_knots,klam_coeffs,taus(i))
c      dtaulam = delz*(h1+h2)/2.0
c      tlam(nz) = tlam(nz-1)+dtaulam

      call spl_def(nz, taus, tlam, tlam_knots, n_tlam_knots,
     .           tlam_coeffs)
      call spl_def(nz, taus, ttot, ttot_knots, n_ttot_knots,
     .           ttot_coeffs)
      call spl_def(ntau, xref, phiQ, phiQ_knots, n_phiQ_knots,
     .           phiQ_coeffs)
      call spl_def(ntau, xref, phiU, phiU_knots, n_phiU_knots,
     .           phiU_coeffs)
      call spl_def(ntau, xref, phiV, phiV_knots, n_phiV_knots,
     .           phiV_coeffs)
      call spl_def(ntau, xref, psiQ, psiQ_knots, n_psiQ_knots,
     .           psiQ_coeffs)
      call spl_def(ntau, xref, psiU, psiU_knots, n_psiU_knots,
     .           psiU_coeffs)
      call spl_def(ntau, xref, psiV, psiV_knots, n_psiV_knots,
     .           psiV_coeffs)

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm


      zd1=spl_ev(z_knots,n_z_knots,z_coeffs,xref(1))
      zd2=spl_ev(z_knots,n_z_knots,z_coeffs,xref(2))
      dbdz = (Planck(t(1)) - Planck(t(2)) )/(zd1-zd2)
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

c      Stokes(1) = Planck(t(1))
c      Stokes(2) = 0.0
c      Stokes(3) = 0.0
c      Stokes(4) = 0.0
c      continuum = Stokes(1)

      delta_tau = -0.05
      emiss_interp(1,1) = emission(ntau)
      emiss_interp(2,1) = emission(ntau)*phiQ(ntau)
      emiss_interp(3,1) = emission(ntau)*phiU(ntau)
      emiss_interp(4,1) = emission(ntau)*phiV(ntau)
      tau_interp(1) = ttot(nz)
      tau_interp_c(1) = tlam(nz)

      call interp_opacities(log10(tauref(ntau)),
     .   kappa_interp, 1, emiss_interp, 2, tau_interp,tau_interp_c,
     .   mu, delta_tau)
      kappa_order(1) = 1
      kappa_order(2) = 2
      emiss_order(1) = 1
      emiss_order(2) = 2
      emiss_order(3) = 3
      do logtau=log10(tauref(ntau))+delta_tau,
     .              log10(tauref(1)),delta_tau
         call interp_opacities(logtau, kappa_interp,
     .        kappa_order(2), emiss_interp, emiss_order(3), tau_interp,
     .        tau_interp_c, mu, delta_tau)
         dtau = (tau_interp(emiss_order(1))-tau_interp(emiss_order(2)))
         etau = dexp(-dtau)

c         write (*,*) tau_interp
c         read (*,*)

         alph = dble(1.0-etau)
         bet =dble((1.0-(1.0+dtau)*etau)/dtau)

         call dcopy(16,ones, 1, matX, 1)
         call dcopy(16,ones, 1, matY, 1)
         call daxpy(16,dble(alph-bet),kappa_interp(:,:,kappa_order(2)),
     .              1,matX,1)
         call dscal(16,etau, matY,1)
         call daxpy(16,dble(-1.0*bet),kappa_interp(:,:,kappa_order(1)),
     .              1,matY,1)

         x = dble(1.0-etau)
         y = dble(dtau - x)
         z = dble(dtau**2.0 - 2.0 * y)
         dtau_i =(tau_interp(emiss_order(2))-tau_interp(emiss_order(3)))
         alph = (z-dtau*y)/((dtau+dtau_i)*dtau_i)
         bet = ((dtau_i+dtau)*y-z)/(dtau*dtau_i)
         gam = x+(z-(dtau_i+2*dtau)*y)/(dtau*(dtau+dtau_i))

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
     .         tau_interp_c(emiss_order(2)))
         etau = dexp(-dtau)
         x = dble(dble(1.0) - etau)
         y = dble(dtau - x)
         z = dble(dtau**2.0 - 2.0*y)
         dtau_i = (tau_interp_c(emiss_order(2))-
     .             tau_interp_c(emiss_order(3)))
         alph = (z-dtau*y)/((dtau+dtau_i)*dtau_i)
         bet = ((dtau_i+dtau)*y - z)/(dtau*dtau_i)
         gam = x+(z-(dtau_i + 2*dtau)*y)/(dtau*(dtau+dtau_i))
         continuum=etau*continuum+alph*emiss_interp(1,emiss_order(3))+
     .              bet*emiss_interp(1,emiss_order(2))+
     .              gam*emiss_interp(1,emiss_order(1))

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
c         write (*,*) logtau, dtau, Stokes(1), continuum
      enddo
c      read (*,*)
      return
      end

      subroutine interp_opacities(logtau, k_interp,
     .   k_ord, e_interp, e_ord, tau_interp, tau_interp_c, mu, deltau)
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
      real*8 t_tot, t_lam, deltau, mu
      real*8 e_1, e_2, e_3, e_4
      real*8 phQ, phU, phV, psQ, psU, psV
      integer k_ord, e_ord

      t_lam=spl_ev(tlam_knots,n_tlam_knots,tlam_coeffs,logtau+deltau)
      t_tot=spl_ev(ttot_knots,n_ttot_knots,ttot_coeffs,logtau+deltau)

      phQ=spl_ev(phiQ_knots,n_phiQ_knots,phiQ_coeffs,logtau)
      phU=spl_ev(phiU_knots,n_phiU_knots,phiU_coeffs,logtau)
      phV=spl_ev(phiV_knots,n_phiV_knots,phiV_coeffs,logtau)
      psQ=spl_ev(psiQ_knots,n_psiQ_knots,psiQ_coeffs,logtau)
      psU=spl_ev(psiU_knots,n_psiU_knots,psiU_coeffs,logtau)
      psV=spl_ev(psiV_knots,n_psiV_knots,psiV_coeffs,logtau)

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

      phQ=spl_ev(phiQ_knots,n_phiQ_knots,phiQ_coeffs,logtau+deltau)
      phU=spl_ev(phiU_knots,n_phiU_knots,phiU_coeffs,logtau+deltau)
      phV=spl_ev(phiV_knots,n_phiV_knots,phiV_coeffs,logtau+deltau)

      emiss = spl_ev(e_knots,n_e_knots,e_coeffs,logtau+deltau)

      e_interp(1,e_ord) = emiss/mu
      e_interp(2,e_ord) = emiss*phQ/mu
      e_interp(3,e_ord) = emiss*phU/mu
      e_interp(4,e_ord) = emiss*phV/mu

      tau_interp(e_ord) = t_tot
      tau_interp_c(e_ord) = t_lam

      return
      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
