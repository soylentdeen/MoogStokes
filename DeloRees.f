
      subroutine traceStokes
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
      real*8 emiss_interp(4,2), kappa_interp(4,4,2)
      real*8 tau_interp(2), tau_interp_c(2), logtau
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V, n3
      real*8 h1, h2, dx, dtau, etau, n1, n2, n
      real*8 matS1(4), matS2(4), bgfl
      integer emiss_order(2), kappa_order(2)
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
c      do i=1,ntau
c         if (i.eq.1) then
c             zdepth(i) = 0.0
c         else
c             h1 = 1.0/(kapref(i-1))
c             h2 = 1.0/(kapref(i))
c             dtau = tauref(i)-tauref(i-1)
c             zdepth(i) = zdepth(i-1)+(h1+h2)/2.0*dtau
c         endif
c      enddo

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
         kaptot(i) = (kaplam(i) + phi_I)

C*****   Assemble the elements of the opacity matrix (K')
         k11(i)=0.0
         k12(i)=phi_Q/kaptot(i)
         k13(i)=phi_U/kaptot(i)
         k14(i)=phi_V/kaptot(i)
         k21(i)=phi_Q/kaptot(i)
         k22(i)=0.0
         k23(i)=psi_V/kaptot(i)
         k24(i)=(-1.0*psi_U)/kaptot(i)
         k31(i)=phi_U/kaptot(i)
         k32(i)=(-1.0*psi_V)/kaptot(i)
         k33(i)=0.0
         k34(i)=psi_Q/kaptot(i)
         k41(i)=phi_V/kaptot(i)
         k42(i)=psi_U/kaptot(i)
         k43(i)=(-1.0*psi_Q)/kaptot(i)
         k44(i)=0.0
c*****  Assemble the Opacity matrix (K')
c         kappa(1,1,i)=0.0
c         kappa(1,2,i)=phi_Q/kaptot(i)
c         kappa(1,3,i)=phi_U/kaptot(i)
c         kappa(1,4,i)=phi_V/kaptot(i)
c         kappa(2,1,i)=phi_Q/kaptot(i)
c         kappa(2,2,i)=0.0
c         kappa(2,3,i)=psi_V/kaptot(i)
c         kappa(2,4,i)=(-1.0*psi_U)/kaptot(i)
c         kappa(3,1,i)=phi_U/kaptot(i)
c         kappa(3,2,i)=(-1.0*psi_V)/kaptot(i)
c         kappa(3,3,i)=0.0
c         kappa(3,4,i)=psi_Q/kaptot(i)
c         kappa(4,1,i)=phi_V/kaptot(i)
c         kappa(4,2,i)=psi_U/kaptot(i)
c         kappa(4,3,i)=(-1.0*psi_Q)/kaptot(i)
c         kappa(4,4,i)=0.0

c*****  Assumes LTE for the Source Function
         source = Planck(t(i))

c*****  Assembles the Emission matrix (J')
         emission(1,i)=source
         emission(2,i)=source*phi_Q/kaptot(i)
         emission(3,i)=source*phi_U/kaptot(i)
         emission(4,i)=source*phi_V/kaptot(i)

         tautot(i) = kaptot(i)/kapref(i)*tauref(i)
         tlam(i) = kaplam(i)/kapref(i)*tauref(i)
      enddo

      call spline(xref, k11, ntau, bgfl, bgfl, dk11)
      call spline(xref, k12, ntau, bgfl, bgfl, dk12)
      call spline(xref, k13, ntau, bgfl, bgfl, dk13)
      call spline(xref, k14, ntau, bgfl, bgfl, dk14)
      call spline(xref, k21, ntau, bgfl, bgfl, dk21)
      call spline(xref, k22, ntau, bgfl, bgfl, dk22)
      call spline(xref, k23, ntau, bgfl, bgfl, dk23)
      call spline(xref, k24, ntau, bgfl, bgfl, dk24)
      call spline(xref, k31, ntau, bgfl, bgfl, dk31)
      call spline(xref, k32, ntau, bgfl, bgfl, dk32)
      call spline(xref, k33, ntau, bgfl, bgfl, dk33)
      call spline(xref, k34, ntau, bgfl, bgfl, dk34)
      call spline(xref, k41, ntau, bgfl, bgfl, dk41)
      call spline(xref, k42, ntau, bgfl, bgfl, dk42)
      call spline(xref, k43, ntau, bgfl, bgfl, dk43)
      call spline(xref, k44, ntau, bgfl, bgfl, dk44)
      call spline(xref, tautot, ntau, bgfl, bgfl, dttot)
      call spline(xref, tlam, ntau, bgfl, bgfl, dtlam)
      call spline(xref, emission(1,:), ntau, bgfl, bgfl, de1)
      call spline(xref, emission(2,:), ntau, bgfl, bgfl, de2)
      call spline(xref, emission(3,:), ntau, bgfl, bgfl, de3)
      call spline(xref, emission(4,:), ntau, bgfl, bgfl, de4)

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm

      Stokes(1) = emission(1,ntau)
      Stokes(2) = dble(0.0)
      Stokes(3) = dble(0.0)
      Stokes(4) = dble(0.0)
      continuum = Stokes(1)

      delta_tau = -0.01
      call dcopy(4, emission(:,ntau), 1, emiss_interp(:,1), 1)
      tau_interp(1) = tauref(ntau)*kaptot(ntau)/kapref(ntau)
      tau_interp_c(1) = tauref(ntau)*kaplam(ntau)/kapref(ntau)

      call interp_opacities(log10(tauref(ntau)),
     .        kappa_interp, 1, emiss_interp, 1, tau_interp,tau_interp_c)
      kappa_order(1) = 1
      kappa_order(2) = 2
      emiss_order(1) = 1
      emiss_order(2) = 2
      do logtau=log10(tauref(ntau))+delta_tau,
     .              log10(tauref(1)),delta_tau
         call interp_opacities(logtau, kappa_interp,
     .        kappa_order(2), emiss_interp, emiss_order(2), tau_interp,
     .        tau_interp_c)
c         write (*,*) logtau, tau_interp(emiss_order(2)), tau_interp_c(
c     .               emiss_order(2))
         dtau = (tau_interp(emiss_order(1))-tau_interp(emiss_order(2)))
     .              *cos(viewing_angle)
         etau = 2.71828183**(-dtau)
c         write (*,*) dtau, etau

         alph = 1.0-etau
         bet =(1.0-(1.0+dtau)*etau)/dtau

         call dcopy(16,ones, 1, matX, 1)
         call dcopy(16,ones, 1, matY, 1)
         call daxpy(16,dble(alph-bet),kappa_interp(:,:,kappa_order(2)),
     .              1,matX,1)
         call dscal(16,etau, matY,1)
         call daxpy(16,dble(-1.0*bet),kappa_interp(:,:,kappa_order(1)),
     .              1,matY,1)

         call dcopy(4, emiss_interp(:,emiss_order(2)), 1, matS1, 1)
         call dcopy(4, emiss_interp(:,emiss_order(1)), 1, matZ, 1)
         call dscal(4, alph-bet, matS1, 1)
         call dscal(4, bet, matZ, 1)
         call daxpy(4, dble(1.0), matS1, 1, matZ, 1)

c****      calculate the RHS of the equation.  Store in matZ
         call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)

c****      Solve the system of differential equations
         call dgesv(4,1,matX,4,IPIV,matZ,4,INFO)

         call dcopy(4, matZ, 1, Stokes, 1)

c****     Now do the same thing for the continuum
         dtau=(tau_interp_c(emiss_order(1))-
     .         tau_interp_c(emiss_order(2)))*cos(viewing_angle)
         etau = 2.71828183**(-dtau)
         alph = 1.0 - etau
         bet =(1.0-(1.0+dtau)*etau)/dtau
c         write (*,*) logtau, dtau, etau, alph, bet,
c     .         emiss_interp(1,emiss_order(2)), 
c     .         emiss_interp(1,emiss_order(1))
         continuum=etau*continuum+(alph-bet)*
     .              emiss_interp(1,emiss_order(2))
     .             +bet*emiss_interp(1,emiss_order(1))

         if (kappa_order(1).eq.1)then
             kappa_order(1) = 2
             kappa_order(2) = 1
         else
             kappa_order(1) = 1
             kappa_order(2) = 2
         endif
         if (emiss_order(1).eq.1) then
             emiss_order(1) = 2
             emiss_order(2) = 1
         elseif (emiss_order(1).eq.2) then
             emiss_order(1) = 1
             emiss_order(2) = 2
         endif
c         write (*,*) logtau, Stokes(1), continuum

      enddo

c      write (*,*) Stokes
      return
      end

      subroutine interp_opacities(logtau, k_interp,
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
      real*8 denom, run, kc_ref, t_tot, t_lam, tref, kref
      real*8 k_ref, k_tot, k_lam, e_1, e_2, e_3, e_4
      real*8 k_11, k_12, k_13, k_14, k_21, k_22, k_23, k_24
      real*8 k_31, k_32, k_33, k_34, k_41, k_42, k_43, k_44
      integer k_ord, e_ord

      call splint(xref, tlam, dtlam, ntau, logtau, t_lam)
      call splint(xref, tautot, dttot, ntau, logtau, t_tot)

      call splint(xref, k11, dk11, ntau, logtau, k_11)
      call splint(xref, k12, dk12, ntau, logtau, k_12)
      call splint(xref, k13, dk13, ntau, logtau, k_13)
      call splint(xref, k14, dk14, ntau, logtau, k_14)
      call splint(xref, k21, dk21, ntau, logtau, k_21)
      call splint(xref, k22, dk22, ntau, logtau, k_22)
      call splint(xref, k23, dk23, ntau, logtau, k_23)
      call splint(xref, k24, dk24, ntau, logtau, k_24)
      call splint(xref, k31, dk31, ntau, logtau, k_31)
      call splint(xref, k32, dk32, ntau, logtau, k_32)
      call splint(xref, k33, dk33, ntau, logtau, k_33)
      call splint(xref, k34, dk34, ntau, logtau, k_34)
      call splint(xref, k41, dk41, ntau, logtau, k_41)
      call splint(xref, k42, dk42, ntau, logtau, k_42)
      call splint(xref, k43, dk43, ntau, logtau, k_43)
      call splint(xref, k44, dk44, ntau, logtau, k_44)

      call splint(xref, emission(1,:), de1, ntau, logtau, e_1)
      call splint(xref, emission(2,:), de2, ntau, logtau, e_2)
      call splint(xref, emission(3,:), de3, ntau, logtau, e_3)
      call splint(xref, emission(4,:), de4, ntau, logtau, e_4)

      k_interp(1,1,k_ord)=k_11
      k_interp(1,2,k_ord)=k_12
      k_interp(1,3,k_ord)=k_13
      k_interp(1,4,k_ord)=k_14
      k_interp(2,1,k_ord)=k_21
      k_interp(2,2,k_ord)=k_22
      k_interp(2,3,k_ord)=k_23
      k_interp(2,4,k_ord)=k_24
      k_interp(3,1,k_ord)=k_31
      k_interp(3,2,k_ord)=k_32
      k_interp(3,3,k_ord)=k_33
      k_interp(3,4,k_ord)=k_34
      k_interp(4,1,k_ord)=k_41
      k_interp(4,2,k_ord)=k_42
      k_interp(4,3,k_ord)=k_43
      k_interp(4,4,k_ord)=k_44

      e_interp(1,e_ord) = e_1
      e_interp(2,e_ord) = e_2
      e_interp(3,e_ord) = e_3
      e_interp(4,e_ord) = e_4

      tau_interp(e_ord) = t_tot
      tau_interp_c(e_ord) = t_lam
c      tau_interp(e_ord) = 10.0**logtau * k_tot/k_ref
c      tau_interp_c(e_ord) = 10.0**logtau * k_lam/k_ref

c      write (*,*) logtau, tau_interp(e_ord), tau_interp_c(e_ord)
      return
      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
