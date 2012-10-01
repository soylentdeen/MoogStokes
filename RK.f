
      subroutine traceStokes (phi_ang, chi_ang, mu)
c**********************************************************************
c     This routine performs the Runge-Kutta integration routine through
c     the photosphere
c**********************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      include 'Stokes.com'
      include 'Angles.com'
      include 'RungeKutta.com'
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V, D_bound
      real*8 tau_start, tau_stop, rtol, atol, rpar, Stokes_c(5), tmp
      real*8 bgfl, phi_ang, chi_ang, mu
      integer itol, iout, lwork, liwork, ipar
      parameter (NDGL=5, NRD=5)
      parameter (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
      DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
      external mat_derivs, Solout

c*****Sets up constants
      bgfl = dble(2e30)


      vewing_angle = acos(mu)

c***** For each layer in the atmosphere, calculate each element of the
c      opacity matrix and emission vector for the DELO algorithm
      do i=1,ntau
         phi_I=(phi_opacity(i,2)*sin(phi_ang)**2.0+
     .                (phi_opacity(i,1)+phi_opacity(i,3))*(1.0+
     .                cos(phi_ang)**2.0)/2.0)/2.0
         phi_Q=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*cos(2.0*chi_ang)/2.0
         phi_U=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*sin(2.0*chi_ang)/2.0
         phi_V=(phi_opacity(i,1)-phi_opacity(i,3))*cos(phi_ang)/2.0
         psi_Q=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*cos(2.0*chi_ang)/2.0
         psi_U=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
     .                /2.0)*sin(phi_ang)**2.0*sin(2.0*chi_ang)/2.0
         psi_V=(psi_opacity(i,1)-psi_opacity(i,3))*cos(phi_ang)/2.0

c*****  The total opacity (line+continuum)
         kaptot(i) = kaplam(i) + phi_I
c         write (*,*) tauref(i), phi_I, kaplam(i), phi_opacity(i,1),
c     .               phi_opacity(i,2), phi_opacity(i,3), kapref(i)
c*****  Assemble the Opacity matrix (K)
         kappa(1,1,i)=kaptot(i)/kapref(i)
         kappa(1,2,i)=phi_Q/kapref(i)
         kappa(1,3,i)=phi_U/kapref(i)
         kappa(1,4,i)=phi_V/kapref(i)
         kappa(2,1,i)=phi_Q/kapref(i)
         kappa(2,2,i)=kaptot(i)/kapref(i)
         kappa(2,3,i)=psi_V/kapref(i)
         kappa(2,4,i)=(-1.0*psi_U)/kapref(i)
         kappa(3,1,i)=phi_U/kapref(i)
         kappa(3,2,i)=(-1.0*psi_V)/kapref(i)
         kappa(3,3,i)=kaptot(i)/kapref(i)
         kappa(3,4,i)=psi_Q/kapref(i)
         kappa(4,1,i)=phi_V/kapref(i)
         kappa(4,2,i)=psi_U/kapref(i)
         kappa(4,3,i)=(-1.0*psi_Q)/kapref(i)
         kappa(4,4,i)=kaptot(i)/kapref(i)

c*****  Assumes LTE for the Source Function
         source = Planck(t(i))

c*****  Assembles the Emission matrix (J')
         emission(1,i)=source!*kaptot(i)/kapref(i)
         emission(2,i)=source!*phi_Q/kapref(i)
         emission(3,i)=source!*phi_U/kapref(i)
         emission(4,i)=source!*phi_V/kapref(i)

         eta0(i) = kaplam(i)/kapref(i)
      enddo

      call spline(xref, kaplam, ntau, bgfl, bgfl, dklam)
      call spline(xref, kaptot, ntau, bgfl, bgfl, dktot)
      call spline(xref, kapref, ntau, bgfl, bgfl, dkref)
      call spline(xref, kappa(1,1,:), ntau, bgfl, bgfl, dk11)
      call spline(xref, kappa(1,2,:), ntau, bgfl, bgfl, dk12)
      call spline(xref, kappa(1,3,:), ntau, bgfl, bgfl, dk13)
      call spline(xref, kappa(1,4,:), ntau, bgfl, bgfl, dk14)
      call spline(xref, kappa(2,1,:), ntau, bgfl, bgfl, dk21)
      call spline(xref, kappa(2,2,:), ntau, bgfl, bgfl, dk22)
      call spline(xref, kappa(2,3,:), ntau, bgfl, bgfl, dk23)
      call spline(xref, kappa(2,4,:), ntau, bgfl, bgfl, dk24)
      call spline(xref, kappa(3,1,:), ntau, bgfl, bgfl, dk31)
      call spline(xref, kappa(3,2,:), ntau, bgfl, bgfl, dk32)
      call spline(xref, kappa(3,3,:), ntau, bgfl, bgfl, dk33)
      call spline(xref, kappa(3,4,:), ntau, bgfl, bgfl, dk34)
      call spline(xref, kappa(4,1,:), ntau, bgfl, bgfl, dk41)
      call spline(xref, kappa(4,2,:), ntau, bgfl, bgfl, dk42)
      call spline(xref, kappa(4,3,:), ntau, bgfl, bgfl, dk43)
      call spline(xref, kappa(4,4,:), ntau, bgfl, bgfl, dk44)
      call spline(xref, emission(1,:), ntau, bgfl, bgfl, de1)
      call spline(xref, emission(2,:), ntau, bgfl, bgfl, de2)
      call spline(xref, emission(3,:), ntau, bgfl, bgfl, de3)
      call spline(xref, emission(4,:), ntau, bgfl, bgfl, de4)
      call spline(xref, eta0, ntau, bgfl, bgfl, deta0)

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm
c      D_bound=kappa(1,1,ntau)**2*(kappa(1,1,ntau)**2-kappa(1,2,ntau)**2-
c     .  kappa(1,3,ntau)**2-kappa(1,4,ntau)**2+kappa(3,4,ntau)**2+
c     .  kappa(4,2,ntau)**2+kappa(2,3,ntau)**2)-(kappa(1,2,ntau)*
c     .  kappa(3,4,ntau)+kappa(1,3,ntau)*kappa(4,2,ntau)+kappa(1,4,ntau)*
c     .  kappa(2,3,ntau))**2
c
c      dB = (emission(1,ntau)-emission(1,ntau-2))/
c     .      (tauref(ntau)-tauref(ntau-2))
c
c      tmp = cos(viewing_angle)*dB/D_bound*
c     .     kappa(1,1,ntau)*(kappa(1,1,ntau)**2+kappa(3,4,ntau)**2+
c     .     kappa(4,2,ntau)**2+kappa(2,3,ntau)**2)
c      Stokes_c(1) = emission(1,ntau)+cos(viewing_angle)*dB/D*
c     .     kappa(1,1,ntau)*(kappa(1,1,ntau)**2+kappa(3,4,ntau)**2+
c     .     kappa(4,2,ntau)**2+kappa(2,3,ntau)**2)
      Stokes_c(1) = emission(1,ntau)
      Stokes_c(2) = dble(0.0)
      Stokes_c(3) = dble(0.0)
      Stokes_c(4) = dble(0.0)
      Stokes_c(5) = emission(1,ntau)
c      Stokes_c(5) = Planck(t(ntau))*kaplam(ntau)/kapref(ntau)

      iout=0
      tau_start = tauref(ntau)
      tau_stop = tauref(1)
      itol = 0
      rtol = 1.0e-14
      atol = 1.0e-5
      do i=1,10
         iwork(i)=0
         work(i)=0.D0
      enddo

      call dop853(ndgl, mat_derivs, tau_start, Stokes_c, tau_stop,
     .            rtol, atol, itol, Solout, iout, work, lwork, iwork,
     .            liwork, rpar, ipar, idid)

      Stokes(1) = Stokes_c(1)!/Stokes_c(5)
      Stokes(2) = Stokes_c(2)!/Stokes_c(5)
      Stokes(3) = Stokes_c(3)!/Stokes_c(5)
      Stokes(4) = Stokes_c(4)!/Stokes_c(5)
      continuum = Stokes_c(5)
c      Stokes(1) = dble(99.9)
c      Stokes(1) = dble(0.0)
c      Stokes(1) = dble(0.0)
c      Stokes(1) = dble(0.0)
c      continuum = dble(100.0)

      return

      end

      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
