
      subroutine rungekutta
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
      real*8 phi_I, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V
      real*8 tau_start, tau_stop, rtol, atol, rpar
      integer itol, iout, lwork, liwork, ipar
      parameter (NDGL=4, NRD=4)
      parameter (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
      DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
      external mat_derivs

c*****Sets up constants

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
         emission(1,i)=source*kaptot(i)/kapref(i)
         emission(2,i)=source*phi_Q/kapref(i)
         emission(3,i)=source*phi_U/kapref(i)
         emission(4,i)=source*phi_V/kapref(i)

      enddo

c*****   Trace the Stokes parameters through the atmosphere
c            via the quadratic DELO algorithm

      Stokes(1) = emission(1,ntau)
      Stokes(2) = 0.0
      Stokes(3) = 0.0
      Stokes(4) = 0.0
      continuum = Stokes(1)

      iout=0
      tau_start = tauref(ntau)
      tau_stop = tauref(1)
      itol = 0
      rtol = 1.0e-9
      atol = 1.0e-9
      do 10 i=1,10
         iwork(i)=0
10       work(i)=0.D0
c      iwork(5) = NDGL

      call dop853(ndgl, mat_derivs, tau_start, Stokes, tau_stop,
     .            rtol, atol, itol, junk, iout, work, lwork, iwork,
     .            liwork, rpar, ipar, idid)

      write (*,*) idid, iwork(17), iwork(18), iwork(19), iwork(20)
      write (*,*) Stokes(1), Stokes(2), Stokes(3), Stokes(4)
c      write (*,*) Stokes
      return
      end


      real*8 function Planck(temperature)
      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      Planck = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .     (dexp(1.43879d+08/(wave*temperature))-1.0d+00))
      return
      end
