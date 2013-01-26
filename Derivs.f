      SUBROUTINE MAT_DERIVS(N,X,Y,DYDX,RPAR,IPAR)
      implicit real*8 (a-h,o-z)
      DIMENSION Y(N), DYDX(N)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'
      include 'Angles.com'
      include 'RungeKutta.com'
      real*8 mu, klam, kref, e_1, e_2, e_3, e_4, logtau, e_0
      real*8 phQ, phU, phV, psQ, psU, psV, kt
      external Planck

      mu = cos(viewing_angle)
      logtau = log10(X)
      call splint(xref, kaptot, dktot, ntau, logtau, kt)
      call splint(xref, phiQ, dphiQ, ntau, logtau, phQ)
      call splint(xref, phiU, dphiU, ntau, logtau, phU)
      call splint(xref, phiV, dphiV, ntau, logtau, phV)
      call splint(xref, psiQ, dpsiQ, ntau, logtau, psQ)
      call splint(xref, psiU, dpsiU, ntau, logtau, psU)
      call splint(xref, psiV, dpsiV, ntau, logtau, psV)
      call splint(xref, emission(1,:), de1, ntau, logtau, e_1)
      call splint(xref, emission(2,:), de2, ntau, logtau, e_2)
      call splint(xref, emission(3,:), de3, ntau, logtau, e_3)
      call splint(xref, emission(4,:), de4, ntau, logtau, e_4)
      call splint(xref, kaplam, dklam, ntau, logtau, klam)
      call splint(xref, kapref, dkref, ntau, logtau, kref)

      DYDX(1)=((kt/kref*Y(1)+phQ*Y(2)+
     .         phU*Y(3)+phV*Y(4)) -
     .         e_1*kt/kref)/mu
      DYDX(2)=((phQ*Y(1)+kt/kref*Y(2)+
     .         psV*Y(3)-psU*Y(4)) -
     .         e_2*phQ)/mu
      DYDX(3)=((phU*Y(1)-psV*Y(2)+
     .         kt/kref*Y(3)+psQ*Y(4)) -
     .         e_3*phU)/mu
      DYDX(4)=((phV*Y(1)+psU*Y(2)-
     .         psQ*Y(3)+kt/kref*Y(4)) -
     .         e_4*phV)/mu
      DYDX(5)=klam/kref*(Y(5)-e_1)/mu
c      DYDX(5)=e_0*(Y(5)-e_1)/mu

      RETURN
      END


      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN,XOUT)
      implicit real*8 (a-h,o-z)
      DIMENSION CON(8*ND),ICOMP(ND)
      real*8 X, XOLD, Y(5), tau_c, tau_tot, kc_ref, ktot, slope
      real*8 denom, run, klam, kref, logtau
      real*8 e_0, e_1, e_2, e_3, e_4
      real*8 phQ, phU, phV, psQ, psU, psV, kt
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'

      logtau = log10(X)
      call splint(xref, kaptot, dktot, ntau, logtau, kt)
      call splint(xref, phiQ, dphiQ, ntau, logtau, phQ)
      call splint(xref, phiU, dphiU, ntau, logtau, phU)
      call splint(xref, phiV, dphiV, ntau, logtau, phV)
      call splint(xref, psiQ, dpsiQ, ntau, logtau, psQ)
      call splint(xref, psiU, dpsiU, ntau, logtau, psU)
      call splint(xref, psiV, dpsiV, ntau, logtau, psV)
      call splint(xref, emission(1,:), de1, ntau, logtau, e_1)
      call splint(xref, emission(2,:), de2, ntau, logtau, e_2)
      call splint(xref, emission(3,:), de3, ntau, logtau, e_3)
      call splint(xref, emission(4,:), de4, ntau, logtau, e_4)
      call splint(xref, kaplam, dklam, ntau, logtau, klam)
      call splint(xref, kapref, dkref, ntau, logtau, kref)
c      write (*,*) logtau, Y(1), Y(2), Y(3), Y(4)
c      write (*,*) logtau, kt, kref, phQ, phU, phV, psQ, psU, psV, e_1
      write (*,*) logtau, Y(1), Y(2), Y(3), Y(4), Y(5)
c      write (*,*) logtau, e_1
      END
