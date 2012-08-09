      SUBROUTINE MAT_DERIVS(N,X,Y,DYDX,RPAR,IPAR)
      implicit real*8 (a-h,o-z)
      DIMENSION Y(N), DYDX(N)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'
      include 'Angles.com'
      real*8 TEFF, k_interp(4,4), j_interp(4), mu, eta_0, kc_lam
      real*8 kc_ref, run
      external Planck

      mu = cos(viewing_angle)
      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.X) THEN
             goto 10
         endif
      enddo
10    denom = (log10(tauref(i+1))-log10(tauref(i)))
      run = (log10(X) - log10(tauref(i)))
      do j=1,4
         do k=1,4
            slope=(kappa(j,k,i+1)-kappa(j,k,i))/denom
            k_interp(j,k)=kappa(j,k,i)+slope*run
         enddo
         slope=(emission(j,i+1)-emission(j,i))/denom
         j_interp(j)=emission(j,i)+slope*run
         slope=(t(i+1)-t(i))/denom
         TEFF=t(i)+slope*run
         slope=(kaplam(i+1)-kaplam(i))/denom
         kc_lam = (kaplam(i)+slope*run)
         slope=(kapref(i+1)-kapref(i))/denom
         kc_ref = (kapref(i)+slope*run)
         eta_0=kc_lam/kc_ref
      enddo

      DYDX(1)=((k_interp(1,1)*Y(1)+k_interp(1,2)*Y(2)+
     .         k_interp(1,3)*Y(3)+k_interp(1,4)*Y(4)) -
     .         j_interp(1)*k_interp(1,1))/mu
      DYDX(2)=((k_interp(2,1)*Y(1)+k_interp(2,2)*Y(2)+
     .         k_interp(2,3)*Y(3)+k_interp(2,4)*Y(4)) -
     .         j_interp(2)*k_interp(1,2))/mu
      DYDX(3)=((k_interp(3,1)*Y(1)+k_interp(3,2)*Y(2)+
     .         k_interp(3,3)*Y(3)+k_interp(3,4)*Y(4)) -
     .         j_interp(3)*k_interp(1,3))/mu
      DYDX(4)=((k_interp(4,1)*Y(1)+k_interp(4,2)*Y(2)+
     .         k_interp(4,3)*Y(3)+k_interp(4,4)*Y(4)) -
     .         j_interp(4)*k_interp(1,4))/mu
c      DYDX(5)=eta_0*(Y(5)-Planck(TEFF))/mu
      DYDX(5)=eta_0*(Y(5)-j_interp(1))/mu

      RETURN
      END


      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN,XOUT)
      implicit real*8 (a-h,o-z)
      DIMENSION CON(8*ND),ICOMP(ND)
      real*8 X, XOLD, Y(5), tau_c, tau_tot, kc_ref, k_tot, slope
      real*8 denom, run
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'

      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.X) THEN
             goto 10
         endif
      enddo
10    denom = (log10(tauref(i+1))-log10(tauref(i)))
      run = (log10(X) - log10(tauref(i)))
      slope=(kaptot(i+1)-kaptot(i))/denom
      k_tot = (kaptot(i)+slope*run)
      slope=(tauref(i+1)*kaptot(i+1)/kapref(i+1)-
     .       tauref(i)*kaptot(i)/kapref(i))/denom
      tau_tot = (tauref(i)*kaptot(i)/kapref(i)+slope*run)
      slope=(tauref(i+1)*kaplam(i+1)/kapref(i+1)-
     .       tauref(i)*kaplam(i)/kapref(i))/denom
      tau_c = (tauref(i)*kaplam(i)/kapref(i)+slope*run)
      
c      write (*,*) log10(X), k_tot, tau_tot, tau_c
      write (*,*) log10(X), Y(1), Y(5)
      END
