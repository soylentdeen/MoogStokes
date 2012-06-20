      SUBROUTINE LINTERPOLATE(Y_OLD, X_NEW, Y_NEW)
      IMPLICIT NONE
      include "Atmos.com"
      include "Linex.com"
      real*8 Y_OLD(100), X_NEW, Y_NEW, SLOPE
      integer I
      
c      DO I=1,NTAU-1
c          IF ((TAUREF(I+1)*KAPLAM(I+1)/(KAPREF(I+1)*MU))+1.0e-10
c     .         .GE.X_NEW) THEN
c              GOTO 10
c          ENDIF
c      ENDDO
c10    SLOPE=(Y_OLD(I+1)-Y_OLD(I))/(TAUREF(I+1)*KAPLAM(I+1)/
c     .   (KAPREF(I+1)*MU)-TAUREF(I)*KAPLAM(I)/
c     .   (KAPREF(I)*MU))
c      Y_NEW = Y_OLD(I)+SLOPE*(X_NEW-TAUREF(I)*KAPLAM(I)/
c     .   (KAPREF(I)*MU))
      DO I=1,NTAU-1
          IF ((TAUREF(I+1)*KAPLAM(I+1)/(KAPREF(I+1)))+1.0e-10
     .         .GE.X_NEW) THEN
              GOTO 10
          ENDIF
      ENDDO
10    SLOPE=(Y_OLD(I+1)-Y_OLD(I))/(TAUREF(I+1)*KAPLAM(I+1)/
     .   (KAPREF(I+1))-TAUREF(I)*KAPLAM(I)/
     .   (KAPREF(I)))
      Y_NEW = Y_OLD(I)+SLOPE*(X_NEW-TAUREF(I)*KAPLAM(I)/
     .   (KAPREF(I)))
      RETURN
      END

      SUBROUTINE MAT_DERIVS(N,X,Y,DYDX,RPAR,IPAR)
      implicit real*8 (a-h,o-z)
      DIMENSION Y(N), DYDX(N)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'
      include 'Angles.com'
      real*8 TEFF, k_interp(4,4), j_interp(4), mu, eta_0, kc_lam, kc_ref
      external Planck

      mu = cos(viewing_angle)
      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.X) THEN
             goto 10
         endif
      enddo
10    denom = (tauref(i+1)-tauref(i))
      do j=1,4
         do k=1,4
            slope=(kappa(j,k,i+1)-kappa(j,k,i))/denom
            k_interp(j,k)=kappa(j,k,i)+slope*
     .         (X-tauref(i))
         enddo
         slope=(emission(j,i+1)-emission(j,i))/denom
         j_interp(j)=emission(j,i)+slope*
     .         (X-tauref(i))
         slope=(t(i+1)-t(i))/denom
         TEFF=t(i)+slope*(X-tauref(i))
         slope=(kaplam(i+1)-kaplam(i))/denom
         kc_lam = (kaplam(i)+slope*(X-tauref(i)))
         slope=(kapref(i+1)-kapref(i))/denom
         kc_ref = (kapref(i)+slope*(X-tauref(i)))
         eta_0=kc_lam/kc_ref
      enddo

      DYDX(1)=((k_interp(1,1)*Y(1)+k_interp(1,2)*Y(2)+
     .         k_interp(1,3)*Y(3)+k_interp(1,4)*Y(4)) -
     .         j_interp(1))/mu
      DYDX(2)=((k_interp(2,1)*Y(1)+k_interp(2,2)*Y(2)+
     .         k_interp(2,3)*Y(3)+k_interp(2,4)*Y(4)) -
     .         j_interp(2))/mu
      DYDX(3)=((k_interp(3,1)*Y(1)+k_interp(3,2)*Y(2)+
     .         k_interp(3,3)*Y(3)+k_interp(3,4)*Y(4)) -
     .         j_interp(3))/mu
      DYDX(4)=((k_interp(4,1)*Y(1)+k_interp(4,2)*Y(2)+
     .         k_interp(4,3)*Y(3)+k_interp(4,4)*Y(4)) -
     .         j_interp(4))/mu
      DYDX(5)=eta_0*(Y(5)-Planck(TEFF))/mu

      RETURN
      END


      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN,XOUT)
      DIMENSION CON(8*ND),ICOMP(ND)
      real*8 X, XOLD
      include 'Atmos.com'
      include 'Linex.com'
      write (*,*) X, XOLD
      END
