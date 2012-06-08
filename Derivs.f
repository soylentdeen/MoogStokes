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
      real*8 TEFF, k_interp(4,4), j_interp(4), mu, eta_0, kc_lam, kc_ref
      external Planck

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
     .         j_interp(1))
      DYDX(2)=((k_interp(2,1)*Y(1)+k_interp(2,2)*Y(2)+
     .         k_interp(2,3)*Y(3)+k_interp(2,4)*Y(4)) -
     .         j_interp(2))
      DYDX(3)=((k_interp(3,1)*Y(1)+k_interp(3,2)*Y(2)+
     .         k_interp(3,3)*Y(3)+k_interp(3,4)*Y(4)) -
     .         j_interp(3))
      DYDX(4)=((k_interp(4,1)*Y(1)+k_interp(4,2)*Y(2)+
     .         k_interp(4,3)*Y(3)+k_interp(4,4)*Y(4)) -
     .         j_interp(4))
      DYDX(5)=eta_0*(Y(5)-Planck(TEFF))

      RETURN
      END


      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN,XOUT)
      DIMENSION CON(8*ND),ICOMP(ND)
      real*8 X, EI, EQ, EV, ZQ, ZV
      real*8 Y(N)
      include 'Atmos.com'
      include 'Linex.com'
      write (*,*) X, Y(1)
      END

      SUBROUTINE DERIVS(N,X,Y,DYDX,RPAR,IPAR)
C*****************************************************************************
C     This subroutine calculates the derivatives of the stokes parameters at
C     Tau = X.
C*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N), DYDX(N)
      include "Atmos.com"
      include "Linex.com"
      real*8 TEFF,EI,EQ,EU,EV,ZQ,ZU,ZV
      
      CALL LINTERPOLATE(ETA_I, X, EI)
      CALL LINTERPOLATE(ETA_Q, X, EQ)
      CALL LINTERPOLATE(ETA_V, X, EV)
      CALL LINTERPOLATE(ZET_Q, X, ZQ)
      CALL LINTERPOLATE(ZET_V, X, ZV)
      CALL LINTERPOLATE(T, X, TEFF)

c      CALL LINTERPOLATE(ETA_I, 10.0**X, EI)
c      CALL LINTERPOLATE(ETA_Q, 10.0**X, EQ)
c      CALL LINTERPOLATE(ETA_V, 10.0**X, EV)
c      CALL LINTERPOLATE(ZET_Q, 10.0**X, ZQ)
c      CALL LINTERPOLATE(ZET_V, 10.0**X, ZV)
c      CALL LINTERPOLATE(T, 10.0**X, TEFF)
      
c      CALL PLANCK(TEFF, B)
c      DYDX(1) = (1.0+EI)*Y(1)+EQ*Y(2)+EU*Y(3)+EV*Y(4) - (1.0+EI)*B
c      DYDX(2) = EQ*Y(1)+(1.0+EI)*Y(2)+ZV*Y(3)-ZU*Y(4) - (EQ)*B
c      DYDX(3) = EU*Y(1)-ZV*Y(2)+(1.0+EI)*Y(3)+ZQ*Y(4) - (EU)*B
c      DYDX(4) = EV*Y(1)+ZU*Y(2)-ZQ*Y(3)+(1.0+EI)*Y(4) - (EV)*B
c      DYDX(1) = ((1.0+EI)*Y(1)+EQ*Y(2)+EV*Y(4) - (1.0+EI)*B)*10**X*2.302
c      DYDX(2) = (EQ*Y(1)+(1.0+EI)*Y(2)+ZV*Y(3) - (EQ)*B)*10**X*2.302
c      DYDX(3) = (-ZV*Y(2)+(1.0+EI)*Y(3)+ZQ*Y(4))*10**X*2.302
c      DYDX(4) = (EV*Y(1)-ZQ*Y(3)+(1.0+EI)*Y(4) - (EV)*B)*10**X*2.302
c      DYDX(5) = (Y(5)-B)*10**X*2.302
      DYDX(1) = ((1.0+EI)*Y(1)+EQ*Y(2)+EV*Y(4) - (1.0+EI)*B)/mu
      DYDX(2) = (EQ*Y(1)+(1.0+EI)*Y(2)+ZV*Y(3) - EQ*B)/mu
      DYDX(3) = (-ZV*Y(2)+(1.0+EI)*Y(3)+ZQ*Y(4))/mu
      DYDX(4) = (EV*Y(1)-ZQ*Y(3)+(1.0+EI)*Y(4) - EV*B)/mu
      DYDX(5) = (Y(5)-B)/mu
      RETURN
      END

c      SUBROUTINE SOLOUT(NR,XOLD,X,Y,N,CON,ICOMP,ND,
c     .                  RPAR, IPAR, IRTRN, XOUT)
c      DIMENSION Y(N),CON(8*ND),ICOMP(ND)
c      DO I=1,N
c         Y(
c      ENDDO
