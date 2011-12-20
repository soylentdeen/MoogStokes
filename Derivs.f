      SUBROUTINE LINTERPOLATE(Y_OLD, X_NEW, Y_NEW)
      IMPLICIT NONE
      include "Atmos.com"
      include "Linex.com"
      real*8 Y_OLD(100), X_NEW, Y_NEW, SLOPE
      integer I
      
      DO I=1,NTAU-1
          IF ((TAUREF(I+1)*KAPLAM(I+1)/(KAPREF(I+1)*MU))+1.0e-10
     .         .GE.X_NEW) THEN
              GOTO 10
          ENDIF
      ENDDO
10    SLOPE=(Y_OLD(I+1)-Y_OLD(I))/(TAUREF(I+1)*KAPLAM(I+1)/
     .   (KAPREF(I+1)*MU)-TAUREF(I)*KAPLAM(I)/
     .   (KAPREF(I)*MU))
      Y_NEW = Y_OLD(I)+SLOPE*(X_NEW-TAUREF(I)*KAPLAM(I)/
     .   (KAPREF(I)*MU))
      RETURN
      END

      SUBROUTINE DERIVS(X,Y,DYDX)
C*****************************************************************************
C     This subroutine calculates the derivatives of the stokes parameters at
C     Tau = X.
C*****************************************************************************
      implicit NONE
      include "Atmos.com"
      include "Linex.com"
      real*8 X, Y(5), DYDX(5), B, TEFF, EI, EQ, EU, EV, ZQ, ZU, ZV
      
      CALL LINTERPOLATE(ETA_I, X, EI)
      CALL LINTERPOLATE(ETA_Q, X, EQ)
      CALL LINTERPOLATE(ETA_V, X, EV)
      CALL LINTERPOLATE(ZET_Q, X, ZQ)
      CALL LINTERPOLATE(ZET_V, X, ZV)
      CALL LINTERPOLATE(T, X, TEFF)

      CALL PLANCK(TEFF, B)
c      DYDX(1) = (1.0+EI)*Y(1)+EQ*Y(2)+EU*Y(3)+EV*Y(4) - (1.0+EI)*B
c      DYDX(2) = EQ*Y(1)+(1.0+EI)*Y(2)+ZV*Y(3)-ZU*Y(4) - (EQ)*B
c      DYDX(3) = EU*Y(1)-ZV*Y(2)+(1.0+EI)*Y(3)+ZQ*Y(4) - (EU)*B
c      DYDX(4) = EV*Y(1)+ZU*Y(2)-ZQ*Y(3)+(1.0+EI)*Y(4) - (EV)*B
      DYDX(1) = (1.0+EI)*Y(1)+EQ*Y(2)+EV*Y(4) - (1.0+EI)*B
      DYDX(2) = EQ*Y(1)+(1.0+EI)*Y(2)+ZV*Y(3) - (EQ)*B
      DYDX(3) = -ZV*Y(2)+(1.0+EI)*Y(3)+ZQ*Y(4)
      DYDX(4) = EV*Y(1)-ZQ*Y(3)+(1.0+EI)*Y(4) - (EV)*B
      DYDX(5) = Y(5)-B
c      write (*,*) "Derivatives: ", DYDX
      RETURN
      END
