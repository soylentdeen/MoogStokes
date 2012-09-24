      function W4(Z)
c******************************************************************************
c     This routine calculates the complex voigt function using the
c     algorithm of Humlicek 1982.
c******************************************************************************

      COMPLEX W4,Z, T, U
      REAL*8 X, Y, S
c      z = cmplx(x, y)
      X=REAL(Z)
      Y=AIMAG(Z)
      T=CMPLX(Y,-X)
      S=ABS(X)+Y
      IF (S.LT.15.) GOTO 1
C********  REGION 1
c      write (*,*) "Region 1"
c      write (*,*) Z, X, Y, T
      W4=T*.5641896/(.5+T*T)
      RETURN
1     IF(S.LT.5.5)GOTO 2
C********  REGION 2
c      write (*,*) "Region 2"
      U=T*T
      W4=T*(1.410474+U*.5641896)/(.75+U*(3.+U))
      RETURN
2     IF(Y.LT..195*ABS(X)-.176)GOTO 3
C********  REGION 3
c      write (*,*) "Region 3"
      W4=(16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*.5642236))))/
     /(16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))))
      RETURN
C********  REGION 4
3     U=T*T
c      write (*,*) "Region 4"
      W4=CEXP(U)-T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313-U*
     *(35.76683-U*(1.320522-U*.56419))))))/(32066.6-U*(24322.84-U*
     *(9022.229-U*(2186.181-U*(364.2191-U*(61.57037-U*(1.841439-U)))))))
      
      RETURN

      END

