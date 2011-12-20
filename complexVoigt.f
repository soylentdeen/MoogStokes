      subroutine complexVoigt(x, y, u, v)
c******************************************************************************
c     This routine calculates the complex voigt function using the
c     algorithm of Z. Shippony and W. G. Read (JQSRT 1993)   
c     Computes the complex Voigt function: Integral from - to + infinity of:
c
c           (y/Pi)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt
c
c     and
c
c           (1/Pi)*(x-t)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt
c
c     Here:
c           x = sqrt(ln 2)* (v-v0)/aD             (x >= 0.0)
c           y = sqrt(ln 2)* (aL/aD)               (y >= 0.0)
c      where v s the wave number, v0 is the line center wave number, aL is
c      the Lorentzian line half-width, and aD is the Doppler line half-width
c******************************************************************************

      implicit real*8 (a-h,o-z) 
      PARAMETER (XF=5.8D0)
      PARAMETER (YS=1.0D0, YF=5.0D0)
      IF(Y.GT.YF.OR.X.GT.XF) THEN
          CALL VOIGTH(X,Y,U,V)                    ! REGION 4
c          write (*,*) "Region 4!"
      ELSE
          IF(Y.GT.YS) THEN
              IF(X.GE.YS) THEN
c                  write (*,*) "Region 3!"
                  CALL VOIGTR3(X,Y,U,V)           ! REGION 3
              ELSE
c                  write (*,*) "Region 2!"
                  CALL VOIGTR2(X,Y,U,V)           ! REGION 2
              ENDIF
          ELSE
c              write (*,*) "Region 1!"
              CALL VOIGTR1(X,Y,U,V)               ! REGION 1
          ENDIF
      ENDIF
      RETURN

      END

      SUBROUTINE VOIGTR1(X,Y,U,V)
C*****************************************************************************
C
C   Computes the Voigt function for the region (xs <= x < xf, y<=ys)
C     (Here xs = 0.0, xf = 5.8, ys = 1.0)
C
C   Based on the algorithm described in:
C
C      "Rapid Computation of the Voigt Profile" by S.R. Drayson,
C         JQSRT, Vol. 16, 611, (1976)
C
C    (Modified by Z. Shippony, JPL 1990)
C
C*****************************************************************************
      IMPLICIT NONE
      REAL*8 XTF, DELX, TINY
      INTEGER MAXD, MAXR, MAXF
      PARAMETER (XTF = 5.0D0, MAXD=11, MAXR=25, MAXF=13)
      PARAMETER (DELX = XTF/(MAXD-1),TINY=1.0D-12)
      REAL*8 DOWSON(MAXD),FCOEF(MAXF,MAXD),HN(MAXD),Q,P,Y2,EPSLN,
     .         DYR,DYI,ZI,ZR,RI(MAXR),V,U,DX,X,Y,TWOVSPI
      LOGICAL FIRST
      INTEGER*4 I,J,N
      SAVE FIRST, RI, HN, FCOEF, EPSLN
      DATA FIRST/.TRUE./
      DATA TWOVSPI /1.1283791670955125739D0/       ! 2.0/SQRT(PI)
C        DAWSON FUNCTION (F(Z)) TABLE, FOR Y=0, X = 0.0 TO: 5.0D DX = 0.5
      DATA DOWSON /  0.0D0 , 4.244363835020223D-1, 5.380795069127685D-1,
     & 4.282490710853986D-1, 3.013403889237919D-1, 2.230837221674355D-1,
     & 1.782710306105583D-1, 1.496215930807565D-1, 1.293480012360052D-1,
     & 1.140886102268250D-1, 1.021340744242768D-1 /
      IF (FIRST) THEN
          FIRST = .FALSE.
          EPSLN = DLOG(1.0D-8)
          DO I = 1, MAXR
              RI(I) = -2.0D0 / DFLOAT(I)
          ENDDO
C     COMPUTE DAWSON'S FUNCTION AT MESH POINTS
          DO I = 1, MAXD
              U = DFLOAT(I-1)*DELX
              HN(I) = U
              P = DOWSON(I)
              Q = 1.0D0 - 2.0D0 * U*P
              FCOEF(1,I) = Q
              DO J = 2, MAXF
                  V = (U*Q+P)*RI(J)
                  FCOEF(J,I) = V
                  P = Q
                  Q = V
              ENDDO
          ENDDO
      ENDIF
C      COMPUTE DAWSON'S FUNCTION AT X FROM TAYLOR SERIES
      Y2 = Y*Y
      J = INT(X/DELX + 0.5) + 1
      N = MIN0(J, MAXD)
      DX = X-HN(N)
      P = DABS(DX)
      IF(P.LE.TINY) THEN
          U = DOWSON(N)                ! F(X)
      ELSE
          U = MAXF
          V = 0.0D0
          IF (P.LT.1.0) U = EPSLN /DLOG(P)
          J = INT(U+0.5)
          I = MIN0(J, MAXF)
          DO J = I, 1, -1
              V = DX*V+FCOEF(J,N)
          ENDDO
          U = DX*V+DOWSON(N)           ! F(X)
      ENDIF
      IF(Y.LE.TINY)THEN
          V = TWOVSPI * U
          U = DEXP(-X*X)
          RETURN
      ENDIF
C  TAYLOR SERIES EXPANSION ABOUT Y = 0.0
C   COMPUTE # OF TERMS TO USE BEORE TRUNCATION:
      IF(X.LT.3.0) THEN
          J = 9
          Q = 16.0
      ELSE IF(X.LT.4.0) THEN
          J = 8
          Q = 17.0
      ELSE
          J = 7
          Q = 14.0
      ENDIF
      P = Q*(Y-0.1)/0.9
      I = J+INT(P+0.5)
      N = MIN0(I, MAXR)            ! NUM OF TERMS IN THE SERIES
      Q = U                        ! F(X)
      V = 1.0D0 - 2.0D0*X*U        ! F'(X)
      ZR = U
      ZI = Y*V
      DYR = 0.0D0
      DYI = Y                      ! DY = COMPLEX(0.0, Y)
      DO I = 2, N
          P = DYR
          DYR = -DYI*Y
          DYI = P*Y                ! = DY^N
          U = (X*V+Q)*RI(I)
          ZR = ZR + U*DYR
          ZI = ZI + U*DYI
          Q = V
          V = U
      ENDDO
      P = -2.0D0 * X*Y
      Q = DEXP(Y2-X*X)
      U = Q*DCOS(P) - TWOVSPI*ZI
      V = Q*DSIN(P) + TWOVSPI*ZR   ! EXP(-Z*Z) + 2*I*F()/SQRT(PI)
      RETURN
      END

      SUBROUTINE VOIGTR2(X, Y, U, V)
C*****************************************************************************
C   Computes the Voigt function in the region:
C             (0.0 <= x < 1.0, 1.0 <= y <= 5.0)
C
C   Hui's method was coded according to:
C        "Rapid computation of the Voigt and Complex Error Function"
C        by A.K. Hui, B.H. Aramstrong and A.A. Wray, JQSRT, Vol 19, 506 (1978)
C
C   Here:
C      x = sqrt(ln 2) * (v - v0)/aD    (x >= 0.0)
C      y = sqrt(ln 2) * aL/aD          (y >= 0.0)
C
C   Where v is the wave number, v0 = is the line center wave number, aL is the
C  Lorentzian line half-width and aD is the Doppler line half width.
C
C*****************************************************************************

      IMPLICIT real*8 (a-h, o-z)
      DIMENSION A(7), B(7)
      COMPLEX*16 Z, AZ, DZ, W
C  HUI'S (P=6) RATIONAL APPROXIMATION COEFFICIENTS:
      DATA A / 122.607931777104326D0, 214.382388694706425D0,
     &         181.928533092181549D0,  93.155580458138441D0,
     &          30.180142196210589D0,   5.912626209773153D0,
     &           0.564189583562615D0 /
      DATA B / 122.607931773875350D0, 352.730625110963558D0,
     &         457.334478783897737D0, 348.703917719495792D0,
     &         170.354001821091472D0,  53.992906912940207D0,
     &          10.479857114260399D0 /
C   REGION: (0.0 <= X < 1.0, 1.0 <= Y < 5.0)
C   HUI (P=6) ALGORITHM
      Z = DCMPLX(Y,-X)
      AZ = A(1)+Z*(A(2)+Z*(A(3)+Z*(A(4)+Z*(A(5)+Z*(A(6)+Z*A(7))))))
      DZ = B(1)+Z*(B(2)+Z*(B(3)+Z*(B(4)+Z*(B(5)+Z*(B(6)+Z*(B(7)+Z))))))
      W = AZ/DZ
      U = DREAL(W)
      V = DIMAG(W)
      RETURN
      END

      SUBROUTINE VOIGTR3(X, Y, U, V)
C*****************************************************************************
C  Computes the Voigt function for the region:
C       (xs <= x, x <= xf, ys <= y <= yf)
C
C (Here xs = 1.0, xf = 5.8, ys = 1.0, yf <= 5.0)
C
C This code was developed by Z. Shippony JPL 1990
C
C  Here:
C            x = sqrt(ln 2)* (v-v0)/aD           (x >= 0.0)
C            y = sqrt(ln 2)* aL/aD               (y >= 0.0)
C
C    Where v is the wave number, v0 is the line center wave number, aL is the
C  Lorentzian line half-width and aD is the Doppler line half-width
C*****************************************************************************

      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER (XS=1.0D0, XF=5.8D0)
      PARAMETER (YS=1.0D0, YF=5.0D0)
      PARAMETER (DX=1.0D0, DY=1.0D0)
      PARAMETER (MAXX=5, MAXY =5, NMAX=13)
      PARAMETER (TINY = 1.0D-12, EPS=1.0D-8)
      LOGICAL FIRST
      DIMENSION WR(MAXX, MAXY),WI(MAXX, MAXY),FN(NMAX)
      COMPLEX*16 WD(NMAX, MAXX, MAXY),Z,W, AZ,DZ,DS,SUM,SRM1
      SAVE FIRST, DXH, DYH, WD
C TABLE FOR W(Z) = W(X,Y)
C   X = 1.0,..., 5.0, STEP = 1.0
C   Y = 1.0,..., 5.0, STEP = 1.0
      DATA ((WR(I,J),J=1,5),I=1,5)/
     . 3.047442052569126D-1, 2.184926152748907D-1, 1.642611363929863D-1,
     . 1.298881599308405D-1, 1.067977383980653D-1, 1.402395813662779D-1,
     . 1.479527595120158D-1, 1.307574696698486D-1, 1.121394779021160D-1,
     . 9.649811260664139D-2, 6.531777728904696D-2, 9.271076642644334D-2,
     . 9.640250558304454D-2, 9.093390419476534D-2, 8.298773797690175D-2,
     . 3.628145648998864D-2, 5.968692961044590D-2, 6.979096164964831D-2,
     . 7.157043342636532D-2, 6.923620958049143D-2, 2.300313259405996D-2,
     . 4.064367633349437D-2, 5.122599656738663D-2, 5.599737714252388D-2,
     . 5.696543988817697D-2 /
      DATA ((WI(I,J),J=1,5),I=1,5)/
     . 2.082189382028316D-1, 9.299780939260188D-2, 5.019713513524858D-2,
     . 3.077886081705883D-2, 2.060408871468425D-2, 2.222134401798991D-1,
     . 1.311797170842179D-1, 8.111265047745665D-2, 5.348899385296694D-2,
     . 3.735165315636876D-2, 1.739183154163490D-1, 1.283169622282616D-1,
     . 9.123632600421876D-2, 6.559233052791429D-2, 4.838936520291309D-2,
     . 1.358389510006551D-1, 1.132100561244882D-1, 8.934000024036491D-2,
     . 6.937451861377146D-2, 5.407022703592907D-2, 1.103328325535800D-1,
     . 9.798731115657190D-2, 8.283691317190719D-2, 6.829488564492278D-2,
     . 5.583874277539103D-2 /
      DATA FIRST/.TRUE./
      DATA SRM1/(0.0D0,1.0D0)/
      DATA SPI/1.7724538509055160273D0/           ! SQRT(PI)
      IF(FIRST) THEN
C   COMPUTE ALL THE NECESSARY DERIVATIVES OF W(Z) AT MESH POINTS:
C   (FOR THE RECURSIVE EXPRESSION OF W(Z) DERIVATIVES, SEE "HANDBOOK OF
C   MATHEMATICAL FUNCTIONS", M. ABRAMOWITZ & A. STEGUN, DOVER PUBLICATIONS,
C   NOV. 1970, PP 298, EQ.: 7.1.12)
          FIRST = .FALSE.
          DXH = 0.5D0*DX
          DYH = 0.5D0*DY
          DO I = 1, NMAX
              FN(I) = -2.0D0/DFLOAT(I+1)
          ENDDO
          XX = XS - DX
          DO IX = 1, MAXX
              XX = XX+DX
              YY = YS-DY
              DO IY = 1, MAXY
                  YY = YY+DY
                  Z = DCMPLX(XX, YY)
                  W = DCMPLX(WR(IX, IY), WI(IX, IY))
                  DS = 2.0D0 * (SRM1 /SPI - Z * W)
                  WD(1,IX,IY) = DS
                  DO I = 1, NMAX-1
                      AZ = FN(I) * (W+Z*DS)
                      WD(I+1,IX,IY) = AZ
                      W = DS
                      DS = AZ
                  ENDDO
              ENDDO
          ENDDO
      ENDIF
C       REGION: (XS <= X <= XF, YS <= Y <= YF)
C         USING TAYLOR SERIES EXPANSION FOR W(Z)
      IX = 1 + INT((X-XS)/DX)
      XX =XS + (IX - 1)*DX
      IF(X-XX.GT.DXH) IX = IX + 1
      IF(IX.GT.MAXX) IX = MAXX
      XX = XS + (IX - 1) *DX
      IY = 1 + INT((Y-YS)/DY)
      YY = YS + (IY - 1) *DY
      IF(Y-YY.GT.DYH) IY = IY+1
      IF(IY.GT.MAXY) IY = MAXY
      YY = YS + (IY - 1) * DY
      U = WR(IX,IY)
      V = WI(IX, IY)
      P = X - XX
      Q = Y - YY
      IF(DABS(P).LT.TINY.AND.DABS(Q).LT.TINY) RETURN
      DZ = DCMPLX(P,Q)
      Q=DSQRT(P*P+Q*Q)
      P=EPS/CDABS(WD(1,IX,IY))
      P=DLOG(P)/DLOG(Q)
      N=1+MIN0(NMAX-1,INT(P))
      SUM=DCMPLX(0.0D0, 0.0D0)
      DO I=N,1,-1
          SUM = DZ*SUM+WD(I,IX,IY)
      ENDDO
      W=DZ*SUM
      U=U+DREAL(W)
      V=V+DIMAG(W)
      RETURN
      END

      SUBROUTINE VOIGTH(X,Y,U,V)
C*****************************************************************************
c  Computes the complex Voigt function: (i/pi)*integral from - to + infinity
c  of: exp(-t*t)/(z-t) dt, where z = (x+iy), x, y >= 0. Using Gauss-Hermite
c  Quadrature (8 points).  The real part is u, the integral of :
c
c      (y/pi)*Exp(-t*t)/(y*y + (x-t)*(x-t)) dt
c
c  The imaginary part is v, the integral of:
c
c      (1/pi)*(x-t)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt
c
c*****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=4)
      DIMENSION GX(N), GW(N)
C  NOTE: THE ROOTS ARE BOTH POSITIVE AND NEGATIVE!
      DATA GX/
     . 3.8118699020732211685D-1, 1.1571937124467801947D0,
     . 1.9816567566958429259D0 , 2.9306374202572440192D0 /
      DATA GW/
     . 6.6114701255824129103D-1, 2.0780232581489187954D-1,
     . 1.7077983007413475456D-2, 1.9960407221136761921D-4 /
      DATA PI / 3.1415926535897932385D0 /         !PI
      Y2 = Y * Y
      SUMU= 0.0D0
      SUMV= 0.0D0
      DO I = 1,N
          T = GX(I)
          XPT = X + T
          XMT = X - T
          FM = 1.0D0/ (Y2 + XMT * XMT)
          FP = 1.0D0/ (Y2 + XPT * XPT)
          SUMU = SUMU + GW(I) * (FM + FP)
          SUMV = SUMV + GW(I) * (XMT * FM + XPT * FP)
      ENDDO
      U = Y * SUMU/PI
      V = SUMV/PI
      RETURN
      END
