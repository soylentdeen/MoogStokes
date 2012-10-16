      SUBROUTINE MAT_DERIVS(N,X,Y,DYDX,RPAR,IPAR)
      implicit real*8 (a-h,o-z)
      DIMENSION Y(N), DYDX(N)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'
      include 'Angles.com'
      include 'RungeKutta.com'
      real*8 mu, klam, kref, e_1, e_2, e_3, e_4, logtau, e_0
      real*8 k_11, k_12, k_13, k_14
      real*8 k_21, k_22, k_23, k_24
      real*8 k_31, k_32, k_33, k_34
      real*8 k_41, k_42, k_43, k_44
      external Planck

      mu = cos(viewing_angle)
      logtau = log10(X)
      call splint(xref, kappa(1,1,:), dk11, ntau, logtau, k_11)
      call splint(xref, kappa(1,2,:), dk12, ntau, logtau, k_12)
      call splint(xref, kappa(1,3,:), dk13, ntau, logtau, k_13)
      call splint(xref, kappa(1,4,:), dk14, ntau, logtau, k_14)
      call splint(xref, kappa(2,1,:), dk21, ntau, logtau, k_21)
      call splint(xref, kappa(2,2,:), dk22, ntau, logtau, k_22)
      call splint(xref, kappa(2,3,:), dk23, ntau, logtau, k_23)
      call splint(xref, kappa(2,4,:), dk24, ntau, logtau, k_24)
      call splint(xref, kappa(3,1,:), dk31, ntau, logtau, k_31)
      call splint(xref, kappa(3,2,:), dk32, ntau, logtau, k_32)
      call splint(xref, kappa(3,3,:), dk33, ntau, logtau, k_33)
      call splint(xref, kappa(3,4,:), dk34, ntau, logtau, k_34)
      call splint(xref, kappa(4,1,:), dk41, ntau, logtau, k_41)
      call splint(xref, kappa(4,2,:), dk42, ntau, logtau, k_42)
      call splint(xref, kappa(4,3,:), dk43, ntau, logtau, k_43)
      call splint(xref, kappa(4,4,:), dk44, ntau, logtau, k_44)
      call splint(xref, emission(1,:), de1, ntau, logtau, e_1)
      call splint(xref, emission(2,:), de2, ntau, logtau, e_2)
      call splint(xref, emission(3,:), de3, ntau, logtau, e_3)
      call splint(xref, emission(4,:), de4, ntau, logtau, e_4)
      call splint(xref, kaplam, dklam, ntau, logtau, klam)
      call splint(xref, kapref, dkref, ntau, logtau, kref)
c      call splint(xref, eta0, deta0, ntau, logtau, e_0)

      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.X) THEN
             goto 10
         endif
      enddo
10    denom =(log10(tauref(i+1))-log10(tauref(i)))
      run =(logtau - log10(tauref(i)))

      slope=(eta0(i+1)-eta0(i))/denom
      e_0 = (eta0(i)+slope*run)

      DYDX(1)=((k_11*Y(1)+k_12*Y(2)+
     .         k_13*Y(3)+k_14*Y(4)) -
     .         e_1*k_11)/mu
      DYDX(2)=((k_21*Y(1)+k_22*Y(2)+
     .         k_23*Y(3)+k_24*Y(4)) -
     .         e_2*k_12)/mu
      DYDX(3)=((k_31*Y(1)+k_32*Y(2)+
     .         k_33*Y(3)+k_34*Y(4)) -
     .         e_3*k_13)/mu
      DYDX(4)=((k_41*Y(1)+k_42*Y(2)+
     .         k_43*Y(3)+k_44*Y(4)) -
     .         e_4*k_14)/mu
c      DYDX(5)=klam/kref*(Y(5)-e_1)/mu
      DYDX(5)=e_0*(Y(5)-e_1)/mu

      RETURN
      END


      subroutine Solout(NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN,XOUT)
      implicit real*8 (a-h,o-z)
      DIMENSION CON(8*ND),ICOMP(ND)
      real*8 X, XOLD, Y(5), tau_c, tau_tot, kc_ref, ktot, slope
      real*8 denom, run, klam, kref, logtau, e_0, k_11, k_11a, e_1
      include 'Atmos.com'
      include 'Linex.com'
      include 'Stokes.com'

      logtau = log10(X)
      do i=1,ntau-1
         if (tauref(i+1)+1.0e-10.ge.X) THEN
             goto 10
         endif
      enddo
10    denom =(log10(tauref(i+1))-log10(tauref(i)))
      run =(logtau - log10(tauref(i)))

      slope=(eta0(i+1)-eta0(i))/denom
      e_0 = (eta0(i)+slope*run)

      call splint(xref, kaptot, dktot, ntau, logtau, ktot)
      call splint(xref, kaplam, dklam, ntau, logtau, klam)
      call splint(xref, kapref, dkref, ntau, logtau, kref)
      call splint(xref, kappa(1,1,:), dk11, ntau, logtau, k_11)
      call splint(xref, emission(1,:), de1, ntau, logtau, e_1)
c      write (*,*) log10(X), klam, kref, e_0, e_1, ktot
c      write (*,*) log10(X), k_11, k_11a
      write (*,*) X, ktot, klam
      END
