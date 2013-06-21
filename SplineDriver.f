      subroutine spl_def (npts, x, y, knots, nknots, coeffs)
      implicit real*8 (a-h,o-z)
      integer npts, k, nest, nknots, lwork
      real*8 x(npts), y(npts), w(npts)
      parameter (k=3)
      real*8 knots(npts+k+1),coeffs(npts+k+1)
      real*8 xe, xb, s, fp, work(npts*(k+1)+(npts+k+1)*(7+3*k))
      integer iwork(npts+k+1), ier
      nest=npts+k+1
      lwork=npts*(k+1)+nest*(7+3*k)

      iopt = 0
      s = 0.0
      xb=x(1)
      xe=x(npts)
      do i = 1, npts
          w(i) = 1.0
      enddo

      call curfit(iopt, npts, x, y, w, xb, xe, k, s, nest,
     .   nknots, knots, coeffs, fp, work, lwork, iwork, ier)

      end

      real*8 function spl_ev(knots, nknots, coeffs, x_new)
      implicit real*8 (a-h,o-z)
      integer nknots, k, ier
      real*8 knots(nknots), coeffs(nknots), x_new, y_new

      k = 3

      call splev(knots, nknots, coeffs, k, x_new, y_new, 1, ier)

      spl_ev = y_new
      return
      end
