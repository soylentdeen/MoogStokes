
      subroutine genwave
c*****************************************************************
c     This routine generates the wavelength grid which will be 
c     used in the synthesis.  It starts with a baseline coarsely
c     sampled grid, and sprinkles in patterns of finely samples of
c     wavelength points at the locations of the spectral lines
c*****************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
