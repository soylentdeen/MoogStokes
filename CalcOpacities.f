
      subroutine calcopacities
c******************************************************************************
c     This routine calculates the line absorption coefficient and the line  
c     opacity at wavelength *wave* for all lines in the spectrum            
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      include 'Stokes.com'
      real*8 voigt, faraday_voigt, voigt_x, voigt_y, gam_L, gam_D,
     .       sqrtpi
      complex Z, hui, w4

c*****compute the total line opacity at each depth    
      
      sqrtpi = 1.772453851
      do i=1,ntau     
         phi_opacity(i,1) = 0.0
         phi_opacity(i,2) = 0.0
         phi_opacity(i,3) = 0.0
         psi_opacity(i,1) = 0.0
         psi_opacity(i,2) = 0.0
         psi_opacity(i,3) = 0.0
         do j=lim1,lim2
            gam_L = a(j,i)
            gam_D = dopp(j,i)/(2.997929d2*wave)
            voigt_x = dabs((1.0d8*(1./wave-1./wave1(j)))/(gam_D))
            voigt_y =(gam_L)
            Z = cmplx(voigt_x, voigt_y)
            hui = w4(Z)
            voigt = REAL(hui)/sqrtpi
            faraday_voigt = AIMAG(hui)/sqrtpi
            phi_opacity(i,int(deltamj(j)+2))=phi_opacity(i,
     .               int(deltamj(j)+2))+ kapnu0(j,i)*voigt
            psi_opacity(i,int(deltamj(j)+2)) = psi_opacity(i,
     .               int(deltamj(j)+2))+ kapnu0(j,i)*faraday_voigt
         enddo                                     

                                                       
c*****do the same for the strong lines
         if (dostrong .gt. 0) then
            do j=nlines+1,nlines+nstrong
               gam_L = a(j,i)
               gam_D = dopp(j,i)/(2.997929d2*wave)
               voigt_x = dabs((1.0d8*(1./wave-1./wave1(j)))/(gam_D))
               voigt_y =(gam_L)
               Z = cmplx(voigt_x, voigt_y)
               hui = w4(Z)
               voigt = REAL(hui)/sqrtpi
               faraday_voigt = AIMAG(hui)/sqrtpi
               phi_opacity(i,int(deltamj(j)+2))=phi_opacity(i,
     .                  int(deltamj(j)+2))+kapnu0(j,i)*voigt
               psi_opacity(i,int(deltamj(j)+2))=psi_opacity(i,
     .                  int(deltamj(j)+2))+kapnu0(j,i)*faraday_voigt
            enddo
         endif
      enddo      

      return                                              
      end                                                


