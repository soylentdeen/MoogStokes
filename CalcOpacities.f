
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
      real*8 voigt_val, faraday_voigt_val,voigt_x,voigt_y,gam_L, gam_D,
     .       sqrtpi
      complex Z, humlicek, w4

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
            voigt_x = (wave-wave1(j))/(wave1(j)*dopp(j,i)/2.997929d10)
            voigt_y = a(j,i)
            Z = cmplx(voigt_x, voigt_y)
            humlicek = w4(Z)
            voigt_val = REAL(humlicek,8)/sqrtpi
            faraday_voigt_val = AIMAG(humlicek)/sqrtpi
            phi_opacity(i,int(deltamj(j)+2))=phi_opacity(i,
     .               int(deltamj(j)+2))+ kapnu0(j,i)*voigt_val
            psi_opacity(i,int(deltamj(j)+2)) = psi_opacity(i,
     .               int(deltamj(j)+2))+ kapnu0(j,i)*faraday_voigt_val
         enddo                                     

                                                       
c*****do the same for the strong lines
         if (dostrong .gt. 0) then
            do j=nlines+1,nlines+nstrong
               voigt_x =(wave-wave1(j))/(wave1(j)*dopp(j,i)/2.997929d10)
               voigt_y = a(j,i)
c               write (*,*) i, wave-wave1(j), voigt_x
               Z = cmplx(voigt_x, voigt_y)
               humlicek = w4(Z)
               voigt_val = REAL(humlicek,8)/sqrtpi
               faraday_voigt_val = REAL(AIMAG(humlicek)/sqrtpi,8)
               phi_opacity(i,int(deltamj(j)+2))=phi_opacity(i,
     .                  int(deltamj(j)+2))+kapnu0(j,i)*voigt_val
               psi_opacity(i,int(deltamj(j)+2))=psi_opacity(i,
     .                  int(deltamj(j)+2))+kapnu0(j,i)*faraday_voigt_val
            enddo
c            j = nlines+1
         endif
      enddo      

      return                                              
      end                                                


