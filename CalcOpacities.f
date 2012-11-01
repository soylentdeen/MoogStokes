
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
     .       sqrtpi, v
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
            gam_L = a(j,i)
            gam_D = dopp(j,i)/(2.997929d2*wave)
            voigt_x = dabs((1.0d8*(1./wave-1./wave1(j)))/(gam_D))
            voigt_y =(gam_L)
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
               gam_L = a(j,i)
               gam_D = dopp(j,i)/(2.997929d2*wave)
               voigt_x = dabs((1.0d8*(1./wave-1./wave1(j)))/(gam_D))
               voigt_y =(gam_L)
               Z = cmplx(voigt_x, voigt_y)
               humlicek = w4(Z)
               voigt_val = REAL(humlicek,8)/sqrtpi
c               v = 2.997929d10*dabs(wave-wave1(j))/(wave1(j)*dopp(j,i))
c               voigt_val = voigt(a(j,i),v
               faraday_voigt_val = REAL(AIMAG(humlicek)/sqrtpi,8)
c               if ((i.eq.30).AND.(j.eq.nlines+1)) then
c                   write (*,*) wave, faraday_voigt_val
c               endif
c               x = real(Z)
c               y = aimag(Z)
c               S = abs(X)+Y
c               if (S.GT.15.) then             ! Region 1
c               if ((S.LT.15.).AND.(S.GT.5.5)) then   ! Region 2
c               if ((S.LT.5.5).AND.(Y.GT..195*abs(X)-.176)) then !Region 3
c               if ((S.LT.5.5).AND.(Y.LT..195*abs(X)-.176)) then !Region 4
c                   write (*,*) wave, i, a(j,i), v,
c     .                    voigt_val, voigt(a(j,i),v)
c               endif
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


