
      subroutine taukap
c******************************************************************************
c     This routine calculates the line absorption coefficient and the line  
c     opacity at wavelength *wave* for all lines in the spectrum            
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      real*8 kapnu_I(100), kapnu_Q(100), kapnu_V(100), kapnu_U(100),
     .       new_voigt, new_fv, voigt_x, voigt_y, gam_L, gam_D,
     .       zetnu_Q(100), zetnu_V(100), zetnu_U(100), sqrtpi, old_voigt
      complex Z, hui, w4

c*****compute the total line opacity at each depth    
      
      sqrtpi = 1.772453851
      do i=1,ntau     
         kapnu_I(i) = 0.0
         kapnu_Q(i) = 0.0
         kapnu_U(i) = 0.0
         kapnu_V(i) = 0.0
         zetnu_Q(i) = 0.0
         zetnu_U(i) = 0.0
         zetnu_V(i) = 0.0
         kapnu(i) = 0.0
         do j=lim1,lim2
c            v = 2.997929d10*dabs(wave-wave1(j))/
c     .             (wave1(j)*dopp(j,i))            
c            old_voigt = voigt(a(j,i),v)
            gam_L = a(j,i)
            gam_D = dopp(j,i)/(2.997929d2*wave)
            voigt_x = !sqrt(log(2.0))*
     .                 dabs((1.0d8*(1./wave-1./wave1(j)))/
     .                 (gam_D))
            voigt_y =(gam_L)!*sqrt(log(2.0))
            Z = cmplx(voigt_x, voigt_y)
            hui = w4(Z)
c            call complexVoigt(voigt_x,voigt_y,
c     .                        new_voigt, new_fv)
            new_voigt = REAL(hui)/(1.772453851)
            new_fv = AIMAG(hui)/(1.772453851)
            if (width(j) .eq. 0.0) then
                kapnu_I(i) = kapnu_I(i) + kapnu0(j,i)*new_voigt*
     .             (sin(phi)**2.0)/2.0
                kapnu_Q(i) = kapnu_Q(i) + kapnu0(j,i)*new_voigt*
     .             (sin(phi)**2.0)/2.0
                zetnu_Q(i) = zetnu_Q(i) + kapnu0(j,i)*new_fv*
     .             (sin(phi)**2.0)/2.0
            else
                kapnu_I(i) = kapnu_I(i) + kapnu0(j,i)*new_voigt*
     .             (1.0+cos(phi)**2.0)/4.0
                kapnu_Q(i) = kapnu_Q(i) - kapnu0(j,i)*new_voigt*
     .             (sin(phi)**2.0)/4.0
                kapnu_V(i) = kapnu_V(i) + kapnu0(j,i)*new_voigt*
     .             (cos(phi))/2.0*width(j)
                zetnu_Q(i) = zetnu_Q(i) - kapnu0(j,i)*new_fv*
     .             (sin(phi)**2.0)/4.0
                zetnu_V(i) = zetnu_V(i) + kapnu0(j,i)*new_fv*
     .             (cos(phi))/2.0*width(j)
            endif
         enddo                                     

                                                       
c*****do the same for the strong lines
         if (dostrong .gt. 0) then
            do j=nlines+1,nlines+nstrong
c               v = 2.997929d10*dabs(wave-wave1(j))/
c     .             (wave1(j)*dopp(j,i)) 
c               old_voigt = voigt(a(j,i),v)
               gam_L = a(j,i)
               gam_D = dopp(j,i)/(2.997929d2*wave)
               voigt_x = !sqrt(log(2.0))*
     .                    dabs((1.0d8*(1./wave-1./wave1(j)))/
     .                    (gam_D))
               voigt_y =(gam_L)!*sqrt(log(2.0))
               Z = cmplx(voigt_x, voigt_y)
               hui = w4(Z)
c               call complexVoigt(voigt_x,voigt_y,
c     .                           new_voigt, new_fv)
               new_voigt = REAL(hui)/(1.772453851)
               new_fv = AIMAG(hui)/(1.772453851)
c               write (*,*) i, old_voigt, new_voigt
               if (width(j) .eq. 0.0) then
                   kapnu_I(i) = kapnu_I(i)+ kapnu0(j,i)*new_voigt*
     .                 (sin(phi)**2.0)/2.0
                   kapnu_Q(i) = kapnu_Q(i)+ kapnu0(j,i)*new_voigt*
     .                 (sin(phi)**2.0)/2.0
c                   kapnu_U(i) = kapnu_U(i)+ kapnu0(j,i)*new_voigt*
c     .                 (sin(phi)**2.0)*sin(2.0*zeta)/2.0
                   zetnu_Q(i) = zetnu_Q(i) + kapnu0(j,i)*new_fv*
     .                 (sin(phi)**2.0)/2.0
c                   zetnu_U(i) = zetnu_U(i) + kapnu0(j,i)*new_fv*
c     .                 (sin(phi)**2.0)*sin(2.0*zeta)/2.0
               else
                   kapnu_I(i) = kapnu_I(i)+ kapnu0(j,i)*new_voigt*
     .                 (1.0+cos(phi)**2.0)/4.0
                   kapnu_Q(i) = kapnu_Q(i)- kapnu0(j,i)*new_voigt*
     .                 (sin(phi)**2.0)/4.0
c                   kapnu_U(i) = kapnu_U(i)- kapnu0(j,i)*new_voigt*
c     .                 (sin(phi)**2.0)*sin(2.0*zeta)/4.0
                   kapnu_V(i) = kapnu_V(i)-kapnu0(j,i)*new_voigt*
     .                 (cos(phi))/2.0*width(j)
                   zetnu_Q(i) = zetnu_Q(i) - kapnu0(j,i)*new_fv*
     .                 (sin(phi)**2.0)/4.0
c                   zetnu_U(i) = zetnu_U(i) - kapnu0(j,i)*new_fv*
c     .                 (sin(phi)**2.0)*sin(2.0*zeta)/4.0
                   zetnu_V(i) = zetnu_V(i) + kapnu0(j,i)*new_fv*
     .                 (cos(phi))/2.0*width(j)
               endif
            enddo
         endif
         eta_I(i) = kapnu_I(i)/(0.4343*kaplam(i))
         eta_Q(i) = kapnu_Q(i)/(0.4343*kaplam(i))
c         eta_U(i) = kapnu_U(i)/(0.4343*kaplam(i))
         eta_V(i) = kapnu_V(i)/(0.4343*kaplam(i))
         zet_Q(i) = zetnu_Q(i)/(0.4343*kaplam(i))
c         zet_U(i) = zetnu_U(i)/(0.4343*kaplam(i))
         zet_V(i) = zetnu_V(i)/(0.4343*kaplam(i))
c         eta_I(i) = kapnu_I(i)/(kaplam(i))!/(0.4343*kaplam(i))
c         eta_Q(i) = kapnu_Q(i)/(kaplam(i))!/(0.4343*kaplam(i))
cc         eta_U(i) = kapnu_U(i)/(0.4343*kaplam(i))
c         eta_V(i) = kapnu_V(i)/(kaplam(i))!/(0.4343*kaplam(i))
c         zet_Q(i) = zetnu_Q(i)/(kaplam(i))!/(0.4343*kaplam(i))
cc         zet_U(i) = zetnu_U(i)/(0.4343*kaplam(i))
c         zet_V(i) = zetnu_V(i)/(kaplam(i))!/(0.4343*kaplam(i))
      enddo      

c*****compute the optical depths                                            
c      first = tauref(1)*kapnu(1)/kapref(1)
c      dummy1(1) = rinteg(xref,dummy1,taunu,ntau,0.)
c      taunu(i) = first
c      do i=2,ntau                                                     
c          taunu(i) = taunu(i-1)+taunu(i)
c      enddo

      return                                              
321   format (f11.3, e11.3, f11.1, 5e11.3)
c322   format (f11.3, 4e11.3)
      end                                                


