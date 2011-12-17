
      subroutine cdcalc (number)                                        
c******************************************************************************
c     this routine calculates continuum and line+continuum                  
c     contribution functions
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      real*8 prev_tau_I, prev_tau_Q, prev_tau_V, prev_tau_lam
      real*8 StokesI, StokesQ, StokesV
      real*8 exptau_I, exptau_Q, exptau_V
      real*8 dTau_I, dTau_Q, dTau_V
      real*8 dQ, dV, junk

c*****continuum "contribution curve" calculation                           
      if (number .eq. 1) then
         do i=1,ntau                                                    
            scont(i) = ((1.19089d+25/wave**2)*1.0d+10)/(wave**3*
     .                 (dexp(1.43879d+08/(wave*t(i)))-1.0d+00))
            if (fluxintopt .eq. 1) then
               cd(i) = kaplam(i)*tauref(i)*scont(i)/mu*
     .                 dexp(-taulam(i)/mu)/(0.4343*kapref(i))
            else
               cd(i) = 2.0*kaplam(i)*tauref(i)*scont(i)*
     .                 expint(taulam(i),2)/(0.4343*kapref(i))
            endif
         enddo
                       
c*****line plus continuum "contribution curve" calculation    
      elseif (number .eq. 2) then
         do i=1,ntau                                                    
            sline(i) = scont(i)   
c            write (*,*) 'Tau nu: ', taunu(i), taulam(i)
            if (fluxintopt .eq. 1) then
               if ((taulam(i)+taunu(i))/mu .le. 50.) then
                  exptau = dexp((-taulam(i)-taunu(i))/mu)
               else
                  exptau = 0.                        
               endif
               cd(i) = tauref(i)*kaplam(i)/(mu*0.4343*flux*
     .                 kapref(i))*(scont(i)*dexp(-taulam(i)/mu) - 
     .                 (1.0+kapnu(i)/kaplam(i))*sline(i)*exptau)    
            else
               cd(i) = 2.0*tauref(i)*kaplam(i)/(0.4343*flux*
     .                 kapref(i))*(scont(i)*expint(taulam(i),2) - 
     .                 (1.0+kapnu(i)/kaplam(i))*sline(i)*
     .                 expint(taulam(i)+taunu(i),2)) 
            endif
         enddo
      endif

c*****normal exit
      return


      end                                                               
