
c******************************************************************************
c     this common block has variables related to the lines.  Most
c     input quantities typically have single dimensions, while the 
c     things that are computed for each line at each atmosphere level
c     have double dimensions.  The variables "a", "dopp", and 
c     "kapnu0" are often over-written with plotting data,
c     so leave them alone or suffer unspeakable programming tortures.
c******************************************************************************

      real*8       deltamj(2500),phi_opacity(100, 3),psi_opacity(100,3),
     .             Stokes(4), continuum, wavelength, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V,
     .             kappa(4,4,100), emission(4,100), kaptot(100),
     .             k11(100), k12(100), k13(100), k14(100),
     .             k21(100), k22(100), k23(100), k24(100),
     .             k31(100), k32(100), k33(100), k34(100),
     .             k41(100), k42(100), k43(100), k44(100),
     .             dk11(100), dk12(100), dk13(100), dk14(100),
     .             dk21(100), dk22(100), dk23(100), dk24(100),
     .             dk31(100), dk32(100), dk33(100), dk34(100),
     .             dk41(100), dk42(100), dk43(100), dk44(100),
     .             de1(100), de2(100), de3(100), de4(100),
     .             eta0(100), deta0(100), tlam(100), dktot(100),
     .             dkref(100),
     .             dttot(100), dtlam(100), tautot(100), dklam(100),
     .             zdepth(100)

      common/stokesparams/deltamj, phi_opacity, psi_opacity,
     .             Stokes, continuum, wavelength, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V,
     .             kappa, emission, kaptot,
     .             k11, k12, k13, k14,
     .             k21, k22, k23, k24,
     .             k31, k32, k33, k34,
     .             k41, k42, k43, k44,
     .             dk11, dk12, dk13, dk14,
     .             dk21, dk22, dk23, dk24,
     .             dk31, dk32, dk33, dk34,
     .             dk41, dk42, dk43, dk44,
     .             de1, de2, de3, de4,
     .             eta0, deta0, tlam, dktot,
     .             dkref,
     .             dttot, dtlam, tautot, dklam,
     .             zdepth

