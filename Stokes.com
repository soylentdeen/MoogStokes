
c******************************************************************************
c     this common block has variables related to the lines.  Most
c     input quantities typically have single dimensions, while the 
c     things that are computed for each line at each atmosphere level
c     have double dimensions.  The variables "a", "dopp", and 
c     "kapnu0" are often over-written with plotting data,
c     so leave them alone or suffer unspeakable programming tortures.
c******************************************************************************

      real*8       deltamj(2500),phi_opacity(100, 3),psi_opacity(100,3),
     .             Stokes(4), continuum, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V,
     .             kappa(4,4,100), emission(4,100), kaptot(100),
     .             phiI(100), phiQ(100), phiU(100), phiV(100),
     .             psiQ(100), psiU(100), psiV(100), dphiI(100),
     .             dphiQ(100), dphiU(100), dphiV(100), dpsiQ(100),
     .             dpsiU(100), dpsiV(100),
     .             de1(100), de2(100), de3(100), de4(100),
     .             eta0(100), deta0(100), tlam(100), dktot(100),
     .             dkref(100),
     .             dttot(100), dtlam(100), tautot(100), dklam(100),
     .             zdepth(100)

      common/stokesparams/deltamj, phi_opacity, psi_opacity,
     .             Stokes, continuum, total_weight,
     .             Stokes_I, Stokes_Q, Stokes_U, Stokes_V,
     .             kappa, emission, kaptot,
     .             phiI, phiQ, phiU, phiV,
     .             psiQ, psiU, psiV, dphiI,
     .             dphiQ, dphiU, dphiV, dpsiQ,
     .             dpsiU, dpsiV,
     .             de1, de2, de3, de4,
     .             eta0, deta0, tlam, dktot,
     .             dkref,
     .             dttot, dtlam, tautot, dklam,
     .             zdepth

