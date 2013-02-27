
c******************************************************************************
c     this common block has variables related to the lines.  Most
c     input quantities typically have single dimensions, while the 
c     things that are computed for each line at each atmosphere level
c     have double dimensions.  The variables "a", "dopp", and 
c     "kapnu0" are often over-written with plotting data,
c     so leave them alone or suffer unspeakable programming tortures.
c******************************************************************************

      real*8       deltamj(2500),crad(2500),c4(2500),phi_opacity(100,3),
     .             psi_opacity(100,3), Stokes(4), continuum,
     .             emission(100), kref_knots(100), kref_coeffs(100),
     .             klam_knots(100), klam_coeffs(100),
     .             kaptot(100), ktot_knots(100), ktot_coeffs(100),
     .             phiI(100), phiI_knots(100), phiI_coeffs(100),
     .             phiQ(100), phiQ_knots(100), phiQ_coeffs(100),
     .             phiU(100), phiU_knots(100), phiU_coeffs(100),
     .             phiV(100), phiV_knots(100), phiV_coeffs(100),
     .             psiQ(100), psiQ_knots(100), psiQ_coeffs(100),
     .             psiU(100), psiU_knots(100), psiU_coeffs(100),
     .             psiV(100), psiV_knots(100), psiV_coeffs(100),
     .             e_knots(100), e_coeffs(100), taus(10000),
     .             tlam(10000), tlam_knots(10000), tlam_coeffs(10000),
     .             ttot(10000), ttot_knots(10000), ttot_coeffs(10000),
     .             zdepth(10000), z_knots(10000), z_coeffs(10000)
      integer      nz, n_phiI_knots, n_phiQ_knots, n_phiU_knots,
     .             n_phiV_knots, n_psiQ_knots, n_psiU_knots,
     .             n_psiV_knots, n_kref_knots, n_klam_knots,
     .             n_ktot_knots, n_e_knots, n_z_knots,
     .             n_tlam_knots, n_ttot_knots

      common/stokesparams/deltamj, crad, c4, phi_opacity,
     .             psi_opacity, Stokes, continuum,
     .             emission, kref_knots, kref_coeffs,
     .             klam_knots, klam_coeffs,
     .             kaptot, ktot_knots, ktot_coeffs,
     .             phiI, phiI_knots, phiI_coeffs,
     .             phiQ, phiQ_knots, phiQ_coeffs,
     .             phiU, phiU_knots, phiU_coeffs,
     .             phiV, phiV_knots, phiV_coeffs,
     .             psiQ, psiQ_knots, psiQ_coeffs,
     .             psiU, psiU_knots, psiU_coeffs,
     .             psiV, psiV_knots, psiV_coeffs,
     .             e_knots, e_coeffs, taus,
     .             tlam, tlam_knots, tlam_coeffs,
     .             ttot, ttot_knots, ttot_coeffs,
     .             zdepth, z_knots, z_coeffs,
     .             nz, n_phiI_knots, n_phiQ_knots, n_phiU_knots,
     .             n_phiV_knots, n_psiQ_knots, n_psiU_knots,
     .             n_psiV_knots, n_kref_knots, n_klam_knots,
     .             n_ktot_knots, n_e_knots, n_z_knots,
     .             n_tlam_knots, n_ttot_knots
