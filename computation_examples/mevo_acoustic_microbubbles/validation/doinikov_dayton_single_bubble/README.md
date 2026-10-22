# Doinikov–Dayton reference, DOI 10.1121/1.2215228

Full original author manuscript read from [PMC3063887](https://pmc.ncbi.nlm.nih.gov/articles/PMC3063887/), including Eqs. 31–40 and the Fig. 8 20-cycle/150-kPa case. Author-code searches (web, GitHub repositories, article supplement links) did not locate a callable implementation. `results/stage2/translation/CODE_SEARCH.json` records the repository searches; absence in these searches is not proof no private code exists. APECSS and PROTEUS do not provide this finite-thickness coupled solid-shell model.

The minimum four-state reference in `src/mevo_mb/validation/doinikov_dayton.py` implements inner radius, its velocity, centre displacement and velocity using SciPy DOP853. Production radial response remains the APECSS wrapper. This reference must not be used as a SonoVue material model.

- Eq. 6: constant shell volume and outer/inner radius relation.
- Eqs. 40 and 33: unstrained radius and equilibrium gas pressure.
- Eq. 31: radial inertia, finite solid-shell elastic stress, viscosity and first-order compressibility.
- Eq. 32: translational inertia, variable added mass, instantaneous acoustic pressure-gradient force, no-slip Stokes drag.
- Eq. 36 eliminates the third derivative. Differentiate the smooth algebraic RHS G by complex step, collect its `G_u Rddot` contribution and the explicit `Rddot` term before division. Units of G are m²/s²; physical-time acceleration is converted to dimensionless time only at the integrator boundary.
- Eq. 63, with Eqs. 44–49, independently checks the mean force in the weak-field limit.

Figure 8's albumin-shell liquid/shell constants are transcribed from Section III. Additional assumptions not explicitly fixed there: air exponent 1.4, gas density 1.2 kg/m³, sinusoidal pulse phase zero, no ambient acoustic particle velocity, no post-pulse displacement. These prevent claiming pixel-level reproduction. The model excludes history force, nonlinear shell rheology/buckling, dissolution and finite-Reynolds drag. No-slip drag is retained exactly as Eq. 32.

`scripts/stage2/translation_validate.py` runs pressure/radius/frequency comparisons and a time-step refinement. The fast comparison uses the *same computed cycle force*, divided by `6*pi*mu*R0`, isolating the translation approximation from radial-model changes. Initial oscillatory transients can dominate the tiny net displacement at very low pressure; a negative short-pulse displacement is not relabeled positive. Measured errors remain in the report.

For monochromatic pressure gradients, `-<V grad(p)>` equals the fundamental-volume formula even for nonlinear V(t), by Fourier orthogonality. Tests also add a second pressure harmonic and recover its additional force explicitly. No spurious nonlinear discrepancy is introduced.

## Contemporary check

Cattaneo & Supponen, **Physical Review Fluids 11, 033606**, published 16 March 2026, DOI [10.1103/p475-xggm](https://doi.org/10.1103/p475-xggm); not Physics of Fluids. Full [author manuscript v2](https://arxiv.org/html/2506.16413v2) was read. The authors use a thermal gas model, no-slip variable-radius history force and finite-Re drag, and report better short-pulse agreement when history is included. Their scaling and stability evidence is specific to their experiments; it does not validate transport in blood or indefinite stable CW exposure. Source data are openly available as [Zenodo 15692129](https://doi.org/10.5281/zenodo.15692129), CC BY 4.0; the record lists video data and a figure-data spreadsheet, not a solver. The classic 2006 comparison does not remove this history-force limitation.
