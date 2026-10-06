# Nonaffine shear boundary localization: figures

Every figure uses the saved numerical data in [verification](../../verification/nonlocal_boundary_nonaffine_shear_localization_v1/). English text; 600 dpi PNG and vector PDF. Reproducible plotter: [plot_nonaffine_shear_boundary.py](../../scripts/plot_nonaffine_shear_boundary.py).

| Figure | PNG | PDF | Evidence / reading |
| --- | --- | --- | --- |
| 1. Residual decomposition | [PNG](figures/nonaffine_shear_residual_decomposition.png) | [PDF](figures/nonaffine_shear_residual_decomposition.pdf) | N19 plane-mean x/z force components; lower panels resolve independently measured model load and corrected residual |
| 2. Naive versus localized conversion | [PNG](figures/shear_naive_vs_localized_conversion.png) | [PDF](figures/shear_naive_vs_localized_conversion.pdf) | Nearest-face mean/max magnitudes on N19, plus independent N27 bulk-depth diagnostic |
| 3. Conversion support | [PNG](figures/shear_boundary_localization_support.png) | [PDF](figures/shear_boundary_localization_support.pdf) | Cumulative nodal force magnitudes and d50/d90/d95/d99; not physical BC thickness |
| 4. N27 boundary decay | [PNG](figures/n27_nonaffine_shear_boundary_decay.png) | [PDF](figures/n27_nonaffine_shear_boundary_decay.pdf) | Z=0 lateral samples, independent N29 reference, 1e-10 gate, explicitly marked fully closed center |
| 5. Perturbed-start recovery | [PNG](figures/nonaffine_shear_manufactured_recovery.png) | [PDF](figures/nonaffine_shear_manufactured_recovery.pdf) | NL and classical loads, same perturbation/constraints/frozen conversion, discrete-target Fbar recovery and solver residuals |
| 6. Kinematic and energy bias | [PNG](figures/nonaffine_shear_kinematic_energy_bias.png) | [PDF](figures/nonaffine_shear_kinematic_energy_bias.pdf) | Fbar at nodes, Ftilde at bond midpoints, unchanged spectral Psi+; max/weighted RMS versus corresponding depth |

The scientific conclusions and scope are in [REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md](../../verification/nonlocal_boundary_nonaffine_shear_localization_v1/REVIEW_NONAFFINE_SHEAR_LOCALIZATION_V1.md). Units, signs, weighting, reference moments, and normalization are documented in the [data dictionary](../../verification/nonlocal_boundary_nonaffine_shear_localization_v1/DATA_DICTIONARY.md).

Boundary-force conversion accuracy and analytic-gradient reconstruction accuracy are distinct diagnostics. This round evaluates intact mechanics and reports Psi+ bias without a fracture acceptance rating. No physical traction-layer thickness is inferred from cumulative conversion support.
