# Nonaffine boundary localization figures

All curves come from the saved numerical arrays of `NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1`. No synthetic data or image generation is used. Each figure is supplied as a 600 dpi PNG and a vector PDF, with English labels.

| Figure | Main question | Files |
| --- | --- | --- |
| Nonaffine residual decomposition | Which residual belongs to boundary conversion, and which remains in the bulk? | [PNG](figures/nonaffine_residual_decomposition.png) · [PDF](figures/nonaffine_residual_decomposition.pdf) |
| Naive versus localized conversion | Why does direct residual cancellation remove interior mechanics? | [PNG](figures/naive_vs_localized_boundary_conversion.png) · [PDF](figures/naive_vs_localized_boundary_conversion.pdf) |
| Cumulative conversion support | How is correction magnitude distributed with physical boundary distance? | [PNG](figures/boundary_localization_cumulative_support.png) · [PDF](figures/boundary_localization_cumulative_support.pdf) |
| N27 depth check | Does correction vanish at a verified fully closed center? | [PNG](figures/n27_nonaffine_boundary_decay.png) · [PDF](figures/n27_nonaffine_boundary_decay.pdf) |
| Manufactured recovery | What happens under nonlocal and classical manufactured body forces? | [PNG](figures/manufactured_solution_recovery.png) · [PDF](figures/manufactured_solution_recovery.pdf) |

The support quantiles are **conversion support, not physical BC thickness**. N19 has no fully closed center stencil. N27 points approach the center from a lateral x face at Z=0; full graph closure is explicitly verified only at its 13h center. Small correction at 12h does not prove closure.

Recovery plots compare the solved and target **discrete reconstructed gradients**. The analytic-gradient reconstruction bias is reported separately. The classical target has a very small nonzero bulk-model residual; solved displacement differences must be interpreted alongside solver precision.

Inputs: [verification directory](../../verification/nonlocal_boundary_nonaffine_localization_v1). Interpretation and limitations: [review](../../verification/nonlocal_boundary_nonaffine_localization_v1/REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md). Plot source: [plot_nonaffine_boundary.py](../../scripts/plot_nonaffine_boundary.py).

No optional shear, mixed traction, curved clot, surface correction, phase, fracture, C5/C7 activation or CFD is run in this round.
