# PD2 intact boundary reintegration figures

Seven figures from recorded numerical data: six required figures and the classical-shear control. English labels, 600 dpi PNG and vector PDF.

[Scientific review](../../verification/higher_order_baqp_boundary_reintegration_v1/REVIEW_HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1.md) · [Data dictionary](../../verification/higher_order_baqp_boundary_reintegration_v1/DATA_DICTIONARY.md)

## 1. Exact integration architecture

Preserved PD2 internal mechanics and the external, frozen boundary conversion. The diagram describes the implementation; it is not a fit to a traction layer.

[PNG](figures/pd2_boundary_reintegration_overview.png) · [PDF](figures/pd2_boundary_reintegration_overview.pdf)

## 2. Affine force, moment, work and solved fields

Both affine solves start from a nonzero perturbation. Errors and predeclared gates use the historical normalization; Psi+ is an unchanged diagnostic.

[PNG](figures/pd2_affine_boundary_consistency.png) · [PDF](figures/pd2_affine_boundary_consistency.pdf)

## 3. Localization and proven closed-center checks

N19 layer maxima/means and four N27 lateral-depth samples are different sampling sets. Only the 13h N27 center has a fully closed dependency graph. The dashed threshold applies to that closed center.

[PNG](figures/pd2_nonaffine_boundary_localization.png) · [PDF](figures/pd2_nonaffine_boundary_localization.pdf)

## 4. C1 versus PD2 conversion support

Magnitude-weighted d50/d90/d95/d99 and cumulative load fractions. A zero PD2 d50 means at least half the magnitude lies at the physical surface. These are conversion-support measures, not physical BC thickness.

[PNG](figures/c1_vs_pd2_boundary_load_support.png) · [PDF](figures/c1_vs_pd2_boundary_load_support.pdf)

## 5. Near-boundary field quality

Accepted historical C1 target-state diagnostics versus newly solved PD2 nonlocal-manufactured states. Node and bond-midpoint references use matched sampling locations. Normalized error in gradients and unchanged Psi+ diagnostic; no C1 campaign was rerun.

[PNG](figures/c1_vs_pd2_boundary_field_quality.png) · [PDF](figures/c1_vs_pd2_boundary_field_quality.pdf)

## 6. Four required perturbed-start recoveries

Target and solved gradients are compared with analytic fields, with actual mechanical residuals. Classical-loaded shear is intentionally separate because its local target has a nonzero bulk-model residual.

[PNG](figures/pd2_manufactured_recovery.png) · [PDF](figures/pd2_manufactured_recovery.pdf)

## 7. Classical shear control

Maximum deviation from the analytic target at each Z plane for displacement, nodal F and nodal Psi+. The classical load retains e_model=q_bulk+b_local. The NL recovery curve exposes numerical resolution; tiny shifts are not all assigned physical meaning.

[PNG](figures/pd2_classical_vs_nonlocal_shear.png) · [PDF](figures/pd2_classical_vs_nonlocal_shear.pdf)
