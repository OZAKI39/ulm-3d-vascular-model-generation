# English figure captions

## fig01 — From bond activation to interpretable clot mechanics

Evidence map for the current clot-model development candidate. Quantities are taken from archived C5/C7 activation and branch audits and the finite-domain BAQP audit. Arrows organize the research argument; they do not represent a validated causal or biological model. Neither bond activation nor a converged positive-phase state establishes macroscopic clot rupture.

[PDF](output/figures/fig01_research_evidence_overview.pdf) · [SVG](output/figures/fig01_research_evidence_overview.svg) · [PNG](output/figures/fig01_research_evidence_overview.png)

## fig02 — Activation shifts with discretization

C5 and C7 activation is compared for the same corrected physical cube and frozen candidate model. Markers show the diagnostic weak-state strain, with numerical brackets shown as intervals rather than statistical confidence intervals; these intervals are smaller than the markers at this scale. Both first sets contain two bonds, so the smaller C7 fraction partly reflects its larger bond denominator. Activation remains boundary-associated, and two grids do not establish convergence or a crack prediction.

[PDF](output/figures/fig02_activation_onset_and_sampling.pdf) · [SVG](output/figures/fig02_activation_onset_and_sampling.svg) · [PNG](output/figures/fig02_activation_onset_and_sampling.png)

## fig03 — End-layer concentration precedes activation

Archived intact equilibria at engineering strain 0.003 show the end-layer displacement increment and the associated centerline bond strain. The affine reference is the prescribed engineering strain, while driver amplification is measured relative to the archived affine energy density. Releasing lateral clamp components does not reduce the recorded peak amplification. These observations support a boundary-sensitive mechanical precursor, but they do not prove a defective operator or identify a physical clot-failure mechanism.

[PDF](output/figures/fig03_preactivation_end_layer_concentration.pdf) · [SVG](output/figures/fig03_preactivation_end_layer_concentration.svg) · [PNG](output/figures/fig03_preactivation_end_layer_concentration.png)

## fig04 — Affine deformation requires distributed nodal forces

Signed internal layer forces are shown for the exact affine state, together with cumulative manufactured external force from the top boundary. All forces are normalized by the matching continuum reaction Pzz L squared. The cumulative L075 curve overshoots before approaching its final half-reaction, so a first threshold crossing is not a confinement thickness. Manufactured force is the negative of the computed internal force and is not an independently validated physical loading law.

[PDF](output/figures/fig04_affine_layer_forces_and_cumulative_loading.pdf) · [SVG](output/figures/fig04_affine_layer_forces_and_cumulative_loading.svg) · [PNG](output/figures/fig04_affine_layer_forces_and_cumulative_loading.png)

## fig05 — First-layer force contains several contributions

Exact algebraic decompositions of the bottom first-interior center-particle force separate direct and Z terms, loaded-surface versus other neighbors, and endpoint-normalization sources. Contributions are divided by the signed net axial force; negative values and values above unity represent cancellation. Z supplies approximately 78 to 82 percent of the net force, while the loaded-surface endpoint normalization is small. All components sum to the original model, and no alternative constitutive law or boundary correction is tested here.

[PDF](output/figures/fig05_first_layer_force_origins.pdf) · [SVG](output/figures/fig05_first_layer_force_origins.svg) · [PNG](output/figures/fig05_first_layer_force_origins.png)

## fig06 — A full family does not isolate the force stencil

Saved affine nodal force magnitudes are plotted against distance to the nearest surface for fixed spacing and horizon. Filled symbols distinguish particles with a geometrically complete family from particles whose family is truncated; the dashed distance marker denotes one horizon. L075 has 27 full-family particles, but its maximum core residual exceeds both declared patch-test criteria and it has no fully closed two- or three-hop core. This failed geometric test does not by itself establish a boundary-independent bulk defect.

[PDF](output/figures/fig06_finite_domain_support_and_residuals.pdf) · [SVG](output/figures/fig06_finite_domain_support_and_residuals.svg) · [PNG](output/figures/fig06_finite_domain_support_and_residuals.png)

## fig07 — Energy consistency and traction recovery differ

Global affine energy errors are compared with the observed half-domain manufactured reaction and the affine virtual-work ratio. The energy and virtual work agree with their continuum references to floating-point accuracy, whereas C5 and C7 half-reactions differ by about 38 and 44 percent. Distributed forces act at different moment arms from surface traction, so energy agreement does not force their half-resultants to match. The one-percent reaction band is the archived audit criterion, not measurement uncertainty.

[PDF](output/figures/fig07_energy_virtual_work_and_half_reaction.pdf) · [SVG](output/figures/fig07_energy_virtual_work_and_half_reaction.svg) · [PNG](output/figures/fig07_energy_virtual_work_and_half_reaction.png)

## fig08 — Positive phase does not establish rupture

Archived positive-phase states show maximum phase, quadrature-weighted mean phase, nominal reaction stress and active-bond counts. Phase and count curves use strain above each grid's diagnostic onset; stress uses absolute engineering strain. Filled C5 markers denote 13 audited stable positive states, while open C7 markers denote 20 accepted stationary states whose stability was not assessed. Connecting segments only guide the eye between saved states; the unequal windows and localized support do not imply convergence or physical crack propagation.

[PDF](output/figures/fig08_positive_phase_branch.pdf) · [SVG](output/figures/fig08_positive_phase_branch.svg) · [PNG](output/figures/fig08_positive_phase_branch.png)

## fig09 — Stationarity and stability require distinct evidence

The C5 mechanical and condensed phase curvatures are shown with the archived curvature threshold, alongside force and projected phase KKT residuals normalized by their acceptance tolerances. All 15 C5 audited states are classified as stationary and stable; the intact state has no movable phase block, and its Schur value is therefore absent. C7 residuals establish stationarity of the 20 accepted states, but no C7 stability calculation exists. No instability was identified in the C5 audited window; this is not a statement about every intervening state or later loads.

[PDF](output/figures/fig09_stationarity_and_stability_evidence.pdf) · [SVG](output/figures/fig09_stationarity_and_stability_evidence.svg) · [PNG](output/figures/fig09_stationarity_and_stability_evidence.png)

## fig10 — Implications for future clot simulation

The current numerical evidence is mapped to interpretation limits and development implications for a future clot model. Activation fractions depend on bond sampling, positive-phase states need not be unstable, and finite-domain loading consistency remains unresolved. These findings justify careful validation of boundary loading and failure observables before interpreting simulated phase as clot rupture. The diagram summarizes development implications only; no clot-specific calibration, biological mechanism or patient-level prediction is claimed.

[PDF](output/figures/fig10_implications_for_clot_simulation.pdf) · [SVG](output/figures/fig10_implications_for_clot_simulation.svg) · [PNG](output/figures/fig10_implications_for_clot_simulation.png)
