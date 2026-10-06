# Figure plan

All numeric values will be loaded from the listed existing files. No solver or new scientific experiment is invoked.

| ID | Title | Scientific message | Type | Clot relevance |
|---|---|---|---|---|
| fig01 | From bond activation to interpretable clot mechanics | Observed activation, numerical admissibility and physical interpretation are separate claims. | Evidence-linked schematic | Places current coupon evidence in the future clot-model validation chain. |
| fig02 | Activation shifts with discretization | The first set contains two bonds on both grids, while its fraction and onset strain change. | Onset intervals, bond fractions and normalized location | A raw activated-bond fraction is not a mesh-independent failure measure. |
| fig03 | End-layer concentration precedes activation | Similar end-layer displacement increments become larger local strain as spacing decreases. | Saved centerline bond increments and loading-control comparison | Identifies a mechanical precursor that must be separated from a clot-specific damage mechanism. |
| fig04 | Affine deformation requires distributed nodal forces | Exact affine kinematics does not imply equilibrium under surface-only constraints. | Signed layer forces and cumulative required force | Shows why prescribed displacement and physically interpreted traction need separate validation. |
| fig05 | First-layer force contains several contributions | Z correction is larger, but neither a large term nor a signed endpoint contribution proves a bias. | Signed algebraic contribution bars | Avoids treating a force-formula component as a physical damage mechanism. |
| fig06 | A full family does not isolate the force stencil | The geometric core gate fails; complete multi-hop force support is absent. | Particle force scatter and support counts | A boundary-free mechanical oracle is needed before interpreting local activation as failure. |
| fig07 | Energy consistency and traction recovery differ | Global affine energy and virtual work agree while small-domain half-reactions differ. | Energy error and reaction/virtual-work comparison | Defines which loading consistency checks are prerequisites for interpretable clot response. |
| fig08 | Positive phase does not establish rupture | C5 has audited stable states; C7 has accepted stationary states with stability unassessed. | Maximum/weighted mean phase, stress and active counts | Separates phase amplitude growth from spatial propagation and macroscopic rupture. |
| fig09 | Stationarity and stability require distinct evidence | All 15 C5 audited states are stable; C7 stationarity alone does not classify stability. | Dimensionless curvature and tolerance-normalized residuals | Prevents solver convergence from being mislabeled as a fracture event. |
| fig10 | Implications for future clot simulation | Mechanical consistency, a physical failure observable and clot calibration remain separate requirements. | Evidence / interpretation limit / development implication matrix | Explicitly states what the current evidence permits and what remains to be established. |

## Exact sources

### fig01 — From bond activation to interpretable clot mechanics

- `verification/pfpd_at1_c7_early_branch_v1/activation/C5_C7_COMPARISON.json`
- `verification/pfpd_at1_c5_stability_v1/SUMMARY.json`
- `verification/pfpd_at1_c7_early_branch_v1/branch/ALL.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/SUMMARY.json`
- `verification/pfpd_at1_c5_activation_stationarity_audit_v1/activation/BRACKET.json`
- `verification/pfpd_at1_c7_early_branch_v1/activation/BRACKET.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/states/S1/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S2/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S3/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S4/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S5/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S6/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S7/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E01/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E02/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E03/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E04/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E05/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E06/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`
- `../PROJECT_ACTIVE_CONTEXT.md`

### fig02 — Activation shifts with discretization

- `verification/pfpd_at1_c7_early_branch_v1/activation/C5_C7_COMPARISON.json`
- `verification/pfpd_at1_c5_activation_stationarity_audit_v1/activation/BRACKET.json`
- `verification/pfpd_at1_c7_early_branch_v1/activation/BRACKET.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`

### fig03 — End-layer concentration precedes activation

- `verification/pfpd_at1_activation_root_cause_and_fix_v1/SUMMARY.json`
- `verification/pfpd_at1_activation_root_cause_and_fix_v1/CONFIG.json`
- `results/pfpd_at1_activation_root_cause_and_fix_v1/tables/C5_FULL_CLAMP_BASELINE_0.0030_CENTERLINE_GAPS.csv`
- `results/pfpd_at1_activation_root_cause_and_fix_v1/tables/C7_FULL_CLAMP_BASELINE_0.0030_CENTERLINE_GAPS.csv`
- `results/pfpd_at1_activation_root_cause_and_fix_v1/tables/C5_FULL_CLAMP_BASELINE_0.0030_AXIAL_BONDS.csv`
- `results/pfpd_at1_activation_root_cause_and_fix_v1/tables/C7_FULL_CLAMP_BASELINE_0.0030_AXIAL_BONDS.csv`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`

### fig04 — Affine deformation requires distributed nodal forces

- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C5_layers.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C7_layers.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C5_cumulative.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C7_cumulative.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_cumulative.csv`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`

### fig05 — First-layer force contains several contributions

- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/NORMALIZATION_CHANNELS.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`

### fig06 — A full family does not isolate the force stencil

- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C5_particles.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L050_particles.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_particles.csv`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/CONFIG.json`

### fig07 — Energy consistency and traction recovery differ

- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/VIRTUAL_WORK_TRACTION.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/CONFIG.json`

### fig08 — Positive phase does not establish rupture

- `verification/pfpd_at1_c7_early_branch_v1/activation/C5_C7_COMPARISON.json`
- `verification/pfpd_at1_c5_stability_v1/SUMMARY.json`
- `verification/pfpd_at1_c7_early_branch_v1/branch/ALL.json`
- `verification/pfpd_at1_c5_stability_v1/states/S1/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S2/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S3/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S4/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S5/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S6/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S7/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E01/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E02/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E03/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E04/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E05/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E06/SUMMARY.json`
- `verification/pfpd_at1_c5_activation_stationarity_audit_v1/activation/BRACKET.json`
- `verification/pfpd_at1_c7_early_branch_v1/activation/BRACKET.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `pd_clot/pfpd_candidate/coupling_stability_audit.py`
- `pd_clot/pfpd_candidate/adr_adapter_at1.py`

### fig09 — Stationarity and stability require distinct evidence

- `verification/pfpd_at1_c5_stability_v1/SUMMARY.json`
- `verification/pfpd_at1_c7_early_branch_v1/branch/ALL.json`
- `verification/pfpd_at1_c5_stability_v1/provenance/PREDECLARED_CONFIG.json`
- `verification/pfpd_at1_petsc_vi_solver_v1/provenance/PREDECLARED_CONFIG.json`
- `verification/pfpd_at1_c7_early_branch_v1/provenance/PREDECLARED_CONFIG.json`
- `verification/pfpd_at1_c7_early_branch_v1/activation/C5_C7_COMPARISON.json`
- `verification/pfpd_at1_c5_stability_v1/states/S1/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S2/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S3/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S4/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S5/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S6/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S7/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E01/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E02/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E03/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E04/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E05/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E06/SUMMARY.json`
- `pd_clot/pfpd_candidate/coupling_stability_audit.py`
- `pd_clot/pfpd_candidate/adr_adapter_at1.py`

### fig10 — Implications for future clot simulation

- `verification/pfpd_at1_c7_early_branch_v1/activation/C5_C7_COMPARISON.json`
- `verification/pfpd_at1_c5_stability_v1/SUMMARY.json`
- `verification/pfpd_at1_c7_early_branch_v1/branch/ALL.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/SUMMARY.json`
- `verification/pfpd_at1_activation_root_cause_and_fix_v1/SUMMARY.json`
- `verification/pfpd_at1_c5_activation_stationarity_audit_v1/activation/BRACKET.json`
- `verification/pfpd_at1_c7_early_branch_v1/activation/BRACKET.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/states/S1/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S2/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S3/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S4/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S5/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S6/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/states/S7/INPUT_VERIFICATION.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E01/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E02/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E03/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E04/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E05/SUMMARY.json`
- `verification/pfpd_at1_c5_stability_v1/continuation/E06/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L050/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`
- `../PROJECT_ACTIVE_CONTEXT.md`

## Exclusions and interpretation decisions

- C3 activation is not pooled with C5/C7: earlier geometry/model and loading differences preclude a controlled convergence plot.
- L050/L075 have affine mechanics data only; no activation or fracture values are inferred.
- No C7 curvature or instability interval is plotted: stability is unassessed.
- No unstable mode, rupture morphology, calibrated clot-failure metric or patient result is fabricated.
- “Artifact-consistent” denotes a possible interpretation; the latest root cause remains INCONCLUSIVE.
- No interpolation, smoothing, fitted curves or resampling is used. Lines between saved points are visual guides.
- The optional timing/performance figure is omitted: it would distract from the scientific evidence.
