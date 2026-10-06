# Development Workflow and Current Progress

## Objective

The core mechanics and solver framework are established. The current bottleneck is finite-domain BAQP boundary consistency, which must be resolved before physical clot-wall failure can be interpreted.

One editable 16:9 slide, based only on archived repository evidence. No new mechanical solve, activation search, stability calculation or parameter sweep was run.

## Source files

- `../PROJECT_ACTIVE_CONTEXT.md`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/affine_internal_force_layers.png`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/first_layer_force_decomposition.png`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/figures/large_domain_affine_patch_test.png`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C5_layers.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/C7_layers.csv`
- `results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_particles.csv`
- `verification/pfpd_at1_c5_stability_v1/REVIEW_C5_STABILITY_V1.md`
- `verification/pfpd_at1_c7_early_branch_v1/REVIEW_C7_EARLY_BRANCH_V1.md`
- `verification/pfpd_at1_petsc_vi_solver_v1/REVIEW_PETSC_VI_SOLVER_V1.md`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/BOUNDARY_THEORY_NOTES.md`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C5/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/C7/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/REVIEW_BAQP_FINITE_DOMAIN_BOUNDARY_CONSISTENCY_V1.md`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/SUMMARY.json`
- `verification/pfpd_baqp_finite_domain_boundary_consistency_v1/author_L075/COMPARISON.json`

## Quantitative statements shown

- Author implementation all-force relative difference: 5.861129936374264e-15; shown as 5.9 × 10⁻¹⁵.
- L075 denotes L=0.75 mm: 27 self-family-complete points; 0 fully three-hop-closed points (also 0 two-hop-closed points).
- 3δ is the conservative geometric coefficient-dependency bound derived in the latest theory notes. It is a proposed isolation-test design, not an observed optimal physical boundary thickness.
- The three plots reproduce actual CSV/JSON coordinates and contributions. Their tick values are axis scales, not additional claimed experimental findings.

## Stage interpretation

### 1. Theory & Discretization Audit — COMPLETED

Nodal quadrature; Kernel / weights; Horizon handling; Nonlocal kinematics defined.

### 2. BAQP Mechanics: Affine, Energy & First Variation — COMPLETED

Affine kinematics; Energy + first variation; Spectral mechanics; Finite-domain equilibrium unresolved.

### 3. Independent Cross-Check: RKCMaterial / Peridynamics.jl — COMPLETED

Local / author match; No clear porting error; Fixed author implementation.

### 4. Finite-Domain Nonlocal Boundary Consistency — CURRENT

Affine ≠ equilibrium; Z / neighbor-Z channels; Cause: INCONCLUSIVE; Boundary influence is a leading hypothesis.

### 5. Fully Closed 3δ Force-Stencil Patch Test — NEXT

Isolate central region; Close force dependencies; Separate boundary vs. bulk; Conservative support bound; test not run.

### 6. Physical Clot Simulation — FUTURE

Adhesion / detachment; Damage / rupture; Flow / microbubble load; Not quantitatively validated.

## Interpretation limits

- Root cause remains formally INCONCLUSIVE. Boundary influence is a leading unresolved mechanism.
- No claim of a confirmed boundary artifact, wrong BAQP operator or incorrect surface BC is made.
- Geometric full family does not imply a closed force-dependency stencil.
- Agreement with the author code and energy gradient is an implementation check, not proof of physical traction correctness.
- C5 audited stable positive-phase states and C7 stationary states are different evidence; C7 stability is unassessed.
- Activation is not a crack prediction; stationarity is not stability.
- The clot drawing is a geometric schematic only: wall, attached mass, near-wall region and flow direction. No simulated contours, cracks or validated clot-wall failure are depicted.
- 3δ test is NEXT, not already run. Physical clot adhesion, detachment, rupture and loading calibration are FUTURE.

## Insets and editability

The same saved result panels are re-typeset from their original data, rather than cropping small raster labels that would be unreadable in projection. No interpolation, smoothing or fabricated trend is introduced. Workflow boxes, arrows, status icons, all narrative text, the progress indicator and the clot schematic are native editable PowerPoint objects; only the three real-data insets are images. PDF/SVG/PNG inset originals are included in assets/.

## Steering and reproduction

Edit `slide_config.py` (central `SLIDE_CONFIG`), then run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python generate_slide.py
```

Controls include title, subtitle, current stage, stage labels/status/body, colors, density, inset selection (up to three), numeric callouts, schematic visibility and bottom messages. The default uses DejaVu Sans, the installed sans-serif font, to keep the editable slide and rendered preview consistent. Low density uses shorter body statements; high density uses more technical body statements. Both wordings are defined centrally in `stage_bodies_by_density`.

The PNG is rendered from the actual PowerPoint via LibreOffice → PDF → PyMuPDF, not a separately reconstructed drawing. The temporary local LibreOffice runtime does not modify the system installation or scientific environment. If LibreOffice is available on PATH it is preferred; otherwise the current task runtime is used.

The earlier 10-figure inventory and plan under `analysis/paper_figures/` were superseded by the user's one-slide steering before their figure generation. They remain as preparation, not as an additional claimed deliverable.
