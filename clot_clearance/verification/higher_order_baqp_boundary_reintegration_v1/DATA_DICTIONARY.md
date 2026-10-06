# Boundary reintegration data

This directory is the record for `HIGHER_ORDER_BAQP_BOUNDARY_REINTEGRATION_V1`. No production or accepted higher-order mechanics code is changed. The configuration was saved before numerical evaluation. `TASK_REQUEST.txt`, the Git snapshot, protected hashes and archived boundary modules retain provenance.

## Coordinate and force conventions

All coordinates/displacements are in metres; nodal volume in m³; nodal force in N; stress and energy density in Pa; total energy and work in J. Reference coordinates are centered Cartesian coordinates. A physical pair appears once, ordered `i < j`; directed gradient weights already include the target volume and have units m⁻¹. F is dimensionless with displacement-output index first and reference-derivative index second.

Primary geometry and PD2 cache are read directly from `../higher_order_baqp_rk_oracle_v1`. The same cache is used throughout all five equilibrium solves. The phase field is absent.

## Main files

| File | Meaning |
| --- | --- |
| PD2_BOUNDARY_CONFIGURATION.json | Frozen physics, geometry, solver settings, gates and execution policy |
| PD2_FORCE_STENCIL_CLOSURE.json | Actual reverse dependency graph; D1 includes the observation point; closure also requires correct nodal volumes |
| PD2_BULK_ORACLE_RAW.json | N27 center and independent shifted N29 evaluations, each with its closure evidence |
| PD2_BULK_ORACLE.csv | 19 matched Z levels per field, q_bulk, local b and q_bulk+b, all density components in N/m³ |
| PD2_TARGET_SUMMARIES.json | Target energy, scales, decomposition errors, affine face data and global integral/work checks |
| PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz | Physical coordinates and volumes, plus arrays under `case__quantity` keys |
| PD2_NORMAL_AFFINE.json etc. | Actual solver termination, residual, displacement/field errors, target and solved energies/work |
| PD2_EQUILIBRIUM_RECOVERY.csv | The five real solve outcomes; classical shear is a control, not an exact-local-recovery requirement |
| PD2_BOUNDARY_SUPPORT.csv | Layer count, depth, mean/max/summed boundary-force magnitudes and cumulative L1 fraction |
| C1_PD2_BOUNDARY_COMPARISON.csv | Matched historical C1 and PD2 load norms and support quantiles; their difference is not a physical corrective force |
| PD2_FORCE_MOMENT_WORK.csv | Affine combined and shear stress-channel comparisons using historical normalization |
| PD2_NONAFFINE_ACCOUNTING.json | Bulk-model, quadrature, internal-force and reference-moment source accounting |
| C1_PD2_FIELD_QUALITY.csv | Saved C1 target diagnostics and new PD2 target/solved field errors by sample depth |
| PD2_CLASSICAL_VS_NONLOCAL_SHEAR.json | Classical/nonlocal displacement difference and force-response resolution relative to measured model load |
| SOURCE_EXTENSION_PROVENANCE.json | Exact before/after hashes for the two allowed boundary-module extensions |

## Nodewise NPZ

Common arrays: `X (6859,3)`, `volume (6859,)`, `depth_over_delta (6859,)`.

Case prefixes: `normal_affine`, `shear_affine`, `nonaffine_axial`, `nonaffine_shear`.

| Suffix after `__` | Meaning |
| --- | --- |
| target | Current nodal coordinates of analytic surrogate, m |
| q_bulk, b_local | Force densities, N/m³ |
| f_FD | Accepted finite-domain internal force at target, N |
| f_bulk_target | q_bulk V, N |
| f_BC_localized | q_bulk V − f_FD, N; the candidate conversion |
| f_BC_naive | −f_FD, N; diagnostic only for nonaffine fields |
| f_body_local | b_local V, N |
| f_model | (q_bulk+b_local)V, N |
| f_external_NL | −q_bulk V + f_BC_localized, N |
| f_external_classical | b_local V + f_BC_localized, N |
| C1_boundary | Matched accepted historical C1 conversion force, N |

## Saved solves and diagnostics

`CASE_SOLVE_INPUT.npz` stores the nonzero perturbed start, fixed-mask and frozen external force. `CASE_RECOVERED_STATE.npz` stores X, solved x, target x, nodal Fbar, internal/external forces, fixed mask and constraint reaction `−(f_int+f_ext)` on fixed nodes only. No full bond tensors or solver scratch columns are retained.

`CASE_FIELD_DIAGNOSTICS.npz` has Fbar and Ftilde prefixes. For each, `F_error` is the Frobenius error to the analytic gradient at the same reference node or bond midpoint; `Psi_plus_discrete`, `Psi_plus_analytic`, and signed `Psi_plus_error` are Pa. Depth is measured at that same sampling location. Ftilde samples count 690020. Nodewise samples count 6859.

Field-profile CSVs use nearest half-h depth bins. Reported maxima use individual sample errors. Weighted RMS uses physical V for nodal F and the accepted energy weights `2*wbar*Vi*Vj` for bonds. Gradient normalization is .003 for affine fields and .004 for the two quadratic fields. Ψ_ref is maximum analytic nodal Ψ⁺ and is common to nodal and bond diagnostics. The old spectral split is diagnostic only.

The unchanged solver's nested summary includes legacy names such as `phase: SCIPY_REFERENCE`, `surface_traction_mean_Pa`, and `maximum_displacement_m`. `phase` here labels the solver stage, not a phase-field evolution; the surface-traction diagnostic is zero because this interface supplies nodal external forces. The legacy maximum-displacement diagnostic refers to its solver geometry/free displacement convention. Use the explicit top-level displacement-shift metric and saved x/X arrays for this round's field comparisons.

## Normalization and interpretation

Affine force/moment/work checks reuse the original total stress-area/stress-volume scales. Mechanical residual is the free-node L2 norm divided by `max(||f_ext,free||, force_floor)`. For affine cases the historical force floor is `||P|| L²`; nonaffine floor is `max(P_ref*h²,b_ref*h³)`. These distinct historical normalization roles are retained explicitly.

Local algebra and closed-center errors use `F_ref = max(P_ref*h²,b_ref*h³)` (affine b_ref=0). A complete closed stencil is established by its graph and coefficients, not merely by observing a small center force. N19 geometric core is not a completely closed bulk region.

The reported external work `sum(f_ext*(x-X))` is work of the frozen nodal load at the specified displacement; it excludes a grip loading path and does not imply equality with stored energy. The classical control retains the nonzero bulk-model difference. No model-error compensation is included.

Support quantiles measure the cumulative sum of nodal force magnitudes versus distance to the nearest face. They are **CONVERSION SUPPORT, NOT PHYSICAL BC THICKNESS**. C1−PD2 loads compare two representations; they cannot uniquely identify a physical correction force.

## Reproduction

Use the existing shared Python environment, `PYTHONPATH=.`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `PYTHONDONTWRITEBYTECODE=1`. The ordered stages of `run_audit.py` are `prepare`, `normal_affine`, `shear_affine`, `nonaffine_axial`, `nonaffine_shear`, `classical_shear`. Recovered-state files guard against accidentally repeating solves in this evidence directory; use a separate output copy for reproduction.

Then run `postprocess.py`, the twelve focused tests, `plot_figures.py`, `write_report.py`, and `finalize_artifacts.py`. Do not run historical full campaigns. Final artifact sealing excludes the self-referential output manifest from its own hash list. Protected input hashes exclude only the two explicitly archived and extended boundary modules; all other protected historical bytes must match.
