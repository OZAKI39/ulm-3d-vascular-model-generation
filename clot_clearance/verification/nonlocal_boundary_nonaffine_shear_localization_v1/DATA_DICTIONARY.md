# Nonaffine shear localization: data dictionary

Task: `NONLOCAL_BOUNDARY_NONAFFINE_SHEAR_LOCALIZATION_V1`. All numerical files use SI units unless a column explicitly says `um`, `depth_h`, or `depth_over_delta`. Reference coordinates are centered on the specimen. Vector component order is x, y, z; tensor indices are spatial row / reference column. No physical phase evolution is represented by these files.

## Frozen inputs and provenance

| File | Meaning |
| --- | --- |
| `TASK_REQUEST.txt` | Exact supplied task text |
| `PREDECLARED_CONFIG.json` | Parameters, gates, normalization, channel and sampling conventions, resources; written before result inspection |
| `STARTING_GIT_STATE.json` | Starting branch, HEAD and working-tree state |
| `ENVIRONMENT.json` | Runtime, CPU and RAM information |
| `PROTECTED_INPUT_HASHES.json` | SHA256 map for 336 historical outputs, configs, manifests, tests and unchanged core sources; paths relative to `clot_clearance/` |
| `APPROVED_NAMESPACE_EXTENSIONS.json` | Three existing diagnostic modules authorized for extension; values are original hashes |
| `source_before/*.py` | Exact pre-extension bytes of those three modules |
| `SOURCE_EXTENSION_PROVENANCE.json` | Before/after hashes and original archive locations; compiled after completion |
| `HISTORICAL_INTEGRITY.json` | Fresh verification of the protected map |
| `PROJECT_ACTIVE_CONTEXT_BEFORE_THIS_ROUND.md` | Original current-context file before updating this round's status |

Old manifests are preserved verbatim. An old output manifest may contain a hash of one of the three extended modules; that entry documents its old version. No assertion is made that such a historical hash equals the new source version. Source archives and extension provenance preserve this distinction.

## Continuum and oracle

| File | Contents and interpretation |
| --- | --- |
| `EXACT_NONAFFINE_SHEAR_REFERENCE.json` | 19 physical Z values and independent closed-form F/C/E/S/P/W/b/J/k; force scales; all-face Gauss integration for combined, normal and tangential channels; 13 continuum virtual modes |
| `NONAFFINE_SHEAR_EXACT_REFERENCE.json` | Byte-identical alias to satisfy the alternate filename in task section 6 |
| `AFFINE_SHEAR_ORACLE.json` | Three beta_s=0 closed-stencil evaluations, performed before nonaffine oracles |
| `AFFINE_SHEAR_LIMIT_RECOVERY.json` | Affine bulk force normalized by F_ref; relative L2 difference from frozen previous B1 load |
| `BULK_STENCIL_CLOSURE.json` | Full dependency and boundary-volume closure for N27 center and shifted N29 reference |
| `BULK_SHEAR_ORACLE_RAW.json` | 57 primary records (19 levels × 3 stress channels), three shifted N29 and three reflected-shear records; force, q=f/h³, observed coordinates and cache counts |
| `BULK_BAQP_SHEAR_REFERENCE.csv` | 57 rows; channel, level, Z_m, q_x/y/z_N_m3, b_x/y/z_N_m3, e_x/y/z_N_m3 with e=q+b; q/b/e_norm_N_m3 are Euclidean vector norms derived from those saved components |
| `BULK_SHEAR_ORACLE_VALIDATION.json` | Independent invariance checks, channel sum, component/norm differences and descriptive model classification |
| `SELECTIVE_VS_FULL_SHEAR.json` | Three unchanged-selective versus full BAQP nodal-force comparisons on the same N19 geometry |

`P_ref_Pa = max ||P||_F`, `b_ref_N_m3 = max ||b||`, and `F_ref_N = max(P_ref h²,b_ref h³)` use the physical N19 Z range. The N27 depth experiment keeps this force scale and the same displacement polynomial length L=18h. N27 is an extended stencil specimen, not a different polynomial benchmark.

The normal channel uses the diagonal of the actual discrete Piola stress; the tangential channel uses its off-diagonal part. On the axis-aligned continuum faces this equals `(t·N)N` and `t−(t·N)N`. The force assembly is linear in supplied P at frozen geometry. This does not assign distributed nonlocal nodal forces to unique faces and does not assert nonlinear material superposition.

## Nodewise decomposition and global balances

`NODEWISE_SHEAR_DECOMPOSITION.npz` has 6859 nodes in the unchanged geometry/family order. Coordinates are in metres, forces in N, force densities in N/m³, volumes in m³, F tensors dimensionless.

| Key | Shape / definition |
| --- | --- |
| `X`, `x_target` | (6859,3), reference coordinates and exact displaced target |
| `volume` | (6859,), physical trapezoid V_i |
| `Fbar_target`, `F_analytic` | (6859,3,3), discrete mean gradient and analytic gradient at the same nodes |
| `q_bulk`, `b_local`, `e_model` | (6859,3), level-mapped bulk oracle, local body force, and q+b |
| `f_FD` | (6859,3), unchanged BAQP internal force at x_target |
| `f_bulk_target` | q_bulk V_i |
| `f_BC_localized` | q_bulk V_i−f_FD |
| `f_BC_naive` | −f_FD; diagnostic control only |
| `f_body_local` | b_local V_i |
| `f_model` | (q_bulk+b_local)V_i |
| `f_external_NL` | −q_bulk V_i+f_BC_localized |
| `f_external_classical` | b_local V_i+f_BC_localized |
| `f_BC_normal`, `f_BC_tangential` | Force conversion from each frozen stress channel |
| `distance_faces` | (6859,6), reference distances to x−, x+, y−, y+, z−, z+ |
| `d_min_over_h`, `d_min_over_delta` | (6859,), nearest-face distance normalized by h or delta |

| File | Meaning |
| --- | --- |
| `SHEAR_NODEWISE_IDENTITIES.json` | Maximum nodewise closure, naive identity, NL target identity and channel/work errors; exact-target BAQP energy |
| `SHEAR_FORCE_BALANCE.json` | Classical surface, exact/nodal body, internal, conversion and model resultants; nodal body quadrature correction |
| `SHEAR_MOMENT_BALANCE.json` | Same terms in reference-coordinate moment about X=0, including the continuum reference-moment source; current-position internal moment reported separately |
| `SHEAR_CHANNEL_ACCOUNTING.json` | Full global force/moment accounting for all three stress channels and nodewise channel sum checks |
| `SHEAR_VIRTUAL_WORK.csv` | 39 rows, 13 variations × 3 channels; internal/localized/body/corrected/model work in J and stable normalized closure |

Force signs: equilibrium is `f_internal + f_external = 0`. `corrected_residual` is `f_FD+f_BC_localized+b_local V`; it equals model load, not zero. Global force bookkeeping keeps the nodal quadrature correction and the computed internal resultant. Global reference-moment bookkeeping additionally subtracts the nonzero continuum source `∫Σ e_a×P_:a dV`; this source must not be omitted in shear.

Virtual modes are `translation_x/y/z` (unit 1 m virtual translations), all `A_ab` (δu_a=X_b), and `quadratic_x` (δu_x=Z²/L). Discrete work scale is `F_ref*max(L,max ||δu||)`; continuum work scale is `P_ref L²*max(L,max ||δu||)`.

## Conversion localization

| File | Meaning |
| --- | --- |
| `SHEAR_BOUNDARY_LOCALIZATION.csv` | Both localized/naive curves: depth_h, depth_m, depth_over_delta, node count, L1 magnitude sum, mean, max, cumulative fraction |
| `SHEAR_SUPPORT_QUANTILES.json` | d50/d90/d95/d99 in m, µm, h and delta; first discrete layer meeting the cumulative magnitude fraction |
| `N27_SHEAR_DEPTH_CHECK.json` | Z=0 x-depth 4h,8h,12h,13h: full force/bulk/correction vectors, magnitudes, naive values and actual dependency closure; independent shifted N29 comparator |

These are conversion-support measures for this benchmark, not physical traction-layer thicknesses. Cumulative magnitudes sum already volume-weighted nodal forces; no additional V_i weighting is applied. N19 center depth is 9h and has no fully closed force stencil. Numerically tiny residual alone does not prove dependency closure; N27 12h is not fully closed, whereas 13h is.

## Mechanical recovery

`SHEAR_NL_SOLVE_INPUT.npz` and `SHEAR_classical_SOLVE_INPUT.npz` store the actual nonzero initial configuration, frozen external force, and Boolean fixed mask. Both z faces have prescribed target positions; all other nodes are free. The same smooth perturbation includes an x component.

| File | Meaning |
| --- | --- |
| `SHEAR_NL_MANUFACTURED_RECOVERY.json` | Unchanged solver diagnostics, force residual, discrete/analytic gradient errors, displacement/energy errors, initial perturbation and fixed-value check |
| `SHEAR_CLASSICAL_MANUFACTURED_RESULT.json` | Same quantities using local body force with unchanged frozen boundary conversion |
| `SHEAR_NL_RECOVERED_STATE.npz`, `SHEAR_classical_RECOVERED_STATE.npz` | X/x/target/initial (m), volume (m³), external/internal/reaction (N), fixed mask, Fbar, nodewise normalized F errors, u_error (m) |
| `SHEAR_NL_RECOVERY_PROFILE.csv`, `SHEAR_classical_RECOVERY_PROFILE.csv` | Nearest-face depth in h/delta, node count, maximum/unweighted nodal RMS discrete-target gradient error, analytic maximum error, maximum displacement error |
| `SHEAR_CLASSICAL_VS_NL_ACCOUNTING.json` | Saved-state comparison of same initial state/constraints, load difference and free-node force response against independently measured model load |

Recovery gradient errors divide Frobenius norms by 0.004. The recovery-profile RMS is an ordinary nodal RMS; the analytic-bias diagnostic below instead uses the explicitly predeclared physical volume/bond weights. Force residual uses free-node L2 divided by max(frozen external free-node L2,F_ref). Fixed-node reactions are stored but do not enter this convergence criterion. Legacy solver `surface_traction_*` zeros are placeholders, not calculated classical traction values.

## Kinematic and positive-energy bias

At exact target positions, Fbar compares to analytic F at reference nodes; Ftilde compares to analytic F at reference bond midpoints. Depth is measured at the corresponding location. The unchanged candidate spectral split supplies Psi+ for both exact and discrete gradients.

| File | Meaning |
| --- | --- |
| `SHEAR_FIELD_BIAS_SUMMARY.json` | Max/weighted RMS F errors (absolute and /0.004), Psi+ errors (Pa and /common reference), signed weighted mean Psi+ error, exact repeated-diagnostic determinism, discrete tensor hashes and descriptive classification |
| `SHEAR_FIELD_BIAS.csv` | Same metrics at each depth, in half-h bins (integer h for nodes, half-h for bond midpoints); sample type, depth_h, depth_over_delta, sample count |
| `SHEAR_FIELD_BIAS_SAMPLES.npz` | `Fbar_*` and `Ftilde_*`: F_error (Frobenius), Psi_plus_discrete/analytic/error (Pa), depth_over_delta, weight; no full large Ftilde tensor archive |

Node weights are V_i; bond weights are `2*wbar*V_i*V_j`, exactly as in the BAQP energy quadrature. RMS normalizes each weight sum before averaging squares. Psi+ relative errors use one stable scale, maximum analytic Psi+ over physical N19 nodes; no component-wise division near zero. Absolute Psi+ error means the maximum of absolute signed errors, not the maximum signed error.

`KINEMATIC_BOUNDARY_BIAS` uses the larger normalized maximum over Fbar and Ftilde: LOW <1e-3, MODERATE ∈[1e-3,1e-2], MATERIAL >1e-2. These are report-only descriptors. `PSI_PLUS_BOUNDARY_BIAS=REPORTED_NOT_ACCEPTANCE_RATED`; there is no fracture-safe threshold.

## Logs, outputs and reproducibility

`REFERENCE_AFFINE_RUN.log`, `ORACLE_RUN.log`, `DECOMPOSITION_RUN.log`, `N27_RUN.log`, `NL_RUN.log`, `CLASSICAL_RUN.log` are the scientific stage logs. Each `RESOURCES_<stage>.json` records wall time, peak process-tree RSS, and minimum system-available RAM sampled every 0.1 s. RSS can count shared pages multiple times. `RESOURCE_SUMMARY.json` sums scientific stage time only; authoring, plots and tests are excluded. Oracle work uses eight CPU workers, equilibrium one serial worker, BLAS one thread. GPU is unused by design.

`FOCUSED_TESTS_PRE_SOLVE.log` records 10 passed with only the solve-dependent test deselected. `FOCUSED_TESTS_FINAL.log` records all 11 requested tests after results exist. No unrelated or full tests are run. `FINAL_STATUS.json` contains all task section 48 labels; `TERMINAL_SUMMARY.json` supplies section 52 values. `OUTPUT_MANIFEST.json` hashes final artifacts and changed/new sources after files are closed, excluding itself.

Entrypoints (relative to project root):

```text
scripts/nonaffine_shear_boundary_audit.py --stage reference_affine
scripts/nonaffine_shear_boundary_audit.py --stage oracle
scripts/nonaffine_shear_boundary_audit.py --stage decomposition
scripts/nonaffine_shear_boundary_audit.py --stage n27
scripts/nonaffine_shear_boundary_audit.py --stage NL
scripts/nonaffine_shear_boundary_audit.py --stage classical
scripts/plot_nonaffine_shear_boundary.py --stage all
scripts/report_nonaffine_shear_boundary.py
tests/test_nonaffine_shear_localization_v1.py
```

Scientific stages are ordered gates, not an instruction to rerun accepted artifacts automatically. Keep their logs and results intact. Use the existing Python environment with `PYTHONPATH=.`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`/`python -B`. All figures read saved data and export English text in 600 dpi PNG plus vector PDF.
