# NONLOCAL_BOUNDARY_FORMULATION_V1 — review

Date: 2026-10-05. Model: `PFPD_AT1_SPECTRAL_CANDIDATE_V1`. **Intact mechanics only.**

Decision: **VERIFIED_FOR_INTACT_AFFINE_MECHANICS**. Scope is the two analytic affine block benchmarks and the exact discretization documented below. This is not validation of general nonlinear clot, fracture, adhesion or CFD boundary conditions. The construction identity is not the acceptance test. B0/B2 virtual-work failures are retained below.

## 1. Why bulk BAQP development is no longer the current blocker

The preceding affine/coefficient audits support the existing variational BAQP operator and explain the small raw Z mismatch through reconstructed-stress floating-point propagation. Their original raw FAIL records remain unchanged. This round changes boundary data, not bulk mechanics. 56 protected input/source files retain their SHA256 values. See [integrity record](HISTORICAL_INTEGRITY.json); historical C5/C7 scientific states are neither changed nor rerun.

## 2. What the four papers actually contribute

The complete source map, equation locations, access limitations and distinctions are in [LITERATURE_TO_BAQP_MAPPING.md](LITERATURE_TO_BAQP_MAPPING.md). D'Elia supplies the conversion architecture; Yu supplies a traction-consistency reference; Scabbia informs an independent ghost-control design; Shojaei/perifit concerns a separate operator correction. Public full texts were read before coding because the uploaded PDF bytes were unavailable locally. The local PDF provenance does not certify identity with those uploads.

## 3. What was deliberately not copied

No scalar diffusion theorem is asserted for BAQP. No Yu LPS traction formula, OSB dilatation/weighted-volume ghost formula, or perifit weight is inserted. No kernel, K/F/Z coefficient rule, material model, phase law or historical assertion is changed. No arbitrary 3-delta physical collar is introduced. [Compatibility notes](PERFIT_BAQP_COMPATIBILITY_NOTES.md) identify future energy/variation obligations.

## 4. Boundary API design

The independent namespace `pd_clot/pfpd_candidate/nonlocal_boundary/` provides Cartesian geometry, face partitions, normals/tangents, distance/projection, physical/interior/exterior interaction masks, analytic local solutions, frozen volume data and diagnostics. [BOUNDARY_SPECIFICATION.md](BOUNDARY_SPECIFICATION.md) gives equations and node roles. Main geometry: N19, 6,859 physical particles, 690,020 bonds, L=1.125 mm, h=62.5 µm, delta=251.25 µm, delta/h=4.02. SVK is explicit: mu=1000 Pa, lambda=0, K=666.6666667 Pa.

Pre-run peak estimates are 1.41 GiB for N19, 2.46 GiB for fixed-delta N21, and 1.85 GiB for fixed-ratio N21. The estimates and all gates preceded results in [PREDECLARED_CONFIG.json](PREDECLARED_CONFIG.json). The geometric core has 729 nodes with d>=delta; it is not a completely closed bulk stencil.

## 5. D'Elia-style BAQP volume-force conversion

B1 applies the existing internal-force evaluator to u_L and freezes its negative on a topology-derived support. Nodal N and volumetric N/m³ units are explicit. Load updates during equilibrium are forbidden. The complete support covers this N19 domain; omitted-support force is 0.000000e+00 N in normal loading and 0.000000e+00 N in shear. The exact cancellation identity is retained as a check, not scientific proof. Analytic continuum targets, all affine work modes, relabeling, perturbed-start recovery and limited refinement provide additional evidence.

This finite-domain in-domain force representation differs from the original exterior nonlocal-flux problem. It is a restricted surrogate-dependent construction, not a theorem that any prescribed traction admits a unique such representation.

## 6. Normal traction benchmark

F=diag(1,1,1.003); Ezz=0.0030045, Szz=6.009 Pa, Pzz=6.027027 Pa, W=0.00902702025 Pa. Every face uses PN. Each z-face resultant has magnitude 7.627956e-06 N; closed-surface target resultant is zero. E/S/P/W and each exact face load are saved in [NORMAL_BENCHMARK.json](NORMAL_BENCHMARK.json) and [CLASSICAL_FACE_TARGETS.csv](CLASSICAL_FACE_TARGETS.csv).

| Candidate | Solver force gate | Force error | Moment error | Max work error | Core F RMS | Boundary F RMS | Max F error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B0 | PASS | 8.682237e-16 | 9.545280e-16 | 1.200251e-01 | 1.193001e-01 | 3.475958e-01 | 9.275087e-01 |
| B1 | PASS | 1.185398e-09 | 7.589290e-10 | 1.171503e-09 | 1.465516e-09 | 2.448481e-09 | 1.006357e-08 |
| B2 | PASS | 2.575594e-15 | 3.429991e-16 | 1.935572e-05 | 5.529249e-03 | 1.969700e-03 | 9.487703e-03 |

F errors are relative to the imposed gradient norm 0.003. Force and moment use stable stress-area and stress-volume scales; nonzero work uses its own target. All three load columns include applied forces plus constraint reactions. Frozen B1 load checks are reported separately in section 9.

## 7. Volume-Dirichlet reference

B2 is **LOCALITY-STRONG REFERENCE**: all physical nodes d<delta on all six faces are prescribed to the local surrogate; no external nodes, ghost energy or extra material are introduced. There are 6130 constrained physical particles and 729 free particles. Constraint matching is exact to 0.000000e+00 m. This is an in-domain support adaptation, not a demonstrated equivalent exterior-domain boundary-value problem. Small constrained-region error alone would not justify selecting it. One-horizon constraints need not remove every multi-hop influence in the free core.

The normal and shear B2 work errors are 1.935572e-05 and 2.387937e-05; both exceed 1e-6 despite mechanical convergence. The prescribed support was not widened after seeing these outcomes.

## 8. Scabbia control if implemented

`SCABBIA_GHOST_CONTROL=DEFERRED` (not implemented in V1). Taylor continuation defines affine ghost positions, but does not define ghost-side BAQP K/F/Z, the physical-only quadrature energy, or the pullback of reactions when fictitious positions depend on physical DOFs. Treating those nodes as ordinary extra material would change the problem. No ghost run or unverifiable independent-control claim is made.

## 9. Normal resultant / moment / virtual work

For the frozen B1 load: force error 8.458372e-18, moment error 7.940165e-18, maximum work error 1.262731e-13. The predeclared gates are 1e-6. All nine affine A basis modes and three translations are saved in [VIRTUAL_WORK_RESULTS.csv](VIRTUAL_WORK_RESULTS.csv). Zero closed-surface force alone is not sufficient: a deliberately zero load fails nonzero-work unit tests.

Normal load relabeling changes the force by 1.543771e-15 relative L2, below 1e-10. The upper-half volume-load versus top-face target discrepancy is 3.834141e-07 on the stress-area scale. That latter check is diagnostic; no unique face assignment is assumed.

## 10. Shear traction benchmark

Shear was run after normal force/work gates passed. F=I+0.003 ex tensor ez. Nonzero P components [Pa] are Pxx=Pzz=0.009, Pxz=3.000027, Pzx=3. The exact energy is 0.00450002025 Pa. Lateral x-face traction is included. All tensors, targets and control states are in [SHEAR_BENCHMARK.json](SHEAR_BENCHMARK.json).

| Candidate | Solver force gate | Force error | Moment error | Max work error | Core F RMS | Boundary F RMS | Max F error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B0 | PASS | 7.810652e-15 | 3.508877e-07 | 1.728892e-01 | 1.262455e-01 | 4.749636e-01 | 1.460653e+00 |
| B1 | PASS | 1.374111e-09 | 6.560412e-10 | 3.610934e-09 | 2.720779e-09 | 6.206027e-09 | 2.936874e-08 |
| B2 | PASS | 8.886873e-15 | 6.744566e-11 | 2.387937e-05 | 4.762073e-03 | 1.696994e-03 | 7.259162e-03 |

## 11. Shear resultant / moment / virtual work

Frozen B1 errors: force 1.306019e-15, reference moment 7.438334e-17, work 1.152770e-11. The reference-coordinate target My is 3.844336e-14 N m, not zero. Relabeling force error is 9.078014e-14.

Moment gates use the stable stress-volume scale specified before the run. The balance CSV also gives error relative to the small nonzero moment itself: shear B0 differs by 5.51% of that target, despite passing the stable moment gate. The separate virtual-work and known-field gates remain necessary.

| Shear channel | Force error | Moment error | Max work error |
| --- | --- | --- | --- |
| normal | 8.724737e-16 | 5.240399e-17 | 3.782262e-16 |
| tangential | 1.940525e-15 | 5.470503e-16 | 1.936518e-16 |

These channels use constant Piola-stress diagonal/off-diagonal components in the unchanged force assembly. They verify Cartesian normal/tangential accounting at frozen geometry; they are not a nonlinear constitutive superposition theorem. Face-wise normal/tangent vectors are retained separately in the target CSV.

## 12. Spatial support of the equivalent nonlocal load

| Load | d50 [µm / h] | d90 [µm / h] | d95 [µm / h] | d99 [µm / h] |
| --- | --- | --- | --- | --- |
| normal | 62.50 / 1 | 125.00 / 2 | 187.50 / 3 | 250.00 / 4 |
| shear | 62.50 / 1 | 125.00 / 2 | 187.50 / 3 | 250.00 / 4 |

The distribution is weighted by summed nodal force magnitudes. Direct-interaction nodes number 6130; the graph support contains all 6,859 nodes. d99 near one horizon does not imply a compact one-horizon physical BC. A small interior tail is retained. [BOUNDARY_FORCE_SUPPORT.csv](BOUNDARY_FORCE_SUPPORT.csv) records both d/h and d/delta.

The final depth bin contains only the center particle. Its near-zero force is consistent with inversion symmetry and must not be used to assert three-hop stencil closure. The favorable N19 half-space diagnostics do not retroactively validate the much smaller historical C5/C7 specimens.

Small predeclared refinement: fixed L, N19 versus N21, with fixed delta or fixed delta/h. Only affine conversion/work/support is reevaluated at N21; these are not additional equilibrium convergence solves.

Both main B1 cases passed before these four N21 evaluations.

| Resolution / load | h [µm] | delta [µm] | Force error | Moment error | Work error | d99 / delta |
| --- | --- | --- | --- | --- | --- | --- |
| base / normal | 6.250000e+01 | 2.512500e+02 | 8.458372e-18 | 7.940165e-18 | 1.262731e-13 | 9.950249e-01 |
| base / shear | 6.250000e+01 | 2.512500e+02 | 1.306019e-15 | 7.438334e-17 | 1.152770e-11 | 9.950249e-01 |
| fixed_delta_refined / normal | 5.625000e+01 | 2.512500e+02 | 2.713728e-17 | 6.307417e-18 | 1.131639e-13 | 8.955224e-01 |
| fixed_delta_refined / shear | 5.625000e+01 | 2.512500e+02 | 2.124066e-15 | 7.285824e-16 | 2.465202e-11 | 8.955224e-01 |
| fixed_ratio_refined / normal | 5.625000e+01 | 2.261250e+02 | 1.547990e-17 | 2.037997e-17 | 1.873855e-13 | 9.950249e-01 |
| fixed_ratio_refined / shear | 5.625000e+01 | 2.261250e+02 | 2.489225e-15 | 9.717531e-16 | 1.148017e-11 | 9.950249e-01 |

Two resolutions cannot establish asymptotic compatibility, mesh objectivity, or fracture convergence. Grid-quantized d99 can change even while the force distribution is close.

## 13. Does the new method reduce the old artificial boundary field?

| Load / method | Psi+ max [Pa] | Psi+ P95 [Pa] | Energy [J] | Top reaction x [N] | Top reaction z [N] | First interior jump norm / affine |
| --- | --- | --- | --- | --- | --- | --- |
| normal / B0 | 8.948699e-02 | 1.954589e-02 | 1.130715e-11 | 1.142724e-21 | 6.712410e-06 | 2.928427e+00 |
| normal / B1 | 9.027020e-03 | 9.027020e-03 | 1.285293e-11 | -4.301820e-15 | 4.409373e-15 | 1.000000e+00 |
| normal / B2 | 9.521126e-03 | 9.115930e-03 | 1.285268e-11 | -3.415576e-23 | 4.219616e-06 | 1.000000e+00 |
| shear / B0 | 7.153618e-02 | 7.566382e-03 | 7.173147e-12 | 3.625969e-06 | 1.335994e-08 | 2.707309e+00 |
| shear / B1 | 2.256760e-03 | 2.256760e-03 | 6.407255e-12 | -4.655011e-15 | -2.773005e-15 | 1.000000e+00 |
| shear / B2 | 2.428590e-03 | 2.284353e-03 | 6.407102e-12 | 2.100429e-06 | 6.301092e-09 | 1.000000e+00 |

Both B0 and B1 meet the mechanical force gate. The normal core gradient RMS changes from 1.193001e-01 to 1.465516e-09; boundary RMS changes from 3.475958e-01 to 2.448481e-09. The positive-energy maximum relative to the affine target changes from 9.913237e+00 to 1.000000e+00.

The first-jump ratio uses the centerline first interior layer relative to its affine displacement increment. Top constraint reaction is separate from already-applied B1 volume force; a nearly zero B1 grip reaction does not mean absent physical loading. Full vectors, reactions and free-node residual profiles are retained in candidate JSON/CSV/NPZ files. These Psi+ changes are intact mechanical diagnostics, not fracture corrections.

Each `*_DISTANCE_PROFILE.csv` records gradient RMS/maximum and nodal energy against the nearest of all six physical faces. `*_FIELD_PROFILE.csv` separately records profiles from the z-grip planes; its original plane-RMS residual column is retained and an explicit free-only RMS/count is supplied. Layers without free nodes have an empty free-only RMS, not an artificially perfect equilibrium metric.

Known-solution recovery is accepted only when the existing solver reaches its 1e-8 force gate and the maximum reconstructed F error is <=1e-5 of the imposed gradient. The final statuses below enforce that additional criterion. All six equilibrium solves meet the force gate, supporting interpretation of the large B0/B2 work and field discrepancies as differences in their discretized boundary treatments at this grid.

| Case | BAQP evaluations | Tangent actions | Wall [s] | Termination |
| --- | --- | --- | --- | --- |
| normal / B0 | 4 | 742 | 7.251947e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |
| normal / B1 | 6 | 989 | 9.036305e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |
| normal / B2 | 3 | 157 | 1.571495e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |
| shear / B0 | 5 | 989 | 7.742080e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |
| shear / B1 | 6 | 989 | 7.657339e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |
| shear / B2 | 3 | 183 | 1.428871e+02 | SCIPY_REFERENCE_FORCE_BALANCE_PASSED |

The adapter runs the existing SciPy solver code, supplying an analytic BAQP tangent action; no new nonlinear iteration is written. Its derivative is checked against force finite differences and reciprocal bilinear products. Matrix-free condition number is explicitly unavailable. A short initial attempt was interrupted to fix the reporting adapter's float conversion; it accepted no state and its raw log/resources remain saved. After the normal cases, an equivalent CSR bond-to-node aggregation was enabled for subsequent cases. [TANGENT_SPARSE_ASSEMBLY_CHECK.json](TANGENT_SPARSE_ASSEMBLY_CHECK.json) compares its action to the archived initial adapter, with relative difference 3.29e-16 on N9. No tolerance, force law or nonlinear algorithm changed; the initial adapter source is retained under `references/`.

## 14. What remains surface/operator correction rather than BC conversion?

`SURFACE_OPERATOR_CORRECTION_STATUS=SEPARATE_ISSUE`. Affine conversion under the unchanged variational operator is the demonstrated object. It does not remove all effects of truncated horizons, prove a boundary law for arbitrary nonaffine fields, or establish any need for a weighted BAQP operator. A manufactured affine load can preserve global work even when local force support is distributed. Future evidence must separate external data conversion from operator changes.

## 15. perifit compatibility assessment

[PERFIT_BAQP_COMPATIBILITY_NOTES.md](PERFIT_BAQP_COMPATIBILITY_NOTES.md) describes candidate BAQP targets: shape invertibility/completeness, homogeneous energy, a common force first variation, and Z-mediated contributions. Symmetric scalar weights alone do not prove all tensor and energy conditions. No fit or kernel modification is performed.

## 16. Implication for future clot-wall mechanics

The API accepts geometric displacement regions and pressure-plus-tangential traction data. Clot-wall adhesion would require an independent interface model. The present Cartesian affine surrogate cannot represent a real clot's curved, heterogeneous or nonaffine solution. Mixed loading, physical clot-wall prototypes, real-clot material calibration, activation, stability, fracture, erosion, microbubble loading and CFD remain NOT_RUN.

## 17. One next development task

After human review only: **NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1** — test a known nonaffine elastic field with compatible body force and surface data, and determine which part of the finite-domain residual is genuinely boundary data. Require spatially resolved traction/work and interior-load accounting before claiming a transferable nonlocal Neumann law. Do not start this task automatically.

## Final status and audit evidence

| Status | Value |
| --- | --- |
| BOUNDARY_API | IMPLEMENTED |
| DELIA_BAQP_VOLUME_FORCE | IMPLEMENTED |
| NORMAL_TRACTION_RESULTANT | PASS |
| NORMAL_TRACTION_MOMENT | PASS |
| NORMAL_VIRTUAL_WORK | PASS |
| NORMAL_BOUNDARY_CANDIDATE | VERIFIED |
| SHEAR_TRACTION_RESULTANT | PASS |
| SHEAR_TRACTION_MOMENT | PASS |
| SHEAR_VIRTUAL_WORK | PASS |
| SHEAR_BOUNDARY_CANDIDATE | VERIFIED |
| VOLUME_DIRICHLET_REFERENCE | IMPLEMENTED |
| SCABBIA_GHOST_CONTROL | DEFERRED |
| PERIFIT | COMPATIBILITY_REVIEW_ONLY |
| BAQP_SPECIFIC_SURFACE_CORRECTION | NOT_IMPLEMENTED |
| BOUNDARY_DATA_CONVERSION_STATUS | VERIFIED_FOR_INTACT_AFFINE_MECHANICS |
| SURFACE_OPERATOR_CORRECTION_STATUS | SEPARATE_ISSUE |
| C5_C7_ACTIVATION | NOT_RUN |
| FRACTURE | NOT_RUN |
| CLOT_WALL_PROTOTYPE | NOT_RUN |

Focused tests only:

```text
...........                                                              [100%]
11 passed in 0.59s
```

An initial test run had two new-test API-name typos; these were fixed without changing scientific assertions, and its log is preserved. No full pytest suite was run. Historical missing artifacts and old failures were not modified.

Resources: one CPU worker, one BLAS thread, GPU_NOT_USED_BY_DESIGN. Summed measured stage wall time 3560.08 s, including the interrupted attempt; peak recorded RSS 1.436 GiB; minimum available system RAM 43.1%. Authoring/literature time is outside this compute timing. [RESOURCE_SUMMARY.json](RESOURCE_SUMMARY.json) preserves per-stage measurements. No memory-heavy solves overlap.

Figures: [results/nonlocal_boundary_formulation_v1/figures](../../results/nonlocal_boundary_formulation_v1/figures), English labels, 600 dpi PNG and vector PDF, generated from saved numerical data by `scripts/plot_nonlocal_boundary.py`.

The observed Python/package versions are in [ENVIRONMENT.json](ENVIRONMENT.json); no dependencies were installed. [SOLVER_PROVENANCE.json](SOLVER_PROVENANCE.json) identifies the adapter version used for each loading family.

**Stopped for human review. No subsequent scientific stage is authorized by this report.**
