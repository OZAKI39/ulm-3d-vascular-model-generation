# PFPD_BAQP_AFFINE_COEFFICIENT_REFERENCE_AUDIT_V1

**Decision: RECONSTRUCTED_STRESS_FLOATING_POINT_PROPAGATION.** The raw excess is explained by reconstructed-stress floating-point propagation within the explicitly measured Z assembly floor. Exact affine stress recovers bulk coefficients. The dominant Z component is well conditioned; cancellation-sensitive summation is not the main source.

**Previous raw reconstructed-Z gate remains FAIL under its original 1e-12 criterion.** This document is a new interpretation layer. No old result, review, tolerance, production formula or accepted state was changed.

## 1. Why the previous audit stopped at NOT_ESTABLISHED

The [previous audit](../pfpd_baqp_fully_closed_force_stencil_audit_v1/REVIEW_FULLY_CLOSED_FORCE_STENCIL_V1.md) measured N27 center force 1.945942391685588e-22 N, local normalized 8.265455792242755e-15, passing 1e-10. Its conjunction also required raw coefficients below 1e-12. Raw Z = 1.1362855097181066e-12 exceeded that threshold. Its `N27_COEFFICIENT_GATE=FAIL`, `BULK_AFFINE_FORCE_PATCH_TEST=NOT_ESTABLISHED` and primary `INCONCLUSIVE` labels remain historical facts. They are not retroactively replaced by this audit.

## 2. What is unchanged from the previous audit

h=6.25e-5 m, delta=0.00025125 m, epsilon=0.003, mu=1000 Pa, lambda=0, bulk=666.6666666666666 Pa; quadratic kernel; phase disabled, s=0. The saved N27 center D1/D2/D3 closure is reused. Saved integer offsets, oracle geometry, weights, F/P and Z are checked bitwise against the unchanged local helper before reuse. Coordinates were not previously saved; the small N27 coordinate array is reconstructed identically. Only 257 Z nodes, 1797 Fbar nodes and 50,356 unique incident physical bonds are evaluated (65,792 directed Z contributions). No full-domain bond state or new force/equilibrium/activation campaign is built.

Historical input/results protection covers 51 directly relevant files, including all previous audit files/figures and the imported production constitutive, BAQP and kernel sources. [HISTORICAL_INTEGRITY.json](HISTORICAL_INTEGRITY.json) records exact SHA256 equality. No historical cleanup manifests or full pytest were run. New computation provenance is in [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json).

## 3. Exact affine F / P reference

Starting from the current SVK energy, with lambda=0,

$$E=\tfrac12(F^T F-I),\quad S=2\mu E,\quad P=FS,\quad W=\mu E:E.$$

For $F=\operatorname{diag}(1,1,1+\epsilon)$,

$$E_{zz}=\epsilon+\epsilon^2/2,\quad S_{zz}=2\mu E_{zz},\quad P_{zz}=(1+\epsilon)S_{zz},\quad W=\mu E_{zz}^2.$$

All other E/S/P components are zero; Fxx=Fyy=1. For intended decimal epsilon=0.003: Ezz=0.0030045, Szz=6.009 Pa, Pzz=6.027027 Pa and W=0.00902702025 Pa. This independent closed form is saved in [EXACT_AFFINE_REFERENCE.json](EXACT_AFFINE_REFERENCE.json).

Finite representation is explicitly separated: stored float64 Fzz = `1.0029999999999998916422327965847216546535491943359375000`; current float64 constitutive gives Pzz=6.027026999999707 Pa. Path B uses exactly that saved P to reproduce the prior diagnostic. Path D promotes exactly the same stored F and evaluates the analytic stress in long double. The intended-decimal reference is an additional diagnostic, never silently substituted for the historical reference. Current P differs from the closed form at the same stored F by 1.229023e-14 of the P scale, and from the intended-decimal reference by 4.857154e-14.

Formula sources: unchanged [energy.py](../../pd_clot/pfpd_candidate/energy.py), [baqp.py](../../pd_clot/pfpd_candidate/baqp.py), and [force_stencil_audit.py](../../pd_clot/pfpd_candidate/force_stencil_audit.py). The reference is independently derived here from those frozen constitutive definitions; no new external physical model is introduced.

## 4. Magnitude of reconstructed F error

All 1797 required Fbar tensors and all 50,356 distinct Ftilde tensors are included, with each physical bond counted once.

| Quantity | Max component | Max Frobenius | Mean Frobenius | RMS Frobenius | Median Frobenius | P95 Frobenius |
|---|---:|---:|---:|---:|---:|---:|
| Fbar | 3.774758284e-15 | 3.774758286e-15 | 3.528618689e-15 | 3.558844117e-15 | 3.774758284e-15 | 3.774758284e-15 |
| Ftilde | 3.774758284e-15 | 3.774758285e-15 | 2.867632020e-15 | 2.995496201e-15 | 3.058002604e-15 | 3.774758284e-15 |
| P [Pa] | 7.594813667e-12 | 7.594813668e-12 | 5.367533959e-12 | 5.701039029e-12 | 5.880667894e-12 | 7.594813667e-12 |

F is dimensionless; P errors are in Pa. Mean signed Fbar zz error = -3.528618688e-15; Ftilde zz = -2.409516439e-15. The error is mainly a systematic small negative axial bias, not zero-mean scatter. Full component means, RMS and p0/p25/p50/p75/p90/p95/p99/p100 appear in [SCIENTIFIC_RESULTS.json](SCIENTIFIC_RESULTS.json) and [F_P_ERROR_DISTRIBUTION.csv](F_P_ERROR_DISTRIBUTION.csv). Per-node/per-bond tensors are preserved in the compact [LOCAL_COEFFICIENT_ARRAYS.npz](LOCAL_COEFFICIENT_ARRAYS.npz); long-double arrays retain their dtype there.

## 5. How F error changes P

Maximum relative P Frobenius error is 1.260126040e-12; RMS is 9.459123095e-13, referenced to ||Pexact||F=6.02702699999971 Pa. Delta Pzz is negative on 48,954 bonds, positive on 1,137, and exactly zero on 265. Mean delta Pzz=-4.847947076e-12 Pa. Orientation and position bins are descriptive and are not independent causal experiments:

| Group | Interval | Bonds | Mean delta Pzz [Pa] | RMS P error [Pa] |
|---|---:|---:|---:|---:|
| abs_bond_cos_z | [0.00, 0.25] | 9134 | -7.426168e-12 | 7.444522e-12 |
| abs_bond_cos_z | [0.25, 0.50] | 14312 | -6.463695e-12 | 6.709591e-12 |
| abs_bond_cos_z | [0.50, 0.75] | 14152 | -4.484458e-12 | 5.210187e-12 |
| abs_bond_cos_z | [0.75, 1.00] | 12758 | -1.592738e-12 | 2.827141e-12 |
| midpoint_z_over_h | [-8.01, -4.00] | 1380 | -3.269072e-12 | 4.311846e-12 |
| midpoint_z_over_h | [-4.00, 0.00] | 24122 | -4.932164e-12 | 5.808062e-12 |
| midpoint_z_over_h | [0.00, 4.00] | 23463 | -4.966006e-12 | 5.751037e-12 |
| midpoint_z_over_h | [4.00, 8.01] | 1391 | -2.962514e-12 | 3.909411e-12 |

The tangent at the stored exact F is

$$\delta E=\operatorname{sym}(F^T\delta F),\quad \delta P_{lin}=\delta F S+2\mu F\delta E.$$

All bonds are evaluated. Max ||delta Praw - delta Plin||F / ||Pexact||F = 3.774741e-15; with high-precision stresses this becomes 8.482751e-18. The raw RMS remainder is 0.281% of the RMS perturbation. This supports first-order sensitivity plus arithmetic roundoff, with no resolvable nonlinear physical discrepancy at this scale.

For the axial mode $P_{zz}=\mu(a^3-a)$, $dP_{zz}/da=\mu(3a^2-1)$, so $|a P'_{zz}/P_{zz}|=335.834082$. Stress measures a small departure from the identity. A roughly 1e-15 F error can therefore produce roughly 1e-12 **relative stress error**; the stress error itself is not merely 1e-15. Transverse bonds retain more axial reconstruction error than nearly axial bonds under the existing bond correction.

## 6. How delta P propagates into Z

At fixed geometry and weights, define $Q_{ij}=I-\xi_{ij}\otimes\xi_{ij}/r_{ij}^2$ and

$$Z_i(P)=\sum_j\bar\omega_{ij}P_{ij}Q_{ij}V_j,\qquad
\Delta Z_{pred,i}=\sum_j\bar\omega_{ij}(P_{raw,ij}-P_*)Q_{ij}V_j.$$

The same float64 endpoint order, projections, neighbor volumes and weights are used for raw, exact-P, and delta-P assembly. No boundary normalization is changed. Every bond P is verified bitwise against the unchanged selective evaluator. The stable scale is ||Zbulk||F = 4.0180179999998 Pa. Ratios never divide by an individual near-zero tensor entry.

## 7. Does delta P reconstruct the observed delta Z?

**PASS**, over all 257 tensors and all nine components. Observed delta Z = Zraw - ZexactP has max Frobenius norm 4.566501658e-12 Pa and max component 4.555467115e-12 Pa. Max observed norm / bulk scale = 1.136506023e-12.

The max causal residual ||delta Zpred - delta Zobs||F is 1.082071187e-14 Pa; max component is 1.082071073e-14 Pa; relative to the stable scale it is **2.693047136e-15 <= 1e-13**. RMS residual/bulk = 1.533896293e-15. Floating-point subtraction and separately rounded assemblies need not be exactly linear. The residual is reported, not forced to zero. Its maximum is 0.237% of the maximum observed shift (different maxima may occur at different nodes).

## 8. Role of floating-point summation

Holding the complete float64 terms fixed, sequential versus balanced summation differs by at most 1.989467988e-15 of the bulk scale (7.993718186e-15 Pa); sequential versus long-double accumulation differs by 1.925220460e-15 (7.735570462e-15 Pa). Classification: **NEGLIGIBLE** relative to the 1e-12 exceedance and predeclared 1e-13 diagnostic criterion.

At the center, zz sum(abs(terms))/abs(sum(terms)) is approximately 1.00. ||sum(abs(terms))||F/||sum(terms)||F=1.060796585; sum(||term||F)/||sum(terms)||F=1.156327199. zx/zy undergo symmetry cancellation, but their final values are approximately 1e-16 Pa and do not dominate the stable tensor norm. Ratios for unsafe denominators are null, with a declared safety threshold of 8.921792e-14 Pa. Thus the proposed cancellation-amplification story is **not supported for the dominant Z norm**. The main amplification occurs in the relative F-to-P constitutive sensitivity.

## 9. High-precision reference result

This platform's long double has 64 significand bits and epsilon=1.084202172e-19; its 128-bit storage does not mean IEEE quad arithmetic. mpmath was not installed; no package was installed. Pure 3x3 algebra and summation use long double, while physical geometry, K/B, weights and projection remain the identical float64 inputs.

| Path | Max Z error / original bulk scale | Worst node ID |
|---|---:|---:|
| RAW_RECONSTRUCTED_FLOAT64 | 1.1362855097e-12 | 9837 |
| EXACT_AFFINE_P_FLOAT64 | 1.7684070888e-15 | 8330 |
| RECONSTRUCTED_F_HIGH_PRECISION_P | 1.1272864724e-12 | 9837 |
| EXACT_AFFINE_HIGH_PRECISION_REFERENCE | 1.3529099860e-14 | 8384 |

All four primary rows compare to the **same saved float64 bulk tensor**. Path C improves constitutive and accumulation precision but retains the reconstructed F error, so its approximately 1e-12 discrepancy persists. A separate C diagnostic casts high-precision P back to float64 before historical Z assembly and gives 1.126803578e-12, confirming the same conclusion without changing accumulation precision.

For D, the high-precision bulk versus historical float64 bulk differs by 1.246853972e-14 of the stable scale. Against the **matched-stress bulk**, D's spatial/weight discrepancy is only **1.060562587e-15**. Both references are reported; the numerical floor is not hidden by a precision-dependent reference offset. The optional intended-decimal F diagnostic has matched bulk error 1.060346723e-15, and historical-reference error 4.981047508e-14. Analytic versus general long-double P at the same stored F differs by 1.734723476e-18 Pa.

## 10. Is the 1.136e-12 mismatch physical or numerical?

**RECONSTRUCTED_STRESS_FLOATING_POINT_PROPAGATION**. Case A is satisfied: raw and exact-P historical values reproduce, exact-P Z is below the original coefficient tolerance, high-precision matched reference is at the fixed-input arithmetic floor, and delta-P reconstructs the raw shift within the declared causal tolerance. No structural coefficient discrepancy is detected. This conclusion is limited to the frozen, closed, intact affine stencil.

[ERROR_BUDGET.json](ERROR_BUDGET.json) and [ERROR_BUDGET.csv](ERROR_BUDGET.csv) record reconstruction, constitutive, Z assembly, summation and reference contributions. Entries use explicitly named relative Frobenius scales. The rows are not an additive error budget. Maximum relative P/F error ratio=578.788612; Z-shift/P ratio=0.901899. This is consistent with constitutive sensitivity followed by a well-conditioned stable-norm Z assembly.

## 11. Updated scientific interpretation

```text
HISTORICAL_RAW_Z_GATE: FAIL_UNCHANGED
RAW_Z_ERROR: 1.1362855097181064e-12
EXACT_AFFINE_P_Z_ERROR: 1.768407088757941e-15
HIGH_PRECISION_Z_ERROR: 1.352909986014041e-14
HIGH_PRECISION_MATCHED_REFERENCE_Z_ERROR: 1.0605625869484475e-15
MAX_F_RECONSTRUCTION_ERROR: 3.774758283725532e-15
MAX_P_ERROR: 7.594813666855771e-12
DELTA_Z_FROM_DELTA_P_RECONSTRUCTION: PASS
SUMMATION_ORDER_EFFECT: NEGLIGIBLE
Z_MISMATCH_SOURCE: RECONSTRUCTED_STRESS_FLOATING_POINT_PROPAGATION
BULK_COEFFICIENT_PHYSICS: CONSISTENT_WITH_EXACT_AFFINE_REFERENCE
BAQP_BULK_OPERATOR_DEFECT: NOT_SUPPORTED
BULK_OPERATOR_CONSISTENCY: SUPPORTED
FINITE_DOMAIN_MULTI_HOP_BOUNDARY_INFLUENCE: SUPPORTED
BOUNDARY_MODEL_DEVELOPMENT_READINESS: READY_FOR_SEPARATE_REVIEW
CURRENT_C5_C7_ACTIVATION_RESULTS: UNCHANGED
BOUNDARY_FIX: NOT_IMPLEMENTED
FRACTURE: NOT_RUN
```

These are this audit's labels. In particular `BULK_OPERATOR_CONSISTENCY=SUPPORTED` does not rewrite the prior raw-gate failure or make a universal operator theorem.

## 12. What this means for finite-domain boundary influence

The previous unchanged N27 noncentral one-hop/two-hop observations have local normalized force approximately 1.53e-2 / 3.53e-7, while the fully closed center passes force and now has an independently explained coefficient reference. Taken together, the evidence supports finite-domain multi-hop influence as the source of prior intact mechanical contamination. This is a scientific inference combining the saved stencil audit with the present reference audit; no new boundary force campaign is claimed. Center inversion symmetry remains a caveat for inferring minimal closure depth.

## 13. What remains unresolved for clot-wall modeling

The correct physical nonlocal loading/interaction-domain condition is unresolved. This audit does not establish that surface Dirichlet is mathematically invalid, that C5/C7 activation is mesh objective, that clot-wall fracture is validated, or that 3delta is the required physical collar thickness. Graph dependency reach is not a physical boundary prescription. Current C5/C7 states and results are unchanged; fracture, stability, phase and CFD were not run.

Six targeted tests only: **6 passed, 0 failed, 0 errors, 0 skipped**; [TARGETED_TESTS.log](TARGETED_TESTS.log). The tests cover analytic SVK, raw reproduction including saved representative tensors, exact-P reproduction, causal assembly, accumulation/reference precision, and unchanged historical bytes. No old test assertions were changed.

Measured numerical computation/export + final figure/test/report phase wall time: 14.364 s. CPU workers=1, BLAS threads=1, peak process-tree RSS=0.570 GiB; GPU=NOT_USED_BY_DESIGN. Resource scope and sampling method are in [RESOURCES.json](RESOURCES.json); authoring/inspection time is not computational runtime.

Scientific figures (English, PNG 600 dpi and vector PDF):

- [error_propagation_chain — PNG](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/error_propagation_chain.png) · [PDF](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/error_propagation_chain.pdf)
- [z_reference_comparison — PNG](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/z_reference_comparison.png) · [PDF](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/z_reference_comparison.pdf)
- [delta_z_closure — PNG](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/delta_z_closure.png) · [PDF](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/delta_z_closure.pdf)
- [z_accumulation_conditioning — PNG](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/z_accumulation_conditioning.png) · [PDF](../../results/pfpd_baqp_affine_coefficient_reference_audit_v1/figures/z_accumulation_conditioning.pdf)

## 14. One next development task

`NONLOCAL_BOUNDARY_FORMULATION_V1` is recommended for **separate human review**, conditional on this Case A interpretation being accepted. Define the physical load and interaction domain, then compare its traction/virtual-work target without tuning collars to match activation thresholds. This task is **not started**. No ghosts, volume constraints, force-density collars, surface corrections or boundary fixes are implemented here.

Reproduce only this local audit from the project root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/clot_affine_reference_mpl /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m scripts.affine_coefficient_reference_audit --recompute
```

Omit `--recompute` to re-export figures/report and run the six focused tests from the saved local tensors. Historical SHA256 protection is checked before and after. Stop for human review.
