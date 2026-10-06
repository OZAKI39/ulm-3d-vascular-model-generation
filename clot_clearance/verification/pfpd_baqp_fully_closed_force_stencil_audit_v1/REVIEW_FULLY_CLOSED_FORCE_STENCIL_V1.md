# PFPD_BAQP_FULLY_CLOSED_FORCE_STENCIL_AUDIT_V1

**The algorithmically closed N27 center recovers zero nodal force at the predeclared force tolerance:** 8.265455792e-15 ≤ 1e-10. **The combined scientific acceptance status is `NOT_ESTABLISHED`**, with primary root cause `INCONCLUSIVE`.

The raw Z bulk mismatch is 1.136285510e-12, exceeding the fixed 1e-12 coefficient gate by 13.629%. This exceedance is retained. It is not relabeled PASS or replaced by a constant-stress diagnostic. The fully closed force result does not support a nonzero bulk-force defect, but the user-specified conjunction needed to upgrade the formal root-cause label has not been met.

## 1. Why d ≥ δ was not enough

An intact nodal force reads incident-neighbor K inverse and Z. Neighbor Z reads endpoint-normalized weights on its own bonds; those weights read omega0 at a second neighbor, whose family reaches a third hop. An isolated own family is not an isolated force stencil. No fracture calculation was used to infer this chain.

## 2. Infinite-lattice bulk oracle

Fixed h=6.25e-05 m, δ=0.00025125 m, δ/h=4.02, ε=0.003, SVK μ=1000 Pa, λ=0, quadratic kernel. The integer-offset set uses the same `0 < |offset| h ≤ δ` inclusion as the existing neighbor construction. The kernel is zero at exactly δ. There are **256** offsets and every volume is h³.

The oracle computes omega0, K, K inverse, wbar, Z, and the incident center force using the existing r² direct formula. Pzz=6.027026999999707 Pa and the fixed local scale is **2.3543074218748854e-08 N**. Oracle force is `[0.0, 0.0, -5.113330854949596e-25]` N; norm/local=2.171904488e-17, PASS. Tensor entries, eigenvalues, offsets and weights are in [BULK_LATTICE_ORACLE.json](BULK_LATTICE_ORACLE.json).

The force scale remains |Pzz|h² for every domain and the zero-strain baseline. No local residual is normalized by Pzz L².

## 3. How the BAQP force dependency propagates

`D0={i}; D1=N(D0); D2=N(D1); D3=N(D2)` are explicit graph-walk sets, not disjoint shells. Infinite-lattice counts are **1, 256, 1,797, 5,821**, with maximum coordinate offsets 0, 4, 8, 12 cells and radii 0, 250, 500, 750 µm. Each finite graph is built separately and compared to its required untruncated counterpart; missing nodes and non-interior volumes are recorded.

The selective evaluator retains local Fbar, Ftilde and P reconstruction so its numerical operation follows `baqp.evaluate()`. This additionally requires K/Fbar on D2; those families still end in D3. Only incident batches of Ftilde/P exist in memory. No full N19/N27 bond state, dense matrix, phase state or equilibrium solve is allocated.

For each unordered bond (a,b), with ξ=Xb−Xa, Bi=ω Ka⁻¹ξ and Bj=−ω Kb⁻¹ξ, the implementation uses

$$
t_a=\bar\omega P\xi/|\xi|^2+Z_aB_i,\qquad
t_b=-\bar\omega P\xi/|\xi|^2+Z_bB_j,
$$
$$f_{ab}=(t_a-t_b)V_aV_b.$$

The lower endpoint receives +fab and the upper endpoint −fab. Both selective and full validation use lexicographic unordered bonds, including the original two endpoint-accumulation passes. No physical bond is duplicated.

## 4. Why N11 / N19 / N27 were chosen

All cubes use fixed h, L=(n−1)h, and tensor-product boundary-node volumes h³, h³/2, h³/4, h³/8. Their center distances are 5h, 9h, 13h, respectively. The largest volume-sum relative error is 2.220e-16.

| Domain | L [mm] | Center d/δ | Center force norm [N] | Norm / local scale | Verified closure with bulk volumes |
|---|---:|---:|---:|---:|---|
| N11 | 0.625 | 1.243781 | 1.904583368e-22 | 8.089781949e-15 | D1:CLOSED/D2:OPEN/D3:OPEN |
| N19 | 1.125 | 2.238806 | 1.311494861e-21 | 5.570618557e-14 | D1:CLOSED/D2:CLOSED/D3:OPEN |
| N27 | 1.625 | 3.233831 | 1.945942392e-22 | 8.265455792e-15 | D1:CLOSED/D2:CLOSED/D3:CLOSED |

## 5. Selective evaluator validation

Eight center/near-center points are selected by lattice offsets, including both signs of all axes and (1,−1,−1). Force, direct and Z correction are compared to `baqp.evaluate()` and its full diagnostic decomposition. Every needed omega0 and every D0/D1 K/Z coefficient is checked, not only the central coefficient. The force error is scaled by the fixed local force scale; coefficient errors by the infinite-lattice scalar or Frobenius norm.

| Grid | Maximum force/direct/Z error | omega0 error | K error | Z error | Gate |
|---|---:|---:|---:|---:|---|
| N11 | 7.026954975e-17 | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 | PASS |
| N13 | 7.026954975e-17 | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 | PASS |

An initial analytic-constant-P prototype exceeded the N13 Z comparison threshold (1.0895332849e-12). It was retired before N19/N27 use. Its records remain as `INITIAL_CONSTANT_P_VALIDATION_N11.json` and `INITIAL_CONSTANT_P_VALIDATION_N13.json`. The final evaluator reconstructs local F and P and passes both grids without changing the 1e-12 threshold. The initial mismatch is not hidden or treated as a physics modification.

## 6. One-hop closure result

N11 has closed D1 with interior volumes but open D2/D3. Its center passes the force-only gate while force-relevant coefficients remain strongly nonuniform. Inversion symmetry at the center can cancel a boundary-contaminated force. This is not evidence that one hop generally suffices.

## 7. Two-hop closure result

N19 closes D2 with bulk volumes, but required D3 dependencies still encounter the finite boundary. Its center again passes the force-only gate. Near-center and geometrically selected off-center results are retained in the observation table rather than using center cancellation to declare translation invariance.

## 8. Three-hop fully closed result

For the N27 center all required D0–D3 nodes exist; every required volume is h³; every D2 neighbor offset set equals the bulk oracle set. `FULL_FORCE_STENCIL_CLOSED=YES` follows from these explicit checks, not just d>3δ. Its force vector is `[-4.366447135679762e-27, -2.531806051174206e-27, 1.9459423910309976e-22]` N, norm **1.945942392e-22 N**, norm/local **8.265455792e-15**.

The nearest 3δ lattice point is at depth 12h, not the center at 13h. At 12h some required D3 nodes lie exactly on the boundary and have half volume; mere node existence does not satisfy the bulk-volume closure test. [DEPENDENCY_CLOSURE.json](DEPENDENCY_CLOSURE.json) records both kinds of closure.

## 9. Coefficient convergence to bulk values

omega0 and K maxima cover D0–D2 (K includes numerical F reconstruction); Z covers D0–D1; wbar covers all bonds incident to D0/D1. Direct-B K restricted to D0/D1 is separately retained in JSON. Matrices use relative Frobenius errors, avoiding division by exact zero entries. Center and representative D1 tensors are also saved.

| Domain | omega0 maximum relative error | K maximum relative error | Z maximum relative error | wbar maximum relative error |
|---|---:|---:|---:|---:|
| N11 | 7.521319756e-01 | 7.635614480e-01 | 1.126462144e-01 | 1.558382712e+00 |
| N19 | 2.373641723e-01 | 3.028387207e-01 | 1.400731970e-03 | 1.556209161e-01 |
| N27 | 1.665334537e-15 | 2.180741146e-15 | 1.136285510e-12 | 1.554312234e-15 |

At the closed N27 center the supplementary diagnostic holds the existing weights and geometry fixed and evaluates Z at the exact analytic affine P. Its maximum Z bulk error is **1.768407089e-15**. The reconstructed-P minus constant-P Z difference is **1.136506023e-12**, while the maximum Fbar component error is **3.774758284e-15**. These observations locate the small raw Z mismatch in floating-point stress reconstruction; they do not waive the original raw coefficient gate. See [NOISE_BASELINE.json](NOISE_BASELINE.json).

## 10. Does affine bulk force vanish?

**Yes, at the prescribed local force tolerance for the fully closed N27 center.** The combined acceptance rule also requires coefficients ≤1e-12 and is therefore **NOT_ESTABLISHED** in this run. These are distinct statements.

`MINIMUM_OBSERVED_CLOSURE_DEPTH_FOR_BULK_FORCE=NOT_ESTABLISHED`. The first center with zero force is already N11, but all centers possess inversion symmetry. Some off-center points at the same nominal closure depth have different residuals. This audit establishes the measured fully closed center force; it does not prove a universal minimum closure depth or prescribe a collar thickness.

Noise references: infinite-lattice norm/local=2.171904488e-17; zero-strain maximum force=0.000000000e+00 N; largest top/bottom paired vector defect/local=1.352600996e-13. The 1e-10 force threshold was fixed before these evaluations.

## 11. Which channel transmits finite-domain influence?

The original force is decomposed as direct + Z correction. Separately, wbar = omega/omega0_inf + (wbar−omega/omega0_inf); Z is independently accumulated for each weight component using the same reconstructed P and K. These are exact summands, not alternative models. Maximum decomposition error/local=4.260091454e-16 ≤1e-12.

N27 observations below use the bottom axis; matching reflected points were also evaluated. The three approximate horizon depths use nearest lattice layers (4h, 8h, 12h). Two extra samples at 5h and 9h provide noncentral one/two-hop closure witnesses. All entries except distance and force norm are signed axial contributions normalized by the same local force scale.

| Observation | Depth [h] | d/δ | Force norm / local | Direct fz / local | Z fz / local | Bulk-weight fz / local | Normalization-excess fz / local |
|---|---:|---:|---:|---:|---:|---:|---:|
| near_1delta_bottom | 4 | 0.995025 | 6.324071447e-02 | -4.647433677e-03 | -5.859328079e-02 | -3.543524787e-02 | -2.780546660e-02 |
| closed_1hop_bottom | 5 | 1.243781 | 1.530761376e-02 | -6.358440922e-04 | -1.467176967e-02 | -9.449101689e-03 | -5.858512070e-03 |
| near_2delta_bottom | 8 | 1.990050 | 9.248977817e-06 | -1.237798119e-13 | -9.248977693e-06 | -7.450434417e-12 | -9.248970366e-06 |
| closed_2hop_bottom | 9 | 2.238806 | 3.529230028e-07 | 2.845038396e-14 | -3.529230311e-07 | 3.592530731e-15 | -3.529230064e-07 |
| near_3delta_bottom | 12 | 2.985075 | 1.335121464e-15 | 0.000000000e+00 | 1.280223360e-15 | 1.405390995e-15 | 8.348108000e-19 |
| center | 13 | 3.233831 | 8.265455792e-15 | 2.635108116e-17 | 8.131504460e-15 | 8.151267771e-15 | -8.065939907e-18 |

Secondary mechanism label: `BOTH`. Nonzero Z and normalization-excess components identify transmission channels; their size does not make either term erroneous.

## 12. Final root cause classification

Primary: **`INCONCLUSIVE`**. Secondary: **`BOTH`**. The strict reason the primary label was not upgraded is the recorded raw Z coefficient exceedance, not a nonzero fully closed bulk force. No tolerance, old scientific result, production formula or boundary condition was changed.

Seven focused tests were run: deterministic lattice stencil; selective/full N11 equivalence; deterministic D1/D2/D3; both force decompositions; zero strain; N27 geometry/volume closure; infinite-lattice symmetry. The [test log](TARGETED_TESTS.log) contains the actual results. Full pytest and historical cleanup QA were not run. An initial new symmetry-test assertion incorrectly used a dimensional 1e-30 bound for exact-zero matrix entries; the final assertion uses the predeclared 1e-12 relative tensor norm. The initial log is retained as `TARGETED_TESTS_INITIAL.log`; no historical test or scientific acceptance limit was changed.

## 13. Implication for existing C5/C7 activation drift

The intact mechanical evidence is consistent with finite-domain multi-hop coefficients contaminating a mechanical driver field. Under the full predeclared acceptance rule, the formal root cause remains unresolved pending the coefficient-reference issue. Existing C5/C7 activation values and accepted states remain **UNCHANGED**. No claim of fracture mesh objectivity, a correct physical nonlocal BC, or an incorrect surface Dirichlet condition follows from this audit. No activation, phase, stability, Gc/Yc or CFD calculation was performed.

## 14. What should be developed next

First obtain human review of the exact-affine versus reconstructed-stress coefficient reference and its fixed tolerance; keep this run's exceedance intact. A focused affine-coefficient-reference audit can settle that numerical acceptance issue. A subsequent, separately approved boundary-development task should compare a physically specified nonlocal loading/interaction-domain formulation against its intended traction and virtual work, at fixed h and δ. Do not select a collar thickness by forcing C5/C7 activation thresholds to agree. No collar, ghost region, force-density loading or volume constraint was implemented here.

CPU workers=1, BLAS threads=1. Recorded audit computation, targeted tests and scientific figure exports: 24.384 s; sampled process-tree peak RSS=0.590 GiB; minimum sampled available RAM fraction=62.717%. GPU=`NOT_USED_BY_DESIGN`. The monitor samples at 50 ms; peaks are reported as sampled measurements, not exact continuous-time maxima.

Reproduce from the `clot_clearance` root using the existing environment:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/clot_closed_stencil_mpl /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m scripts.fully_closed_force_stencil_audit
```

Scientific figures only (English; PNG 600 dpi and PDF):

- [closure_depth_force_residual — PNG](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/closure_depth_force_residual.png) · [PDF](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/closure_depth_force_residual.pdf)
- [coefficient_convergence_to_bulk — PNG](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/coefficient_convergence_to_bulk.png) · [PDF](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/coefficient_convergence_to_bulk.pdf)
- [n27_force_vs_boundary_distance — PNG](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/n27_force_vs_boundary_distance.png) · [PDF](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/n27_force_vs_boundary_distance.pdf)
- [force_stencil_dependency_graph — PNG](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/force_stencil_dependency_graph.png) · [PDF](../../results/pfpd_baqp_fully_closed_force_stencil_audit_v1/figures/force_stencil_dependency_graph.pdf)

```text
INFINITE_LATTICE_BULK_ORACLE: PASS
SELECTIVE_EVALUATOR: VERIFIED
N11_CLOSURE: D1:CLOSED/D2:OPEN/D3:OPEN
N19_CLOSURE: D1:CLOSED/D2:CLOSED/D3:OPEN
N27_CLOSURE: D1:CLOSED/D2:CLOSED/D3:CLOSED
FULL_FORCE_STENCIL_CLOSED: YES
ACTUAL_REQUIRED_GRAPH_DEPTH: 3
MINIMUM_OBSERVED_CLOSURE_DEPTH_FOR_BULK_FORCE: NOT_ESTABLISHED
N27_CENTER_FORCE_N: 1.945942391685588e-22
N27_CENTER_FORCE_VECTOR_N: [-4.366447135679762e-27, -2.531806051174206e-27, 1.9459423910309976e-22]
N27_CENTER_FORCE_LOCAL_NORMALIZED: 8.265455792242755e-15
N27_CENTER_FORCE_ONLY_GATE: PASS
N27_COEFFICIENT_GATE: FAIL
BULK_AFFINE_FORCE_PATCH_TEST: NOT_ESTABLISHED
OMEGA0_BULK_MATCH: 1.6653345369377348e-15
K_BULK_MATCH: 2.180741146465773e-15
Z_BULK_MATCH: 1.1362855097181066e-12
ROOT_CAUSE_PRIMARY: INCONCLUSIVE
ROOT_CAUSE_SECONDARY: BOTH
CURRENT_C5_C7_ACTIVATION_RESULTS: UNCHANGED
BOUNDARY_FIX: NOT_IMPLEMENTED
FRACTURE: NOT_RUN
C7_STABILITY: NOT_RUN
```

Stop for human review. No next development task was started.
