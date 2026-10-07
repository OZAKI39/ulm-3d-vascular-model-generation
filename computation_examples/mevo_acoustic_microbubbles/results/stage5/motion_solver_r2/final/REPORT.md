# Stage 5A-R2 measured outcome

Outcome: **VALIDATION_INCOMPLETE_OR_FAILED**. Promotion: **FAIL**.

Executed source: `7f83f802c6cc435918ea25d69dcbb0f4f4dacce9`. Final archive commit is the commit containing this report; it is not the executed source identity.

Fresh bound tests: 178 passed, 1 retained old xfail. Original-QP regressions: PASS; holdout: PASS; independent final-equation audit: PASS; final-kernel manufactured convergence: PASS.

Native packages: proxsuite/ProxQP 0.7.2, OSQP 1.1.3, Clarabel 0.11.1; float64. The selected solver is bounded component ProxQP, with certified native alternatives. Internal status never replaces the fixed original KKT check.

The final method uses explicit acoustic/flow forcing and endpoint implicit resistance (S1), with bounded start-state S0 recovery at nonsmooth near-wall normal transitions. Every adopted substep has step-doubling and continuous pair/WALL checks. The old tolerance remains 1e-11 m; no projection is applied.

| Case | Kind | Target (ms) | Completed (ms) | Accepted / rejected | Min / median dt (s) | Total wall (s) | Status |
|---|---|---:|---:|---:|---|---:|---|
| B2 | nominal repeat 0 | 1 | 1 | 1988 / 104 | 5e-09 / 7.0247e-08 | 481.754 | COMPLETE |
| B3 | nominal repeat 0 | 1 | 1 | 1872 / 140 | 4.8149e-09 / 8.0561e-08 | 532.820 | COMPLETE |
| B2 | refined repeat 0 | 1 | 1 | 3842 / 169 | 2.1618e-09 / 3.6475e-08 | 778.298 | COMPLETE |
| B3 | refined repeat 0 | 1 | 1 | 3542 / 135 | 2.0254e-09 / 4.4744e-08 | 798.992 | COMPLETE |
| B2 | reference repeat 0 | 1 | 1 | 3242 / — | 3.5597e-09 / 2.9387e-07 | 355.096 | COMPLETE |
| B3 | reference repeat 0 | 1 | 1 | 3984 / — | 1.1847e-09 / 2.3454e-07 | 508.933 | COMPLETE |
| B2 | nominal repeat 0 | 5 | 4.3142345 | 3908 / 202 | 7.6294e-13 / 2.9227e-07 | 775.807 | FAIL |
| B3 | nominal repeat 0 | 5 | 5 | 3478 / 140 | 4.8149e-09 / 2.5e-06 | 794.075 | COMPLETE |
| B2 | refined repeat 0 | 5 | 4.3111491 | 6918 / 248 | 7.6294e-13 / 2.7402e-07 | 1230.813 | FAIL |
| B3 | refined repeat 0 | 5 | 5 | 6730 / 135 | 2.0254e-09 / 1.25e-06 | 1292.447 | COMPLETE |
| B2 | reference repeat 0 | 5 | 3.9264049 | 13854 / — | 1.2128e-09 / 2.6865e-07 | 1800.743 | INCOMPLETE |
| B3 | reference repeat 0 | 5 | 2.8946426 | 11724 / — | 2.3955e-10 / 2.2404e-07 | 1800.798 | INCOMPLETE |
| B2 | nominal repeat 1 | 1 | 1 | 2002 / 104 | 5e-09 / 7.0355e-08 | 482.953 | COMPLETE |
| B3 | nominal repeat 1 | 1 | 1 | 1872 / 140 | 4.8149e-09 / 8.0561e-08 | 530.718 | COMPLETE |
| B2 | reference repeat 1 | 1 | 1 | 3242 / — | 3.5597e-09 / 2.9387e-07 | 351.298 | COMPLETE |
| B3 | reference repeat 1 | 1 | 1 | 3984 / — | 1.1847e-09 / 2.3454e-07 | 506.426 | COMPLETE |
| B2 | nominal repeat 2 | 1 | 1 | 1986 / 104 | 5e-09 / 7.0161e-08 | 477.290 | COMPLETE |
| B3 | nominal repeat 2 | 1 | 1 | 1870 / 127 | 4.8149e-09 / 8.3275e-08 | 520.912 | COMPLETE |
| B2 | reference repeat 2 | 1 | 1 | 3242 / — | 3.5597e-09 / 2.9387e-07 | 353.539 | COMPLETE |
| B3 | reference repeat 2 | 1 | 1 | 3984 / — | 1.1847e-09 / 2.3454e-07 | 505.977 | COMPLETE |

## Accepted-state numerical evidence

| Case / target / kind | QP / original S1 KKT max | Resistance residual max | Force residual max (N) | Pair / WALL swept minimum (m) | Projection (m) |
|---|---|---|---|---|---|
| B2 / 1 ms / nominal | 9.97677e-10 / 9.96187e-10 | 3.58932e-14 | 4.91885e-16 | -9.99757e-12 / -9.99271e-12 | 0 |
| B3 / 1 ms / nominal | 9.9875e-10 / 9.98791e-10 | 5.07033e-14 | 5.33647e-16 | -9.9499e-12 / -9.97209e-12 | 0 |
| B2 / 1 ms / refined | 9.95764e-10 / 9.98957e-10 | 4.15937e-14 | 4.90374e-16 | -9.99757e-12 / -9.97508e-12 | 0 |
| B3 / 1 ms / refined | 9.95214e-10 / 9.99567e-10 | 4.31856e-14 | 5.36839e-16 | -9.9499e-12 / -9.9996e-12 | 0 |
| B2 / 5 ms / nominal | 9.97677e-10 / 9.994e-10 | 3.98674e-14 | 5.00056e-16 | -9.99757e-12 / -9.99999e-12 | 0 |
| B3 / 5 ms / nominal | 9.9875e-10 / 9.98103e-10 | 4.40086e-14 | 5.33192e-16 | -9.9499e-12 / -9.97209e-12 | 0 |
| B2 / 5 ms / refined | 9.95764e-10 / 9.98084e-10 | 3.85068e-14 | 4.89232e-16 | -9.99757e-12 / -9.99999e-12 | 0 |
| B3 / 5 ms / refined | 9.95214e-10 / 9.98044e-10 | 4.64963e-14 | 5.37355e-16 | -9.9499e-12 / -9.9996e-12 | 0 |

These are accepted states only. Complete LOCAL_ERROR rejections and terminal/nonconverged candidates are separately tabulated in RUN_AUDIT.json; failed candidates never contribute to accepted-state maxima. S0, if used, has a start-equation certificate and no endpoint nonlinear KKT. Raw negative gaps inside the unchanged 1e-11 m tolerance are retained.

## Scientific comparison

| Duration / case | Comparison | Rg endpoint / series (%) | Q1 / Q2 / Q3 (%) | Centroid / budget (nm) | Status |
|---|---|---|---|---|---|
| 1 ms / B2 | self_convergence | 0.007709 / 0.0088878 | 0.099097 / 0.061567 / 0.058363 | 5.0195 / 178 | PASS |
| 1 ms / B2 | legacy_comparison | 0.98547 / 0.78805 | 1.7536 / 0.3612 / 1.1698 | 9.6749 / 178 | PASS |
| 1 ms / B2 | refined_legacy_comparison | 0.99325 / 0.79026 | 1.8546 / 0.29945 / 1.2288 | 4.6557 / 178 | PASS |
| 1 ms / B3 | self_convergence | 5.4905e-08 / 0.00691 | 2.1682e-05 / 9.5944e-06 / 0.00038044 | 0.056178 / 178 | PASS |
| 1 ms / B3 | legacy_comparison | 0.011339 / 0.012476 | 0.2239 / 0.013116 / 0.0065219 | 37.741 / 178 | PASS |
| 1 ms / B3 | refined_legacy_comparison | 0.011339 / 0.010777 | 0.22388 / 0.013106 / 0.0069024 | 37.69 / 178 | PASS |
| 5 ms / B3 | self_convergence | 6.2177e-07 / 0.0031066 | 5.2108e-06 / 0.0003223 / 0.00024058 | 0.054715 / 178 | PASS |

The original group membership and all particles are preserved. Rg and centroid have different units/thresholds. Comparisons use actual states at the common clock, with explicit sample-age bounds; no coordinate interpolation. Frozen Stage4L is a numerical reference, not experimental truth.
Passing aggregate group metrics does not mean every individual trajectory is identical. The reports also retain maximum same-ID position difference, rigid-fit residual and neighbor-graph changes; R6/R7 display them without geometric smoothing.

V3 ordinary B2 remains BASELINE_DISAGREEMENT (1 ms Q1 difference 2.20205%). V4/V4A use a documented tighter time-error configuration of the unchanged kernel. V4's rejected entry launches remain recorded; V4A fixes entry routing, with fresh bound tests and explicitly reused V4 holdouts. Earlier results are never relabelled as later configuration runs.
Optional [partial-prefix comparisons](PARTIAL_COMPARISONS.json) describe actual common intervals when a 5 ms run stops early. They retain exact source frames and sampling uncertainty, and never upgrade a complete-duration gate.
[Same-clock prefix refinement analysis](PARTIAL_BASELINE_DISAGREEMENT.json) checks whether a legacy disagreement persists while the new method converges on exactly that comparison clock, including conservative sampling bounds. BASELINE_DISAGREEMENT remains separate from uncompleted full-duration gates; neither reference truth nor new-method invalidity is inferred automatically.
When a fresh reference cannot finish within budget, [archived full-duration science](ARCHIVED_SCIENCE.json) separately compares complete new 5 ms trajectories against the immutable original Stage4L data. These supplementary comparisons do not change the frozen campaign or use historical runtime as a speedup denominator.

### Supplementary actual common prefixes (not complete 5 ms gates)

| Comparison | Actual common duration (ms) | Rg endpoint / series (%) | Q1 / Q2 / Q3 (%) | Centroid / budget (nm) | Scoped status |
|---|---:|---|---|---|---|
| B2_0.005_v4a_partial_legacy | 3.9264049 | 0.9874 / 0.94055 | 1.7579 / 0.3627 / 1.1673 | 275.59 / 178 | FAIL |
| B2_0.005_v4a_partial_self | 4.3111491 | 0.0077235 / 0.0080039 | 0.098822 / 0.061448 / 0.058126 | 1.232 / 178 | PASS |
| B2_0.005_v4a_partial_refined_legacy | 3.9264049 | 0.9952 / 0.94711 | 1.8585 / 0.30107 / 1.2262 | 274.36 / 178 | FAIL |
| B3_0.005_v4a_partial_legacy | 2.8946426 | 0.011702 / 0.010504 | 0.21457 / 0.020041 / 0.01694 | 13.706 / 178 | PASS |
| B3_0.005_v4a_partial_refined_legacy | 2.8946426 | 0.011714 / 0.0098294 | 0.2146 / 0.020237 / 0.016954 | 13.938 / 178 | PASS |

### Supplementary archived-reference science (not fresh timing)

| Comparison | Actual common duration (ms) | Rg endpoint / series (%) | Q1 / Q2 / Q3 (%) | Centroid / budget (nm) | Scoped status |
|---|---:|---|---|---|---|
| B3_0.005_v4a_archived_legacy | 5 | 0.0023936 / 0.0076894 | 0.0026906 / 0.0023102 / 0.0052069 | 0.49413 / 178 | PASS |
| B3_0.005_v4a_archived_refined_legacy | 5 | 0.0023942 / 0.0066462 | 0.0026958 / 0.0026325 / 0.0049664 | 0.45498 / 178 | PASS |

## Performance

[Post-timing scientific repeatability](REPEAT_AUDIT.json): **PASS**. All measured repeats remain in the timing ledger. Adaptive step counts alone are not evidence of identical trajectories; repeated runs are compared against the first same-method run and, for new trajectories, the fresh first legacy reference under the unchanged scientific thresholds.

- B2 1 ms: completion speedup 0.7338586543837424; new/reference repetitions 3/3; three-repeat sufficiency True.
- B3 1 ms: completion speedup 0.954227380396278; new/reference repetitions 3/3; three-repeat sufficiency True.
- B2 5 ms: completion speedup None; new/reference repetitions 0/0; three-repeat sufficiency False.
- B3 5 ms: completion speedup None; new/reference repetitions 1/0; three-repeat sufficiency False.

Only completed equal-duration runs under the same environment enter speedup. Total wall includes process startup, failed trials, nonlinear work, I/O and monitoring. Raw error/wall subcounters overlap; R8 uses a mutually exclusive subset plus the measured remainder. Native solvers do not expose complete factorization/inner iteration histories through all APIs.

## Blockers

- B2 5 ms: numerical=False, science=NOT_EVALUABLE, frozen campaign=INCOMPLETE_OR_FAIL, repeated science=NOT_RUN
- B3 5 ms: numerical=True, science=INCOMPLETE, frozen campaign=INCOMPLETE_OR_FAIL, repeated science=NOT_RUN
- Both 5 ms comparable >=5x timing gates are not established
- B2 BASELINE_DISAGREEMENT on the actual 3.9264049 ms common prefix despite same-clock refinement; separate from uncompleted 5 ms endpoint gates

## Evidence and reproduction

- [Machine-readable gate](MOTION_SOLVER_R2_GATE.json), [all run audits](RUN_AUDIT.json), [performance](PERFORMANCE.json).
- [Figure and data index](../figures/index.md). Actual accepted/rejected ledgers, source snapshots and checkpoint files retain each run ID.
- Numerical code and environment are sealed by VALIDATION_PROTOCOL_V4A.json and TEST_GATE_V4A.json; earlier failed holdouts and campaigns remain unchanged.
- Formal reproduction uses the exact commands in CAMPAIGN_V4A/RESULT.json, in a separate output root with the recorded scientific inputs and dependency lock. Explicitly reused runs keep their earlier commands, run IDs and source commits. Do not rerun the freeze writer over an existing protocol or overwrite immutable run directories.
- Postprocessing: `PYTHONPATH=src:. python -m scripts.stage5a_r2.validation_figures`, `PYTHONPATH=src:. python -m scripts.stage5a_r2.motion_visuals --video --duration 0.001`, `PYTHONPATH=src:. python -m scripts.stage5a_r2.delivery`.
- [Reproduction details](../../../../scripts/stage5a_r2/README.md), including exact source identities, external input files and the immutable-output policy.
- Large scientific input fields, full HDF5 trajectories, source snapshots and MP4 exports are indexed by actual path, bytes and SHA256 in ARTIFACT_MANIFEST.json; they are not all embedded in Git.

Known root causes and bounded revisions are documented in `docs/STAGE5A_R2_ROOT_CAUSE_ANALYSIS.md` and `docs/STAGE5A_R2_MOTION_DEVELOPMENT.md`. Passing numerical checks would not validate the dense classical secondary-force model or a physical aggregation mechanism.
