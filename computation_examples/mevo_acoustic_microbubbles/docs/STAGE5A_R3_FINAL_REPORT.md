# Stage 5A-R3 final development and validation report

Outcome: **CONTINUATION_RECOVERED_ACCURACY_UNRESOLVED / NUMERICAL_RECOVERED_PERFORMANCE_PENDING**. Promotion: **FAIL**.

B2 continuous motion has been repaired without changing forces, radii, the WALL mesh, lubrication coefficients or tolerances. Independent global all-active accuracy remains NOT_EVALUABLE. The measured performance verdict is separate from numerical completion.

## Identity and scope

- Branch: `dev/stage5a-r3-wall-accuracy-performance`.
- Actual R2 archive base: `07db8ed35f3a460a4ba315e42b15d257f8518b0a`.
- V1 executable: `fa2fba5608037ed45143c8ede501cadf415e9180`.
- V2 executable: `b5079ad19be81ffe9e3d4f5b7c548717b3d43d67`.
- Report generation source: `b97da4bd408d60ae623d8ca3c87d3b4e2144dc49`; resolve the containing Git commit for the final archive SHA.
- Server: `/workspace/mevo_acoustic_microbubbles_stage5a_r3`; immutable executable trees are `.execution/<SHA>`.
- Results: `results/stage5/motion_solver_r3`; no old branch, checkpoint, failed gate or stopped service is resumed or overwritten.

V2 changes B3 QP resource organization to component serial, retaining B2 global QP. BLAS and QP workers remain one, and CuPy supplies the pair force. Equivalent B2 V1 executions retain their original identity and cost; no budget reset or relabeling occurred. Formal runs begin at the original 3 s state and extend by at most 5 ms.

## Actual completion and accepted numerical safety

| Case | Duration (ms) | Method | Repeat | Status | Process-inclusive wall (s) | QP / endpoint KKT max | Pair / WALL certified bounds (m) |
|---|---:|---|---:|---|---:|---:|---|
| B2 | 1 | nominal | 0 | COMPLETE | 423.015 | 2.66834e-10 / 9.95134e-10 | -9.99747e-12 / -9.1416e-19 |
| B2 | 5 | nominal | 0 | COMPLETE | 728.034 | 2.66834e-10 / 9.95134e-10 | -9.99747e-12 / -3.78391e-18 |
| B2 | 5 | refined | 0 | COMPLETE | unavailable | 3.98934e-10 / 9.99953e-10 | -9.99747e-12 / -1.14775e-17 |
| B3 | 1 | nominal | 0 | COMPLETE | 506.293 | 9.9875e-10 / 9.99189e-10 | -9.94979e-12 / -4.38532e-13 |
| B3 | 5 | nominal | 0 | COMPLETE | 758.302 | 9.9875e-10 / 9.98511e-10 | -9.94979e-12 / -4.38556e-13 |
| B3 | 5 | refined | 0 | COMPLETE | 1271.87 | 9.90784e-10 / 9.98016e-10 | -9.94979e-12 / -4.38556e-13 |
| B2 | 5 | reference | 0 | COMPLETE | 2472.19 | unavailable / unavailable | unavailable / unavailable |
| B3 | 5 | reference | 0 | COMPLETE | 3572.38 | unavailable / unavailable | unavailable / unavailable |
| B2 | 5 | nominal | 1 | COMPLETE | 726.045 | 2.66834e-10 / 9.95134e-10 | -9.99747e-12 / -3.78391e-18 |
| B3 | 5 | nominal | 1 | COMPLETE | 796.683 | 9.9875e-10 / 9.99638e-10 | -9.94979e-12 / -4.38558e-13 |
| B2 | 5 | reference | 1 | COMPLETE | 2489.18 | unavailable / unavailable | unavailable / unavailable |
| B3 | 5 | reference | 1 | COMPLETE | 3583.34 | unavailable / unavailable | unavailable / unavailable |
| B2 | 5 | nominal | 2 | COMPLETE | 723.288 | 2.66834e-10 / 9.95134e-10 | -9.99747e-12 / -3.78391e-18 |
| B3 | 5 | nominal | 2 | COMPLETE | 749.791 | 9.9875e-10 / 9.98205e-10 | -9.94979e-12 / -4.38558e-13 |
| B2 | 5 | reference | 2 | COMPLETE | 2457.11 | unavailable / unavailable | unavailable / unavailable |
| B3 | 5 | reference | 2 | COMPLETE | 3555.75 | unavailable / unavailable | unavailable / unavailable |
| B2 | 1 | reference | 0 | COMPLETE | 350.136 | unavailable / unavailable | unavailable / unavailable |
| B3 | 1 | reference | 0 | COMPLETE | 503.985 | unavailable / unavailable | unavailable / unavailable |
| B3 | 1 | r2 | 0 | COMPLETE | 521.441 | 9.9875e-10 / 9.9818e-10 | -9.9499e-12 / -9.97209e-12 |
| B3 | 5 | r2 | 0 | COMPLETE | 784.542 | 9.9875e-10 / 9.99839e-10 | -9.9499e-12 / -9.97209e-12 |
| B3 | 5 | r2 | 1 | COMPLETE | 791.006 | 9.9875e-10 / 9.96581e-10 | -9.9499e-12 / -9.97209e-12 |
| B3 | 5 | r2 | 2 | COMPLETE | 785.385 | 9.9875e-10 / 9.99214e-10 | -9.9499e-12 / -9.97209e-12 |

A missing process-inclusive time is not replaced by body timing. In particular, B2 V1 refined completed in 1152.959 s of measured run-body time while its parent scheduler was held. Body time plus 50 s is a conservative budget charge only. Legacy references use different projection/contact semantics, so absent R3 safety fields remain unavailable.

S1 uses endpoint resistance with the frozen forcing/contact discretization. S0 uses start resistance and has no endpoint nonlinear KKT: NOT_APPLICABLE/null. QP KKT ≤1e−9, resistance residual ≤1e−11, pair/WALL tolerance 1e−11 m and zero projection remain unchanged. Rejected candidate residuals are separate from accepted maxima. HDF5 active IDs, finite coordinates, strict time increase and accepted-log correspondence are audited independently of the run summary.

Accepted counts below refer to saved half steps: each successful step-doubling macro step adopts two half steps. The unchanged recovery controller allows eight accepted S0 macro steps per invocation and at most 1024 invocations. Consecutive S0 sequences in the audit are measured from the adopted history, separately from controller invocation counts.

| Case / run | S0 half steps / duration (µs) | S1 half steps / duration (ms) | State integrity |
|---|---:|---:|---|
| `20261007T101533737421Z_B2_S1_0.001s_formal_R3_V1_0.001_nominal_repeat0` | 0 / 0 | 2022 / 1 | PASS |
| `20261007T103209502094Z_B2_S1_0.005s_formal_R3_V1_0.005_nominal_repeat0` | 144 / 3.6 | 3856 / 4.9964 | PASS |
| `20261007T110855321948Z_B2_S1_0.005s_formal_R3_V1_0.005_refined_repeat0` | 64 / 0.8 | 7366 / 4.9992 | PASS |
| `20261007T120312004897Z_B3_S1_0.001s_formal_R3_V2_0.001_nominal_repeat0` | 0 / 0 | 1872 / 1 | PASS |
| `20261007T121138315056Z_B3_S1_0.005s_formal_R3_V2_0.005_nominal_repeat0` | 0 / 0 | 3482 / 5 | PASS |
| `20261007T122416617767Z_B3_S1_0.005s_formal_R3_V2_0.005_refined_repeat0` | 0 / 0 | 6758 / 5 | PASS |
| `20261007T142613106084Z_B2_S1_0.005s_formal_R3_V2_0.005_nominal_repeat1` | 144 / 3.6 | 3856 / 4.9964 | PASS |
| `20261007T143819169614Z_B3_S1_0.005s_formal_R3_V2_0.005_nominal_repeat1` | 0 / 0 | 3494 / 5 | PASS |
| `20261007T163248417997Z_B2_S1_0.005s_formal_R3_V2_0.005_nominal_repeat2` | 144 / 3.6 | 3856 / 4.9964 | PASS |
| `20261007T164451685810Z_B3_S1_0.005s_formal_R3_V2_0.005_nominal_repeat2` | 0 / 0 | 3476 / 5 | PASS |

## Confirmed B2 failure and numerical repair

The responsible bubble is ID 152. Its actual nearest WALL feature switches between finite faces 3178 and 13018, with a normal jump around 0.1343. Cold queries exactly match the stored-query gap/normal, refuting a stale-cache explanation. A 1.526 ps candidate is safe under the exact finite-triangle predicate but fails the old global roundoff guard; 1 ns and 50 ns candidates are genuinely unsafe. Therefore simply removing the guard is insufficient. The minimum-progress failure is an algorithmic 1 ps limit, not an unrepresentable timestamp.

R3 rebuilds finite reachable hard-contact features at the actual position and uses directed-rounding separation bounds with exact binary64-rational segment/triangle fallback. Face/edge/vertex coverage, signs and open ports are tested. The nearest-wall lubrication model is unchanged; extra contact normals do not add extra wall resistance. No wall smoothing, projection, bubble deletion or synthetic accepted frames are used. See the root-cause report, G1–G3 and the full original failed ledgers in the review directory.

## All-active precision and reference lineage

Every comparison preserves original IDs and all 186/189 active particles, with the original 75/91-ID group and its complement also reported. Raw, translation-removed and rigid-fit errors are retained. Independent global trajectory bounds are NOT_ESTABLISHED, so these comparisons cannot promote a per-bubble reference.

| Case / named comparison | Common duration (ms) | Endpoint raw RMS / max (nm) | Sampled maximum (nm) | Worst sampled ID / time (s) | Common-time bound available |
|---|---:|---|---:|---|---|
| B2 / R3_V2_0.001_reference | 1 | 784.95 / 3055.68 | 3085.65 | 187 / 3.000670000 | True |
| B2 / R3_V2_0.001_old_fresh_vs_new_fresh | 1 | 0 / 0 | 0 | 0 / 3.000000000 | True |
| B2 / R3_V2_0.001_R2 | 1 | 0.0170276 / 0.0676669 | 0.756074 | 82 / 3.000810000 | True |
| B2 / R3_V2_0.005_refined | 5 | 31.0879 / 233.979 | 669.608 | 211 / 3.000070000 | True |
| B2 / R3_V2_0.005_reference | 5 | 1057.41 / 3890.95 | 3891.95 | 187 / 3.004920000 | True |
| B2 / R3_V2_0.005_old_fresh_vs_new_fresh | 3.9264 | 0 / 0 | 0 | 0 / 3.000000000 | True |
| B2 / R3_V2_0.005_archived | 5 | 20.8296 / 72.7014 | 718.265 | 211 / 3.000070000 | False |
| B2 / R3_V2_0.005_archived_vs_fresh | 5 | 1051.93 / 3894.18 | 3896 | 187 / 3.004990000 | False |
| B2 / R3_V2_0.005_R2 | 4.31423 | 1.37409 / 8.83457 | 8.83457 | 24 / 3.004314234 | True |
| B2 / R3_V2_0.005_repeat1 | 5 | 0 / 0 | 0 | 0 / 3.000000000 | True |
| B2 / R3_V2_0.005_repeat1_reference | 5 | 1057.41 / 3890.95 | 3891.95 | 187 / 3.004920000 | True |
| B2 / R3_V2_0.005_repeat2 | 5 | 0 / 0 | 0 | 0 / 3.000000000 | True |
| B2 / R3_V2_0.005_repeat2_reference | 5 | 1057.41 / 3890.95 | 3891.95 | 187 / 3.004920000 | True |
| B3 / R3_V2_0.001_reference | 1 | 200.58 / 1161.56 | 1162.24 | 112 / 3.000990000 | True |
| B3 / R3_V2_0.001_protected_R2 | 1 | 2.38201 / 8.97058 | 33.632 | 117 / 3.000250000 | True |
| B3 / R3_V2_0.001_old_fresh_vs_new_fresh | 1 | 0 / 0 | 0 | 1 / 3.000000000 | True |
| B3 / R3_V2_0.001_V1_global_vs_V2_serial | 1 | 1.35366e-05 / 5.30876e-05 | 0.011144 | 33 / 3.000450000 | True |
| B3 / R3_V2_0.001_R2 | 1 | 0.00450307 / 0.0267457 | 0.101456 | 212 / 3.000300000 | True |
| B3 / R3_V2_0.005_refined | 5 | 112.894 / 540.454 | 662.128 | 117 / 3.000210000 | True |
| B3 / R3_V2_0.005_reference | 5 | 209.495 / 1167.47 | 1169.1 | 112 / 3.003000000 | True |
| B3 / R3_V2_0.005_protected_R2 | 5 | 0.0155841 / 0.045547 | 25.6694 | 117 / 3.000250000 | True |
| B3 / R3_V2_0.005_old_fresh_vs_new_fresh | 2.89464 | 0 / 0 | 0 | 1 / 3.000000000 | True |
| B3 / R3_V2_0.005_archived | 5 | 191.993 / 1165.85 | 1167.61 | 112 / 3.003260000 | False |
| B3 / R3_V2_0.005_archived_vs_fresh | 5 | 83.4195 / 379.344 | 577.26 | 117 / 3.000230000 | False |
| B3 / R3_V2_0.005_V1_global_vs_V2_serial | 5 | 2.24842 / 8.67349 | 33.6326 | 117 / 3.000250000 | True |
| B3 / R3_V2_0.005_R2 | 5 | 2.25212 / 8.68884 | 33.6327 | 117 / 3.000250000 | True |
| B3 / R3_V2_0.005_repeat1 | 5 | 2.251 / 8.67497 | 33.6043 | 117 / 3.000250000 | True |
| B3 / R3_V2_0.005_repeat1_reference | 5 | 209.268 / 1167.47 | 1169.11 | 112 / 3.003000000 | True |
| B3 / R3_V2_0.005_repeat1_protected_R2 | 5 | 0.011385 / 0.0289241 | 0.252848 | 212 / 3.003010000 | True |
| B3 / R3_V2_0.005_repeat2 | 5 | 0.00289069 / 0.0115926 | 0.0828548 | 60 / 3.003200000 | True |
| B3 / R3_V2_0.005_repeat2_reference | 5 | 209.496 / 1167.47 | 1169.1 | 112 / 3.003000000 | True |
| B3 / R3_V2_0.005_repeat2_protected_R2 | 5 | 0.0172726 / 0.0455461 | 25.6687 | 117 / 3.000250000 | True |

B2 nominal/refined +5 ms endpoint RMS/max are 31.088/233.979 nm, exceeding frozen 17.8/89 nm targets. The first sampled maximum crossing is bracketed by +30..40 µs (ID 82); the worst sampled error is 669.608 nm at +70 µs (ID 211). This early event is distinct from the +4.31 ms WALL stall. Both selected IDs have no WALL face switch in +30..80 µs; no S0 or fallback occurs there.

Six all-active local +20 µs replays use the actual nominal/refined +30 µs states at three precision levels. Adjacent finest RMS/max differences are 11.004/50.914 nm and 11.430/82.315 nm. Differences between the two actual starts grow from 4.888 nm RMS to 48.655 nm RMS under the finest local settings. Near-contact rearrangement is correlated with this sensitivity, but no unique global-error cause or preferred truth is established. Cold-cache/controller scope is explicit. See B2_EARLY_PRECISION_AUDIT.json and G5b.

Archived Stage4L and fresh frozen Stage4L remain separate numerical lineages. Controlled tests identify early sensitivity to their 500/10 µs clock cutoff; same-input equation/QP checks agree to floating-point precision. This does not resolve the complete archived/fresh trajectory disagreement. At unmatched legacy times, a neighboring-frame chord is only a displacement scale unless a straight recorded path is established; consistency gates are withheld when that sampling bound is unavailable. Actual common stored frames are not interpolated.

B2: all 13855 stored old-reference frames through 3.003926404949 s were compared with the complete new reference prefix; saved-state replay status **PASS**.
B3: all 11725 stored old-reference frames through 3.002894642632 s were compared with the complete new reference prefix; saved-state replay status **PASS**.
This exhaustive saved-frame check covers time, full coordinates, active masks, IDs and radii. It establishes reference replay consistency over the old recorded interval, not physical trajectory accuracy.

## Measured performance and resources

| Case | Complete new / reference 5 ms repeats | New wall samples (s) | Reference wall samples (s) | Median completion speedup | Verdict |
|---|---|---|---|---:|---|
| B2 | 3 / 3 | 728.034, 726.045, 723.288 | 2472.19, 2489.18, 2457.11 | 3.40501 | FAIL |
| B3 | 3 / 3 | 758.302, 796.683, 749.791 | 3572.38, 3583.34, 3555.75 | 4.71103 | FAIL |

| Case | Complete protected R2 5 ms repeats | Protected R2 wall samples (s) | Median R2 / R3 speedup | Verdict |
|---|---:|---|---:|---|
| B2 | 0 |  | unavailable | NOT_EVALUABLE |
| B3 | 3 | 784.542, 791.006, 785.385 | 1.03572 | FAIL |

| Case, +1 ms | R3 wall (s) | Frozen Stage4L wall (s) | Frozen / R3 ratio, non-regression | R2 wall (s) | R2 / R3 ratio, non-regression |
|---|---:|---:|---|---:|---|
| B2 | 423.015 | 350.136 | 0.827717, FAIL | 481.754 | 1.13886, PASS |
| B3 | 506.293 | 503.985 | 0.995441, FAIL | 521.441 | 1.02992, PASS |

The +1 ms table uses one complete sample per method. Non-regression means a measured baseline/R3 ratio ≥1; it does not establish timing variability or statistical significance. A small single-sample difference must not be described as a stable speed advantage or disadvantage.

Per-repeat comparisons against R2 and fresh Stage4L appear in MOTION_SOLVER_R3_GATE.json. Only complete equal-duration runs with verified identical initial arrays enter a speedup. Incomplete R2 B2 +5 ms is never extrapolated. Three repeats are required for the final performance gate; 5× is the target, not a promised result.

The initial 50 ns layout benchmark did not represent B3's full-task cost. Complete 5 µs full/half1/half2 trials at +1.5/+2.5 ms gave global medians 1.618/1.321 s versus serial-component 0.351/0.288 s, with fixed-ID differences ≤2.43e−17 m. Persistent threaded component layouts were slower. The actual full tasks determine the final speedup; layout microbenchmarks do not.

- B2 float64 GPU pair-force inclusive speedup: 2.855×. Transfers and synchronization are included; this is not trajectory speedup.
- B3 float64 GPU pair-force inclusive speedup: 2.972×. Transfers and synchronization are included; this is not trajectory speedup.
- Independent two-job throughput speedup: 1.589×, separately labeled.
- CPU affinity/cgroup, per-process CPU time, actual worker threads, memory and GPU transfer/kernel/synchronization records are retained. Device/host-wide utilization is not evidence of solver parallelism.
- Formal scientific jobs run sequentially with fixed case-specific layout, worker and BLAS limits. Some jobs overlap disclosed read-only transfers of closed artifacts (nice=19, idle I/O, ≤512 KiB/s); this is not perfect I/O isolation. Transfer intervals are retained in diagnosis/*SYNC*.json, and no time is subtracted from measured completion wall. Host scheduling and transfer variability remain timing limitations.
- G7 and DYNAMIC_COST_PROFILE.json partition actual force/geometry/assembly/QP/nonlinear/full-half work without summing overlapping timers. Instrumented diagnostics are separate from formal completion timing.

## Tests, budgets and provenance

Source-bound formal tests: **PASS**, counts `{'tests': 124, 'failures': 0, 'errors': 0, 'skipped': 0}`. Six new numerical-state holdouts: **PASS**. Later reporting/sampling tests have a separate JUnit file; they do not relabel the frozen test executable.

Development conservative budget charge: 7745.227 s / 14400 s. Formal carried charges: `{'new': 9936.771708892658, 'reference': 18984.07657995168, 'r2': 2882.373643346131}`. Overall charged: 39548.449 s / 63000 s. Development allowances and the unmetered exploration reserve are budget accounting, not measured performance. No reset follows a failed candidate.

The original strict bitwise restart failure (one coordinate ULP, 2.1684e−19 m) is retained alongside the separate resolution-based PASS and exact cache serialization PASS. Missing-fixture/import failed test invocations remain recorded. Two early local source ZIPs captured a later deployed driver after workers had loaded V1; original snapshots plus LOCAL_V1_PROVENANCE_CORRECTION.json and a reconstructed loaded-source archive remain. Later formal execution uses immutable trees.

## Deliverables and remaining blockers

- Required audit JSON and gate: `results/stage5/motion_solver_r3/final/`.
- G1–G8 PNG/PDF, underlying actual data and source hashes: `figures/` and `figures/data/`.
- Actual-radius event video and complete old/new videos, decoded timestamp/frame audits: `videos/`.
- Full-system fixtures, candidate arrays, real WALL and complete failed-run ledgers: `review/`; standalone `Stage5A_R3_review.zip`.
- Complete raw/server-only file paths and hashes: `final/ARTIFACT_MANIFEST.json`; exact execution/archive distinctions: `final/DELIVERY_IDENTITY.json`.
- Reproduction commands: `scripts/stage5a_r3/README.md`; isolated dependency lock and license inventory retained.

Remaining blockers are global all-active trajectory accuracy/reference disagreement and any failed measured performance gate. Numerical safety, global mass conservation, GPU agreement or a smooth animation does not establish local physical accuracy. No production/reference promotion and no longer-duration run is authorized by these partial successes.
