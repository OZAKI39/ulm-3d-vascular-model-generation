# Stage5A motion solver development

**FAILED_WITH_EVIDENCE — MOTION_SOLVER_NUMERICAL_GATE_FAILED**

Branch `dev/stage5a-motion-constrained-20261006`, base `40765a4f6252506b9ffdfb35453f62ca537c0797`. Preregistered plan commit `b24fbc0fd`; original active-set implementation `ce1808fd1`; OSQP scientific implementation `7dbc71f19f83687a8ff7e2646bc8402f1d5f634a`. All source/config identities are retained per run. The historical interrupted campaign is preserved separately and was not resumed.

## Bottleneck and mathematical formulation

Historical B2/B3 5 ms half-gap limiting fractions were 0.999425/0.999463. Contact consumed 1092.34/1825.85 s; sparse resistance assembly+solve consumed 804.50/899.01 s. Each old RK2 step uses two mobility/evaluator calls, four contact calls and two projections. Historical L-BFGS iteration counts were unavailable; nine independent saved-state replays per case measured them without overwriting old evidence.

The experimental minimization is 0.5 Uᵀ A U − bᵀ U subject to G U ≥ lower, with unchanged Stage4L A=D+Lpp+Lpw and finite-horizon normal bounds. Pair rows contain at most six nonzeros; wall rows at most three. Vectorized normal outer-product assembly and sparse CSC KKT systems avoid full dense Jacobians/dual matrices. Following explicit user authorization, OSQP 1.1.3 was installed in an isolated dependency overlay. The original NumPy/SciPy/CuPy environment is unchanged, verified by before/after package inventories. OSQP uses builtin CPU sparse algebra; one sparse saddle polish and one cold retry can recover an inaccurate warm candidate only if all original inequality/KKT audits pass. The earlier matrix-free active-set implementation remains separately selectable and its failures remain preserved.

KKT diagnostics use velocities divided by the maximum absolute target velocity and resistance divided by median Stokes drag. They audit scaled stationarity, inequality feasibility, nonnegative multipliers, active slack and complementarity. The separate resistance residual is relative to the physical force right-hand side. Summary KKT/linear maxima include completed solves in rejected local-error trials as well as accepted states; failed solves remain separate failure records.

## Accepted-state reuse, components and integration

Canonical particle/contact IDs carry primal velocity, dimensional dual multipliers, active contacts and sparsity identity only after acceptance. OSQP CSC workspaces are updated when the sparsity pattern matches; every attempted solve explicitly starts from accepted cache values or a cold zero state, not an unaccepted prior trial. Component decomposition includes both lubrication and horizon constraints; walls stay in their owning component. It is exact only for the current local resistance model. It does not generalize to full many-body mobility.

One constrained predictor and one corrected solve refresh forces, lubrication coefficients and exact wall gaps/normals. Both start and predictor tangent inequalities are retained in the corrected displacement problem. Position/gap errors cause binary subdivision of that physical interval only; future nominal intervals reset. Final exact wall/pair checks and all-pair swept-segment checks guard acceptance. The OSQP attempts accepted no new state: zero projection calls are not successful motion-safety evidence. No force caps, damping, contact-tolerance relaxation, ablation or coordinate smoothing was introduced.

## CPU/GPU and tests

123 focused tests passed, including T1–T15, 25 OSQP tests and legacy Stage4L/4M regression; 1 historical active-set real-QP regression remains xfail. That same saved QP passes with OSQP. New dynamic trajectory QPs nevertheless failed, so this test result does not promote the solver. Frozen legacy class-source hashes and force/ablation bytes are unchanged. Archived primary nominal/dt-half transport comparisons were recomputed with input SHA checks; that test does not validate new long transport.

Benchmark: `{"qp": {"components_cold": 0.4306160664418712, "components_parallel": 0.024698270484805107, "components_warm": 0.024696236941963434, "global_cold": 0.49164859391748905, "global_warm": 0.02280270808842033}, "pair": {"auto": 0.0015283729881048203, "cupy": 0.0005266780499368906, "numpy": 0.0015664040111005306}}`. CuPy selected for classical pair forces; CPU sparse QP. Serial components selected because the measured parallel overhead was greater. Identical-state warm-start speedup is 17.437x; this is not a claim about continuous-trajectory speedup. Global warm QP was faster than decomposed warm QP in these small cases. Component decomposition remains explicit and equivalent.

## Dense probe outcomes

| Case | Target (ms) | Macro dt (µs) | Status | Reached (µs) | Wall (s) | Accepted | Rejected | Median dt (µs) | Max KKT |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| B2 | 1 | 20 | FAILED | 0 | 8.768 | 0 | 10 | N/A | 8.914097120390094e-12 (completed trial solves only) |
| B3 | 1 | 20 | FAILED | 0 | 10.812 | 0 | 10 | N/A | 8.631508037590868e-11 (completed trial solves only) |

Complete frozen reference wall times: B2/B3 1 ms 400.37/491.86 s; 5 ms 2216.01/3080.30 s. Ratios involving incomplete new runs are not completion speedups. If the joint 1 ms gate fails, new 5 ms results and their speedups are unavailable, and no longer formation or transport run is launched.

Machine-readable scientific differences (Rg endpoint/series, centroid, pair quartiles, rigid fit, neighbor turnover and sample ages) are in `results/stage5/motion_solver/comparisons/`. Partial-prefix comparisons cannot certify a complete 1 ms or 5 ms gate. Native hard-contact counts remain descriptive; no chatter filtering changes trajectories.

## Runtime interpretation and limitations

This is a coefficient-lagged predictor/corrector with implicit constrained velocities, not a fully implicit nonlinear position solve. The continuous dense predictor/corrector problem is harder than a repeated frozen-state QP. Both OSQP runs rejected ten local-error intervals and failed before accepting a step. Nearly parallel start/predictor tangent constraints challenge both ADMM and active-set polishing. An excellent repeated-state warm-start benchmark is not trajectory acceleration. Accepted dt, new-motion Rg/centroid/quantile differences, and completed-run speedups are unavailable; initial-state identity is not a convergence result.

Tolerance values were frozen before scientific runs and were not loosened after poor performance. Runtime exhaustion is incomplete evidence, not solver convergence or completed science. The 1,800 s budget is checked after an accepted step, so one in-progress bounded solver call can cause an overrun; it is not a hard real-time deadline. Classical near-contact secondary-force validity, normal-only immobile-interface lubrication, absent tangential/history/deformation physics and long-transport accuracy remain outside this numerical development.

## Delivered evidence

Eight publication-size PNG/PDF figures are under `results/stage5/motion_solver/figures/`. Gated-off 5 ms comparisons are explicitly labelled; available shorter/partial data are not relabelled as 5 ms. M3 compares measured original fixed-state replay timing with measured new attempted-run timing; it does not invent unrecorded historical stage timing. M8 uses actual last accepted coordinates at the same clock and records sample age, with a common-origin translation and metre-to-micrometre conversion only.

Authoritative gate: `results/stage5/motion_solver/final/MOTION_SOLVER_GATE.json`. Full per-run source/config/input identities, accepted/rejected ledgers, actual trajectories, last accepted state, resource observations and failures remain in their exclusive run directories. New JSON is strict and finite. No original Stage4L/M result was overwritten.

## Attempt cost decomposition

Completed QP calls below include accepted and local-error-rejected trials. Time inside the terminal failed QP is retained in the failure/campaign evidence and not invented as a completed solve timing.

| Method / case / macro dt (µs) | Rejection reasons | Predictor QP (s) | Corrector QP (s) | Corrector polish (s) | Corrector primal fallback QPs |
|---|---|---:|---:|---:|---:|
| Earlier active-set / B2 / 20 | {'LOCAL_ERROR': 71} | 4.182 | 184.492 | 178.149 | 88 |
| Earlier active-set / B2 / 10 | {'LOCAL_ERROR': 70} | 3.897 | 174.025 | 168.178 | 87 |
| OSQP / B2 / 20 | {'LOCAL_ERROR': 10} | 0.838 | 6.063 | 0.007 | 0 |
| OSQP / B3 / 20 | {'LOCAL_ERROR': 10} | 1.992 | 7.761 | 0.009 | 0 |

## Scientific comparisons on the available prefix

The preregistered comparison clock is 10 µs. A prefix with no nonzero common clock is not evaluable. Passing a short prefix does not pass the requested 1 ms gate.

| Case / comparison | Compared time (µs) | Rg endpoint relative | Rg series L2 | Pair Q1/median/Q3 relative | Centroid difference (m) | Prefix assessment |
|---|---:|---:|---:|---|---:|---|
| B2 / scientific_comparison | 0 | None | None | None | None | NOT_EVALUABLE_NO_NONZERO_COMMON_CLOCK |
| B3 / scientific_comparison | 0 | None | None | None | None | NOT_EVALUABLE_NO_NONZERO_COMMON_CLOCK |

## Observed numerical failure and blockers

Bounded B2 single-interval reconstruction: **REPLAY_REJECTED**. The earlier active-set cache is reconstructed from one preceding accepted interval. No new trajectory is continued. Diagnostic: `results/stage5/motion_solver/diagnostics/20261006T200738743816Z_failure_replay/DIAGNOSTIC.json`.
Failed component: 114 primal unknowns, 192 inequalities; active saddle condition 5.2167e+11. 10 rows have a nearly parallel partner. The feasible fallback exhausted 3940 iterations. Its final working multiplier minimum is -65.56658161459079; this does not satisfy dual feasibility.
The original failure-summary stationarity value was recomputed with zero multipliers after fallback exhaustion; it is not an accepted-state residual. The saved matrix audit additionally reports actual last-working-set stationarity, negative multipliers and active slack. Dense SVD was used only on this diagnostic component, never in production motion.
Bounded B2 single-interval reconstruction: **REPLAY_REJECTED**. The OSQP replay starts at the unchanged initial state with an empty accepted cache. No new trajectory is continued. Diagnostic: `results/stage5/motion_solver/diagnostics/20261006T210101208569Z_B2_osqp/DIAGNOSTIC.json`.
OSQP diagnostic constraint matrix dimensions: [102, 51]; numerical rank 49; 102 rows have a nearly parallel partner. Dense rank/cosine checks are diagnostic-only; all production QP matrices remain sparse.
Bounded B3 single-interval reconstruction: **REPLAY_REJECTED**. The OSQP replay starts at the unchanged initial state with an empty accepted cache. No new trajectory is continued. Diagnostic: `results/stage5/motion_solver/diagnostics/20261006T210102848863Z_B3_osqp/DIAGNOSTIC.json`.
OSQP diagnostic constraint matrix dimensions: [22, 18]; numerical rank 17; 10 rows have a nearly parallel partner. Dense rank/cosine checks are diagnostic-only; all production QP matrices remain sparse.

Blocking issues: robust treatment of nearly dependent start/predictor constraints, continuous-state QP convergence and cost, and stiff integration-error control. Passing isolated or frozen-state tests does not resolve these. No tolerance, activation/floor, secondary force, acoustic field or CFD value was altered to make the gate pass.

## Accepted-state safety and resources

B2, target 1 ms, macro 20 µs: accepted pair/wall gaps are unavailable (zero accepted steps); projection calls 0. Maximum linear relative residual across completed trial solves: 2.126508931901684e-14. The unchanged geometric tolerance is 1e-11 m. Peak RAM 0.865 GiB, GPU allocation 412.0 MiB, mean process CPU 100.1% (100%=one core), sampled GPU utilization 0.27%.
B3, target 1 ms, macro 20 µs: accepted pair/wall gaps are unavailable (zero accepted steps); projection calls 0. Maximum linear relative residual across completed trial solves: 3.6921781615190574e-14. The unchanged geometric tolerance is 1e-11 m. Peak RAM 1.002 GiB, GPU allocation 414.0 MiB, mean process CPU 100.1% (100%=one core), sampled GPU utilization 0.00%.

The earlier FORMAT_AUDIT.json applies to the original active-set revision only. OSQP run identities record source commit 7dbc71f19 and full file hashes. Later analysis scripts are separately hashed and do not rewrite executed provenance. The frozen legacy class/force hashes continue to pass.

## OSQP installation and dynamic failure details

Installation lock: requirements-stage5a-osqp.lock.txt. Invocation: `bash scripts/stage5a/python_osqp.sh ...`. Installation, package hashes, smoke tests, failed diagnostic attempts and combined JUnit report are retained in `results/stage5/motion_solver/environment/osqp_20261006/`. No scientific tolerance was changed after the OSQP trajectory attempts.

B2: failed QP has 51 unknowns and 102 constraints; OSQP returned `solved inaccurate` after 100000 iterations, with independent candidate KKT 1.07494e-09. Sparse polishing did not produce a valid full KKT certificate. This is a rejected trial, not an accepted-state residual.
B3: failed QP has 18 unknowns and 22 constraints; OSQP returned `solved inaccurate` after 100000 iterations, with independent candidate KKT 7.84371e-10. Sparse polishing did not produce a valid full KKT certificate. This is a rejected trial, not an accepted-state residual.

B3's ADMM candidate KKT is below 1e-9, but its `solved inaccurate` status and unsuccessful sparse polish do not satisfy the preregistered adapter acceptance rule. B2 also exceeds the independent KKT bound. Neither failed candidate was used to advance physical time. Installing OSQP resolves the dependency gap and a saved QP failure; it does not by itself resolve dynamic near-dependent constraint handling.
