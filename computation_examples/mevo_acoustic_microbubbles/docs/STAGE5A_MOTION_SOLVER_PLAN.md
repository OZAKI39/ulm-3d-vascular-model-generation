# Stage 5A motion solver preregistration

Status: FROZEN_BEFORE_NEW_SCIENTIFIC_RUNS. Base 40765a4f6252506b9ffdfb35453f62ca537c0797; branch dev/stage5a-motion-constrained-20261006.

## Scope and source audit

Only translation numerics change. Stage4L D + Lpp + Lpw, coefficients, 0.1R activation, 0.001R floor, contact tolerance, acoustic/response/flow and A0 physics remain frozen. No BraVa particle physics is imported. The six requested BraVa files were inspected for interval subdivision, canonical IDs, sparse assembly and acceptance-ledger ideas. Their dense 6N and RBC implementations are not copied.

Read: tracker, translation_backends, contact, forces, secondary_ablation; Stage4L dense/branch/campaign/analysis/native-graph/figure/replay/analytic tooling and rendering entry-point inventory; all three requested tests; model, Stage4L/4M development reports and final gates. No project AGENTS.md exists in the selected tree or its local parents. Server /etc/vast-agents-guide.md applies.

## Measured baseline

Call chain: dense_probe -> GlobalSynchronousIntegrator -> inherited ParticleTracker.evaluate/step; evaluate queries acoustic/response/flow/wall, solves mobility and instantaneous contact. RK2 step applies two horizon contacts and two projections, plus midpoint evaluate; driver evaluates accepted endpoint. Thus steady per-step counts are 2 evaluations/mobility, 4 contact calls and 2 projections; initial evaluation adds one. Counts below are measured from historical logs. No historical L-BFGS iteration count was recorded.

| Probe | Steps | Evaluate/mobility calls | Contact calls | Projection calls | Half-gap limited fraction | Contact seconds | Assembly+solve seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| B2_1ms | 3626 | not recorded | not recorded | 7252 | 0.9994484280198566 | not recorded | not recorded |
| B2_5ms | 17385 | 34771 | 69541 | 34770 | 0.999424791486914 | 1092.3386307933833 | 804.4963788329624 |
| B3_1ms | 3684 | not recorded | not recorded | 7368 | 0.999457111834962 | not recorded | not recorded |
| B3_5ms | 18627 | 37255 | 74509 | 37254 | 0.9994631448971922 | 1825.852679145988 | 899.009501129156 |

Historical dense matrices have 3N=558/567 primal columns and m contact rows (each sparse pair row would need only six entries). Exact constraint distributions and QR refinement counts: `results/stage5/motion_solver/audit_v2/HISTORICAL_AUDIT.json`. 1 ms separate translation logs are absent; their missing totals are null, never fabricated from 5 ms. The first audit failure for that absent log is preserved.

Nine independent original accepted-state replays per B2/B3 (initial through final state, not a new concatenated trajectory) measured response, pair force, wall queries, mobility, contact, projection, L-BFGS iterations and lubrication/horizon component distributions. `FIXED_STATE_PROFILES.json` is authoritative. Timings of nested methods overlap and must not be summed as exclusive totals. These measurements use the final frozen implementation; older historical dense runs predate the final QR fallback repair.

B2 replay totals (9 independently initialized steps): `{"wall_distance": 0.007641170173883438, "wall_gap": 0.0010804911144077778, "mobility": 0.5475127012468874, "contact": 0.5891020956914872, "evaluate": 0.9060647578444332, "wall_project": 0.01310310885310173, "projection": 0.014771735994145274, "response_queries_s": 0.009117994923144579, "pair_force_s": 0.042934359051287174, "flow_interpolation_s": 0.011040763463824987}`.
L-BFGS iterations: `[5, 5, 5, 5, 4, 9, 23, 21, 11, 11, 11, 25, 21, 11, 11, 10, 18, 20, 10, 10, 11, 21, 20, 10, 10, 11, 20, 19, 11, 10, 9, 22, 19, 9, 9, 11, 20, 20, 12, 12, 10, 11, 13, 10, 11]`; initial component sizes: `[75, 17, 5, 12, 38, 11, 3, 1, 8, 9, 1, 1, 3, 1, 1]`.

B3 replay totals (9 independently initialized steps): `{"wall_distance": 0.007058890303596854, "wall_gap": 0.0010741699952632189, "mobility": 0.5533664408139884, "contact": 0.8329040214885026, "evaluate": 0.9510557006578892, "wall_project": 0.013274477794766426, "projection": 0.015486129792407155, "response_queries_s": 0.009161189198493958, "pair_force_s": 0.04524556850083172, "flow_interpolation_s": 0.00991054461337626}`.
L-BFGS iterations: `[5, 4, 4, 5, 4, 11, 24, 23, 11, 11, 11, 27, 23, 10, 11, 11, 44, 30, 11, 11, 11, 28, 31, 10, 11, 13, 25, 26, 11, 12, 11, 24, 21, 11, 12, 12, 28, 24, 12, 12, 13, 17, 17, 14, 14]`; initial component sizes: `[12, 1, 32, 91, 4, 11, 4, 5, 1, 17, 6, 5]`.

## Solver and integration design

Opt-in `constrained_lubrication_v1`; old `lubrication_coupled` calculations remain unchanged. Existing environment contains SciPy 1.15.3 but neither OSQP nor ProxQP. Use the explicitly authorized deterministic sparse active-set fallback: CSC resistance and constraint matrices, sparse saddle-point polishing, matrix-free bounded dual initialization using existing SciPy when required. No production dense G or dense G A^-1 G.T. Small dense references remain test-only. Record full stationarity, primal, dual and complementarity checks. Never drop an original inequality to resolve redundancy.

Accepted-state cache contains canonical particle/contact IDs, primal and dual states, active set and CSC pattern. Failed attempts cannot commit this cache. Components use the exact union of lubrication edges and all horizon constraints. Walls stay within their particle component; isolated particles with wall coupling cannot bypass that coupling. Compare global/components and cold/warm. Benchmark parallel components before enabling them.

Predictor: constrained velocity at x_n. Corrector: refresh forces, coefficients and exact wall geometry at x_p; solve for averaged motion in predictor A metric, with start and predictor tangent constraints mapped to displacement from x_n. One predictor and one corrector; no Newton iteration. Final exact geometry and swept pair checks are mandatory. Residual projection is allowed only below the existing contact tolerance, never a normal advancement mechanism. Reject and bisect only the failed physical interval. Every right child begins at its accepted left endpoint; next nominal interval resets.

## Frozen numerical settings and gates

- Nominal macro interval 20 microseconds; companion dt/2 = 10 microseconds. All runs start from the exact immutable original dense checkpoint. No longer formation run is part of the mandatory campaign.
- Position error: maximum predictor/corrector displacement <=0.002 R + 0.01 times maximum proposed displacement; gap discrepancy <=0.002 R + 0.02 times max(original positive gap, 0.001R) for each local interaction. These are integration-error settings, not changed force/contact parameters. Maximum subdivision depth 24; finite 1,800 s/case runtime budget and 100,000 accepted-step guard.
- Contact normalized KKT <=1e-9; resistance relative residual <=1e-11, finite states and no failed solver state. All inequalities, symmetry and positive self-drag plus nonnegative lubrication coefficient certificate checked. No unresolved penetration beyond the existing config tolerance; no multi-radius correction.
- Same-ID largest-cluster Rg endpoint and relative time-series L2 <=2%; endpoint pair-distance Q1/median/Q3 relative differences <=2%; centroid endpoint difference <=max(0.1R,2% of baseline centroid distance travelled). Other envelope, rigid-fit, neighbor/hard-contact metrics remain descriptive. Hard-contact chatter alone cannot fail the scientific gate.
- Compare on 10 microsecond physical clocks using actual last accepted states and sample age, never interpolated coordinates. Same thresholds apply to new nominal versus dt/2.
- T1–T15, sparse/dense, component/global, warm/cold equivalence must pass. Freeze a source hash regression for legacy class and force/ablation files.
- Both B2 and B3 1 ms numerical/scientific gates must pass before either 5 ms run. Both 5 ms cases must meet numerical/scientific gates AND >=5x total wall speedup versus the immutable Stage4L historical baselines to promote. 10x is a stretch target. Failures close the downstream gate; no tolerance tuning after scientific results.

Existing B2/B3 primary nominal/dt-half evidence is used for transport regression and tolerance context: all IDs agree; crossing maxima 31.10/12.47 ms; early nearest-neighbor relative L2 1.21e-5/5.47e-6. This does not establish dense-secondary pointwise convergence, and therefore does not justify loosening the proposed 2% dense gate. Recompute the available archived identity/outcome comparisons read-only; no new 6 s transport campaign.

## Resources, provenance and deliverables

Fresh server: 8 physical / 16 logical cores, cgroup CPU quota 7.68, RAM limit 45156925440 bytes; RTX4090 24564 MiB. NumPy2.2.6/SciPy1.15.3/CuPy14.2.0/CUDA runtime13020. No installation or system modification. BLAS/OpenMP threads=1. CPU sparse QP; benchmark NumPy/CuPy/auto pair forces with synchronized host outputs, then choose by measured cost. At most one GPU scientific job. No old controller resumed.

All new runs have exclusive directories, strict finite JSON, input/source/config hashes, last accepted states and rejection/failure diagnostics. Raw source paths and SHA checks are in BASELINE_SELECTION/HISTORICAL_AUDIT; absent prior log digests are explicitly distinguished from verified sealed trajectory hashes. Runtime breakdown separates recorded history, replay measurements and new-run exclusive timings.

Commit core code, T1–T15 evidence, benchmark raw rows/summary, eight PNG/PDF figures, development report and final gate in this independent branch. Missing gated scientific results stay explicitly unavailable in figures, never zero or invented curves. Final status: PROMOTED_WITH_SCOPE, RUNTIME_GATE_FAILED or NUMERICAL_GATE_FAILED, with FAILED_WITH_EVIDENCE on failure.
