# Stage 5A-R3 preregistered development plan

Base: the actual R2 archive commit `07db8ed35f3a460a4ba315e42b15d257f8518b0a`.
R2 execution identities `7f83f802c6cc435918ea25d69dcbb0f4f4dacce9` and reused
`5198d3c705c6a3735dcf566f53810274c40a6624` remain distinct. R3 uses a new
worktree, branch, server directory and virtual environment. All previous results
are read-only inputs; no previous campaign is resumed.

## Order and fixed physics

1. Verify baseline/source/input hashes and current resource limits. Read complete
   R2 integration, QP, geometry, force interface, driver and timing implementations.
2. Reproduce the two real B2 terminal events locally in physical time, using every
   active particle and full forcing. Recover nearby actual states and explicitly
   distinguish saved controller/cache from reconstructible or missing metadata.
3. Audit reference lineage and all active same-ID trajectories before selecting
   any new reference. Compare common-state equations before trajectory differences.
4. Test continuous geometry classification, exact/higher precision boundary
   predicates and bounded event handling. No sampled-path substitute for a sweep.
5. Independently refine S0 and S1, protect B3, and profile four dynamic windows.
   Benchmark global/component serial/persistent parallel and CPU/GPU end-to-end.
6. Freeze candidate, accuracy budgets, new holdouts, resources and source-bound
   tests. Run original 3 s starts for +1 ms before +5 ms and refinements.
7. Complete references, three equal-duration timing repeats per case/method where
   preceding gates permit, actual-motion figures/video and an audited Git delivery.

All R2 frozen acoustic, APECSS, flow, A0, radius, particle, mesh, drag and normal
lubrication definitions remain unchanged. Contact tolerance is 1e-11 m, KKT is
1e-9, resistance residual is 1e-11, projection remains zero. Altering a mathematical
geometry/force convention requires MODEL_CHANGE_REQUIRED rather than a numerical
success label. Merely tightening a certificate without a justified bound is not a
repair. Exact predicates apply to the adopted straight physical segment, not an
auxiliary/sampled curve. Rejected candidates never enter accepted history.

## Finite resource and revision budgets

Measured current instance: CPU quota 7.68 cores, affinity 0..15, cgroup memory
42 GiB, RTX 4090 24 GiB. Reserve capacity for existing services; diagnostic
concurrency <=6 CPU workers, each BLAS/OpenMP=1, aggregate RAM <=20 GiB and GPU
<=4 GiB. Formal trajectories run exclusively one at a time. Persistent component
workers compete under these same limits; selection follows actual measurements.

Development budget: at most eight documented numerical revisions and 4 h aggregate
experiment process wall time (parallel tasks charged individually). Local probes
default <=300 s, expanded targeted probes <=900 s, no new long formation runs.
Formal new-method budget <=4 h total, <=2400 s per +5 ms run and <=1200 s per
+1 ms run. Frozen Stage4L complete +5 ms reference budget <=4500 s per run,
<=8 h total; this is chosen before execution from R2's observed 1800 s partial
durations (B2 3.926 ms / B3 2.895 ms), with finite allowance for later stiffness.
R2 protected B3 replays/timing <=1.5 h total. Overall formal + development
allocation <=17.5 h; no budget reset after candidate failures. Calls also have
iteration limits, SIGALRM and supervisor timeout/kill bounds. Hitting a budget
preserves evidence and closes downstream gates; it does not justify relaxation.

## Accuracy protocol before looking at R3 outcomes

All-active errors retain original ID mapping, with raw, translation-removed and
rigid-fit residuals (RMS/p50/p95/p99/max), original 75/91-ID groups and complements.
The provisional target resolution is 0.01R = 17.8 nm in RMS and 0.05R = 89 nm
maximum same-ID error against an independently refined numerical solution. These
are resolution targets, not thresholds inferred from observed discrepancies.
They must be justified using development convergence and sealed before holdout;
without an independent error bound, all-active accuracy is NOT_EVALUABLE.
Unexplained O(R) legacy/new discrepancies prevent reference promotion even when
aggregate cluster metrics pass. Same-clock sampling bounds are included; no
coordinate interpolation or ID rematching. Original cluster thresholds remain.
S0 uses separately halved window ceilings at three development precision levels.

Every hypothesis is CONFIRMED / REFUTED / UNRESOLVED with concrete data. Every
gate identifies case, physical duration, method and repeat; incomplete results
never enter completion speedup. Old archived and fresh references remain separate.
Runtime prediction does not count as a completed measurement. Performance goals
are >=5x for both complete +5 ms cases, not an assertion of achievable speed.

## Deliverables

Versioned R3 backend/certificates and scripts; complete provenance and lightweight
local fixtures/logs; strict JSON audits named by the task; source/environment-bound
JUnit/holdout; G1..G8 PNG/PDF and actual-state video; reproducible commands and
manifest with existing server paths/hashes for large assets. New branch only.
