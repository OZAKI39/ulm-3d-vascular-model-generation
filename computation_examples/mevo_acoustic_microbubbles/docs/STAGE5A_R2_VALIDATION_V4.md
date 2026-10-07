# V4: stricter parameters, identical V6 numerical kernel

V3's ordinary B2 1 ms configuration fails the registered old-reference Q1
threshold: 2.20205% versus 2%. Its already completed refined configuration has
1.75363% and passes. The same-ID self comparison passes, but this does not erase
the ordinary configuration's `BASELINE_DISAGREEMENT`. V3 remains immutable.

One additional **parameter refinement** is registered here. This is not a seventh
integrator revision: every critical source hash remains exactly V3's. Selected
maximum macro dt is 5 us and all original time-error budgets are multiplied by
0.25. Independent refinement uses 2.5 us and factor 0.0625. The original KKT,
resistance, geometry and scientific thresholds are unchanged. No additional
parameter/algorithm tuning is permitted after this V4 validation.

The unchanged V6 recovery function caps S0 macro steps at 50 ns for either
factor below one. Thus this cap is **not** halved again; time-error budgets are
quarter-sized and the actual S0 grid must be reported if recovery occurs. S0
is still labelled start-matrix equilibrium, never endpoint implicit convergence.

New sealed states: floor(0.07/0.31/0.66 times (Nframes-1)) from both immutable
Stage4L 5 ms files. These fractions are distinct from all earlier development
and V1/V2/V3 holdout selections. Freeze source/parameters and run fresh JUnit
before opening the six states at both configurations (12 motion probes). Earlier
failures are retained; no fresh holdout may be used to tune V4.

V3's completed refined 1 ms trajectories are identical to V4's selected
configuration and may be reused, explicitly retaining old run IDs, executed
source identities, measured costs and `original_kind=refined`. Frozen reference
1 ms runs are also reused. Only the extra refinement is newly run. Reuse is
recorded in V4's campaign ledger, never represented as another measurement.

All V1/V3 costs carry forward to the original 10800 s new-method and 10800 s
reference allocations. Per-run 1800 s and all step/retry/resource limits remain
unchanged. The identical-kernel check in the existing run driver still verifies
V3 source hashes; the V4 controller additionally verifies V4 configuration,
fresh test gate and orchestration hashes before each execution. Each actual
driver invocation records its full numerical configuration and executed source.

After both 1 ms cases pass self-refinement and old-reference comparisons, run
5 ms and necessary refinement. A failure or budget boundary is reported without
further retuning. Performance uses only complete equal-duration measured runs.
