# R2 development ledger

Base 7c39828d9; plan 5877e733a. Source before formal validation is development,
not promoted. All frozen force and lubrication classes remain unchanged.

1. QP strategy 1: native ProxQP 0.7.2 at eps_abs=1e-12. All three original
   failed QPs pass unchanged, KKT 1.5e-14 / 3.5e-14 / 9.1e-13 (dense).
   Sparse also passes. Native comparison includes OSQP and Clarabel failures.
   34 new manufactured/interface/real-QP tests passed (tests_dev01.xml).

2. Integrator V1, source 18c7fa38f: converged endpoint resistance AND endpoint
   contact tangent. B2 prefix exposed SWEEP_NOT_CERTIFIED: initial negative gaps
   close to -1e-11 m can pass endpoint conditions but fail intermediate straight
   segment clearance. E.g. dt=1e-8 s minimum swept gap -1.2411e-11 m.
   Steps shrank below 1 ns, not acceptable progress. This diagnostic was stopped
   by the developer at a checkpoint (not resumed old campaign, not a user stop).
   Preserve run `20261006T215835450341Z_B2_S1_1e-06s_prefix_v1`.

3. Integrator V2 hypothesis: a SINGLE start-tangent normal set ensures a convex
   pair distance remains above its initial tangent lower bound along a straight
   segment; retain endpoint self-consistency only in A(x1). Equation:
   A(x1)U = F(x0,t0)+D*u_flow(x0,t0)+J(x0)^T lambda, x1=x0+dt U.
   This is first-order semi-implicit contact/forcing, nonlinear implicit normal
   lubrication resistance. End wall geometry still needs full WALL sweep.
   No new force term, projection, gap floor or tolerance. New source/tests/run
   identity required; V1 evidence is retained in its source snapshot.

4. V2 completed B2 1 us with 60 accepted half steps, 15 nonlinear rejections,
   33.43 s, KKT 9.04e-10, linear residual 3.59e-14, minimum swept pair gap
   -9.99757e-12 m and no projection. Its 10 us continuation was stopped for
   diagnosis; repeated Picard iterations were the bottleneck (1888 QP calls for
   the 1 us prefix). The solver-map velocity difference hovered at 1e-9 while
   independent nonlinear equilibrium can already be substantially better.

5. V3: bounded 6-history Anderson acceleration of the nonlinear solve, with a
   trial-iterate step safeguard; acceptance remains unregularized KKT <=1e-9.
   Velocity-map delta <=1e-7 triggers evaluation of the final original-equation
   residual, it does NOT itself accept a step. Fix the candidate contact set for
   each substep. 34 regressions including whole-step squeeze/wall refinement pass.

6. V3 global sparse execution and CuPy pair force completed both 100 us prefixes.
   Initial-state layout diagnostics measured global 5.07 ms, component serial
   11.92 ms and component parallel 16.67 ms; CuPy pair 0.517 ms vs NumPy 1.553 ms.
   These are development diagnostics, not trajectory completion speedups. Initial
   layout results later proved insufficient to choose a layout for the whole run.

7. V4 retains the V2 discrete equations and V3 bounded Anderson acceleration. It
   adds the positive normal gap-derivative part to the *nonlinear solver metric*
   with the matching right-hand-side correction. The returned step is still
   audited against the original, unmodified endpoint resistance equation. There
   is no added physical force, energy derivative, projection or artificial drag.
   Trial primal/dual caches now stay separate from accepted caches. B2/B3 completed
   1/10/100 us checkpoint continuations. B3's 10->100 us segment used 408 accepted
   half steps and 32 LOCAL_ERROR rejections, 147.52 s, no nonlinear rejections.

8. Manufactured_v2 solves stiff pair squeeze and wall approach for 1 ms against
   an independent scalar ODE. Actual max-dt refinement gives 20/40/80 accepted
   half steps and first-order error reduction. Merely tightening error tolerances
   left the wall case's first two grids identical; that is recorded separately,
   not relabelled as a dt/2 test. A raw 1 ms step can be safe yet inaccurate.

9. Validation V1 was frozen at executed source a22b208f9. Fresh JUnit: 166 PASS,
   one retained old xfail. All three original failed QPs pass without modifying
   their matrices. Six sealed holdout QPs pass global/component and warm/cold
   equivalence; five holdout motion probes complete 1 us. B2 old frame 17384
   stalls at 0.21103515625 us with NONLINEAR_NOT_CONVERGED and minimum dt.
   This is an unresolved robustness failure, not a PASS. Its data are not used
   to tune the subsequent layout performance benchmark.

10. V1 global sparse B2 completes 1 ms from the original 3 s state: 1036 accepted
    half steps, 94 rejections (90 time error, four nonlinear), 559.59 s including
    retries and initialization. Accepted KKT 2.71e-10, nonlinear KKT 9.976e-10,
    resistance residual 3.04e-14; pair/WALL swept bounds stay within 1e-11 m,
    projection is zero. Scientific 1 ms comparison was not yet available at this
    point. About 466.39 s was spent in QP calls: stable late steps remained costly
    despite reaching the 10 us macro cap.

11. The developer safely interrupted the subsequently launched V1 B3 nominal at
    7.2415 us (checkpoint and partial HDF5 retained) to investigate the measured
    late-state cost before consuming the formal budget. The V1 campaign's own
    RESULT.json records this disposition; it is not a scientific failure or a
    complete B3 run. No old/user-stopped campaign was resumed. New V2 holdouts at
    21/59/89% of each frozen 5 ms trajectory were sealed *before* layout tuning.

12. QP strategy 2 benchmark: actual B2 R2 1 ms development states and B3 R2 100 us
    development states, 25/50/75/99% elapsed, two complete 5 us implicit substeps
    per layout. Global sparse, component serial and bounded global dense all pass
    (48 substeps). Maximum cross-layout position difference 1.334e-16 m; relative
    velocity difference 3.41e-9. With six independent CPU workers/BLAS=1, median
    substep times for sparse/component/dense global are B2 0.438/0.135/0.335 s and
    B3 0.518/0.371/0.708 s. This selects component serial for V2 formal validation;
    it is a throughput diagnostic, not an isolated complete-run speedup. The
    bounded dense alternative (<=600 unknowns, <=1200 rows, <=8 MiB dense input)
    remains opt-in and is not selected. No integrator equation/tolerance changes.

13. V2 strengthens infeasibility handling: statuses are checked against known
    feasible witnesses on every original row; unresolved infeasibility stops for
    diagnosis. Three repeats of an identical failed original QP stop ineffective
    interval shrinking. Original matrices and all native failed candidates are
    retained. Native nonconvex status is fatal. V2 receives a fresh source-bound
    JUnit/protocol and new holdout; V1 records remain immutable.

Timing audit: `timings_exclusive_s.error` in early driver summaries includes wall
queries also counted by `wall`. The final mutually exclusive plot uses only
forcing, assembly, QP, I/O and the measured remainder. It never stacks the
overlapping wall/error counters or adds inclusive monitor time to 100% costs.
The label correction is explicit; raw historical counters are retained.
