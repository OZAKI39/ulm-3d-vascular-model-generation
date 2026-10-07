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

14. V2 component execution was frozen at bac62b1d1; 169 tests pass and the same
    old xfail remains. A probe-script configuration error first evaluated the new
    sealed states with global motion again. Its `HOLDOUT_V2.json` is retained as
    a global-layout diagnostic. The corrected `HOLDOUT_V2_COMPONENT.json` uses
    component serial and yields the same result: all six QPs pass, five 1 us
    motion probes pass, and B2 frame 15472 stalls at 0.8810724 us. The earlier B2
    frame 17384 also still stalls. No V2 formal campaign was launched. New V3
    holdouts (13/46/83% of both old trajectories) were sealed before further
    algorithm work; earlier failed states are now explicitly development inputs.

15. Direct failed-state replay confirms a non-smooth numerical geometry trigger:
    particle 152's closest-WALL normal changes by 0.1343 while trial positions
    differ by about 3.9e-19 m. Its nonzero wall lubrication coefficient is
    3.3247e-5 N s/m and the resistance matrix changes by about 1.2%. At 1e-12 s,
    shrinking dt can therefore leave the fixed-point map oscillating. This is
    measured discontinuity, not a proof that the original QP is infeasible or a
    proof that the nonlinear equation has no solution. See `nonlinear_v4/`.

16. Integrator V5 (not selected) updated only wall contact tangents at the endpoint
    while retaining start pair tangents and full endpoint resistance. Six true
    trial probes on the two stalled states fail to converge. The opt-in experiment
    and failures are retained. No normal smoothing or geometry changes were made.

17. Integrator V6 is a **hybrid first-order method**, not an assertion that every
    step is implicit. Its main S1 method is V4. After a failed nonlinear iteration
    with a >0.01 jump of a WALL normal inside the lubrication activation band,
    it tries a burst of at most eight S0 macro steps capped at 100 ns (50 ns in
    refinement), then returns to S1. At most 1024 bursts and all original run,
    accepted/rejected-count and time budgets apply. S0 solves the original
    **start-state** QP; its endpoint nonlinear certificate is NOT_APPLICABLE,
    never fabricated as zero or a converged S1 result. Every S0 trial still
    undergoes independent full-vs-two-half error estimation and continuous pair/
    real WALL sweeps. Failed S1 work remains charged and recorded. No force,
    tolerance, floor, radius, position projection or particle identity changes.
    The numerical recovery state is included in the atomic controller checkpoint.

18. V6 fresh B2/B3 1/10/100 us continuations all complete. B3's 10->100 us segment
    has 408 accepted S1 half steps, 32 time-error rejections, and no S0 recovery;
    its elapsed cost is 157.16 s. Development tests_dev09 has 48 passes, including
    budget interruption preserving the zero-step checkpoint/null certificates.
    V6 exhausts the preregistered six integrator revisions; no further algorithm
    revision will be hidden after the final independent validation.

19. Validation V3 freezes the numerical kernel at 07958c962. Fresh source-bound
    JUnit: **171 PASS, one retained old xfail**. All six new independent V3 states
    pass original-QP and 1 us motion checks. Both previously failed B2 states
    also pass a separate *known-failure regression*, with explicit S0 recovery.
    Those previous failures stay FAIL under their original V1/V2 protocols.
    The real WALL audit confirms exact triangle-subset correspondence to the
    closed geometry, positive triangle areas and three open boundary loops.
    The V3 campaign driver/GPU timing evidence is committed at 5198d3c70; its
    numerical-kernel hashes exactly match the 07958c962 frozen protocol.

20. V3 formal work begins from the original 3 s states, not prefix endpoints.
    The 3 h aggregate new-run allocation includes V1's 559.59 s B2 run and its
    safely interrupted B3 work; the budget is not reset. Each complete accepted
    ledger distinguishes S1 endpoint-matrix and S0 start-matrix equilibrium.
    End-to-end speedup, 1/5 ms scientific results and promotion remain pending
    until their actual runs and comparisons complete.

21. Formal V3 completes both ordinary and refined 1 ms cases. Ordinary B2/B3
    have 1036/996 accepted substeps, 94/63 rejections, and inclusive process-start
    wall times 322.005/333.424 s. Refined cases have 1988/1872 accepted substeps
    (see raw summaries for final values) and wall times 481.754/532.820 s. Frozen
    Stage4L same-environment 1 ms runs complete at 355.096/508.933 s. All new
    accepted intervals pass numerical/path gates with zero projection; no S0
    recovery occurs in these four trajectories. Both self comparisons pass.
    B2 ordinary vs legacy fails Q1 at 2.20205%; refined B2 passes at 1.75363%.
    Both B3 comparisons pass. Therefore V3 closes its own 5 ms gate and keeps
    B2 `BASELINE_DISAGREEMENT`; it is not marked fully successful.

22. V4 registers one stricter configuration of the **unchanged V6 kernel**:
    selected max macro dt 5 us / original error factors 0.25, independent
    refinement 2.5 us / factors 0.0625. All V1/V3 experimental costs carry forward.
    Reuse of completed V3 refined trajectories and references is explicit and
    never counted as reexecution or another repetition. Six fresh states at
    7/31/66% of both old trajectories are sealed; both parameter configurations
    pass on every state (12 probes). Fresh V4 tests: 171 pass, same old xfail.
    V4's new formal launches encounter the still-frozen V3 entry whitelist and
    are rejected before initialization, approximately 0.464 s each. Those logs
    and costs remain under `CAMPAIGN_V4`.

23. V4A fixes only entry routing, with an explicit versioned formal wrapper
    enforcing source, backend, parameters and budget before invoking the unchanged
    integration loop. It preserves the old whitelist. Seven new entry-contract
    tests bring fresh JUnit to 178 pass and one old xfail. The 12 V4 holdout probes
    are reused with their original source identities because the kernel and all
    numerical settings are unchanged. Executed V4A source: 7f83f802c. No further
    algorithm or parameter tuning is allowed after this validation. The new
    deeper-refinement trajectories are running; final 1/5 ms gates remain pending.

24. V4A completes both deeper 1 ms refinements: B2 3842 accepted / 169 rejected,
    B3 3542 / 135. All six selected-vs-refined / selected-vs-legacy /
    refined-vs-legacy comparisons pass. B2's largest registered quartile
    difference is 1.85457%; B3's is approximately 0.22390%. The actual accepted
    grids are finer, with median half-step dt 36.475 / 44.744 ns. Both refined
    runs use S1 only. The campaign therefore opens its 5 ms stage.

25. The selected B2 5 ms run `20261007T004943915783Z` stops at
    **4.314234468079 ms**, with `MINIMUM_PHYSICAL_PROGRESS`; it is **FAIL**, not a
    complete 5 ms trajectory. It has 3908 accepted substeps (3516 S1, 392 S0),
    202 rejected trials, 25 bounded S0 bursts, and approximately 775.08 s driver
    wall time. The last rejected intervals are `SWEEP_NOT_CERTIFIED`, not a
    native QP failure. At the final trial the WALL endpoint gap is
    -9.9999005491e-12 m, but the conservative continuous lower bound is
    -1.0000006237e-11 m after its ~1.05688e-16 m roundoff guard. Shrinking the
    macro interval below 1e-12 s is prohibited. The rejected candidate is never
    committed; original tolerance, geometry and guard are not relaxed. Remaining
    B3/refinement/reference runs continue under the same frozen protocol.

26. Selected B3 `20261007T010239725262Z` completes the entire 5 ms: 3478
    accepted S1 substeps, 140 LOCAL_ERROR rejections, no S0 bursts and zero
    projection. Minimum/median accepted dt are 4.81493 ns / 2.5 us. QP KKT
    <=9.98751e-10, original endpoint KKT <=9.98104e-10, resistance residual
    <=4.401e-14 and force residual <=5.332e-16 N. Continuous pair/WALL minima
    remain inside the unchanged tolerance. Total process-start wall time is
    794.075 s. Refinement and old-reference science are separate pending checks;
    B3 completion does not erase B2's long-run failure.

27. B2 refinement `20261007T011553819251Z` also stops at minimum physical
    progress, elapsed **4.311149113464 ms**, only 3.08535 us before the selected
    configuration. It accepts 6918 substeps (6768 S1, 150 S0), rejects 248 trials
    (179 LOCAL_ERROR, 10 nonlinear, 59 sweep), and uses ten bounded S0 bursts.
    Last accepted half-step dt reaches 7.62939e-13 s; the controller's minimum
    is on the full macro interval, whose accepted half-steps are separately
    recorded. Driver wall time is 1230.04 s. Every adopted interval retains
    its numerical/path checks and zero projection. Tightening the time-error
    budget therefore does not remove the longer-run certification blocker.
28. B3 independently refined 5 ms validation completes with 6730 adopted S1
    substeps and 135 LOCAL_ERROR rejections; no S0 recovery is used. Minimum /
    median adopted dt are 2.02543 ns / 1.25 us, compared with 4.81493 ns /
    2.5 us and 3478 substeps in its nominal 5 ms run. Maximum original S1 KKT
    is 9.98044e-10, resistance residual 4.64963e-14, and physical force residual
    5.37355e-16 N. Minimum continuous pair / WALL bounds are -9.94990e-12 /
    -9.99960e-12 m, with zero projection. Total process wall time is 1292.45 s.
    This establishes a completed finer numerical trajectory; full-duration
    reference science and same-environment timing remain separate checks.
29. Fresh same-environment Stage4L 5 ms references stop at their preregistered
    1800 s per-run budgets, not at a solver exception: B2 completes
    3.926404948744 ms in 13854 accepted steps; B3 completes 2.894642631625 ms
    in 11724 steps. Status is INCOMPLETE with REFERENCE_BUDGET, and process
    wall times are 1800.7434 / 1800.7977 s. Neither is a completion-speedup
    denominator. B2's frozen reference retains one original projection of
    7.61918e-12 m; B3 retains zero. These reference corrections are distinct
    from the R2 trajectories' zero projection.
30. B3's complete 5 ms nominal/refined comparison passes on the actual finer
    grid: endpoint Rg relative difference 6.21772e-9, series relative L2
    3.10658e-5, largest pair-quartile relative difference 3.22298e-6,
    centroid difference 5.47148e-11 m against a 1.78e-7 m budget. Maximum
    same-ID endpoint position difference is 2.45247e-10 m. This self-convergence
    result does not fill the missing fresh full-duration legacy comparison.
31. The repeated 1 ms B2 runs show different adaptive grids (1988, 2002,
    1986 accepted substeps) under the same frozen configuration. Similar timing
    alone is therefore not treated as scientific repeatability. A post-timing
    audit compares each completed repeat to the first same-method trajectory,
    and each new repeat to the first fresh legacy reference, using the unchanged
    registered metrics and exact original initial identities. All timing samples
    remain visible; this additional audit cannot improve the frozen campaign's
    status, and no solver or parameter change follows the observation.
32. Final independent checks pass: twelve reconstructed last-accepted original
    equations (including failed prefixes), ten final-kernel manufactured cases,
    and twelve post-timing repetition comparisons. Manufactured pair orders
    are 0.99417 / 0.99722; wall orders are 0.98180 / 0.99076 on 20/40/80 steps.
    These endpoint audits supplement, rather than replace, online checks on
    every adopted segment.
33. B2 has an additional **BASELINE_DISAGREEMENT** on the actual common
    3.926404948744 ms prefix. Nominal/refined centroid differences from the
    fresh legacy reference are 275.595 / 274.365 nm, above its 178 nm budget.
    Comparing the two new trajectories on exactly that clock gives only
    1.23126 nm centroid difference and passes the original self-comparison
    metrics even with conservative sampling bounds. Both Rg and quartile
    metrics against legacy remain within their 2% limits. This scientific
    disagreement is separate from B2's later continuous-WALL progress failure;
    it is not relabelled as a completed 5 ms test or proof either model is truth.
34. B3's full 5 ms supplementary comparison to the sealed archived Stage4L
    reference passes for both nominal and refined trajectories. Nominal Rg
    endpoint difference is 0.00239359%; pair Q1/Q2/Q3 differences are
    0.00269057 / 0.00231020 / 0.00520692%; centroid difference is 0.494132 nm.
    This does not fill the missing fresh same-environment 5 ms timing sample.
35. Formal timing ends with 9956.1615 s charged to new-method work and
    6182.8102 s to fresh references, including earlier failed/reused work under
    the recorded budget policy. Each case/method has three complete 1 ms
    measurements. No partial 5 ms run is converted into a completion speedup.
