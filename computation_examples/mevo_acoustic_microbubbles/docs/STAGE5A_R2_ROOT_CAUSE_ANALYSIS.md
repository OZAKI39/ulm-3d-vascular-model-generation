# Stage 5A-R2 root-cause analysis

Evidence root: `results/stage5/motion_solver_r2/diagnosis/`.
These findings are new read-only reconstructions, not reinterpretation of old
failures as passes. All three reconstructed A/G matrices match the immutable
failure fixtures exactly (zero entry differences). Row metadata includes stable
IDs, stages, normals, anchors, original/scaled lower bounds and wall exterior
feature indices. Physical velocity/resistance/force scales are recorded.

| Hypothesis | Finding | Evidence and limits |
|---|---|---|
| A: stacked start/predictor near-dependent tangents | CONFIRMED | `reconstruction_v1/REPLAY.json`, spectra and near-row lower differences; identical A/G reconstructions. Dependence is not infeasibility and no near-parallel rows were deleted. |
| B: guessed active set/dual polish fails | CONFIRMED | Re-executed old exceptions retain pre/post-polish residuals; B2/B3 guessed saddle polish has tiny linear residual but fails original KKT. Native ProxQP solves all rows to <=1e-12 without that polish. |
| C: native stopping differs from independent audit | CONFIRMED | Old B3 inaccurate candidate passes 1e-9 KKT but old status rule rejects; new native Clarabel can report Solved with KKT >1e-9. Neither native status alone nor a saddle linear residual is sufficient. |
| D: coefficient-lagged/time error limits old PC | CONFIRMED | All 20 original LOCAL_ERROR trials reproduced; responsible contacts, normals, activation/floor switches and coefficients saved. At B2 39.0625 ns, IDs 110/170 cross activation with gap 190.59->147.16 nm and error 1.7409. B3 IDs 206/221 have coefficient change 5.78e-9->2.70e-8 N s/m, error 2.1831. These are actual changes, not only algebraic failures. |
| E: QP exception bypasses interval recovery | CONFIRMED | Old interval catches only SubdivideInterval. Replays terminate with MotionSolverFailure. R2 uses bounded native alternatives and classifies subsequent interval recovery. |
| F: stage workspace thrashing | CONFIRMED | Old component key contains IDs only; start/corrector G patterns alternate in the same workspace. Reconstructed stages have different row sets. R2 stage+solver+component+pattern identities separate them and numeric updates are explicit. |
| G: known feasible predictor ignored | CONFIRMED | Reconstructed predictor minimum slack: -2.57e-16 / -4.61e-19 / -4.34e-18 on all original rows. Zero violates lower bounds by 1.69e-4 / 2.84e-7 / 2.02e-5 (normalized). Old OSQP ignores its feasible parameter. |
| H: continuous wall certificate absent | CONFIRMED | Old source checks wall endpoints only. New test has endpoint-clear path crossing a triangle, plus finite-radius edge graze. R2 uses full swept finite sphere–WALL certification. This confirms missing coverage, not proof an old accepted path actually crossed WALL. |

Direct evidence: [reconstructed original trials and feasible predictors](../results/stage5/motion_solver_r2/diagnosis/reconstruction_v1/REPLAY.json),
[native solver comparison](../results/stage5/motion_solver_r2/diagnosis/native_v1/SUMMARY.json),
[pre/post polish residuals](../results/stage5/motion_solver_r2/diagnosis/polish_parts/SUMMARY.json),
[old exception and workspace handling](../src/mevo_mb/bubbles/osqp_motion.py),
[old interval driver](../scripts/stage5a/run_probe.py), and
[new continuous geometry regressions](../tests/test_stage5a_r2.py).
The numerical WALL discontinuity is recorded in
[nonlinear reconstruction](../results/stage5/motion_solver_r2/diagnosis/nonlinear_v4/SUMMARY.json);
the independent final validation remains in the
[separated outcome gate](../results/stage5/motion_solver_r2/final/MOTION_SOLVER_R2_GATE.json).

The initial checkpoints have strong transient motion and many near-contact
switches. The old 20 us PC trial changes hundreds of pair states and is far from
the accepted time-error budget. A faster frozen QP therefore cannot alone repair
continuous motion. R2 separates this from algebraic recovery.

The first new endpoint-contact S1 candidate introduced a distinct intermediate
path problem in initial negative-gap tolerance bands. Its rejected sweep is
retained; V2 fixes start tangents and iterates endpoint resistance. V3 uses bounded
Anderson iterations. No initial coordinates, force, radius, floor or geometric
tolerance were changed. Subsequent numerical/scientific/holdout results remain
separate from these root-cause findings.

Later execution evidence distinguishes two remaining issues. V1's six unseen
QP states are all solvable, but one subsequent motion probe stalls at 0.211 us
with a nonlinear map difference that does not vanish as dt shrinks. The precise
trigger is UNRESOLVED; it is not evidence that those original QPs are infeasible.
V1's complete B2 1 ms run spends 466.39 of 559.59 s in QP calls. Late-state
substep benchmarks confirm that the initial global-layout timing did not predict
the full trajectory: component serial is faster on the later development states.
This is an execution-layout finding, not a change to any physical force.

New nonlinear diagnosis (V4/V2 failed states): the earlier trigger is now directly
observed. At original-ID 152, closest-feature WALL normals differ by 0.1343 for
trial coordinates separated by approximately 3.9e-19 m; the nonzero wall
lubrication block changes and the full resistance matrix differs by about 1.2%.
The native QPs remain certified. At tiny dt the nonlinear map can oscillate across
that geometry discontinuity. This confirms a non-smooth coefficient trigger;
mathematical nonexistence or uniqueness of the implicit step is **UNRESOLVED**.
The actual WALL mesh and original gradient evaluator are unchanged.

Endpoint wall contact updates (V5) did not resolve these states. V6 uses explicitly
labelled, bounded S0 start-state recovery only after such a failed nonlinear solve
in the wall-lubrication activation band. This restores the two known short probes
under the new hybrid protocol without claiming that their S1 equation converged.
Every adopted segment still passes original-QP, time-error and continuous geometry
checks. Six newly sealed independent probes also pass. Whole-trajectory scientific
and performance gates are separate and cannot be inferred from these local probes.

The supplementary `diagnosis/polish_parts/` audit applies the exact old saddle
polish equations to the immutable saved OSQP candidates. Before polishing,
B2/B3 KKT are 1.07494e-9 / 7.84371e-10; after polishing they are approximately
0.85746 / 0.13004. The dominant post-polish failure is **dual feasibility**:
the guessed equality set yields negative multipliers. A small saddle residual
therefore does not certify the unilateral QP. Local SciPy pivoting changes a few
last digits in this ill-conditioned solve; the original server replay is kept.

Evidence correction: the early `benchmarks/changing_qp/RESULTS.json` used an
adapter whose default unit scales were 1. Its fields labelled as N or m/s are
normalized quantities, not SI. `AUDITED_SI.json` independently reconstructs the
120 archived candidates with the saved velocity/resistance scales, retains the
original file hash, and compares U, objective and G^T lambda across solvers.
No original record is overwritten, and no multiplier uniqueness is assumed.

The raw summary field `qp_calls` counts whole constrained-stage calls, not every
native component solve. Accepted-substep component solves, workspace hits and
native iterations are counted separately in `final/RUN_AUDIT.json`. Completed
full/half trials include rejected time-error work; unavailable terminal native
iteration histories remain explicitly unavailable rather than fabricated.

Clarabel's native tolerance does not monotonically determine the fixed audit.
At 1e-10, its three fixture statuses are Solved, with original KKT
1.079e-3 / 6.964e-10 / 4.610e-9. At 1e-12 all three return
InsufficientProgress; the large audit term is active slack, despite small
stationarity. R3 deliberately labels its plotted internal tolerance (1e-12),
and the complete two-tolerance table remains in `native_v1/SUMMARY.json`.
This is evidence for selecting certified ProxQP, not a general claim that
Clarabel cannot solve such problems or that tighter requested tolerances always
produce more accurate candidates.

The later B2 5 ms validation exposes a remaining finite-precision progress
blocker, distinct from the original failed QPs. At elapsed 4.314234468079 ms,
repeated S0 candidates have endpoint WALL gaps inside the original tolerance,
but the continuous certificate's conservative roundoff guard places its lower
bound just below -1e-11 m. The last bound misses by about 6.24e-18 m. The driver
correctly refuses the candidate and stops at its minimum physical-progress
limit. This proves an inability to certify continued motion under the frozen
protocol; it does **not** prove an accepted path penetrated WALL beyond tolerance.
Compatibility between displacement constraints, finite-precision safety margin
and the continuous certificate remains unresolved for longer motion. No guard
or tolerance is weakened to make this run pass. The 5 ms run's complete accepted
and rejected ledgers, checkpoint and terminal error are retained.
The independently finer B2 trajectory stops by the same mechanism at
4.311149113464 ms. Its stop differs by only 3.08535 us despite 6918 instead of
3908 accepted substeps. This supplies refinement evidence that the remaining
progress limit is not removed by the current time-error refinement alone.
