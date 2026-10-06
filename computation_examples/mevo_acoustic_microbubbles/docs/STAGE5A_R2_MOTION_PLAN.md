# Stage 5A-R2: recovery, continuous motion and performance

Base: `7c39828d93509fff5551b7a16537ca3aa4ee9c23`, branch
`dev/stage5a-r2-motion-recovery`. Remote source matched this base at creation.
Independent server project: `/workspace/mevo_acoustic_microbubbles_stage5a_r2`.
Independent `.venv`; no old campaigns resumed. Old failures remain failures.

## Preregistered development and resource budget

Development set: the three immutable failed QPs named in the request, original
B2/B3 3 s checkpoints, manufactured problems, and prefixes up to 100 microseconds.
Holdout: 5 ms Stage4L trajectory frames floor(0.37*(Nframes-1)),
floor(0.73*(Nframes-1)), and Nframes-2, for both cases. These are separate from
the earlier 9 evenly spaced replay states and initial/midpoint benchmark states.
Manifest membership, file hashes and selected state hashes must be sealed before
opening these states to solver validation. Do not tune on holdout results.

Total experimental wall allocation: 12 hours (not a promise of success), including
up to 2 h QP/reconstruction/manufactured development, 3 h integrator/prefix
development, 3 h new formal trajectories/refinement, 3 h frozen reference timing,
and 1 h holdout/rendering/final audit. At most six QP strategy revisions and six
integrator revisions; each recorded with hypothesis, source and independent test.
Per scientific run <=1800 s, <=100000 accepted steps, <=10000 rejected trials;
minimum dt 1e-12 s, <=30 consecutive interval reductions. Native solver bounds
and a supervisor subprocess watchdog apply during rejected trials too.
No escalating these budgets after a failure. One failed candidate does not close
the whole development phase; remaining budget is used for traceable repairs.

Fresh quota: 7.68 cores, 45156925440 bytes RAM, RTX4090 24564 MiB. Up to 6 CPU
diagnostic workers, each BLAS/OpenMP=1; RAM budget 16 GiB, GPU budget 4 GiB.
One isolated scientific/timing process at a time. Changing-state diagnostics may
run in partitioned parallel processes, labelled throughput not latency.

## Frozen science and numerical certificate

Unchanged A0 primary/classical secondary force, field/APECSS/flow/real WALL,
Stage4L D+Lpp+Lpw, activation 0.1R, floor 0.001R, radius and contact tolerance
1e-11 m. No friction, adhesion, clipping, locking, deletion or initial projection.
Preserve signed initial gaps and all IDs/active masks.

Original Stage5A fixed-normalization KKT <=1e-9 (all original rows);
unconstrained resistance relative residual <=1e-11. Record all KKT terms,
primal/dual objectives and physical residual mappings. Original matrices remain
immutable during old-QP regression. No near-parallel-row deletion.
Native solver success is insufficient without the certificate. Inaccurate or
AlmostSolved candidates additionally require a second independent solver's
certified solution agreeing in velocity (relative max <=1e-7), objective
(relative <=1e-9), and reaction (relative L2 <=1e-7). Infeasible, nonconvex and
corrupt-input statuses are not eligible. This is a new acceptance protocol;
it never changes old results. Preserve a better pre-polish candidate.

## Discretization development

S0 uses A(x_n), frozen external forcing at (x_n,t_n), one start-anchored finite
horizon constraint set; it is a first-order diagnostic, not stiff implicit motion.
S1 solves x1=x0+dt*U with A(x1) and consistently end-anchored contact linearizations
expressed in the x0 displacement. External acoustic/flow forcing may be explicit
at the start, labelled semi-implicit. No derivative of resistance energy.
Bounded Picard/proximal iterations must converge to a specified residual; do not
stack all historical tangent rows. Frozen forces preserve exact A/G components.

True step-doubling error estimator: full vs two half steps, absolute position,
relative motion, pair gaps and wall gaps. Initial development budget uses
position/gap atol=0.002R, position rtol=0.01, gap rtol=0.02; refinement divides
all error tolerances by four in addition to maximum-dt refinement. Controller
remembers accepted dt/error, does not reset each output interval. Final parameters
and versioned algorithm are sealed before formal 1/5 ms validation.

Every adopted straight subsegment must pass full pair sweep and finite-radius
swept sphere vs real WALL triangle certificate, with open ports excluded.
No residual projection is planned; any later correction needs a new tested method.
QP/nonlinear/time/geometry budgets are separate. Atomic accepted checkpoints
include IDs, active, x/U, time, accepted cache, controller and source/env/config.

## Gates and comparisons

A: three real QPs + manufactured; B: first step, 1/10/100 us prefixes;
C: complete B2/B3 1 ms plus actual refinement; D: both complete 5 ms only after C.
No longer campaign. Failed candidates retain all evidence and close only their
own downstream gate. Keep original regression xfail, no new critical xfails.
Tests bound to actual source/config/environment hashes.

Same-ID largest initial cluster: endpoint Rg and time-series relative L2 <=2%;
pair-distance Q1/median/Q3 endpoint <=2%; centroid endpoint <=max(0.1R,
0.02*baseline centroid travel). Include absolute errors, rigid fit, graph changes,
gaps, IDs and actual sample age. Hard-contact chatter is descriptive, not a sole
gate. No coordinate interpolation. Align new output endpoints at 10 us clocks;
reference last-state sampling error is quantified, not hidden.
Self-convergence inconsistent with Stage4L => BASELINE_DISAGREEMENT, not relaxed
tolerances. Stage4L is a scoped numerical reference, not physical truth.

Performance layers: immutable QPs; changing real QPs; completed same-duration
end-to-end motion. Only the last defines completion speedup. Rerun frozen Stage4L
under same server/env/thread/pair backend conditions. Prefer 3 repetitions,
report fewer as insufficient repeatability within the fixed budget. Include all
failures/fallbacks/nonlinear work, initialization, transfer, IO and monitor cost.
Promotion requires both 5 ms numerical/scientific gates and comparable >=5x.
Numerical/scientific success below target is NUMERICAL_PASS_PERFORMANCE_PENDING.

## Evidence

New root `results/stage5/motion_solver_r2`. Strict finite JSON, immutable run IDs,
source/input/env hashes, lock/install/license records, QP NPZ, accepted/rejected
ledgers, atomic checkpoint, raw trajectories, PNG/PDF R1–R8 and same-clock MP4
after >=1 ms. Missing dynamics stays NOT_EVALUABLE. Distinguish accepted-state,
rejected-trial and terminal-failure statistics. Final gate separates implementation,
QP regression, continuous motion, science, performance and promotion. Push this
new branch only, retaining failures and actual large-file hashes/remote paths.

Official native APIs consulted: https://osqp.org/docs/interfaces/python.html ;
https://simple-robotics.github.io/proxsuite/ ;
https://clarabel.org/stable/python/getting_started_py/ . Installed-version API
introspection and interface tests accompany these references.
