# Stage5A OSQP numerical supplement — 2026-10-06

The user explicitly authorized installation of OSQP after the initial environment
audit. OSQP 1.1.3 and its missing Python dependencies are installed in the isolated
`.solver_deps/osqp_1p1p3` overlay. The original project virtual environment, NumPy,
SciPy, CuPy and frozen scientific inputs are unchanged. The package lock and
installation ledger are retained under `results/stage5/motion_solver/environment/`.

This supplement does not rewrite the original preregistration or failed runs.
The original active-set backend remains reproducible. Explicit configuration
`constrained_solver.qp_solver: osqp` selects `STAGE5A_OSQP_PC_V2`, using OSQP's
builtin CPU sparse direct algebra and CSC workspace updates. Warm starts come
only from accepted canonical particle/constraint caches. A failed warm solve
gets one cold retry, recorded separately with its full cost. Failed iterations
are never accepted positions. All original constraints remain present.

The original physical equations, integration tolerances, KKT <=1e-9, resistance
residual <=1e-11, geometry tolerance, 2% scientific comparisons, 5x runtime gate,
20/10 microsecond nominal intervals and staged 1/5 ms scope remain unchanged.
Before any new trajectory, require manufactured/regression tests and fixed-state
global/component and warm/cold comparisons. Old stopped jobs are not restarted;
any authorized new probe must have its own OSQP source/config/run identity.

Internal OSQP settings: eps_abs=1e-10 (one tenth of the unchanged independent KKT
gate), eps_rel=0, max_iter=100000, polishing=True, polish_refine_iter=10,
adaptive_rho_interval=50, check_termination=25. Status `solved` must also pass
the independent full KKT audit. An ADMM iteration-limit candidate may receive
one sparse saddle-point polish using its estimated active multipliers; this
polish must satisfy <=1e-11 linear residual and <=1e-9 KKT across ALL original
inequalities. Unsuccessful OSQP candidates themselves are never accepted.
The 30-particle manufactured chain exposed the need for this polish; the first
test log is retained. This is a solver repair before any new scientific run.
The first diagnostic used eps_abs=1e-12
and hit iteration limits; this diagnostic is preserved. No scientific run used
either setting when this supplement was written. Cold retry is a solver strategy,
not a relaxation of acceptance tolerances.

The new campaign first attempts B2/B3 nominal 1 ms with the original 1,800 s
per-case bound. Nominal completion is a prerequisite for spending further
runtime on dt/2. Both 1 ms numerical/scientific comparisons remain mandatory
before 5 ms. No stopped active-set run is resumed. OSQP run names carry an
`_osqp` suffix and use separate benchmark/gate/campaign summary files.

References: [OSQP Python interface](https://osqp.org/docs/interfaces/python.html),
[solver settings](https://osqp.org/docs/interfaces/solver_settings.html),
[algorithm](https://osqp.org/docs/solver/index.html).
