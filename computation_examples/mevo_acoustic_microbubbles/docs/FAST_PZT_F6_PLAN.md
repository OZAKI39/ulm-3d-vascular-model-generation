# F6 preregistered validation policy

Baseline c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0 is retained: the current remote
sync branches were inspected on 2026-10-10 and no newer MEVO sync existed.
Independent branch dev/fast-pzt-f6-target-navigation-20261010. F1–F5 are preserved.
New GPU compute uses the already installed CUDA12 Warp 1.18.0; no driver changes.

Admission A: analytical multilinear gradient and CPU/GPU relative agreement
2e-13 on a bounded analytic test, finite differences 2e-9; real primary force
agreement 3e-11 relative / 1e-26 N absolute. Strict domain rejection and explicit
one-sided face policy. Legacy F5 ON/OFF 1 ms prefixes must agree bitwise.
RPY: one-bubble limit, GPU/CPU, reciprocity, symmetry and positive eigenvalues on
one bubble, two bubbles (including overlapping regularization), four bubbles.
These are software checks, not physical validation.

Targets, initial states, horizons and source caps will be frozen in protocol.json
before trajectory integration. Only static field/velocity diagnostics precede that
registration. Each environment receives OFF, unchanged and target-phase controls,
plus at most one bounded optimizer candidate. No controls are tuned to a trajectory.
Numerical overlap <=0.01R and per-pass correction <=0.1R; identical wall budgets,
zero error bits, finite state, no missing IDs for a navigation success.
Arrival radius 2 micrometres, sustained residence 50 ms, largest connected fraction
>=0.9 throughout residence. Initial target distance 10 micrometres (5 tolerances).
Minimum control effect versus flow-only and unchanged: 0.5 micrometres and >3 times
the key timestep centroid discrepancy. A successful flow-only arrival alone is not
PZT navigation. Continuous saved residence is a sampled lower-bound diagnostic,
with 0.5 ms sampling uncertainty at each boundary.

Key dt comparison 10 vs 5 microseconds: maximum common-time centroid/target-error
change <=0.5 micrometres, Rg relative <=5%, largest-fraction change <=0.1; both
must satisfy the unmodified overlap/wall budgets. No bitwise microscopic packing
claim. Reference translation convergence <=2% endpoint change, with absolute
1e-15 m reporting denominator and signed endpoints retained. Failed/reference
nonconvergence is retained; no universal correction factor is inferred.

AI evidence remains EXPLORATORY. Train/validation/test partitioning is by
medium+target scenario; all controls and timestep/mobility variants stay together.
No large random campaign, native CFD rerun, acoustic rerun or response-table rebuild.
