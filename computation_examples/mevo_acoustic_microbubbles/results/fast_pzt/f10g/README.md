# F10G evidence

The exact discrete lubrication reference, frozen-normal implicit prototype, two
32-bubble static full-environment checks and two bounded trajectories passed their
specified numerical tests. Original M0 defaults and budgets were not changed.

SOFT_CONTACT_SHORT_TRAJECTORY_PASS applies only to B2/2kc, 5 ms O/1102 and
2.35 ms preparation1104. There is no new 10 ms or OFF/S/O triplet. No new seed.
Model disagreement is not evidence of physical truth; same-model dt convergence
is untested. Small repulsions after true geometric separation are explicitly
recorded: this frozen-normal spring is not an exact geometric unilateral law.
PRODUCTION_NOT_ADMITTED / PHYSICS_UNCALIBRATED / training_ready=false.

Remote output folders preserve immutable execution identities, configs, logs,
all-step positions/reactions, real accepted edges, budget records and timings.
SOFT_STEPS.json.gz contains the actual soft-contact evolution. The inherited
recorder has some allocation/compatibility labels; EXECUTION_INTERPRETATION.json
clarifies that the new adapter executes zero PBD passes and no CUDA Graph.

lubrication_local_v1/prototype_local_v1 are earlier local CPU reference checks,
not extra independent scientific samples. Final analysis/figures use server v1
reference/prototype and performance_v2 (actual masked/ordered kernel family).
All attempts, including early unmasked timing and deliberate rejection tests,
are retained and charged. Remote hashes and Git snapshots are independently
verified. Four PNG/PDF main groups plus the Figure 3 network supplement are in
figures/. Progress document before/after/backup/diff hashes are preserved separately.
