# Stage 5A-R3: B2 wall event and reference lineage

The actual development base is R2 archive commit
`07db8ed35f3a460a4ba315e42b15d257f8518b0a`. The R2 execution commits
`7f83f802c6cc435918ea25d69dcbb0f4f4dacce9` and reused
`5198d3c705c6a3735dcf566f53810274c40a6624` identify executions, not the archive.
All R2 results remain unchanged. Machine-readable source/input checks are in
`results/stage5/motion_solver_r3/final/R3_BASE_IDENTITY.json`.

## Confirmed B2 mechanism

At nominal absolute time **3.004314234468079 s**, bubble **152** approaches
simultaneous contact with finite WALL triangles **3178 and 13018**. The nearest
points are inside the two finite faces, close to their common edge. The contact
normals differ by approximately 0.1342982 in Euclidean norm. This is a real
nearest-feature transition: independently reconstructed exact distances identify
the change, and cold VTK queries reproduce the cached distances and normals
exactly. A stale wall-query cache is therefore REFUTED for these states.

R2 includes one nearest-wall hard-contact row. At this event that row alone does
not cover both close finite faces. R2's endpoint-resistance S1 iteration encounters
the normal transition; bounded start-resistance S0 recovery then repeatedly
shrinks the physical interval. Two different geometric situations occur:

| Actual full-system nominal candidate | Exact continuous minimum gap | Classification |
|---|---:|---|
| 1.52587890625 ps | −9.999900451221017e−12 m | SAFE |
| 1 ns | −1.001323999589104e−11 m | UNSAFE |
| 50 ns | −1.066787545854893e−11 m | UNSAFE |

The first candidate is safe at the original **1e−11 m** tolerance, but the R2
lower bound is −1.0000006236845907e−11 m. Its global floating-point guard makes
this particular proof too conservative. The larger candidates really exceed the
tolerance. Removing the guard alone cannot repair the motion. Both nominal and
refined terminal fixtures show these classifications.

The 1.52587890625 ps candidate changes the binary64 timestamp by
1.5258905250448151 ps and moves ID 152 by approximately 3.37e−16 m. Thus
`MINIMUM_PHYSICAL_PROGRESS` is an algorithmic request falling below the existing
1 ps floor, not evidence that time or position has no representable update.
`NONLINEAR_NOT_CONVERGED` concerns rejected S1 endpoint iterates; their residuals
are not accepted-state residuals. `SWEEP_NOT_CERTIFIED` includes both a safe but
unproved candidate and genuinely unsafe candidates. The exact audit separates
them rather than treating every rejection as wall penetration.

## Numerical repair and unchanged physics

`motion_r3.py` adds hard-contact supporting rows for close **finite-triangle**
closest points at the actual start position. Each row is a local conservative
supporting constraint for that finite convex feature, with a distance/reach
query; arbitrary distant planes are not made into walls. The normal lubrication
resistance still uses the original **single nearest WALL query**. No second wall
resistance, acoustic force, flow correction, artificial adhesion, radius change,
particle removal, velocity clipping, wall smoothing or projection is introduced.

The continuous certificate first uses complete conservative AABB candidates and
outward-rounded separating lower bounds. Ambiguous cases fall back to exact
rational squared distance between the same adopted straight sphere-center
segment and finite triangles, including face/edge/vertex and degenerate cases.
This is not a sampled-path test. Proof refinement does not advance physical time
or write accepted trajectory frames. Raw negative gaps inside the original
tolerance remain negative in the records.

Signed initial membership remains the original geometry's responsibility. Open
inlet/outlet caps are not WALL triangles. The sphere/WALL barrier argument starts
from the verified interior state; crossing a WALL would require a continuously
unsafe segment. This certificate does not silently turn open ports into closed
surfaces or establish transport semantics beyond the present no-exit probes.
Open-boundary and outside-initial-state cases are separately tested.

Both actual terminal fixtures contain all **186 active B2 particles**, original
IDs/radii, original coordinates and stored accepted caches. Each advances 10 µs
with 16 accepted half steps, no rejection and zero projection. Earlier-window
cache/controller reconstruction is explicitly identified in fixture metadata;
it is not described as a stored exact checkpoint. Original-start campaign results
are separate from these local repairs and appear in the final report and gate.

## S0, S1 and precision limits

S1 evaluates resistance at the endpoint and uses the actual frozen start forcing
and contact protocol. S0 evaluates resistance at the start. The mixed method is
not described as fully implicit. The R2 recovery burst/count limits stay fixed.
S0 has a start-equation certificate and **NOT_APPLICABLE** endpoint nonlinear KKT,
represented by null, never zero.

Local S0 caps were independently halved **50/25/12.5 ns**. S1 local intervals use
three precision levels, all active IDs, and raw/translation-removed/rigid-fit
errors. The finest initial probes hit their finite 300 s budget; their actual
common 10 µs interval is retained without claiming a complete 20 µs result.
Very small nonmonotone late/event differences are reported rather than converted
into a claimed convergence order. Local consistency is not an independent
global 5 ms trajectory error bound.

## Archived and fresh references are different numerical executions

Both legacy lineages use identical original particle arrays, but the archived
macro/output clock is **500 µs** and the fresh comparison clock **10 µs**. Pair
backend settings differ (`auto` versus `cupy`); the actual auto device is not
inferred from its name. The archived contact-solver source differs from the
fresh source by cycle detection and primal/feasibility refinements. The source
diff and hashes remain in `diagnosis/REFERENCE_SOURCE_DIFFERENCE.patch` and the
lineage audit. Neither lineage is selected because it favors R3.

Controlled full-system short replays isolate these differences. Changing only
the archived/fresh contact source at the same clock gives identical endpoints
in the tested 40 µs windows. Changing only 500/10 µs cutoff timing produces B2
RMS/max differences of about **28.65/128.77 nm** and B3 **1.55/5.96 nm** at 40 µs.
Changing the auto/GPU setting alone gives floating-point-scale differences in
these probes. These results identify an early timing sensitivity, but do not
prove that it explains the entire later reference disagreement.

Common-state cross-replays compare forces, flow, original A/b/G/lower bounds and
identical QP horizons before comparing trajectories. Original A and RHS agree
to floating-point precision; legacy/R2/R3 velocities on the same QPs differ by
roughly 1e−12 m/s, while retaining independent certificates. Local three-level
legacy refinement and cross-method refinement are retained. Large full-interval
same-ID discrepancies still require a global accuracy explanation.

All audits cover **186/189 active IDs**, plus the original **75/91-ID groups** and
their complements. No ID matching or coordinate interpolation is used. Actual
frame times and sampling displacement bounds accompany every common-clock
comparison. A small centroid difference does not imply small individual errors.
Without an independent global bound, the predeclared 17.8 nm RMS / 89 nm maximum
resolution targets are descriptive checks and accuracy stays **NOT_EVALUABLE**.
`BASELINE_DISAGREEMENT` blocks automatic reference promotion.

## Performance and provenance constraints

V2 reuses only a same-input already-computed free solution, checks its original
KKT before a fast path, and vectorizes component-row ownership. Changed matrices
still receive numerical updates. Native ProxQP and independently certified
fallback behavior are retained. Isolated changing-state tests compare global,
component serial and persistent component parallel layouts; full-task results
must determine whether their short-window advantage carries over.

Formal execution is source-frozen at
`fa2fba5608037ed45143c8ede501cadf415e9180` in an immutable server execution tree.
Postprocessing/documentation commits are not relabeled as execution commits.
One early local diagnostic campaign imported V1 before two later source ZIPs
captured a deployed V2 driver. Numerical kernels matched, but driver provenance
did not. The original ZIPs and identities are preserved, with a reconstructed
loaded-source ZIP and explicit `LOCAL_V1_PROVENANCE_CORRECTION.json`. Later
benchmarks, attribution and formal runs use immutable execution trees.

The strict bitwise restart check failed by one coordinate ULP (2.1684e−19 m).
Its failure is retained. A separate predeclared-resolution restart check passes,
as does exact accepted-cache serialization. These are distinct claims. All
failed candidates, budget terminations and test invocations remain auditable.
