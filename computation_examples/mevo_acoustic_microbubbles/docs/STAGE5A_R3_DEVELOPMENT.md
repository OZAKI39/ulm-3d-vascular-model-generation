# R3 development evidence

## Revision V1: finite WALL contact coverage and continuous proof

Hypothesis H1 (old certificate rejects a safe near-boundary path): CONFIRMED.
At B2 nominal terminal time 3.004314234468079 s, the actual full-system S0
candidate of 1.52587890625e-12 s is continuously SAFE under exact binary64-input
rational segment/finite-triangle distance: minimum gap -9.999900451221017e-12 m.
The R2 bound is -1.00000062368459e-11 m and rejects it. The original 1e-11 m
tolerance is unchanged. Refined gives the same classification.

Hypothesis H2 (removing the global floating-point guard suffices): REFUTED.
The same full-system states at 1 ns and 50 ns produce exactly UNSAFE candidates.
Their nominal minimum gaps are -1.001323999589104e-11 m and
-1.066787545854893e-11 m. Certificate improvement alone cannot recover motion.

Hypothesis H3 (missing simultaneous finite-feature contact coverage): supported,
pending actual trajectory repair. Particle ID 152 is within the original contact
tolerance of finite WALL triangles 13018 and 3178. Both closest points are in the
face interior (small positive barycentric edge weights); these are two distinct
noncoplanar faces, not an infinite-plane extension. The nearest normal changes
by about 0.1343. Whether the precise switching location includes query numerical
error remains UNRESOLVED until cold-query/independent-distance comparison.

V1 adds start-position supporting contact inequalities for nearby finite triangle
closest points, retaining the original nearest-WALL lubrication resistance and
the original pair contact rows. No extra lubrication term is added. The complete
186-particle coupled B2 fixtures are used; no particle is frozen or removed.
The start supporting plane is a conservative contact linearization of a convex
finite triangle and must be checked by refinement, especially around edges.

Continuous proof uses outward-rounded separating bounds, falling back to exact
rational minimum distance on the SAME straight candidate segment. A proved lower
bound at the tolerance is distinguished from the stored raw exact gap. No path
sampling, projection, position clipping or tolerance enlargement is introduced.
Robust-predicate background: Jonathan Shewchuk,
https://www.cs.cmu.edu/~quake/robust.html . This implementation uses Python exact
rational arithmetic and independently written finite-feature distance code.

Local probes begin with original accepted coordinates/cache and a explicitly
restarted bounded R3 controller. Original controller/cache provenance remains in
the fixture metadata. Such runs are development probes, not original-start gates.
V1 initial tests: analytic and random-bound properties, 13 local tests passed.
Execution results and source-bound server tests are recorded separately.

## Reference lineage before candidate selection

Hypothesis H4 (archived and fresh Stage4L are identical references): REFUTED.
Both cases have identical original particle arrays but archived output/macro
clock is 500 us versus fresh 10 us; pair backend setting is auto versus cupy.
The archived translation_backends.py hash differs: fresh adds contact active-set
cycle detection, primal QR refinement and feasible refinement. Their physics
coefficient formulae match. The full source difference is retained in diagnosis.
These factors require common-state attribution; neither reference is promoted.

All-active same-ID descriptive comparisons cover 186/189 particles, the original
75/91-ID groups and their complements. At their common partial endpoints, B2
archived/fresh RMS is about 1.02 um, much larger than the 1.78 um bubble-radius
resolution target. R2 nominal/refined maxima are about 234 nm (B2) and 540 nm
(B3), so original-group agreement alone does not establish per-bubble accuracy.
Last-actual-frame timing uncertainty is retained; these are not interpolated
trajectories and all-active accuracy remains NOT_EVALUABLE pending refinement.
