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
archived/fresh RMS is about 1.02 um, much larger than the 17.8 nm RMS resolution target
(0.01 of the 1.78 um bubble radius). R2 nominal/refined maxima are about 234 nm (B2) and 540 nm
(B3), so original-group agreement alone does not establish per-bubble accuracy.
Last-actual-frame timing uncertainty is retained; these are not interpolated
trajectories and all-active accuracy remains NOT_EVALUABLE pending refinement.

## Confirmations and original-start V1 campaign

The later exact/cold-query audit resolves H3 as CONFIRMED and the cache-artifact
explanation as REFUTED. See `B2_EVENT_ROOT_CAUSE.json` and the root-cause report.
Both terminal full-system fixtures advance without rejection/projection.
Original-start V1 nominal +1 ms and +5 ms runs complete for both cases and pass
the unchanged numerical gates. Process-inclusive +5 ms times are 728.034 s for
B2 and 1477.759 s for B3. B2 refined +5 ms also completes; its measured run-body
time is 1152.959 s. No independent global per-bubble error bound is claimed.

## Resource revision after measured full-task regression

V1 global QP was selected from isolated short changing-state substeps, where it
was fastest. The full B3 run exposes an unrepresented workload: 29,810 QP calls
cost 1196.99 s, versus 29,795 calls and 560.46 s in R2 component serial. R3/R2
all-active endpoint maximum difference is approximately 15.36 pm; different
physics or radically different motion does not explain this time difference.

The V1 **orchestrator only** was held while its current B2 refined science child
continued unmodified to completion. The scheduler was then closed before any
further tasks launched. Original `RESULT.json` remains untouched; separately
sealed `CLOSED_RESULT.json` records the complete child result and closure. Its
missing orchestrator process-latency measurement is not invented: for budget
accounting the measured run-body time plus a conservative 50 s is charged.
All V1 charges carry into the successor budget, totaling 4404.494 s in the
new-method group. No budget reset occurs.

The first adopted-horizon diagnostic overlapped an artifact transfer and is
retained as `LAYOUT_ADOPTED_HORIZON_WITH_SYNC.json`; a second isolated run resolves
that timing concern. All mathematical checks pass. A further complete-trial
benchmark uses actual formal states, all active bubbles and actual 5 µs
full/half1/half2 work, with three changing macrosteps per layout. It passes all
60 trials and all fixed-ID comparisons. At B3 +1.5/+2.5 ms, global medians are
1.618/1.321 s per trial, versus 0.351/0.288 s for component serial. Parallel
component layouts are slower than serial. Maximum layout trajectory difference
is 2.43e−17 m. Earlier isolated 50 ns tests did not cover these costly states.

V2 therefore keeps B2 global, where initial compression dominates the measured
task time, and changes **only B3 resource layout to component serial**. BLAS and
QP workers remain one; pair forces remain on the GPU. The scientific source
hashes and numerical driver AST are unchanged. The driver modification selects
and verifies the versioned resource protocol. Tolerances, full-system physics,
time-error settings, S0 policy and predeclared accuracy budgets are unchanged.
Equivalent B2 executions keep their actual V1 execution identities and are not
described as newly rerun. New numerical-state holdouts and source-bound tests
precede V2 formal B3 runs and repeated timing.
