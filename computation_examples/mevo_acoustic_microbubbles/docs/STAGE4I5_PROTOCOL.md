# Stage 4I.5 — preserved gold and performance gates

The original Stage 4I B2/B3 cases completed naturally. Their unique full gold
reference is `results/stage4/performance_refactor/gold_baseline/`, with manifest
SHA256 `3d768c763a3f20625b8064039fea5cddd35cc2765810a3511b46eef7a288cef8`.
All 42 preserved files (16,393,783,104 bytes) were independently copied and hashed
on the local workstation. Gold is read-only; do not rerun, overwrite or rename it.

## Short replay protocol

`protocol/PROTOCOL.json` was declared before candidate timing. Both cases use
S1 = 0–0.04 s. S2 selects the largest mean saved contact-pair fraction over two
adjacent saved states, with earliest ties: B2 5.94–5.96 s and B3 4.88–4.90 s.
This choice uses gold contact evidence, never candidate speed or agreement.

P0 retains the original synchronous RK2, global NNLS and the exact per-particle
interior-wall certificate already used in gold. It recomputes its initial fields
instead of injecting checkpoint velocities. Reconstructed initial values and the
endpoint are compared directly with gold. All accepted starts, forces, velocities,
dt values, global IDs, wall flags and contact-event hashes are saved separately.

The implicit `auto` pair-backend input was identified after sealing. Its observed
SHA256 agrees between local and server and its rows agree with the historical
Stage 2 benchmark. It is frozen under `protocol/PAIR_BACKEND_STATE.json`, outside
the immutable gold manifest. This is post-seal corroboration, not a claimed
run-start hash. The crossover is 256: GPU at N >= 256 unless resource fallback,
CPU below that count. Original lifetime GPU time was not instrumented.

## Candidate structure

- `GlobalSynchronousIntegrator` is an explicit name for the retained original
  controller/evaluator/stepper; none of the original 42 production files changed.
- `ExactCoordinateWallCache` uses exact coordinate byte keys, with no bins or
  approximate values. It handles ragged projection subsets and mesh invalidation.
- `PersistentVTKWallBackend` uses the official `vtkStaticCellLocator`, one build,
  reusable generic cells, VTK triangle closest points and geometric pseudonormals.
  The original VTK implementation already builds once and reuses thread-local
  cells; its separate value and gradient calls repeat closest-point work.
- `ContactNeighborCache` caches an enlarged cKDTree broad phase. Displacement,
  shape, radius and gap changes cause conservative rebuilds. The original sorted
  exact narrow phase is unchanged.
- `WarpWallBackend` is optional and experimental. Its install step puts Warp 1.17.0 only in
  `.optional/warp-1.17.0` with `uv --target --no-deps`; its activation is per process.
  The main venv is compared before/after. Float32 geometry, zero-copy CuPy buffers,
  host transfers, mesh winding and projection/trajectory errors are explicit.
- `OrderedContactSweep` is an additional profile-driven CPU candidate: the ordered
  per-pair correction arithmetic is compiled, while the original cKDTree, wall
  projection, tolerance, NNLS and integrator remain unchanged. Its local scalar
  norm first failed bitwise on the server; the corrected implementation reuses
  NumPy's existing ILP64 BLAS norm where required. It records compiler flags,
  binary and BLAS hashes, forbids fast-math/reassociation in corrections, and
  must pass an independent 20-state gold shadow plus the same four trajectories.
  A C compiler is needed only for this optional backend; default Python/VTK works
  without it. No custom CUDA extension is introduced.

The same four S1/S2 windows test each candidate independently. No full baseline is
rerun. P1 requires identical IDs, dt/contact/wall/validity records, identical
assignments and maximum position difference <= 1e-12 m. The 10,000-point geometry
fixture covers lumen, apex, daughters, near-wall/outside and actual gold states.
Neighbor stress tests require zero missed contacts. Passing a query fixture alone
does not authorize production use. The 2x overall engineering target is not a pass
condition; correctness is always required.

An additional official VTK batch option, `vtkSampleImplicitFunctionFilter`, was
screened from the installed-version source. VTK 9.7.1 hardcodes `vtkFloatArray`
for both distance and gradients, so it was not introduced into the float64
equivalence path. This is a rejected feasibility option, not a measured speedup.

`SimulationTerminationPolicy` distinguishes observation objectives from a proven
terminal-only endpoint. The current proximal classification plane fails the
physical terminal-deactivation audit. Existing absorbing-boundary behavior is
preserved for gold comparisons; no additional deletion or plateau stop is enabled.
Formal J checks that its completed observation satisfies its declared objective;
an earlier empty legacy domain without proven terminality is an explicit failure.

## Execution order

Finite Supervisor jobs run sequential baseline profiles, baseline timings, then
candidate timings. Candidate timing never overlaps another candidate. cProfile
runtime is not used as a speedup denominator. Postprocessing computes limiter
categories and slip Reynolds numbers from the frozen flow without inserting
diagnostic wall queries into the scientific stepping loop.

Before adopting a candidate with more than 5% component-timed aggregate gain,
`backend_decision` runs a timing confirmation on the same four declared windows.
`replay --no-component-timers` removes all nested component wrappers and per-point
VTK distance/normal timers. This addresses instrumentation cost that depends on
the query API. The original wall certificate and scientific operations remain;
common accepted-step recording, timestep and classifier timers remain explicit.
Minimal-timer P0 must reproduce instrumented P0 bitwise, and each candidate must
pass both its instrumented counterpart and minimal-timer P0 under the original
strict regression. Adoption also requires more than 5% confirmed aggregate gain.
These serial short timings are recorded separately under `timing_confirmation/`
and `timing_confirmation_report/`; they do not replace the component breakdown
or imply a measured full-six-second production speedup. Closed confirmations are
reused only when the replay sources, protocol and measured environment match.

1. `python -m scripts.stage4i5.prepare` — **one-time only**, already completed.
2. `python -m scripts.stage4i5.replay --profile` — four short P0 profiles.
3. `python -m scripts.stage4i5.replay` — matched P0 timing, without cProfile.
4. `python -m scripts.stage4i5.run_candidates --wait-for <finite-baseline-service>`
   — sequential fixtures, P1 comparisons and optional Warp benchmark.
5. Review limiter evidence; implement P3 only for a supported small stiff subset.
6. `python -m scripts.stage4i5.backend_decision` — requires the complete evidence.
7. Only then run revised `scripts/stage4/convergence.py`, followed by K and audit.

`run_post_candidates` waits for the candidate campaign's natural completion,
executes the additional ordered-sweep gates, selects only an eligible backend,
then runs J, K, audit and report in order. A failed required phase stops the chain
with its evidence intact; it never resumes or reruns the original gold cases.

J includes B2 primary+secondary at half the nominal dt **and** half the displacement
limit, primary B2/B3 checks and one fine acoustic steering field. The observation
remains 6 s because unassigned particles make the plateau nonterminal. Per-ID
outcomes, A/B fractions of the initial population, crossing times and early
cluster metrics accompany conditional exit fractions. No primary-only dt/2 test
is relabelled as secondary timestep validation.

## Interpretation

`FAR_FIELD` is a distance label, not complete physical validity. Pair fractions
weighted by evaluation count, accepted time or wall time differ from the fraction
of computational time whose state contains any near/contact pair. These quantities
are reported separately. Finite slip Re, near-wall hydrodynamics and close-pair
model extrapolation remain warnings regardless of runtime or numerical smoothness.

Implementation references: [VTK static locator](https://vtk.org/doc/nightly/html/classvtkStaticCellLocator.html),
[installed VTK 9.7.1 distance implementation](https://raw.githubusercontent.com/Kitware/VTK/v9.7.1/Filters/Core/vtkImplicitPolyDataDistance.cxx),
[VTK 9.7.1 batch-filter precision](https://raw.githubusercontent.com/Kitware/VTK/v9.7.1/Filters/General/vtkSampleImplicitFunctionFilter.cxx),
[Warp 1.17.0 release](https://github.com/NVIDIA/warp/releases/tag/v1.17.0),
[Warp official mesh-query definitions](https://github.com/NVIDIA/warp/blob/v1.17.0/warp/_src/builtins.py),
[Warp license](https://github.com/NVIDIA/warp/blob/v1.17.0/LICENSE.md).
