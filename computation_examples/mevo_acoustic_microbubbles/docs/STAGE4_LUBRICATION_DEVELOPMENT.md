# Stage 4L development record

## Isolation and preserved legacy state

Local branch: `dev/stage4l-lubrication-20261005`, created from the GitHub snapshot commit `aaf2549fc506933df74229be8e0ec9717f97f7ed` (`sync/mevo-acoustic-microbubbles-stage4-20261005`). All 278 tracked original project source files were byte-identical to that snapshot before development.

Local project: `/home/lzy/projects/computation_examples/_archives/mevo_stage4l_20261005/computation_examples/mevo_acoustic_microbubbles`.
Server project: `/workspace/mevo_acoustic_microbubbles_stage4l`, branch `stage4l-lubrication-20261005`, based on source commit `a6c3153f484113d37ffb5e9b4b57f0c6594ef396` with the same starting source bytes. Its distinct Git ancestry is recorded in run identities; it is not misrepresented as the monorepo snapshot ancestry.

The original server controller remains paused. The already-running Stage J child was observed but not signalled. No old K, final acceptance, multiseed, or BraVa task was resumed. No original gold, configuration, checkpoint or result identity is changed. New evidence lives only below `results/stage4/lubrication_upgrade/` in the isolated project.

## Implemented

- Independent `TranslationBackend`, legacy Stokes and coupled lubrication implementations; default remains legacy.
- Monodisperse normal particle-particle and particle-wall resistance, sparse direct solve, strict residual check and bitwise far-field bypass.
- Exact cached VTK wall gap/normal, predictive closing-gap timestep guard without floor collapse.
- A-metric frictionless instantaneous and finite-horizon contact projection, residual geometric projection diagnostics including RK stages.
- Finite checkpoint runner with input/source hashes, exclusive result paths, copied read-only 3 s inputs and paired unchanged force inputs.
- Same-ID cluster metrics and contact lifetimes at equal physical sampling intervals. Statistics may interpolate accepted states; scientific trajectory output is never filtered.

## Software evidence and experiment protocol

Initial server checks: **47 passed**, covering the requested T1–T10 cases, matrix symmetry/SPD and CG agreement, factor-two convention, predictive crossing, unequal radii rejection, A-metric contact, legacy tracker and ordered contact projection. The local shared Python environment lacks pandas; full tests use the already provisioned server project environment without changing it.

At the sealed t=3 s checkpoints, B2/B3 have 186/189 active particles, 327/344 pairs in lubrication range, and about 150 particles participating in hard contacts. This motivates A-metric contact rather than continued use of a diagonal mobility constraint.

S0 order: paired B2/B3 1 ms, then 5 ms only after stable behavior. Nominal clock, contact tolerance, initial IDs, positions, field, flow, radial response, and acoustic forces are identical between arms. Legacy and new arms use the same ordered geometric fallback implementation. All failures and partial results remain saved as failures.

The first new-model attempts stopped because the diagnostic logger duplicated `time_s`; these are recorder failures, not solver or scientific results. Their directories remain intact. The logger was fixed and new result IDs were generated. The completed legacy 1 ms runs are reused for comparison.

## Gates still pending

Dense rearrangement and severe-projection gate; new primary B0–B3; new B2/B3 dt/2; staged secondary 0.1 s and 0.5–1 s; possible 6 s secondary only after the runtime/behavior gate. No inherited legacy timestep PASS is transferred to this backend. No overall Stage 4 PASS is claimed.

See `LUBRICATION_TRANSLATION_MODEL.md` for equations, contact semantics, physical scope and references. This file will be updated with measured probe outcomes and an explicit final gate decision.
