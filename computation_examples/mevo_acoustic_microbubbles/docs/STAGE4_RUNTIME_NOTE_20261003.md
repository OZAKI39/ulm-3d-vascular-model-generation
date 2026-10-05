# Stage 4 runtime note — 2026-10-03

This is a dated progress/evidence note, not the final Stage 4 acceptance report.

## Completed primary study

All 40 primary-only cases completed: ten paired seeds, four controls, 256 bubbles per case, 6 s observation, 0.5 ms nominal timestep. No failed case was dropped. Same-seed initial arrays match by byte hash. The 20,000-draw bootstrap resamples whole paired seeds.

| Control | Mean conditional f_A | Mean A / initial | Mean B / initial | Mean unassigned / initial |
|---|---:|---:|---:|---:|
| B0 | 51.067% | 49.531% | 47.461% | 3.008% |
| B1 | 49.622% | 12.188% | 12.344% | 75.469% |
| B2 | 99.208% | 34.688% | 0.273% | 65.039% |
| B3 | 0.528% | 0.195% | 35.469% | 64.336% |

B2 and B3 reverse conditional branch preference, but absolute target arrivals are lower than B0. The paired B2−B0 A/initial change is −14.84 percentage points (95% CI −18.75 to −10.63); B3−B0 B/initial is −11.99 points (−15.47 to −8.36). Unassigned means no classified exit by the cutoff, not an outlet success or a numerical failure.

Across matched IDs, B0→B2 includes 838 B→A switches and 1,222 A→unassigned outcomes. B0→B3 includes 858 A→B switches and 1,169 B→unassigned outcomes. Thus rerouting and delayed/non-observed exit both contribute. Pooled transition counts are descriptive; inference uses seeds, not individual particles.

Controls change phase and amplitude together. Pure phase causality was not isolated. Near-wall hydrodynamics and clinical delivery remain unvalidated.

Exact evidence:

- [Paired statistics](../results/stage4/statistics/20261003T115215723091Z/SUMMARY.json)
- [Complete transition table](../results/stage4/statistics/20261003T115215723091Z/BRANCH_SWITCH_TABLE.csv)
- [Per-particle outcomes](../results/stage4/statistics/20261003T115215723091Z/PARTICLE_OUTCOMES.csv)
- [Finite ordered execution](../results/stage4/execution/20261003T100902329567Z/PROGRESS.json)

## Secondary execution and measured limitation

At this note, Stage 4I is running in `results/stage4/primary_secondary/20261003T115216806643Z`. It has not produced a completed 6 s secondary comparison. J (numerical sensitivity), K (selected history references), final audit and final report remain queued in strict order.

A closed independent 0.04 s secondary B2 diagnostic required 1,755 accepted substeps and 63,952 pair-contact events. It measured a maximum saved constrained-slip Reynolds number of 2.5774 at particle 64, t=0.04 s. Re is computed as 2ρR0|v−u_flow|/µ, using velocity after geometric constraints. There were 85 Re≥1 samples among 768 saved active particle snapshots. These are sample counts, not independent particle counts or time fractions; unsaved peaks can be larger.

This exceeds the small-Re assumption of the frozen Stokes transport model. No finite-Re drag correction or other new production physics has been added. Full secondary completion, if obtained, remains model extrapolation and is not physical validation. First-seed primary-only maximum saved constrained-slip Re is 0.066 or below.

- [Reynolds diagnostic](../results/stage4/performance/secondary_slip_reynolds/20261003T121702547297Z/SUMMARY.json)
- [Secondary profile](../results/stage4/performance/profile_B2_secondary/20261003T115538396917Z/SUMMARY.json)

A subsequent saved-state check of all 40 completed primary cases found a maximum constrained-slip Re of 0.08227 and no saved Re≥1 samples. This does not bound unsaved peaks or validate wall hydrodynamics. [All-seed Reynolds evidence](../results/stage4/performance/primary_slip_reynolds/20261003T122958202124Z/SUMMARY.json).

## Runtime investigation and continuation

The profile locates most cost in repeated wall-distance queries and contact projection; force summation accounts for about 2.3% of its instrumented runtime. A diagnostic-only interior-distance certificate reduced 114.8 s to 78.0 s. All 13 trajectory datasets, event bytes and validity summaries matched the reference. The patch was confined to that diagnostic process and was not deployed into production.

- [Certificate trial](../results/stage4/performance/certified_wall_trial/20261003T120728567547Z/SUMMARY.json)
- [Execution and inspection instructions](STAGE4_QUICKSTART.md)

The active supervisor job is `mevo_stage4_ordered_resume_20261003`; it has no automatic restart. Consult `scripts/stage4/status.py` and supervisor before submitting any new job. The existing ordered process continues the full protocol unless explicitly stopped or a real failure occurs. Do not duplicate it or treat the diagnostic 0.04 s run as the requested 6 s secondary result.

All 40 primary cases and their statistics have been synchronized locally. The Stage 3 tag and original evidence remain unchanged. No old BraVa/mouse job or Stage 5 study was restarted.

Note generated at 2026-10-03T12:25:52.101058+00:00.

## Subsequent runtime and numerical investigation

The local advance rate later fell enough that a linear extrapolation suggested 2–3 days for the original secondary calculation; this is a changing runtime estimate, not a completion guarantee.

A component-wise contact NNLS trial was **not accepted for production**. Its 3,519 shadow comparisons against the original solver at identical candidate states differed by at most 8.95e−16 m/s and relative L2 2.76e−15. Nevertheless, by the common t=0.04 s frame the maximum coordinate difference from the original trajectory was 49.1 µm, and contact/validity histories differed. This diagnoses individual-trajectory sensitivity in this secondary-ON calculation, not a measured change in final branch fractions. The original global contact solve is retained. [Exact diagnostic](../results/stage4/performance/component_contact_shadow/20261003T130544461050Z/SUMMARY.json).

At t=0.24 s, stable byte copies of the two still-incomplete secondary files contained 123/133 touching pairs with largest connected components of 13/11 particles. The source HDF5 files were never opened by the reader; only copies with unchanged source file statistics were opened and validated. These are partial diagnostic snapshots, not accepted results or restart checkpoints. [Snapshot provenance](../results/stage4/performance/secondary_contact_snapshot/20261003T130323707624Z/SUMMARY.json).

The whole-batch interior-wall certificate skipped no calls in 64-step replays from these later states because some particles were already close to the wall. A per-particle certificate instead calls the unchanged exact correction on each uncertain subset. Both 64-step replays preserve positions, velocities, contact events, timesteps and validity records bitwise. Measured runtime changed from 13.07 to 6.07 s (B2) and 14.99 to 6.45 s (B3). Four focused tests also cover interior no-ops, actual wall corrections, exact-cache preservation and fallback behavior. [Later equivalence evidence](../results/stage4/performance/later_particle_geometry_certificate/20261003T131349331620Z/SUMMARY.json).

The continuation runner now supports `--secondary-certified-wall`. Each enabled case records the helper hash and equivalence evidence in its input identity. Enabling this option keeps the original physical equations, global contact NNLS, integration parameters and validity guards. It requires a new run from the original initial conditions; a partial HDF5 file is not used as a restart checkpoint.

## Current continuation after the exact optimization

At 13:22:58 UTC, the original ordered job was intentionally stopped with both secondary cases still incomplete (last progress B2 0.29773 s, B3 0.31706 s). This is an operator-initiated performance restart, not a numerical failure. Its partial files and [interruption record](../results/stage4/primary_secondary/20261003T115216806643Z/INTERRUPTED_FOR_EXACT_OPTIMIZATION.json) are preserved and copied locally.

The active supervisor job is now **`mevo_stage4_ordered_certified_wall_20261003`**, with [execution identity](../results/stage4/execution/20261003T132342434994Z/IDENTITY.json) and secondary batch `results/stage4/primary_secondary/20261003T132343245938Z`. It restarts both B2/B3 from the same original particle hash, `42679055d151e4adb067d7efcc97643b205f8802e85dd2dc063fb6ea50a16125`, and then runs J, K, final audit and final report in order. All 42 production source files remain unchanged. Enabled helper SHA-256: `c439d5e17deecc0a66cedf3221611830fb774189d104fbc69579d180364a477f`. Related geometry/contact regression: 22 tests passed.

An initial comparison of the actual restarted runs through the common saved time 0.06 s found all 13 HDF5 datasets bitwise identical to the original runs in both controls. This is a prefix check, not proof that the full 6 s runs have completed. [Prefix evidence](../results/stage4/performance/certificate_prefix_check/20261003T132810756410Z/SUMMARY.json).

The older 0.24 s copied frames measured maximum saved constrained-slip Re of **11.0602 (B2)** and **6.62185 (B3)**, further outside the Stokes small-Re regime. The contact-roundoff diagnostic also distinguishes microtrajectory changes from ensemble statistics: at 0.04 s, mean NN changes from 33.278 to 33.313 µm and the largest cluster fraction remains 0.046875, while contact fraction changes from 0.628906 to 0.566406. Final secondary branch sensitivity remains untested. [Numerical/model diagnostic](../results/stage4/performance/secondary_numerical_diagnostic/20261003T133017582804Z/SUMMARY.json).

## Read-only follow-up at 17:26 UTC / 19:26 Europe/Zurich

The complete saved prefixes of the intentionally interrupted original runs have now been compared: 15 B2 frames through 0.28 s and 16 B3 frames through 0.30 s. All 13 datasets match the current certified-wall runs bitwise. The newer source files were again read only through stable byte copies; these partial copies are neither accepted full runs nor restart checkpoints. [Full old-prefix comparison](../results/stage4/performance/certificate_prefix_check/20261003T165357948621Z/SUMMARY.json).

At this observation, B2 had reached 1.054972 s in 133,652 accepted steps and B3 had reached 1.272341 s in 142,799 steps. Both still had 256 active particles and no recorded worker failure. The two worker processes each reported approximately 100% CPU, with resident memory about 1.3 GB. Progress files continued to refresh.

The recent accepted substeps averaged approximately 7 microseconds, much smaller than the unchanged 0.5 ms nominal step because the displacement guard subdivides the integration. This and geometric contact work explain the long runtime; saved RUNNING status alone would not establish process health. The full 6 s protocol continues under the same supervisor job. J, K, final audit and final report remain pending; no completed secondary branch result or final Stage 4 acceptance is claimed here.

## Follow-up at 2026-10-04 00:09 UTC / 02:09 Europe/Zurich

The same certified-wall secondary batch remains incomplete. At this observation, the live progress was:

| Case | Simulated time [s] / 6 s | Active particles | Confirmed first exits |
| --- | --- | --- | --- |
| B2 | 2.1842 | 206 | 50 through A |
| B3 | 2.5628 | 199 | 57 through B |

These are intermediate counts, not final branch fractions. The append-only event logs were read without opening the active HDF5 files. Within each case, exit IDs are unique; the recorded positions lie on their finite classification planes within tolerance and the recorded velocities point outward. The cumulative exit counts agree with the reduction in active particles. The latest B2 batch contains 21 exits between 2.114007 and 2.115468 s. During that batch, event output briefly preceded the periodically written progress count; rereading after the progress update established agreement. No worker failure has been recorded. J, K, final regression and the final report remain pending.

An 81.15 s read-only resource sample ending at 01:50 Europe/Zurich measured 0.9998/0.9999 CPU cores for the B2/B3 workers, with negligible scheduler wait and resident memory 1.447/1.527 GiB. Container memory was 32.20 GiB against a 42.06 GiB limit. A subsequent placement check found the two workers on distinct physical cores. No quota, priority, affinity or unrelated job was changed. This sample does not indicate a resource stall or establish a fixed completion time.

The final identity audit now explicitly includes the four original formal primary-only cases, in addition to the paired-seed, secondary and convergence runs (commit `1b5ac8e`). A read-only check of the four closed local primary backups verified their initial hashes, identical 42-file production source hashes, particle bookkeeping, complete trajectories and 6 s endpoints. The updated server audit script passes Ruff; the complete final server regression is still queued. This change only expands audit coverage and does not modify the running integration or the certified-wall helper.
