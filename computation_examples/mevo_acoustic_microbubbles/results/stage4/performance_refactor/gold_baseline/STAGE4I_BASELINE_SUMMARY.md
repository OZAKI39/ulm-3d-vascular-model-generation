# Stage 4I gold baseline

Both original cases completed naturally before the orchestration-only boundary stop. J was never started.

| Case | A | B | Observation cutoff | Accepted steps | Runtime (h) | Mean dt (us) | Saved max constrained Re |
|---|---:|---:|---:|---:|---:|---:|---:|
| B2 | 76 | 0 | 180 | 902991 | 33.462174 | 6.644584 | 11.060233 |
| B3 | 0 | 88 | 168 | 925123 | 34.360901 | 6.485624 | 6.899177 |

All 256 IDs are accounted for per case. Both trajectories span 0–6 s and contain 301 saved states. Active fields are finite; inactive force/velocity NaNs are intentional.

The run source commit is `bb6261902dbdc4a4c3a6121144e1cad52181b3f5`; the repository commit at sealing is recorded separately. Each archived source file was checked against its original source manifest.

Full-step maximum dt, rejected-step counter, wall-projection call counter and GPU elapsed time were not recorded. The manifest distinguishes missing measurements, code-derived counters, saved-frame diagnostics and process samples. The minimum dt can be a macro-step remainder and is not by itself a physical stiffness measure.

Contact-event counts include repeated geometric corrections across steps, not distinct colliding pairs or fusion events. Geometry cache misses count point queries, not projection calls. Timing components overlap and must not be summed as independent costs.

The complete source results, including raw event logs, are preserved as independent byte copies. Input snapshots and SHA256 checks prevent accidental reuse of changed geometry, fields or response tables. No 6 s baseline rerun is permitted.

Gold status does not remove model warnings. Small-Re Stokes translation and classical far-field secondary forces are not validated in the observed finite-Re/near-contact region.
