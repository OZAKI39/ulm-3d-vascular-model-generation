# Assigned-particle deactivation audit — Stage 4I.5

Status: **UNSAFE to introduce or justify as a physically proven terminal-outlet optimization.** The gold implementation already deactivates first-exit particles. Stage 4I.5 preserves that existing absorbing-boundary convention for numerical regression and does not claim an additional speedup from it.

## Evidence from the implementation and completed results

`BifurcationExperimentTracker._run` in `src/mevo_mb/coupling/bifurcation_experiment.py` tests classification crossings after every accepted step. It freezes the recorded coordinate at the interpolated first crossing, then filters `active`, positions, and every evaluated state array using the same `keep` mask. All subsequent calls to `evaluate`, contact projection, and the minimum-step calculation receive only the surviving population. Thus assigned particles no longer participate in response lookup, primary/secondary forces, flow/wall queries, contact search, timestep limiting, or integration.

The HDF5 recorder retains all initial IDs for bookkeeping. Inactive coordinates stay at the crossing location; their force/velocity fields are deliberately NaN. This is an output convention, not continuing downstream dynamics.

The completed gold checks verified unique first-exit IDs, frozen final crossing positions, monotonically decreasing active masks, and complete population accounting:

| Case | A exits | B exits | Active at 6 s / observation cutoff | Wall loss |
|---|---:|---:|---:|---:|
| B2 | 76 | 0 | 180 | 0 |
| B3 | 0 | 88 | 168 | 0 |

The source result IDs and immutable checks are in `results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json`, SHA256 `3d768c763a3f20625b8064039fea5cddd35cc2765810a3511b46eef7a288cef8`. The actual run source commit is `bb6261902dbdc4a4c3a6121144e1cad52181b3f5`.

## Why physical terminal safety is not established

The daughter classification planes are 10 mm from the bifurcation; the CFD daughter outlets are at 12 mm. Positive outward velocity at first crossing does not prove that a particle could never return under subsequent acoustic/contact dynamics.

`ClassicalFarFieldBjerknes` evaluates all pairs supplied to it. It does not apply the legacy `secondary_force_cutoff_m` configuration field. No finite interaction cutoff or distance/error bound proves that a just-classified particle can no longer affect a particle immediately upstream of the classification plane. Deactivation removes that interaction by the existing model's population-domain convention, not by a demonstrated physical negligibility result.

The fixed concentration ROI lies 1–3 mm into each daughter, but this alone does not establish terminal safety or validate downstream cluster histories. Gold cluster/contact metrics describe the remaining active population, not the future evolution of all initialized particles.

## Consequences for performance and termination

- Do not add an assigned-particle deletion optimization: it is already present.
- Keep the gold convention unchanged in P0/P1/P2 regression. Restoring downstream particles or moving classification planes would be a separate physical/boundary sensitivity, not an equivalent performance refactor.
- Introduce an explicit objective-based termination policy. Neither current gold case qualifies for all-terminal-assigned early termination: 180 and 168 particles remain censored at 6 s.
- A plateau in cumulative exits does not prove that future branch assignments are determined. Do not shorten a comparison's observation window merely because recent exit counts stopped changing.
- Final Stage 4 interpretation must disclose the absorbing classification-plane convention and the absence of a terminal-interaction safety proof.
