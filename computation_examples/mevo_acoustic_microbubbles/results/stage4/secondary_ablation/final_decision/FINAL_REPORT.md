# Stage 4M final decision

**A2_RUNTIME_GATE_FAILED** — MODEL SENSITIVITY / ABLATION STUDY.

The four requested 5 ms dense diagnostics and the existing-data comparison are evaluated below. No 0.25 s, 1 s or 6 s A2 campaign was launched. L0/A0 and primary-only trajectories were not rerun.

| Case | Policy | Status | Steps | Wall time (s) | Step reduction | Wall speedup |
|---|---|---|---:|---:|---:|---:|
| B2 | A0 | COMPLETE | 17385 | 2216.0 | 1.000× | 1.000× |
| B2 | A1 | COMPLETE | 11960 | 1382.9 | 1.454× | 1.602× |
| B2 | A2 | COMPLETE | 10712 | 1262.6 | 1.623× | 1.755× |
| B3 | A0 | COMPLETE | 18627 | 3080.3 | 1.000× | 1.000× |
| B3 | A1 | COMPLETE | 17279 | 2372.4 | 1.078× | 1.298× |
| B3 | A2 | COMPLETE | 12415 | 1855.2 | 1.500× | 1.660× |

Both A2 cases must individually achieve ≥5× by either measure and pass numerical/envelope checks. This is an engineering gate, not an accuracy threshold. Historical A0 used an earlier contact-solver implementation with the same physical resistance, constraints and tolerance; runtime ratios therefore include a backend-revision and machine-load caveat.

## Q1: EXISTING_EQUAL_TIME_COMPARISON_COMPLETE

At 0.08 s, lubrication retains strong early aggregation. B2: largest cluster 19 to 19; near-contact fraction 0.812500 to 0.843750. B3: largest cluster 19 to 20; near-contact fraction 0.828125 to 0.878906. Some formation metrics change, but a single checkpoint cannot establish a continuous formation-rate change.

## Q2: DENSE_DIAGNOSTIC_COMPLETE

A1 runtime and accepted-step changes are listed per case; it is a hard-switch diagnostic only, with no longer runs.

B2: 1.454× fewer accepted steps, 1.602× wall speedup.
B3: 1.078× fewer accepted steps, 1.298× wall speedup.

## Q3: A2_RUNTIME_GATE_FAILED

A2 did not meet the required fivefold continuation gate. Reduced pair-norm sums do not remove the observed small-step bottleneck sufficiently for a long campaign.

## Q4: NOT_RUN_M0_GATE_FAILED

Whether A2 forms clusters from the original initial population remains unanswered. Legacy-formed dense checkpoints cannot answer formation.

## Q5: NOT_RUN_M0_GATE_FAILED

A2 growth to large clusters by 1 s remains unanswered.

## Q6: NOT_RUN_M0_GATE_FAILED

A2 B2→A and B3→B branch-steering robustness remains unanswered; no full 6 s trajectories were launched.

## Dense-envelope and force-removal observations

B2/A1: same-ID Rg 7.095→6.642 µm; maximum Rg/initial 1.0000; endpoint rigid-fit residual 1.451 µm. Effective/classical pair-norm ratio 0.8448→0.6331.
B2/A2: same-ID Rg 7.095→6.541 µm; maximum Rg/initial 1.0000; endpoint rigid-fit residual 1.827 µm. Effective/classical pair-norm ratio 0.7643→0.6003.
B3/A1: same-ID Rg 7.584→7.062 µm; maximum Rg/initial 1.0000; endpoint rigid-fit residual 1.520 µm. Effective/classical pair-norm ratio 0.8593→0.6815.
B3/A2: same-ID Rg 7.584→7.007 µm; maximum Rg/initial 1.0000; endpoint rigid-fit residual 1.725 µm. Effective/classical pair-norm ratio 0.7907→0.6508.

Native accepted-state graph turnover, censored edge lifetimes, unordered pair-distance quantiles, subfloor counts, contact QP statistics and projection diagnostics are retained under `comparisons/dense/` and each new run. Graphs and positions are never interpolated. Same-ID clusters were formed by legacy L0 evolution; preservation or rearrangement is not evidence of A2 formation.

## Scope and reproducibility

The frozen Stage 4L translation backend hash is `0b3e36eae784f064986caf371a201f1ded0be5e6970eef6f1f856c6dd89c9725`. The only production changes are secondary pair weighting and config forwarding. Each run contains source and input identities, all accepted states and force diagnostics. Two B3 startup attempts failed before integration because a shared compiled cache was read during compilation; each was restarted with an independent cache, using unchanged contact source/compiler flags. The generic smooth-removal banner archived by the first runner applies to A2; A1 is explicitly a hard subfloor switch, as clarified by its POLICY_SCOPE.json.

Figures M1, M2, M3 and M5 are provided in PNG/PDF. M4, M6 and M7 are unavailable because their experiment gates did not open; no synthetic replacement data were created.

[Doinikov & Bouakaz (2015)](https://doi.org/10.1103/PhysRevE.92.043001) and [Chen, Chen & Li (2020)](https://doi.org/10.1016/j.ultsonch.2019.104814) motivate unresolved close-separation and multi-bubble coupling. Neither paper validates a zero near-contact force or provides correction factors for this study.

All dense runs retain the inherited PHYSICS_VALIDITY_WARNING, including far-field separation, acoustic-retardation and near-wall scope limits. Software-test and numerical-stability passes do not validate these physical approximations.

The read-only A0 pair diagnostics cover every saved accepted state. They were reconstructed with NumPy 2.5.3 on the local analysis environment; initial-force relative L2 differences from the archived server arrays are about 1.6–1.7×10⁻¹⁶. These roundoff-scale cross-environment differences are distinct from the same-environment bitwise A0 unit regression. No A0 translation step was taken.
