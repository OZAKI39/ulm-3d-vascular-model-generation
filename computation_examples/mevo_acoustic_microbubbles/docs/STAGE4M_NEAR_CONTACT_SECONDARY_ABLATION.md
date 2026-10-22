# Stage 4M — Near-Contact Secondary Bjerknes Ablation Sensitivity

MODEL SENSITIVITY / ABLATION STUDY. This is not a corrected near-field force model.

> Sensitivity-only smooth removal of unresolved classical secondary interaction inside the lubrication zone.

The final Stage 4L lubrication, contact solver, tolerances, gap-aware timestep and geometric projection are frozen. Acoustic fields, PZT controls, adaptive APECSS response, primary force, flow, geometry, classifier, outlet semantics and particle initial population remain unchanged. Only classical secondary pair weights change.

## Policies

- L0: existing legacy Stokes/full-secondary gold; read-only.
- A0: existing Stage 4L lubrication/full-secondary trajectories; read-only. The default force path preserves the Stage 4L operation order bitwise.
- A1: weight zero for h ≤ 0.001R, otherwise one. Limited to 5 ms dense diagnostics.
- A2: h_floor=0.001R, h_act=0.1R; s=clip((h−h_floor)/(h_act−h_floor),0,1); w=6s⁵−15s⁴+10s³. With R=1.78 µm, the bounds are 1.78 nm and 178 nm. One reciprocal weight multiplies each pair before summation. No particle-level scaling or extra force/velocity clipping.

Classical and effective diagnostics sum magnitudes once over unordered pairs; they are not net force. Their ratio is null when the classical sum is zero. Counts distinguish weight one, transition and weight zero.

## Causality and historical provenance

L0→A0 isolates the translation change. A0→A2 is the incremental ablation comparison; L0→A2 is the combined model difference. The existing 0.08 s HDF5 frames are compared at equal time; final accepted A0 positions at approximately 0.089/0.097 s are not substituted.

The existing dense A0 runs predate final Stage 4L contact QR fallback repairs (same resistance, convex contact objective, constraints and KKT tolerance). Original-initial A0 runs predate final diagnostic additions and explicit nonfinite rejection, with the finite successful solve unchanged. `existing_reference/BACKEND_IDENTITY.json` and the saved source diff document these differences. New experiments use the final frozen Stage 4L backend. Historical wall-time ratios are not an identical-source or isolated-hardware benchmark; accepted-step ratios are also reported. No A0 replay simulation is run.

## Predeclared continuation gates

M0 uses the exact sealed legacy t=3 s dense checkpoint, not a newly formed A2 cluster. A1 and A2 each run B2/B3 for 5 ms, with a 7200 s diagnostic wall budget and 100000 accepted-step safety limit per case. Only A2 determines continuation: each case needs ≥5× accepted-step reduction OR ≥5× wall speedup over the existing complete A0 5 ms run. Completion, finite states, positive steps, converged solves and zero multi-radius projection events are mandatory. A gross expansion guard flags any same-ID initial largest-cluster gyration radius above twice its initial value; this is an engineering pathology guard, not a physical accuracy criterion. The full Rg trace and pair-distance quantiles remain available for inspection. A1 never enters a longer campaign.

If either A2 case fails, final status is `A2_RUNTIME_GATE_FAILED`; no cutoff tuning, new physics or long simulation follows.

If M0 passes, A2 starts from the original population and targets 0.25 s directly with unchanged output semantics and a 7200 s/case budget. Runtime exhaustion is `INCOMPLETE_RUNTIME`. Only the joint endpoint criteria largest≤2, near-contact fraction<0.05 and cluster count≥0.95N define `AGGREGATION_ESSENTIALLY_ABSENT`. Otherwise a stable completed case proceeds to 1 s. The current runner has no exact full solver resume; longer runs restart the identical original population with explicit provenance, never a visualization frame. After 1 s, the measured 0.5–1 s wall/simulated-time rate projects 6 s cost. More than 12 h/case is `FULL_RUN_RUNTIME_GATE_NOT_PASSED`; otherwise stable aggregation permits the full branch-outcome experiment. No primary-only or secondary dt/2 campaign is added.

## Literature scope

[Doinikov & Bouakaz (2015), Physical Review E 92, 043001](https://doi.org/10.1103/PhysRevE.92.043001), “Theoretical model for coupled radial and translational motion of two bubbles at arbitrary separation distances,” motivates the limitations of large-separation coupling: high-order inverse-distance terms, through approximately d⁻⁹ in their close-separation analysis, can matter and attraction/repulsion can change. This does not validate the current shell-coated dense population.

[Chen, Chen & Li (2020), Ultrasonics Sonochemistry 61, 104814](https://doi.org/10.1016/j.ultsonch.2019.104814), “Modulation of the secondary Bjerknes force in multi-bubble systems,” shows that surrounding bubbles alter radial responses and pair forces. Its numerical ratios are not correction factors for this project's frequency, pressure or shell model.

Neither paper implies the true secondary force vanishes below 0.1R. A2 equilibrium gaps cannot be interpreted as physical SonoVue spacing, shell contact, fusion or clinical packing. A2 is not promoted as validated near-field physics.

## Results

The existing 0.08 s comparison is complete: B2 largest cluster 19→19, B3 19→20 (L0→A0); near-contact fractions 0.8125→0.84375 and 0.828125→0.87890625. Lubrication has changed some early metrics but has not removed early aggregation. This single shared time is not an estimate of a continuous formation rate. Detailed distributions and Figure M1 are under `results/stage4/secondary_ablation/l0_vs_a0/`.

M0 is complete; the failed engineering gate prevented downstream formation, growth and full branch experiments.

## Final measured outcome

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
