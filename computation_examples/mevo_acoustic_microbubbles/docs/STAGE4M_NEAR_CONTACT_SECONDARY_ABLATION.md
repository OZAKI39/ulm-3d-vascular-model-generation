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

M0 and downstream results will be populated only after their respective gates.
