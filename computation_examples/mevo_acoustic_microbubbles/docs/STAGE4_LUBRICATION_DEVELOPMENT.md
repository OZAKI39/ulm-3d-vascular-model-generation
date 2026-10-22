# Stage 4L — normal lubrication translation development

## Current status

Final scoped decision: **`LUBRICATION_TRANSLATION_DEVELOPED_WITH_SCOPE`**. The backend, paired 1 ms/5 ms probes, new 6 s primary B0–B3 baseline, and new B2/B3 nominal/dt-half comparison are complete. **60 focused tests pass**, and changed Python modules pass Ruff. This is an engineering development result for branch-control stability with explicit timing and contact limitations, not an overall Stage 4 PASS or validated full near-contact bubble physics.

Both S1 retries stopped at their declared 3,600 s wall budget before 0.1 s: B2 reached **0.089448770 s**, B3 **0.096525398 s**. They remain **INCOMPLETE**, with no numerical solver failure and zero geometric corrections in those retries. No 0.5–1 s secondary extension, full 6 s secondary, multiseed, old K or final audit was launched. The full-secondary runtime gate is **NOT_PASSED**.

Final gate and secondary termination evidence: `results/stage4/lubrication_upgrade/comparison_legacy/final_stage4l_20261005/`.

The final `ARTIFACT_INDEX.json` records SHA-256 values for the curated reports, figures, inputs and source snapshots. Sixteen selected raw HDF5 files (about 205 MB) are preserved locally and on the server outside Git; verbose per-step logs remain on the server. The Git archive contains about 14.4 MB of selected evidence, rather than the raw trajectory/log collection. No original data was deleted.


## Isolation, source and immutable inputs

- Local project: `/home/lzy/projects/computation_examples/_archives/mevo_stage4l_20261005/computation_examples/mevo_acoustic_microbubbles`, branch `dev/stage4l-lubrication-20261005`.
- Baseline snapshot: `aaf2549fc506933df74229be8e0ec9717f97f7ed`, branch `sync/mevo-acoustic-microbubbles-stage4-20261005`. All 278 original tracked source files matched this snapshot byte for byte before development.
- Server project: `/workspace/mevo_acoustic_microbubbles_stage4l`, branch `stage4l-lubrication-20261005`, starting from the byte-identical standalone source baseline `a6c3153f484113d37ffb5e9b4b57f0c6594ef396`. Its Git ancestry differs from the local monorepo and is not misrepresented.
- Final backend code at this update: local `921e989`; server `06af1f5`. B2/B3 nominal and dt/2 scientific reruns use source commit `12331bb`. The later code adds exposed force components, split diagnostics and JSON-safe rejection reporting, without changing successful velocity/contact computations. Every run retains its actual source archive and hashes.
- Original files, gold inputs, configurations and HDF5 files remain untouched. The original legacy controller remains paused; its already-running Stage J child was observed but not signalled. No old BraVa batch or stopped mouse task was resumed.

All results of this task live below `results/stage4/lubrication_upgrade/`. Existing results are interpreted as `LEGACY_STOKES_TRANSLATION_BASELINE`; they are not relabelled as new lubrication physics. Each new experiment has its own exclusive identity, initial population/snapshot, input SHA-256 values and source snapshot. Completed selected runs passed post-run input-hash checks.

## Model and solver evidence

See [LUBRICATION_TRANSLATION_MODEL.md](LUBRICATION_TRANSLATION_MODEL.md) for equations, assumptions, wall-normal convention, activation subtraction, continuum floor, exact wall geometry, timestep and contact semantics, references and reproduction commands.

The implementation retains Stokes self-drag and adds normal particle-particle and particle-wall resistance. It is always on when the configuration selects `lubrication_coupled`, including primary-only and flow-only experiments. R=1.78 µm, activation=0.1R and floor=0.001R. Unequal radii are explicitly rejected. The floor saturates a resistance coefficient; it is not adhesion, a target gap or a cluster lock.

Sparse CSC/LU mobility is checked against 1e-11 relative residual. Contact uses the coupled A metric, with a 1e-9 normalized feasibility/complementarity threshold. Nearly dependent contact rows require primal QR refinement; a feasible active-set fallback resolves cycling working sets. These refinements solve the same constrained optimization with unchanged forces, tolerances and resistance law. Finite-horizon normal constraints protect RK proposals; geometric projection remains a residual correction and is measured at every stage.

Focused checks cover T1–T10, the pair relative-speed factor of two, wall asymptotics, far-field regression, symmetry/SPD and CG agreement, predictive gap guards, floor behavior, unequal-radius rejection, exact VTK cache reuse, always-on lubrication, A-metric contact and unchanged legacy tracking. Four small stored fixtures reproduce actual rejected contact QPs. Further tests check force-component balance, strict-JSON failure evidence and cluster/graph statistics. Current evidence: `unit_tests/final_failures_v8.xml` (60 passed), plus dimensional CSV/JSON tables under `analytic_pair/` and `analytic_wall/`.

## Paired dense-checkpoint probes

Both arms start from the same sealed t=3 s state: 186/189 active B2/B3 particles; initial largest near-neighbor clusters contain 75/91 particles. Initial positions, IDs, radii, masks, field/response/flow samples and primary/secondary force arrays match bitwise. No gold trajectory was overwritten or presented as a new continuation. These mature clusters were formed by legacy translation; the probes do not show that the new model would form the same clusters by t=3 s from the original initial condition.

| Probe | Case | Legacy/new accepted steps | New wall time (s) | Rigid-fit residual, new/legacy | Mean near-neighbor turnover, new/legacy | New maximum projection (m) |
|---|---|---:|---:|---:|---:|---:|
| 1 ms | B2 | 158 / 3,626 | 400.37 | 0.16295 | 0.00917 | 0 |
| 1 ms | B3 | 174 / 3,684 | 491.86 | 0.14110 | 0.00902 | 0 |
| 5 ms | B2 | 791 / 17,385 | 2,216.01 | 0.16597 | 0.001864 | 5.54e-12 |
| 5 ms | B3 | 872 / 18,627 | 3,080.30 | 0.14362 | 0.001806 | 0 |

These ratios use the fixed initial largest-cluster IDs and a 10 µs physical clock. The saved states are actual last-accepted positions, with sample ages recorded; graphs/envelopes do not interpolate coordinates. A separate 2 µs analysis is retained. The earlier interpolated-statistics analysis under `paired_1ms_v3/` is exploratory and superseded for graph/envelope interpretation; its original launch gate remains preserved.

For all initial clusters of at least three particles in the 1 ms probe, particle-weighted rigid-fit RMS falls from 7.02 to 1.18 µm (B2) and 7.49 to 1.19 µm (B3). Some initially quiet small clusters have increased internal displacement, so improvement is not universal per cluster. Largest-cluster gyration radii compact by about 6.8%/6.6% and then stabilize. This is measured behavior, not an imposed rigid aggregate.

Legacy 5 ms maximum projections are about 10 µm, with 28,418/38,281 particle corrections above two radii. The new arms have none above two radii. Median new steps are 0.285/0.254 µs, only about 4.45%/4.36% of legacy medians. The specified closing-gap rule therefore does not improve runtime. Holding the measured dense-state cost over all 6 s would take about **31/43 days per case**. This is a cost scenario, not a prediction of the entire evolving trajectory, and it fails the full-secondary launch gate.

Authoritative comparisons, figures and gate: `comparison_legacy/paired_1ms_actual_states_v5/` and `comparison_legacy/paired_5ms_actual_states_v5/`. The short-probe gate is `PASS_WITH_SCOPE`; it is not long-time physical validation.

### Native-step hard-contact limitation

Every accepted state was also checked without fixed-clock sampling. Near-neighbor edge-change rate ratios remain 0.001996/0.001673, confirming suppression of neighbor replacement. However, strict hard-contact membership at h <= 10 pm changes faster: rate ratios about 4.01/4.57, and observed hard-contact lifetimes are shorter. Near-contact edge lifetime improvements must not be described as an improvement of every contact metric.

At those hard-contact switches, median normal gap changes are 0.044/0.057 nm, 99th percentiles 0.848/0.881 nm, and maxima 1.76/2.22 nm. Legacy medians were 4.45/4.43 µm. This is residual subnanometre contact/gap chatter, not merely floating-point noise or validated physical contact kinetics. No hysteresis, filtering or tolerance change was applied. All reported lifetimes are censored by the 5 ms window, not material persistence estimates.

Evidence: `comparison_legacy/paired_5ms_actual_states_v5/NATIVE_GRAPH_REVIEW.json` and each probe's `NATIVE_GRAPH_METRICS.json`. No tangential-lubrication conclusion is inferred solely from these counts.

## New primary-only baseline

All four new primary cases completed 6 s from the original initial population. The table compares each new run against its legacy counterpart. A/B/unassigned count all 256 particles; unassigned means observation cutoff, not successful outflow. There are no wall-loss assignments. All 256 final particle identities agree in each comparison, but crossing times change.

| Case | A / B / unassigned | ID agreement vs legacy | Maximum crossing-time change (ms) | Accepted steps | Maximum all-stage projection / R |
|---|---:|---:|---:|---:|---:|
| B0 | 128 / 123 / 5 | 100% | 5.329e-12 | 12,000 | 0 |
| B1 | 37 / 28 / 191 | 100% | 4.893 | 12,056 | 0.0034462 |
| B2 | 90 / 0 / 166 | 100% | 97.11 | 17,110 | 0.055533 |
| B3 | 0 / 98 / 158 | 100% | 14.11 | 16,703 | 0.060537 |

B2/B3 steering is retained: assigned B2 particles exit A, and assigned B3 particles exit B. `f_A` and `f_B` in the machine report condition on assigned exits; separate initial-population fractions include all unassigned IDs. Short dense probes having negligible projection must not be generalized to every primary stage: B2/B3 maxima are about 0.056R/0.061R, with RMS corrections about 0.716/0.680 nm and no multi-radius corrections.

Evidence: `comparison_legacy/primary_baseline_v7/PRIMARY_COMPARISON.json`.

## New timestep sensitivity

B2 nominal 0.5 ms and half 0.25 ms both completed. All 256 branch identities agree. Of 90 shared crossings, mean absolute difference is 0.662 ms, median 0.0404 ms, 95th percentile 1.44 ms, and maximum 31.10 ms (particle 188, about 0.948% of its nominal crossing time). Early 0–0.1 s cluster count, near-contact fraction and largest-cluster fraction agree at common recorded times; nearest-neighbor relative L2 difference is 1.21e-5.

This supports the branch assignment result while leaving a measurable timing sensitivity. It does not certify pointwise trajectory convergence or dense-secondary subfloor contact motion. B3 nominal and dt/2 also completed: all 256 branch identities agree; crossing-time differences have mean 0.366 ms, median 0.0244 ms, 95th percentile 1.07 ms and maximum 12.47 ms (particle 25, 0.487% of its nominal time). Its early nearest-neighbor relative L2 difference is 5.47e-6, with the other listed early cluster metrics equal. No legacy dt/2 result is inherited.

The engineering acceptance is limited to branch-control stability with the measured timing bounds. It does not certify every crossing to sub-millisecond accuracy, pointwise trajectory convergence or dense-secondary gap dynamics. The review does not introduce a previously declared accuracy tolerance; observed timing differences remain explicit limitations. Authoritative combined comparison: `comparison_legacy/20261005T232125840331Z_branches/SUMMARY.json`.

## Secondary staging and finite runtime limits

The initial B2 0.1 s attempt exposed a contact solver failure; the initial B3 attempt reached its 1,800 s budget around t=0.074 s. Both remain incomplete. Repaired-solver retries used new identities and the original initial condition, with a declared 3,600 s per-case budget.

| Case | Accepted time (s) | Target (s) | Accepted steps | Median dt (µs) | Maximum projection (m) | Result |
|---|---:|---:|---:|---:|---:|---|
| B2 | 0.089448770 | 0.1 | 16,907 | 0.579 | 0 | INCOMPLETE — runtime budget |
| B3 | 0.096525398 | 0.1 | 30,955 | 1.010 | 0 | INCOMPLETE — runtime budget |

Both retained all 256 active particles. Source input postchecks passed. Maximum linear residuals were 1.36e-14/1.73e-14; contact residuals remained below the unchanged 1e-9 threshold. B2's accumulated contact solves consumed about 2,535 s, compared with 1,640 s for B3. No solver failure or nonpositive timestep occurred in these retries. Stopping at a finite resource limit is not physical completion.

The ordinary visualization HDF5 files end at their last saved 0.08 s frame and have `complete=false`; each final accepted state is separately preserved in `particles/last_accepted_state.npz`. Do not claim the HDF5 contains frames through the checkpoint time. S1 has not passed the requested 0.1 s completion gate, so S2 was not launched. Measured dense-state step behavior and cost also rule out the S3 launch. The task ends at this explicit scoped gate, with no new simulation left running.

The final secondary report is `comparison_legacy/final_stage4l_20261005/SECONDARY_STAGE_REPORT.json`. `secondary_full/NOT_LAUNCHED.json` records staging only, not a successful secondary run. No full secondary-upgrade status is assigned.

## Retained failures and numerical repairs

- The first dense attempts failed because a logger supplied `time_s` twice. New result IDs were created after repair; failed artifacts were not overwritten.
- Some original primary identities predate `effective_acoustic_field`. The loader now reads their historical PZT field path and validates the recorded field SHA-256. Their failed campaign rows remain recorded.
- Longer primary/secondary runs exposed rank-deficient contact rows. Four rejected steps were exactly replayed, retained as fixtures and solved after QR refinements. The original failures remain incomplete; reruns use new IDs and frozen initial conditions.
- An early validation command pointed at the wrong scratch import root. Its collection-error report is retained separately; final tests used the actual isolated project source and passed.
- Runtime-limited partial secondary results are not numerical or scientific completion. The earlier primary dt campaign was allowed to finish rather than being killed; final comparisons select the explicit new nominal/half solver revision.

## Remaining physical scope

Classical secondary Bjerknes remains unchanged despite its known near-contact far-field limitation. Floor/contact occupancy remains substantial. The dense-probe `saved_Reynolds_max` uses total speed, not slip speed, and starts above creeping-flow conditions. Small residuals, same branch counts and smooth envelopes do not establish local physical accuracy.

No tangential/rotational lubrication, friction, shell repulsion, adhesion, fusion, coalescence, deformable-contact/film PDE, history or added-mass production dynamics was added. Rigid-sphere resistance literature is not SonoVue calibration. No artificial damping, speed cap, smoothing or cluster locking was introduced. Full Stage 4 acceptance remains outside this record.

## Recorded-data visualization — 2026-10-06

The Stage 4L renderer reuses the existing `scripts/stage4/render_single_overview.py` v7 scene: black background, one enlarged orthographic view, shaded blue spheres, the same PZT discs and purple pressure slices, upper-left clock, and largest-cluster population label. Camera, field range, physical clustering threshold and 60-fold visibility magnification are unchanged. Only experiment-specific titles, force labels, time precision and a concise scope annotation are added. There is no second view, blue aggregation ring or explanatory footer.

Following the user's title correction, B2/B3 case identifiers are omitted from the video image. Filenames and the viewer retain them to distinguish the datasets.

Four separate 1920 × 1080, 60 fps videos are generated below `results/stage4/lubrication_upgrade/visualization/B2_B3_clean_enlarged_v8_20261006/`:

- B2/B3 dense probes: recorded absolute t=3.000–3.005 s, 1,001 frames, about 16.7 s per video. Every frame takes the actual last-accepted state at a uniform 5 µs display clock. Maximum sampling ages are about 0.44/0.51 µs. No coordinate interpolation or raster blending is used for these dense probes.
- B2/B3 primary-only runs: t=0–6 s, 1,201 frames, about 20 s per video at constant 0.30× playback. All 301 saved states are used; the established v7 raster cross-dissolves supply intermediate display frames without interpolating particle coordinates. Population labels refer to the latest saved state.

These experiments are not spliced into a continuous trajectory. The dense initial clusters were formed by the legacy model; the new short probe does not establish their formation under lubrication. The primary videos contain no secondary force. The incomplete S1 retries and absent full-secondary runs are not represented as completed trajectories.

The renderer verifies raw trajectory hashes against the final Stage 4L artifact index and resolves the original acoustic field and wall mesh by SHA-256. It retains per-frame source/time mappings and per-state active particle identities and connected components. `scripts/stage4l/check_motion.py` checks full video decoding, exact 1/60 s presentation intervals, constant display-time increments, encoded keyframe fidelity, active-ID partitions and unchanged inputs. Any physically stationary or subpixel motion remains visible as such; no artificial motion is introduced to eliminate repeated-looking frames. The delivery report is authoritative for completion of these checks.

All four final videos passed the delivery checker: 1,001/1,201 frames at exact 60 fps, complete decoding, correct source-state and active-ID mappings, and unchanged input hashes. The initial 3D scene and colorbar match v7 pixel for pixel over the compared 824,850-pixel region in both cases; deliberately changed titles and annotations are excluded from that style check. Case labels are absent from the video titles.

The local viewing entry beside the older videos is `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/visualization/B2_B3_lubrication_v8_20261006/`, a link to this delivery. Original v7 exports and all scientific results remain unchanged. This is a visualization-only task; it launches no scientific simulation or controller. An initial render with case labels was superseded after the user's title correction; its partial outputs are retained separately and are not the delivered videos.

### v9 correction: physical centers in complete transport videos

The v8 primary videos passed their encoding checks but had a misleading display transform. A subsequent exact-mesh audit found no exterior physical centers in any of the 301 stored states per case. However, 60-fold cluster-relative coordinate magnification moved displayed centers as far as 2.264/1.991 mm outside the vessel. Regrouping also changed the magnification origin discontinuously: for example, B3 particle 220 moved 6.55 µm physically during t=1.92–1.94 s while its v8 display position jumped 2.543 mm. V8 primary footage is therefore superseded for transport interpretation. Evidence: `visualization_audit/20261006_v8_geometry/AUDIT.json`.

The replacement is `scripts/stage4l/render_transport.py`, exported to `visualization/B2_B3_physical_transport_v9_20261006/`. It retains the camera, field, PZT, black background, shaded blue symbols, clock and count styling, but removes cluster-relative coordinate magnification. Each original stored center is preserved bitwise. Group membership affects annotation only. Enlarged sphere symbols are composited through a silhouette obtained by projecting the actual triangulated vessel mesh; the inflation is a visibility aid, not a physical bubble size, and overlapping glyphs do not establish physical overlap or fusion.

Between 20 ms stored samples, each ID follows a linear display chord, respecting the recorded absorbing branch-crossing event. Raster cross-dissolves are removed. Only unsaved intermediate drawing points are projected against the exact mesh where needed: maximum corrections are 4.926/5.655 µm, about 0.252/0.289 pixels in this camera. This is explicitly display interpolation, not new solver output, recovered contact dynamics or scientific trajectory smoothing. The original 301 states, HDF5 files and branch events are unchanged. Group-count labels still refer to the latest stored state.

The v9 checks verify saved-center identity, every active display center against the same wall mesh, absence of reactivation, the bounded display-only corrections, and glyph-layer clipping on every rendered frame. The validator independently reloads saved data and recomputes mesh distances, as well as checking decoding and uniform presentation times. Maximum ordinary frame displacements are 68.66/77.98 µm; the earlier millimetre jumps from regrouping are absent. Near-wall negative gaps remain no larger than the original approximately 1 pm geometric tolerance. These drawing checks do not establish local physical accuracy.

The viewing entry is `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/visualization/B2_B3_physical_transport_v9_20261006/`. Older exports remain preserved. This correction creates two complete primary-only videos and launches no new scientific simulation.

Both final v9 videos passed the independent delivery validator: 1,201 frames each at exact 60 fps, full decoding, uniform physical-time increments and unchanged source hashes. The mesh check covered 237,547/229,079 active display samples with zero exterior centers; all 301 stored position arrays per case remain bitwise exact. The renderer verified glyph clipping on all 1,201 frames per case, and there were no inactive-to-active transitions. `DELIVERY_VERIFICATION.json`, `DISPLAY_GEOMETRY_CHECKS.json` and `DELIVERY_SHA256.json` retain the checks and artifact identities.

### v10: exact recorded morphology and the complete vessel overview

Following the user's explicit requirement for both unsimplified aggregation and a global vessel view, `scripts/stage4l/render_morphology.py` uses a fixed full-vessel camera. It does not resume the stopped focus-camera prototype. The original exterior mesh extends upstream to X = −14 mm, beyond the former −3.5 mm view limit. The new scene includes that inlet and both branches, renders all 13,978 original WALL triangles without decimation, and leaves the three ports open. WALL vertex coordinates are checked against the authoritative exterior mesh by GlobalNodeID. Black background, shaded blue bubbles, PZT discs, purple field, English axes and population annotations remain; case identifiers remain absent from video images.

The v9 visibility spheres still had 60-fold inflated radii, which obscured physical packing, and its intermediate position chords could overlap. V10 uses every original saved center and radius at a common physical scale, including the solver's original approximately 1 pm wall and 10 pm pair tolerances. No sphere inflation, cluster-relative transform, point filtering, contact projection, synthetic packing, particle interpolation or raster cross-dissolve is used. Every active ID is submitted; true depth occlusion is allowed. Analytic orthographic sphere intersections use sparse depth buffers and 32×32 coverage samples per pixel, without minimum-size markers. Three focused rendering tests cover submission-order-independent depth, subpixel size and preserved separation/radius ratios.

The complete 0–6 s primary-only videos are 3840×2160 at 60 fps and constant 0.30× playback. All 301 stored states are retained, each for four frames except the last endpoint; the clock shows the actual stored time. This is explicitly a 50 Hz scientific record played at 0.30×, not invented intermediate contact dynamics or 60 independent physical states per playback second. Full geometry at this global scale is subpixel: the video cannot resolve all micrometre gaps. Full-scene SVG exports at 0/3/6 s therefore retain one exact projected circle per active ID, with original 3D positions and radii embedded for checking. Their illustrative gradients do not alter geometry. No close-up panel or animated camera is added.

Output: `visualization/B2_B3_full_vessel_morphology_v10_20261006/`. The independent `check_morphology.py` validates source hashes, bitwise saved geometry, all active IDs, raw connected components, a single fixed isotropic projection, vector-circle coordinates/radii, exact frame mapping, full decoding and presentation timestamps. The delivery report records the final execution result; no scientific simulation is launched.

Both final videos passed: 1,201 frames each at exact 60 fps, full decoding, bitwise original centers/radii, all active IDs and exact vector geometry. Projected-center errors were below 5.5e-12 pixels. The true single-bubble diameter is approximately 0.2346 pixels even at 4K; raster footage cannot resolve its internal morphology, and the full-view SVGs preserve the geometry independently of raster resolution. Label leaders start six pixels beyond the projected group envelope to avoid masking that small silhouette. A separate boundary check found exactly three open wall loops. Local focused tests: 3 passed; Ruff and Python compilation passed. The final viewing entry beside the previous exports is `results/stage4/visualization/B2_B3_full_vessel_morphology_v10_20261006/` in the original project, linked to this isolated delivery.

### v11: visible similarity glyphs and reference-style PZT bodies

The user explicitly selected whole-group display magnification: every sphere radius and every within-group displacement share one scale, with that scale visible in the movie. `render_magnified.py` uses a constant 24× for both completed primary-only cases. Each original active ID is retained as an analytic shaded sphere; no point is merged, filtered, randomly repacked, rotated, or replaced by a group surrogate. A connected component is defined by the original 10 µm center-distance threshold. For each component the display map is `y_i = x_anchor + 24 (x_i − x_anchor) + shift_group` and `R_display = 24 R_native`.

Rigid group translations keep every inflated sphere inside the exact lumen and enforce a 5 µm display gap between distinct groups, avoiding new false intersections caused by magnification. The whole-group shift is carried between saved states, with 0.99 relaxation toward the native location and bounded local constrained placement. These are illustrative glyph locations: absolute centers, between-group distances, wall distances, apparent displacements and velocities must not be read from this movie as physical measurements. Native arrays, clustering, outcomes and force models remain unchanged. Regrouping and shape changes between the saved states can still change the display layout; there is no claim that those changes are new physical motion.

The enlarged fixed camera retains the complete inlet to X = −14 mm, both branches and the original 13,978-triangle open WALL mesh. The black background, purple acoustic slices, blue shaded spheres, English axes, upper-left clock and largest-group count remain. The static 3D context now uses VTK depth peeling. PZT front planes and aperture diameters follow the actual configurations; grey-blue front/back/side materials, thin rims and outward reflection lights reproduce the frozen `PZT_control_reference_v6/paper_pzt_reference.py` recipe (SHA256 `207c650fcc2a7d2714d42d6f7cb6aab9d4d26bcb0e29c66d88f91803fcd101ab`). The 0.35 mm backing is a display-only illustration, not a measured transducer thickness. Collision checks and reversed-actor rendering guard PZT intersection and context transparency. Blue glyphs are composited above the translucent context for visibility, with analytic sphere-to-sphere depth occlusion; they are not an optical transmission simulation.

Both primary movies retain the full native 0–6 s record, 301 saved states at 20 ms, four video frames per state and one final endpoint frame: 1,201 frames at 3840×2160/60 fps, constant 0.30× playback. No unsaved sphere coordinates or blended cluster morphologies are invented. “60 fps” describes encoding; the available independent scientific states appear at 15 per playback second. Case identifiers remain absent from the video images. Incomplete secondary runs remain excluded, and no simulation is restarted.

Output: `visualization/B2_B3_visible_morphology_v11_20261006/`. `DISPLAY_CHECKS.json` preserves the layout and scene audit; `*_RECORDED_SPHERES.npz` retains exact native arrays, `*_DISPLAY_SPHERES.npz` the separate glyph arrays, and the compressed group-transform and sphere-projection records expose every submitted identity. `check_magnified.py` independently checks source hashes, original arrays, common radius/relative-vector scaling, complete ID coverage, actual wall gaps, nonintersecting distinct groups, projection, complete decoding and constant timestamps. The delivery verification file is authoritative for the final execution result.

Both final v11 videos passed full validation: all 301 states and 59,544/57,429 active samples retain exact native arrays; every component keeps common 24× relative-vector/radius scaling (maximum vector error below 9e-19 m). All displayed spheres fit the exact lumen and distinct components retain at least 5 µm display clearance. Their projected single-sphere diameter is 9.707 pixels at 4K. Reversed static actor order produces zero RGB difference, and PZT body intersections are zero. Both files fully decode, contain 1,201 frames and have exact 1/60 s timestamp increments. Three encoded keyframes per movie retain more than 99.9% of checked blue sphere interior pixels. Five focused tests passed. The maximum native-step display displacement remains 1.216/1.427 mm, including similarity mapping changes at component merges; this is not a physical displacement claim. The separate native records and display arrays remain available for that distinction.
