# F5 — Original Y vessel, prescribed flow, and acoustic release

The registered 32-bubble experiment completed 600 ms at 10 µs. Continuous flow did not prevent aggregation: all 32 bubbles formed one connected group before the 300 ms switch. After all ultrasound was switched off, the same particles continued downstream, elongated, rearranged and separated into multiple spatial components. The matched no-flow continuation stayed essentially stationary. These are numerical results within `EXPLORATORY / PRESCRIBED_STEADY_FLOW / GEOMETRIC_WALL_ONLY`, not validated blood-flow or bubble-breakup predictions.

## Baseline and execution identity

- Branch: `dev/fast-pzt-f5-flow-release-20261008`.
- Fetched F4 baseline: `e247fa0c61f668a9b0b2e26c909dd687eedde389`; no later F4 commits were found at handoff.
- Isolated local worktree: `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f5_20261008`.
- Server root: `/workspace/mevo_fast_pzt_f5`. Immutable source snapshots are in `.execution/<commit>`; supervisor jobs are bounded, uniquely named and do not auto-restart.
- Scientific trajectories: commit `ee6b35dce83bd3583457a7e85c163b6886e59ce5`, directory `results/fast_pzt/f5/science_v1`.
- Completed short timing: commit `9aba27e3e4bb4a586e0b5a1ecfe0c3b91f769eb9`, `performance_v1`.
- Final GPU checks: commit `2661357dc61031d4af884be5f9c9c8cfd113416c`, `checks_v3`. Subsequent changes concern figures, delivery checks and documentation. The final diagnostics preserve archived ON and OFF 1 ms prefixes bitwise.
- F1–F4 implementation and prior experiment files remain unchanged. New files extend their classes; no legacy batches, CFD solves or QP campaigns were resumed.

## Exact inputs and units

CFD source on the server:

`/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/steady/20261002T141810826990Z/run/8-procs/result_121.vtu`

SHA-256: `131f7cb9b3e4ac9887bb0aa344d0eb2c0a2c8ac1d3b19fafb093ecd9a1045620`.

The original snapshot contains 26,712 nodes and 146,230 linear tetrahedra. Coordinates are metres; `Velocity` is metres/second. Maximum nodal speed is 0.010074017175615313 m/s. Archived steady-flow checks report `PASS_WITH_SCOPE`, inlet flux 2.6354060534101732e-8 m³/s, outlet fractions 0.4999637507 / 0.5000362493, and zero nodal wall velocity. `CFD_STATE.json` and the earlier CPU tetra-P1/VTK validation are included with the exact VTU.

Original labelled Y surface SHA-256:

`e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b`.

Its 50,634 triangles retain the original WALL and open INLET/OUTLET_A/OUTLET_B labels. No vessel remodeling or cap closure is introduced. Exact acoustic basis, response table, prepared wall geometry and initial clouds are included in F5 inputs; `INPUT_MANIFEST.json` records their hashes.

Despite the blood-flow application, this inherited CFD uses **water parameters** (density 1000 kg/m³, viscosity 0.001 Pa·s), consistent with the frozen fast model. It is not a blood-rheology calibration or a mesh-independent near-wall CFD solution.

## Implementation

`flow_gpu.py` loads the original mesh once, promotes its float32 node coordinates to float64, computes P1 inverse maps and cell gradients, and uploads static arrays and a padded float32 AABB BVH. Broad-phase padding is only a search safeguard. Actual membership and interpolation use float64 and the inherited dimensionless barycentric roundoff tolerance of 1e-10. A last-cell cache avoids repeated BVH searches. Invalid/nonfinite queries raise explicit error bit 16384; no nearest-neighbor extrapolation or valid-zero substitution occurs. The exact derivative of the P1 velocity is constant inside each tetrahedron and discontinuous across faces.

`f5_engine.py` adds flow **before** mobility, prediction and contact:

`U_free = u_flow + F_ac / (6*pi*mu*R)`.

The inherited two local Jacobi mobility iterations and component-mean balance approximate `(D + L_pp) U = D u_flow + F_ac`. Common translation is retained while relative motion remains free. This is not a converged global hydrodynamic solve. Central PBD and geometric wall projections remain numerical nonpenetration constraints; their corrections are not physical shear forces.

Each physical step queries flow at the current positions and validates proposed accepted endpoints on GPU. Legal cap crossings retain their exact F4 event time/ID and bypass an unnecessary outside-CFD endpoint lookup. CPU downloads occur at controls and saved frames (every 0.5 ms), never for per-step VTK sampling.

Two independent semantics are implemented and tested:

- `SECONDARY_OFF_PZT_ON`: pair acoustic force is zero, while primary force can remain nonzero.
- `ULTRASOUND_OFF_FLOW_ON`: both acoustic forces are exactly zero; background flow remains active. The main experiment uses this mode.

Changing force switches invalidates captured CUDA graphs without moving particle positions. Ultrasound OFF skips field/response queries and disables the acoustic-grid exit check; real WALL and terminal boundaries remain active. Acoustic ON retains the F4 `FIELD_DOMAIN_EXIT` policy. Invalid flow while flow is active fails the calculation. Disabled-flow diagnostics retain explicit validity flags without restricting the physical F4 domain.

## Registered experiment

The original F4 **Parent** cloud and PZT control were selected by static preflight before scientific integration. No positions or controls were tuned after observing the result. All 32 initial centres were valid in the geometry, CFD and acoustic domains. Initial speeds span approximately 15–126 µm/s, with maximum pairwise background-velocity difference 111.26 µm/s. Minimum initial wall gap is 0.3804 µm.

- 32 bubbles, immutable IDs 0–31, R = 1.78 µm, frequency 2.5 MHz.
- Frozen F2 near field, F1 central PBD, F3 stable reduction, F4 WALL.
- Four original PZT elements, each 291.237687721856 Pa drive amplitude; original Parent phases and geometry retained.
- 0–300 ms: acoustics ON, flow ON.
- 300–600 ms: all PZT excitation OFF, secondary force OFF, flow ON.
- Fixed physical step 10 µs, no adaptive substeps, no velocity clipping or added inertia/repulsion/adhesion.
- Branch A is the exact uninterrupted 300–600 ms main suffix; it is not rerun or double-counted.
- Branch B restarts from the **identical saved switch positions, IDs, radii and active mask**, with ultrasound OFF and flow OFF for 300 ms.

The model is overdamped: velocity, mobility work arrays and lookup caches are recomputed algebraically, not hidden inertial state. `SWITCH_STATE.npz` stores the physical state and saved diagnostics. `SWITCH_STATE.json` explains cache reconstruction. A continued-acoustics Case C was not necessary for this registered flow/no-flow comparison. No 5 µs run was triggered: the prescribed numerical failure/overlap gates did not fail. This is not a time-step convergence claim.

## Verification

`checks_v1` stopped during initialization because the inherited environment lacked PyVista. A separate F5 environment was created; F4 was untouched. PyVista 0.46.5 and VTK 9.5.2 were installed and NumPy pinned back to 2.2.6 to match the inherited numerical stack. The actual Python, packages, driver and GPU are recorded in `ENVIRONMENT_GPU.json`.

The GPU checks passed:

| Check | Result |
|---|---:|
| 44 actual positions: Parent, bifurcation, daughter, near-wall | Max GPU/CPU velocity error 1.735e-18 m/s |
| P1 analytical gradient versus CPU | Zero observed difference |
| Gradient versus CPU finite differences | Max 3.418e-9 s⁻¹ |
| Uniform translation, single and four near-contact bubbles | Max position error 2.386e-17 m |
| Uniform-flow four-bubble Rg change | 4.236e-21 m |
| Simple shear relative motion | 0.320 µm; position error 4.251e-17 m |
| Short real-CFD single-bubble Euler trajectory versus CPU at same dt | Bitwise endpoint agreement |
| Euler integration versus exact local affine P1 solution | 1.012e-11 m after 0.5 ms |
| Ultrasound OFF outside acoustic box | Continues 8.421 µm in 1 ms |
| Acoustic ON grid exit and physical outlet exit | Correct distinct labels and crossing times |
| Disabled-flow F4 compatibility | 100 steps bitwise identical |
| Outside-CFD and nonfinite flow queries | Rejected explicitly |
| Flow OFF with unavailable CFD diagnostics | Motion unrestricted; diagnostic validity false |
| Final source versus completed main ON/OFF prefixes | Bitwise identical at common saved times |

Interpolation error and integration error are separate measurements. Neither validates the underlying CFD or near-contact physical closure. Delivery validation also compares 672 positions from the actual completed trajectory with the independent CPU sampler.

## Scientific outcome

All 32 original particles remain active through 600 ms. There are no physical terminal exits, acoustic-domain exits or flow-domain failures in the main experiment.

| Observable | Initial | At ultrasound switch (300 ms) | Flow-driven endpoint (600 ms) |
|---|---:|---:|---:|
| Spatial components | 32 | 1 | 11 |
| Largest component | 1 | 32 | 7 |
| Fixed-ID Rg (µm) | 12.8555 | 4.9191 | 12.6376 |
| Centroid distance from initial (µm) | 0 | 21.4474 | 46.0935 |
| Shape aspect ratio sqrt(λ1/λ3) | 4.5247 | 1.1128 | 10.8660 |

With flow, the first saved fully connected cloud occurs at 92.5 ms. The unchanged F4 no-flow reference reaches full connectivity at 101.2 ms (different saved cadences: 0.5 and 0.2 ms). Both form 32-bubble groups by 300 ms. With-flow Rg is 4.9191 µm versus 4.9005 µm without flow, about 0.38% larger. These are one matched scene's outcomes, not a general statistical claim that flow accelerates aggregation.

The **release-only** displacement is 24.6621 µm, including 24.6616 µm in X and 0.1606 µm transverse. Mean centroid path speed is 82.207 µm/s. Rg increases by a factor of 2.569. Shape eigenvalues change from `[8.9767, 7.9717, 7.2494]` to `[149.1494, 9.2965, 1.2632]` µm². Fixed-ID RMSD after rigid alignment reaches **9.7643 µm**: this is genuine internal motion in the numerical model, not just translation or rotation.

At release, the mean prescribed gradient includes `∂ux/∂z ≈ 14.0674 s⁻¹`. Mean strain-tensor Frobenius norm is about 9.9630 s⁻¹. Maximum pairwise background-velocity difference ranges from 138.0 to 148.1 µm/s during release. Thus nearby particles at different wall-normal positions have different prescribed velocities; once acoustic attraction is removed, this gradient drives relative motion subject to local mobility and contact. No artificial dispersal force was added.

Connectivity uses the inherited centre distance `<= 2.1R` (surface gap <= 0.1R), not bubble fusion. Accepted-step graph histories show **11 persistent split events and 2 persistent merge events** after release, with persistence threshold 2 ms. One further split at 598.29 ms has only 1.71 ms of observation remaining and is **right-censored**. The final 11 components include this last observed separation; it is not counted as confirmed persistent. These are graph split/merge events, not gas-bubble fragmentation/coalescence.

Mean nearest-neighbor separation increases only from 3.55994 to 3.62309 µm because several close subgroups remain. There are 83 changed graph-neighbor pairs and 27 changed nearest IDs relative to the release state. Do not infer complete particle isolation from the global stretch.

Branch B retains one 32-bubble component, centroid displacement 4.87e-19 m, Rg increase 0.000241 µm (0.00490%), and aligned RMSD 0.342 nm. This tiny residual is contact relaxation of inherited sub-nanometre overlaps, not inertial coasting or physical rearrangement. Six nearest IDs change because near-equal distances at contact are sensitive to tiny corrections; graph-neighbor changes and persistent splits remain zero. B stores the original field/gradients as saved diagnostic values, but `flow_enabled` is false and its applied background velocity is exactly zero.

## Contact and WALL quality

Main-run maximum accepted pair overlap is **2.166 nm** (0.122% of R); maximum pair contact-pass correction is **1.682 nm**. Maximum predicted WALL penetration and one-pass WALL correction are **1.844 nm**. Minimum queried wall gap is -0.00212 nm; minimum saved accepted gap is effectively zero within floating-point error. All inherited numerical gates pass, with zero error bits. This is numerical nonpenetration within the specified tolerance, not exact mathematical nonoverlap or local physical validation.

Summed contact corrections over all particles/passes are 1.43516 mm for the entire main run, including 6.84047 µm during release. Summed wall corrections are 40.0925 µm, including 1.31629 µm during release. These are cumulative numerical path lengths, not net shifts, forces, or shear stresses. The main contact/free-path ratio is 0.2597; corrections are not negligible as a model mechanism. In Branch B the free path is zero, so its inherited floor-denominator ratio is physically undefined; use the absolute 8.214 nm summed correction instead. Raw inherited diagnostics are retained without rewriting history, and `DELIVERY_VALIDATION.json` marks this normalized B ratio null.

## Performance and retained attempts

GPU: NVIDIA GeForce RTX 4090, driver 595.84, Warp 1.18.0, float64 particle state. Main: **60,000 macro steps = 60,000 physical substeps**, 42.9769 s measured loop wall time. Branch B: **30,000 = 30,000**, 12.9358 s. These include host checks, saved-frame downloads and in-run graph work; they are not pure GPU times. A's 30,000 release steps are already included in the main 60,000.

One bounded warmed CUDA-event timing (100 repetitions, 32 particles) reports:

| Operation | Mean GPU time (µs) |
|---|---:|
| Cached flow query | 3.758 |
| BVH query including cache reset | 28.908 |
| Field and radial response | 15.591 |
| Secondary pair force | 23.982 |
| One pair contact pass plus scheduled WALL pass | 25.211 |
| Complete ON step | 544.020 |
| Complete OFF step | 496.876 |

These timings are not additive. Static flow preparation took 0.136 s in the timing job. The main run made 3,840,000 particle flow queries, with 3,839,964 last-cell hits and zero invalid queries; both current and proposed positions are checked. The heavy inherited contact/wall sequence remains the largest part of step cost.

`science_v1` finished and saved both successful trajectories, then its post-run component timer raised an uninitialized contact-pass-counter error. Only instrumentation was fixed. `performance_v1` completed short timing and verified both trajectory hashes were unchanged; the original exception and all logs remain. Later final checks improved disabled-flow diagnostics and nonfinite-input rejection without changing valid physical trajectories. There was no new main run, parameter search, mesh change or CFD solve.

## Figures and video

`results/fast_pzt/f5/visualization/` contains:

- `Fig_F5_1_flow_original_Y.{png,pdf,svg}`: original open Y surface, actual CFD speed on the labelled Z=0 plane, flow-direction arrows, original stored acoustic support and exact initial cloud. Geometry axes keep equal physical scale.
- `Fig_F5_2_attraction_release.{png,pdf,svg}`: actual initial, switch and final positions at one common local scale, with cluster/Rg/displacement metrics.
- `Fig_F5_3_release_comparison.{png,pdf,svg}`: same-state A/B transport, internal shape and connectivity.
- `Video_F5_acoustic_release_flow_60fps.mp4`: 1920×1080, pure black background, reference-v6 PZT materials, open red Y vessel, fixed local camera and true-radius blue spheres. PZT labels and drive appearance change at 300 ms.

All 1,201 saved frames are rendered once at 60 fps: 0.5 ms of simulation per frame, 20.0167 s of playback for 600 ms physical time. The control switch is frame 600 (10 s playback). No position interpolation, recentering, display radius inflation or forced morphology changes are used. Local coordinates are micrometres relative to the initial centroid; the global locator denotes the observed region, not an enlarged bubble. Static paper figures use white backgrounds; the requested movie uses black. Scientific limitations remain in this report instead of irrelevant movie footers.

`FRAME_MAP.npz`, `VIDEO_IDENTITY.json`, full decode timestamps/log, and `VISUALIZATION_VALIDATION.json` preserve source mapping and verify the entire encoded video. A smooth video does not establish physical accuracy.

## Reproduction and deliverables

From the project directory in this branch, install `configs/fast_pzt/f5/requirements-gpu.txt` in an isolated CPython 3.12 environment. Exact small input files and the essential original acoustic basis are included. For rendering, use the shared local visualization environment recorded in `ENVIRONMENT_RENDER.json` (PyVista, Matplotlib, Pillow, imageio-ffmpeg).

```bash
# Local static preflight (uses original archived F4 scene).
PYTHONPATH=src:. python -m scripts.fast_pzt_f5.prepare

# On the GPU host, with the isolated F5 root and committed source snapshot.
# Choose NEW labels: successful or failed runs are never overwritten.
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f5.checks --label checks_repro
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f5.science --label science_repro --checks checks_repro

# Analyze the delivered actual run locally; the F4 reference remains unchanged.
PYTHONPATH=src:. python -m scripts.fast_pzt_f5.analyze --label science_v1 --f4-reference /path/to/F4/aggregation_v1/parent_wall_on
PYTHONPATH=src:. python -m scripts.fast_pzt_f5.render --mode figures
PYTHONPATH=src:. python -m scripts.fast_pzt_f5.render --mode video
PYTHONPATH=src:. python -m scripts.fast_pzt_f5.validate_delivery --label science_v1 --checks checks_v3
```

For bounded remote jobs, use `scripts.fast_pzt_f5.remote` with a new label; it deploys only committed source and records its exact execution SHA. Existing reproduction inputs may be used directly; do not rerun CFD. `ACOUSTIC_RELEASE_PROTOCOL.json` is the pre-registered scene/protocol, while `ON_OFF_FLOW_COMPARISON.json`, `CLUSTER_TRANSPORT_METRICS.json`, `PERFORMANCE.json`, `FLOW_GPU_VALIDATION.json` and `FINAL_GATE.json` are the requested machine-readable outcomes. Full HDF5/NPZ data retain immutable particle IDs, membership, all saved gradients/strain, wall quality, active flags and accepted-step connectivity events.

Remaining scope: one selected 32-bubble near-wall scene, a stored prescribed steady water-parameter CFD snapshot, one-way coupling, approximate local pair lubrication, central PBD contact and geometric wall exclusion. No wall hydrodynamic lubrication, deformable bubble interface, inertia, adhesion, acoustic streaming feedback, two-way coupling, convergence campaign, or experimental validation is claimed. The observed numerical release is supported by actual local velocity gradients and internal displacements; it was not manufactured for the animation.
