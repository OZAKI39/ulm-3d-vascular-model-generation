# F4 — Original Y vessel integration

F4 adds geometric vessel confinement to the existing F3 GPU simulator. The original Stage 4 Y mesh and all four boundary labels are preserved. Two actual-PZT, 32-bubble experiments form connected clusters while respecting the WALL budget; the same-state PZT continuation has a small computed response. This is **EXPLORATORY / GEOMETRIC_WALL_ONLY**, with zero background flow. It does not validate vascular microbubble hydrodynamics or reliable cluster-shape control.

## Identity and reproducibility

Base: `dev/fast-pzt-aggregation-f3-20261008`, commit `3eeb94eb885f5df4f0e214ba9b7904dbdecec744`.
Delivery branch: `dev/fast-pzt-f4-y-geometry-20261008`.
Project: `computation_examples/mevo_acoustic_microbubbles/`.
Server root: `/workspace/mevo_fast_pzt_f4`; each simulation used a committed, immutable `.execution/<SHA>` tree. See `results/fast_pzt/f4/EXECUTION_INDEX.json` for actual execution identities, and `INPUT_MANIFEST.json` / `ARTIFACT_MANIFEST.json` for SHA-256 and full-data locations. An archive commit is not substituted for a simulation commit.

Original mesh: `/workspace/mevo_acoustic_microbubbles/results/stage4/geometry/20261002T135639758501Z/labelled_lumen_m.vtp`.
SHA-256: `e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b`.
The file was read and hashed on the server and copied without regeneration. It contains 25,319 points and 50,634 triangles: WALL 50,068; INLET 254; OUTLET_A 156; OUTLET_B 156. Units are metres. The topology is closed, consistently outward, nondegenerate and manifold; this audit is not a new exhaustive self-intersection proof.

This is the original **parametric smooth Y**, not patient-specific anatomy: parent length/radius 14/1.3 mm; daughter length/radius 12/1.0 mm; ±30° branches. No geometry scale change, remeshing, replacement cylinders or shrinkage was used.

## Field compatibility and initialization

The original field remains 205 × 37 × 37 at 100 µm spacing, spanning X ±10.2 mm and Y/Z ±1.8 mm. The original vessel spans X −14.0000 to +10.8923 mm and Y ±6.8660 mm. Thus it is **not fully covered**.

Valid bubble centres satisfy all three conditions: inside the closed original lumen; distance to actual WALL at least R; and inside the field box inset by R + 20 nm = 1.8 µm. A bounding box alone is not treated as lumen. The field border is a data-domain exit, never a reflecting wall. No acoustic array extension, duplication or extrapolation was performed.

A 2,001-point centreline survey finds support from the junction to 10.192 mm into the parent (7 µm survey spacing), and 3.594 mm along each daughter (6 µm spacing). These are sampled centreline extents, not exact full-cross-section coverage. The complete definition, finite-radius initialization positions and original wall normals are in `GEOMETRY_FIELD_OVERLAP.json`.

Two 32-bubble clouds were initialized by rejection sampling on the actual finite-triangle geometry, with seed 20261008, tangent disk radius 18 µm, positive wall clearance and pair surface gap at least 0.2 µm. Minimum actual initial wall gaps are 0.3804 µm (parent) and 0.2253 µm (junction). Invalid candidates were rejected; initialization did not project bubbles out of solid tissue.

## GPU implementation

`vessel_geometry.py` loads and audits the original VTP, preserves its labels and prepares exact stored coordinates. `GPUVessel` creates two static Warp BVHs once: the closed labelled mesh for inside/outside and open WALL triangles for physical distance. Original coordinates were stored as float32; they are promoted to binary64 without inventing additional mesh accuracy. GPU broad phase uses local float32 coordinates and padded BVH queries; the selected finite triangles are refined in binary64. The broad-phase padding is not a relaxed contact tolerance.

`F4VesselSimulator` extends F3, keeping resident positions, static mesh/work arrays and CUDA graphs. `vessel.enabled=false` delegates to the original F3 path. Enabled kernels add active-ID masks; the existing acoustic, feedback, stable pair reduction, two local-mobility iterations and twelve central-PBD pair passes remain unchanged. WALL passes follow prediction and pair passes 3, 6, 9 and 12. No CPU VTK query occurs in the physical step. GPU geometry observations and CPU transfers occur at save points and are included in measured loop cost.

For a sufficiently distant bubble, signed-distance Lipschitz clearance bounds skip expensive work. Near WALL, finite-triangle endpoint and swept-segment distance tests cover faces, finite edges and vertices. Oversize corrections or unsafe paths reject the candidate before acceptance. Post-correction pair and wall checks precede state acceptance; saved positions are never patched after acceptance.

Physical WALL requires h ≥ −0.01R = −17.8 nm; pair overlap has the same budget. A single pair or wall pass cannot correct more than 0.1R = 178 nm. Wall projection uses only the local normal; there is no tangential arrest, adhesion or rigid-cluster lock.

INLET/OUTLET caps determine closed-mesh sign but are open physical terminals. Finite-cap crossing and rim-wall clearance are distinguished. Accepted exit events retain original ID, label and interpolated crossing time, then deactivate the particle in all acoustic, pair, contact and cluster-link kernels. Field exit is separately labelled `FIELD_DOMAIN_EXIT`. Exited coordinates/IDs remain recorded; fixed-population metrics are withheld after exit. No science trajectory in this delivery exited either domain.

Implementation reference: [NVIDIA Warp 1.18.0 mesh query implementation](https://raw.githubusercontent.com/NVIDIA/warp/v1.18.0/warp/native/mesh.h). The existing project exact-rational finite-triangle implementation was used only in a few independent validation cases, not per simulation step.

## Validation and numerical quality

`checks_v2` passed: 124 actual-mesh points compared with independent VTK signed WALL distances, including near-face/edge/vertex cases. Maximum distance discrepancy is 8.61×10⁻¹⁹ m. This is computational agreement with the stored mesh, not physical surface accuracy. Signs agree away from the 10 nm surface exclusion band. Original float32 mesh discretization remains a scientific limitation.

Prescribed velocities are used only in explicitly labelled **GEOMETRY_UNIT_TEST** cases: approach/contact, tangential sliding (>99% prescribed tangent displacement retained), coupled pair/wall contact, and all three terminal exits. A daughter-to-daughter path through solid tissue is rejected. Cap centres exit at their recorded 20 µs crossing times without reflection. A field-boundary crossing deactivates the ID without a bounce. An inactive test ID cannot influence an active bubble. Disabled F4 and original F3 agree bitwise.

`edge_checks_v1` additionally compared five crossing/parallel/finite-edge/vertex/stationary paths against existing exact-rational squared distances on an actual triangle; maximum distance discrepancy was 3.82×10⁻²¹ m. A simultaneous cap/rim crossing through physical WALL is rejected, without accepting an exit or hiding the particle.

All eight actual-PZT trajectories completed at 10 µs, with one physical substep and no adaptive extra substeps. Maximum accepted pair penetration was 14.38 nm, below 17.8 nm. Maximum WALL pass correction was 24.03 nm, below 178 nm. WALL-enabled minimum accepted gaps remained at floating-point roundoff (~10⁻¹⁹ m); maximum physical-step displacement across science trajectories was 105.64 nm. No NaN/Inf, unrecorded exit or field extrapolation occurred. No difficulty justified a 5 µs comparison, and no time-step scan was performed. Passing these checks is not a proof of local physical accuracy or absence of dense-model rearrangement sensitivity.

## Actual near-wall aggregation

All cases use R = 1.78 µm, 2.5 MHz, eight response harmonics, Secondary ON, F2 near-field feedback, F1 central PBD and F3 stable reduction. No imposed velocity is used in these experiments.

Before any acoustic trajectory, field-only preflight found that the inherited zero-phase controls were weak at the new ROIs (124.2 Pa RMS parent; 50.2 Pa junction). The experiment plan was amended before running trajectories: all four amplitudes at the inherited conservative campaign cap of 291.2377 Pa, phases conjugate to each actual basis contribution at the fixed initial cloud centre, and a fixed 600 ms formation duration. The cap is not a measured hardware rating. The resulting initial-particle RMS fields are 295.8 Pa and 687.9 Pa. See `configs/fast_pzt/f4/controls.json` and `FAST_PZT_F4_PLAN.md` for complete phase vectors and preflight provenance. No trajectory-driven control search or force rescaling was used.

| Region / WALL | First saved 32-ID cluster | Initial → final Rg | Minimum saved WALL gap | Loop / process-body time |
|---|---:|---:|---:|---:|
| Parent OFF | 101.2 ms | 12.855 → 4.864 µm | −1.760 µm | 27.66 / 29.34 s |
| Parent ON | 101.2 ms | 12.855 → 4.900 µm | roundoff near zero | 44.35 / 46.07 s |
| Junction OFF | 19.0 ms | 12.862 → 4.895 µm | −2.982 µm | 27.92 / 30.02 s |
| Junction ON | 19.0 ms | 12.862 → 4.874 µm | roundoff near zero | 33.30 / 35.10 s |

Each pair starts from bitwise-identical coordinates. Both start with 32 separate connected components and end with one component of all 32 IDs. A cluster is a spatial graph with centre distance ≤2.1R; it does not mean fusion or adhesion. First saved wall contact occurs at 63.6 ms parent and 17.6 ms junction; at most four WALL-enabled bubbles meet the saved h ≤0.01R contact criterion simultaneously. Contact begins under actual acoustic forcing. WALL OFF deliberately allows physical wall penetration; its general F3 numerical PASS must not be read as a WALL-validation PASS.

At 600 ms, WALL ON/OFF centroid differences are 1.538 µm parent and 2.327 µm junction. WALL changes constrained geometry; it does not create the underlying attraction or cluster formation. Parent and junction cumulative WALL corrections sum to 49.38 and 75.92 µm over all particles/passes. Dense pair correction is approximately 0.92 and 1.01 times the summed mobility-predicted path. These large cumulative ratios remain a model limitation, despite small per-pass corrections.

## Same-state PZT continuation

The predeclared parent selection rule selected the actual 600 ms parent WALL ON cluster (32 connected IDs). The exact same state was continued for 120 ms in a 2 × 2 WALL ON/OFF × unchanged/switched design. The only switch was Y+ phase +30°; amplitudes, frequency and enable flags stayed fixed. It is not an exactly pressure-matched control (initial-cloud RMS changed by about 2%).

The cluster had detached slightly before continuation: minimum initial wall gap 0.1585 µm. Neither control produced a new wall contact, and each control's WALL ON/OFF trajectories are bitwise identical. Thus this continuation demonstrates correct inactive-wall behaviour, not a new wall-interaction control effect.

At 120 ms, switched versus unchanged centroid separation is **13.91 nm**; relative Rg difference 0.006285%; eigenvalue difference 0.08109%; tensor difference 0.3962%; aligned internal RMSD 16.60 nm. The component remains all 32 IDs, with no component split/merge. Kabsch rotation differences and orientation series are recorded but are not called reliable controlled rotation. Principal-axis angles are reported only when the fixed-ID eigenvalue-gap criterion (>0.1) makes them identifiable.

This is a weak computed PZT response. Robust shape control is **not established**. Prior F3 sensitivity to picometre/nanometre initial perturbations, angular-balance residuals from truncated mobility and the stored-gradient/interpolated-pressure mismatch remain. The older F3 sensitivity envelope is not relabelled as a calibrated uncertainty interval for these different F4 controls and positions.

## Performance

32 bubbles, identical 600 ms physical durations: WALL ON/OFF loop cost ratios are **1.60 parent** and **1.19 junction**. Each trajectory executed 60,000 physical steps, zero additional substeps and five scheduled wall passes per step. Parent/junction physical kernels made 294,117 / 48,968 refined queries and 131,902 / 16,528 swept triangle tests; most passes used the Lipschitz skip. Every counter is in the original run summaries. Initialization/save-point diagnostic queries are separate from these physical-kernel counters.

The bounded N200, 5 ms regression reused the original F3 initial state, uniform 1000 Pa, stable reduction, 500 steps and identical sample cadence. WALL OFF **0.7284 s**, WALL ON **1.1902 s**, ratio **1.634**. Positions agree bitwise because this cloud is far from WALL. It is a performance-only test, not a further aggregation validation.

CUDA events on the resident/warmed GPU graph (100 repetitions) measured: exact 32-point distance query 87.58 µs; cached far-clearance 32-particle wall pass 2.97 µs; single near-contact query plus 1 nm projection including candidate device-copy 70.17 µs. Query and projection share a kernel; the latter is honestly reported as combined cost, not an invented isolated subtraction. These microtimings are state-dependent and not additive. The simulation loop includes saved-state transfers/observations; setup/compilation/BVH/export costs are separately recorded. No repeated BVH construction or CPU VTK is in the physical loop.

## Figures and video

All figure text is English. PNG and PDF outputs are in `results/fast_pzt/f4/figures/`:

- `Fig_F4_1_geometry_and_field`: original open-ended Y surface, labelled rims, current stored-field intersection and both actual initial clouds.
- `Fig_F4_2_near_wall_aggregation`: actual dispersed/intermediate/final states in both ROIs plus cluster size, Rg, wall gap, contacts and correction traces.
- `Fig_F4_3_PZT_control`: same-state unchanged/switched positions and centroid, Rg, eigenvalue, clearance and contact metrics.

`visualization/F4_Y_vessel_aggregation_and_PZT_control_60fps.mp4` shows parent formation and switched continuation over 0–720 ms. It uses 721 actual saved states, 60 fps and 1 ms physical time per video frame (12.017 s duration). Global and local views share physical time. The camera is fixed, every sphere retains R = 1.78 µm, and there is no coordinate interpolation, trajectory smoothing, artificial cluster enlargement or particle hiding. The 0.6 mm global outline marks the ROI location; it is not a sphere or a rescaled cluster. Local coordinates are µm relative to a recorded global origin. The rendered wall is the original finite-triangle surface, clipped only for the local display; its locally flat appearance follows the original mesh facets.

`FRAME_MAP.npz` records every source frame and physical time. Full decoding and timestamp checks are in `VIDEO_VALIDATION.json`. Movie smoothness is not used as a physical-accuracy test.

## Scope and operational record

No CFD/Poiseuille flow, wall lubrication, acoustic scattering/image bubble, adhesion/friction, deformation, coalescence, clot model or QP has been added. Geometry is only confinement. The existing original mesh and 100 µm field grid limit local interpretation. Terminal crossing was tested with prescribed unit velocities because the actual stored field does not reach the terminal caps.

The first basic-check attempt (`checks_v1`) stopped before acoustic experiments because its pair unit-test initialization was too close to a faceted wall edge. The corrected test starts with 0.2 µm clearance. All original thresholds were preserved. Partial outputs and logs are retained.

A cloned launcher initially directed the first F4 check's logs to the F3 `checks_v1` log paths. Only logs were appended. After that F4 process exited, the appended suffixes were copied into F4; the F3 prefixes were restored byte-for-byte using original archive sizes and SHA-256. The launcher was fixed before subsequent jobs. `LOG_RECOVERY.json` records the recovery, and `PRESERVATION.json` verifies every inherited F1/F2/F3 manifest file on both hosts. No F3 scientific output was changed or rerun. An offline VTK annotation API mismatch was fixed during preview generation; it affected no trajectory.

See `scripts/fast_pzt_f4/README.md` for reproduction. Lossless `SAVED_*.npz` copies retain every saved coordinate, ID, metric and diagnostic; full HDF5 trajectories remain indexed at their local/server paths. Large reused basis/response inputs are referenced by exact hash rather than duplicated. Failed attempts are retained separately from accepted evidence.
