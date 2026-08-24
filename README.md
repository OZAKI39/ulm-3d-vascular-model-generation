# ULM 3-D vascular preprocessing

## Included v9 reproducibility snapshot

This repository snapshot includes the source code, configuration, tests, technical notes, and the
minimal inputs used for the formal ROI003274 v9 validation. It intentionally excludes the complete
mouse-brain dataset, historical output runs, scratch files, and intermediate STL collections.

The retained formal run is:

`outputs/model_generate/validation_v9_roi003274_20260823_r1`

Its complete log, A/B/C/D/E diagnostics, same-camera figures, reports, and exported geometries are
included. The strict v9 result is `KEEP_V7_HARD_MIN_V9_ACCEPTANCE_FAILED`: all topology and fidelity
checks passed, but the visible saw-tooth acceptance conditions did not, so the published final lumen
surface is the verified v7 fallback rather than an incorrectly promoted v9 candidate.

The matching test inputs are limited to:

- `outputs/sampling/20260822_171206_radius_plus_structure_k5`
- `outputs/rodent_vasculature/all_run_20260822_171157`

VTP, STL, and PNG artifacts are stored with Git LFS. Run `git lfs pull` after cloning to retrieve
their full contents.

This folder implements the first three preprocessing stages described in the project notes:

1. STL cleanup and selection of one connected main vascular network;
2. uniform voxelization, voxel-island removal, and coarse 3-D skeleton extraction.
3. a hierarchical vascular representation containing branch topology, branch summaries,
   complete centerline sequences, coarse radius sequences, junction geometry, and cycles.

The coarse skeleton is a navigation aid for the later branch-graph stage. It is not final CFD geometry.

## Mouse-brain TIFF/SWC workflow

`main_process_rod.py` processes the dataset under
`vessel_model/T - A high-resolution dataset of mouse brain vasculature`. With no arguments it
processes the first eligible 192 x 192 x 192 block from the 108-block analysis cohort, writes
acceptance evidence, and opens a native PyVista window:

```powershell
python main_process_rod.py
```

Use the mouse to rotate, pan, and zoom. The white/gray raw TIFF volume on black follows the visual
language of Figure 2(a). The left viewport retains the complete block and marks several colored,
clickable ROI cubes. They are real connected/acyclic SWC subtrees ranked by bifurcations, terminals,
depth, centerline length, and tortuosity; ROI 1 is the highest-ranked default. Left-click any cube to
replace the right viewport with that candidate. The right viewport follows Figure 1(b)'s local-cube
idea and suppresses all surrounding vessels, so it remains tree-like rather than dense. Cyan lines
show centerlines, colored points show tree roles, and orange arrows show the saved parent-to-current
relationship. `--tree-roi-branches`, `--tree-roi-depth`, `--tree-roi-padding`, and
`--tree-roi-samples` control this extraction. Both viewports show numerical X/Y/Z axes in
micrometres; the local cube keeps the source model's global coordinates so its physical location
and scale remain explicit. For CI or a headless acceptance run, create the same outputs without
opening the window:

```powershell
python main_process_rod.py --no-show
```

After the directed graph is built, the default command also runs the real connected-ROI sampling
pipeline described in `references/(new) 连通子图plan.md`. This stage uses spatially separated SWC
node anchors, exact line/box intersections, and retains only the connected component containing
each anchor. Original global node IDs and deterministic global edge IDs are retained. A boundary
intersection is exported as `CUT_PORT`; an endpoint that is already degree one in the full graph is
exported separately as `TRUE_TERMINAL`. No physiological inlet/outlet, flow, CFD, or synthetic
vasculature is inferred by this stage.

The default descriptor is `radius_plus_structure`: arc-length-weighted radius P10/P25/P50/P75/P90,
branch count, bifurcation count, total vessel length, and cycle rank. Features are robust-scaled,
clustered with deterministic KMeans, and mapped back to actual candidate ROIs. The same candidate
pool is also evaluated with `radius_only`, providing a direct ablation comparison. Representative
selection supports `coverage_balanced` and `distribution_preserving`, with configurable spatial
overlap control. Use `--no-sampling` to retain the earlier graph-only workflow, or inspect all
sampling options with `--help`.

Each run is stored under `outputs/sampling/<run_id>/` with its config, sampling log, candidate and
selected manifests, global edge mapping, exact cut-port table, feature table, robust scaler,
cluster assignments/centers, per-ROI NPZ geometry, validation metrics, comparison experiments, and
acceptance figures. In the interactive GUI, the sampling layer is display-only: press `A` for all
candidates, `R` or `S` for selected representatives, `C` to cycle clusters, and left-click a box to inspect
that saved connected ROI and its metrics. Both viewports retain physical micrometre coordinates.

### Build CFD-ready lumen surfaces from representative ROIs

`model_generate.py` reads those saved sampling manifests and ROI NPZ files through the existing
sampling loader. It reconstructs the lumen directly from SWC centerlines and source radii; TIFF/Mask
data is never used as geometry. With no `--roi-id`, all selected representative ROIs are processed:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\model_generate.py `
  --sampling-run ulm_3D_vascular\outputs\sampling\<sampling_run> `
  --workers 1 `
  --headless
```

For an initial single-ROI check, pass the exact ID from `manifests/selected_rois.csv`:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\model_generate.py `
  --sampling-run ulm_3D_vascular\outputs\sampling\<sampling_run> `
  --roi-id <roi_id> `
  --workers 1 `
  --headless
```

Each invocation creates a non-overwriting `outputs/model_generate/<run_id>/` tree. The master VTP
surface retains `patch_id`, `patch_type`, and `port_id` face arrays. Micrometre VTP is used for visual
QC, while `*_m.vtp` and `*_m.stl` use SI metre coordinates. `units.json`, source branch/global-edge
mappings, CUT_PORT metadata (`CUT_PORT_UNASSIGNED`), collision reports, watertight/manifold checks,
radius-fidelity sections, 16/24/32/48-side convergence results, and seven acceptance figures are saved
with every successful ROI. No physiological inlet/outlet or blood-flow direction is inferred.

To diagnose a visually abnormal port or junction without changing the source SWC, radii, or formal
reconstruction settings, run the strict v2 diagnostic path serially:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\model_generate.py `
  --diagnose-roi <roi_id> `
  --workers 1 `
  --headless
```

Use `--diagnose-all` for every selected ROI. Results are isolated under
`outputs/model_generate/<run_id>/diagnostics/` and include pre/post-Boolean meshes, port and junction
CSVs, independent boundary/non-manifold/internal/self-intersection checks, triangle quality, synthetic
straight/Y controls, an explicit/implicit comparison, diagnostic overlays, `summary.json`, and
`diagnostic_report.md`.

Every SWC edge is interpreted as `parent_id node -> current node`. Branch arrays retain that
order, and the 3-D, orthogonal, and topology PNGs display it with arrows. This is a structural
direction inferred from annotation topology, not a measured pressure, velocity, or blood-flow
boundary condition. Results are written to `outputs/rodent_vasculature/<stage>_run_*`, including
logs, normalized arrays, CSV/GraphML/NPZ/VTP graph exports, visual evidence (including
`figure2a_interactive_preview.png`), a scene manifest, and HTML acceptance reports. `preprocess`
and `hierarchical-graph --source-run ...` may also be run separately. Use `--help` to select an
exact sample or cohort and to tune arrow count, window size, opacity, and output formats.

## Run the supplied test model

From the repository root:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\main.py
```

The default voxel size is 2 micrometers, which is suitable for a quick full-model check. If the acceptance report warns about a connectivity difference, compare it with a finer version:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\main.py --voxel-size 1.0
```

A finer grid retains more small structures but does not guarantee fewer connected components: voxelization can split thin connections, remove tiny parts, or merge nearby surfaces. Use the saved comparison images and reports when selecting the final resolution.

Every run creates a new timestamped folder below `ulm_3D_vascular/outputs`. Open `acceptance_report.html` in that folder to review statistics and images. `pipeline.log` contains the complete execution trace.

By default, the component with the largest surface area is selected as the main network and every disconnected component is removed from the downstream STL. Removed objects are preserved as separate artifacts: small fragments and larger island networks. The final voxel mask is filtered again so that only its largest 26-connected component reaches skeletonization.

Use `--main-component-id ID` to override the automatic main-network choice after reviewing `mesh_components.csv`. The former behavior can be reproduced with `--component-policy conservative` for comparison, but it does not satisfy the single-network project requirement.

Use `--help` to list all parameters. No input file is modified in place.

## Run Step 3 from a completed Step 2 run

Step 3 deliberately reads an existing timestamped run so that mesh cleanup and voxelization do
not need to be repeated:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\main.py `
  --stage hierarchical-graph `
  --source-run ulm_3D_vascular\outputs\mv03_y.nii.gz_Segment_1\run_20260813_173944
```

This creates a new `hierarchical_graph_run_*` folder beside the source run. The source run is
not modified. Open the new `acceptance_report.html` first, then review:

- `graph_skeleton_overlay.png` and `graph_reconstruction_difference.png` for aggregation fidelity;
- `hierarchical_graph_3d.png` for terminals, junctions, branches, and cycles;
- `branch_as_node_graph.png` for the later GNN-oriented representation;
- `coarse_radius_profiles.png` for navigation-level radius estimates.

`hierarchical_vascular_graph.json` is the complete readable representation. GraphML files carry
portable topology, VTP files carry 3-D lines and points, CSV files support manual inspection, and
`branch_geometry.npz` stores compact lossless numeric sequences.

The raw centerline voxel sequence is retained for every branch. Smoothed coordinates and all
radius/curvature values are derived navigation measurements. They are not final CFD geometry or
final Geometry Generator training truth. Flow direction is not guessed, so parent/daughter,
branch order, depth, and downstream subtree fields remain unavailable until inlet/outlet evidence
is supplied.

## Run the NNE2 Step 1-3 workflow

NNE2 starts from microscope TIFF stacks rather than STL surfaces. Its Step 1 cleans the segmented
image components, Step 2 writes a coarse centerline and radius field, and Step 3 builds the common
branch representation before adding a directed hierarchy from Tree ID, Branching Order, and the
Branching Order 0 diving-trunk anchor.

Run Step 1-2 once:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\main_process_nne2.py `
  --stage preprocess --subject 020514 --tree-id 8
```

Then reuse its timestamped result for Step 3:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_3D_vascular\main_process_nne2.py `
  --stage hierarchical-graph `
  --source-run ulm_3D_vascular\outputs\NNE2_HDbase_v1.0\preprocess_run_TIMESTAMP `
  --subject 020514 --tree-id 8
```

Use `--stage all` to run both parts together. Every stack is segmented only once per run, and
reusable Step 1-2 and branch-graph caches are stored under `outputs/NNE2_HDbase_v1.0/_cache`.
Records missing required source fields remain listed in `skipped_missing_records.csv` and are not
processed. A Step 3 source run is accepted only after its manifest paths, sizes, and SHA-256
checksums have been verified.

The preprocessing run contains candidate, cleaned, removed-island, radius, and centerline volumes,
plus `component_decisions.csv`. The Step 3 run contains JSON, GraphML, NPZ, VTP, branch/node/cycle/
junction tables, directed parent-child files, root-component and excluded-island centerline volumes,
visual review images, logs, and an `acceptance_report.html`. Direction is an anatomical inference
away from the BO0 root, not a measured blood-flow direction. Ambiguous anchors and cross-links are
preserved and marked rather than silently forced into the primary tree.
