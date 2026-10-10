# Aggregation and PZT translation: data-grounded visual stories

Created 2026-10-08. This is an offline visualization delivery plus one explicitly identified, matched Secondary-OFF computation. Existing F3/F4 archives are preserved.

## Deliverables

| Story | Figure | Movie |
|---|---|---|
| Attraction ON / OFF | [Figure 1](Fig_1_attraction_on_off.png) · [PDF](Fig_1_attraction_on_off.pdf) · [SVG](Fig_1_attraction_on_off.svg) | [Video 1](Video_1_attraction_on_off_60fps.mp4) |
| Translation and limits of shape control | [Figure 2](Fig_2_PZT_translation_shape_limit.png) · [PDF](Fig_2_PZT_translation_shape_limit.pdf) · [SVG](Fig_2_PZT_translation_shape_limit.svg) | [Video 2](Video_2_PZT_translation_shape_limit_60fps.mp4) |

Figures are 4800 × 2550 pixels, 300 dpi, with a white background and English text. PDF/SVG contain raster 3D panels and vector text/plots. Videos are 1920 × 1080, H.264, constant 60 fps, 601 frames each. Each is about 10.02 s long. Times displayed in the videos are **simulation time**. Video 1 samples every 1 ms (16.67-fold uniform slow motion); Video 2 samples every 0.2 ms (83.33-fold uniform slow motion). There are no variable-rate introductions, held opening frames, invented intermediate positions, or trajectory smoothing.

## Figure 1 / Video 1: attraction enables aggregation in the matched experiment

The original ON trajectory is `results/fast_pzt/f4/aggregation_v1/parent_wall_on`. Its 32-bubble dispersed cloud lies near the parent-vessel wall. The new OFF trajectory is [data/secondary_off/trajectory.h5](data/secondary_off/trajectory.h5). It uses exactly the same initial positions, four-element PZT drive, time step (10 µs), saved time grid (0.2 ms), radius (1.78 µm), mobility model, pair contact and geometric WALL constraints. Only the addition of the Secondary Bjerknes pair force is disabled. The PZT sources and Primary force remain enabled in both arms.

The OFF run executed on the existing RTX 4090 server under supervisor program `mevo_visual_secondary_off_20261008_v1`, using the frozen F4 engine `e247fa0c61f668a9b0b2e26c909dd687eedde389`. The isolated runner changes the Secondary toggle after the inherited constructor's frozen-model check and overrides only `pairs()`. It does not edit the inherited engine. A pre-run assertion verifies that the Primary force is nonzero and unchanged by the disabled pair-force stage. Runtime identity, engine hashes, configuration and numerical-quality records are preserved beside the HDF5 file. The numerical loop took 18.11 s and completed the requested 600 ms window.

| Quantity | Initial | Attraction ON, 600 ms | Attraction OFF, 600 ms |
|---|---:|---:|---:|
| Largest cluster | 1 / 32 | 32 / 32 | 1 / 32 |
| Radius of gyration | 12.8555 µm | 4.9003 µm | 12.8391 µm |

The ON run first reaches a 32-bubble connected cluster at saved time 101.2 ms. The OFF run remains 32 separate bubbles throughout this 600 ms window. This is a result for the specified initial state, drive and observation window, not a universal statement about all OFF configurations.

Clusters use the inherited spatial connected-component definition: an edge exists when the surface gap is ≤ 0.1 R (centre separation ≤ 2.1 R). Cluster membership does not mean adhesion or bubble fusion. The red local surface is a clipped portion of the actual F4 WALL mesh.

## Figure 2 / Video 2: changed centroid motion, no established robust shape control

These are the existing F3 continuations `repeats_v1/stable_reduction_0` and `controls_v1/amplitude`. Both begin at exactly the same mature 32-bubble state, taken at 102.2 ms of its precursor trajectory. Displayed time starts at zero for the 120 ms continuation. These F3 cases have zero flow and **no wall forces**; the global Y vessel is geometric context, not an assertion that this trajectory was recomputed using the F4 wall solver.

| Setting | Y+ | Y− | Z+ | Z− |
|---|---:|---:|---:|---:|
| Unchanged amplitude (Pa) | 151.3405 | 151.3405 | 151.3405 | 151.3405 |
| Controlled amplitude (Pa) | OFF | 246.2579 | OFF | 246.2579 |

All phases in these two settings are zero. The global panel shows the **controlled** setting. Disabled elements remain visible as faded hardware. Baseline hardware has the same positions, with all four elements enabled.

At 120 ms, the control-minus-baseline centroid distance is **1.3390 µm**. It is a contrast between two continuations, not the controlled cluster's absolute travel distance. Endpoint radius-of-gyration contrast is 1.5935%; eigenvalue contrast is 3.9235%. The inherited comparison classifies the result as `TRANSIENT_SHAPE_RESPONSE` and `NO_ROBUST_ENDPOINT_SHAPE_RESPONSE`. Its longest above-comparison shape-response span is 8.8 ms, below the predefined 10 ms persistence requirement. The inherited model and initial-state comparisons are descriptive sensitivity comparisons, not statistical confidence intervals or physical calibration.

Consequently, the supported wording is **“PZT shifts the cluster; robust shape control is not established.”** It would be inaccurate to say that shape never changes. The actual internal rearrangements, rotation and transient shape changes are retained. Crosses show the current centroid and the common initial centroid using the actual camera projection; their displacement is not magnified separately.

## Geometry, style and display fidelity

The four actual source centres are (0, ±4, 0) mm and (0, 0, ±4) mm, with 6 mm diameter and inward normals. Material colours, transparency, outward reflections and view orientation follow `PZT_control_reference_v6/Fig1a_PZT_control.png`. The 0.35 mm display backing and thin rims are illustrative housing geometry, not measured hardware dimensions. Disc colours do not encode pressure magnitude.

The global red vessel uses the original F4 `inputs/labelled_lumen_m.vtp`, WALL label 1, with all three terminal openings retained. The F4 mesh is used for consistency with the new comparison; this is not the older reference picture's separate mesh tessellation. Transparent surfaces use depth peeling and ordinary geometric occlusion. The global screen-space cross locates the common initial cloud centroid; it is neither a bubble nor a to-scale outline of the magnified region.

Each local view draws all 32 actual particle IDs as individual spheres of true radius 1.78 µm. The fixed local origin is the common initial centroid, and axes are in micrometres. Both arms of each comparison share one fixed camera and scale; the camera does not chase or recenter the cloud. No cluster is replaced by a proxy sphere, ideal lattice or rigid template. Figures use a shared camera sized for their three snapshots; videos use a shared camera sized for both complete histories. Therefore figure and video framing can differ slightly, while data and physical relative sizes are unchanged.

These remain **exploratory model results**. In particular, the inherited near-contact closure, local mobility and contact projection are not independent validation of physical aggregation. The inherited ON case's contact-correction/physical-path ratio is about 0.920. Its numerical quality gate passing and the images being smooth do not establish local physical accuracy.

## Reproduction and verification

Run from `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles` using `/home/lzy/projects/computation_examples/_shared/python_environment/bin/python`:

```bash
python -m scripts.stage4.visual_mechanism_stories
python -m scripts.stage4.visual_mechanism_stories --video aggregation --reuse-evidence
python -m scripts.stage4.visual_mechanism_stories --video translation --reuse-evidence
python -m scripts.stage4.validate_visual_mechanism_stories
```

Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LP_NUM_THREADS=4` for these offline renders. The rendering modules do not start a scientific simulation. The new OFF run has a separate explicit runner, `scripts/stage4/visual_evidence_ablation.py`, whose source and execution identity are retained. Frozen script copies in `source/` are archival copies; the module commands above use the canonical project paths.

`EVIDENCE.json` records actual input paths, HDF5 hashes, configurations and the inherited shape assessment. `data/*.npz` contains compact copies of the exact saved positions, times, radii, IDs and plotted metrics. `*_FRAME_MAP.npz` maps every video frame to its original saved sample. `*_IDENTITY.json` and `RENDER_SOURCE_*.json` record camera, time, controls and rendering source identities. `VALIDATION.json` records full-source metric recomputation, matched configuration checks, archive preservation and complete video decoding with uniform presentation timestamps. `ARTIFACT_MANIFEST.json` records delivery hashes.

## Suggested English captions

**Figure 1.** Matched Secondary-force ablation in a PZT-driven microbubble cloud near the parent-vessel wall. Identical dispersed initial conditions evolve into a 32-bubble connected cluster with the Secondary force enabled, whereas disabling only this force leaves 32 separate bubbles over 600 ms. The PZT drive, Primary force, mobility and contact constraints are retained. Spheres depict individual particles at their actual saved positions and relative physical sizes. This is an exploratory numerical-model comparison.

**Figure 2.** PZT-dependent translation of an existing microbubble cluster and current limits of shape control. Two 120 ms continuations from the same mature state use unchanged or modified source amplitudes. The endpoint centroid contrast is 1.339 µm. Transient shape changes and internal rearrangements occur, but the inherited sensitivity and persistence criteria do not establish robust shape control. The vessel provides geometric context for these F3 trajectories, which do not include wall forces.
