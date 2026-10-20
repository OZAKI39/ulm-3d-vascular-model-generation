# Latest Y-clot discretization: one Z-axis camera orbit

This replaces the three static discretization viewpoints with a single 18 s,
30 fps, 360° camera orbit. The physical state stays fixed at t = 0. This is a
visualization of the latest accepted model's initial configuration, not a new
fracture calculation or a time-dependent loading movie.

## Actual scientific inputs

Read-only source: `verification/y_clot_current_method_v2/20261019T131047Z/`.
The renderer checks the initial coordinates against `provenance/reference.npz`,
checks zero damage, and uses the original atom IDs, reference edges and
`BIPHASIC_FIELD.npz` material labels. Input SHA-256 values are checked again after
rendering. No coordinates, material assignments or connections are generated.

- 1927 clot points: 1045 weak, 882 strong; 1561 fixed anchors.
- Original irregular engineered wall-attached Y-vessel clot, not patient imaging.
- Cartesian sampling with Δx = 50 µm and horizon δ = 150.5 µm.
- Custom variational LPS (`bpm/peri/variational`), equal elastic moduli and
  heterogeneous softening thresholds. Strong/weak colors are material labels,
  not a map of elastic stiffness, damage, loading or preassigned fragments.
- PD quadrature volume 0.240875 mm³; it is distinct from the continuous-shape
  volume of approximately 0.271577 mm³.
- All 82,270 bulk and 22,040 interface reference edges are faintly shown.
  The 5,030 bulk and 290 interface edges at nearest-neighbor distance are
  emphasized for readability. The 53,327 internal fixed-scaffold edges are
  omitted from the picture only; the actual graph is unchanged.

## Display and outputs

Deep/light red points, gray-blue fixed anchors, complete six-face coordinate
grid, orthographic projection, and Z vertical in every frame. The camera
orbits a line parallel to global Z through the clot center, at −16° elevation.
The model remains in its actual global coordinates. Grid lines coincide with
the displayed coordinate ticks. Axis labels fade out and in when changing
silhouette edges. Only coordinates and one shared legend appear in the video.
There are no detached nodes or force arrows in this unloaded initial state.

Results: `results/clot_discretization_current_method_rotation_v1/20261019T205403Z/`.

- `VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4`: 3840×2160 master.
- `VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4`: 1920×1080 PPT copy.
- Both H.264, yuv420p, 540 frames, 30 fps, 18 s, no audio, fast-start MP4.
- `VIEW_*.png`: eight full-resolution viewpoints from the video frames.
- `ROTATION_CONTACT_SHEET.jpg`: eight-view overview.
- `OPEN_RESULTS.html`: local static review index.

Verification is under the matching `verification/` task/run directory:
540-frame projection/label checks, camera path, unchanged input hashes,
complete video decoding, media probes and a final delivery manifest.
The initial preview preceded a grid-only alignment refinement; the final
`VISUALIZATION_VALIDATION.json` identifies the delivered renderer hash.

## Reproduction

Use the existing isolated server environment, `source /workspace/y_clot/env.sh`.
Run `python analysis/clot_discretization_current_method_rotation_v1/render.py
--preview` for static previews. The supplied supervisor configuration and shell
wrapper run the complete export. The renderer refuses to overwrite a finalized
4K video: choose a new run directory before making another final version.
No native LAMMPS library or scientific solver settings are modified.
