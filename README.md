# Vascular Surface Reconstruction

This project converts a saved mouse-brain vascular ROI into a watertight surface with the
published Ultraliser `ultraVessMorpho2Mesh` executable. Ultraliser is the only supported surface
reconstruction backend.

## Input and output

The input is one saved `ROIRecord` from the representative sampling pipeline. It contains source
centerline coordinates, source radii, graph edges, local/global mappings, and `CUT_PORT` metadata.
The source ROI and canonical SWC radii are never modified.

The executable input is a vascular H5 adapter. Its fourth `points` column is diameter:

```text
feed_radius_um = source_radius_um * 0.91
h5_diameter_um = 2 * feed_radius_um
Ultraliser internal radius = h5_diameter_um * 0.5
```

Each successful run writes:

```text
outputs/model_generate/<run_id>/
  input/
    source_swc_stl_model_generate.yaml
    roi_core.swc
    roi_core.h5
    radius_feed_mapping.csv
    cut_port_mapping.csv
    metadata.json
  geometry/
    lumen_surface_um.stl
    lumen_surface_um.vtp
    lumen_surface_m.stl
  qc/
    surface_qc.json
    radius_fidelity.json
    run_summary.json
  report/
    reconstruction_report.md
```

The STL/VTP micrometre files are intended for geometry inspection and QC. The `*_m.stl` file has
coordinates scaled by exactly `1e-6` for later CFD preparation.

## Formal YAML configuration and command

Run from the `ulm_3D_vascular` project directory:

```powershell
D:\anaconda3\envs\pmp\python.exe swc_stl_model_generate.py configs\swc_stl_model_generate.yaml
```

All paths, ROI selectors, Ultraliser settings, and QC thresholds are defined in
`configs/swc_stl_model_generate.yaml`. Exactly one of `selection.roi_anchor` and
`selection.roi_id` must be non-null. Anchor 3274 is the reference validation case, not a code
restriction. The source YAML is copied into every model-generation run before reconstruction.

The existing `Ultraliser/build-wsl` executable is reused. A failed official process raises
`ULTRALISER_EXECUTION_FAILED`; no alternative reconstruction path is attempted.

## Validation

The retained reference evidence is under
`outputs/reference_validation/roi003274_ultraliser_radius091/`. With radius scale 0.91, ROI003274
measured a median signed radius error of `+1.1558%` and P95 absolute radius error of `4.0114%`.
The surface was watertight, single-component, free of detected self-intersections, nonmanifold
edges, and degenerate triangles, and the junction remained visually smooth and continuous.

See `docs/ULTRALISER_PIPELINE.md` for the compact call graph. CFD port preparation, volume meshing,
CFD solving, and microbubble simulation are intentionally outside the current stage.
