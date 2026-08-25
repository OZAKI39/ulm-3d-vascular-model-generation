# CFD surface preparation

## Purpose and scope

`cfd_surface_prepare.py` creates a CFD-specific surface from the already validated Ultraliser
lumen and the frozen, successful CFD preprocessing records. It is a surface-only stage. It does
not rerun ROI generation, Ultraliser, the global 1D flow solver, volume meshing, or a 3D CFD
solver.

The source `lumen_surface_um.stl` is read-only. Its path is taken from the saved
`geometry_reference.json`, and its SHA-256 digest is checked before and after the run. All changes
are made to a new in-memory copy and written under `outputs/cfd_surface_prepare/`.

## Geometry method

The four boundaries are read directly from the existing PASS preprocessing run: one
`ASSUMED_INLET` and three `ASSUMED_OUTLET` boundaries. The saved `CUT_PORT` or `TRUE_TERMINAL`
origin remains attached to each boundary, but both origins receive the same geometric treatment.

For each boundary, the program restricts its search to a short cylindrical neighborhood around
the saved center and outward normal. If more than one unrelated surface patch is present, the run
stops with `LOCAL_PORT_CUT_AMBIGUOUS`. Within the unambiguous patch, every triangle crossing the
saved boundary plane is intersected exactly. Adjacent triangles reuse the same edge intersection
vertex, producing one simple closed cut loop. The outward rounded end is removed while the inward
ROI surface is retained.

The actual, non-circularized cut loop is copied and translated rigidly by the saved extension
length from `port_extension_plan.csv`. Corresponding vertices form a straight side wall, and a
robust planar polygon triangulation closes the distal loop. No taper, flare, contraction, global
voxelization, implicit reconstruction, or whole-surface remeshing is used.

The combined VTP stores per-triangle tags:

- `boundary_type_code`: `0` wall, `1` inlet, `2` outlet;
- `boundary_index`: `-1` for wall, otherwise the saved boundary index;
- `boundary_origin_code`: `0` none/wall, `1` `CUT_PORT`, `2` `TRUE_TERMINAL`;
- `port_id`: empty for wall and the original saved ID for distal-cap triangles.

The side walls remain wall triangles. Only the four flat distal caps are boundary patches.

## Boundary-condition correction

The inlet remains a parabolic `VOLUMETRIC_FLOW_RATE` boundary with exactly the saved flow rate.
No pressure condition is added to it.

For each outlet, the actual distal-cap area defines an equivalent radius. The saved viscosity and
extension length then define the Poiseuille resistance of the added numerical tube:

\[
R_{ext}=\frac{8\mu L_{ext}}{\pi r_{eq}^{4}},\qquad
\Delta P_{ext}=R_{ext}Q_{1D},\qquad
P_{solver}=P_{original}-\Delta P_{ext}.
\]

This is marked `NUMERICAL_ARTIFICIAL_EXTENSION_CORRECTION`. It is only compensation for the
artificial straight tube; it is not a physiological outlet model, a resistance boundary
condition, or a Windkessel model. A negative gauge pressure at the `TRUE_TERMINAL` is permitted
and does not mean a negative absolute pressure.

The original boundary-condition JSON is copied unchanged. The extended solver-plane values are
written to a separate JSON and CSV so the original global-to-ROI pressure and flow remain
traceable.

## Quality controls

The run is accepted only if all four local constructions and the combined surface pass:

- exactly one local cut loop per boundary and the expected one-inlet/three-outlet tags;
- finite positive cap area, saved extension length, straight axis, planar cap, and outward cap
  normal;
- one watertight connected component, zero boundary and nonmanifold edges, zero degenerate
  triangles, and zero detected self-intersections;
- no extension collision with the original vessel or another extension;
- dense original-surface samples outside all four local surgery zones remain within the configured
  maximum and P95 distance tolerances;
- the source STL digest is identical before and after the run;
- the meter STL is the micrometre STL scaled by exactly `1e-6`.

A successful automatic result is deliberately named
`CFD_SURFACE_PREPARE_PASS_PENDING_MANUAL_REVIEW`. It does not mean that the surface is ready for
volume meshing or a solver. The next stage is `MANUAL_CFD_SURFACE_REVIEW`.

## Configuration and command

The strict default configuration is `configs/cfd_surface_prepare.yaml`. Run from the project
directory:

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_surface_prepare.py
```

An alternative strict YAML can be supplied without changing the program:

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_surface_prepare.py configs\my_surface_prepare.yaml
```

## Output layout

Each run writes a new directory:

```text
outputs/cfd_surface_prepare/<run_id>/
  input/
    cfd_surface_prepare.yaml
    source_manifest.json
    original_surface_reference.json
    cfd_preprocess_reference.json
  geometry/
    cfd_surface_extended_um.stl
    cfd_surface_extended_um.vtp
    cfd_surface_extended_m.stl
    cfd_wall_um.stl
  boundaries/
    boundary_00_*.stl
    boundary_01_*.stl
    boundary_02_*.stl
    boundary_03_*.stl
    boundary_manifest.csv
  bc/
    boundary_conditions_original.json
    boundary_conditions_extended.json
    extension_pressure_correction.csv
  qc/
    original_surface_integrity.json
    local_cut_qc.json
    boundary_geometry_qc.json
    core_surface_preservation_qc.json
    surface_qc.json
    extension_collision_qc.json
    meter_scale_qc.json
    run_summary.json
  figures/
    original_vs_cfd_surface.png
    boundary_closeups.png
  report/
    cfd_surface_prepare_report.md
```

The complete micrometre STL and two figures are intended for immediate manual review. The VTP is
the tagged version for preserving solver patch identities, while the four boundary STLs and wall
STL allow later meshing tools to recognize separate patches. This stage stops after producing
those files.
